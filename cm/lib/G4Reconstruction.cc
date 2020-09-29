#include "G4Reconstruction.hh"
#include "DataStructConvert.hh"
#include <TTree.h>
#include <TVector.h>
#include "TF1.h"

#include "CmdLineConfig.hh"

G4Reconstruction::G4Reconstruction(CameraGeometry sim, TH2F* detector)
    : fParams(sim) {
  log->debug(
      "G4Reconstruction::G4Reconstruction(simulationParams, detectorImage)");
  if (detector->GetYaxis()->GetNbins() != sim.detector.binY ||
      detector->GetXaxis()->GetNbins() != sim.detector.binX) {
      log->error("Image bins does not match simulation data");
      throw "bin mismatch";
  }

  fMatrixH.ResizeTo(sim.detector.nBins(), sim.source.nBins());
  fMatrixHTranspose.ResizeTo(sim.source.nBins(), sim.detector.nBins());
  fRecoObject.push_back(TMatrixT<double>(sim.source.nBins(), 1));
  // fRecoObject[0] = 1.0 / sim.source.nBins();
  fImage.ResizeTo(sim.detector.nBins(), 1);

  fImage = SiFi::tools::vectorizeMatrix(
      SiFi::tools::convertHistogramToMatrix(detector));

  fMatrixH = sim.fMatrixHCam;
 if(1){
    S.ResizeTo(fParams.source.nBins(), 1);
    for (int j = 0; j < fParams.source.nBins(); j++){
      S(j, 0) = 0.0;
      for (int i = 0; i < fParams.detector.nBins(); i++)
      {
        // fMatrixH(i,j) = fMatrixH(i,j)/10000.0;
        // fMatrixH(i,j) = fMatrixH(i,j)/(1.5*1e6);
        S(j,0) += fMatrixH(i,j);
      }
      for (int i = 0; i < fParams.detector.nBins(); i++)
      {
        fMatrixH(i,j) = fMatrixH(i,j)/S(j,0);
      }
    }
  }
  fRecoObject[0] = S;
  //  TH2F histoS = SiFi::tools::convertMatrixToHistogram(
    //       "S", "Senesetivity map",
    //       SiFi::tools::unvectorizeMatrix(S, fParams.source.binY,
    //                                     fParams.source.binX),
    //       fParams.source.xRange, fParams.source.yRange);
    // TFile* fileS = new TFile("S.root","RECREATE");
    // fileS->cd();
    // histoS.Write();
    // fileS->Close();
    // exit(0); 
  fMatrixHTranspose.Transpose(fMatrixH);
}

TMatrixT<Double_t> G4Reconstruction::ReadFromTH2F(TH2F* detHist) {
  return SiFi::tools::vectorizeMatrix(
      SiFi::tools::convertHistogramToMatrix(detHist));
}

TMatrixT<Double_t> G4Reconstruction::ReadFromTTree(TBranch* detBranch) {
  TVector3* detPosition = new TVector3();
  detBranch->SetAddress(&detPosition);
  TMatrixT<Double_t> detMatrix(fParams.detector.binY, fParams.detector.binX);
  for (int i = 0; i < detBranch->GetEntries(); i++) {
    detBranch->GetEntry(i);

    auto detectorBin = fParams.detector.getBin(
        detPosition->X(), detPosition->Y(), detPosition->Z());

    if (fParams.detector.isValidBin(detectorBin)) {
      detMatrix(fParams.detector.binY - std::get<1>(detectorBin),
                std::get<0>(detectorBin) - 1) += 1;
    }
  }
  return SiFi::tools::vectorizeMatrix(detMatrix);
}

void G4Reconstruction::RunReconstruction(int nIter) {
  int goon = 1;
  int iter = 0;

  while (iter < nIter && goon == 1) {
    fIter = iter+1;
    goon = SingleIteration();
    iter++;
  }
}

int G4Reconstruction::SingleIteration() {
  log->debug("CMReconstruction::SingleIteration()  iter={}",
             fRecoObject.size());

  // H * f_k
  // vector of correlations
  // i-th element of this vector contains information whether image resulting
  // from source postioned in i-th bin correlates to current iteration.
  // log->debug("Hmatrix = {}, {}, FReco ",
  //              fMatrixH.GetNrows(),fMatrixH.GetNcols() );

  auto hfProduct = fMatrixH * fRecoObject.back();
  log->debug("SingleIteration  H * f_k ({}, {})", hfProduct.GetNrows(),
             hfProduct.GetNcols());

  // Image / (H * f_k)
  // Image vector with inverse of correlations used as weights, the largest
  // weight will be atributed for image element that has the worst corelation
  // with current iteration.
  TMatrixT<Double_t> weightedImage(fParams.detector.nBins(), 1);
  for (int i = 0; i < fParams.detector.nBins(); i++) {
    // if (fRecoObject.size() < 50){
    //   weightedImage(i, 0) = 1 / hfProduct(i, 0);
    // } else {
    //   weightedImage(i, 0) = fImage(i, 0) / hfProduct(i, 0);
    // }
    weightedImage(i, 0) = fImage(i, 0) / hfProduct(i, 0);
  }
  log->debug("SingleIteration Image / (H * f_k) ({}, {})",
             weightedImage.GetNrows(), weightedImage.GetNcols());

  TMatrixT<Double_t> nextIteration = fMatrixHTranspose * weightedImage;
  for (int i = 0; i < fParams.source.nBins(); i++) {
    nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0);
  }
  // if(1){
  //   for (int i = 0; i < fParams.source.nBins(); i++) {
  //     if (fRecoObject.size() < 11)
  //     {
  //       nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0);
  //     } else {
  //       nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0)/S(i,0);
  //     }
  //   }
  // } else {

  //   for (int i = 0; i < fParams.source.nBins(); i++) {
  //     // nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0)/S(i,0);
  //     nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0);
  //   }
  // }

  fRecoObject.push_back(nextIteration);

  if(CmdLineOption::GetFlagValue("Autoiter") && fRecoObject.size() % 25 == 0 && fRecoObject.size() > 98){
    Double_t sigma = CheckConvergence(SiFi::tools::convertMatrixToHistogram(
          "reco", TString::Format("iteration %d", (int)fRecoObject.size()).Data(),
          SiFi::tools::unvectorizeMatrix(fRecoObject[fRecoObject.size()-1], fParams.source.binY,
                                         fParams.source.binX),
          fParams.source.xRange, fParams.source.yRange));
    log->info("iter = {}, sigma = {}, relative = {}",fRecoObject.size(),sigma,abs(sigma-fSigma)/fSigma);
    if (abs(sigma-fSigma)/fSigma > 0.01){
      fSigma = sigma;
    } else {
      return 0;
    }
  }

  log->debug("end CMReconstruction::SingleIteration()  iter={}",
             fRecoObject.size() - 1);
  return 1;
}

void G4Reconstruction::Write(TString filename) const {
  log->info("CMReconstruction::Write({})", filename.Data());
  TFile file(filename, "RECREATE");
  file.cd();

  SiFi::tools::convertMatrixToHistogram(
      "histH", "histogram of matrix H(probability matrix)", fMatrixH)
      .Write();

  SiFi::tools::convertMatrixToHistogram(
      "image", "image on detector",
      SiFi::tools::unvectorizeMatrix(fImage, fParams.detector.binX, fParams.detector.binY))
      .Write();

  int nIterations = fRecoObject.size() - 1;

  log->debug("Save {} iterations", nIterations);
  for (int i = 0; i < nIterations; i++) {
    log->debug("saving iteration {}", i);

    TH2F recoIteration = SiFi::tools::convertMatrixToHistogram(
        "reco", TString::Format("iteration %d", i).Data(),
        SiFi::tools::unvectorizeMatrix(fRecoObject[i], fParams.source.binY,
                                       fParams.source.binX),
        fParams.source.xRange, fParams.source.yRange);
    if((i+1)%100  == 0){
    // if(i == 500){
      TH2F* smoothed = SmoothGauss(&recoIteration, 1.3);
      smoothed->Write();
    }
    recoIteration.Write();
  }
  // log->info("Sigma {}",fSigma);
  TVector sig(1);
  sig[0] = fSigma;
  sig.Write("sigma");

  TVector iter(1);
  iter[0] = fIter;
  iter.Write("maxIter");

  file.Close();

  log->debug("end CMReconstruction::Write({})", filename.Data());
}

//Returns the average sigma value of X and Y Projections of a given 2D histogram
Double_t G4Reconstruction::CheckConvergence(TH2F reco){
  TH1D *hx, *hy;
  TF1 *fSignal;
  Double_t sigmaX, sigmaY;


  Double_t xmin = reco.GetXaxis()->GetXmin();
  Double_t xmax = reco.GetXaxis()->GetXmax();
  Double_t ymin = reco.GetYaxis()->GetXmin();
  Double_t ymax = reco.GetYaxis()->GetXmax();

  hx = reco.ProjectionX();
  hx->SetTitle("ProjectionX");
  hy = reco.ProjectionY();
  hy->SetTitle("ProjectionY");

  Double_t sx = hx->GetXaxis()->GetBinCenter( hx->GetMaximumBin() );
  Double_t sy = hy->GetYaxis()->GetBinCenter( hy->GetMaximumBin() );

  fSignal = new TF1("fSignal","gaus",sx-5,sx+5  );

  fSignal->SetParameters(hx->GetMaximum(), sx,0.5,10);
  hx->Fit("fSignal", "Q", "",sx-5, sx+5);

  sigmaX = fSignal->GetParameter(2);


  fSignal = new TF1("fSignal","gaus",sy-5,sy+5);

  fSignal->SetParameters(hx->GetMaximum(), sy,0.5,10);
  hy->Fit("fSignal", "Q", "",sy-5, sy+5);

  sigmaY = fSignal->GetParameter(2);



  // TFile file("out2.root", "UPDATE");
  // file.cd();

  // hx->Write("hx");
  // hy->Write("hy");

  // file.Close();

  return 0.5*(sigmaX+sigmaY);
}


TH2F* G4Reconstruction::SmoothGauss(TH2F* hin, double sigma) const{

  if(sigma <= 0){
    std::cout << "Smearing with sigma = " << sigma 
              << " will not work, provide a positive number here..." << std::endl;
    return nullptr;
  }
  
  TH2F* hout = dynamic_cast<TH2F*>(hin->Clone(Form("%s_smooth", hin->GetName())));
  hout->Reset();
  const int nbinsx = hin->GetNbinsX();
  const int nbinsy = hin->GetNbinsY();
  double binwx = hin->GetXaxis()->GetBinWidth(1);
  double binwy = hin->GetYaxis()->GetBinWidth(1);
  
  double kernelx[nbinsx];
  double kernely[nbinsy];
  
  TF1* gaus = new TF1("gaus", "gaus");
  gaus->SetParameters(1./sqrt(TMath::TwoPi())/sigma, 0, sigma);

  for(int i=0; i<nbinsx; i++)
    kernelx[i] = gaus->Eval(binwx*i);
  for(int i=0; i<nbinsy; i++)
    kernely[i] = gaus->Eval(binwy*i);

  int deltabin = 0;
  double z = 0;

  //smearing in rows
  for(int biny=1; biny<nbinsy; biny++){
    for(int binx=1; binx<nbinsx; binx++){
      z = 0;
      for(int binxp=1; binxp<nbinsx; binxp++){
        deltabin = abs(binxp-binx);
        z += kernelx[deltabin]*hin->GetBinContent(binxp, biny);
      }
      hout->SetBinContent(binx, biny, z);      
    }
  }
  TH2F* htmp = dynamic_cast<TH2F*>(hout->Clone());
  hout->Reset();
  
  //smearing in columns
  for(int binx=1; binx<nbinsx; binx++){
    for(int biny=1; biny<nbinsy; biny++){
      z = 0;
      for(int binyp=1; binyp<nbinsy; binyp++){
        deltabin = abs(binyp-biny);
        z += kernely[deltabin]*htmp->GetBinContent(binx, binyp);        
      }
      hout->SetBinContent(binx, biny, z);      
    }
  }

  return hout;
}