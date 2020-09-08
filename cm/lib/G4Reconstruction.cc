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

  //OSEM
  // Ns = CmdLineOption::GetIntValue("Subset");
  // T = sim.detector.nBins()/Ns;
  // nums = new int*[Ns];
  // length = new int[Ns];
  // product = new Double_t[T];
  // for (int i = 0; i < Ns; i++)
  // {
  //   nums[i]=new int[T];
  // }
  
  // for (int i = 0; i < Ns; i++)
  // {
  //   length[i]=0;
  // }
  
  // int randomNumber;
  // int i = 0;
  // while (i < sim.detector.nBins())
  // {
  //   randomNumber = rand() % Ns;
  //   // log->info("random {}", randomNumber);
  //   if(length[randomNumber]<=T-1){
  //     nums[randomNumber][length[randomNumber]]=i;
  //     length[randomNumber]++;
  //   } else{
  //     i--;
  //   }
  //   i++;
  // }
  //OSEM



  fMatrixH.ResizeTo(sim.detector.nBins(), sim.source.nBins());
  fMatrixHTranspose.ResizeTo(sim.source.nBins(), sim.detector.nBins());
  fRecoObject.push_back(TMatrixT<double>(sim.source.nBins(), 1));
  // fRecoObject[0] = 1.0 / sim.source.nBins();
  fRecoObject[0] = 1.0;
  fImage.ResizeTo(sim.detector.nBins(), 1);

  fImage = SiFi::tools::vectorizeMatrix(
      SiFi::tools::convertHistogramToMatrix(detector));

  // if (CmdLineOption::GetFlagValue("Hmatrix")){
    fMatrixH = sim.fMatrixHCam;
    fMatrixHTranspose.Transpose(fMatrixH);
    TMatrixT<Double_t>* psf = new TMatrixT<Double_t>(fMatrixH.GetNcols(),4);
    // TMatrixT<Double_t>* psf;
    TMatrixT<Double_t> fMatrixH2;
    TFile* file2 = new TFile("matr220_170_n15e5_nowallpet1cm_4lay_mask31_70mm_v2.root","READ");
    file2->cd();
    fMatrixH2.Read("matrixH");
    file2->Close();

    psf = GetPSF(fMatrixH);
    exit(0);
    // S.ResizeTo(fParams.source.nBins(), 1);
  // for (int i = 0; i < sim.detector.nBins(); i++)
  // {
  //   fImage(i,0) = fMatrixH(i,70 * 100 + 13)*1000;
  // }
  
    // for (int j = 0; j < fParams.source.nBins(); j++){
    //   S(j, 0) = 0.0;
    //   for (int i = 0; i < fParams.detector.nBins(); i++)
    //   {
    //     // fMatrixH(i,j) = fMatrixH(i,j)/1e7;
    //     S(j,0) += fMatrixH(i,j);
    //   }
    //     log->info("S({}) = {}",j,S(j,0));
    // }
  // fMatrixHTranspose.Transpose(fMatrixH);
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

  // OSEM
  // if(Ns > 1){ 
  //   TMatrixT<Double_t> nextSubset(fMatrixH.GetNcols(),1);
  //   TMatrixT<Double_t> prevSubset(fMatrixH.GetNcols(),1);
  //   prevSubset = fRecoObject.back();
  //   TMatrixT<Double_t> weightedImage(fParams.detector.nBins(), 1);
  //   for (int i = 0; i < Ns; i++)
  //   {
  //     auto hfProduct = fMatrixH * prevSubset;    
  //     for (int j = 0; j < fMatrixH.GetNcols(); j++)
  //     {
  //       sumH = 0;
  //       sum2 = 0;
  //       for (int t = 0; t < T; t++)
  //       {
  //         sum2 += fMatrixH(nums[i][t],j);
  //         sumH += fMatrixH(nums[i][t],j)*fImage(nums[i][t], 0)/hfProduct(nums[i][t],0);
  //       }
  //       nextSubset(j,0) = prevSubset(j,0)*sumH/sum2;
  //     }
  //     prevSubset = nextSubset;
  //   }
  //   // exit(0);
  //   fRecoObject.push_back(nextSubset);
  
  // }else{
  // OSEM
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
    weightedImage(i, 0) = fImage(i, 0) / hfProduct(i, 0);
  }
  log->debug("SingleIteration Image / (H * f_k) ({}, {})",
             weightedImage.GetNrows(), weightedImage.GetNcols());

  TMatrixT<Double_t> nextIteration = fMatrixHTranspose * weightedImage;
  for (int i = 0; i < fParams.source.nBins(); i++) {
    // nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0)/S(i,0);
    nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0);
  }
  // if(fRecoObject.size()>100){
  // if(0){
  // if(fRecoObject.size()>100 && fRecoObject.size()%20==0){
  //   TH2F recoIterationTMP = SiFi::tools::convertMatrixToHistogram(
  //         "reco", "iteration",
  //         SiFi::tools::unvectorizeMatrix(nextIteration, fParams.source.binY,
  //                                       fParams.source.binX),
  //         fParams.source.xRange, fParams.source.yRange);
  //   TH2F* smoothed = SmoothGauss(&recoIterationTMP, 0.4);
  //   TMatrixT<double> nextiterationsmoothed = SiFi::tools::vectorizeMatrix(
  //       SiFi::tools::convertHistogramToMatrix(smoothed));

  //   fRecoObject.push_back(nextiterationsmoothed);
  // } else {
  fRecoObject.push_back(nextIteration);
  // }
// }
  if(CmdLineOption::GetFlagValue("Autoiter") && fRecoObject.size() % 50 == 0 && fRecoObject.size() > 98){
    Double_t sigma = CheckConvergence(SiFi::tools::convertMatrixToHistogram(
          "reco", TString::Format("iteration %d", (int)fRecoObject.size()).Data(),
          SiFi::tools::unvectorizeMatrix(fRecoObject[fRecoObject.size()-1], fParams.source.binY,
                                         fParams.source.binX),
          fParams.source.xRange, fParams.source.yRange));
    log->info("iter = {}, sigma = {}, relative = {}",fRecoObject.size(),sigma,abs(sigma-fSigma)/fSigma);
    if (abs(sigma-fSigma)/fSigma > 0.1){
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

TH2F* G4Reconstruction::SmoothGauss(TH2F* hin, double sigma){

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


 TMatrixT<Double_t>* G4Reconstruction::GetPSF(TMatrixT<Double_t>  fMatrixH2){
  
  
  int detbins = fMatrixH2.GetNrows();
  int sourcebins = fMatrixH2.GetNcols();

  TH1D* hx;
  Double_t xmin, xmax, ymin, ymax, histmax;
  TF1* fGaus;
  TH2F recoIteration;

  TFile file("test.root", "RECREATE");
  file.cd();

  TMatrixT<Double_t>* psf = new  TMatrixT<Double_t>(sourcebins,4);
  // fImage.ResizeTo(detbins, 1);
  // psf->ResizeTo(sourcebins,4);

    // for (int j = 0; j < sourcebins; j++)
    for (int j = 5020; j < 5035; j++)
    {
      // int j = 5080;
      fRecoObject.clear();
      fRecoObject.push_back(TMatrixT<double>(sourcebins, 1));
      fRecoObject[0] = 1.0 / sourcebins;
      for (int i = 0; i < detbins; i++)
      {
        fImage(i,0)=fMatrixH2(i,j)*1000;          
      }
      int row = j%fParams.source.binY;
      int col = (j-row)/fParams.source.binY;
      RunReconstruction(100);
      recoIteration = SiFi::tools::convertMatrixToHistogram(
          "reco", TString::Format("bin   %d", j).Data(),
          SiFi::tools::unvectorizeMatrix(fRecoObject[fRecoObject.size()-1], fParams.source.binY,
                                        fParams.source.binX),
          fParams.source.xRange, fParams.source.yRange);
      double sx = recoIteration.GetXaxis()->GetBinCenter(col);
      double sy = -recoIteration.GetYaxis()->GetBinCenter(row);
      log->info("x = {}, y = {}",sx,sy);
      recoIteration.Write();
    
      hx = recoIteration.ProjectionX();
      // TH1F hy = recoIteration.ProjectionY();

      // exit(0);
      // double sy = 25;

      xmin = recoIteration.GetXaxis()->GetXmin();
      xmax = recoIteration.GetXaxis()->GetXmax();
      histmax = recoIteration.GetMaximum();
      // ymin = Image->GetYaxis()->GetXmin();
      // ymax = Image->GetYaxis()->GetXmax();
      
      // TF1* fSignal = new TF1("fSignal","[1]*exp(-0.5*((x-[2])/[3])^2)+[4]",sx-7,sx+7);
      fGaus = new TF1("fGaus","gaus",sx-7,sx+7);
      hx->SetTitle("ProjectionX");
      // fSignal->SetParameters(hx->GetMaximum(), sx,0.5,7);
      // fSignal->SetParameters(hx->GetMaximum(), sx,0.5,2);
      fGaus->SetParameters(histmax, sx,0.2,7);
      hx->Fit("fGaus", "", "",sx-5, sx+5);
      
      hx->Write();
      fRecoObject.clear();
    }

    file.Close();

    return psf;
}