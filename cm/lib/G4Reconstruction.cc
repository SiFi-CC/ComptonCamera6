#include "G4Reconstruction.hh"
#include "DataStructConvert.hh"
#include <TTree.h>
#include <TVector.h>
#include "TF1.h"
#include "TF2.h"

#include "CmdLineConfig.hh"
using namespace std;

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
  if(1){
    S.ResizeTo(fParams.source.nBins(), 1);
    for (int j = 0; j < fParams.source.nBins(); j++){
      S(j, 0) = 0.0;
      for (int i = 0; i < fParams.detector.nBins(); i++)
      {
        fMatrixH(i,j) = fMatrixH(i,j)/10000.0;
        // fMatrixH(i,j) = fMatrixH(i,j)/(1.5*1e6);
        S(j,0) += fMatrixH(i,j);
      }
      // for (int i = 0; i < fParams.detector.nBins(); i++)
      // {
        // fMatrixH(i,j) = fMatrixH(i,j)/S(j,0);
      // }
    }
  }
    fMatrixHTranspose.Transpose(fMatrixH);
    ReadFit("Fit1D_matr220_170_n15e5_nowallpet1cm_4lay_mask31_70mm_noNormalized.txt");
    log->info("S1 = {}", S(1,0));
    S = ImageSpaceConvolute(S);
    
    // TMatrixT<Double_t> fMatrixH2;
    // TFile* file2 = new TFile("matr220_170_n15e5_nowallpet1cm_4lay_mask31_70mm_noNormalized.root","READ");
    // file2->cd();
    // fMatrixH2.Read("matrixH");
    // file2->Close();

    // GetPSF(fMatrixH2);
    // exit(0);
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

  auto hfProduct = fMatrixH * ImageSpaceConvolute(fRecoObject.back());
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

  TMatrixT<Double_t> nextIteration = ImageSpaceConvolute(fMatrixHTranspose * weightedImage);
  
  if(1){
    for (int i = 0; i < fParams.source.nBins(); i++) {
      if (fRecoObject.size() < 11)
      {
        nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0);
      } else {
        nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0)/S(i,0);
      }
    }
  } else {

    for (int i = 0; i < fParams.source.nBins(); i++) {
      // nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0)/S(i,0);
      nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0);
    }
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

#include <fstream>
 void G4Reconstruction::GetPSF(TMatrixT<Double_t>  fMatrixH2){
  
  
  int detbins = fMatrixH2.GetNrows();
  int sourcebins = fMatrixH2.GetNcols();

  TH1D* hx;
  TH1D* hy;
  Double_t xmin, xmax, ymin, ymax, histmax;
  TF1* fGaus;
  TH2F recoIteration;

  // TFile file("test.root", "RECREATE");
  // file.cd();

  	ofstream myfile;
	myfile.open("./data.txt");
	myfile<< left << setw(15) << "x"<< left << setw(15)<<"y"<< left << setw(15)
	<< left << setw(15) << "xFit"<< left << setw(15)<<"yFit"<< left << setw(15)
	<<"sigmax"<< left << setw(15) <<"sigmay"<< left << setw(15) <<"height"<<endl;
  Double_t param[5];

  // fImage.ResizeTo(detbins, 1);
  // psf->ResizeTo(sourcebins,4);

    // for (int j = 0; j < sourcebins; j++)
    for (int j =0; j < sourcebins; j++)
    {
      // int j = 5020;
      fRecoObject.clear();
      fRecoObject.push_back(TMatrixT<double>(sourcebins, 1));
      fRecoObject[0] = 1.0 / sourcebins;
      for (int i = 0; i < detbins; i++)
      {
        fImage(i,0)=fMatrixH2(i,j);          
      }
      int row = j%fParams.source.binY;
      int col = (j-row)/fParams.source.binY;
      RunReconstruction(100);
      recoIteration = SiFi::tools::convertMatrixToHistogram(
          "reco", TString::Format("bin   %d", j).Data(),
          SiFi::tools::unvectorizeMatrix(fRecoObject[fRecoObject.size()-1], fParams.source.binY,
                                        fParams.source.binX),
          fParams.source.xRange, fParams.source.yRange);
      double sx = recoIteration.GetXaxis()->GetBinCenter(col+1);
      double sy = -recoIteration.GetYaxis()->GetBinCenter(row+1);
      // log->info("x = {}, y = {}",sx,sy);
      // recoIteration.Write();
    
      
      TF2* fGaus2D = new TF2("fGaus2D","[0]*exp(-0.5*((x-[1])^2/[2]^2+(y-[3])^2/[4]^2))",sx-7,sx+7,sy-7,sy+7);
      // TF2* fGaus2D = new TF2("fGaus2D","[0]*exp(-0.5*((x-sx)^2/[1]^2+(y-sy)^2/[2]^2))",sx-7,sx+7,sy-7,sy+7);
      fGaus2D->SetParameters(recoIteration.GetMaximum(), sx,0.5,sy,0.5);
      // log->info("sx = {}, sy = {}, Maximum = {}",sx, sy,recoIteration.GetMaximum());

      recoIteration.Fit("fGaus2D");
      // recoIteration.Write();
      // recoIteration.ProjectionX()->Write();
      // recoIteration.ProjectionY()->Write();
     
      fGaus2D->GetParameters(param);
      
      // hy->Write();

      myfile << left << setw(15) << sx << left << setw(15)<< sy
      << left << setw(15) << param[1] << left << setw(15)<< param[3]
			<< left << setw(15) << abs(param[2])
			<< left << setw(15) << abs(param[4])
			<< left << setw(15) << param[0] << endl;

      fRecoObject.clear();
    }

    // file.Close();

}

void G4Reconstruction::ReadFit(TString  filename){

  Double_t temp1;
  Double_t temp2;
  Double_t temp3;
  Double_t temp4;
  Double_t temp5;

  sx.ResizeTo(fParams.source.nBins(), 1);
  sy.ResizeTo(fParams.source.nBins(), 1);
  sigmax.ResizeTo(fParams.source.nBins(), 1);
  sigmay.ResizeTo(fParams.source.nBins(), 1);
  histomax.ResizeTo(fParams.source.nBins(), 1);

  ifstream file;
  file.open(filename);

  // while (!file.eof())
  int n = 0;
  while(file >> sx(n,0) >> sy(n,0) >> sigmax(n,0) >> sigmay(n,0) >> histomax(n,0))
  {
    // sx(n,0)=temp1;
    // sy(n,0) = temp2;
    // sigmax(n,0) = temp3;
    // sigmay(n,0) = temp4;
    // histomax(n,0) = temp5;
    // log->info("{},  {}, {}, {}, {}",sx(n,0), sy(n,0),sigmax(n,0),sigmay(n,0),histomax(n,0));
    n++;
  }
  log->info("total number {}",n);
  file.close();
}

TMatrixT<Double_t> G4Reconstruction::ImageSpaceConvolute(TMatrixT<Double_t> image){

  // return image;

  int detbins = fMatrixH.GetNrows();
  int sourcebins = fMatrixH.GetNcols();

  int minrow, maxrow, mincol, maxcol, rad;
  int row, col, n;
  int nRows = fParams.source.binY;
  int nCols = fParams.source.binX;
  double k, maxEl;

  maxEl = 1;

  rad = 1;
  TMatrixT<Double_t> convoluted;

  convoluted.ResizeTo(sourcebins,1);

  for (int j = 0; j < sourcebins; j++)
  {
    convoluted(j,0) = 0.0;
    row = j%fParams.source.binY;
    col = (j-row)/fParams.source.binY;
    minrow = row - rad < 0 ? 0 : row - rad;
    maxrow = row + rad > nRows - 1 ? nRows - 1 : row + rad;
    mincol = col - rad < 0 ? 0 : col - rad;
    maxcol = col + rad > nCols - 1 ? nCols - 1 : col + rad;


    // if (row > 3 && row < 97 && col > 5 && col < 95){
      for (int irow = minrow; irow <= maxrow; irow++)
      {
        for (int icol = mincol; icol <= maxcol; icol++)
        {
          n = icol * nRows + irow;
          // k= histomax(j,0)*exp(-0.5*(pow(((sx(n,0)-sx(j,0))/sigmax(j,0)),2)+
          //                         pow(((sy(n,0)-sy(j,0))/sigmay(j,0)),2)));
          k= exp(-0.5*(pow(((sx(n,0)-sx(j,0))/sigmax(j,0)),2)+
                                  pow(((sy(n,0)-sy(j,0))/sigmay(j,0)),2)));
          convoluted(j,0) += image(n,0)*k;
        }
      }
      if (convoluted(j,0) > maxEl) maxEl = convoluted(j,0);
      if (convoluted(j,0) == 0) convoluted(j,0)=1e-9;
    // }
    // log->info("image({}) = {}, convoluted({}) = {}",j,image(j,0), j,convoluted(j,0));
    // log->info("image({}) = {}",j,image(j,0));
  }
    // for (int j = 0; j < sourcebins; j++)
    // {
    //   convoluted(j,0) /= maxEl;
    //   if (convoluted(j,0) == 0) convoluted(j,0) = 1;
    // }
    return convoluted;
}