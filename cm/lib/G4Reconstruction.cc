#include "G4Reconstruction.hh"
#include "DataStructConvert.hh"
#include "Smoothing.hh"
#include <TTree.h>
#include <TVector.h>
#include "TF1.h"
#include <TGraph.h>
#include <TCanvas.h>
#include <TMultiGraph.h>
#include <TLegend.h>

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
  fRecoObject[0] = 1.0;
  fImage.ResizeTo(sim.detector.nBins(), 1);

  fImage = SiFi::tools::vectorizeMatrix(
      SiFi::tools::convertHistogramToMatrix(detector));

  fMatrixH = sim.fMatrixHCam;
  if(1){ // Normalization of Hmatrix
    S.ResizeTo(fParams.source.nBins(), 1);
    for (int j = 0; j < fParams.source.nBins(); j++){
      S(j, 0) = 0.0;
      for (int i = 0; i < fParams.detector.nBins(); i++)
      {
        S(j,0) += fMatrixH(i,j);
      }
      for (int i = 0; i < fParams.detector.nBins(); i++)
      {
        fMatrixH(i,j) = fMatrixH(i,j)/S(j,0);
      }
    }
  }
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
    weightedImage(i, 0) = fImage(i, 0) / hfProduct(i, 0);
  }
  log->debug("SingleIteration Image / (H * f_k) ({}, {})",
             weightedImage.GetNrows(), weightedImage.GetNcols());

  TMatrixT<Double_t> nextIteration = fMatrixHTranspose * weightedImage;
  for (int i = 0; i < fParams.source.nBins(); i++) {
    nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0);
    // nextIteration(i, 0) = nextIteration(i, 0) * fRecoObject.back()(i, 0)/S(i,0);
  }

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

void G4Reconstruction::Write(TString filename, TH2F* simHist) const {
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
  std::vector <double> mse_list, uqi_list, iter_list;
  for (int i = 0; i < nIterations; i++) {
    log->debug("saving iteration {}", i);

    TH2F recoIteration = SiFi::tools::convertMatrixToHistogram(
        "reco", TString::Format("iteration %d", i).Data(),
        SiFi::tools::unvectorizeMatrix(fRecoObject[i], fParams.source.binY,
                                       fParams.source.binX),
        fParams.source.xRange, fParams.source.yRange);
    if((i+1)%50  == 0){ // Add smoothed histo each 50 iterations
      std::vector<Double_t> params = SiFi::tools::UQI_MSE(simHist,&recoIteration);
      recoIteration.SetTitle(Form("#splitline{MSE = %f}{UQI = %f}", params[0],params[1]));
      log->info("mse = {}", params[0]);
      log->info("uqi = {}", params[1]);
      mse_list.push_back(params[0]);
      uqi_list.push_back(params[1]);
      iter_list.push_back(i+1);
      TH2F* smoothed = SiFi::tools::SmoothGauss(&recoIteration, 1.5);
      smoothed->SetName(TString::Format("smoothed_%d",i+1));
      smoothed->Write();
    }
    recoIteration.Write();
  }
  TCanvas* c1 = new TCanvas("mas_uqi","mse_uqi",1);
  TGraph* gr = new TGraph(mse_list.size(),&iter_list[0], &mse_list[0]);
  TGraph* gr2 = new TGraph(mse_list.size(),&iter_list[0], &uqi_list[0]);
  TMultiGraph *mg = new TMultiGraph();
  
  gr->SetMarkerColor(4);
  gr->SetMarkerSize(1.5);
  gr->SetMarkerStyle(21);
  gr2->SetMarkerColor(3);
  gr2->SetMarkerSize(1.5);
  gr2->SetMarkerStyle(5);
  mg -> Add(gr);
  mg -> Add(gr2);
  mg -> Draw("ACP");
  TLegend *legend = new TLegend();
  legend->AddEntry(gr,"MSE","p");
  legend->AddEntry(gr2,"UQI","p");
  legend->Draw();
  c1->Write();

  file.WriteObject(&mse_list,"mse");
  file.WriteObject(&uqi_list,"uqi");
  file.WriteObject(&iter_list,"iter");

  TH2F histoS = SiFi::tools::convertMatrixToHistogram(
      "S", "Senesitivity map",
      SiFi::tools::unvectorizeMatrix(S, fParams.source.binY,
                                    fParams.source.binX),
      fParams.source.xRange, fParams.source.yRange);
  histoS.Write("SensitivityMap");

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


  return 0.5*(sigmaX+sigmaY);
}

