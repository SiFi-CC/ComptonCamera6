#include "G4Reconstruction.hh"
#include "DataStructConvert.hh"
#include <TTree.h>

void SimulationParams::initSimulationMetadata() {
  source.xRange = std::make_pair(-105, 105);
  source.yRange = std::make_pair(-105, 105);
  source.zRange = std::make_pair(-0.1, 0.1);
  source.binX = 21;
  source.binY = 21;
  source.binZ = 1;

  mask.xRange = std::make_pair(-100, 100);
  mask.yRange = std::make_pair(-100, 100);
  mask.zRange = std::make_pair(485, 515);
  mask.binX = 21;
  mask.binY = 21;
  mask.binZ = 1;

  detector.xRange = std::make_pair(-150, 150);
  detector.yRange = std::make_pair(-150, 150);
  detector.zRange = std::make_pair(635, 685);
  detector.binX = 100;
  detector.binY = 100;
  detector.binZ = 1; // TODO: very temporary
  spdlog::info("setup hardcoded defaults for simulation options");

  // TODO temporary (this should be extracted from metadata inside adapter)
  //
  //   for (auto data : recoData) {
  //     log::debug("reading data from simulation {}", data->GetName())
  //   }
}

G4Reconstruction::G4Reconstruction(SimulationParams sim, TH2F* detector)
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
  fRecoObject[0] = 1.0 / sim.source.nBins();
  fImage.ResizeTo(sim.detector.nBins(), 1);

  fImage = SiFi::tools::vectorizeMatrix(
      SiFi::tools::convertHistogramToMatrix(detector));

  log->info("Building H matrix");
  for (auto file : sim.recoData) {
    auto srcBranch =
        static_cast<TTree*>(file->Get("source"))->GetBranch("position");
    auto detBranch =
        static_cast<TTree*>(file->Get("deposits"))->GetBranch("position");
    TVector3* srcPosition = new TVector3();
    srcBranch->SetAddress(&srcPosition);
    TVector3* detPosition = new TVector3();
    detBranch->SetAddress(&detPosition);

    srcBranch->GetEntry(0); // seting value on position variable
    // check if std dev is small enough

    auto sourceHistBin =
        sim.source.getBin(srcPosition->X(), srcPosition->Y(), srcPosition->Z());
    auto sourceMatBin =
        std::make_tuple<int, int>(sim.source.binY - std::get<1>(sourceHistBin),
                                  std::get<0>(sourceHistBin) - 1);
    int colIndexMatrixH =
        std::get<1>(sourceMatBin) * sim.source.binY + std::get<0>(sourceMatBin);

    log->debug("processing point source at {}, {} in histogram bin({}, {}) = "
               "matrix bin ({}, {})",
               srcPosition->X(), srcPosition->Y(), std::get<0>(sourceHistBin),
               std::get<1>(sourceHistBin), std::get<0>(sourceMatBin),
               std::get<1>(sourceMatBin));

    TMatrixT<Double_t> detMatrix(sim.detector.binY, sim.detector.binX);
    for (int i = 0; i < detBranch->GetEntries(); i++) {
      detBranch->GetEntry(i);

      auto detectorBin = sim.detector.getBin(detPosition->X(), detPosition->Y(),
                                             detPosition->Z());

      if (sim.detector.isValidBin(detectorBin)) {
        detMatrix(sim.detector.binY - std::get<1>(detectorBin),
                  std::get<0>(detectorBin) - 1) += 1;
      }
    }
    auto column = SiFi::tools::vectorizeMatrix(detMatrix);
    double normFactor = 0;
    for (int row = 0; row < sim.detector.nBins(); row++) {
      normFactor += column(row, 0);
    }
    for (int row = 0; row < sim.detector.nBins(); row++) {
      fMatrixH(row, colIndexMatrixH) =
          column(row, 0) == 0 ? 1e-9 : (column(row, 0) / normFactor);
    }
  }
  fMatrixHTranspose.Transpose(fMatrixH);
}

void G4Reconstruction::RunReconstruction(int nIter) {
  for (int iter = 0; iter < nIter; iter++) {
    SingleIteration();
  }
}

void G4Reconstruction::SingleIteration() {
  log->debug("CMReconstruction::SingleIteration()  iter={}",
             fRecoObject.size());

  // H * f_k
  // vector of correlations
  // i-th element of this vector contains information whether image resulting
  // from source postioned in i-th bin correlates to current iteration.
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
  }

  fRecoObject.push_back(nextIteration);

  log->debug("end CMReconstruction::SingleIteration()  iter={}",
             fRecoObject.size() - 1);
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
      SiFi::tools::unvectorizeMatrix(fImage, 100, 100))
      .Write();

  int nIterations = fRecoObject.size() - 1;

  log->debug("Save {} iterations", nIterations);
  for (int i = 0; i < nIterations; i++) {
    log->debug("saving iteration {}", i);

    auto recoIteration = SiFi::tools::convertMatrixToHistogram(
        "reco", TString::Format("iteration %d", i).Data(),
        SiFi::tools::unvectorizeMatrix(fRecoObject[i], fParams.source.binY,
                                       fParams.source.binX));
    recoIteration.Write();
  }

  file.Close();
  log->debug("end CMReconstruction::Write({})", filename.Data());
}
