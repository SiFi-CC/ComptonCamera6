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
  mask.zRange = std::make_pair(400, 420);
  mask.binX = 21;
  mask.binY = 21;
  mask.binZ = 1;

  detector.xRange = std::make_pair(-100, 100);
  detector.yRange = std::make_pair(-100, 100);
  detector.zRange = std::make_pair(550, 600);
  mask.binX = 21;
  mask.binY = 21;
  mask.binZ = 1; // TODO: very temporary

  // TODO temporary (this should be extracted from metadata
  //
  //   for (auto data : recoData) {
  //     log::debug("reading data from simulation {}", data->GetName())
  //   }
}

G4Reconstruction::G4Reconstruction(SimulationParams sim, TH2F* detector)
    : fParams(sim) {
  fMatrixH = TMatrixT<Double_t>(sim.source.nBins(), sim.detector.nBins());
  fMatrixHTranspose =
      TMatrixT<Double_t>(sim.detector.nBins(), sim.source.nBins());
  fImage = SiFi::tools::convertHistogramToMatrix(detector);
  for (auto file : sim.recoData) {
    TTree* srcTree = static_cast<TTree*>(file->Get("source"));
    TTree* detTree = static_cast<TTree*>(file->Get("deposits"));
    TVector3 srcPosition;
    srcTree->Branch("position", &srcPosition);
    TVector3 detPosition;
    detTree->Branch("position", &detPosition);

    if (srcTree->GetEntries() < 100) {
      log->error("too small number of events");
    }
    srcTree->GetEntry(0); // seting value on position variable

    auto sourceBin =
        sim.source.getBin(srcPosition.X(), srcPosition.Y(), srcPosition.Z());

    TMatrixT<Double_t> detMatrix(sim.detector.binY, sim.detector.binX);
    for (int i = 0; i < detTree->GetEntries(); i++) {
      detTree->GetEntry(i);

      auto detectorBin = sim.detector.getBin(detPosition.X(), detPosition.Y(),
                                             detPosition.Z());

      if (sim.detector.isValidBin(detectorBin)) {
        detMatrix(std::get<1>(detectorBin), std::get<0>(detectorBin)) += 1;
      }
    }
    auto column = SiFi::tools::vectorizeMatrix(detMatrix);
    double normFactor = 0;
    for (int row = 0; row < sim.detector.nBins(); row++) {
      normFactor += column(row, 0);
    }
    for (int row = 0; row < sim.detector.nBins(); row++) {
      int col =
          std::get<0>(sourceBin) * sim.source.binY + std::get<1>(sourceBin);
      fMatrixH(row, col) = column(row, 0) / normFactor;
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
    nextIteration(i, 0) =
        nextIteration(i, 0) * fRecoObject.back()(i, 0);
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

  int nIterations = fRecoObject.size() - 1;
  std::vector<TH2F*> histReco(nIterations);
  std::vector<TH1D*> histProjX(nIterations);
  std::vector<TH1D*> histProjY(nIterations);
  
  log->debug("Save {} iterations", nIterations);
  for (int i = 0; i < nIterations; i++) {
    log->debug("saving iteration {}", i);

    auto recoIteration = SiFi::tools::convertMatrixToHistogram(
        "reco", TString::Format("iteration %d", i).Data(),
        SiFi::tools::unvectorizeMatrix(fRecoObject[i], fParams.source.binX,
                          fParams.source.binY));
    recoIteration.Write();

  }

  file.Close();
  log->debug("end CMReconstruction::Write({})", filename.Data());
}
