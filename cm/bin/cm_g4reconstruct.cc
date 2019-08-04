#include "CLog.hh"
#include "CMReconstruction.hh"
#include "G4Reconstruction.hh"
#include "G4SimulationAdapter.hh"
#include <TStopwatch.h>

int main(int argc, char** argv) {
  spdlog::set_level(spdlog::level::debug);
  if (argc != 4) {
    spdlog::info(
        "type: './cm_g4reconstruct [DETECTOR] [SIMULATIONS] [ITERATIONS]' to "
        "start:\n\n"
        "where:\n\n"
        "DETECTOR - is a result of single simulation, detector image from "
        "this file will be used to reconstruct source\n\n"
        "SIMULATIONS - is a set of simulation results for difrent source "
        "position\n\n"
        "ITERATIONS - is the numer of iterations to be processed.\n\n");
    return 1;
  }

  TString simFile(argv[1]);
  TString dataFile(argv[2]);
  TString reconstructFile("reconstruct.root");
  Int_t iterations = TString(argv[3]).Atoi();

  G4SimulationAdapter adapter(dataFile);

  // for now I'm assuming only one set of data is available in file
  CameraGeometry geometryData = adapter.GetFirstReconstructData();
  spdlog::info("Extracted {} inputs from data file",
               geometryData.recoData.size());
  geometryData.Print();

  // TODO: switch to using collecton of hits instead of histogram
  TFile simulationFile(simFile, "READ");
  if (!simulationFile.IsOpen()) {
    spdlog::error("Unable to open file {}", simFile);
    return -1;
  }
  auto detectorImage = static_cast<TH2F*>(simulationFile.Get("energyDeposits"));
  adapter.VerifyForReconstruct(&simulationFile);

  spdlog::info("Prepare reconstruction");
  G4Reconstruction reconstruction(geometryData, detectorImage);
  spdlog::info("Run reconstruction");
  reconstruction.RunReconstruction(iterations);

  spdlog::info("Save reconstruction to file {}.", reconstructFile);
  reconstruction.Write(reconstructFile);

  spdlog::info("Finished simulation");

  return 0;
}
