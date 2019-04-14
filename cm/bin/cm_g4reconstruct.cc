#include "CLog.hh"
#include "CMReconstruction.hh"
#include "G4Reconstruction.hh"
#include "G4SimulationAdapter.hh"
#include <TStopwatch.h>

int main(int argc, char** argv) {
  spdlog::set_level(spdlog::level::info);
  if (argc != 3) {
    spdlog::info(
        "type: './cm_g4reconstruct [FILENAME] [ITERATIONS]' to start:\n\n"
        "where:\n\n"
        "FILE - is an input file from simulations\n\n"
        "ITERATIONS - is the numer of iterations to be processed.\n\n");
    return 1;
  }
  TString filename(argv[1]);
  Int_t iterations = TString(argv[2]).Atoi();


  // TODO: add cmd arg
  G4SimulationAdapter adapter("../../g4sim/build/grid_simulation.root");
  // for now I'm assuming only one set of data is available in file
  auto input = adapter.GetFirstReconstructData();
  spdlog::info("Extracted {} inputs from data file", input.size());

  // TODO: switch to usng collecton of hits instead of histogram
  TFile simulationFile(filename, "READ");
  auto detectorImage = static_cast<TH2F*>(simulationFile.Get("energyDepositions"));

  SimulationParams params;
  params.recoData = input;
  params.initSimulationMetadata();

  G4Reconstruction reconstruction(params, detectorImage);
  reconstruction.RunReconstruction(iterations);
  reconstruction.Write(filename.ReplaceAll(".root", "_reconstruct.root"));

  spdlog::info("Finished simulation");

  return 0;
}
