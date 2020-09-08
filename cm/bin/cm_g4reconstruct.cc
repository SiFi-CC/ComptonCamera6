#include "CLog.hh"

#include <TStopwatch.h>

#include <CmdLineConfig.hh>

#include "CMReconstruction.hh"
#include "G4Reconstruction.hh"
#include "G4SimulationAdapter.hh"

int main(int argc, char** argv) {
  CmdLineOption cmdopt_output("Output", "-o",
                              "Output file (string), default: reconstruct.root",
                              "reconstruct.root");
  CmdLineOption cmdopt_iter("Iterations", "-n",
                            "Number of iterations, default: 20 (integer)", 20);
  // CmdLineOption cmdopt_hmat("Hmatrix", "-hmat",
  //                               "Data file provides H matrix, default: NO");
  CmdLineOption cmdopt_autoiter("Autoiter", "-autoiter",
                                "Dynamic number of iterations, but smaller then 'n', default: NO");
  // CmdLineOption cmdopt_subset("Subset", "-subset",
  //                               "Dynamic number of iterations, but smaller then 'n', default: NO",1);



  CmdLineArg cmdarg_simf("simfile", "Simulation file", CmdLineArg::kString);
  CmdLineArg cmdarg_dataf("datafile", "Data file", CmdLineArg::kString);

  CmdLineConfig::instance()->ReadCmdLine(argc, argv);

  spdlog::set_level(spdlog::level::debug);

  const Positional& pargs = CmdLineConfig::GetPositionalArguments();
  if (pargs.size() < 2) {
    spdlog::error("Not enough arguments, {} are required", 2);

    spdlog::info(
        "type: './cm_g4reconstruct [DETECTOR] [SIMULATIONS] [-n ITERATIONS] "
        "[-o OUTPUTFILE]' to "
        "start,\n"
        " where:\n"
        "  DETECTOR - is a result of single simulation, detector image from "
        "this file will be used to reconstruct source\n"
        "  SIMULATIONS - is a set of simulation results for difrent source "
        "position\n"
        "  ITERATIONS - is the numer of iterations to be processed, (default: "
        "20)\n"
        "  OUTPUTFILE - name of the output file (default: reconstruct.root\n");

    abort();
  }
  TString simFile(pargs.at("simfile")->GetStringValue());
  TString dataFile(pargs.at("datafile")->GetStringValue());

  TString reconstructFile = CmdLineOption::GetStringValue("Output");
  Int_t iterations = CmdLineOption::GetIntValue("Iterations");

  G4SimulationAdapter adapter(dataFile);

  // for now I'm assuming only one set of data is available in file
  CameraGeometry geometryData = adapter.GetFirstReconstructData();
  spdlog::info("Extracted {} inputs from data file",
               geometryData.recoData.size());
  geometryData.Print();

  // TODO: switch to using collecton of hits instead of histogram
  TFile simulationFile(simFile, "READ");
  if (!simulationFile.IsOpen()) {
    spdlog::error("Unable to open file {}", simFile.Data());
    return -1;
  }


  // TTree* fTree = (TTree*)simulationFile.Get("deposits");
  // TTree* fTreeSource = (TTree*)simulationFile.Get("source");
  // TVector3* fPosition = 0;
  // Double_t energy = 0;
  // Double_t energysource = 0;
  // Double_t sumenergy = 0;
  // Double_t sumenergysource = 0;
  // TH2F* weightedhisto = new TH2F("weighted","weighted histogram",
  //                               16,-8*1.3,8*1.3,
  //                               16,-8*1.3,8*1.3);
  
  // fTree->SetBranchAddress("position",&fPosition);
  // fTree->SetBranchAddress("energy",&energy);
  // fTreeSource->SetBranchAddress("energy",&energysource);
  // int size = fTree->GetEntries();
  // spdlog::info("entries = {}", size);

  // for (int i = 0; i < size; i++)
  // {
  //   fTree->GetEntry(i);
  //   // spdlog::info("x = {}, y = {}, energy = {}",fPosition->X(),fPosition->Y(),energy);
  //   // spdlog::info("x = {}",fPosition->X());
  //   weightedhisto->Fill(fPosition->X(),fPosition->Y(),energy);
  //   sumenergy += energy;
  //   fTreeSource->GetEntry(i);
  //   sumenergysource += energysource;
  // }
  // // TFile out("test.root", "RECREATE");
  // // out.cd();
  // // weightedhisto.Write();
  // // out.Close();
  // spdlog::info("totalenergyDeposited: {}", sumenergy);
  // spdlog::info("totalenergySource: {}", sumenergysource);
  // exit(0);








  auto detectorImage = static_cast<TH2F*>(simulationFile.Get("energyDeposits"));
  adapter.VerifyForReconstruct(&simulationFile);

  spdlog::info("Prepare reconstruction");
  G4Reconstruction reconstruction(geometryData, detectorImage);
  // G4Reconstruction reconstruction(geometryData, weightedhisto);
  spdlog::info("Run reconstruction");
  reconstruction.RunReconstruction(iterations);

  spdlog::info("Save reconstruction to file {}.", reconstructFile.Data());
  reconstruction.Write(reconstructFile);

  spdlog::info("Finished simulation");

  return 0;
}
