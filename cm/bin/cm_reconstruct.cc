#include "CLog.hh"
#include "CMReconstruction.hh"
#include "G4SimulationAdapter.hh"
#include <TStopwatch.h>
#include "CmdLineConfig.hh"

int main(int argc, char** argv) {
  spdlog::set_level(spdlog::level::info);


  CmdLineOption opt_hmatrix("Hmatrix", "-hmat",
                           "File with H matrix, default: Calculate","");

  CmdLineConfig::instance()->ReadCmdLine(argc, argv);
  PositionalArgs args = CmdLineOption::GetPositionalArguments();


  if (args.size() != 2) {
    spdlog::info(
        "type: './cm_reconstruct [FILENAME] [ITERATIONS]' to start:\n\n"
        "where:\n\n"
        "FILE - is an input file from simulations\n\n"
        "ITERATIONS - is the numer of iterations to be processed.\n\n"); 
    return 1;
  }

  if(opt_hmatrix.GetStringValue()){ 
    spdlog::info("Hmatrix file: {}",opt_hmatrix.GetStringValue());
  } else {
    spdlog::info("Hmatrix will be calculated");
  }
  // return 1;

  TString filename(args[0]);
  Int_t iterations = TString(args[1]).Atoi();

  CMReconstruction reconstruction(filename);
  reconstruction.RunReconstruction(iterations);
  reconstruction.Write(filename.ReplaceAll(".root", "_reconstruct.root"));

  spdlog::info("Finished simulation");

  return 0;
}
