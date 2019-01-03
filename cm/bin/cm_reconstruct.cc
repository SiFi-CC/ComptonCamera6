#include "CLog.hh"
#include "CMReconstruction.hh"
#include <TStopwatch.h>

namespace log = SiFi::log;

int main(int argc, char** argv) {
  log::setLevel(log::level::info);
  if (argc != 3) {
    log::info("type: './cm_reconstruct [FILENAME] [ITERATIONS]' to start:\n\n"
              "where:\n\n"
              "FILE - is an input file from simulations\n\n"
              "ITERATIONS - is the numer of iterations to be processed.\n\n");
    return 1;
  }

  TString filename(argv[1]);
  Int_t iterations = TString(argv[2]).Atoi();

  CMReconstruction reconstruction(filename);
  reconstruction.RunReconstruction(iterations);
  reconstruction.Write(filename.ReplaceAll(".root", "_reconstruct.root"));

  log::info("Finished simulation");

  return 0;
}
