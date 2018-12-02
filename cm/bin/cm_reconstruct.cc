#include "CMReconstruction.hh"
#include "TStopwatch.h"
#include <iostream>
#include <spdlog/spdlog.h>

int main(int argc, char** argv) {
  spdlog::set_level(spdlog::level::debug);
  auto log = spdlog::stdout_logger_mt("cm_reconstruct_main", true);

  if (argc != 3) {
    log->info("type: './cm_reconstruct fSimuFilename fNIter' to start:\n\n"
              " where:\n\n"
              " fSimuFilename is input filefrom simulations stored in "
              "the ../results/ directory\n\n"
              " fNIter is the numer of iterations to be processed.\n\n");
    return 1;
  }

  TStopwatch t;
  t.Start();

  TString inputfile(argv[1]);
  TString niterstring(argv[2]);
  Int_t niter = niterstring.Atoi();

  CMReconstruction reco(inputfile, 1);
  reco.MLEMIterate(niter);
  log->info("Finished simulation");

  t.Stop();
  t.Print();

  return 0;
}
