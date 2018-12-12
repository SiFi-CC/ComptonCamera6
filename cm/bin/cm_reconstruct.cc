#include "CMReconstruction.hh"
#include "TStopwatch.h"
#include <iostream>

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cout << "type: './cm_reconstruct fSimuFilename fNIter' to start:\n\n"
                 " where:\n\n"
                 " fSimuFilename is input filefrom simulations stored in "
                 "the ../results/ directory\n\n"
                 " fNIter is the numer of iterations to be processed."
              << std::endl;
    return 1;
  }

  TStopwatch t;
  t.Start();

  TString inputfile(argv[1]);
  TString niterstring(argv[2]);
  Int_t niter = niterstring.Atoi();

  CMReconstruction reco(inputfile, 1);
  reco.MLEMIterate(niter);
  std::cout << "Finished simulation" << std::endl;

  t.Stop();
  t.Print();

  return 0;
}
