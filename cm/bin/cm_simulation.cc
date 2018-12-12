#include "CMSimulation.hh"
#include "Sources/PointSource.hh"
#include <TSystem.h>
#include <iostream>

int main(int argc, char** argv) {
  TString maskFilename(argv[2]);
  TString path =
      TString(gSystem->Getenv("CC6DIR")) + "/share/ComptonCamera6/masks/";
  TFile* maskfile = new TFile(path + "hMURA" + maskFilename + ".root", "READ");
  TH2F* h = TString(maskfile->GetName()).Contains("2d")
                ? (TH2F*)maskfile->Get("hMURA2d")
                : (TH2F*)maskfile->Get("hMURA1d");

  PointSource source(TVector3(0, 0, 0), 1);
  DetPlane detector(1, 0, 0, 600, 300, 300, "detector");
  Mask mask(1, 0, 0, 500, 300, 300, h, "mask");

  CMSimulation sim(&source, &mask, &detector);
  sim.RunSimulation(10000);
  sim.Write("results/simuation" + maskFilename);

  return 0;
}
