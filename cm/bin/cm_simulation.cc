#include "CLog.hh"
#include "CMSimulation.hh"
#include "Sources/PointSource.hh"
#include <TSystem.h>
#include <iostream>

int main(int argc, char** argv) {
  spdlog::set_level(spdlog::level::info);
  TString maskFilename(argv[1]);
  TString path =
      TString(gSystem->Getenv("CC6DIR")) + "/share/ComptonCamera6/masks/";
  TFile* maskfile = new TFile(path + "hMURA" + maskFilename + ".root", "READ");
  TH2F* h = TString(maskfile->GetName()).Contains("2d")
                ? (TH2F*)maskfile->Get("hMURA2d")
                : (TH2F*)maskfile->Get("hMURA1d");

  MultiPointSource source(TVector3(0, 0, 0));
  source.AddSourceElement(PointSource(TVector3(0, 0, 0), 1));
  source.AddSourceElement(PointSource(TVector3(0, 0, 100), 1));
  source.AddSourceElement(PointSource(TVector3(0, 100, 0), 1));
  source.AddSourceElement(PointSource(TVector3(0, -100, -100), 1));

  DetPlane detector(1, 0, 0, 600, 300, 300, "detector");
  Mask mask(1, 0, 0, 500, 300, 300, h, "mask");

  CMSimulation sim(&source, &mask, &detector);
  sim.RunSimulation(100000);
  sim.Write("results/simuation" + maskFilename + ".root");

  return 0;
}
