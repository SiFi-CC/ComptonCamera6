#include "CLog.hh"
#include "CMSimulation.hh"
#include "CmdLineConfig.hh"
#include "Sources/MultiPointSource.hh"
#include "Sources/PlanarSource.hh"
#include "Sources/PointSource.hh"
#include <TSystem.h>
#include <iostream>

int main() {

  CmdLineOption _mask("mask", "-mask", "mask type, like 2d_29 or 1d_11",
                      "2d_23");

  spdlog::set_level(spdlog::level::info);
  TString maskFilename = CmdLineOption::GetStringValue("mask");
  TString path =
      TString(gSystem->Getenv("CC6DIR")) + "/share/ComptonCamera6/masks/";
  TString fullname = path + "hMURA" + maskFilename + ".root";
  TFile* maskfile = new TFile(fullname, "READ");
  if (maskfile == NULL) {
    std::cout << "File with mask: " << fullname << " does not exist, exit...\n"
              << std::endl;
    return 0;
  }
  TH2F* h = TString(maskfile->GetName()).Contains("2d")
                ? (TH2F*)maskfile->Get("hMURA2d")
                : (TH2F*)maskfile->Get("hMURA1d");

  // PointSource(TVector3(0,0,0), 4.4);
  // PointSource source("SinglePoint.mac");
  PlanarSource source("Planar.mac");
  // MultiPointSource source("Multipoint.mac");
  source.Print();

  // DetPlane detector(1, 0, 0, 600, 300, 300, "detector");
  // Mask mask(1, 0, 0, 500, 300, 300, h, "mask");
  DetPlane detector(1, 0, 0, 600, 300, 300, "detector");
  Mask mask(1, 0, 0, 500, 300, 300, h, "mask");

  CMSimulation sim(&source, &mask, &detector);
  sim.RunSimulation(10000);
  sim.Write("results/simulation" + maskFilename + ".root");

  return 1;
}
