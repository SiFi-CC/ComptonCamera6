#include "CMSimulation.hh"
#include <iostream>
using namespace std;

int main(int argc, char** argv) {

  if (argc != 3) {
    cout << "type: './cm_simulation fGenVersion fMask' to start:\n\n 1 - "
            "pointlike source in (0,0,0)\n 2 - uniform distribution along z, "
            "on beam axis, right-angled to it\n 3 - isotropic two point-like "
            "sources, 20mm gap\n 4 - isotropic two point-like sources, 40mm "
            "gap\n\n and fMask is a mask code, e.g. 1d_47 (for more options "
            "look into ../sources/macros/masks/\n\n"
         << endl;
    return 1;
  }
  TString mask(argv[2]);
  TString path = TString(gSystem->Getenv("CC6DIR"))+"/share/ComptonCamera6/masks/";
  TFile* maskfile = new TFile(path+"hMURA" + mask + ".root", "READ");
  TH2F* h = 0;
  TString fname = maskfile->GetName();
  if (fname.Contains("2d"))
    h = (TH2F*)maskfile->Get("hMURA2d");
  else
    h = (TH2F*)maskfile->Get("hMURA1d");

  TString gen(argv[1]);

  CMSimulation sim("CMSimulation_" + gen + "_MURA" + mask, 1);
  // maskdist,maskZsize,maskYsize,detdist,detZsize,detYsize, all mm
  sim.BuildSetup(500, 300, 300, 600, 300, 300);
  // 1 - point-like source at z=0, isotropic,
  // 2 - uniform along z,
  // 3 - two pointlike sources with 20mm gap,
  // 4 - two pointlike sources with 40mm gap
  sim.SetGenVersion(gen.Atoi());
  sim.SetPattern(h);
  sim.Print();
  sim.SetupSpectra();
  sim.Loop(1000);

  return 0;
}
