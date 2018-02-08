#include "CMSimulation.hh"
#include <iostream>
using namespace std;

int main(int argc, char **argv){
  
  if(argc!=2){
    cout<<"type: './mask fGenVersion' to start..."<<endl;
    return 0;
  }
  TFile* maskfile = new TFile("macros/masks/hMURA1d_47.root","READ");
  TH2F* h = (TH2F*)maskfile->Get("hMURA1d");
  
  TString gen(argv[1]);

  CMSimulation* sim = new CMSimulation("CMSimulation_"+gen, 1);
  //maskdist,maskZsize,maskYsize,detdist,detZsize,detYsize, all mm
  sim->BuildSetup(500,300,300,600,300,300);
  //1 - point-like source at z=0, isotropic, 
  //2 - uniform along z, 
  //3 - two pointlike sources with 20mm gap, 
  //4 - two pointlike sources with 40mm gap
  sim->SetGenVersion(gen.Atoi());
  sim->SetPattern(h);
  sim->Print();
  sim->SetupSpectra();
  sim->Loop(1000);
  delete sim;
  
  return 1; 
}
