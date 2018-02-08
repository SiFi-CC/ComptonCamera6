#include "CMReconstruction.hh"
#include "TStopwatch.h"
#include <iostream>
using namespace std;

TStopwatch t; 

int main(void){
  t.Start();
  CMReconstruction *reco = new CMReconstruction("CMsimulation_1.root",1);
  //  reco->FillHMatrix();
  reco->MLEMIterate(30);
  reco->Write();
  t.Stop();
  t.Print();
  
  //delete reco;
  return 1;
}
