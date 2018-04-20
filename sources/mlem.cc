 
#include "CCMLEM.hh"
#include "IsectionPoint.hh"
#include <iostream>
using namespace std;

int main(void){
  
  Int_t nev = 10;
  Int_t gen = 4;
  
  CCMLEM *rec = new CCMLEM(Form("../sources/results/CCSimulation_gen%i.root",gen),
					       Form("CCMLEM_gen%i",gen),1,kFALSE,80,80);
 
  rec->Reconstruct(7,8);
  
  
  
  delete rec;
 
  return 1;
}
  
  