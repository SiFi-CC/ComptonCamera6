#include "CCReconstruction.hh"
#include <iostream>
using namespace std;

int main(void){
 
  Int_t nev = 100000;
  Int_t gen = 5;
  
  CCReconstruction *rec;
  
  try{
    rec =  new CCReconstruction(Form("../sources/results/CCSimulation_gen%i.root",gen),
			        Form("CCReconstruction_gen%i",gen),kFALSE);
  }
  catch(const char *message){
    cout << message << endl;
    return 0;
  }
  
  rec->RebuildSetupTxt();
  rec->ReconstructImage(0,nev);
  delete rec;
  
  return 1;
}
