 
#include "CCMLEM.hh"
#include "IsectionPoint.hh"
#include <iostream>
using namespace std;

int main(int argc, char *argv[]){
  
  //expected arguments: startevent, stopevent

  TString tmps;
  Int_t nstart, nstop;
  Int_t verbose;
  cout<<"argc = "<<argc<<endl;
  if(argc<3){
    cout<<"I will analyze events 0-1000"<<endl;
    nstart = 0; 
    nstop = 1000;
  }
  else if(argc>3){
    tmps = argv[3];
    verbose = tmps.Atoi();
  }
  if(argc>2){
    tmps = argv[1];
    nstart = tmps.Atoi();
    tmps = argv[2];
    nstop = tmps.Atoi();
  }
  
  cout<<"I will analyze events "<<nstart<<"-"<<nstop<<endl;
     
  Int_t nev = 10;
  Int_t gen = 4;
  
  CCMLEM *rec = new CCMLEM(Form("../sources/results/CCSimulation_gen%i.root",gen),
			   Form("CCMLEM_gen%i",gen),1,verbose,80,80,80,80);
 
  rec->Reconstruct(nstart,nstop);
  
  
  
  delete rec;
 
  return 1;
}
  
  
