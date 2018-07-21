 
#include "CCMLEM.hh"
#include "IsectionPoint.hh"
#include <iostream>
using namespace std;

int main(int argc, char *argv[]){
  
  //expected arguments: startevent, stopevent, verbose, niter

  TString tmps;
  
  Int_t nstart, nstop, niter, verbose;
  cout<<"argc = "<<argc<<endl;
  if(argc!=5){
    cout<<"Please run the program like this:\n ./mlem <nstart> <nstop> <verbose> <niter>"<<
	"\n where:\n\t nstart - number of first event to process (int),\n\t nstop - number of last event to process (int),\n\t verbose - verbose flag (1/0),\n\t niter - number of iterations in MLEM (int)"<< endl;
    return 0;
  }
  tmps = argv[3];
  verbose = tmps.Atoi();
  
  tmps = argv[1];
  nstart = tmps.Atoi();
  tmps = argv[2];
  nstop = tmps.Atoi();
  tmps = argv[4];
  niter = tmps.Atoi();
  
  cout<<"I will analyze events "<<nstart<<"-"<<nstop<<" with "<<
    niter<<" MLEM  iterations..."<<endl;
     
  
  
  CCMLEM *rec = new CCMLEM();
 
  rec->Reconstruct(nstart,nstop);
  
  
  
  delete rec; 
 
  return 1;
}