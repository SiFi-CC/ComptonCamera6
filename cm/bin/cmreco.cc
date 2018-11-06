#include "CMReconstruction.hh"
#include "TStopwatch.h"
#include <iostream>
using namespace std;

TStopwatch t; 

int main(int argc, char **argv){
  t.Start();
  if(argc!=3){
    cout<<"type: './cmreco fSimuFilename fNIter' to start:\n\n where:\n\n fSimuFilename is input filefrom simulations stored in the ../results/ directory\n\n fNIter is the numer of iterations to be processed.\n\n"<<endl;
    return 0;
  }
  TString inputfile(argv[1]);
  TString niterstring(argv[2]);
  Int_t niter = niterstring.Atoi();
  CMReconstruction *reco = new CMReconstruction(inputfile,1);
  reco->MLEMIterate(niter);
  reco->Write();
  t.Stop();
  t.Print();
  
  //delete reco;
  return 1;
}
