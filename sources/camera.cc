#include "CCSimulation.hh"
#include <iostream>
using namespace std;

int main(void){

  Int_t nev = 100000;
  Int_t gen = 5;
  
  CCSimulation *sim = new CCSimulation(Form("CCSimulation_gen%i",gen),kFALSE);
  sim->SetGenVersion(gen);
  sim->BuildSetup(50,80,80,100,300,300);
  sim->Loop(nev);
  delete sim;
   
  return 1; 
}
