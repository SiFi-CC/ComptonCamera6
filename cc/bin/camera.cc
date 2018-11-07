#include "CCSimulation.hh"
#include <iostream>
using namespace std;

int main(void) {

  Int_t nev = 100000;
  Int_t gen = 1;

  CCSimulation* sim = new CCSimulation(Form("CCSimulation_gen%i", gen), kFALSE);
  sim->SetGenVersion(gen);
  sim->BuildSetup(230, 80, 80, 460, 300, 300);
  sim->Loop(nev);
  delete sim;

  return 1;
}
