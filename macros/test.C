#include <iostream>
#include <fstream>
#include "/home/kasia/ComptonCamera/include/ComptonCone.hh"

using namespace std;

int test(void){

 gSystem->Load("/home/kasia/ComptonCamera/libComptonCamera.so");

 TString name = "cone";
 TVector3 apex(1,1,1);
 TVector3 axis(3,3,3);
 Double_t angle = TMath::Pi()/180.;
 ComptonCone *cone = new ComptonCone(name,apex,axis,angle);

 cone->Print();
  
  return 1;
}

