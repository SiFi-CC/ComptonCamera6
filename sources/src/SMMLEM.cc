#include "SMMLEM.hh"
#include <iostream>

using namespace std;

ClassImp(SMMLEM);

//------------------------------------------------------------------
SMMLEM::SMMLEM(){
 feventno= 0;
 fGlobalBin = 0;
 fdist = 0;
}
//------------------------------------------------------------------
SMMLEM::SMMLEM(Int_t eventno, Int_t bin, Double_t dist){
  feventno=eventno;
  fGlobalBin=bin;
  fdist=dist;
}
//------------------------------------------------------------------
SMMLEM::~SMMLEM(){
  
}
//------------------------------------------------------------------
void SMMLEM::SetEvBinDist(Int_t eventno, Int_t bin, Double_t dist){
  feventno=eventno;
  fGlobalBin=bin;
  fdist=dist;
}
//------------------------------------------------------------------
void SMMLEM::SetEvBin(Int_t eventno, Int_t bin){
  feventno=eventno;
  fGlobalBin=bin;
  
}
//------------------------------------------------------------------
void SMMLEM::SetEvent(Int_t i){
  feventno=i;
}
//------------------------------------------------------------------
void SMMLEM::SetBin(Int_t b){
  fGlobalBin=b;
}
//------------------------------------------------------------------
void SMMLEM::SetDist(Double_t d){
  fdist=d;
}
//------------------------------------------------------------------
Int_t SMMLEM::GetEvent(void) const{
  return feventno;
}
//------------------------------------------------------------------
Int_t SMMLEM::GetBin(void) const{
  
  return fGlobalBin;
}
//------------------------------------------------------------------
Double_t SMMLEM::GetDist(void){
  
  return fdist;
}
//------------------------------------------------------------------
void SMMLEM::Print(void){
 cout << "\nSMMLEM::Print()" <<endl;
 cout << "\tEvent number: \t"<<feventno<<endl;
 cout << "\tGlobal bin number: \t"<<fGlobalBin<<endl;
 cout << "\tLength of track: "<<fdist<< endl;
}
//------------------------------------------------------------------

