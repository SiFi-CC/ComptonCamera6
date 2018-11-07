#include "SMElement.hh"
#include <iostream>

using namespace std;

ClassImp(SMElement);

//------------------------------------------------------------------
SMElement::SMElement() {
  feventno = 0;
  fGlobalBin = 0;
  fdist = 0;
}
//------------------------------------------------------------------
SMElement::SMElement(Int_t eventno, Int_t bin, Double_t dist) {
  feventno = eventno;
  fGlobalBin = bin;
  fdist = dist;
}
//------------------------------------------------------------------
SMElement::~SMElement() {}
//------------------------------------------------------------------
void SMElement::SetEvBinDist(Int_t eventno, Int_t bin, Double_t dist) {
  feventno = eventno;
  fGlobalBin = bin;
  fdist = dist;
}
//------------------------------------------------------------------
void SMElement::SetEvent(Int_t i) { feventno = i; }
//------------------------------------------------------------------
void SMElement::SetBin(Int_t b) { fGlobalBin = b; }
//------------------------------------------------------------------
void SMElement::SetDist(Double_t d) { fdist = d; }
//------------------------------------------------------------------
Int_t SMElement::GetEvent(void) const { return feventno; }
//------------------------------------------------------------------
Int_t SMElement::GetBin(void) const { return fGlobalBin; }
//------------------------------------------------------------------
Double_t SMElement::GetDist(void) { return fdist; }
//------------------------------------------------------------------
void SMElement::Print(void) {
  cout << "\nSMElement::Print()" << endl;
  cout << "\tEvent number: \t" << feventno << endl;
  cout << "\tGlobal bin number: \t" << fGlobalBin << endl;
  cout << "\tLength of track: " << fdist << endl;
}
//------------------------------------------------------------------
