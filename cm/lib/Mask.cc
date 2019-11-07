#include "Mask.hh"
#include <iostream>

using namespace std;

ClassImp(Mask);

//------------------------------------------------------------------
Mask::Mask() {
  SetPlane(0, 0, 1, 0);
  SetDimensions(10, 10);
  fPattern = NULL;
  SetName("mask");
}
//------------------------------------------------------------------
Mask::~Mask() {
  // if(fPattern) delete fPattern;
}
//------------------------------------------------------------------
Mask::Mask(Double_t a, Double_t b, Double_t c, Double_t d, Double_t dimZ,
           Double_t dimY, TH2F* h, TString name) {

  SetPlane(a, b, c, d);
  SetDimensions(dimZ, dimY);
  SetName(name);
  SetPattern(h);
}
//------------------------------------------------------------------
Int_t Mask::IsOpaque(TVector3 point) {
  // at the moment I assume fPattern->Xaxis() is along Z axis and
  // fPattern->Yaxis() is along y
  if (fPattern == 0) return 0;
  Int_t binx = fPattern->GetXaxis()->FindBin(point.Z());
  Int_t biny = fPattern->GetYaxis()->FindBin(point.Y());
  if (fPattern->GetBinContent(binx, biny) < 0.001)
    return 1;
  else
    return 0;
};
//------------------------------------------------------------------
void Mask::Print(void) {
  cout << "\nMask::Print() for object " << GetName() << endl;
  if (fPattern)
    cout << "\tPattern from histogram " << fPattern->GetName() << "\t"
         << fPattern << endl;
  else
    cout << "\tNo pattern defined yet..." << endl;
  ((DetPlane*)this)->Print();
}
//------------------------------------------------------------------
void Mask::SetPattern(TH2F* h) {
  fPattern = NULL;
  if (h) fPattern = (TH2F*)h->Clone("mask");
  fPattern->GetXaxis()->SetLimits(-GetDimZ() / 2, GetDimZ() / 2);
  fPattern->GetYaxis()->SetLimits(-GetDimY() / 2, GetDimY() / 2);
}
//------------------------------------------------------------------
