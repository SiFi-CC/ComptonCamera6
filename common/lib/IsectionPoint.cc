#include "IsectionPoint.hh"
#include <iostream>

using namespace std;

ClassImp(IsectionPoint);

//------------------------------------------------------------------
/// Default constructor.
IsectionPoint::IsectionPoint() {
  fGlobalBin = 0;
  fPoint = new TVector3(0, 0, 0);
}
//------------------------------------------------------------------
/// Standard constructor
///\param bin (Int_t) - global bin number.
///\param pos (TVector3) - coordinates of intersection point.
IsectionPoint::IsectionPoint(Int_t bin, TVector3* pos) {
  fGlobalBin = bin;
  fPoint = new TVector3(pos->X(), pos->Y(), pos->Z());
}
//------------------------------------------------------------------
/// Standard constructor
///\param bin (Int_t) - global bin number.
///\param x (Double_t) - x-component of intersection point.
///\param y (Double_t) - y-component of intersection point.
///\param z (Double_t) - z-component of intersection point.
IsectionPoint::IsectionPoint(Int_t bin, Double_t x, Double_t y, Double_t z) {
  fGlobalBin = bin;
  fPoint = new TVector3(x, y, z);
}
//------------------------------------------------------------------
/// Deafault destructor
IsectionPoint::~IsectionPoint() { delete fPoint; }
//------------------------------------------------------------------
/// Sets coordinates of intersection point and global bin number.
///\param bin (Int_t) - global bin number.
///\param x (Double_t) - x-component of intersection point.
///\param y (Double_t) - y-component of intersection point.
///\param z (Double_t) - z-component of intersection point.
void IsectionPoint::SetBinPoint(Int_t bin, Double_t x, Double_t y, Double_t z) {
  fGlobalBin = bin;
  fPoint->SetXYZ(x, y, z);
}
//------------------------------------------------------------------
/// Sets coordinates of intersection point.
///\param x (Double_t) - x-component of intersection point.
///\param y (Double_t) - y-component of intersection point.
///\param z (Double_t) - z-component of intersection point.
void IsectionPoint::SetPointCoordinates(Double_t x, Double_t y, Double_t z) {
  fPoint->SetXYZ(x, y, z);
}
//------------------------------------------------------------------
/// Sets global bin number.
///\param bin (Int_t) - global bin number.
void IsectionPoint::SetBin(Int_t b) { fGlobalBin = b; }
//------------------------------------------------------------------
/// Returns coordinates of intersection point.
TVector3* IsectionPoint::GetPointCoordinates(void) { return fPoint; }
//------------------------------------------------------------------
/// Returns global bin number.
Int_t IsectionPoint::GetBin(void) const { return fGlobalBin; }
//------------------------------------------------------------------
/// Compares global bin numbers to sort them.
Int_t IsectionPoint::Compare(const TObject* run2o) const {
  if (!run2o) {
    cout << "IsectionPoint::operator> : you attempt to compare to a "
            "non-existing object!"
         << endl;
  }
  const IsectionPoint* run2 = static_cast<const IsectionPoint*>(run2o);
  if (this->GetBin() < run2->GetBin())
    return -1;
  else if (this->GetBin() > run2->GetBin())
    return 1;
  else
    return 0;
}
//------------------------------------------------------------------
/// Prints details of the IsectionPoint class object.
void IsectionPoint::Print(void) {
  cout << "\nIsectionPoint::Print()" << endl;
  cout << "\tGlobal bin number: \t" << fGlobalBin << endl;
  cout << "\tPoint position: " << fPoint->X() << ", " << fPoint->Y() << ", "
       << fPoint->Z() << endl;
}
//------------------------------------------------------------------
