#include "DetPlane.hh"
#include <iostream>

using namespace std;

ClassImp(DetPlane);

//------------------------------------------------------------------
/// Default constructor.
DetPlane::DetPlane() {
  SetPlane(0, 0, 1, 0);
  SetDimensions(10, 10);
  SetName("plane");
}
//------------------------------------------------------------------
/// Standard constructor.
///\param a (Double_t) coefficient A of the cartesian plane equation
///\param b (Double_t) coefficient B
///\param c (Double_t) coefficient C
///\param d (Double_t) coefficient D
///\param dimZ (Double_t) full length along Z axis [mm]
///\param dimY (Double_t) full length along Y axis [mm]
///\param name (TString) object name
DetPlane::DetPlane(Double_t a, Double_t b, Double_t c, Double_t d,
                   Double_t dimZ, Double_t dimY, TString name) {

  SetPlane(a, b, c, d);
  SetDimensions(dimZ, dimY);
  SetName(name);
}
//------------------------------------------------------------------
/// Deafult destructor.
DetPlane::~DetPlane() {}
//------------------------------------------------------------------
/// Sets coefficients of the cartesian plane equation.
///\param a (Double_t) coefficient A of the cartesian plane equation
///\param b (Double_t) coefficient B
///\param c (Double_t) coefficient C
///\param d (Double_t) coefficient D
void DetPlane::SetPlane(Double_t a, Double_t b, Double_t c, Double_t d) {

  fA = a;
  fB = b;
  fC = c;
  fD = d;

  if (fA == 0 && fB == 0 && fC == 0) {
    cout << "##### Parameters incorrect! A, B and C cannot be equal 0 at the "
            "same time!"
         << endl;
    cout << "##### Please correct!" << endl;
    return;
  }
}
//------------------------------------------------------------------
/// Sets size of the detector plane.
///\param dimZ (Double_t) full length along Z axis [mm]
///\param dimY (Double_t) full length along Y axis [mm]
void DetPlane::SetDimensions(Double_t dimZ, Double_t dimY) {
  fDimZ = dimZ;
  fDimY = dimY;
}
//------------------------------------------------------------------
/// Returns coordinates of versor perpendicular to the given detector plane.
TVector3 DetPlane::GetNormal(void) {

  TVector3 norm;

  Double_t x = fA;
  Double_t y = fB;
  Double_t z = fC;
  Double_t n = sqrt(fA * fA + fB * fB + fC * fC);
  norm.SetXYZ(x / n, y / n, z / n);

  return norm;
}
//------------------------------------------------------------------
/// Checks whether given point belongs to the plane or not. Returns kTURE if
/// it does, and kFALSE if it doesn't.
///\param point (TVector3) - coordinates of the point.
Bool_t DetPlane::CheckPoint(TVector3 point) {

  if (fabs(point.Z()) > 0.5 * fDimZ) return kFALSE;
  if (fabs(point.Y()) > 0.5 * fDimY) return kFALSE;

  Double_t num = fabs(fA * point.X() + fB * point.Y() + fC * point.Z() + fD);
  Double_t denom = sqrt(fA * fA + fB * fB + fC * fC);
  Double_t dist = num / denom;

  if (dist != 0) return kFALSE;

  return kTRUE;
}
//------------------------------------------------------------------
/// Prints details of the DetPlane class object.
void DetPlane::Print(void) {
  cout << "\nDetPlane::Print() for object " << GetName() << endl;
  cout << "\tA = " << GetA() << endl;
  cout << "\tB = " << GetB() << endl;
  cout << "\tC = " << GetC() << endl;
  cout << "\tD = " << GetD() << endl;
  cout << "\tZ full size = " << GetDimZ() << endl;
  cout << "\tY full size = " << GetDimY() << endl;
  cout << "\tDimensions of the detector plane: " << fDimZ << "mm x " << fDimY
       << "mm" << endl;
  cout << "\tCartesian equation of the plane: " << fA << "x + " << fB << "y + "
       << fC << "z + " << fD << " = 0" << endl;
}
//------------------------------------------------------------------
