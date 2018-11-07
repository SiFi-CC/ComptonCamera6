#include "ComptonCone.hh"
#include <iostream>

using namespace std;

ClassImp(ComptonCone);

//------------------------------------------------------------------
/// Default constructor.
ComptonCone::ComptonCone() {
  TVector3 apex(0, 0, 0);
  TVector3 axis(-1, 0, 0);
  SetName("cone");
  SetApex(apex);
  SetAxis(axis);
  SetAngle(0);
}
//------------------------------------------------------------------
/// Standard constructor.
///\param name (TString) - name of the object
///\param apex (TVector3) - coordinates of the apex
///\param axis (TVector3) - cone axis
///\param angle (Double_t) - apex angle
ComptonCone::ComptonCone(TString name, TVector3 apex, TVector3 axis,
                         Double_t angle) {
  SetName(name);
  SetApex(apex);
  SetAxis(axis);
  SetAngle(angle);
}
//------------------------------------------------------------------
/// Standard constructor. Allows to build ComptonCone based on required
/// positions in the detector and registered energies.
///\param posScat (TVector3*) - position of Compton scattering, i.e. interaction
/// point in the scatterer
///\param posAbs (TVector3*) - position of absorption, i.e. interaction point in
/// the absorber
///\param enPrimary (Double_t) - energy of primary gamma [MeV]
///\param enScat (Double_t) - energy of gamma after Compton scattering [MeV]
ComptonCone::ComptonCone(TVector3* posScat, TVector3* posAbs,
                         Double_t enPrimary, Double_t enScat) {

  //----- Calculating angle theta
  const double me = 0.510999; // MeV/c2
  Double_t costheta = 1. - me * (1. / enScat - 1. / enPrimary);
  Double_t theta = TMath::ACos(costheta); // rad

  //----- Calculating cone axis
  TVector3 axis;
  axis.SetX(posScat->X() - posAbs->X());
  axis.SetY(posScat->Y() - posAbs->Y());
  axis.SetZ(posScat->Z() - posAbs->Z());
  axis.SetMag(1);

  //----- Setting values of private class members
  SetName("cone");
  SetAngle(theta);
  SetAxis(axis);
  SetApex(*posScat);
}
//------------------------------------------------------------------
/// Default destructor.
ComptonCone::~ComptonCone() {}
//------------------------------------------------------------------
/// Prints details of the ComptonCone class object.
void ComptonCone::Print(void) {
  cout << "\nComptonCone::Print() for object " << GetName() << endl;
  cout << "\tCoordinates of the apex: \n\t";
  fApex.Print();
  cout << "\tCoordinates of the leading versor of the axis: \n\t";
  fAxis.Print();
  cout << "\tAperture angle: " << fAngle << " radians\n" << endl;
}
//------------------------------------------------------------------
