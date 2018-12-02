#include "PhysicsBase.hh"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <iostream>

using namespace std;

const double kR0 = 2.817E-15; // m
const double kMe = 0.510999;  // MeV/c2

//------------------------------------------------------------------
/// Klein-Nishina formula. Returns probablility of Compton scattering in m2.
///\param x - table of x, where x[0] is theta [deg]
///\param par - table of parameters, where par[0] is energy of incident photon
///[MeV].
Double_t KleinNishina(Double_t* x, Double_t* par) {

  Double_t costheta = cos(x[0] * TMath::DegToRad());
  Double_t alpha = par[0] / kMe;
  Double_t factor1 = (1. + costheta * costheta) / 2.;
  Double_t factor2 = 1. / (1. + alpha * (1. - costheta));
  Double_t factor3 =
      1 + (alpha * alpha * (1 - costheta) * (1 - costheta)) /
              ((1 + alpha * (1 - costheta)) * (1 + costheta * costheta));

  Double_t prob = kR0 * kR0 * factor1 * factor2 * factor2 * factor3; // m2
  // prob = prob*1E31;						// mb

  return prob;
}
//------------------------------------------------------------------

ClassImp(PhysicsBase);

//------------------------------------------------------------------
/// Deafault constructor
PhysicsBase::PhysicsBase() {
  SetName("physics");
  fTheta = 0.;
  fPhi = 0.;
}
//------------------------------------------------------------------
/// Standard constructor.
///\param name - object name.
PhysicsBase::PhysicsBase(TString name) {
  SetName(name);
  fTheta = 0.;
  fPhi = 0.;
}
//------------------------------------------------------------------
/// Deafault destructor.
PhysicsBase::~PhysicsBase() {
  if (fFunction) delete fFunction;
}
//------------------------------------------------------------------
/// Finds value of phi angle [rad].
Double_t PhysicsBase::FindPhi(void) {
  Double_t phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi()); // rad
  return phi;
}
//------------------------------------------------------------------
/// Generates Klein-Nishina function for requested energy of gamma quantum
/// and finds value of theta scattering angle [rad].
///\param energy (Double_t) energy of incident gamma quantum [MeV].
Double_t PhysicsBase::FindTheta(Double_t energy) {
  fFunction = new TF1("fFunction", KleinNishina, 0, 180, 1);
  fFunction->SetParameter(0, energy);
  Double_t tmp = fFunction->GetRandom();    // deg
  Double_t theta = tmp * TMath::DegToRad(); // rad
  return theta;
}
//------------------------------------------------------------------
/// Returns energy [MeV] of gamma quantum after Compton scattering.
///\param theta (Double_t) theta scattering angle (must be given in radians)
///\param initE (Double_t) initial energy of gamma quantum [MeV]
Double_t PhysicsBase::NewEnergy(Double_t theta, Double_t initE) {
  Double_t costheta = cos(theta); // theta must be in rad
  Double_t alpha = initE / kMe;
  Double_t finE = initE / (1 + alpha * (1 - costheta)); // MeV
  return finE;
}
//------------------------------------------------------------------
/// Method for Compton scattering process.
///\param initTrack (*Track) - track representing incident gamma quantum
///\param plane (*DetPlane) - plane of the detector
///
/// 1. Checks whether incident Track hits the detector plane. If not returns
/// NULL.
/// 2. Finds theta and phi scattering angles.
/// 3. Calculates gamma energy after scattering.
/// 4. Calculates coordinates of the new Track, representing scattered gamma
/// quantum.
/// 5. Assignes values to the scattered Track object and returns it.
Track* PhysicsBase::ComptonScatter(Track* initTrack, DetPlane* plane) {
  TVector3 finVersor;
  Double_t initE, finE;
  Double_t epsilon = 1.E-8;

  auto [crossPoint, crossFlag] = plane->FindCrossPoint(*initTrack);
  if (crossFlag == kFALSE) return NULL;

  initE = initTrack->GetEnergy();
  fTheta = FindTheta(initE); // rad
  // fTheta = TMath::Pi()/10.;
  fPhi = FindPhi(); // rad

  if (fabs(fTheta) < epsilon)
    cout << "##### Warning! Theta angle after scattering still equals 0!"
         << endl;
  if (fabs(fPhi) < epsilon)
    cout << "##### Warning! Phi angle after scattering still equals 0!" << endl;

  finE = NewEnergy(fTheta, initE);

  // finVersor.SetXYZ(-1,0,0);
  //----- scattering
  TVector3 xPrim, yPrim, zPrim;
  TVector3 xVersor, yVersor;
  TVector3 xComp, yComp, zComp;
  xVersor.SetXYZ(1, 0, 0);
  yVersor.SetXYZ(0, 1, 0);
  zPrim = initTrack->GetVersor();
  yPrim = zPrim.Cross(yVersor);
  xPrim = yPrim.Cross(zPrim);
  yPrim.SetMag(1.);
  xPrim.SetMag(1.);
  xComp = cos(fPhi) * sqrt(1 - pow(cos(fTheta), 2)) * xPrim;
  yComp = sin(fPhi) * sqrt(1 - pow(cos(fTheta), 2)) * yPrim;
  zComp = cos(fTheta) * zPrim;
  finVersor = xComp + yComp + zComp;
  //-----

  Track* finTrack = new Track();
  finTrack->SetPoint(crossPoint);
  finTrack->SetEnergy(finE);
  finTrack->SetVersor(finVersor);

  //----- theta check
  Double_t ang = initTrack->GetVersor().Angle(finVersor);
  if (fabs(ang - fTheta) > epsilon) {
    cout << "##### Error in PhysicsBase::ComptonScatter!" << endl;
    cout << "##### Incorrect theta angle! PLease check!" << endl;
    cout << "Chosen theta = " << fTheta * TMath::RadToDeg() << " deg \t "
         << "Set theta = " << ang * TMath::RadToDeg() << " deg" << endl;
    return NULL;
  }
  //----- end of the theta check

  return finTrack;
}
//------------------------------------------------------------------
/// Prints details of the PhysicsBase class object.
void PhysicsBase::Print(void) {
  cout << "\nPhysicsBase::Print() for object " << GetName() << endl;
  cout << "\tTheta scattering angle: " << fTheta * TMath::RadToDeg() << " deg"
       << endl;
  cout << "\tPhi scattering angle: " << fPhi * TMath::RadToDeg() << " deg"
       << endl;
}
//------------------------------------------------------------------
