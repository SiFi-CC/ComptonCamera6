#include "PhysicsBase.hh"

#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"

#include <cmath>
#include <iostream>
#include <memory>

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

//------------------------------------------------------------------
/// Prints details of the PhysicsBase class object.
/*void PhysicsBase::Print(void) {
  cout << "\nPhysicsBase::Print() for object " << GetName() << endl;
  cout << "\tTheta scattering angle: " << fTheta * TMath::RadToDeg() << " deg"
       << endl;
  cout << "\tPhi scattering angle: " << fPhi * TMath::RadToDeg() << " deg"
       << endl;
}*/
//------------------------------------------------------------------

namespace CC6 {

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
auto ComptonScatter(Double_t E, const TVector3& p_versor,
                    const TVector3& cross_point)
    -> std::pair<Double_t, TVector3> {
  Double_t epsilon = 1.E-8;

  auto theta = CC6::RandomKleinNishinaTheta(E); // rad
  // fTheta = TMath::Pi()/10.;
  auto phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi()); // rad

  if (fabs(theta) < epsilon)
    cout << "##### Warning! Theta angle after scattering still equals 0!"
         << endl;
  if (fabs(phi) < epsilon)
    cout << "##### Warning! Phi angle after scattering still equals 0!" << endl;

  auto finE = CC6::ComptonScatteringGammaE(theta, E);

  // finVersor.SetXYZ(-1,0,0);
  //----- scattering
  TVector3 xVersor, yVersor;
  xVersor.SetXYZ(1, 0, 0);
  yVersor.SetXYZ(0, 1, 0);
  auto yPrim = p_versor.Cross(yVersor);
  auto xPrim = yPrim.Cross(p_versor);
  yPrim.SetMag(1.);
  xPrim.SetMag(1.);
  auto xComp = cos(phi) * sqrt(1 - pow(cos(theta), 2)) * xPrim;
  auto yComp = sin(phi) * sqrt(1 - pow(cos(theta), 2)) * yPrim;
  auto zComp = cos(theta) * p_versor;
  auto finVersor = xComp + yComp + zComp;
  //-----

  //----- theta check
  Double_t ang = p_versor.Angle(finVersor);
  if (fabs(ang - theta) > epsilon) {
    cout << "##### Error in PhysicsBase::ComptonScatter!" << endl;
    cout << "##### Incorrect theta angle! PLease check!" << endl;
    cout << "Chosen theta = " << theta * TMath::RadToDeg() << " deg \t "
         << "Set theta = " << ang * TMath::RadToDeg() << " deg" << endl;
    return {};
  }
  //----- end of the theta check

  return {finE, finVersor};
}

//------------------------------------------------------------------
/// Generates Klein-Nishina function for requested energy of gamma quantum
/// and finds value of theta scattering angle [rad].
///\param energy (Double_t) energy of incident gamma quantum [MeV].
Double_t RandomKleinNishinaTheta(Double_t energy) {
  auto f = std::make_unique<TF1>("KNFunction", KleinNishina, 0, 180, 1);
  f->SetParameter(0, energy);
  return f->GetRandom() * TMath::DegToRad(); // deg
}

/// Returns energy [MeV] of gamma quantum after Compton scattering.
///\param theta (Double_t) theta scattering angle (must be given in radians)
///\param initE (Double_t) initial energy of gamma quantum [MeV]
Double_t ComptonScatteringGammaE(Double_t theta, Double_t initE) {
  Double_t costheta = cos(theta); // theta must be in rad
  Double_t alpha = initE / kMe;
  Double_t finE = initE / (1 + alpha * (1 - costheta)); // MeV
  return finE;
}
//------------------------------------------------------------------

inline namespace literals {

long double operator"" _rad(long double deg) { return deg * M_PI / 180; }

long double operator"" _rad(unsigned long long deg) { return deg * M_PI / 180; }

} // namespace literals

} // namespace CC6
