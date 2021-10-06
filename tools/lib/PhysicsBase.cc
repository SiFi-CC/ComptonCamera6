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

/**
 * Klein-Nishina formula. Returns probablility of Compton scattering in m2.
 * \param x - table of x, where x[0] is theta [deg]
 * \param par - table of parameters, where par[0] is energy of incident photon [MeV].
 * \return probability of scattering into energy and angle phasespace element
 */
Double_t KleinNishina(Double_t* x, Double_t* par)
{
    Double_t costheta = cos(x[0] * TMath::DegToRad());
    Double_t alpha = par[0] / kMe;
    Double_t factor1 = (1. + costheta * costheta) / 2.;
    Double_t factor2 = 1. / (1. + alpha * (1. - costheta));
    Double_t factor3 = 1 + (alpha * alpha * (1 - costheta) * (1 - costheta)) /
                               ((1 + alpha * (1 - costheta)) * (1 + costheta * costheta));

    Double_t prob = kR0 * kR0 * factor1 * factor2 * factor2 * factor3; // m2
    // prob = prob*1E31;						// mb

    return prob;
}

/// Prints details of the PhysicsBase class object.
/*void PhysicsBase::Print(void) {
  cout << "\nPhysicsBase::Print() for object " << GetName() << endl;
  cout << "\tTheta scattering angle: " << fTheta * TMath::RadToDeg() << " deg"
       << endl;
  cout << "\tPhi scattering angle: " << fPhi * TMath::RadToDeg() << " deg"
       << endl;
}*/

namespace CC6
{

auto ComptonScatter(Double_t E, const TVector3& p_versor) -> std::tuple<Double_t, TVector3>
{
    Double_t epsilon = 1.E-8;

    auto theta = CC6::RandomKleinNishinaTheta(E); // rad
    // fTheta = TMath::Pi()/10.;
    auto phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi()); // rad

    if (fabs(theta) < epsilon)
        cout << "##### Warning! Theta angle after scattering still equals 0!" << endl;
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
    if (fabs(ang - theta) > epsilon)
    {
        cout << "##### Error in PhysicsBase::ComptonScatter!" << endl;
        cout << "##### Incorrect theta angle! PLease check!" << endl;
        cout << "Chosen theta = " << theta * TMath::RadToDeg() << " deg \t "
             << "Set theta = " << ang * TMath::RadToDeg() << " deg" << endl;
        return {};
    }
    //----- end of the theta check

    return {finE, finVersor};
}

Double_t RandomKleinNishinaTheta(Double_t energy)
{
    static auto f = std::make_unique<TF1>("KNFunction", KleinNishina, 0, 180, 1);
    f->SetParameter(0, energy);
    return f->GetRandom() * TMath::DegToRad(); // deg
}

Double_t ComptonScatteringGammaE(Double_t theta, Double_t initE)
{
    Double_t costheta = cos(theta); // theta must be in rad
    Double_t alpha = initE / kMe;
    Double_t finE = initE / (1 + alpha * (1 - costheta)); // MeV
    return finE;
}

} // namespace CC6
