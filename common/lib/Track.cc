#include "Track.hh"
#include "DetPlane.hh"
#include "TMath.h"
#include <iostream>

using namespace std;

ClassImp(Track);

//------------------------------------------------------------------
/// Standard constructor.
///\param point (TVector3) - starting point.
///\param vec (TVector3) - leading versor.
///\param energy (Double_t) - energy (MeV).
///\param name (TString) - object name.
Track::Track(TVector3 point, TVector3 vec, Double_t energy)
{
    SetPoint(point);
    SetVersor(vec);
    SetEnergy(energy);
}
//------------------------------------------------------------------
/// Sets coordinates of the starting point of the Track.
///\param point (TVector3) - starting point.
void Track::SetPoint(TVector3 point)
{
    fPoint.SetX(point.X());
    fPoint.SetY(point.Y());
    fPoint.SetZ(point.Z());
}
//------------------------------------------------------------------
/// Sets coordinates of the leading versor of the Track.
///\param vec (TVector3) - leading vector, later normalized into versor.
void Track::SetVersor(TVector3 vec)
{
    Double_t n = sqrt(vec.X() * vec.X() + vec.Y() * vec.Y() + vec.Z() * vec.Z());
    fVersor.SetX(vec.X() / n);
    fVersor.SetY(vec.Y() / n);
    fVersor.SetZ(vec.Z() / n);
}
//------------------------------------------------------------------
/// Sets Track energy.
///\param energy - energy in MeV
void Track::SetEnergy(Double_t energy) { fEnergy = energy; }

//------------------------------------------------------------------
/// Prints details of the Track class object.
void Track::Print() const
{
    cout << "\nTrack::Print() for object " << GetName() << endl;
    cout << "\tCoordinates of the starting point of the track: \n\t";
    fPoint.Print();
    cout << "\tCoordinates of the the track versor: \n\t";
    fVersor.Print();
    cout << "\tEnergy connected with this track: " << fEnergy << " MeV\n\n";
}
//------------------------------------------------------------------
