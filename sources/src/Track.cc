#include "Track.hh"
#include "DetPlane.hh"
#include "TMath.h"
#include <iostream>

using namespace std;

ClassImp(Track);

//------------------------------------------------------------------
///Default constructor. 
Track::Track(){
  SetPoint(TVector3(0,0,0));
  SetVersor(TVector3(-1,0,0));
  SetEnergy(4.44);
  SetName("track");
}
//------------------------------------------------------------------
///Standard constructor.
///\param point (TVector3) - starting point.
///\param vec (TVector3) - leading versor.
///\param energy (Double_t) - energy (MeV).
///\param name (TString) - object name.
Track::Track(TVector3 point, TVector3 vec, Double_t energy, TString name){
  SetPoint(point);
  SetVersor(vec);
  SetEnergy(energy);
  SetName(name);
}
//------------------------------------------------------------------
/// Standard destructor.
Track::~Track(){
}
//------------------------------------------------------------------
///Sets coordinates of the starting point of the Track.
///\param point (TVector3) - starting point.
void Track::SetPoint(TVector3 point){ 
  fPoint.SetX(point.X());
  fPoint.SetY(point.Y());
  fPoint.SetZ(point.Z());
}
//------------------------------------------------------------------
///Sets coordinates of the leading versor of the track.
///\param vec (TVector3) - leading vector, later normalized into versor.
void Track::SetVersor(TVector3 vec){
  Double_t n = sqrt(vec.X()*vec.X() + 
                    vec.Y()*vec.Y() +
		    vec.Z()*vec.Z());
  fVersor.SetX(vec.X()/n);
  fVersor.SetY(vec.Y()/n);
  fVersor.SetZ(vec.Z()/n);
}
//------------------------------------------------------------------
///Sets Track's energy.
///\param energy - energy in MeV
void Track::SetEnergy(Double_t energy){
 fEnergy = energy;
}
//------------------------------------------------------------------
///Calculates coordinates of cross point of the Track and given plane. 
///Function returns kTRUE when there is single cross point and kFALSE in the following 
///cases: (1) when Track is parallel to the plane and there are no cross points, (2) 
///when plane includes the Track and all points are common, (3) when cross point is outside 
///of the plane (size of the plane is specified in DetPlane object).
///\param plane (DetPlane*) - plane wich intersects with the given Track.
///\param position (TVector3) - coordinates of the cross point returned as reference.
Bool_t Track::FindCrossPoint(DetPlane* plane, TVector3 &position){
  
  Double_t a    = plane->GetA();
  Double_t b    = plane->GetB();
  Double_t c    = plane->GetC();
  Double_t d    = plane->GetD();
  Double_t dimZ = plane->GetDimZ();
  Double_t dimY = plane->GetDimY();
  Double_t l    = fVersor.X();
  Double_t m    = fVersor.Y();
  Double_t n    = fVersor.Z();
  Double_t x    = fPoint.X();
  Double_t y    = fPoint.Y();
  Double_t z    = fPoint.Z();
    
  Double_t num   = a*x + b*y + c*z + d;
  Double_t denom = a*l + b*m + c*n;
  
  if(fabs(denom) < 1.E-8){
    cout << "##### The track is parallel to the plane! No cross points!" << endl;
    return kFALSE;
  }
  if(fabs(num) < 1.E-8){
    cout << "##### The plane includes the track! All points are common!" << endl;
    return kFALSE;
  }
  
  Double_t rho = num/denom;
  Double_t crossX = x - l*rho;
  Double_t crossY = y - m*rho;
  Double_t crossZ = z - n*rho;
  
  position.SetXYZ(crossX,crossY,crossZ);
  
  if(l>0 || fabs(crossZ)>(0.5*dimZ) || fabs(crossY)>(0.5*dimY)){
    //cout << "\n\tCrossing point outside of the detector plane..." << endl;
    return kFALSE;
  }
  
  return kTRUE;
}
//------------------------------------------------------------------
///Prints details of the Track class object.
void Track::Print(void){
  cout << "\nTrack::Print() for object " << GetName() << endl;
  cout << "\tCoordinates of the starting point of the track: \n\t";
  fPoint.Print();
  cout << "\tCoordinates of the the track versor: \n\t";
  fVersor.Print();
  cout << "\tEnergy connected with this track: " << fEnergy << " MeV\n\n";
}
//------------------------------------------------------------------