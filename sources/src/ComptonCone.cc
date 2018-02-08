#include "ComptonCone.hh"
#include <iostream>

using namespace std;

ClassImp(ComptonCone);

//------------------------------------------------------------------
ComptonCone::ComptonCone(){
 TVector3 apex(0,0,0);
 TVector3 axis(-1,0,0);
 SetName("cone");
 SetApex(apex);
 SetAxis(axis);
 SetAngle(0); 
}
//------------------------------------------------------------------
ComptonCone::ComptonCone(TString name, TVector3 apex, TVector3 axis, Double_t angle){
  SetName(name);
  SetApex(apex);
  SetAxis(axis);
  SetAngle(angle);
}
//------------------------------------------------------------------
ComptonCone::~ComptonCone(){
}
//------------------------------------------------------------------
void ComptonCone::Print(void){
 cout << "\nComptonCone::Print() for object " << GetName() << endl;
 cout << "\tCoordinates of the apex: \n\t";
 fApex.Print();
 cout << "\tCoordinates of the leading versor of the axis: \n\t";
 fAxis.Print();
 cout << "\tAperture angle: " << fAngle << " radians\n" << endl; 
}
//------------------------------------------------------------------
