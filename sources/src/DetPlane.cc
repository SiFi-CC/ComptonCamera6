#include "DetPlane.hh"
#include <iostream>

using namespace std;

ClassImp(DetPlane);

//------------------------------------------------------------------
DetPlane::DetPlane(){
  SetPlane(0,0,1,0);
  SetDimensions(10,10);
  SetName("plane");
}
//------------------------------------------------------------------
DetPlane::DetPlane(Double_t a, Double_t b, Double_t c, Double_t d, 
		   Double_t dimZ, Double_t dimY, TString name){
  
  SetPlane(a,b,c,d);
  SetDimensions(dimZ,dimY);
  SetName(name);
}
//------------------------------------------------------------------
DetPlane::~DetPlane(){
}
//------------------------------------------------------------------
void DetPlane::SetPlane(Double_t a, Double_t b, Double_t c, Double_t d){
  
  fA = a;
  fB = b;
  fC = c;
  fD = d;
  
  if(fA==0 && fB==0 && fC==0){
   cout << "##### Parameters incorrect! A, B and C cannot be equal 0 at the same time!" << endl;
   cout << "##### Please correct!" << endl;
   return;
  }
}
//------------------------------------------------------------------
void DetPlane::SetDimensions(Double_t dimZ, Double_t dimY){
 fDimZ = dimZ;
 fDimY = dimY;
}
//------------------------------------------------------------------
TVector3 DetPlane::GetNormal(void){
  
  TVector3 norm;
  
  Double_t x = fA;
  Double_t y = fB;
  Double_t z = fC;
  Double_t n = sqrt(fA*fA+fB*fB+fC*fC);
  norm.SetXYZ(x/n,y/n,z/n);
  
  return norm;
}
//------------------------------------------------------------------
Bool_t DetPlane::CheckPoint(TVector3 point){
 
  if(fabs(point.Z())>0.5*fDimZ) return kFALSE;
  if(fabs(point.Y())>0.5*fDimY) return kFALSE;
  
  Double_t num   = fabs(fA*point.X() + 
                   fB*point.Y() + fC*point.Z() + fD);
  Double_t denom = sqrt(fA*fA + fB*fB + fC*fC);
  Double_t dist  = num/denom;
  
  if(dist!=0) return kFALSE;
  
  return kTRUE;
}
//------------------------------------------------------------------
void DetPlane::Print(void){
  cout << "\nDetPlane::Print() for object " << GetName() << endl;
  cout << "\tA = " << GetA() << endl;
  cout << "\tB = " << GetB() << endl;
  cout << "\tC = " << GetC() << endl;
  cout << "\tD = " << GetD() << endl;
  cout << "\tZ full size = " << GetDimZ() << endl;
  cout << "\tY full size = " << GetDimY() << endl;
  cout << "\tDimensions of the detector plane: "
       << fDimZ << "mm x " << fDimY << "mm" << endl;
  cout << "\tCartesian equation of the plane: "
       << fA << "x + " << fB << "y + " << fC << "z + " 
       << fD << " = 0" << endl; 
}
//------------------------------------------------------------------
