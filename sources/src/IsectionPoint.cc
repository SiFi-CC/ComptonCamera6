#include "IsectionPoint.hh"
#include <iostream>

using namespace std;

ClassImp(IsectionPoint);

//------------------------------------------------------------------
IsectionPoint::IsectionPoint(){
 fGlobalBin = 0;
 fPoint = new TVector3(0,0,0);
}
//------------------------------------------------------------------
IsectionPoint::IsectionPoint(Int_t bin, TVector3 *pos){
  fGlobalBin=bin;
  fPoint = new TVector3(pos->X(), pos->Y(), pos->Z());
}
//------------------------------------------------------------------
IsectionPoint::IsectionPoint(Int_t bin, Double_t x, Double_t y, Double_t z){
  fGlobalBin=bin;
  fPoint = new TVector3(x, y, z);
}
//------------------------------------------------------------------
IsectionPoint::~IsectionPoint(){
  delete fPoint;
}
//------------------------------------------------------------------
void IsectionPoint::SetBinPoint(Int_t bin, Double_t x, Double_t y, Double_t z){
  fGlobalBin=bin;
  fPoint->SetXYZ(x, y, z);
}
//------------------------------------------------------------------
void IsectionPoint::SetPointCoordinates(Double_t x, Double_t y, Double_t z){
  fPoint->SetXYZ(x, y, z);
}
//------------------------------------------------------------------
void IsectionPoint::SetBin(Int_t b){
  fGlobalBin=b;
}
//------------------------------------------------------------------
TVector3* IsectionPoint::GetPointCoordinates(void){
  return fPoint;
}
//------------------------------------------------------------------
Int_t IsectionPoint::GetBin(void) const{
  
  return fGlobalBin;
}
//------------------------------------------------------------------
Int_t IsectionPoint::Compare(const TObject* run2o) const{
  if(!run2o){
    cout<<"IsectionPoint::operator> : you attempt to compare to a non-existing object!"<<endl;
  }
  const IsectionPoint*  run2 = static_cast<const IsectionPoint*>(run2o);
  if(this->GetBin() < run2->GetBin())
    return -1;
  else if(this->GetBin() > run2->GetBin())
    return 1;
  else
    return 0;
}
//------------------------------------------------------------------
void IsectionPoint::Print(void){
 cout << "\nIsectionPoint::Print()" <<endl;
 cout << "\tGlobal bin number: \t"<<fGlobalBin<<endl;
 cout << "\tPoint position: "<<fPoint->X() << ", " << fPoint->Y() <<", " << fPoint->Z() << endl;
}
//------------------------------------------------------------------
