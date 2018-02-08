#ifndef __PhysicsBase_H_
#define __PhysicsBase_H_ 1
#include "TObject.h"
#include "TF1.h"
#include "Track.hh"
#include "DetPlane.hh"

class PhysicsBase : public TObject{
  
public:
  PhysicsBase();
  PhysicsBase(TString name);
  ~PhysicsBase();
  
  Double_t NewEnergy(Double_t theta, Double_t initE);
  Double_t FindPhi(void);
  Double_t FindTheta(Double_t energy);
  Track*   ComptonScatter(Track *initTrack, DetPlane *plane);
  void     Print(void);
  void     SetName(TString name){ fName = name; };
  const char*  GetName() const { return fName.Data(); };
  Double_t GetTheta(void){ return fTheta; };
  Double_t GetPhi(void){ return fPhi; };
  
private:
  
  TString  fName;
  Double_t fTheta;
  Double_t fPhi;
  TF1*     fFunction;
  
  ClassDef(PhysicsBase,0)
};

#endif