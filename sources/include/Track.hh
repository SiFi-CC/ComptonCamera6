#ifndef __Track_H_
#define __Track_H_ 1
#include "TVector3.h"
#include "TString.h"
#include "TObject.h"

class DetPlane;

class Track : public TObject{
  
public:
  Track();
  Track(TVector3 point, TVector3 vec, Double_t energy, TString name);
  ~Track();
  
  void     SetPoint(TVector3 point);
  void     SetVersor(TVector3 vec);
  void     SetEnergy(Double_t energy);
  Bool_t   FindCrossPoint(DetPlane* plane, TVector3 &position);
  void     Print(void);
  Double_t GetEnergy(void) { return fEnergy; };
  TVector3 GetVersor(void){ return fVersor; };
  TVector3 GetPoint(void){ return fPoint; };
  void     SetName(TString name){fName = name;};
  const char*  GetName() const { return fName.Data(); }
  
private:
  
  TVector3 fPoint;
  TVector3 fVersor;
  Double_t fEnergy;
  TString  fName;
  
  ClassDef(Track,0)
  
};

#endif
