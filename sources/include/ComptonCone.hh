#ifndef __ComptonCone_H_
#define __ComptonCone_H_ 1 
#include "TObject.h"
#include "TVector3.h"
#include "TString.h"

class ComptonCone : public TObject{
  
public:
  
  ComptonCone();
  ComptonCone(TString name, TVector3 apex, TVector3 axis, Double_t angle);
  ~ComptonCone();
  
  void     SetName(TString name) { fName = name; };
  void     SetApex(TVector3 apex){ fApex = apex; };
  void     SetAxis(TVector3 axis){ fAxis = axis; };
  void     SetAngle(Double_t angle){ fAngle = angle; };
  void     Print(void);
  const char*  GetName() const { return fName.Data(); }
  TVector3 GetAxis(void){ return fAxis; };
  TVector3 GetApex(void){ return fApex; };
  Double_t GetAngle(void){ return fAngle; };
  
private:
  TString  fName;
  Double_t fAngle; //in radians
  TVector3 fApex;
  TVector3 fAxis;
	   
  ClassDef(ComptonCone,0)
  
};

#endif
