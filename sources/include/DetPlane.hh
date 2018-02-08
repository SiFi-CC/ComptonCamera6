#ifndef __DetPlane_H_
#define __DetPlane_H_ 1
#include "TVector3.h"
#include "TObject.h"
#include "TString.h"

class DetPlane : public TObject{
 
public:
  DetPlane();
  DetPlane(Double_t a, Double_t b, Double_t c, Double_t d, 
	   Double_t dimZ, Double_t dimY, TString name);
  ~DetPlane();
  
  void     SetPlane(Double_t a, Double_t b, Double_t c, Double_t d);
  void     SetDimensions(Double_t dimZ, Double_t dimY);
  void     SetName(TString name){fName = name;};
  TVector3 GetNormal(void);
  Bool_t   CheckPoint(TVector3 point);
  void     Print(void);
  Double_t GetA(void){ return fA; };
  Double_t GetB(void){ return fB; };
  Double_t GetC(void){ return fC; };
  Double_t GetD(void){ return fD; };
  Double_t GetDimZ(void){ return fDimZ; };
  Double_t GetDimY(void){ return fDimY; };
  const char*  GetName() const { return fName.Data(); }
  
private:
  Double_t fA;
  Double_t fB;
  Double_t fC;
  Double_t fD;
  Double_t fDimZ; //full length in mm
  Double_t fDimY; //full length in mm
  TString  fName;
  
  ClassDef(DetPlane,0)
  
};

#endif
