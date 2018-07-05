#ifndef __SMMLEM_H_
#define __SMMLEM_H_ 1 
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

class SMMLEM : public TObject{
  
public:
  
  SMMLEM();
  SMMLEM(Int_t eventno, Int_t bin, Double_t dist);
  ~SMMLEM();
  
  //Double_t Distance(TVector3 point1, TVector3 point2);
  void SetEvBinDist(Int_t eventno, Int_t bin, Double_t dist);
  void SetEvBin(Int_t eventno, Int_t bin);
  void SetEvent(Int_t i);
  void SetBin(Int_t b);
  void SetDist(Double_t d);
  void Print(void);
  
  Int_t GetEvent(void) const; 
  Int_t GetBin(void) const;
  Double_t GetDist(void);
  
private:
  
 
  Int_t    feventno;
  Int_t    fGlobalBin;
  Double_t fdist;
  
  
  
  
  
  
  ClassDef(SMMLEM,0)
  
};

#endif 
