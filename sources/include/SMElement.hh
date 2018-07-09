#ifndef __SMElement_H_
#define __SMElement_H_ 1 
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

class SMElement : public TObject{
  
public:
  
  SMElement();
  SMElement(Int_t eventno, Int_t bin, Double_t dist);
  ~SMElement();
  
 
  void SetEvBinDist(Int_t eventno, Int_t bin, Double_t dist);
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
  
  
  
  
  
  
  ClassDef(SMElement,0)
  
};

#endif 
