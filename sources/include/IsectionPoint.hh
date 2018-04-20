#ifndef __IsectionPoint_H_
#define __IsectionPoint_H_ 1 
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

class IsectionPoint : public TObject{
  
public:
  
  IsectionPoint();
  IsectionPoint(Int_t bin, TVector3 *pos);
  IsectionPoint(Int_t bin, Double_t x, Double_t y, Double_t z);
  ~IsectionPoint();
  
  void SetBinPoint(Int_t bin, Double_t x, Double_t y, Double_t z);
  void SetPointCoordinates(Double_t x, Double_t y, Double_t z);
  void SetBin(Int_t b);
  void Print(void);
  
  TVector3* GetPointCoordinates(void); //{ return fPoint; };
  Int_t GetBin(void) const; //{ return fGlobalBin; };
  
private:
  
 
  TVector3 *fPoint;
  Int_t    fGlobalBin;
  
  Int_t  Compare(const TObject* run2) const;
  Bool_t IsSortable() const {return kTRUE;};
  
  
  
  
  ClassDef(IsectionPoint,0)
  
};

#endif