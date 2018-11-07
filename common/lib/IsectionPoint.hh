#ifndef __IsectionPoint_H_
#define __IsectionPoint_H_ 1 
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"
///Class which gets intersection point coordinates and global bin number to be used in class CCMLEM.  
class IsectionPoint : public TObject{
  
public:
  
  IsectionPoint();
  IsectionPoint(Int_t bin, TVector3 *pos);
  IsectionPoint(Int_t bin, Double_t x, Double_t y, Double_t z);
  ~IsectionPoint();
  
  ///Sets coordinates of intesection point and global bin number
  void SetBinPoint(Int_t bin, Double_t x, Double_t y, Double_t z);
  ///Sets coordinates of intesection point.
  void SetPointCoordinates(Double_t x, Double_t y, Double_t z);
  ///Sets global bin number.
  void SetBin(Int_t b);
  
  void Print(void);
  
  ///Returns coordinates of intesection point.
  TVector3* GetPointCoordinates(void);
  ///Returns global bin number.
  Int_t GetBin(void) const; 
  
private:
  
 
  TVector3 *fPoint;				///< Coordinates of the interaction point
  Int_t    fGlobalBin;				///< Global bin number
  
  Int_t  Compare(const TObject* run2) const;	///< Compares sortable objects (global bin numbers)
  Bool_t IsSortable() const {return kTRUE;};	///< Returns kTRUE for all sortable objects
  
  
  
  
  ClassDef(IsectionPoint,0)
  
};

#endif