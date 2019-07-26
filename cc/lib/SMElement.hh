#ifndef __SMElement_H_
#define __SMElement_H_ 1
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"
/// Class which gets event number and global bin number and length of track
/// used in class CCMLEM.
class SMElement : public TObject {

public:
  SMElement();
  SMElement(Int_t eventno, Int_t bin, Double_t dist);
  ~SMElement();
  
  /// Sets event number and global bin number and length of track.
  void SetEvBinDist(Int_t eventno, Int_t bin, Double_t dist);
  /// Sets event number.
  void SetEvent(Int_t i);
  /// Sets global bin number.
  void SetBin(Int_t b);
  /// Sets length of track.
  void SetDist(Double_t d);
  
  void Print(void);
  /// Returns event number.
  Int_t GetEvent(void) const;
  /// Returns global bin number.
  Int_t GetBin(void) const;
  /// Returns length of track.
  Double_t GetDist(void);

private:
  Int_t feventno;       ///< event number
  Int_t fGlobalBin;     ///< Global bin number
  Double_t fdist;       ///< length of track

  ClassDef(SMElement, 0)
};

#endif
