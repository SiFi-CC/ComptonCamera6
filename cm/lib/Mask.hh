#ifndef __Mask_H_
#define __Mask_H_ 1
#include "DetPlane.hh"
#include "TH2F.h"
#include "TVector3.h"

class Mask : public DetPlane {

public:
  Mask();
  Mask(Double_t a, Double_t b, Double_t c, Double_t d, Double_t dimZ,
       Double_t dimY, TH2F* h, TString name);
  ~Mask();
  void SetPattern(TH2F* h);
  TH2F* GetPattern(void) { return fPattern; };
  Bool_t IsOpaque(TVector3 point);
  void Print(void);

private:
  TH2F* fPattern;

  ClassDef(Mask, 1)
};

#endif
