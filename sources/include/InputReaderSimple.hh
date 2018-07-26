#ifndef __InputReaderSimple_H_
#define __InputReaderSimple_H_ 1
#include "TString.h"
#include "TVector3.h"
#include "InputReader.hh"
#include <iostream>

using namespace std;

class InputReaderSimple : public InputReader{
  
public:
  InputReaderSimple();
  InputReaderSimple(TString path);
  ~InputReaderSimple();
 
  TVector3 *GetPositionPrimary(void);
  TVector3 *GetPositionScattering(void);
  TVector3 *GetPositionAbsorption(void);
  TVector3 *GetGammaDirPrimary(void);
  TVector3 *GetGammaDirScattered(void);
  double   GetEnergyPrimary(void);
  double   GetEnergyLoss(void);
  double   GetEnergyScattered(void);
  
private:
  TVector3 *fPoint0;
  TVector3 *fPoint1;
  TVector3 *fPoint2;
  TVector3 *fVersor1;
  TVector3 *fVersor2;
  double   fEnergy0;
  double   fEnergy1;
  double   fEnergy2;
  
  bool     AccessTree(TString name);
  
  ClassDef(InputReaderSimple,0)
};

#endif