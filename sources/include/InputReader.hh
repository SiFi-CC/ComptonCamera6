#ifndef __InputReader_H_
#define __InputReader_H_ 1
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TVector3.h"
#include <iostream>

using namespace std;

class InputReader : public TObject{
  
public:
   InputReader();
   InputReader(TString path);
  ~InputReader();
  
  bool     LoadEvent(int i);
  void     Clear(void);
  void     Print(void);
  TVector3 virtual *GetPositionPrimary(void);
  TVector3 virtual *GetPositionScattering(void);
  TVector3 virtual *GetPositionAbsorption(void);
  TVector3 virtual *GetGammaDirPrimary(void);
  TVector3 virtual *GetGammaDirScattered(void);
  double   virtual GetEnergyPrimary(void);
  double   virtual GetEnergyLoss(void);
  double   virtual GetEnergyScattered(void);
  
protected:
  TFile   *fFile;
  TTree   *fTree;
  
  bool    virtual AccessTree(TString name);
  bool    SetInputFile(TString path);
  
  ClassDef(InputReader,0)
};

#endif
