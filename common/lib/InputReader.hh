#ifndef __InputReader_H_
#define __InputReader_H_ 1
#include "TFile.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"
#include <iostream>

using namespace std;

/// Base class for InputReaderSimple and InputReaderGeant. It allows
/// for reading simulation data saved in different forms and passing it
/// to reconstruction classes via set of getter functions.

class InputReader : public TObject {

public:
  InputReader();
  InputReader(TString path);
  ~InputReader();

  void Print(void);
  void virtual Clear(void);
  bool virtual LoadEvent(int i);
  TVector3 virtual* GetPositionPrimary(void);
  TVector3 virtual* GetPositionScattering(void);
  TVector3 virtual* GetPositionAbsorption(void);
  TVector3 virtual* GetGammaDirPrimary(void);
  TVector3 virtual* GetGammaDirScattered(void);
  TVector3 virtual* GetScattererPosition(void);
  TVector3 virtual* GetAbsorberPosition(void);
  
  double virtual GetEnergyPrimary(void);
  double virtual GetEnergyLoss(void);
  double virtual GetEnergyScattered(void);
  double virtual GetScatThickx(void);
  double virtual GetScatThicky(void);
  double virtual GetScatThickz(void);
  double virtual GetAbsThickx(void);
  double virtual GetAbsThicky(void);
  double virtual GetAbsThickz(void);

protected:
  TFile* fFile; ///< Input file
  TTree* fTree; ///< Tree containing simulation results
  bool virtual AccessTree(TString name);
  bool SetInputFile(TString path);
  
  ClassDef(InputReader, 0)
};

#endif
