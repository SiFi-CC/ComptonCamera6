#ifndef __InputReaderEI_H_
#define __InputReaderEI_H_ 1
#include "InputReader.hh"
#include <iostream>
#include <vector>
#include "DR_GenerallStructs.hh"
#include "G4Input.hh"
using namespace std;

/// This is class derived from InputReader class. It opens requested
/// ROOT file containing tree with simulation results and via set of
/// getter function passes information to reconstruction classes, i.e.
/// CCREconstruction and CCMLEM.

class InputReaderEI : public InputReader {

public:
  InputReaderEI();
  InputReaderEI(TString path);
  ~InputReaderEI();

  bool LoadEvent(int i);
  void Clear(void);
  
  TVector3* GetPositionScattering(void);
  
  TVector3* GetPositionAbsorption(void);
  
  double GetEnergyLoss(void);
  double GetEnergyScattered(void);
  double GetTotalEnergy(void);
  
  int GetScatSize(void);
  int GetAbsSize(void);

  

private:
  
  TVector3* fPositionScatReco;
  TVector3* fPositionAbsReco;
  
  PhysicVar* fEnergyRecoScat2;
  PhysicVar* fEnergyRecoAbs2;
  
  PhysicVec* fPointRecoScat2;
  PhysicVec* fPointRecoAbs2;
  
  Int_t fSizeScat2;
  Int_t fSizeAbs2;
  
  Double_t fTotalEnergy2;
    
  bool AccessTree(TString name, TString name1);
  TTree* fTree;
  TTree* fTree1;

//  ClassDef(InputReaderEI, 0)
};

#endif

