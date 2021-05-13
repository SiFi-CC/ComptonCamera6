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
  vector<TVector3> virtual* GetElectronPosition(void);
  int virtual GetRealPosPSize(void);
  vector<int> virtual* GetRealInteractionP(void);
  int virtual GetRealPosESize(void);
  vector<int> virtual* GetRealInteractionE(void);
  TVector3 virtual* GetPositionScatteringReco(void);
  vector<TVector3> virtual* GetPhotonPosition(void);
  TVector3 virtual* GetPositionAbsorption(void);
  TVector3 virtual* GetPositionAbsorptionReco(void);
  TVector3 virtual* GetGammaDirPrimary(void);
  TVector3 virtual* GetGammaDirScattered(void);
  TVector3 virtual* GetGammaDirScatteredReco(void);
  TVector3 virtual* GetScattererPosition(void);
  TVector3 virtual* GetAbsorberPosition(void);
  
  int virtual GetIdentified(void);
  int virtual GetMultiplicityNum(void);
  int virtual GetClassID(void);
  int virtual GetRecoClusterPosSize(void);
  double virtual GetEP(void);
  double virtual GetReES(void);
  double virtual GetES(void);
  double virtual GetReEP(void);
  double virtual GetEnergyPrimary(void);
  double virtual GetEnergyPrimaryReco(void);
  double virtual GetEnergyLoss(void);
  double virtual GetEnergyLossReco(void);
  double virtual GetEnergyScattered(void);
  double virtual GetEnergyScatteredReco(void);
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
  //vector<PhysicVec*>* GetRecoClusterPosSize(void);
  
  ClassDef(InputReader, 0)
};

#endif
