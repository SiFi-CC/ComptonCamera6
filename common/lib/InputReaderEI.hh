#ifndef __InputReaderEI_H_
#define __InputReaderEI_H_ 1
#include "InputReader.hh"
#include <iostream>
#include <vector>
#include "DR_GenerallStructs.hh"
#include "G4Input.hh"
using namespace std;

/// This is class derived from InputReader class. It opens requested
/// ROOT file containing tree with reconstruction results by Geant4 
/// and via set of getter function passes information to reconstruction 
/// classes, i.e. CCMLEM.

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
  //double GetTotalEnergy(void);
  
  

  

private:
  int fEventNumber;     ///< Event number
  TVector3* fPositionScatReco;      ///< Position in scatterer module
  TVector3* fPositionAbsReco;       ///< Position in absorber module
  
  PhysicVar* fEnergyReco0;      ///< Deposited energy in scatterer module
  PhysicVar* fEnergyReco1;      ///< Deposited energy in absorber module
  
  PhysicVec* fPosScatClus;
  PhysicVec* fPosAbsClus;
  
  Double_t fEnergyRe0;
  Double_t fEnergyRe1;
  
//   TVector3* fPosScatReco;
//   TVector3* fPosAbsReco;
//   
//   Double_t fRecoEnergy_e;
//   Double_t fRecoEnergy_p;
//   
//   TVector3* fPosScatReal;
//   TVector3* fPosAbsReal;
//   
//   Double_t fRealEnergy_e;
//   Double_t fRealEnergy_p;
  
  
  bool AccessTree(TString name/*, TString name1, TString name2*/);
  TTree* fTree;
//   TTree* fTree1;
//   TTree* fTree2;

//  ClassDef(InputReaderEI, 0)
};

#endif

