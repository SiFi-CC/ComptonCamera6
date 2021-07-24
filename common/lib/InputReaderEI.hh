#ifndef __InputReaderEI_H_
#define __InputReaderEI_H_ 1
#include "InputReader.hh"
#include <iostream>
#include <vector>
#include "DR_GenerallStructs.hh"
#include "G4Input.hh"
using namespace std;

/// This is class derived from InputReader class. It opens requested
/// ROOT file containing tree with reconstruction results by for example, Geant4
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
  double GetEP(void);
  double GetReES(void);
//  double GetReEP(void);
  double GetES(void);
  
  int GetMultiplicityNum(void);
  int GetClassID(void);

  

private:
    
  int fEventNumber;     ///< Event number
  int fS;
  int fClassID;
  
  TVector3* fPositionElectron;      ///< Position in scatterer module
  TVector3* fPositionPhoton;       ///< Position in absorber module
  
 
  
  TVector3* fPos_Scat;
  TVector3* fPos_Abs;
   
  Double_t fEnergy_Scat;
  Double_t fEnergy_Abs;

  Double_t fEnergy_Primary;
  Double_t fReEnS;
  Double_t fEnergyS;
//  Double_t fReEnP;
  
  bool AccessTree(TString name/*, TString name1, TString name2*/);
  TTree* fTree;
  

  ClassDef(InputReaderEI, 0)
};

#endif

