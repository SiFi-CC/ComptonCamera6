#ifndef __InputReaderGeant_H_
#define __InputReaderGeant_H_ 1
#include "DR_GenerallStructs.hh"
#include "G4Input.hh"
#include "InputReader.hh"
#include <iostream>

using namespace std;

/// Class for accessing data from Geant4 simulation performed by JK.
/// This is class derived from InputReader class. It opens requested
/// ROOT file containing tree with simulation results and via set of
/// getter function passes information to reconstruction classes, i.e.
/// CCREconstruction and CCMLEM.

class InputReaderGeant : public InputReader {

public:
  InputReaderGeant();
  InputReaderGeant(TString path);
  ~InputReaderGeant();

  bool LoadEvent(int i);
  void Clear(void);
  TVector3* GetPositionPrimary(void);
  TVector3* GetPositionScattering(void);
  TVector3* GetPositionAbsorption(void);
  TVector3* GetGammaDirPrimary(void);
  TVector3* GetGammaDirScattered(void);
  double GetEnergyPrimary(void);
  double GetEnergyLoss(void);
  double GetEnergyScattered(void);

  TVector3* GetScattererPosition(void);
  TVector3* GetAbsorberPosition(void);

  double GetScatThickx(void);
  double GetScatThicky(void);
  double GetScatThickz(void);
  double GetAbsThickx(void);
  double GetAbsThicky(void);
  double GetAbsThickz(void);

private:
  int fEventNumber; ///< Event number
  bool fIdentified; ///< Flag indicating whether the event was labeled or not
  Double_t fEnergy_Primary;
  Double_t fRealEnergy_e; ///< Electron energy + uncertainty [MeV]
  Double_t fRealEnergy_p; ///< Photon energy + uncertainty [MeV]
  TVector3* fRealPosition_source;
  TVector3* fRealDirection_source;
  TVector3* fRealPosition_e; ///< Electron creation position + uncertainty
  TVector3*
      fRealPosition_p; ///< Photon energy deposition position + incertainty
  TVector3* fRealDirection_scatter; ///< Direction of the scattered photon +
                                    ///< uncertainty
  // vector<PhysicVec*>*
  // fRecoClusterPositions; ///< Positions cluster with uncertainties
  // vector<PhysicVar*>*
  // fRecoClusterEnergies; ///< Energies cluster with uncertainties

  TVector3* fPositionScat;  ///< Position of interaction in scatterer
  TVector3* fPositionAbs;   ///< Position of interaction in absorber
  TVector3* fDirectionScat; ///< Direction of scattered gamma
  TVector3* fPositionSource;
  TVector3* fDirectionSource;

  TVector3* fScattererPosition;
  TVector3* fAbsorberPosition;

  TVector3* fScatPlanePos;
  TVector3* fAbsPlanePos;

  Double_t fScattererThickness_x;
  Double_t fScattererThickness_y;
  Double_t fScattererThickness_z;
  Double_t fAbsorberThickness_x;
  Double_t fAbsorberThickness_y;
  Double_t fAbsorberThickness_z;
  Double_t fNumberOfSimulatedEvents;

  bool AccessTree(TString name, TString name1);
  TTree* fTree1;

  ClassDef(InputReaderGeant, 0)
};

#endif
