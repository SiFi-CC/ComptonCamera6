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

private:
  int fEventNumber; ///< Event number
  bool fIdentified; ///< Flag indicating whether the event was labeled or not
  PhysicVar* fRecoEnergy_e;   ///< Electron energy + uncertainty [MeV]
  PhysicVar* fRecoEnergy_p;   ///< Photon energy + uncertainty [MeV]
  PhysicVec* fRecoPosition_e; ///< Electron creation position + uncertainty
  PhysicVec*
      fRecoPosition_p; ///< Photon energy deposition position + incertainty
  PhysicVec* fRecoDirection_scatter; ///< Direction of the scattered photon +
                                     ///< uncertainty
  vector<PhysicVec*>*
      fRecoClusterPositions; ///< Positions cluster with uncertainties
  vector<PhysicVar*>*
      fRecoClusterEnergies; ///< Energies cluster with uncertainties

  TVector3* fPositionScat;  ///< Position of interaction in scatterer
  TVector3* fPositionAbs;   ///< Position of interaction in absorber
  TVector3* fDirectionScat; ///< Direction of scattered gamma

  bool AccessTree(TString name);

  ClassDef(InputReaderGeant, 0)
};

#endif
