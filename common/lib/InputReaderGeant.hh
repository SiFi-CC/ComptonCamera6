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
/// Here, we have two types information: 
/// 1. the real information from the simulation 
/// 2. the information of the reconstructed events 
class InputReaderGeant : public InputReader {

public:
  InputReaderGeant();
  InputReaderGeant(TString path);
  ~InputReaderGeant();

  bool LoadEvent(int i);
  void Clear(void);
  TVector3* GetPositionPrimary(void);
  TVector3* GetPositionScattering(void);
  TVector3* GetPositionScatteringReal(void);
  TVector3* GetPositionScatteringReco(void);
  TVector3* GetPositionAbsorption(void);
  TVector3* GetPositionAbsorptionReal(void);
  TVector3* GetPositionAbsorptionReco(void);
  TVector3* GetGammaDirPrimary(void);
  TVector3* GetGammaDirScattered(void);
  TVector3* GetGammaDirScatteredReal(void);
  TVector3* GetGammaDirScatteredReco(void);
  int GetRecoClusterPosSize(void);
   
  double GetEnergyPrimary(void);
  double GetEnergyLoss(void);
  double GetEnergyLossReal(void);
  double GetEnergyLossReco(void);
  double GetEnergyScattered(void);
  double GetEnergyScatteredReal(void);
  double GetEnergyScatteredReco(void);
  
  
  TVector3* GetScattererPosition(void);
  TVector3* GetAbsorberPosition(void);

  double GetScatThickx(void);
  double GetScatThicky(void);
  double GetScatThickz(void);
  double GetAbsThickx(void);
  double GetAbsThicky(void);
  double GetAbsThickz(void);

  /// Set if real or reconstructed information are to be read from tree 
  ///\see fUseRealInformation
  void SetUseRealInformation(bool b){fUseRealInformation = b;};

private:
  int fEventNumber;     ///< Event number
  int fEventNumberReco;     ///< Event number from tree with reconstructed events, should be the same as fEventNumber
  bool fBIdentified;      ///< Flag if reconstruction found one cluster in absorber und one cluster in scatterer (used until summer 2019)
  /// Identifier which event configuration the reconstruction found (used from summer 2019)
  /// 0: did not match one of the classes below
  /// +-1: S1A1-EP
  /// +-2: S1A1-PE
  /// +-3: S1AX
  /// A positive value means that the classification was correct, a negative value means, it was false 
  Int_t fIIdentified;      
  Double_t fEnergy_Primary;				///< Primary photon energy from simulation [MeV]
  Double_t fRealEnergy_e;				///< Electron energy [MeV]
  Double_t fRealEnergy_p;				///< Photon energy [MeV]
  TVector3* fRealPosition_source;
  TVector3* fRealDirection_source;
  TVector3* fComptonPosition;			///< Position of the Compton effect [mm]
  TVector3* fRealPosition_e;			///< Electron creation position [mm]
  TVector3* fRealPosition_p;			///< Photon energy deposition position [mm]
  TVector3* fRealDirection_scatter;		///< Direction of the scattered photon [mm]
  PhysicVar* fRecoEnergy_e;				///< Electron energy + uncertainty [MeV]
  PhysicVar* fRecoEnergy_p;				///< Photon energy + uncertainty [MeV]
  PhysicVec* fRecoPosition_e;			///< Electron creation position + uncertainty
  PhysicVec* fRecoPosition_p;			///< Photon energy deposition position + uncertainty
  PhysicVec* fRecoDirection_scatter;	///< Direction of the scattered photon + uncertainty
  vector<PhysicVec*>*
      fRecoClusterPositions; ///< Positions cluster with uncertainties
  vector<PhysicVar*>*
      fRecoClusterEnergies; ///< Energies cluster with uncertainties


  TVector3* fPositionScat;  ///< Position of interaction in scatterer
  TVector3* fPositionAbs;   ///< Position of interaction in absorber
  TVector3* fDirectionScat; ///< Direction of scattered gamma
  
  TVector3* fPositionScatReco;
  TVector3* fPositionAbsReco;
  TVector3* fDirectionScatReco;
  
  
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

  ///\brief Switch to use either real or reconstructed information from simulation
  ///\details The status of this variable defines, which events are returned by
  /// the standard getter functions like GetEnergyLoss().
  /// Both real and recontructed events can still be accessed by the
  /// specific functions like GetEnergyLossReal() or GetEnergyLossReco().
  /// The default behavior is to return the real information.
  bool fUseRealInformation;

  bool AccessTree();
  TTree* fTreeSetup;	///< Tree containing placement of scatterer and absorber
  TTree* fTreeReco;		///< Separate tree with reconstructed events (only used in structure before summer 2019)

  bool fJointTree;		///< Internal flag, if structure with two trees for real and reconstructed events is used or if they are both stored in one common tree

  ClassDef(InputReaderGeant, 0)
};

#endif
