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
/// Two Types of information can be loaded:
/// 1. the real information from the simulation
/// 2. the information of the reconstructed events (only identified ones are
/// passed) Which one is loaded is toggled with SetLoadMonteCarlo()
class InputReaderGeant : public InputReader
{

public:
    /// Default constructor.
    InputReaderGeant();
    /// Standard constructor.
    ///\param path (TString) - path to the input file.
    InputReaderGeant(TString path);
    /// Default destructor.
    ~InputReaderGeant();

    /// loads events from trees to analyze them
    ///\param i (int) - number of events
    bool LoadEvent(int i);
    int GetNumberOfEventsInFile(void);
    void Clear(void);
    TVector3* GetPositionPrimary(void);
    TVector3* GetPositionScattering(void);
    TVector3* GetPositionAbsorption(void);
    TVector3* GetGammaDirPrimary(void);
    TVector3* GetGammaDirScattered(void);

    double GetEnergyPrimary(void);
    double GetEnergyLoss(void);
    double GetEnergyScattered(void);

    void SetLoadMonteCarlo(void);
    void SetLoadOnlyCorrect(void);

private:
    // TREEVARIABLES
    int fEventNumber; ///< Event number
    int fIdentified;  ///< Number of events were labeled
    bool fPurCrossed;
    Double_t fEnergy_Primary; ///< Primary photon energy
    Double_t fRealEnergy_e;   ///< Electron energy + uncertainty [MeV]
    Double_t fRealEnergy_p;   ///< Photon energy + uncertainty [MeV]
    TVector3* fRealPosition_source;
    TVector3* fRealDirection_source;
    // TVector3* fRealPosition_e;        ///< Electron creation position +
    // uncertainty
    TVector3* fRealComptonPosition;
    // TVector3* fRealPosition_p;      ///< Photon energy deposition position +
    // uncertainty
    TVector3* fRealDirection_scatter;         ///< Direction of the scattered photon +
                                              ///< uncertainty
    PhysicVar* fRecoEnergy_e;                 ///< Electron energy + uncertainty [MeV]
    PhysicVar* fRecoEnergy_p;                 ///< Photon energy + uncertainty [MeV]
    PhysicVec* fRecoPosition_e;               ///< Electron creation position + uncertainty
    PhysicVec* fRecoPosition_p;               ///< Photon energy deposition position + uncertainty
    PhysicVec* fRecoDirection_scatter;        ///< Direction of the scattered photon +
                                              ///< uncertainty
    vector<PhysicVec>* fRecoClusterPositions; ///< Positions cluster with uncertainties
    vector<PhysicVar>* fRecoClusterEnergies;  ///< Energies cluster with uncertainties

    vector<TVector3>* fRealPosition_e;
    vector<TVector3>* fRealPosition_p;

    vector<int>* fRealInteractions_e;
    vector<int>* fRealInteractions_p;

    Long64_t fNumberOfSimulatedEvents;
    ////////////////////////////////////

    // CLASSVARIABLES

    bool fLoadReal;
    bool fCorrectOnly;

    TVector3* fPositionSource;
    TVector3* fDirectionSource;

    TVector3* fScattererPosition;
    TVector3* fAbsorberPosition;

    TVector3* fScattererDim;
    TVector3* fAbsorberDim;

    TVector3* fPositionScat;  ///< Position of interaction in scatterer
    TVector3* fPositionAbs;   ///< Position of interaction in absorber
    TVector3* fDirectionScat; ///< Direction of scattered gamma
    double fPrimaryEnergy;
    double fEnergyLoss;
    double fEnergyScattered;
    /// Accesses data of trees'branches in ROOT file.
    ///\param name (TString) - name of tree.
    ///\param namesetup (TString) - name of tree holding setp information.
    bool AccessTree(TString name, TString namesetup);
    TTree* fTree;
    TTree* fTreeSetup;

    ClassDef(InputReaderGeant, 0)
};

inline int InputReaderGeant::GetNumberOfEventsInFile(void) { return fTree->GetEntries(); }
inline void InputReaderGeant::SetLoadMonteCarlo(void) { fLoadReal = true; }
inline void InputReaderGeant::SetLoadOnlyCorrect(void) { fCorrectOnly = true; }
inline TVector3* InputReaderGeant::GetPositionScattering(void) { return fPositionScat; }
inline TVector3* InputReaderGeant::GetPositionAbsorption(void) { return fPositionAbs; }
inline TVector3* InputReaderGeant::GetGammaDirScattered(void) { return fDirectionScat; }
inline double InputReaderGeant::GetEnergyPrimary(void) { return fPrimaryEnergy; }
inline double InputReaderGeant::GetEnergyLoss(void) { return fEnergyLoss; }
inline double InputReaderGeant::GetEnergyScattered(void) { return fEnergyScattered; }

#endif
