#ifndef __InputReaderNN_H_
#define __InputReaderNN_H_ 1
#include "InputReader.hh"
#include <iostream>

using namespace std;

/// Class for accessing data from the selection done by the NN.
/// This is class derived from InputReader class.
/// A ROOTfile is opened a tree containing the selected events and a tree
/// containing information about the SiFi-CC setup are loaded.
class InputReaderNN : public InputReader
{

public:
    InputReaderNN();
    InputReaderNN(TString path);
    ~InputReaderNN();

  bool LoadEvent(int i);
  int GetNumberOfEventsInFile(void);
  void SelectEvents(void);  
  list<int> GetSelectedEvents(void);
  void Clear(void);
  
  TVector3* GetPositionScattering(void);
  TVector3* GetPositionAbsorption(void);
  TVector3* GetGammaDirScattered(void);
  
  double GetEnergyPrimary(void);
  double GetEnergyLoss(void);
  double GetEnergyScattered(void);
  
  void SetLoadOnlyCorrect(int value);

private:
  int fCorrectOnly; // 0 all results from the NN, 1 loading eventtype 1 and 2, 2 loading eventtype 2, 3 only eventtype 0, 4 only eventtype 1 
  list<int> fSelectedEvents;
  
  float fX1;///< x dimension of electron interaction 
  float fY1;///< y dimension of electron interaction
  float fZ1;///< z dimension of electron interaction
  float fX2;///< x dimension of photon interaction
  float fY2;///< y dimension of photon interaction
  float fZ2;///< z dimension of photon interaction
  float fX3;///< x dimension of third interaction in the detector, artifact from the file form from Luebeck here always =0 
  float fY3;///< y dimension of third interaction in the detector, artifact from the file form from Luebeck here always =0
  float fZ3;///< z dimension of third interaction in the detector, artifact from the file form from Luebeck here always =0
  float fVX;///< x dimension of electron interaction, the same as x1 
  float fVY;///< y dimension of electron interaction, the same as y1
  float fVZ;///< z dimension of electron interaction, the same as z1
  float fVUncX;///< uncertainty of the x dimension of electron interaction, in the case of the current NN always =0 
  float fVUncY;///< uncertainty of the y dimension of electron interaction, in the case of the current NN always =0
  float fVUncZ;///< uncertainty of the z dimension of electron interaction, in the case of the current NN always =0
  float fPX;///< x dimension of photon interaction, the same as x2
  float fPY;///< y dimension of photon interaction, the same as y2
  float fPZ;///< z dimension of photon interaction, the same as z2
  float fPUncX;///< uncertainty of the x dimension of photon interaction, in the case of the current NN always =0 
  float fPUncY;///< uncertainty of the y dimension of photon interaction, in the case of the current NN always =0
  float fPUncZ;///< uncertainty of the z dimension of photon interaction, in the case of the current NN always =0
  float fE0Calc; ///< primary energy calculated from electron energy+photon energy (leaf is artifact from 3 interaction model from Luebeck)
  float fE0CalcUnc; ///< uncertainty to E0Calc, here always =0
  float fArc;///< scattering angle calculated from electron energy and photon energy
  float fArcUnc;///< uncertainty to Arc, here always =0
  float fE1;///< energy of the electron
  float fE1Unc;///< uncertainty to E1, here always =0
  float fE2;///< energy of the photon 
  float fE2Unc;///< uncertainty to E2, here always =0
  float fE3;///< energy of the third interaction in the detector (3 interaction model to recover full energy), artefact from Luebeck 
  float fE3Unc;///< uncertainty to E3, here always =0
  int fClassID;///< Artefact from Luebeck, not meaningfull here 
  int fEventType;///< Quality of the event: =0 nether energy not position correct, =1 position correct energy not, =2 energy and position correct 
  int fEnergyBinID;///< Artefact from Luebeck, not meaningfull here 
  
  TVector3* fPositionScat;  ///< Position of interaction in scatterer
  TVector3* fPositionAbs;   ///< Position of interaction in absorber
  TVector3* fDirectionScat; ///< Direction of scattered gamma
  double fPrimaryEnergy;
  double fEnergyLoss;
  double fEnergyScattered;

    bool AccessTree(TString name);
    TTree* fTree;

    void Rotate(TVector3* vec);

    ClassDef(InputReaderNN, 0)
};

inline int InputReaderNN::GetNumberOfEventsInFile(void){return fTree->GetEntries();}
inline list<int> InputReaderNN::GetSelectedEvents(void){return fSelectedEvents;}
inline void InputReaderNN::SetLoadOnlyCorrect(int value){fCorrectOnly=value;}
inline TVector3* InputReaderNN::GetPositionScattering(void) {return fPositionScat;}
inline TVector3* InputReaderNN::GetPositionAbsorption(void) {  return fPositionAbs;}
inline TVector3* InputReaderNN::GetGammaDirScattered(void) {return fDirectionScat;}
inline double InputReaderNN::GetEnergyPrimary(void) {return fPrimaryEnergy;}
inline double InputReaderNN::GetEnergyLoss(void) { return fEnergyLoss; }
inline double InputReaderNN::GetEnergyScattered(void) { return fEnergyScattered; }

#endif
