#ifndef __InputReaderNN_H_
#define __InputReaderNN_H_ 1
#include "InputReader.hh"
#include <iostream>

using namespace std;

/// Class for accessing data from the selection done by the NN.
/// This is class derived from InputReader class.
/// A ROOTfile is opened a tree containing the selected events and a tree containing information about the SiFi-CC setup are loaded. 
class InputReaderNN: public InputReader {

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
  int fCorrectOnly;
  list<int> fSelectedEvents;
  
  float fX1;
  float fY1;
  float fZ1;
  float fX2;
  float fY2;
  float fZ2;
  float fX3;
  float fY3;
  float fZ3;
  float fVX;
  float fVY;
  float fVZ;
  float fVUncX;
  float fVUncY;
  float fVUncZ;
  float fPX;
  float fPY;
  float fPZ;
  float fPUncX;
  float fPUncY;
  float fPUncZ;
  float fE0Calc;
  float fE0CalcUnc;
  float fArc;
  float fArcUnc;
  float fE1;
  float fE1Unc;
  float fE2;
  float fE2Unc;
  float fE3;
  float fE3Unc;
  int fClassID;
  int fEventType;
  int fEnergyBinID;
  
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
