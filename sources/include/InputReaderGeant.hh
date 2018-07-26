#ifndef __InputReaderGeant_H_
#define __InputReaderGeant_H_ 1
#include "InputReader.hh"
#include "DR_GenerallStructs.hh"
#include "G4Input.hh"
#include <iostream>

using namespace std;

class InputReaderGeant : public InputReader{
  
public:
  InputReaderGeant();
  InputReaderGeant(TString path);
  ~InputReaderGeant();
  
  TVector3 *GetSourcePosition(void);
  TVector3 *GetScatPosition(void);
  TVector3 *GetAbsPosition(void);
  TVector3 *GetPrimaryGammaDir(void);
  TVector3 *GetScatGammaDir(void);
  double    GetEnSource(void);
  double    GetEnScat(void);
  double    GetEnAbs(void);
  
private:
  int        fEventNumber;
  bool       fIdentified;
  PhysicVar *fRecoEnergy_e;
  PhysicVar *fRecoEnergy_p;
  PhysicVec *fRecoPosition_e;
  PhysicVec *fRecoPosition_p;
  PhysicVec *fRecoDirection_scatter;
  vector <PhysicVec*> *fRecoClusterPositions;
  vector <PhysicVar*> *fRecoClusterEnergies;
  
  double     fScatDimX;
  double     fScatDimY;
  double     fScatDimZ;
  double     fAbsDimX;
  double     fAbsDimY;
  double     fAbsDimZ;
  TVector3   fScatPosition;
  TVector3   fAbsPosition;
  
  bool AccessTree(TString name);
  bool AccessSetup(void);
  
  ClassDef(InputReaderGeant,0)
};

#endif
