#ifndef __InputReaderGeant_H_
#define __InputReaderGeant_H_ 1
#include "InputReader.hh"
#include "DR_GenerallStructs.hh"
#include <iostream>

using namespace std;

class InputReaderGeant : public InputReader{
  
public:
  InputReaderGeant();
  InputReaderGeant(TString path);
  ~InputReaderGeant();
  
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
  
  bool AccessTree(TString name);
  
  ClassDef(InputReaderGeant,0)
};

#endif
