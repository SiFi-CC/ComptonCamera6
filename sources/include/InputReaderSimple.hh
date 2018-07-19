#ifndef __InputReaderSimple_H_
#define __InputReaderSimple_H_ 1
#include "TString.h"
#include "TVector3.h"
#include "InputReader.hh"
#include <iostream>

using namespace std;

class InputReaderSimple : public InputReader{
  
public:
  InputReaderSimple();
  InputReaderSimple(TString path);
  ~InputReaderSimple();
  
  bool AccessTree(void);
  bool ReadEvent(int i, TVector3 &point1, TVector3 &point2, double &energy1, double &energy2);
  TVector3 *GetSourcePosition(void);
  TVector3 *GetScatPosition(void);
  TVector3 *GetAbsPosition(void);
  TVector3 *GetPrimaryGammaDir(void);
  TVector3 *GetScatGammaDir(void);
  double   GetEnSource(void);
  double   GetEnScat(void);
  double   GetEnAbs(void);
  
private:
  TVector3 *fPoint0;
  TVector3 *fPoint1;
  TVector3 *fPoint2;
  TVector3 *fVersor1;
  TVector3 *fVersor2;
  double   fEnergy0;
  double   fEnergy1;
  double   fEnergy2;
  
  ClassDef(InputReaderSimple,0)
};

#endif