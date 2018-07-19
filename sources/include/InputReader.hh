#ifndef __InputReader_H_
#define __InputReader_H_ 1
#include "TObject.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TVector3.h"
#include <iostream>

using namespace std;

class InputReader : public TObject{
  
public:
   InputReader();
   InputReader(TString path);
  ~InputReader();
  
  bool     LoadEvent(int i);
  void     Clear(void);
  void     Print(void);
  TVector3 virtual *GetSourcePosition(void);
  TVector3 virtual *GetScatPosition(void);
  TVector3 virtual *GetAbsPosition(void);
  TVector3 virtual *GetPrimaryGammaDir(void);
  TVector3 virtual *GetScatGammaDir(void);
  double   virtual GetEnSource(void);
  double   virtual GetEnScat(void);
  double   virtual GetEnAbs(void);
  
protected:
  TFile   *fFile;
  TTree   *fTree;
  
  bool    virtual AccessTree(TString name);
  bool    SetInputFile(TString path);
  
  ClassDef(InputReader,0)
};

#endif
