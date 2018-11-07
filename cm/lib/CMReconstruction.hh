#ifndef __CMReconstruction_H_
#define __CMReconstruction_H_ 1
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "Mask.hh"
#include "Track.hh"
#include <vector>

class CMReconstruction : public TObject{

public:
  CMReconstruction();
  CMReconstruction(TString name, Int_t vLevel);
  ~CMReconstruction();
  void   RebuildSetupTxt(void);
  Bool_t MLEMIterate(Int_t ni);
  void   Print(void);
  void   Write(void);
  void   SetupSpectra(void);
  void   SetName(TString name){fName = name;};
  const char*  GetName() const { return fName.Data(); };

private:
  Bool_t SingleToDoubleIdx(TString which, int i, int &binx, int& biny);
  Bool_t DoubleToSingleIdx(TString which, int binx, int biny, int &i);
  Double_t Image(Int_t i);
  Double_t Object(Int_t i);
  Double_t RecoObject(Int_t i);
  Double_t H(Int_t i, Int_t j);
  Double_t Hprime(Int_t i, Int_t j);
  Bool_t CalculateS(void);
  Bool_t FillHMatrix(void);
  Bool_t SingleIteration(void);

private:
  TString   fInputName;
  TString   fName;
  Bool_t    fVerbose;
  TFile     *fFileIn;
  TFile     *fFileOut;
  TH2F      *fImage;
  TH2F      *fRecoObject[100];
  TH2F      *fObject;
  TH2D      *fHmatrix;
  DetPlane  fDetPlane;
  Mask      fMask;

  Int_t     fNvoxelsI;
  Int_t     fNvoxelsO;
  Int_t     fNbinsxI;
  Int_t     fNbinsyI;
  Int_t     fNbinsxO;
  Int_t     fNbinsyO;

  Int_t fNiter;
  Int_t fThisIter;
  //Double_t* S;
  std::vector <Double_t> S;
  
  ClassDef(CMReconstruction,1)
};

#endif
