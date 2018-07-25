#ifndef __CCMLEM_H_
#define __CCMLEM_H_ 1 
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include "ComptonCone.hh"
#include "CCReconstruction.hh"
#include "TGraph.h"
#include "TStopwatch.h"
class CCMLEM : public TObject{
  
public:
  CCMLEM();
  CCMLEM(TString path);
  ~CCMLEM();
  Bool_t Iterate(Int_t nstop, Int_t iter);
  Bool_t Reconstruct(void);

  //Bool_t SaveToFile(TGraph *h);
  Bool_t SaveToFile(TObject *ob);
  Bool_t Config(TString path);
  Bool_t DrawHisto(void);
  Double_t Smear(double val, double sigma);
  bool SetInputReader(void);
  
  void SetName(TString name){ fName = name; };
  void SetInputName(TString inputName){ fInputName = inputName; };
  void SetIter(Int_t iter){ fIter = iter; };
  void Print(void);
  
  Int_t AddIsectionPoint(TString dir, Double_t x, Double_t y, Double_t z);
  
private:
  TString   fInputName;
  TString   fName;
  Bool_t    fVerbose;
  Bool_t    fFreshOutput;
  Bool_t    fSmear;
  Int_t     fIter;
  TFile     *fOutputFile;
  TH2F      *fImage[100];
  TGraph    *fGraph;
  TGraph    *g;
  Double_t  fDimZ;
  Double_t  fDimY;
  Int_t     fNbinsZ;
  Int_t     fNbinsY;
  Int_t     fNIpoints;
  Int_t     fpoints;
  Int_t     fStart;
  Int_t     fStop;
  Double_t  fPixelSizeZ;
  Double_t  fPixelSizeY;
  Double_t  fXofRecoPlane;
  Double_t  fYofRecoPlane;
  Double_t  fZofRecoPlane;
  Double_t  fSigmaE;
  Double_t  fResolutionX;
  Double_t  fResolutionY;
  Double_t  fResolutionZ;
  TClonesArray*  fArray;
  TClonesArray* fSM;
  InputReader *fReader;
  
  ClassDef(CCMLEM,0)
};

#endif
