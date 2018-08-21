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
#include "TGraph.h"
#include "TStopwatch.h"
#include "InputReader.hh"
#include "InputReaderSimple.hh"
#include "InputReaderGeant.hh"

class CCMLEM : public TObject{
  
public:
  CCMLEM();
  CCMLEM(TString path);
  ~CCMLEM();
  
  Bool_t   Iterate(Int_t nstop, Int_t iter);
  Bool_t   Reconstruct(void);
  Int_t    AddIsectionPoint(TString dir, Double_t x, Double_t y, Double_t z);
  Double_t Smear(double val, double sigma);
  Bool_t   ReadConfig(TString path);
  Bool_t   SetInputReader(void);
  Bool_t   DrawHisto(void);
  Bool_t   SaveToFile(TObject *ob);
  void     Print(void);
  void     Clear(void);
  
private:
  TString   fInputName;
  Double_t  fXofRecoPlane;
  Double_t  fYofRecoPlane;
  Double_t  fZofRecoPlane;
  Double_t  fDimZ;
  Double_t  fDimY;
  Int_t     fNbinsZ;
  Int_t     fNbinsY;
  Bool_t    fSmear;
  Double_t  fResolutionX;
  Double_t  fResolutionY;
  Double_t  fResolutionZ;
  Double_t  fSigmaE;
  Int_t     fIter;
  Bool_t    fFreshOutput;
  Int_t     fStart;
  Int_t     fStop;
  Bool_t    fVerbose;
  
  Int_t         fNIpoints;
  Int_t         fPoints;
  Double_t      fPixelSizeZ;
  Double_t      fPixelSizeY;
  TFile        *fOutputFile;
  TH2F         *fImage[100];
  TClonesArray *fArray;
  TClonesArray *fSM;
  InputReader  *fReader;
  TGraph       *fGraph;
  
  ClassDef(CCMLEM,0)
};

#endif
