#ifndef __CCMLEM_H_
#define __CCMLEM_H_ 1
#include "ComptonCone.hh"
#include "InputReader.hh"
#include "InputReaderGeant.hh"
#include "InputReaderSimple.hh"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TObject.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"
#include "Track.hh"

class CCMLEM : public TObject {

public:
  CCMLEM();
  CCMLEM(TString path);
  ~CCMLEM();

  Bool_t Iterate(Int_t nstop, Int_t iter);
  Bool_t Reconstruct(void);
  Bool_t GetSigmaError(void);
  Int_t AddIsectionPoint(TString dir, Double_t x, Double_t y, Double_t z);
  Double_t SmearGaus(double val, double sigma);
  Double_t SmearBoxX(double x);
  Double_t SmearBoxZ(double z);
  Double_t GetSigmaE(double energy);
  Bool_t ReadConfig(TString path);
  Bool_t SetInputReader(void);
  Bool_t DrawHisto(void);
  Bool_t SaveToFile(TObject* ob);
  void Print(void);
  void Clear(void);

private:
  TString fInputName;
  Double_t fDistScatter;
  Double_t fDistAbsorber;
  Double_t fXofRecoPlane;
  Double_t fYofRecoPlane;
  Double_t fZofRecoPlane;
  Double_t fDimZ;
  Double_t fDimY;
  Double_t fDimX;
  Int_t fNbinsZ;
  Int_t fNbinsY;
  Int_t fNbinsX;
  Bool_t fSmear;
  Double_t fResolutionX;
  Double_t fResolutionY;
  Double_t fResolutionZ;
  Int_t fIter;
  Bool_t fFreshOutput;
  Int_t fStart;
  Int_t fStop;
  Bool_t fVerbose;

  Int_t fNIpoints;
  Int_t fPoints;
  Double_t fPixelSizeZ;
  Double_t fPixelSizeY;
  Double_t fPixelSizeX;
  Double_t fSigma[100];
  TTree* fTree;
  TH1D* fHisto;
  TFile* fOutputFile;
  TH2F* fImage[150];
  TH1F* fAngDiff;
  TH1F* fGeoAng;
  TH1F* fScatAng;
  TClonesArray* fArray;
  TClonesArray* fSM;
  InputReader* fReader;
  TGraph* fGraph;

  ClassDef(CCMLEM, 0)
};

#endif
