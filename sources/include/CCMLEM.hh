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
  ~CCMLEM();
  Bool_t Iterate(Int_t nstart, Int_t nstop, Int_t iter);
  Bool_t Reconstruct(Int_t iStart, Int_t iStop);

  Bool_t SaveToFile(TGraph *h);
  Bool_t SaveToFile(TObject *h);
  Bool_t Config(void);
  Bool_t Smear();
  Bool_t Freshoutput();
  Bool_t Drawhisto(void);
  ///Sets name of the CCReconstruction object.
  void SetName(TString name){ fName = name; };
  ///Sets path to the file containing simulation results.
  void SetInputName(TString inputName){ fInputName = inputName; };
  void SetIter(Int_t iter){ fIter = iter; };
  
  Int_t AddIsectionPoint(TString dir, Double_t x, Double_t y, Double_t z);
  
private:
  TString   fInputName;
  TString   fName;
  TString   fPath;
  Int_t     fNev;
  Bool_t    fVerbose;
  Int_t     fIter;
  TFile     *fFile;
  TFile     *file;
  TTree     *fTree;
  TH2F      *fImage[100];
  TGraph    *fGraph;
  TGraph    *g;
  TVector3  *fPoint0;
  TVector3  *fPoint1;
  TVector3  *fPoint2;
  TVector3  *fVersor1;
  TVector3  *fVersor2;
  Double_t  fEnergy0;
  Double_t  fEnergy1;
  Double_t  fEnergy2;
  Double_t  fDimZ;
  Double_t  fDimY;
  Int_t     fNbinsZ;
  Int_t     fNbinsY;
  Int_t     fNIpoints;
  Int_t     fpoints;
  Double_t  fPixelSizeZ;
  Double_t  fPixelSizeY;
  Double_t  fXofRecoPlane;
  Double_t  fYofRecoPlane;
  Double_t  fZofRecoPlane;
  Double_t  fSigmaE;
  Double_t  fSigmaX;
  Double_t  fSigmaY;
  Double_t  fSigmaZ;
  TClonesArray*  fArray;
  TClonesArray* fSM;
  
  
  ClassDef(CCMLEM,0)
};

#endif
