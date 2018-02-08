#ifndef __CCReconstruction_H_
#define __CCReconstruction_H_ 1 
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include "ComptonCone.hh"
#include "DetPlane.hh"

class CCReconstruction : public TObject{
  
public:
  CCReconstruction(TString inputName, TString name, Int_t iter, Bool_t verbose);
  ~CCReconstruction();
  
  TVector3 ConnectPoints(TVector3 point1, TVector3 point2);
  Double_t CalculateTheta(Double_t e1, Double_t e2);
  void RebuildSetupTxt(void);
  ComptonCone* ReconstructCone(Int_t i);
  Bool_t ReconstructImage(Int_t iStart, Int_t iStop);
  void   Clear(void);
  Bool_t SaveHistogram(TH1F *h);
  Bool_t SaveHistogram(TH2F *h);
  void Print(void);
  void SetName(TString name){ fName = name; }; 
  void SetInputName(TString inputName){ fInputName = inputName; };
  void SetIter(Int_t iter){ fIter = iter; };
  
private:
  TString   fInputName;
  TString   fName;
  Bool_t    fVerbose;
  Int_t     fNev;
  Int_t     fIter;
  Int_t     fGenVersion;
  TFile     *fFile;
  TTree     *fTree;
  TH2F      *fImage;
  TH1F      *fNpixels;
  TVector3  *fPoint0;
  TVector3  *fPoint1;
  TVector3  *fPoint2;
  TVector3  *fVersor1;
  TVector3  *fVersor2;
  Double_t  fEnergy0;
  Double_t  fEnergy1;
  Double_t  fEnergy2;
  DetPlane  fScatterer, fAbsorber;
  
  ClassDef(CCReconstruction,0)
};

#endif
