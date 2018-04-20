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
//#include "DetPlane.hh"
#include "CCReconstruction.hh"
#include "TGraph.h"
#include "IsectionPoint.hh"

class CCMLEM : public TObject{
  
public:
  CCMLEM(TString inputName, TString name, Int_t iter, Bool_t verbose, Double_t dimZ, Double_t dimY);
  ~CCMLEM();
 
  //TVector3 ConnectPoints(TVector3 *point1, TVector3 *point2);
  Bool_t Reconstruct(Int_t iStart, Int_t iStop);
  Bool_t SaveHistogram(TH2F *h);
 // Bool_t Sort(Int_t iFirst, Int_t iEnd);
  Bool_t SaveToFile(TGraph *h);
  //void SetClass(const char *classname, Int_t s = 1000);
  void SetName(TString name){ fName = name; }; 
  void SetInputName(TString inputName){ fInputName = inputName; };
  void SetIter(Int_t iter){ fIter = iter; };
  //TH2F* GetImage(void) { return fImage; };
  
private:
  TString   fInputName;
  TString   fName;
  Int_t     fNev;
  Bool_t    fVerbose;
  Int_t     fIter;
  TFile     *fFile;
  TTree     *fTree;
  TH2F      *fImage;
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
  Double_t  fXofRecoPlane;
  TClonesArray*  fArray;
  //DetPlane  fScatterer, fAbsorber;
  
  ClassDef(CCMLEM,0)
};

#endif
