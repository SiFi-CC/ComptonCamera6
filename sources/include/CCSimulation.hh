#ifndef __CCSimulation_H_
#define __CCSimulation_H_ 1
#include "TObject.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "PhysicsBase.hh"
#include "Track.hh"
#include "TGeoManager.h"

class CCSimulation : public TObject{
  
public:
  
  CCSimulation();
  CCSimulation(TString name, Bool_t verbose);
  ~CCSimulation();
  
  void   BuildSetup(Double_t scatDist, Double_t scatZ, Double_t scatY,
                    Double_t absDist, Double_t absZ, Double_t absY);
  Bool_t GenerateRay(void);
  Bool_t ProcessEvent(void);
  void   Loop(Int_t nev);
  void   Clear(void);
  void   SaveGeometryTxt(void);
  void   BuildTGeometry(void);
  void   Print(void);
  void   SetGenVersion(Int_t gen){ fGenVersion = gen; };
  void   SetName(TString name){ fName = name; };
  const char*  GetName() const { return fName.Data(); };
  
private:
  
  TString     fName;
  Bool_t      fVerbose;
  Int_t       fGenVersion;
  Int_t       fNev;
  TTree       *fTree;
  TFile       *fFile;
  PhysicsBase fPhysics;
  DetPlane    fScatterer, fAbsorber;
  Track       fTrack1;
  Track       *fTrack2;
  TVector3    fPoint0, fPoint1, fPoint2;
  TVector3    fVersor1, fVersor2;
  Double_t    fEnergy0, fEnergy1, fEnergy2;
  Double_t    fYgap, fZgap;
  Double_t    fRadius;
  
  TH2F *hSource;
  TH2F *hScat;
  TH2F *hAbs;
  TH1F *hEnergy;
  
  ClassDef(CCSimulation,1)
};

#endif