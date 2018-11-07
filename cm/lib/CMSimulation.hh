#ifndef __CMSimulation_H_
#define __CMSimulation_H_ 1
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "Mask.hh"
#include "Track.hh"
#include "TGeoManager.h"

class CMSimulation : public TObject{

public:
  CMSimulation();
  CMSimulation(TString name, Int_t vLevel);
  ~CMSimulation();
  void   BuildSetup(double maskdist,double maskZsize, double maskYsize,
		    double detdist, double detZsize, double detYsize);
  void   BuildSetup(DetPlane plane, Mask mask){fDetPlane=plane; fMask=mask;};
  void   Print(void);
  void   Clear(void);
  void   Loop(Int_t nev);
  Bool_t GenerateRay(void);
  Bool_t ProcessEvent(void);
  void   SaveGeometryTxt(void);
  void   BuildTGeometry(void);
  void   SetupSpectra(void);
  void   ClearSpectra(void);
  void   SetName(TString name){fName = name;};
  Bool_t  SetPattern(TH2F* h);
  void   SetGenVersion(Int_t i){fGenVersion = i;};
  void   SetSourcePosition(Double_t x, Double_t y, Double_t z){fSourcePos.SetXYZ(x,y,z); };  
  TVector3   GetSourcePosition(void){return fSourcePos; };
  const char*  GetName() const { return fName.Data(); };
  TH2F*  GetImage() { return hYZdetected;};

private:
  TVector3 fPoint0, fPoint1, fPoint2, fDir;
  Int_t    fOpaque;
  TFile*   fFile;
  TTree*   fTree;
  TString  fName;
  Mask     fMask;
  DetPlane fDetPlane;
  Track    fTrack;
  Int_t    fVerbLevel;
  Int_t    fNev;
  Int_t    fGenVersion;
  TVector3 fSourcePos;

  TH2F* hYZ;
  TH2F* hYZdetected;
  TH1F* hCosTheta;
  TH1F* hPhi;
  
  ClassDef(CMSimulation,1)
};

#endif
