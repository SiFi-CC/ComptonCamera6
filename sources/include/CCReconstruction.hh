#ifndef __CCReconstruction_H_
#define __CCReconstruction_H_ 1 
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"
#include "TFile.h"
#include "TMath.h"
#include "TH2F.h"
#include "TH1F.h"
#include "ComptonCone.hh"
#include "DetPlane.hh"
#include "InputReader.hh"
#include "InputReaderSimple.hh"
#include "InputReaderGeant.hh"

///Class for image reconstruction from events simulated in class CCSimulations.
///This class usues simple back projection algorithm based on ComptonCone
///class objects. Detailed description of the algorithm can be found in 
///presentation by KR on wiki [LINK](http://bragg.if.uj.edu.pl/gccbwiki/images/0/0e/KR_20170222_CCandCarbonLine.pdf)

class CCReconstruction : public TObject{
  
public:
  CCReconstruction(TString inputName, TString name, Bool_t verbose);
  ~CCReconstruction();
  
  bool         SetInputReader(TString inputName);
  TVector3     ConnectPoints(TVector3 point1, TVector3 point2);
  Double_t     CalculateTheta(Double_t e1, Double_t e2);
  void         RebuildSetupTxt(void);
  ComptonCone* ReconstructCone(Int_t i);
  Bool_t       ReconstructImage(Int_t iStart, Int_t iStop);
  void         Clear(void);
  Bool_t       SaveHistogram(TH1F *h);
  Bool_t       SaveHistogram(TH2F *h);
  void         Print(void);
  
  ///Sets name of the CCReconstruction object. 
  void SetName(TString name)          { fName = name; }; 
  ///Sets path to the file containing simulation results.
  void SetInputName(TString inputName){ fInputName = inputName; };
  
private:
  TString   fInputName;		///< Path to the file with simulation data
  TString   fName;		///< Name of the CCReconstruction object
  Bool_t    fVerbose;		///< Verbose level for print-outs on screen
  Int_t     fGenVersion;	///< Version of gamma generator 
  TH2F      *fImage;		///< Reconstructed image histogramgram
  TH1F      *fNpixels;		///< Distribution of filled pixels in the reconstructed image
  TVector3  *fPoint0;		///< Coordinates of the gamma origin
  TVector3  *fPoint1;		///< Coordinates of the interaction point in scatterer
  TVector3  *fPoint2;		///< Coordinates of the interaction point in absorber
  TVector3  *fVersor1;		///< Direction of gamma ray coming from the source
  TVector3  *fVersor2;		///< Direction of the scattered gamma ray
  Double_t  fEnergy0;		///< Initial gamma energy
  Double_t  fEnergy1;		///< Energy loss due to Compton scattering
  Double_t  fEnergy2;		///< Energy after Compton scattering
  DetPlane  fScatterer;		///< Scatterer detection plane
  DetPlane  fAbsorber;		///< Absorber detection plane
  InputReader *fReader;
  
  ClassDef(CCReconstruction,0)
};

#endif
