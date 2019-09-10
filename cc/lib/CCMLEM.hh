#ifndef __CCMLEM_H_
#define __CCMLEM_H_ 1
#include "ComptonCone.hh"
#include "InputReader.hh"
#include "InputReaderGeant.hh"
#include "InputReaderSimple.hh"
#include "InputReaderEI.hh"
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
#include "TStyle.h"
///Class for image reconstruction from events simulated in class CCSimulations
///and using MLEM method based on ComptonCone class objects. More details about this class
/// is available on wiki [LINK] (http://bragg.if.uj.edu.pl/gccbwiki/index.php/File:MK_20180513_Image_Reconstruction_Analysis_For_CC_Toy_Model.pptx)

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
  Double_t SmearBox(double x, double resolution);
  Double_t GetSigmaE(double energy);
  static Double_t FitProjection(Double_t *x,Double_t *par);
  Bool_t ReadConfig(TString path);
  Bool_t SetInputReader(void);
  Bool_t DrawHisto(void);
  Bool_t SaveToFile(TObject* ob);
  void Print(void);
  void Clear(void);

private:
  TString fInputName;       ///< Path to the file with simulation data
  Double_t fXofRecoPlane;       ///< x-component of image plane coordinate
  Double_t fYofRecoPlane;       ///< y-component of image plane coordinate
  Double_t fZofRecoPlane;       ///< z-component of image plane coordinate
  Double_t fDimZ;       ///< Size of image plane in direction z-axis
  Double_t fDimY;       ///< Size of image plane in direction y-axis
  Double_t fDimX;       ///< Size of image plane in direction x-axis
  Int_t fNbinsZ;        ///< Numbers of bins of image plane in direction of z-axis
  Int_t fNbinsY;        ///< Numbers of bins of image plane in direction of y-axis
  Int_t fNbinsX;        ///< Numbers of bins of image plane in direction of x-axis
  Bool_t fSmear;        ///< Smear flag for applying to camera performance 
  Double_t fResolutionX;        ///< Position resolution in direction x-axis
  Double_t fResolutionY;        ///< Position resolution in direction y-axis
  Double_t fResolutionZ;        ///< Position resolution in direction z-axis
  Double_t fP0;     ///< Fitting parameter
  Double_t fP1;     ///< Fitting parameter
  Double_t fP2;     ///< Fitting parameter
  Int_t fIter;      ///< Number of iterations for MLEM
  Bool_t fFreshOutput;      ///< FreshOutput flag to recreate or update output file
  Int_t fStart;     ///< first event number
  Int_t fStop;      ///< last event number
  Int_t fGeantFilter;	///< Select which events to use from Geant4 simulation @see InputReaderGeant.fFilter
  Bool_t fVerbose;      ///< Verbose level for print-outs on screen

  Int_t fNIpoints;      ///< Numbers of intersection points for each Compton cone
  Int_t fPoints;        ///< Numbers of intersection points for all Compton cones
  Double_t fPixelSizeZ;     ///< Size of pixel in z-axis direction
  Double_t fPixelSizeY;     ///< Size of pixel in y-axis direction
  Double_t fPixelSizeX;     ///< Size of pixel in x-axis direction
  Double_t fSigma[100];     ///< Relative sigma value to compare different iterations
  TH1D* fHisto;     ///< Histogram containing energy resolution obtained by Geant4
  TFile* fOutputFile;       ///< ROOT file containing reconstruction results 
  TH2F* fImage[150];        ///< Reconstructed image histogram
  //TH1F* fAngDiff;       
  TClonesArray* fArray;     ///< Array of information for intersection of Compton cone with image plane
  TClonesArray* fSM;        ///< Array of information for intersection of all Compton cones with image plane
  InputReader* fReader;     ///< InputReader to read different given input simulation files
  //TGraph* fGraph;

  ClassDef(CCMLEM, 0)
};

#endif
