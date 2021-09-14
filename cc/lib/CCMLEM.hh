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
#include <TMatrixT.h>
#include "Coordinates.hh"
#include <vector>
///To be acquainted with Class for image reconstruction from events simulated in class CCSimulations
///and using MLEM method based on ComptonCone class objects. More details about this class
/// is available on wiki [LINK] (http://bragg.if.uj.edu.pl/gccbwiki/index.php/File:MK_20180513_Image_Reconstruction_Analysis_For_CC_Toy_Model.pptx)
/// Moreover, it can be used for reconstruction of Geant4 and machine learning outputs
class CCMLEM : public TObject {

public:
  CCMLEM();
  CCMLEM(TString path);
  ~CCMLEM();

  Bool_t Iterate(Int_t nstop, Int_t iter);
  Bool_t Reconstruct();
  Bool_t GetSigmaError(void);
  Int_t AddIsectionPoint(TString dir, Double_t x, Double_t y, Double_t z);
  //Bool_t Sensitivity(void);
  TH2 *SmoothGauss(TH2 *hin, double sigma);
  Double_t SmearGaus(double val, double sigma);
  Double_t SmearBox(double x, double resolution);
  Double_t GetSigmaE(double energy);
  static Double_t FitProjection(Double_t *x,Double_t *par);
  static Double_t FitProjection1(Double_t *x,Double_t *par);
  Bool_t ReadConfig(TString path);
  Bool_t SetInputReader(void);
  Bool_t DrawCanvas(void);
  Bool_t DrawHisto(void);
  Bool_t SaveToFile(TObject* ob);
  void Print(void);
  void Clear(void);
  /** Calculate H matrix and save it in file */
  void SmatrixToFile(const TString& filename);
  Bool_t DrawREGraph(void);

private:
  TString fInputName;       ///< Path to the file with simulation data
  Double_t fXofRecoPlane;       ///< x-component of image plane coordinate
  Double_t fYofRecoPlane;       ///< y-component of image plane coordinate
  Double_t fZofRecoPlane;       ///< z-component of image plane coordinate
  Double_t fScatthick_x;
  Double_t fScatthick_y;
  Double_t fScatthick_z;
  
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
  Bool_t fVerbose;      ///< Verbose level for print-outs on screen
  Int_t fEvent[1000000];
  Int_t fSubFirst[1000000];
  Int_t fSubSecond[1000000];
  Int_t fNIpoints;      ///< Numbers of intersection points for each Compton cone
  Int_t fPoints;        ///< Numbers of intersection points for all Compton cones
  Double_t fPixelSizeZ;     ///< Size of pixel in z-axis direction
  Double_t fPixelSizeY;     ///< Size of pixel in y-axis direction
  Double_t fPixelSizeX;     ///< Size of pixel in x-axis direction
  Double_t fSigma[250];     ///< Relative sigma value to compare different iterations
  Double_t sigma[250];
  Double_t fDenominator[10000000];
  std::vector< std::pair<Int_t, Double_t>> fRErr;
  TH2F* fSensitivity;
  TH1D* fHisto;     ///< Histogram containing energy resolution obtained by Geant4
  TFile* fOutputFile;       ///< ROOT file containing reconstruction results 

  TH2F* fImage[250];        ///< Reconstructed image histogram
  
  TH2F* fSH[250];
  TH2F* fSmatrix;
  TTree* fTree;
  TTree* fTree1;
  TH2F* fSenHisto[250];
  //TH1F* fAngDiff; 
  TH1D* fProX[250];
  TH1D* fProZ[250];
  TH1D* fProY[250];
  TClonesArray* fArray;     ///< Array of information for intersection of Compton cone with image plane
  TClonesArray* fSM;        ///< Array of information for intersection of all Compton cones with image plane
  
  InputReader* fReader;     ///< InputReader to read different given input simulation files
  //TGraph* fGraph;
  
  ClassDef(CCMLEM, 0)
};

#endif
