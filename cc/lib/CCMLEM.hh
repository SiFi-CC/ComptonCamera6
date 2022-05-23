#ifndef __CCMLEM_H_
#define __CCMLEM_H_ 1
#include "ComptonCone.hh"
#include "Coordinates.hh"
#include "InputReader.hh"
#include "InputReaderEI.hh"
#include "InputReaderGeant.hh"
#include "InputReaderNN.hh"
#include "InputReaderSimple.hh"
#include "InputReaderPMI.hh"
#include "TFile.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TObject.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVector3.h"
#include "Track.hh"
#include <TMatrixT.h>
#include <vector>
#include <map>
#include <list>
///To be acquainted with Class for image reconstruction from events simulated in class CCSimulations
///and using MLEM method based on ComptonCone class objects. More details about this class
/// is available on wiki [LINK] (http://bragg.if.uj.edu.pl/gccbwiki/index.php/File:MK_20180513_Image_Reconstruction_Analysis_For_CC_Toy_Model.pptx)
/// Moreover, it can be used for reconstruction of Geant4 and machine learning outputs
/// The usage is mostly commmented in the files now. 
class CCMLEM : public TObject {

public:
  CCMLEM();
  CCMLEM(TString path);
  ~CCMLEM();
  Bool_t Iterate(Int_t nstop, Int_t iter);
  Bool_t Reconstruct();

  Bool_t DetermineConvergence(Int_t iter);
  Bool_t GetSigmaError(void);
  Int_t AddIsectionPoint(TString dir, Double_t x, Double_t y, Double_t z);
  //Bool_t Sensitivity(void);
  TH2 *SmoothGauss(TH2 *hin, double sigma);
  static Double_t FitProjection(Double_t *x,Double_t *par);
  static Double_t FitProjection1(Double_t *x,Double_t *par);
  Bool_t ReadConfig(TString path);
  Bool_t SetInputReader(void);
  void DrawAtConvergence(int iter);
  Double_t DrawSubSetResults(int iter, int subsetnumber);
  void DrawREGraph(int iter);
  Bool_t DrawCanvas(void);
  Bool_t DrawHisto(void);
  Bool_t SaveToFile(TObject* ob);
  void Print(void);
  void Clear(void);
  /** Calculate H matrix and save it in file */
  void SmatrixToFile(const TString& filename);
  void DrawAllIterations(void);

private:
  TString fInputName;       ///< Path to the file with simulation data
  TString fOutputName;       ///< Path to the file with the results 
//  Double_t fScatthick_x;
//  Double_t fScatthick_y;
//  Double_t fScatthick_z;
//  Double_t fAbsthick_x;
//  Double_t fAbsthick_y;
//  Double_t fAbsthick_z;
//  TVector3* fScatposition; 
//  TVector3* fAbsposition;

    ///SIMPLE SIMULATION INPUT
    Bool_t fSmear; ///< Smear flag for smearing energy and position only used for simple Simulation results
    Double_t fScatResolutionX; ///< Position resolution of the scatterer in direction x-axis only used for simple Simulation results
    Double_t fScatResolutionY; ///< Position resolution of the scatterer in direction y-axis only used for simple Simulation results
    Double_t fScatResolutionZ; ///< Position resolution of the scatterer in direction z-axis only used for simple Simulation results
    Double_t fAbsResolutionX; ///< Position resolution of the absorber in direction x-axis only used for simple Simulation results
    Double_t fAbsResolutionY; ///< Position resolution of the absorber in direction y-axis only used for simple Simulation results
    Double_t fAbsResolutionZ; ///< Position resolution of the absorber in direction z-axis only used for simple Simulation results
    Double_t fP0; ///< Fitting parameter used for the energy smearing in simple simulation
    Double_t fP1; ///< Fitting parameter used for the energy smearing in simple simulation
    Double_t fP2; ///< Fitting parameter used for the energy smearing in simple simulation

    // Convergence
    Bool_t fSimpleConvergence;      ///< Flag or function or pixel convergence
    Double_t fROIX;                 ///< x-component of ROI for convergence criterium
    Double_t fROIY;                 ///< y-component of ROI for convergence criterium
    Double_t fROIZ;                 ///< z-component of ROI for convergence criterium
    Double_t fConvergenceCriterium; ///< value that leads to stop of iterations
    Double_t fSigma[251];           ///< Relative sigma values between adjustent iterations

///GEANT4 SIMULATION INPUT
  Bool_t fLoadReal; ///< Flag if Real or Reco data is loaded from the input
  Int_t fCorrectIdentified; ///< Indicates if only correct identified events from the Reco are loaded

//INTERNAL VARIABLES
  Int_t fIter;      ///< Number of iterations for MLEM
  Bool_t fFreshOutput;      ///< FreshOutput flag to recreate or update output file
  Bool_t fSaveAll;      ///<SaveAll flag to save the intermediate steps 
  Int_t fStart;     ///< first event number
  Int_t fStop;      ///< last event number
  Bool_t fVerbose;      ///< Verbose level for print-outs on screen
  Int_t fNIpoints;      ///< Numbers of intersection points for each Compton cone
  Int_t fPoints;        ///< Numbers of intersection points for all Compton cones
  //Image parametes 
  Double_t fDimZ;       ///< Size of image plane in direction z-axis
  Double_t fDimY;       ///< Size of image plane in direction y-axis
  Double_t fDimX;       ///< Size of image plane in direction x-axis
  Int_t fNbinsZ;        ///< Numbers of bins of image plane in direction of z-axis
  Int_t fNbinsY;        ///< Numbers of bins of image plane in direction of y-axis
  Int_t fNbinsX;        ///< Numbers of bins of image plane in direction of x-axis
  Double_t fPixelSizeZ;     ///< Size of pixel in z-axis direction
  Double_t fPixelSizeY;     ///< Size of pixel in y-axis direction
  Double_t fPixelSizeX;     ///< Size of pixel in x-axis direction

  //Reconstructionparameters
  Double_t fDenominator[10000000];
  TH2F* fSensitivity;
  TFile* fOutputFile;       ///< ROOT file containing reconstruction results 
  std::list<int> fEvents;
  TRandom3* fRandomGen;
  double fGausFilter; 
  TH2F* fImage[251];        ///< Reconstructed image histogram
  
  TH2F* fSH[251];
  TH2F* fSmatrix;

  list<Double_t> fBraggPeakPosition;     ///< determined Braggpeak position with a ERR function fit
  TH1F* fResolutionDistribution;

  TTree* fTree;
  TClonesArray* fArray;     ///< Array of information for intersection of Compton cone with image plane
  TClonesArray* fSM;        ///< Array of information for intersection of all Compton cones with image plane
  
  InputReader* fReader;     ///< InputReader to read different given input simulation files
  
  ClassDef(CCMLEM, 0)
};

#endif
