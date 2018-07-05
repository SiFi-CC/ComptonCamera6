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
#include "IsectionPoint.hh"
///Class for image reconstruction from events simulated in class CCSimulations
///and using MLEM method based on ComptonCone class objects. More details about this class
/// is available on wiki [LINK] (http://bragg.if.uj.edu.pl/gccbwiki/index.php/File:MK_20180513_Image_Reconstruction_Analysis_For_CC_Toy_Model.pptx)

class CCMLEM : public TObject{
  
public:
  CCMLEM(TString inputName, TString name, Int_t iter, Bool_t verbose, 
	 Double_t dimZ, Double_t dimY, Int_t nbinsz, Int_t nbinsy);
  ~CCMLEM();
 
  
  Bool_t Reconstruct(Int_t iStart, Int_t iStop);
  Bool_t SaveHistogram(TH2F *h);
 
  Bool_t SaveToFile(TGraph *h);
  Bool_t SaveToFile(TObject *h);
  ///Sets name of the CCReconstruction object.
  void SetName(TString name){ fName = name; };
  ///Sets path to the file containing simulation results.
  void SetInputName(TString inputName){ fInputName = inputName; };
  ///Sets number or iterations for MLEM (not used yet).
  void SetIter(Int_t iter){ fIter = iter; };
 
  Int_t AddIsectionPoint(TString dir, Double_t x, Double_t y, Double_t z);
  
private:
  TString   fInputName;		 ///< Path to the file with simulation data         
  TString   fName;		 ///< Name of the CCMLEM object
  Int_t     fNev;		 ///< Number of events to be reconstructed
  Bool_t    fVerbose;		 ///< Verbose level for print-outs on screen
  Int_t     fIter;		 ///< Number of iterations for MLEM (not used yet)
  TFile     *fFile;		 ///< ROOT file containing simulation data
  TTree     *fTree;		 ///< ROOT tree containing simulation data
  TH2F      *fImage;		 ///< Reconstructed image histogram
  TGraph    *fGraph;		 ///< Non-pixelized reconstructed image
  TGraph    *g;			 ///< Non-pixelized reconstructed image graph
  TVector3  *fPoint0;		 ///< Coordinates of the gamma origin
  TVector3  *fPoint1;		 ///< Coordinates of the interaction point in scatterer
  TVector3  *fPoint2;		 ///< Coordinates of the interaction point in absorber
  TVector3  *fVersor1;		 ///< Direction of gamma ray coming from the source
  TVector3  *fVersor2;		 ///< Direction of the scattered gamma ray
  Double_t  fEnergy0;		 ///< Initial gamma energy
  Double_t  fEnergy1;		 ///< Energy loss due to Compton scattering
  Double_t  fEnergy2;		 ///< Energy after Compton scattering
  Double_t  fDimZ;		 ///< Size of image plane in direction z-axis
  Double_t  fDimY;		 ///< Size of image plane in direction y-axis
  Int_t fNbinsZ;		 ///< Numbers of bins of image plane in direction of z-axis
  Int_t fNbinsY;		 ///< Numbers of bins of image plane in direction of y-axis
  Int_t fNIpoints;		 ///< Numbers of intersection points for each ComptonCone
  Double_t  fPixelSizeZ;	 ///< Size of pixel in z-axis direction 
  Double_t  fPixelSizeY;	 ///< Size of pixel in y-axis direction 
  Double_t  fXofRecoPlane;	 ///< x-component of image plane coordinate
  TClonesArray*  fArray;	 ///< Array of information for intersection of ComptonCone with image plane 
  
  
  ClassDef(CCMLEM,0)
};

#endif
