#ifndef __InputReaderSimple_H_
#define __InputReaderSimple_H_ 1
#include "InputReader.hh"
#include "TF1.h"
#include "TH1D.h"
#include "TRandom.h"
#include "TString.h"
#include "TSystem.h"
#include "TVector3.h"
#include <iostream>
#include <fstream>

using namespace std;

/// Class for accessing data from simple simulation performed by KR.
/// This is class derived from InputReader class. It opens requested
/// ROOT file containing tree with simulation results and via set of
/// getter function passes information to reconstruction classes, i.e.
/// CCREconstruction and CCMLEM.

class InputReaderSimple : public InputReader
{

public:
    InputReaderSimple();
    InputReaderSimple(TString path);
    ~InputReaderSimple();

    void Clear(void);
    bool LoadEvent(int i);
    TVector3* GetPositionPrimary(void);
    TVector3* GetPositionScattering(void);
    TVector3* GetPositionAbsorption(void);
    TVector3* GetGammaDirPrimary(void);
    TVector3* GetGammaDirScattered(void);
    double GetEnergyPrimary(void);
    double GetEnergyLoss(void);
    double GetEnergyScattered(void);

    void SetSmearing(bool smear, Double_t posScatX, Double_t posAbsX, Double_t posScatY, Double_t posAbsY, Double_t posScatZ, Double_t posAbsZ);
    
    void SelectEvents(void);
    list<int> GetSelectedEvents(void);
    
    TString fInputFile; ///< Path to the root file
    
    void ReadGeometry(void); ///< Reads scatterer's and absorber's width and height
    /* Scatterer's and absorber's width and height read from geometry file */
    Double_t fScatWidth;
    Double_t fScatHeight;
    Double_t fAbsWidth;
    Double_t fAbsHeight;
    void PrintEvent();

private:
    // LOADING FROM THE TREE
    TVector3* fPoint0;        ///< Coordinates of the gamma source
    TVector3* fPoint1;        ///< Coordinates of the interaction in the scatterer
                              ///< (Comptin scattering)
    TVector3* fPoint2;        ///< Coordinates of the interaction in the absorber (absorption)
    TVector3* fVersor1;       ///< Direction of the primary gamma
    TVector3* fVersor2;       ///< Direction of the scattered gamma
    double fEnergy0;          ///< Energy of the primary gamma [MeV]
    double fEnergy1;          ///< Energy deposited in the scatterer [MeV]
    double fEnergy2;          ///< Energy of the scattered gamma [MeV]
                              // PARAMETERS PASSED ON
    TVector3* fPositionScat;  ///< Position of interaction in scatterer
    TVector3* fPositionAbs;   ///< Position of interaction in absorber
    TVector3* fDirectionScat; ///< Direction of scattered gamma
    double fPrimaryEnergy;
    double fEnergyLoss;
    double fEnergyScattered;
    // INTERN PARAMETES FOR SMEARING
    Bool_t fSmear;
    TH1D* fHisto; ///< Histogram containing energy resolution obtained by Geant4

    Double_t fScatResolutionX; ///< Position resolution in direction x-axis of the scatterer
    Double_t fScatResolutionY; ///< Position resolution in direction y-axis of the scatterer
    Double_t fScatResolutionZ; ///< Position resolution in direction z-axis of the scatterer
    Double_t fAbsResolutionX; ///< Position resolution in direction x-axis of the absorber
    Double_t fAbsResolutionY; ///< Position resolution in direction y-axis of the absorber
    Double_t fAbsResolutionZ; ///< Position resolution in direction z-axis of the absorber

    Double_t SmearGaus(double val, double sigma);
    Double_t SmearBox(double x, double resolution);
    Double_t GetSigmaEScat(double energy);
    Double_t GetSigmaEAbs(double energy);

    bool AccessTree(TString name);
    
    list<int> fSelectedEvents;
    bool SelectSingleEvent(void);

    ClassDef(InputReaderSimple, 0)
};
//------------------------------------------------------------------
inline TVector3* InputReaderSimple::GetPositionPrimary(void) { return fPoint0; }
//------------------------------------------------------------------
inline TVector3* InputReaderSimple::GetPositionScattering(void) { return fPositionScat; }
//------------------------------------------------------------------
inline TVector3* InputReaderSimple::GetPositionAbsorption(void) { return fPositionAbs; }
//------------------------------------------------------------------
inline TVector3* InputReaderSimple::GetGammaDirPrimary(void) { return fVersor1; }
//------------------------------------------------------------------
inline TVector3* InputReaderSimple::GetGammaDirScattered(void) { return fDirectionScat; }
//------------------------------------------------------------------
inline double InputReaderSimple::GetEnergyPrimary(void) { return fPrimaryEnergy; }
//------------------------------------------------------------------
inline double InputReaderSimple::GetEnergyLoss(void) { return fEnergyLoss; }
//------------------------------------------------------------------
inline double InputReaderSimple::GetEnergyScattered(void) { return fEnergyScattered; }
//------------------------------------------------------------------
inline list<int> InputReaderSimple::GetSelectedEvents(void){return fSelectedEvents; }

#endif
