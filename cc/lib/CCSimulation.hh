#ifndef __CCSimulation_H_
#define __CCSimulation_H_

#include "DetPlane.hh"

#include <TGeoManager.h>

class TFile;
class TH1F;
class TH2F;
class TTree;

/// Class which performs simulations of Compton Camera.
///
/// SIMULATIONS:
///
/// Several types of gamma rays sources (generators) are available:
/// 1. Isotropic point-like source located at (0,0,0)
/// 2. Beam-like uniform distribution along Z axis
/// 3. Two isotropic point-like sources located along Z axis and separated by
/// distance fZgap.
/// fZgap by default is set to 20 mm. It can be changed via SetGapZ(Double_t
/// gapZ) function.
/// 4. Two isotropic point-like sources located along Y axis and separated by
/// distance fYgap.
/// fYgap by default is set to 20mm. It can be changed via SetGapY(Double_t
/// gapY) function.
/// 5. Circular source of given radius fRadius. By defauls fRadius is 15 mm and
/// it can be changed
/// via SetRadius(Double_t radius) function. Center of the source is at (0,0,0).
///
/// RESULTS:
///
/// Simulated data are saved in the TTree and subsequently in the ROOT file. The
/// tree has the
/// following branches:
///- fPoint0  (TVector3) - position of the gamma source
///- fPoint1  (TVector3) - intersection point in the scatterer
///- fPoint2  (TVector3) - intersection point in the absorber
///- fVersor1 (TVector3) - direction of gamma ray comming form the source
///- fVersor2 (TVector3) - direction of the scattered gamma ray
///- fEnergy0 (Double_t) - initial gamma energy (set to be 4.44 MeV)
///- fEnergy1 (Double_t) - energy loss due to Compton scattering
///- fEnergy2 (Double_t) - energy of scattered gamma after scattering.
///
/// Additionally 4 histograms are filled:
///- hSource (TH2D*) - distribution of the gamma source
///- hScat (TH2D*) - distribution of events on scatterer
///- hAbs (TH2D*) - distribution of events on absorber
///- hEnergy (TH1D*) - energy of gammas after scaterring.
///
/// Details of the simulated setup are saved in text file and as TGeometry in
/// ROOT file.

class CCSimulation
{

public:
    CCSimulation() = delete;
    CCSimulation(const TString& name, const TString& outputPath, Int_t genVer, Double_t negAngle, Double_t posAngle, Bool_t verbose);
    ~CCSimulation();

    void BuildSetup(Double_t scatDist, Double_t scatZ, Double_t scatY, Double_t absDist,
                    Double_t absZ, Double_t absY);
    Bool_t GenerateRay();

    Bool_t ProcessEvent();
    void Loop(Int_t nev);
    void Clear();
    void SaveGeometryTxt();
    void BuildTGeometry();
    void SetCoordinate(Double_t x, Double_t y, Double_t z);
    /// Sets radius of the gamma source if generator #5 is used.
    ///\param radius (Double_t) - source radius.
    void SetRadius(Double_t radius) { fRadius = radius; };
    /// Sets distance between two point-like sources along Y axis
    /// if generator #4 is used.
    ///\param gapY (Double_t) - gap between gamma sources (in mm).
    void SetGapY(Double_t gapY) { fYgap = gapY; };
    /// Sets distance between two point-like sources along Z axis
    /// if generator #3 is used.
    ///\param gapZ (Double_t) - gap between gamma sources (in mm).
    void SetGapZ(Double_t gapZ) { fZgap = gapZ; };
    /// Sets version of the generator
    ///\param gen (Int_t) - generator version (numbers from 1 to 5).
    void SetGenVersion(Int_t gen) { fGenVersion = gen; };

private:
    // TString fOutputName;
    Bool_t fVerbose;       ///< Verbose level
    TString fOutputPath;   ///< Output path
    Int_t fGenVersion;     ///< Version of generator (numbers from 1 to 5)
    Int_t fNev;            ///< Counter of events currently in the acceptance of the absorber
    TTree* fTree{nullptr}; ///< Tree containing results of the simulations
    TFile* fFile{nullptr}; ///< Results ROOT file
    DetPlane fScatterer;   ///< Scatterer plane
    DetPlane fAbsorber;    ///< Absorber plane
    TVector3 fPoint0;      ///< Source of emitted gamma ray
    TVector3 fPoint1;      ///< Interaction point in the scatterer
    TVector3 fPoint2;      ///< Interaction point in the absorber
    TVector3 fVersor1;     ///< Direction of gamma ray comming from the source
    TVector3 fVersor2;     ///< Direction of scattered gamma ray
    Double_t fEnergy0;     ///< Initial gamma energy
    Double_t fEnergy1;     ///< Energy loss due to Compton scattering
    Double_t fEnergy2;     ///< Energy after Compton scattering
    Double_t fXofSource;   ///< x-component of source coordinate
    Double_t fYofSource;   ///< y-component of source coordinate
    Double_t fZofSource;   ///< z-component of source coordinate
    Double_t fYgap;        ///< For generator #4 - distance between two point-like sources
                           /// along Y axis
    Double_t fZgap;        ///< For generator #3 - distance between two point-like sources
                           /// along Z axis
    Double_t fRadius;      ///< For generator #5 - radius of the gamma source
    Double_t fNegAngle;    ///< Negative constraint on generated photons
    Double_t fPosAngle;    ///< Positive constraint on generated photons

    TH2F* hSource{nullptr};     ///< 2D histogram of distribution of the gamma source
    TH2F* hScat{nullptr};       ///< 2D histogram of distribution of evens on scatterer plane
    TH2F* hAbs{nullptr};        ///< 2D histogram of distribution of events on absorber plane
    TH1F* hEnergyAbs{nullptr};  ///< Histogram of scattered gammas energy
    TH1F* hEnergyLoss{nullptr}; ///< Histogram of energy loss
};

#endif
