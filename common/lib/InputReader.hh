#ifndef __InputReader_H_
#define __InputReader_H_ 1
#include "TFile.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"
#include <iostream>

using namespace std;

/// Base class for InputReaderSimple and InputReaderGeant. It allows
/// for reading simulation data saved in different forms and passing it
/// to reconstruction classes via set of getter functions.

class InputReader : public TObject
{

public:
    InputReader();
    InputReader(TString path);
    ~InputReader();

    void Print(void);
    //------------------------------------------------------------------
    /// Sets default values of the protected class members.
    void virtual Clear(void);
    //------------------------------------------------------------------
    /// Loads requested event from the opened tree with simulations results.
    ///\param i (int) - number of the requested event.
    bool virtual LoadEvent(int i);
    int virtual GetNumberOfEventsInFile(void);
    TVector3 virtual* GetPositionPrimary(void);
    TVector3 virtual* GetPositionScattering(void);
    TVector3 virtual* GetPositionAbsorption(void);
    TVector3 virtual* GetGammaDirPrimary(void);
    TVector3 virtual* GetGammaDirScattered(void);

    int virtual GetMultiplicityNum(void);
    int virtual GetClassID(void);
    // REMOVE/CHANGE IN EI
    double virtual GetEP(void);
    double virtual GetReES(void);
    double virtual GetES(void);
    //  double virtual GetReEP(void);
    double virtual GetEnergyPrimary(void);
    double virtual GetEnergyLoss(void);
    double virtual GetEnergyScattered(void);

    // SIMPLE SPECIFIC

    void virtual SetSmearing(bool smear, Double_t posX, Double_t posY, Double_t posZ);
    Double_t virtual SmearGaus(double val, double sigma);
    Double_t virtual SmearBox(double x, double resolution);
    Double_t virtual GetSigmaE(double energy);
    // EI SPECIFIC

    // GEANT4 SPECIFIC
    void virtual SetLoadMonteCarlo(void);
    void virtual SetLoadOnlyCorrect(void);

    // NN SPECIFIC

protected:
    TFile* fFile; ///< Input file
    TTree* fTree; ///< Tree containing simulation results
    /// Opens tree containing simulations results. Sets branches addresses
    ///\param name (TString) - name of the tree
    bool virtual AccessTree(TString name);
    bool SetInputFile(TString path);
    // vector<PhysicVec*>* GetRecoClusterPosSize(void);

    ClassDef(InputReader, 0)
};
inline int InputReader::GetNumberOfEventsInFile(void) { return 0; }
/// GEANT4 SPECIFIC
inline void InputReader::SetLoadMonteCarlo(void) {}
inline void InputReader::SetLoadOnlyCorrect(void) {}
/// SIMPLE SPECIFIC
inline void InputReader::SetSmearing(bool smear, Double_t posX, Double_t posY, Double_t posZ) {}
inline Double_t InputReader::SmearGaus(double val, double sigma) { return 0; }
inline Double_t InputReader::SmearBox(double x, double resolution) { return 0; }
inline Double_t InputReader::GetSigmaE(double energy) { return 0; }
#endif
