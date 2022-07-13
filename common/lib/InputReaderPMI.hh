#ifndef __InputReaderPMI_H_
#define __InputReaderPMI_H_ 1
#include "InputReader.hh"
#include "TString.h"
#include "TVector3.h"
#include <iostream>
using namespace std;

/// This is class derived from InputReader class. It opens requested
/// ROOT file containing tree with reconstruction results by for example, Geant4
/// and via set of getter function passes information to reconstruction 
/// classes, i.e. CCMLEM.

class InputReaderPMI : public InputReader {

public:
    InputReaderPMI();
    InputReaderPMI(TString path);
    ~InputReaderPMI();

    bool LoadEvent(int i);
    void Clear(void);
    int GetNumberOfEventsInFile(void);
    void SelectEvents(void);
    list<int> GetSelectedEvents(void);

    TVector3* GetPositionScattering(void);
    TVector3* GetPositionAbsorption(void);
    TVector3* GetGammaDirScattered(void);

    double GetEnergyLoss(void);
    double GetEnergyScattered(void);

private:

    //LOADING FROM THE TREE
    TVector3* fPoint1; ///< Coordinates of the interaction in the scatterer
    TVector3* fPoint2; ///< Coordinates of the interaction in the absorber (absorption)
    double fEnergy1;    ///< Energy deposited in the scatterer [keV]
    double fEnergy2;    ///< Energy of the scattered gamma [keV]
    //PARAMETERS PASSED ON
    TVector3* fPositionScat;  ///< Position of interaction in scatterer
    TVector3* fPositionAbs;   ///< Position of interaction in absorber
    TVector3* fDirectionScat; ///< Direction of the scattered gamma
    double fEnergyLoss;
    double fEnergyScattered;
    int fabsclustersize;

    bool AccessTree(TString name);
    //TTree* fTree;
    list<int> fSelectedEvents;
    bool SelectSingleEvent(void);

    ClassDef(InputReaderPMI, 0)
};

inline int InputReaderPMI::GetNumberOfEventsInFile(void){return fTree->GetEntries();}
inline list<int> InputReaderPMI::GetSelectedEvents(void){return fSelectedEvents;}

#endif