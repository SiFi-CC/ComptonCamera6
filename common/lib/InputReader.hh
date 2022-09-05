#ifndef __InputReader_H_
#define __InputReader_H_ 1
#include "TFile.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"
#include <iostream>
#include <list>

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
  /// Returns the number of events in the file.
  int virtual GetNumberOfEventsInFile(void);
  /// Applies the cuts to all events in the files and determine the ones passing the cuts. The event numbers are stored in a std::list<int>.
  void virtual SelectEvents(void); 
  /// Applies the cut to a single event in the file. Used in SelectEvents 
  bool virtual SelectSingleEvent(void);
  /// Returns the list of event numbers that fullfilled the cuts. These were determined in SelectEvents 
  list<int> virtual GetSelectedEvents(void);
  TVector3 virtual* GetPositionPrimary(void);
  TVector3 virtual* GetPositionScattering(void);
  TVector3 virtual* GetPositionAbsorption(void);
  TVector3 virtual* GetGammaDirPrimary(void);
  TVector3 virtual* GetGammaDirScattered(void);
  
  int virtual GetMultiplicityNum(void);
  int virtual GetClassID(void);
//EI Input Reader is at the moment not normalized to the other ones since it was a real mess
  double virtual GetEP(void);
  double virtual GetReES(void);
  double virtual GetES(void);
//  double virtual GetReEP(void);
  double virtual GetEnergyPrimary(void);
  double virtual GetEnergyLoss(void);
  double virtual GetEnergyScattered(void);

    // SIMPLE SPECIFIC

    void virtual SetSmearing(bool smear, Double_t posScatX, Double_t posAbsX, Double_t posScatY, Double_t posAbsY, Double_t posScatZ, Double_t posAbsZ);
    Double_t virtual SmearGaus(double val, double sigma);
    Double_t virtual SmearBox(double x, double resolution);
    Double_t virtual GetSigmaEScat(double energy);
    Double_t virtual GetSigmaEAbs(double energy);
    void virtual ReadGeometry(void) { }; ///< Reads scatterer's and absorber's width and height
    void virtual PrintEvent(int i) { };
    // EI SPECIFIC

//GEANT4 SPECIFIC
  void virtual SetLoadMonteCarlo(void);
//Geant4/NN SPECIFIC
  void virtual SetLoadOnlyCorrect(int value);

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

inline void InputReader::SelectEvents(void){}
inline bool InputReader::SelectSingleEvent(void){return false;}
inline int InputReader::GetNumberOfEventsInFile(void){return 0;}
inline list<int> InputReader::GetSelectedEvents(void){return std::list<int>();}
/// GEANT4 SPECIFIC
inline void InputReader::SetLoadMonteCarlo(void){}
inline void InputReader::SetLoadOnlyCorrect(int value){}
/// SIMPLE SPECIFIC
inline void InputReader::SetSmearing(bool smear, Double_t posScatX, Double_t posAbsX, Double_t posScatY, Double_t posAbsY, Double_t posScatZ, Double_t posAbsZ) {}
inline Double_t InputReader::SmearGaus(double val, double sigma) { return 0; }
inline Double_t InputReader::SmearBox(double x, double resolution) { return 0; }
inline Double_t InputReader::GetSigmaEScat(double energy) { return 0; }
inline Double_t InputReader::GetSigmaEAbs(double energy) { return 0; }
#endif
