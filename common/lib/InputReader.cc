#include "InputReader.hh"
#include <vector>
ClassImp(InputReader);

//------------------------------------------------------------------
/// Deafault constructor.
InputReader::InputReader() {
  cout << "##### Warning in InputReader constructor!" << endl;
  cout << "You are using default constructor. Set the input file!" << endl;
  Clear();
}
//------------------------------------------------------------------
/// Deafault destructor.
InputReader::~InputReader() {
  if (fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
/// Standard constructor (recommended).
///\param path (TString) - path to the input file.
InputReader::InputReader(TString path) {
  Clear();
  if (!SetInputFile(path)) {
    throw "##### Exception in InputReader constructor!";
  }
}
//------------------------------------------------------------------
/// Opens input file contatinig tree with simulations results.
///\param path (TString) - path to the input file.
bool InputReader::SetInputFile(TString path) {

  fFile = new TFile(path, "READ");
  if (!fFile->IsOpen()) {
    cout << "##### Error in InputReader::SetInputFile()!" << endl;
    cout << "Could not open the file!" << endl;
    return false;
  }

  return true;
}
//------------------------------------------------------------------
/// Opens tree containing simulations results. Sets branches addresses
///\param name (TString) - name of the tree
bool InputReader::AccessTree(TString name) {

  fTree = (TTree*)fFile->Get(name);
  if (fTree == NULL) {
    cout << "##### Error in InputReader::AccessTree()!" << endl;
    cout << "Could not access tree: " << name << endl;
    return false;
  }
  
  return true;
}
//------------------------------------------------------------------
/// Loads requested event from the opened tree with simulations results.
///\param i (int) - number of the requested event.
bool InputReader::LoadEvent(int i) {

  int imax = fTree->GetEntries();

  if (i > imax) {
    cout << "##### Error in InputReader::LoadEvent()!" << endl;
    cout << "Requested event number larger than number of events in the tree!"
         << endl;
    return false;
  }

  fTree->GetEntry(i);
  return true;
}
//------------------------------------------------------------------
/// Returns pointer to the vector representing souurce of the gamma.
TVector3* InputReader::GetPositionPrimary(void) { return NULL; }
//------------------------------------------------------------------
/// Returns pointer to the vector representing place of interaction
/// in the scatterer (Compton scattering).
TVector3* InputReader::GetPositionScattering(void) { return NULL; }
//------------------------------------------------------------------
vector<TVector3>* InputReader::GetElectronPosition(void) { return NULL; }
//------------------------------------------------------------------
int InputReader::GetRealPosESize(void) { return -100; }    
//------------------------------------------------------------------

vector<int>* InputReader::GetRealInteractionE(void) { 
    
    return NULL; }
//------------------------------------------------------------------
TVector3* InputReader::GetPositionScatteringReco(void) { return NULL; }
//------------------------------------------------------------------
/// Returns pointer to the vector representing place of interaction
/// in the absorber (absorption).
TVector3* InputReader::GetPositionAbsorption(void) { return NULL; }
//------------------------------------------------------------------
vector<TVector3>* InputReader::GetPhotonPosition(void) { return NULL; }
//------------------------------------------------------------------
int InputReader::GetRealPosPSize(void) { return -100; }      
//------------------------------------------------------------------
vector<int>* InputReader::GetRealInteractionP(void) { 
    
    return NULL; }
//------------------------------------------------------------------
TVector3* InputReader::GetPositionAbsorptionReco(void) { return NULL; }
//------------------------------------------------------------------
/// Returns pointer to the vector representing direction of the primary
/// gamma
TVector3* InputReader::GetGammaDirPrimary(void) { return NULL; }
//------------------------------------------------------------------
/// Returns pointer to the vector representing direction of the
/// scattered gamma.
TVector3* InputReader::GetGammaDirScattered(void) { return NULL; }
//------------------------------------------------------------------
TVector3* InputReader::GetGammaDirScatteredReco(void) { return NULL; }
//------------------------------------------------------------------
int InputReader::GetRecoClusterPosSize(void) { return -100; }
//------------------------------------------------------------------
int InputReader::GetIdentified(void) {return -1000;}

//------------------------------------------------------------------
int InputReader::GetMultiplicityNum(void) {return -1000;}
//------------------------------------------------------------------
int InputReader::GetClassID(void) {return -1000;}
//------------------------------------------------------------------
double InputReader::GetEP(void) { return -100; }
//------------------------------------------------------------------
double InputReader::GetReES(void) { return -100; }
//------------------------------------------------------------------
double InputReader::GetES(void) { return -100; }
//------------------------------------------------------------------
double InputReader::GetReEP(void) { return -100; }
//------------------------------------------------------------------
/// Returns energy of the gamma emitted from the source [MeV].
double InputReader::GetEnergyPrimary(void) { return -100; }
//------------------------------------------------------------------
double InputReader::GetEnergyPrimaryReco(void) { return -100; }
//------------------------------------------------------------------
/// Returns energy deposited in the scatterer [MeV].
double InputReader::GetEnergyLoss(void) { return -100; }
//------------------------------------------------------------------
double InputReader::GetEnergyLossReco(void) { return -100; }
//------------------------------------------------------------------
/// Returns energy of the scattered gamma [MeV].
double InputReader::GetEnergyScattered(void) { return -100; }
//------------------------------------------------------------------
double InputReader::GetEnergyScatteredReco(void) { return -100; }
//------------------------------------------------------------------
TVector3* InputReader::GetScattererPosition(void) { return NULL; }
//------------------------------------------------------------------
TVector3* InputReader::GetAbsorberPosition(void) { return NULL; }
//------------------------------------------------------------------
double InputReader::GetScatThickx(void) { return -100; }
//------------------------------------------------------------------
double InputReader::GetScatThicky(void) { return -100; }
//------------------------------------------------------------------
double InputReader::GetScatThickz(void) { return -100; }
//------------------------------------------------------------------
double InputReader::GetAbsThickx(void) { return -100; }
//------------------------------------------------------------------
double InputReader::GetAbsThicky(void) { return -100; }
//------------------------------------------------------------------
double InputReader::GetAbsThickz(void) { return -100; }
//------------------------------------------------------------------
/// Sets default values of the protected class members.
void InputReader::Clear(void) {
  fFile = NULL;
  fTree = NULL;
  
}
//------------------------------------------------------------------
/// Prints details of the InputReader class object.
void InputReader::Print(void) {
  cout << "\n-------------------------------------------------------" << endl;
  cout << "This is Print() for InputReader class object" << endl;
  if (fFile != NULL) {
    cout << "Opened file: " << fFile->GetName() << endl;
  } else {
    cout << "It's empty!" << endl;
  }
  cout << "-------------------------------------------------------\n" << endl;
}
//------------------------------------------------------------------
