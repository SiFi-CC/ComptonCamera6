#include "InputReader.hh"
//#include <vector>
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
///\brief Get number of entries in tree
///\details InputReader::AccessTree has to be called before.
///\return number of events in fTree
Int_t InputReader::GetEntries(void){
  if (fTree == NULL) {
    cout << "##### Error in InputReader::GetEntries()" << endl;
    cout << "No tree loaded to read entries from." << endl;
    return 0;
  }
  return fTree->GetEntries();
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
/// Returns pointer to the vector representing place of interaction
/// in the absorber (absorption).
TVector3* InputReader::GetPositionAbsorption(void) { return NULL; }
//------------------------------------------------------------------
/// Returns pointer to the vector representing direction of the primary
/// gamma
TVector3* InputReader::GetGammaDirPrimary(void) { return NULL; }
//------------------------------------------------------------------
/// Returns pointer to the vector representing direction of the
/// scattered gamma.
TVector3* InputReader::GetGammaDirScattered(void) { return NULL; }
//------------------------------------------------------------------
/// Returns energy of the gamma emitted from the source [MeV].
double InputReader::GetEnergyPrimary(void) { return -100; }
//------------------------------------------------------------------
/// Returns energy deposited in the scatterer [MeV].
double InputReader::GetEnergyLoss(void) { return -100; }
//------------------------------------------------------------------
/// Returns energy of the scattered gamma [MeV].
double InputReader::GetEnergyScattered(void) { return -100; }
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
///\brief Check if a single event fulfils a given selection cut.
///\details The cut is passes as a ROOT TCut element using the variable
/// names from the tree. 
/// For the Geant4 tree, a valid cut would for example be\n
/// RecoEnergy_e.value+RecoEnergy_p.value>4.3\n\n
/// Several cuts can be added as cut1+cut2.
///\param i event number
///\param cut selection cut
///\return true if event i passes selection cut, false if it doesn't pass or in case of error
bool InputReader::ApplySelectionCut(int i, TCut cut){
  if (fTree == NULL) {
    cout << "##### Error in InputReader::ApplySelectionCut()" << endl;
    cout << "No tree loaded." << endl;
    return false;
  }
  // no = 
  // 1: event passed cut
  // 0: event did not pass cut
  //-1: error
  int no = fTree->Draw("",cut,"goff",1,i);
  if (1==no) return true;
  else if (-1==no){
    cout << "##### Error in InputReader::ApplySelectionCut" <<endl;
    cout << "Applying selection " << cut.GetTitle() << " to event " << i 
			<< " of tree "<< fTree->GetName() << " failed." <<endl;
  }
  return false;
}
//------------------------------------------------------------------
