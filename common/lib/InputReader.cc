#include "InputReader.hh"

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
bool InputReader::AccessTree(TString name1) {

  fTree = (TTree*)fFile->Get(name1);
  if (fTree == NULL) {
    cout << "##### Error in InputReader::AccessTree()!" << endl;
    cout << "Could not access tree: " << name1 << endl;
    return false;
  }
  /*fTree1 = (TTree*)fFile->Get(name2);
  if (fTree1 == NULL) {
    cout << "##### Error in InputReader::AccessTree()!" << endl;
    cout << "Could not access tree: " << name2 << endl;
    return false;
  }
*/
  return true;
}
//------------------------------------------------------------------
/// Loads requested event from the opened tree with simulations results.
///\param i (int) - number of the requested event.
bool InputReader::LoadEvent(int i) {

  int imax = fTree->GetEntries();
  // int imax1 = fTree1->GetEntries();

  if (i > imax) {
    cout << "##### Error in InputReader::LoadEvent()!" << endl;
    cout << "Requested event number larger than number of events in the tree!"
         << endl;
    return false;
  }

  fTree->GetEntry(i);
  // fTree1->GetEntry(i);
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
  fTree1 = NULL;
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
