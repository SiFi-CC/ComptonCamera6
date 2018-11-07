#include "InputReaderSimple.hh"

ClassImp(InputReaderSimple);

//------------------------------------------------------------------
/// Default constructor.
InputReaderSimple::InputReaderSimple() : InputReader() {
  Clear();
  cout << "##### Warning in InputReaderSimple constructor!" << endl;
  cout << "You are using default constructor." << endl;
}
//------------------------------------------------------------------
/// Standard constructor (recommended).
///\param path (TString) - path to the input file.
InputReaderSimple::InputReaderSimple(TString path) : InputReader(path) {

  if (!AccessTree("data")) {
    throw "##### Exception in InputReaderSimple constructor!";
  }
}
//------------------------------------------------------------------
/// Default destructor.
InputReaderSimple::~InputReaderSimple() {
  if (fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
bool InputReaderSimple::AccessTree(TString name) {

  fTree = (TTree*)fFile->Get(name);

  if (fTree == NULL) {
    cout << "##### Error in InputReaderSimple::AccessTree()!" << endl;
    cout << "Could not access the tree!" << endl;
    return false;
  }

  fPoint0 = new TVector3();
  fPoint1 = new TVector3();
  fPoint2 = new TVector3();
  fVersor1 = new TVector3();
  fVersor2 = new TVector3();

  fTree->SetBranchAddress("point0", &fPoint0);
  fTree->SetBranchAddress("point1", &fPoint1);
  fTree->SetBranchAddress("point2", &fPoint2);
  fTree->SetBranchAddress("versor1", &fVersor1);
  fTree->SetBranchAddress("versor2", &fVersor2);
  fTree->SetBranchAddress("energy0", &fEnergy0);
  fTree->SetBranchAddress("energy1", &fEnergy1);
  fTree->SetBranchAddress("energy2", &fEnergy2);

  cout << "\n\nIn InputReaderSimple::AccessTree()." << endl;
  cout << fTree->GetName() << " tree accessed.\n" << endl;

  return true;
}
//------------------------------------------------------------------
TVector3* InputReaderSimple::GetPositionPrimary(void) { return fPoint0; }
//------------------------------------------------------------------
TVector3* InputReaderSimple::GetPositionScattering(void) { return fPoint1; }
//------------------------------------------------------------------
TVector3* InputReaderSimple::GetPositionAbsorption(void) { return fPoint2; }
//------------------------------------------------------------------
TVector3* InputReaderSimple::GetGammaDirPrimary(void) { return fVersor1; }
//------------------------------------------------------------------
TVector3* InputReaderSimple::GetGammaDirScattered(void) { return fVersor2; }
//------------------------------------------------------------------
double InputReaderSimple::GetEnergyPrimary(void) { return fEnergy0; }
//------------------------------------------------------------------
double InputReaderSimple::GetEnergyLoss(void) { return fEnergy1; }
//------------------------------------------------------------------
double InputReaderSimple::GetEnergyScattered(void) { return fEnergy2; }
//------------------------------------------------------------------
void InputReaderSimple::Clear(void) {
  fPoint0 = NULL;
  fPoint1 = NULL;
  fPoint2 = NULL;
  fVersor1 = NULL;
  fVersor2 = NULL;
  fFile = NULL;
  fTree = NULL;
  fEnergy0 = -100;
  fEnergy1 = -100;
  fEnergy2 = -100;
  return;
}
//------------------------------------------------------------------