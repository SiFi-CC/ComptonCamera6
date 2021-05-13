#include "InputReaderEI.hh"
#include <vector>
ClassImp(InputReaderEI);

//------------------------------------------------------------------
/// Default constructor.
InputReaderEI::InputReaderEI() : InputReader() {
  Clear();
  cout << "##### Warning in InputReaderEI constructor!" << endl;
  cout << "You are using default constructor." << endl;
} 
//------------------------------------------------------------------
/// Standard constructor.
///\param path (TString) - path to the input file.
InputReaderEI::InputReaderEI(TString path) : InputReader(path) {

  if (!AccessTree("TreeSB")) {
    throw "##### Exception in InputReaderEI constructor!";
  }
  fPositionElectron = new TVector3();
  fPositionPhoton = new TVector3();
}
//------------------------------------------------------------------
/// Default destructor.
InputReaderEI::~InputReaderEI() {
  if (fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
/// Accesses data of trees'branches in ROOT file.
///\param name (TString) - name of tree.
bool InputReaderEI::AccessTree(TString name/*, TString name1, TString name2*/) {

  fTree = (TTree*)fFile->Get(name);

  
  if (fTree == NULL) {
    cout << "##### Error in InputReaderEI::AccessTree()!" << endl;
    cout << "Could not access the tree!" << endl;
    return false;
  }
  
  fPos_Scat = new TVector3();
  fPos_Abs = new TVector3();

  fTree->SetBranchAddress("EventNumber", &fEventNumber);
  fTree->SetBranchAddress("Energy_Primary", &fEnergy_Primary);
  fTree->SetBranchAddress("ReEnergy_Sum", &fReEnS);
  fTree->SetBranchAddress("ReEnergy_Primary", &fReEnP);
//    fTree->SetBranchAddress("PosX_Scat", &fPosX_Scat);
//    fTree->SetBranchAddress("PosY_Scat", &fPosY_Scat);
//    fTree->SetBranchAddress("PosZ_Scat", &fPosZ_Scat);
  fTree->SetBranchAddress("Pos_Scat", &fPos_Scat);
  fTree->SetBranchAddress("Energy_Scat", &fEnergy_Scat);
//    fTree->SetBranchAddress("PosX_Abs", &fPosX_Abs);
//    fTree->SetBranchAddress("PosY_Abs", &fPosY_Abs);
//    fTree->SetBranchAddress("PosZ_Abs", &fPosZ_Abs);
  fTree->SetBranchAddress("Pos_Abs", &fPos_Abs);
  fTree->SetBranchAddress("Energy_Abs", &fEnergy_Abs);
  fTree->SetBranchAddress("Energy_Sum", &fEnergyS);
  //fTree->SetBranchAddress("Multiplicity", &fS);
  fTree->SetBranchAddress("classID", &fClassID);
  

  
  cout << "\n\nIn InputReaderEI::AccessTree()." << endl;
  cout << fTree->GetName() << " tree accessed.\n" << endl;


 
  
  return true;
}
//------------------------------------------------------------------
/// loads events from trees to analyze them in CCMLEM class.
///\param i (int) - number of events
bool InputReaderEI::LoadEvent(int i) {

  int imax = fTree->GetEntries();
  if (i > imax) {
    cout << "##### Error in InputReaderEI::LoadEvent() in cluster tree!" << endl;
    cout << "Requested event number larger than number of events in the tree!"
         << endl;
    return false;
  }
  fTree->GetEntry(i);
  

  return true;
}
//------------------------------------------------------------------
TVector3* InputReaderEI::GetPositionScattering(void)
{
  fPositionElectron->SetX(fPos_Scat->X());
  fPositionElectron->SetY(fPos_Scat->Y());
  fPositionElectron->SetZ(fPos_Scat->Z());
  return fPositionElectron;
}
//------------------------------------------------------------------
TVector3* InputReaderEI::GetPositionAbsorption(void) {
  fPositionPhoton->SetX(fPos_Abs->X());
  fPositionPhoton->SetY(fPos_Abs->Y());
  fPositionPhoton->SetZ(fPos_Abs->Z());
  return fPositionPhoton;
}
//------------------------------------------------------------------
double InputReaderEI::GetEP(void) { return fEnergy_Primary; }
//------------------------------------------------------------------
double InputReaderEI::GetReES(void) { return fReEnS; }
//-----------------------------------------------------------------
double InputReaderEI::GetES(void) { return fEnergyS; }
//------------------------------------------------------------------
double InputReaderEI::GetReEP(void) { return fReEnP; }
//------------------------------------------------------------------
double InputReaderEI::GetEnergyLoss(void) { return fEnergy_Scat; }
//------------------------------------------------------------------
double InputReaderEI::GetEnergyScattered(void) { return fEnergy_Abs; }
//------------------------------------------------------------------
int InputReaderEI::GetMultiplicityNum(void) { return fS; }
//------------------------------------------------------------------
int InputReaderEI::GetClassID(void) { return fClassID; }
//------------------------------------------------------------------
void InputReaderEI::Clear(void) {
  
 
  fEventNumber = -1;
  
  fS = -1;
  fClassID = -1;
  fPos_Scat = NULL;
  fPos_Abs = NULL;
  fEnergy_Primary = -1000;
  fReEnS = -1000;
  fReEnP = -1000;
  fEnergyS = -1000;
  fEnergy_Scat = -1000;
  fEnergy_Abs = -1000;
  
  fPositionElectron = NULL;
  fPositionPhoton = NULL;

  fTree = NULL;

  
  fFile = NULL;

  return;
}
//------------------------------------------------------------------
