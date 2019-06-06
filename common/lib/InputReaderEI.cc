#include "InputReaderEI.hh"
#include <vector>
//ClassImp(InputReaderEI);

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

  if (!AccessTree("Cluster21", "Cluster22")) {
    throw "##### Exception in InputReaderEI constructor!";
  }
  fPositionScatReco = new TVector3();
  fPositionAbsReco = new TVector3();
}
//------------------------------------------------------------------
/// Default destructor.
InputReaderEI::~InputReaderEI() {
  if (fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
bool InputReaderEI::AccessTree(TString name, TString name1) {

  fTree = (TTree*)fFile->Get(name);
  fTree1 = (TTree*)fFile->Get(name1);
  
  if (fTree == NULL) {
    cout << "##### Error in InputReaderEI::AccessTree()!" << endl;
    cout << "Could not access the tree!" << endl;
    return false;
  }
  if (fTree1 == NULL) {
    cout << "##### Error in InputReaderEI::AccessTree()!" << endl;
    cout << "Could not access the tree!" << endl;
    return false;
  }
  
  fEnergyRecoScat2 = new PhysicVar();
  fEnergyRecoAbs2 = new PhysicVar();
  fPointRecoScat2 = new PhysicVec();
  fPointRecoAbs2 = new PhysicVec();


  fTree->SetBranchAddress("PointRecoScat2", &fPointRecoScat2);
  fTree->SetBranchAddress("EnergyRecoScat2", &fEnergyRecoScat2);
  fTree->SetBranchAddress("SizeScat2", &fSizeScat2);
  
  
  fTree1->SetBranchAddress("PointRecoAbs2", &fPointRecoAbs2);
  fTree1->SetBranchAddress("EnergyRecoAbs2", &fEnergyRecoAbs2);
  fTree1->SetBranchAddress("SizeAbs2", &fSizeAbs2);
  fTree1->SetBranchAddress("TotalEnergy'2'", &fTotalEnergy2);
  
  cout << "\n\nIn InputReaderEI::AccessTree()." << endl;
  cout << fTree->GetName() << " tree accessed.\n" << endl;

  cout << "\n\nIn InputReaderEI::AccessTree()." << endl;
  cout << fTree1->GetName() << " tree accessed.\n" << endl;
  
 
  
  return true;
}
//------------------------------------------------------------------
bool InputReaderEI::LoadEvent(int i) {

  int imax = fTree->GetEntries();
  fTree->GetEntry(i);
  
  int imax1 = fTree1->GetEntries();
  if (i > imax1) {
    cout << "##### Error in InputReaderEI::LoadEvent() in reconstruction tree!" << endl;
    cout << "Requested event number larger than number of events in the tree!"
         << endl;
    return false;
  }
  fTree1->GetEntry(i);
  

  return true;
}
//------------------------------------------------------------------
TVector3* InputReaderEI::GetPositionScattering(void)
{
  fPositionScatReco->SetX(fPointRecoScat2->position.X());
  fPositionScatReco->SetY(fPointRecoScat2->position.Y());
  fPositionScatReco->SetZ(fPointRecoScat2->position.Z());
  return fPositionScatReco;
}
//------------------------------------------------------------------
TVector3* InputReaderEI::GetPositionAbsorption(void) {
  fPositionAbsReco->SetX(fPointRecoAbs2->position.X());
  fPositionAbsReco->SetY(fPointRecoAbs2->position.Y());
  fPositionAbsReco->SetZ(fPointRecoAbs2->position.Z());
  return fPositionAbsReco;
}
//------------------------------------------------------------------
double InputReaderEI::GetEnergyLoss(void) { return fEnergyRecoScat2->value; }
//------------------------------------------------------------------
double InputReaderEI::GetEnergyScattered(void) { return fEnergyRecoAbs2->value; }
//------------------------------------------------------------------
int InputReaderEI::GetScatSize(void) { return fSizeScat2; }
//------------------------------------------------------------------
int InputReaderEI::GetAbsSize(void) { return fSizeAbs2; }
//------------------------------------------------------------------
double InputReaderEI::GetTotalEnergy(void) { return fTotalEnergy2; }
//------------------------------------------------------------------
void InputReaderEI::Clear(void) {
  
  fTotalEnergy2 = -1000;
  fSizeScat2 = -1;
  fSizeAbs2 = -1;
  fPositionScatReco = NULL;
  fPositionAbsReco = NULL;
  fEnergyRecoScat2 = NULL;
  fEnergyRecoAbs2 = NULL;
  fPointRecoScat2 = NULL;
  fPointRecoAbs2 = NULL;
  fTree = NULL;
  fTree1 = NULL;
  
  fFile = NULL;

  return;
}
//------------------------------------------------------------------
