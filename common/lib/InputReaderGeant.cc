#include "InputReaderGeant.hh"

ClassImp(InputReaderGeant);

//------------------------------------------------------------------
/// Default constructor.
InputReaderGeant::InputReaderGeant() : InputReader() {
  Clear();
  cout << "##### Warning in InputReaderGeant constructor!" << endl;
  cout << "You are usinf default constructor." << endl;
}
//------------------------------------------------------------------
/// Standard constructor.
///\param path (TString) - path to the input file.
InputReaderGeant::InputReaderGeant(TString path) : InputReader(path) {

  if (!AccessTree("G4SimulationData_Reconstruction")) {
    throw "##### Exception in InputReaderGeant constructor!";
  }

  fPositionScat = new TVector3();
  fPositionAbs = new TVector3();
  fDirectionScat = new TVector3();
}
//------------------------------------------------------------------
/// Default destructor.
InputReaderGeant::~InputReaderGeant() {
  if (fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
bool InputReaderGeant::AccessTree(TString name) {

  fTree = (TTree*)fFile->Get(name);

  if (fTree == NULL) {
    cout << "##### Error in InputReaderGeant::AccessTree()!" << endl;
    cout << "Could not access the tree!" << endl;
    return false;
  }

  fRecoEnergy_e = new PhysicVar();
  fRecoEnergy_p = new PhysicVar();
  fRecoPosition_e = new PhysicVec();
  fRecoPosition_p = new PhysicVec();
  fRecoDirection_scatter = new PhysicVec();
  fRecoClusterPositions = 0;
  fRecoClusterEnergies = 0;

  fTree->SetBranchAddress("EventNumber", &fEventNumber);
  fTree->SetBranchAddress("Identified", &fIdentified);
  fTree->SetBranchAddress("RecoEnergy_e", &fRecoEnergy_e);
  fTree->SetBranchAddress("RecoEnergy_p", &fRecoEnergy_p);
  fTree->SetBranchAddress("RecoPosition_e", &fRecoPosition_e);
  fTree->SetBranchAddress("RecoPosition_p", &fRecoPosition_p);
  fTree->SetBranchAddress("RecoDirection_scatter", &fRecoDirection_scatter);
  fTree->SetBranchAddress("RecoClusterPositions", &fRecoClusterPositions);
  fTree->SetBranchAddress("RecoClusterEnergies", &fRecoClusterEnergies);

  cout << "\n\nIn InputReaderGeant::AccessTree()." << endl;
  cout << fTree->GetName() << " tree accessed.\n" << endl;

  return true;
}
//------------------------------------------------------------------
bool InputReaderGeant::LoadEvent(int i) {

  int imax = fTree->GetEntries();
  if (i > imax) {
    cout << "##### Error in InputReaderGeant::LoadEvent()!" << endl;
    cout << "Requested event number larger than number of events in the tree!"
         << endl;
    return false;
  }

  fTree->GetEntry(i);

  return fIdentified;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionPrimary(void) {
  cout << "##### Warning in InputReaderGeant::GetPositionPrimary()!" << endl;
  cout << "\t Position of gamma source is unknown!" << endl;
  return NULL;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionScattering(void) {
  fPositionScat->SetX(fRecoPosition_e->position.X());
  fPositionScat->SetY(fRecoPosition_e->position.Y());
  fPositionScat->SetZ(fRecoPosition_e->position.Z());
  return fPositionScat;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionAbsorption(void) {
  fPositionAbs->SetX(fRecoPosition_p->position.X());
  fPositionAbs->SetY(fRecoPosition_p->position.Y());
  fPositionAbs->SetZ(fRecoPosition_p->position.Z());
  return fPositionAbs;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetGammaDirPrimary(void) {
  cout << "##### Warning in InputReaderGeant::GetGammaDirPrimary()!" << endl;
  cout << "\t Direction of primary gamma is unknown!" << endl;
  return NULL;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetGammaDirScattered(void) {
  TVector3* temp = new TVector3();
  fDirectionScat->SetX(fRecoDirection_scatter->position.X());
  fDirectionScat->SetY(fRecoDirection_scatter->position.Y());
  fDirectionScat->SetZ(fRecoDirection_scatter->position.Z());
  return fDirectionScat;
}
//------------------------------------------------------------------
double InputReaderGeant::GetEnergyPrimary(void) {
  double sum = fRecoEnergy_e->value + fRecoEnergy_p->value;
  return sum;
}
//------------------------------------------------------------------
double InputReaderGeant::GetEnergyLoss(void) { return fRecoEnergy_e->value; }
//------------------------------------------------------------------
double InputReaderGeant::GetEnergyScattered(void) {
  return fRecoEnergy_p->value;
}
//------------------------------------------------------------------
void InputReaderGeant::Clear(void) {
  fEventNumber = -1;
  fIdentified = false;
  fRecoEnergy_e = NULL;
  fRecoEnergy_p = NULL;
  fRecoPosition_e = NULL;
  fRecoPosition_p = NULL;
  fPositionScat = NULL;
  fPositionAbs = NULL;
  fDirectionScat = NULL;
  fTree = NULL;
  fFile = NULL;
  fRecoDirection_scatter = NULL;
  if (!fRecoClusterPositions->empty()) fRecoClusterPositions->clear();
  if (!fRecoClusterEnergies->empty()) fRecoClusterEnergies->clear();
  return;
}
//------------------------------------------------------------------