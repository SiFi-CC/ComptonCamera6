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

  if (!AccessTree("G4SimulationData_Real", "G4SimulationData_Setup")) {
    throw "##### Exception in InputReaderGeant constructor!";
  }

  fPositionScat = new TVector3();
  fPositionAbs = new TVector3();
  fDirectionScat = new TVector3();
  fPositionSource = new TVector3();
  fDirectionSource = new TVector3();

  fScatPlanePos = new TVector3();
  fAbsPlanePos = new TVector3();
}
//------------------------------------------------------------------
/// Default destructor.
InputReaderGeant::~InputReaderGeant() {
  if (fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
bool InputReaderGeant::AccessTree(TString name, TString name1) {

  fTree = (TTree*)fFile->Get(name);
  fTree1 = (TTree*)fFile->Get(name1);
  if (fTree == NULL) {
    cout << "##### Error in InputReaderGeant::AccessTree()!" << endl;
    cout << "Could not access the tree!" << endl;
    return false;
  }
  if (fTree1 == NULL) {
    cout << "##### Error in InputReaderGeant::AccessTree()!" << endl;
    cout << "Could not access the tree!" << endl;
    return false;
  }
  // fEnergy_Primary = new PhysicVar();
  // fRealEnergy_e = new PhysicVar();
  // fRealEnergy_p = new PhysicVar();
  fRealPosition_source = new TVector3();
  fRealDirection_source = new TVector3();
  fRealPosition_e = new TVector3();
  fRealPosition_p = new TVector3();
  fRealDirection_scatter = new TVector3();
  // fRecoClusterPositions = 0;
  // fRecoClusterEnergies = 0;
  fScattererPosition = new TVector3();
  fAbsorberPosition = new TVector3();

  fTree->SetBranchAddress("EventNumber", &fEventNumber);
  fTree->SetBranchAddress("Energy_Primary", &fEnergy_Primary);
  fTree->SetBranchAddress("RealEnergy_e", &fRealEnergy_e);
  fTree->SetBranchAddress("RealEnergy_p", &fRealEnergy_p);
  fTree->SetBranchAddress("RealPosition_source", &fRealPosition_source);
  fTree->SetBranchAddress("RealDirection_source", &fRealDirection_source);
  fTree->SetBranchAddress("RealPosition_e", &fRealPosition_e);
  fTree->SetBranchAddress("RealPosition_p", &fRealPosition_p);
  fTree->SetBranchAddress("RealDirection_scatter", &fRealDirection_scatter);
  // fTree->SetBranchAddress("RecoClusterPositions", &fRecoClusterPositions);
  // fTree->SetBranchAddress("RecoClusterEnergies", &fRecoClusterEnergies);
  fTree1->SetBranchAddress("ScattererThickness_x", &fScattererThickness_x);
  fTree1->SetBranchAddress("ScattererThickness_y", &fScattererThickness_y);
  fTree1->SetBranchAddress("ScattererThickness_z", &fScattererThickness_z);
  fTree1->SetBranchAddress("AbsorberThickness_x", &fAbsorberThickness_x);
  fTree1->SetBranchAddress("AbsorberThickness_y", &fAbsorberThickness_y);
  fTree1->SetBranchAddress("AbsorberThickness_z", &fAbsorberThickness_z);
  fTree1->SetBranchAddress("ScattererPosition", &fScattererPosition);
  fTree1->SetBranchAddress("AbsorberPosition", &fAbsorberPosition);
  fTree1->SetBranchAddress("NumberOfSimulatedEvents",
                           &fNumberOfSimulatedEvents);

  cout << "\n\nIn InputReaderGeant::AccessTree()." << endl;
  cout << fTree->GetName() << " tree accessed.\n" << endl;

  cout << "\n\nIn InputReaderGeant::AccessTree()." << endl;
  cout << fTree1->GetName() << " tree accessed.\n" << endl;
  return true;
}
//------------------------------------------------------------------
bool InputReaderGeant::LoadEvent(int i) {

  int imax1 = fTree1->GetEntries();
  fTree1->GetEntry(i);
  int imax = fTree->GetEntries();
  if (i > imax) {
    cout << "##### Error in InputReaderGeant::LoadEvent()!" << endl;
    cout << "Requested event number larger than number of events in the tree!"
         << endl;
    return false;
  }

  fTree->GetEntry(i);

  if (0 == fRealPosition_e->X() || 0 == fRealPosition_p->X()) return false;

  return true;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionPrimary(void) {
  // cout << "##### Warning in InputReaderGeant::GetPositionPrimary()!" << endl;
  // cout << "\t Position of gamma source is unknown!" << endl;
  fPositionSource->SetX(fRealPosition_source->X());
  fPositionSource->SetY(fRealPosition_source->Y());
  fPositionSource->SetZ(fRealPosition_source->Z());
  return fPositionSource;
  // return NULL;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionScattering(void) {
  fPositionScat->SetX(fRealPosition_e->X());
  fPositionScat->SetY(fRealPosition_e->Y());
  fPositionScat->SetZ(fRealPosition_e->Z());
  return fPositionScat;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionAbsorption(void) {
  fPositionAbs->SetX(fRealPosition_p->X());
  fPositionAbs->SetY(fRealPosition_p->Y());
  fPositionAbs->SetZ(fRealPosition_p->Z());
  return fPositionAbs;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetGammaDirPrimary(void) {
  // cout << "##### Warning in InputReaderGeant::GetGammaDirPrimary()!" << endl;
  // cout << "\t Direction of primary gamma is unknown!" << endl;
  fDirectionSource->SetX(fRealDirection_source->X());
  fDirectionSource->SetY(fRealDirection_source->Y());
  fDirectionSource->SetZ(fRealDirection_source->Z());
  return fDirectionSource;
  // return NULL;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetGammaDirScattered(void) {
  // TVector3* temp = new TVector3();
  fDirectionScat->SetX(fRealDirection_scatter->X());
  fDirectionScat->SetY(fRealDirection_scatter->Y());
  fDirectionScat->SetZ(fRealDirection_scatter->Z());
  return fDirectionScat;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetScattererPosition(void) {
  fScatPlanePos->SetX(fScattererPosition->X());
  fScatPlanePos->SetY(fScattererPosition->Y());
  fScatPlanePos->SetZ(fScattererPosition->Z());
  return fScatPlanePos;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetAbsorberPosition(void) {
  fAbsPlanePos->SetX(fAbsorberPosition->X());
  fAbsPlanePos->SetY(fAbsorberPosition->Y());
  fAbsPlanePos->SetZ(fAbsorberPosition->Z());
  return fAbsPlanePos;
}
//------------------------------------------------------------------
double InputReaderGeant::GetEnergyPrimary(void) {
  double sum = fRealEnergy_e + fRealEnergy_p;
  return sum;
}
//------------------------------------------------------------------
double InputReaderGeant::GetEnergyLoss(void) { return fRealEnergy_e; }
//------------------------------------------------------------------
double InputReaderGeant::GetEnergyScattered(void) { return fRealEnergy_p; }
//------------------------------------------------------------------
double InputReaderGeant::GetScatThickx(void) { return fScattererThickness_x; }
//------------------------------------------------------------------
double InputReaderGeant::GetScatThicky(void) { return fScattererThickness_y; }
//------------------------------------------------------------------
double InputReaderGeant::GetScatThickz(void) { return fScattererThickness_z; }
//------------------------------------------------------------------
double InputReaderGeant::GetAbsThickx(void) { return fAbsorberThickness_x; }
//------------------------------------------------------------------
double InputReaderGeant::GetAbsThicky(void) { return fAbsorberThickness_y; }
//------------------------------------------------------------------
double InputReaderGeant::GetAbsThickz(void) { return fAbsorberThickness_z; }
//------------------------------------------------------------------
void InputReaderGeant::Clear(void) {
  fEventNumber = -1;
  // fIdentified = false;
  fEnergy_Primary = -1000;
  fRealEnergy_e = -1000;
  fRealEnergy_p = -1000;
  fRealPosition_e = NULL;
  fRealPosition_p = NULL;
  fPositionScat = NULL;
  fPositionAbs = NULL;
  fDirectionScat = NULL;
  fDirectionSource = NULL;
  fPositionSource = NULL;
  fTree = NULL;
  fTree1 = NULL;
  fFile = NULL;
  fRealDirection_scatter = NULL;
  fRealDirection_source = NULL;
  fRealPosition_source = NULL;
  // if (!fRecoClusterPositions->empty()) fRecoClusterPositions->clear();
  // if (!fRecoClusterEnergies->empty()) fRecoClusterEnergies->clear();
  fNumberOfSimulatedEvents = -1;
  fScattererThickness_x = -1000;
  fScattererThickness_y = -1000;
  fScattererThickness_z = -1000;
  fAbsorberThickness_x = -1000;
  fAbsorberThickness_y = -1000;
  fAbsorberThickness_z = -1000;

  fScattererPosition = NULL;
  fAbsorberPosition = NULL;
  fScatPlanePos = NULL;
  fAbsPlanePos = NULL;

  return;
}
//------------------------------------------------------------------
