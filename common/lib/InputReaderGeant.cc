#include "InputReaderGeant.hh"
#include <vector>
ClassImp(InputReaderGeant);

//------------------------------------------------------------------
/// Default constructor.
InputReaderGeant::InputReaderGeant() : InputReader() {
  Clear();
  cout << "##### Warning in InputReaderGeant constructor!" << endl;
  cout << "You are using default constructor." << endl;
}
//------------------------------------------------------------------
/// Standard constructor.
///\param path (TString) - path to the input file.
InputReaderGeant::InputReaderGeant(TString path) : InputReader(path) {

  if (!AccessTree()) {
    throw "##### Exception in InputReaderGeant constructor!";
  }

  // Load setup
  fTreeSetup->GetEntry(0);

  fPositionScat = new TVector3();
  fPositionAbs = new TVector3();
  fDirectionScat = new TVector3();
  
  fPositionScatReco = new TVector3();
  fPositionAbsReco = new TVector3();
  fDirectionScatReco = new TVector3();
  
  
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
/// Accesses data of trees' branches in ROOT file.
bool InputReaderGeant::AccessTree() {

  // Check structure of file: 
  // 3 trees: G4SimulationData_Setup+G4SimulationData_Real+G4SimulationData_Reconstruction
  // or 2 trees: Setup+Events
  if (fFile->Get("G4SimulationData_Setup"))  fJointTree = false;
  else if (fFile->Get("Setup")) fJointTree = true;
  else {
    cout << "##### Error in InputReaderGeant::AccessTree()!" << endl;
    return false;
  }

  if(fJointTree){ 
    fTree = (TTree*)fFile->Get("Events");
    fTreeSetup = (TTree*)fFile->Get("Setup");
    if (fTree == NULL || fTreeSetup == NULL) {
      cout << "##### Error in InputReaderGeant::AccessTree()!" << endl;
      cout << "Could not access tree!" << endl;
      return false;
    }
  } else {
    fTree = (TTree*)fFile->Get("G4SimulationData_Real");
    fTreeReco = (TTree*)fFile->Get("G4SimulationData_Reconstruction");
    fTreeSetup = (TTree*)fFile->Get("G4SimulationData_Setup");
    if (fTree == NULL || fTreeReco == NULL || fTreeSetup == NULL) {
      cout << "##### Error in InputReaderGeant::AccessTree()!" << endl;
      cout << "Could not access tree!" << endl;
      return false;
    }
  }
  
  fRecoEnergy_e = new PhysicVar();
  fRecoEnergy_p = new PhysicVar();
  fRecoPosition_e = new PhysicVec();
  fRecoPosition_p = new PhysicVec();
  fRecoDirection_scatter = new PhysicVec();
  fRecoClusterPositions = 0;
  fRecoClusterEnergies = 0;
  
  fRealPosition_source = new TVector3();
  fRealDirection_source = new TVector3();
  fRealPosition_e = new TVector3();
  fRealPosition_p = new TVector3();
  fRealDirection_scatter = new TVector3();
  
  fScattererPosition = new TVector3();
  fAbsorberPosition = new TVector3();

  fTreeSetup->SetBranchAddress("ScattererThickness_x", &fScattererThickness_x);
  fTreeSetup->SetBranchAddress("ScattererThickness_y", &fScattererThickness_y);
  fTreeSetup->SetBranchAddress("ScattererThickness_z", &fScattererThickness_z);
  fTreeSetup->SetBranchAddress("AbsorberThickness_x", &fAbsorberThickness_x);
  fTreeSetup->SetBranchAddress("AbsorberThickness_y", &fAbsorberThickness_y);
  fTreeSetup->SetBranchAddress("AbsorberThickness_z", &fAbsorberThickness_z);
  fTreeSetup->SetBranchAddress("ScattererPosition", &fScattererPosition);
  fTreeSetup->SetBranchAddress("AbsorberPosition", &fAbsorberPosition);
  fTreeSetup->SetBranchAddress("NumberOfSimulatedEvents",
                           &fNumberOfSimulatedEvents);

  fTree->SetBranchAddress("EventNumber", &fEventNumber);
  fTree->SetBranchAddress("Energy_Primary", &fEnergy_Primary);
  fTree->SetBranchAddress("RealEnergy_e", &fRealEnergy_e);
  fTree->SetBranchAddress("RealEnergy_p", &fRealEnergy_p);
  fTree->SetBranchAddress("RealPosition_source", &fRealPosition_source);
  fTree->SetBranchAddress("RealDirection_source", &fRealDirection_source);
  fTree->SetBranchAddress("RealPosition_e", &fRealPosition_e);
  fTree->SetBranchAddress("RealPosition_p", &fRealPosition_p);
  fTree->SetBranchAddress("RealDirection_scatter", &fRealDirection_scatter);

 
  if(fJointTree){ 
    fTree->SetBranchAddress("Identified", &fIIdentified);
    fTree->SetBranchAddress("RecoEnergy_e", &fRecoEnergy_e);
    fTree->SetBranchAddress("RecoEnergy_p", &fRecoEnergy_p);
    fTree->SetBranchAddress("RecoPosition_e", &fRecoPosition_e);
    fTree->SetBranchAddress("RecoPosition_p", &fRecoPosition_p);
    fTree->SetBranchAddress("RecoDirection_scatter", &fRecoDirection_scatter);
    fTree->SetBranchAddress("RecoClusterPositions", &fRecoClusterPositions);
    fTree->SetBranchAddress("RecoClusterEnergies", &fRecoClusterEnergies);

  } else {
    fTreeReco->SetBranchAddress("EventNumber", &fEventNumberReco);
    fTreeReco->SetBranchAddress("Identified", &fBIdentified);
    fTreeReco->SetBranchAddress("RecoEnergy_e", &fRecoEnergy_e);
    fTreeReco->SetBranchAddress("RecoEnergy_p", &fRecoEnergy_p);
    fTreeReco->SetBranchAddress("RecoPosition_e", &fRecoPosition_e);
    fTreeReco->SetBranchAddress("RecoPosition_p", &fRecoPosition_p);
    fTreeReco->SetBranchAddress("RecoDirection_scatter", &fRecoDirection_scatter);
    fTreeReco->SetBranchAddress("RecoClusterPositions", &fRecoClusterPositions);
    fTreeReco->SetBranchAddress("RecoClusterEnergies", &fRecoClusterEnergies);
  }

  return true;
}
//------------------------------------------------------------------
/// loads events from trees to analyze them in CCMLEM class.
///\param i (int) - index of event
bool InputReaderGeant::LoadEvent(int i) {

  if (i > fTree->GetEntries()) {
    cout << "##### Error in InputReaderGeant::LoadEvent()" << endl;
    cout << "Requested event number larger than number of events in the tree "
         << fTree->GetName() << endl;
    return false;
  }

  fTree->GetEntry(i);
  
  if (!fJointTree) {
    if (i > fTreeReco->GetEntries()) {
      cout << "##### Error in InputReaderGeant::LoadEvent()" << endl;
      cout << "Requested event number larger than number of events in the tree "
           << fTreeReco->GetName() << endl;
      return false;
    }

    fTreeReco->GetEntry(i);
  }

  // Check if event contains valid data
  if ( fUseRealInformation && 0 == fRealEnergy_e ) return false;
  if ( !fUseRealInformation && (-11 == fIIdentified || 0 == fIIdentified )) return false;

  return true;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionPrimary(void) {
  fPositionSource->SetX(fRealPosition_source->X());
  fPositionSource->SetY(fRealPosition_source->Y());
  fPositionSource->SetZ(fRealPosition_source->Z());
  return fPositionSource;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionScattering(void) {
  if(fUseRealInformation){
    fPositionScat->SetX(fRealPosition_e->X());
    fPositionScat->SetY(fRealPosition_e->Y());
    fPositionScat->SetZ(fRealPosition_e->Z());
  } else {
    fPositionScat->SetX(fRecoPosition_e->position.X());
    fPositionScat->SetY(fRecoPosition_e->position.Y());
    fPositionScat->SetZ(fRecoPosition_e->position.Z());
  }
  return fPositionScat;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionScatteringReal(void) {
  fPositionScat->SetX(fRealPosition_e->X());
  fPositionScat->SetY(fRealPosition_e->Y());
  fPositionScat->SetZ(fRealPosition_e->Z());
  return fPositionScat;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionScatteringReco(void) {
  fPositionScatReco->SetX(fRecoPosition_e->position.X());
  fPositionScatReco->SetY(fRecoPosition_e->position.Y());
  fPositionScatReco->SetZ(fRecoPosition_e->position.Z());
  return fPositionScatReco;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionAbsorption(void) {
  if(fUseRealInformation){
    fPositionAbs->SetX(fRealPosition_p->X());
    fPositionAbs->SetY(fRealPosition_p->Y());
    fPositionAbs->SetZ(fRealPosition_p->Z());
  } else {
    fPositionAbs->SetX(fRecoPosition_p->position.X());
    fPositionAbs->SetY(fRecoPosition_p->position.Y());
    fPositionAbs->SetZ(fRecoPosition_p->position.Z());
  }
  return fPositionAbs;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionAbsorptionReal(void) {
  fPositionAbs->SetX(fRealPosition_p->X());
  fPositionAbs->SetY(fRealPosition_p->Y());
  fPositionAbs->SetZ(fRealPosition_p->Z());
  return fPositionAbs;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPositionAbsorptionReco(void) {
  fPositionAbsReco->SetX(fRecoPosition_p->position.X());
  fPositionAbsReco->SetY(fRecoPosition_p->position.Y());
  fPositionAbsReco->SetZ(fRecoPosition_p->position.Z());
  return fPositionAbsReco;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetGammaDirPrimary(void) {
  fDirectionSource->SetX(fRealDirection_source->X());
  fDirectionSource->SetY(fRealDirection_source->Y());
  fDirectionSource->SetZ(fRealDirection_source->Z());
  return fDirectionSource;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetGammaDirScattered(void) {
  if(fUseRealInformation){
    fDirectionScat->SetX(fRealDirection_scatter->X());
    fDirectionScat->SetY(fRealDirection_scatter->Y());
    fDirectionScat->SetZ(fRealDirection_scatter->Z());
  } else {
    fDirectionScat->SetX(fRecoDirection_scatter->position.X());
    fDirectionScat->SetY(fRecoDirection_scatter->position.Y());
    fDirectionScat->SetZ(fRecoDirection_scatter->position.Z());
  }
  return fDirectionScat;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetGammaDirScatteredReal(void) {
  fDirectionScat->SetX(fRealDirection_scatter->X());
  fDirectionScat->SetY(fRealDirection_scatter->Y());
  fDirectionScat->SetZ(fRealDirection_scatter->Z());
  return fDirectionScat;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetGammaDirScatteredReco(void) {
  fDirectionScatReco->SetX(fRecoDirection_scatter->position.X());
  fDirectionScatReco->SetY(fRecoDirection_scatter->position.Y());
  fDirectionScatReco->SetZ(fRecoDirection_scatter->position.Z());
  return fDirectionScatReco;
}
//------------------------------------------------------------------
int InputReaderGeant::GetRecoClusterPosSize(void) { return fRecoClusterPositions->size(); }
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
  return fEnergy_Primary;
}
//------------------------------------------------------------------
double InputReaderGeant::GetEnergyLoss(void) { 
  if(fUseRealInformation)
    return fRealEnergy_e; 
  else
    return fRecoEnergy_e->value;
}
//------------------------------------------------------------------
double InputReaderGeant::GetEnergyLossReal(void) { return fRealEnergy_e; }
//------------------------------------------------------------------
double InputReaderGeant::GetEnergyLossReco(void) { return fRecoEnergy_e->value; }
//------------------------------------------------------------------
double InputReaderGeant::GetEnergyScattered(void) {
  if(fUseRealInformation)
    return fRealEnergy_p; 
  else
    return fRecoEnergy_p->value;
}
//------------------------------------------------------------------
double InputReaderGeant::GetEnergyScatteredReal(void) { return fRealEnergy_p; }
//------------------------------------------------------------------
double InputReaderGeant::GetEnergyScatteredReco(void) { return fRecoEnergy_p->value; }
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
  fEventNumberReco = -1;
  fBIdentified = false;
  fIIdentified = 0;
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
  fTreeReco = NULL;
  fTreeSetup = NULL;
  fFile = NULL;
  fRealDirection_scatter = NULL;
  fRealDirection_source = NULL;
  fRealPosition_source = NULL;
  
  fRecoEnergy_e = NULL;
  fRecoEnergy_p = NULL;
  fRecoPosition_e = NULL;
  fRecoPosition_p = NULL;
  fPositionScatReco = NULL;
  fPositionAbsReco = NULL;
  fDirectionScatReco = NULL;
  
  if (!fRecoClusterPositions->empty()) fRecoClusterPositions->clear();
  if (!fRecoClusterEnergies->empty()) fRecoClusterEnergies->clear();
  
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

  fUseRealInformation = true;

  return;
}
//------------------------------------------------------------------
