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

  if (!AccessTree("Pos&EnergyRecoClus"/*, "Reco", "Real"*/)) {
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
/// Accesses data of trees'branches in ROOT file.
///\param name (TString) - name of tree.
bool InputReaderEI::AccessTree(TString name/*, TString name1, TString name2*/) {

  fTree = (TTree*)fFile->Get(name);
//   fTree1 = (TTree*)fFile->Get(name1);
//   fTree2 = (TTree*)fFile->Get(name2);
  
  if (fTree == NULL) {
    cout << "##### Error in InputReaderEI::AccessTree()!" << endl;
    cout << "Could not access the tree!" << endl;
    return false;
  }
//   if (fTree1 == NULL) {
//     cout << "##### Error in InputReaderEI::AccessTree()!" << endl;
//     cout << "Could not access the tree!" << endl;
//     return false;
//   }
//   if (fTree2 == NULL) {
//     cout << "##### Error in InputReaderEI::AccessTree()!" << endl;
//     cout << "Could not access the tree!" << endl;
//     return false;
//   }
  
  fEnergyReco0 = new PhysicVar();
  fEnergyReco1 = new PhysicVar();
  fPosScatClus = new PhysicVec();
  fPosAbsClus = new PhysicVec();
  
//   fPosScatReco = new TVector3();
//   fPosAbsReco = new TVector3();
//   
//   fPosScatReal = new TVector3();
//   fPosAbsReal = new TVector3();
  fTree->SetBranchAddress("EventNumber", &fEventNumber);
  fTree->SetBranchAddress("PosScat", &fPosScatClus);
  fTree->SetBranchAddress("EnergyReco0", &fEnergyReco0);
  fTree->SetBranchAddress("EnergyRe0", &fEnergyRe0);
  fTree->SetBranchAddress("PosAbs", &fPosAbsClus);
  fTree->SetBranchAddress("EnergyReco1", &fEnergyReco1);
  fTree->SetBranchAddress("EnergyRe1", &fEnergyRe1);
  
//   fTree1 = new TTree("Reco", "Reco");
//   fTree1->Branch("PosScatReco", &fPosScatReco);
//   fTree1->Branch("PosAbsReco", &fPosAbsReco);
//   fTree1->Branch("RecoEnergy_e", &fRecoEnergy_e);
//   fTree1->Branch("RecoEnergy_p", &fRecoEnergy_p);
//      
//      
//   fTree2 = new TTree("Real", "Real");
//   fTree2->Branch("PosScatReal", &fPosScatReal);
//   fTree2->Branch("PosAbsReal", &fPosAbsReal);
//   fTree2->Branch("RealEnergy_e", &fRealEnergy_e);
//   fTree2->Branch("RealEnergy_p", &fRealEnergy_p);
  
  cout << "\n\nIn InputReaderEI::AccessTree()." << endl;
  cout << fTree->GetName() << " tree accessed.\n" << endl;

//   cout << "\n\nIn InputReaderEI::AccessTree()." << endl;
//   cout << fTree1->GetName() << " tree accessed.\n" << endl;
//   
//   cout << "\n\nIn InputReaderEI::AccessTree()." << endl;
//   cout << fTree2->GetName() << " tree accessed.\n" << endl;
 
  
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
  
//   int imax1 = fTree1->GetEntries();
//   if (i > imax1) {
//     cout << "##### Error in InputReaderEI::LoadEvent() in reconstruction tree!" << endl;
//     cout << "Requested event number larger than number of events in the tree!"
//          << endl;
//     return false;
//   }
//   fTree1->GetEntry(i);
//   
//   int imax2 = fTree2->GetEntries();
//   if (i > imax2) {
//     cout << "##### Error in InputReaderEI::LoadEvent() in reconstruction tree!" << endl;
//     cout << "Requested event number larger than number of events in the tree!"
//          << endl;
//     return false;
//   }
//   fTree2->GetEntry(i);

  return true;
}
//------------------------------------------------------------------
TVector3* InputReaderEI::GetPositionScattering(void)
{
  fPositionScatReco->SetX(fPosScatClus->position.X());
  fPositionScatReco->SetY(fPosScatClus->position.Y());
  fPositionScatReco->SetZ(fPosScatClus->position.Z());
  return fPositionScatReco;
}
//------------------------------------------------------------------
TVector3* InputReaderEI::GetPositionAbsorption(void) {
  fPositionAbsReco->SetX(fPosAbsClus->position.X());
  fPositionAbsReco->SetY(fPosAbsClus->position.Y());
  fPositionAbsReco->SetZ(fPosAbsClus->position.Z());
  return fPositionAbsReco;
}
//------------------------------------------------------------------
double InputReaderEI::GetEnergyLoss(void) { return fEnergyReco0->value; }
//------------------------------------------------------------------
double InputReaderEI::GetEnergyScattered(void) { return fEnergyReco1->value; }
//------------------------------------------------------------------
void InputReaderEI::Clear(void) {
  
  //fTotalEnergy2 = -1000;
  fEventNumber = -1;
  fPositionScatReco = NULL;
  fPositionAbsReco = NULL;
  
  fPosScatClus = NULL;
  fPosAbsClus = NULL;
  fEnergyReco0 = NULL;
  fEnergyReco1 = NULL;
  fEnergyRe0 = -1000;
  fEnergyRe1 = -1000;
  
//   fPosScatReco = NULL;
//   fPosAbsReco = NULL;
//   fRecoEnergy_e = -1000;
//   fRecoEnergy_p = -1000;
//   
//   fPosScatReal = NULL;
//   fPosAbsReal = NULL;
//   fRealEnergy_e = -1000;
//   fRealEnergy_p = -1000;
//   
  fTree = NULL;
//   fTree1 = NULL;
//   fTree2 = NULL;
  
  fFile = NULL;

  return;
}
//------------------------------------------------------------------
