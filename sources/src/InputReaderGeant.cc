#include "InputReaderGeant.hh" 

ClassImp(InputReaderGeant);

//------------------------------------------------------------------
InputReaderGeant::InputReaderGeant()
		 :InputReader(){
  cout << "##### Warning in InputReaderGeant constructor!" << endl;
  cout << "You are usinf default constructor." << endl;
}
//------------------------------------------------------------------
InputReaderGeant::InputReaderGeant(TString path)
		 :InputReader(path){

  if(!AccessTree("G4SimulationData_Reconstruction")){
    throw "##### Exception in InputReaderGeant constructor!";
  }
  AccessSetup();
  
}
//------------------------------------------------------------------
InputReaderGeant::~InputReaderGeant(){
  if(fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
bool InputReaderGeant::AccessTree(TString name){
  
  fTree = (TTree*)fFile->Get(name);
  
  if(fTree==NULL){
    cout << "##### Error in InputReaderGeant::AccessTree()!" << endl;
    cout << "Could not access the tree!" << endl;
    return false;
  }
  
  fTree->SetBranchAddress("EventNumber",&fEventNumber);
  fTree->SetBranchAddress("Identified",&fIdentified);
  fTree->SetBranchAddress("RecoEnergy_e",&fRecoEnergy_e);
  fTree->SetBranchAddress("RecoEnergy_p",&fRecoEnergy_p);
  fTree->SetBranchAddress("RecoPosition_e",&fRecoPosition_e);
  fTree->SetBranchAddress("RecoPosition_p",&fRecoPosition_p);
  fTree->SetBranchAddress("RecoDirection_scatter",&fRecoDirection_scatter);
  fTree->SetBranchAddress("RecoClusterPositions",&fRecoClusterPositions);
  fTree->SetBranchAddress("RecoClusterEnergies",&fRecoClusterEnergies);

  cout << "\n\nIn InputReaderGeant::AccessTree()." << endl;
  cout << fTree->GetName() << " tree accessed.\n" << endl;
  
  return true;
}
//------------------------------------------------------------------
bool InputReaderGeant::AccessSetup(void){

  TString fname = fFile->GetName();
  G4Input *input = new G4Input(fname,false);
  
  fScatDimX = input->GetScattererXLength();
  fScatDimY = input->GetScattererYLength();
  fScatDimZ = input->GetScattererZLength();
  fAbsDimX  = input->GetAbsorberXLength();
  fAbsDimY  = input->GetAbsorberYLength();
  fAbsDimZ  = input->GetAbsorberZLength();
  fScatPosition = input->GetScattererPosition();
  fAbsPosition  = input->GetAbsorberPosition();
  
  cout << "\n\n----- In InputReaderGeant::AccessSetup()" << endl;
  cout << "\t Scatterer dimensions: " << fScatDimX << "\t" 
       << fScatDimY << "\t" << fScatDimZ << endl;
  cout << "\t Absorber dimensions:  " << fAbsDimX << "\t"
       << fAbsDimY << "\t" << fAbsDimZ << endl;
  cout << "\t Scatterer position: " << fScatPosition.X() << "\t" 
       <<  fScatPosition.Y() << "\t" << fScatPosition.Z() << endl; 
  cout << "\t Absorber position:  " << fAbsPosition.X() << "\t" 
       <<  fAbsPosition.Y() << "\t" << fAbsPosition.Z() << endl << endl;
  
  delete input;
  
  return true;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetSourcePosition(void){
 return NULL; 
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetScatPosition(void){
  return NULL;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetAbsPosition(void){
  return NULL;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetPrimaryGammaDir(void){
  return NULL;
}
//------------------------------------------------------------------
TVector3* InputReaderGeant::GetScatGammaDir(void){
  return NULL;
}
//------------------------------------------------------------------
double InputReaderGeant::GetEnSource(void){
  return -100;
}
//------------------------------------------------------------------
double InputReaderGeant::GetEnScat(void){
  return -100;
}
//------------------------------------------------------------------
double InputReaderGeant::GetEnAbs(void){
  return -100;
}
//------------------------------------------------------------------