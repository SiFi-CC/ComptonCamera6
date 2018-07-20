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