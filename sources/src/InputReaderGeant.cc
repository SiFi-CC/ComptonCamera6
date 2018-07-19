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
}
//------------------------------------------------------------------