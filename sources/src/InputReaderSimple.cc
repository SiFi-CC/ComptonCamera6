#include "InputReaderSimple.hh"

ClassImp(InputReaderSimple);

//------------------------------------------------------------------
InputReaderSimple::InputReaderSimple()
                     :InputReader(){
  cout << "##### Warning in InputReaderSimple constructor!" << endl;
  cout << "You are using default constructor." << endl;
}
//------------------------------------------------------------------
InputReaderSimple::InputReaderSimple(TString path)
                     :InputReader(path,true){

  if(!AccessTree()){
   throw "##### Error in InputReaderSimple constructor!"; 
  }
}
//------------------------------------------------------------------
InputReaderSimple::~InputReaderSimple(){
  if(fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
bool InputReaderSimple::AccessTree(void){
 
  fTree = (TTree*)fFile->Get("data");
  
  if(fTree==NULL){
    cout << "##### Error in InputReaderSimple::AccessTree()!" << endl;
    cout << "Could not access the tree!" << endl;
    return false;
  }
  
  fPoint0  = new TVector3();
  fPoint1  = new TVector3();
  fPoint2  = new TVector3();
  fVersor1 = new TVector3();
  fVersor2 = new TVector3();
  
  fTree->SetBranchAddress("point0",&fPoint0);
  fTree->SetBranchAddress("point1",&fPoint1);
  fTree->SetBranchAddress("point2",&fPoint2);
  fTree->SetBranchAddress("versor1",&fVersor1);
  fTree->SetBranchAddress("versor2",&fVersor2);
  fTree->SetBranchAddress("energy0",&fEnergy0);
  fTree->SetBranchAddress("energy1",&fEnergy1);
  fTree->SetBranchAddress("energy2",&fEnergy2);
  
  int Nentries = fTree->GetEntries();
 
  cout << "\n\nIn InputReaderSimple::AccessTree()." << endl;
  cout << fTree->GetName() << " tree accessed.\n" << endl;
  
  return true;
}
//------------------------------------------------------------------
TVector3* InputReaderSimple::GetSourcePosition(void){
  return fPoint0;
}
//------------------------------------------------------------------
TVector3* InputReaderSimple::GetScatPosition(void){
  return fPoint1;
}
//------------------------------------------------------------------
TVector3* InputReaderSimple::GetAbsPosition(void){
  return fPoint2;
}
//------------------------------------------------------------------
TVector3* InputReaderSimple::GetPrimaryGammaDir(void){
  return fVersor1;
}
//------------------------------------------------------------------
TVector3* InputReaderSimple::GetScatGammaDir(void){
  return fVersor2;
}
//------------------------------------------------------------------
double InputReaderSimple::GetEnSource(void){
  return fEnergy0;
}
//------------------------------------------------------------------
double InputReaderSimple::GetEnScat(void){
  return fEnergy1;
}
//------------------------------------------------------------------
double InputReaderSimple::GetEnAbs(void){
  return fEnergy2;
}
//------------------------------------------------------------------