#include "InputReader.hh"

ClassImp(InputReader);

//------------------------------------------------------------------
InputReader::InputReader(){
  cout << "##### Warning in InputReader constructor!" << endl;
  cout << "You are using default constructor. Set the input file!" << endl;
  Clear();
}
//------------------------------------------------------------------
InputReader::~InputReader(){
  if(fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
InputReader::InputReader(TString path, bool verbose){
  if(!SetInputFile(path)){
    throw "##### Error in InputReader constructor!";
  }
  SetVerbose(verbose);
}
//------------------------------------------------------------------
bool InputReader::SetInputFile(TString path){
  
  fFile = new TFile(path,"READ");
  if(!fFile->IsOpen()){
   cout << "##### Error in InputReader::SetInputFile()!" << endl;
   cout << "Could not open the file!" << endl;
   return false;
  }
  
  return true;
}
//------------------------------------------------------------------
bool InputReader::AccessTree(TString name){
  
  fTree = (TTree*)fFile->Get(name);
  if(fTree==NULL){
   cout << "##### Error in InputReader::AccessTree()!" << endl;
   cout << "Could not access tree: " << name << endl;
   return false;
  }
  
  return true;
}
//------------------------------------------------------------------
bool InputReader::LoadEvent(int i){
  
  if(fVerbose)
    cout << "Reading event number " << i << endl;
  
  int imax = fTree->GetEntries();
  if(i>imax){
   cout << "##### Error in InputReader::LoadEvent()!" << endl;
   cout << "Requested event number larger than number of events in the tree!" << endl;
   return false;
  }
  
  fTree->GetEntry(i);
  return true; 
}
//------------------------------------------------------------------
TVector3* InputReader::GetSourcePosition(void){
  TVector3 *vec = new TVector3(-100,-100,-100);
  return vec;
}
//------------------------------------------------------------------
TVector3* InputReader::GetScatPosition(void){
  TVector3 *vec = new TVector3(-100,-100,-100);
  return vec;
}
//------------------------------------------------------------------
TVector3* InputReader::GetAbsPosition(void){
  TVector3 *vec = new TVector3(-100,-100,-100);
  return vec;
}
//------------------------------------------------------------------
TVector3* InputReader::GetPrimaryGammaDir(void){
  TVector3 *vec = new TVector3(-100,-100,-100);
  return vec;
}
//------------------------------------------------------------------
TVector3* InputReader::GetScatGammaDir(void){
  TVector3 *vec = new TVector3(-100,-100,-100);
  return vec;
}
//------------------------------------------------------------------
double InputReader::GetEnSource(void){
  return -100;
}
//------------------------------------------------------------------
double InputReader::GetEnScat(void){
  return -100;
}
//------------------------------------------------------------------
double InputReader::GetEnAbs(void){
  return -100;
}
//------------------------------------------------------------------
void InputReader::Clear(void){
 fVerbose  = false;
 fFile     = NULL;
 fTree     = NULL;
}
//------------------------------------------------------------------
void InputReader::Print(void){
  cout << "\n-------------------------------------------------------" << endl;
  cout << "This is Print() for InputReader class object" << endl;
  if(fFile!=NULL){
    cout << "Opened file: " << fFile->GetName() << endl;
    cout << "Verbose level: " << fVerbose << endl;
  }
  else{
    cout << "It's empty!" << endl;
  }
  cout << "-------------------------------------------------------\n" << endl;
} 
//------------------------------------------------------------------