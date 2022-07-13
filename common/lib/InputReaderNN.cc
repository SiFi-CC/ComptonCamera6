#include "InputReaderNN.hh"
#include <vector>
ClassImp(InputReaderNN);

//------------------------------------------------------------------
/// Default constructor.
InputReaderNN::InputReaderNN() : InputReader() {
  Clear();
  cout << "##### Warning in InputReaderNN constructor!" << endl;
  cout << "You are usinf default constructor." << endl;
}
//------------------------------------------------------------------
/// Standard constructor.
///\param path (TString) - path to the input file.
InputReaderNN::InputReaderNN(TString path):
InputReader(path),
fCorrectOnly(false),
fX1(0),
fY1(0),
fZ1(0),
fX2(0),
fY2(0),
fZ2(0),
fX3(0),
fY3(0),
fZ3(0),
fVX(0),
fVY(0),
fVZ(0),
fVUncX(0),
fVUncY(0),
fVUncZ(0),
fPX(0),
fPY(0),
fPZ(0),
fPUncX(0),
fPUncY(0),
fPUncZ(0),
fE0Calc(0),
fE0CalcUnc(0),
fArc(0),
fArcUnc(0),
fE1(0),
fE1Unc(0),
fE2(0),
fE2Unc(0),
fE3(0),
fE3Unc(0),
fClassID(0),
fEventType(0),
fEnergyBinID(0),
fPrimaryEnergy(0),
fEnergyLoss(0), 
fEnergyScattered(0) 
{

  if (!AccessTree("ConeList")) {
    throw "##### Exception in InputReaderNN constructor!";
  }

  fPositionScat = new TVector3();
  fPositionAbs = new TVector3();
  fDirectionScat = new TVector3();
  
}
//------------------------------------------------------------------
/// Default destructor.
InputReaderNN::~InputReaderNN() {
  if (fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
/// Accesses data of trees'branches in ROOT file.
///\param name (TString) - name of tree.
bool InputReaderNN::AccessTree(TString name) {

  fTree = (TTree*)fFile->Get(name);
  if (fTree == NULL) {
    cout << "##### Error in InputReaderNN::AccessTree()!" << endl;
    cout << "Could not access the tree with requested name" << name.Data() << endl;
    return false;
  }
  fTree->SetBranchAddress("x_1",&fX1);
  fTree->SetBranchAddress("y_1",&fY1);
  fTree->SetBranchAddress("z_1",&fZ1);
  fTree->SetBranchAddress("x_2",&fX2);
  fTree->SetBranchAddress("y_2",&fY2);
  fTree->SetBranchAddress("z_2",&fZ2);
  fTree->SetBranchAddress("x_3",&fX3);
  fTree->SetBranchAddress("y_3",&fY3);
  fTree->SetBranchAddress("z_3",&fZ3);
  fTree->SetBranchAddress("v_x",&fVX);
  fTree->SetBranchAddress("v_y",&fVY);
  fTree->SetBranchAddress("v_z",&fVZ);
  fTree->SetBranchAddress("v_unc_x",&fVUncX);
  fTree->SetBranchAddress("v_unc_y",&fVUncY);
  fTree->SetBranchAddress("v_unc_z",&fVUncZ);
  fTree->SetBranchAddress("p_x",&fPX);
  fTree->SetBranchAddress("p_y",&fPY);
  fTree->SetBranchAddress("p_z",&fPZ);
  fTree->SetBranchAddress("p_unc_x",&fPUncX);
  fTree->SetBranchAddress("p_unc_y",&fPUncY);
  fTree->SetBranchAddress("p_unc_z",&fPUncZ);
  fTree->SetBranchAddress("E0Calc",&fE0Calc);
  fTree->SetBranchAddress("E0Calc_unc",&fE0CalcUnc);
  fTree->SetBranchAddress("arc",&fArc);
  fTree->SetBranchAddress("arc_unc",&fArcUnc);
  fTree->SetBranchAddress("E1",&fE1);
  fTree->SetBranchAddress("E1_unc",&fE1Unc);
  fTree->SetBranchAddress("E2",&fE2);
  fTree->SetBranchAddress("E2_unc",&fE2Unc);
  fTree->SetBranchAddress("E3",&fE3);
  fTree->SetBranchAddress("E3_unc",&fE3Unc);
  fTree->SetBranchAddress("ClassID",&fClassID);
  fTree->SetBranchAddress("EventType",&fEventType);
  fTree->SetBranchAddress("EnergyBinID",&fEnergyBinID);

  cout << "\n\nIn InputReaderNN::AccessTree()." << endl;
  cout << fTree->GetName() << " tree accessed.\n" << endl;

  return true;
}

void InputReaderNN::SelectEvents(){  
  for(int i=0;i<fTree->GetEntries();i++) {
  	fTree->GetEntry(i);
	if(fCorrectOnly==1 && fEventType==0)continue;
	else if(fCorrectOnly==2 && fEventType!=2)continue;
	else if(fCorrectOnly==3 && fEventType!=0)continue;
	else if(fCorrectOnly==4 && fEventType!=1)continue;
 	fSelectedEvents.push_back(i);
  }
}
//------------------------------------------------------------------
/// loads events from trees to analyze them in CCMLEM class.
///\param i (int) - number of events
bool InputReaderNN::LoadEvent(int i) {

  
  int imax = fTree->GetEntries();
  if (i > imax) {
    cout << "##### Error in InputReaderNN::LoadEvent() in Event tree!" << endl;
    cout << "Requested event number larger than number of events in the tree!"
         << endl;
    return false;
  }

  fTree->GetEntry(i);
  if(fCorrectOnly==1 && fEventType==0)return false;
  else if(fCorrectOnly==2 && fEventType!=2)return false;
  else if(fCorrectOnly==3 && fEventType!=0)return false;
  else if(fCorrectOnly==4 && fEventType!=1)return false;
  fPositionScat->SetXYZ(fX1,fY1,fZ1);
  fPositionAbs->SetXYZ(fX2,fY2,fZ2);
  fDirectionScat->SetXYZ(fPX,fPY,fPZ);
  //ROTATION BACK TO KRAKOW FRAME
  //INPUTDATA IS IN CS OF LUEBECK
  Rotate(fPositionScat);
  Rotate(fPositionAbs);
  Rotate(fDirectionScat);
  fPrimaryEnergy = fE0Calc;
  fEnergyLoss=fE1; 
  fEnergyScattered=fE2; 

  return true;
}
//------------------------------------------------------------------
void InputReaderNN::Clear(void) {
  fX1=-1000;
  fY1=-1000;
  fZ1=-1000;
  fX2=-1000;
  fY2=-1000;
  fZ2=-1000;
  fX3=-1000;
  fY3=-1000;
  fZ3=-1000;
  fVX=-1000;
  fVY=-1000;
  fVZ=-1000;
  fVUncX=-1000;
  fVUncY=-1000;
  fVUncZ=-1000;
  fPX=-1000;
  fPY=-1000;
  fPZ=-1000;
  fPUncX=-1000;
  fPUncY=-1000;
  fPUncZ=-1000;
  fE0Calc=-1000,
  fE0CalcUnc=-1000;
  fArc=-1000;
  fArcUnc=-1000;
  fE1=-1000;
  fE1Unc=-1000;
  fE2=-1000;
  fE2Unc=-1000;
  fE3=-1000;
  fE3Unc=-1000;
  fClassID=-1000;
  fEventType=-1000;
  fEnergyBinID=-1000;

  fPrimaryEnergy=-1000;
  fEnergyLoss=-1000;
  fEnergyScattered=-1000;

  fPositionScat = NULL;
  fPositionAbs = NULL;
  fDirectionScat = NULL;

  return;
}
//------------------------------------------------------------------
void InputReaderNN::Rotate(TVector3* vec){
	vec->RotateZ(-TMath::Pi()/2);
	vec->RotateY(TMath::Pi());
	vec->RotateY(-TMath::Pi()/2);
	vec->RotateZ(-TMath::Pi());
}
