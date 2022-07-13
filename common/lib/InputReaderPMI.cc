#include "InputReaderPMI.hh"
#include <math.h>

ClassImp(InputReaderPMI);

//------------------------------------------------------------------
/// Default constructor.
InputReaderPMI::InputReaderPMI() : InputReader() {
    Clear();
    cout << "##### Warning in InputReaderPMI constructor!" << endl;
    cout << "You are using default constructor." << endl;
} 
//------------------------------------------------------------------
/// Standard constructor.
///\param path (TString) - path to the input file.
InputReaderPMI::InputReaderPMI(TString path) : InputReader(path) {

    if (!AccessTree("CalibratedEvents")) {
        throw "##### Exception in InputReaderPMI constructor!";
    }
}
//------------------------------------------------------------------
/// Default destructor.
InputReaderPMI::~InputReaderPMI() {
  if (fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
/// Accesses data of the tree's branches in the ROOT file.
///\param name (TString) - name of the tree.
bool InputReaderPMI::AccessTree(TString name) {

    fTree = (TTree*)fFile->Get(name);

    if (fTree == NULL) {
        cout << "##### Error in InputReaderPMI::AccessTree()!" << endl;
        cout << "Could not access the tree!" << endl;
        return false;
    }

    fPoint1 = new TVector3();
    fPoint2 = new TVector3();

    int fscafibernumber, fabsfibernumber;
    long long int ftimestampR, ftimestampL, ftimestampAbs;

    fTree->SetBranchAddress("ScaPosition", &fPoint1);
    fTree->SetBranchAddress("AbsPosition", &fPoint2);
    fTree->SetBranchAddress("ScaE", &fEnergy1);
    fTree->SetBranchAddress("AbsE", &fEnergy2);
    fTree->SetBranchAddress("ScaFiberNumber", &fscafibernumber);
    fTree->SetBranchAddress("AbsFiberNumber", &fabsfibernumber);
    fTree->SetBranchAddress("AbsTimeStamp", &ftimestampAbs);
    fTree->SetBranchAddress("ScaTimeStampL", &ftimestampL);
    fTree->SetBranchAddress("ScaTimeStampR", &ftimestampR);
    fTree->SetBranchAddress("AbsClusterSize", &fabsclustersize);
  
    cout << "\n\nIn InputReaderPMI::AccessTree()." << endl;
    cout << fTree->GetName() << " tree accessed.\n" << endl;

    return true;
}
//------------------------------------------------------------------
/// Loads events from the tree to analyze them.
bool InputReaderPMI::LoadEvent(int i) {

    fPositionScat = new TVector3();
    fPositionAbs  = new TVector3();
    fDirectionScat = new TVector3();

    int imax = fTree->GetEntries();

    if (i > imax)
    {
        cout << "##### Error in InputReader::LoadEvent()!" << endl;
        cout << "Requested event number larger than number of events in the tree!" << endl;
        return false;
    }

    fTree->GetEntry(i);

    return true;

}
//------------------------------------------------------------------
TVector3* InputReaderPMI::GetPositionScattering(void) {
    fPositionScat->SetX(fPoint1->X());
    fPositionScat->SetY(fPoint1->Y());
    fPositionScat->SetZ(fPoint1->Z());
    return fPositionScat;
}
//------------------------------------------------------------------
TVector3* InputReaderPMI::GetPositionAbsorption(void) {
    fPositionAbs->SetX(fPoint2->X());
    fPositionAbs->SetY(fPoint2->Y());
    fPositionAbs->SetZ(fPoint2->Z());
    return fPositionAbs;
}
//------------------------------------------------------------------
TVector3* InputReaderPMI::GetGammaDirScattered(void) {
    fDirectionScat->SetX((fPoint2->X()-fPoint1->X())/sqrt(pow(fPoint1->X(), 2)+pow(fPoint2->X(),2)));
    fDirectionScat->SetY((fPoint2->Y()-fPoint1->Y())/sqrt(pow(fPoint1->Y(), 2)+pow(fPoint2->Y(),2)));
    fDirectionScat->SetZ((fPoint2->Z()-fPoint1->Z())/sqrt(pow(fPoint1->Z(), 2)+pow(fPoint2->Z(),2)));
    return fDirectionScat;
}
//------------------------------------------------------------------
double InputReaderPMI::GetEnergyLoss(void) {
    fEnergyLoss = fEnergy1;
    return fEnergyLoss;
}
//------------------------------------------------------------------
double InputReaderPMI::GetEnergyScattered(void) {
    fEnergyScattered = fEnergy2;
    return fEnergyScattered;
}
//------------------------------------------------------------------
void InputReaderPMI::SelectEvents() {  
  fSelectedEvents.clear(); 
  for(int i = 0; i < fTree->GetEntries(); i++) {
  	fTree->GetEntry(i);
 	if(SelectSingleEvent()) fSelectedEvents.push_back(i);
  }
}
//------------------------------------------------------------------
bool InputReaderPMI::SelectSingleEvent() {
    return fabsclustersize < 3;
}
//------------------------------------------------------------------
void InputReaderPMI::Clear(void) {
    fPoint1 = NULL;
    fPoint2 = NULL;
    fFile = NULL;
    fTree = NULL;
    fEnergy1 = -100;
    fEnergy2 = -100;
    return;
}
//------------------------------------------------------------------
