#include "InputReaderPMIDec2021.hh"
#include <math.h>

ClassImp(InputReaderPMIDec2021);

//------------------------------------------------------------------
/// Default constructor.
InputReaderPMIDec2021::InputReaderPMIDec2021() : InputReader() {
    Clear();
    cout << "##### Warning in InputReaderPMIDec2021 constructor!" << endl;
    cout << "You are using default constructor." << endl;
} 
//------------------------------------------------------------------
/// Standard constructor.
///\param path (TString) - path to the input file.
InputReaderPMIDec2021::InputReaderPMIDec2021(TString path) : InputReader(path) {

    if (!AccessTree("CalibratedEventsDec2021") ) {
        throw "##### Exception in InputReaderPMIDec2021 constructor!";
    }
}
//------------------------------------------------------------------
/// Default destructor.
InputReaderPMIDec2021::~InputReaderPMIDec2021() {
  if (fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
/// Accesses data of the tree's branches in the ROOT file.
///\param name (TString) - name of the tree.
bool InputReaderPMIDec2021::AccessTree(TString name) {

    fTree = (TTree*)fFile->Get(name);

    if (fTree == NULL) {
        cout << "##### Error in InputReaderPMIDec2021::AccessTree()!" << endl;
        cout << "Could not access the tree!" << endl;
        return false;
    }

    fPoint1 = new TVector3();
    fPoint2 = new TVector3();


    fTree->SetBranchAddress("ScaPosition", &fPoint1);
    fTree->SetBranchAddress("AbsPosition", &fPoint2);
    fTree->SetBranchAddress("ScaE", &fEnergy1);
    fTree->SetBranchAddress("AbsE", &fEnergy2);
    //fTree->SetBranchAddress("ScaFiberNumber", &fscafibernumber);
    //fTree->SetBranchAddress("AbsNeedleNumber", &fabsfibernumber);
    fTree->SetBranchAddress("AbsTimeStamp", &fTimeStampAbs);
    fTree->SetBranchAddress("ScaTimeStamp", &fTimeStampSca);
//    fTree->SetBranchAddress("AbsClusterSize", &fabsclustersize);
  
    cout << "\n\nIn InputReaderPMIDec2021::AccessTree()." << endl;
    cout << fTree->GetName() << " tree accessed.\n" << endl;

    return true;
}
//------------------------------------------------------------------
/// Loads events from the tree to analyze them.
bool InputReaderPMIDec2021::LoadEvent(int i) {

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
TVector3* InputReaderPMIDec2021::GetPositionScattering(void) {
    fPositionScat->SetX(fPoint1->Z());
    fPositionScat->SetY(fPoint1->Y());
    fPositionScat->SetZ(-fPoint1->X());
    return fPositionScat;
}
//------------------------------------------------------------------
TVector3* InputReaderPMIDec2021::GetPositionAbsorption(void) {
    fPositionAbs->SetX(fPoint2->Z());
    fPositionAbs->SetY(fPoint2->Y());
    fPositionAbs->SetZ(-fPoint2->X());
    return fPositionAbs;
}
//------------------------------------------------------------------
TVector3* InputReaderPMIDec2021::GetGammaDirScattered(void) {
    fDirectionScat->SetX((-fPoint2->Z()-(-fPoint1->Z()))/sqrt(pow(-fPoint1->Z(), 2)+pow(-fPoint2->Z(),2)));
    fDirectionScat->SetY((fPoint2->Y()-fPoint1->Y())/sqrt(pow(fPoint1->Y(), 2)+pow(fPoint2->Y(),2)));
    fDirectionScat->SetZ((fPoint2->X()-fPoint1->X())/sqrt(pow(fPoint1->X(), 2)+pow(fPoint2->X(),2)));
    return fDirectionScat;
}
//------------------------------------------------------------------
double InputReaderPMIDec2021::GetEnergyLoss(void) {
    fEnergyLoss = fEnergy1;
    return fEnergyLoss;
}
//------------------------------------------------------------------
double InputReaderPMIDec2021::GetEnergyScattered(void) {
    fEnergyScattered = fEnergy2;
    return fEnergyScattered;
}
//------------------------------------------------------------------
void InputReaderPMIDec2021::SelectEvents() {  
  fSelectedEvents.clear();
  for(int i = 0; i < fTree->GetEntries(); i++) {
	fTree->GetEntry(i);
 	if(SelectSingleEvent()) fSelectedEvents.push_back(i);
  }
}
//------------------------------------------------------------------
bool InputReaderPMIDec2021::SelectSingleEvent() {
	if(3000<(fTimeStampAbs-fTimeStampSca) && (fTimeStampSca-fTimeStampAbs)<7000 && fEnergy1+fEnergy2>1100 && fEnergy1+fEnergy2<1450 && fEnergy1<450 && fPoint1->y()!=50 && fPoint1->y()!=-50) return true;
	//if(1100 <( fEnergy1+fEnergy2) && (fEnergy1+fEnergy2)<1300) return true;
	else return false;    
	//return true;
}
//------------------------------------------------------------------
void InputReaderPMIDec2021::Clear(void) {
    fPoint1 = NULL;
    fPoint2 = NULL;
    fFile = NULL;
    fTree = NULL;
    fEnergy1 = -100;
    fEnergy2 = -100;
    return;
}
//------------------------------------------------------------------
