#include "InputReaderSimple.hh"

ClassImp(InputReaderSimple);

//------------------------------------------------------------------
/// Default constructor.
InputReaderSimple::InputReaderSimple() : InputReader() {
  Clear();
  cout << "##### Warning in InputReaderSimple constructor!" << endl;
  cout << "You are using default constructor." << endl;
}
//------------------------------------------------------------------
/// Standard constructor (recommended).
///\param path (TString) - path to the input file.
InputReaderSimple::InputReaderSimple(TString path) : InputReader(path) {

  if (!AccessTree("data")) {
    throw "##### Exception in InputReaderSimple constructor!";
  }
}
//------------------------------------------------------------------
/// Default destructor.
InputReaderSimple::~InputReaderSimple() {
  if (fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
bool InputReaderSimple::AccessTree(TString name) {

  fTree = (TTree*)fFile->Get(name);

  if (fTree == NULL) {
    cout << "##### Error in InputReaderSimple::AccessTree()!" << endl;
    cout << "Could not access the tree!" << endl;
    return false;
  }

  fPoint0 = new TVector3();
  fPoint1 = new TVector3();
  fPoint2 = new TVector3();
  fVersor1 = new TVector3();
  fVersor2 = new TVector3();

  fTree->SetBranchAddress("point0", &fPoint0);
  fTree->SetBranchAddress("point1", &fPoint1);
  fTree->SetBranchAddress("point2", &fPoint2);
  fTree->SetBranchAddress("versor1", &fVersor1);
  fTree->SetBranchAddress("versor2", &fVersor2);
  fTree->SetBranchAddress("energy0", &fEnergy0);
  fTree->SetBranchAddress("energy1", &fEnergy1);
  fTree->SetBranchAddress("energy2", &fEnergy2);

  cout << "\n\nIn InputReaderSimple::AccessTree()." << endl;
  cout << fTree->GetName() << " tree accessed.\n" << endl;

  return true;
}
bool InputReaderSimple::LoadEvent(int i) {

  int imax = fTree->GetEntries();

  if (i > imax) {
    cout << "##### Error in InputReader::LoadEvent()!" << endl;
    cout << "Requested event number larger than number of events in the tree!"
         << endl;
    return false;
  }

  fTree->GetEntry(i);
  fPositionScat->SetXYZ(fPoint1->X(), fPoint1->Y(), fPoint1->Z());
  fPositionAbs->SetXYZ(fPoint2->X(), fPoint2->Y(), fPoint2->Z());
  fDirectionScat->SetXYZ(fVersor2->X(), fVersor2->Y(), fVersor2->Z());
  fPrimaryEnergy = fEnergy0;
  fEnergyLoss = fEnergy1;
  fEnergyScattered = fEnergy2;
  if (fSmear) {
    fPositionScat->SetXYZ(SmearBox(fPositionScat->X(), fResolutionX),
                          SmearGaus(fPositionScat->Y(), fResolutionY),
                          SmearBox(fPositionScat->Z(), fResolutionZ));
    fPositionAbs->SetXYZ(SmearBox(fPositionAbs->X(), fResolutionX),
                         SmearGaus(fPositionAbs->Y(), fResolutionY),
                         SmearBox(fPositionAbs->Z(), fResolutionZ));
    fEnergyLoss = SmearGaus(fEnergyLoss, GetSigmaE(fEnergyLoss));
    fEnergyScattered = SmearGaus(fEnergyScattered, GetSigmaE(fEnergyScattered));
  }
  return true;
}

//------------------------------------------------------------------
void InputReaderSimple::Clear(void) {
  fPoint0 = NULL;
  fPoint1 = NULL;
  fPoint2 = NULL;
  fVersor1 = NULL;
  fVersor2 = NULL;
  fFile = NULL;
  fTree = NULL;
  fEnergy0 = -100;
  fEnergy1 = -100;
  fEnergy2 = -100;
  return;
}
void InputReaderSimple::SetSmearing(bool smear, Double_t posX, Double_t posY,
                                    Double_t posZ) {
  /// To smear camera performance(simple simulation), this file is called for
  /// energy resolution. It shows a function of energy deposited in a 10-cm-long
  /// LuAG(Ce) fiber
  ///  with a square cross-section from Geant4 simulation
  fSmear = smear;
  TString path = gSystem->Getenv("CC6DIR");
  TString name =
      path + "/share/ComptonCamera6/mlem_reco/EnergyResolutionExample.root";
  TFile* file = new TFile(name, "READ");
  fHisto = (TH1D*)file->Get("Scintillator_0ResolutionCombiEnergy");
  TF1* func = new TF1("fit1", "[0]+[1]/sqrt(x)+[2]/x^(3/2)", 0, 4.5);

  fHisto->Fit("fit1", "r");

  fResolutionX = posX;
  fResolutionY = posY;
  fResolutionZ = posZ;
}
//------------------------------------------------------------------
/// Gausian function to smear position resolution along y axis.
/// returns a double value from gausian distribution with the given mean and
/// sigma. \param val (double) - the given position value. \param sigma (double)
///- the given sigma value.
Double_t InputReaderSimple::SmearGaus(double val, double sigma) {
  return gRandom->Gaus(val, sigma);
}
//------------------------------------
/// Returns a double value from Uniform function with respect to position
/// resolution value. \param x (double) - the given position value.
Double_t InputReaderSimple::SmearBox(double x, double resolution) {
  return gRandom->Uniform(x - (resolution / 2), x + (resolution / 2));
}
//------------------------------------
/// Returns the sigma value from the fitting function of deposited energy plot.
///\param energy (double) - the given energy value.
Double_t InputReaderSimple::GetSigmaE(double energy) {
  double sigma = fHisto->GetFunction("fit1")->Eval(energy) * energy;
  return sigma;
}
