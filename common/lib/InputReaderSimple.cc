#include "InputReaderSimple.hh"

ClassImp(InputReaderSimple);

//------------------------------------------------------------------
/// Default constructor.
InputReaderSimple::InputReaderSimple() : InputReader()
{
    Clear();
    cout << "##### Warning in InputReaderSimple constructor!" << endl;
    cout << "You are using default constructor." << endl;
}
//------------------------------------------------------------------
/// Standard constructor (recommended).
///\param path (TString) - path to the input file.
InputReaderSimple::InputReaderSimple(TString path) : InputReader(path)
{
    fInputFile = path;
    if (!AccessTree("data")) { throw "##### Exception in InputReaderSimple constructor!"; }
}
//------------------------------------------------------------------
/// Default destructor.
InputReaderSimple::~InputReaderSimple()
{
    if (fFile->IsOpen()) fFile->Close();
}
//------------------------------------------------------------------
bool InputReaderSimple::AccessTree(TString name)
{

    fTree = (TTree*)fFile->Get(name);

    if (fTree == NULL)
    {
        cout << "##### Error in InputReaderSimple::AccessTree()!" << endl;
        cout << "Could not access the tree!" << endl;
        return false;
    }

    fPoint0 = new TVector3();
    fPoint1 = new TVector3();
    fPoint2 = new TVector3();
    fVersor1 = new TVector3();
    fVersor2 = new TVector3();
    fPositionScat = new TVector3();
    fPositionAbs  = new TVector3();
    fDirectionScat = new TVector3();

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
bool InputReaderSimple::LoadEvent(int i)
{
    int imax = fTree->GetEntries();

    if (i > imax)
    {
        cout << "##### Error in InputReader::LoadEvent()!" << endl;
        cout << "Requested event number larger than number of events in the tree!" << endl;
        return false;
    }

    fTree->GetEntry(i);

    fPositionScat->SetXYZ(fPoint1->X(), fPoint1->Y(), fPoint1->Z());
    fPositionAbs->SetXYZ(fPoint2->X(), fPoint2->Y(), fPoint2->Z());
    fDirectionScat->SetXYZ(fVersor2->X(), fVersor2->Y(), fVersor2->Z());

    fPrimaryEnergy = fEnergy0;
    fEnergyLoss = fEnergy1;
    fEnergyScattered = fEnergy2;

    if (fSmear)
    {
        fPositionScat->SetXYZ(SmearBox(fPositionScat->X(), fScatResolutionX),
                              SmearGaus(fPositionScat->Y(), fScatResolutionY),
                              SmearBox(fPositionScat->Z(), fScatResolutionZ));
        fPositionAbs->SetXYZ(SmearGaus(fPositionAbs->X(), fAbsResolutionX),
                             SmearBox(fPositionAbs->Y(), fAbsResolutionY),
                             SmearBox(fPositionAbs->Z(), fAbsResolutionZ));
        fEnergyLoss = SmearGaus(fEnergyLoss, GetSigmaEScat(fEnergyLoss));
        fEnergyScattered = SmearGaus(fEnergyScattered, GetSigmaEAbs(fEnergyScattered));

        
        ///Makes sure that events are not out of the detector volumes.
         if (fPositionScat->Y() < -fScatHeight/2) fPositionScat->SetY(-fScatHeight/2);
         if (fPositionScat->Y() > fScatHeight/2) fPositionScat->SetY(fScatHeight/2);
         if (fPositionScat->Z() < -fScatWidth/2) fPositionScat->SetZ(-fScatWidth/2);
         if (fPositionScat->Z() > fScatWidth/2) fPositionScat->SetZ(fScatWidth/2);
         
         if (fPositionAbs->Y() < -fAbsHeight/2) fPositionAbs->SetY(-fAbsHeight/2);
         if (fPositionAbs->Y() > fAbsHeight/2) fPositionAbs->SetY(fAbsHeight/2);
         if (fPositionAbs->Z() < -fAbsWidth/2) fPositionAbs->SetZ(-fAbsWidth/2);
         if (fPositionAbs->Z() > fAbsWidth/2) fPositionAbs->SetZ(fAbsWidth/2);
    }
    return true;
}

//------------------------------------------------------------------
void InputReaderSimple::SelectEvents() {  
  fSelectedEvents.clear(); 
  for(int i = 0; i < fTree->GetEntries(); i++) {
  	fTree->GetEntry(i);
 	if(SelectSingleEvent()) fSelectedEvents.push_back(i);
  }
  cout << "Number of selected events: " << fSelectedEvents.size() << endl;
}
//------------------------------------------------------------------
bool InputReaderSimple::SelectSingleEvent() {
    if (fEnergy1 > 0.000015) return true;
    else return false;
}
//------------------------------------------------------------------
void InputReaderSimple::PrintEvent() {
    cout << "position in the scatterer: " << fPositionScat->X() << " " << fPositionScat->Y() << " " << fPositionScat->Z() << endl;
    cout << "position in the absorber: " << fPositionAbs->X() << " " << fPositionAbs->Y() << " " << fPositionAbs->Z() << endl;
    cout << "energy in the scatterer: " << fEnergyLoss << endl;
    cout << "energy in the absorber: " << fEnergyScattered << endl;
}
//------------------------------------------------------------------
/// Reads scatterer's and absorber's dimensions from Simple Simulation.
/// Works only if the files are in certain directories!
void InputReaderSimple::ReadGeometry(void) {
    string inputFile(fInputFile.Data());
    string geometryPath = inputFile;
    geometryPath.replace(0, 21, "results/CCSimulation_geometry_");
    geometryPath = geometryPath.substr(0, geometryPath.find("_scat")) + ".txt";
    ifstream geometry(geometryPath);
    
    if (geometry.is_open()) {
        string line;
        bool scattererSection = true;
        while (getline(geometry, line)) {
            if (line.rfind("fDimZ", 0) == 0) {
                double value = stod(line.substr(8));
                if (scattererSection) {
                    fScatWidth = value;
                } else {
                    fAbsWidth = value;
                }
            } else if (line.rfind("fDimY", 0) == 0) {
                double value = stod(line.substr(8));
                if (scattererSection) {
                    fScatHeight = value;
                    scattererSection = false;
                } else {
                    fAbsHeight = value;
                }
            }
        }
        geometry.close();
    }
}
//------------------------------------------------------------------
void InputReaderSimple::Clear(void)
{
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
void InputReaderSimple::SetSmearing(bool smear, Double_t posScatX, Double_t posAbsX, Double_t posScatY, Double_t posAbsY, Double_t posScatZ, Double_t posAbsZ)
{
    /// To smear camera performance(simple simulation), this file is called for
    /// energy resolution. It shows a function of energy deposited in a 10-cm-long
    /// LuAG(Ce) fiber
    ///  with a square cross-section from Geant4 simulation
    fSmear = smear;
    TString path = gSystem->Getenv("CC6DIR");
    TString name = path + "/share/ComptonCamera6/mlem_reco/EnergyResolutionExample.root";
    TFile* file = new TFile(name, "READ");
    fHisto = (TH1D*)file->Get("Scintillator_0ResolutionCombiEnergy");
    TF1* func = new TF1("fit1", "[0]+[1]/sqrt(x)+[2]/x^(3/2)", 0, 4.5);

    fHisto->Fit("fit1", "r");

    fScatResolutionX = posScatX;
    fScatResolutionY = posScatY;
    fScatResolutionZ = posScatZ;
    fAbsResolutionX = posAbsX;
    fAbsResolutionY = posAbsY;
    fAbsResolutionZ = posAbsZ;
}
//------------------------------------------------------------------
/// Gausian function to smear position resolution along y axis.
/// returns a double value from gausian distribution with the given mean and
/// sigma. \param val (double) - the given position value. \param sigma (double)
///- the given sigma value.
Double_t InputReaderSimple::SmearGaus(double val, double sigma)
{
    return gRandom->Gaus(val, sigma);
}
//------------------------------------
/// Returns a double value from Uniform function with respect to position
/// resolution value. \param x (double) - the given position value.
Double_t InputReaderSimple::SmearBox(double x, double resolution)
{
    return gRandom->Uniform(x - (resolution / 2), x + (resolution / 2));
}
//------------------------------------
/// Returns the sigma value from the fitting function of deposited energy plot.
///\param energy (double) - the given energy value.
Double_t InputReaderSimple::GetSigmaEScat(double energy)
{
    double sigma = fHisto->GetFunction("fit1")->Eval(energy) * energy * 7.7/7.2;
    return sigma;
}
//------------------------------------
Double_t InputReaderSimple::GetSigmaEAbs(double energy)
{
    double sigma = (-0.00488 * energy + 0.04688) * energy;
    return sigma;
}
