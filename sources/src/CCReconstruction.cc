#include "CCReconstruction.hh"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

ClassImp(CCReconstruction);

const double kMe  = 0.510999;	// MeV/c2

//------------------------------------------------------------------
///Standard constructor. Opens given ROOT file containing simulation data
///and accesses tree with events.
///\param inputName (TString) - path to the ROOT file with the simulation data
///\param name (TString) - name of the object
///\param iter (Int_t) - number of iterations for MLEM (not used for now)
///\param verbose (Bool_t) - verbose level for print-outs on the screen.
CCReconstruction::CCReconstruction(TString inputName, TString name, Bool_t verbose){
  
  SetInputName(inputName);
  SetName(name);
  fVerbose = verbose;
  TString tmp = name(20,22);
  fGenVersion = tmp.Atoi();
  
  if(!SetInputReader(inputName)){
    throw "##### Exception in CCReconstruction constructor!";
  }
  
  fPoint0 = new TVector3;
  fPoint1 = new TVector3;
  fPoint2 = new TVector3;
  fVersor1 = new TVector3;
  fVersor2 = new TVector3;
}
//------------------------------------------------------------------
///Default destructor.
CCReconstruction::~CCReconstruction(){
  delete fReader;
}
//------------------------------------------------------------------
bool CCReconstruction::SetInputReader(TString inputName){
  
  TFile *file = new TFile(inputName,"READ");
  if(!file->IsOpen()){
   cout << "##### Error in CCReconstruction::SetInputReader!" << endl;
   cout << "Could not open requested file" << endl;
   return false;
  }
  
  if(file->Get("data")){
   file->Close();
   fReader = new InputReaderSimple(inputName);
  }
  else if(file->Get("G4SimulationData_Reconstruction")){
    file->Close();
    fReader = new InputReaderGeant(inputName);
  }
  else{
    cout << "##### Error in CCReconstruction::SetInputReader()!" << endl;
    cout << "Unknown data format" << endl;
    return false;
  }
  
  return true;
}
//------------------------------------------------------------------
///Calculates coordinates of versor connecting two given points.
///\param point1 (TVector3) - coordinates of first point
///\param point2 (TVector3) - coordinates of second point
///
///In this reconstruction algorithm this function is used to find ComptonCone
///axis. Then point1 is interaction point in the scatterer and point2 is 
///interaction point in the absorber.
TVector3 CCReconstruction::ConnectPoints(TVector3 point1, TVector3 point2){
  TVector3 vec;
  vec.SetX(point2.X() - point1.X());
  vec.SetY(point2.Y() - point1.Y());
  vec.SetZ(point2.Z() - point1.Z());
  vec.SetMag(1.);
  return vec;	//versor
}
//------------------------------------------------------------------
///Calculates scattering angle theta. Angle is returned in radians.
///\param e1 (Double_t) - initail gamma energy
///\param e2 (Double_t) - gamma energy after Compton scattering.
Double_t CCReconstruction::CalculateTheta(Double_t e1, Double_t e2){
  Double_t costheta = 1. - kMe*(1./e2 - 1./e1);
  Double_t theta = TMath::ACos(costheta);	//rad
  return theta;
}
//------------------------------------------------------------------
///Opens text file containing details of simulated Compton Camera setup.
///Based on the information form the file rebuilds the setup for image ,
///reconstruction i.e. stes parameters and dimensions of fScatterer and 
///fAbsorber.
void CCReconstruction::RebuildSetupTxt(void){
  
  cout << "\n----- Rebuilding setup from the txt file \n" << endl;
  
  ifstream input(Form("../sources/results/CCSimulation_geometry_gen%i.txt",fGenVersion), std::ios::in);
  if(!(input.is_open())){
    cout << "##### Could not open CCSimulation_geometry.txt file! " << endl;
    cout << "##### Please check!" << endl;
    return;
  }
  
  Double_t scatPar[4];
  Double_t absPar[4];
  Double_t scatDim[2];
  Double_t absDim[2];
  TString scatName, absName;
  TString dummy;
  char line[100];
  
  input.getline(line,100);
  input.getline(line,100);
   
  if(fVerbose) cout << "\n----- Rebuilding the scatterer" << endl;
  input >> dummy >> dummy >> scatName;
  input >> dummy >> dummy >> scatPar[0];
  input >> dummy >> dummy >> scatPar[1];
  input >> dummy >> dummy >> scatPar[2];
  input >> dummy >> dummy >> scatPar[3];
  input >> dummy >> dummy >> scatDim[0];
  input >> dummy >> dummy >> scatDim[1];
  
  fScatterer.SetPlane(scatPar[0],scatPar[1],scatPar[2],scatPar[3]);
  fScatterer.SetDimensions(scatDim[0],scatDim[1]);
  fScatterer.SetName(scatName);
  if(fVerbose) fScatterer.Print();
  
  if(fVerbose) cout << "\n----- Rebuilding the scatterer" << endl;
  input >> dummy >> dummy >> absName;
  input >> dummy >> dummy >> absPar[0];
  input >> dummy >> dummy >> absPar[1];
  input >> dummy >> dummy >> absPar[2];
  input >> dummy >> dummy >> absPar[3];
  input >> dummy >> dummy >> absDim[0];
  input >> dummy >> dummy >> absDim[1];
  
  fAbsorber.SetPlane(absPar[0],absPar[1],absPar[2],absPar[3]);
  fAbsorber.SetDimensions(absDim[0],absDim[1]);
  fAbsorber.SetName(absName);
  if(fVerbose) fAbsorber.Print();
 
  input.close();
}
//------------------------------------------------------------------
///Based on the simulated data accessed from the tree reconstructs
///single ComptonCone class object. 
///\param i (Int_t) - number of event (cone) to be reconstructed.
ComptonCone* CCReconstruction::ReconstructCone(Int_t i){
  
  Double_t epsilon = 1.E-8;
  
  Clear();
  fReader->LoadEvent(i);
  fPoint0  = fReader->GetPositionPrimary();
  fPoint1  = fReader->GetPositionScattering();
  fPoint2  = fReader->GetPositionAbsorption();
  fVersor1 = fReader->GetGammaDirPrimary();
  fVersor2 = fReader->GetGammaDirScattered();
  fEnergy0 = fReader->GetEnergyPrimary();
  fEnergy1 = fReader->GetEnergyLoss();
  fEnergy2 = fReader->GetEnergyScattered();
  
  //----- energy check
  Double_t en = fEnergy1 + fEnergy2;
  if(fabs(en-fEnergy0) > epsilon){
    cout << "##### Error in CCReconstruction::ReconstructCones!" << endl;
    cout << "##### Energies do not sum correctly! Please check!" << endl;
    cout << "Event " << i << ": E0 = " << fEnergy0 << "\t E1 = " << fEnergy1 
         << "\t E2 = " << fEnergy2 << endl;
    return NULL;
  }
  //----- end of the energy check
  
  Double_t theta = CalculateTheta(fEnergy0,fEnergy2);
  
  //----- angle check
  Double_t ang = fVersor2->Angle(*fVersor1);
  if((TMath::Abs(ang-theta) > epsilon)){
    cout << "##### Error in CCReconstruction::ReconstructCones!" << endl;
    cout << "##### Theta angle is incorrect! Please check!" << endl;
    cout << "Event " << i << ": ang = " << ang*TMath::RadToDeg() 
         << " deg \t theta = " << theta*TMath::RadToDeg() << " deg \t delta = " 
         << ang*TMath::RadToDeg() - theta*TMath::RadToDeg() << " deg" << endl;
    return NULL;
  }
  //----- end of the angle check
  
  TVector3 axis = ConnectPoints(*fPoint1,*fPoint2);
  
  //----- axis check
  if(TMath::Abs(axis.X()-fVersor2->X()) > epsilon &&
     TMath::Abs(axis.Y()-fVersor2->Y()) > epsilon &&
     TMath::Abs(axis.Z()-fVersor2->Z()) > epsilon){
    cout << "##### Error in CCReconstruction::ReconstructCones!" << endl;
    cout << "##### Incorrect axis of the cone! Please check!" << endl;
    cout << "Event " << i << ": \nAxis \t";
    axis.Print();
    cout << "\nfVersor2 \t";
    fVersor2->Print();
    return NULL;
  }
  //----- end of the axis check
  
  ComptonCone *cone = new ComptonCone();
  cone->SetName(Form("cone_%i",i));
  cone->SetAngle(theta);
  cone->SetApex(*fPoint1);
  cone->SetAxis(axis);
  
  if(fVerbose){
    cout << "\n\tpoint1:\t";
    fPoint1->Print();
    cout << "\tpoint2:\t";
    fPoint2->Print(); 
    cout << "\taxis:\t";
    axis.Print();
    cout << "\tenergy0 = " << fEnergy0 << "\tenergy1 = " << fEnergy1 << "\tenergy2 = " << fEnergy2 << "\n";
    cout << "\ttheta = " << theta << " radians\t" << theta*TMath::RadToDeg()<<" deg\n\t";
    }
    
  //if(i%1000==0) cout << i << " events out of " << fNev << " processed\n" << endl;
  
  return cone;
}
//------------------------------------------------------------------
///Function for image reconstruction. Sets 2D image histogram and 
///1D histogram of filled pixels distribution. In the loop reconstructs
///requested numbers of ComptonCones and looks for their intersection with 
///the pixels of the image plane. Details of this simple back projection 
///algorithm can be found in presentation linked in the description of this class.
///\param iStart (Int_t) - number of first event for image reconstruction
///\param iStop (Int_t) - number of last event for image reconstruction
Bool_t CCReconstruction::ReconstructImage(Int_t iStart, Int_t iStop){
  
  cout << "\n\n----- Image reconstruction ----- \n\n" << endl;
  
  //image histogram
  Int_t nbinsz = (Int_t)fScatterer.GetDimZ();
  Int_t nbinsy = (Int_t)fScatterer.GetDimY();
  Double_t zlimit = fScatterer.GetDimZ()/2.;
  Double_t ylimit = fScatterer.GetDimY()/2.;
  fImage = new TH2F("image","image",nbinsz,-zlimit,zlimit,nbinsy,-ylimit,ylimit);
  fImage->GetXaxis()->SetTitle("z [mm]");
  fImage->GetYaxis()->SetTitle("y [mm]");
  fNpixels = new TH1F("npixels","npixels",5*nbinsz,-0.5,5*nbinsz-0.5);
  fNpixels->GetXaxis()->SetTitle("number of filled pixels");
  fNpixels->GetYaxis()->SetTitle("counts");

  Double_t pixelX = 0.;
  Double_t pixelY, pixelZ;
  Int_t bin;
  Double_t pixelSize = zlimit*2/nbinsz;	//mm
  Double_t resolution = 0.;
  Double_t npixels = 0;
  TVector3 pixelCenter;
  TVector3 interactionPoint;
  TVector3 coneAxis;
  TVector3 linkingVector;
  Double_t coneTheta, angle;
  
  for(Int_t i=iStart; i<iStop; i++){
  npixels = 0;
    if(fVerbose) cout << "----- processing event " << i << endl;
    ComptonCone *cone = ReconstructCone(i);
    interactionPoint = cone->GetApex();
    coneAxis = cone->GetAxis();
    coneTheta = cone->GetAngle();
    for(Int_t z=0; z<nbinsz; z++){
      for(Int_t y=0; y<nbinsy; y++){
        pixelY = fImage->GetYaxis()->GetBinCenter(y);
        pixelZ = fImage->GetXaxis()->GetBinCenter(z);
        pixelCenter.SetXYZ(pixelX,pixelY,pixelZ);
        linkingVector = ConnectPoints(pixelCenter,interactionPoint);
        angle = coneAxis.Angle(linkingVector);
        bin = fImage->GetBin(z,y);
	resolution = TMath::ATan(0.5*pixelSize*(sqrt(2))/
	             (fScatterer.GetD()/fScatterer.GetA()));
        if(TMath::Abs(coneTheta - angle) <= resolution){
	  fImage->Fill(pixelZ,pixelY);
	  npixels++;
	}
      }
    }
    fNpixels->Fill(npixels);
    delete cone;
  }
  
  SaveHistogram(fImage);
  SaveHistogram(fNpixels);
  
  return kTRUE;
}
//------------------------------------------------------------------ 
///Resests values of chosen private members of the class to their default values.
///Varaibles reset in this function: fPoint0, fPoint1, fPoint2, fVersor1, fVersor2,
///fEnergy0, fEnergy1, fEnergy2.
void CCReconstruction::Clear(void){
  fPoint0->SetXYZ(-1000,-1000,-1000);
  fPoint1->SetXYZ(-1000,-1000,-1000);
  fPoint2->SetXYZ(-1000,-1000,-1000);
  fVersor1->SetXYZ(-1000,-1000,-1000);
  fVersor2->SetXYZ(-1000,-1000,-1000);
  fEnergy0 = -1000;
  fEnergy1 = -1000;
  fEnergy2 = -1000;
}
//------------------------------------------------------------------
///Saves 1D histogram in the ROOT file. Name of the file is based on the name 
///of the CCReconstruction object.
Bool_t CCReconstruction::SaveHistogram(TH1F *h){
  TString name = "../sources/results/" + fName + ".root";
  TFile *file = new TFile(name,"UPDATE");
  h->Write();
  file->Close();
  if(fVerbose) cout << "\nHistogram " << h->GetName() << 
                       " saved in the file " << name << endl;
  return kTRUE;
}
//------------------------------------------------------------------
///Saves 2D histogram in the ROOT file. Name of the file is based on the name 
///of the CCReconstruction object.
Bool_t CCReconstruction::SaveHistogram(TH2F *h){
  TString name = "../sources/results/" + fName + ".root";
  TFile *file = new TFile(name,"UPDATE");
  h->Write();
  file->Close();
  if(fVerbose) cout << "\nHistogram " << h->GetName() << 
                       " saved in the file " << name << endl;
  return kTRUE;
}
//------------------------------------------------------------------
///Prints details of the CCReconstruction class object.
void CCReconstruction::Print(void){
  cout << "\nCCReconstruction::Print() for the object: " << fName << endl;
  cout << "\tData reconstruction form the file: " << fInputName << endl;
}
//------------------------------------------------------------------
