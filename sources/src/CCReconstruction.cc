#include "CCReconstruction.hh"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

ClassImp(CCReconstruction);

const double kMe  = 0.510999;	// MeV/c2

//------------------------------------------------------------------
CCReconstruction::CCReconstruction(TString inputName, TString name, Int_t iter, Bool_t verbose){
  
  SetInputName(inputName);
  SetName(name);
  SetIter(iter);
  fVerbose = verbose;
  TString tmp = name(20,22);
  fGenVersion = tmp.Atoi();
  
  fPoint0 = new TVector3;
  fPoint1 = new TVector3;
  fPoint2 = new TVector3;
  fVersor1 = new TVector3;
  fVersor2 = new TVector3;
  
  fFile = new TFile(fInputName,"READ");
  fTree = (TTree*)fFile->Get("data");
  
  fTree->SetBranchAddress("point0",&fPoint0);
  fTree->SetBranchAddress("point1",&fPoint1);
  fTree->SetBranchAddress("point2",&fPoint2);
  fTree->SetBranchAddress("versor1",&fVersor1);
  fTree->SetBranchAddress("versor2",&fVersor2);
  fTree->SetBranchAddress("energy0",&fEnergy0);
  fTree->SetBranchAddress("energy1",&fEnergy1);
  fTree->SetBranchAddress("energy2",&fEnergy2);

  fNev = fTree->GetEntries();
}
//------------------------------------------------------------------
CCReconstruction::~CCReconstruction(){
  if(fFile) fFile->Close();
}
//------------------------------------------------------------------
TVector3 CCReconstruction::ConnectPoints(TVector3 point1, TVector3 point2){
  TVector3 vec;
  vec.SetX(point2.X() - point1.X());
  vec.SetY(point2.Y() - point1.Y());
  vec.SetZ(point2.Z() - point1.Z());
  vec.SetMag(1.);
  return vec;	//versor
}
//------------------------------------------------------------------
Double_t CCReconstruction::CalculateTheta(Double_t e1, Double_t e2){
  Double_t costheta = 1. - kMe*(1./e2 - 1./e1);
  Double_t theta = TMath::ACos(costheta);	//rad
  return theta;
}
//------------------------------------------------------------------
void CCReconstruction::RebuildSetupTxt(void){
  
  cout << "\n----- Rebuilding setup from the txt file \n" << endl;
  
  ifstream input(Form("results/CCSimulation_geometry_gen%i.txt",fGenVersion), std::ios::in);
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
ComptonCone* CCReconstruction::ReconstructCone(Int_t i){
  
  Double_t epsilon = 1.E-8;
  
  Clear();
  fTree->GetEntry(i);
  
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
    
  if(i%1000==0) cout << i << " events out of " << fNev << " processed\n" << endl;
  
  return cone;
}
//------------------------------------------------------------------
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
Bool_t CCReconstruction::SaveHistogram(TH1F *h){
  TString name = "results/" + fName + ".root";
  TFile *file = new TFile(name,"UPDATE");
  h->Write();
  file->Close();
  if(fVerbose) cout << "\nHistogram " << h->GetName() << 
                       " saved in the file " << name << endl;
  return kTRUE;
}
//------------------------------------------------------------------
Bool_t CCReconstruction::SaveHistogram(TH2F *h){
  TString name = "results/" + fName + ".root";
  TFile *file = new TFile(name,"UPDATE");
  h->Write();
  file->Close();
  if(fVerbose) cout << "\nHistogram " << h->GetName() << 
                       " saved in the file " << name << endl;
  return kTRUE;
}
//------------------------------------------------------------------
void CCReconstruction::Print(void){
  cout << "\nCCReconstruction::Print() for the object: " << fName << endl;
  cout << "\tData reconstruction form the file: " << fInputName << endl;
  cout << "\tNumber of iterationf for MLEM: " << fIter << endl;
}
//------------------------------------------------------------------
