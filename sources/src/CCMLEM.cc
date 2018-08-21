#include "CCMLEM.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "IsectionPoint.hh"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "SMElement.hh"
#include "TRandom.h"

using namespace std;

ClassImp(CCMLEM);

//--------------------
///Standard constructor (recommended).
///\param path (TString) - full path to the configuration file.
CCMLEM::CCMLEM(TString path){
  
  Clear();
  
  bool stat_config = ReadConfig(path);
  bool stat_reader = SetInputReader();
  if(!stat_config || !stat_reader){
    throw "##### Exception in CCMLEM constructor!";
  }
  
  Size_t len = strlen(fInputName);
  TString fname = fInputName;
  fname.Insert(len-5,"_MLEM");
  TString outputName = "../sources/results/" + fname;
  outputName.ReplaceAll("CCSimulation","CCReconstruction");
  TString option = (fFreshOutput ? "RECREATE" : "UPDATE");
  fOutputFile = new TFile(outputName,option);
  
  fArray = new TClonesArray("IsectionPoint",1000);
  fSM    = new TClonesArray("SMElement",1000000);
  fNIpoints = 0;
  fPoints   = 0;
}
//--------------------
///Default constructor.
CCMLEM::CCMLEM(){
  cout << "##### Warning in CCMLEM constructor!" << endl;
  cout << "You are using default constructor!" << endl;
  Clear();
}
//----------------------------------------
///Default destructor.
CCMLEM::~CCMLEM(){
  if(fOutputFile) fOutputFile->Close();
}
//------------------------------------------
///Based on the given input simulation file creates suitable
///InputReader object - either for simple input or Geant4 input. 
Bool_t CCMLEM::SetInputReader(void){
  
  TString fullName = "../sources/results/"+fInputName;
  
  TFile *file = new TFile(fullName,"READ");
  if(!file->IsOpen()){
   cout << "##### Error in CCMLEM::SetInputReader!" << endl;
   cout << "Could not open requested file" << endl;
   return false;
  }
  
  if(file->Get("data")){
   file->Close();
   fReader = new InputReaderSimple(fullName);
  }
  else if(file->Get("G4SimulationData_Reconstruction")){
    file->Close();
    fReader = new InputReaderGeant(fullName);
  }
  else{
    cout << "##### Error in CCMLEM::SetInputReader()!" << endl;
    cout << "Unknown data format" << endl;
    return false;
  }
  
  if(fVerbose){
    fReader->Print();
  }
  
  return true;
}
//------------------------------------------
Bool_t CCMLEM::Reconstruct(void){
  
  //image histogram
  fImage[0] = new TH2F("image","image",fNbinsZ,-fDimZ/2.,fDimZ/2.,fNbinsY,-fDimY/2.,fDimY/2.);
  fImage[0]->GetXaxis()->SetTitle("z [mm]");
  fImage[0]->GetYaxis()->SetTitle("y [mm]");

 /*fGraph = new TGraph();
  fGraph->SetName("g");
  fGraph->SetTitle("non-pixelized reco image");
  fGraph->SetMarkerStyle(7);
  fGraph->SetMarkerColor(kRed);
  TH1F* htmp = new TH1F("g","non-pixelized reco image",fNbinsZ,-fDimZ/2.,fDimZ/2.);
  htmp->GetYaxis()->SetRangeUser(-fDimY/2.,fDimY/2.);
  htmp->GetXaxis()->SetTitle("z [mm]");
  htmp->GetYaxis()->SetTitle("y [mm]");
  fGraph->SetHistogram(htmp);
 */
 
  fPixelSizeZ = fDimZ/fNbinsZ;
  fPixelSizeY = fDimY/fNbinsY;
  Double_t A = fDimY/2.;  
  Double_t B = fDimZ/2.;
  Double_t F = fXofRecoPlane;
  
  Int_t    j,k;
  Double_t z,y;
  
  TVector3 interactionPoint;
  TVector3 coneAxis;
  Double_t coneTheta;
  Double_t energy1, energy2;
  TVector3 *point1, *point2;
  
  fNIpoints = 0;
  fPoints   = 0;
  IsectionPoint* tmppoint1;
  IsectionPoint* tmppoint2;

  const Double_t maxdist = sqrt(pow(fPixelSizeY,2)+pow(fPixelSizeZ,2));
  
  bool status;
  int  counter = 0;
  
  TStopwatch t;
  t.Start(); 
  
  //This loop will always analyze fStop-fStart events starting with 
  //event number fStart. If some of the events are not valid they will 
  //be skipped, but still fStop-fStart events will be analyzed i.e.
  //last analyzed event will have number fStop+n, where n is number of
  //skipped events. If you want to change this - remove 'counter' variable.
  
  for(Int_t i=fStart; i<fStop; i++){
    
    fNIpoints = 0;
     
    if(fVerbose)  cout<<"CCMLEM::Reconstruct(...) event "<< i<<endl<<endl;
    
    status = fReader->LoadEvent(counter+i);
    
    if(status==false){
      counter++;
      continue;
    }
    
    energy1 = fReader->GetEnergyPrimary();
    energy2 = fReader->GetEnergyScattered();
    point1 =  fReader->GetPositionScattering();
    point2 =  fReader->GetPositionAbsorption();
    
    ComptonCone *cone = new ComptonCone(point1,point2,energy1,energy2);
    interactionPoint  = cone->GetApex();
    coneAxis          = cone->GetAxis();
    
    if(fSmear){
    interactionPoint.SetXYZ(Smear(interactionPoint.X(), fResolutionX),
		      Smear(interactionPoint.Y(), fResolutionY),Smear(interactionPoint.Z(), fResolutionZ));
    coneAxis.SetXYZ(Smear(coneAxis.X(), fResolutionX),
		      Smear(coneAxis.Y(), fResolutionY),Smear(coneAxis.Z(), fResolutionZ));
    }
   
    coneTheta = cone->GetAngle();
    Double_t K = cos(coneTheta);
    Double_t a, b, c1, c2, c3, c4, c, z1, z2;
    Double_t d, e, f1, f2, f3, f4, f, y1, y2;  
    
    y = -A;
    if(fVerbose) cout<<"Loop over horizontal lines..."<<endl;
    
    for(j=0; j<=fNbinsY; j++){
      
      a = 2*(pow(-coneAxis.Z(),2) - pow(K,2));
      
      b = 2*((-coneAxis.Z())*((interactionPoint.Z())*(-coneAxis.Z()) - y*((-coneAxis.Y())) + (interactionPoint.Y())*
	  (-coneAxis.Y()) + (interactionPoint.X())*(-coneAxis.X()) + F - (-coneAxis.X())*F) - (interactionPoint.Z())*pow(K,2));
      
      c1 = pow(interactionPoint.Y(),2) + pow(interactionPoint.X(),2) + ((-2)*pow(interactionPoint.Y(),2)*pow(-coneAxis.Z(),2)) +  
	   ((-2)*pow(interactionPoint.X(),2)*pow(-coneAxis.Z(),2)) + ((-2)*pow(interactionPoint.Y(),2)*pow(-coneAxis.Y(),2)) + 
	   ((-4)*(interactionPoint.Y())*(interactionPoint.X())*(-coneAxis.Y())*(-coneAxis.X())) + ((-2)*pow(interactionPoint.X(),2)*
	   pow(-coneAxis.X(),2));
      
      c2 = (-2)*(interactionPoint.X())*F + 4*(interactionPoint.X())*F*
	   pow(-coneAxis.Z(),2) - 4*(interactionPoint.Y())*(-coneAxis.Y())*F - 4*(interactionPoint.X())*
	   (-coneAxis.X())*F + 4*(interactionPoint.Y())*(-coneAxis.Y())*(-coneAxis.X())*F + 
	   4*(interactionPoint.X())*pow(-coneAxis.X(),2)*F - pow(F,2) - 2*pow(-coneAxis.Z(),2)*
	   pow(F,2) + 4*(-coneAxis.X())*pow(F,2) - 2*pow(-coneAxis.X(),2)*pow(F,2);
      
      c3 = pow(y,2)*(1 - 2*pow(-coneAxis.Z(),2) - 2*pow(-coneAxis.Y(),2)) + y*
	  ((interactionPoint.Y())*((-2) + 4*pow(-coneAxis.Z(),2) + 4*pow(-coneAxis.Y(),2)) + 4*(-coneAxis.Y())*
	  ((interactionPoint.X())*(-coneAxis.X()) +  F - (-coneAxis.X())*F));
      
      c4 = (pow(y,2) - 2*y*(interactionPoint.Y()) + pow(interactionPoint.Y(),2) + pow((interactionPoint.X() - F),2))*
	   ((-1) + 2*pow(K,2));
      
      
      c = (-2)*pow(K,2)*(c1 + c2 + c3 + c4); 
      
      if(c >= 0){
	z1 = (b - sqrt(c))/a;
        z2 = (b + sqrt(c))/a;
	if(fabs(z1)<1e-14) z1 = z1*1e6;
	if(fabs(z2)<1e-14) z2 = z2*1e6;
	AddIsectionPoint("hor", fXofRecoPlane, y, z1);
	AddIsectionPoint("hor", fXofRecoPlane, y, z2);
	
      }
      y=y+fPixelSizeY;
      
    } //end of loop over horizontal lines
    
    
    z = -B; 
    if(fVerbose) cout<<"Loop over vertical lines..."<<endl;
    for(k=0; k<=fNbinsZ; k++){
      
      d = 2*(pow(-coneAxis.Y(),2) - pow(K,2));
      
      e = 2*((-coneAxis.Y())*((interactionPoint.Z())*(-coneAxis.Z()) - z*((-coneAxis.Z())) + (interactionPoint.Y())*
	  (-coneAxis.Y()) + (interactionPoint.X())*(-coneAxis.X()) + F - (-coneAxis.X())*F) - (interactionPoint.Y())*pow(K,2));
      
      f1 = pow(interactionPoint.Z(),2) + pow(interactionPoint.X(),2) + ((-2)*pow(interactionPoint.Z(),2)*pow(-coneAxis.Z(),2)) +  
	   ((-2)*pow(interactionPoint.Z(),2)*pow(-coneAxis.Y(),2)) + ((-2)*pow(interactionPoint.X(),2)*pow(-coneAxis.Y(),2)) + 
	   ((-4)*(interactionPoint.Z())*(interactionPoint.X())*(-coneAxis.Z())*(-coneAxis.X())) + ((-2)*pow(interactionPoint.X(),2)*
	   pow(-coneAxis.X(),2));
      
      f2 = (-2)*(interactionPoint.X())*F + 4*(interactionPoint.X())*F*
	   pow(-coneAxis.Y(),2) - 4*(interactionPoint.Z())*(-coneAxis.Z())*F - 4*(interactionPoint.X())*
	   (-coneAxis.X())*F + 4*(interactionPoint.Z())*(-coneAxis.Z())*(-coneAxis.X())*F + 
	   4*(interactionPoint.X())*pow(-coneAxis.X(),2)*F - pow(F,2) - 2*pow(-coneAxis.Y(),2)*
	   pow(F,2) + 4*(-coneAxis.X())*pow(F,2) - 2*pow(-coneAxis.X(),2)*pow(F,2);
      
      f3 = pow(z,2)*(1 - 2*pow(-coneAxis.Z(),2) - 2*pow(-coneAxis.Y(),2)) + z*
	   ((interactionPoint.Z())*((-2) + 4*pow(-coneAxis.Z(),2) + 4*pow(-coneAxis.Y(),2)) + 4*(-coneAxis.Z())*
	   ((interactionPoint.X())*(-coneAxis.X()) + F - (-coneAxis.X())*F));
      
      f4 = (pow(z,2) - 2*z*(interactionPoint.Z()) + pow(interactionPoint.Z(),2) + pow((interactionPoint.X() - F),2))*
	   ((-1) + 2*pow(K,2));
      
      f = (-2)*pow(K,2)*(f1 + f2 + f3 + f4);
      
      if(f >= 0){
	y1 = (e - sqrt(f))/d;
	y2 = (e + sqrt(f))/d;
	if(fabs(y1)<1e-14) y1 = y1*1e6;
	if(fabs(y2)<1e-14) y2 = y2*1e6;
	AddIsectionPoint("ver", fXofRecoPlane, y1, z);
	AddIsectionPoint("ver", fXofRecoPlane, y2, z);
	
      }      
      z=z+fPixelSizeZ;
      
    } //end of loop over vertical lines
    
    Int_t index[fNIpoints]; 
    Int_t fA[fNIpoints];
    IsectionPoint *tempp;
    
    for(Int_t i=0; i<fNIpoints; i++){
       tempp = (IsectionPoint*)fArray->At(i);
       fA[i] = tempp->GetBin();
    }

    TMath::Sort(fNIpoints, fA, index, kFALSE);
    
    TVector3 *tmpvec1;
    TVector3 *tmpvec2;
    Double_t dist;
    Int_t binno1, binno2;
    SMElement* temp;
    
    for(int h=0; h<fNIpoints; h=h+2) {
     
      if(fVerbose)  cout<<" index["<<h<<"]="<<index[h]<<", index["<<h+1<<"]="<<index[h+1]<<endl;
      tmppoint1 = (IsectionPoint*)fArray->At(index[h]);
      tmppoint2 = (IsectionPoint*)fArray->At(index[h+1]);
      if(tmppoint1==NULL || tmppoint2==NULL){
	cout<<"Something went wrong"<<tmppoint1<<"\t"<<tmppoint2<<endl;
      }
      tmpvec1 = tmppoint1->GetPointCoordinates();
      binno1 = tmppoint1->GetBin();
      tmpvec2 = tmppoint2->GetPointCoordinates();
      binno2 = tmppoint2->GetBin();
      if(fVerbose)  cout<<" binno1="<<binno1<<", binno2="<<binno2<<endl<<endl;
      dist = ((*tmpvec1)-(*tmpvec2)).Mag();
      if(dist > maxdist){
	//cout<<"Event "<<h<<": distance exceeds pixel diagonal "<<dist/maxdist<<" times"<<endl;
	continue;
      }
      if(binno1!=binno2){
	if(fVerbose) cout<<binno1<<"!="<<binno2<<" ->Bin numbers are different when they should not!"<<endl;
	//h--;
      }
       
      fImage[0]->SetBinContent(binno1, fImage[0]->GetBinContent(binno1) + dist);
      if(fVerbose) cout<<"fPoints = "<<fPoints<<endl;
      temp = (SMElement*)fSM->ConstructedAt(fPoints++);
      temp->SetEvBinDist(i, binno1, dist);
      if(fVerbose) temp->Print();

    }
    //if(fVerbose) cout<<"end of loop"<<endl;
   
    delete cone;
   
    if(fVerbose) cout<<"----------------------------------------------------------------"<<endl;
     
  }// end of loop over events
  
  fArray->Clear("C");
  
  SaveToFile(fImage[0]);
  
  TH1D* hProZ[100];
  TH1D* hProY[100];

  for(int iter=1; iter<fIter+1; iter++){
    Iterate(fStop, iter);
  }
 
  t.Stop(); 
  t.Print();
   
  return kTRUE;
}

//------------------------------------
Int_t CCMLEM::AddIsectionPoint(TString dir, Double_t x, Double_t y, Double_t z){
  if(fVerbose) cout<<dir<<"\t"<<x<<"\t"<<y<<"\t"<<z<<endl;
  dir.ToLower();
  if(dir!="hor" && dir!="ver"){
    if(fVerbose) cout<<"Unknown direction of intrsecting line: "<<dir<<endl;
    return 0;
  }
 
  if(fabs(y)>fDimY/2. || fabs(z)>fDimZ/2.){
   // if(fVerbose) cout<<"point outside of image range..."<<endl;
    return 0;
  }
 
  Int_t i;
  IsectionPoint* tmppoint;
  //SMElement* temp;
  Int_t added = 0;
  Int_t pixelZ, pixelY;
  Double_t yplus, yminus;
  Double_t zplus, zminus;
  if(dir=="hor"){ // adding point from intersections with horizontal lines
    pixelZ = fImage[0]->GetXaxis()->FindBin(z);
    if(pixelZ>fNbinsZ) pixelZ=fNbinsZ; //inclusion of upper edges of histo
    yplus = y +0.5*fPixelSizeY;
    yminus = y -0.5*fPixelSizeY;
   
    if(fabs(yplus)<=fDimY/2){
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);
      pixelY = fImage[0]->GetYaxis()->FindBin(yplus);
      tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ,pixelY), x, y, z);
     //if(fVerbose) tmppoint->Print();
     //fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
     
    }
    if(fabs(yminus)<=fDimY/2){
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);
      pixelY = fImage[0]->GetYaxis()->FindBin(yminus);
      tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ,pixelY), x, y, z);
      //if(fVerbose) tmppoint->Print();
      //fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
      
    }
      
  }
  if(dir=="ver"){ // adding point from intersections with vertical lines
    pixelY = fImage[0]->GetYaxis()->FindBin(y);
    if(pixelY>fNbinsY) pixelY=fNbinsY; //inclusion of upper edges of histo
    zplus = z +0.5*fPixelSizeZ;
    zminus = z -0.5*fPixelSizeZ;
   
    if(fabs(zplus)<=fDimZ/2){
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);
      
      pixelZ = fImage[0]->GetXaxis()->FindBin(zplus);
      tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ,pixelY), x, y, z);
     //if(fVerbose) tmppoint->Print();
     //fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
      
    }
    if(fabs(zminus)<=fDimZ/2){
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);
      
      pixelZ = fImage[0]->GetXaxis()->FindBin(zminus);
      tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ,pixelY), x, y, z); 
      //if(fVerbose) tmppoint->Print();
      //fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
      
    }
    
  }
  //if(fVerbose)  cout<<added<<" points added..."<<endl<<endl;;
  return added;
  
}

//------------------------------------
Bool_t CCMLEM::Iterate(Int_t nstop, Int_t iter){
 
  int lastiter = iter-1;

  if(fImage[lastiter]==NULL){
    cout << "Error in CCMELM::Iterate(). Last iteration NULL" << endl;
    return kFALSE;
  }
  
  TH2F* hlastiter = (TH2F*)fImage[lastiter];
  fImage[lastiter+1]=(TH2F*)hlastiter->Clone();
  TH2F* hthisiter = fImage[lastiter+1];
  hthisiter->Reset();
  hthisiter->SetName(Form("%s_iter%i",fImage[0]->GetName(), lastiter+1));
  hthisiter->SetTitle(Form("%s_iter%i",fImage[0]->GetTitle(), lastiter+1));
  Int_t eventno;
  Int_t eventno_prev=0;
  Int_t entry;
  Int_t binno;
  Double_t dist, addvalue;
  SMElement* temp;
  Double_t denominator[nstop+1];
  for(int i=0; i<nstop+1; i++)
    denominator[i]=0;
  Int_t nSMentries = fSM->GetEntries();
  //cout<<"nSMentries = "<<nSMentries <<endl;
  for(entry=0; entry<nSMentries; entry++){
    temp = (SMElement*)fSM->At(entry);
    binno=temp->GetBin();
    eventno=temp->GetEvent();
    dist=temp->GetDist();
    denominator[eventno]= denominator[eventno]+dist*hlastiter->GetBinContent(binno);
  }
  
  for(entry=0; entry<nSMentries; entry++){
    addvalue = 0;
    temp = (SMElement*)fSM->At(entry);
    binno=temp->GetBin();
    eventno=temp->GetEvent();
    dist=temp->GetDist();
    addvalue=dist*hlastiter->GetBinContent(binno)/denominator[eventno];
    hthisiter->SetBinContent(binno,hthisiter->GetBinContent(binno)+addvalue);
  }
  
  
  SaveToFile(hthisiter);
  
  return kTRUE;
}
//------------------------------------
///Reads configuration file and sets values of private class
///members according to read information.
///\param path (TString) - full path to the configuration file.
Bool_t CCMLEM::ReadConfig(TString path){

  ifstream config(path);

  if(!(config.is_open())){
    cout << "##### Could not open configuration file!" << endl;
    return kFALSE;
  }
  
  cout << "\n\nIn CCMLEM::Config(). Reading config file..." << endl;
  
  TString comment;
  
  while(!config.eof()){
    comment.ReadLine(config);
    if(comment.Contains("Name of the input file")){
      config >> fInputName;
      if(!fInputName.Contains(".root")){
	cout << "##### Error in CCMLEM::Config()! Unknown file type!" << endl;
	return false;
      }
    }
    else if(comment.Contains("Center of reco plane")){
      config >> fXofRecoPlane >> fYofRecoPlane >> fZofRecoPlane;
    }
    else if(comment.Contains("Size of image")){
      config >> fDimZ >> fDimY;
      if(fDimZ<1 || fDimY<1){
	cout << "##### Error in CCMLEM::Config()! Image size incorrect!" << endl;
	return false;
      }
    }
    else if(comment.Contains("No. of bins")){
      config >> fNbinsZ >> fNbinsY;
      if(fNbinsZ<1 || fNbinsY<1){
	cout << "##### Error in CCMLEM::Config()! Number of bins incorrect!" << endl;
	return false;
      }
    }
    else if(comment.Contains("Smear")){
      config >> fSmear;
    }
    else if(comment.Contains("Position resolution")){
      config >> fResolutionX >> fResolutionY >> fResolutionZ;
    }
    else if(comment.Contains("Energy resolution")){
      config >> fSigmaE;
    }
    else if(comment.Contains("No. of MLEM iterations")){
      config >> fIter; 
      if(fIter<0){
	cout << "##### Error in CCMLEM::Config()! Number of iterations incorrect!" << endl;
	return false;
      }
    }
    else if(comment.Contains("Fresh output")){
      config >> fFreshOutput;
    }
    else if(comment.Contains("No. of first and last event")){
      config >> fStart >> fStop;
      if(fStart<0 || fStop<0 || fStop<fStart){
	cout << "##### Error in CCMLEM::Config()! Number of first or last event incorrect!" << endl;
	return false;
      }
    }
    else if(comment.Contains("Verbose flag")){
      config >> fVerbose;
    }
    else{
      cout << "##### Warning in CCMLEM::Config()! Unknown syntax!" << endl;
      cout << comment << endl;
    }
  }
  
  if(fVerbose) Print();
  
  config.close();

  return true;
}
//------------------------------------
Bool_t CCMLEM::DrawHisto(void){
  
  TH1D* hProZ[100];
  TH1D* hProY[100];
  TCanvas* can  = new TCanvas("MLEM2D","MLEM2D",1000,1000);
  TCanvas* canz = new TCanvas("MLEM1DZ","MLEM1DZ",1000,1000);
  TCanvas* cany = new TCanvas("MLEM1DY","MLEM1DY",1000,1000);
  can->Divide((int)sqrt(fIter)+1, (int)sqrt(fIter)+1);
  canz->Divide((int)sqrt(fIter)+1, (int)sqrt(fIter)+1);
  cany->Divide((int)sqrt(fIter)+1, (int)sqrt(fIter)+1);
  for(int iter=0; iter<fIter+1; iter++){
    can->cd(iter+1);
    gPad->SetLogz(1);
    fImage[iter]->Draw("colz");
    //fImage[iter]->SetMinimum(200);
    hProZ[iter]=fImage[iter]->ProjectionX();
    hProY[iter]=fImage[iter]->ProjectionY();
    canz->cd(iter+1);
    hProZ[iter]->Draw();
    cany->cd(iter+1);
    hProY[iter]->Draw();
  }
  
  SaveToFile(can);
  SaveToFile(canz);
  SaveToFile(cany);
  return kTRUE;
}
//------------------------------------
Double_t CCMLEM::Smear(double val, double sigma){ 
  return gRandom->Gaus(val,sigma);
}
//-------------------------------------------
///Prints details of the CCMLEM class obejct.
void CCMLEM::Print(void){
 cout << "\nCCMLEM::Print()" << endl;
 cout << setw(35) << "Name of input file: \t" << "../sources/results/" + fInputName << endl;
 cout << setw(35) << "Center of reco plane: \t" << fXofRecoPlane << ", " 
      << fYofRecoPlane << ", " << fZofRecoPlane << endl;
 cout << setw(35) << "Size of image: \t" << fDimZ << ", " << fDimY << endl;
 cout << setw(35) << "No. of bins: \t" << fNbinsZ << ", " << fNbinsY << endl;
 cout << setw(35) << "Smear level: \t" << fSmear << endl;
 cout << setw(35) << "Pos resolution: \t" << fResolutionX << ", " 
      << fResolutionY << ", " << fResolutionZ << endl;
 cout << setw(35) << "Energy resolution: \t" << fSigmaE << endl;
 cout << setw(35) << "No. of MLEM iterations: \t" << fIter << endl;
 cout << setw(35) << "FreshOutput level: \t" << fFreshOutput << endl;
 cout << setw(35) << "No. of first and last event: \t" << fStart << ", " << fStop << endl;
 cout << setw(35) << "Verbose level: \t" << fVerbose << endl << endl;
}
//--------------------------------------
///Sets default values of the private class members.
void CCMLEM::Clear(void){
  fInputName    = "dummy";
  fXofRecoPlane = -1000;
  fYofRecoPlane = -1000;
  fZofRecoPlane = -1000;
  fDimZ         = -1000;
  fDimY         = -1000;
  fNbinsZ	= -1000;
  fNbinsY       = -1000;
  fResolutionX  = -1000;
  fResolutionY  = -1000;
  fResolutionZ  = -1000;
  fSigmaE       = -1000;
  fIter         = -1000;
  fStart        = -1000;
  fStop         = -1000;
  fSmear        = kFALSE;
  fFreshOutput  = kFALSE;
  fVerbose      = kFALSE;
  fNIpoints     = -1000;
  fPoints       = -1000;
  fPixelSizeZ   = -1000;
  fPixelSizeY   = -1000;
  fReader       = NULL;
  fArray        = NULL;
  fSM           = NULL;
  fOutputFile   = NULL;
  fGraph        = NULL;
}
//--------------------------------------
///Saves object in the output file.
///\param ob (TObject*) - saved object.
Bool_t CCMLEM::SaveToFile(TObject *ob){
  fOutputFile->cd();
  ob->Write();
  cout << ob->ClassName() << " " << ob->GetName() 
       << " saved in the file " << fOutputFile->GetName() << endl;
  return kTRUE;
} 
//--------------------------------------