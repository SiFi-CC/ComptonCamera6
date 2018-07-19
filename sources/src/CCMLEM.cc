#include "CCMLEM.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "CCReconstruction.hh"
#include "IsectionPoint.hh"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TMath.h"
#include "SMElement.hh"
#include "TRandom.h"
using namespace std;

ClassImp(CCMLEM);



//--------------------

CCMLEM::CCMLEM(){
  
  
  
  fPoint0 = new TVector3;
  fPoint1 = new TVector3;
  fPoint2 = new TVector3;
  fVersor1 = new TVector3;
  fVersor2 = new TVector3;
  
  
  
  Config();
  TString name = "../sources/results/" + fName + ".root";
  TFile *file = new TFile(name,"RECREATE");
  //ofstream outfile ("../sources/results/" + fOutputName + ".root");
  Freshoutput(h);
  
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

  
  
  fArray = new TClonesArray("IsectionPoint",1000);
  fNIpoints = 0;
  fSM = new TClonesArray("SMElement",1000000);
  fpoints = 0;
  cout<<"CCMLEM: I will work in fVerbose="<<fVerbose<<" mode"<<endl;

}
//----------------------------------------
CCMLEM::~CCMLEM(){
  if(fFile) fFile->Close();
  if(file) file->Close();
}
//------------------------------------------
Bool_t CCMLEM::Reconstruct(Int_t iStart,Int_t iStop){

  CCReconstruction *reco = new CCReconstruction(fInputName, fName, fVerbose);

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
  Double_t pixelX = 0.;
  Double_t pixelY, pixelZ;
 

  Int_t j,k;
  
  TVector3 interactionPoint;
  TVector3 coneAxis;
  Double_t coneTheta;
  TVector3 uvec;
 
  Double_t F = fXofRecoPlane;
  Double_t z,y;
  Double_t A = fDimY/2.;  
  Double_t B = fDimZ/2.;
  
  
  fNIpoints = 0;
  fpoints = 0;
  IsectionPoint* tmppoint1;
  IsectionPoint* tmppoint2;

  const Double_t maxdist = sqrt(pow(fPixelSizeY,2)+pow(fPixelSizeZ,2));
  TStopwatch t;
  t.Start(); 
  
  for(Int_t i=iStart; i<iStop; i++){
    fNIpoints = 0;
     
    if(fVerbose)  cout<<"CCMLEM::Reconstruct(...) event "<< i<<endl<<endl;;
    ComptonCone *cone = reco->ReconstructCone(i);
    interactionPoint = cone->GetApex();
    coneAxis = cone->GetAxis();
    uvec = coneAxis.Unit();
    coneTheta = cone->GetAngle();
    Double_t K = cos(coneTheta);
    Double_t a, b, c1, c2, c3, c4, c, z1, z2;
    Double_t d, e, f1, f2, f3, f4, f, y1, y2;  
    
    y = -A;
    if(fVerbose) cout<<"Loop over horizontal lines..."<<endl;
    for (j=0; j<=fNbinsY; j++){
      
      a = 2*(pow(-uvec.Z(),2) - pow(K,2));
      
      b = 2*((-uvec.Z())*((interactionPoint.Z())*(-uvec.Z()) - y*((-uvec.Y())) + (interactionPoint.Y())*
			  (-uvec.Y()) + (interactionPoint.X())*(-uvec.X()) + F - (-uvec.X())*F) - (interactionPoint.Z())*pow(K,2));
      
      c1 = pow(interactionPoint.Y(),2) + pow(interactionPoint.X(),2) + ((-2)*pow(interactionPoint.Y(),2)*pow(-uvec.Z(),2)) +  
	((-2)*pow(interactionPoint.X(),2)*pow(-uvec.Z(),2)) + ((-2)*pow(interactionPoint.Y(),2)*pow(-uvec.Y(),2)) + 
	((-4)*(interactionPoint.Y())*(interactionPoint.X())*(-uvec.Y())*(-uvec.X())) + ((-2)*pow(interactionPoint.X(),2)*
											pow(-uvec.X(),2));
      
      c2 = (-2)*(interactionPoint.X())*F + 4*(interactionPoint.X())*F*
	pow(-uvec.Z(),2) - 4*(interactionPoint.Y())*(-uvec.Y())*F - 4*(interactionPoint.X())*
	(-uvec.X())*F + 4*(interactionPoint.Y())*(-uvec.Y())*(-uvec.X())*F + 
	4*(interactionPoint.X())*pow(-uvec.X(),2)*F - pow(F,2) - 2*pow(-uvec.Z(),2)*
	pow(F,2) + 4*(-uvec.X())*pow(F,2) - 2*pow(-uvec.X(),2)*pow(F,2);
      
      c3 = pow(y,2)*(1 - 2*pow(-uvec.Z(),2) - 2*pow(-uvec.Y(),2)) + y*
	((interactionPoint.Y())*((-2) + 4*pow(-uvec.Z(),2) + 4*pow(-uvec.Y(),2)) + 4*(-uvec.Y())*((interactionPoint.X())*(-uvec.X()) + 
												  F - (-uvec.X())*F));
      
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
    for (k=0; k<=fNbinsZ; k++){
      
      d = 2*(pow(-uvec.Y(),2) - pow(K,2));
      
      e = 2*((-uvec.Y())*((interactionPoint.Z())*(-uvec.Z()) - z*((-uvec.Z())) + (interactionPoint.Y())*
			  (-uvec.Y()) + (interactionPoint.X())*(-uvec.X()) + F - (-uvec.X())*F) - (interactionPoint.Y())*pow(K,2));
      
      f1 = pow(interactionPoint.Z(),2) + pow(interactionPoint.X(),2) + ((-2)*pow(interactionPoint.Z(),2)*pow(-uvec.Z(),2)) +  
	((-2)*pow(interactionPoint.Z(),2)*pow(-uvec.Y(),2)) + ((-2)*pow(interactionPoint.X(),2)*pow(-uvec.Y(),2)) + 
	((-4)*(interactionPoint.Z())*(interactionPoint.X())*(-uvec.Z())*(-uvec.X())) + ((-2)*pow(interactionPoint.X(),2)*
											pow(-uvec.X(),2));
      
      f2 = (-2)*(interactionPoint.X())*F + 4*(interactionPoint.X())*F*
	pow(-uvec.Y(),2) - 4*(interactionPoint.Z())*(-uvec.Z())*F - 4*(interactionPoint.X())*
	(-uvec.X())*F + 4*(interactionPoint.Z())*(-uvec.Z())*(-uvec.X())*F + 
	4*(interactionPoint.X())*pow(-uvec.X(),2)*F - pow(F,2) - 2*pow(-uvec.Y(),2)*
	pow(F,2) + 4*(-uvec.X())*pow(F,2) - 2*pow(-uvec.X(),2)*pow(F,2);
      
      f3 = pow(z,2)*(1 - 2*pow(-uvec.Z(),2) - 2*pow(-uvec.Y(),2)) + z*
	((interactionPoint.Z())*((-2) + 4*pow(-uvec.Z(),2) + 4*pow(-uvec.Y(),2)) + 4*(-uvec.Z())*((interactionPoint.X())*(-uvec.X()) + 
												  F - (-uvec.X())*F));
      
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
      if(fVerbose) cout<<"fpoints = "<<fpoints<<endl;
      temp = (SMElement*)fSM->ConstructedAt(fpoints++);
      temp->SetEvBinDist(i, binno1, dist);
      if(fVerbose) temp->Print();

    }
    //if(fVerbose) cout<<"end of loop"<<endl;
   
    delete cone;
   
    if(fVerbose) cout<<"----------------------------------------------------------------"<<endl;
     
  }// end of loop over events
  
  fArray->Clear("C");
  file->cd();
  fImage[0]->Write();
  //SaveToFile(fImage[0]);
  
  TH1D* hProZ[100];
  TH1D* hProY[100];

  for(int iter=1; iter<fIter+1; iter++){
    Iterate(iStop, iter);
  }
 
  t.Stop(); 
  t.Print();

  delete reco;
   
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
  //Int_t nstop;
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
  
  file->cd();
  hthisiter->Write();
 //SaveToFile(hthisiter);
  
  return kTRUE;
}
//------------------------------------
Bool_t CCMLEM::Config(void){
  
  
  string nextline;
  
  TString config_name = "../sources/results/config.txt";
  ifstream config(config_name);
  if(!(config.is_open())){

    cout << "##### Could not open configuration file!" << endl;
    return kFALSE;
  
    TString ReadLine(std::istream& config);  //Read a line from stream upto newline skipping any whitespace.
    getline(config, nextline);
 
    config >> fInputName;
   
    while(fInputName){
      
      config >> fXofRecoPlane >> fYofRecoPlane >> fZofRecoPlane;
      if(fXofRecoPlane<0 || fYofRecoPlane<0 || fZofRecoPlane<0){
        cout << "##### Please check the center of image plane!" << endl;
        return kFALSE;
      }
  
      config >> fDimZ >> fDimY;
      if(fDimZ<1 || fDimY<1){
        cout << "##### impossible size dimensions!" << endl;
        return kFALSE;
      }
  
      config >> fNbinsZ >> fNbinsY;
      if(fNbinsZ<2 || fNbinsY<2){
        cout << "##### impossible no. of bins!" << endl;
        return kFALSE;
      }
  
      config >> fSigmaX >> fSigmaY >> fSigmaZ;
      if(fSigmaX<-0.1 || fSigmaY<-0.1 || fSigmaZ<-0.1){
        cout << "##### Unexpected values of sigma" << endl;
        return kFALSE;
      }
  
      config >> fSigmaE;
      if(fSigmaE<0){
        cout << "##### Unexpected value of Energy Sigma!" << endl;
        return kFALSE;
      }
  
      config >> fIter;
      if(fIter<0){
        cout << "##### impossible no. of iterations!" << endl;
        return kFALSE;
      }

     config >> fStart >> fStop;
     if(fStart<0 || fStop<1){
        cout << "##### impossible no. of events!" << endl;
        return kFALSE;
      }
      
     config >> fVerbose;
     if(fVerbose<-1){
        cout << "##### Please check verbose!" << endl;
        return kFALSE;
      }
    
    }
  }
  config.close();

  return kTRUE;
}
//------------------------------------
Bool_t CCMLEM::Drawhisto(void){
  
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
  file->cd();
  can->Write();
  canz->Write();
  cany->Write();
  /*SaveToFile(can);
  SaveToFile(canz);
  SaveToFile(cany);*/
  return kTRUE;
}
//------------------------------------
Double_t Smear(double val, double sigma){ 
 
  return gRandom->Gaus(val, sigma);
}
//-------------------------------------------
void CCMLEM::Print(void){
 cout << "\nCCMLEM::Print()" <<endl;
 cout << "\tName of input file: \t"<<"../sources/results/" + fInputName <<endl;
 cout << "\tCenter of reco plane: "<<fXofRecoPlane << ", " << fYofRecoPlane <<", " << fZofRecoPlane << endl;
 cout << "\tSize of image: \t"<<fDimZ<< ", "<<fDimY<<endl;
 cout << "\tNo. of bins: \t"<<fNbinsZ<< ", "<<fNbinsY<<endl;
 cout << "\tPos resolution: \t"<<fSigmaX<< ", "<<fSigmaY<<", "<<fSigmaZ<<endl;
 cout << "\tEnergy resolution: \t"<<fSigmaE<<endl;
 cout << "\tNo. of MLEM iterations: \t"<<fIter<<endl;
 cout << "\tNo. of first and last event: \t"<<fStart<<", "<<fStop<<endl;
 cout << "\tVerbose flag: \t"<<fVerbose<<endl;
}
//-------------------------------------------
Bool_t CCMLEM::Freshoutput(TObject *ob){
  TString name = "../sources/results/" + fName + ".root";
  TFile *file = new TFile(name,"RECREATE");
  ob->Write();
  file->Close();
  cout << ob->ClassName()<<" " << ob->GetName() << 
                       " saved in the file " << name << endl;
  return kTRUE;
} 
//------------------------------------
/*Bool_t CCMLEM::SaveToFile(TGraph *h){
  TString name = "../sources/results/" + fName + ".root";
  TFile *file = new TFile(name,"UPDATE");
  h->Write();
  file->Close();
  if(fVerbose) cout << "\nObject " << h->GetName() << 
                       " saved in the file " << name << endl;
  return kTRUE;
}
//--------------------------------------
Bool_t CCMLEM::SaveToFile(TObject *ob){
  TString name = "../sources/results/" + fName + ".root";
  TFile *file = new TFile(name,"UPDATE");
  ob->Write();
  file->Close();
  cout << ob->ClassName()<<" " << ob->GetName() << 
                       " saved in the file " << name << endl;
  return kTRUE;
} 
*/

