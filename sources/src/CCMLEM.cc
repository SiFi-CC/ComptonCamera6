#include "CCMLEM.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include "CCReconstruction.hh"
#include "IsectionPoint.hh"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TMath.h"
using namespace std;

ClassImp(CCMLEM);



//--------------------
///Standard constructor. Opens given ROOT file containing simulation data
///and accesses tree with events.
///\param inputName (TString) - path to the ROOT file with the simulation data
///\param name (TString) - name of the object
///\param iter (Int_t) - number of iterations for MLEM (not used yet)
///\param verbose (Bool_t) - verbose level for print-outs on the screen.
///\param dimZ (Double_t) - size of image plane in z-axis direction 
///\param dimY (Double_t) - size of image plane in y-axis direction 
///\param nbinsz (Int_t) - numbers of bins of image plane in z-axis direction 
///\param nbinsy (Int_t) - numbers of bins of image plane in y-axis direction 
CCMLEM::CCMLEM(TString inputName, TString name, Int_t iter, Bool_t verbose, Double_t dimZ, Double_t dimY, Int_t nbinsz, Int_t nbinsy){
  
  SetInputName(inputName);
  SetName(name);
  SetIter(iter);
  fVerbose = verbose;
  fDimZ = dimZ;
  fDimY = dimY;
  fXofRecoPlane = 0;
  fNbinsZ = nbinsz;
  fNbinsY = nbinsy;
  fPixelSizeZ = fDimZ/fNbinsZ;
  fPixelSizeY = fDimY/fNbinsY;
  
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
///TClonesArray is used to be filled by 
///intersection points and their global bin number of ComptonCones with the image plane
  fArray = new TClonesArray("IsectionPoint",10000);
  fNIpoints = 0;
  
  cout<<"CCMLEM: I will work in fVerbose="<<fVerbose<<" mode"<<endl;

}
//----------------------------------------
///Default destructor
CCMLEM::~CCMLEM(){
  if(fFile) fFile->Close();
}
//------------------------------------------
///Function for image reconstruction and returns kTRUE when a reconstructed cone 
///is bulit.Sets 2D image histogram and 
///1D histogram (graph) of non-pixelized reconstructed image. 
///By calling ComptonCones from CCReconstruction class, In the loop it reconstructs
///requested numbers of ComptonCones and uses the numerical method 
///to get intersection of ComptonCone with the image plane.
/// For each ComptonCone, the following processes will be done: 
///cases: (1) To fill TClonesArray by intersection points and global bin number,
///AddIsectionPoint() is used (in the following I will describe this function).
///(2) The global bin numbers are sorted, then for each touched pixel
///the length of track is calculated to estimate the reconstructed cone.
///\param iStart (Int_t) - number of first event for image reconstruction
///\param iStop (Int_t) - number of last event for image reconstruction
Bool_t CCMLEM::Reconstruct(Int_t iStart,Int_t iStop){

  CCReconstruction *reco = new CCReconstruction(fInputName, fName, fIter, fVerbose);

//image histogram
  fImage = new TH2F("image","image",fNbinsZ,-fDimZ/2.,fDimZ/2.,fNbinsY,-fDimY/2.,fDimY/2.);
  fImage->GetXaxis()->SetTitle("z [mm]");
  fImage->GetYaxis()->SetTitle("y [mm]");

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
  Double_t abit = 1.E-3;

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
  //IsectionPoint* tmppoint;
  IsectionPoint* tmppoint1;
  IsectionPoint* tmppoint2;

  const Double_t maxdist = sqrt(pow(fPixelSizeY,2)+pow(fPixelSizeZ,2));
  TStopwatch t;
  t.Start(); 
  for(Int_t i=iStart; i<iStop; i++){
    fNIpoints = 0; 
    //if(fVerbose) 
    cout<<"CCMLEM::Reconstruct(...) event "<< i<<endl<<endl;;
    ComptonCone *cone = reco->ReconstructCone(i);
    interactionPoint = cone->GetApex();
    coneAxis = cone->GetAxis();
    uvec = coneAxis.Unit();
    coneTheta = cone->GetAngle();
    Double_t a, b, c1, c2, c3, c4, c, z1, z2;
    Double_t d, e, f1, f2, f3, f4, f, y1, y2;  
    
    y = -A;
    if(fVerbose) cout<<"Loop over horizontal lines..."<<endl;
    for (j=0; j<=fNbinsY; j++){
      
      a = 2*(pow(-uvec.Z(),2) - pow(cos(coneTheta),2));
      
      b = 2*((-uvec.Z())*((interactionPoint.Z())*(-uvec.Z()) - y*((-uvec.Y())) + (interactionPoint.Y())*
			  (-uvec.Y()) + (interactionPoint.X())*(-uvec.X()) + F - (-uvec.X())*F) - (interactionPoint.Z())*pow(cos(coneTheta),2));
      
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
	((-1) + 2*pow(cos(coneTheta),2));
      
      
      c = (-2)*pow(cos(coneTheta),2)*(c1 + c2 + c3 + c4); 
      
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
      
      d = 2*(pow(-uvec.Y(),2) - pow(cos(coneTheta),2));
      
      e = 2*((-uvec.Y())*((interactionPoint.Z())*(-uvec.Z()) - z*((-uvec.Z())) + (interactionPoint.Y())*
			  (-uvec.Y()) + (interactionPoint.X())*(-uvec.X()) + F - (-uvec.X())*F) - (interactionPoint.Y())*pow(cos(coneTheta),2));
      
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
	((-1) + 2*pow(cos(coneTheta),2));
      
      f = (-2)*pow(cos(coneTheta),2)*(f1 + f2 + f3 + f4);
      
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
    IsectionPoint *temp;
    
    for(Int_t i=0; i<fNIpoints; i++){
       temp = (IsectionPoint*)fArray->At(i);
       fA[i] = temp->GetBin();
    }

    TMath::Sort(fNIpoints, fA, index, kFALSE);
    
    TVector3 *tmpvec1;
    TVector3 *tmpvec2;
    Double_t dist;
    Int_t binno1, binno2;
    
    for(int i=0; i<fNIpoints-1; i=i+2) {
     
      cout<<" index["<<i<<"]="<<index[i]<<", index["<<i+1<<"]="<<index[i+1]<<endl;
      tmppoint1 = (IsectionPoint*)fArray->At(index[i]);
      tmppoint2 = (IsectionPoint*)fArray->At(index[i+1]);
      if(tmppoint1==NULL || tmppoint2==NULL){
	cout<<"Something went wrong"<<tmppoint1<<"\t"<<tmppoint2<<endl;
      }
      tmpvec1 = tmppoint1->GetPointCoordinates();
      binno1 = tmppoint1->GetBin();
      tmpvec2 = tmppoint2->GetPointCoordinates();
      binno2 = tmppoint2->GetBin();
      cout<<" binno1="<<binno1<<", binno2="<<binno2<<endl<<endl;
      dist = ((*tmpvec1)-(*tmpvec2)).Mag();

      if(dist > maxdist){
        cout<<"Event "<<i<<": distance exceeds pixel diagonal "<<dist/maxdist<<" times"<<endl;
	continue;
      }
      if(binno1!=binno2){
        cout<<binno1<<"!="<<binno2<<" ->Bin numbers are different when they should not!"<<endl;
	i--;
      }
      
      fImage->SetBinContent(binno1, fImage->GetBinContent(binno1) + dist);
    }
    //if(fVerbose) cout<<"end of loop"<<endl;
    
    delete cone;
    cout<<"---------------------"<<endl;
  }// end of loop over events
  
  fArray->Clear("C");
  t.Stop(); 
  t.Print();
  SaveHistogram(fImage);
  //SaveToFile(fGraph);
  delete reco;
   
  return kTRUE;
}

//------------------------------------
///Function to check direction of intersecting lines and get coordinate of intersection points.
///TClonesArray is filled by coordinate of intersection points and their global bin number.
///To get intersection points and global bin number, adding and subtracting 
///a bit value (half value of pixel size) make us sure each point belongs to which two of pixels. 
/// So this function returns two points added for each pixel if all conditions are met.
///\param dir (TString) - direction of intersecting lines marched to 
///get intersection points on pixel borders
///\param x (Double_t) - x-component of intersection point on the image plane
///\param y (Double_t) - y-component of intersection point on the image plane
///\param z (Double_t) - z-component of intersection point on the image plane
Int_t CCMLEM::AddIsectionPoint(TString dir, Double_t x, Double_t y, Double_t z){
  if(fVerbose) cout<<dir<<"\t"<<x<<"\t"<<y<<"\t"<<z<<endl;
  dir.ToLower();
  if(dir!="hor" && dir!="ver"){
    if(fVerbose) cout<<"Unknown direction of intrsecting line: "<<dir<<endl;
    return 0;
  }
  if(fabs(y)>fDimY/2. || fabs(z)>fDimZ/2.){
    //if(fVerbose) cout<<"point outside of image range..."<<endl;
    return 0;
  }
    
  IsectionPoint* tmppoint;
  Int_t added = 0;
  Int_t pixelZ, pixelY;
  Double_t yplus, yminus;
  Double_t zplus, zminus;
  if(dir=="hor"){ // adding point from intersections with horizontal lines
    pixelZ = fImage->GetXaxis()->FindBin(z);
    if(pixelZ>fNbinsZ) pixelZ=fNbinsZ; //inclusion of upper edges of histo
    yplus = y +0.5*fPixelSizeY;
    yminus = y -0.5*fPixelSizeY;

    if(fabs(yplus)<=fDimY/2){
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);
      pixelY = fImage->GetYaxis()->FindBin(yplus);
      tmppoint->SetBinPoint(fImage->GetBin(pixelZ,pixelY), x, y, z);
      //if(fVerbose) tmppoint->Print();
      //fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
    }
    if(fabs(yminus)<=fDimY/2){
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);
      pixelY = fImage->GetYaxis()->FindBin(yminus);
      tmppoint->SetBinPoint(fImage->GetBin(pixelZ,pixelY), x, y, z);
      //if(fVerbose) tmppoint->Print();
      //fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
    }
  }
  if(dir=="ver"){ // adding point from intersections with vertical lines
    pixelY = fImage->GetYaxis()->FindBin(y);
    if(pixelY>fNbinsY) pixelY=fNbinsY; //inclusion of upper edges of histo
    zplus = z +0.5*fPixelSizeZ;
    zminus = z -0.5*fPixelSizeZ;
    
    if(fabs(zplus)<=fDimZ/2){
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);
      pixelZ = fImage->GetXaxis()->FindBin(zplus);
      tmppoint->SetBinPoint(fImage->GetBin(pixelZ,pixelY), x, y, z);
      //if(fVerbose) tmppoint->Print();
      //fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
    }
    if(fabs(zminus)<=fDimZ/2){
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);
      pixelZ = fImage->GetXaxis()->FindBin(zminus);
      tmppoint->SetBinPoint(fImage->GetBin(pixelZ,pixelY), x, y, z);
      //if(fVerbose) tmppoint->Print();
      //fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
    }
  }
  //if(fVerbose)  cout<<added<<" points added..."<<endl<<endl;;
  return added;
}

//------------------------------------
///Saves 2D histogram in the ROOT file. Name of the file is based on the name 
///of the CCMLEM object.
Bool_t CCMLEM::SaveHistogram(TH2F *h){
  TString name = "../sources/results/" + fName + ".root";
  TFile *file = new TFile(name,"UPDATE");
  h->Write();
  file->Close();
  if(fVerbose) cout << "\nHistogram " << h->GetName() << 
                       " saved in the file " << name << endl;
  return kTRUE;
} 
//------------------------------------
///Saves 1D histogram (graph) in the ROOT file. Name of the file is 
///based on the name of the CCMLEM object.
Bool_t CCMLEM::SaveToFile(TGraph *h){
  TString name = "../sources/results/" + fName + ".root";
  TFile *file = new TFile(name,"UPDATE");
  h->Write();
  file->Close();
  if(fVerbose) cout << "\nObject " << h->GetName() << 
                       " saved in the file " << name << endl;
  return kTRUE;
} 



