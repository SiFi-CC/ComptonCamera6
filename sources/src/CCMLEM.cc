#include "CCMLEM.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include "CCReconstruction.hh"
#include "IsectionPoint.hh"
#include "TClonesArray.h"

using namespace std;

ClassImp(CCMLEM);



//--------------------

CCMLEM::CCMLEM(TString inputName, TString name, Int_t iter, Bool_t verbose, Double_t dimZ, Double_t dimY){
  
  SetInputName(inputName);
  SetName(name);
  SetIter(iter);
  fVerbose = verbose;
  fDimZ = dimZ;
  fDimY = dimY;
  fXofRecoPlane = 0;
  
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
  
  fArray = new TClonesArray("IsectionPoint",10000);

}
//----------------------------------------
CCMLEM::~CCMLEM(){
  if(fFile) fFile->Close();
}
//------------------------------------------

Bool_t CCMLEM::Reconstruct(Int_t iStart,Int_t iStop){

  CCReconstruction *reco = new CCReconstruction(fInputName, fName, fIter, fVerbose);

//image histogram
  Int_t nbinsz = (Int_t)fDimZ/10;
  Int_t nbinsy = (Int_t)fDimY/10;
  Double_t zlimit = fDimZ/2.;
  Double_t ylimit = fDimY/2.;
  fImage = new TH2F("image","image",nbinsz,-zlimit,zlimit,nbinsy,-ylimit,ylimit);
  g = new TGraph();
  g->SetName("g");
  g->SetTitle("non-pixelized reco image");
  g->SetMarkerStyle(7);
  g->SetMarkerColor(kRed);
  fImage->GetXaxis()->SetTitle("z [mm]");
  fImage->GetYaxis()->SetTitle("y [mm]");
 
  
  //TFile *file = new TFile("myTree.root","RECREATE");
  //TTree *tree = new TTree("intersection","points");
  //tree->Branch("intersectionPoints",&intersectionPoints,"x/D:y/D:z/D");
  Double_t pixelX = 0.;
  Double_t pixelY, pixelZ;
  Double_t abit = 1.E-3;
  Double_t pixelSizeZ = zlimit*2/nbinsz;	//mm
  Double_t pixelSizeY = ylimit*2/nbinsy;
  
  Int_t j,k;
  
  TVector3 interactionPoint;
  TVector3 coneAxis;
  Double_t coneTheta;
  TVector3 uvec;
 
  Double_t F = fXofRecoPlane;
  Double_t z,y;
  Double_t A = fDimY/2.;  
  Double_t B = fDimZ/2.;
  
  
  fArray->Clear("C");
  Int_t npoints = 0;
  IsectionPoint* tmppoint;
  IsectionPoint* tmppoint1;
  IsectionPoint* tmppoint2;

  for(Int_t i=iStart; i<iStop; i++){
 
    cout<<"Recontruction for event "<< i<<endl<<endl;;
    ComptonCone *cone = reco->ReconstructCone(i);
    interactionPoint = cone->GetApex();
    coneAxis = cone->GetAxis();
    uvec = coneAxis.Unit();
    coneTheta = cone->GetAngle();
      
      
    y = -A;
    
   // cout<<"pixelSizeY = "<<pixelSizeY<<", pixelSizeZ = "<<pixelSizeZ<<endl;
   //cout<<"nbinsy = "<<nbinsy<<", nbinsz = "<<nbinsz<<endl;
  cout<<"Loop over horizontal lines..."<<endl;
    for (j=0; j<=nbinsy; j++){

      Double_t a = 2*(pow(-uvec.Z(),2) - pow(cos(coneTheta),2));

      Double_t b = 2*((-uvec.Z())*((interactionPoint.Z())*(-uvec.Z()) - y*((-uvec.Y())) + (interactionPoint.Y())*
                   (-uvec.Y()) + (interactionPoint.X())*(-uvec.X()) + F - (-uvec.X())*F) - (interactionPoint.Z())*pow(cos(coneTheta),2));
	
      Double_t c1 = pow(interactionPoint.Y(),2) + pow(interactionPoint.X(),2) + ((-2)*pow(interactionPoint.Y(),2)*pow(-uvec.Z(),2)) +  
                    ((-2)*pow(interactionPoint.X(),2)*pow(-uvec.Z(),2)) + ((-2)*pow(interactionPoint.Y(),2)*pow(-uvec.Y(),2)) + 
                    ((-4)*(interactionPoint.Y())*(interactionPoint.X())*(-uvec.Y())*(-uvec.X())) + ((-2)*pow(interactionPoint.X(),2)*
		    pow(-uvec.X(),2));
		      
      Double_t c2 = (-2)*(interactionPoint.X())*F + 4*(interactionPoint.X())*F*
		    pow(-uvec.Z(),2) - 4*(interactionPoint.Y())*(-uvec.Y())*F - 4*(interactionPoint.X())*
		    (-uvec.X())*F + 4*(interactionPoint.Y())*(-uvec.Y())*(-uvec.X())*F + 
		    4*(interactionPoint.X())*pow(-uvec.X(),2)*F - pow(F,2) - 2*pow(-uvec.Z(),2)*
		    pow(F,2) + 4*(-uvec.X())*pow(F,2) - 2*pow(-uvec.X(),2)*pow(F,2);
	
      Double_t c3 = pow(y,2)*(1 - 2*pow(-uvec.Z(),2) - 2*pow(-uvec.Y(),2)) + y*
		    ((interactionPoint.Y())*((-2) + 4*pow(-uvec.Z(),2) + 4*pow(-uvec.Y(),2)) + 4*(-uvec.Y())*((interactionPoint.X())*(-uvec.X()) + 
                    F - (-uvec.X())*F));
	
      Double_t c4 = (pow(y,2) - 2*y*(interactionPoint.Y()) + pow(interactionPoint.Y(),2) + pow((interactionPoint.X() - F),2))*
                    ((-1) + 2*pow(cos(coneTheta),2));
	

      Double_t c = (-2)*pow(cos(coneTheta),2)*(c1 + c2 + c3 + c4); 
	
        if(c >= 0){
      
          Double_t z1 = (b - sqrt(c))/a;
          Double_t z2 = (b + sqrt(c))/a;
        
          
	  
	  
	  
	  
	  
	  if(fabs(z1) <= zlimit){
	    
	    if(fabs(y+abit) <= ylimit){
	    
	      tmppoint = (IsectionPoint*)fArray->ConstructedAt(npoints++);
	      pixelZ = fImage->GetXaxis()->FindBin(z1);
	      
              pixelY = fImage->GetYaxis()->FindBin(y+abit);
	      tmppoint->SetBinPoint(fImage->GetBin(pixelZ,pixelY), fXofRecoPlane, y, z1);
	    }
	   
	    if(fabs(y-abit) <= ylimit){
	      tmppoint = (IsectionPoint*)fArray->ConstructedAt(npoints++);
	      pixelZ = fImage->GetXaxis()->FindBin(z1);
	      
              pixelY = fImage->GetYaxis()->FindBin(y-abit);
	      tmppoint->SetBinPoint(fImage->GetBin(pixelZ,pixelY), fXofRecoPlane, y, z1);
	    }
	    g->SetPoint(g->GetN(), z1, y);
	    //fImage->Fill(z1,y);
	    cout<<"j="<<j<<", y="<<y<<", z="<<z1<<endl;
	  }
	  if(fabs(z2) <= zlimit){
	    
	    if(fabs(y+abit) <= ylimit){
	    
	      tmppoint = (IsectionPoint*)fArray->ConstructedAt(npoints++);
	      pixelZ = fImage->GetXaxis()->FindBin(z2);
	    
              pixelY = fImage->GetYaxis()->FindBin(y+abit);
	      tmppoint->SetBinPoint(fImage->GetBin(pixelZ,pixelY), fXofRecoPlane, y, z2);
	    }
	    
	    if(fabs(y-abit) <= ylimit){
	      tmppoint = (IsectionPoint*)fArray->ConstructedAt(npoints++);
	      pixelZ = fImage->GetXaxis()->FindBin(z2);
	      
              pixelY = fImage->GetYaxis()->FindBin(y-abit);
	      tmppoint->SetBinPoint(fImage->GetBin(pixelZ,pixelY), fXofRecoPlane, y, z2);
	    }
	    g->SetPoint(g->GetN(), z2, y);
	    //fImage->Fill(z2,y);
	    cout<<"j="<<j<<", y="<<y<<", z="<<z2<<endl;
	  }
      }

      y=y+pixelSizeY;
     
     
    // tree->Fill();

    } //end of loop over horizontal lines
   
  
   z = -B;
    
   //cout<<"Loop over vertical lines..."<<endl;
    for (k=0; k<=nbinsz; k++){

      Double_t d = 2*(pow(-uvec.Y(),2) - pow(cos(coneTheta),2));

      Double_t e = 2*((-uvec.Y())*((interactionPoint.Z())*(-uvec.Z()) - z*((-uvec.Z())) + (interactionPoint.Y())*
                    (-uvec.Y()) + (interactionPoint.X())*(-uvec.X()) + F - (-uvec.X())*F) - (interactionPoint.Y())*pow(cos(coneTheta),2));

      Double_t f1 = pow(interactionPoint.Z(),2) + pow(interactionPoint.X(),2) + ((-2)*pow(interactionPoint.Z(),2)*pow(-uvec.Z(),2)) +  
                    ((-2)*pow(interactionPoint.Z(),2)*pow(-uvec.Y(),2)) + ((-2)*pow(interactionPoint.X(),2)*pow(-uvec.Y(),2)) + 
                    ((-4)*(interactionPoint.Z())*(interactionPoint.X())*(-uvec.Z())*(-uvec.X())) + ((-2)*pow(interactionPoint.X(),2)*
		    pow(-uvec.X(),2));
		      
      Double_t f2 = (-2)*(interactionPoint.X())*F + 4*(interactionPoint.X())*F*
		    pow(-uvec.Y(),2) - 4*(interactionPoint.Z())*(-uvec.Z())*F - 4*(interactionPoint.X())*
		    (-uvec.X())*F + 4*(interactionPoint.Z())*(-uvec.Z())*(-uvec.X())*F + 
		    4*(interactionPoint.X())*pow(-uvec.X(),2)*F - pow(F,2) - 2*pow(-uvec.Y(),2)*
		    pow(F,2) + 4*(-uvec.X())*pow(F,2) - 2*pow(-uvec.X(),2)*pow(F,2);
	
      Double_t f3 = pow(z,2)*(1 - 2*pow(-uvec.Z(),2) - 2*pow(-uvec.Y(),2)) + z*
		    ((interactionPoint.Z())*((-2) + 4*pow(-uvec.Z(),2) + 4*pow(-uvec.Y(),2)) + 4*(-uvec.Z())*((interactionPoint.X())*(-uvec.X()) + 
                    F - (-uvec.X())*F));
	
      Double_t f4 = (pow(z,2) - 2*z*(interactionPoint.Z()) + pow(interactionPoint.Z(),2) + pow((interactionPoint.X() - F),2))*
                    ((-1) + 2*pow(cos(coneTheta),2));
       
      Double_t f = (-2)*pow(cos(coneTheta),2)*(f1 + f2 + f3 + f4);

          if(f >= 0){

            Double_t y1 = (e - sqrt(f))/d;
            Double_t y2 = (e + sqrt(f))/d;
	   
	      
	   
	    
	    
	 
	  
	  
	   
	   
	  
	   if(fabs(y1) <= ylimit){
	     
	     if(fabs(z+abit) <= zlimit){
	     
	       tmppoint = (IsectionPoint*)fArray->ConstructedAt(npoints++);
	       pixelY = fImage->GetYaxis()->FindBin(y1);
	       
	       pixelZ = fImage->GetXaxis()->FindBin(z+abit);
	       tmppoint->SetBinPoint(fImage->GetBin(pixelZ,pixelY), fXofRecoPlane, y1, z);
	     }
	     
	     if(fabs(z-abit) <= zlimit){
	       tmppoint = (IsectionPoint*)fArray->ConstructedAt(npoints++);
	       pixelY = fImage->GetYaxis()->FindBin(y1);
	       
	       pixelZ = fImage->GetXaxis()->FindBin(z-abit);
               
	       tmppoint->SetBinPoint(fImage->GetBin(pixelZ,pixelY), fXofRecoPlane, y1, z);
	     }
	     g->SetPoint(g->GetN(), z, y1);
	     //fImage->Fill(z,y1);
	     cout<<"k="<<k<<", y="<<y1<<", z="<<z<<endl;
	   }
	   
	  if(fabs(y2) <= ylimit){
	    
	    if(fabs(z+abit) <= zlimit){
	    
	      tmppoint = (IsectionPoint*)fArray->ConstructedAt(npoints++);
	      
              pixelY = fImage->GetYaxis()->FindBin(y2);
	      
	      pixelZ = fImage->GetXaxis()->FindBin(z+abit);
	      tmppoint->SetBinPoint(fImage->GetBin(pixelZ,pixelY), fXofRecoPlane, y2, z);
	    }
	    
	    if(fabs(z-abit) <= zlimit){
	      tmppoint = (IsectionPoint*)fArray->ConstructedAt(npoints++);
	      
              pixelY = fImage->GetYaxis()->FindBin(y2);
	      
	      pixelZ = fImage->GetXaxis()->FindBin(z-abit);
	      tmppoint->SetBinPoint(fImage->GetBin(pixelZ,pixelY), fXofRecoPlane, y2, z);
	    }
	    g->SetPoint(g->GetN(), z, y2);
	    //fImage->Fill(z,y2);
	    cout<<"k="<<k<<", y="<<y2<<", z="<<z<<endl;
	  }
	  
      }
      
      z=z+pixelSizeZ;
      
   
  // tree->Fill();

    } //end of loop over vertical lines
    

    fArray->Sort();
    
   cout<<"Print fArray after sorting..."<<endl;
    for(int i=0; i<npoints; i++){
      cout<<"i="<<i<<endl;
      ((IsectionPoint*)fArray->At(i))->Print();
    }
    cout<<"End of print fArray after sorting..."<<endl;

    TVector3 *tmpvec1;
    TVector3 *tmpvec2;
    Double_t dist;
    Int_t binno1, binno2;
    //cout<<npoints<<endl;
    for(int i=0; i<npoints; i=i+2){
      
      tmppoint1 = (IsectionPoint*)fArray->At(i);
      //cout<<"tmppoint1="<<tmppoint1<<endl;
      tmpvec1 = tmppoint1->GetPointCoordinates();
      binno1 = tmppoint1->GetBin();
      //cout<<"binno1="<<binno1<<endl;
      tmppoint2 = (IsectionPoint*)fArray->At(i+1);
     // cout<<"tmppoint2="<<tmppoint2<<endl;
      tmpvec2 = tmppoint2->GetPointCoordinates();
      binno2 = tmppoint2->GetBin();
      //cout<<"binno2="<<binno2<<endl;
      
      if(dist <= sqrt(pow(pixelSizeY,2)+pow(pixelSizeZ,2))){
	dist = ((*tmpvec1)-(*tmpvec2)).Mag();
        //cout<<i<<"  "<<dist<<endl;
      }
      if(binno1!=binno2){
        cout<<"Bin numbers are different when they should not!"<<endl;
        //return kFALSE;
      }
      
      
      //cout<<"\n\n\n"<<endl;
      fImage->SetBinContent(binno1, fImage->GetBinContent(binno1) + dist);
      //fImage->SetBinContent(binno2, fImage->GetBinContent(binno1) + dist);
    }
    cout<<"end of loop"<<endl;
    
    delete cone;
  }// end of loop over events
   
  SaveHistogram(fImage);
  SaveToFile(g);
  //tree->Write();
  delete reco;
   
  return kTRUE;
}


//------------------------------------
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
Bool_t CCMLEM::SaveToFile(TGraph *h){
  TString name = "../sources/results/" + fName + ".root";
  TFile *file = new TFile(name,"UPDATE");
  h->Write();
  file->Close();
  if(fVerbose) cout << "\nObject " << h->GetName() << 
                       " saved in the file " << name << endl;
  return kTRUE;
} 



