#include "CCSimulation.hh" 
#include "TMath.h"
#include "TRandom.h"
#include <iostream>
#include <fstream>

using namespace std;

ClassImp(CCSimulation);

//------------------------------------------------------------------
CCSimulation::CCSimulation(){
  SetName("simulation");
  fVerbose = kTRUE;
  fFile = 0;
  fTree = 0;
  fNev = 0;
}
//------------------------------------------------------------------
CCSimulation::CCSimulation(TString name, Bool_t verbose){
  
  SetName(name);
  fVerbose = verbose;
  fFile = new TFile("../sources/results/"+name+".root","RECREATE");
  fTree = new TTree("data","data");
  fTree->Branch("point0",&fPoint0);	//source position
  fTree->Branch("point1",&fPoint1);	//interaction point on the scaterrer
  fTree->Branch("point2",&fPoint2);	//interaction point on the absorber
  fTree->Branch("versor1",&fVersor1);	//direction of gamma quantum comming from the source
  fTree->Branch("versor2",&fVersor2);	//direction of scattered gamma quantum
  fTree->Branch("energy0",&fEnergy0);	//energy of gamma from the source
  fTree->Branch("energy1",&fEnergy1);	//energy loss
  fTree->Branch("energy2",&fEnergy2);	//energy after scattering
  fNev = 0;
  fGenVersion = -1000;
  
  Double_t size = 200.;
  Int_t nbins = 100;
  hSource = new TH2F("hSource","hSource",nbins,-size,size,nbins,-size,size);
  hSource->GetXaxis()->SetTitle("z [mm]");
  hSource->GetYaxis()->SetTitle("y [mm]");
  hScat = new TH2F("hScat","hScat",nbins,-size,size,nbins,-size,size);
  hScat->GetXaxis()->SetTitle("z [mm]");
  hScat->GetYaxis()->SetTitle("y [mm]");
  hAbs = new TH2F("hAbs","hAbs",nbins,-size,size,nbins,-size,size);
  hAbs->GetXaxis()->SetTitle("z [mm]");
  hAbs->GetYaxis()->SetTitle("y [mm]");
  hEnergy = new TH1F("hEnergy","energy of scattered gamma",nbins,0.,5.);
  hEnergy->GetXaxis()->SetTitle("energy [MeV]");
  
}
//------------------------------------------------------------------
CCSimulation::~CCSimulation(){
  if(fVerbose) cout << "Inside destructor of the CCSimulation class" << endl;
  if(fFile && fTree){
    fTree->Write();
    hSource->Write();
    hScat->Write();
    hAbs->Write();
    hEnergy->Write();
    fFile->Close(); 
  }
  SaveGeometryTxt();
  BuildTGeometry();
}
//------------------------------------------------------------------
void CCSimulation::BuildSetup(Double_t scatDist = 50, //mm
			      Double_t scatZ    = 80,
			      Double_t scatY    = 80,
                              Double_t absDist  = 60,
			      Double_t absZ     = 300,
			      Double_t absY     = 300){
  
  if(fVerbose) cout << "\n----- Building the scatterer..." << endl;
  fScatterer.SetPlane(1,0,0,scatDist);
  fScatterer.SetDimensions(scatZ,scatY);
  fScatterer.SetName("scatterer");
  if(fVerbose) fScatterer.Print();
  
  if(fVerbose) cout << "\n----- Building the absorber..." << endl;
  fAbsorber.SetPlane(1,0,0,absDist);
  fAbsorber.SetDimensions(absZ,absY);
  fAbsorber.SetName("absorber");
  if(fVerbose) fAbsorber.Print();
}
//------------------------------------------------------------------
Bool_t CCSimulation::GenerateRay(void){
  
  Double_t theta, phi, maxz;
  fYgap = -20.;
  fZgap = -20.;
  fRadius = 15.;
  Double_t theta_circ, phi_circ;
  Double_t x,y,z,rad;
  TVector3 circ;
  
  switch(fGenVersion){
    case 1:	//isotropic point-like source
      fPoint0.SetXYZ(0,0,0);
      theta = acos(gRandom->Uniform(-1,1));			//rad
      phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());	//rad
      //Double_t theta = TMath::Pi()/2.;	//rad
      //Double_t phi = TMath::Pi();		//rad
      break;
    case 2:	//uniform distribution along z axis (beam)
      maxz = fScatterer.GetDimZ()/2.;
      fPoint0.SetXYZ(0,0,gRandom->Uniform(-maxz,maxz));
      theta = acos(gRandom->Uniform(-1,1));			//rad
      phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());	//rad
      break;
    case 3:	//two isotropic point-like sources along z
      if(gRandom->Uniform(0,1)<0.5)
	fPoint0.SetXYZ(0,0,0);
      else
	fPoint0.SetXYZ(0,0,fZgap);
      theta = acos(gRandom->Uniform(-1,1));			//rad
      phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());	//rad
      break;
    case 4:	//two isotropic point-like sources along y
      if(gRandom->Uniform(0,1)<0.5)
	fPoint0.SetXYZ(0,0,0);
      else 
	fPoint0.SetXYZ(0,fYgap,0);
      theta = acos(gRandom->Uniform(-1,1));			//rad
      phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());	//rad
      break;
    case 5:	//circular source of given radius
      phi_circ = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
      theta_circ = acos(gRandom->Uniform(-1,1));
      rad = gRandom->Uniform(0,fRadius);
      fPoint0.SetMagThetaPhi(rad,theta_circ,phi_circ);
      theta = acos(gRandom->Uniform(-1,1));			//rad
      phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());	//rad
      break;
    default:
      cout << "##### Please choose correcr version of the generaror!" << endl;
      break;
  }
  
  fVersor1.SetMagThetaPhi(1,theta,phi);
  
  return kTRUE;
}
//------------------------------------------------------------------
Bool_t CCSimulation::ProcessEvent(void){
 
  Clear();
  GenerateRay();

  fEnergy0 = 4.44;
  fTrack1.SetPoint(fPoint0);
  fTrack1.SetVersor(fVersor1);
  fTrack1.SetEnergy(fEnergy0);
  
  Bool_t scatFlag = fTrack1.FindCrossPoint(&fScatterer,fPoint1);
  if(scatFlag == kFALSE){
    if(fVerbose) cout << "\tNo cross points with the scatterer\n" << endl;
    return kFALSE;
  }
  
  fTrack2 = fPhysics.ComptonScatter(&fTrack1,&fScatterer);
  fVersor2 = fTrack2->GetVersor();
  fEnergy2 = fTrack2->GetEnergy();
  fEnergy1 = fEnergy0 - fEnergy2;
  
  //----- energy check
  Double_t en = fEnergy1 + fEnergy2;
  if(fabs(fEnergy0-en)>1.E-8){
    cout << "\t##### Error in CCSimulation::ProcessEvent!" << endl;
    cout << "\t##### Energies of tracks do not sum correctly! Please check!" << endl;
    return kFALSE;
  }
  //----- end of the energy check
  
  //----- position check
  TVector3 point = fTrack2->GetPoint();
  if(point!=fPoint1){
    cout << "\t##### Error in CCSimulation::ProcessEvent!" << endl;
    cout << "\t##### Cross point with the scaterrer is different than starting point of the new track" << endl;
    return kFALSE;
  }
  //----- end of the position check
  
  Bool_t absFlag = fTrack2->FindCrossPoint(&fAbsorber,fPoint2);
  if(absFlag == kFALSE){
     if(fVerbose) cout << "\tNo cross point with the absorber\n" << endl;
     return kFALSE;
  }

  fTree->Fill();
  hSource->Fill(fPoint0.Z(),fPoint0.Y());
  hScat->Fill(fPoint1.Z(),fPoint1.Y());
  hAbs->Fill(fPoint2.Z(),fPoint2.Y());
  hEnergy->Fill(fEnergy2);
  fNev++;
  
  return kTRUE;
}
//------------------------------------------------------------------
void CCSimulation::Loop(Int_t nev){
  Int_t i=0;
  while(fNev<nev){
    ProcessEvent();
    i++;
    if(i%1000==0) cout << "i = " << i << " events simulated, " 
		       << fNev << " events in acceptance" << endl;
  }
} 
//------------------------------------------------------------------ 
void CCSimulation::Clear(void){
  fPoint0.SetXYZ(-1000,-1000,-1000);
  fPoint1.SetXYZ(-1000,-1000,-1000);
  fPoint2.SetXYZ(-1000,-1000,-1000);
  fVersor1.SetXYZ(-1000,-1000,-1000);
  fVersor2.SetXYZ(-1000,-1000,-1000);
  fEnergy0 = -1000;
  fEnergy1 = -1000;
  fEnergy2 = -1000;
}
//------------------------------------------------------------------ 
void CCSimulation::SaveGeometryTxt(void){
  ofstream output(Form("../sources/results/CCSimulation_geometry_gen%i.txt",fGenVersion), std::ios::out | std::ios::trunc);
  output << "Generator version: " << fGenVersion << endl;
  if(fGenVersion==1)
    output << fPoint0.X() << " " << fPoint0.Y() << " " << fPoint0.Z() << endl;
  else if(fGenVersion==2)
    output << "-" << fScatterer.GetDimZ()/2. << " " << fScatterer.GetDimZ()/2. << endl;
  else if(fGenVersion==3)
    output << "0 0 0\t0 0 " << fZgap << endl;
  else if(fGenVersion==4)
    output << "0 0 0\t0 " << fYgap << " 0" << endl;
  else if(fGenVersion==5)
    output << "Center: 0 0 0 \t Radius: " << fRadius << endl;
  else
    cout << "##### Please choose correcr version of the generaror!" << endl;
  
  output << "DetPlane object: " << fScatterer.GetName() << endl;
  output << "fA = " << fScatterer.GetA() << endl;
  output << "fB = " << fScatterer.GetB() << endl;
  output << "fC = " << fScatterer.GetC() << endl;
  output << "fD = " << fScatterer.GetD() << endl;
  output << "fDimZ = " << fScatterer.GetDimZ() << endl;
  output << "fDimY = " << fScatterer.GetDimY() << endl;
  output << "DetPlane object: " << fAbsorber.GetName() << endl;
  output << "fA = " << fAbsorber.GetA() << endl;
  output << "fB = " << fAbsorber.GetB() << endl;
  output << "fC = " << fAbsorber.GetC() << endl;
  output << "fD = " << fAbsorber.GetD() << endl;
  output << "fDimZ = " << fAbsorber.GetDimZ() << endl;
  output << "fDimY = " << fAbsorber.GetDimY() << endl;
  output.close();
  if(fVerbose) cout << "Geometry of the setup saved in the text file" << endl;
}
//------------------------------------------------------------------ 
void CCSimulation::BuildTGeometry(void){
  
   TGeoManager *geom = new TGeoManager("world", "CCSimulation - world");
   
   //----- define materials
   TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum",0,0,0);
   matVacuum->SetTransparency(0);
   
   //----- define media
   TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1,matVacuum);
   
   //----- define translations
   Double_t scatA = fScatterer.GetA();
   Double_t scatD = fScatterer.GetD();
   Double_t absA = fAbsorber.GetA();
   Double_t absD = fAbsorber.GetD();
   Double_t scatTrans = -scatD/scatA;
   Double_t absTrans = -absD/absA;

   TGeoTranslation *trScat = new TGeoTranslation(scatTrans,0.,0.);  
   TGeoTranslation *trAbs = new TGeoTranslation(absTrans,0.,0.);
   
   //----- define dimensions of objects
   Double_t scatDimZ = fScatterer.GetDimZ()/2.;	//half length
   Double_t scatDimY = fScatterer.GetDimY()/2.;
   Double_t absDimZ = fAbsorber.GetDimZ()/2.;
   Double_t absDimY = fAbsorber.GetDimY()/2.;
   
   //----- make the top container volume
   Double_t worldX = -absTrans + 50.; //cm half length 
   Double_t worldY = absDimY + 50.;
   Double_t worldZ = absDimZ + 50.;
   TGeoVolume *top = geom->MakeBox("top",Vacuum,worldX,worldZ,worldY);
   geom->SetTopVolume(top);
   
   //----- make other objects
   TGeoVolume *scat = geom->MakeBox("scat",Vacuum,0.001,scatDimZ,scatDimY);
   TGeoVolume *abs = geom->MakeBox("abs",Vacuum,0.001,absDimZ,absDimY);
   top->AddNode(scat,1,trScat);
   top->AddNode(abs,1,trAbs);
   
   //----- colouring 
   scat->SetLineColor(kMagenta);
   abs->SetLineColor(kOrange+7);
   
   //-----setting the source (depending on the version of the generator)
   if(fGenVersion==1){			//isotropic point-like source 
       TGeoTranslation *tr1 = new TGeoTranslation(fPoint0.X(),fPoint0.Z(),fPoint0.Y());
       TGeoVolume *source1 = geom->MakeBox("source",Vacuum,0.1,0.1,0.1);
       top->AddNode(source1,1,tr1);
       source1->SetLineColor(kRed);
  }
     else if(fGenVersion==2){		//uniform distribution along z axis (beam)
       TGeoTranslation *tr2 = new TGeoTranslation(0.,0.,0.);
       TGeoVolume *source2 = geom->MakeBox("source",Vacuum,0.1,scatDimZ,0.1);
       top->AddNode(source2,1,tr2);
       source2->SetLineColor(kRed);
    }
     else if(fGenVersion==3){		//two isotropic point-like sources along z 
       TGeoTranslation *tr3 = new TGeoTranslation(0.,0.,0.);
       TGeoTranslation *tr4 = new TGeoTranslation(0.,fZgap,0.);
       TGeoVolume *source3 = geom->MakeBox("source",Vacuum,0.1,0.1,0.1);
       TGeoVolume *source4 = geom->MakeBox("source",Vacuum,0.1,0.1,0.1);
       top->AddNode(source3,1,tr3);
       top->AddNode(source4,1,tr4);
       source3->SetLineColor(kRed);
       source4->SetLineColor(kRed);
    }
     else if(fGenVersion==4){		//two isotropic point-like sources along y
       TGeoTranslation *tr5 = new TGeoTranslation(0.,0.,0.);
       TGeoTranslation *tr6 = new TGeoTranslation(0.,0.,fYgap);
       TGeoVolume *source5 = geom->MakeBox("source",Vacuum,0.1,0.1,0.1);
       TGeoVolume *source6 = geom->MakeBox("source",Vacuum,0.1,0.1,0.1);
       top->AddNode(source5,1,tr5);
       top->AddNode(source6,1,tr6);
       source5->SetLineColor(kRed);
       source6->SetLineColor(kRed);
    }
    else if(fGenVersion==5){
       TGeoTranslation *tr7 = new TGeoTranslation(0.,0.,0.);
       TGeoVolume *source7 = geom->MakeSphere("source",Vacuum,0.1,fRadius);
       top->AddNode(source7,1,tr7);
       source7->SetLineColor(kRed);
    }
     else
      cout << "##### Please choose correcr version of the generaror!" << endl;
   
   //----- close geometry and save
   geom->CloseGeometry();
   geom->SetVisLevel(4);
   geom->Export(Form("../sources/results/CCSimulation_TGeometry_gen%i.root",fGenVersion));
}
//------------------------------------------------------------------ 
void CCSimulation::Print(void){
 cout << "\nCCSimulation::Print() for object " << GetName() << endl; 
}
//------------------------------------------------------------------ 
