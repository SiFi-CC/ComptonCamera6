#include "CCMLEM.hh"
#include "IsectionPoint.hh"
#include "SMElement.hh"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TRandom.h"
#include "TStyle.h"
#include "TSystem.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "Math/GenVector/Plane3D.h"
#include "Math/Point3D.h"
#include <CmdLineConfig.hh>
using namespace std;

ClassImp(CCMLEM);

//--------------------
/// Standard constructor (recommended).
///\param path (TString) - full path to the configuration file.
CCMLEM::CCMLEM(TString path) {

  Clear();

  bool stat_config = ReadConfig(path);
  bool stat_reader = SetInputReader();
  if (!stat_config || !stat_reader) {
    throw "##### Exception in CCMLEM constructor!";
  }

  Size_t len = strlen(fInputName);
  TString fname = fInputName;
  fname.Insert(len - 5, "_MLEM");
  TString outputName = fname;
  outputName.ReplaceAll("CCSimulation", "CCReconstruction");
  TString option = (fFreshOutput ? "RECREATE" : "UPDATE");
  fOutputFile = new TFile(outputName, option);

  fArray = new TClonesArray("IsectionPoint", 10000);
  fSM = new TClonesArray("SMElement", 1000000);
  fNIpoints = 0;
  fPoints = 0;
/// To smear camera performance, this file is called for energy resolution.
/// It shows a function of energy deposited in a 10-cm-long LuAG(Ce) fiber 
///  with a square cross-section from Geant4 simulation   
  if (fSmear) {
    TString path = gSystem->Getenv("CC6DIR");
    TString name = path+"/share/ComptonCamera6/mlem_reco/EnergyResolutionExample.root";
     //TString name = path+"/share/ComptonCamera6/mlem_reco/Fiber_1_0.root";
    // TString name = path+"/share/ComptonCamera6/mlem_reco/Fiber_1_1.root";
    TFile* file = new TFile(name, "READ");
    fHisto = (TH1D*)file->Get("Scintillator_0ResolutionCombiEnergy");
     //fHisto = (TH1D*)file->Get("FiberEnergyUnc_Type3");
    // fHisto = (TH1D*)file->Get("FiberEnergyUnc_Type5");
    TF1* func = new TF1("fit1", "[0]+[1]/sqrt(x)+[2]/x^(3/2)", 0, 4.5);

//     func->FixParameter(0, -0.0178);
//     func->FixParameter(1, 0.05);
//     func->FixParameter(2, 0.01);
    fHisto->Fit("fit1","r");
//     std::cout << "par0 " << func->GetParameter(0) << std::endl;
//     std::cout << "par1 " << func->GetParameter(1) << std::endl;
//     std::cout << "par2 " << func->GetParameter(2) << std::endl;
  }
}
//--------------------
/// Default constructor.
CCMLEM::CCMLEM() {
  cout << "##### Warning in CCMLEM constructor!" << endl;
  cout << "You are using default constructor!" << endl;
  Clear();
}
//----------------------------------------
/// Default destructor.
CCMLEM::~CCMLEM() {
  if (fOutputFile) fOutputFile->Close();
}
//------------------------------------------
/// Based on the given input simulation file creates suitable
/// InputReader object - either for simple input or Geant4 input.
Bool_t CCMLEM::SetInputReader(void) {

  TString fullName = fInputName;

  TFile* file = new TFile(fullName, "READ");
  if (!file->IsOpen()) {
    cout << "##### Error in CCMLEM::SetInputReader!" << endl;
    cout << "Could not open requested file" << endl;
    return false;
  }

  if (file->Get("data")) {
    file->Close();
    fReader = new InputReaderSimple(fullName);
  } else if (file->Get("Setup") &&
             file->Get("Events")/*&&
             file->Get("G4SimulationData_Reconstruction")*/) {
    file->Close();
    fReader = new InputReaderGeant(fullName);
  } else if (file->Get("Pos&EnergyRecoClus") /*&&
             file->Get("Reco") && 
             file->Get("Real")*/) {
    file->Close();
    fReader = new InputReaderEI(fullName);
  } else {
    cout << "##### Error in CCMLEM::SetInputReader()!" << endl;
    cout << "Unknown data format" << endl;
    return false;
  }

  if (fVerbose) { fReader->Print(); }

  return true;
}
//------------------------------------------
///Function for image reconstruction and returns kTRUE when the last iteration of   
///the reconstructed image is done.Sets 2D image histogram and 
///1D histogram (graph) of non-pixelized reconstructed image.
/// Even it can produce 3D image histogram but commented for this moment.
///By calling cones from ComptonCone class, it uses the necessary input data from InputReader,
///In the loop, it reconstructs requested numbers of ComptonCones and uses the numerical method 
///to get intersection of ComptonCone with the image plane.
/// For each ComptonCone, the following processes will be done: 
///cases: (1) To fill TClonesArray by intersection points and global bin number,
///AddIsectionPoint() is used (in the following I will describe this function).
///(2) The global bin numbers are sorted, then for each touched pixel
///the length of track is calculated to estimate the reconstructed cone.
Bool_t CCMLEM::Reconstruct(void) {

  // image histogram
  /*fImage[0] = new TH3F("image", "image", fNbinsZ, -fDimZ / 2. , fDimZ / 2.,
                       fNbinsY, -fDimY / 2., fDimY / 2., fNbinsX, -fDimX/2., fDimX/2.);*/
  fImage[0] = new TH2F("image", "image", fNbinsZ, -fDimZ / 2. , fDimZ / 2.,
                       fNbinsY, -fDimY / 2., fDimY / 2.);
  fImage[0]->GetXaxis()->SetTitle("z [mm]");
  fImage[0]->GetYaxis()->SetTitle("y [mm]");
  //fImage[0]->GetZaxis()->SetTitle("x [mm]");
  
  //fAngDiff = new TH1F("Angle Difference", "Angle Difference", 90, 0, 90);
  //fAngDiff->GetXaxis()->SetTitle("Angle Diff.(deg)");
  //fAngDiff->GetYaxis()->SetTitle("events");
  TH1F *Energy = new TH1F ("Deposited Energy", " Deposited Energy",1000, 0, 10);
  Energy->GetXaxis()->SetTitle("Deposited_Energy");
  Energy->GetYaxis()->SetTitle("Counts");
  
//   TH1F *EnergyPri = new TH1F ("Primary Energy", " Primary Energy",1000, 0, 15);
//   EnergyPri->GetXaxis()->SetTitle("Primary_Energy");
//   EnergyPri->GetYaxis()->SetTitle("Counts");
  
  TH2* h2 = new TH2F("EScat. vs EAbs.", "EScat. vs EAbs.", 50, 0, 10, 50, 0, 10);
  h2->GetXaxis()->SetTitle("Energy_Absorber");
  h2->GetYaxis()->SetTitle("Energy_Scatterer");
  
  /*fGraph = new TGraph();
   fGraph->SetName("g");
   fGraph->SetTitle("non-pixelized reco image");
   fGraph->SetMarkerStyle(7);
   fGraph->SetMarkerColor(kRed);
   TH1F* htmp = new TH1F("g","non-pixelized reco
   image",fNbinsZ,-fDimZ/2.,fDimZ/2.);
   htmp->GetYaxis()->SetRangeUser(-fDimY/2.,fDimY/2.);
   htmp->GetXaxis()->SetTitle("z [mm]");
   htmp->GetYaxis()->SetTitle("y [mm]");
   fGraph->SetHistogram(htmp);
  */
   //Int_t ncells = fImage[0]->GetSize();
   //cout<<"total number of bins :"<< ncells<<endl;
  fPixelSizeZ = fDimZ / fNbinsZ;
  fPixelSizeY = fDimY / fNbinsY;
  fPixelSizeX = fDimX / fNbinsX;
  Double_t A = fDimY / 2.;
  Double_t B = fDimZ / 2.;
  Double_t F = fDimX / 2.;

  Int_t m, n, j, k;
  Double_t z, y, x;
  Double_t sum;
  TVector3 interactionPoint;
  TVector3 coneAxis;
  Double_t coneTheta;
  Double_t energy1, energy2, energy0, energy4, energyreco1, energyreco2;
  Double_t energy3 = 4.44;
  
  Double_t Absthick_z, Absthick_y, Absthick_x;
  TVector3 *point_e, *point_p, *point_dir, *point_abs_sc, *pointreco_e, *pointreco_p;
  TVector3 *Scatposition, *Absposition;
  point_abs_sc = new TVector3();
  fNIpoints = 0;
  fPoints = 0;
  IsectionPoint* tmppoint1;
  IsectionPoint* tmppoint2;

  const Double_t maxdist = sqrt(pow(fPixelSizeY, 2) + pow(fPixelSizeZ, 2));
  int identified;
  bool status;
  int counter = 0;
  int count_before = 0;
  int count = 0;
  int countback = 0;
  int countforward = 0;
  int count1 = 0;
  int count2 = 0;
  int count3 = 0;
  int count40 = 0;
  int count_other = 0;
  TStopwatch t;
  t.Start();
  
//   for (Int_t i = fStart; i < fStop/2; i++) {
//   //if (coneTheta <= M_PI/18.) {
//         
//         fSubFirst[i] = fEvent[i];
//         cout<< fSubFirst[i] << endl; 
//   }
//   for (Int_t i = fStop/2; i < fStop; i++) {
//         //} else {
//             
//             fSubSecond[i - fStop/2] = fEvent[i];
//             cout<< fSubSecond[i] << endl;
//         //}
//   }

/// This loop will always analyze fStop-fStart events starting with
/// event number fStart. If some of the events are not valid they will
/// be skipped, but still fStop-fStart events will be analyzed i.e.
/// last analyzed event will have number fStop+n, where n is number of
/// skipped events. If you want to change this - remove 'counter' variable.

  for (Int_t i = fStart; i < fStop; i++) {
      

      fNIpoints = 0;
      //count= 0;
      if (fVerbose)
      cout << "CCMLEM::Reconstruct(...) event " << i << endl << endl;

      status = fReader->LoadEvent(counter + i);

      if (status == false) {
         counter++;
         continue;
      }
      
      
      energy1 = fReader->GetEnergyLoss();
      energy2 = fReader->GetEnergyScattered();
      energy0 = fReader->GetEP();
      energy4 = fReader->GetEnergyPrimary();
      point_e = fReader->GetPositionScattering();
      //point_e->Print();
      point_p = fReader->GetPositionAbsorption();
      identified = fReader->GetIdentified();
      //point_p->Print();
      //point_dir = fReader->GetGammaDirScattered();
      //point_dir->Print();
      //point_abs_sc->SetXYZ(point_p->X()-point_e->X(),
      //point_p->Y()-point_e->Y(), point_p->Z()-point_e->Z());
      
      energyreco1 = fReader->GetEnergyLossReco();
      energyreco2 = fReader->GetEnergyScatteredReco();
      
      pointreco_e = fReader->GetPositionScatteringReco();
      pointreco_p = fReader->GetPositionAbsorptionReco();
      
      fScatthick_z = fReader->GetScatThickz();
      fScatthick_y = fReader->GetScatThicky();
      fScatthick_x = fReader->GetScatThickx();
      Absthick_z = fReader->GetAbsThickz();
      Absthick_y = fReader->GetAbsThicky();
      Absthick_x = fReader->GetAbsThickx();
      Scatposition = fReader->GetScattererPosition();
      Absposition = fReader->GetAbsorberPosition();
   /*   
      sum = energy1 + energy2;
      if (sum < 4.2 || sum > 4.6) {
          
          continue;
      }
      cout << sum << endl;*/
      
      energy2 = 4.4 - energy1;
/*      
      if (identified == 1 || identified == 3) {
          
      
      
       sum = energyreco1 + energyreco2;
       //sum = energy1 + energy2;
     
       if (sum >= 3 && sum < 4.6) {
           //count++;
           count1++;
         
            
       } else if (sum >= 3.8 && sum < 4.1) {
           
           count2++;
           
       } else if (sum >= 4.2 && sum < 4.6) {
           
           count3++;
         
       } else {
           
           continue;
       }
       
      //count++;
      //point_e->Print();
      //point_p->Print();
      
      Energy->Fill(sum);
      h2->Fill(energyreco2,energyreco1);
      //EnergyPri->Fill(energy0);
      //h2->Fill(energy2,energy1);
*/      
      
///To analyze Geant4 results (point source (Real info.)), set cuts for energies and positions. 
///Otherwise, comment them from line:242 to line:277. 
//       fA = 1, fB = 0, fC = 0;
//       fD = Absposition->Z()/2.;
//       
//       Double_t num = fabs(fA * point_dir->x() + fB * point_dir->y() + fC * point_dir->z() + fD);
//       Double_t denom = sqrt(fA * fA + fB * fB + fC * fC);
//       Double_t distance = num / denom;
//     
//       
//       point_dir->SetXYZ(point_e->x()+point_dir->x()*distance, point_e->y()+point_dir->y()*distance, point_e->z()+point_dir->z()*distance);
//       
//               energy0 = 4.40;
//               Double_t en = energy1 + energy2;
//               if (fabs(energy0 - en) > 0.5) {
// // //                  //cout << "energy_e : "<< energy1 << "\t" << "+" << "\t" << "energy_p :" << energy2 << "\t" << " = " << "\t" << "total energy : " << en << endl;
// // //                //cout<< "X of scatterer : " <<point1->X() <<endl;
// // //                //cout<< "----------------------------------------------------"<<endl;
//                   continue;
//                }
//     
//            if(point_e->X() >= (Scatposition->X() - (Scatthick_x/2)) && point_e->X() <=
//            (Scatposition->X() + (Scatthick_x/2))){
//               //cout<< "X of scatterer : " <<point1->X() <<endl;
//               if(point_p->X() >= (Absposition->X() - (Absthick_x/2)) && point_p->X() <=
//               (Absposition->X() + (Absthick_x/2))){
//                 //cout<< "X of absorber : " <<point2->X() <<endl;
//                   if(point_p->Y() >= (-Absthick_y) && point_p->Y() <= Absthick_y){
//                     //cout<< "Y of absorber : " <<point2->Y() <<endl;
//                   }
//                   if(point_p->Z() >= (-Absthick_z) && point_p->Z() <= Absthick_z){
//                      //cout<< "Z of absorber : " <<point2->Z() <<endl;
//                   }
//                  //count++;
//               } else {
//                    //cout<< "out of scatterer : " << point_e->X() << "\t"<< "out of absorber : " <<       point_p->X()<< "\t" << point_p->Y() << "\t" << point_p->Z() <<endl; 
//                    //cout<< "----------------------------------------------------"<<endl; 
//                    continue;
//               }
//   
//   //          cout<< "----------------------------------------------------"<<endl;
//   
//           } else {
//   //           cout<< "----------------------------------------------------"<<endl;
//                continue;
//    
//          }
//    cout<<"---------------------"<<endl;
    //fAngDiff->Fill(point_abs_sc->Angle(*point_dir)*TMath::RadToDeg());
   
    if (fSmear) {
          point_e->SetXYZ(SmearBox(point_e->X(), fResolutionX),
                         SmearGaus(point_e->Y(), fResolutionY),
                         SmearBox(point_e->Z(), fResolutionZ));
          point_p->SetXYZ(SmearBox(point_p->X(), fResolutionX),
                         SmearGaus(point_p->Y(), fResolutionY),
                         SmearBox(point_p->Z(), fResolutionZ));
      energy1 = SmearGaus(energy1, GetSigmaE(energy1));
      energy2 = SmearGaus(energy2, GetSigmaE(energy2));
    }
    
    //ComptonCone* cone = new ComptonCone(pointreco_e, pointreco_p, energyreco1+energyreco2, energyreco2);
    ComptonCone* cone =
        new ComptonCone(point_e, point_p, energy1 + energy2, energy2);    
    interactionPoint = cone->GetApex();
    coneAxis = cone->GetAxis();
    coneTheta = cone->GetAngle();
    //cout<< "theta :" << "\t"<< coneTheta << endl;
   

    
    

   
    Double_t K = cos(coneTheta);
    Double_t a, b, c1, c2, c3, c4, c, z1, z2;
    Double_t d, e, f1, f2, f3, f4, f, y1, y2;
    
///For 3D image, when No. of bins is odd along third component of volume,
///user should add a half value of pixel size to it.

    //x = -F + fPixelSizeX/2;
    x = 0;
    //for (m = 1; m < fNbinsX + 1; m++) {
        
        //fNIpoints = 0;
    
    y = -A;
      if(fVerbose) cout<<"Loop over horizontal lines..."<<endl;
    
      for(j=0; j< fNbinsY + 1; j++){
      
        a = 2*(pow(-coneAxis.Z(),2) - pow(K,2));
      
        b = 2*((-coneAxis.Z())*((interactionPoint.Z())*(-coneAxis.Z()) - y*((-coneAxis.Y())) + (interactionPoint.Y())*
	  (-coneAxis.Y()) + (interactionPoint.X())*(-coneAxis.X()) + x - (-coneAxis.X())*x) - (interactionPoint.Z())*pow(K,2));
      
        c1 = pow(interactionPoint.Y(),2) + pow(interactionPoint.X(),2) + ((-2)*pow(interactionPoint.Y(),2)*pow(-coneAxis.Z(),2)) +  
	   ((-2)*pow(interactionPoint.X(),2)*pow(-coneAxis.Z(),2)) + ((-2)*pow(interactionPoint.Y(),2)*pow(-coneAxis.Y(),2)) + 
	   ((-4)*(interactionPoint.Y())*(interactionPoint.X())*(-coneAxis.Y())*(-coneAxis.X())) + ((-2)*pow(interactionPoint.X(),2)*
	   pow(-coneAxis.X(),2));
      
        c2 = (-2)*(interactionPoint.X())*x + 4*(interactionPoint.X())*x*
	   pow(-coneAxis.Z(),2) - 4*(interactionPoint.Y())*(-coneAxis.Y())*x - 4*(interactionPoint.X())*
	   (-coneAxis.X())*x + 4*(interactionPoint.Y())*(-coneAxis.Y())*(-coneAxis.X())*x + 
	   4*(interactionPoint.X())*pow(-coneAxis.X(),2)*x - pow(x,2) - 2*pow(-coneAxis.Z(),2)*
	   pow(x,2) + 4*(-coneAxis.X())*pow(x,2) - 2*pow(-coneAxis.X(),2)*pow(x,2);
      
        c3 = pow(y,2)*(1 - 2*pow(-coneAxis.Z(),2) - 2*pow(-coneAxis.Y(),2)) + y*
	  ((interactionPoint.Y())*((-2) + 4*pow(-coneAxis.Z(),2) + 4*pow(-coneAxis.Y(),2)) + 4*(-coneAxis.Y())*
	  ((interactionPoint.X())*(-coneAxis.X()) +  x - (-coneAxis.X())*x));
      
        c4 = (pow(y,2) - 2*y*(interactionPoint.Y()) + pow(interactionPoint.Y(),2) + pow((interactionPoint.X() - x),2))*
	   ((-1) + 2*pow(K,2));
      
      
        c = (-2)*pow(K,2)*(c1 + c2 + c3 + c4); 
      
        if(c >= 0){
	  z1 = (b - sqrt(c))/a;
          z2 = (b + sqrt(c))/a;
	  if(fabs(z1)<1e-14) z1 = z1*1e6;
	  if(fabs(z2)<1e-14) z2 = z2*1e6;
	  AddIsectionPoint("hor", x, y, z1);
	  AddIsectionPoint("hor", x, y, z2);
	
        }
        y=y+fPixelSizeY;
      
      } //end of loop over horizontal lines

    z = -B; 
      if(fVerbose) cout<<"Loop over vertical lines..."<<endl;
      for(k=0; k< fNbinsZ + 1; k++){
      
        d = 2*(pow(-coneAxis.Y(),2) - pow(K,2));
      
        e = 2*((-coneAxis.Y())*((interactionPoint.Z())*(-coneAxis.Z()) - z*((-coneAxis.Z())) + (interactionPoint.Y())*
	  (-coneAxis.Y()) + (interactionPoint.X())*(-coneAxis.X()) + x - (-coneAxis.X())*x) - (interactionPoint.Y())*pow(K,2));
      
        f1 = pow(interactionPoint.Z(),2) + pow(interactionPoint.X(),2) + ((-2)*pow(interactionPoint.Z(),2)*pow(-coneAxis.Z(),2)) +  
	   ((-2)*pow(interactionPoint.Z(),2)*pow(-coneAxis.Y(),2)) + ((-2)*pow(interactionPoint.X(),2)*pow(-coneAxis.Y(),2)) + 
	   ((-4)*(interactionPoint.Z())*(interactionPoint.X())*(-coneAxis.Z())*(-coneAxis.X())) + ((-2)*pow(interactionPoint.X(),2)*
	   pow(-coneAxis.X(),2));
      
        f2 = (-2)*(interactionPoint.X())*x + 4*(interactionPoint.X())*x*
	   pow(-coneAxis.Y(),2) - 4*(interactionPoint.Z())*(-coneAxis.Z())*x - 4*(interactionPoint.X())*
	   (-coneAxis.X())*x + 4*(interactionPoint.Z())*(-coneAxis.Z())*(-coneAxis.X())*x + 
	   4*(interactionPoint.X())*pow(-coneAxis.X(),2)*x - pow(x,2) - 2*pow(-coneAxis.Y(),2)*
	   pow(x,2) + 4*(-coneAxis.X())*pow(x,2) - 2*pow(-coneAxis.X(),2)*pow(x,2);
      
        f3 = pow(z,2)*(1 - 2*pow(-coneAxis.Z(),2) - 2*pow(-coneAxis.Y(),2)) + z*
	   ((interactionPoint.Z())*((-2) + 4*pow(-coneAxis.Z(),2) + 4*pow(-coneAxis.Y(),2)) + 4*(-coneAxis.Z())*
	   ((interactionPoint.X())*(-coneAxis.X()) + x - (-coneAxis.X())*x));
      
        f4 = (pow(z,2) - 2*z*(interactionPoint.Z()) + pow(interactionPoint.Z(),2) + pow((interactionPoint.X() - x),2))*
	   ((-1) + 2*pow(K,2));
      
        f = (-2)*pow(K,2)*(f1 + f2 + f3 + f4);
      
        if(f >= 0){
	  y1 = (e - sqrt(f))/d;
	  y2 = (e + sqrt(f))/d;
	  if(fabs(y1)<1e-14) y1 = y1*1e6;
	  if(fabs(y2)<1e-14) y2 = y2*1e6;
	  AddIsectionPoint("ver",x , y1, z);
	  AddIsectionPoint("ver",x , y2, z);
	
        }      
        z=z+fPixelSizeZ;
      
      } //end of loop over vertical lines

    Int_t index[fNIpoints];
    Int_t fA[fNIpoints];
    IsectionPoint* tempp;

    for (Int_t i = 0; i < fNIpoints; i++) {
      tempp = (IsectionPoint*)fArray->At(i);
      fA[i] = tempp->GetBin();
      // if(fVerbose) tempp->Print();
    }

    TMath::Sort(fNIpoints, fA, index, kFALSE);

    TVector3* tmpvec1;
    TVector3* tmpvec2;
    Double_t dist;
    Int_t binno1, binno2, binno;
    SMElement* temp;
    

    for (int h = 0; h < fNIpoints; h = h + 2) {

      if (fVerbose)
        cout << " index[" << h << "]=" << index[h] << ", index[" << h + 1
             << "]=" << index[h + 1] << endl;
      tmppoint1 = (IsectionPoint*)fArray->At(index[h]);
      tmppoint2 = (IsectionPoint*)fArray->At(index[h + 1]);
      if (tmppoint1 == NULL || tmppoint2 == NULL) {
        if (fVerbose)
          cout << "Something went wrong"
               << "\t" << tmppoint1 << "\t" << tmppoint2 << endl;
         continue;
      }
      tmpvec1 = tmppoint1->GetPointCoordinates();
      binno1 = tmppoint1->GetBin();
      tmpvec2 = tmppoint2->GetPointCoordinates();
      binno2 = tmppoint2->GetBin();
      if (fVerbose)
        cout << " binno1=" << binno1 << ", binno2=" << binno2 << endl << endl;
      dist = ((*tmpvec1) - (*tmpvec2)).Mag();
      if (dist > maxdist) {
        // cout<<"Event "<<h<<": distance exceeds pixel
        // diagonal"<<dist/maxdist<<" times"<<endl;
        continue;
      }

      if (binno1 != binno2) {
         //count++;
        cout << "eventno :" << i << endl;
        cout << binno1 << "!=" << binno2
             << " ->Bin numbers are different when they should not!" << endl;
        cout<< "Position 1 : " << tmpvec1->X() << "," << "\t" << tmpvec1->Y() << "," << "\t" << tmpvec1->Z() <<endl; 
        cout<< "Position 2 : " << tmpvec2->X() << "," << "\t" << tmpvec2->Y() << "," << "\t" << tmpvec2->Z() <<endl; 
             continue;
      }

      fImage[0]->SetBinContent(binno1, fImage[0]->GetBinContent(binno1) + dist);
      
      temp = (SMElement*)fSM->ConstructedAt(fPoints++);
      temp->SetEvBinDist(i, binno1, dist);

      if (fVerbose) temp->Print();
    }
    //cout << "event = " << i << endl;
    //cout << "fPoints = " << fPoints << endl;
    //cout<< "no. of Ipoints : " << fNIpoints<<endl; 
   // x = x + fPixelSizeX;
  //}
    
    //cout<< " number of bins : " << count <<endl;
    delete cone;

    if (fVerbose)
      cout << "----------------------------------------------------------------"
           << endl;
           

   // } 
    
  } // end of loop over events
  
  //Int_t ncells = fImage[0]->GetSize();
  //cout<<"total number of bins :"<< ncells<<endl;
  //cout<< "no. of Points : " << fPoints<<endl;
  //SaveToFile(fAngDiff);
  //SaveToFile(EnergyPri);
  //SaveToFile(Energy);
  //SaveToFile(h2);
  
  fArray->Clear("C");
  
  //count= count1 + count2 + count3;
  
    //cout<< " number of events(3.4) : " << count1 <<endl;
    //cout<< " number of events(3.9) : " << count2 << endl;
    //cout<< " number of events(4.4) : " << count3 << endl;
    //cout<< " number of events : " << count <<endl;
//    cout<< " number of events(30 < x < 40) : " << count40 << endl;
//    cout<< " number of events(others) : " << count_other << endl;
//    
  SaveToFile(fImage[0]);
  
/// Sensitivity map calculation
/// with respect to our geometry in Geant4 simulation
  
  Int_t binz, biny, binx;
  Double_t pixelX;
  Double_t pixelY, pixelZ;
  Double_t pixelx, pixely, pixelz;
  Double_t VecSize;
  TVector3 pixelCenter;
  TVector3 DetCenter;
  TVector3 UnitDetCenter;
  TVector3 Vector;
  Double_t angle, costheta;
  Double_t factor,factor1,factor2;
  Int_t nbinZ = fScatthick_z;
  //cout<<fScatthick_z<<endl;
  Int_t nbinY = fScatthick_y;
  Int_t nbinX = fScatthick_x;
  fSensitivity = new TH2F("Sensitivity map", "Sensitivity map", fNbinsZ, -fDimZ / 2. , fDimZ / 2.,
                       fNbinsY, -fDimY / 2., fDimY / 2./*, fNbinsX, -fDimX / 2. , fDimX / 2.*/);
  fSensitivity->GetXaxis()->SetTitle("z [mm]");
  fSensitivity->GetYaxis()->SetTitle("y [mm]");
  //fSensitivity->GetZaxis()->SetTitle("x [mm]");
  TH3F *Scat = new TH3F("scatterer", "scatterer", 99, -49.5 , 49.5,
                       100, -50, 50, 13, 193.5, 206.5);
  /*TH2F *Scat = new TH2F("scatterer", "scatterer", nbinZ, -fScatthick_z/2. , fScatthick_z/2., nbinY, -fScatthick_y/2. , fScatthick_y/2./*, nbinX, Scatposition->X() - fScatthick_x/2. , Scatposition->X() + fScatthick_x/2.);*/  
  //for (Int_t k = 1; k < fNbinsX + 1; k++) {
    for (Int_t i = 0; i < fNbinsZ + 1; i++) {
     
      for (Int_t j = 0; j < fNbinsY + 1; j++) {
          factor = 0;
          factor1 = 0;
          factor2 = 0;
          //pixelX = fSensitivity->GetZaxis()->GetBinCenter(k);
          pixelZ = fSensitivity->GetXaxis()->GetBinCenter(i) + 0.5;
          pixelY = fSensitivity->GetYaxis()->GetBinCenter(j) + 0.5;
          
         // cout<< pixelZ << "," << pixelY << "," << pixelX << endl;
          pixelCenter.SetXYZ(pixelX, pixelY, pixelZ);
          //pixelCenter.Print();
          //double x = Scatposition->X() - fScatthick_x/2.;
          double x = 193.5;
          //while (x <= Scatposition->X() + fScatthick_x/2.) {  
          while (x <= 206.5) {     
          //for (int n = 1; n <= nbinZ; n++) {
          for (int n = 1; n <= 99; n++) {    
              //for ( int m = 1; m <= nbinY; m++) {
               pixelz = Scat->GetXaxis()->GetBinCenter(n);
               //pixely = Scat->GetYaxis()->GetBinCenter(m);
               //cout<< pixelx << "," << pixely << "," << pixelz << endl;
               DetCenter.SetXYZ(x, 0, pixelz);
               //DetCenter.Print();
               
               Vector.SetXYZ(pixelCenter.X() - DetCenter.X(), pixelCenter.Y() - DetCenter.Y(), pixelCenter.Z() - DetCenter.Z());
               //Vector.Print();
               UnitDetCenter = DetCenter.Unit();
               //UnitDetCenter.Print();
               VecSize = Vector.Mag();
               //angle = Vector.Angle(UnitDetCenter)*TMath::RadToDeg();
               costheta = ((Vector.X()*UnitDetCenter.X()) + (Vector.Y()*UnitDetCenter.Y()) + (Vector.Z()*UnitDetCenter.Z()))/VecSize;
               //costheta = TMath::Cos(angle);
               //cout << "angle :" << angle << " cos : " << costheta << endl;
               //cout<< pixelZ << "\t" << pixelY << "\t" << fabs(costheta)/(VecSize) << endl;
               //fSensitivity->Fill(pixelZ,pixelY,pixelX,fabs(costheta)/VecSize);
               factor1 = factor1 + fabs(costheta)/(VecSize*VecSize);
               //cout<< factor1 << " \t z" << endl;
            }
          //for (int m = 1; m <= nbinY; m++) {
          for (int m = 1; m <= 100; m++) {    
              
              pixely = Scat->GetYaxis()->GetBinCenter(m);
                  
              //cout<< pixelx << "," << pixely << "," << pixelz << endl;
              DetCenter.SetXYZ(x + 1, pixely, 0);
              //DetCenter.Print();
              
              Vector.SetXYZ(pixelCenter.X() - DetCenter.X(), pixelCenter.Y() - DetCenter.Y(), pixelCenter.Z() - DetCenter.Z());
              //Vector.Print();
              UnitDetCenter = DetCenter.Unit();
              //UnitDetCenter.Print();
              VecSize = Vector.Mag();
              //angle = Vector.Angle(UnitDetCenter)*TMath::RadToDeg();
              costheta = ((Vector.X()*UnitDetCenter.X()) + (Vector.Y()*UnitDetCenter.Y()) + (Vector.Z()*UnitDetCenter.Z()))/VecSize;
              //costheta = TMath::Cos(angle);
              //cout << "angle :" << angle << " cos : " << costheta << endl;
              //cout<< pixelZ << "\t" << pixelY << "\t" << fabs(costheta)/(VecSize) << endl;
              //fSensitivity->Fill(pixelZ,pixelY,pixelX,fabs(costheta)/VecSize);
              factor2 = factor2 + fabs(costheta)/(VecSize*VecSize);
          
              //cout<< factor2 << " \t y" << endl;
            }    
       
          
          x = x + 2;
        }
          factor = factor1 + factor2;
          //cout<< factor << endl;
          fSensitivity->Fill(pixelZ,pixelY,factor); 
          //fSensitivity->Fill(pixelZ,pixelY,pixelX,factor); 
      }
    //}
  }
  SaveToFile(fSensitivity);
  
/// End of Sensitivity map calculation
  
  fSigma [150] = 1.;
  for (int iter = 1; iter < fIter + 1; iter++) {
      
      Iterate(fStop, iter);
      
      
      if (fSigma[150] < 0.01) {
          t.Stop();
          t.Print();
          //cout<< "sigma value" << fSigma[150] << endl;
          //DrawCanvas();
          return 0;
      }

      

  }
  
  //DrawCanvas();
  //fSM->Clear("C");
  
  //DrawHisto();
  //GetSigmaError();
   //t.Stop();
   //t.Print();

  return kTRUE;
}

//------------------------------------
///Function to check direction of intersecting lines and get coordinate of intersection points.
///TClonesArray is filled by coordinate of intersection points and their global bin number.
///To get intersection points and global bin number, adding and subtracting 
///a bit value make us sure each point belongs to which two of pixels. 
/// So this function returns two points added for each pixel if all conditions are met.
///\param dir (TString) - direction of intersecting lines marched to 
///get intersection points on pixel borders
///\param x (Double_t) - x-component of intersection point on the image plane
///\param y (Double_t) - y-component of intersection point on the image plane
///\param z (Double_t) - z-component of intersection point on the image plane
Int_t CCMLEM::AddIsectionPoint(TString dir, Double_t x, Double_t y,
                               Double_t z) {
  if (fVerbose) cout << dir << "\t" << x << "\t" << y << "\t" << z << endl;
  dir.ToLower();
  if (dir != "hor" && dir != "ver") {
    // if(fVerbose) cout<<"Unknown direction of intrsecting line: "<<dir<<endl;
    return 0;
  }

  if (fabs(y) > fDimY / 2. || fabs(z) > fDimZ / 2./*|| fabs(x)>fDimX/2.*/) {
    // if(fVerbose) cout<<"point outside of image range..."<<endl;
    return 0;
  }

  IsectionPoint* tmppoint;
  // SMElement* temp;
  Int_t added = 0;
  Int_t pixelZ, pixelY, pixelX;
  Double_t yplus, yminus;
  Double_t zplus, zminus;
  if (dir == "hor") { // adding point from intersections with horizontal lines
    pixelZ = fImage[0]->GetXaxis()->FindBin(z);
    //pixelX = fImage[0]->GetZaxis()->FindBin(x);
    if (pixelZ > fNbinsZ) pixelZ = fNbinsZ; // inclusion of upper edges of histo
    // cout<<"pixelZ : " <<pixelZ<< endl;
    yplus = y + 0.005 * fPixelSizeY;
    yminus = y - 0.005 * fPixelSizeY;

    if (fabs(yplus) <= fDimY / 2) {
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);
      pixelY = fImage[0]->GetYaxis()->FindBin(yplus);
      //tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ, pixelY,pixelX), x, y, z);
      tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ, pixelY), x, y, z);
      // if(fVerbose) tmppoint->Print();
      // fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
    }
    if (fabs(yminus) <= fDimY / 2) {
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);
      pixelY = fImage[0]->GetYaxis()->FindBin(yminus);
      //tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ, pixelY,pixelX), x, y, z);
      tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ, pixelY), x, y, z);
      // if(fVerbose) tmppoint->Print();
      // fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
    }
  }
  if (dir == "ver") { // adding point from intersections with vertical lines
    pixelY = fImage[0]->GetYaxis()->FindBin(y);
    //pixelX = fImage[0]->GetZaxis()->FindBin(x);
    if (pixelY > fNbinsY) pixelY = fNbinsY; // inclusion of upper edges of histo
    // cout<<"pixelY : " <<pixelY<< endl;
    zplus = z + 0.005 * fPixelSizeZ;
    zminus = z - 0.005 * fPixelSizeZ;

    if (fabs(zplus) <= fDimZ / 2) {
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);

      pixelZ = fImage[0]->GetXaxis()->FindBin(zplus);
      //tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ, pixelY,pixelX), x, y, z);
      tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ, pixelY), x, y, z);
      // if(fVerbose) tmppoint->Print();
      // fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
    }
    if (fabs(zminus) <= fDimZ / 2) {
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);

      pixelZ = fImage[0]->GetXaxis()->FindBin(zminus);
      //tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ, pixelY,pixelX), x, y, z);
      tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ, pixelY), x, y, z);
      // if(fVerbose) tmppoint->Print();
      // fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
    }
  }
  // if(fVerbose)  cout<<added<<" points added..."<<endl<<endl;;
  return added;
}

//------------------------------------
/// Applies the iteration formula (multiplicative correction) 
///and update the image of reconstructed cones after each iteration.
///then returns kTRUE when the last reconstructed image is bulit.
///\param nstop (Int_t) - number of last event for image reconstruction
///\param iter (Int_t) - number of iteration for image reconstruction
Bool_t CCMLEM::Iterate(Int_t nstop, Int_t iter) {

  //gStyle->SetOptStat(1101);
  //gStyle->SetOptFit(1011);
    
  int lastiter = iter - 1;
  
  double sigma[lastiter + 1];
  //TH1D* ProZ[150];
  //TH1D* ProY[150];
  //TH1D* ProX[150];
  
///For the unknown source distribution, user should define an appropriate fitting
///function to compare reasonablely sigma value of the new 10th iteration with previous one.
/*  
  TF1* func_z = new TF1("func_z", "[2]*TMath::Erfc((x-[0])/sqrt(2.)/[1])",-10.5, 20.5);
          func_z->SetParName(0, "X_{0}");
          func_z->SetParName(1, "#sigma");
          func_z->SetParName(2, "A");
*/          
/*          
  TF1* func_z = new TF1(
      "func", "[0]/2*TMath::Erf((x-[1])/[3])-[0]/2*TMath::Erf((x-[2])/[3])",
      -30.5, 30.5);
   func_z->SetParameter(0, 2000);
   func_z->SetParameter(1, 0);
   func_z->SetParameter(2, 0);
   //func_z->SetParameter(2, 200);
   func_z->SetParameter(3, 5);
  func_z->SetParNames("Constant_L","Mean_value_L","Mean_value_R","Sigma_z");
  */
/*  
   TF1* func_y = new TF1("FitProjection", FitProjection, -10.5, 10.5, 3);
   func_y->SetParameters(5,-5,5);
   func_y->SetParNames("Constant_y","Mean_value_y","Sigma_y");
   */
/*
  TF1 *func_z = new TF1("func_z","[0]/2*TMath::Erf((x-[1])/[3])-[0]/2*TMath::Erf((x-[2])/[3])",-30, 30);
  func_z->SetParameter(0, 700);
  func_z->SetParameter(1, -20);
  func_z->SetParameter(2, 20);
  func_z->SetParameter(3, 5);
  func_z->SetParNames("Constant","Mean_value_L","Mean_value_R","Sigma_z");
  
  TF1 *func = new TF1("func","[0]/2*TMath::Erf((x-[1])/[3])-[0]/2*TMath::Erf((x-[2])/[3])",-30, 30);
  func->SetParameter(0, 70000);
  func->SetParameter(1, -20);
  func->SetParameter(2, 20);
  func->SetParameter(3, 5);
  func->SetParNames("Constant","Mean_value_L","Mean_value_R","Sigma_z");
 */  

   TF1* func_z = new TF1("FitProjection", FitProjection, -12.5, 12.5, 3);
   func_z->SetParameters(fP0,fP1,fP2);
   func_z->SetParNames("Constant","Mean_value","Sigma_z");
   
   TF1* func = new TF1("FitProjection", FitProjection, -12.5, 12.5, 3);
   func->SetParameters(fP0,fP1,fP2);
   func->SetParNames("Constant","Mean_value","Sigma_z");
/*   
   TF1* func = new TF1("func","[0]*TMath::Exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))", -14.5, 14.5);
   //func->SetParameters(fP0,fP1,fP2);
   func->SetParameter(0, 400000);
   func->SetParameter(1, 0.05);
   func->SetParameter(2, 2);
   func->SetParNames("Constant","Mean_value","Sigma_z");
   */
//   func->SetParameter(0, 18000);
//   func->SetParameter(1, 0.01);
//   func->SetParameter(2, 0.1);
  // func->SetParameter(3,8.1);

  if (fImage[lastiter] == NULL) {
    cout << "Error in CCMELM::Iterate(). Last iteration NULL" << endl;
    return kFALSE;
  }

//   TH3F* hlastiter = (TH3F*)fImage[lastiter];
//   fImage[lastiter + 1] = (TH3F*)hlastiter->Clone();
//   TH3F* hthisiter = fImage[lastiter + 1];
  
  TH2F* hlastiter = (TH2F*)fImage[lastiter];
  fImage[lastiter + 1] = (TH2F*)hlastiter->Clone();
  TH2F* hthisiter = fImage[lastiter + 1];
  
  //TH2F* Sensitivity = (TH2F*)hthisiter->Clone(Form("Sensitivity_%i", lastiter  ));
  //Sensitivity->SetTitle(Form("Sensitivity_iter%i", lastiter  ));
  
  hthisiter->Reset();
  hthisiter->SetName(Form("%s_iter%i", fImage[0]->GetName(), lastiter + 1));
  hthisiter->SetTitle(Form("%s_iter%i", fImage[0]->GetTitle(), lastiter + 1));
  Int_t eventno;
  Int_t eventno_prev = 0;
  Int_t entry;
  Int_t binno;
  //Int_t nz = fNbinsZ + 2.;
  //Int_t ny = fNbinsY + 2.;
  Double_t dist, addvalue;
  SMElement* temp;
  fDenominator[nstop + 1];
  Int_t binz, biny, binx;
  Double_t pixelX = 0;
  Double_t pixelY, pixelZ;
  Double_t VecSize;
  TVector3 pixelCenter;
  TVector3 DetCenter;
  TVector3 UnitDetCenter;
  TVector3 connectingVector;
  Double_t angle, costheta;
  Double_t factor;
  Int_t count = 0;
  for (int i = 0; i < nstop + 1; i++) {
    
    fDenominator[i] = 0;
    //fSensitivity[i] = 0;
  }
  Int_t nSMentries = fSM->GetEntries();
  //cout<<"nSMentries = "<<nSMentries <<endl;
  
  for (entry = 0; entry < nSMentries; entry++) {
    temp = (SMElement*)fSM->At(entry);
    binno = temp->GetBin();
    eventno = temp->GetEvent();
    dist = temp->GetDist();
    fDenominator[eventno] =
        fDenominator[eventno] + dist * hlastiter->GetBinContent(binno);
        //cout << " event : " << eventno << endl;
        //cout << " fac : " << dist * hlastiter->GetBinContent(binno) << endl;
        //cout << " denom : " << fDenominator[eventno] << endl;
        
  }
  

  //cout<< "--------------------------------------" << endl; 
  for (entry = 0; entry < nSMentries; entry++) {
    addvalue = 0;
    temp = (SMElement*)fSM->At(entry);
    binno = temp->GetBin();
    //binz = binno%(fNbinsZ + 2);
    //biny = ((binno-binz)/(fNbinsZ + 2))%(fNbinsY + 2);
    eventno = temp->GetEvent();
    dist = temp->GetDist();
    //cout<< " sen 123 : " << fSensitivity[binno] << endl;
    
    addvalue = dist  * hlastiter->GetBinContent(binno) / fDenominator[eventno];
    //cout<< "bin : "<< binno << " binz : " << binz << " biny : " << biny << " sen123 : "<< fSensitivity->GetBinContent(binz,biny) << endl;
    hthisiter->SetBinContent(binno, hthisiter->GetBinContent(binno) + addvalue);
  }
  
  //Sensitivity->Divide(fSensitivity);
  
    
  fSenHisto[lastiter + 1] = (TH2F*)hthisiter->Clone(Form("%s_SenHisto_iter%i", fImage[0]->GetName(), lastiter + 1 ));
  fSenHisto[lastiter + 1]->SetTitle(Form("%s_SenHisto_iter%i", fImage[0]->GetTitle(), lastiter + 1 ));
  
  //int nx = hlastiter->GetZaxis()->GetNbins();
  
  int nz = hlastiter->GetXaxis()->GetNbins();

  int ny = hlastiter->GetYaxis()->GetNbins();
  //cout<< "nx : " << nz << "ny : " << ny <<endl;
  //for (int k = 1; k < nx + 1; k++) {
      
    for (int i = 0; i < nz + 1; i++) {
      
      for (int j = 0; j < ny + 1; j++) {
          
          double c1 = fSensitivity->GetBinContent(i,j);
          double c2 = hthisiter->GetBinContent(i,j);
          if (c1 != 0 && c2 != 0) {
          fSenHisto[lastiter + 1]->SetBinContent(i,j,c2/c1);
          //cout<< i << j << c2/c1 << endl;
          //count++;
          }
      }
    //}
  }
  //cout << count << endl;
  //TH1D* ProX[150];
  TH1D* ProZ[150];
  TH1D* ProY[150];

  for (int i = 10; i <= lastiter ; i = i + 10) {
    
    //fProX[i] = fImage[i]->ProjectionZ();
    
    //ProX[i] = fSenHisto[i]->ProjectionZ();
    
    fProY[i] = fImage[i]->ProjectionY();
    //ProY[i]->Fit(func_y,"r");
    ProY[i] = fSenHisto[i]->ProjectionY();
    
    
    ProZ[i] = fSenHisto[i]->ProjectionX();
    ProZ[i]->Fit(func, "r");
    
    fProZ[i] = fImage[i]->ProjectionX();
    fProZ[i]->Fit(func_z, "r");
    
    sigma[i] = func_z->GetParameter(2);
    // cout<< "i : \t" << i << "\t" << "sigma : \t " << sigma[i]<< endl;
  }
///Comparing the sigma_z of new 10th iteration with previous one
///to check if it continues iterating more or not. 
  for (int j = 20; j <= lastiter ; j = j + 10) {
       //cout << fSigma[100] << endl;
      if (fabs(sigma[j]) !=0) fSigma[150] = (fabs(sigma[j] - sigma[j - 10]))/fabs(sigma[j]);
       //cout << fSigma[100] << endl;
      if (fSigma[150] < 0.01) {
          
          cout << "SigmaIter" << j << " - " << "SigmaIter" << j - 10 << " = " << "SigmaPercentageErr : " << fSigma[150] << endl;
 
          TCanvas* can = new TCanvas("can", "MLEM2D", 1000, 1000);
          TCanvas* canz = new TCanvas("MLEM1DZ","MLEM1DZ",1000,1000);
          TCanvas* cany = new TCanvas("MLEM1DY","MLEM1DY",1000,1000);
          //TCanvas* canx = new TCanvas("MLEM1DX","MLEM1DX",1000,1000);
          can->Divide(2, 2);
          can->cd(1);
          fImage[j - 10]->Draw("colz");
          can->cd(2);
          fSenHisto[j - 10]->Draw("colz");
          can->cd(3);
          fImage[j]->Draw("colz");
          can->cd(4);
          fSenHisto[j]->Draw("colz");
          
          canz->Divide(2, 2);
          canz->cd(1);
          fProZ[j - 10]->Draw();
          canz->cd(2);
          ProZ[j - 10]->Draw();
          canz->cd(3);
          fProZ[j ]->Draw();
          canz->cd(4);
          ProZ[j]->Draw();
          
          cany->Divide(2, 2);
          cany->cd(1);
          fProY[j - 10]->Draw();
          cany->cd(2);
          ProY[j - 10]->Draw();
          cany->cd(3);
          fProY[j]->Draw();
          cany->cd(4);
          ProY[j]->Draw();
        /*  
          canx->Divide(2, 2);
          canx->cd(1);
          fProX[j - 10]->Draw();
          canx->cd(2);
          ProX[j - 10]->Draw();
          canx->cd(3);
          fProX[j]->Draw();
          canx->cd(4);
          ProX[j]->Draw();
          */
          
          SaveToFile(can);
          SaveToFile(canz);
          SaveToFile(cany);
          //SaveToFile(canx);
 
          return 0;
      }
  }
            
           
  SaveToFile(hthisiter);
  //SaveToFile(Sensitivity);
  SaveToFile(fSenHisto[lastiter + 1]);
  
  return kTRUE;
}
//------------------------------------
Bool_t CCMLEM::DrawCanvas(void){
  
  //int iter;
/*  
  TF1 *func = new TF1("func","[0]/2*TMath::Erf((x-[1])/[4])-[2]/2*TMath::Erf((x-[3])/[4])",-100,200);
  func->SetParameter(0, 650000);
  func->SetParameter(1, -40);
  func->SetParameter(3, 40);
  func->SetParameter(2, 450000);
  func->SetParameter(4, 100);
  func->SetParNames("Constant_L","Mean_value_L","Constant_R","Mean_value_R","Sigma_z");
  
 /* 
  TF1 *func_z = new TF1("func_z","[0]/2*TMath::Erf((x-[1])/[4])-[2]/2*TMath::Erf((x-[3])/[4])",-80, 100);
  func_z->SetParameter(0, 700);
  func_z->SetParameter(1, -20);
  func_z->SetParameter(3, 10);
  func_z->SetParameter(2, 100);
  func_z->SetParameter(4, 5);
  func_z->SetParNames("Constant_L","Mean_value_L","Constant_R","Mean_value_R","Sigma_z");
 
 
  TF1* func_z = new TF1("func_z", "[0]*TMath::Erfc((x-[1])/sqrt(2.)/[2])",-80, 100);
  
  func_z->SetParameter(0, 800);
  func_z->SetParameter(1, 20);
  func_z->SetParameter(2, 5);
  
  func_z->SetParName(0, "A");
  func_z->SetParName(1, "X_{0}");
  func_z->SetParName(2, "#sigma");
  */
  TF1* func_z = new TF1("FitProjection", FitProjection, -14.5, 14.5, 3);
  //func_z->SetParameters(fP0,fP1,fP2);
  func_z->SetParameter(0, 18000);
  func_z->SetParameter(1, 0.1);
  func_z->SetParameter(2, 5);
  func_z->SetParNames("Constant","Mean_value","Sigma_z");
  
  
  TF1* func = new TF1("FitProjection", FitProjection, -14.5, 14.5, 3);
  //func_z->SetParameters(fP0,fP1,fP2);
  func->SetParameter(0, 1800000);
  func->SetParameter(1, 0.1);
  func->SetParameter(2, 10);
  func->SetParNames("Constant","Mean_value","Sigma_z");
  
  TH1D* ProZ[150];
  TH1D* ProY[150];
  TCanvas* can  = new TCanvas("MLEM2D","MLEM2D",1000,1000);
  TCanvas* canz = new TCanvas("MLEM1DZ","MLEM1DZ",1000,1000);
  TCanvas* cany = new TCanvas("MLEM1DY","MLEM1DY",1000,1000);
  can->Divide(4, 5);
  canz->Divide(4, 5);
  cany->Divide(4, 5);
  
  for(int i = 1; i <= fIter; i++) {
//       if (i = 0) {
//           continue;
//       }
    /*  if ( i = 1) {
         can->cd(1);
          fImage[i]->Draw("colz");
          can->cd(2);
          fSenHisto[i]->Draw("colz");
          
          canz->cd(1);
          fProZ[i]=fImage[i]->ProjectionX();
          fProZ[i]->Fit(func_z, "r");
          fProZ[i]->Draw();
          canz->cd(2);
          ProZ[i]=fSenHisto[i]->ProjectionX();
          ProZ[i]->Fit(func_z, "r");
          ProZ[i]->Draw();
          
          cany->cd(1);
          fProY[i]=fImage[i]->ProjectionY();
          fProY[i]->Draw();
          cany->cd(2);
          ProY[i]=fSenHisto[i]->ProjectionY();
          ProY[i]->Draw(); 
      }*/ 
      //if ( i > 0) {
          can->cd(2*i - 1);
          fImage[i]->Draw("colz");
          can->cd(2*i);
          fSenHisto[i]->Draw("colz");
          
          canz->cd(2*i - 1);
          fProZ[i]=fImage[i]->ProjectionX();
          fProZ[i]->Fit(func_z, "r");
          fProZ[i]->Draw();
          canz->cd(2*i);
          ProZ[i]=fSenHisto[i]->ProjectionX();
          ProZ[i]->Fit(func, "r");
          ProZ[i]->Draw();
          
          cany->cd(2*i - 1);
          fProY[i]=fImage[i]->ProjectionY();
          fProY[i]->Draw();
          cany->cd(2*i);
          ProY[i]=fSenHisto[i]->ProjectionY();
          ProY[i]->Draw();
    // }
      
  
  }
/*  
  TCanvas* can  = new TCanvas("MLEM2D","MLEM2D",1000,1000);
  TCanvas* canz = new TCanvas("MLEM1DZ","MLEM1DZ",1000,1000);
  TCanvas* cany = new TCanvas("MLEM1DY","MLEM1DY",1000,1000);
  can->Divide(4, 3);
  canz->Divide(4, 3);
  cany->Divide(4, 3);
  
  for(iter=0; iter<fIter+1; iter=iter + 2) {
      
    //cout<< " number of events :" <<endl;
      if (iter < 2) {
          
          //cout<< " number = 0 : " <<endl; 
          can->cd(1);
    
          fImage[0]->Draw("colz");
    
          canz->cd(1);
          fProZ[0]=fImage[0]->ProjectionX();
          fProZ[0]->Fit(func_z, "r");
          fProZ[0]->Draw();
          
          cany->cd(1);
          fProY[0]=fImage[0]->ProjectionY();
          fProY[0]->Draw();
      
      } 
      
      if (iter >= 2 && iter <= 20) {
          
        //cout<< " number > 0 : " <<endl;
          can->cd(iter-(iter/2.)+1);
    //gPad->SetLogz(1);
          fImage[iter]->Draw("colz");
    //fImage[iter]->SetMinimum(50);
          canz->cd(iter-(iter/2.)+1);
          fProZ[iter]=fImage[iter]->ProjectionX();
          fProZ[iter]->Fit(func_z, "r");
          fProZ[iter]->Draw();
          cany->cd(iter-(iter/2.)+1);
          fProY[iter]=fImage[iter]->ProjectionY();
          fProY[iter]->Draw();
      }
/*      
      if (iter == fIter) {
          
        //cout<< " number > 0 : " <<endl;
          can->cd(12);
    //gPad->SetLogz(1);
          fImage[fIter]->Draw("colz");
    //fImage[iter]->SetMinimum(50);
          canz->cd(12);
          fProZ[fIter]=fImage[fIter]->ProjectionX();
          fProZ[fIter]->Fit(func_z, "r");
          fProZ[fIter]->Draw();
          cany->cd(12);
          fProY[fIter]=fImage[fIter]->ProjectionY();
          fProY[fIter]->Draw();
      }
     
    //cout<< " end " <<endl;
  }
 */ 
  SaveToFile(can);
  SaveToFile(canz);
  SaveToFile(cany);
  return kTRUE;
}
//------------------------------------
/// Reads configuration file and sets values of private class
/// members according to read information.
///\param path (TString) - full path to the configuration file.
Bool_t CCMLEM::ReadConfig(TString path) {

  ifstream config(path);

  if (!(config.is_open())) {
    cout << "##### Could not open configuration file!" << endl;
    return kFALSE;
  }

  cout << "\n\nIn CCMLEM::Config(). Reading config file..." << endl;

  TString comment;

  while (!config.eof()) {
    comment.ReadLine(config);
    if (comment.Contains("Name of the input file")) {
      config >> fInputName;
      if (!fInputName.Contains(".root")) {
        cout << "##### Error in CCMLEM::Config()! Unknown file type!" << endl;
        return false;
      }
    }
    /*else if(comment.Contains("Center of reco plane")){
      config >> fXofRecoPlane >> fYofRecoPlane >> fZofRecoPlane;
    }*/
    else if (comment.Contains("Size of image volume")) {
      config >> fDimZ >> fDimY >> fDimX;
      //config >> fDimZ >> fDimY;
      //if (fDimZ < 1 || fDimY < 1 || fDimX<1) {
      if (fDimZ < 1 || fDimY < 1) {    
        cout << "##### Error in CCMLEM::Config()! Image size incorrect!"
             << endl;
        return false;
      }
    } else if (comment.Contains("No. of bins")) {
      config >> fNbinsZ >> fNbinsY >> fNbinsX;
      //config >> fNbinsZ >> fNbinsY;
      //if (fNbinsZ < 1 || fNbinsY < 1 || fNbinsX<1) {
      if (fNbinsZ < 1 || fNbinsY < 1) {
        cout << "##### Error in CCMLEM::Config()! Number of bins incorrect!"
             << endl;
        return false;
      }
    } else if (comment.Contains("Smear")) {
      config >> fSmear;
    } else if (comment.Contains("Position resolution")) {
      config >> fResolutionX >> fResolutionY >> fResolutionZ;
    } else if (comment.Contains("Fitting parameters")) {
      config >> fP0 >> fP1 >> fP2;
    } else if (comment.Contains("No. of MLEM iterations")) {
      config >> fIter;
      if (fIter < 0) {
        cout << "##### Error in CCMLEM::Config()! Number of iterations "
                "incorrect!"
             << endl;
        return false;
      }
    } else if (comment.Contains("Fresh output")) {
      config >> fFreshOutput;
    } else if (comment.Contains("No. of first and last event")) {
      config >> fStart >> fStop;
      if (fStart < 0 || fStop < 0 || fStop < fStart) {
        cout << "##### Error in CCMLEM::Config()! Number of first or last "
                "event incorrect!"
             << endl;
        return false;
      }
    } else if (comment.Contains("Verbose flag")) {
      config >> fVerbose;
    } else {
      cout << "##### Warning in CCMLEM::Config()! Unknown syntax!" << endl;
      cout << comment << endl;
    }
  }

  if (fVerbose) Print();

  config.close();

  return true;
}
//------------------------------------
///Saves 2D histogram and the projections of the final reconstructed image in the ROOT file. 
Bool_t CCMLEM::DrawHisto(void) {

  int lastiter = fIter;

  TH1D* hProZ[150];
  TH1D* hProY[150];
  //TH1D* hProX[150];
  TCanvas* can = new TCanvas("MLEM2D", "MLEM2D", 1000, 1000);
  
  can->Divide(2, 2);
  can->cd(1);
  // gPad->SetLogz(1);
  
  fImage[lastiter]->Draw("colz");
  // fImage[iter]->SetMinimum(200);
  hProZ[lastiter] = fImage[lastiter]->ProjectionX();
  hProY[lastiter] = fImage[lastiter]->ProjectionY();
  //hProX[lastiter] = fImage[lastiter]->ProjectionZ();
  can->cd(2);
  hProZ[lastiter]->Draw();
  can->cd(3);
  hProY[lastiter]->Draw();
  //can->cd(4);
  //hProX[lastiter]->Draw();
  

  SaveToFile(can);
 
  return kTRUE;
}
//------------------------------------
///Comparing the sigma_z of new 10th iteration with previous one
///to check if iteration continues more or not.
///user should define an appropriate fitting function
///with respect to source distribution.
Bool_t CCMLEM::GetSigmaError(void) {

  double sigma[fIter + 1];
  TH1D* ProZ[150];
  TH1D* ProY[150];

  TF1* func = new TF1(
      "func", "[0]/2*TMath::Erf((x-[1])/[3])-[0]/2*TMath::Erf((x-[2])/[3])",
      -40.5, 40.5);
  func->SetParameter(0, 2500);
  func->SetParameter(1, 0.01);
  func->SetParameter(2, 0.01);
  // func->SetParameter(3,200);
  func->SetParameter(3, 8.1);
  // func->SetParameter(5,8);
  for (int i = 10; i <= fIter + 1; i = i + 10) {
    
    //ProX[i] = fImage[i]->ProjectionZ();  
      ProY[i] = fImage[i]->ProjectionY();
      ProZ[i] = fImage[i]->ProjectionX();
      ProZ[i]->Fit(func, "r");
    //ProY[i]->Fit(func_y,"r");
      sigma[i] = func->GetParameter(2);
    // cout<< "i : \t" << i << "\t" << "sigma : \t " << sigma[i]<< endl;
  }
  
  for (int j = 20; j <= fIter + 1; j = j + 10) {
       //cout << fSigma[100] << endl;
      if (fabs(sigma[j]) !=0) fSigma[150] = (fabs(sigma[j] - sigma[j - 10]))/fabs(sigma[j]);
       //cout << fSigma[100] << endl;
      if (fSigma[150] < 0.01) {
 
          cout << "SigmaIter" << j << " - " << "SigmaIter" << j - 10 << " = " << "SigmaPercentageErr : " << fSigma[150] << endl;
 
          TCanvas* can = new TCanvas("can", "MLEM2D", 1000, 1000);
          can->Divide(2, 2);
          can->cd(1);
          fImage[j - 10]->Draw("colz");
          can->cd(2);
          ProZ[j - 10]->Draw();
          can->cd(3);
          ProY[j - 10]->Draw();
          //can->cd(4);
          //ProX[j - 10]->Draw();
                
          SaveToFile(can);
          //can->Update();
          //delete can;
 
          return 0;
            }
     }
 

  return kTRUE;
}
//-------------------------------------
///Gausian function to smear position resolution along y axis.
///returns a double value from gausian distribution with the given mean and sigma.
///\param val (double) - the given position value.
///\param sigma (double) - the given sigma value.
Double_t CCMLEM::SmearGaus(double val, double sigma) {
  return gRandom->Gaus(val, sigma);
}
//------------------------------------
///Returns a double value from Uniform function with respect to position resolution value.
///\param x (double) - the given position value.
Double_t CCMLEM::SmearBox(double x, double resolution) {
  return gRandom->Uniform(x - (resolution / 2), x + (resolution / 2));
}
//------------------------------------
///Returns the sigma value from the fitting function of deposited energy plot.
///\param energy (double) - the given energy value.
Double_t CCMLEM::GetSigmaE(double energy) {
  double sigma = fHisto->GetFunction("fit1")->Eval(energy) * energy;
  return sigma;
}
//-------------------------------------------
///Returns a double value from fitting function(gausian) used 
///for comparing sigma_z values of different iterations.
///\param x (Double_t*) - unknown value.
///\param p (Double_t*) - fitting parameter.
Double_t CCMLEM::FitProjection(Double_t *x,Double_t *p) {
      Double_t a = 0;
      if (p[2]!=0) a = (x[0] - p[1])/p[2];
      Double_t fitvalue = p[0]*TMath::Exp(-0.5*a*a);
      return fitvalue;
}
//-------------------------------------------
///Returns a double value from fitting function(Breit-Wiegner) used 
///for comparing sigma_z values of different iterations.
///\param x (Double_t*) - unknown value.
///\param p (Double_t*) - fitting parameter.
Double_t CCMLEM::FitProjection1(Double_t *x,Double_t *p) {
      Double_t a = 0;
      a = (x[0] - p[1])*(x[0] - p[1]) + 0.5*0.5*p[2]*p[2];
      Double_t fitvalue = p[0]*0.5*p[2]/a;
      return fitvalue;
}
//-------------------------------------------
///Prints details of the CCMLEM class obejct.
void CCMLEM::Print(void) {
  cout << "\nCCMLEM::Print()" << endl;
  cout << setw(35) << "Name of input file: \t"
       << fInputName << endl;
   //cout << setw(35) << "Center of reco plane: \t" << fXofRecoPlane << ", "
   //<< fYofRecoPlane << ", " << fZofRecoPlane << endl;
  cout << setw(35) << "Size of image volume: \t" << fDimZ << ", " << fDimY
       << ", " << fDimX << endl;
  cout << setw(35) << "No. of bins: \t" << fNbinsZ << ", " << fNbinsY << ", "
       << fNbinsX << endl;
  cout << setw(35) << "Smear level: \t" << fSmear << endl;
  cout << setw(35) << "Pos resolution: \t" << fResolutionX << ", "
       << fResolutionY << ", " << fResolutionZ << endl;
  cout << setw(35) << "No. of MLEM iterations: \t" << fIter << endl;
  cout << setw(35) << "FreshOutput level: \t" << fFreshOutput << endl;
  cout << setw(35) << "No. of first and last event: \t" << fStart << ", "
       << fStop << endl;
  cout << setw(35) << "Verbose level: \t" << fVerbose << endl << endl;
}
//--------------------------------------
///Sets default values of the private class members.
void CCMLEM::Clear(void) {
  fInputName = "dummy";
  //fSigma[100] = 1;
  fEvent[1000000] = 1;
  fSubFirst[1000000] = 1;
  fSubSecond[1000000] = 1;
  fXofRecoPlane = -1000;
  fYofRecoPlane = -1000;
  fZofRecoPlane = -1000;
  fDimZ = -1000;
  fDimY = -1000;
  fDimX = -1000;
  fNbinsZ = -1000;
  fNbinsY = -1000;
  fNbinsX = -1000;
  fResolutionX = -1000;
  fResolutionY = -1000;
  fResolutionZ = -1000;
  fP0 = -1000;
  fP1 = -1000;
  fP2 = -1000;
  fIter = -1000;
  fStart = -1000;
  fStop = -1000;
  fSmear = kFALSE;
  fFreshOutput = kFALSE;
  fVerbose = kFALSE;
  fNIpoints = -1000;
  fPoints = -1000;
  fPixelSizeZ = -1000;
  fPixelSizeY = -1000;
  fPixelSizeX = -1000;
  fReader = NULL;
  fArray = NULL;
  fSM = NULL;
  fOutputFile = NULL;
  //fGraph = NULL;
}
//--------------------------------------
///Saves object in the output file.
///\param ob (TObject*) - saved object.
Bool_t CCMLEM::SaveToFile(TObject* ob) {
  fOutputFile->cd();
  ob->Write();
  cout << ob->ClassName() << " " << ob->GetName() << " saved in the file "
       << fOutputFile->GetName() << endl;
  return kTRUE;
}
//--------------------------------------
