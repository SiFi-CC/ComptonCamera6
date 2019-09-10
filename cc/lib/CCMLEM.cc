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

  fArray = new TClonesArray("IsectionPoint", 1000);
  fSM = new TClonesArray("SMElement", 1000000);
  fNIpoints = 0;
  fPoints = 0;
/// To smear camera performance, this file is called for energy resolution.
/// It shows a function of energy deposited in a 10-cm-long LuAG(Ce) fiber 
///  with a square cross-section from Geant4 simulation   
  if (fSmear) {
    TString path = gSystem->Getenv("CC6DIR");
    TString name = path+"/share/ComptonCamera6/mlem_reco/EnergyResolutionExample.root";
     //TString name = "../work/results/Fiber_1_0.root";
    // TString name = "../work/results/Fiber_1_1.root";
    TFile* file = new TFile(name, "READ");
    fHisto = (TH1D*)file->Get("Scintillator_0ResolutionCombiEnergy");
    // fHisto = (TH1D*)file->Get("FiberEnergyUnc_Type3");
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
  } else if ((file->Get("G4SimulationData_Real") &&
             file->Get("G4SimulationData_Setup")&&
             file->Get("G4SimulationData_Reconstruction")) ||	// structure of simulation output before summer 2019
			 (file->Get("Setup") && file->Get("Events")) ) {	// structure of simulation output after summer 2019
    file->Close();
    fReader = new InputReaderGeant(fullName);
    dynamic_cast<InputReaderGeant*>(fReader)->SetFilter(fGeantFilter);
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
  fImage[0] = new TH2F("image", "image", fNbinsZ, -fDimZ / 2., fDimZ / 2.,
                       fNbinsY, -fDimY / 2., fDimY / 2.);
  fImage[0]->GetXaxis()->SetTitle("z [mm]");
  fImage[0]->GetYaxis()->SetTitle("y [mm]");
  //fImage[0]->GetZaxis()->SetTitle("x [mm]");

  //fAngDiff = new TH1F("Angle Difference", "Angle Difference", 90, 0, 90);
  //fAngDiff->GetXaxis()->SetTitle("Angle Diff.(deg)");
  //fAngDiff->GetYaxis()->SetTitle("events");

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
  // Int_t ncells = fImage[0]->GetSize();
  // cout<<"total number of bins :"<< ncells<<endl;
  fPixelSizeZ = fDimZ / fNbinsZ;
  fPixelSizeY = fDimY / fNbinsY;
  // fPixelSizeX = fDimX / fNbinsX;
  Double_t A = fDimY / 2.;
  Double_t B = fDimZ / 2.;
  Double_t F = fDimX / 2.;

  Int_t m, n, j, k;
  Double_t z, y, x;

  TVector3 interactionPoint;
  TVector3 coneAxis;
  Double_t coneTheta;
  Double_t energy1, energy2, energy0;
  Double_t Scatthick_z, Scatthick_x, Absthick_z, Absthick_y, Absthick_x;
  TVector3 *point_e, *point_p, *point_dir, *point_abs_sc;
  TVector3 *Scatposition, *Absposition;
  point_abs_sc = new TVector3();
  fNIpoints = 0;
  fPoints = 0;
  IsectionPoint* tmppoint1;
  IsectionPoint* tmppoint2;

  const Double_t maxdist = sqrt(pow(fPixelSizeY, 2) + pow(fPixelSizeZ, 2));

  bool status;
  int counter = 0;
  int count = 0;
  TStopwatch t;
  t.Start();

/// This loop will always analyze fStop-fStart events starting with
/// event number fStart. If some of the events are not valid they will
/// be skipped, but still fStop-fStart events will be analyzed i.e.
/// last analyzed event will have number fStop+n, where n is number of
/// skipped events. If you want to change this - remove 'counter' variable.

  for (Int_t i = fStart; i < fStop; i++) {
      

      fNIpoints = 0;

      if (fVerbose)
      cout << "CCMLEM::Reconstruct(...) event " << i << endl << endl;

      status = fReader->LoadEvent(counter + i);

      if (status == false) {
         counter++;
         continue;
      }
 
      energy1 = fReader->GetEnergyLoss();
      energy2 = fReader->GetEnergyScattered();
      point_e = fReader->GetPositionScattering();
      //point_e->Print();
      point_p = fReader->GetPositionAbsorption();
      //point_p->Print();
      //point_dir = fReader->GetGammaDirScattered();
      //point_dir->Print();
      //point_abs_sc->SetXYZ(point_p->X()-point_e->X(),
      //point_p->Y()-point_e->Y(), point_p->Z()-point_e->Z());
    
      Scatthick_z = fReader->GetScatThickz();
      Scatthick_x = fReader->GetScatThickx();
      Absthick_z = fReader->GetAbsThickz();
      Absthick_y = fReader->GetAbsThicky();
      Absthick_x = fReader->GetAbsThickx();
      Scatposition = fReader->GetScattererPosition();
      Absposition = fReader->GetAbsorberPosition();
      
///To analyze Geant4 results, set cuts for energies and positions. 
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
       energy0 = 4.40;
/*
       Double_t en = energy1 + energy2;
       if (fabs(energy0 - en) > 1.E-1) {
            //cout << "energy_e : "<< energy1 << "\t" << "+" << "\t" << "energy_p :" << energy2 << "\t" << " = " << "\t" << "total energy : " << en << endl;
//           //cout<< "X of scatterer : " <<point1->X() <<endl;
//           //cout<< "----------------------------------------------------"<<endl;
           continue;
        }
//  
        if(point_e->X() >= (Scatposition->X() - (Scatthick_x/2)) && point_e->X() <=
        (Scatposition->X() + (Scatthick_x/2))){
//           //cout<< "X of scatterer : " <<point1->X() <<endl;
           if(point_p->X() >= (Absposition->X() - (Absthick_x/2)) && point_p->X() <=
           (Absposition->X() + (Absthick_x/2))){
//               //cout<< "X of absorber : " <<point2->X() <<endl;
               if(point_p->Y() >= (-Absthick_y) && point_p->Y() <= Absthick_y){
//                   //cout<< "Y of absorber : " <<point2->Y() <<endl;
               }
               if(point_p->Z() >= (-Absthick_z) && point_p->Z() <= Absthick_z){
//                    //cout<< "Z of absorber : " <<point2->Z() <<endl;
               }
               count++;
           } else {
                //cout<< "out of scatterer : " << point_e->X() << "\t"<< "out of absorber : " <<       point_p->X()<< "\t" << point_p->Y() << "\t" << point_p->Z() <<endl; 
                //cout<< "----------------------------------------------------"<<endl; 
                continue;
           }
// 
// //          cout<< "----------------------------------------------------"<<endl;
// // 
       } else {
// //           cout<< "----------------------------------------------------"<<endl;
            continue;
//  
      }
//    cout<<"---------------------"<<endl;
    //fAngDiff->Fill(point_abs_sc->Angle(*point_dir)*TMath::RadToDeg());
*/   
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
    
    ComptonCone* cone =
        new ComptonCone(point_e, point_p, energy1 + energy2, energy2);
    interactionPoint = cone->GetApex();
    coneAxis = cone->GetAxis();
    coneTheta = cone->GetAngle();
    
    Double_t K = cos(coneTheta);
    Double_t a, b, c1, c2, c3, c4, c, z1, z2;
    Double_t d, e, f1, f2, f3, f4, f, y1, y2;
    x = 0;

     //for (m = 1; m <= fNbinsX; m++) {
    
    y = -A;
    if (fVerbose) cout << "Loop over horizontal lines..." << endl;

    for (j = 0; j < fNbinsY + 1; j++) {

      a = 2 * (pow(-coneAxis.Z(), 2) - pow(K, 2));

      b = 2 * ((-coneAxis.Z()) * ((interactionPoint.Z()) * (-coneAxis.Z()) -
                                  y * ((-coneAxis.Y())) +
                                  (interactionPoint.Y()) * (-coneAxis.Y()) +
                                  (interactionPoint.X()) * (-coneAxis.X()) + x -
                                  (-coneAxis.X()) * x) -
               (interactionPoint.Z()) * pow(K, 2));

      c1 = pow(interactionPoint.Y(), 2) + pow(interactionPoint.X(), 2) +
           ((-2) * pow(interactionPoint.Y(), 2) * pow(-coneAxis.Z(), 2)) +
           ((-2) * pow(interactionPoint.X(), 2) * pow(-coneAxis.Z(), 2)) +
           ((-2) * pow(interactionPoint.Y(), 2) * pow(-coneAxis.Y(), 2)) +
           ((-4) * (interactionPoint.Y()) * (interactionPoint.X()) *
            (-coneAxis.Y()) * (-coneAxis.X())) +
           ((-2) * pow(interactionPoint.X(), 2) * pow(-coneAxis.X(), 2));

      c2 = (-2) * (interactionPoint.X()) * x +
           4 * (interactionPoint.X()) * x * pow(-coneAxis.Z(), 2) -
           4 * (interactionPoint.Y()) * (-coneAxis.Y()) * x -
           4 * (interactionPoint.X()) * (-coneAxis.X()) * x +
           4 * (interactionPoint.Y()) * (-coneAxis.Y()) * (-coneAxis.X()) * x +
           4 * (interactionPoint.X()) * pow(-coneAxis.X(), 2) * x - pow(x, 2) -
           2 * pow(-coneAxis.Z(), 2) * pow(x, 2) +
           4 * (-coneAxis.X()) * pow(x, 2) -
           2 * pow(-coneAxis.X(), 2) * pow(x, 2);

      c3 = pow(y, 2) *
               (1 - 2 * pow(-coneAxis.Z(), 2) - 2 * pow(-coneAxis.Y(), 2)) +
           y * ((interactionPoint.Y()) * ((-2) + 4 * pow(-coneAxis.Z(), 2) +
                                          4 * pow(-coneAxis.Y(), 2)) +
                4 * (-coneAxis.Y()) *
                    ((interactionPoint.X()) * (-coneAxis.X()) + x -
                     (-coneAxis.X()) * x));

      c4 = (pow(y, 2) - 2 * y * (interactionPoint.Y()) +
            pow(interactionPoint.Y(), 2) + pow((interactionPoint.X() - x), 2)) *
           ((-1) + 2 * pow(K, 2));

      c = (-2) * pow(K, 2) * (c1 + c2 + c3 + c4);

      if (c >= 0) {
        z1 = (b - sqrt(c)) / a;
        z2 = (b + sqrt(c)) / a;
        if (fabs(z1) < 1e-14) z1 = z1 * 1e6;
        if (fabs(z2) < 1e-14) z2 = z2 * 1e6;
        AddIsectionPoint("hor", x, y, z1);
        AddIsectionPoint("hor", x, y, z2);
      }
      y = y + fPixelSizeY;

    } // end of loop over horizontal lines

    z = -B;
    if (fVerbose) cout << "Loop over vertical lines..." << endl;
    for (k = 0; k < fNbinsZ + 1; k++) {

      d = 2 * (pow(-coneAxis.Y(), 2) - pow(K, 2));

      e = 2 * ((-coneAxis.Y()) * ((interactionPoint.Z()) * (-coneAxis.Z()) -
                                  z * ((-coneAxis.Z())) +
                                  (interactionPoint.Y()) * (-coneAxis.Y()) +
                                  (interactionPoint.X()) * (-coneAxis.X()) + x -
                                  (-coneAxis.X()) * x) -
               (interactionPoint.Y()) * pow(K, 2));

      f1 = pow(interactionPoint.Z(), 2) + pow(interactionPoint.X(), 2) +
           ((-2) * pow(interactionPoint.Z(), 2) * pow(-coneAxis.Z(), 2)) +
           ((-2) * pow(interactionPoint.Z(), 2) * pow(-coneAxis.Y(), 2)) +
           ((-2) * pow(interactionPoint.X(), 2) * pow(-coneAxis.Y(), 2)) +
           ((-4) * (interactionPoint.Z()) * (interactionPoint.X()) *
            (-coneAxis.Z()) * (-coneAxis.X())) +
           ((-2) * pow(interactionPoint.X(), 2) * pow(-coneAxis.X(), 2));

      f2 = (-2) * (interactionPoint.X()) * x +
           4 * (interactionPoint.X()) * x * pow(-coneAxis.Y(), 2) -
           4 * (interactionPoint.Z()) * (-coneAxis.Z()) * x -
           4 * (interactionPoint.X()) * (-coneAxis.X()) * x +
           4 * (interactionPoint.Z()) * (-coneAxis.Z()) * (-coneAxis.X()) * x +
           4 * (interactionPoint.X()) * pow(-coneAxis.X(), 2) * x - pow(x, 2) -
           2 * pow(-coneAxis.Y(), 2) * pow(x, 2) +
           4 * (-coneAxis.X()) * pow(x, 2) -
           2 * pow(-coneAxis.X(), 2) * pow(x, 2);

      f3 = pow(z, 2) *
               (1 - 2 * pow(-coneAxis.Z(), 2) - 2 * pow(-coneAxis.Y(), 2)) +
           z * ((interactionPoint.Z()) * ((-2) + 4 * pow(-coneAxis.Z(), 2) +
                                          4 * pow(-coneAxis.Y(), 2)) +
                4 * (-coneAxis.Z()) *
                    ((interactionPoint.X()) * (-coneAxis.X()) + x -
                     (-coneAxis.X()) * x));

      f4 = (pow(z, 2) - 2 * z * (interactionPoint.Z()) +
            pow(interactionPoint.Z(), 2) + pow((interactionPoint.X() - x), 2)) *
           ((-1) + 2 * pow(K, 2));

      f = (-2) * pow(K, 2) * (f1 + f2 + f3 + f4);

      if (f >= 0) {
        y1 = (e - sqrt(f)) / d;
        y2 = (e + sqrt(f)) / d;
        if (fabs(y1) < 1e-14) y1 = y1 * 1e6;
        if (fabs(y2) < 1e-14) y2 = y2 * 1e6;
        AddIsectionPoint("ver", x, y1, z);
        AddIsectionPoint("ver", x, y2, z);
      }
      z = z + fPixelSizeZ;

    } // end of loop over vertical lines

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
        // continue;
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

        cout << "eventno :" << i << endl;
        cout << binno1 << "!=" << binno2
             << " ->Bin numbers are different when they should not!" << endl;
      }

      fImage[0]->SetBinContent(binno1, fImage[0]->GetBinContent(binno1) + dist);

      temp = (SMElement*)fSM->ConstructedAt(fPoints++);
      temp->SetEvBinDist(i, binno1, dist);

      if (fVerbose) temp->Print();
    }
    //cout<< "no. of Ipoints : " << fNIpoints<<endl; 
    //x = x + fPixelSizeX;
    //}
    

    delete cone;

    if (fVerbose)
      cout << "----------------------------------------------------------------"
           << endl;

  } // end of loop over events
  
  
  //cout<< "no. of Points : " << fPoints<<endl;
  //SaveToFile(fAngDiff);
  fArray->Clear("C");
  //cout<< "number of events : " << count <<endl;
  SaveToFile(fImage[0]);
  
  fSigma [100] = 1.;
  for (int iter = 1; iter < fIter + 1; iter++) {
      
      Iterate(fStop, iter);
      
      if (fSigma[100] < 0.01) {
          t.Stop();
          t.Print();
          //cout<< "sigma value" << fSigma[100] << endl;
          return 0;
      }


  }

  //DrawHisto();
  //GetSigmaError();
   //t.Stop();
  // t.Print();

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

  if (fabs(y) > fDimY / 2. || fabs(z) > fDimZ / 2.) {
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
    yplus = y + 0.0005 * fPixelSizeY;
    yminus = y - 0.0005 * fPixelSizeY;

    if (fabs(yplus) <= fDimY / 2) {
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);
      pixelY = fImage[0]->GetYaxis()->FindBin(yplus);
      tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ, pixelY), x, y, z);
      // if(fVerbose) tmppoint->Print();
      // fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
    }
    if (fabs(yminus) <= fDimY / 2) {
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);
      pixelY = fImage[0]->GetYaxis()->FindBin(yminus);
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
    zplus = z + 0.0005 * fPixelSizeZ;
    zminus = z - 0.0005 * fPixelSizeZ;

    if (fabs(zplus) <= fDimZ / 2) {
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);

      pixelZ = fImage[0]->GetXaxis()->FindBin(zplus);
      tmppoint->SetBinPoint(fImage[0]->GetBin(pixelZ, pixelY), x, y, z);
      // if(fVerbose) tmppoint->Print();
      // fGraph->SetPoint(fGraph->GetN(), z, y);
      added++;
    }
    if (fabs(zminus) <= fDimZ / 2) {
      tmppoint = (IsectionPoint*)fArray->ConstructedAt(fNIpoints++);

      pixelZ = fImage[0]->GetXaxis()->FindBin(zminus);
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
  gStyle->SetOptFit(1011);
    
  int lastiter = iter - 1;

  double sigma[lastiter + 1];
  TH1D* ProZ[150];
  TH1D* ProY[150];
  //TH1D* ProX[150];
  
///For the unknown source distribution, user should define an appropriate fitting
///function to compare reasonablely sigma value of the new 10th iteration with previous one. 
  
   TF1* func_z = new TF1("FitProjection", FitProjection, -14.5, 14.5, 3);
   func_z->SetParameters(fP0,fP1,fP2);
   func_z->SetParNames("Constant","Mean_value","Sigma_z");
   //TF1* func_y = new TF1("FitProjection", FitProjection, -14.5, 14.5, 3);
   //func_y->SetParameters(fP0,fP1,fP2);
   //func_y->SetParNames("Constant_y","Mean_value_y","Sigma_y");
//   func->SetParameter(0, 18000);
//   func->SetParameter(1, 0.01);
//   func->SetParameter(2, 0.1);
  // func->SetParameter(3,8.1);

  if (fImage[lastiter] == NULL) {
    cout << "Error in CCMELM::Iterate(). Last iteration NULL" << endl;
    return kFALSE;
  }

  TH2F* hlastiter = (TH2F*)fImage[lastiter];
  fImage[lastiter + 1] = (TH2F*)hlastiter->Clone();
  TH2F* hthisiter = fImage[lastiter + 1];
  hthisiter->Reset();
  hthisiter->SetName(Form("%s_iter%i", fImage[0]->GetName(), lastiter + 1));
  hthisiter->SetTitle(Form("%s_iter%i", fImage[0]->GetTitle(), lastiter + 1));
  Int_t eventno;
  Int_t eventno_prev = 0;
  Int_t entry;
  Int_t binno;
  Double_t dist, addvalue;
  SMElement* temp;
  Double_t denominator[nstop + 1];
  for (int i = 0; i < nstop + 1; i++)
    denominator[i] = 0;
  Int_t nSMentries = fSM->GetEntries();
  // cout<<"nSMentries = "<<nSMentries <<endl;
  for (entry = 0; entry < nSMentries; entry++) {
    temp = (SMElement*)fSM->At(entry);
    binno = temp->GetBin();
    eventno = temp->GetEvent();
    dist = temp->GetDist();
    denominator[eventno] =
        denominator[eventno] + dist * hlastiter->GetBinContent(binno);
  }

  for (entry = 0; entry < nSMentries; entry++) {
    addvalue = 0;
    temp = (SMElement*)fSM->At(entry);
    binno = temp->GetBin();
    eventno = temp->GetEvent();
    dist = temp->GetDist();
    addvalue = dist * hlastiter->GetBinContent(binno) / denominator[eventno];
    hthisiter->SetBinContent(binno, hthisiter->GetBinContent(binno) + addvalue);
  }
  
  

  for (int i = 10; i <= lastiter + 1; i = i + 10) {
    
    //ProX[i] = fImage[i]->ProjectionZ();  
    ProY[i] = fImage[i]->ProjectionY();
    ProZ[i] = fImage[i]->ProjectionX();
    ProZ[i]->Fit(func_z, "r");
    //ProY[i]->Fit(func_y,"r");
    sigma[i] = func_z->GetParameter(2);
    // cout<< "i : \t" << i << "\t" << "sigma : \t " << sigma[i]<< endl;
  }
///Comparing the sigma_z of new 10th iteration with previous one
///to check if it continues iterating more or not. 
  for (int j = 20; j <= lastiter + 1; j = j + 10) {
       //cout << fSigma[100] << endl;
      if (fabs(sigma[j]) !=0) fSigma[100] = (fabs(sigma[j] - sigma[j - 10]))/fabs(sigma[j]);
       //cout << fSigma[100] << endl;
      if (fSigma[100] < 0.01) {
          
          cout << "SigmaIter" << j << " - " << "SigmaIter" << j - 10 << " = " << "SigmaPercentageErr : " << fSigma[100] << endl;
 
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
            
            
  SaveToFile(hthisiter);
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

  while (!config.eof()) { // TODO: this reads past the last line, so will always
						  //give a warning about unknown syntax of an empty comment
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
      if (fDimZ < 1 || fDimY < 1) {
        cout << "##### Error in CCMLEM::Config()! Image size incorrect!"
             << endl;
        return false;
      }
    } else if (comment.Contains("No. of bins")) {
      config >> fNbinsZ >> fNbinsY >> fNbinsX;
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
    } else if (comment.Contains("Event filter Geant4")) {
      config >> fGeantFilter;
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
      if (fabs(sigma[j]) !=0) fSigma[100] = (fabs(sigma[j] - sigma[j - 10]))/fabs(sigma[j]);
       //cout << fSigma[100] << endl;
      if (fSigma[100] < 0.01) {
 
          cout << "SigmaIter" << j << " - " << "SigmaIter" << j - 10 << " = " << "SigmaPercentageErr : " << fSigma[100] << endl;
 
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
///Prints details of the CCMLEM class obejct.
void CCMLEM::Print(void) {
  cout << "\nCCMLEM::Print()" << endl;
  cout << setw(35) << "Name of input file: \t"
       << fInputName << endl;
  // cout << setw(35) << "Center of reco plane: \t" << fXofRecoPlane << ", "
  // << fYofRecoPlane << ", " << fZofRecoPlane << endl;
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
  cout << setw(35) << "Event filter Geant4: \t" << fGeantFilter << endl;
  cout << setw(35) << "Verbose level: \t" << fVerbose << endl << endl;
}
//--------------------------------------
///Sets default values of the private class members.
void CCMLEM::Clear(void) {
  fInputName = "dummy";
  //fSigma[100] = 1;
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
  fGeantFilter = 0;
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
