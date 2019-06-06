#include "CCMLEM.hh"
#include "IsectionPoint.hh"
#include "SMElement.hh"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TRandom.h"
#include "TStyle.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
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
  TString outputName = "../work/results/" + fname;
  outputName.ReplaceAll("CCSimulation", "CCReconstruction");
  TString option = (fFreshOutput ? "RECREATE" : "UPDATE");
  fOutputFile = new TFile(outputName, option);

  fArray = new TClonesArray("IsectionPoint", 1000);
  fSM = new TClonesArray("SMElement", 1000000);
  fNIpoints = 0;
  fPoints = 0;

  if (fSmear) {
     TString name = "../work/results/EnergyResolutionExample.root";
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
    std::cout << "par0 " << func->GetParameter(0) << std::endl;
    std::cout << "par1 " << func->GetParameter(1) << std::endl;
    std::cout << "par2 " << func->GetParameter(2) << std::endl;
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

  TString fullName = "../work/results/" + fInputName;

  TFile* file = new TFile(fullName, "READ");
  if (!file->IsOpen()) {
    cout << "##### Error in CCMLEM::SetInputReader!" << endl;
    cout << "Could not open requested file" << endl;
    return false;
  }

  if (file->Get("data")) {
    file->Close();
    fReader = new InputReaderSimple(fullName);
  } else if (file->Get("G4SimulationData_Real") &&
             file->Get("G4SimulationData_Setup")&&
             file->Get("G4SimulationData_Reconstruction")) {
    file->Close();
    fReader = new InputReaderGeant(fullName);
  } else if (file->Get("Cluster21") && file->Get("Cluster22")) {
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
Bool_t CCMLEM::Reconstruct(void) {

  // image histogram
  fImage[0] = new TH2F("image", "image", fNbinsZ, -fDimZ / 2., fDimZ / 2.,
                       fNbinsY, -fDimY / 2., fDimY / 2.);
  fImage[0]->GetXaxis()->SetTitle("z [mm]");
  fImage[0]->GetYaxis()->SetTitle("y [mm]");

  //   fAngDiff = new TH1F("Angle Difference", "Angle Difference", 20, -5, 5);
  //   fAngDiff->GetXaxis()->SetTitle("Angle Diff.(deg)");
  //   fAngDiff->GetYaxis()->SetTitle("events");

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
  Double_t Scatthick_z, Absthick_z, Absthick_y, Absthick_x;
  TVector3 *point1, *point2;
  TVector3 *Scatposition, *Absposition;

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

  // This loop will always analyze fStop-fStart events starting with
  // event number fStart. If some of the events are not valid they will
  // be skipped, but still fStop-fStart events will be analyzed i.e.
  // last analyzed event will have number fStop+n, where n is number of
  // skipped events. If you want to change this - remove 'counter' variable.

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
    point1 = fReader->GetPositionScattering();
    // point1->Print();
    point2 = fReader->GetPositionAbsorption();
    /* point2 = fReader->GetGammaDirScattered();
      //point2->Print();
     Scatthick_z = fReader->GetScatThickz();
     //cout<< Scatthick_z<<endl;
     Absthick_z = fReader->GetAbsThickz();
     Absthick_y = fReader->GetAbsThicky();
     Absthick_x = fReader->GetAbsThickx();
     //cout<<  Absthick_z<<endl;
     Scatposition = fReader->GetScattererPosition();
     Absposition = fReader->GetAbsorberPosition();

     Double_t DistSA = Absposition->Z() - Scatposition->Z();

     point2->SetXYZ(point1->x()+(point2->x()*DistSA),
     point1->y()+(point2->y()*DistSA), point1->z()+(point2->z()*DistSA));

     energy0 = 4.40;
     Double_t en = energy1 + energy2;
     if (fabs(energy0 - en) > 1.E-1) {
         //cout << "energy_e : "<< energy1 << "\t" << "+" << "\t" << "energy_p :
     " << energy2 << "\t" << " = " << "\t" << "total energy : " << en << endl;
         //cout<< "X of scatterer : " <<point1->X() <<endl;
         //cout<< "----------------------------------------------------"<<endl;
         continue;
     }

     if(point1->X() >= (Scatposition->Z() - (Scatthick_z/2)) && point1->X() <=
     (Scatposition->Z() + (Scatthick_z/2))){
         //cout<< "X of scatterer : " <<point1->X() <<endl;
         if(point2->X() >= (Absposition->Z() - (Absthick_z/2)) && point2->X() <=
     (Absposition->Z() + (Absthick_z/2))){
             //cout<< "X of absorber : " <<point2->X() <<endl;
             if(point2->Y() >= (-Absthick_y) && point2->Y() <= Absthick_y){
                //cout<< "Y of absorber : " <<point2->Y() <<endl;
             }
             if(point2->Z() >= (-Absthick_x) && point2->Z() <= Absthick_x){
                  //cout<< "Z of absorber : " <<point2->Z() <<endl;
             }
             count++;
         } else {
              //cout<< "out of scatterer : " << point1->X() << "\t"<< "out of
     absorber : " << point2->X()<< "\t" << point2->Y() << "\t" << point2->Z()
     <<endl;

              //cout<<
     "----------------------------------------------------"<<endl; continue;
         }


         //cout<< "----------------------------------------------------"<<endl;

     } else {
          //cout<< "----------------------------------------------------"<<endl;
          continue;

     }*/
    //
    //     Double_t mag1 = point1->X()*point2->X() + point1->Y()*point2->Y() +
    //     point1->Z()*point2->Z(); Double_t mag2 = point1->Mag() *
    //     point2->Mag(); Double_t costheta = mag1 / mag2; Double_t theta =
    //     (TMath::ACos(costheta))* 180/M_PI;
    //     //fGeoAng->Fill(theta);
    //     //cout<< "theta = "<< theta <<endl;
    // if (energy1>=0.6){
    // point2->SetXYZ(point1->x()+point2->x(), point1->y()+point2->y(),
    // point1->z()+point2->z()); point2->Print();
    if (fSmear) {
          point1->SetXYZ(SmearBoxX(point1->X()),
                         SmearGaus(point1->Y(), fResolutionY),
                         SmearBoxZ(point1->Z()));
          point2->SetXYZ(SmearBoxX(point2->X()),
                         SmearGaus(point2->Y(), fResolutionY),
                         SmearBoxZ(point2->Z()));
      energy1 = SmearGaus(energy1, GetSigmaE(energy1));
      energy2 = SmearGaus(energy2, GetSigmaE(energy2));
    }

    ComptonCone* cone =
        new ComptonCone(point1, point2, energy1 + energy2, energy2);
    interactionPoint = cone->GetApex();
    coneAxis = cone->GetAxis();
    coneTheta = cone->GetAngle();

    Double_t K = cos(coneTheta);
    Double_t a, b, c1, c2, c3, c4, c, z1, z2;
    Double_t d, e, f1, f2, f3, f4, f, y1, y2;
    x = 0;

    // for (m = 1; m <= fNbinsX; m++) {
    // fNIpoints = 0;
    // fPoints   = 0;
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
      Int_t binno1, binno2;
      SMElement* temp;

      for (int h = 0; h < fNIpoints; h = h + 2) {

        if (fVerbose)
          cout << " index[" << h << "]=" << index[h] << ", index[" << h + 1
               << "]=" << index[h + 1] << endl;
        tmppoint1 = (IsectionPoint*)fArray->At(index[h]);
        tmppoint2 = (IsectionPoint*)fArray->At(index[h + 1]);
        if (tmppoint1 == NULL || tmppoint2 == NULL) {
          cout << "Something went wrong" << tmppoint1 << "\t" << tmppoint2
               << endl;
        }
        tmpvec1 = tmppoint1->GetPointCoordinates();
        binno1 = tmppoint1->GetBin();
        tmpvec2 = tmppoint2->GetPointCoordinates();
        binno2 = tmppoint2->GetBin();
        // if(fVerbose)  cout<<" binno1="<<binno1<<",
        // binno2="<<binno2<<endl<<endl;
        dist = ((*tmpvec1) - (*tmpvec2)).Mag();
        if (dist > maxdist) {
          // cout<<"Event "<<h<<": distance exceeds pixel diagonal
          // "<<dist/maxdist<<" times"<<endl;
          continue;
        }
        if (binno1 != binno2) {
          // cout<<binno1<<"!="<<binno2<<" ->Bin numbers are different when they
          // should not!"<<endl;
        }

        fImage[0]->SetBinContent(binno1,
                                 fImage[0]->GetBinContent(binno1) + dist);
        if (fVerbose) cout << "fPoints = " << fPoints << endl;
        temp = (SMElement*)fSM->ConstructedAt(fPoints++);
        temp->SetEvBinDist(i, binno1, dist);
        if (fVerbose) temp->Print();
      }
      // cout<< fNIpoints<<endl;
      x = x + fPixelSizeX;
    }

    // if(fVerbose) cout<<"end of loop"<<endl;

    delete cone;

    if (fVerbose)
      cout << "----------------------------------------------------------------"
           << endl;

  } // end of loop over events
  //}

  // SaveToFile(fAngDiff);
  fArray->Clear("C");

  SaveToFile(fImage[0]);

  for (int iter = 1; iter < fIter + 1; iter++) {
      
      Iterate(fStop, iter);
      if (fSigma[100] < 0.01) {
          
          t.Stop();
          t.Print();
// //       // cout<< "sigma value" <<fSigma[100] << endl;
          return 0;
      }
  }

  // DrawHisto();
  // GetSigmaError();
  // t.Stop();
  // t.Print();

  return kTRUE;
}

//------------------------------------
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
    // pixelX = fImage[0]->GetZaxis()->FindBin(x);
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
    // pixelX = fImage[0]->GetZaxis()->FindBin(x);
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
Bool_t CCMLEM::Iterate(Int_t nstop, Int_t iter) {

  //gStyle->SetOptStat(1101);
  gStyle->SetOptFit(1011);
    
  int lastiter = iter - 1;

  double sigma[lastiter + 1];
  TH1D* ProZ[150];
  TH1D* ProY[150];
  
  // For the unknown source distribution, user should define an appropriate fitting
  // function to compare reasonablely sigma value of the new 10th iteration with previous one. 
  
  TF1* func = new TF1("FitProjection", FitProjection, -14.5, 14.5, 3);
  func->SetParameters(fP0,fP1,fP2);
  func->SetParNames("Constant","Mean_value","Sigma_z");
//   TF1* func_y = new TF1("FitProjection", FitProjection, -14.5, 14.5, 3);
//   func_y->SetParameters(fP0,fP1,fP2);
//   func_y->SetParNames("Constant_y","Mean_value_y","Sigma_y");
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
    // sigma[j]=0;
    ProY[i] = fImage[i]->ProjectionY();
    ProZ[i] = fImage[i]->ProjectionX();
    ProZ[i]->Fit(func,"r");
    //ProY[i]->Fit(func_y,"r");
    sigma[i] = func->GetParameter(2);
    // cout<< "i : \t" << i << "\t" << "sigma : \t " << sigma[i]<< endl;
  }

    for (int j = 20; j <= lastiter + 1; j = j + 10) {
      //cout << fSigma[100] << endl;
      if (fabs(sigma[j]) !=0) fSigma[100] = (fabs(sigma[j] - sigma[j - 10]))/fabs(sigma[j]);
      //cout << fSigma[100] << endl;
      if (fSigma[100] < 0.01) {

        cout << "SigmaIter" << j << " - "
           << "SigmaIter" << j - 10 << " = "
           << "SigmaPercentageErr : " << fSigma[100] << endl;

        TCanvas* can = new TCanvas("can", "MLEM2D", 1000, 1000);
        can->Divide(1, 2);
        can->cd(1);
        //fImage[j - 10]->Draw("colz");
        //can->cd(2);
        ProZ[j - 10]->Draw();
        can->cd(2);
        ProY[j - 10]->Draw();

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
Bool_t CCMLEM::DrawHisto(void) {

  int lastiter = fIter;
  // int iter = fIter - 10;

  TH1D* hProZ[150];
  TH1D* hProY[150];
  // TH1D* hProX[100];
  TCanvas* can = new TCanvas("MLEM2D", "MLEM2D", 1000, 1000);
  // TCanvas* canz = new TCanvas("MLEM1DZ", "MLEM1DZ", 1000, 1000);
  // TCanvas* cany = new TCanvas("MLEM1DY", "MLEM1DY", 1000, 1000);
  // TCanvas* canx = new TCanvas("MLEM1DX", "MLEM1DX", 1000, 1000);
  can->Divide(2, 2);
  // canz->Divide(1,1);
  // cany->Divide(1,1);
  // canx->Divide(1,1);
  // canz->Divide((int)sqrt(fIter) + 1, (int)sqrt(fIter) + 1);
  // cany->Divide((int)sqrt(fIter) + 1, (int)sqrt(fIter) + 1);
  // for (int iter = 0; iter < fIter + 1; iter++) {
  can->cd(1);
  // gPad->SetLogz(1);
  // fImage[iter]->Draw("colz");
  // can->cd(2);
  fImage[lastiter]->Draw("colz");
  // fImage[iter]->SetMinimum(200);
  // hProZ[iter] = fImage[iter]->ProjectionX();
  // hProY[iter] = fImage[iter]->ProjectionY();
  hProZ[lastiter] = fImage[lastiter]->ProjectionX();
  hProY[lastiter] = fImage[lastiter]->ProjectionY();
  // hProX[lastiter] = fImage[lastiter]->ProjectionZ();
  // can->cd(2);
  // hProZ[iter]->Draw();
  can->cd(2);
  hProZ[lastiter]->Draw();
  // can->cd(5);
  // hProY[iter]->Draw();
  can->cd(3);
  hProY[lastiter]->Draw();
  // can->cd(4);
  // hProX[lastiter]->Draw();
  //}

  SaveToFile(can);
  // SaveToFile(canz);
  // SaveToFile(cany);
  // SaveToFile(canx);
  return kTRUE;
}
//------------------------------------
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
  //     for(int i=10; i < fIter + 1; i=i+10){
  //
  //         //sigma[i]=0;
  //         ProY[i] = fImage[i]->ProjectionY();
  //         ProZ[i] = fImage[i]->ProjectionX();
  //         ProZ[i]->Fit("func");
  //         sigma[i] = func->GetParameter(3);
  //
  //         cout<< "i : \t" << i << "\t" << "sigma : \t " << sigma[i]<< endl;
  //
  //     }
  for (int j = 10; j < fIter + 1; j = j + 10) {
    sigma[j] = 0;
    ProY[j] = fImage[j]->ProjectionY();
    ProZ[j] = fImage[j]->ProjectionX();
    ProZ[j]->Fit("func");
    sigma[j] = func->GetParameter(3);
    // cout<< "j : \t" << j << "\t" << "sigma : \t " << sigma[j]<< endl;

    if (fabs(sigma[j] - sigma[j - 10]) <= 0.5) {

      cout << "SigmaIter" << j << " - "
           << "SigmaIter" << j - 10 << " = "
           << "SigmaErr : " << fabs(sigma[j] - sigma[j - 10]) << endl;
      TCanvas* can = new TCanvas("MLEM2D", "MLEM2D", 1000, 1000);
      can->Divide(2, 2);
      can->cd(1);
      fImage[j - 10]->Draw("colz");
      can->cd(2);
      ProZ[j - 10]->Draw();
      can->cd(3);
      ProY[j - 10]->Draw();

      SaveToFile(can);

      return 0;
    }
  }

  return kTRUE;
}
//-------------------------------------
Double_t CCMLEM::SmearGaus(double val, double sigma) {
  return gRandom->Gaus(val, sigma);
}
//------------------------------------
Double_t CCMLEM::SmearBoxX(double x) {
  return gRandom->Uniform(x - (fResolutionX / 2), x + (fResolutionX / 2));
}
//------------------------------------
Double_t CCMLEM::SmearBoxZ(double z) {
  return gRandom->Uniform(z - (fResolutionZ / 2), z + (fResolutionZ / 2));
}
//------------------------------------
Double_t CCMLEM::GetSigmaE(double energy) {
  double sigma = fHisto->GetFunction("fit1")->Eval(energy) * energy;
  return sigma;
}
//-------------------------------------------
Double_t CCMLEM::FitProjection(Double_t *x,Double_t *p) {
      Double_t a = 0;
      if (p[2]!=0) a = (x[0] - p[1])/p[2];
      Double_t fitvalue = p[0]*TMath::Exp(-0.5*a*a);
      return fitvalue;
}
//-------------------------------------------
/// Prints details of the CCMLEM class obejct.
void CCMLEM::Print(void) {
  cout << "\nCCMLEM::Print()" << endl;
  cout << setw(35) << "Name of input file: \t"
       << "../work/results/" + fInputName << endl;
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
  cout << setw(35) << "Verbose level: \t" << fVerbose << endl << endl;
}
//--------------------------------------
/// Sets default values of the private class members.
void CCMLEM::Clear(void) {
  fInputName = "dummy";
  fSigma[100] = 1;
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
  fGraph = NULL;
}
//--------------------------------------
/// Saves object in the output file.
///\param ob (TObject*) - saved object.
Bool_t CCMLEM::SaveToFile(TObject* ob) {
  fOutputFile->cd();
  ob->Write();
  cout << ob->ClassName() << " " << ob->GetName() << " saved in the file "
       << fOutputFile->GetName() << endl;
  return kTRUE;
}
//--------------------------------------
