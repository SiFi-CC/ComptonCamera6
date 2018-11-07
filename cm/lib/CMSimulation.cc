#include "CMSimulation.hh"
#include "TCanvas.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"
#include <fstream>
#include <iostream>
using namespace std;

ClassImp(CMSimulation);
//------------------------------------------------------------------
CMSimulation::CMSimulation() {
  fFile = 0;
  fTree = 0;
  fVerbLevel = 0;
  fNev = 0;
}
//------------------------------------------------------------------
CMSimulation::CMSimulation(TString name, Int_t vlevel) {
  SetName(name);
  fVerbLevel = vlevel;
  fFile = new TFile("results/" + name + ".root", "RECREATE");
  fTree = new TTree("data", "data");
  fTree->Branch("point0", &fPoint0);
  fTree->Branch("point1", &fPoint1);
  fTree->Branch("point2", &fPoint2);
  fTree->Branch("dir", &fDir);
  fTree->Branch("opaque", &fOpaque, "opaque/I");
  fNev = 0;
  fGenVersion = -1000;
}
//------------------------------------------------------------------
CMSimulation::~CMSimulation() {
  if (fVerbLevel) cout << "Inside destructor of CMsimulation class" << endl;
  if (fTree && fFile) {
    fTree->Write();
    fMask.GetPattern()->Write();
    hYZ->Write();
    hYZdetected->Write();
    hCosTheta->Write();
    hPhi->Write();
    TH1D* h1 = hYZ->ProjectionX("hZ");
    h1->SetTitle("hZ of vertex");
    h1->SetFillColor(kRed);
    h1->SetLineColor(kRed);
    h1->SetFillStyle(3004);
    h1->Write();
    TH1D* h2 = hYZdetected->ProjectionX("hZdetected");
    h2->SetTitle("hZ detected");
    h2->SetFillColor(kBlue);
    h2->SetLineColor(kBlue);
    h2->SetFillStyle(3005);
    h2->Write();

    TCanvas* a = new TCanvas("summary", "summary");
    a->Divide(3, 2);
    a->cd(1);
    hCosTheta->Draw();
    a->cd(4);
    hPhi->Draw();
    a->cd(2);
    hYZ->Draw("colz");
    a->cd(5);
    hYZdetected->Draw("colz");
    a->cd(3);
    h1->Draw("");
    a->cd(6);
    h2->Draw("");
    a->Write();
    fFile->Close();
  }
  SaveGeometryTxt();
  BuildTGeometry();
}
//------------------------------------------------------------------
void CMSimulation::BuildSetup(double maskdist = 50, double maskZsize = 300,
                              double maskYsize = 300, double detdist = 70,
                              double detZsize = 300, double detYsize = 300) {

  fMask.SetPlane(1, 0, 0, maskdist);
  fMask.SetDimensions(maskZsize, maskYsize);
  fMask.SetName("mask");
  //  TH2F* h = new TH2F("h","pattern",10,0,1,10,0,1);
  //  SetPattern(h);
  //  delete h;

  fDetPlane.SetPlane(1, 0, 0, detdist);
  fDetPlane.SetDimensions(detZsize, detYsize);
  fDetPlane.SetName("detector");
}
//------------------------------------------------------------------
void CMSimulation::SetupSpectra(void) {

  double maskZdim = fMask.GetDimZ() / 2;
  double maskYdim = fMask.GetDimY() / 2;
  double detZdim = fDetPlane.GetDimZ() / 2;
  double detYdim = fDetPlane.GetDimY() / 2;
  int nbinsz = fMask.GetPattern()->GetXaxis()->GetNbins();
  int nbinsy = fMask.GetPattern()->GetYaxis()->GetNbins();
  hYZ = new TH2F("hYZ", "Y vs Z for vertex", nbinsz, -maskZdim, maskZdim,
                 nbinsy, -maskYdim, maskYdim);
  hYZdetected = new TH2F("hYZdetected", "Y vs Z detected", nbinsz, -maskZdim,
                         maskZdim, nbinsy, -maskYdim, maskYdim);
  hCosTheta = new TH1F("hCosTheta", "cos(theta) at vertex", 100, -1, 1);
  hPhi = new TH1F("hPhi", "phi at vertex", 100, -TMath::Pi(), TMath::Pi());
}
//------------------------------------------------------------------
void CMSimulation::ClearSpectra(void) {
  hYZ->Reset();
  hYZdetected->Reset();
  hCosTheta->Reset();
  hPhi->Reset();
  fNev = 0;
}
//------------------------------------------------------------------
Bool_t CMSimulation::GenerateRay(void) {
  Double_t theta, phi, maxz;
  switch (fGenVersion) {
    case 1: // isotropic point-like source
      fPoint0 = fSourcePos;
      theta = acos(gRandom->Uniform(-1, 1));
      phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
      // cout<<"fPoint0 ";
      // fPoint0.Print();
      break;
    case 2: // uniform distribution along z, on beam axis, right-angled to it
      maxz = fDetPlane.GetDimZ() / 2.;
      fPoint0.SetXYZ(0, 0, gRandom->Uniform(-maxz, maxz));
      phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
      theta = TMath::Pi() / 2.;
      break;
    case 3: // isotropic two point-like sources, 20mm gap
      if (gRandom->Uniform(0, 1) < 0.5)
        fPoint0.SetXYZ(0, 0, 0);
      else
        fPoint0.SetXYZ(0, 0, -20);
      theta = acos(gRandom->Uniform(-1, 1));
      phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
      break;
    case 4: // isotropic two point-like sources, 40mm gap
      if (gRandom->Uniform(0, 1) < 0.5)
        fPoint0.SetXYZ(0, 0, 0);
      else
        fPoint0.SetXYZ(0, 0, -40);
      theta = acos(gRandom->Uniform(-1, 1));
      phi = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
      break;
    default:
      cout << "You have to choose properly the version of event generator"
           << endl;
      break;
  }
  fDir.SetMagThetaPhi(1., theta, phi);
  return kTRUE;
}
//------------------------------------------------------------------
Bool_t CMSimulation::ProcessEvent() {
  Clear();
  GenerateRay();
  if (fabs(fDir.Phi()) < TMath::Pi() / 2.) {
    if (fVerbLevel) cout << "Phi out of range" << endl;
    return kFALSE;
  }
  fTrack.SetPoint(fPoint0);
  fTrack.SetVersor(fDir);
  if (!fTrack.FindCrossPoint(&fMask, fPoint1)) {
    if (fVerbLevel) cout << "No cross point with Mask" << endl;
    return kFALSE;
  }
  if (!fTrack.FindCrossPoint(&fDetPlane, fPoint2)) {
    if (fVerbLevel) cout << "No cross point with DetPlane" << endl;
    return kFALSE;
  }
  fOpaque = fMask.IsOpaque(fPoint1);
  fTree->Fill();
  hYZ->Fill(fPoint0.Z(), fPoint0.Y());
  hCosTheta->Fill(fDir.CosTheta());
  hPhi->Fill(fDir.Phi());
  if (fOpaque == 0) hYZdetected->Fill(fPoint2.Z(), fPoint2.Y());
  fNev++;
  return kTRUE;
}
//------------------------------------------------------------------
void CMSimulation::Loop(Int_t nev) {
  Int_t i = 0;
  while (fNev < nev) {
    ProcessEvent();
    i++;
    if (fVerbLevel)
      cout << "i=" << i << " events simulated, " << fNev << " in acceptance"
           << endl;
  }
}
//------------------------------------------------------------------
void CMSimulation::Clear(void) {
  fPoint0.SetXYZ(-1000, -1000, -1000);
  fPoint1.SetXYZ(-1000, -1000, -1000);
  fPoint2.SetXYZ(-1000, -1000, -1000);
  fDir.SetXYZ(-1000, -1000, -1000);
  fOpaque = -1000;
}
//------------------------------------------------------------------
void CMSimulation::SaveGeometryTxt(void) {
  ofstream output("results/" + fName + "_geometry.txt",
                  std::ios::out | std::ios::trunc);
  output << "TVector3 object: fPoint0" << endl;
  output << fPoint0.X() << " " << fPoint0.Y() << " " << fPoint0.Z() << endl;
  output << "Mask object: " << fMask.GetName() << endl;
  output << "fA = " << fMask.GetA() << endl;
  output << "fB = " << fMask.GetB() << endl;
  output << "fC = " << fMask.GetC() << endl;
  output << "fD = " << fMask.GetD() << endl;
  output << "fDimZ = " << fMask.GetDimZ() << endl;
  output << "fDimY = " << fMask.GetDimY() << endl;
  output << "DetPlane object: " << fDetPlane.GetName() << endl;
  output << "fA = " << fDetPlane.GetA() << endl;
  output << "fB = " << fDetPlane.GetB() << endl;
  output << "fC = " << fDetPlane.GetC() << endl;
  output << "fD = " << fDetPlane.GetD() << endl;
  output << "fDimZ = " << fDetPlane.GetDimZ() << endl;
  output << "fDimY = " << fDetPlane.GetDimY() << endl;
  output << "Generator Version: = " << fGenVersion << endl;
  if (fVerbLevel)
    cout << "Geometry of the setup saved in the text file" << endl;
}
//------------------------------------------------------------------
void CMSimulation::BuildTGeometry(void) {

  TGeoManager* geom = new TGeoManager("world", "CMSimulation - world");

  //----- define materials
  TGeoMaterial* matVacuum = new TGeoMaterial("Vacuum", 0, 0, 0);

  //----- define media
  TGeoMedium* Vacuum = new TGeoMedium("Vacuum", 1, matVacuum);

  //----- define translations
  Double_t maskA = fMask.GetA();
  Double_t maskD = fMask.GetD();
  Double_t detA = fDetPlane.GetA();
  Double_t detD = fDetPlane.GetD();
  Double_t maskTrans = -maskD / maskA;
  Double_t detTrans = -detD / detA;
  Double_t x = fPoint0.X();
  Double_t y = fPoint0.Y();
  Double_t z = fPoint0.Z();

  TGeoTranslation* tr1 = new TGeoTranslation(x, y, z);
  TGeoTranslation* tr2 = new TGeoTranslation(maskTrans, 0., 0.);
  TGeoTranslation* tr3 = new TGeoTranslation(detTrans, 0., 0.);

  //----- define dimensions of objects
  Double_t maskDimZ = fMask.GetDimZ() / 2.; // half length
  Double_t maskDimY = fMask.GetDimY() / 2.;
  Double_t detDimZ = fDetPlane.GetDimZ() / 2.;
  Double_t detDimY = fDetPlane.GetDimY() / 2.;

  //----- make the top container volume
  Double_t worldX = -detTrans + 50; // cm hafl length
  Double_t worldY = detDimY + 50.;
  Double_t worldZ = detDimZ + 50;
  TGeoVolume* top = geom->MakeBox("top", Vacuum, worldX, worldY, worldZ);
  geom->SetTopVolume(top);

  //----- make other objects
  TGeoVolume* source = geom->MakeBox("source", Vacuum, 0.1, 0.1, 0.1);
  TGeoVolume* mask = geom->MakeBox("mask", Vacuum, 0.001, maskDimY, maskDimZ);
  TGeoVolume* detector =
      geom->MakeBox("detector", Vacuum, 0.001, detDimY, detDimZ);
  top->AddNode(source, 1, tr1);
  top->AddNode(mask, 1, tr2);
  top->AddNode(detector, 1, tr3);

  //----- colouring
  source->SetLineColor(kRed);
  source->SetFillColor(kRed);
  mask->SetLineColor(kMagenta);
  mask->SetFillColor(kMagenta);
  detector->SetLineColor(kOrange + 7);
  detector->SetFillColor(kOrange + 7);

  //----- close geometry and save
  geom->CloseGeometry();
  geom->SetVisLevel(4);
  geom->Export("results/" + fName + "_geometry.root");
}
//------------------------------------------------------------------
void CMSimulation::Print(void) {
  cout << "\nCMSimulation::Print() for object " << GetName() << endl;
}
//------------------------------------------------------------------
Bool_t CMSimulation::SetPattern(TH2F* h) {
  if (h) {
    fMask.SetPattern(h);
    return kTRUE;
  } else {
    cout << "No such pattern for mask available, please check..." << endl;
    return kFALSE;
  }
}
//------------------------------------------------------------------
