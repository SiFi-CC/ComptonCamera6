#include "CMReconstruction.hh"
#include "CMSimulation.hh"
#include "TCanvas.h"
#include "TF2.h"
#include "TMath.h"
#include "TRandom.h"
#include <fstream>
#include <iostream>
#include <spdlog/spdlog.h>
using namespace std;

ClassImp(CMReconstruction);
//------------------------------------------------------------------
CMReconstruction::CMReconstruction() {
  fFileIn = 0;
  fFileOut = 0;
  fImage = 0;
  fObject = 0;
  fVerbose = 0;
  fNiter = 0;
  for (int i = 0; i < 100; i++)
    fRecoObject[i] = 0;
}
//------------------------------------------------------------------
CMReconstruction::CMReconstruction(TString filename, Int_t vlevel) {
  fInputName = filename;
  filename.ReplaceAll(".root", "");
  SetName(filename);
  fVerbose = vlevel;
  fFileIn = new TFile(fInputName, "READ");
  fFileIn->Print();
  if (!fFileIn) {
    cout << "Input file " << fInputName << " not opened correctly..." << endl;
    // exit(0);
  }
  fMask = (Mask*)fFileIn->Get("mask");
  fDetPlane = (DetPlane*)fFileIn->Get("detector");
  fSource = (Source*)fFileIn->Get("source");
  fImage = (TH2F*)fFileIn->Get("detectedYZ");
  fObject = (TH2F*)fFileIn->Get("sourceYZ");
  fNbinsxO = fObject->GetXaxis()->GetNbins();
  fNbinsyO = fObject->GetYaxis()->GetNbins();
  fNvoxelsO = fNbinsxO * fNbinsyO;
  fNbinsxI = fImage->GetXaxis()->GetNbins();
  fNbinsyI = fImage->GetYaxis()->GetNbins();
  fNvoxelsI = fNbinsxI * fNbinsyI;

  fFileOut = new TFile(filename + "_reco.root", "RECREATE");
  fHmatrix = new TH2D("hHmatrix", "H matrix", fNvoxelsO, 0.5, fNvoxelsO + 0.5,
                      fNvoxelsI, 0.5, fNvoxelsI + 0.5);
  fHmatrix->GetXaxis()->SetTitle("Nr of Object voxel");
  fHmatrix->GetYaxis()->SetTitle("Nr of Image voxel");
  fNiter = 0;
  fThisIter = 0;
  for (int i = 0; i < 100; i++)
    fRecoObject[i] = 0;
}
//------------------------------------------------------------------
CMReconstruction::~CMReconstruction() {
  if (fVerbose) cout << "Inside destructor of CMsimulation class" << endl;
  Write();
}
//------------------------------------------------------------------
void CMReconstruction::Write(void) {
  if (fVerbose) cout << "Inside CMReconstruction::Write()..." << endl;
  if (fFileOut) {
    fFileOut->cd();
    fMask->GetPattern()->Write();
    fObject->Write();
    fImage->Write();
    fMask->GetPattern()->Write();
    fHmatrix->Write();
    int i = 0;
    while (fRecoObject[i] != 0) {
      fRecoObject[i]->Write();
      i++;
    }

    // TCanvas* a = new TCanvas("summary","summary");
    // a->Divide(3,2);
    // a->cd(1);
    // fObject->Draw("col");
    // a->cd(3);
    // fMask->GetPattern()->Draw();
    // a->cd(2);
    // fImage->Draw("colz");
    // a->cd(4);
    // fObject->Draw("colz");
    // a->Write();

    fFileOut->Close();
    fFileIn->Close();
  }
}
void CMReconstruction::Print(void) {
  cout << "\nCMReconstruction::Print() for object " << GetName() << endl;
}
//------------------------------------------------------------------
Bool_t CMReconstruction::FillHMatrix(void) {
  spdlog::set_level(spdlog::level::warn);
  CMSimulation* sim = new CMSimulation(fSource, fMask, fDetPlane);
  int bx, by, bz;
  double x, y, prob;
  TH2F* tmpimg = sim->GetImage();
  int nev = 10000;
  for (int i = 1; i <= fNvoxelsO; i++) { // loop over object voxels
    sim->ResetSimulation();
    SingleToDoubleIdx("O", i, bx, by);
    x = fObject->GetXaxis()->GetBinCenter(bx);
    y = fObject->GetYaxis()->GetBinCenter(by);
    // sim->SetSourcePosition(0, y, x);
    sim->RunSimulation(nev);
    for (int j = 1; j <= fNvoxelsI; j++) { // loop over image voxels
      SingleToDoubleIdx("I", j, bx, by);
      prob = tmpimg->GetBinContent(bx, by) / nev;
      fHmatrix->SetBinContent(i, j, prob);
    }
  }
  spdlog::set_level(spdlog::level::debug);
  // delete sim;
  return kTRUE;
}
//------------------------------------------------------------------
Bool_t CMReconstruction::SingleToDoubleIdx(TString which, int i, int& binx,
                                           int& biny) {
  int nbinsx, nbinsy;
  if (which == "I") {
    nbinsx = fNbinsxI;
    nbinsy = fNbinsyI;
  } else if (which == "O") {
    nbinsx = fNbinsxO;
    nbinsy = fNbinsyO;
  } else
    return kFALSE;

  binx = i % nbinsx;
  if (binx == 0) binx = nbinsx;
  biny = i / nbinsx + 1;
  if (biny > nbinsy) biny = nbinsy;
  return kTRUE;
}
//------------------------------------------------------------------
Bool_t CMReconstruction::DoubleToSingleIdx(TString which, int binx, int biny,
                                           int& i) {
  int nbinsx, nbinsy;
  if (which == "I") {
    nbinsx = fNbinsxI;
    nbinsy = fNbinsyI;
  } else if (which == "O") {
    nbinsx = fNbinsxO;
    nbinsy = fNbinsyO;
  } else
    return kFALSE;

  i = (biny - 1) * nbinsx + binx;
  return kTRUE;
}
//------------------------------------------------------------------
Double_t CMReconstruction::Image(Int_t i) {
  int bx, by;
  SingleToDoubleIdx("I", i, bx, by);
  return fImage->GetBinContent(bx, by);
}
//------------------------------------------------------------------
Double_t CMReconstruction::Object(Int_t i) {
  int bx, by;
  SingleToDoubleIdx("O", i, bx, by);
  return fObject->GetBinContent(bx, by);
}
//------------------------------------------------------------------
Double_t CMReconstruction::RecoObject(Int_t i) {
  int bx, by;
  SingleToDoubleIdx("O", i, bx, by);
  return fRecoObject[fThisIter - 1]->GetBinContent(bx, by);
}
//------------------------------------------------------------------
Double_t CMReconstruction::H(Int_t i, int j) {
  return fHmatrix->GetBinContent(i, j);
}
//------------------------------------------------------------------
Double_t CMReconstruction::Hprime(Int_t i, int j) {
  return fHmatrix->GetBinContent(j, i);
}
//------------------------------------------------------------------
Bool_t CMReconstruction::CalculateS(void) {
  // S = new Double_t[fNvoxelsO+1];

  S.reserve(fNvoxelsO + 1);
  S[0] = 0;
  for (int j = 1; j < fNvoxelsO + 1; j++) {
    S[j] = 0;
    for (int i = 1; i < fNvoxelsI + 1; i++) {
      S[j] += H(j, i);
    }
  }
  cout << fObject << endl;
  fObject->Print();
  fRecoObject[0] = (TH2F*)fObject->Clone("hYZreco00");
  fRecoObject[0]->SetTitle("reco YZ, iter 00");
  fRecoObject[0]->Reset();
  for (int i = 0; i < fNbinsxO; i++) {
    for (int j = 0; j < fNbinsyO; j++) {
      fRecoObject[0]->SetBinContent(i + 1, j + 1, 1);
    }
  }
  // TF2* f = new
  // TF2("gaus2d","exp(-(x*x)/(2.*[0]*[0]))*exp(-(y*y)/(2.*[1]*[1]))",-150,150,-150,150);
  // f->SetParameters(10.,10);
  // fRecoObject[0]->FillRandom("gaus2d",100000);
  // delete f;

  return kTRUE;
}
//------------------------------------------------------------------
Bool_t CMReconstruction::SingleIteration(void) {
  if (fThisIter == 0) {
    FillHMatrix();
    CalculateS();
  }
  fThisIter++;
  fRecoObject[fThisIter] =
      (TH2F*)fObject->Clone(Form("hYZreco%02i", fThisIter));
  fRecoObject[fThisIter]->SetTitle(Form("reco YZ, iter %02i", fThisIter));
  fRecoObject[fThisIter]->Reset();
  Double_t den, P, fj_new;
  Double_t R[fNvoxelsI + 1];
  Int_t bx, by;

  // calculating I/(H*fk)
  for (int i = 1; i < fNbinsxO + 1; i++) {
    R[i] = 0;
    den = 0;
    for (int l = 1; l < fNvoxelsO + 1; l++) // single element of denominator
      den += (H(l, i) * RecoObject(l));
    R[i] = Image(i) / den;
  }

  for (int j = 1; j < fNvoxelsO + 1; j++) {
    fj_new = 0;
    // P
    P = 0;
    for (int k = 1; k < fNvoxelsI + 1; k++)
      P += (Hprime(k, j) * R[k]);
    // full expression
    fj_new = RecoObject(j) * P / S[j];
    SingleToDoubleIdx("O", j, bx, by);
    fRecoObject[fThisIter]->SetBinContent(bx, by, fj_new);
  }
  return kTRUE;
}
//------------------------------------------------------------------
Bool_t CMReconstruction::MLEMIterate(Int_t ni) {
  fNiter = ni;
  if (fNiter > 100) {
    cout << "Too many iterations requested. Currently <100 feasible. "
         << "\nFor more please adjust the code." << endl;
    return kFALSE;
  }
  for (int i = 0; i < fNiter; i++) {
    cout << "Before " << i + 1 << "th iteration..." << endl;
    SingleIteration();
    cout << "\tdone!" << endl;
  }

  TH1D* hProZ[100];
  TH1D* hProY[100];

  fLogger->info("canvas created");
  TCanvas* can = new TCanvas("MLEM2D", "MLEM2D", 1000, 1000);
  TCanvas* canz = new TCanvas("MLEM1DZ", "MLEM1DZ", 1000, 1000);
  TCanvas* cany = new TCanvas("MLEM1DY", "MLEM1DY", 1000, 1000);

  can->Divide((int)sqrt(fNiter) + 1, (int)sqrt(fNiter) + 1);
  canz->Divide((int)sqrt(fNiter) + 1, (int)sqrt(fNiter) + 1);
  cany->Divide((int)sqrt(fNiter) + 1, (int)sqrt(fNiter) + 1);

  fLogger->info("save iterations");
  for (int iter = 0; iter < fNiter + 1; iter++) {
    fLogger->info("saving iteration %d", iter);
    can->cd(iter + 1);
    gPad->SetLogz(1);
    fRecoObject[iter]->Draw("colz");
    hProZ[iter] = fRecoObject[iter]->ProjectionX();
    hProY[iter] = fRecoObject[iter]->ProjectionY();
    canz->cd(iter + 1);
    hProZ[iter]->Draw();
    cany->cd(iter + 1);
    hProY[iter]->Draw();
  }
  fFileOut->cd();
  fLogger->info("canvas write");
  can->Write();
  fLogger->info("canvas z write");
  canz->Write();
  fLogger->info("canvas y write");
  cany->Write();
  return true;
}
