#include "CMReconstruction.hh"
#include "CMSimulation.hh"
#include "TRandom.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TF2.h"
#include <iostream>
#include <fstream>
using namespace std;

ClassImp(CMReconstruction);
//------------------------------------------------------------------
CMReconstruction::CMReconstruction(){
  fFileIn=0;
  fFileOut=0;
  fImage = 0;
  fObject = 0;
  fVerbose = 0;
  fNiter = 0;
  for(int i=0; i<100; i++)
    fRecoObject[i] = 0;
}
//------------------------------------------------------------------
CMReconstruction::CMReconstruction(TString filename, Int_t vlevel){
  fInputName = filename;
  filename.ReplaceAll(".root","");
  SetName(filename);
  fVerbose = vlevel;
  fFileIn = new TFile("results/"+fInputName,"READ");
  fFileIn->Print();
  if(!fFileIn){
    cout<<"Input file "<<fInputName<<" not opened correctly..."<<endl;
    //exit(0);
  }
  fImage = (TH2F*)fFileIn->Get("hYZdetected");
  fObject = (TH2F*)fFileIn->Get("hYZ");
  fNbinsxO = fObject->GetXaxis()->GetNbins();
  fNbinsyO = fObject->GetYaxis()->GetNbins();
  fNvoxelsO = fNbinsxO * fNbinsyO;
  fNbinsxI = fImage->GetXaxis()->GetNbins();
  fNbinsyI = fImage->GetYaxis()->GetNbins();
  fNvoxelsI = fNbinsxI * fNbinsyI;
  RebuildSetupTxt();
  fMask.SetPattern((TH2F*)fFileIn->Get("mask"));
  fFileOut = new TFile("results/"+filename+"_reco.root","RECREATE");
  fHmatrix = new TH2D("hHmatrix","H matrix",
		     fNvoxelsO, 0.5, fNvoxelsO+0.5,
		     fNvoxelsI, 0.5, fNvoxelsI+0.5);
  fHmatrix->GetXaxis()->SetTitle("Nr of Object voxel");
  fHmatrix->GetYaxis()->SetTitle("Nr of Image voxel");
  fNiter = 0;
  fThisIter = 0;
  for(int i=0; i<100; i++)
    fRecoObject[i] = 0;
}
//------------------------------------------------------------------
CMReconstruction::~CMReconstruction(){
  if(fVerbose) cout<<"Inside destructor of CMsimulation class"<<endl;
  Write();
}
//------------------------------------------------------------------
void CMReconstruction::Write(void){
  if(fVerbose) cout<<"Inside CMReconstruction::Write()..."<<endl;
  if(fFileOut){
    fFileOut->cd();
    fMask.GetPattern()->Write();
    fObject->Write();
    fImage->Write();
    fMask.GetPattern()->Write();
    fHmatrix->Write();
    int i=0;
    while(fRecoObject[i] != 0){
      fRecoObject[i]->Write();
      i++;
    }

    // TCanvas* a = new TCanvas("summary","summary");
    // a->Divide(3,2);
    // a->cd(1);
    // fObject->Draw("col");
    // a->cd(3);
    // fMask.GetPattern()->Draw();
    // a->cd(2);
    // fImage->Draw("colz");
    // a->cd(4);
    // fObject->Draw("colz");
    // a->Write();

    fFileOut->Close();
    fFileIn->Close();
  }
}
//------------------------------------------------------------------
void CMReconstruction::RebuildSetupTxt(void){
  
  if(fVerbose) cout << "\n----- Rebuilding setup from the txt file \n" << endl;
  TString fname = fFileIn->GetName();
  fname.ReplaceAll(".root","_geometry.txt");
  ifstream input(fname.Data(), std::ios::in);
  if(!(input.is_open())){
    cout << "##### Could not open "<<fname<<" file! " << endl;
    cout << "##### Please check!" << endl;
    return;
  }

  Double_t maskPar[4];
  Double_t absPar[4];
  Double_t maskDim[2];
  Double_t absDim[2];
  TString maskName, absName;
  TString dummy;

  input >> dummy >> dummy >> dummy;
  input >> dummy >> dummy >> dummy;

  if(fVerbose) cout << "\n----- Rebuilding the mask" << endl;
  input >> dummy >> dummy >> maskName;
  input >> dummy >> dummy >> maskPar[0];
  input >> dummy >> dummy >> maskPar[1];
  input >> dummy >> dummy >> maskPar[2];
  input >> dummy >> dummy >> maskPar[3];
  input >> dummy >> dummy >> maskDim[0];
  input >> dummy >> dummy >> maskDim[1];

  fMask.SetPlane(maskPar[0],maskPar[1],maskPar[2],maskPar[3]);
  fMask.SetDimensions(maskDim[0],maskDim[1]);
  fMask.SetName(maskName);
  if(fVerbose) fMask.Print();

  if(fVerbose) cout << "\n----- Rebuilding the scatterer" << endl;
  input >> dummy >> dummy >> absName;
  input >> dummy >> dummy >> absPar[0];
  input >> dummy >> dummy >> absPar[1];
  input >> dummy >> dummy >> absPar[2];
  input >> dummy >> dummy >> absPar[3];
  input >> dummy >> dummy >> absDim[0];
  input >> dummy >> dummy >> absDim[1];

  fDetPlane.SetPlane(absPar[0],absPar[1],absPar[2],absPar[3]);
  fDetPlane.SetDimensions(absDim[0],absDim[1]);
  fDetPlane.SetName(absName);
  if(fVerbose) fDetPlane.Print();

  input.close();
}
//------------------------------------------------------------------
void CMReconstruction::SetupSpectra(void){

  double maskZdim = fMask.GetDimZ()/2;
  double maskYdim = fMask.GetDimY()/2;
  double detZdim = fDetPlane.GetDimZ()/2;
  double detYdim = fDetPlane.GetDimY()/2;
  int nbinsz = fMask.GetPattern()->GetXaxis()->GetNbins();
  int nbinsy = fMask.GetPattern()->GetYaxis()->GetNbins();
}
//------------------------------------------------------------------

void CMReconstruction::Print(void){
  cout<<"\nCMReconstruction::Print() for object "<<GetName()<<endl;
}
//------------------------------------------------------------------
Bool_t CMReconstruction::FillHMatrix(void){
  CMSimulation* sim = new CMSimulation("tmpsim", 0);
  sim->BuildSetup(fDetPlane, fMask);
  sim->SetGenVersion(1);
  sim->SetupSpectra();
  int bx,by,bz;
  double x,y, prob;
  TH2F* tmpimg = sim->GetImage();
  int nev = 10000;
  for(int i=1; i<=fNvoxelsO; i++){ //loop over object voxels
    sim->ClearSpectra();
    SingleToDoubleIdx("O",i, bx,by);
    x=fObject->GetXaxis()->GetBinCenter(bx);
    y=fObject->GetYaxis()->GetBinCenter(by);
    sim->SetSourcePosition(0, y, x);
    sim->Loop(nev);
    for(int j=1; j<=fNvoxelsI; j++){//loop over image voxels
      SingleToDoubleIdx("I",j, bx,by);
      prob = tmpimg->GetBinContent(bx,by)/nev;
      fHmatrix->SetBinContent(i,j,prob);
    }
  }
  //delete sim;
  return kTRUE;
}
//------------------------------------------------------------------
Bool_t CMReconstruction::SingleToDoubleIdx(TString which, int i, 
					   int &binx, int& biny){
  int nbinsx, nbinsy;
  if(which=="I"){
    nbinsx=fNbinsxI; 
    nbinsy=fNbinsyI;
  }
  else if(which=="O"){
    nbinsx=fNbinsxO; 
    nbinsy=fNbinsyO;
  }
  else return kFALSE;
  
  binx = i%nbinsx;
  if(binx==0) binx=nbinsx;
  biny = i/nbinsx + 1;
  if(biny>nbinsy) biny=nbinsy;
  return kTRUE;
}
//------------------------------------------------------------------
Bool_t CMReconstruction::DoubleToSingleIdx(TString which, 
					   int binx, int biny, int &i){
  int nbinsx, nbinsy;
  if(which=="I"){
    nbinsx=fNbinsxI; 
    nbinsy=fNbinsyI;
  }
  else if(which=="O"){
    nbinsx=fNbinsxO; 
    nbinsy=fNbinsyO;
  }
  else return kFALSE;

  i=(biny-1)*nbinsx+binx;
  return kTRUE;
}
//------------------------------------------------------------------
Double_t CMReconstruction::Image(Int_t i){
  int bx, by;
  SingleToDoubleIdx("I",i, bx, by);
  return fImage->GetBinContent(bx,by);
}
//------------------------------------------------------------------
Double_t CMReconstruction::Object(Int_t i){
  int bx, by;
  SingleToDoubleIdx("O",i, bx, by);
  return fObject->GetBinContent(bx,by);
}
//------------------------------------------------------------------
Double_t CMReconstruction::RecoObject(Int_t i){
  int bx, by;
  SingleToDoubleIdx("O",i, bx, by);
  return fRecoObject[fThisIter-1]->GetBinContent(bx,by);
}
//------------------------------------------------------------------
Double_t CMReconstruction::H(Int_t i, int j){
  return fHmatrix->GetBinContent(i,j);
}
//------------------------------------------------------------------
Double_t CMReconstruction::Hprime(Int_t i, int j){
  return fHmatrix->GetBinContent(j,i);
}
//------------------------------------------------------------------
Bool_t CMReconstruction::CalculateS(void){
  S = new Double_t[fNvoxelsO+1];
  
  for(int j=1; j<fNvoxelsO+1; j++){
    S[j] = 0;
    for(int i=1; i<fNvoxelsI+1; i++){
      S[j] += H(j,i);
    }
  }
  cout<<fObject<<endl;
  fObject->Print();
  fRecoObject[0] = (TH2F*)fObject->Clone("hYZreco00");
  fRecoObject[0]->SetTitle("reco YZ, iter 00");
  fRecoObject[0]->Reset();
  for(int i=0; i<fNbinsxO; i++){
    for(int j=0; j<fNbinsyO; j++){
      fRecoObject[0]->SetBinContent(i+1,j+1,1);
    }
  }
  TF2* f = new TF2("gaus2d","exp(-(x*x)/(2.*[0]*[0]))",-150,150,-150,150);
  f->SetParameters(10.,0.001);
  fRecoObject[0]->FillRandom("gaus2d",100000);
  delete f;
  
  return kTRUE;
}
//------------------------------------------------------------------
Bool_t CMReconstruction::SingleIteration(void){
  if(fThisIter==0){
    FillHMatrix();
    CalculateS();
  }
  fThisIter++;
  fRecoObject[fThisIter] = (TH2F*)fObject->Clone(Form("hYZreco%02i",fThisIter));
  fRecoObject[fThisIter]->SetTitle(Form("reco YZ, iter %02i",fThisIter));
  fRecoObject[fThisIter]->Reset();
  Double_t den, P, fj_new;
  Double_t R[fNvoxelsI+1];
  Int_t bx, by;

  //calculating I/(H*fk)
  for(int i=1; i<fNbinsxO+1; i++){
    R[i] = 0;
    den = 0;
    for(int l=1; l<fNvoxelsO+1; l++) //single element of denominator
      den+= (H(l,i)*RecoObject(l));
    R[i] = Image(i)/den;
  }

  for(int j=1; j<fNvoxelsO+1; j++){
    fj_new = 0;
    //P
    P = 0;
    for(int k=1; k<fNvoxelsI+1; k++)
      P += (Hprime(k,j)*R[k]);
    //full expression
    fj_new = RecoObject(j)*P/S[j];
    SingleToDoubleIdx("O",j, bx, by);
    fRecoObject[fThisIter]->SetBinContent(bx,by,fj_new);
  }
  return kTRUE;
}
//------------------------------------------------------------------
Bool_t CMReconstruction::MLEMIterate(Int_t ni){
  fNiter = ni;
  if(fNiter>100){
    cout<<"Too many iterations requested. Currently <100 feasible. "<<
      "\nFor more please adjust the code."<<endl;
    return kFALSE;
  }
  for(int i=0; i<fNiter; i++){
    cout<<"Before "<<i+1<<"th iteration..."<<endl;
    SingleIteration();
    cout<<"\tdone!"<<endl;
  }
}
