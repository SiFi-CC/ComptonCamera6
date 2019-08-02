#include "Source.hh"


void Source::CreateHistograms(){
    fhTheta = new TH1F("hThetaAll", "theta angles of all generated particles; #theta / rad; ", 100, 0,TMath::Pi());
    fhTheta->SetLineColor(kRed);
    fhPhi =   new TH1F("hPhiAll", "phi angles of all generated particles; #phi / rad; ", 100, -TMath::Pi(), TMath::Pi());
    fhPhi->SetLineColor(kRed);
    fhZY = new TH2F("hZYAll", "Y vs Z for vertices of all generated particles; z / mm; y / mm", 300, -300, 300,300, -300, 300);
}
