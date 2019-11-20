#include "Source.hh"

Source::Source() { CreateHistograms(); }

Source::Source(const TVector3& position, const Bool_t histos)
    : fPosition(position), fhTheta(nullptr), fhPhi(nullptr), fhZY(nullptr) {
  if (histos) CreateHistograms();
}

Source::Source(const TString& fname) : fInFileName(fname) {
  CreateHistograms();
}

Source::~Source() { DeleteHistograms(); }

void Source::CreateHistograms() {
  fhTheta = new TH1F("hThetaAll",
                     "theta angles of all generated particles; #theta / rad; ",
                     100, 0, TMath::Pi());
  fhTheta->SetLineColor(kRed);
  fhPhi =
      new TH1F("hPhiAll", "phi angles of all generated particles; #phi / rad; ",
               100, -TMath::Pi(), TMath::Pi());
  fhPhi->SetLineColor(kRed);
  fhZY =
      new TH2F("hZYAll",
               "Y vs Z for vertices of all generated particles; z / mm; y / mm",
               300, -300, 300, 300, -300, 300);
}

void Source::DeleteHistograms() {
  if (fhTheta) delete fhTheta;
  if (fhPhi) delete fhPhi;
  if (fhZY) delete fhZY;
};
