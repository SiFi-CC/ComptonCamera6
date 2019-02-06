#include "CMSimulation.hh"
#include "CLog.hh"
#include <TCanvas.h>
#include <TH2F.h>
#include <TMath.h>
#include <TRandom.h>

CMSimulation::~CMSimulation() {
  log->info("delete CMSimulation");
  log->debug("delete fTree");
  delete fTree;
  log->debug("delete fH2Source");
  delete fH2Source;
  log->debug("delete fH2Detector");
  delete fH2Detector;
  log->debug("delete fH1Theta");
  delete fH1Theta;
  log->debug("delete fH1Phi");
  delete fH1Phi;
}

void CMSimulation::Init() {
  TDirectory::AddDirectory(false);
  fTree = new TTree("track", "tracks");
  fTree->SetDirectory(nullptr);
  fTree->Branch("sourceTrack", &fPersist.sourceTrack);
  fTree->Branch("maskTrack", &fPersist.maskTrack);
  fTree->Branch("detectorTrack", &fPersist.detectorTrack);
  fTree->Branch("absorbedInMask", &fPersist.absorbed);

  double maskZdim = fMask->GetDimZ() / 2;
  double maskYdim = fMask->GetDimY() / 2;
  int nbinsz = fMask->GetPattern()->GetXaxis()->GetNbins();
  int nbinsy = fMask->GetPattern()->GetYaxis()->GetNbins();
  fH2Source = new TH2F("sourceYZ", "Y vs Z source", nbinsz, -maskZdim, maskZdim,
                       nbinsy, -maskYdim, maskYdim);
  fH2Detector = new TH2F("detectedYZ", "Y vs Z detected", nbinsz, -maskZdim,
                         maskZdim, nbinsy, -maskYdim, maskYdim);
  fH1Theta = new TH1F("hTheta", "theta angle of generated particles", 100, 0,
                      TMath::Pi());
  fH1Phi = new TH1F("hPhi", "phi angle of generated particles", 100,
                    -TMath::Pi(), TMath::Pi());
}

Bool_t CMSimulation::ProcessEvent() {
  Track sourceTrack = fSource->GenerateEvent();

  auto maskCross = fMask->FindCrossPoint(sourceTrack);
  if (!maskCross.second) {
    log->debug("No cross point with Mask");
    return false;
  }
  auto detCross = fDetPlane->FindCrossPoint(sourceTrack);
  if (!detCross.second) {
    log->debug("No cross point with DetPlane");
    return false;
  }

  // commit data
  fPersist.sourceTrack = sourceTrack;
  fPersist.maskTrack =
      Track(maskCross.first, sourceTrack.GetVersor(), sourceTrack.GetEnergy());
  fPersist.detectorTrack =
      Track(detCross.first, sourceTrack.GetVersor(), sourceTrack.GetEnergy());
  fPersist.absorbed = !fMask->IsOpaque(maskCross.first);
  fTree->Fill();

  fH2Source->Fill(sourceTrack.GetPoint().Z(), sourceTrack.GetPoint().Y());
  fH1Theta->Fill(sourceTrack.GetVersor().Theta());
  fH1Phi->Fill(sourceTrack.GetVersor().Phi());
  if (!fPersist.absorbed) {
    fH2Detector->Fill(detCross.first.Z(), detCross.first.Y());
  }
  return true;
}

void CMSimulation::RunSimulation(Int_t nEvents) {
  Int_t eventsProcessed = 0;
  Int_t eventsAccepted = 0;
  while (nEvents > eventsAccepted) {
    if (ProcessEvent()) { eventsAccepted++; }
    eventsProcessed++;
    log->debug("{} events simulated, {} in acceptance", eventsProcessed,
               eventsAccepted);
    if (eventsProcessed % 100000 == 0) {
      log->info("{} processed events", eventsProcessed);
    }
  }
  log->info("Finished simulation: {} events simulated, {} in acceptance",
            eventsProcessed, eventsAccepted);
}

void CMSimulation::Write(TString name) const {
  log->info("Saving results of simulation to file");

  log->debug("Save raw data");
  TFile file(name, "RECREATE");
  fTree->SetDirectory(&file);
  fTree->Write();
  fTree->SetDirectory(nullptr);

  log->debug("Save histograms");
  fH2Source->Write("H2Source");
  fH2Detector->Write("H2Detector");
  fH1Theta->Write("H1Theta");
  fH1Phi->Write("H1Phi");

  log->debug("Save simulation configuration (source, mask, detector)");
  fSource->Write("source");
  fMask->Write("mask");
  fDetPlane->Write("detector");

  TH1D* h1 = fH2Source->ProjectionX("sourceProjectionZ");
  h1->SetTitle("source along Z axis");
  h1->SetFillColor(kRed);
  h1->SetLineColor(kRed);
  h1->SetFillStyle(3004);
  h1->Write();
  TH1D* h2 = fH2Detector->ProjectionX("detectorProjectionZ");
  h2->SetTitle("detector along Z axis");
  h2->SetFillColor(kBlue);
  h2->SetLineColor(kBlue);
  h2->SetFillStyle(3005);
  h2->Write();

  log->debug("Prepare simulation summary");
  TCanvas canvas("summary", "summary");
  canvas.Divide(3, 2);

  canvas.cd(1);
  fH1Theta->Draw();

  canvas.cd(4);
  fH1Phi->Draw();

  canvas.cd(2);
  fH2Source->Draw("COLZ");

  canvas.cd(5);
  fH2Detector->Draw("COLZ");

  canvas.cd(3);
  h1->Draw();

  canvas.cd(6);
  h2->Draw();
  canvas.Write();

  log->debug("Closing file");
  file.Close();

  // BuildTGeometry(name);
}

void CMSimulation::ResetSimulation() {
  fH2Source->Reset();
  fH2Detector->Reset();
  fH1Theta->Reset();
  fH1Phi->Reset();
  fTree->Reset();
}

void CMSimulation::Print() const {
  log->info("Mask object: {}\n"
            "fA = {}\nfB = {}\nfC = {}\nfD = {}\n"
            "fDimZ = {}\nfDimY = {}\n",
            "DetPlane object: {}\n", "fA = {}\nfB = {}\nfC = {}\nfD = {}\n"
                                     "fDimZ = {}\nfDimY = {}\n",
            fMask->GetName(), fMask->GetA(), fMask->GetB(), fMask->GetC(),
            fMask->GetD(), fMask->GetDimZ(), fMask->GetDimY(),
            fDetPlane->GetName(), fDetPlane->GetA(), fDetPlane->GetB(),
            fDetPlane->GetC(), fDetPlane->GetD(), fDetPlane->GetDimZ(),
            fDetPlane->GetDimY());

  log->info("Geometry of the setup saved in the text file");
}

void CMSimulation::BuildTGeometry(TString name) const {
  TGeoManager geom("world", "CMSimulation - world");

  //----- define materials
  TGeoMaterial matVacuum("Vacuum", 0, 0, 0);

  //----- define media
  TGeoMedium Vacuum("Vacuum", 1, &matVacuum);

  //----- define translations
  Double_t maskA = fMask->GetA();
  Double_t maskD = fMask->GetD();
  Double_t detA = fDetPlane->GetA();
  Double_t detD = fDetPlane->GetD();
  Double_t maskTrans = -maskD / maskA;
  Double_t detTrans = -detD / detA;
  Double_t x = 0;
  Double_t y = 0;
  Double_t z = 0;

  TGeoTranslation tr1(x, y, z);
  TGeoTranslation tr2(maskTrans, 0., 0.);
  TGeoTranslation tr3(detTrans, 0., 0.);

  //----- define dimensions of objects
  Double_t maskDimZ = fMask->GetDimZ() / 2.; // half length
  Double_t maskDimY = fMask->GetDimY() / 2.;
  Double_t detDimZ = fDetPlane->GetDimZ() / 2.;
  Double_t detDimY = fDetPlane->GetDimY() / 2.;

  //----- make the top container volume
  Double_t worldX = -detTrans + 50; // cm hafl length
  Double_t worldY = detDimY + 50.;
  Double_t worldZ = detDimZ + 50;
  TGeoVolume* top = geom.MakeBox("top", &Vacuum, worldX, worldY, worldZ);
  geom.SetTopVolume(top);

  //----- make other objects
  TGeoVolume* source = geom.MakeBox("source", &Vacuum, 0.1, 0.1, 0.1);
  TGeoVolume* mask = geom.MakeBox("mask", &Vacuum, 0.001, maskDimY, maskDimZ);
  TGeoVolume* detector =
      geom.MakeBox("detector", &Vacuum, 0.001, detDimY, detDimZ);
  top->AddNode(source, 1, &tr1);
  top->AddNode(mask, 1, &tr2);
  top->AddNode(detector, 1, &tr3);

  //----- colouring
  source->SetLineColor(kRed);
  source->SetFillColor(kRed);
  mask->SetLineColor(kMagenta);
  mask->SetFillColor(kMagenta);
  detector->SetLineColor(kOrange + 7);
  detector->SetFillColor(kOrange + 7);

  //----- close geometry and save
  geom.CloseGeometry();
  geom.SetVisLevel(4);
  geom.Export("results/" + name + "_geometry.root");
}
