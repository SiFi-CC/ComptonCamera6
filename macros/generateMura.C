TH2F* generateMura1d(Int_t L = 101) {
  gROOT->Reset();
  TString title = Form("1d MURA mask %i", L);
  TH2F* h = new TH2F("hMURA1d", title, L, -0.5, L - 0.5, 1, -1, 1);
  h->SetStats(0);
  h->SetFillColor(kBlue);
  h->SetBinContent(1, 1, 0);
  for (Int_t i = 1; i < L; i++) {
    if (IsQuaResidue(i, L) == 1) h->SetBinContent(i, 1, 1);
  }
  h->Draw("col");
  TString fname = Form("masks/%s_%i.png", h->GetName(), L);
  gPad->SaveAs(fname);
  TString fname = Form("masks/%s_%i.root", h->GetName(), L);
  TFile* fout = new TFile(fname, "RECREATE");
  h->Write();
  fout->Close();
  return h;
}

TH2F* generateMura2d(Int_t L = 101) {
  gROOT->Reset();
  TString title = Form("2d MURA mask %i", L);
  TH2F* h = new TH2F("hMURA2d", title, L, -0.5, L - 0.5, L, -0.5, L - 0.5);
  h->SetStats(0);
  for (int bin = 1; bin < L + 1; bin++)
    h->SetBinContent(1, bin, 0);
  for (int bin = 2; bin < L + 1; bin++)
    h->SetBinContent(bin, 1, 1);
  int binx, biny;
  for (Int_t i = 1; i < L; i++) {
    for (Int_t j = 1; j < L; j++) {
      binx = i + 1;
      biny = j + 1;
      if (IsQuaResidue(i, L) * IsQuaResidue(j, L) == -1)
        h->SetBinContent(binx, biny, 1);
    }
  }
  h->Draw("col");
  TString fname = Form("masks/%s_%i.png", h->GetName(), L);
  gPad->SaveAs(fname);
  TString fname = Form("masks/%s_%i.root", h->GetName(), L);
  TFile* fout = new TFile(fname, "RECREATE");
  h->Write();
  fout->Close();
  return h;
}

Int_t IsQuaResidue(Int_t q, Int_t p) {
  Int_t result = -1;
  for (int i = 1; i < p; i++) {
    if ((i * i) % p == q) {
      result = 1;
      break;
    }
  }
  return result;
}

void generateMura(void) {
  gStyle->SetPalette(8, 0);
  int kL[23] = {7,  11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
                53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101};
  for (int i = 0; i < 23; i++) {
    generateMura1d(kL[i]);
    generateMura2d(kL[i]);
  }
}
