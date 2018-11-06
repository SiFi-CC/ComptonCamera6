Bool_t DrawData(void){
 
  //Setup histograms
  //gStyle->SetOptStat(111111);
  gStyle->SetOptStat(0);
  Double_t size = 300.;
  TH2F* hStart = new TH2F("start","start",100,-size,size,100,-size,size);
  hStart->GetXaxis()->SetTitle("z");
  hStart->GetYaxis()->SetTitle("y");
  TH2F* hScat = new TH2F("scat","scat",100,-size,size,100,-size,size);
  hScat->GetXaxis()->SetTitle("z");
  hScat->GetYaxis()->SetTitle("y");
  TH2F* hAbs = new TH2F("abs","abs",100,-size,size,100,-size,size);
  hAbs->GetXaxis()->SetTitle("z");
  hAbs->GetYaxis()->SetTitle("y");
  TH1F* hEnergy = new TH1F("energy","energy",100,0,5);
  hEnergy->GetYaxis()->SetTitle(" ");
  hEnergy->GetXaxis()->SetTitle("energy [MeV]");
  
  
  //Connect tree file
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../results/CCSimulation_gen4.root");
  if (!f) {
     f = new TFile("../results/CCSimulation_gen4.root");
  }
  TTree* tree = (TTree*)f->Get("data");
  
  //Declaration of leaves types
  TVector3* point0 = new TVector3(0,0,0);
  TVector3* point1 = new TVector3(0,0,0);
  TVector3* point2 = new TVector3(0,0,0);
  TVector3* versor1 = new TVector3(0,0,0);
  TVector3* versor2 = new TVector3(0,0,0);
  Double_t energy0;
  Double_t energy1;
  Double_t energy2;
  
  //Set branch addresses
  data->SetBranchAddress("point0",&point0);
  data->SetBranchAddress("point1",&point1);
  data->SetBranchAddress("point2",&point2);
  data->SetBranchAddress("versor1",&versor1);
  data->SetBranchAddress("versor2",&versor2);
  data->SetBranchAddress("energy0",&energy0);
  data->SetBranchAddress("energy1",&energy1);
  data->SetBranchAddress("energy2",&energy2);
  
  //Fill histograms
  Int_t nentries = data->GetEntries();
  for(Int_t ii=0; ii<nentries; ii++) {
     data->GetEntry(ii);
     hStart->Fill(point0->Z(), point0->Y());
     hScat->Fill(point1->Z(), point1->Y());
     hAbs->Fill(point2->Z(), point2->Y());
     hEnergy->Fill(energy2);
   }
  
  //Setup canvas & draw
  TCanvas *can = new TCanvas("can","can",1200,1200);
  can->Divide(2,2);
  can->cd(1);
  //  hStart->SetStats(0);
  hStart->Draw("colz");
  can->cd(2);
  //hScat->SetStats(0);
  hScat->Draw("colz");
  can->cd(3);
  //hAbs->SetStats(0);
  hAbs->Draw("colz");
  can->cd(4);
  //hEnergy->SetStats(0);
  hEnergy->Draw();
  
  return kTRUE; 
}
