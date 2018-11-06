{
//////////////////////////////////////////////////////////
//   This file has been automatically generated 
//     (Mon Nov 21 14:45:57 2016 by ROOT version5.34/36)
//   from TTree data/data
//   found on file: CMsimulation.root
//////////////////////////////////////////////////////////

  gStyle->SetOptStat(111111);
  TH2F* hInit = new TH2F("init","init",100,-200,200,100,-200,200);
  hInit->GetXaxis()->SetTitle("z");
  hInit->GetYaxis()->SetTitle("y");
  TH2F* hMask = new TH2F("mask","mask",100,-200,200,100,-200,200);
  hMask->GetXaxis()->SetTitle("z");
  hMask->GetYaxis()->SetTitle("y");
  TH2F* hDet = new TH2F("det","det",100,-200,200,100,-200,200);
  hDet->GetXaxis()->SetTitle("z");
  hDet->GetYaxis()->SetTitle("y");

//Reset ROOT and connect tree file
//   gROOT->Reset();
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("results/CMsimulation.root");
   if (!f) {
      f = new TFile("results/CMsimulation.root");
   }
   TTree* tree = (TTree*)f->Get("data");

//Declaration of leaves types
   TVector3* point0 = new TVector3(0,0,0);
   TVector3* point1 = new TVector3(0,0,0);
   TVector3* point2 = new TVector3(0,0,0);
   TVector3* dir    = new TVector3(0,0,0);
   
   Int_t           opaque;

   // Set branch addresses.
   data->SetBranchAddress("point0",&point0);
   data->SetBranchAddress("point1",&point1);
   data->SetBranchAddress("point2",&point2);
   data->SetBranchAddress("dir",&dir);
   data->SetBranchAddress("opaque",&opaque);

   Int_t nentries = data->GetEntries();
   for (Long64_t ii=0; ii<nentries; ii++) {
     data->GetEntry(ii);
     cout<<hMask<<endl;
     hMask->Fill(point1->Z(), point1->Y());
     hDet->Fill(point2->Z(), point2->Y());
     hInit->Fill(point0->Z(), point0->Y());
   }
   TCanvas* can = new TCanvas();
   can->Divide(2,2);
   can->cd(1);
   hInit->Draw("col");
   can->cd(2);
   hMask->Draw("col");
   can->cd(3);
   hDet->Draw("col");
}
