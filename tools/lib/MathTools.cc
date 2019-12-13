#include "MathTools.hh"

//------------------------------------------------------------------
void MathTools::ScaleGraphAndMove(TGraphErrors* g, Double_t factor, Double_t offset){
    
  Double_t x, y;
  for(Int_t i=0; i<g->GetN(); i++){
    g->GetPoint(i, x, y);
    g->SetPoint(i, x+offset, y*factor);
    g->SetPointError(i, g->GetErrorX(i), g->GetErrorY(i)*factor); // TODO do errors scale the same?
  }
  return;
  
}
//------------------------------------------------------------------
void MathTools::ScaleGraphAndMove(TGraph* g, Double_t factor, Double_t offset){
    
  Double_t x, y;
  for(Int_t i=0; i<g->GetN(); i++){
    g->GetPoint(i, x, y);
    g->SetPoint(i, x+offset, y*factor);
  }
  return;
  
}
//------------------------------------------------------------------
TGraphErrors* MathTools::SubtractBackground(TGraphErrors* gorig, TF1* bg){
 
  Double_t x, y;
  TString newname  = gorig->GetName();
  TString newtitle = gorig->GetTitle();
  newname  += "_bgs";
  newtitle += " bgs";
  TGraphErrors* g = (TGraphErrors*)gorig->Clone(newname);
  g->SetTitle(newtitle);
  
  for(Int_t i = 0; i<g->GetN(); i++){
    g->GetPoint(i, x, y);
    g->SetPoint(i, x, y-bg->Eval(x));
    g->SetPointError(i, g->GetErrorX(i), g->GetErrorY(i));
  }
  return g;
  
}
//------------------------------------------------------------------
TGraphErrors* MathTools::SubtractBackground(TGraphErrors* gorig, TSpline3* bg){
    
  Double_t x, y;
  TString newname  = gorig->GetName();
  TString newtitle = gorig->GetTitle();
  newname  += "_bgs";
  newtitle += " bgs";
  TGraphErrors* g = (TGraphErrors*)gorig->Clone(newname);
  g->SetTitle(newtitle);
  
  for(Int_t i=0; i<g->GetN(); i++){
    g->GetPoint(i, x, y);
    g->SetPoint(i, x, y-bg->Eval(x));
    g->SetPointError(i, g->GetErrorX(i), g->GetErrorY(i));
  }
  return g;
  
}
//------------------------------------------------------------------
TGraphErrors* MathTools::SubtractBackground(TGraphErrors* gorig, TGraphErrors* bg){
    
  Double_t x, y;
  TString newname  = gorig->GetName();
  TString newtitle = gorig->GetTitle();
  newname  += "_bgs";
  newtitle += " bgs";
  TGraphErrors* g = (TGraphErrors*)gorig->Clone(newname);
  g->SetTitle(newtitle);
  
  for(Int_t i=0; i<g->GetN(); i++){
    g->GetPoint(i, x, y);
    g->SetPoint(i, x, y-bg->Eval(x));
    g->SetPointError(i, g->GetErrorX(i), g->GetErrorY(i));
  }
  return g;
  
}
//------------------------------------------------------------------
/// The function implements gaussian smearing/blur in two passes, first x and then y,
/// based on the statements from https://en.wikipedia.org/wiki/Gaussian_blur
///  that this is equivalent to a single 2D pass, but requires less computing time.
/// \param hin is the pointer to the histogram which should be smeared
/// \param sigma is the standard deviation of the gaussian kernel
///
/// Function returns a new object of a class TH2, therefore to get the desired type
/// one needs to cast
/// Example of usage: TH2F* h = dynamic_cast<TH2F*>(SmearGauss(hin, 3.0));

TH2* MathTools::SmoothGauss(TH2* hin, double sigma){

  if(sigma <= 0){
    std::cout << "Smearing with sigma = " << sigma 
              << " will not work, provide a positive number here..." << std::endl;
    return nullptr;
  }
  
  TH2* hout = dynamic_cast<TH2*>(hin->Clone(Form("%s_smooth", hin->GetName())));
  hout->Reset();
  const int nbinsx = hin->GetNbinsX();
  const int nbinsy = hin->GetNbinsY();
  double binwx = hin->GetXaxis()->GetBinWidth(1);
  double binwy = hin->GetYaxis()->GetBinWidth(1);
  
  double kernelx[nbinsx];
  double kernely[nbinsy];
  
  TF1* gaus = new TF1("gaus", "gaus");
  gaus->SetParameters(1./sqrt(TMath::TwoPi())/sigma, 0, sigma);

  for(int i=0; i<nbinsx; i++)
    kernelx[i] = gaus->Eval(binwx*i);
  for(int i=0; i<nbinsy; i++)
    kernely[i] = gaus->Eval(binwy*i);

  int deltabin = 0;
  double z = 0;

  //smearing in rows
  for(int biny=1; biny<nbinsy; biny++){
    for(int binx=1; binx<nbinsx; binx++){
      z = 0;
      for(int binxp=1; binxp<nbinsx; binxp++){
        deltabin = abs(binxp-binx);
        z += kernelx[deltabin]*hin->GetBinContent(binxp, biny);
      }
      hout->SetBinContent(binx, biny, z);      
    }
  }
  TH2* htmp = dynamic_cast<TH2*>(hout->Clone());
  hout->Reset();
  
  //smearing in columns
  for(int binx=1; binx<nbinsx; binx++){
    for(int biny=1; biny<nbinsy; biny++){
      z = 0;
      for(int binyp=1; binyp<nbinsy; binyp++){
        deltabin = abs(binyp-biny);
        z += kernely[deltabin]*htmp->GetBinContent(binx, binyp);	
      }
      hout->SetBinContent(binx, biny, z);      
    }
  }

  return hout;
}
//------------------------------------------------------------------
