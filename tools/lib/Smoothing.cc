#include "Smoothing.hh"

namespace SiFi {
namespace tools {

TH2F* SmoothGauss(TH2F* hin, double sigma){

  if(sigma <= 0){
    std::cout << "Smearing with sigma = " << sigma 
              << " will not work, provide a positive number here..." << std::endl;
    return nullptr;
  }
  
  TH2F* hout = dynamic_cast<TH2F*>(hin->Clone(Form("%s_smooth", hin->GetName())));
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
  TH2F* htmp = dynamic_cast<TH2F*>(hout->Clone());
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


} // namespace tools
} // namespace SiFi
