#include "PointSource.hh"
#include "CLog.hh"
#include "TCollection.h"
#include "TObjArray.h"
#include "TRandom.h"
#include "TObjString.h"
#include <iostream>
#include <fstream>

// -----------------------
// PointSource
// -----------------------

PointSource::PointSource(const TString fname) {
    SetName("PointSource");
    fInFileName = fname;    
    if(!Init()){
       TString message;
       message.Form("Exception in %s constructor!", this->ClassName());
       throw message.Data();
    }
    SetName("PointSource");
  }

Track PointSource::GenerateEvent() {
  Double_t theta = gRandom->Uniform(fMinAngle, fMaxAngle);
  TVector3 versor;
  versor.SetMagThetaPhi(1.,theta, gRandom->Uniform(0, TMath::TwoPi()));
  versor.RotateY(-TMath::Pi()/2);
  return Track(fPosition, versor, fEnergy);
}

Bool_t PointSource::Init(){
  if(!fInFileName.EndsWith(".mac")){
    std::cout<<"PointSource::Init(): not appropriate input file format, .mac file expected..."<<std::endl;
    return kFALSE;
  }
  std::ifstream infile;
  infile.open(fInFileName.Data(), std::ios::in); 
  if(!infile.is_open() || infile.fail() || infile.bad()) {
    std::cout<<"Input file "<<fInFileName<<" not opened correctly or corrupted"<<std::endl;
    std::cout<<infile.is_open() <<" "<< infile.fail() <<" "<< infile.bad()<<std::endl;
    return kFALSE;
  }
  
  TString line;
  while(!infile.eof()){
    line.ReadLine(infile);
    if(line.BeginsWith("#") || line.IsWhitespace())
      continue;
    
    TObjArray* words = line.Tokenize(" ");
    if(words->GetEntries()==0) continue;
    TString sitem;
    Double_t tmp;
    sitem = ((TObjString*)words->At(0))->GetString();
    if(sitem.Contains("/gps/pos/type")){
      sitem = ((TObjString*)words->At(1))->GetString();
      if(!sitem.Contains("Point")){
	std::cout<<"PointSource::Init(): Inconsistency in source types. "<<
	  "Point expected, "<<sitem<<" found in source config file "<<
	  fInFileName<<std::endl;
      }
    }
    else if(sitem.Contains("/gps/pos/centre")){
      Double_t pos[3];
      for(int i=0; i<3; i++){   
	pos[i] = ((TObjString*)words->At(i+1))->GetString().Atof();
      }
      SetPosition(TVector3(pos));
      sitem = ((TObjString*)words->At(4))->GetString();
      if(sitem.Contains("cm"))
	fPosition.SetMag(fPosition.Mag()*10.); //default units in sim. are mm
    }
    else if(sitem.Contains("/gps/ang/mintheta")){
      tmp =  ((TObjString*)words->At(1))->GetString().Atof();
      sitem = ((TObjString*)words->At(2))->GetString();
      if(sitem=="deg") tmp*=TMath::DegToRad();
      fMaxAngle = TMath::Pi()-tmp;
    }
    else if(sitem.Contains("/gps/ang/maxtheta")){
      tmp =  ((TObjString*)words->At(1))->GetString().Atof();
      sitem = ((TObjString*)words->At(2))->GetString();
      if(sitem=="deg") tmp*=TMath::DegToRad();
      fMinAngle = TMath::Pi()-tmp;
    }
    else if(sitem.Contains("/gps/ene/type")){
      sitem = ((TObjString*)words->At(1))->GetString();
      if(!sitem.Contains("Mono")){
	std::cout<<"PointSource::Init(): /gps/ene/type "<<sitem<<" not implemented, only Mono available"<<std::endl;
      }
    }
    else if(sitem.Contains("/gps/energy")){
      fEnergy = ((TObjString*)words->At(1))->GetString().Atof();
      sitem = ((TObjString*)words->At(2))->GetString();
      if(sitem=="keV") fEnergy*=0.001;
      else if(sitem=="GeV") fEnergy*=1000;
    }
    else{
      std::cout<<"Command "<<sitem<<" unknown, thus ignored..."<<std::endl;
    }
    
    delete words;
  }
  return kTRUE;
}

void PointSource::Print(){
  std::cout<<"\n\nPointSource::Print():\nName =\t"<<GetName()<<"\nEnergy =\t"<<fEnergy<<" MeV"<<"\nPosition=\t("<<fPosition.X()<<" , "<<fPosition.Y()<<" , "<<fPosition.Z()<<")\nAngularRange =\t("<<fMinAngle<<" , "<<fMaxAngle<<") rad\n"<<std::endl; 
}
