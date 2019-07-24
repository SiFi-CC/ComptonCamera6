#include "PointSource.hh"
#include "CLog.hh"
#include <TRandom.h>
#include <iostream>
#include <fstream>

// -----------------------
// PointSource
// -----------------------

PointSource(const TString& fname) : fInFileName(fname){
  if(!Init())
    throw "Exception in PointSource constructor!";
};

Track PointSource::GenerateEvent() {
  Double_t angleY = gRandom->Uniform(-fAngleY, fAngleY);
  Double_t angleZ = gRandom->Uniform(-fAngleZ, fAngleZ);
  TVector3 versor = TVector3(-1, tan(angleY), tan(angleZ)).Unit();
  return Track(fPosition, versor, fEnergy);
}

Bool_t PointSource::Init(){
  if(!fInFileName.EndsWith(".mac")){
    std::cout<<"PointSource::Init(): not appropriate input file format, .mac file expected..."<<std::endl;
    return kFALSE;
  }
  ifstream infile();
  infile.open(fInFileName.Data(), ios::in); 
  if(!infile.is_open() || infile.fail() || infile.bad()) {
    cout<<"Input file "<<fMcaFile<<" not opened correctly or corrupted"<<endl;
    cout<<infile.is_open() <<" "<< infile.fail() <<" "<< infile.bad()<<endl;
    return kFALSE;
  }
  
  TString line;
  while(!infile.eof()){
    line.ReadLine(infile);
    if(line.BeginsWith("#"))
      continue;
    
    TObjArray* words = mcafile.Tokenize("/");
    TIter iWord(words);
    TString sitem;
    Double_t tmp;
    TObjString* os = iWord.Begin();
    sitem = ((TObjString*)iWord())->GetString();
    if(sitem.Contains("/gps/pos/type")){
      iWord++;
      sitem = ((TObjString*)iWord())->GetString();
      if(!sitem.Contains("Point")){
	std::cout<<"PointSource::Init(): Inconsistency in source types. "<<
	  "Point expected, "<<sitem<<" found in source config file "<<
	  fInFileName<<std::endl;
      }
    }
    else if(sitem.Contains("/gps/pos/centre")){
      Double_t pos[3];
      for(int i=0; i<3; i++){
	iWord++;
	pos[i] = ((TObjString*)iWord())->GetString()->Atof();
      }
      SetPosition(TVector3(pos));
      iWord++;
      sitem = ((TObjString*)iWord())->GetString();
      if(sitem.Contains("cm"))
	fPosition.SetMag(fPosition.Mag()*10.); //default units in sim. are mm
    }
    else if(sitem.Contains("/gps/pos/mintheta")){
      iWord++;
      tmp =  ((TObjString*)iWord())->GetString()->Atof();
      iWord++;
      sitem = ((TObjString*)iWord())->GetString();
      if(sitem=="deg") tmp*=TMath::DegToRad();
      fMaxAngle = TMath::Pi()-tmp;
    }
    else if(sitem.Contains("/gps/pos/maxtheta")){
      iWord++;
      tmp =  ((TObjString*)iWord())->GetString()->Atof();
      iWord++;
      sitem = ((TObjString*)iWord())->GetString();
      if(sitem=="deg") tmp*=TMath::DegToRad();
      fMinAngle = TMath::Pi()-tmp;
    }
    else if(sitem.Contains("/gps/ene/type")){
      iWord++;
      sitem = ((TObjString*)iWord())->GetString();
      if(!sitem.Contains("Mono")){
	std::cout<<"PointSource::Init(): /gps/ene/type "<<sitem<<" not implemented, only Mono available"<<std::endl;
      }
    }
    else if(sitem.Contains("/gps/energy")){
      iWord++;
      fEnergy = ((TObjString*)iWord())->GetString()->Atof();
      iWord++;
      sitem = ((TObjString*)iWord())->GetString();
      if(sitem=="keV") fEnergy*=0.001;
      else if(sitem=="GeV") fEnergy*=1000;
    }
    else{
      std::cout<<"Command "<<sitem<<" unknown, thus ignored..."<<endl;
    }
    
    delete words;
  }
}

void PointSource:Print(){
  std::cout<<"\n\nPointSource::Print():\nName =\t"<<fName<<"\nEnergy =\t"<<fEnergy<<" MeV"<<"\nPosition=\t("<<fPosition.X()<<" , "<<fPosition.Y()<<" , "<<fPosition.Z()<<")\nAngularRange =\t("<<fMinAngle<<" , "<<fMaxAngle<<"\n"<<std::endl; 
}
