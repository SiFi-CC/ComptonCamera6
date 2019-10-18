#include "PlanarSource.hh"
#include "CLog.hh"
#include "TCollection.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"
#include <TRandom.h>
#include <fstream>
#include <iostream>

PlanarSource::PlanarSource(const TString fname) {
  SetName("PlanarSource");
  fInFileName = fname;
  if (!Init()) {
    TString message;
    message.Form("Exception in %s constructor!", this->ClassName());
    throw message.Data();
  }
}

Track PlanarSource::GenerateEvent() {
  TVector3 offset(0, gRandom->Uniform(-fHalfY, fHalfY),
                  gRandom->Uniform(-fHalfZ, fHalfZ));
  Double_t theta = gRandom->Uniform(fMinAngle, fMaxAngle);
  TVector3 versor;
  versor.SetMagThetaPhi(1., theta, gRandom->Uniform(0, TMath::TwoPi()));
  versor.RotateY(-TMath::Pi() / 2);
  Track track;
  if (fEnergyDist)
    track = Track(fPosition + offset, versor,
                  fEnergyDist->GetRandom(fMinEnergy, fMaxEnergy));
  else
    track = Track(fPosition + offset, versor, fEnergy);
  fhTheta->Fill(track.GetVersor().Theta());
  fhPhi->Fill(track.GetVersor().Phi());
  fhZY->Fill(track.GetPoint().Z(), track.GetPoint().Y());
  return track;
}

Bool_t PlanarSource::Init() {
  if (!fInFileName.EndsWith(".mac")) {
    std::cout << "PlanarSource::Init(): not appropriate input file format, "
                 ".mac file expected..."
              << std::endl;
    return kFALSE;
  }
  std::ifstream infile;
  infile.open(fInFileName.Data(), std::ios::in);
  if (!infile.is_open() || infile.fail() || infile.bad()) {
    std::cout << "Input file " << fInFileName
              << " not opened correctly or corrupted" << std::endl;
    std::cout << infile.is_open() << " " << infile.fail() << " " << infile.bad()
              << std::endl;
    return kFALSE;
  }

  TString line;
  TObjArray* words;
  while (!infile.eof()) {
    // if(words) words->Delete();
    line.ReadLine(infile);
    // std::cout<<"Read Line : "<<line<<std::endl;
    if (line.BeginsWith("#") || line.IsWhitespace()) continue;

    words = line.Tokenize(" ");
    // for(int i=0; i<words->GetEntries(); i++)
    //    std::cout<<i<<"="<<((TObjString*)words->At(i))->GetString()<<",\t";
    // std::cout<<std::endl;
    if (words->GetEntries() == 0) continue;
    TString sitem = "";
    Double_t tmp;
    sitem = ((TObjString*)words->At(0))->GetString();
    if (sitem.Contains("/gps/pos/type")) {
      sitem = ((TObjString*)words->At(1))->GetString();
      if (!sitem.Contains("Plane")) {
        std::cout << "PlanarSource::Init(): Inconsistency in source types. "
                  << "Plane expected, " << sitem
                  << " found in source config file " << fInFileName
                  << std::endl;
      }
    }
    if (sitem.Contains("/gps/pos/shape")) {
      sitem = ((TObjString*)words->At(1))->GetString();
      if (!sitem.Contains("Rectangle")) {
        std::cout << "PlanarSource::Init(): Inconsistency in source shapes. "
                  << "Rectangle expected, " << sitem
                  << " found in source config file " << fInFileName
                  << std::endl;
      }
    } else if (sitem.Contains("/gps/pos/centre")) {
      Double_t pos[3];
      for (int i = 0; i < 3; i++) {
        pos[i] = ((TObjString*)words->At(i + 1))->GetString().Atof();
      }
      SetPosition(TVector3(pos));
      sitem = ((TObjString*)words->At(4))->GetString();
      if (sitem.Contains("cm"))
        if (fPosition.Mag() > 1e-14)
          fPosition.SetMag(fPosition.Mag() *
                           10.); // default units in sim. are mm
    } else if (sitem.Contains("/gps/ang/mintheta")) {
      tmp = ((TObjString*)words->At(1))->GetString().Atof();
      sitem = ((TObjString*)words->At(2))->GetString();
      if (sitem == "deg") tmp *= TMath::DegToRad();
      fMaxAngle = TMath::Pi() - tmp;
    } else if (sitem.Contains("/gps/ang/maxtheta")) {
      tmp = ((TObjString*)words->At(1))->GetString().Atof();
      sitem = ((TObjString*)words->At(2))->GetString();
      if (sitem == "deg") tmp *= TMath::DegToRad();
      fMinAngle = TMath::Pi() - tmp;
    } else if (sitem.Contains("/gps/ene/type")) {
      sitem = ((TObjString*)words->At(1))->GetString();
      if (!sitem.Contains("Mono") && !sitem.Contains("Lin")) {
        std::cout << "PlanarSource::Init(): /gps/ene/type " << sitem
                  << " not implemented, only Mono and Lin available"
                  << std::endl;
        if (sitem.Contains("Lin"))
          if (fEnergyDist == NULL)
            fEnergyDist = new TF1("fEnergyDist", "pol1", 0, 10000);
      }
    } else if (sitem.Contains("/gps/energy")) {
      fEnergy = ((TObjString*)words->At(1))->GetString().Atof();
      sitem = ((TObjString*)words->At(2))->GetString();
      if (sitem == "keV")
        fEnergy *= 0.001;
      else if (sitem == "GeV")
        fEnergy *= 1000;
    } else if (sitem.Contains("/gps/ene/gradient")) {
      if (fEnergyDist == NULL)
        fEnergyDist = new TF1("fEnergyDist", "pol1", 0, 10000);
      fEnergyDist->SetParameter(
          1, ((TObjString*)words->At(1))->GetString().Atof());
    } else if (sitem.Contains("/gps/ene/intercept")) {
      if (fEnergyDist == NULL)
        fEnergyDist = new TF1("fEnergyDist", "pol1", 0, 10000);
      fEnergyDist->SetParameter(
          0, ((TObjString*)words->At(1))->GetString().Atof());
    } else if (sitem.Contains("/gps/ene/min")) {
      fMinEnergy = ((TObjString*)words->At(1))->GetString().Atof();
      sitem = ((TObjString*)words->At(2))->GetString();
      if (sitem == "keV")
        fMinEnergy *= 0.001;
      else if (sitem == "GeV")
        fMinEnergy *= 1000;
    } else if (sitem.Contains("/gps/ene/max")) {
      fMaxEnergy = ((TObjString*)words->At(1))->GetString().Atof();
      sitem = ((TObjString*)words->At(2))->GetString();
      if (sitem == "keV")
        fMaxEnergy *= 0.001;
      else if (sitem == "GeV")
        fMaxEnergy *= 1000;
    } else if (sitem.Contains("/gps/pos/halfy")) {
      fHalfY = ((TObjString*)words->At(1))->GetString().Atof();
      sitem = ((TObjString*)words->At(2))->GetString();
      if (sitem.Contains("cm")) fHalfY *= 10; // default units in sim. are mm
    } else if (sitem.Contains("/gps/pos/halfx")) {
      fHalfZ = ((TObjString*)words->At(1))->GetString().Atof();
      sitem = ((TObjString*)words->At(2))->GetString();
      if (sitem.Contains("cm")) fHalfZ *= 10; // default units in sim. are mm
    } else {
      std::cout << "Command " << sitem << " unknown, thus ignored..."
                << std::endl;
    }

    delete words;
  }
  return kTRUE;
}

void PlanarSource::Print() {
  std::cout << "\n\nPlanarSource::Print():\nName =\t" << GetName()
            << "\nEnergy =\t" << fEnergy << " MeV"
            << "\nPosition=\t(" << fPosition.X() << " , " << fPosition.Y()
            << " , " << fPosition.Z() << ") mm"
            << "\nSize=\t(" << 0 << " , " << fHalfY * 2 << " , " << fHalfZ * 2
            << ") mm\nAngularRange =\t(" << fMinAngle << " , " << fMaxAngle
            << ") rad wrt -X axis" << std::endl;
  if (fEnergyDist)
    std::cout << "Energy generated according to distribution f(E) = "
              << fEnergyDist->GetParameter(0) << " + E*"
              << fEnergyDist->GetParameter(1) << " in the range (" << fMinEnergy
              << " , " << fMaxEnergy << ") MeV\n"
              << std::endl;
}
