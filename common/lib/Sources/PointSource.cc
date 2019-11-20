#include "PointSource.hh"
#include "CLog.hh"
#include "TCollection.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom.h"

// -----------------------
// PointSource
// -----------------------

PointSource::PointSource(const TString fname) {
  SetName("PointSource");
  fInFileName = fname;
  if (!Init()) {
    TString message;
    message.Form("Exception in %s constructor!", this->ClassName());
    throw message.Data();
  }
  CreateHistograms();
}

Track PointSource::GenerateEvent() {
  Double_t theta = gRandom->Uniform(fMinAngle, fMaxAngle);
  TVector3 versor;
  versor.SetMagThetaPhi(1., theta, gRandom->Uniform(0, TMath::TwoPi()));
  versor.RotateY(-TMath::Pi() / 2);
  Track track(fPosition, versor, fEnergy);
  if (fhTheta) fhTheta->Fill(track.GetVersor().Theta());
  if (fhPhi) fhPhi->Fill(track.GetVersor().Phi());
  if (fhZY) fhZY->Fill(track.GetPoint().Z(), track.GetPoint().Y());
  return track;
}

Bool_t PointSource::Init() {
  if (!fInFileName.EndsWith(".mac")) {
    std::cout << "PointSource::Init(): not appropriate input file format, .mac "
                 "file expected..."
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
  Bool_t result = Init(infile);
  infile.close();
  return result;
}

Bool_t PointSource::Init(std::ifstream& infile) {
  SetName("PointSource");
  TString line = "";
  while (!infile.eof() && !(line.Contains("####") && line.Contains("Point"))) {
    line.ReadLine(infile);
    // std::cout<<"Read Line : "<<line<<std::endl;
    if (line.BeginsWith("#") || line.IsWhitespace()) continue;

    TObjArray* words = line.Tokenize(" ");
    // for(int i=0; i<words->GetEntries(); i++)
    //    std::cout<<i<<"="<<((TObjString*)words->At(i))->GetString()<<",\t";
    // std::cout<<std::endl;
    if (words->GetEntries() == 0) continue;
    TString sitem;
    Double_t tmp;
    sitem = ((TObjString*)words->At(0))->GetString();
    if (sitem.Contains("/gps/pos/type")) {
      sitem = ((TObjString*)words->At(1))->GetString();
      if (!sitem.Contains("Point")) {
        std::cout << "PointSource::Init(): Inconsistency in source types. "
                  << "Point expected, " << sitem
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
    } else if (sitem.Contains("/gps/source/intensity") ||
               sitem.Contains("/gps/source/add")) {
      fIntensity = ((TObjString*)words->At(1))->GetString().Atof();
    } else if (sitem.Contains("/gps/ang/maxtheta")) {
      tmp = ((TObjString*)words->At(1))->GetString().Atof();
      sitem = ((TObjString*)words->At(2))->GetString();
      if (sitem == "deg") tmp *= TMath::DegToRad();
      fMinAngle = TMath::Pi() - tmp;
    } else if (sitem.Contains("/gps/ene/type")) {
      sitem = ((TObjString*)words->At(1))->GetString();
      if (!sitem.Contains("Mono")) {
        std::cout << "PointSource::Init(): /gps/ene/type " << sitem
                  << " not implemented, only Mono available" << std::endl;
      }
    } else if (sitem.Contains("/gps/energy")) {
      fEnergy = ((TObjString*)words->At(1))->GetString().Atof();
      sitem = ((TObjString*)words->At(2))->GetString();
      if (sitem == "keV")
        fEnergy *= 0.001;
      else if (sitem == "GeV")
        fEnergy *= 1000;
    } else {
      std::cout << "Command " << sitem << " unknown, thus ignored..."
                << std::endl;
    }
    delete words;
  }
  return kTRUE;
}

void PointSource::Print() {
  std::cout << "\n\nPointSource::Print():\nName =\t" << GetName()
            << "\nEnergy =\t" << fEnergy << " MeV"
            << "\nPosition=\t(" << fPosition.X() << " , " << fPosition.Y()
            << " , " << fPosition.Z() << ")\nAngularRange =\t(" << fMinAngle
            << " , " << fMaxAngle << ") rad wrt -X axis\nIntensity = \t"
            << fIntensity << std::endl;
}
