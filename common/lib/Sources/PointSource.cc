#include "PointSource.hh"
#include <TRandom.h>

// -----------------------
// PointSource
// -----------------------

Track PointSource::GenerateEvent() {
  Double_t angleY = gRandom->Uniform(-fAngleY, fAngleY);
  Double_t angleZ = gRandom->Uniform(-fAngleZ, fAngleZ);
  TVector3 versor = TVector3(-1, tan(angleY), tan(angleZ)).Unit();
  return Track(fPosition, versor, fEnergy);
}

// -----------------------
// MultiPointSource
// -----------------------

Track MultiPointSource::GenerateEvent() {
  Track track = fIterator->GenerateEvent();
  track.SetPoint(track.GetPoint() + fPosition);

  if (fIterator == fSources.end()) {
    fIterator = fSources.begin(); // on last element go to start
  } else {
    ++fIterator;
  }
  return track;
}

void MultiPointSource::AddSourceElement(const PointSource& source) {
  fSources.push_back(source);
  fIterator = fSources.begin();
}
