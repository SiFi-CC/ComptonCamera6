#include "MultiPointSource.hh"
#include "CLog.hh"
#include <TRandom.h>

// -----------------------
// MultiPointSource
// -----------------------

Track MultiPointSource::GenerateEvent() {
  if (fIterator == fSources.end()) {
    fIterator = fSources.begin(); // on last element go to start
  }
  Track track = fIterator->GenerateEvent();
  track.SetPoint(track.GetPoint() + fPosition);

  ++fIterator;
  return track;
}

void MultiPointSource::AddSourceElement(const PointSource& source) {
  fSources.push_back(source);
  fIterator = fSources.begin();
}
