#include "LineSource.hh"
#include <TRandom.h>

using ROOT::Math::Interpolation::kCSPLINE;
using ROOT::Math::Interpolation::kLINEAR;

Interpolator*
LineSource::CreateEnergyInterpolator(std::vector<Double_t> energies) {
  auto interpolationType = energies.size() <= 2 ? kLINEAR : kCSPLINE;

  // ensure at least 2 points
  if (energies.size() == 1) { energies.push_back(energies[0]); }
  if (energies.size() == 0) { energies = {0, 0}; }

  Double_t lineLength = (fEndPosition - fStartPosition).Mag();

  // generate equidistant points
  std::vector<Double_t> indicies(energies.size());
  for (int i = 0; i < energies.size(); i++) {
    indicies[i] = (lineLength / energies.size() - 1) * i;
  }

  return new Interpolator(indicies, energies, interpolationType);
}

Track LineSource::GenerateEvent() {
  TVector3 line = fEndPosition - fStartPosition;
  Double_t length = line.Mag();
  TVector3 pathVersor = line.Unit();
  Double_t random = gRandom->Uniform(0, length);

  Double_t angleY = gRandom->Uniform(-fAngleY, fAngleY);
  Double_t angleZ = gRandom->Uniform(-fAngleZ, fAngleZ);

  TVector3 point = fStartPosition + (pathVersor * random);
  TVector3 versor = TVector3(-1, tan(angleY), tan(angleZ)).Unit();
  Double_t energy = fEnergyInterpolator->Eval(random);
  return Track(point, versor, energy);
}
