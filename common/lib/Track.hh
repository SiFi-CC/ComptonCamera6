#ifndef __Track_H_
#define __Track_H_ 1
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

class DetPlane;

/// Class, which represents single gamma ray. The Track object consists of:
///- starting point i.e. place of origin of gamma ray,
///- leading versor i.e. direction of the gamma ray,
///- energy - gamma energy in MeV,
///- name.

class Track : public TObject {

public:
  Track() = default;
  Track(TVector3 point, TVector3 vec, Double_t energy);

  void SetPoint(TVector3 point);
  void SetVersor(TVector3 vec);
  void SetEnergy(Double_t energy);
  void Print() const;

  /// Returns energy assigned to the track [MeV].
  Double_t GetEnergy() const { return fEnergy; };
  /// Returns Track leading versor coordinates.
  TVector3 GetVersor() const { return fVersor; };
  /// Returns Track starting point coordinates.
  TVector3 GetPoint() const { return fPoint; };

private:
  /** Coordinates of the starting point of the Track */
  TVector3 fPoint = TVector3(0, 0, 0);
  TVector3 fVersor = TVector3(-1, 0, 0); ///< Leading versor of the track
  Double_t fEnergy = 4.44;               ///< Energy assigned to the Track [MeV]

  ClassDef(Track, 1)
};

#endif
