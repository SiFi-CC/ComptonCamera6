#ifndef __Source_H_
#define __Source_H_ 1
#include "Track.hh"
#include <TObject.h>
#include <TString.h>
#include <TVector3.h>

/** Abstract representation of gamma ray source.  */
class Source : public TObject {
public:
  Source() = default;

  /** Create source at certain position.
   * \param position source position
   */
  Source(const TVector3& position) : fPosition(position){};

  /** Generates particle represented by Track object. */
  virtual Track GenerateEvent() = 0;

  /** Set max values of angles between generated particles and Y and Z axis.
   *  \param angleY max angle between axis Y and direction of generated
   *  particle.
   *  \param angleZ max angle between axis Z and direction of generated
   *  particle.
   */
  void SetAngleRanges(Double_t angleY, Double_t angleZ) {
    fAngleY = angleY;
    fAngleZ = angleZ;
  }

  /** Get name */
  const char* GetName() const override { return fName.Data(); }
  /** Set name */
  void SetName(TString name) { fName = name; }

protected:
  /** \brief Position of the source.
   *
   *  Exact meaning of this variable might differ between implementations,
   *  in most cases it will define where the center is located.
   */
  TVector3 fPosition;

  /** Max angle between axis Y and direction of generated particle. */
  Double_t fAngleY = TMath::PiOver4();

  /** Max angle between axis Z and direction of generated particle. */
  Double_t fAngleZ = TMath::PiOver4();

private:
  TString fName = "generic_source"; ///< Object name

  ClassDef(Source, 1)
};

#endif
