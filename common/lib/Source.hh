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

  /** Create source based on the input file
      \param fname input file name
   */
  Source(const TString fname) {
    fInFileName = fname;
    fPosition.SetXYZ(0,0,0);
  };

  /** Generates particle represented by Track object. */
  virtual Track GenerateEvent() = 0;

  virtual void Print() = 0;

  /** Set min and max values of angles between generated particles and -X axis.
   *  \param minAngle min angle between negative part of axis -X and 
   *  direction of generated particle.
   *  \param maxAngle max angle between negative part of axis -X and 
   *  direction of generated particle.
   *  In .mac files Geant4 interprets these angles as given with respect to
   *  the negative Z axis, and out Geant4 simulations really work in 
   *  the coordinate frame such, that gammas fly mostly along z axis.
   *  In simple simulations, though, and in reco, gammas should fly mostly
   *  along negatice x axis. A rotation from one to the other coordinate system
   *  is taken care of in the source implementation classes.
   */
  void SetAngleRange(Double_t minAngle, Double_t maxAngle) {
    fMinAngle = minAngle;
    fMaxAngle = maxAngle;
  }
  /** Set  source position   */
  void SetPosition(const TVector3& position) {
    fPosition = position;
  }

  /** Get name */
  const char* GetName() const override { return fName.Data(); }
  /** Set name */
  void SetName(TString name) { fName = name; }

protected:
  
  /** Init() initializes source properties based on configuration file */
  virtual Bool_t Init() = 0;
  
  /** \brief Position of the source.
   *
   *  Exact meaning of this variable might differ between implementations,
   *  in most cases it will define where the center is located.
   */
  TVector3 fPosition;
  
  /** File name of input file containing source configuration. */
  TString fInFileName;

  /** Min angle between -X axis and direction of generated particle. */
  Double_t fMinAngle = 0;

  /** Max angle between -X axis and direction of generated particle. */
  Double_t fMaxAngle = TMath::PiOver2();
  
private:
  TString fName = "generic_source"; ///< Object name

  //ClassDef(Source, 2)
};

#endif
