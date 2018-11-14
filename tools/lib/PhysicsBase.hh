#ifndef __PhysicsBase_H_
#define __PhysicsBase_H_ 1
#include "DetPlane.hh"
#include "TF1.h"
#include "TObject.h"
#include "Track.hh"

/// Class containing methods, which are necessary in order to
/// simulate Compton scattering of gamma quanta in Compon Camera.
/// Details of this class are described in presentation by KR
/// available on wiki:
///[LINK](http://bragg.if.uj.edu.pl/gccbwiki/index.php/File:KR_20170222_CCandCarbonLine.pdf)

class PhysicsBase : public TObject {

public:
  PhysicsBase();
  PhysicsBase(TString name);
  ~PhysicsBase();

  Double_t NewEnergy(Double_t theta, Double_t initE);
  Double_t FindPhi(void);
  Double_t FindTheta(Double_t energy);
  Track* ComptonScatter(Track* initTrack, DetPlane* plane);
  void Print(void);

  /// Sets object name.
  void SetName(TString name) { fName = name; };
  /// Returns object name.
  const char* GetName() const { return fName.Data(); };
  /// Returns theta angle [rad]
  Double_t GetTheta(void) { return fTheta; };
  /// Returns phi angle [rad]
  Double_t GetPhi(void) { return fPhi; };

private:
  TString fName;   ///< Object name
  Double_t fTheta; ///< Theta angle [rad]
  Double_t fPhi;   ///< Phi angle [rad]
  TF1* fFunction;  ///< Klein-Nishina function for requested energy

  ClassDef(PhysicsBase, 0)
};

#endif
