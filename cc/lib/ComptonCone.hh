#ifndef __ComptonCone_H_
#define __ComptonCone_H_ 1
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

/// Class representing Compton Cone necessary for image reconstruction
/// from simulated data using back projection algorithm. ComptonCone
/// object consists of
///- apex (fApex) - point of interaction on the scatterer plane. It's
///coordiantes are stored in TVector3
///- axis (fAxis) - it is calculated as vector connecting interaction points in
///the scatterer and absorber. Coordiantes stored in TVector3
///- apex angle (fAngle) - calculated based on gamma energy loss during Compton
/// scattering.
///- name (fName) - name of the object.

class ComptonCone : public TObject {

public:
  ComptonCone();
  ComptonCone(TString name, TVector3 apex, TVector3 axis, Double_t angle);
  ComptonCone(TVector3* posScat, TVector3* posAbs, Double_t enPrimary,
              Double_t enScat);
  ~ComptonCone();

  void Print(void);

  /// Sets name of the object.
  void SetName(TString name) { fName = name; };
  /// Sets coordinates of the cone apex.
  void SetApex(TVector3 apex) { fApex = apex; };
  /// Sets cone axis.
  void SetAxis(TVector3 axis) { fAxis = axis; };
  /// Sets apex angle.
  void SetAngle(Double_t angle) { fAngle = angle; };
  /// Returns cone axis.
  TVector3 GetAxis(void) { return fAxis; };
  /// Returns coordinates of the apex.
  TVector3 GetApex(void) { return fApex; };
  /// Returns apex angle.
  Double_t GetAngle(void) { return fAngle; };
  /// Returns object name.
  const char* GetName() const { return fName.Data(); }

private:
  TString fName;   ///< Name of the object
  Double_t fAngle; ///< Apex angle [rad]
  TVector3 fApex;  ///< Coordinates of the apex
  TVector3 fAxis;  ///< Cone axis

  ClassDef(ComptonCone, 0)
};

#endif
