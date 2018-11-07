#ifndef __DetPlane_H_
#define __DetPlane_H_ 1
#include "TObject.h"
#include "TString.h"
#include "TVector3.h"

/// Class which represents plane of the detector. The DetPlane
/// object consists of:
///- 4 coefficients of cartesian plane equation,
///- dimensions of the detector plane in Y and Z directions,
///- name.
///
/// All dimensions and distances are in mm.
/// Cartesian plane equation: Ax + By + Cy + D = 0, where D corresponds
/// to the distance (along X axis) between center of the coordinate system and
/// the detector plane.

class DetPlane : public TObject {

public:
  DetPlane();
  DetPlane(Double_t a, Double_t b, Double_t c, Double_t d, Double_t dimZ,
           Double_t dimY, TString name);
  ~DetPlane();

  void SetPlane(Double_t a, Double_t b, Double_t c, Double_t d);
  void SetDimensions(Double_t dimZ, Double_t dimY);
  TVector3 GetNormal(void);
  Bool_t CheckPoint(TVector3 point);
  void Print(void);

  /// Sets object name.
  void SetName(TString name) { fName = name; };
  /// Returns coefficient A of cartesian plane equation.
  Double_t GetA(void) { return fA; };
  /// Returns coefficient B of cartesian plane equation.
  Double_t GetB(void) { return fB; };
  /// Returns coefficient C of cartesian plane equation.
  Double_t GetC(void) { return fC; };
  /// Returns coefficient D of cartesian plane equation.
  Double_t GetD(void) { return fD; };
  /// Returns size of the detector plane along Z axis.
  Double_t GetDimZ(void) { return fDimZ; };
  /// Returns size of the detector plane along Y axis
  Double_t GetDimY(void) { return fDimY; };
  /// Returns object name.
  const char* GetName() const { return fName.Data(); }

private:
  Double_t fA; ///< Coefficient A of the cartesian plane equation
  Double_t fB; ///< Coefficient B of the cartesian plane equation
  Double_t fC; ///< Coefficient C of the cartesian plane equation
  Double_t fD; ///< Coefficient D of the cartesian plane equation
  Double_t
      fDimZ; ///< Size of the detector plane along Z axis (full length in mm)
  Double_t fDimY; ///< Size of the detector plane along Y (full length in mm)
  TString fName;  ///< Object name

  ClassDef(DetPlane, 0)
};

#endif
