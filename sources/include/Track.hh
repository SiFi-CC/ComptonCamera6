#ifndef __Track_H_
#define __Track_H_ 1
#include "TVector3.h"
#include "TString.h"
#include "TObject.h"

class DetPlane;

///Class, which represents single gamma ray. The track object consists of:
///- starting point i.e. place of origin of gamma ray, 
///- leading versor i.e. direction of the gamma ray,
///- energy - gamma energy in MeV,
///- name.

class Track : public TObject{
  
public:
  Track();
  Track(TVector3 point, TVector3 vec, Double_t energy, TString name);
  ~Track();
  
  void     SetPoint(TVector3 point);
  void     SetVersor(TVector3 vec);
  void     SetEnergy(Double_t energy);
  Bool_t   FindCrossPoint(DetPlane* plane, TVector3 &position);
  void     Print(void);
  
  ///Returns energy assigned to the track [MeV].
  Double_t GetEnergy(void) { return fEnergy; };
  ///Returns Track's leading versor coordinates.
  TVector3 GetVersor(void){ return fVersor; };
  ///Returns Track's starting point coordinates.
  TVector3 GetPoint(void){ return fPoint; };
  ///Sets name of the Track.
  void     SetName(TString name){fName = name;};
  ///Returns name of the Track.
  const char*  GetName() const { return fName.Data(); }
  
private:
  
  TVector3 fPoint;	///< Coordinates of the starting point of the Track
  TVector3 fVersor;	///< Leading versor of the track
  Double_t fEnergy;	///< Energy assigned to the Track [MeV]
  TString  fName;	///< Name of the object
  
  ClassDef(Track,0)
  
};

#endif
