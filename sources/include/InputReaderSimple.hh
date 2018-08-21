#ifndef __InputReaderSimple_H_
#define __InputReaderSimple_H_ 1
#include "TString.h"
#include "TVector3.h"
#include "InputReader.hh"
#include <iostream>

using namespace std;

///Class for accessing data from simple simulation performed by KR.
///This is class derived from InputReader class. It opens requested
///ROOT file containing tree with simulation results and via set of
///getter function passes information to reconstruction classes, i.e.
///CCREconstruction and CCMLEM.

class InputReaderSimple : public InputReader{
  
public:
  InputReaderSimple();
  InputReaderSimple(TString path);
  ~InputReaderSimple();
 
  void     Clear(void);
  TVector3 *GetPositionPrimary(void);
  TVector3 *GetPositionScattering(void);
  TVector3 *GetPositionAbsorption(void);
  TVector3 *GetGammaDirPrimary(void);
  TVector3 *GetGammaDirScattered(void);
  double   GetEnergyPrimary(void);
  double   GetEnergyLoss(void);
  double   GetEnergyScattered(void);
  
private:
  TVector3 *fPoint0;	///< Coordinates of the gamma source
  TVector3 *fPoint1;	///< Coordinates of the interaction in the scatterer (Comptin scattering)
  TVector3 *fPoint2;	///< Coordinates of the interaction in the absorber (absorption)
  TVector3 *fVersor1;	///< Direction of the primary gamma
  TVector3 *fVersor2;	///< Direction of the scattered gamma
  double   fEnergy0;	///< Energy of the primary gamma [MeV]
  double   fEnergy1;	///< Energy deposited in the scatterer [MeV]
  double   fEnergy2;	///< Energy of the scattered gamma [MeV]
  
  bool     AccessTree(TString name);
  
  ClassDef(InputReaderSimple,0)
};

#endif