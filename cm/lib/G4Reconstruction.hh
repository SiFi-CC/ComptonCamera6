#ifndef __G4Reconstruction_H_
#define __G4Reconstruction_H_ 1

#include "BinnedGeometry.hh"
#include "CLog.hh"
#include <TFile.h>
#include <TH2F.h>
#include <TObject.h>
#include <TVector3.h>

struct SimulationParams {
  BinnedGeometry source;
  BinnedGeometry mask;
  BinnedGeometry detector;
  std::vector<TFile*> recoData;

  void initSimulationMetadata();
};

class G4Reconstruction : public TObject {
public:
  G4Reconstruction() = default;
  G4Reconstruction(SimulationParams sim, TH2F* detector);
  virtual ~G4Reconstruction() = default;

  void RunReconstruction(int nIter);
  void Write(TString filename) const;

private:
  void SingleIteration();

  std::vector<TMatrixT<Double_t>> fRecoObject;
  TMatrixT<Double_t> fImage;
  TMatrixT<Double_t> fMatrixH;
  TMatrixT<Double_t> fMatrixHTranspose;

  SimulationParams fParams;

  SiFi::logger log = SiFi::createLogger("G4Reconstruction");

  ClassDef(G4Reconstruction, 0)
};

#endif
