#ifndef __G4Reconstruction_H_
#define __G4Reconstruction_H_ 1

#include "BinnedGeometry.hh"
#include "CLog.hh"
#include "G4SimulationAdapter.hh"
#include <TBranch.h>
#include <TFile.h>
#include <TH2F.h>
#include <TObject.h>
#include <TVector3.h>

class G4Reconstruction : public TObject {
public:
  G4Reconstruction() = default;
  G4Reconstruction(CameraGeometry sim, TH2F* detector);
  virtual ~G4Reconstruction() = default;

  void RunReconstruction(int nIter);
  void Write(TString filename) const;

private:
  void SingleIteration();
  TMatrixT<Double_t> ReadFromTTree(TBranch* branch);
  TMatrixT<Double_t> ReadFromTH2F(TH2F* hist);

  std::vector<TMatrixT<Double_t>> fRecoObject;
  // vectorized column matrix representing detector image
  TMatrixT<Double_t> fImage;
  TMatrixT<Double_t> fMatrixH;
  TMatrixT<Double_t> fMatrixHTranspose;

  CameraGeometry fParams;

  SiFi::logger log = SiFi::createLogger("G4Reconstruction");

  ClassDef(G4Reconstruction, 0)
};

#endif
