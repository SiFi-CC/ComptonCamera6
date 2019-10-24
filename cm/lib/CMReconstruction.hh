#ifndef __CMReconstruction_H_
#define __CMReconstruction_H_ 1
#include "Coordinates.hh"
#include "Mask.hh"
#include "Source.hh"
#include "Sources/PointSource.hh"
#include "Track.hh"
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMatrixT.h>
#include <TTree.h>
#include <vector>

/**
 *  Reconstruction algorithm:
 *  - calculate H matrix
 *
 *   Matrix represent probability of j-th bin of detector to be hit by source at
 * i-th bin of source. Matrix is genereted by positioning source at different
 * position in XY grid and simulating resutlts for that source.
 *
 *  - calculate S normalization vector
 *
 *    S[i] - sum of probabilites for i-th source postion
 *
 *  - MLEM iterations
 *
 *    f_{k+1} = f_k/S * H' (Image / H * f_k)
 *
 *    based on epjconf-I106-05003-2016
 */
class CMReconstruction : public TObject {
public:
  CMReconstruction() = default;
  CMReconstruction(TString simulationFile);
  virtual ~CMReconstruction();

  /** Start reconstruction */
  void RunReconstruction(Int_t nIterations);

  /** Save to file */
  void Write(TString filename) const;

  /** Calculate H matrix and save it in file */
  void HmatrixToFile(TString filename, Int_t t);

private:
  /** calulate probability matrix */
  void FillHMatrix();
  /** run single iteration of simulation */
  void SingleIteration();
  /** calculate vector s to normalize probabilities */
  Bool_t CalculateS();

private:
  /** List of reconstructed objects for every iteration. Every object is kept as
   * column vector.
   */
  std::vector<TMatrixT<Double_t>> fRecoObject;

  /** Probability matrix
   *  - row represent i-th detector pixel
   *  - column represent j-th position of point source
   *
   *  H[i, j] represents probability of hitting i-th pixel of detector from j-th
   * segment of source.
   */
  TMatrixT<Double_t> fMatrixH;

  /** Transposition of  H matrix */
  TMatrixT<Double_t> fMatrixHPrime;

  /** Normalization terms for reconstruction itertions */
  std::vector<Double_t> fNormS;

  /** detector image vectorized into column matrix */
  TMatrixT<Double_t> fImageMat;

  /** image from simimulation result file */
  TH2F fImage;
  /** object that stores coordinates/dimensions of source object */
  H2Coords fObjectCoords;
  /** object that stores coordinates/dimensions of detector image */
  H2Coords fImageCoords;

  /** object from simimulation result file, only data relating to dimensions are
   * used.
   * TODO: remove reference after switching to angular segments in reconstructed
   * object.
   */
  TH2F fObject;

  /** Detector plane used in simulation */
  DetPlane fDetPlane;
  /** Mask used in simulation */
  Mask fMask;

  SiFi::logger log = SiFi::createLogger("CMReconstruction");

  ClassDef(CMReconstruction, 0)
};

#endif
