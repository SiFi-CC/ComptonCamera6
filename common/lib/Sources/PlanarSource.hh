#ifndef __PlanarSource_H_
#define __PlanarSource_H_ 1

#include "Source.hh"
#include "TF1.h"

/** 2 dimensional, rectangular source, where energy is either fixed,
 * or energy dist. is linear.
 */
class PlanarSource : public Source {
public:
  PlanarSource() = default;

  /** Create planar, rectangular source from configuration file, examples in
   * data/sources/  */
  PlanarSource(const TString fname);

  /** Generate particle */
  Track GenerateEvent() override;

  /** Printing PlanarSource configuration */
  void Print() override;

private:
  /** Init() initializes source properties based on configuration file */
  Bool_t Init() override;

  /** energy of gammas emitted from this source, if discrete */
  Double_t fEnergy;

  /** linear energy distribution of gammas emitted from this source, if this
   * option chosen */
  TF1* fEnergyDist = NULL;

  /** lower limit of energy of gammas emitted from this source, if continuous
   * dist. */
  Double_t fMinEnergy;

  /** upper limit of energy of gammas emitted from this source, if continuous
   * dist. */
  Double_t fMaxEnergy;

  /** half size of the rectangle measured along horizontal axis*/
  Double_t fHalfZ;

  /** half size of the rectangle measured along vertical axis*/
  Double_t fHalfY;

  ClassDefOverride(PlanarSource, 1)
};

#endif
