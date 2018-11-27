#ifndef __LineSource_H_
#define __LineSource_H_ 1

#include "Source.hh"
#include <Math/Interpolator.h>

using ROOT::Math::Interpolator;

/** 1 dimensional Source positioned on line between two points, where energy is
 * interpolated based on provided list of values or Interpolator object.
 */
class LineSource : public Source {
public:
  LineSource() = default;

  /** Create line source
   * \param start start position
   * \param end end position
   * \param energies values of energy
   */
  LineSource(const TVector3& start, const TVector3& end,
             std::vector<Double_t> energies)
      : Source((start + end) * 0.5), fStartPosition(start), fEndPosition(end),
        fEnergyInterpolator(CreateEnergyInterpolator(energies)) {}

  /** Create line source
   * \param start start position
   * \param end end position
   * \param interpolator energy interpolator
   */
  LineSource(const TVector3& start, const TVector3& end,
             Interpolator* interpolator)
      : Source((start + end) * 0.5), fStartPosition(start), fEndPosition(end),
        fEnergyInterpolator(interpolator) {}

  virtual ~LineSource() {
    if (fEnergyInterpolator != nullptr) { delete fEnergyInterpolator; };
  }

  /** Create energy interpolator by passing list of values.
   * Energies will be interpolated between fStartPosition and fEndPosition,
   * using cubic splines(third-order polynomials) with exception of 2-element or
   * shorter lists where interpolation will be linear.
   * \param energies list of energy values
   */
  Interpolator* CreateEnergyInterpolator(std::vector<Double_t> energies);

  /** Set energy distribution by passing Interpolator. */
  void AdoptEnergyInterpolator(Interpolator* interpolator) {
    fEnergyInterpolator = interpolator;
  }

  /** Generates Track object that is positioned between fStartPosition and
   *  fEndPosition with energy specified by fEnergyInterpolator.
   */
  Track* GenerateEvent() override;

private:
  /** Line start position*/
  TVector3 fStartPosition;
  /** Line end position*/
  TVector3 fEndPosition;
  /** Interpolator define energy values between fStartPosition and
   * fEndPosition.
   */
  Interpolator* fEnergyInterpolator;

  ClassDef(LineSource, 0)
};

#endif
