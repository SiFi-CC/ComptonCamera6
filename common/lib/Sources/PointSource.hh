#ifndef __PointSource_H_
#define __PointSource_H_ 1

#include "Source.hh"

/** Source with defined position and enrgy that generates gammas only in that
 * point with that energy
 */
class PointSource : public Source {
public:
  PointSource() = default;
  /** Create point surce
   * \param position source position
   * \param energy initial energy of particles generated by this source.
   */
  PointSource(const TVector3& position, Double_t energy)
      : Source(position), fEnergy(energy){};

  /** Generte particle */
  Track* GenerateEvent() override;

private:
  /** energy of gammas emited from this source */
  Double_t fEnergy;

  ClassDef(PointSource, 0)
};

/** Source constructed from multiple point sources */
class MultiPointSource : public Source {
public:
  MultiPointSource() = default;

  /** Create multipoint source
   * \param position position of source (all point source positions are treated
   * as coordinates relative to this.
   */
  MultiPointSource(const TVector3& position) : Source(position){};

  /** Add point source
   *  Positions of all point sources is assumed to be speciied relative to
   *  fPosition of this source.
   */
  void AddSourceElement(const PointSource& source);
  /** Generate particle */
  Track* GenerateEvent() override;

private:
  /** List of point sources */
  std::vector<PointSource> fSources;
  /** Iterator pointing to one of list elements
   *
   * Iterator is not cyclic so GenerateEvent needs to ccontain logic
   * responsible to reseting interator at the end.
   */
  std::vector<PointSource>::iterator fIterator = fSources.begin();

  ClassDef(MultiPointSource, 0)
};

#endif
