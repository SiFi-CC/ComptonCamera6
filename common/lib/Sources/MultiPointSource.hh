#ifndef __MultiPointSource_H_
#define __MultiPointSource_H_ 1

#include "PointSource.hh"

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
   *  Positions of all point sources is assumed to be specified relative to
   *  fPosition of this source.
   */
  void AddSourceElement(const PointSource& source);
  /** Generate particle */
  Track GenerateEvent() override;

private:
  /** List of point sources */
  std::vector<PointSource> fSources;
  /** Iterator pointing to one of list elements
   *
   * Iterator is not cyclic so GenerateEvent needs to contain logic
   * responsible to reseting interator at the end.
   */
  std::vector<PointSource>::iterator fIterator = fSources.begin();

  ClassDef(MultiPointSource, 1)
};

#endif
