#ifndef __MultiPointSource_H_
#define __MultiPointSource_H_ 1

#include "PointSource.hh"

/** Source constructed from multiple point sources */
class MultiPointSource : public Source
{
public:
    MultiPointSource() = default;

    /** Create multipoint source
     * \param fname filename of input file in style of Geant4 gps macro
     * individual point-like sources must be separated by "#####...#  Point i"
     * lines, example in data/sources
     */
    MultiPointSource(const TString fname);

    /** Add point source
     */
    void AddSourceElement(const PointSource& source);

    /** Generate particle */
    Track GenerateEvent() override;

    /** Printing MultiPointSource configuration */
    void Print() override;

private:
    /** Init method, empty for this class*/
    Bool_t Init() override { return kTRUE; };

    /** Normalization of source intensities, such that they sum up to 1 */
    Bool_t NormalizeIntensities();

    /** Dice source number, based on relative intensities*/
    Int_t DiceSource();

    /** List of point sources */
    std::vector<PointSource> fSources;

    /** List of source relative intensities */
    std::vector<Double_t> fIntensities;

    /** Number of sources in this MultiPointSource*/
    Int_t fNumOfSources = 0;

    /** Iterator pointing to one of list elements*/
    std::vector<PointSource>::iterator fIterator = fSources.begin(); //!

    ClassDefOverride(MultiPointSource, 1)
};

#endif
