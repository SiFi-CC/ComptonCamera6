#ifndef __Coordinates_H_
#define __Coordinates_H_ 1

#include "CLog.hh"
#include <TH2.h>
#include <sstream>

class H2Coords {
public:
  H2Coords() = default;
  H2Coords(TH2* hist)
      : fNBinsXY(std::make_pair(hist->GetXaxis()->GetNbins(),
                                hist->GetYaxis()->GetNbins())),
        fNBins(hist->GetXaxis()->GetNbins() * (hist->GetYaxis()->GetNbins())) {}
  Int_t NBinsX() const { return fNBinsXY.first; };
  Int_t NBinsY() const { return fNBinsXY.second; };
  Int_t NBins() const { return fNBins; };
  std::pair<Int_t, Int_t> BinXY(Int_t bin) const {
    if (bin >= fNBins) {
      spdlog::warn("invalid bin coordinate");
      return std::make_pair(-1, -1);
    }
    return std::make_pair(bin / fNBinsXY.first, fNBinsXY.second - bin % fNBinsXY.first);//VU
  }

  std::string String() const {
    std::stringstream ss;
    ss << "number of bins: " << fNBins << ", number of bins along x axis "
       << fNBinsXY.first << ", number of bins along y axis " << fNBinsXY.second;
    return ss.str();
  }

private:
  std::pair<Int_t, Int_t> fNBinsXY;
  Int_t fNBins;
};

#endif
