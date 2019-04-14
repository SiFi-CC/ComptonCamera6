#ifndef __BinnedGeometry_H_
#define __BinnedGeometry_H_ 1

#include <TFile.h>
#include <TH2F.h>
#include <TObject.h>
#include <TVector3.h>

struct BinnedGeometry {
  std::pair<double, double> xRange;
  std::pair<double, double> yRange;
  std::pair<double, double> zRange;

  int binX;
  int binY;
  int binZ;

  int nBins() { return binX * binY * binZ; };

  std::tuple<int, int, int> getBin(double x, double y, double z);
  std::tuple<double, double, double> getBinCenter(int x, int y, int z);
  bool isValidBin(std::tuple<int, int, int> bin);
  bool isInside(std::tuple<double, double, double> pos);
};

#endif
