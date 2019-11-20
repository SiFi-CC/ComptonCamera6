#include "BinnedGeometry.hh"

std::tuple<int, int, int> BinnedGeometry::getBin(double x, double y, double z) {
  auto segXSize = (xRange.second - xRange.first) / binX;
  auto segYSize = (yRange.second - yRange.first) / binY;
  auto segZSize = (zRange.second - zRange.first) / binZ;

  auto nBinX = static_cast<int>((x - xRange.first) / segXSize) + 1;
  auto nBinY = static_cast<int>((y - yRange.first) / segYSize) + 1;
  auto nBinZ = static_cast<int>((z - zRange.first) / segZSize) + 1;

  nBinX = nBinX < 0 ? 0 : nBinX;
  nBinX = nBinX > binX ? binX + 1 : nBinX;
  nBinY = nBinY < 0 ? 0 : nBinY;
  nBinY = nBinY > binY ? binY + 1 : nBinY;
  nBinZ = nBinZ < 0 ? 0 : nBinZ;
  nBinZ = nBinZ > binZ ? binZ + 1 : nBinZ;

  return std::make_tuple(nBinX, nBinY, nBinZ);
}

std::tuple<double, double, double> BinnedGeometry::getBinCenter(int x, int y,
                                                                int z) {
  auto segXSize = (xRange.second - xRange.first) / binX;
  auto segYSize = (yRange.second - yRange.first) / binY;
  auto segZSize = (zRange.second - zRange.first) / binZ;

  return std::make_tuple(xRange.first + segXSize * (x - 0.5),
                         yRange.first + segYSize * (y - 0.5),
                         zRange.first + segZSize * (z - 0.5));
}

bool BinnedGeometry::isValidBin(std::tuple<int, int, int> bin) {
  return std::get<0>(bin) > 0 && std::get<0>(bin) <= binX &&
         std::get<1>(bin) > 0 && std::get<1>(bin) <= binY &&
         std::get<2>(bin) > 0 && std::get<2>(bin) <= binZ;
}

bool BinnedGeometry::isInside(std::tuple<double, double, double> pos) {
  return std::get<0>(pos) >= xRange.first &&
         std::get<0>(pos) <= xRange.second &&
         std::get<1>(pos) >= yRange.first &&
         std::get<1>(pos) <= yRange.second &&
         std::get<2>(pos) >= zRange.first && std::get<2>(pos) <= zRange.second;
}
