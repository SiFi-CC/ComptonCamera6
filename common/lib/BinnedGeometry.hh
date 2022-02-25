#ifndef __BinnedGeometry_H_
#define __BinnedGeometry_H_ 1

#include "CLog.hh"
#include <TFile.h>
#include <TH2F.h>
#include <TObject.h>
#include <TVector3.h>

// Describes dimensions of 3d object that can be divided into segemnts along any
// of three axis
struct BinnedGeometry
{
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

    void Print()
    {
        spdlog::info("geometry info\n"
                     "range - x ({}, {}), y ({}, {}), z({}, {})\n"
                     "bins - ({}, {}, {})",
                     xRange.first, xRange.second, yRange.first, yRange.second, zRange.first,
                     zRange.second, binX, binY, binZ);
    }
};

#endif
