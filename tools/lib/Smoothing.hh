#include <TH2F.h>
#include <TF1.h>

#include "CLog.hh"
#ifndef __Smoothing_H_
#define __Smoothing_H_ 1

namespace SiFi {
namespace tools {


TH2F* SmoothGauss(TH2F* hin, double sigma);

std::vector<Double_t> UQI_MSE(TH2F* sourceImage, TH2F* recoImage);


} // namespace tools
} // namespace SiFi

#endif
