#include <TH2F.h>
#include <TF1.h>

#include "CLog.hh"
#ifndef __Smoothing_H_
#define __Smoothing_H_ 1

namespace SiFi {
namespace tools {


TH2F* SmoothGauss(TH2F* hin, double sigma);


} // namespace tools
} // namespace SiFi

#endif
