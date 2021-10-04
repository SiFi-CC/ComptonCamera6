#ifndef __PhysicsBase_H_
#define __PhysicsBase_H_ 1

#include "TF1.h"
#include "TObject.h"
#include "TVector3.h"

/// Class containing methods, which are necessary in order to
/// simulate Compton scattering of gamma quanta in Compon Camera.
/// Details of this class are described in presentation by KR
/// available on wiki:
///[LINK](http://bragg.if.uj.edu.pl/gccbwiki/index.php/File:KR_20170222_CCandCarbonLine.pdf)
namespace CC6
{
auto ComptonScatter(Double_t E, const TVector3& p_versor, const TVector3& cross_point)
    -> std::pair<Double_t, TVector3>;

Double_t RandomKleinNishinaTheta(Double_t energy);

Double_t ComptonScatteringGammaE(Double_t theta, Double_t initE);

inline namespace literals
{
long double operator"" _rad(long double deg);
long double operator"" _rad(unsigned long long);
} // namespace literals

} // namespace CC6

#endif
