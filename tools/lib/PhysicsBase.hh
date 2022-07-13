#ifndef __PhysicsBase_H_
#define __PhysicsBase_H_ 1

#include "TVector3.h"

#include <tuple>

/*
 * Namespace containing methods, which are necessary in order to simulate Compton scattering of
 * gamma quanta in Compon Camera. Details of this class are described in presentation by KR
 * available on wiki:
 * [LINK](http://bragg.if.uj.edu.pl/gccbwiki/index.php/File:KR_20170222_CCandCarbonLine.pdf)
 */
namespace CC6
{
/**
 * Method for Compton scattering process.
 * 1. Checks whether incident Track hits the detector plane. If not returns NULL.
 * 2. Finds theta and phi scattering angles.
 * 3. Calculates gamma energy after scattering.
 * 4. Calculates coordinates of the new Track, representing scattered gamma quantum.
 * 5. Assignes values to the scattered Track object and returns it.

 * \param E incident gamma particle energy
 * \param p_versor incident gamma particle momentum versor
 * \return tuple of energy and track verson for scattered particle, or empty if angle incorrect.
*/
auto ComptonScatter(Double_t E, const TVector3& p_versor) -> std::tuple<Double_t, TVector3>;

/**
 * Generates Klein-Nishina function for requested energy of gamma quantum and finds value of theta
 * scattering angle. \param energy (Double_t) energy of incident gamma quantum [MeV]. \return
 * scattering angle [rad]
 * \return
 */
Double_t RandomKleinNishinaTheta(Double_t energy);

/**
 * Returns energy [MeV] of gamma quantum after Compton scattering.
 * \param theta theta scattering angle (must be given in radians)
 * \param initE (Double_t) initial energy of gamma quantum [MeV]
 * \return energy of scattered gamma
 */
Double_t ComptonScatteringGammaE(Double_t theta, Double_t initE);

/**
 * Namespace with helper literals. Use it with
 * ```using namespace CC6::literals;
 */
inline namespace literals
{
/**
 * Convert X_rad wher X is a int-type angle in degrees into angle in radians.
 * \param deh angle in degrees
 * \return angle in radians
 */
constexpr long double operator"" _rad(long double deg) { return deg * M_PI / 180; }
/**
 * Convert X_rad wher X is a double-type angle in degrees into angle in radians.
 * \param deg angle in degrees
 * \return angle in radians
 */
constexpr long double operator"" _rad(unsigned long long deg) { return deg * M_PI / 180; };
} // namespace literals

} // namespace CC6

#endif
