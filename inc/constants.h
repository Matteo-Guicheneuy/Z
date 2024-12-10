#ifndef CONSTANTS_H
#define CONSTANTS_H
// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
#include <cmath>
// ------------------EW Constants ----------------------//
//Constants found in the MadGraph configuration
const double MZ = 91.188; //Z boson mass in GeV
const double GF = 1.16639e-5; //Fermi constant
//double alpha_b = 7.29735257e-3;  //alpha
const double alpha_s = 0.118; //alpha_s(MZ)
//double alpha_s = 0.05; //alpha_s(MZ)
const double alphaM1 = 1.32507e+2; //1/alpha
const double alpha = 1/alphaM1;//7.5467711139788835e-3;  //alpha(MZ)
//double tw = (1.+sqrt(1.-4.*M_PI*alpha/(sqrt(2.)*std::pow(MZ,2.)*GF)))/2.; //tw=cos(theta_w)^2=mW^2/mZ^2
const double mw2 = MZ*MZ/2. + sqrt(std::pow(MZ,4)/4. - alpha*M_PI*MZ*MZ/(GF*sqrt(2)));
const double tw = mw2/(MZ*MZ);
const double tw2 = std::pow(0.88190334743339216,2.); // cos(theta_w)^2=MW^2/MZ^2 at scale MZ

const double gu=(8./9.*std::pow(1-tw,2.)+1./4.-2./3.*(1-tw))/(tw*(1-tw)); //g_L^2+g_R^2 for u type quark
const double gd=(std::pow(1-tw,2.)*2./9.+1./4.-(1-tw)/3.)/(tw*(1-tw)); //g_L^2+g_R^2 for d type quark

extern double muF;
extern double muR;

const int Nflav=5;

// Precision for numerical derivative
const double eps = 0.0001;

// Mellin space constants
const double beta0 = 23./6., beta1 = 29./3.;
const double Gamma_0 = 16./3., Gamma_1 = 16./3.* (67./3. - std::pow(M_PI, 2) - 50./9.);
const double gamma = -4.;
const double conversion_factor = 0.38937966e9; // GeV^-2 to pb

#endif