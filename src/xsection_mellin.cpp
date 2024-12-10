// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <complex>                                      //
// -------- Classes ----------------------------------- //
#include "constants.h"    // Constants                  //
#include "process.h"      // Process definition         //
#include "settings.h"     // Setting routines           //
#include <iostream>


std::complex<double> set_hat_cross_section(const std::complex<double> &N, const int &id)
{
  // Numerical stability safeguard
  //if (abs(N) > 140.) return 0.;

  // Variables
  std::complex<double> Nb=N*exp(-Psi(1.)), lambda0=alpha_s/(2.*M_PI)*beta0*log(Nb);
  const double sigma0 = 2./3.*alpha*pow(M_PI/MZ,2);// eq (A.6): Cf*alpha*pi²/(2*Mz²), the rest in the PDFs

  // Determine g1, g2, and H based on Col_imp
  const bool collinear = (id==6);
  const std::complex<double> sudakov_log = log(1.-2.*lambda0);
  const std::complex<double> sudakov_logterm = collinear ? 1.0 : (1.-2.*lambda0);
  const std::complex<double> gamma_term = collinear ? (1./beta0*gamma_0*sudakov_log) : 0.0;

  std::complex<double> g1 = 1./(2.*lambda0*beta0)*Gamma_0*(2.*lambda0 + sudakov_logterm*sudakov_log);
  std::complex<double> g2 = beta1*pow(beta0,-3)/4.*Gamma_0*(4.*lambda0+2.*sudakov_log+pow(sudakov_log, 2))
      - pow(beta0,-2)/4.*Gamma_1*(2.*lambda0 + sudakov_log)
      + 1./beta0*Gamma_0*sudakov_log*log(MZ/muR)
      - 2./beta0*Gamma_0*lambda0*log(muR/muF)
      + gamma_term;
  double H = conversion_factor*sigma0*pow(MZ,2)/(2.*M_PI)*(1.+alpha_s/M_PI*4./3.* ((collinear ? (2.*pow(M_PI,2.)/3.) : (-2.*pow(log(MZ/muR),2) + 3.*log(MZ/muF) + 4.*pow(M_PI,2.)/3.)) - 4.));

  // Compute Sudakov factor and result
  return H*std::exp(g1*log(Nb)+g2);

  //test
  //g1 = Gamma_0/beta0*lambda0;
  //return H*std::exp(g1*log(Nb));
  
}

std::complex<double> set_hat_NLO_cross_section (const std::complex<double> N, void *prm)// dsigma_hat/(dM² dphi_1) 
{
  double sigma0=0, Cf=0, dZ=0, tau=0;
  std::complex<double> Phat=0, Pbar=0;
  Process *p = (Process *)prm;
  tau = pow(MZ,2.)/p->sh;
  Cf=4./3.;
  sigma0=2./3.*alpha*pow(M_PI/MZ,2); 
  dZ=Cf*(-4.+pow(M_PI,2)/3.)+3.*Cf*log(MZ/muF)+4.*Cf*log(1.-tau)*log(MZ/muF)+2.*Cf*pow(log(1.-tau),2);
  Pbar=Cf*(-pow(2. + N,-2) - pow(3. + N,-2) + 2.*PolyGamma(1, 2. + N))*2.*log(MZ/muF);
  Phat=2.*Cf*(-pow(2. + N,-2) - pow(3. + N,-2) + 2.*PolyGamma(1, 2. + N));
  return sigma0*(1.+alpha_s/M_PI*(dZ+Pbar+Phat));
}