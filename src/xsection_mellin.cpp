// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <complex>                                      //
// -------- Classes ----------------------------------- //
#include "constants.h"    // Constants                  //

std::complex<double> Psi(std::complex<double> N)
{
  std::complex<double> res(0.,0.);
  while(real(N)<10.) { res=res-1./N; N=N+1.; }
  res=res+log(N)-1./(2.*N)-pow(N,-2.)/2520.*(210.+pow(N,-2.)*(-21.+pow(N,-2.)*10.));
  return res;
}


std::complex<double> set_hat_cross_section(const std::complex<double> &N, const int &id)
{
  // Numerical stability safeguard
  if (abs(N) > 140.) return 0.;

  // Variables
  std::complex<double> Nb=N*exp(-Psi(1.)), lambda0=alpha_s/(2.*M_PI)*beta0*log(Nb);
  const double prefactor = 2./3.*alpha*pow(M_PI/MZ,2);// eq (A.6): Cf*alpha*pi²/(2*Mz²), the rest in the PDFs

  // Determine g1, g2, and H based on Col_imp
  const bool collinear = (id==6);
  const std::complex<double> sudakov_log = log(1.-2.*lambda0);
  const std::complex<double> sudakov_logterm = collinear ? 1.0 : (1.-2.*lambda0);
  const std::complex<double> gamma_term = collinear ? (1./beta0*gamma*sudakov_log) : 0.0;

  std::complex<double> g1 = 1./(2.*lambda0*beta0)*Gamma_0*(2.*lambda0 + sudakov_logterm*sudakov_log);
  std::complex<double> g2 = beta1*pow(beta0,-3)/4.*Gamma_0*(4.*lambda0+2.*sudakov_log+pow(sudakov_log, 2))
      - pow(beta0,-2)/4.*Gamma_1*(2.*lambda0 + sudakov_log)
      + 1./beta0*Gamma_0*sudakov_log*log(MZ/muR)
      - 2./beta0*Gamma_0*lambda0*log(muR/muF)
      + gamma_term;

  double H = conversion_factor*prefactor*pow(MZ,2)/(2.*M_PI)*(1.+alpha_s/M_PI*4./3.* ((collinear ? (2.*pow(M_PI,2.)/3.) : (-2.*pow(log(MZ/muR),2) + 3.*log(MZ/muF) + 4.*pow(M_PI,2.)/3.)) - 4.));

  // Compute Sudakov factor and result
  return H*std::exp(g1*log(Nb)+g2);
}



