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
  if (abs(N) > 10.) return 0.;
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
  //std::cout << "g1'= " << 1./beta0*Gamma_0 << " g1''=" << 1./beta0*Gamma_0* (1./(2.*lambda0)-1.)*sudakov_log << " g1=" << 1./beta0*Gamma_0 + 1./beta0*Gamma_0* (1./(2.*lambda0)-1.)*sudakov_log << std::endl;
  //g1 = 1./beta0*Gamma_0;
  //std::complex<double> g1bis = 1./beta0*Gamma_0* (1./(2.*lambda0)-1.)*sudakov_log;

  //return H*std::exp(g1*log(Nb))*std::exp(g1bis*log(Nb)+g2);
  //return H*std::pow(g1*log(Nb),3)/6.;
  
  
}


std::complex<double> set_hat_NLO_cross_section (const std::complex<double> N, void *prm, bool qq)// dsigma_hat/(dM² dphi_1) 
{
  double sigma0=2./3.*alpha*pow(M_PI/MZ,2), Cf=4./3., Tf=1., dZ=0;
  std::complex<double> EulerGamma=-Psi(1.), Phat=0, Pbar=0, sigmaqg=0., sigmaqq=0.;
  Process *p = (Process *)prm;
  double tau = pow(MZ,2.)/p->sh;
  Cf=4./3.;
  dZ=Cf*(-4.+pow(M_PI,2)/3.)+3.*Cf*log(MZ/muF)+4.*Cf*log(1.-tau)*log(MZ/muF)+2.*Cf*pow(log(1.-tau),2);
  Pbar=Cf*(-pow(2. + N,-2) - pow(3. + N,-2) + 2.*PolyGamma(1, 2. + N))*2.*log(MZ/muF);
  Phat=2.*Cf*(-(EulerGamma*(-EulerGamma - 1./(1. + N) - 1./(2. + N))) - (-7. - 3.*N*(3. + N))/pow(2. + 3.*N + pow(N,2.),2) + 
     pow(M_PI,2)/6. - pow(log(1 - tau),2) + ((3. + 2.*N + 2.*EulerGamma*(1. + N)*(2. + N))*PolyGamma(0,1. + N))/
      ((1. + N)*(2. + N)) - PolyGamma(1,1. + N) + pow(PolyGamma(0,1. + N),2));
  sigmaqq=conversion_factor*sigma0*(1.+alpha_s/M_PI*(dZ+Pbar+Phat));

  sigmaqg=conversion_factor*sigma0*alpha_s/M_PI*Tf*(-(pow(N,4)+pow(N,3)-11.*pow(N,11)-19.*pow(N,4)-4.)*pow((N+1.)*(N+2.)*(N+3.),-2)+(2.*log(MZ/muF)-2.*Psi(N+1.)-2.*EulerGamma)*(pow(N,2)+3.*N+4.)/(2.*(N+1.)*(N+2.)*(N+3.)));
  
  return qq ? sigmaqq: sigmaqg;
}