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
  if (abs(N) > 300.) return 0.;
  // Variables
  std::complex<double> Nb=N*exp(-Psi(1.)), lambda0=AlphaS/(2.*M_PI)*beta0*log(Nb);
  const double sigma0 = 2./3.*Alpha*pow(M_PI/MZ,2);// eq (A.6): Cf*alpha*pi²/(2*Mz²), the rest in the PDFs

  // Determine g1, g2, and H based on Col_imps
  const bool collinear = (id==7);
  const std::complex<double> sudakov_log = log(1.-2.*lambda0);
  const std::complex<double> sudakov_logterm = collinear ? 1.0 : (1.-2.*lambda0);
  const std::complex<double> gamma_term = collinear ? (1./beta0*gamma_0*sudakov_log) : 0.0;

  std::complex<double> g1 = 1.e+00/(2.*lambda0*beta0)*Gamma0*(2.*lambda0 + sudakov_logterm*sudakov_log);
  std::complex<double> g2 = beta1*pow(beta0,-3)/4.*Gamma0*(4.*lambda0+2.*sudakov_log+pow(sudakov_log, 2))
      - pow(beta0,-2)/4.*Gamma1*(2.*lambda0 + sudakov_log)
      + 1./beta0*Gamma0*sudakov_log*log(MZ/muR)
      - 2./beta0*Gamma0*lambda0*log(muR/muF)
      + gamma_term;
  double H = conversion_factor*sigma0*pow(MZ,2)/(2.*M_PI)*(1.+AlphaS/M_PI*4./3.* ((collinear ? (2.*pow(M_PI,2.)/3.) : (-2.*pow(log(MZ/muR),2) + 3.*log(MZ/muF) + 4.*pow(M_PI,2.)/3.)) - 4.));

  // Compute Sudakov factor and result
  return H*std::exp(g1*log(Nb)+g2);
  //Only g1
  ////incertitude = 11.6561 res1 = 44.9039 res2 = 44.4586 err1 = 0.00315261 err2 = 0.0380738////
  ///////     g1 lambda0 expansion   //////
  // Expansion order 6
  //g1=(gamma*lambda0)/beta0+(2.*gamma*pow(lambda0,2))/(3.*beta0)+(2.*gamma*pow(lambda0,3))/(3.*beta0)+(4.*gamma*pow(lambda0,4))/(5.*beta0)+(16.*gamma*pow(lambda0,5))/(15.*beta0)+(32.*gamma*pow(lambda0,6))/(21.*beta0);
////incertitude = 1.10555 res1 = 20.4549 res2 = 20.4453 err1 = 0.00371664 err2 = 0.00781316////
  
  // Expansion order 10
  //g1+=(16.*gamma*pow(lambda0,7))/(7.*beta0)+(32.*gamma*pow(lambda0,8))/(9.*beta0)+(256.*gamma*pow(lambda0,9))/(45.*beta0)+(512.*gamma*pow(lambda0,10))/(55.*beta0);
  ////incertitude = 1.3201 res1 = 20.4595 res2 = 20.471 err1 = 0.00374644 err2 = 0.00783814////

  ///////     g1 AlphaS expansion   ///////
  //std::complex<double> M=N;
  //g1=0.;
  //Order 1
  //g1=(AlphaS*Gamma0*(EulerGamma + log(M)))/(2.*Pi); 
  //// E-3 incertitude = 0.142926 res1 = 43.1407 res2 = 43.1352 err1 = 0.00297702 err2 = 0.0383873////
  // Not stable E-4 
  // With cut E-4 ////incertitude = 16.6963 res1 = 43.14 res2 = 43.0631 err1 = 0.0029665 err2 = 0.00351966//// (Non-physical work)
    //Order 2
  //g1+=(pow(AlphaS,2)*beta0*Gamma0*pow(EulerGamma +log(M),2))/(6.*pow(Pi,2));
  ////Not StableCumulativ E-3 incertitude = 5.47952 res1 = 44.6723 res2 = 44.4362 err1 = 0.00309703 err2 = 0.0429841////
  // With cut E-3 ////incertitude = 7.53888 res1 = 44.6653 res2 = 44.3805 err1 = 0.00312356 err2 = 0.0376536//// (Non-physical work)
    //Order 3
  //g1+=(pow(AlphaS,3)*pow(beta0,2)*Gamma0*pow(EulerGamma +log(M),3))/(12.*pow(Pi,3));
  ////Cumulativ E-3 incertitude = nan res1 = 44.865 res2 = -nan err1 = 0.00317653 err2 = 0.0993999////
    //Order 4
  //g1+=(pow(AlphaS,4)*pow(beta0,3)*Gamma0*pow(EulerGamma +log(M),4))/(20.*pow(Pi,4));
  ////Cumulativ E-3 incertitude = nan res1 = 44.8892 res2 = -nan err1 = 0.00292435 err2 = 268.664////
    //Order 5
  //g1+=(pow(AlphaS,5)*pow(beta0,4)*Gamma0*pow(EulerGamma +log(M),5))/(30.*pow(Pi,5));
  ////Not Stable Cumulativ E-3 incertitude = 11.3134 res1 = 44.898 res2 = 44.4495 err1 = 0.00285925 err2 = 0.0395372////
    //Order 6
  //g1+=(pow(AlphaS,6)*pow(beta0,5)*Gamma0*pow(EulerGamma +log(M),6))/(42.*pow(Pi,6));
  ////Not Stable Cumulativ E-3 incertitude = 12.5653 res1 = 44.9066 res2 = 44.3963 err1 = 0.00284556 err2 = 0.0405168////
    //Order 7
  //g1+=(pow(AlphaS,7)*pow(beta0,6)*Gamma0*pow(EulerGamma +log(M),7))/(56.*pow(Pi,7)); 
    //Order 8
  //g1+=(pow(AlphaS,8)*pow(beta0,7)*Gamma0*pow(EulerGamma+log(M),8))/(72.*pow(Pi,8)); 
    //Order 9
  //g1+=(pow(AlphaS,9)*pow(beta0,8)*Gamma0*pow(EulerGamma+log(M),9))/(90.*pow(Pi,9)); 
    //Order 10
  //g1+=(pow(AlphaS,10)*pow(beta0,9)*Gamma0*pow(EulerGamma +log(M),10))/(110.*pow(Pi,10));
  ////incertitude = 101.533 res1 = 44.9063 res2 = 44.4407 err1 = 0.00292429 err2 = 0.0035322////
  // cut 200////incertitude = 96.582 res1 = 44.903 res2 = 44.4492 err1 = 0.0030457 err2 = 0.00357741////

  //return H*std::exp(g1*log(Nb));

}



std::complex<double> set_hat_NLO_cross_section (const std::complex<double> N, void *prm, bool qq)// dsigma_hat/(dM² dphi_1) 
{
  double sigma0=2./3.*Alpha*pow(M_PI/MZ,2);
  if(qq)
  {
  std::complex<double> EulerGamma=-Psi(1.), Phat=0, Pbar=0;
  Process *p = (Process *)prm;
  double tau = pow(MZ,2.)/p->sh;
  double dZ=Cf*(-4.+pow(M_PI,2)/3.)+3.*Cf*log(MZ/muF)+4.*Cf*log(1.-tau)*log(MZ/muF)+2.*Cf*pow(log(1.-tau),2);
  Pbar=Cf*(-pow(2. + N,-2) - pow(3. + N,-2) + 2.*PolyGamma(1, 2. + N))*2.*log(MZ/muF);
  Phat=2.*Cf*(-(EulerGamma*(-EulerGamma - 1./(1. + N) - 1./(2. + N))) - (-7. - 3.*N*(3. + N))/pow(2. + 3.*N + pow(N,2.),2) + 
     pow(M_PI,2)/6. - pow(log(1 - tau),2) + ((3. + 2.*N + 2.*EulerGamma*(1. + N)*(2. + N))*PolyGamma(0,1. + N))/
      ((1. + N)*(2. + N)) - PolyGamma(1,1. + N) + pow(PolyGamma(0,1. + N),2));

  return conversion_factor*sigma0*(1.+AlphaS/M_PI*(dZ+Pbar+Phat));
  }

  else // qg case
  {
  std::complex<double> fqg=0.,Phatqg1, Phatqg2;
  fqg=(5. + 4.*N)/(2.*(6. + 11.*N + 6.*pow(N,2) + pow(N,3)));
  Phatqg1=-0.5*(26. + N*(60. + N*(63. + 5.*N*(6. + N))) + 2.*(1. + N)*(2. + N)*(3. + N)*(4. + N*(3. + N))*HarmonicNumber(N))/
    (pow(1. + N,2)*pow(2. + N,2)*pow(3. + N,2));
  Phatqg2=-(1./(2. + N)) + 1./(3. + N) + 1./(2. + 2.*N);

  return conversion_factor*sigma0*AlphaS/M_PI*Tf*(Phatqg1+Phatqg2*2.*log(MZ/muF)+fqg);
  }
}