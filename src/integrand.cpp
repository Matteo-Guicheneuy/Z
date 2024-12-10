// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <complex>                                      //
// -------- Classes ----------------------------------- //
#include "messages.h"     // Message services           //
#include "process.h"      // Process definition         //
#include "constants.h"    // Constants                  //
#include "luminosity.h"
// -------- Functions --------------------------------- //

std::complex<double> set_hat_cross_section(const std::complex<double>&, const int&);
std::complex<double> set_hat_NLO_cross_section (const std::complex<double> N, void *prm);// dsigma_hat/(dM² dphi_1) 

//void EvolvePDF(const std::complex<double>&, std::complex<double>*, std::complex<double>*, std::complex<double>&);


//------------------------------------------------------------//
//----------------Resummed computation------------------------//
//------------------------------------------------------------//
double Resummed(double *x, size_t dim, void *prm,  bool usePDFtrick = false)
{
  // Parameters and constants
  const double C = 2.015, Phi = 1.7*M_PI/2., a=-1.;
  const std::complex<double> i(0., 1.);
  double xa, xb;
  std::complex<double> N, sigma, jac, g_12, prefactor;

  // Process information
  Process *p = (Process *)prm;

  // Tau computation and change of variables in Mellin space
  double tau = pow(MZ, 2) / p->sh;
  N = C + (cos(Phi)+i*sin(Phi))*tan(M_PI*x[0]/2.);
  jac = 0.5*(cos(Phi)+i*sin(Phi))*(1.+pow(tan(M_PI*x[0]/2.), 2));

  // PDF initialisation
  if (usePDFtrick)
  {
    // Kinematics
    xa=tau/pow(tau,x[1]);
    xb=tau/pow(tau/xa,x[2])/xa;
    jac *= xa*xb*pow(log(tau),2.)*x[1];

    // If using PDF tricks, we use a different method for setting PDFs -> luminosity
    std::complex<double> fAB = Luminosity(5-p->xs_id, xa, xb, a, p->usePDFAnzat);  // Use PDF fit
    // Adjust N for second PDF trick (xs_id 4 uses an offset) + prefactor
    if( p->xs_id==4) { N+=a+1; prefactor=pow(tau/xa/xb,-N-a)/xa/xb; }
    else prefactor=tau*pow(tau/xa/xb,-N);

    //Fake PDFs for the PDF trick
    g_12=pow(N+a,-2.);

    // Definiton of sigma_hat
    if(p->xs_id==1) sigma=set_hat_NLO_cross_section(N, prm);
    if(p->xs_id==4 || p->xs_id==5) sigma = 2.*M_PI*pow(MZ,-2)*set_hat_cross_section(N, p->xs_id);

    // Fake PDFs in Mellin space for the PDF trick
    return std::imag(jac*prefactor*fAB*g_12*sigma);
  }
  else
  {
    // If not using PDF tricks, use normal PDF setting
    std::complex<double> q[5], g, qbar[5], fAB;
    set_pdfN(N+1., q, qbar, g, p->usePDFAnzat);  // Set PDFs with N

    // Evaluate sigma based on NLO conditions
    if (p->xs_id==2) sigma = 2.*M_PI*pow(MZ,-2)*set_hat_cross_section(N+1., p->xs_id);
    else if (p->xs_id==6)
    {
      // EvolvePDF(N+1., q, qbar, g);
      sigma = 2.*M_PI*pow(MZ,-2)*set_hat_cross_section(N+1., p->xs_id);
    }

    // Calculate fAB the luminosity
    for (int i=0; i<Nflav; i++) fAB+= (i%2==0 ? gd : gu)*q[i]*qbar[i]*2.;
    //remove charm quark contribution
    fAB-= gu*q[3]*qbar[3]*2.;

    // PDF Anzat implementation
    if(p->usePDFAnzat) fAB=q[0]*qbar[0];
    // Calculate sigma and final result
    return std::imag(jac*pow(tau,-N)*fAB*sigma);  // Integrand
  }
}

// Integrand (only needed ones)
double Tot(double *x, size_t dim, void *prm)
{
  Process *p = (Process *)prm;
  // Dispatch computation based on xs_id
  switch (p->xs_id)
  {
    case 2: return Resummed(x, dim, prm); // No collinear improvement
    case 1: case 4: case 5: return Resummed(x, dim, prm, true); // NLO+PDF tricks // Second PDF trick // First PDF trick
    default:
        error("Invalid Process ID: ID = " + std::to_string(p->xs_id));
        return 0.; // Unreachable, but avoids compiler warnings
  }
}
