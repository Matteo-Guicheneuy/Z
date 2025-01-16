// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <complex>                                      //
#include <vector>                                       //*
// -------- Classes ----------------------------------- //
#include "messages.h"     // Message services           //
#include "process.h"      // Process definition         //
#include "constants.h"    // Constants                  //
#include "luminosity.h"
#include "settings.h"
// -------- Functions --------------------------------- //

std::complex<double> set_hat_cross_section(const std::complex<double>&, const int&);
std::complex<double> set_hat_NLO_cross_section (const std::complex<double> N, void *prm, bool qq);// dsigma_hat/(dM² dphi_1) 

//void EvolvePDF(const std::complex<double>&, std::complex<double>*, std::complex<double>*, std::complex<double>&);


//------------------------------------------------------------//
//----------------Resummed computation------------------------//
//------------------------------------------------------------//
double Resummed(double *x, size_t dim, void *prm,  bool usePDFtrick = false, int k=1, int icut = -1 )
{   // Process information
    Process *p = (Process *)prm;
  // For ic over icut we remove one part of the integral and use the default contour.
  // Parameters and constants
  const double C = 2.015, Phi = ((usePDFtrick && p->ic<=icut) ? 1.0 : 1.7)*M_PI/2., a=(usePDFtrick && p->ic<=icut) ? 0.7 : -1.;//0.99999999
  const std::complex<double> i(0., 1.);
  std::complex<double> N, sigma, jac, g_12, prefactor;

  // Tau computation and change of variables in Mellin space
  double tau = pow(MZ, 2) / p->sh;
  N = C + (cos(Phi)+i*sin(Phi))*tan(M_PI*x[0]/2.);
  jac = 0.5*(cos(Phi)+i*sin(Phi))*(1.+pow(tan(M_PI*x[0]/2.), 2));
  
  if (usePDFtrick)
  { 
    double xa=(p->ic<=icut) ? x[1]: tau/pow(tau,x[1]), xb=(p->ic<=icut) ? x[2] : tau/pow(tau/xa,x[2])/xa;
    if(p->ic>icut) jac *= xa*xb*pow(log(tau),2.)*x[1];

    // If using PDF tricks, we use a different method for setting PDFs -> luminosity
    double fAB = Luminosity(p->xs_id, xa, xb, a, p->usePDFAnzat, true, k);  // Use PDF fit
    //double fABqg = Luminosity(p->xs_id, xa, xb, a, p->usePDFAnzat, false, k);
    // Adjust N for second PDF trick (xs_id 4 uses an offset) + prefactor
        
    //if( p->xs_id==4) { prefactor=pow(tau/xa/xb,-N-a)/xa/xb; }
    if( p->xs_id==4 || p->xs_id==5) { prefactor=pow(tau/xa/xb,-N-a)/xa/xb;}
    else prefactor=tau*pow(tau/xa/xb,-N);

    //Fake PDFs for the PDF trick
    g_12=pow(N+a,-2.*k);
    // Definiton of sigma_hat
    //if(p->xs_id==1) sigma=set_hat_NLO_cross_section(N, prm,true);
    if(p->xs_id==4 || p->xs_id==5 || p->xs_id==6) sigma = 2.*M_PI*pow(MZ,-2)*set_hat_cross_section(N+a+1., p->xs_id);

    // Fake PDFs in Mellin space for the PDF trick
    double res= std::imag(jac*prefactor*fAB*g_12*sigma); 
    //if(abs(res)>1e5) std::cout<<"N= " << N << std::setprecision(10) << " xa=" << xa << " xb=" << xb << " x[0]= " << x[0] << " prefactor = " << prefactor << " fAB = " << fAB << " g_12 = " << g_12 << " sigma = " << sigma << " jac = " << jac <<" res = " << res << std::endl;
    return std::imag(jac*prefactor*fAB*g_12*sigma);//+(p->xs_id==1 ? std::imag(jac*prefactor*fABqg*g_12*set_hat_NLO_cross_section(N, prm,false)): 0.);
  }
  else 
  {
    if(p->xs_id==2 || p->xs_id==3 || p->xs_id==6)
    {
    // If not using PDF tricks, use normal PDF setting
    // Evaluate sigma based on NLO conditions
    if (p->xs_id==2) sigma = 2.*M_PI*pow(MZ,-2)*set_hat_cross_section(N+1., p->xs_id);
    else if (p->xs_id==3) sigma =set_hat_NLO_cross_section(N, prm, true);
    else if (p->xs_id==6)
    {
      //EvolvePDF(N+1., q, qbar, g);
      sigma = 2.*M_PI*pow(MZ,-2)*set_hat_cross_section(N+1., p->xs_id);
    }
    // Calculate fAB the luminosity

    // PDF Anzat implementation
    std::complex<double> fAB= MLuminosity(N+1., p->usePDFAnzat);
    // Calculate sigma and final result
    //double res= std::imag(jac*pow(tau,-N)*fAB*sigma);
    return std::imag(jac*pow(tau,-N)*fAB*sigma);  // Integrand
  }
  else if(p->xs_id==1)
  {
     // Kinematics
    std::complex<double> fAB = Luminosity(p->xs_id, x[1], x[2], a, p->usePDFAnzat, true, k);  // Use PDF fit
    prefactor=pow(tau/x[1]/x[2],-N)/x[1]/x[2]; 
    g_12=pow(N,-2.);
    sigma = 2.*M_PI*pow(MZ,-2)*set_hat_cross_section(N+1., p->xs_id);
    return std::imag(jac*prefactor*fAB*g_12*sigma);
  }
  }
}
// Integrand (only needed ones)
double Tot(double *x, size_t dim, void *prm)
{
  Process *p = (Process *)prm;
  // Dispatch computation based on xs_id
  switch (p->xs_id)
  {
    case 1: case 2: case 3: return Resummed(x, dim, prm); // No collinear improvement
    case 4:  return Resummed(x, dim, prm, true, 2); // NLO+PDF tricks // Second PDF trick // First PDF trick
    case 5: return Resummed(x, dim, prm, true, 2);
    default:
        error("Invalid Process ID: ID = " + std::to_string(p->xs_id));
        return 0.; // Unreachable, but avoids compiler warnings
  }
}
