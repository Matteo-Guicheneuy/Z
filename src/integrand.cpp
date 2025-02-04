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
std::complex<double> set_hat_cross_section_Exp(const std::complex<double>&, const int&);
std::complex<double> set_hat_NLO_cross_section (const std::complex<double> N, void *prm, bool qq);// dsigma_hat/(dM² dphi_1) 

//void EvolvePDF(const std::complex<double>&, std::complex<double>*, std::complex<double>*, std::complex<double>&);


//------------------------------------------------------------//
//----------------Resummed computation------------------------//
//------------------------------------------------------------//
double Resummed(double *x, size_t dim, void *prm,  bool usePDFtrick = false, int k=1, int icut = -1 )
{ // Process information
  Process *p = (Process *)prm;
  //(Not Useful) For ic under icut we add an unphysical part to the integral and use another countour.I do not why it works
  //(Not Useful) all the (Not Useful) following have to be remove to get this implementation
  // Parameters and constants
  double C = 2.015, Phi = 1.7*M_PI/2., a=-1.;
  //(Not Useful) if(usePDFtrick && p->ic<=icut)Phi= 1.0*M_PI/2., a=0.7;

  const std::complex<double> i(0., 1.);
  std::complex<double> N, sigma, jac, g_12, prefactor;

  // Tau definition
  double tau = pow(MZ, 2) / p->sh;
  // Mellin variable:
  N = C + (cos(Phi)+i*sin(Phi))*tan(M_PI*x[0]/2.);
  // jacobian associated to this contour in complexe plane
  jac = 0.5*(cos(Phi)+i*sin(Phi))*(1.+pow(tan(M_PI*x[0]/2.), 2));
  
  if (usePDFtrick) // PDF trick
  { 
    double xa=tau/pow(tau,x[1]), xb=tau/pow(tau/xa,x[2])/xa;
    jac *= xa*xb*pow(log(tau),2.)*x[1];
    //(Not Useful)if(usePDFtrick && p->ic<=icut) xa=x[1], xb=x[2], jac = 0.5*(cos(Phi)+i*sin(Phi))*(1.+pow(tan(M_PI*x[0]/2.), 2)); 
    

    // If using PDF tricks, we use a different method for setting PDFs -> luminosity
    double fAB = Luminosity(p->xs_id, xa, xb, a, p->usePDFAnzat, true, k);
    double fABqg = Luminosity(p->xs_id, xa, xb, a, p->usePDFAnzat, false, k);

    // Prefactor
    if( p->xs_id==2 || p->xs_id==3 || p->xs_id==5) { prefactor=pow(tau/xa/xb,-N-a)/xa/xb;}
    else prefactor=tau*pow(tau/xa/xb,-N);

    //Fake PDFs for the PDF trick
    g_12=pow(N+a,-2.*k);

    // Definiton of sigma_hat
    if(p->xs_id==2 || p->xs_id==5) sigma = 2.*M_PI*pow(MZ,-2)*set_hat_cross_section(N+a+1., p->xs_id);
    if(p->xs_id==3) sigma = 2.*M_PI*pow(MZ,-2)*set_hat_cross_section_Exp(N+a+1., p->xs_id);

    return std::imag(jac*prefactor*fAB*g_12*sigma) + (p->xs_id==1 ? std::imag(jac*prefactor*fABqg*g_12*set_hat_NLO_cross_section(N, prm,false)): 0.);
  }
  else // Resum with Mellin PDF
  {
    // If not using PDF tricks, use normal PDF setting in Mellin space
    // Evaluate sigma based on NLO conditions
    if (p->xs_id==4) sigma = 2.*M_PI*pow(MZ,-2)*set_hat_cross_section(N+1., p->xs_id);
    // Calculate fAB the luminosity
    std::complex<double> fAB= MLuminosity(N+1., p->usePDFAnzat);
    // Calculate sigma and final result
    return std::imag(jac*pow(tau,-N)*fAB*sigma);  // Integrand
  }
}

double NLO_cross_section(double *x, size_t dim, void *prm)
{   
  Process *p = (Process *)prm;
  double xa, xb, s, jac, res=0.;
  // PDFs
   double fAB=0;

  ////      BORN      ////
   // Computation of the kinematical variables
   Kinematics(xa,xb,jac,s,x,p);
   SetCouplings(0,p,fAB,xa,xb); //Sets PDFs with their appropriate u/d coupling, summed and symmetrized
   
   res+=pow(M_PI,2.)*2./3.*Alpha*jac*fAB/xa/xb/p->sh*conversion_factor;

  ////       REAL      ////
   // Computation of the kinematical variables
   Kinematics2(xa,xb,jac,s,x,p);
   double z=pow(MZ,2.)/p->sh/xa/xb; //z variable to express functions more easily
   SetCouplings(3,p,fAB,xa,xb); //Sets PDFs with their appropriate u/d coupling, summed and symmetrized

   // Couplings and scale setting
   double prefqg=Alpha*AlphaS*M_PI/3./p->sh; //prefactors for real (qg->zq) NLO subprocesses
   double integqg=(1+7*z)*(1-z)/4.+(pow(z,2.)+pow(1-z,2.))/2.*log((1-z)*MZ/muF/sqrt(z))*2; //Function to integrate for real NLO subprocesses needing 2 variables
 
   res+=prefqg*integqg*fAB*jac/xa/xb*conversion_factor;
   
   ////    VIRTUAL 1D    ////
   // Computation of the kinematical variables
   Kinematics(xa,xb,jac,s,x,p);
   SetCouplings(0,p,fAB,xa,xb); //Sets PDFs with their appropriate u/d coupling, summed and symmetrized

   // Couplings and scale setting
   double prefqq=Alpha*AlphaS*M_PI*8./9./p->sh; //Prefactors for virtual (qqbar->z) NLO subprocesses
   double integqq=3.*log(MZ/muF)+pow(M_PI,2.)/3.-4; //Function to integrate for virtual NLO subprocesses needing 2 variables

   res+=prefqq*integqq*fAB*jac/xa/xb*conversion_factor;

  ////    VIRTUAL 2D     ////
   // Computation of the kinematical variables
   Kinematics2(xa,xb,jac,s,x,p);
   SetCouplings(0,p,fAB,xa,xb); //Sets PDFs with their appropriate u/d coupling, summed and symmetrized

   // Couplings and scale setting
   prefqq=Alpha*AlphaS*M_PI*8./9./p->sh; //Prefactors for virtual (qqbar->z) NLO subprocesses
   integqq=-(1+z*z)*log(z)/(1-z); //Function to integrate for virtual NLO subprocesses needing 2 variables

   res+=prefqq*integqq*fAB*jac/xa/xb*conversion_factor;

  ////     VIRTUAL PLUS      ////
   // Computation of the kinematical variables
   double jac1, y, jac2;
   Kinematics2Plus(y,z,jac1,s,x,p);
   jac2=-log(pow(MZ,2.)/p->sh); // jacobian for 1D +distribution leftovers

   // Couplings and scale setting
   double fAB1=0, fAB2=0;
   prefqq=Alpha*AlphaS*M_PI*8./9./p->sh; // Prefactor for q qbar > z NLO subprocess

   SetCouplings(1,p,fAB1,y,z); //Sets PDFs with their appropriate u/d coupling, summed and symmetrized
   SetCouplings(2,p,fAB2,y,z); //Sets PDFs with their appropriate u/d coupling, summed and symmetrized

   res+=prefqq*conversion_factor*(jac1*fAB1+jac2*fAB2);

  return res;

}
// Integrand (only needed ones)
double Tot(double *x, size_t dim, void *prm)
{
  Process *p = (Process *)prm;
  // Dispatch computation based on xs_id
  switch (p->xs_id)
  {
    case 1: return NLO_cross_section(x, dim, prm);
    case 4: return Resummed(x, dim, prm); 
    case 2: case 3: case 5: return Resummed(x, dim, prm, true, 2); 
    default:
        error("Invalid Process ID: ID = " + std::to_string(p->xs_id));
        return 0.; // Unreachable, but avoids compiler warnings
  }
}
