// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <fstream>        // File streams               //
#include <cstring>        // In/Out streams             //
#include <cmath>      	// Mathematical functions       //
#include <math.h>
#include <complex>
#include <gsl/gsl_sf_psi.h>


// -------- Classes ----------------------------------- //
#include "messages.h"     // Message services           //
#include "LHAPDF/LHAPDF.h"// LHAPDF
#include "constants.h"   
#include "main.h"
#include "polygamma.h"
//extern double AlphaS;




void PDFInit()
{
  // Init
  const std::string pathstr = "/home/matteo/LHAPDF/share/LHAPDF";
  LHAPDF::setPaths(pathstr);

  info("PDF Paths:");
  for (const std::string& p : LHAPDF::paths()) info("  " + p);
  info("Available PDFs:");
  for (const std::string& s : LHAPDF::availablePDFSets()) info("  " + s);
}



// Function to initialise the PDF coefficients
void InitialiseCoefficients(double A[9][8])
  //const int idx_q[6] = {1, 2, 4, 7, 5, 8};  
  //const   int idx_qb[6] = {3, 6, 4, 7, 5, 8}; 
{
  // Only for CT18NLO
  // if(Proc->NLO==2 || Proc->NLO==3) pdfFit(A1min,A,muF,-1.6,-1.6,-1.6,Proc); 
  //Fitting gluon PDF...#chisq/dof = 5.08633e-05
  //Fitting valence down quark PDF...#chisq/dof = 1.03978e-07
  //Fitting valence up quark PDF...#chisq/dof = 1.29063e-07
  //Fitting sea down quark PDF...#chisq/dof = 7.17778e-08
  //Fitting bottom quark PDF...#chisq/dof = 1.7599e-08
  //Fitting sea up quark PDF...#chisq/dof = 8.35594e-09
  //Fitting charm quark PDF...#chisq/dof = 1.61266e-08

  const double coeffs[9][8] = {
    {3.36311, -1.35686, 12.01851, -9.06557, 59.82888, -209.61905, 358.78810, -214.30364},
    {0.17088, -1.32599, 4.86054, -2.74299, 27.12173, -51.56600, 58.71166, -29.80953},
    {0.15891, -1.33284, 4.04498, -1.28483, 25.16123, -21.03617, 37.34074, -32.12262},
    {0.18060, -1.32147, 4.80454, -4.01851, 19.31558, -55.90824, 71.69672, -33.26465},
    {0.19715, -1.31408, 12.42845, -7.36475, 44.88323, -146.11711, 243.91655, -145.10893},
    {0.10914, -1.34332, 16.06492, -9.85651, 68.81877, -260.94287, 489.90364, -321.74473},
    {0.19950, -1.31377, 10.66072, -5.55560, 37.16200, -128.66354, 217.17546, -124.50014},
    {0.14115, -1.33602, -0.27091, -6.04545, 15.21004, -19.76661, 13.13392, -3.53289},
    {1.0, -1.32147, 16.0, 0.0, 0.0, 0.0, 0.0, 0.0}
  };
  memcpy(A, coeffs, sizeof(coeffs));
}

// Function to handle file output streams
void OpenFiles(std::vector<std::ofstream>& files, const std::vector<std::string>& filenames)
{
  for (size_t i = 0; i < filenames.size(); ++i)
  {
    files[i].open(filenames[i]);
    if (!files[i].is_open()) error("Error opening file: " + filenames[i]);
  }
}

void CloseFiles(std::vector<std::ofstream>& files)
  { for (auto& file : files) file.close(); }

std::complex<double> Psi(std::complex<double> N)
{
  const std::complex<double> i(0., 1.);
  gsl_sf_result resR,resI;
  gsl_sf_complex_psi_e(std::real(N),std::imag(N),&resR,&resI);
  //std::cout << resR.val<< "," << resI.val << std::endl;
  return resR.val+i*resI.val;
}

std::complex<double> HarmonicNumber(std::complex<double> N)
{
  std::complex<double> res(0.,0.);
  res=(Psi(N+1.-Psi(1.)))/N;
  return res;
}

// Complex Digamma function using series expansion


int coefBinomial(int n, int k){
 
    if (k > n)
        return 0;
    if (k == 0 || k == n)
        return 1;
  
    return coefBinomial(n-1, k-1) + coefBinomial(n-1, k);
}

double derivkbis(double (* f)(double x, int i0), double x, int i0, double Eps, int k, int n)//(Not used)
{ // wiki finite difference
  //n=-1 backwar 
  //n=0 central 
  //n=1 forward
  double res=0.;
  for(int j=0;j<=k;j++)
  {
    switch (n)
    {
    case -1: res+=pow(-1, k-j)*coefBinomial(k,j)*f(x-j*Eps,i0);
      break;
    case 0: res+=pow(-1, j)*coefBinomial(k,j)*f(x+(k/2.-j)*Eps,i0);
      break;
    case 1: res+=pow(-1, j)*coefBinomial(k,j)*f(x+j*Eps,i0);
      break;
    default:
      break;
    }
  }
  return res*pow(Eps,-k);
}

double derivk(double (* f)(double x, int i0, const LHAPDF::PDF* F), double x, int i0, const LHAPDF::PDF* F ,double eps, int k)
{
    switch(k){
      case 1: return (3.*f(x-4.*eps,i0,F)-32.*f(x-3.*eps,i0,F)+168.*f(x-2.*eps,i0,F)-672.*f(x-eps,i0,F)+672.*f(x+eps,i0,F)-168.*f(x+2.*eps,i0,F)+32.*f(x+3.*eps,i0,F)-3.*f(x+4.*eps,i0,F))*pow(840.*eps,-1);
        break;
      case 2: return (-9.*f(x-4.*eps,i0,F)+128.*f(x-3.*eps,i0,F)-1008.*f(x-2.*eps,i0,F)+8064.*f(x-eps,i0,F)-14350*f(x,i0,F)+8064.*f(x+eps,i0,F)-1008.*f(x+2.*eps,i0,F)+128.*f(x+3.*eps,i0,F)-9.*f(x+4.*eps,i0,F))/5040.*pow(eps,-2);
        break;
      case 3: return (-7.*f(x-4.*eps,i0,F)+72.*f(x-3.*eps,i0,F)-338.*f(x-2.*eps,i0,F)+488.*f(x-eps,i0,F)-488.*f(x+eps,i0,F)+338.*f(x+2.*eps,i0,F)-72.*f(x+3.*eps,i0,F)+7.*f(x+4.*eps,i0,F))/240.*pow(eps,-3);
        break;
    }
}

double derive_x_k(double (* f)(double x, int i0, const LHAPDF::PDF* F), double x, int i0, const LHAPDF::PDF* F, int k)
{ //k: order of derivative
  //For derivativ:
  //n=-1 backwar 
  //n=0 central 
  //n=1 forward
  double eps=x*Eps;
  switch(k)
  {
    case 1: return f(x, i0, F)+x*derivk(f, x, i0, F,eps, 1);
      break;
    case 2: return f(x, i0, F)+3.*x*derivk(f, x, i0, F,eps, 1)+pow(x,2)*derivk(f, x, i0, F, eps, 2);
      break;
    case 3: return f(x, i0, F)+7.*x*derivk(f, x, i0, F, eps, 1)+6.*pow(x,2)*derivk(f, x, i0, F, eps, 2)+pow(x,3)*derivk(f, x, i0, F, eps, 3);
      break;
    default: return 0.;
      break;
  }

}

int indice(double x, double xmin, double xmax, int np)
{
  if(x>=xmin){
    if(x<=xmax){
      return (int) np*(x-xmin)/(xmax-xmin);
    }
  }
  else return -1;

}
long long factorial(int n) {
  if (n <= 1) {
      return 1;
  }
  return n * factorial(n - 1);
}

std::complex<double> pgamma(int m, std::complex<double> z){
  std::complex<double> res=0.;
  while(std::real(z)<14){
      res+= - pow(-1.,m)*factorial(m)*pow(z,-m-1);
      z=z+1.;
  }
  res+= pGamma(m, z);
  return res;
}