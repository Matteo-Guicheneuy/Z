// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <string>     // Strings                        //
#include <vector>
#include <complex>
// -------- Classes ----------------------------------- //
#include "messages.h"     // Message services           //
#include "luminosity.h"   // Luminosities and PDFs      //
#include "constants.h"    // Constants                  //
#include "coefficients.h" //
#include "settings.h"     // Setting                    //
#include "process.h"
#include "main.h"
#include <iostream>


double f_LHAPDF(double x, int i0, const LHAPDF::PDF* F)
{ 
  double res=0.;
  if(x > 0 && x < 1) res=F->xfxQ2(i0,x,muF*muF)/x;
  return res;
}

double Set_f_LHAPDF(double &x, int i0, const LHAPDF::PDF* F, int k)
{ 
  if(k==0) return f_LHAPDF(x, i0, F);
  else return x*derive_x_k(f_LHAPDF,x,i0,F ,k);
}

void set_pdf_LHAPDF(const LHAPDF::PDF* F , double &x, std::vector<double>& q, std::vector<double>& qbar, double &g, const bool usePDFAnzat, int k, bool num)
{
  // Quarks PDG list
  std::vector<int>  idx_q = usePDFAnzat ? std::vector<int>{8} : std::vector<int>{1, 2, 3, 4, 5};// d , u , s , c , b
  std::vector<int>  idx_qb = usePDFAnzat ? std::vector<int>{8} : std::vector<int>{-1, -2, -3, -4, -5};// dbar , ubar , sbar , cbar , bbar
  //PDFs  
  g = Set_f_LHAPDF(x, 0, F, k);
  for (int j=0; j<Nflav; j++)
  {
    qbar.push_back( Set_f_LHAPDF(x, idx_qb[j], F, k));  // qbar for b, c, s, u
    q.push_back( Set_f_LHAPDF(x, idx_q[j], F, k));    // q for b, c, s 
  }
}
//PDF analitic function
double f(double x, int i0, const LHAPDF::PDF* F)
{
  return A[i0][0]*pow(x,A[i0][1])*pow(1.-x,A[i0][2])*(1.+A[i0][3]*sqrt(x)+A[i0][4]*x+A[i0][5]*pow(x,1.5)+A[i0][6]*pow(x,2)+A[i0][7]*pow(x, 2.5));
}
//PDF Anzat analitic function
double fAnzat(double x, int i0, const LHAPDF::PDF* F)
{
  return A[i0][0]*pow(x,A[i0][1])*pow(1.-x,A[i0][2]);
}

// Derivative of the PDF function at order k numericaly
double F_real(const double &x, int i0, const LHAPDF::PDF* F, const bool usePDFAnzat, int k)
{
  if (x>1-5.*Eps*x) return 0;
  else if(x<5.*Eps*x) return 0;
  if(k==0) return usePDFAnzat ? fAnzat(x, i0, F) : f(x, i0, F);
  else return usePDFAnzat ? x*derive_x_k(fAnzat,x,i0,F,k) : x*derive_x_k(f,x,i0,F,k);
}

// Derivative of the PDF function at order k analiticaly
double F_real_bis(const double &x, int i0, const LHAPDF::PDF* F, const bool usePDFAnzat, int k)
{
  // Derivative of PDF in real space
  if(usePDFAnzat){
    switch (k){
       case 0 :return A[i0][0]*pow(x,A[i0][1])*pow(1.-x,A[i0][2]);
       break;
       case 1 :return A[i0][0]*pow(x,A[i0][1]+1.)*pow(1.-x,A[i0][2])*(A[i0][1]+1.-x*A[i0][2]/(1.-x));
       break;
       //case 2 :return x*(pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])-1.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]);
       case 2 :return x*(pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])-1.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]+x*(pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])-2.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]+pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]));
       break;
       case 3 :return x*(pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])-1.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]+x*(pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])-2.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]+pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2])+x*(2.*pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])-4.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]+2.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]+x*(pow(1.-1.*x,A[i0][2])*pow(x,-2.+A[i0][1])*A[i0][0]*(-1.+A[i0][1])*A[i0][1]*(1.+A[i0][1])-3.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*A[i0][2]+3.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(-1.+A[i0][2])*A[i0][2]-1.*pow(1.-1.*x,-3.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-2.+A[i0][2])*(-1.+A[i0][2])*A[i0][2])));
       break;
       case 4: return x*(pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])-1.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]+x*(pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])-2.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]+pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2])+x*(2.*pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])-4.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]+2.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]+x*(pow(1.-1.*x,A[i0][2])*pow(x,-2.+A[i0][1])*A[i0][0]*(-1.+A[i0][1])*A[i0][1]*(1.+A[i0][1])
-3.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*A[i0][2]+3.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(-1.+A[i0][2])*A[i0][2]-1.*pow(1.-1.*x,-3.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-2.+A[i0][2])*(-1.+A[i0][2])*A[i0][2]))+x*(4.*pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])-8.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]+4.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]+2.*x*(pow(1.-1.*x,A[i0][2])*pow(x,-2.+A[i0][1])*A[i0][0]*(-1.+A[i0][1])*A[i0][1]*(1.+A[i0][1])-3.*pow(1.-1.*
x,-1.+A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*A[i0][2]+3.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(-1.+A[i0][2])*A[i0][2]-1.*pow(1.-1.*x,-3.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-2.+A[i0][2])*(-1.+A[i0][2])*A[i0][2])+x*(3.*pow(1.-1.*x,A[i0][2])*pow(x,-2.+A[i0][1])*A[i0][0]*(-1.+A[i0][1])*A[i0][1]*(1.+A[i0][1])-9.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*A[i0][2]+9.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(-1.+A[i0][2])*A[i0][2]-3.*pow(1.-1.*x,-3.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-2.+A[i0][2])*(-1.+A[i0][2])*A[i0][2]+x*(pow(1.-1.*x,A[i0][2])*pow(x,-3.+A[i0][1])*A[i0][0]*(-2.+
A[i0][1])*(-1.+A[i0][1])*A[i0][1]*(1.+A[i0][1])-4.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,-2.+A[i0][1])*A[i0][0]*(-1.+A[i0][1])*A[i0][1]*(1.+A[i0][1])*A[i0][2]+6.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*(-1.+A[i0][2])*A[i0][2]-4.*pow(1.-1.*x,-3.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(-2.+A[i0][2])*(-1.+A[i0][2])*A[i0][2]+pow(1.-1.*x,-4.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-3.+A[i0][2])*(-2.+A[i0][2])*(-1.+A[i0][2])*A[i0][2]))));
       break;
    }
  }
  else{
    switch (k)
    {
    case 0: A[i0][0]*pow(x,A[i0][1])*pow(1.-x,A[i0][2])*(1.+A[i0][3]*sqrt(x)+A[i0][4]*x+A[i0][5]*pow(x,1.5)+A[i0][6]*pow(x,2)+A[i0][7]*pow(x, 2.5));
      break;
    case 1: return A[i0][0]*pow(x,A[i0][1]+1.)*pow(1.-x,A[i0][2])*((A[i0][1]+1.-x*A[i0][2]/(1.-x))*(1.+A[i0][3]*sqrt(x)+A[i0][4]*x+A[i0][5]*pow(x,1.5)+A[i0][6]*pow(x,2.)+A[i0][7]*pow(x,2.5))+0.5*A[i0][3]*sqrt(x)+A[i0][4]*x+1.5*A[i0][5]*pow(x,1.5)+2.*A[i0][6]*pow(x,2.)+2.5*A[i0][7]*pow(x,2.5));
      break;
    case 2: return x*(pow(1.-x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(A[i0][3]/(2.*sqrt(x))+A[i0][4]+(3.*sqrt(x)*A[i0][5])/2.+2.*x*A[i0][6]+(5.*pow(x,1.5)*A[i0][7])/2.)+pow(1.-x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(1.+sqrt(
x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-pow(1.-x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+x*(pow(1.-x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(
-0.25*A[i0][3]/pow(x,1.5)+(3.*A[i0][5])/(4.*sqrt(x))+2*A[i0][6]+(15.*sqrt(x)*A[i0][7])/4.)+2*pow(1.-x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(A[i0][3]/(2.*sqrt(x))+A[i0][4]+(3.*sqrt(x)*A[i0][5])/2.+2*x*A[i0][6]+(5*pow(x,1.5)*A[i0][7])/2.)-2.*pow(1.-x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*(A[i0][3]/(2.*sqrt(x))+A[i0][4]+(3.*sqrt(x)*
A[i0][5])/2.+2*x*A[i0][6]+(5.*pow(x,1.5)*A[i0][7])/2.)+pow(1.-x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-2.*pow(1-x,-1+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0]
[2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+pow(1.-x,-2+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])));
      break;
    case 3: return x*(pow(1.-1.*x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-1.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+x*(pow(1.-1.*x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])+2.*pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])-2.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-2.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7]))+x*(2.*pow(1.-1.*x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])+4.*pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])-4.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+2.*pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-4.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+2.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+x*(pow(1.-1.*x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*((0.375*A[i0][3])/pow(x,2.5)-(0.375*A[i0][5])/pow(x,1.5)+(1.875*A[i0][7])/sqrt(x))+3.*pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])-3.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])+3.*pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])-6.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+3.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+pow(1.-1.*x,A[i0][2])*pow(x,-2.+A[i0][1])*A[i0][0]*(-1.+A[i0][1])*A[i0][1]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-3.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+3.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-1.*pow(1.-1.*x,-3.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7]))));
      break;
      case 4: return x*(pow(1.-1.*x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-1.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+x*(pow(1.-1.*x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])+2.*pow(1.-1.*x,A[i0][2])*
pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])-2.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-2.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7]))+x*(2.*pow(1.-1.*x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+
2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])+4.*pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])-4.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+2.*pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-4.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+
2.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+x*(pow(1.-1.*x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*((0.375*A[i0][3])/pow(x,2.5)-(0.375*A[i0][5])/pow(x,1.5)+(1.875*A[i0][7])/sqrt(x))+3.*pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])-3.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])+3.*pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])-6.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+3.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*
pow(x,1.5)*A[i0][7])+pow(1.-1.*x,A[i0][2])*pow(x,-2.+A[i0][1])*A[i0][0]*(-1.+A[i0][1])*A[i0][1]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-3.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+3.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-1.*pow(1.-1.*x,-3.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-2.+A[i0][2])*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])))+x*(4.*pow(1.-1.*x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])+
8.*pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])-8.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+4.*pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-8.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+4.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+
pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+2.*x*(pow(1.-1.*x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*((0.375*A[i0][3])/pow(x,2.5)-(0.375*A[i0][5])/pow(x,1.5)+(1.875*A[i0][7])/sqrt(x))+3.*pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])-3.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])+3.*pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])-6.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*
A[i0][7])+3.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+pow(1.-1.*x,A[i0][2])*pow(x,-2.+A[i0][1])*A[i0][0]*(-1.+A[i0][1])*A[i0][1]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-3.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+3.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-1.*pow(1.-1.*x,-3.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-2.+A[i0][2])*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+
pow(x,2.5)*A[i0][7]))+x*(3.*pow(1.-1.*x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*((0.375*A[i0][3])/pow(x,2.5)-(0.375*A[i0][5])/pow(x,1.5)+(1.875*A[i0][7])/sqrt(x))+9.*pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])-9.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])+9.*pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])-18.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+9.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+3.*pow(1.-1.*x,A[i0][2])*pow(x,-2.+A[i0][1])*A[i0][0]*(-1.+
A[i0][1])*A[i0][1]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-9.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+9.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-3.*pow(1.-1.*x,-3.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-2.+A[i0][2])*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+x*(pow(1.-1.*x,A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*((-0.9375*A[i0][3])/pow(x,3.5)+(0.5625*A[i0][5])/pow(x,2.5)-(0.9375*A[i0][7])/pow(x,1.5))+4.*pow(1.-1.*x,A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*((0.375*A[i0][3])/pow(x,2.5)-(0.375*A[i0][5])/pow(x,1.5)+(1.875*A[i0][7])/sqrt(x))-4.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*A[i0][2]*((0.375*A[i0][3])/pow(x,2.5)-(0.375*A[i0][5])/pow(x,1.5)+(1.875*A[i0][7])/sqrt(x))+6.*pow(1.-1.*x,A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+
2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])-12.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*A[i0][2]*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])+6.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-1.+A[i0][2])*A[i0][2]*((-0.25*A[i0][3])/pow(x,1.5)+(0.75*A[i0][5])/sqrt(x)+2.*A[i0][6]+3.75*sqrt(x)*A[i0][7])+4.*pow(1.-1.*x,A[i0][2])*pow(x,-2.+A[i0][1])*A[i0][0]*(-1.+A[i0][1])*A[i0][1]*(1.+A[i0][1])*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])-12.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+12.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(-1.+A[i0][2])*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+
2.5*pow(x,1.5)*A[i0][7])-4.*pow(1.-1.*x,-3.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-2.+A[i0][2])*(-1.+A[i0][2])*A[i0][2]*((0.5*A[i0][3])/sqrt(x)+A[i0][4]+1.5*sqrt(x)*A[i0][5]+2.*x*A[i0][6]+2.5*pow(x,1.5)*A[i0][7])+pow(1.-1.*x,A[i0][2])*pow(x,-3.+A[i0][1])*A[i0][0]*(-2.+A[i0][1])*(-1.+A[i0][1])*A[i0][1]*(1.+A[i0][1])*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-4.*pow(1.-1.*x,-1.+A[i0][2])*pow(x,-2.+A[i0][1])*A[i0][0]*(-1.+A[i0][1])*A[i0][1]*(1.+A[i0][1])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+6.*pow(1.-1.*x,-2.+A[i0][2])*pow(x,-1.+A[i0][1])*A[i0][0]*A[i0][1]*(1.+A[i0][1])*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])-4.*pow(1.-1.*x,-3.+A[i0][2])*pow(x,A[i0][1])*A[i0][0]*(1.+A[i0][1])*(-2.+A[i0][2])*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])+pow(1.-1.*x,-4.+A[i0][2])*pow(x,1.+A[i0][1])*A[i0][0]*(-3.+A[i0][2])*(-2.+A[i0][2])*(-1.+A[i0][2])*A[i0][2]*(1.+sqrt(x)*A[i0][3]+x*A[i0][4]+pow(x,1.5)*A[i0][5]+pow(x,2)*A[i0][6]+pow(x,2.5)*A[i0][7])))));
      break;
    default:
      break;
    }
  }
}

void set_pdf_fit(const LHAPDF::PDF* F, double &x,std::vector<double>& q, std::vector<double>& qbar, double &g, const bool usePDFAnzat, int k, bool num)
{ 
  // Quarks PDG list
  std::vector<int>  idx_q = usePDFAnzat ? std::vector<int>{8} : std::vector<int>{1, 2, 4, 7, 5};// d , u , s , c , b
  std::vector<int>  idx_qb = usePDFAnzat ? std::vector<int>{8} : std::vector<int>{3, 6, 4, 7, 5};// dbar , ubar , sbar , cbar , bbar
  double (*func)(const double&, int, const LHAPDF::PDF* F,const bool, int);
  if (num)
    func = F_real;
  else
    func = F_real_bis;
  //PDFs  
  g = func(x, 0, F,usePDFAnzat, k);
  for (int j=0; j<Nflav; j++)
  {
    qbar.push_back( func(x, idx_qb[j], F,usePDFAnzat, k));  // qbar for b, c, s, u
    q.push_back( func(x, idx_q[j], F,usePDFAnzat, k));    // q for b, c, s 
  }
}


void Mset_pdf(const std::complex<double>& N, std::vector<std::complex<double>>& q, std::vector<std::complex<double>>& qbar, std::complex<double>& g, const bool usePDFAnzat)
{
  // Quarks PDG list
  std::vector<int>  idx_q = usePDFAnzat ? std::vector<int>{8} : std::vector<int>{1, 2, 4, 7, 5};// d , u , s , c , b
  std::vector<int>  idx_qb = usePDFAnzat ? std::vector<int>{8} : std::vector<int>{3, 6, 4, 7, 5};// dbar , ubar , sbar , cbar , bbar
  std::complex<double> (*func)(const std::complex<double> &,const int &, const bool);
  
  if (std::abs(N)<140.)
    func = Fbis;
  else
    func = F3;
  //PDFs
  g=func(N,0,usePDFAnzat);
  
  for (int j=0; j<Nflav; j++)
  {
    qbar.push_back(func(N, idx_qb[j], usePDFAnzat));  // qbar for b, c, s, u
    q.push_back(func(N, idx_q[j], usePDFAnzat));    // q for b, c, s
  }
}



double Luminosity(const int &i, double &xa, double &xb, const LHAPDF::PDF* F,const double &a, const bool usePDFAnzat, bool qq, int k)
{

  void (*set_pdf)( const LHAPDF::PDF* F, double &x, std::vector<double>& q, std::vector<double>& qbar, double &g, const bool usePDFAnzat, int k, bool num);
  if(use_LHAPDF)  set_pdf = set_pdf_LHAPDF;
  else set_pdf = set_pdf_fit;
  // Variables for PDFs
  double fAB = 0.;
  std::vector<double> qa_plus, qb_plus, qbara_plus, qbarb_plus;
  std::vector<double> qa_minus, qb_minus, qbara_minus, qbarb_minus;
  double ga, gb;
  double ha = 2.*Eps*xa;
  double hb = 2.*Eps*xb;

  // Adjust limits for xa_plus and xb_plus
  double xa_plus = std::min(xa*(1.+Eps), 1.0);
  double xb_plus = std::min(xb*(1.+Eps), 1.0);
  double xa_minus = xa*(1.-Eps);
  double xb_minus = xb*(1.-Eps);

  // Compute PDFs for variations
  if(i==7)  // First PDF trick
  {
    set_pdf(F, xa_plus, qa_plus, qbara_plus, ga, usePDFAnzat, k, false);   // PDFs with xa+eps
    set_pdf(F, xa_minus, qa_minus, qbara_minus, ga, usePDFAnzat, k, false); // PDFs with xa-eps
    set_pdf(F, xb_plus, qb_plus, qbarb_plus, gb, usePDFAnzat, k, false);   // PDFs with xb+eps
    set_pdf(F, xb_minus, qb_minus, qbarb_minus, gb, usePDFAnzat, k, false); // PDFs with xb-eps

    // Compute contributions to fAB
    //for (int fl = 0; fl < (usePDFAnzat ? 1 : Nflav); fl++)
    for (int fl = 0; (fl < (usePDFAnzat ? 1 : Nflav) && fl!=3) ; fl++)
    {
      double xa_der = (std::pow(xa_plus,-a)*qa_plus[fl]    - std::pow(xa_minus,-a)*qa_minus[fl])/ha;
      double xb_der = (std::pow(xb_plus,-a)*qbarb_plus[fl] - std::pow(xb_minus,-a)*qbarb_minus[fl])/hb;
      double xa_der_qbar = (std::pow(xa_plus,-a)*qbara_plus[fl] - std::pow(xa_minus,-a)*qbara_minus[fl])/ha;
      double xb_der_q = (std::pow(xb_plus,-a)*qb_plus[fl]       - std::pow(xb_minus,-a)*qb_minus[fl])/hb;
      fAB += (fl%2==0 ? gd : gu)*std::pow(xa,a)*std::pow(xb,a)*(xa_der*xb_der+xa_der_qbar*xb_der_q);

      if(usePDFAnzat) fAB = std::pow(xa,a)*std::pow(xb,a)*xa_der*xb_der; 
    }

    
    
  }
  else if (i==2 || i == 3 || i== 4 || i==5) // Second PDF trick
  {
    std::vector<double> qa, qb, qbara, qbarb;
    bool num=(i==2||i==3);
    set_pdf(F, xa, qa, qbara, ga, usePDFAnzat, k, num);
    set_pdf(F, xb, qb, qbarb, gb, usePDFAnzat, k, num);
    // Compute fAB contributions
    for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*(qq ? qa[fl]*qbarb[fl] + qb[fl]*qbara[fl] : qa[fl]*gb+qb[fl]*ga+qbara[fl]*gb+qbarb[fl]*ga) ;
    
    // Remove charm quark
    if(remove_charm_quark) fAB-= gu*(qq ? qa[3]*qbarb[3] + qb[3]*qbara[3] : qa[3]*gb+qb[3]*ga+qbara[3]*gb+qbarb[3]*ga);

    //PDF Anzat
    if(usePDFAnzat) fAB=qa[0]*qbarb[0];

    if (std::abs(fAB) > 1e10)
     info("xa = " + std::to_string(xa)  + "; xb = " + std::to_string(xb) + "; fAB = " + std::to_string(fAB));
  }
  // Output
  return fAB;
}
std::complex<double> EvolOperator(const std::complex<double> &N, std::complex<double> &EigenV)
{
  std::complex<double> res(0.0,0.0), Nb=N*exp(-Psi(1.));
  double beta0=23./6.;
  
  res=pow((1.+beta0/M_PI*AlphaS*log(MZ/Nb/muR))/(1.+beta0/M_PI*AlphaS*log(muF/muR)),EigenV/beta0);
  return res;
}

void EvolvePDF(const std::complex<double> &N, std::vector<std::complex<double>>& q, std::vector<std::complex<double>>& qbar, std::complex<double> &g)
{
  double beta0=23./6., nf=(double)Nflav, cmp=0;
  std::complex<double> V[5], NS[5], Sigma(0.,0.), tmp(0.0,0.0), Gqq, Ggg, Gqg, Ggq, rp, rm;
  std::complex<double> VE[5], NSE[5], SigmaE(0.,0.);

  Gqq=4./3.*(3./2.+1./(N*(N+1.))-2.*(Psi(N+1.)-Psi(1.)));
  Gqg=(2.+N+N*N)/(2.*N*(N+1.)*(N+2.));
  Ggq=4./3.*(2.+N+N*N)/(N*(N*N-1.));
  Ggg=beta0+2.*3.*(1./(N*(N-1.))+1./((N+1.)*(N+2.))+Psi(1.)-Psi(N+1.));
		   
  rp=0.5*(Ggg+Gqq+pow(pow(Gqq-Ggg,2.)+8.*nf*Gqg*Ggq,0.5));
  rm=0.5*(Ggg+Gqq-pow(pow(Gqq-Ggg,2.)+8.*nf*Gqg*Ggq,0.5));
  
  
  for(int i0=0; i0<Nflav; i0++){
    cmp+=1.;
    
    V[i0]=q[i0]-qbar[i0];
    Sigma+=q[i0]+qbar[i0];
    NS[i0]=Sigma-cmp*(q[i0]+qbar[i0]);
    
    VE[i0]=V[i0]*EvolOperator(N,Gqq);
    NSE[i0]=NS[i0]*EvolOperator(N,Gqq);
  }

  SigmaE=EvolOperator(N,rp)*((Gqq-rm)*Sigma+2.*nf*Gqg*g)/(rp-rm)-EvolOperator(N,rm)*((Gqq-rp)*Sigma+g*Gqg*2.*nf)/(rp-rm);

  //Evolved gluon PDF
  g=EvolOperator(N,rp)*((Ggg-rm)*g+Ggq*Sigma)/(rp-rm)-EvolOperator(N,rm)*((Ggg-rp)*g+Ggq*Sigma)/(rp-rm);

  //Evolved quark PDF
  for (int i0 = Nflav-1; i0 >= 0; i0--) {  
     q[i0] = (1.0 / nf * SigmaE - 1.0 / cmp * NSE[i0] + tmp + VE[i0]) * 0.5;
     qbar[i0] = (1.0 / nf * SigmaE - 1.0 / cmp * NSE[i0] + tmp - VE[i0]) * 0.5;
     tmp += 1.0 / cmp / (cmp - 1.0) * NSE[i0];
     cmp -= 1.0;
   }
  
}

std::complex<double> MLuminosity(const std::complex<double> &N, const bool usePDFAnzat,std::string A,bool collinear)
{  
  std::vector<std::complex<double>> q, qbar;
  std::complex<double>  g, fAB=0.;
  Mset_pdf(N, q, qbar, g, usePDFAnzat);  // Set PDFs with N
  if(collinear) EvolvePDF(N, q,qbar, g);
  // qq case
  if(A=="qq"){
  // Compute fAB contributions
  for (int fl=0; fl<Nflav; fl++) fAB+= (fl%2==0 ? gd : gu)*q[fl]*qbar[fl]*2.;

  if(usePDFAnzat) fAB=2.*q[0]*qbar[0];
  //remove charm  quark contribution 
  if(remove_charm_quark) fAB-= gu*q[3]*qbar[3]*2.;
  }
  //
  if(A=="qg"){
  for (int fl=0; fl<Nflav; fl++) fAB+= (fl%2==0 ? gd : gu)*2.*g*(q[fl]+qbar[fl]);
  //remove charm  quark contribution 
  if(remove_charm_quark) fAB-= gu*g*2.*(q[3]+qbar[3]);
  }
 

  return fAB;
}

std::complex<double> Gamma(const std::complex<double> x) {
    const int intx = (int)real(x);
    const double n = -(double)(intx < 0 ? -intx : intx);
    if (real(x) == n && imag(x) == 0.0) error("Gamma("  + std::to_string(real(x))+ "i*" + std::to_string(imag(x)) + ") undefined");

    // Works with Re(xx) > 0
    const std::complex<double> xx = (real(x) > 0.0 ? x : -x);

    // Magic numbers for Gamma function
    const double q0 = 75122.6331530;
    const double q1 = 80916.6278952;
    const double q2 = 36308.2951477;
    const double q3 = 8687.24529705;
    const double q4 = 1168.92649479;
    const double q5 = 83.8676043424;
    const double q6 = 2.50662827511;

    std::complex<double> gamma = exp((xx + .5) * log(xx + 5.5) - xx - 5.5) *
                            (q0 + q1 * xx + q2 * pow(xx, 2) + q3 * pow(xx, 3) + q4 * pow(xx, 4) +
                             q5 * pow(xx, 5) + q6 * pow(xx, 6))
                            / xx / (xx + 1.0) / (xx + 2.0) / (xx + 3.0) / (xx + 4.0) / (xx + 5.0) / (xx + 6.0);

    return (x == xx ? gamma : -M_PI / xx / sin(M_PI * xx) / gamma);
}

std::complex<double> B2(const std::complex<double> a, double b){ return Gamma(a)*Gamma(b)/Gamma(a+b);}

std::complex<double> Fbis(const std::complex<double> &N, const int &i0, const bool usePDFAnzat)
{
  return A[i0][0]*(B2(A[i0][1]+N,A[i0][2]+1.)+ (usePDFAnzat ? 0.0 : A[i0][3]*B2(A[i0][1]+N+0.5,A[i0][2]+1.)+A[i0][4]*B2(A[i0][1]+N+1.,A[i0][2]+1.)+A[i0][5]*B2(A[i0][1]+N+1.5,A[i0][2]+1.)+A[i0][6]*B2(A[i0][1]+N+2.,A[i0][2]+1.)+A[i0][7]*B2(A[i0][1]+N+2.5,A[i0][2]+1.)));
}

std::complex<double> B(const std::complex<double> &n, double a, double b)// Series expansion of B2 function
{
  return (1. - (b*(-1.+ + 2*a + b))/(2.*n) + (b*(1. + b)*(2. + 12*pow(a,2) + 12*a*(-1. + b) - 5*b + 3*pow(b,2)))/(24.*pow(n,2)) - (b*(2. + 3*b + pow(b,2))*(8*pow(a,3) + 12*pow(a,2)*(-1. + b) + pow(-1. + b,2)*b + 2*a*(2. - 5*b + 3*pow(b,2))))/(48.*pow(n,3)) + (b*(6. + 11*b + 6*pow(b,2) + pow(b,3))*(-8. + 240*pow(a,4) + 480*pow(a,3)*(-1. + b) + 18*b + 120*a*pow(-1. + b,2)*b + 5*pow(b,2) -  30*pow(b,3) + 15*pow(b,4) + 120*pow(a,2)*(2. - 5*b + 3*pow(b,2))))/(5760.*pow(n,4)) -   (b*(24. + 50*b + 35*pow(b,2) + 10*pow(b,3) + pow(b,4))* (96*pow(a,5) + 240*pow(a,4)*(-1. + b) + 120*pow(a,2)*pow(-1. + b,2)*b + 80*pow(a,3)*(2. - 5*b + 3*pow(b,2)) +      pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) + 2*a*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4))))/(11520.*pow(n,5)) +      (b*(120. + 274*b + 225*pow(b,2) + 85*pow(b,3) + 15*pow(b,4) + pow(b,5))*        (96. + 4032*pow(a,6) + 12096*pow(a,5)*(-1. + b) - 236*b + 10080*pow(a,3)*pow(-1. + b,2)*b - 84*pow(b,2) + 539*pow(b,3) -           315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6) + 5040*pow(a,4)*(2. - 5*b + 3*pow(b,2)) +           252*a*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) + 252*pow(a,2)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4))))/   (2.90304e6*pow(n,6)) - (b*(720. + 1764*b + 1624*pow(b,2) + 735*pow(b,3) + 175*pow(b,4) + 21*pow(b,5) + pow(b,6))*
        (1152*pow(a,7) + 4032*pow(a,6)*(-1. + b) + 5040*pow(a,4)*pow(-1. + b,2)*b + 2016*pow(a,5)*(2. - 5*b + 3*pow(b,2)) +    252*pow(a,2)*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) +           pow(-1. + b,2)*b*(80. - 34*b - 57*pow(b,2) + 18*pow(b,3) + 9*pow(b,4)) +       168*pow(a,3)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4)) +        2*a*(96. - 236*b - 84*pow(b,2) + 539*pow(b,3) - 315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6))))/(5.80608e6*pow(n,7)) +    (b*(5040. + 13068*b + 13132*pow(b,2) + 6769*pow(b,3) + 1960*pow(b,4) + 322*pow(b,5) + 28*pow(b,6) + pow(b,7))*    (-1152. + 34560*pow(a,8) + 138240*pow(a,7)*(-1. + b) + 3088*b + 241920*pow(a,5)*pow(-1. + b,2)*b + 884*pow(b,2) -      8140*pow(b,3) + 6055*pow(b,4) + 840*pow(b,5) - 1890*pow(b,6) + 180*pow(b,7) + 135*pow(b,8) +       80640*pow(a,6)*(2. - 5*b + 3*pow(b,2)) + 20160*pow(a,3)*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) +        240*a*pow(-1. + b,2)*b*(80. - 34*b - 57*pow(b,2) + 18*pow(b,3) + 9*pow(b,4)) +         10080*pow(a,4)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4)) +      240*pow(a,2)*(96. - 236*b - 84*pow(b,2) + 539*pow(b,3) - 315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6))))/   (1.3934592e9*pow(n,8)) - (b*(40320. + 109584*b + 118124*pow(b,2) + 67284*pow(b,3) + 22449*pow(b,4) + 4536*pow(b,5) +        546*pow(b,6) + 36*pow(b,7) + pow(b,8))*(7680*pow(a,9) + 34560*pow(a,8)*(-1. + b) + 80640*pow(a,6)*pow(-1. + b,2)*b +        23040*pow(a,7)*(2. - 5*b + 3*pow(b,2)) + 10080*pow(a,4)*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) +         240*pow(a,2)*pow(-1. + b,2)*b*(80. - 34*b - 57*pow(b,2) + 18*pow(b,3) + 9*pow(b,4)) +         4032*pow(a,5)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4)) +       pow(-1. + b,2)*b*(-1008. + 668*b + 768*pow(b,2) - 527*pow(b,3) - 135*pow(b,4) + 75*pow(b,5) + 15*pow(b,6)) +         160*pow(a,3)*(96. - 236*b - 84*pow(b,2) + 539*pow(b,3) - 315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6)) +         2*a*(-1152. + 3088*b + 884*pow(b,2) - 8140*pow(b,3) + 6055*pow(b,4) + 840*pow(b,5) - 1890*pow(b,6) + 180*pow(b,7) +            135*pow(b,8))))/(2.7869184e9*pow(n,9)) + (b*(362880. + 1026576*b + 1172700*pow(b,2) + 723680*pow(b,3) + 269325*pow(b,4) + 63273*pow(b,5) + 9450*pow(b,6) +         870*pow(b,7) + 45*pow(b,8) + pow(b,9))*(7680. + 101376*pow(a,10) + 506880*pow(a,9)*(-1. + b) - 22112*b +      1520640*pow(a,7)*pow(-1. + b,2)*b - 3960*pow(b,2) + 62524*pow(b,3) - 56958*pow(b,4) - 1265*pow(b,5) + 20559*pow(b,6) -       5082*pow(b,7) - 1980*pow(b,8) + 495*pow(b,9) + 99*pow(b,10) + 380160*pow(a,8)*(2. - 5*b + 3*pow(b,2)) +         266112*pow(a,5)*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) +          10560*pow(a,3)*pow(-1. + b,2)*b*(80. - 34*b - 57*pow(b,2) + 18*pow(b,3) + 9*pow(b,4)) +        88704*pow(a,6)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4)) +           132*a*pow(-1. + b,2)*b*(-1008. + 668*b + 768*pow(b,2) - 527*pow(b,3) - 135*pow(b,4) + 75*pow(b,5) + 15*pow(b,6)) +           5280*pow(a,4)*(96. - 236*b - 84*pow(b,2) + 539*pow(b,3) - 315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6)) +   132*pow(a,2)*(-1152. + 3088*b + 884*pow(b,2) - 8140*pow(b,3) + 6055*pow(b,4) + 840*pow(b,5) - 1890*pow(b,6) + 180*pow(b,7) + 135*pow(b,8))))/(3.678732288e11*pow(n,10)))/pow(n,b);
}


std::complex<double> F3(const std::complex<double> &N, const int &i0, const bool usePDFAnzat) // Calcul simplifié dans le cas |N| > 1e2
{
  return A[i0][0]*Gamma(A[i0][2]+1.)*(B(N,A[i0][1],1.+A[i0][2])+(usePDFAnzat ? 0.0 : B(N,A[i0][1]+0.5,1.+A[i0][2])*A[i0][3]+B(N,A[i0][1]+1.,1.+A[i0][2])*A[i0][4]+B(N,A[i0][1]+1.5,1.+A[i0][2])*A[i0][5]+B(N,A[i0][1]+2.,1.+A[i0][2])*A[i0][6]+B(N,A[i0][1]+2.5,1.+A[i0][2])*A[i0][7]));
}


///**********************************************************///
///                   NLO CALCULATION                        ///  
///**********************************************************///

void Kinematics(double &xa, double &xb, double &jac, double &s, double *x, Process *p)
{
  // Estimate tau
  double tau = pow(MZ,2.)/p->sh;

  // Computes xa + integration from 0->1 becomes from tau to 1
  xa=tau/pow(tau,x[0]); //Only 1 integration here

  // Computes xb
  xb=tau/xa;
  
  // Computes shat
  s = xa*xb*p->sh;
  
  // Jacobian + delta function integration
  jac = -xa*log(tau)*xb; //-tau*log(tau); ?
}

// NLO Kinematics with 2 variables: xa and xb
void Kinematics2(double &xa, double &xb, double &jac, double &s, double *x, Process *p)
{
  // Estimate tau
  double tau = pow(MZ,2.)/p->sh;

  // Computes xa + integration from 0->1 becomes from tau to 1
  xa=tau/pow(tau,x[0]);

  // Computes xb + integration from 0->1 becomes from tau/xa to 1 , ensures that xa*xb >= tau
  xb=tau/pow(tau/xa,x[1])/xa;

  // Computes shat
  s = xa*xb*p->sh;

  // Jacobian of xa, xb -> x0, x1
  jac = xa*xb*pow(log(tau),2.)*x[0];
}

// NLO Kinematics with 2 variables: y and z
void Kinematics2Plus(double &y, double &z, double &jac, double &s, double *x, Process *p)
{
  // Estimate tau
  double tau = pow(MZ,2)/p->sh;

  // Computes y + integration from log(tau)/2 -> -log(tau)/2 becomes from tau/xa to 1 , ensures that xa*xb >= tau
  y=log(tau)*(0.5-x[0]);

  // Computes z + integration from 0->1 becomes from tau to 1
  z=x[1]+tau*exp(2*abs(y))*(1-x[1]);
  //z=pow(tau,x[0]*x[1]);

  // Computes shat
  s = tau/z*p->sh;

  // Jacobian of dy dz/z -> dx0 dx1 !! Carrefull a 1/z term was ommited as it has to be absorbed in the fAB factor to be integrated
  jac = -log(tau)*(1-tau*exp(2*abs(y)));
}

void SetLHAPDF(const LHAPDF::PDF* F, double &x, double q[5], double qb[5], double &g)
{

  for(int i0=-5; i0<6; i0++){
    if(i0<0) qb[-1-i0]=F->xfxQ2(i0,x,muF*muF)/x;
    else if(i0==0) g=F->xfxQ2(i0,x,muF*muF)/x;
    else q[i0-1]=F->xfxQ2(i0,x,muF*muF)/x;
  }
  
}
void SetCouplings(const int &Subproc, Process *p, double &fAB, double &xa, double &xb)
//void SetCouplings(const int &Subproc, Process *p, const PDF *F double &fAB, double &xa, double &xb)
{ 
  void (*set_pdf)( const LHAPDF::PDF* F, double &x, std::vector<double>& q, std::vector<double>& qbar, double &g, const bool usePDFAnzat, int k, bool num);
  if(use_LHAPDF)  set_pdf = set_pdf_LHAPDF;
  else set_pdf = set_pdf_fit;
  fAB=0.;
	double  ga, gb;
  std::vector<double> qa, qb, qbara, qbarb;
	double tau = pow(MZ,2.)/p->sh;

	if(Subproc==0) // For Born and regular virtual corrections (not +-distribution)
	{
	  set_pdf(p->F,xa,qa,qbara,ga, usePDFAnzat, 0, true); // PDFs with xa
	  set_pdf(p->F,xb,qb,qbarb,gb, usePDFAnzat, 0, true); // PDFs with xb
    for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*(qa[fl]*qbarb[fl]+qb[fl]*qbara[fl]);
    if(remove_charm_quark){fAB -= gu*(qa[3]*qbarb[3]+qb[3]*qbara[3]);}
    if(usePDFAnzat) fAB=qa[0]*qbarb[0]+qb[0]*qbara[0];
	}

	else if(Subproc==1) // For +-distribution part of virtual corrections. !!Carefull xa and xb variables are replaced by z and y variables!!
	{
		double z=xb; // Renaming variable
		double y=xa; // Renaming variable
		double xapdf=sqrt(tau/z)*exp(y); // xa value in the PDFs
		double xbpdf=sqrt(tau/z)*exp(-y); // xb value in the PDFs
		double xaplus=sqrt(tau)*exp(y); // xa for z=1
		double xbplus=sqrt(tau)*exp(-y); // xb for z=1
   
		double  gaplus, gbplus;
    std::vector<double>  qaplus, qbplus, qbaraplus, qbarbplus;

		set_pdf(p->F,xapdf,qa,qbara,ga,usePDFAnzat, 0, true); // xa(y,z) and xb(y,z)
		set_pdf(p->F,xbpdf,qb,qbarb,gb,usePDFAnzat, 0, true);

		set_pdf(p->F,xaplus,qaplus,qbaraplus,gaplus,usePDFAnzat, 0, true); // Evaluation of PDFs in z=1
    set_pdf(p->F,xbplus,qbplus,qbarbplus,gbplus,usePDFAnzat, 0, true);

    for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*(2.*(log(MZ/muF)+log(1-z))*((1+z*z)*qa[fl]*qbarb[fl]/z-2.*qaplus[fl]*qbarbplus[fl])/(1-z));
    for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*(2.*(log(MZ/muF)+log(1-z))*((1+z*z)*qb[fl]*qbara[fl]/z-2.*qbaraplus[fl]*qbplus[fl])/(1-z)); //dbar quark
    

    //for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*(2.*(log(MZ/muF))*((1+z*z)*qa[fl]*qbarb[fl]/z-2.*qaplus[fl]*qbarbplus[fl])/(1-z));
    //for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*(2.*(log(MZ/muF))*((1+z*z)*qb[fl]*qbara[fl]/z-2.*qbaraplus[fl]*qbplus[fl])/(1-z)); //dbar quark

    //for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*2.*(log(1-z))*((1+z*z)*qa[fl]*qbarb[fl]/z-2.*qaplus[fl]*qbarbplus[fl])/(1-z);
    //for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*2.*(log(1-z))*((1+z*z)*qb[fl]*qbara[fl]/z-2.*qbaraplus[fl]*qbplus[fl])/(1-z); 
    if(remove_charm_quark){
      fAB -= gu*(2.*(log(MZ/muF)+log(1-z))*((1+z*z)*qa[3]*qbarb[3]/z-2.*qaplus[3]*qbarbplus[3])/(1-z));
      fAB -= gu*(2.*(log(MZ/muF)+log(1-z))*((1+z*z)*qb[3]*qbara[3]/z-2.*qbaraplus[3]*qbplus[3])/(1-z));
    }
    if(usePDFAnzat) fAB=(2.*(log(MZ/muF)+log(1-z))*((1+z*z)*qa[0]*qbarb[0]/z-2.*qaplus[0]*qbarbplus[0])/(1-z));
	}

	else if(Subproc==2) // For +-distribution 1D part of virtual corrections. !!Carefull xa and xb variables are replaced by z and y variables!!
  {
    double z0=log(1-tau*exp(2.*abs(xa))); // Leftover factor to take into account the mismatch between the distribution and the integral range
    double y=xa; // Renaming variable
    double xaplus=sqrt(tau)*exp(y); // xa for z=1
    double xbplus=sqrt(tau)*exp(-y); // xb for z=1

    set_pdf(p->F,xaplus,qa,qbara,ga,usePDFAnzat, 0, true); // Evaluation of PDFs in z=1
    set_pdf(p->F,xbplus,qb,qbarb,gb,usePDFAnzat, 0, true);
		
    for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*z0*(4.*log(MZ/muF)+2.*z0)*(qa[fl]*qbarb[fl]+qb[fl]*qbara[fl]);

    //for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*z0*(4.*log(MZ/muF))*(qa[fl]*qbarb[fl]+qb[fl]*qbara[fl]);
    //for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*z0*(2.*z0)*(qa[fl]*qbarb[fl]+qb[fl]*qbara[fl]);

    if(remove_charm_quark) fAB -= gu*z0*(4*log(MZ/muF)+z0)*(qa[3]*qbarb[3]+qb[3]*qbara[3]);
    if(usePDFAnzat) fAB=z0*(4.*log(MZ/muF)+2.*z0)*(qa[0]*qbarb[0]);
  }


	else if(Subproc==3) // For real corrections
	{
	  set_pdf(p->F,xa,qa,qbara,ga,usePDFAnzat, 0, true); // PDFs with xa
	  set_pdf(p->F,xb,qb,qbarb,gb,usePDFAnzat, 0, true); // PDFs with xb
   
    for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*(qa[fl]*gb+qb[fl]*ga+qbara[fl]*gb+qbarb[fl]*ga);
    if(remove_charm_quark) fAB -= gu*(qa[3]*gb+qb[3]*ga+qbara[3]*gb+qbarb[3]*ga);
    if(usePDFAnzat) fAB=qa[0]*gb+qb[0]*ga;

	}

	else std::cout << "Label Error, please specify Born or Virtual (0), Virtual+ (1/2), Real (3)" << std::endl;
}