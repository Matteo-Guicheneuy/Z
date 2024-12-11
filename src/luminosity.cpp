// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <string>     // Strings                        //
// -------- Classes ----------------------------------- //
#include "messages.h"     // Message services           //
#include "luminosity.h"   // Luminosities and PDFs      //
#include "constants.h"    // Constants                  //
#include "coefficients.h"
#include <iostream>

double F_real(const double &x, int i0, const bool usePDFAnzat)
{
  // PDF in real space from MSTW-like parameterisatio (eq 2.70 of Yehudi's thesis)
  return A[i0][0] * pow(x,A[i0][1]) * pow(1.-x,A[i0][2]) *
    (usePDFAnzat ? 1. : (1.+A[i0][3]*sqrt(x)+A[i0][4]*x+A[i0][5]*pow(x,1.5)+A[i0][6]*pow(x,2)+A[i0][7]*pow(x, 2.5)));

}

double F_real_bis(const double &x, int i0, const bool usePDFAnzat)
{
  // Derivative of PDF in real space
  return (usePDFAnzat ?  A[i0][0]*pow(x,A[i0][1]+1.)*pow(1.-x,A[i0][2])*(A[i0][1]+1.-x*A[i0][2]/(1.-x)) :
  A[i0][0]*pow(x,A[i0][1]+1.)*pow(1.-x,A[i0][2])*((A[i0][1]+1.-x*A[i0][2]/(1.-x))*(1.+A[i0][3]*sqrt(x)+A[i0][4]*x+A[i0][5]*pow(x,1.5)+A[i0][6]*pow(x,2.)+A[i0][7]*pow(x,2.5))+0.5*A[i0][3]*sqrt(x)+A[i0][4]*x+1.5*A[i0][5]*pow(x,1.5)+2.*A[i0][6]*pow(x,2.)+2.5*A[i0][7]*pow(x,2.5)));
}

void set_pdf(int i, const double &x, double q[5], double qbar[5], double &g, const bool usePDFAnzat)
{
  const int idx_q[5] = {1, 2, 4, 7, 5};  // d , u , s , c , b
  const int idx_qb[5] = {3, 6, 4, 7, 5}; // dbar , ubar , sbar , cbar , bbar
  //auto func = (i == 0) ? F_real : F_real_bis;
  double (*func)(const double&, int, const bool);
  if (i==5)
    func = F_real;
  else
    func = F_real_bis;

  for (int j=0; j<5; j++)
  {
    qbar[j] = func(x, idx_qb[j], usePDFAnzat);  // qbar for b, c, s, u
    q[j] = func(x, idx_q[j], usePDFAnzat);    // q for b, c, s
    g = func(x, 0, usePDFAnzat);              // g
  }
}
void set_pdfN(const std::complex<double> &N, std::complex<double> q[5], std::complex<double> qbar[5], std::complex<double> &g, const bool usePDFAnzat)
{
  const int idx_q[5] = {1, 2, 4, 7, 5};
  const int idx_qb[5] = {3, 6, 4, 7, 5};
  std::complex<double> (*func)(const std::complex<double> &,const int &, const bool);
  
  if (std::abs(N)<140.)
    func = Fbis;
  else
    func = F3;
  //PDFs
  g=func(N,0,usePDFAnzat);
  for (int j=0; j<Nflav; j++)
  {
    qbar[j] = func(N, idx_qb[j], usePDFAnzat);  // qbar for b, c, s, u
    q[j] = func(N, idx_q[j], usePDFAnzat);    // q for b, c, s
  }
}



double Luminosity(const int &i, const double &xa, const double &xb, const double &a, const bool usePDFAnzat)
{
  // Variables for PDFs
  double fAB = 0.;
  double qa_plus[Nflav], qb_plus[Nflav], qbara_plus[Nflav], qbarb_plus[Nflav];
  double qa_minus[Nflav], qb_minus[Nflav], qbara_minus[Nflav], qbarb_minus[Nflav];
  double g;
  double ha = 2.*eps*xa;
  double hb = 2.*eps*xb;

  // Adjust limits for xa_plus and xb_plus
  double xa_plus = std::min(xa*(1.+eps), 1.0);
  double xb_plus = std::min(xb*(1.+eps), 1.0);
  double xa_minus = xa*(1.-eps);
  double xb_minus = xb*(1.-eps);

  // Compute PDFs for variations
  if(i==5)  // First PDF trick
  {
    set_pdf(i, xa_plus, qa_plus, qbara_plus, g, usePDFAnzat);   // PDFs with xa+eps
    set_pdf(i, xa_minus, qa_minus, qbara_minus, g, usePDFAnzat); // PDFs with xa-eps
    set_pdf(i, xb_plus, qb_plus, qbarb_plus, g, usePDFAnzat);   // PDFs with xb+eps
    set_pdf(i, xb_minus, qb_minus, qbarb_minus, g, usePDFAnzat); // PDFs with xb-eps

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
  else if (i==1 || i == 4) // Second PDF trick
  {
    double qa[Nflav], qb[Nflav], qbara[Nflav], qbarb[Nflav];
    set_pdf(i, xa, qa, qbara, g, usePDFAnzat);
    set_pdf(i, xb, qb, qbarb, g, usePDFAnzat);
    // Compute fAB contributions
    for (int fl = 0; fl < Nflav; ++fl) fAB += (fl%2==0 ? gd : gu)*(qa[fl]*qbarb[fl] + qb[fl]*qbara[fl]) ;
    
    // Remove charm quark
    fAB-= gu*(qa[3]*qbarb[3] + qb[3]*qbara[3]);

    //PDF Anzat
    if(usePDFAnzat) fAB=qa[0]*qbarb[0];

    if (std::abs(fAB) > 1e4)
     info("xa = " + std::to_string(xa)  + "; xb = " + std::to_string(xb) + "; fAB = " + std::to_string(fAB));
  }

  // Output
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

std::complex<double> B(const std::complex<double> &n, double a, double b)
{
  return (1. - (b*(-1.+ + 2*a + b))/(2.*n) + (b*(1. + b)*(2. + 12*pow(a,2) + 12*a*(-1. + b) - 5*b + 3*pow(b,2)))/(24.*pow(n,2)) - (b*(2. + 3*b + pow(b,2))*(8*pow(a,3) + 12*pow(a,2)*(-1. + b) + pow(-1. + b,2)*b + 2*a*(2. - 5*b + 3*pow(b,2))))/(48.*pow(n,3)) + (b*(6. + 11*b + 6*pow(b,2) + pow(b,3))*(-8. + 240*pow(a,4) + 480*pow(a,3)*(-1. + b) + 18*b + 120*a*pow(-1. + b,2)*b + 5*pow(b,2) -  30*pow(b,3) + 15*pow(b,4) + 120*pow(a,2)*(2. - 5*b + 3*pow(b,2))))/(5760.*pow(n,4)) -   (b*(24. + 50*b + 35*pow(b,2) + 10*pow(b,3) + pow(b,4))* (96*pow(a,5) + 240*pow(a,4)*(-1. + b) + 120*pow(a,2)*pow(-1. + b,2)*b + 80*pow(a,3)*(2. - 5*b + 3*pow(b,2)) +      pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) + 2*a*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4))))/(11520.*pow(n,5)) +      (b*(120. + 274*b + 225*pow(b,2) + 85*pow(b,3) + 15*pow(b,4) + pow(b,5))*        (96. + 4032*pow(a,6) + 12096*pow(a,5)*(-1. + b) - 236*b + 10080*pow(a,3)*pow(-1. + b,2)*b - 84*pow(b,2) + 539*pow(b,3) -           315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6) + 5040*pow(a,4)*(2. - 5*b + 3*pow(b,2)) +           252*a*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) + 252*pow(a,2)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4))))/   (2.90304e6*pow(n,6)) - (b*(720. + 1764*b + 1624*pow(b,2) + 735*pow(b,3) + 175*pow(b,4) + 21*pow(b,5) + pow(b,6))*
        (1152*pow(a,7) + 4032*pow(a,6)*(-1. + b) + 5040*pow(a,4)*pow(-1. + b,2)*b + 2016*pow(a,5)*(2. - 5*b + 3*pow(b,2)) +    252*pow(a,2)*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) +           pow(-1. + b,2)*b*(80. - 34*b - 57*pow(b,2) + 18*pow(b,3) + 9*pow(b,4)) +       168*pow(a,3)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4)) +        2*a*(96. - 236*b - 84*pow(b,2) + 539*pow(b,3) - 315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6))))/(5.80608e6*pow(n,7)) +    (b*(5040. + 13068*b + 13132*pow(b,2) + 6769*pow(b,3) + 1960*pow(b,4) + 322*pow(b,5) + 28*pow(b,6) + pow(b,7))*    (-1152. + 34560*pow(a,8) + 138240*pow(a,7)*(-1. + b) + 3088*b + 241920*pow(a,5)*pow(-1. + b,2)*b + 884*pow(b,2) -      8140*pow(b,3) + 6055*pow(b,4) + 840*pow(b,5) - 1890*pow(b,6) + 180*pow(b,7) + 135*pow(b,8) +       80640*pow(a,6)*(2. - 5*b + 3*pow(b,2)) + 20160*pow(a,3)*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) +        240*a*pow(-1. + b,2)*b*(80. - 34*b - 57*pow(b,2) + 18*pow(b,3) + 9*pow(b,4)) +         10080*pow(a,4)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4)) +      240*pow(a,2)*(96. - 236*b - 84*pow(b,2) + 539*pow(b,3) - 315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6))))/   (1.3934592e9*pow(n,8)) - (b*(40320. + 109584*b + 118124*pow(b,2) + 67284*pow(b,3) + 22449*pow(b,4) + 4536*pow(b,5) +        546*pow(b,6) + 36*pow(b,7) + pow(b,8))*(7680*pow(a,9) + 34560*pow(a,8)*(-1. + b) + 80640*pow(a,6)*pow(-1. + b,2)*b +        23040*pow(a,7)*(2. - 5*b + 3*pow(b,2)) + 10080*pow(a,4)*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) +         240*pow(a,2)*pow(-1. + b,2)*b*(80. - 34*b - 57*pow(b,2) + 18*pow(b,3) + 9*pow(b,4)) +         4032*pow(a,5)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4)) +       pow(-1. + b,2)*b*(-1008. + 668*b + 768*pow(b,2) - 527*pow(b,3) - 135*pow(b,4) + 75*pow(b,5) + 15*pow(b,6)) +         160*pow(a,3)*(96. - 236*b - 84*pow(b,2) + 539*pow(b,3) - 315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6)) +         2*a*(-1152. + 3088*b + 884*pow(b,2) - 8140*pow(b,3) + 6055*pow(b,4) + 840*pow(b,5) - 1890*pow(b,6) + 180*pow(b,7) +            135*pow(b,8))))/(2.7869184e9*pow(n,9)) + (b*(362880. + 1026576*b + 1172700*pow(b,2) + 723680*pow(b,3) + 269325*pow(b,4) + 63273*pow(b,5) + 9450*pow(b,6) +         870*pow(b,7) + 45*pow(b,8) + pow(b,9))*(7680. + 101376*pow(a,10) + 506880*pow(a,9)*(-1. + b) - 22112*b +      1520640*pow(a,7)*pow(-1. + b,2)*b - 3960*pow(b,2) + 62524*pow(b,3) - 56958*pow(b,4) - 1265*pow(b,5) + 20559*pow(b,6) -       5082*pow(b,7) - 1980*pow(b,8) + 495*pow(b,9) + 99*pow(b,10) + 380160*pow(a,8)*(2. - 5*b + 3*pow(b,2)) +         266112*pow(a,5)*pow(-1. + b,2)*b*(-6. + b + 3*pow(b,2)) +          10560*pow(a,3)*pow(-1. + b,2)*b*(80. - 34*b - 57*pow(b,2) + 18*pow(b,3) + 9*pow(b,4)) +        88704*pow(a,6)*(-8. + 18*b + 5*pow(b,2) - 30*pow(b,3) + 15*pow(b,4)) +           132*a*pow(-1. + b,2)*b*(-1008. + 668*b + 768*pow(b,2) - 527*pow(b,3) - 135*pow(b,4) + 75*pow(b,5) + 15*pow(b,6)) +           5280*pow(a,4)*(96. - 236*b - 84*pow(b,2) + 539*pow(b,3) - 315*pow(b,4) - 63*pow(b,5) + 63*pow(b,6)) +   132*pow(a,2)*(-1152. + 3088*b + 884*pow(b,2) - 8140*pow(b,3) + 6055*pow(b,4) + 840*pow(b,5) - 1890*pow(b,6) + 180*pow(b,7) + 135*pow(b,8))))/(3.678732288e11*pow(n,10)))/pow(n,b);
}


std::complex<double> F3(const std::complex<double> &N, const int &i0, const bool usePDFAnzat) // Calcul simplifié dans le cas |N| > 1e2
{
  return A[i0][0]*Gamma(A[i0][2]+1.)*(B(N,A[i0][1],1.+A[i0][2])+(usePDFAnzat ? 0.0 : B(N,A[i0][1]+0.5,1.+A[i0][2])*A[i0][3]+B(N,A[i0][1]+1.,1.+A[i0][2])*A[i0][4]+B(N,A[i0][1]+1.5,1.+A[i0][2])*A[i0][5]+B(N,A[i0][1]+2.,1.+A[i0][2])*A[i0][6]+B(N,A[i0][1]+2.5,1.+A[i0][2])*A[i0][7]));
}



