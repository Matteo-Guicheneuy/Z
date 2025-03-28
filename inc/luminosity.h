#ifndef LUMINOSITY_H
#define LUMINOSITY_H
// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <complex>                                      //
#include <vector>   
#include "process.h"                                    //

// -------- Functions --------------------------------- //
double F_real(const double &, int,const bool, int);
double F_real_bis(const double &, int,const bool, int);
void set_pdf(int i, const double &x, std::vector<double>& q, std::vector<double>& qbar, double &g, const bool usePDFAnzat, int k, bool num);
void set_pdfN(const std::complex<double>& N, std::vector<std::complex<double>>& q, std::vector<std::complex<double>>& qbar, std::complex<double>& g, const bool usePDFAnzat);
double Luminosity(const int &, double &, double &, const LHAPDF::PDF* F,const double &,const bool, bool, int);
std::complex<double> MLuminosity(const std::complex<double> &N, const bool usePDFAnzat,std::string A,bool collinear=false);
std::complex<double> Gamma(const std::complex<double>);
std::complex<double> B2(const std::complex<double>, double);
std::complex<double> Fbis(const std::complex<double> &, const int &,const bool);
std::complex<double> B(const std::complex<double> &, double, double);
std::complex<double> F3(const std::complex<double> &, const int &,const bool);
void Kinematics(double &xa, double &xb, double &jac, double &s, double *x, Process *p);
void Kinematics2(double &xa, double &xb, double &jac, double &s, double *x, Process *p);
void Kinematics2Plus(double &y, double &z, double &jac, double &s, double *x, Process *p);
void SetCouplings(const int &Subproc, Process *p, double &fAB, double &xa, double &xb);
void SetLHAPDF(const LHAPDF::PDF* F, double &x, double q[5], double qb[5], double &g);
void Mset_pdf(const std::complex<double>& N, std::vector<std::complex<double>>& q, std::vector<std::complex<double>>& qbar, std::complex<double>& g, const bool usePDFAnzat);


#endif
