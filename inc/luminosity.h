#ifndef LUMINOSITY_H
#define LUMINOSITY_H
// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <complex>                                      //


// -------- Functions --------------------------------- //
double F_real(const double &, int,const bool);
double F_real_bis(const double &, int,const bool);
void set_pdf(int, const double &, double q[5], double qbar[5],const bool usePDFAnzat);
void set_pdfN(const std::complex<double> &, std::complex<double> q[5], std::complex<double> qbar[5], std::complex<double> &g,const bool usePDFAnzat);
double Luminosity(const int &, const double &, const double &, const double &,const bool);
std::complex<double> Gamma(const std::complex<double>);
std::complex<double> B2(const std::complex<double>, double);
std::complex<double> Fbis(const std::complex<double> &, const int &,const bool);
std::complex<double> B(const std::complex<double> &, double, double);
std::complex<double> F3(const std::complex<double> &, const int &,const bool);
#endif
