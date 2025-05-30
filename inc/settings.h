#ifndef SETTINGS_H
#define SETTINGS_H
// ************************************************************************* //
// Declaration of setting routines                                           //
//                                                                           //
// By Benjamin Fuks - 02.12.2024.                                            //
// ************************************************************************* //

// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //

// ************************************************************************* //
// Definition of the message services class                                  //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <complex>

void PDFInit();
void InitialiseCoefficients(double A[9][8]);
void OpenFiles(std::vector<std::ofstream>&, const std::vector<std::string>&);
void CloseFiles(std::vector<std::ofstream>&);
std::complex<double> Psi(std::complex<double> );
std::complex<double> HarmonicNumber(std::complex<double> );
int coefBinomial(int n, int k);
double derive_x_k(double (* f)(double x, int i0, const LHAPDF::PDF* F), double x, int i0, const LHAPDF::PDF* F, double eps,int k);
int indice(double x, double xmin, double xmax, int np);
std::complex<double> pgamma(int m, std::complex<double> z);

#endif

