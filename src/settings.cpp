// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <fstream>        // File streams               //
// -------- Classes ----------------------------------- //
#include "messages.h"     // Message services           //
#include "LHAPDF/LHAPDF.h"// LHAPDF

void PDFInit()
{
  // Init
  const std::string pathstr = "/usr/local/share/LHAPDF";
  LHAPDF::setPaths(pathstr);

  info("PDF Paths:");
  for (const std::string& p : LHAPDF::paths()) info("  " + p);
  info("Available PDFs:");
  for (const std::string& s : LHAPDF::availablePDFSets()) info("  " + s);
}



// Function to initialise the PDF coefficients
void InitialiseCoefficients(double A[8][8])
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

  const double coeffs[8][8] = {
    {3.36311, -1.35686, 12.01851, -9.06557, 59.82888, -209.61905, 358.78810, -214.30364},
    {0.17088, -1.32599, 4.86054, -2.74299, 27.12173, -51.56600, 58.71166, -29.80953},
    {0.15891, -1.33284, 4.04498, -1.28483, 25.16123, -21.03617, 37.34074, -32.12262},
    {0.18060, -1.32147, 4.80454, -4.01851, 19.31558, -55.90824, 71.69672, -33.26465},
    {0.19715, -1.31408, 12.42845, -7.36475, 44.88323, -146.11711, 243.91655, -145.10893},
    {0.10914, -1.34332, 16.06492, -9.85651, 68.81877, -260.94287, 489.90364, -321.74473},
    {0.19950, -1.31377, 10.66072, -5.55560, 37.16200, -128.66354, 217.17546, -124.50014},
    {0.14115, -1.33602, -0.27091, -6.04545, 15.21004, -19.76661, 13.13392, -3.53289}
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

