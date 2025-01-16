// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <chrono>     // Timing service                 //
#include <fstream>    // File streams                   //
#include <iostream>   // In/Out streams                 //

#include <complex>
// -------- Classes ----------------------------------- //
#include "main.h"
#include "messages.h" // Message services               //
#include "process.h"  // Process setting                //
#include "settings.h" // Setting routines               //
#include "coefficients.h"

#include "luminosity.h"
// -------- Functions --------------------------------- //
void PerformIntegration(double&, double&, double&, Process*, size_t);
double A[8][8];

// ************************************************************************* //
//  Main code                                                                //
// ************************************************************************* //
int main(int argc, char* argv[])
{
  // Checking the arguments
  if(argc!=2)
  {
    std::cout << "Usage: program <PDF ID (1=CT18 NLO, 2=CT14LN, 3=CT14LL, 4=CTEQ6M)>" << std::endl;
    return 0;
  }

  // Initializing message services
  InitMessages();

  // PDF
  PDFInit();

  // Setup of the collision setup, pdfs, scales, etc...
  const int nsteps = 10;
  const double s0 = 300.0, sf = 1e6;
  //const double s0 = 200.0, sf = 600.0; < 1.6 10 steps
  InitialiseCoefficients(A) ; // For CT18NLO

  // Output file handling
  std::vector<std::string> filenames =
  {
    "LOXsec.txt", "NLOXsec.txt", "ResumXsec_C.txt", "ResumXsec_dyn.txt",
    "ResumXsec_dyn_1.txt", "ResumXsec_NC.txt", "ExpandedXsec.txt",
    "ResumXsecHS.txt", "ExpandedXsecHS.txt", "ResumXsecHSConvert.txt",
    "ExpandedXsecHSConvert.txt", "test.txt", "int_res.txt"
  };
  std::vector<std::ofstream> files(filenames.size());
  OpenFiles(files, filenames);

// Use Pdf Anzat or not: //
  //const bool usePDFAnzat=true;

  // Main computation loop
  for (int ic=7; ic<=10; ++ic)
  {
    double sc = s0 * pow(sf / s0, 1.0 * ic / nsteps);
    double res1 = 0, err1 = 0, res2 = 0, err2 = 0, res3 = 0, err3 = 0;
    for (int kID=2; kID<=5; ++kID)
    {
      // Initialisation
      auto start = std::chrono::high_resolution_clock::now();
      Process* proc = new Process(argv[1], sc, kID, ic, usePDFAnzat);
      info("Processing kID = " + std::to_string(kID) + ", ic = " + std::to_string(ic));

      // Results
      double res = 0, err = 0,  chi = 0;
      if (kID==1 || kID == 4)PerformIntegration(res2, err2, chi, proc, 3); 
      else if (kID == 5) PerformIntegration(res3, err3, chi, proc, 3);
      else if (kID == 2) PerformIntegration(res1, err1, chi, proc, 1);
      
      // Output file handling
      files[11] << "kID=" << kID << " sigma=" << res << "+/-" << err << " sc=" << sc << std::endl;
      switch (kID)
      {
        case 0: files[0] << res << "," << sc << std::endl; break;
        case 1: files[1] << res << "," << sc << std::endl; break;
        case 2: files[5] << res1 << "+/-" << err1 << "," << sc << std::endl; break;
        case 3: files[6] << res << "," << sc << std::endl; break;
        case 4: files[3] << res2 << "+/-" << err2 << "," << sc << std::endl; break;
        case 5: files[4] << res2 << "+/-" << err2 << ","  << std::abs(res1-res2)/std::sqrt(err1*err1+err2*err2)<< sc << std::endl; break;
        case 6: files[4] << res2 << "+/-" << err2 << "," << sc << std::endl; break;

        default: error("Invalid Process ID: kID=" + std::to_string(kID));
      }
      
      // Timer
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      info("kID = " + std::to_string(kID) + ", duration=" + std::to_string(duration.count() / 1e6) + "s, sc = " +  std::to_string(sc));
      delete proc;
    }

          std::cout << "////" << "incertitude = " << std::abs(res1-res2)/std::sqrt(err1*err1+err2*err2)  << "  "<< std::abs(res1-res3)/std::sqrt(err1*err1+err3*err3) << " res1 = " << res1 << " res2 = " << res2 << " res3 = " << res3 <<" err1 = " << err1 << " err2 = " << err2 << " err3 = " << err3 << "////" <<   std::endl;
    
  }

  // End of the program
  CloseFiles(files);
  return 0;
}
