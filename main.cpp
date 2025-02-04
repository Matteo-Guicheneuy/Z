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
void PerformIntegration(double&, double&, double&, Process*, size_t, double);
double A[8][8];

// ************************************************************************* //
//  Main code                                                                //
// ************************************************************************* //
// To Run do: ./Run 1
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
  // To use PDF Anzat  see main.h

  // Setup square root of the center of mass energy
  const int nsteps = 10;
  const double s0 = 300.0, sf = 1e6;
  //Setup PDF fits
  InitialiseCoefficients(A) ; // For CT18NLO

  // Output file handling
  std::vector<std::string> filenames =
  {
    "LOXsec.txt", "NLOXsec.txt", "ResumXsec_C.txt", "ExpandedXsec.txt",
    "ResumXsec_Mellin.txt", "ResumXsec_Analitic.txt", "test.txt", "int_res.txt"
  };
  std::vector<std::ofstream> files(filenames.size());
  OpenFiles(files, filenames);

  // Integration precision
  double precision_target = 1e-2;
// Use Pdf Anzat or not:look at main.h;

  // Main computation loop
  for (int ic=0; ic<=3; ++ic)// For all the different center of mass energy
  {
    double sc = s0 * pow(sf / s0, 1.0 * ic / nsteps);// square root of the center of mass energy
    double res1 = 0, err1 = 0, res2 = 0, err2 = 0, res3 = 0, err3 = 0, res4=0, err4=0;
    for (int kID=1; kID<=1; ++kID)//For all the different processes
    { // kiD = 1 NLO Calculation
      // kiD = 2 Calculation at NLO+NLL with PDF Trick numerical seconde derivative
      // kiD = 3 Calculation at Expended NLO+NLL with PDF Trick numerical seconde derivative
      // kiD = 4 Calculation at NLO+NLL with PDF Trick in Mellin space
      // kiD = 5 Calculation at NLO+NLL with PDF Trick analytic seconde derivative
      // Initialisation
      auto start = std::chrono::high_resolution_clock::now();
      Process* proc = new Process(argv[1], sc, kID, ic, usePDFAnzat);
      info("Processing kID = " + std::to_string(kID) + ", ic = " + std::to_string(ic));
      std::cout << 1 << std::endl; 
      // Results
      double res = 0, err = 0,  chi = 0;
      if (kID==1)PerformIntegration(res1, err1, chi, proc, 2, precision_target); 
      else if (kID == 2) PerformIntegration(res2, err2, chi, proc, 3, precision_target);
      else if (kID == 3) PerformIntegration(res3, err3, chi, proc, 3, precision_target);
      else if (kID == 4) PerformIntegration(res4, err4, chi, proc, 1, precision_target);
      else if (kID == 5) PerformIntegration(res, err, chi, proc, 1, precision_target);
      
      // Output file handling
      switch (kID)
      {
        case 0: files[0] << res << "," << sc << std::endl; break;
        case 1: files[1] << res1 << "," << sc << std::endl; break;
        case 2: files[2] << res2 << "," << sc << std::endl; break;
        case 3: files[3] << res3 << "," << sc << std::endl; break;
        case 4: files[4] << res4 << "," << sc << std::endl; break;
        case 5: files[5] << res << "+/-" << err << "," << sc << std::endl; break;

        default: error("Invalid Process ID: kID=" + std::to_string(kID));
      }
      
      // Timer
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      info("kID = " + std::to_string(kID) + ", duration=" + std::to_string(duration.count() / 1e6) + "s, sc = " +  std::to_string(sc));
      delete proc;
    }
          //res1 = 2.00531e+06,err1 = 876.833;
          //std::cout << "////" << "incertitude = " << std::abs(res1-res2)/std::sqrt(err1*err1+err2*err2)  << "  "<< std::abs(res1-res3)/std::sqrt(err1*err1+err3*err3) << " res1 = " << res1 << " res2 = " << res2 << " res3 = " << res3 <<" err1 = " << err1 << " err2 = " << err2 << " err3 = " << err3 << "////" <<   std::endl;
      std::cout << "NLO= " << res1 << " NLL_num= " << res2 << " NLL_exp= " << res3 << " NLL_Mellin= " << res4 << std::endl;
  }

  // End of the program
  CloseFiles(files);
  return 0;
}
