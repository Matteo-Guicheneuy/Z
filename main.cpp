// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <chrono>     // Timing service                 //
#include <fstream>    // File streams                   //
#include <iostream>   // In/Out streams                 //

#include <complex>
// -------- Classes ----------------------------------- //
//#include "main.h"
#include "messages.h" // Message services               //
#include "process.h"  // Process setting                //
#include "settings.h" // Setting routines               //
#include "coefficients.h"
#include "constants.h"
#include "luminosity.h"
#include "main.h"
#include "LHAPDF/AlphaS.h"


double AlphaS=0.118;
double U[20][20];
void PerformIntegration(double&, double&, double&, Process*, size_t, double);
double A[9][8];

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

  // PDF170.0, 200.0,
  PDFInit();
  // To use PDF Anzat  see main.h
  for (int i=0; i<20;i++){for (int j=0; j<20;j++){U[i][j]=0.;}};


  // Setup square root of the center of mass energy
  const int nsteps = 10;
  const double s0 = 170.0, sf = 1000.0;
  //const double s0 = 300., sf = 1000.;
  double Sc[]={110.0,140.0, 150.0, 160.0, 170.0, 200.0, 240.0, 290.0, 345.0, 410.0, 490.0, 580.0, 700.0, 840.0, 1000.0, 1500.0, 2000.0, 3000.0, 5000.0}; 
  //Setup PDF fits
  InitialiseCoefficients(A) ; // For CT18NLO

  // Output file handling
  std::vector<std::string> filenames =
  {
    "LOXsec.txt", "NLOXsec.txt", "ResumXsec_C.txt", "ExpandedXsec.txt",
    "ResumXsec_Mellin.txt", "ResumXsec_Analitic.txt", "test.txt", "test2.txt","all.txt", "int_res.txt"
  };
  std::vector<std::ofstream> files(filenames.size());
  OpenFiles(files, filenames);
  // Integration precision
  double precision_target = 5.e-4;
// Use Pdf Anzat or not:look at main.h;
// AlphaS test :
  double Alpha_S[]={0.000118, 0.00118, 0.0118, 0.118};
  //LHAPDF::AlphaS_Analytic alphaS;  
  double sc;
  double Nlo[]={0.0000009175, 0.005267,0.02203, 0.06837,0.1722, 1.262, 6.662, 24.89, 65.00, 142.3, 279.3, 479.3, 803.7, 1240.0, 1790.0, 3694.0, 5744.0, 10050.0, 19060.0};
  //files[6] <<  "[ ";
  //std::flush(files[6]);
  //files[7] <<  "[ ";
  //std::flush(files[7]);
  // Main computation loop
  for (int ic=9; ic<=18; ++ic)// For all the different center of mass energy
  {     
    //for(int j=0; j<=3; j++){
      //files[6] <<  "[ ";
      //std::flush(files[6]);
      //AlphaS=Alpha_S[j];
      //std::cout << "alpha_S= " << AlphaS << std::endl;
    //sc = s0 * pow(sf / s0, 1.0 * ic / nsteps);// square root of the center of mass energy
    sc = Sc[ic];
    double res0=0., err0=0., res1 = 0, err1 = 0, res2 = 0, err2 = 0, res3 = 0, err3 = 0, res4=0, err4=0, res5=0., err5=0., res6=0., err6=0., res7=0., err7=0.;
    for (int kID=0; kID<=7; ++kID)//For all the different processes
    { // kiD = 0 LO Calculation
      // kiD = 1 NLO Calculation
      // kiD = 2 Calculation at NLO+NLL with PDF Trick numerical seconde derivative
      // kiD = 3 Calculation at Expended NLO+NLL with PDF Trick numerical seconde derivative
      // kiD = 4 Calculation at NLO+NLL with PDF Trick in Mellin space
      // kiD = 5 Calculation at NLO+NLL with PDF Trick analytic seconde derivative
      // Initialisation
      auto start = std::chrono::high_resolution_clock::now();
      Process* proc = new Process(argv[1], sc, kID, ic, usePDFAnzat);
      // Define AlphaS from LHAPDF
      Process *p = (Process *)proc;
      AlphaS = p->F->alphasQ(muR);
      info("Processing kID = " + std::to_string(kID) + ", ic = " + std::to_string(ic));
      // Results
      double res = 0, err = 0,  chi = 0;
      if (kID==0) PerformIntegration(res0, err0, chi, proc, 1, precision_target);
      else if (kID==1)PerformIntegration(res1, err1, chi, proc, 2, precision_target); 
      //else if (kID == 2) PerformIntegration(res2, err2, chi, proc, 3, precision_target);
      //else if (kID == 3) PerformIntegration(res3, err3, chi, proc, 3, precision_target);
      //else if (kID == 4) PerformIntegration(res4, err4, chi, proc, 1, precision_target);
      //else if (kID == 5) PerformIntegration(res5, err5, chi, proc, 1, precision_target);
      //else if (kID == 6) PerformIntegration(res6, err6, chi, proc, 1, precision_target);
      else if (kID==7) PerformIntegration(res7, err7, chi, proc, 1, precision_target);
      // Output file handling
      switch (kID)
      {
        case 0: files[0] << res << "," << sc << std::endl; break;
        case 1: files[1] << res1 << "," << sc << std::endl; break;
        case 2: files[2] << res2 << "," << sc << std::endl; break;
        case 3: break;//files[3] << res3 << "," << sc << std::endl; break;
        case 4: break;//files[4] << res4 << "," << sc << std::endl; break;
        case 5: files[5] << res << "+/-" << err << "," << sc << std::endl; break;
        case 6: case 7: break;
        default: error("Invalid Process ID: kID=" + std::to_string(kID));
      }

      // Timer
      auto stop = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      info("kID = " + std::to_string(kID) + ", duration=" + std::to_string(duration.count() / 1e6) + "s, sc = " +  std::to_string(sc));
      delete proc;

    }
        //std::cout << "////" << "incertitude = " << std::abs(res1-res2)/std::sqrt(err1*err1+err2*err2)  << "  "<< std::abs(res1-res3)/std::sqrt(err1*err1+err3*err3) << " res1 = " << res1 << " res2 = " << res2 << " res3 = " << res3 <<" err1 = " << err1 << " err2 = " << err2 << " err3 = " << err3 << "////" <<   std::endl;
      std::cout << "sc= " << sc << " AlphaS= " << AlphaS << std::endl;
      std::cout << "Born= " << res0 << " err= " << err0 << std::endl;
      std::cout << "NLO_xspace= " << res1 << " err= " << err1 << std::endl;
      std::cout << "NLO_Mellin= " << res7 << " err= " << err7 << std::endl;
      std::cout << "Diff= " << res7-res1 << " err= " << err7 << std::endl;
      std::cout << "NLL_num= " << res2 << " err= " << err2 << std::endl;
      std::cout << "NLL_exp= " << res3 << " err= " << err3 << std::endl;
      std::cout << "NLL_Mel= " << res5 << " err= " << err5 << std::endl;
      std::cout << "Exp_Mel= " << res6 << " err= "  << err6 << std::endl;
      std::cout << std::abs((res6-res1)/res0)<< std::endl;
      std::cout << std::abs((res3-res1)/res0)  << std::endl;
      files[6] << "sc= " << sc << " AlphaS= " << AlphaS << std::endl;
      files[6] << "Born= " << res0 << " err= " << err0 << std::endl;
      files[6] << "NLO= " << res1 << " err= " << err1 << std::endl;
      files[6] << "NLO_Mellin= " << res7 << " err= " << err7 << std::endl;
      //files[6] << "NLL_num= " << res2 << " err= " << err2 << std::endl;
      //files[6] << "NLL_exp= " << res3 << " err= " << err3 << std::endl;
      files[6] << "NLL_Mel= " << res5 << " err= " << err5 << std::endl;
      files[6] << "Exp_Mel= " << res6 << " err= "  << err6 << std::endl;
      //files[6] << "(Exp_Mel-NLO)/LO= " <<std::abs((res3-Nlo[ic])/(res3+Nlo[ic]))<< std::endl;
      files[6] << "(Exp-NLO)/LO= " <<std::abs((res3-res1)/(res3+res1)) << std::endl;
      //files[6] << "(Exp_Mel-NLO)/LO= " <<std::abs((res6-Nlo[ic])/(res6+Nlo[ic]))<< std::endl;
      //files[6] << "(Exp_Mel-NLO)/LO= " <<std::abs((res6-res1)/(res6+res1))<< std::endl;
      files[3] << ((res3+5.0*err3+5.0*err1-res1)/(res3+res1)) << "," << std::endl;
      files[4] << ((res3-5.0*err3-5.0*err1-res1)/(res3+res1)) << "," << std::endl;
 
  }


  // End of the program
  CloseFiles(files);
  return 0;
}
