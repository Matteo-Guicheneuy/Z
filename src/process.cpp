// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <unordered_map> // Maps                        //
// -------- Classes ----------------------------------- //
#include "process.h"  // Process                        //
#include "messages.h" // Message services               //
// ---------------------------------------------------- //

// ************************************************************************* //
//  Constructor                                                              //
// ************************************************************************* //

Process::Process(cst &s1, double s2, int s3, int ic,const bool usePDFAnzat): sh(s2*s2), xs_id(s3),ic(ic),usePDFAnzat(usePDFAnzat) 
{
   // Printing information
   info("Initializing the collision setup");
   // Parse PDF set from string
   // PDF set
   //pdf=s1;
   std::istringstream is1(s1); is1 >> pdf;

  
   //Loading the PDF grids
   std::string name;
   if(pdf==1) name = "CT18NLO";
   else if(pdf==2) name = "CT14lo";
   else if (pdf==3) name = "CT14llo";
   else if (pdf==4) name = "cteq6";

   F=LHAPDF::mkPDF(name,0);

}

