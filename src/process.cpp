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
Process::Process(const std::string&s1, double s2, int s3, const bool usePDFAnzat): sh(s2*s2), xs_id(s3),usePDFAnzat(usePDFAnzat) 
{
   // Printing information
   info("Initializing the collision setup");

   // Parse PDF set from string
   std::istringstream is1(s1);
   if (!(is1 >> pdf)) throw std::invalid_argument("Invalid PDF set string: " + s1);

    // Map PDF identifiers to names
    static const std::unordered_map<std::string, std::string> pdf_map =
      { {"1", "CT18NLO"}, {"2", "CT14lo"}, {"3", "CT14llo"}, {"4", "cteq6"} };
    auto it = pdf_map.find(pdf);
    if (it == pdf_map.end()) throw std::invalid_argument("Unsupported PDF set identifier: " + pdf);

    // Load the PDF grids
    F = LHAPDF::mkPDF(it->second, 0);
}

