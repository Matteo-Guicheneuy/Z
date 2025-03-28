#ifndef PROCESS_H
#define PROCESS_H
// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <string>    // String streams                  //
// -------- Classes ----------------------------------- //
#include "LHAPDF/LHAPDF.h"
// -------- Processor variables ----------------------- //
#define cst const std::string                           //
// ---------------------------------------------------- //

// ************************************************************************* //
// Definition of the class Process                                           //
// ************************************************************************* //
class Process
{
public:
  Process() : usePDFAnzat(false) { };
  Process(cst &, double, int, int, const bool);
  ~Process() { };

  // Kinematical quantities and scales
  double sh;            // Hadronic center of mass energy
  int xs_id;            // Cross section ID 0=LO, 1=NLO, 2=Resummed, 3=Expanded 
  int ic;
  int pdf;              // PDF identifier
  const bool usePDFAnzat;     // Use PDF Anzat or not

  // PDF handler
  const LHAPDF::PDF* F;
};
#endif
