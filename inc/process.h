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
  Process()  { };
  Process(const std::string&, double, int);
  ~Process() { };

  // Kinematical quantities and scales
  double sh;            // Hadronic center of mass energy
  int xs_id;            // Cross section ID 0=LO, 1=NLO, 2=Resummed, 3=Expanded 
  std::string pdf;      // PDF identifier

  // PDF handler
  const LHAPDF::PDF* F;
};
#endif
