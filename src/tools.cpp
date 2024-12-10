// ************************************************************************* //
// Tools for std::string manipulations                                            //
//                                                                           //
// By Benjamin Fuks - 31.05.2012.                                            //
// ************************************************************************* //

// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cmath>      	// Mathematical functions       //
#include <string>     // Strings                        //
#include <sstream>    // String streams                 //
// -------- Classes ----------------------------------- //
#include "messages.h"        // Message services        //
// ---------------------------------------------------- //

// ************************************************************************* //
// Printing a xsection                                                       //
// ************************************************************************* //
void DisplayXsec(const double &r, const double &e, const std::string &tag)
{
  std::ostringstream osr; osr << r; std::string rstr= osr.str();
  std::ostringstream ose; ose << e; std::string estr= ose.str();
  double prec=0.;
  if(r!=0.) prec=std::abs(100.*e/r);
  std::ostringstream opr; opr.precision(2); 
  opr << std::fixed << prec; std::string prstr= opr.str();
  info(tag + " results: " + rstr + " +/- " + estr + " pb  (@ " + prstr + "%)");
}

