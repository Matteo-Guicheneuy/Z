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
void PDFInit();
void InitialiseCoefficients(double A[8][8]);
void OpenFiles(std::vector<std::ofstream>&, const std::vector<std::string>&);
void CloseFiles(std::vector<std::ofstream>&);

#endif

