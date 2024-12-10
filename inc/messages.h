#ifndef MESSAGES_H
#define MESSAGES_H
// ************************************************************************* //
// Declaration of message services                                           //
//                                                                           //
// By Benjamin Fuks - 18.04.2012.                                            //
// ************************************************************************* //

// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <string>     // Strings                        //
// ---------------------------------------------------- //

// ************************************************************************* //
// Definition of the message services class                                  //
// ************************************************************************* //
void InitMessages();
void info(const std::string&);
void warning(const std::string&);
void error(const std::string&);
void print(const std::string&, const std::string&);
void printerr(const std::string&);


#endif

