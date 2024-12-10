// ************************************************************************* //
// Message services                                                          //
//                                                                           //
// By Benjamin Fuks - 18.04.2012.                                            //
// ************************************************************************* //

// ************************************************************************* //
//  Includes                                                                 //
// ************************************************************************* //
// -------- Standard headers -------------------------- //
#include <cstdlib>    // C Standard General Utilities   //
#include <fstream>    // Read/Write files               //
#include <iostream>   // In/Out streams                 //
// -------- Classes ----------------------------------- //
#include "messages.h"	// Message services		//
// ---------------------------------------------------- //

// ************************************************************************* //
//  Initializing the log file                                                //
// ************************************************************************* //
void InitMessages()
{
   std::ofstream log("messages.log"); log.close();
}

// ************************************************************************* //
//  Printing information                                                     //
// ************************************************************************* //
void info(const std::string &str) { print(str,"INFO"); }
void warning(const std::string &str) { print(str,"WARNING");  }
void error(const std::string &str) { printerr(str); }

void print(const std::string &str, const std::string &tag)
{
   std::ofstream log("messages.log", std::ios::app);
   log << " ** " << tag << ": " << str.c_str() <<"."<< std::endl;
   std::cout << " ** " << tag << ": " << str.c_str() <<"."<< std::endl;
   log.close();
}

void printerr(const std::string &str)
{
   std::cout << " ** ERROR: " << str.c_str() <<"."<< std::endl;
   std::ofstream log("messages.log", std::ios::app);
   log << " ** ERROR: " << str.c_str() <<"."<< std::endl;
   log.close();
   exit(0);
}
