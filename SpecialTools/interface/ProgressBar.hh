#ifndef PROGRESSBAR_HH
#define PROGRESSBAR_HH

///////////////////////////////////////////////////////////////////
// ProgressBar.hh
// --------------
//
//            03/19/2014 Alexx Perloff  <aperloff@physics.tamu.edu>
///////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <string>

using std::cout;
using std::setw;
using std::flush;
using std::string;

class ProgressBar {
public:

  // C Versions
  static void loadbar(unsigned int x, unsigned int n, unsigned int w = 50); //static inline?
  static void loadbar2(unsigned int x, unsigned int n, unsigned int w = 50); //static inline?
  static void loadbar3(string prefix, unsigned int x, unsigned int n, unsigned int w = 50); //static inline?
  // C++ Version
  static void loadBar(int x, int n, int r, int w); //static inline?

//ClassDef (ProgressBar,1)

};
#endif
