// Name: PhysicsProcessNEW.cc
// Author: Travis Lamb

// This class will eventually replace an older Physics class. It is a base class that simply holds data members.

#include "TFile.h"

#include "TAMUWW/Tools/interface/PhysicsProcessNEW.hh"

using namespace std;

// ##################################################
// ################ PHYSICS PROCESS #################
// ##################################################

PhysicsProcessNEW::PhysicsProcessNEW (string procName,
                                      string fileName,
                                      double cross_section,
                                      double lum, 
                                      unsigned int in_ev,
                                      string treeName):
   name(procName),
   chain (0),
   sigma(cross_section),
   intLum(lum),
   initial_events(in_ev)
{
  
   TFile * file = TFile::Open(fileName.c_str());
   if (!file->IsOpen())
   {
      cout << "ERROR proc::proc() could not open file " << fileName << endl;
      return;
   }

   if (!file->cd("PS"))
   {
      cout << "ERROR proc::proc() could not CD into directory PS in file " << fileName << endl;
      return;
   }
  
   chain = (TChain*) file->Get(treeName.c_str());
   if (chain == 0)
   {
      cout << "ERROR proc::proc() could not find tree named " << treeName << " in file " << fileName << endl;
      return;
   }
}

// ##################################################
// ############ COLORED PHYSICS PROCESS #############
// ##################################################

ColoredPhysicsProcessNEW::ColoredPhysicsProcessNEW (string procName,
                                                    string fileName,
                                                    double cross_section,
                                                    double lum, 
                                                    unsigned int in_ev,
                                                    int col,
                                                    string treeName):
   PhysicsProcessNEW(procName, fileName, cross_section, lum, in_ev, treeName),
   color(col)
{
   
}