#ifndef PLOTFILLER_HH
#define PLOTFILLER_HH

// Title:     PlotFiller.hh
// Author:    Travis Lamb
// Started:   22nd of June, 2012

// Uses Plot and PhysicsProcess to return a plot for a given process and yield.

// Our libraries
#include "TAMUWW/Tools/interface/Plots.hh"
//#include "TAMUWW/Tools/interface/PlotMap.hh"
//#include "TAMUWW/SpecialTools/interface/ProtectedMap.hh"
#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"
#include "TAMUWW/MEPATNtuple/interface/METree.hh"
#include "TAMUWW/MEPATNtuple/interface/MicroNtuple.hh"
#include "TAMUWW/SpecialTools/interface/MVAVar.hh"
#include "TAMUWW/SpecialTools/interface/ProgressBar.hh"
#include "TAMUWW/Tools/interface/PUreweight.hh"
#include "TAMUWW/Tools/interface/TriggerEfficiency.hh"

//ROOT libraries
#include "TBenchmark.h"

//C++ libraries
#include <iostream>
#include <string>
#include <vector>
#include <map>

// PlotFiller is a modular solution to extracting plots from root files.
// It allows the user to customize how he/she wants to produce the plots.
// The user can set various functions to change how the class calculates cut, weights, etc.
class PlotFiller{

public:

  typedef map<DEFS::LeptonCat, map<std::string,  Plot * > >  MapOfPlots;

   // NOTE That the user must provide the fill function at construction (because it is a required function).
   PlotFiller(MapOfPlots &plotsTemp,
              std::vector<PhysicsProcess*> &procsTemp,
              void (*userFillFuncTemp) (MapOfPlots &, TString, EventNtuple*, METree*, MicroNtuple*, vector<TString>, std::map<TString,MVAVar>&, double) );
   ~PlotFiller();
   
   // Simple functions to change the functionality of the code.
   void setWeightFunction(double (*userWeightFuncTemp) (EventNtuple*, MicroNtuple*, const PhysicsProcess*));
   void setCutFunction(bool (*userCutFuncTemp) (EventNtuple*, MicroNtuple*, const PhysicsProcess*, map<TString,MVAVar>&, vector<TString>));
   void setProcessFunction(void (*userProcessFuncTemp) (EventNtuple*, const PhysicsProcess*));
   void setInitializeEventFunction(void (*userInitEventFuncTemp) (EventNtuple*, const PhysicsProcess*));
   void setMVAWeightDir(TString MVAWD);
   void setMVAMethods(vector<TString> MVAM);
   void setMVAVar(vector<TString> mvav);
   void setMVASpec(vector<TString> mvas);
   void setLimitBranches(int lb);
   
   // Debug functions
   void setMaximumEventsDEBUG(unsigned int maxEvts);
   
   // Runs the events and produces the output.
   void run();
   
private:
   // NOTE That the plots and processes are references.
   MapOfPlots &plots;
   std::vector<PhysicsProcess*> &processes;
   unsigned int numberOfEvents;
   unsigned int debugNumberOfEvents;
   bool debug;
   int limitBranches;
   TString MVAWeightDir;
      vector<TString> MVAMethods, MVAV, MVAS;
   std::map<TString,MVAVar> MVAVars;
   TBenchmark* event_benchmark;
   TBenchmark* func_benchmark;

   // These are the custom functions.
   // Fills the plots
   void (*userFillFunc) (MapOfPlots &, TString, EventNtuple*, METree*, MicroNtuple*, vector<TString>, std::map<TString,MVAVar>&, double);
   // Returns a double that will multiply the weight
   double (*userWeightFunc) (EventNtuple*, MicroNtuple*, const PhysicsProcess*);
   // Returns true if the event passes the cut
   bool (*userCutFunc) (EventNtuple*, MicroNtuple*, const PhysicsProcess*, map<TString,MVAVar>&, vector<TString>);
   // This function is called once for each process before the events are run
   void (*userProcessFunc) (EventNtuple*, const PhysicsProcess*);
   // This function is called once for each event before any cuts are made
   void (*userInitEventFunc) (EventNtuple*, const PhysicsProcess*);
   
   // These default functions allow the user to only have to create functions for weights etc that he/she wants to add.
   static bool defaultCutFunc(EventNtuple*, MicroNtuple*, const PhysicsProcess*, map<TString,MVAVar>&, vector<TString>)
   {
      return true;
   }
   static double defaultWeightFunc(EventNtuple*, MicroNtuple*, const PhysicsProcess*)
   {
      return 1.0;
   }
   static void defaultProcessFunc (EventNtuple*, const PhysicsProcess*)
   {
      return;
   }
   static void defaultInitEventFunc (EventNtuple*, const PhysicsProcess*)
   {
      return;
   }
};

#endif
