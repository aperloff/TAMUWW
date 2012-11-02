// This Program Takes input Monte Carlo histograms of several processes
// and fits them to a data histogram. It only scales any two processes.
// 
// Written by Alexx Perloff
// Started on the 16th of October, 2012

// C++ Libraries
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <utility>

// Root Libraries
#include "TROOT.h"
#include "TString.h"
#include "TKey.h"
#include "TPad.h"
#include "TH1.h"
#include "THStack.h"
#include "TFile.h"
#include "TMath.h"
#include "TMinuit.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

// Root Math Core headers
#include <Math/SpecFuncMathCore.h>
#include <Math/PdfFuncMathCore.h>
#include <Math/ProbFuncMathCore.h>

// TAMU Libraries
#include "TAMUWW/Tools/interface/Plots.hh"
#include "TAMUWW/SpecialTools/interface/FigureOfMerit.hh"

namespace HistogramsFitter
{
   std::map<std::string, TH1D*> monteCarloHistograms;
   TH1D* dataHistogram;
   //names of the processes that will vary
   std::vector<std::string> fProcessNames;
   TH1D* signalHistogram;
   TH1D* backgroundHistogram;
}

//##################################################
//########## BACKGROUND ESTIMATION CLASS ###########
//##################################################

class Fitter
{
private:
   //read location
   std::string inRootFileLocation;
   //write location
   std::string outRootFileLocation;
   // The lepton and object we will perform on
   std::string leptonName;
   std::string objectName;
   
   //names of the histograms
   std::vector<std::string> histNames;
   
   //output histograms
   THStack* resultStack;
   
   //Debug options
   unsigned int rebinSizeDEBUG;
   bool debug;
   
   //std::pair<double, double> scaleParameters;
   std::vector<double> scaleParameters;
   double reducedChiSquared;
   
   //These are used to clone the input canvas
   TCanvas* canvas;
   TPad* pad;
   THStack* mcStack;
   THStack* dataStack;
   TList* mcList;
   TList* dataList;
   
public:
   Fitter(std::string lepton, std::string object, std::vector<std::string> fproc, std::string inFileLoc, std::string outFileLoc = "");
   ~Fitter();
   
   //reads histograms from inRootFileLocation
   void readHistograms();
   //performs the fitting
   void fitHistograms();
   //writes the results to outRootFileLocation
   void writeHistograms();
   //adds all the signal histograms together and adds all of the background histograms together
   void addSigBkgHistograms(std::vector<std::string> sig);
   
   //These are the functions that return chi^2
   static double fitFunc(const double *par);
   
   //Changes the write location
   void setWriteLocation(std::string outFileLoc);
   //Changes the write location
   void setReadLocation(std::string inFileLoc);
   //A debug option that rebins the histograms as it reads them in.
   void setRebinSizeDEBUG(unsigned int rebinSize);
   
   // Returns a pair with the parameters
   // first is process1 scale factor, second is process2 scale factor
   //std::pair<double, double> getParameters();
   std::vector<double> getParameters();
   // Returns the Chi2/NDF
   double getReducedChiSquared();
   // Returns a FOM for the specified signal and background
   double getFOM(double FOM);
   // Finds if a specific string exists in a vector
   int vfind(std::vector<std::string> a, std::string b);
   
private:
   void initializeFileLocations(std::string inFileLoc, std::string outFileLoc);
   void initializeHistNames();
   //Returns the numbers that we scale by to fit data
   //double[0] is the scale factor for process1, double[1] is for process2
   const double* fitAndReturnParameters();
   //Colors each monteCarlo Process
   void colorMonteCarlo();
   //This fills resultStack
   void produceResultStacks(const double* pars);
   //Returns a vector filled with the different plot names (including MET)
   // e.g. "AngleJ1J2", "MET", "j1Pt_Mjj"
   std::vector<std::string> getPlotNames();
};
