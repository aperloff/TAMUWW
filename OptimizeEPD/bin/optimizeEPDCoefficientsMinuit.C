
// This macro is intended to calculate the most discriminating EPD coefficients
// Author: Ricardo Eusebi, November 10th, 2008
//       : Charlie Zhang, July 2013.

//ROOT libraries
#include <TStyle.h>
#include "TDirectory.h"
#include "TH1F.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TMath.h"
#include "TGraph.h"

//C++ libraries
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

//Other libraries
#include "TAMUWW/OptimizeEPD/interface/PhysicsProcessForOpt.hh"
#include "TAMUWW/OptimizeEPD/interface/EPDOptimizer.hh"
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"
#include "TAMUWW/SpecialTools/interface/FigureOfMerit.hh"
#include "TAMUWW/MEPATNtuple/interface/ProbsForEPD.hh"
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"



using std::cout;
using std::endl;
using std::string;
using std::map;
using std::vector;
using std::ofstream;
using std::stringstream;

bool debug = false;
// ===================================================================
double WWEPD(const ProbsForEPD & probs){
  return (probs.ww + probs.wz)/(probs.schan + probs.tchan + probs.tchan2 + probs.wbb +
				probs.wc + probs.qcd + probs.tt + probs.ww + probs.wz);
}

// ===================================================================
double HiggsEPD(const ProbsForEPD & probs){
  return (probs.wh + probs.hww)/(probs.schan + probs.tchan + probs.tchan2 + probs.wbb +
				 probs.wc + probs.qcd + probs.tt + probs.ww + probs.wz +
				 probs.wh + probs.hww);
}


// ===================================================================
// This method gets all the PhysicsProcesses
// and loads the relevant info into memory. 
vector<PhysicsProcessForOpt*> loadProcessesIntoMemory(DEFS::Ana::Type anaType, DEFS::JetBin jetBin, DEFS::TagCat tagcat){
  
  // Get the physical processes for WW. The false means "data" is not one of the processes. 
  vector<PhysicsProcess*> phy_processes;
  if(debug){ //needs to be set to HiggsAnalysis
    //phy_processes = DefaultValues::getProcessesTest(jetBin, tagcat, false, false);
  } else if(anaType == DEFS::Ana::WWAnalysis){
    phy_processes = DefaultValues::getProcessesWW(jetBin, tagcat, false, false, DEFS::MicroNtuple);
  } else if(anaType == DEFS::Ana::HiggsAnalysis){
    phy_processes = DefaultValues::getProcessesHiggs(jetBin, tagcat, false, false, false, DEFS::MicroNtuple);
  }

  // Now construct a set of PhysicsProcessForOpt out of them
  vector<PhysicsProcessForOpt*> processes ;
  for (unsigned int pp = 0 ; pp < phy_processes.size(); pp++){
    PhysicsProcessForOpt* ppo = new PhysicsProcessForOpt(PhysicsProcessMemory(*phy_processes[pp]));
    //    ppo->fillMETreeIndexMap();
    if(anaType == DEFS::Ana::WWAnalysis){
      ppo->setEPDFunction(WWEPD);
    } else if(anaType == DEFS::Ana::HiggsAnalysis){
      ppo->setEPDFunction(HiggsEPD);
    }
    processes.push_back(ppo);
  }
 
  return  processes ;

}//loadProcessesIntoMemory


int main(int argc,char**argv){
  //
  // evaluate command-line / configuration file options
  //
  CommandLine cl;
  if (!cl.parse(argc,argv)) return 0;
  string          ofilename        = cl.getValue<string>       ("ofilename",        "EPDCoefficients.txt");
  //string          lepCat           = cl.getValue<string>       ("lep",              "both");
  //DEFS::LeptonCat leptonCat        = DEFS::getLeptonCat(lepCat);
  string          jetBinS          = cl.getValue<string>       ("jetBin",           "jets2");
  DEFS::JetBin    jetBin           = DEFS::getJetBin(jetBinS);
  string          tagcatS          = cl.getValue<string>       ("tagcat",           "pretag");
  DEFS::TagCat    tagcat           = DEFS::getTagCat(tagcatS);
  string          anaTypeS         = cl.getValue<string>       ("anaType",          "HiggsAnalysis");
  DEFS::Ana::Type anaType          = DEFS::Ana::getAnaType(anaTypeS);

  if (!cl.check()) return 0;
  cl.print();

  // Select the analysis type, only to report to the screen...
  //DEFS::Ana::Type anaType = DEFS::Ana::HiggsAnalysis;

  cout<<endl;
  cout<<"----------------------------------------------"<<endl;
  cout<<" OPTIMIZING  coefficients for "<<DEFS::Ana::getTypeString(anaType)<<endl;
  cout<<"----------------------------------------------"<<endl<<endl; 


  // jetbin we are optimizing over
  //DEFS::JetBin jetBin = DEFS::jets2;

  // our category is pretag
  //DEFS::TagCat tagcat = DEFS::TagCat::pretag;

  // load all the processes into memory
  vector<PhysicsProcessForOpt*>  processes  = loadProcessesIntoMemory(anaType, jetBin, tagcat);
  cout<<"loadProcessIntoMemory is done"<<endl;

  // do the actual computation of best EPD coefficients for the given processes
  // last parameter is debug flag. set to false.
  EPDOptimizer opt(anaType, tagcat, processes, 125, DEFS::jets2, false);  //125 = mhiggs
  opt.optimize();

  // Optimize the coefficients for the pretag
  ProbsForEPD coeffs0 = opt.getEPD();
  // The resulting figure of merit of the optimization for PreTag,=1T and >= 2 tags
  double figOfMerit0 = opt.getFOM();

  // Print the best EPD coefficientsnorm
  cout<<endl<<" BEST Coefficients with fig of merit="<< figOfMerit0<<" follow"<<endl;
  cout<<coeffs0.getStringConstructor()<<endl;


  // - - - - - - - - - - - - - - - - - - - - - - - - - -
  // The rest just create and print the optimized 
  //     coefficients to the output file
  // - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Create the outfile; no need to close it as it will be destructed automatically.
  ofstream outFile(ofilename);

  // Print Instructions
  outFile <<" You can copy and paste from this to MicroNtuple.cc"<<endl;
  outFile << endl << endl << endl;
  outFile << coeffs0.getStringConstructor()<<" // DEFS::TagCat="<<DEFS::pretag<<" figOfMerit="<<figOfMerit0<<endl;
  
  // close the outFile
  outFile.close();

  // Free memory
  for (unsigned p=0;p<processes.size();p++)
    delete processes[p];

 
} // main



  
