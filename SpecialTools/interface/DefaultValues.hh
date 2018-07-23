#ifndef DefaultValues_HH
#define DefaultValues_HH

//ROOT libraries
#include "TROOT.h"
#include "TSystem.h"
#include "TError.h"
#include "TMath.h"
#include "TBenchmark.h"
#include "TList.h"
#include "TString.h"
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"

//C++ libraries
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <utility>
#include <assert.h>
#include <sstream>


//Our libraries
#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/SpecialTools/interface/Value.hh"
#include "TAMUWW/SpecialTools/interface/Table.hh"
#include "TAMUWW/SpecialTools/interface/TableRow.hh"
#include "TAMUWW/SpecialTools/interface/TableCellVal.hh"
#include "TAMUWW/SpecialTools/interface/TableCellInt.hh"
#include "TAMUWW/SpecialTools/interface/TableCellText.hh"
#include "TAMUWW/SpecialTools/interface/PhysicsProcess.hh"
#include "TAMUWW/SpecialTools/interface/FileLocationTable.hh"

using std::string;
using std::vector;
using std::map;
using std::cout;
using std::endl;
using std::pair;
using std::make_pair;

class DefaultValues{

public:

   // Returns the path to the config files
   // basically returns $CMSSW_BASE+"/src/TAMUWW/Config/Official/"
   static std::string getConfigPath();
  
   static std::string getBDTLocation(DEFS::JetBin jetBin, DEFS::TagCat tagcat, DEFS::University univ, bool filenameOnly = false);

   static std::string getWeightForCategory(DEFS::TagCat tagcat,  
                                           DEFS::PhysicsProcessType, int );

   // Return a table with normalization of process for the given tag category
   static Table getNormTable(DEFS::LeptonCat evtcat, DEFS::TagCat tagcat);

   // Return a table with location of files for the given tag category
   static Table getFileLocationTable(DEFS::TagCat tagcat);

   // Return a vector of PhysicsProcess's for the requested processName's
   // Each PhysicsProces has norm obtained from getNormTable, and 
   // ntuples chained from files in getFileLocationTable
   // Calls the getSingleProcess recursively
   static std::vector < PhysicsProcess * > getProcesses(std::vector<DEFS::PhysicsProcessType> processName,
                                                        DEFS::JetBin jetBin, 
                                                        DEFS::TagCat tagcat,
                                                        bool for_plots,
                                                        DEFS::NtupleType ntuple_type);

   
   // Return a single PhysicsProces with norm obtained from normTable, and 
   // ntuples chained from files in fileTable
   static PhysicsProcess * getSingleProcess(DEFS::PhysicsProcessType process,
                                            DEFS::JetBin jetBin,
                                            std::map<DEFS::LeptonCat,Table> normTable,
                                            Table fileTable,
                                            bool for_plots,
                                            DEFS::NtupleType ntuple_type);



  
   //GIVE THE OFFICIAL LIST OF PROCESSES FOR THE GIVEN ANALYSIS
   // DON'T MESS WITH THIS METHOD, CREATE ANOTHER ONE FOR TEST PURPOSES
   static std::vector < PhysicsProcess * > getProcessesWW(DEFS::JetBin,
                                                          DEFS::TagCat,
                                                          bool include_data,
                                                          bool for_plots,
                                                          DEFS::NtupleType ntuple_type);

   static std::vector < PhysicsProcess * > getProcessesHiggs(DEFS::JetBin,
                                                             DEFS::TagCat,
                                                             bool include_data,
                                                             bool include_systematics,
                                                             bool for_plots,
                                                             DEFS::NtupleType ntuple_type,
                                                             DEFS::LeptonCat lepton = DEFS::both);

   static std::vector < PhysicsProcess * > getProcessesTest(DEFS::JetBin jetBin,
                                                            DEFS::TagCat tagcat, 
                                                            bool include_data,
                                                            bool forPlots,
                                                            DEFS::NtupleType ntuple_type,
                                                            DEFS::LeptonCat lepton = DEFS::both);

   // Returns the cross section for the given process
   static pair<double,double> getCrossSectionAndError(TString channelName);
   
   // Returns the branching ratio for the given process
   static double getBranchingRatio(TString channelName);
   
   // Returns the number of Monte Carlo events for the given process
   static double getNumMCEvts(TString channelName);

   // Returns the scale factor for the given process
   static double getScaleFactor(TString channelName);
   static double getScaleFactor(TString channelName, DEFS::LeptonCat leptonCat);

   // Returns the maximum tEventProb for a given ME found in that MEs matching 
   // MC sample (i.e. the WW ME in the WW MC sample).
   //This is not a unique identifier. It will mix up STopT and STopS in STopT_T/STopT_Tbar and STopS_T/STopS_Tbar respectively.
   //Takes the highest value.
   static pair<double,double> getMaxEventProbAndError(int probStatIndex, Table* inputTable = 0);
   //Not unique. Will mix up WLg and WLgSub. Takes the highest value.
   static pair<double,double> getMaxEventProbAndError(DEFS::PhysicsProcessType ppType, string meType, Table* inputTable = 0);
   //Alias for the previous member function so it can be used in CINT
   static pair<double,double> getMaxEventProbAndError(string ppType, string meType, Table* inputTable = 0);
   //Not unique. Will mix up WLg/WLgSub from WJets, STopT from STopT_T and STopT_Tbar, and STopS from STopS_T and STopS_Tbar.
   //Takes the highest value.
   static pair<double,double> getMaxEventProbAndError(string meType, Table* inputTable = 0);

   // Returns the median BDT value for a given jet bin and lepton category
   static pair<double,double> getMedianPurity(DEFS::JetBin jetBin, DEFS::LeptonCat leptonCat, string BDTType, Table* inputTable = 0);

   // Sets two vectors with the variables and spectators for a given MVA training
   static void getMVAVar(TString filename, vector<TString>& MVAV, vector<TString>& MVAS);

   // Returns the index of the location of b within a, or -1 if b is not found in a
   static int vfind(vector<string> a, string b);
   static int vfind(vector<TString> a, TString b);

   // Remove a substring p from a given string s
   template<typename T>
   void removeSubstrs(std::basic_string<T>& s, const std::basic_string<T>& p);

   // Destroy all open canvases
   static void DestroyCanvases();

   // Returns an object from any ROOT file in TAMUWW/ConfigFiles/Official/
   static TObject* getConfigTObject(TString histoFile, TString hname, TString newName);
   static TH1* getConfigTH1(TString histoFile, TString hname, TString newName);
   static TH2* getConfigTH2(TString histoFile, TString hname, TString newName);

   // Returns a rebinned TH2 object
   // There is still a problem with rebinning the underflow and overflow bins
   // Everything else works though
   TH2* Rebin2D(TH2* old, Int_t nxgroup, Int_t nygroup, const char*newname, const Double_t *xbins, const Double_t *ybins, TString Options = "");

   // Tests the Rebin2D function
   void Rebin2DTest(TString Options = "");

   //prints a TBenchmark summary with variable precision and a list of timers
   static void printSummary(TBenchmark* bench, int precision, Float_t &rt, Float_t &ct, vector<string> timers);

};
 
#endif
