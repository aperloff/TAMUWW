//Our libraries
#include "TAMUWW/Tools/interface/MicroNtupleMaker.hh"
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
 int main(int argc,char**argv) {
   // evaluate command-line / configuration file options
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;

   vector<TString> inputPaths                = cl.getVector<TString> ("inputPaths",   "/uscms_data/d3/ilyao/Winter12to13ME8TeV/MEResults/rootOutput/");
   vector<TString> addDir                    = cl.getVector<TString> ("addDir",                                                                    "");
   bool            useXROOTD                 = cl.getValue<bool>     ("useXROOTD",                                                              false);
   TString         outputPath                = cl.getValue<TString>  ("outputPath",  "/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples/");
   int             largeProcessCase          = cl.getValue<int>      ("largeProcessCase",                                                           0);
   //TString         smallProcessLabel         = cl.getValue<TString>  ("smallProcessLabel",                                                  "ggH170");
   vector<TString> processes                 = cl.getVector<TString> ("processes", "ggH125:::qqH125:::WH125:::WW:::WZ:::ZZ:::ZJets:::TTbar:::STopS_T:::STopS_Tbar:::STopT_T:::STopT_Tbar:::STopTW_T:::STopTW_Tbar:::WJets");
   TString         mergeNewEventNtuple       = cl.getValue<TString>  ("mergeNewEventNtuple",                                                       "");
   bool            saveMissingEvents         = cl.getValue<bool>     ("saveMissingEvents",                                                      false);
   bool            fillBDT                   = cl.getValue<bool>     ("fillBDT",                                                                false);
   bool            debug                     = cl.getValue<bool>     ("debug",                                                                  false);

   if (!cl.check()) return 0;
   cl.print();

   // Trying to speed up the code
   gEnv->SetValue("TFile.AsyncPrefetching", 1);
   
   TBenchmark* m_benchmark = new TBenchmark();
   m_benchmark->Reset();
   m_benchmark->Start("event");

   // This loads the TMVA library
   TMVA::Tools::Instance();

   if(inputPaths.size()>1 && inputPaths.size()!=processes.size()) {
      cout << "WARNING::makeMicroNtuple_x The size of the inputPaths vector and the processes vector are not the same" << endl
           << "Make sure you meant to merge the ROOT files in multiple folders into one microNtuple." << endl;
   }
   else if(inputPaths.size()==1 && inputPaths.size()!=processes.size()) {
      cout << "WARNING::makeMicroNtuple_x The size of the inputPaths vector is 1 while the number of processes is greater." << endl
           << "The inputPath will be clones to match the number of processes and all ROOT files will come from the same folder." << endl
           << "Make bloody sure this is what you want!" << endl;
      vector<TString> tmp(processes.size(),inputPaths.front());
      inputPaths = tmp;
   }
   else if(addDir.size()>0 && (inputPaths.size()>1 || processes.size()>1)) {
      cout << "WARNING::makeMicroNtuple_x The size of addDir is greater than 0." << endl
           << "Make bloody sure this is what you want!" << endl;
   }

   if(!outputPath.EndsWith("/")) outputPath += "/";
   for(unsigned int ip=0; ip<inputPaths.size(); ip++)
      if(!inputPaths[ip].EndsWith("/")) inputPaths[ip] += "/";

   MicroNtupleMaker* mnm;
   for(unsigned int i=0; i<processes.size(); i++) {
      mnm = new MicroNtupleMaker(); 
      mnm->setEventNtuplePath(mergeNewEventNtuple);
      mnm->setProcess(processes[i]);
      mnm->setAddDir(addDir);
      mnm->setXROOTD(useXROOTD);
      mnm->setOutputPath(outputPath);
      mnm->setMissingEventsFlag(saveMissingEvents);
      mnm->setFillBDT(fillBDT);
      mnm->setDebug(debug);
      mnm->createMicroNtuple(inputPaths[i],outputPath,largeProcessCase,processes[i]);
      delete mnm;
   }      

   m_benchmark->Stop("event");
   cout << "makeMicroNtuple_x" << endl << "\tCPU time = " << m_benchmark->GetCpuTime("event") << " s" << endl
        << "\tReal time = " << m_benchmark->GetRealTime("event") << " s" << endl;
   delete m_benchmark;

   return 0;

}//main
