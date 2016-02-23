//Our libraries
#include "TAMUWW/Tools/interface/CreateCondorScriptME.hh"
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
//// Calls CreateCondorScript_runME for each process, with the option of running on JES shifted files as well (runOnJESShifted=true).
//// ScriptName=LME+PLabel+suffix, OutputName=PLabel+ME+suffix; where suffix="","JESp1s","JESm1s"
//// rootInputName=inputmode+PLabel+filenameSuffix
//// Set checkCompletedJobsInstead=true to check completed jobs (in the directory completedOutputDir) and set createCompletionScripts=true to create scripts to rerun the missing jobs.
 int main(int argc,char**argv) {
   // evaluate command-line / configuration file options
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;

   TString         inputmode                  = cl.getValue<TString>  ("inputmode");
   TString         scriptDir                  = cl.getValue<TString>  ("scriptDir");
   TString         rootInputDir               = cl.getValue<TString>  ("rootInputDir");
   int             nEventsPerJob              = cl.getValue<int>      ("nEventsPerJob");
   const TString   globalRunMESuffix          = cl.getValue<TString>  ("globalRunMESuffix",            "");// Use "Wgg" when running over additional MEs, set to "" by default
   bool            checkCompletedJobsInstead  = cl.getValue<bool>     ("checkCompletedJobsInstead", false);
   TString         completedOutputDir         = cl.getValue<TString>  ("completedOutputDir",         "./");
   bool            createCompletionScripts    = cl.getValue<bool>     ("createCompletionScripts",   false);
   bool            verbose                    = cl.getValue<bool>     ("verbose",                   false);
   int             nMaxJobs                   = cl.getValue<int>      ("nMaxJobs",                     -1);
   TString         filenameSuffix             = cl.getValue<TString>  ("filenameSuffix",          ".root");
   bool            pbsScripts                 = cl.getValue<bool>     ("pbsScripts",                false);
   TString         pbsQueue                   = cl.getValue<TString>  ("pbsQueue",               "hepxrt");
   TString         pbsWalltimeStr             = cl.getValue<TString>  ("pbsWalltimeStr",       "48:00:00");
   TString         pbsUser                    = cl.getValue<TString>  ("pbsUser",                      "");
   TString         TDirectoryFilename         = cl.getValue<TString>  ("TDirectoryFilename",         "PS");
   bool            useSpecificEvents          = cl.getValue<bool>     ("useSpecificEvents",         false);
   TString         specificEventsFileLocation = cl.getValue<TString>  ("specificEventsFileLocation",  "1");
   vector<TString> processes                  = cl.getVector<TString> ("processes", "ggH125:::qqH125:::WH125:::WW:::WZ:::ZZ:::ZJets:::TTbar:::STopS_T:::STopS_Tbar:::STopT_T:::STopT_Tbar:::STopTW_T:::STopTW_Tbar:::WJets");

   if (!cl.check()) return 0;
   cl.print();
   
   TBenchmark* m_benchmark = new TBenchmark();
   m_benchmark->Reset();
   m_benchmark->Start("event");

   if(!scriptDir.EndsWith("/")) scriptDir += "/";
   if(!rootInputDir.EndsWith("/")) rootInputDir += "/";
   if(!completedOutputDir.EndsWith("/")) completedOutputDir += "/";


   bool insideTDirectoryFile = true;
   TString scriptBaseName,rootinputBaseName,outputBaseName;
   int nJobs;
   double nEntries;

   for (UInt_t np=0; np<processes.size();np++) {
      cout << "Process=" << processes[np] << endl;
      scriptBaseName=globalRunMESuffix;
      scriptBaseName=processes[np]+scriptBaseName;
      scriptBaseName=inputmode+scriptBaseName;
      rootinputBaseName=inputmode;
      rootinputBaseName=processes[np]+rootinputBaseName;
      outputBaseName=inputmode;
      outputBaseName=globalRunMESuffix+outputBaseName;
      outputBaseName=processes[np]+outputBaseName;

      TFile* infile;
      if(pbsScripts)
         infile = TFile::Open(rootInputDir+rootinputBaseName+filenameSuffix);
         //infile = new TFile(rootInputDir+processes[np]+"/"+rootinputBaseName+filenameSuffix);
      else
         infile = TFile::Open(rootInputDir+rootinputBaseName+filenameSuffix);
      TTree* InTree;

      if ( !insideTDirectoryFile ) {
         InTree = (TTree*)infile->Get("jets2p");
      } else {
         TDirectoryFile *tdf = (TDirectoryFile*)infile->Get(TDirectoryFilename);
         InTree = (TTree*)tdf->Get("jets2p");
      }

      nEntries=InTree->GetEntries();
      nJobs=ceil(nEntries/nEventsPerJob);
      nEventsPerJob=nEventsPerJob;
      if ( (nMaxJobs>0) && (nJobs>nMaxJobs) ) {
         // Only work with the first nMaxJobs jobs and nMaxJobs*nEventsPerJob entries
         cout << "Restricting the number of jobs (and entries) which will be considered" << endl;
         nJobs=nMaxJobs;
         nEntries=nJobs*nEventsPerJob;
      }

      cout << "nEntries=" << nEntries ;
      cout << ", nJobs=" << nJobs << endl;


      CreateCondorScriptME* ccs = new CreateCondorScriptME();
      ccs->setTotalJobs(nJobs);
      ccs->setTotalEvents();
      ccs->setEventsPerJob(nEventsPerJob);
      ccs->setFirstJob(0);
      ccs->setTreeName();
      ccs->setScriptDir(scriptDir);
      ccs->setScriptNameSuffix(scriptBaseName);
      ccs->setRootInputDir(rootInputDir);
      ccs->setRootInputName(rootinputBaseName+filenameSuffix);
      ccs->setOutputName(outputBaseName);
      ccs->setDisjointJob();
      ccs->setDisjointJobString();
      ccs->setGlobalRunMESuffix(globalRunMESuffix);
      ccs->setCompletedOutputDir(string(completedOutputDir));
      ccs->setSpecificEventsFile(useSpecificEvents,specificEventsFileLocation);
      if ( !checkCompletedJobsInstead ) {
         if (pbsScripts)
            ccs->writePBSScripts(pbsQueue, pbsWalltimeStr, pbsUser);
         else
            ccs->writeCondorScripts();
      }
      else {
         cout << "Checking Completed Jobs" << endl;
         TString infileprefix=globalRunMESuffix;
         infileprefix=inputmode+infileprefix;
         infileprefix=processes[np]+infileprefix;
         infileprefix=completedOutputDir+infileprefix;
         
         //cout << "infileprefix=" << infileprefix << endl;
         //CheckCONDOROutput(infileprefix,invalidFileStr,nInvalidFiles,0,nJobs-1,true,NEvtsPerJob,nEntries,false);
         ccs->checkForBadFiles(verbose);
         ccs->checkForMissingFiles();
         pair<TString,int> invalidFiles = ccs->getInvalidFiles();
         if ( createCompletionScripts && (invalidFiles.second>0) ) {
            cout << "Creating the completion script" << endl;
            CreateCondorScriptME* ccs_completion = new CreateCondorScriptME();
            ccs_completion->setTotalJobs(invalidFiles.second);
            ccs_completion->setTotalEvents();
            ccs_completion->setEventsPerJob(nEventsPerJob);
            ccs_completion->setFirstJob(0);
            ccs_completion->setTreeName();
            ccs_completion->setScriptDir(scriptDir);
            ccs_completion->setScriptNameSuffix(scriptBaseName+"Completion");
            ccs_completion->setRootInputDir(rootInputDir);
            ccs_completion->setRootInputName(rootinputBaseName+filenameSuffix);
            ccs_completion->setOutputName(outputBaseName);
            ccs_completion->setDisjointJob(true);
            ccs_completion->setDisjointJobString(invalidFiles.first);
            ccs_completion->setGlobalRunMESuffix(globalRunMESuffix);
            ccs_completion->setSpecificEventsFile(useSpecificEvents,specificEventsFileLocation);
            ccs->setCompletedOutputDir(string(completedOutputDir));
            if (pbsScripts)
               ccs_completion->writePBSScripts(pbsQueue, pbsWalltimeStr, pbsUser);
            else
               ccs_completion->writeCondorScripts();

         }
      }

   }

   m_benchmark->Stop("event");
   cout << "CreateCondorScript_runME_x" << endl << "\tCPU time = " << m_benchmark->GetCpuTime("event") << " s" << endl
        << "\tReal time = " << m_benchmark->GetRealTime("event") << " s" << endl;
   delete m_benchmark;

   return 0;

}//main

