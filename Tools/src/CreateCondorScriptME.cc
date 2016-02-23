// This class libraries
#include "TAMUWW/Tools/interface/CreateCondorScriptME.hh"

//______________________________________________________________________________
CreateCondorScriptME::CreateCondorScriptME() {

}//C'tor

//______________________________________________________________________________
CreateCondorScriptME::~CreateCondorScriptME() {

}//D'tor

//______________________________________________________________________________
// read_directory()
//   Return an ASCII-sorted vector of filename entries in a given directory.
//   If no path is specified, the current working directory is used.
//
//   Always check the value of the global 'errno' variable after using this
//   function to see if anything went wrong. (It will be zero if all is well.)
//
vector<TString> CreateCondorScriptME::read_directory( const string& path ) {
   vector <TString> result;
   dirent* de;
   DIR* dp;
   errno = 0;
   dp = opendir( path.empty() ? "." : path.c_str() );
   if (dp)
   {
      while (true)
      {
         errno = 0;
         //readdir() gives you the directory entries in the order you get them. 
         //If you want to sort them you must do that yourself once you get the names.
         de = readdir( dp );
         if (de == NULL) break;
         result.push_back( TString( de->d_name ) );
      }
      closedir( dp );
      std::sort( result.begin(), result.end() );
   }
   else
   {
      cout << "\tERROR::CreateCondorScriptME::read_directory could not open directory " << path << endl; 
   }
   return result;
}

//______________________________________________________________________________
int CreateCondorScriptME::get_file_number(TString filename) {
   char* str = strdup(filename.Data()), *p = str;
   vector<int> numbers;
   
   while (*p) { // While there are more characters to process...
      if (isdigit(*p)) { // Upon finding a digit, ...
         long val = strtol(p, &p, 10); // Read a number, ...
         numbers.push_back(val);
      }
      else { // Otherwise, move on to the next character.
         p++;
      }
   }

   free(str);

   if(numbers.size()>0)
      return numbers.back();
   else
      return -1;
}

//______________________________________________________________________________
void CreateCondorScriptME::checkForBadFiles(bool verbose) {
   cout << "Starting checkFiles..." << endl;

   string directory_to_check = completedOutputDir;

   if(!TString(directory_to_check).EndsWith("/")) directory_to_check += "/";
   cout << "\tChecking directory " << directory_to_check << endl;

   Int_t eil = gErrorIgnoreLevel;
   gErrorIgnoreLevel=kError;
   files = read_directory(directory_to_check);
   TFile* current_file_being_checked;
   TTree* METree;
   good_files.reserve(files.size());   
   bad_files.reserve(files.size());
   bad = false;
   counter = 0;

   cout << "\t" << files.size() << " files to check" << endl;
   for(unsigned int i=0; i<files.size(); i++) {
   //for(unsigned int i=11000; i<files.size(); i++) {
      bad = false;
      if(verbose)
         cout << "\tChecking file " << files[i] << " ... " << endl;
      else {
         if(i%500==0)
            cout << "\t" << i << " : Checking file " << files[i] << " ... " << endl;
      }

      void* dirp = gSystem->OpenDirectory(TString(directory_to_check)+files[i]);
      if(dirp) {
         if(verbose)
            cout << "\tWARNING::" << TString(directory_to_check)+files[i] << " is a directory and not a file." << endl
                 << "\tIt will be skipped." << endl;
         continue;
      }
      
      current_file_being_checked = new TFile(TString(directory_to_check)+files[i],"READ");
      
      if(!bad && current_file_being_checked->IsZombie()) {
         if(verbose)
            cout << "\t\tWARNING:: File " << current_file_being_checked->GetName() << " is corrupt." << endl;
   
         if(get_file_number(files[i])>=0) {
            bad_files.push_back(get_file_number(files[i]));
         }
         bad = true;
      }
      if(!bad && current_file_being_checked->TestBit(TFile::kRecovered)) {
         //the Recover procedure has been run when opening the file
         //and the Recover was successful
         if(verbose)
            cout << "\t\tWARNING::File " << current_file_being_checked->GetName() << " needed to have keys recovered." << endl;

         if(get_file_number(files[i])>=0) {
            bad_files.push_back(get_file_number(files[i]));
         }
         bad = true;
      }
      
      METree = (TTree*)current_file_being_checked->Get(tree_name);
      if(!bad && !METree) {
         if(verbose)
            cout << "\t\tWARNING::METree does not exist in " << current_file_being_checked->GetName() << endl;

         if(get_file_number(files[i])>=0) {
            bad_files.push_back(get_file_number(files[i]));
         }
         bad = true;
      }
      if(!bad && METree->GetEntries()!=nEventsPerJob) {
         if((nEventsPerJob>0)&&(nEventsTot>0)&&(get_file_number(files[i])>0)&&((METree->GetEntries())==(nEventsTot-(get_file_number(files[i])*nEventsPerJob)))&&(get_file_number(files[i])==nJobsTot-1)) {
            bad = false;
         }
         else {
            if(verbose)
               cout << "\t\tWARNING::METree contains " << METree->GetEntries() << "/" << nEventsPerJob << "entries." <<  endl
                    << "\t\t" << current_file_being_checked->GetName() << " is bad." << endl;
            
            if(get_file_number(files[i])>=0) {
               bad_files.push_back(get_file_number(files[i]));
            }
            bad = true;
         }
      }
      
      if(!bad) {
         if(verbose)
            cout << "\t\tFile " << current_file_being_checked->GetName() << " is okay." << endl;
         good_files.push_back(get_file_number(files[i]));
      }
      else {
         counter++;
      }
      
      current_file_being_checked->Close();
      delete current_file_being_checked;
      
      gSystem->FreeDirectory(dirp);
   }

   gErrorIgnoreLevel = eil;

   for(unsigned int i=0; i<bad_files.size(); i++) {
      cout << bad_files[i];
      if(i<bad_files.size()-1)
         cout << ",";
      if(i==bad_files.size()-1)
         cout << endl;
   }

   cout << counter << " bad files found." << endl;
   cout << bad_files.size() << " bad files printed." << endl;
   cout << "DONE" << endl;
}

//______________________________________________________________________________
void CreateCondorScriptME::checkForMissingFiles() {
   missing_files = vector<int> (nJobsTot+100);

   if(nJobsTot>0) {
      cout << "\nChecking for missing files..." << endl;
      all_files.insert(all_files.end(),good_files.begin(),good_files.end());
      all_files.insert(all_files.end(),bad_files.begin(),bad_files.end());
      sort(all_files.begin(),all_files.end());
      vector<int> temp;
      for(int i=0; i<nJobsTot-1; i++)
         temp.push_back(i);
      vector<int>::iterator it;
      it=set_difference(temp.begin(),temp.end(),all_files.begin(),all_files.end(),missing_files.begin());
      missing_files.resize(it-missing_files.begin());

      for(unsigned int i=0; i<missing_files.size(); i++) {
         cout <<  missing_files[i];
         if(i<missing_files.size()-1)
            cout << ",";
         if(i==missing_files.size()-1)
            cout << endl;
      }

      cout << missing_files.size() << " missing files." << endl;
      
      cout << "DONE" << endl;
   }
}

//______________________________________________________________________________
pair<TString, int> CreateCondorScriptME::getInvalidFiles() {
   int retInt = missing_files.size()+bad_files.size();
   TString retStr;
   vector<int> tempVec(missing_files);
   tempVec.insert(tempVec.end(),bad_files.begin(),bad_files.end());
   sort(tempVec.begin(),tempVec.end());
   for(unsigned int i=0; i<tempVec.size(); i++) {
      retStr +=  tempVec[i];
      if(i<tempVec.size()-1)
         retStr += ",";
   }

   return make_pair(retStr,retInt);

}

//______________________________________________________________________________
//// Creates a condor launcher (ScriptDir/CondorLauncher_ScriptNameSuffix) script and a .csh script used by it (ScriptDir/CondorScript_ScriptNameSuffix.csh).
//// These can be used from any directory, provided it has cteq5l.tbl & cteq6l.tbl libraries (normally located in TAMUWW/PDFs/Pdfdata/), a /log subdirectory as well as run_MatrixElement macro in it. Note: rootInputDir is the location of the .root input file as seen from this directory.
//// Running the CondorLauncher will generate root files of the form OutputName#.root (with #=FirstJob,...,FirstJob+NJobs-1).
//// To create a script which will run over disjoint jobs set createDisjointJobScript=false and input disjointJobsString= JobA , JobB , ... , JobZ. Note, the job numbers & corresponding events are the same as in the continuos string, you are simply specifying explicitly which jobs to run over (through a series of if statements).
void CreateCondorScriptME::writeCondorScripts() {
   ofstream outlauncher, outscript;
   TString tempStr;
   tempStr=ScriptNameSuffix;
   tempStr="CondorLauncher_" + tempStr;
   tempStr=ScriptDir + tempStr;
   outlauncher.open(tempStr,ios::out);
   tempStr=".csh";
   tempStr=ScriptNameSuffix + tempStr;
   tempStr="CondorScript_" + tempStr;
   tempStr=ScriptDir + tempStr;
   outscript.open(tempStr,ios::out);

   //Temporary Variables
//   char NEvtsPerJob_C[6], FirstJob_C[6], NJobs_C[6];
//   //sprintf(FirstEvt_C,"%i",FirstEvt);
//   sprintf(NEvtsPerJob_C,"%i",NEvtsPerJob);
//   sprintf(FirstJob_C,"%i",FirstJob);
//   sprintf(NJobs_C,"%i",NJobs);

   /// Make the launcher
   outlauncher << "CB = /uscms_data/d2/aperloff/YOURWORKINGAREA/MatrixElement/gitty/CMSSW_5_3_2_patch5/" << endl
               << "SD = $(CB)/src/" << endl
               << "TCD = $(SD)/TAMUWW/ConfigFiles" << endl
               << "LD = $(CB)/lib/slc5_amd64_gcc462/" << endl
               << "BD = $(CB)/bin/slc5_amd64_gcc462/" << endl << endl;
   outlauncher << "universe = vanilla" << endl;
   outlauncher << "Executable = " << tempStr << endl;
   //outlauncher << "Transfer_Input_Files = " << rootInputDir << rootInputName << ", run_MatrixElement" << globalRunMESuffix 
   //            << ", cteq5l.tbl, cteq6l.tbl, " << specificEventsFileLocation << "/micro" << ScriptNameSuffix << "missingEvents.txt" << endl;
   outlauncher << "transfer_input_files = " << "$(TCD), $(BD)/run_MatrixElement" << globalRunMESuffix 
               << ", cteq5l.tbl, cteq6l.tbl"
               <<  ", $(LD)/libTAMUWWMatrixElement.so, $(LD)/libTAMUWWTools.so, $(LD)/libTAMUWWSpecialTools.so, $(LD)/libTAMUWWMEPATNtuple.so, $(LD)/libElectroWeakAnalysisVPlusJets.so, $(LD)/libEGammaEGammaAnalysisTools.so, $(LD)/libCondFormatsEgammaObjects.so, $(LD)/libDataFormatsPatCandidates.so, $(LD)/libJetMETAnalysisJetUtilities.so, $(LD)/libTAMUWWAuxFunctions.so, $(LD)/libTAMUWWDHELAS.so, $(LD)/libTAMUWWPDFs.so, $(LD)/libSHyFTTemplateMakers.so, $(LD)/libTAMUWWIntegrator.so, $(LD)/libCondFormatsEgammaObjects.so, $(LD)/libTAMUWWCubaCuhre.so, $(LD)/libTAMUWWCubaDivonne.so, $(LD)/libTAMUWWCubaSuave.so, $(LD)/libTAMUWWCubaVegas.so" << endl;
   outlauncher << "Requirements = OpSys == \"LINUX\" && (Arch != \"DUMMY\" )" << endl;
   outlauncher << "request_disk = 1000000" << endl;
   outlauncher << "request_memory = 1945" << endl;
   outlauncher << "+BigMemoryJob = TRUE" << endl;
   outlauncher << "Should_Transfer_Files = ALWAYS" << endl;
   outlauncher << "WhenToTransferOutput = ON_EXIT_OR_EVICT" << endl;
   outlauncher << "x509userproxy = $ENV(X509_USER_PROXY)" << endl;

   outlauncher << "Output = log/" << ScriptNameSuffix << "/CondorME_" << ScriptNameSuffix << "_C$(Cluster)_$(Process).stdout" << endl;
   outlauncher << "Error = log/" << ScriptNameSuffix << "/CondorME_" << ScriptNameSuffix << "_C$(Cluster)_$(Process).stderr" << endl;
   outlauncher << "Log = log/" << ScriptNameSuffix << "/CondorME_" << ScriptNameSuffix << "_C$(Cluster)_$(Process).log" << endl;
   outlauncher << "notify_user = none" << endl;
   outlauncher << "Arguments = $(Process) " << nEventsPerJob <<  endl;
   outlauncher << "Queue " << nJobsTot << endl;

   /// Make the internal script
   outscript << "#! /bin/csh" << endl;
   outscript << "echo \"Starting\" " << endl;
   outscript << "setenv SCRAM_ARCH slc5_amd64_gcc462" << endl
             << "source /cvmfs/cms.cern.ch/cmsset_default.csh" << endl
             << "set CONDORDIR=$PWD" << endl << endl
             << "scram p CMSSW_5_3_2_patch5" << endl
             << "cd CMSSW_5_3_2_patch5" << endl
             << "mkdir -p src/TAMUWW" << endl
             << "mkdir -p lib/${SCRAM_ARCH}" << endl << endl
             << "mv ../ConfigFiles src/TAMUWW" << endl
             << "mv ../lib*.so lib/slc5_amd64_gcc462" << endl << endl;
   outscript << "scram b" << endl;
   outscript << "eval `scram runtime -csh` " << endl
             << "cd .." << endl << endl;
   outscript << "@ JobNumber = $argv[1]" << endl;
   outscript << "echo \"JobNumber=$JobNumber\"" << endl;
   outscript << "@ NEvtsPerJob = $argv[2]" << endl;
   outscript << "echo \"NEvtsPerJob=$NEvtsPerJob\"" << endl;
   outscript << "@ TheJob = " << firstJob << " + $JobNumber" << endl;

   if ( createDisjointJobScript ) {
      outscript << "###### START: Explicitly Setting Starting Job Number ######" << endl;
      int nJob=0;
      TString nj;
      TString temp="temp";
      istrstream djstr(disjointJobsString);
      std::string token;
      //while ( temp!="" ) {
      while (std::getline(djstr, token, ',')) {
         //djstr >> nj >> temp;
         outscript << "if ( $JobNumber == " << nJob << " ) then" << endl;
         //outscript << "   @ TheJob = " << nj << endl;
         outscript << "   @ TheJob = " << token << endl;
         outscript << "endif" << endl;
         nJob++;
      }

      if ( nJob!=nJobsTot ) { cerr << "Error: the job count in the disjointJobsString does not correspond to NJobs" << endl; }
      outscript << "###### FINISH: Explicitly Setting Starting Job Number ######" << endl;
   }
   //  outscript << "@ StartEvt = $argv[2] * $argv[1] + " << FirstEvt_C << endl;
   // outscript << "echo \"StartEvt=NEvts*Process+" << FirstEvt_C << "=$StartEvt\"" << endl;
   // outscript << "@ StartJob = $argv[1] + " << FirstJob_C << endl;

   outscript << "echo \"TheJob=$TheJob\"" << endl;
   outscript << "@ StartEvt = $TheJob * $NEvtsPerJob" << endl;
   outscript << "echo \"StartEvt=$StartEvt\"" << endl;


   //  outscript << "run_MatrixElement " << rootInputName << " " << OutputName << "$StartJob.root EvtTree " << NEvtsPerJob_C << " $StartEvt 1" << endl;

   outscript << "./run_MatrixElement " << rootInputDir << rootInputName << " " << OutputName << "$TheJob.root jets2p $NEvtsPerJob $StartEvt 1 1 PS "
             << useSpecificEvents << " "  << specificEventsFileLocation << "/micro" << ScriptNameSuffix << "missingEvents.txt" << endl;
   outscript << "echo \"Finished\"" << endl;

   tempStr="chmod +x " + tempStr;
   system(tempStr);
}

//______________________________________________________________________________
void CreateCondorScriptME::writePBSScripts(TString pbsQueue, TString pbsWalltimeStr, TString pbsUser) {
   ofstream outlauncher, outscript;
   TString tempStr;
   tempStr=".pbs";
   tempStr=ScriptNameSuffix+tempStr;
   tempStr="CondorLauncher_" + tempStr;
   tempStr=ScriptDir + ScriptNameSuffix + "/" + tempStr;

   DIR* dp;
   errno = 0;
   dp = opendir((ScriptDir+ScriptNameSuffix+"/").Data());
   if (dp) {
      closedir(dp);
   }
   else {
      mkdir((ScriptDir+ScriptNameSuffix+"/").Data(), S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH);
   }

   outlauncher.open(tempStr,ios::out);
   tempStr=".csh";
   tempStr=ScriptNameSuffix + tempStr;
   tempStr="CondorScript_" + tempStr;
   tempStr=ScriptDir + ScriptNameSuffix + "/" + tempStr;
   outscript.open(tempStr,ios::out);

   /// Make the launcher
   outlauncher << "#PBS -l nodes=1:ppn=1" << endl;
   outlauncher << "#PBS -l walltime=" << pbsWalltimeStr << endl;
   outlauncher << "#PBS -q " << pbsQueue << endl;
   outlauncher << "#PBS -j oe" << endl;
   if ( createDisjointJobScript ) {
      outlauncher << "#PBS -t " << disjointJobsString << endl;
   }
   else{
      outlauncher << "#PBS -t 0-" << nJobsTot << endl;
   }
///home/aperloff/MatrixElement/CMSSW_5_3_2_patch4/src/MECONDORJobs_Alexx/WJets_part2/CondorScript_WJets_part2.csh ${PBS_ARRAYID} 400 > /fdata/hepx/store/user/aperloff/ME532/log/CondorME_WJets_part2_${PBS_JOBID}.stdout
   outlauncher << endl << tempStr << " ${PBS_ARRAYID} " << nEventsPerJob << " > " << completedOutputDir << "/log/" << ScriptNameSuffix << "/CondorME_" 
               << ScriptNameSuffix << "_${PBS_JOBID}.stdout" << endl; // FIX ME 
   outlauncher << endl << "exit 0" << endl;


   /// Make the internal script
   outscript << "#! /bin/csh" << endl;
   outscript << "echo \"Starting\" " << endl;
   outscript << "cd $HOME/MatrixElement/gitty/CMSSW_5_3_2_patch4/src" << endl;//***SET DEPENDENT ON THE CMSSW RELEASE***
   outscript << "setenv SCRAM_ARCH slc5_amd64_gcc462" << endl;
   outscript << "pwd" << endl;
   outscript << "eval `scram runtime -csh` " << endl;
   //outscript << "cd MECONDORJobs_" << pbsUser << endl;// FIX ME !!!
   outscript << "cd /fdata/hepx/store/user/aperloff/ME532_Completion/" << endl;
   outscript << "@ JobNumber = $argv[1]" << endl;
   outscript << "echo \"JobNumber=$JobNumber\"" << endl;
   outscript << "@ NEvtsPerJob = $argv[2]" << endl;
   outscript << "echo \"NEvtsPerJob=$NEvtsPerJob\"" << endl;
   outscript << "@ TheJob = " << firstJob << " + $JobNumber" << endl;
   outscript << "echo \"TheJob=$TheJob\"" << endl;
   outscript << "@ StartEvt = $TheJob * $NEvtsPerJob" << endl;
   outscript << "echo \"StartEvt=$StartEvt\"" << endl;
   //outscript << "./run_MatrixElement " << completedOutputDir << "/" << ScriptNameSuffix << "/" << rootInputName << " " << completedOutputDir << "/" 
   outscript << "./run_MatrixElement " << rootInputDir << "/" << rootInputName << " " << completedOutputDir << "/" 
             << ScriptNameSuffix << "/" << OutputName << "$TheJob.root jets2p $NEvtsPerJob $StartEvt 1 1 PS " << useSpecificEvents << " "  
             << specificEventsFileLocation << "/micro" << ScriptNameSuffix << "_missingEvents.txt" << endl;
   outscript << "echo \"Finished\"" << endl;

   tempStr="chmod +x " + tempStr;
   system(tempStr);
}
