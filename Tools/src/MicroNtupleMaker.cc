// This class libraries
#include "TAMUWW/Tools/interface/MicroNtupleMaker.hh"

//______________________________________________________________________________
MicroNtupleMaker::MicroNtupleMaker() {

}//C'tor

//______________________________________________________________________________
MicroNtupleMaker::~MicroNtupleMaker() {

   if(mergeChain)
      delete mergeChain;
   if(index)
      delete index;
   if(mergeEventNtuple)
      delete mergeEventNtuple;

}//D'tor

//______________________________________________________________________________
void MicroNtupleMaker::createMicroNtuple(TString inputPath, TString outputPath, TString outputSuffix, int largeProcessCase, TString smallProcessLabel){
//// Large processes need to be split up into parts 'by hand'
//// largeProcessCase=0 corresponds to small processes
//// The files with similar names, but which need to be ran on separately should be put into separate directories (e.g. STopT_T and STopT_Tbar)

   //string basePath="/uscms_data/d3/ilyao/Winter12to13ME8TeV/MEResults/";

   //Input files for 2-jet bin
//   string inputPath=basePath+"rootOutput/";
//   string inputPathBeta="/eos/uscms/store/user/eusebi/Winter12to13ME8TeV/rootOutput/";
   unsigned nJet = 2;

   //Input files for 3-jet bin
   //string inputPath=basePath+"Event_Probs_Summer09_3jets_StndTF/";
   //unsigned nJet = 3;

   //Output files
   //  string outputPath=basePath+"microNtuples/";

   //Create the list of MicroNtuples
   vector<MyStr> listOfMicroNtuples;

   /// Test Ntuples:
//   inputPath=basePath+"test_rootOutput/";
//   outputPath=basePath+"test_microNtuples/";
//   listOfMicroNtuples.push_back(MyStr("WZ*","testmicroWZ",false,false,false));
   TString microOutputSuffix="_BaseCuts"; //Set a suffix entire microNtuple production
   TString microOutputName;

   switch (largeProcessCase) {
      case 0:
      {
         ///Small Process:
         TString inputName="*";
         inputName=smallProcessLabel+inputName;
         microOutputName=microOutputSuffix;
         microOutputName=smallProcessLabel+microOutputName;
         microOutputName="micro"+microOutputName;
         bool addFilesInBaseDir;
         if(useXROOTD)
            addFilesInBaseDir = false;
         else
            addFilesInBaseDir = true;
         MyStr fileLocations(inputName,microOutputName,false,false,false,addFilesInBaseDir);
         
         if(useXROOTD) {
            fileLocations.AddDir("/");
         }
         else {
            for(unsigned int idir=0; idir<addDir.size(); idir++) {
               if(addDir[idir].CompareTo("/")==0)
                  continue;
               fileLocations.AddDir((addDir[idir]+inputName).Data());
            }
         }
         listOfMicroNtuples.push_back(fileLocations);
      }
      break;
      case -1:
         ///Multiple Small Processes
         break;
         ///W+3Jets
      case 1001:
      {
         microOutputName="microW3Jets_pt1"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("W3Jets1*",microOutputName,false,false,false));
      }
      break;
      case 1002:
      {
         microOutputName="microW3Jets_pt2"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("W3Jets2*",microOutputName,false,false,false));
      }
      break;
      case 1003:
      {
         microOutputName="microW3Jets_pt3"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("W3Jets3*",microOutputName,false,false,false));
      }
      break;
      case 1004:
      {
         microOutputName="microW3Jets_pt4"+microOutputSuffix;
         MyStr fileLocations_W3Jets_pt4("W3Jets4*",microOutputName,false,false,false);
         fileLocations_W3Jets_pt4.AddDir("W3Jets5*");
         fileLocations_W3Jets_pt4.AddDir("W3Jets6*");
         fileLocations_W3Jets_pt4.AddDir("W3Jets7*");
         fileLocations_W3Jets_pt4.AddDir("W3Jets8*");
         fileLocations_W3Jets_pt4.AddDir("W3Jets9*");
         fileLocations_W3Jets_pt4.AddDir("W3Jets0*");
         listOfMicroNtuples.push_back(fileLocations_W3Jets_pt4);
      }
      break;


      ///Data: Muon
      case 1011:
      {
         microOutputName="microSingleMu_Data_19p279fb_pt1"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("SingleMu_Data_19p279fb10*",microOutputName,false,false,false));
      }
      break;
      case 1012:
      {
         microOutputName="microSingleMu_Data_19p279fb_pt2"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("SingleMu_Data_19p279fb2*",microOutputName,false,false,false));
      }
      break;
      case 1013:
      {
         microOutputName="microSingleMu_Data_19p279fb_pt3"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("SingleMu_Data_19p279fb3*",microOutputName,false,false,false));
      }
      break;
      case 1014:
      {
         microOutputName="microSingleMu_Data_19p279fb_pt4"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("SingleMu_Data_19p279fb4*",microOutputName,false,false,false));
      }
      break;
      case 1015:
      {
         microOutputName="microSingleMu_Data_19p279fb_pt5"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("SingleMu_Data_19p279fb5*",microOutputName,false,false,false));
      }
      break;
      case 1016:
      {
         microOutputName="microSingleMu_Data_19p279fb_pt6"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("SingleMu_Data_19p279fb6*",microOutputName,false,false,false));
      }
      break;
      case 1017:
      {
         microOutputName="microSingleMu_Data_19p279fb_pt7"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("SingleMu_Data_19p279fb7*",microOutputName,false,false,false));
      }
      break;
      case 1018:
      {
         microOutputName="microSingleMu_Data_19p279fb_pt8"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("SingleMu_Data_19p279fb8*",microOutputName,false,false,false));
      }
      break;
      case 1019:
      {
         microOutputName="microSingleMu_Data_19p279fb_pt9"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("SingleMu_Data_19p279fb9*",microOutputName,false,false,false));
      }
      break;
      case 1020:
      {
         microOutputName="microSingleMu_Data_19p279fb_pt10"+microOutputSuffix;
         listOfMicroNtuples.push_back(MyStr("SingleMu_Data_19p279fb11*",microOutputName,false,false,false));
      }
      break;
      case 1021:
      {
         microOutputName="microSingleMu_Data_19p279fb_pt11"+microOutputSuffix;
         MyStr fileLocations_SingleMu_Data_19p279fb_pt11("SingleMu_Data_19p279fb12*",microOutputName,false,false,false);
         fileLocations_SingleMu_Data_19p279fb_pt11.AddDir("SingleMu_Data_19p279fb13*");
         fileLocations_SingleMu_Data_19p279fb_pt11.AddDir("SingleMu_Data_19p279fb14*");
         fileLocations_SingleMu_Data_19p279fb_pt11.AddDir("SingleMu_Data_19p279fb15*");
         fileLocations_SingleMu_Data_19p279fb_pt11.AddDir("SingleMu_Data_19p279fb16*");
         fileLocations_SingleMu_Data_19p279fb_pt11.AddDir("SingleMu_Data_19p279fb17*");
         fileLocations_SingleMu_Data_19p279fb_pt11.AddDir("SingleMu_Data_19p279fb18*");
         fileLocations_SingleMu_Data_19p279fb_pt11.AddDir("SingleMu_Data_19p279fb19*");
         fileLocations_SingleMu_Data_19p279fb_pt11.AddDir("SingleMu_Data_19p279fb0*");
         listOfMicroNtuples.push_back(fileLocations_SingleMu_Data_19p279fb_pt11);
      }
      break;
   }

   //Loop over all of them
   for ( vector<MyStr>::iterator it = listOfMicroNtuples.begin();
         it != listOfMicroNtuples.end(); it++){
    
      makeMicroNtuple(it->GetLocationVector(inputPath), outputPath+ it->name + outputSuffix + ".root", nJet, it->mistag, it->nonw, it->untag);
    
   }//for


}//createAllMicroNtuples()

//______________________________________________________________________________
void MicroNtupleMaker::makeMicroNtuple(TString location, TString output, unsigned nJets, bool doLight, bool doNonW, bool doUntag){

   vector<TString> locations;
   locations.push_back(location);
   makeMicroNtuple(locations, output, nJets, doLight, doNonW, doUntag);
  
}

//______________________________________________________________________________
void MicroNtupleMaker::makeMicroNtuple(vector<TString> locations, TString output, unsigned nJets, bool doLight, bool doNonW, bool doUntag){

   cout<<endl;
   cout<<"  Generating MicroNtuple "<<output<<endl;

   TChain chain("METree");
   int file_count =0;
   for (unsigned e=0;e<locations.size(); e++) {
      if(useXROOTD) {
         //cout << "Using Xrootd to open files from " << locations[e] << endl;
         ////Get list of files at location
         //TSystemDirectory dir(locations[e], locations[e]);
         //TString pwd(gSystem->pwd());
         //TList *files = dir.GetListOfFiles();
         //gSystem->cd(pwd.Data());
         //TFileCollection fc("fc","fc");
         //if (files) {
         //   cout << "Got list of files from " << locations[e] << endl;
         //   //files->Print();
         //   TSystemFile *file;
         //   TString fname;
         //   TIter next(files);
         //   while ((file=(TSystemFile*)next())) {
         //      fname = file->GetName();
         //      if (!file->IsDirectory() && fname.EndsWith(".root")) {
         //         //create file info object
         //         cout << "filename (Hash):" << locations[e].ReplaceAll("/eos/uscms","root://cmsxrootd.fnal.gov//") + fname.Data()
         //              << " (" << TString(locations[e].ReplaceAll("/eos/uscms","root://cmsxrootd.fnal.gov//") + fname.Data()).Hash()
         //              <<")"<< endl;
         //         TFileInfo* file_info = new TFileInfo(locations[e].ReplaceAll("/eos/uscms","root://cmsxrootd.fnal.gov//") + fname.Data());
         //         //cout << fname.Data() << endl;
         //         //Add file info object to file collection object
         //         //usleep(500);
         //         //file_count += fc.Add(&file_info);
         //         file_count += fc.Add(file_info);
         //         delete file_info;
         //         cout << "file_count = " << file_count << endl;
         //      }
         //   }
         //}
         cout << "Using Xrootd to open files from " << locations[e] << endl;
         TFileCollection fc("fc","fc",locations[e]+"/"+currentProcess+"_fileList.txt");
         file_count = fc.GetNFiles();
         chain.AddFileInfoList(fc.GetList());
      }
      else {
         cout<<"\tAdding "<<(locations[e] + "*.root")<<endl;   
         file_count += chain.Add((locations[e] + "*.root"));
      }
   }

   if (file_count==0){
      cout << "\tNo files found!  Aborting.\n";
      return;
   }

   if (mergeNewEventNtuple.CompareTo("")!=0) {
      mergeChain = new TChain("jets2p","Output tree for matrix element");
      if(currentProcess=="ZJets") currentProcess="DYJets";
      cout << "\tAdding " << mergeNewEventNtuple << currentProcess << ".root/PS/jets2p" << endl;
      mergeChain->Add(mergeNewEventNtuple+currentProcess+".root/PS/jets2p");
      mergeChain->LoadTree(0);
      mergeTree = mergeChain->GetTree();
      mergeTree->SetCacheSize(10000000);
      mergeTree->AddBranchToCache("*");

   }

   chain.SetBranchStatus("*",0);
   chain.SetBranchStatus("event*",1);
   chain.LoadTree(0);
   Int_t treeEntries = chain.GetTree()->GetEntries(); treeEntries=treeEntries;
   Long64_t* offsets = chain.GetTreeOffset();
   TTree* tmpTree;
   index = new TTreeIndex();
   int sumIndex = 0, sumTreeEntries = 0;
   if(mergeNewEventNtuple.CompareTo("")!=0) {
      for(Int_t t=0; t<chain.GetNtrees(); t++) {
      //for(Int_t t=0; t<2; t++) {
         //chain.LoadTree((t*treeEntries)+1);
         chain.LoadTree(offsets[t]);
         cout << "MicroNtupleMaker::makeMicroNtuple Building index for tree " << chain.GetTreeNumber() << " of " << chain.GetNtrees()-1 << " ... ";
         tmpTree = chain.GetTree();
         //tmpTree->SetBranchStatus("*",0);
         //tmpTree->SetBranchStatus("m_event",1);
         tmpTree->BuildIndex("event");
         index->Append((TTreeIndex*)tmpTree->GetTreeIndex());
         cout << "DONE" << endl;
         sumIndex += ((TTreeIndex*)tmpTree->GetTreeIndex())->GetN();
         sumTreeEntries += tmpTree->GetEntries();
         if(debug) {
            cout << "\tAdded " << ((TTreeIndex*)tmpTree->GetTreeIndex())->GetN() << " entries to make " << index->GetN() << " total (" << sumTreeEntries << " summed and " << flush;
            cout << chain.GetEntries() << " in chain)" << endl;
         }
         else
            cout << "\tAdded " << ((TTreeIndex*)tmpTree->GetTreeIndex())->GetN() << " entries to make " << index->GetN() << " total (" << sumTreeEntries << " in chain)" << endl;
         if (((TTreeIndex*)tmpTree->GetTreeIndex())->GetN()!=tmpTree->GetEntries()) {
            cout << "\t\tERROR::makeMicroNtuple Added TTreeIndex has a different number of entries ("
                 << ((TTreeIndex*)tmpTree->GetTreeIndex())->GetN() << ") than the TTree is comes from ("
                 << tmpTree->GetEntries() << ")" << endl;
            assert(((TTreeIndex*)tmpTree->GetTreeIndex())->GetN()==tmpTree->GetEntries());
         }
      }
   }
   if(index->GetN()>chain.GetEntries()) {
      cout << "ERROR::makeMicroNtuple The TTreeIndex for the chain of METree files has more entries ("
           << index->GetN() << ") than the chain itself (" << chain.GetEntries() << ")" << endl;
      cout << "Compare this to the summed values of " << sumIndex << " and " << sumTreeEntries << " respectively" << endl;
      assert(index->GetN()==chain.GetEntries());
   }
   chain.SetBranchStatus("*",1);
   chain.SetCacheSize(10000000);
   chain.AddBranchToCache("*");

   cout << endl << endl << "makeMicroNtuple(TChain& chain, TString output,"
        << " unsigned nJets, bool doLight, bool doNonW, bool doUntag) ... " << endl;

   makeMicroNtuple(chain, output, nJets, doLight, doNonW, doUntag);

}

//______________________________________________________________________________
//void MicroNtupleMaker::makeMicroNtuple(TChain & chain, TString output, unsigned nJets, 
void MicroNtupleMaker::makeMicroNtuple(TChain& chain, TString output, unsigned nJets, 
                                       bool doLight, bool doNonW, bool doUntag)
{
   cout << "MicroNtupleMaker::makeMicroNtuple Make ntuples ... ";
   // Create the base objects
   EventNtuple * eventNtuple      = new EventNtuple();
   METree      * meNtuple         = new METree();
   METree      * meNtupleBlank    = new METree();
   MicroNtuple * microNtuple      = new MicroNtuple(nJets);
   MicroNtuple * microNtupleBlank = new MicroNtuple(nJets);
   cout << "DONE" << endl;

   cout << "MicroNtupleMaker::makeMicroNtuple Set METree and EvtTree branches ... ";
   // and set the meNtuple to read from it
   chain.SetBranchAddress("METree", &meNtuple);
   chain.SetBranchAddress("EvtTree", &eventNtuple);
   cout << "DONE" << endl;

   //sort chains by event number
   if(mergeNewEventNtuple.CompareTo("")!=0) {
      
      cout << "MicroNtupleMaker::makeMicroNtuple Set merge EvtNtuple branch ... ";
      mergeEventNtuple = new EventNtuple();
      mergeChain->GetFile()->cd();
      mergeTree->SetBranchAddress("EvtNtuple", &mergeEventNtuple);
      cout << "DONE" << endl;

      Long64_t local = -1;
      chain.SetBranchStatus("*",0);
      chain.SetBranchStatus("event*",1);
      int nIndex = index->GetN();      
      cout << "MicroNtupleMaker::makeMicroNtuple Making event index map for " << nIndex << " entries ... " << endl;
      cout.flush();
      for( int i = 0; i < nIndex ; ++i ) {
         loadbar2(i+1,nIndex);
         local = index->GetIndex()[i];
         //chain.GetEntry(local,0);
         //eventIndex[eventNtuple->event] = local;
         //Warn if duplicate event index
         if(eventIndex.find(index->GetIndexValues()[i]>>31)!=eventIndex.end()) {
            duplicateIndex[local] = index->GetIndexValues()[i]>>31;
            if(debug) cout << "WARNING::makeMicroNtuple::There could be duplicate events in the chain" << endl;
         }
         if(local>chain.GetEntries()) {
            cout << "WARNING::makeMicroNtuple::The local index is greater than the number of entries in the chain" << endl;
         }
         eventIndex[index->GetIndexValues()[i]>>31] = local;
      }
      chain.SetBranchStatus("*",1);
      cout << endl <<"DONE" << endl;
   }

   cout << "MicroNtupleMaker::makeMicroNtuple Create output file and setup output tree ... ";
   // Create and output file and clone the tree that will be in the output and set microNtuple that fills it
   TFile outputFile(output, "RECREATE");
   TTree* outputTree = new TTree(chain.GetName(),chain.GetTitle());
   outputTree->Branch("METree", "METree", &meNtuple);
   if(mergeNewEventNtuple.CompareTo("")!=0)
      outputTree->Branch("EvtTree", "EventNtuple", &mergeEventNtuple);
   else
      outputTree->Branch("EvtTree", "EventNtuple", &eventNtuple);
   outputTree->Branch("mnt", "MicroNtuple", &microNtuple);
   cout << "DONE" << endl << endl << endl;

   // Holds the list of events to avoid duplication
   set<Signature> sigSet;

   // The index map has not been filled yet
   imFilled = false;

   //
   // If the paths for the BDT weights are set, store them in the tree
   //
   if (fillBDT) {
      setBDTReadersFromTable();
   }

   // Get the entries and report if zero
   mergeNewEventNtuple.CompareTo("")!=0 ? nentries = static_cast<unsigned>(mergeTree->GetEntries()) : nentries = static_cast<unsigned>(chain.GetEntries());
   cout << endl << endl << "Original chain has " << nentries << " entries" << endl;
   if (nentries == 0){
      cout << "\t\tNo entries found!  Aborting.\n";
      return;
   }

   // If the start and end entries are not 0 and -1 respectively then only loop through the requested entries
   startEndEntries.first!=0 ? startEntry = startEndEntries.first : startEntry = 0;
   startEndEntries.second!=-1 ? endEntry = startEndEntries.second : endEntry = nentries;
   // Check if the user made a mistake and the endEntry requested is greater than the last entry in the TChain or TTree
   endEntry>nentries ? endEntry = nentries : endEntry = endEntry;
   cout << "MicroNtupleMaker::makeMicroNtuple Starting with entry " << startEntry << " and ending with entry " << endEntry << endl;

   //Loop over all the entries in the original chain. The idea is to copy everything to microNtuple
   cout << "MicroNtupleMaker::makeMicroNtuple Filling MicroNtuple ..." << endl;
   if(missingEventsFlag){
      mergeTree->SetBranchStatus("*",0);
      mergeTree->SetBranchStatus("event",1);
   }
   for (unsigned ientry = startEntry; ientry < endEntry; ++ientry){

      loadbar2(ientry+1,nentries);

      // get the entry
      if (mergeNewEventNtuple.CompareTo("")!=0) {
         mergeTree->GetEntry(ientry);
         if(eventIndex.find(mergeEventNtuple->event)!=eventIndex.end()) {
            if(!missingEventsFlag)
               chain.GetEntry(eventIndex[mergeEventNtuple->event]);
            else
               continue;
         }
         else {
            missingME[ientry] = mergeEventNtuple->event;
            if(!missingEventsFlag) {
               microNtuple->clear();
               meNtuple->clear();
               outputTree->Fill();
            }
            continue;
         }
      }
      else {
         chain.GetEntry(ientry);
         //Check and remove duplicates here
         if (addDir.size()!=0 && runEventSet.alreadySeen(eventNtuple->run, eventNtuple->event)) {
            //cout << "WARNING, event Run = " << eventNtuple->run
            //     << ", Event = " << eventNtuple->event
            //     <<" is duplicated" << endl
            //     << "We will skip this event." << endl;
            duplicateIndex[ientry] = eventNtuple->event;
            if(!duplicateEventsFlag)
               continue;
         }
      }

      if (!imFilled) {
         microNtuple->indexMap = meNtuple->fillIndexMap();
         imFilled = true;
      }

      // clear the microNtuple
      microNtuple->clear();

      // First copy all the event probabilities
      if(meNtuple->getNProbStat()<1) {
         cout << endl << endl << "ientry = " << ientry << "\tThis is eventNtuple->event " << eventNtuple->event
              << ", mergeEventNtuple->event " << mergeEventNtuple->event 
              << ", meNtuple->getEvent() " << meNtuple->getEvent() << endl << endl;
         cout << endl << endl<< "There is a problem with getNProbStat (size = " << meNtuple->getNProbStat()
                                          << ") in entry " << meNtuple->getEvent() << " or " << eventNtuple->event << endl << endl; 
         cout << endl << endl << "eventIndex[mergeEventNtuple->event] = " << eventIndex[mergeEventNtuple->event] << endl << endl;
         cout << endl << endl << chain.GetFile()->GetName() << endl << endl;
         cout << endl << endl << endl;
      }
      bool passCopy = copyEventProbs(mergeEventNtuple,eventNtuple,meNtuple,microNtuple);
      if(!passCopy) {
         microNtuple->clear();
         continue;
      }
      setSizeRunEvent(eventNtuple,meNtuple,microNtuple);
      setEPDs(meNtuple, microNtuple);

      // TMVA variables
      if(mergeNewEventNtuple.CompareTo("")!=0) {
         setTMVAInformation(mergeEventNtuple,microNtuple);
      }
      else {
         setTMVAInformation(eventNtuple,microNtuple);
      }
      microNtuple->reader = 0;

      //
      // Finalize the output tree
      //

      outputTree->Fill();

      microNtuple->clear();
      meNtuple->clear();
   }//for entries

   //
   // Report some results
   //
   cout << endl << endl << "Wrote " << output << " with " << outputTree->GetEntries() << " entries" << endl;
   outputFile.Write();
  
   delete meNtuple;
   delete microNtuple;
   delete eventNtuple;    
   delete meNtupleBlank;
   delete microNtupleBlank;

   if(duplicateIndex.size()>0) {
      cout << endl << endl << "Duplicate entries in the event index ..." << endl;
      stringstream duplicateIndexStream;
      stringstream duplicateEventStream;
      map<int,int>::iterator finalIter = duplicateIndex.end();
      if(duplicateIndex.size()>0)
         --finalIter;
      for (map<int,int>::iterator it=duplicateIndex.begin(); it!=duplicateIndex.end(); it++) {
         if(it!=finalIter) {
            duplicateIndexStream << it->first << ",";
            duplicateEventStream << it->second << ",";
         }
         else {
            duplicateIndexStream << it->first << endl;
            duplicateEventStream << it->second << endl;
         }
      }
      cout << "The indecies for the duplicate events in the chain (" << duplicateIndex.size() << " total):" << endl << endl 
           << duplicateIndexStream.str() << endl
           << "The event numbers for the duplicate events in the chain (" << duplicateIndex.size() << " total):" << endl << endl
           << duplicateEventStream.str();
   }

   if(mergeNewEventNtuple.CompareTo("")!=0 && !missingEventsFlag) {
      cout << endl << endl <<"Missing Information ..." << endl;
      stringstream missingIndex;
      stringstream missingEvents;
      map<int,int>::iterator finalIter = missingME.end();
      if(missingME.size()>0)
         --finalIter;
      for (map<int,int>::iterator it=missingME.begin(); it!=missingME.end(); ++it) {
         if(it!=finalIter) {
            missingIndex << it->first << ",";
            missingEvents << it->second << ",";
         }
         else {
            missingIndex << it->first << endl;
            missingEvents << it->second << endl;
         }
      }
      cout << "The indecies for the entries missing ME (" << missingME.size() << " total):" << endl << endl << missingIndex.str() << endl
           << "The event numbers for the entries missing ME (" << missingME.size() << " total):" << endl << endl << missingEvents.str();
   }
   else if(mergeNewEventNtuple.CompareTo("")!=0 && missingEventsFlag) {
      writeMissingEventsFile(missingME);
   }

}// makeMicroNtuple

/*
void MicroNtupleMaker::writeMissingEventsFile(map<int,int> &missingME) {
   string filename = string(outputPath+"micro"+currentProcess+"_missingEvents.txt");
   cout << "\t\tWriting missing events file to " << filename << " ... ";
   fstream ofile;
   ofile.open(filename,std::fstream::out | std::fstream::trunc);
   int len = getIntLength(missingME.size());

   if (ofile.is_open()) {
      ofile << setw(34) << "Missing Index In New EventNtuple (" << setw(len) << missingME.size() << setw(7) << " total)" << setw(10) << " " << setw(23) << "Missing Event Numbers (" <<  setw(len) << missingME.size() << setw(7) << " total)" << endl;
      ofile << setfill('=') << setw(34+len+7) << "" << setfill(' ') << setw(10) << " " << setfill('=') << setw(23+len+7) << ""  << setfill(' ') << endl;

      for (map<int,int>::iterator it=missingME.begin(); it!=missingME.end(); ++it) {
         ofile << setw(34+len+7) << it->first << setw(10) << " " << setw(23+len+7) << it->second << endl; 
      }

      ofile.close();
      cout << "DONE" << endl;
   }
   else {
      cout << endl << "\t\tERROR::Could not open file " << filename << endl;
   }
}//writeMissingEventsFile
*/
void MicroNtupleMaker::writeMissingEventsFile(map<int,int> &missingME) {
   string filename = string(outputPath+"micro"+currentProcess+"_missingEvents.txt");
   cout << "\t\tWriting missing events file to " << filename << " ... ";

   Table* table = new Table("MissingEvents");
   TableRow* tableRow;
   TableCellInt* tableCellInt;
   //int value;
   int counter = 0;

   //for (Int_t index=0; index<missingME.size();index++) {
   for (map<int,int>::iterator it=missingME.begin(); it!=missingME.end(); ++it) {
      counter++;

      std::string s;
      std::stringstream ss;
      ss << counter;
      s = ss.str();
      ss.str("");
      ss << missingME.size();

      tableRow = new TableRow(s);

      for (Int_t col=0; col<2; col++) {
         string title;
         if (col==0)
            title = string("Missing Index In New EventNtuple (" + ss.str() + " total)");
         else
            title = string("Missing Event Numbers (" + ss.str() + " total)");

         tableCellInt = new TableCellInt(title);

         if (col==0)
            tableCellInt->val = it->first;
         else
            tableCellInt->val = it->second;

         tableRow->addCellEntries(tableCellInt);
      }
      table->addRow(*tableRow);
      delete tableRow;
   }

   table->printToFile(filename);

   cout << "DONE" << endl;

}//writeMissingEventsFile

int MicroNtupleMaker::getIntLength(int i) {
   int power = -1;
   int result = -1;
   while (result!=0) {
      power++;
      result = i/int(pow(10,power));
   }
   return power;
}//getIntLength

void MicroNtupleMaker::setBDTReadersFromTable() {
   map<TString,MVAVar> empty;
   for(int i=0; i<9; i++) {
      BDTReaders.push_back(make_pair(empty,new TMVA::Reader( "!Color:!Silent" )));
      MVAVar::setVarMap(BDTReaders.back().first);
   }
   if(debug) cout << "Setting BDT variables for " << DefaultValues::getBDTLocation(DEFS::jets2,DEFS::eq0tag,DEFS::UVa,true) << " ... " << endl;
   MVAVar::setUseFromTable(BDTReaders[0].first, DefaultValues::getBDTLocation(DEFS::jets2,DEFS::eq0tag,DEFS::UVa,true));
   if(debug) cout << "Setting BDT variables for " << DefaultValues::getBDTLocation(DEFS::jets3,DEFS::eq0tag,DEFS::UVa,true) << " ... " << endl;
   MVAVar::setUseFromTable(BDTReaders[1].first, DefaultValues::getBDTLocation(DEFS::jets3,DEFS::eq0tag,DEFS::UVa,true));
   if(debug) cout << "Setting BDT variables for " << DefaultValues::getBDTLocation(DEFS::jets4,DEFS::eq0tag,DEFS::UVa,true) << " ... " << endl;
   MVAVar::setUseFromTable(BDTReaders[2].first, DefaultValues::getBDTLocation(DEFS::jets4,DEFS::eq0tag,DEFS::UVa,true));
   if(debug) cout << "Setting BDT variables for " << DefaultValues::getBDTLocation(DEFS::jets2,DEFS::eq0tag,DEFS::TAMU,true) << " ... " << endl;
   MVAVar::setUseFromTable(BDTReaders[3].first, DefaultValues::getBDTLocation(DEFS::jets2,DEFS::eq0tag,DEFS::TAMU,true));
   if(debug) cout << "Setting BDT variables for " << DefaultValues::getBDTLocation(DEFS::jets3,DEFS::eq0tag,DEFS::TAMU,true) << " ... " << endl;
   MVAVar::setUseFromTable(BDTReaders[4].first, DefaultValues::getBDTLocation(DEFS::jets3,DEFS::eq0tag,DEFS::TAMU,true));
   if(debug) cout << "Setting BDT variables for " << DefaultValues::getBDTLocation(DEFS::jets4,DEFS::eq0tag,DEFS::TAMU,true) << " ... " << endl;
   MVAVar::setUseFromTable(BDTReaders[5].first, DefaultValues::getBDTLocation(DEFS::jets4,DEFS::eq0tag,DEFS::TAMU,true));
   if(debug) cout << "Setting BDT variables for " << DefaultValues::getBDTLocation(DEFS::jets2,DEFS::eq0tag,DEFS::Combined,true) << " ... " << endl;
   MVAVar::setUseFromTable(BDTReaders[6].first, DefaultValues::getBDTLocation(DEFS::jets2,DEFS::eq0tag,DEFS::Combined,true));
   if(debug) cout << "Setting BDT variables for " << DefaultValues::getBDTLocation(DEFS::jets3,DEFS::eq0tag,DEFS::Combined,true) << " ... " << endl;
   MVAVar::setUseFromTable(BDTReaders[7].first, DefaultValues::getBDTLocation(DEFS::jets3,DEFS::eq0tag,DEFS::Combined,true));
   if(debug) cout << "Setting BDT variables for " << DefaultValues::getBDTLocation(DEFS::jets4,DEFS::eq0tag,DEFS::Combined,true) << " ... " << endl;
   MVAVar::setUseFromTable(BDTReaders[8].first, DefaultValues::getBDTLocation(DEFS::jets4,DEFS::eq0tag,DEFS::Combined,true));
   setAllBDTReaders();
}

void MicroNtupleMaker::setAllBDTReaders() {
   vector<TString> weightfiles;
   weightfiles.push_back(DefaultValues::getBDTLocation(DEFS::jets2,DEFS::eq0tag,DEFS::UVa));
   weightfiles.push_back(DefaultValues::getBDTLocation(DEFS::jets3,DEFS::eq0tag,DEFS::UVa));
   weightfiles.push_back(DefaultValues::getBDTLocation(DEFS::jets4,DEFS::eq0tag,DEFS::UVa));
   weightfiles.push_back(DefaultValues::getBDTLocation(DEFS::jets2,DEFS::eq0tag,DEFS::TAMU));
   weightfiles.push_back(DefaultValues::getBDTLocation(DEFS::jets3,DEFS::eq0tag,DEFS::TAMU));
   weightfiles.push_back(DefaultValues::getBDTLocation(DEFS::jets4,DEFS::eq0tag,DEFS::TAMU));
   weightfiles.push_back(DefaultValues::getBDTLocation(DEFS::jets2,DEFS::eq0tag,DEFS::Combined));
   weightfiles.push_back(DefaultValues::getBDTLocation(DEFS::jets3,DEFS::eq0tag,DEFS::Combined));
   weightfiles.push_back(DefaultValues::getBDTLocation(DEFS::jets4,DEFS::eq0tag,DEFS::Combined));
   for(unsigned int i=0; i<BDTReaders.size(); i++) {
      cout << "Doing weightfile " << weightfiles[i] << " ... " << endl;
      TString var;
      int size = 0;
      MVAVarMap::iterator firstUsed;
      MVAVarMap::iterator firstUsedSpectator;
      //for(firstUsed = BDTReaders[i].first.begin(); firstUsed!=BDTReaders[i].first.end(); ++firstUsed) {
      //   if(firstUsed->second.use && !firstUsed->second.isSpectator) break;
      //}
      //for(firstUsedSpectator = BDTReaders[i].first.begin(); firstUsedSpectator!=BDTReaders[i].first.end(); ++firstUsedSpectator) {
      //   if(firstUsedSpectator->second.use && firstUsedSpectator->second.isSpectator) break;
      //}
      firstUsed = MVAVar::getFirstUsed(BDTReaders[i].first,false);
      firstUsedSpectator = MVAVar::getFirstUsed(BDTReaders[i].first,true);

      while(size!=firstUsed->second.maxIndex) {
         for(map<TString,MVAVar>::iterator it=BDTReaders[i].first.begin(); it!=BDTReaders[i].first.end(); ++it) {
            if(it->second.use && !it->second.isSpectator && it->second.index == size) {
               cout << "Using variable " << it->second.name << " in the TMVA::Reader" << endl;
               if(i<3)
                  BDTReaders[i].second->AddVariable(it->second.name,&it->second.value);
               else
                  BDTReaders[i].second->AddVariable(it->second.definition,&it->second.value);
               size++;
            }
         }
         //cout << "size = " << size << "\tfirstUsed->second.maxIndex = " << firstUsed->second.maxIndex << endl;
      }
      size = 0;
      while(size!=firstUsedSpectator->second.maxIndex && firstUsedSpectator!=BDTReaders[i].first.end()) {
         for(map<TString,MVAVar>::iterator it=BDTReaders[i].first.begin(); it!=BDTReaders[i].first.end(); ++it) {
            if(it->second.use && it->second.isSpectator && it->second.index == size) {
               cout << "Using spectator " << it->second.name << " in the TMVA::Reader" << endl;
               BDTReaders[i].second->AddSpectator(it->second.name,&it->second.value);
               size++;
            }
         }
         //cout << "size = " << size << "\tfirstUsedSpectator->second.maxIndex = " << firstUsedSpectator->second.maxIndex << endl;
      }

      cout << "\tBooking the reader ... ";
      BDTReaders[i].second->BookMVA("BDT method", weightfiles[i]);
      cout << "DONE" << endl;
   }

}//setAllBDTReaders

pair<int,int> MicroNtupleMaker::setMVAVar(EventNtuple * ntuple, MicroNtuple * mnt) {
   pair<int,int> ret;

   if(ntuple->jLV.size()==2 && ntuple->getNBTags()==0) {
      ret = make_pair(0,3);

      //UVa
      BDTReaders[0].first["leptonPt"].value = ntuple->lLV[0].Pt();
      BDTReaders[0].first["lepMT"].value = TMath::Sign(TMath::Sqrt(TMath::Abs(TMath::Power(ntuple->lLV[0].Pt()+ntuple->METLV[0].Pt(),2)-(TMath::Power(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+TMath::Power(0,2)))),TMath::Power(ntuple->lLV[0].Pt()+ntuple->METLV[0].Pt(),2)-(TMath::Power(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+TMath::Power(0,2)));
      BDTReaders[0].first["jet1dRLep"].value = ntuple->jLV[0].DRlj;
      BDTReaders[0].first["jet2dRLep"].value = ntuple->jLV[1].DRlj;
      double tmpHT = 0.0;
      for(unsigned int i=0; i<ntuple->jLV.size(); i++) {
         tmpHT+=ntuple->jLV[i].Et();
      }
      tmpHT+=ntuple->lLV[0].Pt();
      BDTReaders[0].first["ht"].value = tmpHT;
      BDTReaders[0].first["Ptlnujj"].value = TMath::Sqrt(TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2));
      BDTReaders[0].first["dRlepjj"].value = mnt->dRlepjj;
      BDTReaders[0].first["dPhiMETJet"].value = mnt->dPhiMETJet;
      BDTReaders[0].first["dPhiJetJet"].value = mnt->dPhiJetJet;
      BDTReaders[0].first["CosTheta_l"].value = mnt->CosTheta_l;
      BDTReaders[0].first["CosTheta_WH"].value = mnt->CosTheta_WH;

      //Combined
      BDTReaders[6].first["leptonPt"].value = ntuple->lLV[0].Pt();
      BDTReaders[6].first["lepMT"].value = TMath::Sign(TMath::Sqrt(TMath::Abs(TMath::Power(ntuple->lLV[0].Pt()+ntuple->METLV[0].Pt(),2)-(TMath::Power(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+TMath::Power(0,2)))),TMath::Power(ntuple->lLV[0].Pt()+ntuple->METLV[0].Pt(),2)-(TMath::Power(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+TMath::Power(0,2)));
      BDTReaders[6].first["jet1dRLep"].value = ntuple->jLV[0].DRlj;
      BDTReaders[6].first["jet2dRLep"].value = ntuple->jLV[1].DRlj;
      tmpHT = 0.0;
      for(unsigned int i=0; i<ntuple->jLV.size(); i++) {
         tmpHT+=ntuple->jLV[i].Et();
      }
      tmpHT+=ntuple->lLV[0].Pt();
      BDTReaders[6].first["ht"].value = tmpHT;
      BDTReaders[6].first["Ptlnujj"].value = TMath::Sqrt(TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2));
      BDTReaders[6].first["dRlepjj"].value = mnt->dRlepjj;
      BDTReaders[6].first["dPhiMETJet"].value = mnt->dPhiMETJet;
      BDTReaders[6].first["dPhiJetJet"].value = mnt->dPhiJetJet;
      BDTReaders[6].first["CosTheta_l"].value = mnt->CosTheta_l;
      BDTReaders[6].first["CosTheta_WH"].value = mnt->CosTheta_WH;
   }
   else if(ntuple->jLV.size()==3 && ntuple->getNBTags()==0) {
      ret = make_pair(1,4);

      //UVa
      BDTReaders[1].first["leptonPt"].value = ntuple->lLV[0].Pt();
      BDTReaders[1].first["leptonEtaCharge"].value = ntuple->lLV[0].Eta()*ntuple->lLV[0].lQ;
      BDTReaders[1].first["jet2dRLep"].value = ntuple->jLV[1].DRlj;
      BDTReaders[1].first["jet3dRLep"].value = ntuple->jLV[2].DRlj;
      double tmpHT = 0.0;
      for(unsigned int i=0; i<ntuple->jLV.size(); i++) {
         tmpHT+=ntuple->jLV[i].Et();
      }
      tmpHT+=ntuple->lLV[0].Pt();
      BDTReaders[1].first["ht"].value = tmpHT;
      BDTReaders[1].first["Mlnujj"].value = TMath::Sign(TMath::Sqrt(TMath::Abs(TMath::Power(ntuple->jLV[0].E()+ntuple->jLV[1].E()+ntuple->lLV[0].E()+ntuple->METLV[0].E(),2)-(TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+TMath::Power(ntuple->jLV[0].Pz()+ntuple->jLV[1].Pz()+ntuple->lLV[0].Pz()+ntuple->METLV[0].Pz(),2)))),TMath::Power(ntuple->jLV[0].E()+ntuple->jLV[1].E()+ntuple->lLV[0].E()+ntuple->METLV[0].E(),2)-(TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+TMath::Power(ntuple->jLV[0].Pz()+ntuple->jLV[1].Pz()+ntuple->lLV[0].Pz()+ntuple->METLV[0].Pz(),2)));
      BDTReaders[1].first["dRlepjj"].value = mnt->dRlepjj;
      BDTReaders[1].first["dPhiMETJet"].value = mnt->dPhiMETJet;
      BDTReaders[1].first["dEtaJetJet"].value = TMath::Abs(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta());
      BDTReaders[1].first["minDPhiLepJet"].value = mnt->minDPhiLepJet;
      BDTReaders[1].first["CosTheta_l"].value = mnt->CosTheta_l;
      BDTReaders[1].first["CosTheta_j"].value = mnt->CosTheta_j;
      BDTReaders[1].first["CosTheta_WH"].value = mnt->CosTheta_WH;

      //Combined
      BDTReaders[7].first["leptonPt"].value = ntuple->lLV[0].Pt();
      BDTReaders[7].first["leptonEtaCharge"].value = ntuple->lLV[0].Eta()*ntuple->lLV[0].lQ;
      BDTReaders[7].first["jet2dRLep"].value = ntuple->jLV[1].DRlj;
      BDTReaders[7].first["jet3dRLep"].value = ntuple->jLV[2].DRlj;
      tmpHT = 0.0;
      for(unsigned int i=0; i<ntuple->jLV.size(); i++) {
         tmpHT+=ntuple->jLV[i].Et();
      }
      tmpHT+=ntuple->lLV[0].Pt();
      BDTReaders[7].first["ht"].value = tmpHT;
      BDTReaders[7].first["Mlnujj"].value = TMath::Sign(TMath::Sqrt(TMath::Abs(TMath::Power(ntuple->jLV[0].E()+ntuple->jLV[1].E()+ntuple->lLV[0].E()+ntuple->METLV[0].E(),2)-(TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+TMath::Power(ntuple->jLV[0].Pz()+ntuple->jLV[1].Pz()+ntuple->lLV[0].Pz()+ntuple->METLV[0].Pz(),2)))),TMath::Power(ntuple->jLV[0].E()+ntuple->jLV[1].E()+ntuple->lLV[0].E()+ntuple->METLV[0].E(),2)-(TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+TMath::Power(ntuple->jLV[0].Pz()+ntuple->jLV[1].Pz()+ntuple->lLV[0].Pz()+ntuple->METLV[0].Pz(),2)));
      BDTReaders[7].first["dRlepjj"].value = mnt->dRlepjj;
      BDTReaders[7].first["dPhiMETJet"].value = mnt->dPhiMETJet;
      BDTReaders[7].first["dEtaJetJet"].value = TMath::Abs(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta());
      BDTReaders[7].first["minDPhiLepJet"].value = mnt->minDPhiLepJet;
      BDTReaders[7].first["CosTheta_l"].value = mnt->CosTheta_l;
      BDTReaders[7].first["CosTheta_j"].value = mnt->CosTheta_j;
      BDTReaders[7].first["CosTheta_WH"].value = mnt->CosTheta_WH;
   }
   else if(ntuple->jLV.size()>=4 && ntuple->getNBTags()==0) {
      ret = make_pair(2,5);

      //UVa
      BDTReaders[2].first["leptonEtaCharge"].value = ntuple->lLV[0].Eta()*ntuple->lLV[0].lQ;
      BDTReaders[2].first["jet2dRLep"].value = ntuple->jLV[1].DRlj;
      BDTReaders[2].first["jet3dRLep"].value = ntuple->jLV[2].DRlj;
      double tmpHT = 0.0;
      for(unsigned int i=0; i<ntuple->jLV.size(); i++) {
         tmpHT+=ntuple->jLV[i].Et();
      }
      tmpHT+=ntuple->lLV[0].Pt();
      BDTReaders[2].first["ht"].value = tmpHT;
      BDTReaders[2].first["Mlnujj"].value = TMath::Sign(TMath::Sqrt(TMath::Abs(TMath::Power(ntuple->jLV[0].E()+ntuple->jLV[1].E()+ntuple->lLV[0].E()+ntuple->METLV[0].E(),2)-(TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+TMath::Power(ntuple->jLV[0].Pz()+ntuple->jLV[1].Pz()+ntuple->lLV[0].Pz()+ntuple->METLV[0].Pz(),2)))),TMath::Power(ntuple->jLV[0].E()+ntuple->jLV[1].E()+ntuple->lLV[0].E()+ntuple->METLV[0].E(),2)-(TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+TMath::Power(ntuple->jLV[0].Pz()+ntuple->jLV[1].Pz()+ntuple->lLV[0].Pz()+ntuple->METLV[0].Pz(),2)));
      BDTReaders[2].first["dPhiMETJet"].value = mnt->dPhiMETJet;
      BDTReaders[2].first["dPhiMETLep"].value = TMath::Abs(ntuple->METLV[0].Phi()-ntuple->lLV[0].Phi());

      //Combined
      BDTReaders[8].first["leptonEtaCharge"].value = ntuple->lLV[0].Eta()*ntuple->lLV[0].lQ;
      BDTReaders[8].first["jet2dRLep"].value = ntuple->jLV[1].DRlj;
      BDTReaders[8].first["jet3dRLep"].value = ntuple->jLV[2].DRlj;
      tmpHT = 0.0;
      for(unsigned int i=0; i<ntuple->jLV.size(); i++) {
         tmpHT+=ntuple->jLV[i].Et();
      }
      tmpHT+=ntuple->lLV[0].Pt();
      BDTReaders[8].first["ht"].value = tmpHT;
      BDTReaders[8].first["Mlnujj"].value = TMath::Sign(TMath::Sqrt(TMath::Abs(TMath::Power(ntuple->jLV[0].E()+ntuple->jLV[1].E()+ntuple->lLV[0].E()+ntuple->METLV[0].E(),2)-(TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+TMath::Power(ntuple->jLV[0].Pz()+ntuple->jLV[1].Pz()+ntuple->lLV[0].Pz()+ntuple->METLV[0].Pz(),2)))),TMath::Power(ntuple->jLV[0].E()+ntuple->jLV[1].E()+ntuple->lLV[0].E()+ntuple->METLV[0].E(),2)-(TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+TMath::Power(ntuple->jLV[0].Pz()+ntuple->jLV[1].Pz()+ntuple->lLV[0].Pz()+ntuple->METLV[0].Pz(),2)));
      BDTReaders[8].first["dPhiMETJet"].value = mnt->dPhiMETJet;
      BDTReaders[8].first["dPhiMETLep"].value = TMath::Abs(ntuple->METLV[0].Phi()-ntuple->lLV[0].Phi());
   }

   //TAMU+Combined
   if(ntuple->jLV.size()>=2) {
      for(unsigned int ireader=3; ireader<BDTReaders.size(); ireader++) {
         for(int i=0; i<13; i++)
            BDTReaders[ireader].first[Form("logEventProb%i",i)].value = TMath::Log(mnt->eventProb[i]);
         BDTReaders[ireader].first["logEventProb19"].value = TMath::Log(mnt->eventProb[19]);
         BDTReaders[ireader].first["logEventProb54"].value = TMath::Log(mnt->eventProb[54]);
         BDTReaders[ireader].first["event"].value = ntuple->event;
         BDTReaders[ireader].first["lumi"].value = ntuple->lumi;
         BDTReaders[ireader].first["run"].value = ntuple->run;
      }
   }

   return ret;

}//setMVAVar

bool MicroNtupleMaker::copyEventProbs(EventNtuple* mergeEventNtuple, EventNtuple* eventNtuple, METree* meNtuple, MicroNtuple* microNtuple) {
   //assert(meNtuple->getNProbStat()>0);
   if(meNtuple->getNProbStat()==0) {
      cout << "WARNING::MicroNtupleMaker::copyEventProbs Skipping event because NProbStat==0" << endl
           << "\tProbably a TTbar event" << endl;
      return false;
   }
   for (int i = 0; i < meNtuple->getNProbStat(); ++i){
      microNtuple->eventProb[i]    = meNtuple->getProbStat(i)->tEventProb;
      if(currentProcess.Contains("TTbar") && meNtuple->getProbStat(i)->tEventProb==0) {
         cout << "WARNING::MicroNtupleMaker::copyEventProbs Skipping event because tEventProb for ProbStat(" << i << ")==0" << endl
              << "\tThis is a TTbar event." << endl;
         return false;
      }
      else if(currentProcess.CompareTo("TTbar")!=0) {
         assert(meNtuple->getProbStat(i)->tEventProb!=0);
         assert(microNtuple->eventProb[i]!=0);
      }
      microNtuple->eventMaxProb[i] = meNtuple->getProbStat(i)->tEventMaxProb;  
   }
   return true;   
}

void MicroNtupleMaker::setSizeRunEvent(EventNtuple* eventNtuple, METree* meNtuple, MicroNtuple* microNtuple) {
   microNtuple->size = meNtuple->getNProbStat();
   microNtuple->run = eventNtuple->run;//meNtuple->getRun();
   microNtuple->event = eventNtuple->event;//meNtuple->getEvent();
}

void MicroNtupleMaker::setEPDs(METree* meNtuple, MicroNtuple* microNtuple) {
   microNtuple->epdPretagWWandWZ = microNtuple->calcWZEPD(DEFS::pretag);
   microNtuple->epd0tagWWandWZ = microNtuple->calcWZEPD(DEFS::eq0tag);
   microNtuple->epd1tagWWandWZ = microNtuple->calcWZEPD(DEFS::eq1tag);
   microNtuple->epd2tagWWandWZ = microNtuple->calcWZEPD(DEFS::eq2tag);
   int counterHWW = 0;
   int counterWH = 0;
   int counterH = 0;
   //If using TArrayD
   //microNtuple->epdPretagHiggs.Set(MicroNtuple::nHiggsMasses);
   for(int i = 0; i < meNtuple->getNProbStat(); ++i) {
      if(meNtuple->getProbStat(i)->tmeType == DEFS::EP::HWW) {
         microNtuple->epd1tagHWW[counterHWW] = microNtuple->calcHWWEPD(DEFS::eq1tag,meNtuple->getProbStat(i)->tmeParam);
         microNtuple->epd2tagHWW[counterHWW++] = microNtuple->calcHWWEPD(DEFS::eq2tag,meNtuple->getProbStat(i)->tmeParam);
      }
      if(meNtuple->getProbStat(i)->tmeType == DEFS::EP::WH) {
         microNtuple->epd1tagWH[counterWH] = microNtuple->calcWHEPD(DEFS::eq1tag,meNtuple->getProbStat(i)->tmeParam);
         microNtuple->epd2tagWH[counterWH++] = microNtuple->calcWHEPD(DEFS::eq2tag,meNtuple->getProbStat(i)->tmeParam);
      }
      if(meNtuple->getProbStat(i)->tmeType == DEFS::EP::HWW || meNtuple->getProbStat(i)->tmeType == DEFS::EP::WH) {
         //Used to absorb the error of the missing first entry
         microNtuple->absorbError[counterH] = -1;
         microNtuple->epdPretagHiggs[counterH] = microNtuple->calcHiggsEPD(DEFS::pretag,meNtuple->getProbStat(i)->tmeParam);
         microNtuple->epd1tagHiggs[counterH] = microNtuple->calcHiggsEPD(DEFS::eq1tag,meNtuple->getProbStat(i)->tmeParam);
         microNtuple->epd2tagHiggs[counterH++] = microNtuple->calcHiggsEPD(DEFS::eq2tag,meNtuple->getProbStat(i)->tmeParam);
      }
   }
}

void MicroNtupleMaker::setTMVAInformation(EventNtuple* eventNtuple, MicroNtuple* microNtuple) {
   microNtuple->dPhiJetJet = eventNtuple->getDeltaPhiJetJet();
   microNtuple->dPhiMETJet = eventNtuple->getDeltaPhiMETJet();
   microNtuple->minDPhiMETJet = eventNtuple->getMinDeltaPhiMETJet();
   microNtuple->minDPhiLepJet = eventNtuple->getMinDPhiLepJet();
   eventNtuple->getAngularVariables(microNtuple->Cos_dPhiWW,microNtuple->Cos_dPhiWH,microNtuple->CosTheta_l,
                                    microNtuple->CosTheta_j,microNtuple->CosTheta_WH,microNtuple->JacksonAngle);
   microNtuple->JacobePeak = eventNtuple->getJacobePeak();
   microNtuple->dRlepjj = eventNtuple->getDeltaRlepjj();
   microNtuple->sumJetEt = eventNtuple->getSumJetEt();

   if(fillBDT) {
      pair<int,int> readerIndex = setMVAVar(eventNtuple,microNtuple);
      //if(eventNtuple->jLV.size()>1 && eventNtuple->jLV.size()<5 && eventNtuple->getNBTags()==0) {
      if(eventNtuple->jLV.size()>1 && eventNtuple->getNBTags()==0) {
         microNtuple->KinBDT = BDTReaders[readerIndex.first].second->EvaluateMVA("BDT method");
         microNtuple->MEBDT = BDTReaders[readerIndex.second].second->EvaluateMVA("BDT method");
         microNtuple->KinMEBDT = BDTReaders[readerIndex.second+3].second->EvaluateMVA("BDT method");
      }
      else {
         microNtuple->KinBDT = -9999;
         microNtuple->MEBDT = -9999;
         microNtuple->KinMEBDT = -9999;
      }
   }
   else {
      microNtuple->KinBDT = -9999;
      microNtuple->MEBDT = -9999;
      microNtuple->KinMEBDT = -9999;
   }
}

void MicroNtupleMaker::updateMicroNtuple(TString outputPath, TString updateMicroNtuple) {

   TFile* file_to_update = TFile::Open(updateMicroNtuple,"READ");
   TTree* tree_to_update = (TTree*)gDirectory->Get("METree");
   if(!tree_to_update) {
      cout << "ERROR::MicroNtupleMaker::updateMicroNtuple Could not find the tree METree in the file " << updateMicroNtuple << endl;
      assert(tree_to_update); 
   }

   cout << "MicroNtupleMaker::updateMicroNtuple Make ntuples ... ";
   // Create the base objects
   EventNtuple * eventNtuple      = new EventNtuple();
   METree      * meNtuple         = new METree();
   MicroNtuple * microNtuple      = new MicroNtuple(2);
   cout << "DONE" << endl;

   cout << "MicroNtupleMaker::updateMicroNtuple Set EvtTree, METree, and MicroNtuple branches ... ";
   tree_to_update->SetBranchStatus("mnt*",0);
   tree_to_update->SetBranchAddress("EvtTree", &eventNtuple);
   tree_to_update->SetBranchAddress("METree", &meNtuple);
   cout << "DONE" << endl;   

   cout << "MicroNtupleMaker::updateMicroNtuple Create output file and setup output tree ... ";
   // Create and output file and clone the tree that will be in the output and set microNtuple that fills it
   TString output = outputPath+"micro"+currentProcess+"_BaseCuts.root";
   TFile outputFile(output, "RECREATE");
   TTree* outputTree = new TTree(tree_to_update->GetName(),tree_to_update->GetTitle());
   outputTree->Branch("EvtTree", "EventNtuple", &eventNtuple);
   outputTree->Branch("METree", "METree", &meNtuple);
   outputTree->Branch("mnt", "MicroNtuple", &microNtuple);
   cout << "DONE" << endl << endl << endl;

   // The index map has not been filled yet
   imFilled = false;

   //
   // If the paths for the BDT weights are set, store them in the tree
   //
   if (fillBDT) {
      setBDTReadersFromTable();
   }

   // Get the entries and report if zero
   nentries = static_cast<unsigned>(tree_to_update->GetEntries());
   cout << endl << endl << "Original chain has " << nentries << " entries" << endl;
   if (nentries == 0){
      cout << "\t\tNo entries found!  Aborting.\n";
      return;
   }

   // If the start and end entries are not 0 and -1 respectively then only loop through the requested entries
   startEndEntries.first!=0 ? startEntry = startEndEntries.first : startEntry = 0;
   startEndEntries.second!=-1 ? endEntry = startEndEntries.second : endEntry = nentries;
   // Check if the user made a mistake and the endEntry requested is greater than the last entry in the TChain or TTree
   endEntry>nentries ? endEntry = nentries : endEntry = endEntry;
   cout << "MicroNtupleMaker::updateMicroNtuple Starting with entry " << startEntry << " and ending with entry " << endEntry << endl;

   //Loop over all the entries in the original chain. The idea is to copy everything to microNtuple
   cout << "MicroNtupleMaker::updateMicroNtuple Filling MicroNtuple ..." << endl;

   for (unsigned ientry = startEntry; ientry < endEntry; ++ientry){

      loadbar2(ientry+1,nentries);

      // get the entry
      tree_to_update->GetEntry(ientry);

      if (!imFilled) {
         microNtuple->indexMap = meNtuple->fillIndexMap();
         imFilled = true;
      }

      // clear the microNtuple
      microNtuple->clear();

      // First copy all the event probabilities
      if(meNtuple->getNProbStat()<1) {
         cout << endl << endl << "ientry = " << ientry << "\tThis is eventNtuple->event " << eventNtuple->event
              << ", meNtuple->getEvent() " << meNtuple->getEvent() << endl << endl;
         cout << endl << endl<< "There is a problem with getNProbStat (size = " << meNtuple->getNProbStat()
                                          << ") in entry " << meNtuple->getEvent() << " or " << eventNtuple->event << endl << endl; 
         cout << endl << endl << endl;
         continue;
      }
      bool passCopy = copyEventProbs(0,eventNtuple,meNtuple,microNtuple);
      if(!passCopy) {
         microNtuple->clear();
         continue;
      }
      setSizeRunEvent(eventNtuple,meNtuple,microNtuple);
      setEPDs(meNtuple, microNtuple);
      setTMVAInformation(eventNtuple,microNtuple);
      microNtuple->reader = 0;

      //
      // Finalize the output tree
      //

      outputTree->Fill();

      microNtuple->clear();
      meNtuple->clear();
   }//for entries

   //
   // Report some results
   //
   cout << endl << endl << "Wrote " << output << " with " << outputTree->GetEntries() << " entries" << endl;
   outputFile.Write();
   outputFile.Close();
   file_to_update->Close();
   delete meNtuple;
   delete microNtuple;
   delete eventNtuple;    
}



