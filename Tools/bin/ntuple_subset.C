//STD include
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <utility>
#include <cmath>
//Root Files
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TTreeIndex.h"
//TAMUWW
#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"
#include "TAMUWW/MEPATNtuple/interface/METree.hh"
#include "TAMUWW/MEPATNtuple/interface/MicroNtuple.hh"
#include "TAMUWW/Selection/bin/RunEventSet.h"


using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::time_t;
using std::ios;
using std::pair;

bool cuts(EventNtuple* ntuple) {
   if(ntuple->METLV[0].Pt() < 25.0)
      return false;
   if((ntuple->lLV[0].lpfIso < 0.2 || ntuple->lLV[0].lpfIso > 2.0))
      return false;
   if(ntuple->lLV[0].leptonCat == DEFS::electron && ntuple->lLV[0].Pt() < 30)
     return false;
   if(ntuple->lLV[0].leptonCat == DEFS::muon && ntuple->lLV[0].Pt() < 25)
     return false;
   if (ntuple->jLV[0].Pt() < 30 || ntuple->jLV[1].Pt() < 25)
      return false;

   return true;
}

/*****************************************************************************/
static inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50) {
   if ( (x != n) && (x % (n/100) != 0) ) return;
 
   float ratio  =  x/(float)n;
   int   c      =  ratio * w;
 
   cout << setw(3) << (int)(ratio*100) << "% [";
   for (int x=0; x<c; x++) cout << "=";
   for (unsigned int x=c; x<w; x++) cout << " ";
   cout << "]\r" << flush;
}

/*****************************************************************************/
void ntuple_subset(TString inputRootFile, TString baseFilename, bool useCuts = false, int evtsPerFile = -1) {
   TFile* ifile = TFile::Open(inputRootFile);
   ifile->cd("PS");

   //Set up date and Root Chain
   TChain *t1 = (TChain*)gDirectory->Get("jets2p");
   EventNtuple *ntuple = new EventNtuple();
   t1->SetBranchAddress("EvtNtuple", &ntuple);
   int nEventNtuple = static_cast<int>(t1->GetEntries());

   int numFiles = 0;
   if (evtsPerFile==-1)
      numFiles = 1;
   else
      numFiles = std::ceil(double(nEventNtuple)/evtsPerFile);

   for(int f=0; f<numFiles; f++) {
      //Opens the .ROOT file and set it so that new graph will be added, but nothing will be deleted
      TString index = Form("%d",f);
      TFile myfile(TString(baseFilename+index+".root"), "RECREATE");
      myfile.mkdir("PS");
      myfile.cd("PS");
      cout << "File " << f << " of " << numFiles << endl;

      TChain *newntuple = (TChain*)t1->CloneTree(0);
      if(numFiles == 1) evtsPerFile = nEventNtuple;
      for(int i=0; i<evtsPerFile; i++)
      {
         t1->GetEntry((f*evtsPerFile)+i);
         if(!useCuts || (useCuts && cuts(ntuple)))
            newntuple->Fill();
         if (i%10000==0) cout << "\tEntry " << (f*evtsPerFile)+i << " of " << nEventNtuple << endl;
      }
      
      newntuple->AutoSave();
      newntuple->Write();
      myfile.Close();
   }
}


/*****************************************************************************/
void microNtuple_subtraction(TString inputRootFile, TString compareTo, TString baseFilename) {

   cout << "Opening files ... ";
   TFile* ifile = TFile::Open(inputRootFile);
   TFile* compareToFile = TFile::Open(compareTo);
   cout << "DONE" << endl;

   //Set up the EventNtuple
   cout << "Setting up EventNtuple ... " << flush;
   compareToFile->cd("PS");
   TChain *t1 = (TChain*)gDirectory->Get("jets2p");
   EventNtuple *ntuple = new EventNtuple();
   t1->SetBranchAddress("EvtNtuple", &ntuple);
   t1->SetBranchStatus("*",0);
   t1->SetBranchStatus("event",1);
   int nEventNtuple = static_cast<int>(t1->GetEntries());
   cout << "DONE" << endl;

   //Set up the MicroNtuple
   cout << "Setting up original MicroNtuple ... " << flush;
   ifile->cd();
   TChain *t2 = (TChain*)gDirectory->Get("METree");
   MicroNtuple *mnt = new MicroNtuple();
   t2->SetBranchAddress("mnt", &mnt);
   t2->GetTree()->BuildIndex("m_event");
   TTreeIndex* microNtupleIndex = new TTreeIndex();
   microNtupleIndex->Append((TTreeIndex*)t2->GetTree()->GetTreeIndex());
   //int nMicroNtuple = static_cast<int>(t2->GetEntries());
   cout << "DONE" << endl;

   cout << "Building the eventIndex ... " << endl;
   Long64_t local = -1;
   t2->SetBranchStatus("*",0);
   t2->SetBranchStatus("event*",1);
   map<int,int> eventIndex;
   map<int,int> duplicateIndex;
   int nIndex = microNtupleIndex->GetN();
   for( int i = 0; i < nIndex ; ++i ) {
      loadbar(i+1,nIndex);
      local = microNtupleIndex->GetIndex()[i];
      if(eventIndex.find(microNtupleIndex->GetIndexValues()[i]>>31)!=eventIndex.end()) {
         duplicateIndex[local] = microNtupleIndex->GetIndexValues()[i]>>31;
         cout << "WARNING::makeMicroNtuple::There could be duplicate events in the chain" << endl;
      }
      if(local>t2->GetEntries()) {
         cout << "WARNING::makeMicroNtuple::The local index is greater than the number of entries in the chain" << endl;
      }
      eventIndex[microNtupleIndex->GetIndexValues()[i]>>31] = local;
   }
   t2->SetBranchStatus("*",1);

   //Opens the .ROOT file and set it so that new graph will be added, but nothing will be deleted
   cout << endl << "Setting up the output file and chain ... ";
   TFile myfile(TString(baseFilename+".root"), "RECREATE");   
   TChain *newntuple = (TChain*)t2->CloneTree(0);
   cout << "DONE" << endl;

   //counters
   int eventsRemoved = 0;
   int duplicates = 0;
   // Holds the list of events to avoid duplication
   RunEventSet runEventSet;

   cout << "Looping through the EventNtuple and filling subtracted MicroNtuple ... " << endl;
   for(int ientry=0; ientry<nEventNtuple; ientry++)
   {
      t1->GetEntry(ientry);
      loadbar(ientry+1,nEventNtuple);
      if (runEventSet.alreadySeen(ntuple->run, ntuple->event)) {
         cout << "WARNING, event Run = " << ntuple->run << ", Event = " << ntuple->event
              <<" is duplicated" << endl << "We will skip this event." << endl;
         duplicates++;
         eventsRemoved++;
      }
      else if(eventIndex.find(ntuple->event)!=eventIndex.end()) {
         t2->GetEntry(eventIndex[ntuple->event]);
         newntuple->Fill();
      }
      else
         eventsRemoved++;
   }
   
   cout << endl << "Events Removed: " << eventsRemoved << endl;
   cout << "Events Remaining: " << newntuple->GetEntries() << endl;
   cout << "Duplicate Events: " << duplicateIndex.size() << " (" << duplicates << ")" << endl;

   newntuple->AutoSave();
   newntuple->Write();
   myfile.Close();
}

void microNtuple_addition(TString inputRootFile, TString compareTo, TString addFrom, TString baseFilename) {

   cout << "Opening files ... ";
   //TFile* ifile = TFile::Open(inputRootFile);
   //TFile* addFromFile = TFile::Open(addFrom);
   //TFile* compareToFile = TFile::Open(compareTo);
   cout << "DONE" << endl;

   cout << "Making all of the ntuples ... ";
   EventNtuple *ntuple      = new EventNtuple();
   EventNtuple *ntupleFinal = new EventNtuple();
   METree *metree           = new METree();
   MicroNtuple *mnt         = new MicroNtuple();
   cout << "DONE" << endl;

   //Set up the EventNtuple
   cout << "Setting up EventNtuple ... " << flush;
   //compareToFile->cd("PS");
   TChain *t1 = new TChain("jets2p","Output tree for matrix element");//(TChain*)gDirectory->Get("jets2p");
   t1->Add(compareTo+"/PS/jets2p");
   t1->LoadTree(0);
   TTree* t1t = t1->GetTree();
   t1t->SetCacheSize(10000000);
   t1t->AddBranchToCache("*");
   t1->SetBranchAddress("EvtNtuple", &ntuple);
   t1->SetBranchStatus("*",0);
   t1->SetBranchStatus("event",1);
   t1->GetTree()->BuildIndex("event");
   TTreeIndex* eventNtupleIndex = new TTreeIndex();
   eventNtupleIndex->Append((TTreeIndex*)t1->GetTree()->GetTreeIndex());
   int nEventNtuple = static_cast<int>(t1->GetEntries());
   t1->SetBranchStatus("*",1);
   cout << "DONE" << endl;

   //Set up the final MicroNtuple
   cout << "Setting up original/final MicroNtuple ... ";
   //ifile->cd();
   TChain *t2 = new TChain("METree");//(TChain*)gDirectory->Get("METree");
   t2->Add(inputRootFile);
   t2->SetBranchAddress("EvtTree",&ntupleFinal);
   t2->SetBranchAddress("METree",&metree);
   t2->SetBranchAddress("mnt", &mnt);
   cout << "DONE" << endl;

   //Set up the MicroNtuple you will be adding events from
   cout << "Setting up MicroNtuple that has the events to be added ... ";
   //addFromFile->cd();
   TChain* t3 = new TChain("METree");//(TChain*)gDirectory->Get("METree");
   t3->Add(addFrom);
   t3->SetBranchAddress("EvtTree",&ntupleFinal);
   t3->SetBranchAddress("METree",&metree);
   t3->SetBranchAddress("mnt", &mnt);
   //t3->GetTree()->BuildIndex("m_event");
   //TTreeIndex* microNtupleIndex = new TTreeIndex();
   //microNtupleIndex->Append((TTreeIndex*)t3->GetTree()->GetTreeIndex());
   int nMicroNtuple = static_cast<int>(t3->GetEntries());
   cout << "DONE" << endl;

   cout << "Building the eventIndex ... " << endl;
   Long64_t local = -1;
   t3->SetBranchStatus("*",0);
   t3->SetBranchStatus("event*",1);
   map<int,int> eventIndex;
   int nIndex = eventNtupleIndex->GetN();
   for( int i = 0; i < nIndex ; ++i ) {
      loadbar(i+1,nIndex);
      local = eventNtupleIndex->GetIndex()[i];
      eventIndex[eventNtupleIndex->GetIndexValues()[i]>>31] = local;
   }
   t3->SetBranchStatus("*",1);

   //counters
   int eventsSkipped = 0;
   int eventsAdded = 0;
   
   cout << endl << "Looping through the additional MicroNtuple and checking the EventNtuple ... " << endl;
   //ifile->cd();
   for(int ientry=0; ientry<nMicroNtuple; ientry++)
   {
      //cout << "sfsg1.1" << endl; 
      t3->GetEntry(ientry);
      //cout << "sfsg1.2" << endl; 
      loadbar(ientry+1,nMicroNtuple);
      if(eventIndex.find(mnt->event)!=eventIndex.end()) {
         //cout << "sfsg1.3" << endl; 
         //t3->GetEntry(eventIndex[ntuple->event]);
         //cout << "sfsg1.4" << endl; 
         //t2->Fill();
         //cout << "sfsg1.5" << endl; 
         eventsAdded++;
      }
      else
         eventsSkipped++;
   }

   //Opens the .ROOT file and set it so that new graph will be added, but nothing will be deleted
   //cout << endl << "Setting up the output file ... ";
   //TFile myfile(TString(baseFilename+".root"), "RECREATE");   
   //cout << "DONE" << endl;

   if(eventsSkipped == 0 && eventsAdded == t3->GetEntries()) {
      cout << endl << "Merging the original/final MicroNtuple with the additional MicroNtuple ... " << flush;
      t2->Add(addFrom);
      t2->Merge(TString(baseFilename+".root"));
      cout << "DONE" << endl;
   }
   else {
      cout << endl << "WARNING::No merge because either eventsSkipped [" << eventsSkipped << "] is greater than 0 or eventsAdded("
           << eventsAdded << ") does not equal t3->GetEntries() [" << t3->GetEntries() << "]" << endl;
   }
   
   cout << endl << "Events Skipped: " << eventsSkipped << endl;
   cout << "Events Added: " << eventsAdded << endl;
   cout << "Events Total: " << t2->GetEntries() << endl;
   cout << "Difference from EventNtuple: " << t1->GetEntries() - t2->GetEntries() << endl;

   //t2->AutoSave();
   //t2->Write();
   //myfile.Close();
}

void splitNtupleByGenInfo(TString inputFilepath = "/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/", TString outputDir = "./", TString hTo = "WW", bool eventNtupleOnly = true, TString postfix = "") {

   TString folderName = "";
   if(eventNtupleOnly)
      folderName = "PS/";
   TString inputTreeName = "";
   if(eventNtupleOnly)
      inputTreeName = "jets2p";
   else
      inputTreeName = "METree";

   cout << "Opening input files ... ";
   TFile* ifile;
   TString ifilename;
   if(eventNtupleOnly) {
      ifilename = inputFilepath+"WH_ZH_TTH_HTo"+hTo+"_M125"+postfix+".root";
      ifile = new TFile(ifilename.Data(),"READ");
   }
   else {
      ifilename = inputFilepath+"microWH_ZH_TTH_HTo"+hTo+"_M125"+postfix+"_BaseCuts.root";
      ifile = new TFile(ifilename.Data(),"READ");
   }
   cout << "DONE" << endl;
   cout << "\t" << ifilename << endl;

   cout << "Making all of the input ntuples ... ";
   EventNtuple *ntuple      = new EventNtuple();
   METree *metree           = new METree();
   MicroNtuple *mnt         = new MicroNtuple();
   cout << "DONE" << endl;

   cout << "Setting up EventNtuple ... " << flush;
   TTree *t1 = (TTree*)gDirectory->Get(folderName+inputTreeName);
   if(eventNtupleOnly)
      t1->SetBranchAddress("EvtNtuple", &ntuple);
   else
      t1->SetBranchAddress("EvtTree", &ntuple);
   int nEventNtuple = static_cast<int>(t1->GetEntries());
   cout << "DONE" << endl;

   cout << endl << "Setting up the output file and chain ... " << flush;
   TString prefix = "";
   if(!eventNtupleOnly) {
      prefix = "micro";
      postfix += "_BaseCuts";
   }
   TFile WHFile(outputDir+prefix+"WH_HTo"+hTo+"_M125"+postfix+".root", "RECREATE");
   if(eventNtupleOnly) {
      WHFile.mkdir("PS");
      WHFile.cd("PS");
   }
   TChain *WHEventNtuple = (TChain*)t1->CloneTree(0);

   TFile ZHFile(outputDir+prefix+"ZH_HTo"+hTo+"_M125"+postfix+".root", "RECREATE");
   if(eventNtupleOnly) {
      ZHFile.mkdir("PS");
      ZHFile.cd("PS");
   }
   TChain *ZHEventNtuple = (TChain*)t1->CloneTree(0);

   TFile TTHFile(outputDir+prefix+"TTH_HTo"+hTo+"_M125"+postfix+".root", "RECREATE");
   if(eventNtupleOnly) {
      TTHFile.mkdir("PS");
      TTHFile.cd("PS");
   }
   TChain *TTHEventNtuple = (TChain*)t1->CloneTree(0);
   cout << "DONE" << endl;

   cout << endl << "Looping through the input EventNtuple/MicroNtuple ... " << endl;
   int nMatchedToWH = 0, nMatchedToZH = 0, nMatchedToTTH = 0, nUnmatched = 0, nMultiMatch = 0;
   vector<pair<int,int> > MultiMatchSampleEvents, UnmatchedSampleEvents;
   for(int ientry=0; ientry<nEventNtuple; ientry++) {
      t1->GetEntry(ientry);
      loadbar(ientry+1,nEventNtuple);

      int nW = 0, nZ = 0, nT = 0, nH = 0;
      for (unsigned int i=0; i<ntuple->genParticleCollection.size() && i< 30; i++) {
         if(abs(ntuple->genParticleCollection[i].pdgId)==24) {
            nW++;
            for (int j=0; j<ntuple->genParticleCollection[i].numberOfMothers; j++) {
               if (abs(ntuple->genParticleCollection[ntuple->genParticleCollection[i].motherPositions[j]].pdgId)==6)
                  nW--;
            }
         }
         if(abs(ntuple->genParticleCollection[i].pdgId)==23)
            nZ++;
         if(abs(ntuple->genParticleCollection[i].pdgId)==6)
            nT++;
         if(abs(ntuple->genParticleCollection[i].pdgId)==25)
            nH++;
      }

      if(hTo.CompareTo("WW")==0) {
         if(nW >= 3 && nZ == 0 && nT < 2 && nH == 1) {
            WHEventNtuple->Fill();
            nMatchedToWH++;
         }
         else if(nW == 2 && nZ == 1 && nT == 0 && nH == 1) {
            ZHEventNtuple->Fill();
            nMatchedToZH++;
         }
         else if(nW == 2 && nZ == 0 && nT == 2 && nH == 1) {
            TTHEventNtuple->Fill();
            nMatchedToTTH++;
         }
         else if((nW > 2 && nZ > 0 && nT == 0 && nH > 0) ||
                 (nW > 2 && nZ == 0 && nT > 0 && nH > 0) ||
                 (nW == 2 && nZ > 0 && nT > 0 && nH > 0)){
            nMultiMatch++;
            MultiMatchSampleEvents.push_back(make_pair(ientry,ntuple->event));
         }
         else {
            nUnmatched++;
            UnmatchedSampleEvents.push_back(make_pair(ientry,ntuple->event));
         }
      }
      else if(hTo.CompareTo("ZZ")==0) {
         if(ntuple->event == 201693) {
            ZHEventNtuple->Fill();
            nMatchedToZH++;
         }
         else if(nW == 1 && nZ == 2 && nT < 2 && nH == 1) {
            WHEventNtuple->Fill();
            nMatchedToWH++;
         }
         else if(nW == 0 && nZ >= 3 && nT == 0 && nH == 1) {
            ZHEventNtuple->Fill();
            nMatchedToZH++;
         }
         else if(nW == 0 && nZ == 2 && nT == 2 && nH == 1) {
            TTHEventNtuple->Fill();
            nMatchedToTTH++;
         }
         else if((nW > 0 && nZ > 2 && nT == 0 && nH > 0) ||
                 (nW > 0 && nZ == 2 && nT > 0 && nH > 0) ||
                 (nW == 0 && nZ > 2 && nT > 0 && nH > 0)){
            nMultiMatch++;
            MultiMatchSampleEvents.push_back(make_pair(ientry,ntuple->event));
         }
         else {
            nUnmatched++;
            UnmatchedSampleEvents.push_back(make_pair(ientry,ntuple->event));
         }
      }
   }

   cout << endl << "Printing the matching information ... " << endl;
   cout << "\t" << setw(12) << "Sample"     << setw(12) << "nMatched"    << endl
        << "\t" << setw(12) << "------"     << setw(12) << "--------"    << endl
        << "\t" << setw(12) << "Original"   << setw(12) << nEventNtuple  << endl
        << "\t" << setw(12) << "------"     << setw(12) << "--------"    << endl
        << "\t" << setw(12) << "WH"         << setw(12) << nMatchedToWH  << endl
        << "\t" << setw(12) << "ZH"         << setw(12) << nMatchedToZH  << endl
        << "\t" << setw(12) << "TTH"        << setw(12) << nMatchedToTTH << endl
        << "\t" << setw(12) << "MultiMatch" << setw(12) << nMultiMatch   << endl
        << "\t" << setw(12) << "Unmatched"  << setw(12) << nUnmatched    << endl
        << "\t" << setw(12) << "------"     << setw(12) << "--------"    << endl
        << "\t" << setw(12) << "Total"      << setw(12) << nMatchedToWH+nMatchedToZH+nMatchedToTTH+nMultiMatch+nUnmatched << endl;

   cout << endl << "If MultiMatch or Unmatched events, print the first 10 sample (entry,event) for each ... " << endl;
   cout << "\tMultiMatch: ";
   for(unsigned int i=0; i<MultiMatchSampleEvents.size() && i<10; i++) {
      cout << "(" << MultiMatchSampleEvents[i].first << "," << MultiMatchSampleEvents[i].second << ")";
      if(i==9 || i==MultiMatchSampleEvents.size()-1)
         cout << endl;
      else
         cout << ", ";
   }
   if(MultiMatchSampleEvents.size()==0)
      cout << endl;
   cout << "\tUnmatched: ";
   for(unsigned int i=0; i<UnmatchedSampleEvents.size() && i<10; i++) {
      cout << "(" << UnmatchedSampleEvents[i].first << "," << UnmatchedSampleEvents[i].second << ")";
      if(i==9 || i==UnmatchedSampleEvents.size()-1)
         cout << endl;
      else
         cout << ", ";
   }
   if(UnmatchedSampleEvents.size()==0)
      cout << endl;

   cout << "Write the output TTrees in the appropriate files ... " << flush;
   if(eventNtupleOnly)
      WHFile.cd("PS");
   WHEventNtuple->AutoSave();
   WHEventNtuple->Write();
   if(eventNtupleOnly)
      ZHFile.cd("PS");
   ZHEventNtuple->AutoSave();
   ZHEventNtuple->Write();
   if(eventNtupleOnly)
      TTHFile.cd("PS");
   TTHEventNtuple->AutoSave();
   TTHEventNtuple->Write();
   cout << "DONE" << endl;

   cout << "Close all the files ... " << flush;
   ifile->Close();
   WHFile.Close();
   ZHFile.Close();
   TTHFile.Close();
   cout << "DONE" << endl;

}



