//Our libraries
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"
#include "TAMUWW/SpecialTools/interface/Table.hh"
#include "TAMUWW/SpecialTools/interface/TableRow.hh"
#include "TAMUWW/SpecialTools/interface/TableCellVal.hh"
#include "TAMUWW/SpecialTools/interface/TableCellText.hh"
#include "TAMUWW/SpecialTools/interface/Value.hh"
#include "TAMUWW/SpecialTools/interface/ProgressBar.hh"
#include "TAMUWW/MEPATNtuple/interface/METree.hh"
#include "TAMUWW/MatrixElement/interface/EventProbDefs.hh"

// ROOT libraries
#include "TROOT.h"
#include "TSystem.h"
#include "TBenchmark.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2D.h"
#include "TChain.h"

// C++ libraries
#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <utility>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//  Local Functions
////////////////////////////////////////////////////////////////////////////////

/// prints out the average CPU times and real times for each of the ME processes
void average_times_per_process(string ifile);

/// prints out the average CPU times and real times for each of the ME processes,
///  but uses the logs stored in the ROOT files instead of the text file logs
void average_times_per_process_root(string ifile);

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

/*
Examples:
   txt/stdout file: readMELog_x -ifile CondorME_WlnuJets_M-10To50_lowerSubLeadingJet_HighPtLeptonKept050mu061el_13413158_45.stdout -mcSample WlnuJets_M-10To50_lowerSubLeadingJet_HighPtLeptonKept050mu061el
   root file: readMELog_x -ifile root://cmseos.fnal.gov//store/user/eusebi/Winter12to13ME8TeV/rootOutput/WlnuJets_M-10To50_lowerSubLeadingJet_HighPtLeptonKept050mu061el/WlnuJets_M-10To50_lowerSubLeadingJet_HighPtLeptonKept050mu061el45.root -mcSample WlnuJets_M-10To50_lowerSubLeadingJet_HighPtLeptonKept050mu061el
*/

//______________________________________________________________________________
int main(int argc,char**argv) {
 
   // evaluate command-line / configuration file options
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;
   
   string ifile    = cl.getValue<string> ("ifile");
   string mcSample = cl.getValue<string> ("mcSample", "");
   
   if (!cl.check()) return 0;
   cl.print();
   
   // benchmark this executable
   TBenchmark* m_benchmark = new TBenchmark();
   m_benchmark->Reset();
   m_benchmark->Start("event");

   // Analyze the log file
   cout << "readMELog::Average times for " << mcSample << " sample:" << endl;
   if(ifile.find(".txt")!=string::npos || ifile.find(".stdout")!=string::npos)
      average_times_per_process(ifile);
   else if(ifile.find(".root")!=string::npos)
      average_times_per_process_root(ifile);
   else {
      cout << "ERROR::readMELog_x Can't process the ME logs because the input file type is not a known type [txt or root]" << endl;
      std::terminate();
   }

   
   // benchmark this executable
   m_benchmark->Stop("event");
   cout << "readMELog_x" << endl << "\tCPU time = " << m_benchmark->GetCpuTime("event") << " s" << endl
        << "\tReal time = " << m_benchmark->GetRealTime("event") << " s" << endl;
   delete m_benchmark;

   return 0;

}//main

////////////////////////////////////////////////////////////////////////////////
//  Implement Local Functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void average_times_per_process(string ifile) {

   // open the log file
   fstream fin;
   fin.open(ifile,fstream::in);

   map<string,vector<double> > CPUTimes;
   map<string,vector<double> > realTimes;
   char line[1024];
   string currentME;
   double currentCPUTime = 0;
   double currentRealTime = 0;
   size_t CPUPos;
   size_t realPos;
   size_t sPos;
   int nevents;

   // clear header
   //cout << endl << endl << "First while" << endl << endl;
   //fin.getline(line,1024);
   while(TString(line).Contains("Starting")!=1) {
      fin.getline(line,1024);
      //cout << TString(line) << endl;
   }
   //cout << endl << endl << "First for" << endl << endl;
   //50 to first event
   //56 to first process
   for(int i=0; i<50; i++) {
      fin.getline(line,1024);
      //cout << TString(line) << endl;
   }

   fin.getline(line,1024);
   while(TString(line).Contains("Finished")!=1) {
      if(TString(line).Contains("Processing")) {
         nevents++;
         //Can't use a ProgressBar because we don't know how many total events there are
         if(string(line).substr(19,string(line).find("with")-19-1)=="1")
            cout << "\tEvent " << string(line).substr(19,string(line).find("with")-20) << "...";
         else
            cout << "DONE" << endl
                 << "\tEvent " << string(line).substr(19,string(line).find("with")-20) << "...";  
      }
      while(TString(line).Contains("Starting")!=1) {
         fin.getline(line,1024);
      }
      currentME = string(line).substr(11,string(line).size()-24);
      while(TString(line).Contains("CPU")!=1) {
         fin.getline(line,1024);
      }
      CPUPos = string(line).find("CPU");
      sPos = string(line).find("s");
      currentCPUTime = atof(string(line).substr(CPUPos+10,sPos-15).c_str());
      realPos = string(line).find("Real");
      sPos = string(line).rfind("s");
      currentRealTime = atof(string(line).substr(realPos+11,sPos-realPos-12).c_str());

      CPUTimes[currentME].push_back(currentCPUTime);
      realTimes[currentME].push_back(currentRealTime);
      fin.getline(line,1024);
   }
   cout << "DONE" << endl;

   double jobCPU = 0;
   double jobReal = 0;
   map<string,vector<double> >::iterator it_real=realTimes.begin();
   //cout << "sfsg1 " << (*it_real).second.size()<< endl;

   for(map<string,vector<double> >::iterator it_cpu=CPUTimes.begin(); it_cpu!=CPUTimes.end(); it_cpu++) {
      cout << endl << "\tME Process: " << (*it_cpu).first << endl;
      double averageCPU = 0;
      double averageReal = 0;
      double sumCPU = 0;
      double sumReal = 0;
      //cout << "\tAveraging CPU times " << endl;
      for(unsigned int i=0; i<(*it_cpu).second.size(); i++) {
         ProgressBar::loadbar3("\tAveraging CPU times",i+1,(*it_cpu).second.size());
         //if(i%10 == 0)
         //   cout << ".";
         sumCPU+=(*it_cpu).second[i];
         jobCPU+=(*it_cpu).second[i];
      }
      cout << endl;
      //cout << endl << "\tAveraging real times " << endl;
      for(unsigned int i=0; i<(*it_real).second.size(); i++) {
         ProgressBar::loadbar3("\tAveraging real times",i+1,(*it_real).second.size());
         //cout << ".";
         sumReal+=(*it_real).second[i];
         jobReal+=(*it_real).second[i];
      }
      averageCPU = sumCPU/(*it_cpu).second.size();
      averageReal = sumReal/(*it_real).second.size();

      
      cout << endl << "\t\tCPU: " << averageCPU << " s" << endl
           << "\t\tReal: " << averageReal << " s" << endl;

      it_real++;
   }

   cout << endl << "\t\tBatch Job: " << endl
        << "\t\tCPU: " << jobCPU/3600. << " hour(s)" << endl
        << "\t\tReal: " << jobReal/3600. << " hour(s)" << endl
        << "\t\t~CPU/event: " << jobCPU/nevents/60. << " min" << endl
        << "\t\t~Real/event: " << jobReal/nevents/60. << " min" << endl;

   // close the log file
   fin.close();

}//average_times_per_process

//______________________________________________________________________________
void average_times_per_process_root(string ifile) {
   TFile* inf = TFile::Open(ifile.c_str(),"READ");
   TTree* t = (TTree*)inf->Get("METree");
   METree* meNtuple = new METree();
   t->SetBranchAddress("METree",&meNtuple);
   int nentries = t->GetEntries();

   //int=event number from MC (not place in file)
   //pair<double,double> = CPU and real time
   map<int,pair<double,double> > sum_over_processes_per_event;

   Table sum_over_events_per_process("Averages Per Physics Process");
   TableRow* tableRow;
   pair<TableCellVal*,TableCellVal*> tableCellValues;

   string tmeTypeName;
   double nHWW=0, nWH=0;
   //nentries = 1;

   for (int ientry = 0; ientry < nentries; ++ientry){
      ProgressBar::loadbar3("Processing Events:",ientry+1,nentries);
      t->GetEntry(ientry);

      if(sum_over_processes_per_event.find(meNtuple->getEvent())==sum_over_processes_per_event.end()) {
         sum_over_processes_per_event[meNtuple->getEvent()] = make_pair(0,0);
      }

      for(int i = 0; i < meNtuple->getNProbStat(); ++i) {
         sum_over_processes_per_event[meNtuple->getEvent()].first+=meNtuple->getProbStat(i)->tCpuTime;
         sum_over_processes_per_event[meNtuple->getEvent()].second+=meNtuple->getProbStat(i)->tRealTime;

         tmeTypeName = DEFS::EP::getTypeString(meNtuple->getProbStat(i)->tmeType);
         if(ientry==0) {
            if(tmeTypeName=="HWW") nHWW++;
            else if(tmeTypeName=="WH") nWH++;

            if(tmeTypeName=="HWW" && sum_over_events_per_process(tmeTypeName,"CPU Time (s)")) continue;
            if(tmeTypeName=="WH"  && sum_over_events_per_process(tmeTypeName,"CPU Time (s)")) continue;
            if(tmeTypeName=="WLg" && sum_over_events_per_process(tmeTypeName,"CPU Time (s)")) continue;
            tableRow = new TableRow(tmeTypeName);
            tableCellValues.first = new TableCellVal("CPU Time (s)");
            tableCellValues.second = new TableCellVal("Real Time (s)");
            tableCellValues.first->val.value = 0;
            tableCellValues.second->val.value = 0;
            tableRow->addCellEntries(tableCellValues.first);
            tableRow->addCellEntries(tableCellValues.second);
            sum_over_events_per_process.addRow(*tableRow);
            delete tableRow;
         }
         dynamic_cast<TableCellVal*>(sum_over_events_per_process(tmeTypeName,"CPU Time (s)"))->val.value+=meNtuple->getProbStat(i)->tCpuTime;
         dynamic_cast<TableCellVal*>(sum_over_events_per_process(tmeTypeName,"Real Time (s)"))->val.value+=meNtuple->getProbStat(i)->tRealTime;
      }
   }
   //Divide sums by the number of entries
   bool HWW_div = false, WH_div = false, WLg_div = false;
   for(int i = 0; i < meNtuple->getNProbStat(); ++i) {
      tmeTypeName = DEFS::EP::getTypeString(meNtuple->getProbStat(i)->tmeType);
      if(tmeTypeName!="HWW" && tmeTypeName!="WH" && tmeTypeName!="WLg") {
         dynamic_cast<TableCellVal*>(sum_over_events_per_process(tmeTypeName,"CPU Time (s)"))->val.value/=nentries;
         dynamic_cast<TableCellVal*>(sum_over_events_per_process(tmeTypeName,"Real Time (s)"))->val.value/=nentries;
      }
      else if(tmeTypeName=="WLg" && !WLg_div) {
         dynamic_cast<TableCellVal*>(sum_over_events_per_process(tmeTypeName,"CPU Time (s)"))->val.value/=(nentries*2.);
         dynamic_cast<TableCellVal*>(sum_over_events_per_process(tmeTypeName,"Real Time (s)"))->val.value/=(nentries*2.);
         WLg_div = true;
      }
      else if(tmeTypeName=="HWW" && !HWW_div) {
         dynamic_cast<TableCellVal*>(sum_over_events_per_process(tmeTypeName,"CPU Time (s)"))->val.value/=(nentries*nHWW);
         dynamic_cast<TableCellVal*>(sum_over_events_per_process(tmeTypeName,"Real Time (s)"))->val.value/=(nentries*nHWW);
         HWW_div = true;
      }
      else if(tmeTypeName=="WH" && !WH_div) {
         dynamic_cast<TableCellVal*>(sum_over_events_per_process(tmeTypeName,"CPU Time (s)"))->val.value/=(nentries*nWH);
         dynamic_cast<TableCellVal*>(sum_over_events_per_process(tmeTypeName,"Real Time (s)"))->val.value/=(nentries*nWH);
         WH_div = true;
      }
   }
   cout << endl;

   //Total Time for the job
   double jobCPU = 0, jobReal = 0;
   for(map<int,pair<double,double> >::iterator it_evt=sum_over_processes_per_event.begin(); it_evt!=sum_over_processes_per_event.end(); it_evt++) {
      jobCPU+=it_evt->second.first;
      jobReal+=it_evt->second.second;
   }

   //
   // Print results
   //

   //Average over all events per process
   sum_over_events_per_process.printTable(cout);

   //Average time per event (all processes summed)
   cout << "Average per event: " << endl
        << "CPU/event: " << jobCPU/nentries/60. << " min" << endl
        << "Real/event: " << jobReal/nentries/60. << " min" << endl;

   cout << endl << "Batch Job: " << endl
        << "CPU: " << jobCPU/3600. << " hour(s)" << endl
        << "Real: " << jobReal/3600. << " hour(s)" << endl << endl;

}//average_times_per_process_root
