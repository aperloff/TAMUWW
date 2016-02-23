///////////////////////////////////////////////////////////////////
//
// MakeMaxEventProbTable_x
// -----------------------
//
//            03/18/2014 Alexx Perloff  <aperloff@physics.tamu.edu>
///////////////////////////////////////////////////////////////////

// C++ libraries
#include <cmath>
#include <iostream>
#include <iomanip>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <sstream>
#include <fstream>

// ROOT libraries
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TBranch.h"
#include "TString.h"
#include "TBenchmark.h"
#include "TTreeIndex.h"

// Our libraries
#include "TAMUWW/MEPATNtuple/interface/METree.hh"
#include "TAMUWW/MEPATNtuple/interface/MicroNtuple.hh"
#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"
#include "TAMUWW/SpecialTools/interface/Table.hh"
#include "TAMUWW/SpecialTools/interface/TableRow.hh"
#include "TAMUWW/SpecialTools/interface/TableCellInt.hh"
#include "TAMUWW/SpecialTools/interface/TableCellVal.hh"
#include "TAMUWW/SpecialTools/interface/TableCellText.hh"
#include "TAMUWW/SpecialTools/interface/Value.hh"
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"
#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/SpecialTools/interface/PhysicsProcess.hh"
#include "TAMUWW/SpecialTools/interface/ProgressBar.hh"
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"

using namespace std;

class eventProbInfo {
//private:
public:
   DEFS::EP::Type meType;
   DEFS::PhysicsProcess::Type ppType;
   double tmeParam;
   int probStatIndex;
   int nentries;
public :
   eventProbInfo() : meType(DEFS::EP::UNKNOWN), ppType(DEFS::PhysicsProcess::UNKNOWN), tmeParam(-1.0), probStatIndex(-1), nentries(0) {}
   eventProbInfo(DEFS::EP::Type meType_, DEFS::PhysicsProcess::Type ppType_, double tmeParam_, int probStatIndex_, int nentries_=0) :
      meType(meType_), ppType(ppType_), tmeParam(tmeParam_), probStatIndex(probStatIndex_), nentries(nentries_) {}

   bool operator<(const eventProbInfo& rhs) const {
      //if (meType < rhs.meType) {cout << endl << "op1 lt " << (int)meType-(int)rhs.meType << endl; return true;}
      //if (ppType < rhs.ppType) {cout << endl << "op2 lt " << (int)ppType-(int)rhs.ppType << endl; return true;}
      //if (tmeParam < rhs.tmeParam) {cout << endl << "op3 lt " << tmeParam-rhs.tmeParam << endl; return true;}
      double resLHS = 1.0e9*probStatIndex + 1000000.0*(int)meType + 1000.0*(int)ppType + 1.0*tmeParam;
      double resRHS = 1.0e9*rhs.probStatIndex + 1000000.0*(int)rhs.meType + 1000.0*(int)rhs.ppType + 1.0*rhs.tmeParam;
      //if(resLHS<resRHS) {cout << endl << "op lt " << resLHS-resRHS << endl; return true;}
      if(resLHS<resRHS) return true;
      else return false;
      //else {cout << endl << "gtet" << endl; return false;}
   }
   bool operator==(const eventProbInfo& rhs) const {
      if (meType == rhs.meType && ppType == rhs.ppType && tmeParam == rhs.tmeParam && probStatIndex == rhs.probStatIndex)
         return true;
      else
         return false;
   }
};
typedef map<eventProbInfo, Value> eventProbMap;
typedef map<eventProbInfo, vector<double> > eventProbMapVec;

void checkAndSetEventProb(int probStatIndex, DEFS::PhysicsProcess::Type ppType, DEFS::EP::Type meType,
                          double tmeParam, double eventProb, double eventProbError, eventProbMap& epm) {
   if(epm.find(eventProbInfo(meType,ppType,tmeParam,probStatIndex))!=epm.end()) {
   //if(epm.count(eventProbInfo(meType,ppType,tmeParam))>0) {
      if(epm[eventProbInfo(meType,ppType,tmeParam,probStatIndex)].value < eventProb) {
        epm[eventProbInfo(meType,ppType,tmeParam,probStatIndex)] = Value(eventProb,eventProbError);
      }
   }
   else {
      epm[eventProbInfo(meType,ppType,tmeParam,probStatIndex)] = Value(eventProb,eventProbError);
   }
}

void addEventProb(int probStatIndex, DEFS::PhysicsProcess::Type ppType, DEFS::EP::Type meType,
                  double tmeParam, double eventProb, double eventProbError, eventProbMap& epm, int nentries) {
   if(epm.find(eventProbInfo(meType,ppType,tmeParam,probStatIndex))!=epm.end()) {
      epm[eventProbInfo(meType,ppType,tmeParam,probStatIndex)] = Value(epm[eventProbInfo(meType,ppType,tmeParam,probStatIndex)].value+eventProb,
                                                                       epm[eventProbInfo(meType,ppType,tmeParam,probStatIndex)].error+TMath::Power(eventProbError,2));
   }
   else {
      eventProbInfo first(meType,ppType,tmeParam,probStatIndex,nentries);
      epm[first] = Value(eventProb,eventProbError);
   }
}

void storeEventProb(int probStatIndex, DEFS::PhysicsProcess::Type ppType, DEFS::EP::Type meType,
                    double tmeParam, double eventProb, double eventProbError, eventProbMapVec& epmv, int nentries) {
   if(epmv.find(eventProbInfo(meType,ppType,tmeParam,probStatIndex))!=epmv.end()) {
      epmv[eventProbInfo(meType,ppType,tmeParam,probStatIndex)].push_back(eventProb);
   }
   else {
      eventProbInfo first(meType,ppType,tmeParam,probStatIndex,nentries);
      epmv[first].push_back(eventProb);
   }
}

Table* makeTable(vector<eventProbMap>& epm) {

   Table* table = new Table("EventProbs");
   TableRow* tableRow;
   TableCellVal* tableCellVal;
   TableCellText* tableCellText;
   Value val;
   string text;

   eventProbMap::iterator it2=epm[1].begin();
   eventProbMap::iterator it3=epm[2].begin();
   for (eventProbMap::iterator it=epm[0].begin(); it!=epm[0].end(); ++it) {

      std::string s;
      std::stringstream ss;
      ss << it->first.probStatIndex;
      s = ss.str();
      ss.str("");
      tableRow = new TableRow(s);

      for (Int_t col=0; col<6; col++) {
         string title;
         if (col==0) {
            title = string("MatrixElementType");
            tableCellText = new TableCellText(title);
            tableCellText->text = DEFS::EP::getTypeString(it->first.meType);
            tableRow->addCellEntries(tableCellText);
         }
         else if (col==1) {
            title = string("PhysicsProcessType");
            tableCellText = new TableCellText(title);
            tableCellText->text = DEFS::PhysicsProcess::getTypeString(it->first.ppType);
            tableRow->addCellEntries(tableCellText);
         }
         else if (col==2) {
            title = string("tmeParam");
            tableCellVal = new TableCellVal(title);
            tableCellVal->val = Value(it->first.tmeParam,0.0);
            tableRow->addCellEntries(tableCellVal);
         }
         else if (col==3){
            title = string("MaxEventProb");
            tableCellVal = new TableCellVal(title);
            tableCellVal->val = it->second;
            tableRow->addCellEntries(tableCellVal);
         }
         else if (col==4){
            title = string("MeanEventProb");
            tableCellVal = new TableCellVal(title);
            tableCellVal->val = it2->second;
            tableRow->addCellEntries(tableCellVal);
         }
         else{
            title = string("MedianEventProb");
            tableCellVal = new TableCellVal(title);
            tableCellVal->val = it3->second;
            tableRow->addCellEntries(tableCellVal);
         }
      }
      table->addRow(*tableRow);
      delete tableRow;
      it2++;
      it3++;
   }
   return table;
}

////////////////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
int main(int argc, char**argv) {
   //
   // evaluate command-line / configuration file options
   //
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;
   string          ofile     = cl.getValue<string>       ("ofile",         "MaxMeanMedianEventProbs.txt");
   //string          lepCat    = cl.getValue<string>       ("lep",              "both");
   //DEFS::LeptonCat leptonCat = DEFS::getLeptonCat(lepCat);
   string          jetBinS   = cl.getValue<string>       ("jetBin",           "jets2");
   DEFS::JetBin    jetBin    = DEFS::getJetBin(jetBinS);
   string          tagcatS   = cl.getValue<string>       ("tagcat",           "pretag");
   DEFS::TagCat    tagcat    = DEFS::getTagCat(tagcatS);
   bool            debug     = cl.getValue<bool>         ("debug",            false);
   bool            batch     = cl.getValue<bool>         ("batch",            false);

   if (!cl.check()) return 0;
   cl.print();

   TBenchmark* m_benchmark = new TBenchmark();
   m_benchmark->Reset();
   m_benchmark->Start("event");

   vector<PhysicsProcess*> processes = DefaultValues::getProcessesHiggs(jetBin, tagcat, false, true, false, DEFS::MicroNtuple);
   if(debug) {
      processes.erase(processes.begin(),processes.begin()+11);
   }

   METree      * meNtuple    = 0;
   MicroNtuple * microNtuple = 0;

   eventProbMap    maxEventProbs;
   eventProbMap    meanEventProbs;
   eventProbMap    medianEventProbs;
   eventProbMapVec medianEventProbsStore;

   for(unsigned int iprocess=0; iprocess<processes.size(); iprocess++) {

      TString processName = processes[iprocess]->name;

      //Skip Data
      if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::SingleEl_Data) {
         continue;
      }
      else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::SingleMu_Data) {
         continue;
      }
      //Don't have any qqH or STopTW matrix elements
      else if(processName.Contains("qqH")) {
         continue;
      }
      else if(processName.Contains("STopTW")) {
         continue;
      }

      cout << endl << "MakeMaxEventProbTable::Doing physics process " << processName << endl;
      TChain* c = processes[iprocess]->chain;
      c->SetBranchStatus("lLV*",0);
      c->SetBranchStatus("jLV*",0);
      c->SetBranchStatus("vLV*",0);
      c->SetBranchStatus("METLV*",0);
      c->SetBranchStatus("genParticleCollection*",0);
      c->SetBranchStatus("epd*",0);

      if (c->GetBranch("METree")) {
         meNtuple = new METree();
         c->SetBranchAddress("METree",&meNtuple);
      }
      else {
         cout << "\tMakeMaxEventProbTable::METree branch not found." << endl
              << "\tSetting METree pointer to 0x0." << endl;
      }
      if (c->GetBranch("mnt")) {
         microNtuple = new MicroNtuple((int)jetBin);
         c->SetBranchAddress("mnt",&microNtuple);
      }
      else {
         cout << "\tMakeMaxEventProbTable::mnt branch not found." << endl
              << "\tSetting MicroNtuple pointer to 0x0." << endl;
      }

      unsigned int nentries = c->GetEntries();
      if(debug) nentries = 100;
      for (unsigned int ientry = 0 ; ientry < nentries ; ientry++) {
         if(!debug)
            ProgressBar::loadbar2(ientry+1,nentries);
         c->GetEntry(ientry);

         for (int i = 0; i < meNtuple->getNProbStat(); ++i){
            if(DEFS::PhysicsProcess::isHiggs(DEFS::PhysicsProcess::getProcessType(string(processName)))) {
               if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::ggH125 &&
                  meNtuple->getProbStat(i)->tmeType == DEFS::EP::HWW &&
                  meNtuple->getProbStat(i)->tmeParam == 125.0) {
                  checkAndSetEventProb(i,DEFS::PhysicsProcess::ggH125, DEFS::EP::HWW, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
                  addEventProb(i,DEFS::PhysicsProcess::ggH125, DEFS::EP::HWW, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
                  storeEventProb(i,DEFS::PhysicsProcess::ggH125, DEFS::EP::HWW, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
               }
               else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::WH125 &&
                  meNtuple->getProbStat(i)->tmeType == DEFS::EP::WH &&
                  meNtuple->getProbStat(i)->tmeParam == 125.0) {
                  checkAndSetEventProb(i,DEFS::PhysicsProcess::WH125, DEFS::EP::WH, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
                  addEventProb(i,DEFS::PhysicsProcess::WH125, DEFS::EP::WH, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
                  storeEventProb(i,DEFS::PhysicsProcess::WH125, DEFS::EP::WH, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
               }
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::WW &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::WW) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::WW, DEFS::EP::WW, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::WW, DEFS::EP::WW, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::WW, DEFS::EP::WW, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::WZ &&
                    (meNtuple->getProbStat(i)->tmeType == DEFS::EP::WZ ||
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::WZbb)) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::WZ, meNtuple->getProbStat(i)->tmeType, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::WZ, meNtuple->getProbStat(i)->tmeType, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::WZ, meNtuple->getProbStat(i)->tmeType, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::ZZ &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::ZZ) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::ZZ, DEFS::EP::ZZ, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::ZZ, DEFS::EP::ZZ, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::ZZ, DEFS::EP::ZZ, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::ZJets &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::ZLight) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::ZJets, DEFS::EP::ZLight, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::ZJets, DEFS::EP::ZLight, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::ZJets, DEFS::EP::ZLight, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::WJets &&
                    (meNtuple->getProbStat(i)->tmeType == DEFS::EP::WJets ||
                     meNtuple->getProbStat(i)->tmeType == DEFS::EP::WLg ||
                     meNtuple->getProbStat(i)->tmeType == DEFS::EP::Wgg || 
                     meNtuple->getProbStat(i)->tmeType == DEFS::EP::WLL ||
                     meNtuple->getProbStat(i)->tmeType == DEFS::EP::WLb ||
                     meNtuple->getProbStat(i)->tmeType == DEFS::EP::Wbb ||
                     meNtuple->getProbStat(i)->tmeType == DEFS::EP::Wc  ||
                     meNtuple->getProbStat(i)->tmeType == DEFS::EP::WLight)) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::WJets, meNtuple->getProbStat(i)->tmeType, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::WJets, meNtuple->getProbStat(i)->tmeType, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::WJets, meNtuple->getProbStat(i)->tmeType, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::TTbar &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::TTbar) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::TTbar, DEFS::EP::TTbar, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::TTbar, DEFS::EP::TTbar, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::TTbar, DEFS::EP::TTbar, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::STopS_T &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::STopS) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::STopS_T, DEFS::EP::STopS, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::STopS_T, DEFS::EP::STopS, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::STopS_T, DEFS::EP::STopS, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::STopS_Tbar &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::STopS) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::STopS_Tbar, DEFS::EP::STopS, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::STopS_Tbar, DEFS::EP::STopS, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::STopS_Tbar, DEFS::EP::STopS, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::STopT_T &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::STopT) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::STopT_T, DEFS::EP::STopT, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::STopT_T, DEFS::EP::STopT, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::STopT_T, DEFS::EP::STopT, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::STopT_Tbar &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::STopT) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::STopT_Tbar, DEFS::EP::STopT, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::STopT_Tbar, DEFS::EP::STopT, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::STopT_Tbar, DEFS::EP::STopT, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::STopTW_T &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::STopTW) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::STopTW_T, DEFS::EP::STopTW, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::STopTW_T, DEFS::EP::STopTW, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::STopTW_T, DEFS::EP::STopTW, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::STopTW_Tbar &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::STopTW) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::STopTW_Tbar, DEFS::EP::STopTW, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::STopTW_Tbar, DEFS::EP::STopTW, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::STopTW_Tbar, DEFS::EP::STopTW, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::QCD_ElEnriched &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::QCD) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::QCD_ElEnriched, DEFS::EP::QCD, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::QCD_ElEnriched, DEFS::EP::QCD, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::QCD_ElEnriched, DEFS::EP::QCD, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::QCD_ElFULL &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::QCD) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::QCD_ElFULL, DEFS::EP::QCD, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::QCD_ElFULL, DEFS::EP::QCD, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::QCD_ElFULL, DEFS::EP::QCD, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::QCD_MuFULL &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::QCD) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::QCD_MuFULL, DEFS::EP::QCD, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::QCD_MuFULL, DEFS::EP::QCD, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::QCD_MuFULL, DEFS::EP::QCD, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
            else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::QCD_MuEnriched &&
                    meNtuple->getProbStat(i)->tmeType == DEFS::EP::QCD) {
               checkAndSetEventProb(i,DEFS::PhysicsProcess::QCD_MuEnriched, DEFS::EP::QCD, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, maxEventProbs);
               addEventProb(i,DEFS::PhysicsProcess::QCD_MuEnriched, DEFS::EP::QCD, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, meanEventProbs, nentries);
               storeEventProb(i,DEFS::PhysicsProcess::QCD_MuEnriched, DEFS::EP::QCD, meNtuple->getProbStat(i)->tmeParam, meNtuple->getProbStat(i)->tEventProb, meNtuple->getProbStat(i)->tEventErr, medianEventProbsStore, nentries);
            }
         }
      }
      delete meNtuple;
      delete microNtuple;
   }

   //Divide the sum of the eventprobs by the number of entries to get the mean
   for (eventProbMap::iterator it=meanEventProbs.begin(); it!=meanEventProbs.end(); ++it) {
      it->second.value /= it->first.nentries;
      it->second.error = TMath::Sqrt(it->second.error)/it->first.nentries;
   }

   //Find the median event probs
   for (eventProbMapVec::iterator it=medianEventProbsStore.begin(); it!=medianEventProbsStore.end(); ++it) {
      medianEventProbs[it->first] = Value(TMath::Median(it->first.nentries,&(it->second.at(0))),0.0);
   }

   //Make a table out of the eventProbMap
   vector<eventProbMap> epm;
   epm.push_back(maxEventProbs);
   epm.push_back(meanEventProbs);
   epm.push_back(medianEventProbs);
   Table* table = makeTable(epm);

   //Print the map to a file
   if (batch) {
      table->printToFile("./"+ofile);
   }
   else {
      table->printToFile(DefaultValues::getConfigPath()+ofile);
   }

   m_benchmark->Stop("event");
   cout << endl << "MakeMaxEventProbTable_x" << endl << "\tCPU time = " << m_benchmark->GetCpuTime("event") << " s" << endl
        << "\tReal time = " << m_benchmark->GetRealTime("event") << " s" << endl;
   delete m_benchmark;

   return 0;
}
