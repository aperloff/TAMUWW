///////////////////////////////////////////////////////////////////
//
// MakeMedianPurityTable_x
// -----------------------
//
//            09/19/2017 Alexx Perloff  <aperloff@physics.tamu.edu>
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

class BDTInfo {
public:
   DEFS::JetBin jetBin;
   DEFS::LeptonCat leptonCat;
   string BDTType;
   int nentries;
public :
   BDTInfo() : jetBin(DEFS::JetBin::jets0), leptonCat(DEFS::LeptonCat::none), BDTType("KinBDT"), nentries(0) {}
   BDTInfo(DEFS::JetBin jetBin_, DEFS::LeptonCat leptonCat_, string BDTType_, int nentries_=0) :
      jetBin(jetBin_), leptonCat(leptonCat_), BDTType(BDTType_), nentries(nentries_) {}

   bool operator<(const BDTInfo& rhs) const {
      double resLHS = 1.0e3*(int)jetBin + 1.0e0*(int)leptonCat;
      double resRHS = 1.0e3*(int)rhs.jetBin + 1.0e0*(int)rhs.leptonCat;
      if(resLHS<resRHS) return true;
      else return false;
   }
   bool operator==(const BDTInfo& rhs) const {
      if (leptonCat == rhs.leptonCat && jetBin == rhs.jetBin && BDTType == rhs.BDTType)
         return true;
      else
         return false;
   }
};
typedef map<BDTInfo, Value> BDTMap;
typedef map<BDTInfo, vector<double> > BDTMapVec;

void addBDT(DEFS::JetBin jetBin, DEFS::LeptonCat leptonCat, string BDTType, 
            double BDT, double BDTError, BDTMap& bmap, int nentries) {
   if(bmap.find(BDTInfo(jetBin,leptonCat,BDTType))!=bmap.end()) {
      bmap[BDTInfo(jetBin,leptonCat,BDTType)] = Value(bmap[BDTInfo(jetBin,leptonCat,BDTType)].value+BDT,
                                                      bmap[BDTInfo(jetBin,leptonCat,BDTType)].error+TMath::Power(BDTError,2));
   }
   else {
      BDTInfo first(jetBin,leptonCat,BDTType,nentries);
      bmap[first] = Value(BDT,BDTError);
   }
}

void storeBDT(DEFS::JetBin jetBin, DEFS::LeptonCat leptonCat, string BDTType,
              double BDT, double BDTError, BDTMapVec& bmapv, int nentries) {
   if(bmapv.find(BDTInfo(jetBin,leptonCat,BDTType))!=bmapv.end()) {
      bmapv[BDTInfo(jetBin,leptonCat,BDTType)].push_back(BDT);
   }
   else {
      BDTInfo first(jetBin,leptonCat,BDTType,nentries);
      bmapv[first].push_back(BDT);
   }
}

Table* makeTable(vector<BDTMap>& bmap) {

   Table* table = new Table("JetBin");
   TableRow* tableRow;
   TableCellVal* tableCellVal;
   TableCellText* tableCellText;
   Value val;
   string text;

   BDTMap::iterator it2=bmap[1].begin();
   for (BDTMap::iterator it=bmap[0].begin(); it!=bmap[0].end(); ++it) {

      std::string s;
      std::stringstream ss;
      ss << DEFS::getJetBinString(it->first.jetBin);
      s = ss.str();
      ss.str("");
      tableRow = new TableRow(s);

      string title;

      title = string("LeptonCat");
      tableCellText = new TableCellText(title);
      tableCellText->text = DEFS::getLeptonCatString(it->first.leptonCat);
      tableRow->addCellEntries(tableCellText);

      title = string("BDTType");
      tableCellText = new TableCellText(title);
      tableCellText->text = it->first.BDTType;
      tableRow->addCellEntries(tableCellText);

      title = string("Mean");
      tableCellVal = new TableCellVal(title);
      tableCellVal->val = it->second;
      tableRow->addCellEntries(tableCellVal);

      title = string("Median");
      tableCellVal = new TableCellVal(title);
      tableCellVal->val = it2->second;
      tableRow->addCellEntries(tableCellVal);

      title = string("Comment");
      tableCellText = new TableCellText(title);
      tableCellText->text = "Found using all of the Higgs signal ntuples";
      tableRow->addCellEntries(tableCellText);

      table->addRow(*tableRow);
      delete tableRow;
      it2++;
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
   string          ofile     = cl.getValue<string>       ("ofile",         "MedianPurity.txt");
   string          tagcatS   = cl.getValue<string>       ("tagcat",        "pretag");
   DEFS::TagCat    tagcat    = DEFS::getTagCat(tagcatS);
   bool            debug     = cl.getValue<bool>         ("debug",         false);
   bool            batch     = cl.getValue<bool>         ("batch",         false);

   if (!cl.check()) return 0;
   cl.print();

   TBenchmark* m_benchmark = new TBenchmark();
   m_benchmark->Reset();
   m_benchmark->Start("event");

   vector<PhysicsProcess*> processes = DefaultValues::getProcessesHiggs(DEFS::JetBin::jets2, tagcat, false, false, false, DEFS::MicroNtuple);
   for(unsigned int i=0; i<processes.size(); i++) {
        if(!DEFS::PhysicsProcess::isHiggs(DEFS::PhysicsProcess::getProcessType(string(processes[i]->name)))) {
           cout << "Erasing process " << processes[i]->name << endl;
           processes.erase(processes.begin()+i);
           i--;
        }
    }
    cout << "Processes left to run on:" << endl;
   	for(auto proc : processes) {
   		cout << "\t" << proc->name << endl;
   	}

   EventNtuple * ntuple      = 0;
   MicroNtuple * microNtuple = 0;

   BDTMap    meanBDT;
   BDTMap    medianBDT;
   BDTMapVec medianBDTStore;

   for(unsigned int iprocess=0; iprocess<processes.size(); iprocess++) {

      TString processName = processes[iprocess]->name;

      cout << endl << "MakeMedianPurityTable::Doing physics process " << processName << endl;
      TChain* c = processes[iprocess]->chain;
      c->SetBranchStatus("*",0);
      c->SetBranchStatus("KinBDT*",1);
      c->SetBranchStatus("lLV.leptonCat",1);
      c->SetBranchStatus("jLV",1);

      if (c->GetBranch("EvtTree")) {
         ntuple = new EventNtuple();
         c->SetBranchAddress("EvtTree",&ntuple);
      }
      else {
         cout << "\tMakeMedianPurityTable::EvtTree branch not found." << endl
              << "\tSetting ntuple pointer to 0x0." << endl;
      }
      if (c->GetBranch("mnt")) {
         microNtuple = new MicroNtuple(2);
         c->SetBranchAddress("mnt",&microNtuple);
      }
      else {
         cout << "\tMakeMedianPurityTable::mnt branch not found." << endl
              << "\tSetting MicroNtuple pointer to 0x0." << endl;
      }

      unsigned int nentries = c->GetEntries();
      if(debug) nentries = 100;
      DEFS::LeptonCat leptonCat;
      DEFS::JetBin    jetBin;
      for (unsigned int ientry = 0 ; ientry < nentries ; ientry++) {
         if(!debug)
            ProgressBar::loadbar2(ientry+1,nentries);
         c->GetEntry(ientry);

        jetBin = (ntuple->jLV.size()>=4) ? DEFS::JetBin::jets4 : (ntuple->jLV.size()==3) ? DEFS::JetBin::jets3 :
        		 (ntuple->jLV.size()==2) ? DEFS::JetBin::jets2 : DEFS::JetBin::jets0;
        if((int)jetBin<2) continue;
        leptonCat = ntuple->lLV[0].leptonCat;
        if(leptonCat!=DEFS::LeptonCat::electron && leptonCat!=DEFS::LeptonCat::muon) continue;
        if(microNtuple->KinBDT<-10) continue;

        addBDT(jetBin,leptonCat,"KinBDT",microNtuple->KinBDT,0.0,meanBDT,nentries);
        storeBDT(jetBin,leptonCat,"KinBDT",microNtuple->KinBDT,0.0,medianBDTStore,nentries);
      }
      delete microNtuple;
   }

   //Divide the sum of the eventprobs by the number of entries to get the mean
   for (BDTMap::iterator it=meanBDT.begin(); it!=meanBDT.end(); ++it) {
      it->second.value /= it->first.nentries;
      it->second.error = TMath::Sqrt(it->second.error)/it->first.nentries;
   }

   //Find the median event probs
   for (BDTMapVec::iterator it=medianBDTStore.begin(); it!=medianBDTStore.end(); ++it) {
      medianBDT[it->first] = Value(TMath::Median(it->first.nentries,&(it->second.at(0))),0.0);
   }

   //Make a table out of the eventProbMap
   vector<BDTMap> bmap;
   bmap.push_back(meanBDT);
   bmap.push_back(medianBDT);
   Table* table = makeTable(bmap);

   //Print the map to a file
   if (batch) {
      table->printToFile("./"+ofile);
   }
   else {
      table->printToFile(DefaultValues::getConfigPath()+ofile);
   }

   m_benchmark->Stop("event");
   cout << endl << "MakeMedianPurityTable_x" << endl << "\tCPU time = " << m_benchmark->GetCpuTime("event") << " s" << endl
        << "\tReal time = " << m_benchmark->GetRealTime("event") << " s" << endl;
   delete m_benchmark;

   return 0;
}