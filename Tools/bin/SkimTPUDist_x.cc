//Our libraries
#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"
#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"
#include "TAMUWW/SpecialTools/interface/Value.hh"
#include "TAMUWW/SpecialTools/interface/Table.hh"
#include "TAMUWW/SpecialTools/interface/TableRow.hh"
#include "TAMUWW/SpecialTools/interface/TableCellVal.hh"
#include "TAMUWW/SpecialTools/interface/TableCellText.hh"
#include "TAMUWW/SpecialTools/interface/FileLocationTable.hh"

#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"

// ROOT libraries
#include "TROOT.h"
#include "TSystem.h"
#include "TBenchmark.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"
#include "TMath.h"
#include "TError.h"

// C++ libraries
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <stdlib.h>

using namespace std;


////////////////////////////////////////////////////////////////////////////////
//  Local Functions
////////////////////////////////////////////////////////////////////////////////

///  fills the histograms and controls the output canvas and file for the rest of the program
vector<TH1D*> skimHistograms(vector<PhysicsProcess*> procs, Table& fileTable, string histogramName, bool verbose);

///  saves the resultant histograms to a user defined file
void saveHistograms(TString ofilename, vector<TH1D*> hists);

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
int main(int argc,char**argv) {
 
   // evaluate command-line / configuration file options
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;

   string  histogramName = cl.getValue<string>  ("histogramName",               "PS/TPUDist");
   bool    verbose       = cl.getValue<bool>    ("verbose",                            false);
   TString ofilename     = cl.getValue<TString> ("ofilename",        "TPUDistributions.root");

   if (!cl.check()) return 0;
   cl.print();

   TBenchmark* m_benchmark = new TBenchmark();
   m_benchmark->Reset();
   m_benchmark->Start("event");

   // The vector holding all processes.
   vector <PhysicsProcess*> procs = DefaultValues::getProcessesHiggs(DEFS::jets2, DEFS::pretag, true, true, false, DEFS::EventNtuple);
      
   Table fileTable = DefaultValues::getFileLocationTable(DEFS::pretag);

   vector<TH1D*> hists = skimHistograms(procs, fileTable, histogramName, verbose);

   saveHistograms(DefaultValues::getConfigPath()+ofilename, hists);

   m_benchmark->Stop("event");
   cout << "SkimTPUDist_x" << endl << "\tCPU time = " << m_benchmark->GetCpuTime("event") << " s" << endl
        << "\tReal time = " << m_benchmark->GetRealTime("event") << " s" << endl;
   delete m_benchmark;

   return 0;

}//main


////////////////////////////////////////////////////////////////////////////////
//  Implement Local Functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
vector<TH1D*> skimHistograms(vector<PhysicsProcess*> procs, Table& fileTable, string histogramName, bool verbose) {

   vector<TH1D*> ret;

   for(unsigned int p=0; p<procs.size(); p++) {
      
      cout << "Doing Process " << procs[p]->name << " ... " << endl;

      // find the file location for that process
      TableCellText * cellFile = (TableCellText *) fileTable.getCellRowColumn(string(procs[p]->name),"FilePath_"+DEFS::getNtupleTypeString(DEFS::EventNtuple));
      
      // make sure we found the cell
      if (cellFile == 0){
         cout<<"ERROR SkimTPUDist_x::main Table "<<fileTable.getTableOrigin()
             <<" does not have row "<<procs[p]->name
             <<" for column FilePath_" << DEFS::getNtupleTypeString(DEFS::EventNtuple) <<endl;
         cout<<" SKIPPING PROCESS "<<procs[p]->name<<endl;
         return ret;
      }

      TFile * file = TFile::Open(TString(cellFile->text),"READ");
      if (!file->IsOpen())
      {
         cout << "ERROR SkimTPUDist_x::skimHistograms() could not open file " << cellFile->text << endl;
         return ret;
      }

      /*
   if (!file->cd("PS"))
   {
      cout << "ERROR PhysicsProcess::PhysicsProcess() could not CD into directory PS in file " << fileName << endl;
      return;
   }
      */
      TH1D* temp;
      if(procs[p]->name.Contains("QCD")) { //&& procs[p]->name.Contains("Enriched")) {
         temp = new TH1D(Form("TPUDist_%s",procs[p]->name.Data()),Form("TPUDist_%s (jet1+jets2p)",procs[p]->name.Data()),600,0,60);
         if (!file->cd("PS"))
         {
            cout << "ERROR SkimTPUDist_x::skimHistograms() could not CD into directory PS in file " << cellFile->text << endl;
            assert(file->cd("PS"));
         }
         TTree* jet1 = (TTree*)gDirectory->Get("jet1");
         TTree* jets2p = (TTree*)gDirectory->Get("jets2p");
         if(jet1)
            jet1->Draw("vLV[0].npv>>j1(600,0,60)","","goff");
         else {
            cout << "ERROR SkimTPUDist_x::skimHistograms() TTree jet1 is NULL" << endl;
            assert(jet1 != 0);
         }
         if(jets2p)
            jets2p->Draw("vLV[0].npv>>j2p(600,0,60)","","goff");
         else {
            cout << "ERROR SkimTPUDist_x::skimHistograms() TTree jets2p is NULL" << endl;
            assert(jets2p != 0);
         }
         TH1D* j1 = (TH1D*) gDirectory->Get("j1");
         TH1D* j2p = (TH1D*) gDirectory->Get("j2p");
         temp->Add(j1);
         temp->Add(j2p);
      }
      else{
         temp = (TH1D*) file->Get(histogramName.c_str());
      }

      if (temp == 0)
      {
         cout << "ERROR SkimTPUDist_x::skimHistograms() could not find histogram named " << histogramName << " in file " << cellFile->text << endl;
         if(procs[p]->name.Contains("QCD")!=0 && procs[p]->name.Contains("Data")!=0)
            assert(temp != 0);
      }
      else {
         ret.push_back((TH1D*)temp->Clone(TString("TPUDist_")+procs[p]->name));
         if(ret.back() == 0) {
            cout << "WARNING SkimTPUDist_x::skimHistograms() the clone of histogram named " << histogramName << " in file " << cellFile->text << " was a NULL pointer." << endl;
            ret.pop_back();
         }
         else {
            ret.back()->SetDirectory(0);
         }
      }

      //delete temp;
      file->Close();
      //delete file;
      
   }

   return ret;

}//skimHistograms

//______________________________________________________________________________
void saveHistograms(TString ofilename, vector<TH1D*> hists) {

   TFile * file = new TFile(ofilename,"RECREATE");
   if (!file->IsOpen())
   {
      cout << "ERROR SkimTPUDist_x::saveHistograms() could not open file " << ofilename << endl;
      return;
   }

   for(unsigned int h=0; h<hists.size(); h++) {
      hists[h]->Write();
   }

   file->Close();


}//saveHistograms


//Method to save QCD histos. To be implemented.
//QCD_ElEnriched
/*
TFile* fQCD = new TFile("/uscms/home/aperloff/nobackup/Summer12ME8TeV/MEInput/QCD_Electron_Dp6p7.root","READ")
TFile* fTPU = new TFile("TPUDistributions.root","UPDATE")
fQCD->cd("PS")
TH1D* TPUDist_QCD_ElEnriched = new TH1D("TPUDist_QCD_ElEnriched","TPUDist_QCD_ElEnriched",600,0,60)
jets2p->Draw("vLV[0].npv>>TPUDist_QCD_ElEnriched")
fTPU->cd()
TPUDist_QCD_ElEnriched->Write()
.q
*/
//QCD_ElFULL
/*
TFile* fQCD = new TFile("/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/SingleEl_Full.root","READ")
TFile* fTPU = new TFile("TPUDistributions.root","UPDATE")
fQCD->cd("PS")
TH1D* temp = new TH1D("TPUDist_QCD_ElFULL","TPUDist_QCD_ElFULL (jet1+jets2p)",600,0,60)
TH1D* j1 = new TH1D("j1","pileup_QCD_ElFULL jet1",600,0,60)
TH1D* j2p = new TH1D("j2p","pileup_QCD_ElFULL jets2p",600,0,60)
TTree* j1t = (TTree*)gDirectory->Get("jet1")
TTree* j2pt = (TTree*)gDirectory->Get("jets2p")
j1t->Draw("vLV[0].npv>>j1")
j2pt->Draw("vLV[0].npv>>j2p")
temp->Add(j1)
temp->Add(j2p)
fTPU->cd()
fTPU->Write()
.q
*/
//Save jet1 + jet2p
/*
TFile* _file0 = new TFile("pileup12_noTrig.root","UPDATE")
TFile* fQCD = new TFile("/uscms/home/aperloff/nobackup/Summer12ME8TeV/MEInput/SingleEl_Full.root","READ")
fQCD->cd("PS")
TH1D* temp = new TH1D("pileup_QCD","pileup_QCD (jet1+jets2p)",600,0,60)
TH1D* j1 = new TH1D("j1","pileup_QCD jet1",600,0,60)
TH1D* j2p = new TH1D("j2p","pileup_QCD jets2p",600,0,60)
TTree* j1t = (TTree*)gDirectory->Get("jet1")
TTree* j2pt = (TTree*)gDirectory->Get("jets2p")
j1t->Draw("vLV[0].npv>>j1")
j2pt->Draw("vLV[0].npv>>j2p")
temp->Add(j1)
temp->Add(j2p)
_file0->cd()
temp->Write()
_file0->Write()
 */
//Split the WH ZH and TTH contributions (i.e. copy same distribution 3 times
/*
TFile* fTPU = new TFile("TPUDistributions.root","UPDATE")
TH1D* orig = (TH1D*) fTPU->Get("TPUDist_WH_ZH_TTH_HToWW_M125")
TH1D* TPUDist_WH_HToWW_M125 = (TH1D*)orig->Clone("TPUDist_WH_HToWW_M125")
TH1D* TPUDist_ZH_HToWW_M125 = (TH1D*)orig->Clone("TPUDist_ZH_HToWW_M125")
TH1D* TPUDist_TTH_HToWW_M125 = (TH1D*)orig->Clone("TPUDist_TTH_HToWW_M125")
TPUDist_WH_HToWW_M125->Write()
TPUDist_ZH_HToWW_M125->Write()
TPUDist_TTH_HToWW_M125->Write()
TH1D* orig2 = (TH1D*) fTPU->Get("TPUDist_WH_ZH_TTH_HToZZ_M125")
TH1D* TPUDist_WH_HToZZ_M125 = (TH1D*)orig2->Clone("TPUDist_WH_HToZZ_M125")
TH1D* TPUDist_ZH_HToZZ_M125 = (TH1D*)orig2->Clone("TPUDist_ZH_HToZZ_M125")
TH1D* TPUDist_TTH_HToZZ_M125 = (TH1D*)orig2->Clone("TPUDist_TTH_HToZZ_M125")
TPUDist_WH_HToZZ_M125->Write()
TPUDist_ZH_HToZZ_M125->Write()
TPUDist_TTH_HToZZ_M125->Write()
fTPU->Close()
*/
//Add TPUDist to the split WH ZH and TTH files
/*
TFile* fMerged = new TFile("/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/WH_ZH_TTH_HToZZ_M125.root","READ")
TFile* fWH = new TFile("/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/WH_HToZZ_M125.root","UPDATE")
TFile* fZH = new TFile("/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/ZH_HToZZ_M125.root","UPDATE")
TFile* fTTH = new TFile("/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/TTH_HToZZ_M125.root","UPDATE")
fMerged->cd("PS")
TH1D* orig = (TH1D*) gDirectory->Get("TPUDist")
fWH->cd("PS")
orig->Write()
fZH->cd("PS")
orig->Write()
fTTH->cd("PS")
orig->Write()
fWH->Close()
fZH->Close()
fTTH->Close()
fMerged->Close()
*/
//Add npv histograms from data samples to the data pileup distributions
/*
TFile* fpu = new TFile("pileup12_noTrig_minBiasXsec69400_coarseBinning_withAdditionalNPVHist.root","UPDATE")
TFile* fel = new TFile("/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/SingleEl_Data_19p148fb.root","READ")
TFile* fmu = new TFile("/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/SingleMu_Data_19p279fb.root","READ")
fel->cd("PS")
TH1D* hel = new TH1D("pileup_SingleEl","pileup_SingleEl",60,0,60)
jets2p->Draw("vLV[0].npv>>pileup_SingleEl")
fmu->cd("PS")
TH1D* hmu = new TH1D("pileup_SingleMu","pileup_SingleMu",60,0,60)
jets2p->Draw("vLV[0].npv>>pileup_SingleMu")
fpu->cd()
hel->Write()
hmu->Write()
fel->Close()
fmu->Close()
fpu->Close()
*/
