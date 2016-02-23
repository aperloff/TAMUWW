//ROOT libraries
#include "TROOT.h"
#include "TSystem.h"
#include "TBenchmark.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"

//Standard libraries
#include <iostream>
#include <vector>
#include <string>
#include <map>

//Our libraries
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"
#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
int main(int argc,char**argv) {
   // evaluate command-line / configuration file options
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;

   TString             inputElectronChannel = cl.getValue<TString> ("inputElectronChannel",  "histos_electron.root");
   TString             inputMuonChannel     = cl.getValue<TString> ("inputMuonChannel",          "histos_muon.root");
   TString             outputFilename       = cl.getValue<TString> ("outputFilename",          "./histos_both.root");
   TString             distribution         = cl.getValue<TString> ("distribution",                           "MET");
   string              lepCat               = cl.getValue<string>  ("lep",                                   "both");
   DEFS::LeptonCat     leptonCat            = DEFS::getLeptonCat(lepCat);
   string              ntype                = cl.getValue<string>  ("ntype",                          "MicroNtuple");
   DEFS::NtupleType    ntupleType           = DEFS::getNtupleType(ntype);
   string              jBin                 = cl.getValue<string>  ("jBin",                                 "jets2");
   DEFS::JetBin        jetBin               = DEFS::getJetBin(jBin);
   string              tcat                 = cl.getValue<string>  ("tcat",                                "pretag");
   DEFS::TagCat        tagCat               = DEFS::getTagCat(tcat);

   if (!cl.check()) return 0;
   cl.print();
   
   TBenchmark* m_benchmark = new TBenchmark();
   m_benchmark->Reset();
   m_benchmark->Start("event");

   // Check that the leptcat actually exists
   if (leptonCat == DEFS::none) {
      cout << "ERROR::combineChannels_x called with unknown lepton category " << lepCat << endl;
      return 1;
   }

   //if(!outputPath.EndsWith("/")) outputPath += "/";

   cout << "Getting processes ... ";
   vector <PhysicsProcess*> procs = DefaultValues::getProcessesHiggs(jetBin, tagCat, true, true, true, ntupleType, leptonCat);
   cout << "DONE" << endl;

   cout << "Opening the histogram files ... ";
   TFile* histosElectron = new TFile(inputElectronChannel,"READ");
   TFile* histosMuon     = new TFile(inputMuonChannel,"READ");
   cout << "DONE" << endl;

   map<TString,TH1D*> outputHistograms;

   cout << "Looping through the processes ... " << endl;
   for(unsigned int iproc=0; iproc<procs.size(); iproc++) {
      cout << "\tDoing process " << procs[iproc]->name << " (" << iproc+1 << "/" << procs.size() << ")... ";
      TString distNameEle = Form("%s_%s_electron",distribution.Data(),procs[iproc]->name.Data());
      TString distNameMu  = Form("%s_%s_muon",distribution.Data(),procs[iproc]->name.Data());
      TString distNameBoth;
      if(procs[iproc]->name.Contains("QCD",TString::kIgnoreCase)) {
         distNameBoth = Form("%s_QCD_both",distribution.Data());
      }
      else if(procs[iproc]->name.Contains("Data",TString::kIgnoreCase)) {
         distNameBoth = Form("%s_Data_both",distribution.Data());
      }
      else {
         distNameBoth = Form("%s_%s_both",distribution.Data(),procs[iproc]->name.Data());
      }
      double intEle = 0, intMu = 0, intBoth = 0;
      TH1D* tempHist;

      histosElectron->cd();
      tempHist = (TH1D*)gDirectory->Get(distNameEle);
      if(tempHist) {
         outputHistograms[distNameBoth] = (TH1D*)(tempHist)->Clone(distNameBoth);
         intEle = outputHistograms[distNameBoth]->Integral();
         tempHist = 0;
      }
      histosMuon->cd();
      tempHist = (TH1D*)gDirectory->Get(distNameMu);
      if(tempHist) {
         if(outputHistograms.find(distNameBoth)!=outputHistograms.end())
            outputHistograms[distNameBoth]->Add(tempHist,1.0);
         else
            outputHistograms[distNameBoth] = (TH1D*)(tempHist)->Clone(distNameBoth);
         intMu = ((TH1D*)gDirectory->Get(distNameMu))->Integral();
      }
      intBoth = outputHistograms[distNameBoth]->Integral();

      if((intEle+intMu)-intBoth > 0.0001) {
         cout << endl << "\t\tWARNING::The sum of the integrals of the electron (" << distNameEle << " = " << intEle
              << ") and muon (" << distNameMu << " = " << intMu
              << ") channels does not equal the integral of the combined channel (" << distNameBoth << " = " << intBoth
              << ")." << endl;
      }
      cout << "DONE" << endl;
   }

   cout << "Opening the output file ... ";
   TFile* ofile = new TFile(outputFilename,"RECREATE");
   cout << "DONE" << endl;

   cout << "Writing the histograms to the output file ... ";
   for (map<TString,TH1D*>::iterator it=outputHistograms.begin(); it!=outputHistograms.end(); ++it) {
      it->second->Write();
   }
   ofile->Close();
   cout << "DONE" << endl;

   m_benchmark->Stop("event");
   cout << "combineChannels_x" << endl << "\tCPU time = " << m_benchmark->GetCpuTime("event") << " s" << endl
        << "\tReal time = " << m_benchmark->GetRealTime("event") << " s" << endl;
   delete m_benchmark;

   return 0;

}//main
