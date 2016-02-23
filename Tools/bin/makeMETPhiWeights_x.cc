///////////////////////////////////////////////////////////////////////////////
//
// makeMETPhiWeights_x
// -------------------
//
//               09/09/2013 Alexx Perloff            <aperloff@physics.tamu.edu>
////////////////////////////////////////////////////////////////////////////////

#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TLine.h"
#include "TString.h"
#include "TPaveText.h"
#include "TLatex.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// define local functions
////////////////////////////////////////////////////////////////////////////////
void getHistosFromPlotterOutput(TFile* ifile, TH2D *& hd, TH2D *& hmc, bool full);

////////////////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
int main(int argc,char**argv) {
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;

   string lepton    = cl.getValue<string> ("lepton",          "");
   string ifilePath = cl.getValue<string> ("ifilePath",     "./");
   string ifileName = cl.getValue<string> ("ifile",           "");
   string ofilePath = cl.getValue<string> ("ofilePath",     "./");
   string ofileName = cl.getValue<string> ("ofile",           "");
   string objects   = cl.getValue<string> ("object",          "");
   bool   full      = cl.getValue<bool>   ("full",          true);
   
   if (!cl.check()) return 0;
   cl.print();

   if(!TString(ifilePath).EndsWith("/")) ifilePath+="/";
   if(!TString(ofilePath).EndsWith("/")) ofilePath+="/";
   if(!gSystem->OpenDirectory(TString(ofilePath))) gSystem->mkdir(TString(ofilePath));

   // Check that lepton is set to either muon or electron
   if(lepton != "muon" && lepton != "electron")
   {
      cout << "makeMETPhiWeights::ERROR Lepton was not properly set. Options are electron and muon." << endl;
      return 1;
   }

   if(ifileName == "")
      ifileName = "histos_" + lepton + ".root";
   
   if(ofileName == "")
      ofileName = "METPhiWeights_" + lepton + ".root";
   
   // Open the input file and get all histos
   TFile* ifile = new TFile((ifilePath+ifileName).c_str(),"READ");
   if (!ifile->IsOpen()) { cout<<"Can't open "<<ifilePath<<ifileName<<endl; return 0; }
   
   TH1D* hd, hmc;
   getHistosFromPlotterOutput(ifile,hd,hmc,full);
      
   // set the styles
   hd->UseCurrentStyle();
   hmc->UseCurrentStyle();
   
   TH1D* hdw = new TH1D("hdw","MET Phi Weights (Data)",hd->GetNbinsX(),hd->GetXaxis()->GetXmin(),hd->GetXaxis()->GetXmax());
   TH1D* hmcw = new TH1D("hmcw","MET Phi Weights (MC)",hmc->GetNbinsX(),hmc->GetXaxis()->GetXmin(),hmc->GetXaxis()->GetXmax());
   

   TFile* ofile = new TFile((outFileLoc+outFile).c_str(),"RECREATE");
   if (!ofile->IsOpen()) { cout<<"Can't create "<<ofilePath<<ofileName<<endl; return 0; }

   return 0;
}


////////////////////////////////////////////////////////////////////////////////
// implement local functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void getHistosFromPlotterOutput(TFile* ifile, TH1D *& hd, TH1D  *& hmc, bool full){
   bool removesw2 = false;// true;

   // list of background processes. WJets and QCD are handled separately
   TString hn = "METPhi_";

   // Name of the data histo
   TString data = hn+"SingleEl_Data_electron"; // 19148 pb-1

   // Vector containing all other histos. Proper normalization
   // to the luminosity of the signal data is assumed
   vector<TString> backs;
   backs.push_back(hn+"STopS_T_electron");
   backs.push_back(hn+"STopS_Tbar_electron");
   backs.push_back(hn+"STopT_T_electron");
   backs.push_back(hn+"STopT_Tbar_electron");
   backs.push_back(hn+"STopTW_T_electron");
   backs.push_back(hn+"STopTW_Tbar_electron");
   backs.push_back(hn+"TTbar_electron");
   backs.push_back(hn+"ZJets_electron");
   backs.push_back(hn+"WW_electron");
   backs.push_back(hn+"WZ_electron");
   backs.push_back(hn+"ggH125_electron");
   backs.push_back(hn+"qqH125_electron");
   backs.push_back(hn+"WH125_electron");
   backs.push_back(hn+"WJets_electron");
   // Name of the QCD histo
   TString qcd;
   if(full)
      qcd = hn+"QCD_ElFULL_electron"; // 19148 pb-1
   else
      qcd = hn+"QCD_ElEnriched_electron"; // 1169 pb-1
   backs.push_back(qcd);

   // Data
   hd =  (TH2D*) _file0->Get(data);
   if (removesw2) hd->GetSumw2()->Set(0);
   //hd->Rebin2D(rebinEta, rebinMet);

   //WJets
   hw = (TH2D*) _file0->Get(wj);
   if (removesw2) hw->GetSumw2()->Set(0);
   //hw->Rebin2D(rebinEta, rebinMet);

   //QCD
   hq = (TH2D*) _file0->Get(qcd);
   if (removesw2) hq->GetSumw2()->Set(0);
   //hq->Rebin2D(rebinEta, rebinMet);

   // Remainding Processes
   // Substract all other backgrounds from the data
   hr = 0;
   for (unsigned int b = 0 ; b < backs.size();b++){
      TH2D * aux = (TH2D*) _file0->Get(backs[b]);
      //aux->Rebin2D(rebinEta,rebinMet);
      if (hr==0) hr = aux;
      else hr->Add(aux);
   }//for

}// getHistosFromPlotterOutput
