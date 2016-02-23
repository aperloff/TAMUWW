#include "TROOT.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TH1.h"
#include "TH1D.h"
#include "TError.h"
#include "TString.h"
#include "TAxis.h"

#include <iostream>
#include <string>
#include <vector>

#include "TAMUWW/Tools/interface/Style.hh"

using namespace std;

void MakeZJetsToLLPtWeightFile() {
   // Trying to speed up the code
   gEnv->SetValue("TFile.AsyncPrefetching", 1);

   Style* st = new Style();
   st->setTDRStyle();
   //TGaxis::SetMaxDigits(3);

   cout << "Opening output file ... " << flush;
   TFile* ofile = TFile::Open("ZllPtWeights.root","RECREATE");
   ofile->mkdir("validation");
   cout << "DONE" << endl;
   /*
   const int nJetBins = 3;
   const int nLeptonBins = 2;
   string jetBins[nJetBins] = {"jets2","jets3","jets4"};
   string leptonBins[nLeptonBins] = {"electron","muon"};
   */
   const int nJetBins = 1;
   const int nLeptonBins = 1;
   string jetBins[nJetBins] = {"jets2"};
   string leptonBins[nLeptonBins] = {"electron"};

   string object = "LeptPt_";

   for(int jbin=0; jbin<nJetBins; jbin++) {
      for(int lbin=0; lbin<nLeptonBins; lbin++) {
         cout << "Opening input file ... " << flush;
         TFile* ifile = TFile::Open((jetBins[jbin]+"/"+leptonBins[lbin]+"/noMETCut/histos_"+leptonBins[lbin]+"_"+jetBins[jbin]+".root").c_str(),"READ");
         cout << "DONE" << endl;

         cout << "Getting input histograms ... " << flush;
         string hname = object+"ZJetsToLL_M50_"+leptonBins[lbin];
         TH1D* ZJetsToLL = (TH1D*)ifile->Get(hname.c_str())->Clone("ZJetsToLL");
         hname = object+"ZJetsToLL_M10To50_"+leptonBins[lbin];
         ZJetsToLL->Add((TH1D*)ifile->Get(hname.c_str()));
         hname = object+"WJets_"+leptonBins[lbin];
         TH1D* WJets = (TH1D*)ifile->Get(hname.c_str())->Clone("WJets");
         cout << "DONE" << endl;

         cout << "Scaling ZJetToLL to WJets ... " << flush;
         ZJetsToLL->Scale(WJets->Integral()/ZJetsToLL->Integral());
         cout << "DONE" << endl;

         cout << "Making weights ... " << flush;
         hname = "ZllWeight_"+jetBins[jbin]+"_"+leptonBins[lbin];
         TH1D* ZllWeight = (TH1D*)WJets->Clone(hname.c_str());
         ZllWeight->Divide(ZJetsToLL);
         cout << "DONE" << endl;

         cout << "Writing weights ... " << flush;
         ofile->cd();
         ZllWeight->Write();
         cout << "DONE" << endl;

         cout << "Making and writing a validation canvas ... " << flush;
         ofile->cd("validation");
         string cname = "WeightValidation_"+jetBins[jbin]+"_"+leptonBins[lbin];
         TH1D* hup = new TH1D("hup","hup",10000,0,1000);
         hup->GetXaxis()->SetLimits(20,1000);
         hup->GetYaxis()->SetRangeUser(0.2,ZJetsToLL->GetMaximum()*100.0);
         TH1D* hdown = new TH1D("hdown","hdown",10000,0,1000);
         hdown->GetXaxis()->SetLimits(20,1000);
         hdown->GetXaxis()->SetTitle("p_{T}^{lepton}");
         hdown->GetYaxis()->SetRangeUser(0.1,10);
         hup->GetYaxis()->SetTitle("Events");
         hdown->GetYaxis()->SetTitle("Weight");
         TCanvas* c = st->tdrDiCanvas(cname.c_str(),hup,hdown);
         WJets->SetLineColor(kGreen);
         ZJetsToLL->SetLineColor(kRed);
         ZllWeight->SetLineColor(kBlue);
         TVirtualPad* p1 = c->cd(1);
         p1->SetLogx();
         p1->SetLogy();
         st->tdrDraw(WJets,"",kFullCircle,kGreen);
         st->tdrDraw(ZJetsToLL,"",kFullCircle,kRed);
         TVirtualPad* p2 = c->cd(2);
         p2->SetLogx();
         p2->SetLogy();
         st->tdrDraw(ZllWeight,"",kFullCircle,kBlue);
         c->Write();
         cout << "DONE" << endl;

         cout << "Closing the ifile ... " << flush;
         ifile->Close();
         cout << "DONE" << endl;
      }
   }
   ofile->Close();
}
