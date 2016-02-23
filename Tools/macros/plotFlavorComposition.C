#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveText.h"

#include <iostream>
#include <vector>
#include <utility>

using namespace std;

//
// Usage
//
//.L $CMSSW_BASE/src/TAMUWW/Tools/macros/plotFlavorComposition.C+
//plotFlavorCompositionLoop()

vector<TString> getWJetFlavors(bool doOther) {
   vector<TString> ret;
   ret.push_back("");
   ret.push_back("WpInclusive");
   ret.push_back("WpB");
   ret.push_back("WpC");
   ret.push_back("WpLight");
   if(doOther) ret.push_back("WpOther");
   return ret;
}

pair<TString,TH1D*> getDataMinusMCHistogram(TFile* inf, TString lepton) {
   TString short_lepton = (lepton=="electron") ? "El" : "Mu";
   TString data_name = Form("KinMEBDT_Single%s_Data_%s",short_lepton.Data(),lepton.Data());
   TH1D* data = (TH1D*)inf->Get(data_name);

   int nprocesses = 22;
   TString processes[] = {"WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar",
                          "TTbar","WW",Form("QCD_%sFULL",short_lepton.Data()),"ZJets","WJets","WH_HToZZ_M125","ZH_HToZZ_M125",
                          "TTH_HToZZ_M125","WH125_HToBB","TTH_HToBB_M125","ggH125","qqH125","WH_HToWW_M125",
                          "ZH_HToWW_M125","TTH_HToWW_M125"};
   vector<TH1D*> mc;
   TString hname;
   TH1D* htmp;
   for(int iproc=0; iproc<nprocesses; iproc++) {
      if(processes[iproc]!="WJets") {
         hname = Form("KinMEBDT_%s_%s",processes[iproc].Data(),lepton.Data());
         htmp = (TH1D*)inf->Get(hname);
         data->Add(htmp,-1.0);
      }
   }
   return make_pair(data_name,data);
}

void plotFlavorComposition(TString jetBin, TString lepton, bool scaleToInclusive = false, 
                           bool doLogY = false, bool doGrid = true, bool doOther = false,
                           bool doData = false) {
   cout << "Opening file ... " << flush;
   TFile* inf = TFile::Open(Form("%s/%s/histos_%s.root",jetBin.Data(),lepton.Data(),lepton.Data()),"READ");
   cout << "DONE" << endl;

   cout << "Getting WJets flavors ... " << flush;
   vector<TString> wjets_flavors = getWJetFlavors(doOther);
   cout << "DONE" << endl;

   cout << "Getting histograms from file ... " << flush;
   vector<pair<TString,TH1D*> > wjets_histograms;
   TString hname;
   TH1D* htmp;
   for(unsigned int iflavor=0; iflavor<wjets_flavors.size(); iflavor++) {
      hname = (wjets_flavors[iflavor]=="") ? Form("KinMEBDT_WJets_%s",lepton.Data())
                                           : Form("KinMEBDT_%s_WJets_%s",wjets_flavors[iflavor].Data(),lepton.Data());
      htmp = (TH1D*)inf->Get(hname);
      wjets_histograms.push_back(make_pair(hname,htmp));
   }
   cout << "DONE" << endl;

   if(doData) {
      cout << "Getting data-(MC!WJets) ... " << flush;
      wjets_histograms.push_back(getDataMinusMCHistogram(inf,lepton));
      wjets_flavors.push_back("data");
      cout << "DONE" << endl;
   }

   cout << "Checking that histograms were successfully retrieved ... " << endl;
   for(unsigned int iflavor=0; iflavor<wjets_histograms.size(); iflavor++) {
      cout << "\t" << wjets_histograms[iflavor].first << ": " << wjets_histograms[iflavor].second << endl;
   }
   cout << "DONE" << endl;

   cout << "Calculating ratios ... " << flush;
   vector<double> wjets_ratios;
   for(unsigned int iflavor=0; iflavor<wjets_histograms.size(); iflavor++) {
      wjets_ratios.push_back(wjets_histograms[iflavor].second->Integral()/wjets_histograms[0].second->Integral());
   }
   cout << "DONE" << endl;

   if(scaleToInclusive) {
      cout << "Scaling the histograms ... " << flush;
      for(unsigned int iflavor=0; iflavor<wjets_histograms.size(); iflavor++) {
         wjets_histograms[iflavor].second->Scale(1.0/wjets_ratios[iflavor]);
      }
      cout << "DONE" << endl;
   }

   cout << "Formatting the histograms ... " << flush;
   int colors[] = {kBlack, kBlack, kRed, kGreen, kBlue, kMagenta};
   for(unsigned int iflavor=0; iflavor<wjets_histograms.size(); iflavor++) {
      wjets_histograms[iflavor].second->SetMarkerColor(colors[iflavor]);
      wjets_histograms[iflavor].second->SetLineColor(colors[iflavor]);
      if(iflavor==0) {
         wjets_histograms[iflavor].second->SetLineStyle(2);
         wjets_histograms[iflavor].second->SetMarkerStyle(kOpenCircle);
         wjets_histograms[iflavor].second->GetYaxis()->SetRangeUser(0.1,130000.0);
         wjets_histograms[iflavor].second->GetYaxis()->SetTitle("events (normalized to inclusive)");
         wjets_histograms[iflavor].second->GetXaxis()->SetTitle("KinMEBDT");
      }
      if(wjets_histograms[iflavor].first.Contains("data",TString::kIgnoreCase)) {
         wjets_histograms[iflavor].second->SetMarkerStyle(kFullCircle);
      }
   }
   cout << "DONE" << endl;

   TCanvas* c1 = new TCanvas("c","c",600,600);

   cout << "Drawing the histograms ... " << flush;
   for(unsigned int iflavor=0; iflavor<wjets_histograms.size(); iflavor++) {
      if(iflavor==0)
         wjets_histograms[iflavor].second->Draw();
      else
         wjets_histograms[iflavor].second->Draw("same");
   }
   cout << "DONE" << endl;
   
   cout << "Formatting the canvas ... " << flush;
   if(doLogY) c1->SetLogy();
   if(doGrid) {
      gPad->SetGridx();
      gPad->SetGridy();
   }
   cout << "DONE" << endl;

   cout << "Creating the legend ... " << flush;
   TLegend * leg = new TLegend(0.7,0.7,1.0,1.0);
   leg->AddEntry(wjets_histograms[0].second, "Inclusive (original)", "lep");
   leg->AddEntry(wjets_histograms[1].second, "W+Inclusive", "lep");
   leg->AddEntry(wjets_histograms[2].second, "W+b", "lep");
   leg->AddEntry(wjets_histograms[3].second, "W+c (no b)", "lep");
   leg->AddEntry(wjets_histograms[4].second, "W+light (no b or c)", "lep");
   if(doOther) leg->AddEntry(wjets_histograms[5].second, "W+Other", "lep");
   if(doData) leg->AddEntry(wjets_histograms[5].second, "data-MC!WJets", "lep");
   leg->Draw("same");
   cout << "DONE" << endl;

   cout << "Creating the ratio text ... " << flush;
   TPaveText* text = new TPaveText(0.7,0.4,1.0,0.7,"brNDC");
   text->AddText("Ratio To Inclusive");
   text->AddText(Form("Inclusive: %e",wjets_ratios[1]));
   text->AddText(Form("W+b: %e",wjets_ratios[2]));
   text->AddText(Form("W+c: %e",wjets_ratios[3]));
   text->AddText(Form("W+light: %e",wjets_ratios[4]));
   if(doOther) text->AddText(Form("W+Other: %e",wjets_ratios[5]));
   if(doData) text->AddText(Form("data-MC: %e",wjets_ratios[5]));
   text->Draw("same");
   cout << "DONE" << endl;

   cout << "Saving the canvases ... " << flush;
   TString save_name = "FlavorComposition_"+jetBin+"_"+lepton;
   if(doLogY) save_name+="_LogY";
   if(scaleToInclusive) save_name+="_ScaledToInclusive";
   c1->SaveAs(jetBin+"/"+lepton+"/"+save_name+".root");
   c1->SaveAs(jetBin+"/"+lepton+"/"+save_name+".png");
   cout << "DONE" << endl;

   //cout << "Closing input file ... " << flush;
   //inf->Close();
   //cout << "DONE" << endl;
}

void plotFlavorCompositionLoop() {
   vector<TString> jetBins;
   jetBins.push_back("Jets2");
   jetBins.push_back("Jets3");
   jetBins.push_back("Jets4");

   vector<TString> leptons;
   leptons.push_back("electron");
   leptons.push_back("muon");

   for(unsigned int iJetBin=0; iJetBin<jetBins.size(); iJetBin++) {
      for(unsigned int iLepton=0; iLepton<leptons.size(); iLepton++) {
         cout << "DOING JETBIN=" << jetBins[iJetBin] << " AND LEPTON=" << leptons[iLepton] << " ... " << endl;
         plotFlavorComposition(jetBins[iJetBin],leptons[iLepton],false,false,true,false,true);
         plotFlavorComposition(jetBins[iJetBin],leptons[iLepton],false,true,true,false,true);
         plotFlavorComposition(jetBins[iJetBin],leptons[iLepton],true,false,true,false,true);
      }
   }
}
