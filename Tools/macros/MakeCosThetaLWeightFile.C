#include "TROOT.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TError.h"
#include "TString.h"

#include <iostream>
#include <vector>

#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"

using namespace std;

void MakeCosThetaLWeightFile(TString lepton) {
   // Trying to speed up the code
   gEnv->SetValue("TFile.AsyncPrefetching", 1);

   Float_t cosdPhiWW = 0;
   Float_t cosdPhiWH = 0;
   Float_t costhetal = 0;
   Float_t costhetaj = 0;
   Float_t costhetaWH = 0;
   Float_t jacksonAngle = 0;
   int nevts = 0;
   int origIgnoreLevel = gErrorIgnoreLevel;
   gErrorIgnoreLevel = kBreak;
   int leptonNumber = 0;

   TFile* _file0;
   TFile* _file1;
   TTree* jet1Tree;
   EventNtuple* ntuple = new EventNtuple();
   TH1F* Data1j;
   TH1F* WJets1j;

   TFile* ofile = TFile::Open("CosThetaLWeight_jet1.root","UPDATE");
   if(lepton.CompareTo("electron")==0) {
      leptonNumber=2;
      Data1j = new TH1F("CosThetaL_Data_jet1_electron","CosThetaL_Data_jet1_electron",25,-1,1);
      WJets1j = new TH1F("CosThetaL_WJets_jet1_electron","CosThetaL_WJets_jet1_electron",25,-1,1);
      _file0 = new TFile("/eos/uscms/store/user/aperloff//MatrixElement/Summer12ME8TeV/MEInput/SingleEl_Data_19p148fb.root","READ");
   }
   else {
      leptonNumber = 1;
      Data1j = new TH1F("CosThetaL_Data_jet1_muon","CosThetaL_Data_jet1_muon",25,-1,1);
      WJets1j = new TH1F("CosThetaL_WJets_jet1_muon","CosThetaL_WJets_jet1_muon",25,-1,1);
      _file0 = new TFile("/eos/uscms/store/user/aperloff//MatrixElement/Summer12ME8TeV/MEInput/SingleMu_Data_19p279fb.root","READ");
   }
   _file1 = new TFile("/eos/uscms/store/user/aperloff//MatrixElement/Summer12ME8TeV/MEInput/WJets.root","READ");

   ofile->cd();

   jet1Tree = (TTree*)_file0->Get("PS/jet1");
   jet1Tree->SetBranchAddress("EvtNtuple",&ntuple);
   nevts = 1000000;//jet1Tree->GetEntries();
   for(int i=0; i<nevts; i++) {
      if(i%10000==0){
         cout<<"Doing event "<<i<<"/"<<nevts<<endl;
      }
      jet1Tree->GetEntry(i);
      if(ntuple->jLV.size()==1){
         if(ntuple->jLV[0].Pt()<=0||ntuple->lLV[0].Pt()<=0) {
            cout << "WARNING!!!" << endl;
            continue;
         }
         ntuple->getAngularVariables(cosdPhiWW,cosdPhiWH,costhetal,costhetaj,costhetaWH,jacksonAngle);
         if(!TMath::IsNaN(costhetal)) {
            Data1j->Fill(costhetal);
         }
         else{
            cout << "WARNING::CosTheta_l is a NaN" << endl;
         }
      }
   }
   //jet1->Draw("CosTheta_l>>CosThetaL_Data_jet1_electron","@jLV.size()==1","goff")

   jet1Tree = (TTree*)_file1->Get("PS/jet1");
   jet1Tree->SetBranchAddress("EvtNtuple",&ntuple);
   nevts = 1000000;//jet1Tree->GetEntries();
   for(int i=0; i<nevts; i++) {
      if(i%10000==0){
         cout<<"Doing event "<<i<<"/"<<nevts<<endl;
      }
      jet1Tree->GetEntry(i);
      if(ntuple->jLV.size()==1&&ntuple->lLV[0].leptonCat==leptonNumber){
         if(ntuple->jLV[0].Pt()<=0||ntuple->lLV[0].Pt()<=0) {
            cout << "WARNING!!!" << endl;
            continue;
         }
         ntuple->getAngularVariables(cosdPhiWW,cosdPhiWH,costhetal,costhetaj,costhetaWH,jacksonAngle);
         if(!TMath::IsNaN(costhetal)) {
            WJets1j->Fill(costhetal);
         }
         else{
            cout << "WARNING::CosTheta_l is a NaN" << endl;
         }
      }
   }
   //jet1->Draw("CosTheta_l>>CosThetaL_WJets_jet1_electron","leptonCat==2&&@jLV.size()==1","goff")
   gErrorIgnoreLevel = origIgnoreLevel;

   Data1j->Scale(1.0/Data1j->Integral());
   WJets1j->Scale(1.0/WJets1j->Integral());
   TH1F* weights1j = (TH1F*)Data1j->Clone(Form("CosThetaL_WeightHist_jet1_%s",lepton));
   weights1j->Divide(WJets1j);
   Data1j->SetLineColor(kBlack);
   WJets1j->SetLineColor(kGreen);
   weights1j->SetLineColor(kRed);
   TCanvas* c = new TCanvas(Form("CosThetaL_Canvas_jet1_%s",lepton),Form("CosThetaL_Canvas_jet1_%s",lepton),800,400);
   c->Divide(2,1);
   c->cd(1);
   Data1j->Draw();
   WJets1j->Draw("same");
   c->cd(2);
   weights1j->Draw();

   ofile->cd();
   Data1j->Write();
   WJets1j->Write();
   weights1j->Write();
   c->Write();
   ofile->Close();
   _file0->Close();
   _file1->Close();
}

/*
TFile* _file0 = new TFile("/eos/uscms/store/user/aperloff//MatrixElement/Summer12ME8TeV/MEResults/2015_06_05_microNtuples_optimized_BDT/rootFiles/microSingleEl_Data_19p148fb_BaseCuts.root","READ")
_file0->cd()
TH1F* DataEl2j = new TH1F("CosThetaL_Data_jets2_electron","CosThetaL_Data_jets2_electron",25,-1,1)
TH1F* DataEl3j = new TH1F("CosThetaL_Data_jets3_electron","CosThetaL_Data_jets3_electron",25,-1,1)
TH1F* DataEl4j = new TH1F("CosThetaL_Data_jets4_electron","CosThetaL_Data_jets4_electron",25,-1,1)
METree->Draw("CosTheta_l>>CosThetaL_Data_jets2_electron","@jLV.size()==2","goff")
METree->Draw("CosTheta_l>>CosThetaL_Data_jets3_electron","@jLV.size()==3","goff")
METree->Draw("CosTheta_l>>CosThetaL_Data_jets4_electron","@jLV.size()>=4","goff")
TFile* _file1 = new TFile("/eos/uscms/store/user/aperloff//MatrixElement/Summer12ME8TeV/MEResults/2015_06_05_microNtuples_optimized_BDT/rootFiles/microWJets_BaseCuts.root","READ")
_file1->cd()
TH1F* WJets2j = new TH1F("CosThetaL_WJets_jets2_electron","CosThetaL_WJets_jets2_electron",25,-1,1)
TH1F* WJets3j = new TH1F("CosThetaL_WJets_jets3_electron","CosThetaL_WJets_jets3_electron",25,-1,1)
TH1F* WJets4j = new TH1F("CosThetaL_WJets_jets4_electron","CosThetaL_WJets_jets4_electron",25,-1,1)
METree->Draw("CosTheta_l>>CosThetaL_WJets_jets2_electron","leptonCat==2&&@jLV.size()==2","goff")
METree->Draw("CosTheta_l>>CosThetaL_WJets_jets3_electron","leptonCat==2&&@jLV.size()==3","goff")
METree->Draw("CosTheta_l>>CosThetaL_WJets_jets4_electron","leptonCat==2&&@jLV.size()>=4","goff")
DataEl2j->Scale(1.0/DataEl2j->Integral())
DataEl3j->Scale(1.0/DataEl3j->Integral())
DataEl4j->Scale(1.0/DataEl4j->Integral())
WJets2j->Scale(1.0/WJets2j->Integral())
WJets3j->Scale(1.0/WJets3j->Integral())
WJets4j->Scale(1.0/WJets4j->Integral())
TH1F* weights2j = (TH1F*)DataEl2j->Clone("CosThetaL_WeightHist_jets2_electron")
TH1F* weights3j = (TH1F*)DataEl3j->Clone("CosThetaL_WeightHist_jets3_electron")
TH1F* weights4j = (TH1F*)DataEl4j->Clone("CosThetaL_WeightHist_jets4_electron")
weights2j->Divide(WJets2j)
weights3j->Divide(WJets3j)
weights4j->Divide(WJets4j)
DataEl2j->SetLineColor(kBlack)
DataEl3j->SetLineColor(kBlack)
DataEl4j->SetLineColor(kBlack)
WJets2j->SetLineColor(kGreen)
WJets3j->SetLineColor(kGreen)
WJets4j->SetLineColor(kGreen)
weights2j->SetLineColor(kRed)
weights3j->SetLineColor(kRed)
weights4j->SetLineColor(kRed)
TCanvas* c = new TCanvas("c","c",800,1000)
c->Divide(2,3)
c->cd(1)
DataEl2j->Draw()
WJets2j->Draw("same")
c->cd(2)
weights2j->Draw()
c->cd(3)
DataEl3j->Draw()
WJets3j->Draw("same")
c->cd(4)
weights3j->Draw()
c->cd(5)
DataEl4j->Draw()
WJets4j->Draw("same")
c->cd(6)
weights4j->Draw()
*/
