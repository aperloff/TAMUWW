////////////////////////////////////////////////////////////////////////////////
//
// CompareSystematicShapes_x
// -------------------------
//
//                          08/03/2015 Alexx Perloff <aperloff@Physics.tamu.edu>
////////////////////////////////////////////////////////////////////////////////

// ROOT libraries
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TString.h"
#include "TLegend.h"

// C++ libraries
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <stdlib.h>

// My libraries
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"
#include "TAMUWW/SpecialTools/interface/Defs.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////////////////

int main(int argc,char**argv)
{
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;

   TString basepath = cl.getValue<TString> ("basepath", "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_04_27_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/");
   
   if (!cl.check()) 
      return 0;
   cl.print();
   
   vector<TString> jetBins;
   jetBins.push_back("Jets2");
   jetBins.push_back("Jets3");
   jetBins.push_back("Jets4");

   vector<TString> leptons;
   leptons.push_back("electron");
   leptons.push_back("muon");

   vector<TString> shapes;
   shapes.push_back("KinBDT");
   shapes.push_back("MEBDT");
   shapes.push_back("KinMEBDT");

   vector<TString> systematics;
   systematics.push_back("_CMS_hww_lnujj_matching");
   systematics.push_back("_CMS_hww_lnujj_scale");
   systematics.push_back("_CMS_hww_lnujj_QCDEtaWeight");
   systematics.push_back("_CMS_scale_j");
   systematics.push_back("_CMS_hww_lnujj_topPtWeight");

   vector<TString> sysLevels;
   sysLevels.push_back("");
   sysLevels.push_back("_shapeUp");
   sysLevels.push_back("_shapeDown");

   TString processNames[] = {"ggH125","qqH125","WH125_HToBB","WH125_HToBB_M125","WH_ZH_TTH_HToWW",
                             "WH_ZH_TTH_HToZZ","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125",
                             "WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","TTH_HToBB","TTH_HToBB_M125",
                             "WJets","ZJets","TTbar","Diboson","singleTop","QCD_ElFULL","QCD_MuFULL",
                             "WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","WW"};
   vector<TString> processes(processNames, processNames + sizeof(processNames) / sizeof(TString));

   TFile* ofile = TFile::Open(basepath+"/CompareSystematicShapes.root","RECREATE");

   for(unsigned int j = 0; j<jetBins.size(); j++) {
      for(unsigned int l = 0; l<leptons.size(); l++) {
         TFile* ifile = TFile::Open(basepath+"/"+jetBins[j]+"/"+leptons[l]+"/histos_"+leptons[l]+"_SysNames.root","READ");
         for(unsigned int s = 0; s<shapes.size(); s++) {
            for(unsigned int p = 0; p<processes.size(); p++) {
               for(unsigned int isys = 0; isys<systematics.size(); isys++) {
                  if((systematics[isys].Contains("matching") || systematics[isys].Contains("scale"))&&processes[p].CompareTo("WJets")!=0)
                     continue;
                  if(systematics[isys].Contains("topPtWeight") && processes[p].CompareTo("TTbar")!=0)
                     continue;
                  if(systematics[isys].Contains("QCDEtaWeight") && (processes[p].CompareTo("WJets")!=0 && !processes[p].Contains("QCD")))
                     continue;
                  if(leptons[l].CompareTo("electron")==0 && processes[p].CompareTo("QCD_MuFULL")==0)
                     continue;
                  if(leptons[l].CompareTo("muon")==0 && processes[p].CompareTo("QCD_ElFULL")==0)                  
                     continue;
                  TString canName = Form("CompareSystematicShapes_%s_%s_%s_%s%s",jetBins[j].Data(),leptons[l].Data(),
                                         shapes[s].Data(),processes[p].Data(),systematics[isys].Data());
                  cout << "Making canvas " << canName << " ..." << endl;
                  TCanvas * c = new TCanvas(canName,canName,600,600);
                  TLegend * leg = new TLegend(0.6,0.6,0.9,0.9);
                  for(unsigned int ilevel = 0; ilevel<sysLevels.size(); ilevel++) {
                     //Get input histograms
                     TString ihname;
                     if(ilevel==0)
                        ihname = Form("%s_%s_%s",shapes[s].Data(),processes[p].Data(),leptons[l].Data());
                     else
                        ihname = Form("%s_%s_%s%s%s",shapes[s].Data(),processes[p].Data(),leptons[l].Data(),
                                      systematics[isys].Data(),sysLevels[ilevel].Data());
                     cout << "\tGetting histogram " << ihname << " ..." << endl;
                     ifile->cd();
                     assert(ifile->Get(ihname));
                     TH1D* tmp = (TH1D*)ifile->Get(ihname);
                     if(tmp) {
                        if(ilevel==0) {
                           tmp->SetLineColor(kBlack);
                           tmp->SetMarkerColor(kBlack);
                           tmp->SetMarkerStyle(20);
                           tmp->Draw();
                           leg->AddEntry(tmp,"nominal","lp");
                        }
                        else if(ilevel==1) {
                           tmp->SetLineColor(kRed);
                           tmp->SetMarkerColor(kRed);
                           tmp->SetMarkerStyle(21);
                           tmp->Draw("same");
                           leg->AddEntry(tmp,Form("%s_%s",systematics[isys].Data(),sysLevels[ilevel].Data()),"lp");
                        }
                        else {
                           tmp->SetLineColor(kBlue);
                           tmp->SetMarkerColor(kBlue);
                           tmp->SetMarkerStyle(22);
                           tmp->Draw("same");
                           leg->AddEntry(tmp,Form("%s_%s",systematics[isys].Data(),sysLevels[ilevel].Data()),"lp");
                        }
                        leg->Draw("same");
                     }
                  }
                  ofile->cd();
                  c->Write();
               }
            }
         }
         ifile->Close();
      }
   }

   ofile->Close();

   return 0;
}
