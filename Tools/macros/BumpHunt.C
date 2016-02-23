#include "TROOT.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TError.h"
#include "TString.h"
#include "TLegend.h"

#include <iostream>
#include <vector>

using namespace std;

void BumpHunt(TString jetBin, TString lepton) {
	gROOT->SetBatch(kTRUE);
	//argc = (batch) ? 2 : 1; if (batch) argv[1] = (char*)"-b";
	//TApplication* app = new TApplication("jet_l3_correction_x",&argc,argv);

	cout << "Doing jet bin " << jetBin << " and lepton " << lepton << " ... " << endl;

	TString distributions[19] = {"MVADiscriminator","MEBDT","leptonEtaCharge","LeptPt","WmT","Mlvjj","jet1dRLep","jet2dRLep","jet3dRLep","ht","Ptlnujj","dRlepjj","dPhiMETJet","minDPhiLepJet","DeltaEtaJ1J2","DeltaPhi_J1J2","CosTheta_l","CosTheta_j","CosTheta_WH"};
	TString processes[10] = {"WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","WH125_HToBB","TTH_HToBB_M125","ggH125","qqH125","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125"};

	TFile *bumpFile = new TFile(Form("/uscms_data/d2/aperloff/Summer12ME8TeV/2015_04_27_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/MVA_Bump_SignalOnly/%s/%s/histos_%s.root",jetBin.Data(),lepton.Data(),lepton.Data()),"READ");
	TFile *antiBumpFile = new TFile(Form("/uscms_data/d2/aperloff/Summer12ME8TeV/2015_04_27_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/MVA_AntiBump_signalOnly/%s/%s/histos_%s.root",jetBin.Data(),lepton.Data(),lepton.Data()),"READ");
	TFile *nominalFile = new TFile(Form("/uscms_data/d2/aperloff/Summer12ME8TeV/2015_04_27_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/%s/%s/histos_%s.root",jetBin.Data(),lepton.Data(),lepton.Data()),"READ");

	TCanvas* c;
	for(unsigned int i=0; i<19; i++) {
		for(unsigned int p=0; p<10; p++) {

			TString name = Form("%s_%s_%s_%s",distributions[i].Data(),lepton.Data(),jetBin.Data(),processes[p].Data());
			cout << "\tDoing canvas " << name << " ... " << endl;
			c = new TCanvas(name,name);

			TString hname = Form("%s_%s_%s",distributions[i].Data(),processes[p].Data(),lepton.Data());
			cout << "\t\tDoing histogram " << hname << " ... " << endl;
			TH1D* h_bump = (TH1D*)bumpFile->Get(hname.Data());
			TH1D* h_antiBump = (TH1D*)antiBumpFile->Get(hname.Data());
			TH1D* h_nominal = (TH1D*)nominalFile->Get(hname.Data());

			if(!h_bump) {
				cout << "ERROR::BumpHunt Could not find histogram " << hname << " in file " << bumpFile->GetName() << endl;
				continue; 
			}
			if(!h_antiBump) {
				cout << "ERROR::BumpHunt Could not find histogram " << hname << " in file " << antiBumpFile->GetName() << endl;
				continue; 
			}
			if(!h_nominal && (distributions[i].Contains("MVA")||distributions[i].Contains("BDT"))) {
				cout << "ERROR::BumpHunt Could not find histogram " << hname << " in file " << nominalFile->GetName() << endl;
				continue; 
			}

			h_bump->Scale(1.0/h_bump->Integral());
			h_antiBump->Scale(1.0/h_antiBump->Integral());
			if(h_nominal) h_nominal->Scale(1.0/h_nominal->Integral());

			if(distributions[i].CompareTo("MVADiscriminator")==0 && jetBin.CompareTo("Jets4")==0 && lepton.CompareTo("electron")==0)
				h_bump->GetYaxis()->SetRangeUser(0,1.1);

			h_bump->Draw();
			h_antiBump->Draw("same");
			if(h_nominal) h_nominal->Draw("same");

			h_bump->SetLineColor(kRed);
			h_antiBump->SetLineColor(kBlue);
			if(h_nominal) h_nominal->SetLineColor(kBlack);

			TLegend* leg = new TLegend(0.8,0.8,1.0,1.0);
			leg->AddEntry(h_bump,"MVADiscriminator>0.2","l");
			leg->AddEntry(h_antiBump,"MVADiscriminator<0","l");
			if(h_nominal) leg->AddEntry(h_nominal,"Nominal","l");
			leg->Draw("same");

			c->SaveAs(Form("%s.png",name.Data()));
		}
	}
}