#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH1F.h"
#include "TAxis.h"
#include "TString.h"

#include <iostream>
#include <vector>

using namespace std;

vector<TH1F*> getInputHistograms(TFile* inputFile) {
	inputFile->cd();
	vector<TH1F*> ret;

	TIter next(inputFile->GetListOfKeys());
	TKey *key;
   	while ((key = (TKey*)next())) {
   		//cout << key->GetName() << endl;
   		TClass *cl = gROOT->GetClass(key->GetClassName());
      	if (!cl->InheritsFrom("TH1F")) continue;
      	if (!TString(key->GetName()).Contains("CosThetaL_WeightHist_")) continue;
      	TH1F *h = (TH1F*)key->ReadObj();
      	ret.push_back(h);
   }
   return ret;
}

void changeHistogramNamesAndTitles(vector<TH1F*> hists) {
	//CosThetaL_WeightHist_2Jets_ele
	//CosThetaL_WeightHist_2Jets_mu
	//CosThetaL_WeightHist_2Jets_lep
	for(unsigned int ihist = 0; ihist<hists.size(); ihist++) {
		TString name = hists[ihist]->GetName();
		if(TString(hists[ihist]->GetName()).Contains("2Jets")) {
			name.ReplaceAll("2Jets","jets2");
		}
		else if(TString(hists[ihist]->GetName()).Contains("3Jets")) {
			name.ReplaceAll("3Jets","jets3");
		}
		else if(TString(hists[ihist]->GetName()).Contains("4Jets")) {
			name.ReplaceAll("4Jets","jets4");
		}

		if(TString(hists[ihist]->GetName()).Contains("_ele")) {
			name.ReplaceAll("_ele","_electron");
		}
		else if(TString(hists[ihist]->GetName()).Contains("_mu")) {
			name.ReplaceAll("_mu","_muon");
		}
		else if(TString(hists[ihist]->GetName()).Contains("_lep")) {
			name.ReplaceAll("_lep","_both");
		}
		hists[ihist]->SetNameTitle(name,name);
	}
}

void setHistogramXaxisTitles(vector<TH1F*> hists) {
	for(unsigned int ihist = 0; ihist<hists.size(); ihist++) {
		hists[ihist]->GetXaxis()->SetTitle("CosTheta_l");
	}
}

void saveHistogramsToNewFile(TFile* outFile, vector<TH1F*> hists) {
	outFile->cd();
	for(unsigned int ihist = 0; ihist<hists.size(); ihist++) {
		hists[ihist]->Write();
	}
}

void ChangeCosThetaLWeightfileNames() {
	TFile* inputFile = TFile::Open("/uscms_data/d2/aperloff/YOURWORKINGAREA/MatrixElement/gitty/CMSSW_5_3_2_patch5/src/TAMUWW/ConfigFiles/Official/CosThetaLWeight_eq1tag_TTbarSub.root","READ");
	TFile* outFile = TFile::Open("/uscms_data/d2/aperloff/YOURWORKINGAREA/MatrixElement/gitty/CMSSW_5_3_2_patch5/src/TAMUWW/ConfigFiles/Official/CosThetaLWeight_eq1tag_TTbarSub_changedNames.root","RECREATE");

	cout << "Getting initial histograms ... ";
	vector<TH1F*> hists = getInputHistograms(inputFile);
	cout << "DONE" << endl
		 << "\tRetrieved " << hists.size() << " histograms" << endl;
	cout << "Changing the histogram names ... ";
	changeHistogramNamesAndTitles(hists);
	cout << "DONE" << endl;
	cout << "Setting the histogram x-axis titles ... ";
	setHistogramXaxisTitles(hists);
	cout << "DONE" << endl;
	cout << "Saving the histograms to a new file ... ";
	saveHistogramsToNewFile(outFile, hists);
	cout << "DONE" << endl;

	outFile->Close();
	inputFile->Close();
}

