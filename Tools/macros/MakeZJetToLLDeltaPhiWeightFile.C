#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TEnv.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLorentzVector.h"
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

////////////////////////////////////////////////////////////////////////////////
//  Declare Global Variables
////////////////////////////////////////////////////////////////////////////////
const double ZMASS = 91.1876; //GeV
const double WMASS = 80.3980; //GeV

////////////////////////////////////////////////////////////////////////////////
//  Local Functions
////////////////////////////////////////////////////////////////////////////////

// Get the values for the jetBins, leptonBins, and processes vectors
void getBinValues(vector<string>& jetBins, vector<string>& leptonBins, vector<string>& processes);

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void MakeZJetsToLLDeltaPhiWeightFile(string input_file_basepath) {
    // Trying to speed up the code
    gEnv->SetValue("TFile.AsyncPrefetching", 1);

    Style* st = new Style();
    st->setTDRStyle();

    cout << "Getting the binning values ... " << flush;
    vector<string> jetBins;
    vector<string> leptonBins;
    vector<string> processes;
    getBinValues(jetBins,leptonBins,processes);
    cout << "DONE" << endl;

    cout << "Opening input files, get the histograms, and closing the input files ... " << flush;
    TFile* current_input_file;
    map<string, TH1D*> histograms;
    string hname, htitle, objectName = "DeltaPhi_GenLNuOrRecoLL_", current_input_file_name;
    for(unsigned int jbin=0; jbin<jetBins.size(); jbin++) {
        for(unsigned int lbin=0; lbin<leptonBins.size(); lbin++) {
            current_input_file_name = input_file_basepath+"/"+jetBins[jbin]+"/"+leptonBins[lbin]+
                                      "/histos_"+leptonBins[lbin]+"_"+jetBins[jbin]+".root";
            current_input_file = TFile::Open(current_input_file_name.c_str(),"READ");
            if(!current_input_file) {
                cout << "ERROR::MakeZJetsToLLDeltaPhiWeightFile Couldn't open the input file " << current_input_file_name << endl;
                std::terminate();
            }
            for(unsigned int iprocess=0; iprocess<processes.size(); iprocess++) {
                hname = processes[iprocess]+"_"+jetBins[jbin]+"_"+leptonBins[lbin];
                htitle = "#Delta#phi(l,#nu) (";
                htitle+= processes[iprocess]=="WJets"?"GEN)":"RECO)";
                histograms[hname] = dynamic_cast<TH1D*>(current_input_file->Get((objectName+processes[iprocess]+
                                                       "_"+leptonBins[lbin]).c_str()));
                histograms[hname]->SetNameTitle(hname.c_str(),htitle.c_str());
                histograms[hname]->SetDirectory(0);
            }
            current_input_file->Close();
        }
    }
    cout << "DONE" << endl;

    cout << "Opening output file ... " << flush;
    TFile* ofile = TFile::Open("ZllDeltaPhiWeights.root","RECREATE");
    ofile->mkdir("validation");
    cout << "DONE" << endl;

    cout << "Checking the input histograms ... " << flush;
    for(map<string,TH1D*>::iterator it=histograms.begin(); it!=histograms.end(); it++) {
        if(it->second->GetEntries()==0) {
            cout << endl << "ERROR::MakeZJetsToLLDeltaPhiWeightFile The histogtam " << it->first << " has zero entries." << endl;
            std::terminate();
        }
    }
    cout << "DONE" << endl;

    cout << "Add ZJetsToLL_M10To50 to ZJetsToLL_M50 and remove ZJetsToLL_M10To50 from map ... " << flush;
    for(unsigned int jbin=0; jbin<jetBins.size(); jbin++) {
        for(unsigned int lbin=0; lbin<leptonBins.size(); lbin++) {
            string M10To50_name = "ZJetsToLL_M10To50_"+jetBins[jbin]+"_"+leptonBins[lbin];
            string M50_name = "ZJetsToLL_M50_"+jetBins[jbin]+"_"+leptonBins[lbin];
            string combined_name = "ZJetsToLL_"+jetBins[jbin]+"_"+leptonBins[lbin];
            if(histograms.find(M10To50_name)!=histograms.end() && histograms.find(M50_name)!=histograms.end()) {
                histograms[M50_name]->Add(histograms[M10To50_name]);
                histograms.erase(M10To50_name);
                histograms[M50_name]->SetName(combined_name.c_str());
                histograms[combined_name] = histograms[M50_name];
                histograms.erase(M50_name);
            }
            else if(histograms.find(M10To50_name)!=histograms.end() && histograms.find(M50_name)==histograms.end()) {
                cout << endl << "ERROR::MakeZJetsToLLDeltaPhiWeightFile Found the histogram " << M10To50_name
                     << ", but not the histogram " << M50_name << endl
                     << "\tContinuing now would be inaccurate" << endl;
                std::terminate();                    
            }
            else if(histograms.find(M10To50_name)==histograms.end() && histograms.find(M50_name)!=histograms.end()) {
                cout << endl << "ERROR::MakeZJetsToLLDeltaPhiWeightFile Found the histogram " << M50_name
                     << ", but not the histogram " << M10To50_name << endl
                     << "\tContinuing now would be inaccurate" << endl;
                std::terminate();         
            }
            else {
                cout << endl << "ERROR::MakeZJetsToLLDeltaPhiWeightFile Couldn't fine either the histogram "
                     << M50_name << " or the histogram " << M10To50_name << endl
                     << "\tContinuing now would make no sense" << endl;
                std::terminate();         
            }
        }
    }
    cout << "DONE" << endl;

    cout << "Scaling ZJetToLL to WJets, making the weights, and writing the weights ... " << flush;
    for(unsigned int jbin=0; jbin<jetBins.size(); jbin++) {
        for(unsigned int lbin=0; lbin<leptonBins.size(); lbin++) {
            string zll_name = "ZJetsToLL_"+jetBins[jbin]+"_"+leptonBins[lbin];
            string w_name = "WJets_"+jetBins[jbin]+"_"+leptonBins[lbin];
            histograms[zll_name]->Scale(histograms[w_name]->Integral()/histograms[zll_name]->Integral());

            hname = "ZllWeight_"+jetBins[jbin]+"_"+leptonBins[lbin];
            histograms[hname] = (TH1D*)histograms[w_name]->Clone(hname.c_str());
            histograms[hname]->SetTitle(hname.c_str());
            histograms[hname]->Divide(histograms[zll_name]);

            ofile->cd();
            histograms[hname]->Write();
        }
    }
    cout << "DONE" << endl;

    cout << "Making and writing a validation canvas ... " << endl;
    map<string, TCanvas*> canvases;
    for(unsigned int jbin=0; jbin<jetBins.size(); jbin++) {
        for(unsigned int lbin=0; lbin<leptonBins.size(); lbin++) {
            string zll_name = "ZJetsToLL_"+jetBins[jbin]+"_"+leptonBins[lbin];
            string w_name = "WJets_"+jetBins[jbin]+"_"+leptonBins[lbin];
            hname = "ZllWeight_"+jetBins[jbin]+"_"+leptonBins[lbin];
            ofile->cd("validation");
            string cname = "WeightValidation_"+jetBins[jbin]+"_"+leptonBins[lbin];
            TH1D* hup = new TH1D("hup","hup",62,-TMath::Pi(),TMath::Pi());
            hup->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
            hup->GetYaxis()->SetRangeUser(0.0,histograms[zll_name]->GetMaximum()*1.2);
            TH1D* hdown = new TH1D("hdown","hdown",62,-TMath::Pi(),TMath::Pi());
            hdown->GetXaxis()->SetLimits(-TMath::Pi(),TMath::Pi());
            hdown->GetXaxis()->SetTitle("#Delta#phi(l,#nu(l))");
            hdown->GetYaxis()->SetRangeUser(0,2);
            hup->GetYaxis()->SetTitle("Events");
            hdown->GetYaxis()->SetTitle("Weight");
            canvases[cname] = st->tdrDiCanvas(cname.c_str(),hup,hdown);
            canvases[cname]->cd(1);
            st->tdrDraw(histograms[w_name],"",kFullCircle,kGreen,kSolid,kGreen);
            st->tdrDraw(histograms[zll_name],"",kOpenCircle,kRed,kSolid,kRed);
            TLegend* l = st->tdrLeg(0.5,0.8,0.85,0.9);
            l->AddEntry(histograms[w_name],"WJets (GEN)", "lep");
            l->AddEntry(histograms[zll_name],"ZJetsToLL (RECO)", "lep");
            l->Draw("same");
            canvases[cname]->cd(2);
            st->tdrDraw(histograms[hname],"",kFullCircle,kBlue,kSolid,kBlue);
            canvases[cname]->Write();
        }
    }
    cout << "DONE" << endl;

    cout << "Closing the output file ... " << flush;
    ofile->Close();
    cout << "DONE" << endl;
}

////////////////////////////////////////////////////////////////////////////////
//  Implement Local Functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void getBinValues(vector<string>& jetBins, vector<string>& leptonBins, vector<string>& processes) {
    jetBins.push_back("jets2");
    jetBins.push_back("jets3");
    jetBins.push_back("jets4");

    leptonBins.push_back("electron");
    leptonBins.push_back("muon");

    processes.push_back("WJets");
    processes.push_back("ZJetsToLL_M10To50");
    processes.push_back("ZJetsToLL_M50");

    return;
}
