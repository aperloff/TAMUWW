#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
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
#include "TH2D.h"
#include "TError.h"
#include "TString.h"
#include "TAxis.h"
#include "TLine.h"

#include <iostream>
#include <string>
#include <vector>

#include "TAMUWW/Tools/interface/Style.hh"
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//  Local Functions
////////////////////////////////////////////////////////////////////////////////

// Get the values for the jetBins, leptonBins, and processes vectors
void getBinValues(vector<string>& jetBins, vector<string>& leptonBins,
                  vector<string>& processes, vector<double>& MET_binning);

// Plot the simplest necessary canvas
void PlotBasic();

// Once the final lepton pT split has been decided for electrons and muons combined two files together
void combineTwoFiles(string electronPercentage = "050", string muonPercentage = "050", bool saveToConfig = true);

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void MakeZJetsToLLGenNeutrinoWeightFile(string input_file_basepath, string highPtLeptonPercentage = "050",
                                        bool batch = true, bool close = false, bool save = false,
                                        bool saveToConfig = true, bool validation = false) {
    if(batch)
        gROOT->SetBatch(kTRUE);

    // Trying to speed up the code
    gEnv->SetValue("TFile.AsyncPrefetching", 1);

    Style* st = new Style();
    st->setTDRStyle();

    cout << "Getting the binning values ... " << flush;
    vector<string> jetBins;
    vector<string> leptonBins;
    vector<string> processes;
    vector<double> MET_binning;
    getBinValues(jetBins,leptonBins,processes,MET_binning);
    cout << "DONE" << endl;

    cout << "Opening input files, get the histograms, and closing the input files ... " << flush;
    TFile* current_input_file;
    map<string, TH1D*> histograms;
    string hname, htitle, objectName = "Pt_GenNuOrRecoMETLept_", current_input_file_name;
    for(unsigned int jbin=0; jbin<jetBins.size(); jbin++) {
        for(unsigned int lbin=0; lbin<leptonBins.size(); lbin++) {
            current_input_file_name = input_file_basepath+"/HighPtLeptonKept"+highPtLeptonPercentage+
                                      "/"+jetBins[jbin]+"/"+leptonBins[lbin]+"/histos_"+leptonBins[lbin]+
                                      "_"+jetBins[jbin]+".root";
            current_input_file = TFile::Open((current_input_file_name).c_str(),"READ");
            for(unsigned int iprocess=0; iprocess<processes.size(); iprocess++) {
                hname = processes[iprocess]+"_"+jetBins[jbin]+"_"+leptonBins[lbin];
                htitle = "p_{T}^{#nu} ("+processes[iprocess]+")";
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
    string ofile_name = (saveToConfig) ? DefaultValues::getConfigPath()+"ZllGenNeutrinoPtWeightFiles/" : "./";
    ofile_name += "ZllGenNeutrinoPtWeightFile_HighPtLeptonKept"+highPtLeptonPercentage+"_WGenZRec.root";
    TFile* ofile = TFile::Open(ofile_name.c_str(),"RECREATE");
    ofile->mkdir("validation");
    cout << "DONE" << endl;

    cout << "Checking the input histograms ... " << flush;
    for(map<string,TH1D*>::iterator it=histograms.begin(); it!=histograms.end(); it++) {
        if(it->second->GetEntries()==0) {
            cout << endl << "ERROR::MakeZJetsToLLGenNeutrinoWeightFile The histogtam " << it->first << " has zero entries." << endl;
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
                cout << endl << "ERROR::MakeZJetsToLLGenNeutrinoWeightFile Found the histogram " << M10To50_name
                     << ", but not the histogram " << M50_name << endl
                     << "\tContinuing now would be inaccurate" << endl;
                std::terminate();                    
            }
            else if(histograms.find(M10To50_name)==histograms.end() && histograms.find(M50_name)!=histograms.end()) {
                cout << endl << "ERROR::MakeZJetsToLLGenNeutrinoWeightFile Found the histogram " << M50_name
                     << ", but not the histogram " << M10To50_name << endl
                     << "\tContinuing now would be inaccurate" << endl;
                std::terminate();         
            }
            else {
                cout << endl << "ERROR::MakeZJetsToLLGenNeutrinoWeightFile Couldn't fine either the histogram "
                     << M50_name << " or the histogram " << M10To50_name << endl
                     << "\tContinuing now would make no sense" << endl;
                std::terminate();         
            }
        }
    }
    cout << "DONE" << endl;

    cout << "Sacaling ZJetsToLL to WJets and making the ratio histogram ... " << flush;
    for(unsigned int jbin=0; jbin<jetBins.size(); jbin++) {
        for(unsigned int lbin=0; lbin<leptonBins.size(); lbin++) {
            string wjets_name = "WJets_"+jetBins[jbin]+"_"+leptonBins[lbin];
            string zjets_name = "ZJetsToLL_"+jetBins[jbin]+"_"+leptonBins[lbin];
            hname = "ZllWeight_"+jetBins[jbin]+"_"+leptonBins[lbin];
            histograms[zjets_name]->Scale(dynamic_cast<TH1D*>(histograms[wjets_name])->Integral()/histograms[zjets_name]->Integral());
            // Needed because WJets has three more bins of gen neutrino pt <15.0 GeV and ZJetsToLL does not have those bins
            // ((WJets Integral-Sum of content of 3 bins)/WJets Integral) = scale factor
            // From canvas: ((1.48201e+06-96741.5)/1.48201e+06)=0.9347227751
            // From command line: ((1.48200519787897007e+06-9.67414995618765242e+04)/1.48200519787897007e+06)=9.34722563928701566e-01
            // Measured: 9.34690809464712835e-01
            if(validation) {
                histograms[zjets_name]->Scale(9.34690809464712835e-01);
            }
            histograms[hname] = (TH1D*)histograms[wjets_name]->Clone(hname.c_str());
            histograms[hname]->Divide(dynamic_cast<TH1D*>(histograms[zjets_name]));
        }
    }
    cout << "DONE" << endl;

    cout << "Making and writing a validation canvas ... " << flush;
    map<string, TCanvas*> canvases;
    map<string, TLegend*> legends;
    TH1D *hup(0), *hdown(0);
    TLine *line(0);
    for(unsigned int jbin=0; jbin<jetBins.size(); jbin++) {
        for(unsigned int lbin=0; lbin<leptonBins.size(); lbin++) {
            string w_name = "WJets_"+jetBins[jbin]+"_"+leptonBins[lbin];
            string zll_name = "ZJetsToLL_"+jetBins[jbin]+"_"+leptonBins[lbin];
            string r_name = "ZllWeight_"+jetBins[jbin]+"_"+leptonBins[lbin];
            double max_y = max(histograms[w_name]->GetMaximum(),
                               histograms[zll_name]->GetMaximum());
            ofile->cd("validation");
            string cname = "Canvas_"+jetBins[jbin]+"_"+leptonBins[lbin];

            hup = new TH1D("hup","hup",10000,MET_binning.front(),MET_binning.back());
            hup->GetXaxis()->SetLimits(1,MET_binning.back());
            hup->GetYaxis()->SetRangeUser(1,max_y*1.2);
            hup->GetYaxis()->SetTitle("Entries");
            hdown = new TH1D("hdown","hdown",10000,MET_binning.front(),MET_binning.back());
            hdown->GetXaxis()->SetLimits(1,MET_binning.back());
            hdown->GetYaxis()->SetRangeUser(0.1,histograms[r_name]->GetMaximum()*1.2);
            hdown->GetXaxis()->SetTitle("p_{T}^{#nu(l)} (GeV)");
            hdown->GetYaxis()->SetTitle("Ratio");
            canvases[cname] = st->tdrDiCanvas(cname.c_str(),hup,hdown,2,11); //8TeV, inside the frame

            canvases[cname]->cd(1)->SetLogx();
            st->tdrDraw(histograms[w_name],"E1",kFullCircle,kGreen);
            st->tdrDraw(histograms[zll_name],"E1",kOpenCircle,kRed);
            string l_name = "legend_"+jetBins[jbin]+"_"+leptonBins[lbin];
            legends[l_name] = st->tdrLeg(0.60,0.75,0.89,0.89);
            legends[l_name]->AddEntry(histograms[w_name],"WJets (GEN)","ep");
            legends[l_name]->AddEntry(histograms[zll_name],"ZJetsToLL (RECO)","p");
            legends[l_name]->Draw("same");

            canvases[cname]->cd(2);
            gPad->SetLogx();
            gPad->SetLogy();
            line = new TLine(hdown->GetXaxis()->GetXmin(),1,hdown->GetXaxis()->GetXmax(),1);
            line->SetLineStyle(2);
            line->Draw("same");
            st->tdrDraw(histograms[r_name],"E1",kFullCircle,kBlack);

            canvases[cname]->Write();
        }
    }
    cout << "DONE" << endl;

    cout << "Writing the weight histograms ... " << flush;
    ofile->cd();
    for(unsigned int jbin=0; jbin<jetBins.size(); jbin++) {
        for(unsigned int lbin=0; lbin<leptonBins.size(); lbin++) {
            hname = "ZllWeight_"+jetBins[jbin]+"_"+leptonBins[lbin];
            histograms[hname]->Write();
        }
    }
    cout << "DONE" << endl;

    if(save) {
        cout << "Saving all of the canvases ... " << flush;
        for (std::map<string,TCanvas*>::iterator it=canvases.begin(); it!=canvases.end(); ++it) {
            it->second->SaveAs((it->first+".eps").c_str());
            it->second->SaveAs((it->first+".png").c_str());
        }
        cout << "DONE" << endl;
    }

    cout << "Closing the output file ... " << flush;
    ofile->Close();
    cout << "DONE" << endl;

    if(batch && close)
        gApplication->Terminate();
    return;
}

////////////////////////////////////////////////////////////////////////////////
//  Implement Local Functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void getBinValues(vector<string>& jetBins, vector<string>& leptonBins,
                  vector<string>& processes, vector<double>& MET_binning) {
    jetBins.push_back("jets2");
    jetBins.push_back("jets3");
    jetBins.push_back("jets4");

    leptonBins.push_back("electron");
    leptonBins.push_back("muon");

    processes.push_back("WJets");
    processes.push_back("ZJetsToLL_M10To50");
    processes.push_back("ZJetsToLL_M50");

    for(int i=0; i<=300; i+=5) {
        MET_binning.push_back(i);
    }
    MET_binning.push_back(1000);

    return;
}

//______________________________________________________________________________
void PlotBasic() {
    TFile* ofile = TFile::Open((DefaultValues::getConfigPath()+"ZllGenNeutrinoWeightFile.root").c_str(),"READ");
    ((TH1D*)ofile->Get("ZllWeight_jets2_electron"))->Draw("P");
    ((TCanvas*)ofile->Get("validation/Canvas_jets2_electron"))->Draw();
}

//______________________________________________________________________________
void combineTwoFiles(string electronPercentage, string muonPercentage, bool saveToConfig) {
    gROOT->SetBatch(kTRUE);

    Style* st = new Style();
    st->setTDRStyle();

    cout << "Getting the binning values ... " << flush;
    vector<string> jetBins;
    vector<string> leptonBins;
    vector<string> processes;
    vector<double> MET_binning;
    getBinValues(jetBins,leptonBins,processes,MET_binning);
    cout << "DONE" << endl;

    cout << "Opening input files, get the histograms and canvases, and closing the input files ... " << endl;
    TFile* current_input_file;
    map<string, TH1D*> histograms;
    map<string, TCanvas*> canvases;
    string hname, cname, dir, current_input_file_name;
    for(unsigned int lbin=0; lbin<leptonBins.size(); lbin++) {
        current_input_file_name = DefaultValues::getConfigPath()+"ZllGenNeutrinoPtWeightFiles/ZllGenNeutrinoPtWeightFile_HighPtLeptonKept"+
                          ((leptonBins[lbin]=="electron") ? electronPercentage : muonPercentage) + "_WGenZRec.root";
        cout << "\tOpening the file " << current_input_file_name << endl;
        current_input_file = TFile::Open((current_input_file_name).c_str(),"READ");
        for(unsigned int jbin=0; jbin<jetBins.size(); jbin++) {
            current_input_file->cd("");
            hname = "ZllWeight_"+jetBins[jbin]+"_"+leptonBins[lbin];
            cout << "\t\tGetting histogram " << hname << " ... " << flush;
            if(current_input_file->Get(hname.c_str())) {
                histograms[hname] = dynamic_cast<TH1D*>(current_input_file->Get(hname.c_str()));
                histograms[hname]->SetDirectory(0);
                cout << "DONE" << endl;
            }
            else {
                cout << "FAIL" << endl;
            }
            current_input_file->cd("validation");
            dir = "validation/";
            cname = "Canvas_"+jetBins[jbin]+"_"+leptonBins[lbin];
            cout << "\t\tGetting canvas " << dir << cname << "... " << flush;
            if(current_input_file->Get((dir+cname).c_str())) {
                canvases[cname] = dynamic_cast<TCanvas*>(current_input_file->Get((dir+cname).c_str()));
                cout << "DONE" << endl;
            }
            else {
                cout << "FAIL" << endl;
            }
        }
        current_input_file->Close();
    }
    cout << "DONE" << endl;

    cout << "Opening output file ... " << flush;
    string ofile_name = (saveToConfig) ? DefaultValues::getConfigPath()+"ZllGenNeutrinoPtWeightFiles/" : "./";
    ofile_name += "ZllGenNeutrinoPtWeightFile_HighPtLeptonKept"+muonPercentage+"mu"+electronPercentage+"el"+"_WGenZRec.root";
    TFile* ofile = TFile::Open(ofile_name.c_str(),"RECREATE");
    ofile->mkdir("validation");
    cout << "DONE" << endl;

    cout << "Checking the input histograms ... " << flush;
    unsigned int totalHistograms = jetBins.size()*leptonBins.size();
    if(histograms.size()!=totalHistograms) {
        cout << endl << "ERROR::combineTwoFiles Something is wrong. We should have found " << totalHistograms
             << " histograms, but instead we only found " << histograms.size() << endl;
        std::terminate();
    }
    if(canvases.size()!=totalHistograms) {
        cout << endl << "ERROR::combineTwoFiles Something is wrong. We should have found " << totalHistograms
             << " canvases, but instead we only found " << canvases.size() << endl;
        std::terminate();
    }
    cout << "DONE" << endl;

    cout << "Writing the weight histograms ... " << flush;
    ofile->cd();
    for(unsigned int jbin=0; jbin<jetBins.size(); jbin++) {
        for(unsigned int lbin=0; lbin<leptonBins.size(); lbin++) {
            hname = "ZllWeight_"+jetBins[jbin]+"_"+leptonBins[lbin];
            histograms[hname]->Write();
        }
    }
    cout << "DONE" << endl;

    cout << "Writing the validation canvases ... " << flush;
    ofile->cd("validation");
    for(unsigned int jbin=0; jbin<jetBins.size(); jbin++) {
        for(unsigned int lbin=0; lbin<leptonBins.size(); lbin++) {
            cname = "Canvas_"+jetBins[jbin]+"_"+leptonBins[lbin];
            canvases[cname]->Write();
        }
    }
    cout << "DONE" << endl;    

    cout << "Closing the output file ... " << flush;
    ofile->Close();
    cout << "DONE" << endl;

    return;
}