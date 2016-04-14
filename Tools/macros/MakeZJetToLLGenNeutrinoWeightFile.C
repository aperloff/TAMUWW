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
#include "TAMUWW/SpecialTools/interface/ProgressBar.hh"
#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"
#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void MakeZJetsToLLGenNeutrinoWeightFile(string highPtLeptonPercentage = "050", int maxEvts = 0,
                                        bool batch = true, bool close = false, bool save = false) {
    if(batch)
        gROOT->SetBatch(kTRUE);

    // Trying to speed up the code
    gEnv->SetValue("TFile.AsyncPrefetching", 1);

    Style* st = new Style();
    st->setTDRStyle();
    //TGaxis::SetMaxDigits(3);

    cout << "Opening output file ... " << flush;
    TFile* ofile = TFile::Open((DefaultValues::getConfigPath()+"ZllGenNeutrinoPtWeightFiles/ZllGenNeutrinoPtWeightFile_HighPtLeptonKept"+
                               highPtLeptonPercentage+"_WGenZRec.root").c_str(),"RECREATE");
    ofile->mkdir("validation");
    cout << "DONE" << endl;

    const int nJetBins = 3;
    const int nLeptonBins = 2;
    const int nProcesses = 2;
    string jetBins[nJetBins] = {"jets2","jets3","jets4"};
    string leptonBins[nLeptonBins] = {"electron","muon"};
    string processes[nProcesses] = {"WJets","ZJetsToLL"};

    TTree *jets2p(0);
    string hname, htitle, gname;
    map<string, TH1*> histograms;
    map<string, TCanvas*> canvases;
    map<string, TLegend*> legends;
    TLorentzVector lepton;
    TLorentzVector neutrino;
    int *nu_pdgid = new int(0);
    vector<TFile*> ifiles;
    EventNtuple* ntuple(0);

    cout << "Creating the histograms ... " << flush;
    vector<double> MET_binning;
    for(int i=0; i<=300; i+=5) {
        MET_binning.push_back(i);
    }
   MET_binning.push_back(1000);
    for(int jbin=0; jbin<nJetBins; jbin++) {
        for(int lbin=0; lbin<nLeptonBins; lbin++) {
            for(int iprocess=0; iprocess<nProcesses; iprocess++) {
                hname = processes[iprocess]+"_"+jetBins[jbin]+"_"+leptonBins[lbin];
                htitle = "p_{T}^{#nu} ("+processes[iprocess]+")";
                histograms[hname] = new TH1D(hname.c_str(),htitle.c_str(),MET_binning.size()-1,&MET_binning[0]);
                histograms[hname]->Sumw2();
            }
        }
    }
    cout << "DONE" << endl;

    cout << "Opening input file ... " << flush;
    ifiles.push_back(TFile::Open("$SMEInput/WJets.root","READ"));
    ifiles.push_back(TFile::Open(("$SMEInput/WlnuJetsTest/WlnuJets_M-10To50_HighPtLeptonKept"+highPtLeptonPercentage+".root").c_str(),"READ"));
    ifiles.push_back(TFile::Open(("$SMEInput/WlnuJetsTest/WlnuJets_M-50_HighPtLeptonKept"+highPtLeptonPercentage+".root").c_str(),"READ"));
    cout << "DONE" << endl;

    cout << "Using files ... " << endl;
    for(unsigned int ifile=0; ifile<ifiles.size(); ifile++) {
        cout << "\t" << ifiles[ifile]->GetName() << endl;
    }

    cout << "Filling the input histograms ... " << endl;
    for(unsigned int ifile=0; ifile<ifiles.size(); ifile++) {
        int ifile_for_process = ifile<2 ? ifile : ifile-1;
        ntuple = new EventNtuple();
        jets2p = (TTree*)ifiles[ifile]->Get("PS/jets2p");
        jets2p->SetBranchAddress("EvtNtuple", &ntuple);
        jets2p->SetBranchStatus("*",0);
        jets2p->SetBranchStatus("METLV.refLV*",1);
        jets2p->SetBranchStatus("jLV",1);
        jets2p->SetBranchStatus("lLV.leptonCat",1);
        jets2p->SetBranchStatus("genParticleCollection*",1);
        unsigned int nEntries = (maxEvts>0 && maxEvts<jets2p->GetEntries()) ? maxEvts : jets2p->GetEntries();
        for(unsigned int ientry=0; ientry<nEntries; ientry++) {
            ProgressBar::loadbar2(ientry+1,nEntries);
            jets2p->GetEntry(ientry);
            hname = processes[ifile_for_process]+"_"+DEFS::getJetBinString(DEFS::getJetBin(ntuple->jLV.size()))+"_"+DEFS::getLeptonCatString(ntuple->lLV[0].leptonCat);
            if(histograms.find(hname)==histograms.end()) {
                cout << "ERROR::MakeZJetsToLLGenNeutrinoWeightFile Cannot find the histogram named " << hname << endl;
                std::terminate();
            }
            if(processes[ifile_for_process]=="ZJetsToLL") {
                histograms.find(hname)->second->Fill(ntuple->METLV[0].refLV.Pt());
            }
            else{
                neutrino = ntuple->getGenVorDaughter(EventNtuple::NEUTRINO, 24, 0, nu_pdgid, false);
                if(neutrino.Pt()>0) {
                    histograms.find(hname)->second->Fill(neutrino.Pt());
                }
            }
        }
        cout << endl;
        delete ntuple;
    }

    cout << "Checking filled histograms ... " << flush;
    for(map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); it++) {
        if(it->second->GetEntries()==0) {
            cout << endl << "ERROR::MakeZJetsToLLGenNeutrinoWeightFile The histogtam " << it->first << " has zero entries." << endl;
            std::terminate();
        }
    }
    cout << "DONE" << endl;

    cout << "Make the ratio histogram ... " << flush;
    for(int jbin=0; jbin<nJetBins; jbin++) {
        for(int lbin=0; lbin<nLeptonBins; lbin++) {
            string wjets_name = "WJets_"+jetBins[jbin]+"_"+leptonBins[lbin];
            string zjets_name = "ZJetsToLL_"+jetBins[jbin]+"_"+leptonBins[lbin];
            hname = "ZllWeight_"+jetBins[jbin]+"_"+leptonBins[lbin];
            histograms[zjets_name]->Scale(dynamic_cast<TH1D*>(histograms[wjets_name])->Integral()/histograms[zjets_name]->Integral());
            histograms[hname] = (TH1D*)histograms[wjets_name]->Clone(hname.c_str());
            histograms[hname]->Divide(dynamic_cast<TH1D*>(histograms[zjets_name]));
        }
    }
    cout << "DONE" << endl;

    cout << "Making and writing a validation canvas ... " << flush;
    TH1D *hup(0), *hdown(0);
    TLine *line(0);
    for(int jbin=0; jbin<nJetBins; jbin++) {
        for(int lbin=0; lbin<nLeptonBins; lbin++) {
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
            st->tdrDraw(histograms[zll_name],"E1",kFullCircle,kRed);
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
    for(int jbin=0; jbin<nJetBins; jbin++) {
        for(int lbin=0; lbin<nLeptonBins; lbin++) {
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

    cout << "Closing the input files ... " << flush;
    for(unsigned int ifile=0; ifile<ifiles.size(); ifile++) {
        ifiles[ifile]->Close();
    }
    cout << "DONE" << endl;

    cout << "Closing the output file ... " << flush;
    ofile->Close();
    cout << "DONE" << endl;

    if(batch && close)
        gApplication->Terminate();
    return;
}

void PlotBasic() {
    TFile* ofile = TFile::Open((DefaultValues::getConfigPath()+"ZllGenNeutrinoWeightFile.root").c_str(),"READ");
    ((TH1D*)ofile->Get("ZllWeight_jets2_electron"))->Draw("P");
    ((TCanvas*)ofile->Get("validation/Canvas_jets2_electron"))->Draw();
}