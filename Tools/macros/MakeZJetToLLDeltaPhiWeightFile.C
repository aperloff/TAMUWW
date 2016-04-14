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
#include "TAMUWW/SpecialTools/interface/ProgressBar.hh"
#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//  Declare Global Variables
////////////////////////////////////////////////////////////////////////////////
const double ZMASS = 91.1876; //GeV
const double WMASS = 80.3980; //GeV

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void MakeZJetsToLLDeltaPhiWeightFile() {
    // Trying to speed up the code
    gEnv->SetValue("TFile.AsyncPrefetching", 1);

    Style* st = new Style();
    st->setTDRStyle();
    //TGaxis::SetMaxDigits(3);

    cout << "Opening output file ... " << flush;
    TFile* ofile = TFile::Open("ZllDeltaPhiWeights.root","RECREATE");
    ofile->mkdir("validation");
    cout << "DONE" << endl;

    const int nJetBins = 3;
    const int nLeptonBins = 2;
    const int nProcesses = 2;
    string jetBins[nJetBins] = {"jets2","jets3","jets4"};
    string leptonBins[nLeptonBins] = {"electron","muon"};
    string processes[nProcesses] = {"WJets","ZJetsToLL"};

    TTree *jets2p(0);
    string hname, htitle;
    map<string, TH1D*> histograms;
    map<string, TCanvas*> canvases;
    TLorentzVector lepton;
    TLorentzVector neutrino;
    int *l_pdgid = new int(0);
    int *nu_pdgid = new int(0);

    cout << "Creating the histograms ... " << flush;
    for(int jbin=0; jbin<nJetBins; jbin++) {
        for(int lbin=0; lbin<nLeptonBins; lbin++) {
            for(int iprocess=0; iprocess<nProcesses; iprocess++) {
                hname = processes[iprocess]+"_"+jetBins[jbin]+"_"+leptonBins[lbin];
                htitle = "#Delta#phi(l,#nu) (";
                htitle+= processes[iprocess]=="WJets"?"GEN)":"RECO)";
                histograms[hname] = new TH1D(hname.c_str(),htitle.c_str(),62,-TMath::Pi(),TMath::Pi());
                histograms[hname]->Sumw2();
            }
        }
    }
    cout << "DONE" << endl;

    cout << "Opening input file ... " << flush;
    TFile* _file0 = TFile::Open("$SMEInput/WJets.root","READ");
    //TFile* _file1 = TFile::Open("$SMEInput/DYJetsToLL_M-50.root","READ");
    TFile* _file1 = TFile::Open("$SMEInput/DYJetsToLL_M-50_lowerSubLeadingJet.root","READ");
    cout << "DONE" << endl;

    cout << "Filling the input histograms ... " << endl;
    EventNtuple* ntuple = new EventNtuple();
    jets2p = (TTree*)_file0->Get("PS/jets2p");
    jets2p->SetBranchAddress("EvtNtuple", &ntuple);
    jets2p->SetBranchStatus("*",0);
    jets2p->SetBranchStatus("jLV*",1);
    jets2p->SetBranchStatus("lLV*",1);
    jets2p->SetBranchStatus("genParticleCollection*",1);
    int nEntries = jets2p->GetEntries();
    for(unsigned int ientry=0; ientry<jets2p->GetEntries(); ientry++) {
        ProgressBar::loadbar2(ientry+1,nEntries);
        jets2p->GetEntry(ientry);
        lepton = ntuple->getGenVorDaughter(EventNtuple::LEPTON, 24, 0, l_pdgid, false);
        neutrino = ntuple->getGenVorDaughter(EventNtuple::NEUTRINO, 24, 0, nu_pdgid, false);
        //if(lepton.Pt()<28 || abs(lepton.Eta())>2.0 || abs(*l_pdgid)!=11){
        //if(DEFS::getLeptonCatString(ntuple->lLV[0].leptonCat)!=leptonBins[lbin] ||
        //   (ntuple->jLV.size() != 1 && jetBin == DEFS::jet1) ||
        //   (ntuple->jLV.size() != 2 && jetBin == DEFS::jets2) ||
        //   (ntuple->jLV.size() != 3 && jetBin == DEFS::jets3) ||
        //   (ntuple->jLV.size() < 4 && jetBin == DEFS::jets4)) {
        //    continue;
        //}
        hname = "WJets_"+DEFS::getJetBinString(DEFS::getJetBin(ntuple->jLV.size()))+"_"+DEFS::getLeptonCatString(ntuple->lLV[0].leptonCat);
        if(histograms.find(hname)==histograms.end()) {
            cout << "ERROR::MakeZJetsToLLDeltaPhiWeightFile Cannot find the histogram named " << hname << endl;
            std::terminate();
        }
        if(lepton.Pt()>0 && neutrino.Pt()>0) {
            histograms.find(hname)->second->Fill(lepton.DeltaPhi(neutrino));
        }
    }
    cout << endl;

    jets2p = (TTree*)_file1->Get("PS/jets2p");
    jets2p->SetBranchAddress("EvtNtuple", &ntuple);
    jets2p->SetBranchStatus("*",0);
    jets2p->SetBranchStatus("jLV*",1);
    jets2p->SetBranchStatus("lLV*",1);
    nEntries = jets2p->GetEntries();
    for(unsigned int ientry=0; ientry<jets2p->GetEntries(); ientry++) {
        ProgressBar::loadbar2(ientry+1,nEntries);
        jets2p->GetEntry(ientry);
        //lepton = ntuple->getGenVorDaughter(EventNtuple::LEPTON, 23, 0, l_pdgid, false);
        //neutrino = ntuple->getGenVorDaughter(EventNtuple::LEPTON, 23, 1, nu_pdgid, false);
        //if(lepton.Pt()<28 || abs(lepton.Eta())>2.0 || abs(*l_pdgid)!=11){
        //if(DEFS::getLeptonCatString(ntuple->lLV[0].leptonCat)!=leptonBins[lbin] ||
        //   (ntuple->jLV.size() != 1 && jetBin == DEFS::jet1) ||
        //   (ntuple->jLV.size() != 2 && jetBin == DEFS::jets2) ||
        //   (ntuple->jLV.size() != 3 && jetBin == DEFS::jets3) ||
        //   (ntuple->jLV.size() < 4 && jetBin == DEFS::jets4)) {
        //    continue;
        //}
        //ZJetsToLL->Fill(lepton.DeltaPhi(neutrino));
        hname = "ZJetsToLL_"+DEFS::getJetBinString(DEFS::getJetBin(ntuple->jLV.size()))+"_"+DEFS::getLeptonCatString(ntuple->lLV[0].leptonCat);
        if(histograms.find(hname)==histograms.end()) {
            cout << "ERROR::MakeZJetsToLLDeltaPhiWeightFile Cannot find the histogram named " << hname << endl;
            std::terminate();
        }
        ntuple->lLV[0]*=(WMASS/ZMASS);
        ntuple->lLV[1]*=(WMASS/ZMASS);
        histograms.find(hname)->second->Fill(ntuple->lLV[0].DeltaPhi(ntuple->lLV[1]));
    }

    cout << "Checking filled histograms ... " << flush;
    for(map<string,TH1D*>::iterator it=histograms.begin(); it!=histograms.end(); it++) {
        if(it->second->GetEntries()==0) {
            cout << endl << "ERROR::MakeZJetsToLLDeltaPhiWeightFile The histogtam " << it->first << " has zero entries." << endl;
            std::terminate();
        }
    }
    cout << "DONE" << endl;

    cout << "Scaling ZJetToLL to WJets, making the weights, and writing the weights ... " << flush;
    for(int jbin=0; jbin<nJetBins; jbin++) {
        for(int lbin=0; lbin<nLeptonBins; lbin++) {
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

    cout << "Making and writing a validation canvas ... " << flush;
    for(int jbin=0; jbin<nJetBins; jbin++) {
        for(int lbin=0; lbin<nLeptonBins; lbin++) {
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
            hdown->GetYaxis()->SetRangeUser(0,4);
            hup->GetYaxis()->SetTitle("Events");
            hdown->GetYaxis()->SetTitle("Weight");
            canvases[cname] = st->tdrDiCanvas(cname.c_str(),hup,hdown);
            canvases[cname]->cd(1);
            st->tdrDraw(histograms[w_name],"",kFullCircle,kGreen,kSolid,kGreen);
            st->tdrDraw(histograms[zll_name],"",kFullCircle,kRed,kSolid,kRed);
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

    cout << "Closing the input files ... " << flush;
    _file0->Close();
    _file1->Close();
    cout << "DONE" << endl;

    cout << "Closing the output file ... " << flush;
    ofile->Close();
    cout << "DONE" << endl;
}
