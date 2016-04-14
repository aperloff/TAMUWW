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
#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void MakeZJetsToLLLepResSmearingFactor(int maxEvts = 0, bool batch = true, bool close = false, bool save = false) {
    if(batch)
        gROOT->SetBatch(kTRUE);

    // Trying to speed up the code
    gEnv->SetValue("TFile.AsyncPrefetching", 1);

    Style* st = new Style();
    st->setTDRStyle();
    //TGaxis::SetMaxDigits(3);

    cout << "Opening output file ... " << flush;
    TFile* ofile = TFile::Open("ZllLepResSmearing.root","RECREATE");
    ofile->mkdir("validation");
    cout << "DONE" << endl;

    const int nJetBins = 3;
    const int nLeptonBins = 2;
    const int nProcesses = 2;
    const int nAxis = 2;
    string jetBins[nJetBins] = {"jets2","jets3","jets4"};
    string leptonBins[nLeptonBins] = {"electron","muon"};
    string processes[nProcesses] = {"WJets","ZJetsToLL"};
    string processes_extra[nProcesses+1] = {"WJets","ZJetsToLL","Difference"};
    string axis[nAxis] = {"x","y"};

    TTree *jets2p(0);
    string hname, htitle, gname;
    map<string, TH1*> histograms;
    map<string, TGraph*> graphs;
    map<string, TCanvas*> canvases;
    map<string, TLegend*> legends;
    TLorentzVector lepton;
    TLorentzVector neutrino;
    int *nu_pdgid = new int(0);

    cout << "Creating the histograms ... " << flush;
    for(int jbin=0; jbin<nJetBins; jbin++) {
        for(int lbin=0; lbin<nLeptonBins; lbin++) {
            for(int iprocess=0; iprocess<nProcesses; iprocess++) {
                for(int iAxis=0; iAxis<nAxis; iAxis++) {
                    hname = processes[iprocess]+"_"+jetBins[jbin]+"_"+leptonBins[lbin]+"_p"+axis[iAxis];
                    htitle = "P_{"+axis[iAxis]+"}^{MET}/P_{"+axis[iAxis]+"}^{#nu} vs. P_{"+axis[iAxis]+"}^{#nu}";
                    histograms[hname] = new TH2D(hname.c_str(),htitle.c_str(),100,-500,500,10000,-500,500);
                    histograms[hname]->Sumw2();
                }
            }
        }
    }
    cout << "DONE" << endl;

    cout << "Opening input file ... " << flush;
    TFile* _file0 = TFile::Open("$SMEInput/WJets.root","READ");
    TFile* _file1 = TFile::Open("$SMEInput/WlnuJets_M-50_HighPtLeptonKept062.root","READ");
    cout << "DONE" << endl;

    cout << "Filling the input histograms ... " << endl;
    EventNtuple* ntuple = new EventNtuple();
    jets2p = (TTree*)_file0->Get("PS/jets2p");
    jets2p->SetBranchAddress("EvtNtuple", &ntuple);
    jets2p->SetBranchStatus("*",0);
    jets2p->SetBranchStatus("METLV*",1);
    jets2p->SetBranchStatus("jLV",1);
    jets2p->SetBranchStatus("lLV.leptonCat",1);
    jets2p->SetBranchStatus("genParticleCollection*",1);
    unsigned int nEntries = (maxEvts>0 && maxEvts<jets2p->GetEntries()) ? maxEvts : jets2p->GetEntries();
    for(unsigned int ientry=0; ientry<nEntries; ientry++) {
        ProgressBar::loadbar2(ientry+1,nEntries);
        jets2p->GetEntry(ientry);
        neutrino = ntuple->getGenVorDaughter(EventNtuple::NEUTRINO, 24, 0, nu_pdgid, false);

        for(int iAxis=0; iAxis<nAxis; iAxis++) {
            hname = "WJets_"+DEFS::getJetBinString(DEFS::getJetBin(ntuple->jLV.size()))+"_"+DEFS::getLeptonCatString(ntuple->lLV[0].leptonCat)+"_p"+axis[iAxis];
            if(histograms.find(hname)==histograms.end()) {
                cout << "ERROR::MakeZJetsToLLLepResSmearingFactor Cannot find the histogram named " << hname << endl;
                std::terminate();
            }
            if(neutrino.Pt()>0) {
                if(axis[iAxis]=="x")
                    //histograms.find(hname)->second->Fill(neutrino.Px(),ntuple->METLV[0].Px()/neutrino.Px());
                    histograms.find(hname)->second->Fill(neutrino.Px(),ntuple->METLV[0].Px());
                else if(axis[iAxis]=="y")
                    //histograms.find(hname)->second->Fill(neutrino.Py(),ntuple->METLV[0].Py()/neutrino.Py());
                    histograms.find(hname)->second->Fill(neutrino.Py(),ntuple->METLV[0].Py());
            }
        }
    }
    cout << endl;

    jets2p = (TTree*)_file1->Get("PS/jets2p");
    jets2p->SetBranchAddress("EvtNtuple", &ntuple);
    jets2p->SetBranchStatus("*",0);
    jets2p->SetBranchStatus("METLV*",1);
    jets2p->SetBranchStatus("jLV",1);
    jets2p->SetBranchStatus("lLV.leptonCat",1);
    //jets2p->SetBranchStatus("genParticleCollection*",1);
    nEntries = (maxEvts>0 && maxEvts<jets2p->GetEntries()) ? maxEvts : jets2p->GetEntries();
    for(unsigned int ientry=0; ientry<nEntries; ientry++) {
        ProgressBar::loadbar2(ientry+1,nEntries);
        jets2p->GetEntry(ientry);
        //neutrino = ntuple->getGenVorDaughter(EventNtuple::NEUTRINO, 24, 0, nu_pdgid, false);

        for(int iAxis=0; iAxis<nAxis; iAxis++) {
            hname = "ZJetsToLL_"+DEFS::getJetBinString(DEFS::getJetBin(ntuple->jLV.size()))+"_"+DEFS::getLeptonCatString(ntuple->lLV[0].leptonCat)+"_p"+axis[iAxis];
            if(histograms.find(hname)==histograms.end()) {
                cout << "ERROR::MakeZJetsToLLLepResSmearingFactor Cannot find the histogram named " << hname << endl;
                std::terminate();
            }
            //if(neutrino.Pt()>0) {
                if(axis[iAxis]=="x")
                    //histograms.find(hname)->second->Fill(neutrino.Px(),ntuple->METLV[0].refLV.Px()/neutrino.Px());
                    //histograms.find(hname)->second->Fill(ntuple->METLV[0].refLV.Px(),ntuple->METLV[0].rawLV.Px()/ntuple->METLV[0].refLV.Px());
                    //histograms.find(hname)->second->Fill(ntuple->METLV[0].refLV.Px(),ntuple->METLV[0].rawLV.Px());
                    histograms.find(hname)->second->Fill(ntuple->METLV[0].refLV.Px(),ntuple->METLV[0].Px());
                else if(axis[iAxis]=="y")
                    //histograms.find(hname)->second->Fill(neutrino.Py(),ntuple->METLV[0].refLV.Py()/neutrino.Py());
                    //histograms.find(hname)->second->Fill(ntuple->METLV[0].refLV.Py(),ntuple->METLV[0].rawLV.Py()/ntuple->METLV[0].refLV.Py());
                    //histograms.find(hname)->second->Fill(ntuple->METLV[0].refLV.Py(),ntuple->METLV[0].rawLV.Py());
                    histograms.find(hname)->second->Fill(ntuple->METLV[0].refLV.Py(),ntuple->METLV[0].Py());
            //}
        }
    }

    cout << "Checking filled histograms ... " << flush;
    for(map<string,TH1*>::iterator it=histograms.begin(); it!=histograms.end(); it++) {
        if(it->second->GetEntries()==0) {
            cout << endl << "ERROR::MakeZJetsToLLLepResSmearingFactor The histogtam " << it->first << " has zero entries." << endl;
            std::terminate();
        }
    }
    cout << "DONE" << endl;

    cout << "Make the difference histogram ... " << flush;
    for(int jbin=0; jbin<nJetBins; jbin++) {
        for(int lbin=0; lbin<nLeptonBins; lbin++) {
            for(int iAxis=0; iAxis<nAxis; iAxis++) {
                string wjets_name = "WJets_"+jetBins[jbin]+"_"+leptonBins[lbin]+"_p"+axis[iAxis];
                string zjets_name = "ZJetsToLL_"+jetBins[jbin]+"_"+leptonBins[lbin]+"_p"+axis[iAxis];
                hname = "Difference_"+jetBins[jbin]+"_"+leptonBins[lbin]+"_p"+axis[iAxis];
                histograms[hname] = (TH2D*)histograms[wjets_name]->Clone(hname.c_str());
                histograms[hname]->Scale(dynamic_cast<TH2D*>(histograms[zjets_name])->Integral()/histograms[hname]->Integral());
                histograms[hname]->Add(dynamic_cast<TH2D*>(histograms[zjets_name]),-1.0);
            }
        }
    }
    cout << "DONE" << endl;

    cout << "Making and writing a validation canvas ... " << flush;
    TH1D *hup(0), *hdown(0);
    TLine *line(0);
    for(int jbin=0; jbin<nJetBins; jbin++) {
        for(int lbin=0; lbin<nLeptonBins; lbin++) {
            for(int iprocess=0; iprocess<nProcesses+1; iprocess++) {
                for(int iAxis=0; iAxis<nAxis; iAxis++) {
                    hname = processes_extra[iprocess]+"_"+jetBins[jbin]+"_"+leptonBins[lbin]+"_p"+axis[iAxis];
                    ofile->cd("validation");
                    string cname = "Canvas_"+hname;

                    hup = new TH1D("hup","hup",100,-500,500);
                    hup->GetXaxis()->SetLimits(-500.0,500.0);
                    hup->GetYaxis()->SetRangeUser(-500.0,500.0);
                    if(processes_extra[iprocess].find("ZJets")!=string::npos)
                        hup->GetYaxis()->SetTitle(("P_{"+axis[iAxis]+"}^{MET} (GeV)").c_str());
                        //hup->GetYaxis()->SetTitle(("P_{"+axis[iAxis]+"}^{MET}/P_{"+axis[iAxis]+"}^{#nu}").c_str());
                    else if(processes_extra[iprocess].find("Difference")!=string::npos)
                        hup->GetYaxis()->SetTitle(("P_{"+axis[iAxis]+",MET}^{WJets} - P_{"+axis[iAxis]+",MET}^{ZJetsToLL} (GeV)").c_str());
                    else
                        hup->GetYaxis()->SetTitle(("P_{"+axis[iAxis]+"}^{MET} (GeV)").c_str());
                        //hup->GetYaxis()->SetTitle(("P_{"+axis[iAxis]+"}^{MET}/P_{"+axis[iAxis]+"}^{#nu}").c_str());
                    hdown = new TH1D("hdown","hdown",100,-500,500);
                    hdown->GetXaxis()->SetLimits(-500.0,500.0);
                    hdown->GetYaxis()->SetRangeUser(-500.0,500.0);
                    if(processes_extra[iprocess].find("ZJets")!=string::npos)
                        hdown->GetXaxis()->SetTitle(("P_{"+axis[iAxis]+"}^{l (RECO)} (GeV)").c_str());
                    else if(processes_extra[iprocess].find("Difference")!=string::npos)
                        hdown->GetXaxis()->SetTitle(("P_{"+axis[iAxis]+"}^{Mixed} (GeV)").c_str());
                    else
                        hdown->GetXaxis()->SetTitle(("P_{"+axis[iAxis]+"}^{#nu (GEN)} (GeV)").c_str());
                    hdown->GetYaxis()->SetTitle("Prof+Res");
                    canvases[cname] = st->tdrDiCanvas(cname.c_str(),hup,hdown,2,0); //8TeV, out-of-frame

                    //Setting a smoother palette
                    const Int_t NRGBs = 5;
                    const Int_t NCont = 104;
                    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
                    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
                    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
                    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
                    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
                    gStyle->SetNumberContours(NCont);

                    canvases[cname]->cd(1);
                    st->tdrDraw(histograms[hname],"colz",kNone,kNone);
                    canvases[cname]->cd(2);
                    line = new TLine(hdown->GetXaxis()->GetXmin(),1,hdown->GetXaxis()->GetXmax(),1);
                    line->SetLineStyle(2);
                    //ProfileX options:
                    //s: standard deviation of Y, S(Y)
                    //i: standard error on the mean, S(Y)/SQRT(N). When the standard deviation in Y is zero the error a
                    //   standard deviation = 1/SQRT(12) is assumed and the error is 1./SQRT(12*N). This approximation
                    //   assumes that the Y values are integer (e.g. ADC counts) and have an implicit uncertainty of y +/- 0.5.
                    //   With the assumption that the probability that y takes any value between y-0.5 and y+0.5 is uniform,
                    //   its standard error is 1/SQRT(12)
                    histograms[hname+"_pfx"] = dynamic_cast<TH2D*>(histograms[hname])->ProfileX((hname+"_pfx").c_str(),1,-1,"s");
                    st->tdrDraw(histograms[hname+"_pfx"],"",kFullCircle,kBlack,kSolid,kBlue);

                    gname = "Resolution_"+hname;
                    htitle = "Resolution vs. P_{"+axis[iAxis]+"}^{#nu}";
                    histograms[gname] = new TH1D(gname.c_str(),htitle.c_str(),histograms[hname+"_pfx"]->GetNbinsX(),
                                                 histograms[hname+"_pfx"]->GetXaxis()->GetXmin(),histograms[hname+"_pfx"]->GetXaxis()->GetXmax());
                    for(int ibin=1; ibin<=histograms[hname+"_pfx"]->GetNbinsX(); ibin++) {
                        histograms[gname]->SetBinContent(ibin,histograms[hname+"_pfx"]->GetBinError(ibin));
                    }
                    st->tdrDraw(histograms[gname],"P",kFullCircle,kRed,kNone,kNone,kNone,kNone);

                    legends["legend_"+hname] = st->tdrLeg(0.75,0.75,0.89,0.89);
                    legends["legend_"+hname]->AddEntry(histograms[hname+"_pfx"],"ProfileX","ep");
                    legends["legend_"+hname]->AddEntry(histograms[gname],"Resolution","p");
                    legends["legend_"+hname]->Draw("same");

                    canvases[cname]->Write();
                }
            }
        }
    }
    cout << "DONE" << endl;

    cout << "Making the smearing factors ... " << flush;
    ofile->cd();
    for(int jbin=0; jbin<nJetBins; jbin++) {
        for(int lbin=0; lbin<nLeptonBins; lbin++) {
            for(int iAxis=0; iAxis<nAxis; iAxis++) {
                string wjets_name = "Resolution_WJets_"+jetBins[jbin]+"_"+leptonBins[lbin]+"_p"+axis[iAxis];
                string zjets_name = "Resolution_ZJetsToLL_"+jetBins[jbin]+"_"+leptonBins[lbin]+"_p"+axis[iAxis];
                gname = "Resolution_"+jetBins[jbin]+"_"+leptonBins[lbin]+"_p"+axis[iAxis];
                htitle = "Resolution (#nu) vs. P_{"+axis[iAxis]+"}^{#nu}";

                histograms[gname] = new TH1D(gname.c_str(),htitle.c_str(),histograms[wjets_name]->GetNbinsX(),
                                             histograms[wjets_name]->GetXaxis()->GetXmin(),histograms[wjets_name]->GetXaxis()->GetXmax());

                double wjets_y=0, zjets_y=0, smearing_factor=0;
                for(int ibin=1; ibin<=histograms[wjets_name]->GetNbinsX(); ibin++) {
                    wjets_y = histograms[wjets_name]->GetBinContent(ibin);
                    zjets_y = histograms[zjets_name]->GetBinContent(ibin);
                    smearing_factor = sqrt(pow(wjets_y,2)-pow(zjets_y,2));
                    if(!TMath::IsNaN(smearing_factor))
                        histograms[gname]->SetBinContent(ibin,smearing_factor);
                }
                histograms[gname]->Write();
            }
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
    _file0->Close();
    _file1->Close();
    cout << "DONE" << endl;

    cout << "Closing the output file ... " << flush;
    ofile->Close();
    cout << "DONE" << endl;

    if(batch && close)
        gApplication->Terminate();
    return;
}

void PlotBasic() {
    TFile* ofile = TFile::Open("ZllLepResSmearing.root","READ");
    ((TH1D*)ofile->Get("Resolution_jets2_electron_px"))->Draw("P");
    ((TCanvas*)ofile->Get("validation/Canvas_WJets_jets2_electron_px"))->Draw();
    ((TCanvas*)ofile->Get("validation/Canvas_ZJetsToLL_jets2_electron_px"))->Draw();
    ((TCanvas*)ofile->Get("validation/Canvas_Difference_jets2_electron_px"))->Draw();
}