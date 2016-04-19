//
// ROOT includes
//
#include "TROOT.h"
#include "TSystem.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TProfile3D.h"
#include "TRandom3.h"
#include "TString.h"
#include "TF1.h"
#include "TMath.h"
#include "TBenchmark.h"
#include "TRegexp.h"
#include "TKey.h"
#include "TList.h"
#include "TLegend.h"
#include "TLegendEntry.h"

//
// Standard Library Includes
//
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

TCanvas* draw_all_hist(TString objectName = "", TString scaleToHistPattern = "", bool removeSys = true,
                   bool doLogX = false, bool doLogY = false, bool doGrid = true, double xmin = 20,
                   double xmax = 1000) {

    TCanvas* c = new TCanvas(objectName,objectName);
    c->cd();

    TH1D* scaleToHist = 0;
    double scaleToHistInt = 0;
    if(!scaleToHistPattern.IsNull()) {
        TString reg = objectName+scaleToHistPattern;
        TRegexp re(reg,kTRUE);
        TIter next(gDirectory->GetListOfKeys());
        TKey* key(0);
        while ((key= (TKey*)next())) {
            TString st = key->GetName();
            if (st.Index(re) == kNPOS) continue;
            if(scaleToHist!=0) {
                cout << "ERROR::draw_all_hist More than one histogram matches the scaleToHistPattern (" << scaleToHistPattern << ")" << endl
                     << "\tAre you trying to match to " << key->GetName() << endl;
            }
            scaleToHist = (TH1D*)key->ReadObj();
        }
        if(scaleToHist==0) {
            cout << "ERROR::draw_all_hist A scaleToHistPattern was specified (" << scaleToHistPattern << "), but no histogram was found." << endl;
            return 0;
        }
        scaleToHistInt = scaleToHist->Integral();
    }

    TIter next = gDirectory->GetListOfKeys();
    TKey* key(0);
    int count = 0;
    TLegend* leg = new TLegend(0.65,0.6,0.8,1.0);
    leg->SetName("leg");
    TLegend* leg_ks = new TLegend(0.8,0.6,0.9,1.0);
    leg_ks->SetName("leg_ks");
    TLegend* leg_chi2NDF = new TLegend(0.9,0.6,1.0,1.0);
    leg_chi2NDF->SetName("leg_chi2NDF");
    leg->AddEntry((TObject*)0,"Sample Name","");
    leg_ks->AddEntry((TObject*)0,"K-S","");
    leg_chi2NDF->AddEntry((TObject*)0,"Chi2/NDF","");
    TH1D* h_zll_merge = 0;
    while(key = (TKey*)next()) {
        if(TString(key->GetName()).Contains("shape")) continue;
        if(TString(key->GetName()).Contains("HToZZ")) continue;
        if(TString(key->GetName()).Contains("HToWW")) continue;
        if(TString(key->GetName()).Contains("HToBB")) continue;
        //if(TString(key->GetName()).Contains("ZJets") && !TString(key->GetName()).Contains("ZJetsToLL")) continue;
        if(!objectName.IsNull() && !TString(key->GetName()).BeginsWith(objectName)) continue;
        if(removeSys && TString(key->GetName()).Contains("JESUp"))                  continue;
        if(removeSys && TString(key->GetName()).Contains("JESDown"))                continue;
        if(removeSys && TString(key->GetName()).Contains("matchingup"))             continue;
        if(removeSys && TString(key->GetName()).Contains("matchingdown"))           continue;
        if(removeSys && TString(key->GetName()).Contains("scaleup"))                continue;
        if(removeSys && TString(key->GetName()).Contains("scaledown"))              continue;
        cout << key->GetName() << endl;
        TH1D* h = (TH1D*)key->ReadObj();

        h->SetTitle(objectName);

        if(TString(key->GetName()).Contains("ZJetsToLL_M50")) {
            h->SetMarkerColor(kRed);
            h->SetLineColor(kRed);
            if(h_zll_merge) {
                //xs*lumi*BR*SF/init. ev.
                h->Scale((3387.60*19213.5*1.0*1.0/30209426.0));
                cout << "WARNING::draw_all_hist Adding the ZJetsToLL_M10To50 histogram to the ZJetsToLL_M50 histogram" << endl;
                h->Add(h_zll_merge);
            }
        }
        else if(TString(key->GetName()).Contains("ZJetsToLL_M10To50")) {
            h->SetMarkerColor(kMagenta);
            h->SetLineColor(kMagenta);
            h->Scale((5900.10*19213.5*1.0*1.5468132424)/37835275.0);
            h_zll_merge = h;
            count++;
            continue;
        }
        else {
            h->SetLineColor(h->GetMarkerColor());
        }

        if(scaleToHist!=0) {
            h->Scale(scaleToHistInt/h->Integral());
        }

        if(count==0 || (count==1 && h_zll_merge!=0)) {
            h->Draw();
            h->GetXaxis()->SetRangeUser(xmin,xmax);
            h->GetXaxis()->SetTitle(objectName);
            if(scaleToHistPattern.IsNull()) h->GetYaxis()->SetRangeUser(1,200000);
            else h->GetYaxis()->SetRangeUser(1,scaleToHist->GetMaximum()*1.1);
            h->GetYaxis()->SetTitle(TString("Scaled to histogram matching ")+scaleToHistPattern);
            if(doLogX) c->SetLogx(1);
            if(doLogY) c->SetLogy(1);
            if(doGrid) {
                gPad->SetGridx();
                gPad->SetGridy();
            }
        }
        else{
            h->Draw("same");
        }

        string name = h->GetName();
        int pos = name.find(string(objectName));
        int n = string(objectName).length();
        name.erase(pos,n);
        pos = name.find("_electron");
        if(pos!=std::string::npos) name.erase(pos,9);
        pos = name.find("_muon");
        if(pos!=std::string::npos) name.erase(pos,5);

        if(h_zll_merge!=0 && TString(h->GetName()).Contains("ZJetsToLL_M50"))
            leg->AddEntry(h,"ZJetsToLL","lep");
        else
            leg->AddEntry(h,name.c_str(),"lep");

        double chi2 = 0;
        int NDF = 0;
        int igood;
        double chi2NDF = scaleToHist->Chi2TestX(h,chi2,NDF,igood,"CHI2/NDF");
        leg_ks->AddEntry((TObject*)0,Form("%5.4g",scaleToHist->KolmogorovTest(h)),"");
        leg_chi2NDF->AddEntry((TObject*)0,Form("%5.4g",chi2/NDF),"");        

        count++;
    }
    leg->Draw("same");
    leg_ks->Draw("same");
    leg_chi2NDF->Draw("same");

    return c;
}

void save_all_canvases(string folder = "./", string format = ".eps") {
    TList* loc = (TList*)gROOT->GetListOfCanvases();
    TListIter itc(loc);
    TObject *o(0);
    while ((o = itc())) {
        TCanvas* c = (TCanvas*)o;
        string name = c->GetName();
        if(name[name.length()-1]=='_')
            name = name.substr(0,name.length()-1);
        c->SaveAs((folder+"/"+name+format).c_str());
    }
}

void DestroyCanvases() {
    TList* loc = (TList*)gROOT->GetListOfCanvases();
    TListIter itc(loc);
    TObject *o(0);
    while ((o = itc())) delete o;
}

void draw_all_standard(bool save = false, bool batch = false, bool doLogY = true) {
    if(batch)
        gROOT->SetBatch(kTRUE);
    
    draw_all_hist("LeptPt_","WJets_*",true,true,doLogY,true);
    draw_all_hist("MET_","WJets_*",true,true,doLogY,true);
    draw_all_hist("WmT_","WJets_*",true,true,doLogY,true);
    draw_all_hist("Ptjj_","WJets_*",true,true,doLogY,true);
    draw_all_hist("Jet1Pt_","WJets_*",true,true,doLogY,true);
    draw_all_hist("Jet2Pt_","WJets_*",true,true,doLogY,true);
    draw_all_hist("Ptlv_","WJets_*",true,true,doLogY,true);
    draw_all_hist("Mlv_","WJets_*",true,true,doLogY,true);
    draw_all_hist("Mlvjj_","WJets_*",true,true,doLogY,true);
    draw_all_hist("LeptPhi_","WJets_*",true,false,false,true,-TMath::Pi(),TMath::Pi());
    draw_all_hist("LeptEta_","WJets_*",true,false,false,true,-TMath::Pi(),TMath::Pi());
    draw_all_hist("Jet1Eta_","WJets_*",true,false,false,true,-TMath::Pi(),TMath::Pi());
    draw_all_hist("Jet1Phi_","WJets_*",true,false,false,true,-TMath::Pi(),TMath::Pi());
    draw_all_hist("DeltaPhi_LMET_","WJets_*",true,false,false,true,-TMath::Pi(),TMath::Pi());

    if(save) {
        save_all_canvases("./",".png");
        save_all_canvases("./",".eps");
    }
}

void makeSetOfStandardCanvases(bool save = true, bool batch = true) {
    if(batch)
      gROOT->SetBatch(kTRUE);

    TFile* ifile(0);
    TCanvas* c(0);
    TString jet[3] = {"jets2","jets3","jets4"};
    TString lepton[2] = {"electron","muon"};
    stringstream ss;

    for(int j=0; j<3; j++) {
        cout << "\t" << jet[j] << endl;
        for(int k=0; k<2; k++) {
            cout << "\t\t" << lepton[k] << endl;
            ifile = TFile::Open(Form("%s/%s/histos_%s_%s.root",jet[j].Data(),lepton[k].Data(),lepton[k].Data(),jet[j].Data()),"READ");
            draw_all_standard(false,true,false);
            save_all_canvases(Form("%s/%s/",jet[j].Data(),lepton[k].Data()),".png");
            DestroyCanvases();
        }
    }
}

void makeSetOfLeptonPtCanvases(bool batch = true) {
    if(batch)
        gROOT->SetBatch(kTRUE);

    TFile* ifile(0);
    TCanvas* c(0);
    TString jet[3] = {"jets2","jets3","jets4"};
    TString lepton[2] = {"electron","muon"};
    stringstream ss;

    for(int i=40; i<=70; i++) {
        cout << "HighPtLeptonPercentage = " << i << "%:" << endl;
        vector<TString> stats;
        for(int j=0; j<3; j++) {
            cout << "\t" << jet[j] << endl;
            for(int k=0; k<2; k++) {
                cout << "\t\t" << lepton[k] << endl;
                ifile = TFile::Open(Form("HighPtLeptonKept0%i/%s/%s/histos_%s_%s.root",i,jet[j].Data(),lepton[k].Data(),lepton[k].Data(),jet[j].Data()),"READ");
                c=draw_all_hist("LeptPt_","WJets_*",true,true,false,true);
                stats.push_back(((TLegendEntry*)((TList*)((TLegend*)c->FindObject("leg_ks"))->GetListOfPrimitives())->At(1))->GetLabel());
                stats.push_back(((TLegendEntry*)((TList*)((TLegend*)c->FindObject("leg_chi2NDF"))->GetListOfPrimitives())->At(1))->GetLabel());
                c=draw_all_hist("MET_","WJets_*",true,true,false,true);
                save_all_canvases(Form("HighPtLeptonKept0%i/%s/%s/",i,jet[j].Data(),lepton[k].Data()),".png");
                DestroyCanvases();
            }
        }
        ss << i << ", ";
        for(unsigned int s=0; s<stats.size(); s++) {
            ss << stats[s];
            if(s<stats.size()-1) ss << ", ";
        }
        ss << endl;
    }
    cout << ss.str() << endl;
}


