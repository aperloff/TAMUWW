//
// ROOT includes
//
#include "TROOT.h"
#include "TSystem.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile3D.h"
#include "TRandom3.h"
#include "TString.h"
#include "TF1.h"
#include "TMath.h"
#include "TBenchmark.h"
#include "TRegexp.h"

//
// Standard Library Includes
//
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>


void draw_all_hist(TString objectName = "", TString scaleToHistPattern = "", bool removeSys = true,
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
            return;
        }
        scaleToHistInt = scaleToHist->Integral();
    }

    next = gDirectory->GetListOfKeys();
    TKey* key(0);
    int count = 0;
    TLegend* leg = new TLegend(0.65,0.6,0.8,1.0);
    TLegend* leg_ks = new TLegend(0.8,0.6,0.9,1.0);
    TLegend* leg_chi2NDF = new TLegend(0.9,0.6,1.0,1.0);
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
