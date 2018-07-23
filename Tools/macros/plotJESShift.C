#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGaxis.h"
#include "TH1D.h"
#include "TLegend.h"

#include <vector>
#include <string>
#include <algorithm>

#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"
#include "TAMUWW/Tools/interface/Style.hh"
#include "TAMUWW/Tools/interface/PUreweight.hh"

using namespace std;

/*
Use:
.L /uscms_data/d2/aperloff/YOURWORKINGAREA/MatrixElement/gitty/CMSSW_5_3_22_patch1/src/TAMUWW/Tools/macros/plotJESShift.C+
*/

void plotJESShift(string ifilename="", string process="WJets", string variable="KinMEBDT",
                  string lepton="electron",string systematic="CMS_scale_j_shape",bool normalize = true){
    Style* st = new Style();
    st->setTDRStyle();
    TGaxis::SetMaxDigits(3);

    //Get the input histograms
    if (ifilename.empty())
        ifilename = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/jets2/electron/histos_electron_jets2_SysNames.root";
    TFile* ifile = TFile::Open(ifilename.c_str());
    ifile->cd();
    string nominal_name = variable+"_"+process+"_"+lepton;
    TH1D* nominal = (TH1D*)((TH1D*)gDirectory->Get(nominal_name.c_str()))->Clone("nominal");
    string shape_name = nominal_name+"_"+systematic;
    TH1D* up = (TH1D*)((TH1D*)gDirectory->Get((shape_name+"Up").c_str()))->Clone("up");
    TH1D* down = (TH1D*)((TH1D*)gDirectory->Get((shape_name+"Down").c_str()))->Clone("down");
    nominal->SetDirectory(0);
    up->SetDirectory(0);
    down->SetDirectory(0);

    //Scale the histograms to the nominal integral
    float nominal_int = nominal->Integral();
    if(normalize) {
        up->Scale(nominal_int/up->Integral());
        down->Scale(nominal_int/down->Integral());
    }

    TH1D* frame;
    TCanvas* can;
    TPad* padMain;
    TLegend* leg;

   //Setup the first frame and create the canvas
    frame = new TH1D();
    frame->GetXaxis()->SetLimits(-1, 1);
    vector<double> maximums;
    maximums.push_back(nominal->GetBinContent(nominal->GetMaximumBin()));
    maximums.push_back(up->GetBinContent(up->GetMaximumBin()));
    maximums.push_back(down->GetBinContent(down->GetMaximumBin()));
    double multiplier = 1.25;
    frame->GetYaxis()->SetRangeUser(0,*std::max_element(maximums.begin(),maximums.end())*multiplier);
    frame->GetXaxis()->SetTitle(variable.c_str());
    frame->GetYaxis()->SetTitle("Events");
    can = st->tdrCanvas("JES",frame,12,11,true);
    can->cd();
    padMain = (TPad*)can->GetPad(0);
    //padMain->SetLogy();
    //Draw the pileup distributions
    st->tdrDraw(nominal,"hist",kFullCircle,kBlack,kSolid,kBlack,0,kNone);
    st->tdrDraw(up,"hist",kFullCircle,kRed,kSolid,kRed,0,kNone);
    st->tdrDraw(down,"hist",kFullCircle,kBlue,kSolid,kBlue,0,kNone);
    //Create and draw the legend
    leg = st->tdrLeg(0.65,0.74,0.88,0.91);
    leg->AddEntry(nominal,"Nominal","l");
    leg->AddEntry(up,"+1#sigma","l");
    leg->AddEntry(down,"-1#sigma","l");
    leg->Draw("same");
    //Save the canvas
    can->SaveAs(("Shape_"+variable+"_"+systematic+"_"+process+".eps").c_str());
    can->SaveAs(("Shape_"+variable+"_"+systematic+"_"+process+".pdf").c_str());
    can->SaveAs(("Shape_"+variable+"_"+systematic+"_"+process+".png").c_str());
    can->SaveAs(("Shape_"+variable+"_"+systematic+"_"+process+".C").c_str());
    can->SaveAs(("Shape_"+variable+"_"+systematic+"_"+process+".root").c_str());
}