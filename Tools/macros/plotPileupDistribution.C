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

void plotPileupDistribution(bool drawNPI = true, bool drawWeights = true){
   Style* st = new Style();
   st->setTDRStyle();
   TGaxis::SetMaxDigits(3);

   //Get the input histograms
   string dataname = DefaultValues::getConfigPath()+"pileup12_noTrig_minBiasXsec69400_coarseBinning_withAdditionalNPVHist.root";
   TFile* ifile_data = TFile::Open(dataname.c_str());
   ifile_data->cd();
   TH1D* data = (TH1D*)((TH1D*)gDirectory->Get("pileup_noTrig"))->Clone("data");
   data->SetDirectory(0);
   string MCname = DefaultValues::getConfigPath()+"TPUDistributions.root";
   TFile* ifile_mc = TFile::Open(MCname.c_str());
   ifile_mc->cd();
   TH1D* mc = (TH1D*)((TH1D*)gDirectory->Get("TPUDist_WJets"))->Clone("mc");
   mc->SetDirectory(0);

   //Setup the pileup weights to be used later
   PUreweight* puweight = new PUreweight(dataname,MCname,"pileup_noTrig","TPUDist_WJets",make_pair(1,10));

   //Scale the histograms to be drawn
   float data_int = data->Integral();
   data->Scale(1.0/data_int);
   mc->Rebin(10);
   float mc_int = mc->Integral();
   mc->Scale(1.0/mc_int);

   TH1D* frame;
   TCanvas* can;
   TPad* padMain;
   TLegend* leg;
   if(drawNPI) {
      //Setup the first frame and create the canvas
      frame = new TH1D();
      frame->GetXaxis()->SetLimits(0, 60);
      vector<double> maximums;
      maximums.push_back(data->GetBinContent(data->GetMaximumBin()));
      maximums.push_back(mc->GetBinContent(mc->GetMaximumBin()));
      double multiplier = 100.0;//1.45;//1.25;
      frame->GetYaxis()->SetRangeUser(1e-9,*std::max_element(maximums.begin(),maximums.end())*multiplier);
      frame->GetXaxis()->SetTitle("Number of pileup interactions");
      frame->GetYaxis()->SetTitle("A.U.");
      can = st->tdrCanvas("NPU",frame,2,11,true);
      can->cd();
      padMain = (TPad*)can->GetPad(0);
      padMain->SetLogy();

      //Draw the pileup distributions
      st->tdrDraw(data,"P",kFullCircle,kBlack,kNone,kNone,0,kNone);
      st->tdrDraw(mc,"P",kFullTriangleUp,kRed+2,kNone,kNone,0,kNone);

      //Create and draw the legend
      leg = st->tdrLeg(0.65,0.75,0.88,0.92);
      leg->AddEntry(data,"Data","p");
      leg->AddEntry(mc,"Simulation","p");
      leg->Draw("same");

      //Save the canvas
      can->SaveAs((DefaultValues::getConfigPath()+"PileupDistribution.eps").c_str());
      can->SaveAs((DefaultValues::getConfigPath()+"PileupDistribution.pdf").c_str());
      can->SaveAs((DefaultValues::getConfigPath()+"PileupDistribution.png").c_str());
      can->SaveAs((DefaultValues::getConfigPath()+"PileupDistribution.C").c_str());
      can->SaveAs((DefaultValues::getConfigPath()+"PileupDistribution.root").c_str());
   }
   if(drawWeights){
      //Make the weights histogram
      TH1D* weights = new TH1D("weights","weights",60,0,60);
      for(unsigned int i=0; i<61; i++) {
         weights->Fill(i+0.5,puweight->getWeight(i));
      }

      //Setup the first frame and create the canvas
      frame = new TH1D();
      frame->GetXaxis()->SetLimits(0, 60);
      vector<double> maximums;
      maximums.push_back(weights->GetBinContent(weights->GetMaximumBin()));
      double multiplier = 1.1;
      frame->GetYaxis()->SetRangeUser(0,*std::max_element(maximums.begin(),maximums.end())*multiplier);
      frame->GetXaxis()->SetTitle("Number of pileup interactions");
      frame->GetYaxis()->SetTitle("Pileup weight");
      can = st->tdrCanvas("PileupWeight",frame,2,11,true);
      can->cd();

      //Draw the weight histogram
      weights->SetLineWidth(2);
      st->tdrDraw(weights,"hist",kNone,kNone,kSolid,kRed+2,1001,kRed-7);

      //Save the canvas
      can->SaveAs((DefaultValues::getConfigPath()+"PileupWeights.eps").c_str());
      can->SaveAs((DefaultValues::getConfigPath()+"PileupWeights.pdf").c_str());
      can->SaveAs((DefaultValues::getConfigPath()+"PileupWeights.png").c_str());
      can->SaveAs((DefaultValues::getConfigPath()+"PileupWeights.C").c_str());
      can->SaveAs((DefaultValues::getConfigPath()+"PileupWeights.root").c_str());
   }
}

void plotPileupWeight() {



   
}
