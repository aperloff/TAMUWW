// Our libraries
#include "TAMUWW/Tools/interface/Plots.hh"

// ROOT libraries
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TVirtualPad.h"
#include "TH1.h"
#include "THStack.h"
#include "TChain.h"
#include "TLegend.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TMath.h"

// C++ libraries
#include <iostream>
#include <sstream>
#include <map>
#include <string>
#include <vector>

using namespace std;

// ##################################################
// ################### PLOT CLASS ###################
// ##################################################

// Create a new histo here
void Plot::prepareToFillProcess(TString suffix, TString groupName)
{
   TString n = templateHisto->GetTitle();
   TH1 * clone = (TH1*) templateHisto->Clone(n+suffix);
   clone->SetTitle(groupName);
   clone->Sumw2();
   histos.push_back(clone);
}

// Fill the last histo here
void Plot::Fill(double x, double w){
   if (histos.size()>0) 
      histos.back()->Fill(x,w);
   else
      cout<<"ERROR in Plot named "<<templateHisto->GetName()
	  <<" Fill(..) called without calling prepareToFillProcess first"<<endl; 
}

// Fill the last histo here
void Plot::Fill(double x, double y, double w){
   if (histos.size()>0) 
      ((TH2*)histos.back())->Fill(x,y,w);
   else
      cout<<"ERROR in Plot named "<<templateHisto->GetName()
	  <<" Fill(..) called without calling prepareToFillProcess first"<<endl; 
}

// Fill the last histo here
void Plot::Fill(vector<Double_t> coord, double v, double w){
   if (histos.size()>0) 
      ((TProfileMDF*)histos.back())->Fill(coord,v,w);
   else
      cout<<"ERROR in Plot named "<<templateHisto->GetName()
	  <<" Fill(..) called without calling prepareToFillProcess first"<<endl; 
}

// Do the scaling to luminosity or data.
void Plot::scaleToData(vector<PhysicsProcess*> procs, DEFS::LeptonCat lepCat)
{
   // Make sure we can't do the scaling twice
   if (scaled)
      return;
  
   // Make sure that the number of histos and processes is the same
   if (procs.size() != histos.size())
   {
      cout<<"ERROR Plot::scaleToData() wrong sizes. Something is wrong"<<endl;
      return;
   }

   // First normalize to luminosity, and in case we might need to normalize to data
   // obtain also the total MC and total data
   double totalData = 0;
   double totalMC = 0;
   for (unsigned int p = 0 ; p < procs.size() ; p++)
   {
      //cout << "Normalizing to luminosity and cross section..." << endl;
      
      // This works for MC and data as well.
      histos[p]->Scale(procs[p]->getScaleFactor(lepCat));
      //cout << " histo="<<histos[p]->GetTitle()<<"\tscaleFactor="<<procs[p]->getScaleFactor(lepCat) << "\tleptonCat=" << DEFS::getLeptonCatString(lepCat) << endl;
      //cout<<" histo="<<histos[p]->GetTitle()<<", has Integral()="<<histos[p]->Integral()<<endl;

      TString hName = histos[p]->GetTitle();
      hName.ToUpper();
      if ( hName.Contains("DATA") )
         totalData += histos[p]->Integral();
      else 
         totalMC += histos[p]->Integral();
   }
   
   // If norm to data do the scaling of the all the MC to the data
   if (normToData)
   {
      //cout<<"Normalizing to data..."<<endl;
      
      if (totalData == 0)
      {
         cout<<"ERROR  Plot::scaleToData() integral of data is zero, cannot normalize to data."
             <<"WILL NOT NORMALIZE"<<endl;
         return;
      }
      
      //cout<<" totalData = "<<totalData<<"  totalMC = "<<totalMC<<endl;
      
      for (unsigned int p = 0 ; p < procs.size() ; p++)
      {
         // This works for MC and data as well.
         //histos[p]->Scale(procs[p]->getScaleFactor());
      
         TString hName = histos[p]->GetTitle();
         hName.ToUpper();
         if (!hName.Contains("DATA") )
            histos[p]->Scale(totalData/totalMC);
         
         //cout<<"AfterNorm histo="<<histos[p]->GetTitle()<<", has Integral()="<<histos[p]->Integral()<<endl; 
      }
    
   }
   
   //set this flag to true so we won't do the scaling again.
   scaled = true;
}

/* 
// ------------------------------------------------------------
// This groups all the histograms that need grouping.
// Typically all the Single tops and Dibosons
vector<TH1*> Plot::doGrouping(vector<PhysicsProcess*> procs)
{
   // The returning vector
   vector<TH1*> groups;
   if (procs.size()  != histos.size() )
   {
      cout<<"ERROR Plot::doGrouping procs.size()="<<procs.size()
          <<" is different than histos.size()="<<histos.size()<<endl;
      cout<<"\t Returning NO groups!"<<endl;
      return groups;
   }

   // The grouped ones
   TH1 * stop = (TH1*) templateHisto->Clone(TString(templateHisto->GetName())+"_STop");
   stop->SetTitle("STop");
   stop->Sumw2();
   TH1 * dibo = (TH1*) templateHisto->Clone(TString(templateHisto->GetName())+"_Diboson");
   dibo->SetTitle("WW+WZ+ZZ");
   dibo->Sumw2();
   TH1 * higgs = (TH1*) templateHisto->Clone(TString(templateHisto->GetName())+"_Higgs");
   higgs->SetTitle("ggH+WH+VBF");
   higgs->Sumw2();
   
   // Loop over histos grouping around
   for (unsigned int h=0; h < histos.size(); h ++)
   {
      TString hName = histos[h]->GetTitle();
      
      if (hName.Contains("STopT_Tbar")  ||
          hName.Contains("STopT_T")     ||
          hName.Contains("STopTW_Tbar") ||
          hName.Contains("STopTW_T")    ||
          hName.Contains("STopS_Tbar")  ||
          hName.Contains("STopS_T")     )
      {
         // if first time set the attributes
         if (stop->GetEntries()==0)
         {
            stop->SetLineColor(histos[h]->GetLineColor());
            stop->SetFillColor(histos[h]->GetFillColor());
            stop->SetMarkerColor(histos[h]->GetMarkerColor());
         }
         
         // add to single top
         stop->Add(histos[h]);
         
      }
      else if (hName.Contains("WW") ||
               hName.Contains("WZ") ||
               hName.Contains("ZZ") )
      {

         // if first time set the attributes
         if (dibo->GetEntries()==0)
         {
            dibo->SetLineColor(histos[h]->GetLineColor());
            dibo->SetFillColor(histos[h]->GetFillColor());
            dibo->SetMarkerColor(histos[h]->GetMarkerColor());
         }

         // add to diboson
         dibo->Add(histos[h]);

      }
      else if (hName.Contains("ggH") || 
               hName.Contains("qqH") ||
               hName.Contains("WH")  )
      {
         // if first time set the attributes
         if (higgs->GetEntries()==0)
         {
            higgs->SetLineColor(histos[h]->GetLineColor());
            higgs->SetFillColor(histos[h]->GetFillColor());
            higgs->SetMarkerColor(histos[h]->GetMarkerColor());
         }

         // add to higgs
         higgs->Add(histos[h]);
      }
      else
      {
         // add as an independent process.
         groups.push_back(histos[h]);
      }
   }

   // Add the higgs, diboson, and single top if they were present
   if (stop->GetEntries() >0)  groups.insert(groups.begin(), stop);
   if (dibo->GetEntries() >0)  groups.insert(groups.begin(), dibo);
   if (higgs->GetEntries() >0)  groups.insert(groups.begin(), higgs);

   // return the groupedHistos
   return groups;

}
*/
 // This is wrong, someone needs to fix it;
 // ** the output->close before the ProfileMDF::WriteToFile call will 
 // not allow the regular histo to be save after a TProfileMDF has been saved
void Plot::saveHistogramsToFile(TString histoFile, TString option)
{
   TFile output(histoFile, option);
   for (unsigned int i = 0; i < histos.size(); i++){
     if(TString(histos[i]->ClassName()).Contains("TProfileMDF")==1) {
       output.Close();
       ((TProfileMDF*)histos[i])->WriteToFile(histoFile,"UPDATE");
     }
     else if(histos[i]->GetEntries()>0){
       histos[i]->Write();
     }
     else {
       cout << "WARNING::Plot Not saving histogram " << histos[i]->GetName() << " because it doesn't have any entries." << endl;
       continue;
     }
   }
   if(output.IsOpen())
     output.Close();
}
void Plot::loadHistogramsFromFile(TDirectoryFile* dir, std::vector<PhysicsProcess*> procs, DEFS::LeptonCat lepCat, bool procNameOnly, TDirectoryFile* data_dir)
{
  int loaded = 0;

  cout << "Plot::loadHistogramsFromFile Loading histograms ... ";
  for(unsigned int i = 0; i < procs.size(); i++) {
    TString name;
    if(procNameOnly)
      name = Form("%s",procs[i]->name.Data());
    else
      name = Form("%s_%s_%s",templateHisto->GetTitle(),procs[i]->name.Data(),
                        DEFS::getLeptonCatString(lepCat).c_str());

    if(procs[i]->name.Contains("data",TString::kIgnoreCase) && data_dir!=0) {
      TString nameAlt = Form("%s_%s_%s",templateHisto->GetTitle(),procs[i]->name.Data(),
                        DEFS::getLeptonCatString(lepCat).c_str());
      TH1* data_hist(0); 

      if(data_dir->Get(name)) {
        data_hist = (TH1*)data_dir->Get(name);
      }
      else if(data_dir->Get(nameAlt)) {
        data_hist = (TH1*)data_dir->Get(nameAlt);
      }
      else {
        cout << "\tWARNING::Plot::loadHistogramsFromFile did not find the data histogram in the file provided" << endl;
      }

      histos.push_back((TH1D*)templateHisto->Clone(data_hist->GetName()));
      histos.back()->Reset();
      histos.back()->SetTitle(procs[i]->groupName);
      for(int ibin=1; ibin<=templateHisto->GetXaxis()->GetNbins(); ibin++) {
        histos.back()->SetBinContent(ibin,data_hist->GetBinContent(ibin));
        histos.back()->SetBinError(ibin,data_hist->GetBinError(ibin));
      }
      loaded++;
    }
    else if(dir->Get(name)) {
      histos.push_back((TH1*)dir->Get(name));
      histos.back()->SetTitle(procs[i]->groupName);
      loaded++;
    }
    else {
      TH1* tmp = (TH1*)histos.back()->Clone("empty");
      tmp->Reset();
      histos.push_back(tmp);
    }

    if(i==0 && (histos.back()->GetXaxis()->GetXmin()!=templateHisto->GetXaxis()->GetXmin() ||
       histos.back()->GetXaxis()->GetXmax()!=templateHisto->GetXaxis()->GetXmax())) {
      TString title = templateHisto->GetTitle();
      templateHisto = (TH1D*)histos.back()->Clone(templateHisto->GetName());
      templateHisto->SetTitle(title);
      templateHisto->Reset();
    }
  }
  cout << "DONE" << endl;
  cout << "\tPlot::loadHistogramsFromFile " << loaded << " histograms loaded for " << histos.size() << " processes" << endl; 
}

void Plot::loadHistogramsFromCombinedFile(TDirectoryFile* dir, DEFS::LeptonCat lepCat)
{
  int loaded = 0;
  bool found_data = false;

  cout << "Plot::loadHistogramsFromCombinedFile Loading histograms ... ";
  TIter nextHist(dir->GetListOfKeys());
  TKey* histKey(0);
  while ((histKey=(TKey*)nextHist())) {
    if (strcmp(histKey->GetClassName(),"TH1D")!=0) continue;
    TH1* tmp = (TH1*)histKey->ReadObj();
    if(tmp!=nullptr && TString(histKey->GetName()).BeginsWith(templateHisto->GetTitle())) {
      string name(tmp->GetName());
      if(TString(tmp->GetName()).Contains("data",TString::kIgnoreCase)) {
        found_data=true;
        histos.push_back((TH1D*)templateHisto->Clone(tmp->GetName()));
        histos.back()->Reset();
        int pos = name.find("Single");
        string processName = name.substr(pos,name.find("_"+DEFS::getLeptonCatString(lepCat))-pos);
        histos.back()->SetTitle(DEFS::PhysicsProcess::getTypeTitle(DEFS::PhysicsProcess::getProcessType(processName)).c_str());
        for(int ibin=1; ibin<=templateHisto->GetXaxis()->GetNbins(); ibin++) {
          histos.back()->SetBinContent(ibin,tmp->GetBinContent(ibin));
          histos.back()->SetBinError(ibin,tmp->GetBinError(ibin));
        }
      }
      else {
        histos.push_back(tmp);
      }
      loaded++;
    }
    if(histos.size()==1 && (histos.back()->GetXaxis()->GetXmin()!=templateHisto->GetXaxis()->GetXmin() ||
       histos.back()->GetXaxis()->GetXmax()!=templateHisto->GetXaxis()->GetXmax())) {
      TString title = templateHisto->GetTitle();
      templateHisto = (TH1D*)histos.back()->Clone(templateHisto->GetName());
      templateHisto->SetTitle(title);
      templateHisto->Reset();
    }
  }
  cout << "DONE" << endl;
  if(!found_data)
    cout << "\tWARNING::Plot::loadHistogramsFromCombinedFile did not find the data histogram in the file provided" << endl;
  cout << "\tPlot::loadHistogramsFromCombinedFile " << loaded << " histograms loaded" << endl; 
}

void Plot::printList() {
  cout << "LIST OF HISTOGRAMS FOR PLOT " << templateHisto->GetName() << ":" << endl;
  for(unsigned int i=0; i<histos.size(); i++) {
    cout << "\tname=" << histos[i]->GetName() << "\ttitle=" << histos[i]->GetTitle() << endl;
  }
}

// ##################################################
// ############## FORMATTED PLOT CLASS ##############
// ##################################################

FormattedPlot::FormattedPlot()
{
   overlaySignalFactor = 1.0;
   overlaySignalName = "";
   logxy = make_pair(false,false);
}

// ------------------------------------------------------------
// Make the Canvas here ala Ricardo Eusebi
// ------------------------------------------------------------
TCanvas* FormattedPlot::getCanvas(vector<PhysicsProcess*> procs)
{
   // These flags are to determine if this FormattedPlot contains data or Monte Carlo histograms
   // With these flags, it is possible for the user to create canvases with only data or only Monte Carlo
   bool areMCHists = false;
   bool areDataHists = false;
   
   // Format the colors
   formatColors(procs);
   
   // Do the scaling of the histos to lum or to data
   scaleToData(procs,leptonCat);
  
   // Do the grouping of the histos and return the list of histograms to be plotted.
   // This takes care of putting all the processes with the same groupName together
   // typically SingleTop Histos together, or diboson together etc.
   std::vector<TH1*> groupedHistos = doGrouping(procs);
   //for(unsigned int i = 0; i < procs.size(); i++) {cout << "procs[i]->groupName = " << procs[i]->groupName << endl;}
   // Create the total Data and MC histos and Stacks
   TString tempName = templateHisto->GetName();
   TH1     * tMC = (TH1*) templateHisto->Clone(tempName+"_TotalMC");   tMC->Sumw2();
   TH1     * tDa = (TH1*) templateHisto->Clone(tempName+"_TotalData"); tDa->Sumw2();
   TH1     * sys_JESUp = (TH1*) templateHisto->Clone(tempName+"_Systematic_JESUp"); sys_JESUp->Sumw2();
   TH1     * sys_JESDown = (TH1*) templateHisto->Clone(tempName+"_Systematic_JESDown"); sys_JESDown->Sumw2();
   TH1     * sys_matchingup = (TH1*) templateHisto->Clone(tempName+"_Systematic_matchingup"); sys_matchingup->Sumw2();
   TH1     * sys_matchingdown = (TH1*) templateHisto->Clone(tempName+"_Systematic_matchingdown"); sys_matchingdown->Sumw2();
   TH1     * sys_scaleup = (TH1*) templateHisto->Clone(tempName+"_Systematic_scaleup"); sys_scaleup->Sumw2();
   TH1     * sys_scaledown = (TH1*) templateHisto->Clone(tempName+"_Systematic_scaledown"); sys_scaledown->Sumw2();
   TH1     * signal = 0;
   THStack * sMC = new THStack(tempName+"_stackMC",tempName+"_stackMC");
   THStack * sDa = new THStack(tempName+"_stackData",tempName+"_stackData");

   // Make the legend
   TLegend * l = new TLegend(0.8,0.4,0.96,0.89);
   l->SetName(tempName+"_legend");
   l->SetBorderSize(0);
   l->SetFillColor(0);
   l->SetFillStyle(0); 
  
   // Create the totalHistos and Stacks for the MC and Data processes
   for (unsigned int h = 0 ; h < groupedHistos.size() ; h++)
   {
      TString hname = groupedHistos[h]->GetTitle();
      hname.ToUpper();
      //cout << "groupedHistos[h]->GetTitle() = " << hname << endl; 
      if (hname.Contains("DATA"))
      {
         sDa->Add(groupedHistos[h],"hist");
         tDa->Add(groupedHistos[h]);
         l->AddEntry(groupedHistos[h],groupedHistos[h]->GetTitle(),"lpe");
         
         areDataHists = true;
      }
      else if (hname.Contains("SYSTEMATIC")) {
        if(hname.Contains("JES")) {
          if(hname.Contains("UP"))
            sys_JESUp->Add(groupedHistos[h]);
          else if(hname.Contains("DOWN"))
            sys_JESDown->Add(groupedHistos[h]);
        }
      }
      else
      {
         sMC->Add(groupedHistos[h],"hist");
         tMC->Add(groupedHistos[h]);

         if(hname.CompareTo("WJets")!=0 && !hname.Contains("JES")) {
           if(hname.Contains("MATCHING")) {
             if(hname.Contains("UP"))
               sys_matchingup->Add(groupedHistos[h]);
             else if(hname.Contains("DOWN"))
               sys_matchingdown->Add(groupedHistos[h]);
           }
           else if(hname.Contains("SCALE")) {
             if(hname.Contains("UP"))
               sys_scaleup->Add(groupedHistos[h]);
             else if(hname.Contains("DOWN"))
               sys_scaledown->Add(groupedHistos[h]);
           }
           else {
             sys_matchingup->Add(groupedHistos[h]);
             sys_matchingdown->Add(groupedHistos[h]);
             sys_scaleup->Add(groupedHistos[h]);
             sys_scaledown->Add(groupedHistos[h]);
           }
         }

         l->AddEntry(groupedHistos[h],groupedHistos[h]->GetTitle(),"f");
         
         areMCHists = true;
      }
      
      // Check if this is the signal we want to overlay
      if (overlaySignalName.Contains(groupedHistos[h]->GetTitle())){
	
	//TString className = tMC->ClassName(); 
	//if (!className.Contains("TH2")){
	  
	  ostringstream signalNameStream;
	  signalNameStream <<  groupedHistos[h]->GetName() << "_x" << overlaySignalFactor;
	  signal = (TH1*)groupedHistos[h]->Clone(TString(signalNameStream.str()));
	  signal->SetFillColor(0);
    signal->SetLineColor(kRed-4);
	  signal->Scale(overlaySignalFactor);
          
	  ostringstream legendTitleStream;
	  legendTitleStream << signal->GetTitle() << "(x" << overlaySignalFactor << ")";
	  string legendTitle = legendTitleStream.str();
	  l->AddEntry(signal,legendTitle.c_str(),"lp0");

	  //}// skip TH2 for some reasons.

      } // if this is the signal histo
      
   }// for grouped histos
   
   // General Style Formatting
   gStyle->SetOptTitle(0);
   //gStyle->SetOptStat(2211);

   // Create and Draw the Canvas
   TCanvas * can = new TCanvas(templateHisto->GetTitle());
   can->SetFillColor(0);
   can->Divide(1,2);
   
   // Define the two different pads that go on the canvas
   // canMain: this is the pad with the actual histogram on it
   TVirtualPad *canMain = can->GetPad(1);
   canMain->SetPad(0.01,0.32,0.99,0.99);
   canMain->SetObjectStat(true);
   canMain->SetBottomMargin(0.16); 
   canMain->SetLeftMargin(0.115); 
   canMain->SetRightMargin(0.03); 
   canMain->cd();
   TPad* padMain = (TPad*)canMain->GetPad(0);
   padMain->SetLogx(logxy.first);
   padMain->SetLogy(logxy.second);
   
   //Define the graphics option
   TString gOption = stacked ? "": "nostack";
   if (areDataHists)
      sMC->Draw(gOption);
   else
      sMC->Draw(gOption+"colz");
   
   if (areMCHists)
      sDa->Draw(gOption+"ep SAME");
   else
      sDa->Draw(gOption+"ep");

   //Draw the MC error bars
   if(tMC != 0)
   {
      //Loop through each of the beins in the tMC histogram and reset the errors to be stats+systematics
      TGraphAsymmErrors* errorBand = makeSystematicErrors(tMC);
      addSystematicErrors(errorBand, tMC,sys_JESUp,sys_JESDown);
    //addSystematicErrors(errorBand, tMC,sys_matchingup,sys_matchingdown);
    //addSystematicErrors(errorBand, tMC,sys_scaleup,sys_scaledown);
      errorBand->SetMarkerColor(kGray+2);
      errorBand->SetLineColor(kGray+2);
      errorBand->SetFillColor(kGray+2);
      errorBand->SetFillStyle(3154);
      errorBand->Draw("2 SAME");
      //tMC->SetFillColor(kGray+2);
      //tMC->SetFillStyle(3005);
      //tMC->Draw("E2 SAME");
   }
   
   //Draw the signal overlay
   if (signal != 0)
   {
      signal->SetLineWidth(2);
      signal->Draw("hist,SAME");
   }
   
   // Set display range
   if (areMCHists)
      sMC->GetXaxis()->SetRangeUser(range.first, range.second);
   if (areDataHists)
      sDa->GetXaxis()->SetRangeUser(range.first, range.second);
   if (signal != 0)
      signal->GetXaxis()->SetRangeUser(range.first, range.second);

   // Format the stacks, make sure to do this after the Draw command.
   // Set the maximum range for the Y-axis
   double maxi = max(sMC->GetMaximum(),sDa->GetMaximum())*1.1;
   if (areMCHists)
      formatStack(sMC,maxi);
   //formatStack(sDa,maxi);
   
   // Add the KS and Chi2 info in the active canvas
   if (areMCHists && areDataHists)
      drawKSandChi2Tests(tDa, tMC, range);

   // Add the Legend
   l->Draw("same");

   // Add the luminosity
   //drawLumi(3.6);
   drawLumi(procs[0]->intLum[leptonCat]/1000);

   // canRatio: this is the pad with the Data/MC on it
   TVirtualPad * canRatio = can->GetPad(2);
   canRatio->SetPad(0.01,0.01,0.99,0.32);
   canRatio->SetBottomMargin(0.16);
   canRatio->SetLeftMargin(0.115); 
   canRatio->SetRightMargin(0.03);  
   canRatio->SetGridx();
   canRatio->SetGridy();
   canRatio->cd();
   TPad* padRatio = (TPad*)canRatio->GetPad(0);
   padRatio->SetLogx(logxy.first);
   //padRatio->SetLogy(logxy.second);
  
   // Create the Ratio plot
   TH1* hRatio = (TH1*) tDa->Clone(tempName+"_Ratio");
   hRatio->SetTitle("#frac{Data - MC}{MC}");
   hRatio->Add(tMC,-1);
   hRatio->Divide(tMC);

   // Format the ratio plot and Draw it
   formatRatio(hRatio);
   hRatio->Draw();

   // Return the Canvas
   return can;
}

// ------------------------------------------------------------
// Just like the getCanvas member function, but in the TDR style
TCanvas* FormattedPlot::getCanvasTDR(vector<PhysicsProcess*> procs) {
  Style* st = new Style();
  st->setTDRStyle();
  TGaxis::SetMaxDigits(3);

  // These flags are to determine if this FormattedPlot contains data or Monte Carlo histograms
  // With these flags, it is possible for the user to create canvases with only data or only Monte Carlo
  bool areMCHists = false;
  bool areDataHists = false;

  // Format the colors
  formatColors(procs);

  // Do the scaling of the histos to lum or to data
  scaleToData(procs,leptonCat);

  // Do the grouping of the histos and return the list of histograms to be plotted.
  // This takes care of putting all the processes with the same groupName together
  // typically SingleTop Histos together, or diboson together etc.
  std::vector<TH1*> groupedHistos = doGrouping(procs);
  // Create the total Data and MC histos and Stacks
  TString tempName = templateHisto->GetName();
  TH1     * tMC = (TH1*) templateHisto->Clone(tempName+"_TotalMC");   tMC->Sumw2();
  TH1     * tWJets = (TH1*) templateHisto->Clone(tempName+"_TotalWJets");   tWJets->Sumw2();
  TH1     * tTTbar = (TH1*) templateHisto->Clone(tempName+"_TotalTTbar");   tTTbar->Sumw2();
  TH1     * tQCDandWJets = (TH1*) templateHisto->Clone(tempName+"_TotalQCDandWJets");   tQCDandWJets->Sumw2();
  TH1     * tDa = (TH1*) templateHisto->Clone(tempName+"_TotalData"); tDa->Sumw2();
  TH1     * sys_JESUp = (TH1*) templateHisto->Clone(tempName+"_Systematic_JESUp"); sys_JESUp->Sumw2();
  TH1     * sys_JESDown = (TH1*) templateHisto->Clone(tempName+"_Systematic_JESDown"); sys_JESDown->Sumw2();
  TH1     * sys_matchingup = (TH1*) templateHisto->Clone(tempName+"_Systematic_matchingup"); sys_matchingup->Sumw2();
  TH1     * sys_matchingdown = (TH1*) templateHisto->Clone(tempName+"_Systematic_matchingdown"); sys_matchingdown->Sumw2();
  TH1     * sys_scaleup = (TH1*) templateHisto->Clone(tempName+"_Systematic_scaleup"); sys_scaleup->Sumw2();
  TH1     * sys_scaledown = (TH1*) templateHisto->Clone(tempName+"_Systematic_scaledown"); sys_scaledown->Sumw2();
  TH1     * sys_csvup = (TH1*) templateHisto->Clone(tempName+"_Systematic_CSVWeightUp"); sys_csvup->Sumw2();
  TH1     * sys_csvdown = (TH1*) templateHisto->Clone(tempName+"_Systematic_CSVWeightDown"); sys_csvdown->Sumw2();
  TH1     * sys_puup = (TH1*) templateHisto->Clone(tempName+"_Systematic_PUWeightUp"); sys_puup->Sumw2();
  TH1     * sys_pudown = (TH1*) templateHisto->Clone(tempName+"_Systematic_PUWeightDown"); sys_pudown->Sumw2();
  TH1     * sys_topptup = (TH1*) templateHisto->Clone(tempName+"_Systematic_topPtWeightUp"); sys_topptup->Sumw2();
  TH1     * sys_topptdown = (TH1*) templateHisto->Clone(tempName+"_Systematic_topPtWeightDown"); sys_topptdown->Sumw2();
  TH1     * sys_qcdetaup = (TH1*) templateHisto->Clone(tempName+"_Systematic_QCDEtaWeightUp"); sys_qcdetaup->Sumw2();
  TH1     * sys_qcdetadown = (TH1*) templateHisto->Clone(tempName+"_Systematic_QCDEtaWeightDown"); sys_qcdetadown->Sumw2();
  TH1     * sys_statup = (TH1*) templateHisto->Clone(tempName+"_Systematic_StatUp"); sys_statup->Sumw2();
  TH1     * sys_statdown = (TH1*) templateHisto->Clone(tempName+"_Systematic_StatDown"); sys_statdown->Sumw2();
  TH1     * signal = 0;
  THStack * sMC = new THStack(tempName+"_stackMC",tempName+"_stackMC");
  THStack * sDa = new THStack(tempName+"_stackData",tempName+"_stackData");

  // Make the legend
  TLegend * l = st->tdrLeg(0.44,0.64,0.62,0.89);
  l->SetName(tempName+"_legend");
  TLegend * lsig = st->tdrLeg(0.66,0.64,0.86,0.89);
  lsig->SetName(tempName+"_signal_legend");

  // Create the totalHistos and Stacks for the MC and Data processes
  for (unsigned int h = 0 ; h < groupedHistos.size() ; h++) {
    TString hname = groupedHistos[h]->GetTitle();
    hname.ToUpper();

    if (hname.Contains("DATA")) {
      sDa->Add(groupedHistos[h],"hist");
      tDa->Add(groupedHistos[h]);
      lsig->AddEntry(groupedHistos[h],groupedHistos[h]->GetTitle(),"lpe");

      areDataHists = true;
    }
    else if (hname.Contains("SYSTEMATIC") || hname.Contains("CMS_")) {
      if(hname.Contains("JES") || hname.Contains("SCALE_J_SHAPE")) {
        if(hname.Contains("UP"))
          sys_JESUp->Add(groupedHistos[h]);
        else if(hname.Contains("DOWN"))
          sys_JESDown->Add(groupedHistos[h]);
      }
      else if(hname.Contains("CSVWEIGHT_SHAPE")) {
        if(hname.Contains("UP"))
          sys_csvup->Add(groupedHistos[h]);
        else if(hname.Contains("DOWN"))
          sys_csvdown->Add(groupedHistos[h]); 
      }
      else if(hname.Contains("PUWEIGHT_SHAPE")) {
        if(hname.Contains("UP"))
          sys_puup->Add(groupedHistos[h]);
        else if(hname.Contains("DOWN"))
          sys_pudown->Add(groupedHistos[h]); 
      }
      else if(hname.Contains("MATCHING")) {
        if(hname.Contains("UP"))
          sys_matchingup->Add(groupedHistos[h]);
        else if(hname.Contains("DOWN"))
          sys_matchingdown->Add(groupedHistos[h]);
      }
      else if(hname.Contains("SCALE")) {
        if(hname.Contains("UP"))
          sys_scaleup->Add(groupedHistos[h]);
        else if(hname.Contains("DOWN"))
          sys_scaledown->Add(groupedHistos[h]);
      }
      else if(hname.Contains("TOPPTWEIGHT_SHAPE")) {
        if(hname.Contains("UP"))
          sys_topptup->Add(groupedHistos[h]);
        else if(hname.Contains("DOWN"))
          sys_topptdown->Add(groupedHistos[h]);
      }
      else if(hname.Contains("QCDETAWEIGHT_SHAPE")) {
        if(hname.Contains("UP"))
          sys_qcdetaup->Add(groupedHistos[h]);
        else if(hname.Contains("DOWN"))
          sys_qcdetadown->Add(groupedHistos[h]);
      }
      else if (hname.Contains("STATERROR")) {
        if(hname.Contains("UP"))
          sys_statup->Add(groupedHistos[h]);
        else if(hname.Contains("DOWN"))
          sys_statdown->Add(groupedHistos[h]); 
      }

    }
    else {
      sMC->Add(groupedHistos[h],"hist");
      tMC->Add(groupedHistos[h]);

      if(hname.Contains("W+JETS")) {
        tWJets->Add(groupedHistos[h]);
        tQCDandWJets->Add(groupedHistos[h]);
      }
      else if(hname.CompareTo("QCD")==0) {
        tQCDandWJets->Add(groupedHistos[h]);
      }
      else if(hname.Contains("T#BAR{T}")) {
        tTTbar->Add(groupedHistos[h]);
      }

      if(hname.Contains("H(125)"))
        lsig->AddEntry(groupedHistos[h],groupedHistos[h]->GetTitle(),"f");
      else
        l->AddEntry(groupedHistos[h],groupedHistos[h]->GetTitle(),"f");

      areMCHists = true;
    }

    // Check if this is the signal we want to overlay
    if (overlaySignalName.Contains(groupedHistos[h]->GetTitle())){

      if(titleContainedInTH1Vector("DATA",groupedHistos)) {
        overlaySignalFactor =  titleContainedInTH1Vector("DATA",groupedHistos)->Integral()/groupedHistos[h]->Integral();
        //cout << "titleContainedInTH1Vector(\"DATA\",groupedHistos)->Integral() = " << titleContainedInTH1Vector("DATA",groupedHistos)->Integral() << endl;
        //cout << "groupedHistos[h]->Integral() = " << groupedHistos[h]->Integral() << endl;
        //cout << "original overlaySignalFactor = " << overlaySignalFactor << endl;
        //cout << "(double)RoundToNearest(int(overlaySignalFactor),100) = " << (double)RoundToNearest(int(overlaySignalFactor),100) << endl; 
        overlaySignalFactor = max(3100.0,(double)RoundToNearest(int(overlaySignalFactor),100));
      }

      ostringstream signalNameStream;
      signalNameStream <<  groupedHistos[h]->GetName() << "_x" << overlaySignalFactor;
      signal = (TH1*)groupedHistos[h]->Clone(TString(signalNameStream.str()));
      signal->Scale(overlaySignalFactor);

      //ostringstream legendTitleStream;
      //legendTitleStream << signal->GetTitle() << "(x" << overlaySignalFactor << ")";
      //string legendTitle = legendTitleStream.str();
      //lsig->AddEntry(signal,legendTitle.c_str(),"lp0");
      //Form("%s(x%.0d)",signal->GetTitle(),overlaySignalFactor);
      lsig->AddEntry(signal,signal->GetTitle(),"lp0");
      lsig->AddEntry((TObject*)0,Form("(x%.0f)",overlaySignalFactor),"");

    } // if this is the signal histo
      
  }// for grouped histos

  //Make the error bands
  TGraphAsymmErrors* errorBand;
  if(tMC != 0) {
    //Setup the normalization errors
    double lumi_8TeV = 0.026;
    double CMS_eff_l = 0.02;
    double CMS_scale_met = 0.002;
    //double CMS_QCDscale_pdf_ttbar = 0.057;
    //double CMS_QCDscale_pdf_zjets = 0.034;
    //double CMS_QCDscale_pdf_singlet = 0.05;
    //double CMS_QCDscale_pdf_diboson = 0.03;
    double totalNorm = sqrt(pow(lumi_8TeV,2)+pow(CMS_eff_l,2)+pow(CMS_scale_met,2));

    //Loop through each of the bins in the tMC histogram and reset the errors to be stats+systematics
    errorBand = makeSystematicErrors(tMC,totalNorm);
    addSystematicErrors(errorBand, tMC, sys_JESUp, sys_JESDown, "JES");
    addSystematicErrors(errorBand, tWJets, sys_matchingup, sys_matchingdown, "Matching");
    addSystematicErrors(errorBand, tWJets, sys_scaleup, sys_scaledown, "Scale");
    addSystematicErrors(errorBand, tMC, sys_csvup, sys_csvdown, "CSV");
    addSystematicErrors(errorBand, tMC, sys_puup, sys_pudown, "PU");
    addSystematicErrors(errorBand, tTTbar, sys_topptup, sys_topptdown, "TopPt");
    addSystematicErrors(errorBand, tQCDandWJets, sys_qcdetaup, sys_qcdetadown, "QCDEta");
    addSystematicErrors(errorBand, tMC, sys_statup, sys_statdown, "stat");
  }

  // Create and Draw the Canvas
  TH1D* frame_up = new TH1D();
  frame_up->GetXaxis()->SetLimits(range.first, range.second);
  frame_up->GetXaxis()->SetMoreLogLabels();
  frame_up->GetXaxis()->SetNoExponent();
  vector<double> maximums;
  if(tMC) maximums.push_back(tMC->GetBinContent(tMC->GetMaximumBin()));
  if(tDa) maximums.push_back(tDa->GetBinContent(tDa->GetMaximumBin()));
  if(signal) maximums.push_back(signal->GetBinContent(signal->GetMaximumBin()));
  if(errorBand) {
    double maxError = 0;
    for(int ipoint=0; ipoint<errorBand->GetN(); ipoint++) {
      if((errorBand->GetY()[ipoint]+errorBand->GetErrorYhigh(ipoint))>maxError)
        maxError = errorBand->GetY()[ipoint]+errorBand->GetErrorYhigh(ipoint);
    }
    maximums.push_back(maxError);
  }
  double multiplier = 1.45;//1.25;
  frame_up->GetYaxis()->SetRangeUser(0.0001,*std::max_element(maximums.begin(),maximums.end())*multiplier);
  frame_up->GetXaxis()->SetTitle(axisTitles[0].c_str());
  frame_up->GetYaxis()->SetTitle(axisTitles[1].c_str());
  TH1D* frame_down = new TH1D();
  frame_down->GetXaxis()->SetLimits(range.first, range.second);
  frame_down->GetXaxis()->SetMoreLogLabels();
  frame_down->GetXaxis()->SetNoExponent();
  frame_down->GetYaxis()->SetRangeUser(0.0,2.0); //(-0.5,0.5) <-- (-1.0,1.0)
  frame_down->GetXaxis()->SetTitle(axisTitles[0].c_str());
  frame_down->GetYaxis()->SetTitle("Data/MC");
  if (signal != 0) {
    signal->GetXaxis()->SetRangeUser(range.first, range.second);
  }
  TCanvas * can = st->tdrDiCanvas(templateHisto->GetTitle(),frame_up,frame_down,2,11,procs[0]->intLum[leptonCat]/1000);
  can->cd(1);

  TPad* padMain = (TPad*)can->GetPad(0);
  padMain->SetLogx(logxy.first);
  padMain->SetLogy(logxy.second);

   //Define the graphics option
  string gOption = stacked ? "": "nostack";
  if (areDataHists)
    st->tdrDraw(sMC,gOption);
  else
    st->tdrDraw(sMC,gOption+"colz");

  //Draw the signal overlay
  if (signal != 0) {
    signal->SetLineWidth(2);
    st->tdrDraw(signal,"hist",kNone,kRed-4,kSolid,kRed-4,kNone,0);
  }

  if (areMCHists)
    st->tdrDraw(sDa,gOption+"ep");
  else
    st->tdrDraw(sDa,gOption+"ep");

  //Draw the MC error bars
  if(errorBand!=nullptr) {
    cout << "Drawing systematic error band..." << flush;
    gStyle->SetHatchesSpacing(0.03);   
    st->tdrDraw(errorBand,"2",kNone,kGray+2,kNone,kGray+2,3004,kGray+2);
    //st->tdrDraw(tMC,"E2",kNone,kGray+2,kNone,kGray+2,3005,kGray+2);
    cout << "DONE" << endl;
  }
  else {
    cout << "Can't draw the systematic error band because the object \"errorBand\" is null." << endl;
  }

  // Add the Legend
  l->Draw("same");
  lsig->Draw("same");

  // Draw some more label information
  drawBinInfo();

  // canRatio: this is the pad with the Data/MC on it
  TVirtualPad * canRatio = can->GetPad(2);
  canRatio->SetGridx();
  canRatio->SetGridy();
  canRatio->cd();
  TPad* padRatio = (TPad*)canRatio->GetPad(0);
  padRatio->SetLogx(logxy.first);
  st->tdrGrid(true);
  
  //Draw a line at 0
  TLine* line = new TLine(range.first,1,range.second,1);
  line->SetLineColor(kBlack);
  line->SetLineWidth(2);
  line->SetLineStyle(2);
  line->Draw("same");

  // Add the KS and Chi2 info in the active canvas
  //if (areMCHists && areDataHists) {
  //  drawKSandChi2Tests(tDa, tMC, range, true);
  //}

  // Create the Ratio plot
  TH1* hRatio = (TH1*) tDa->Clone(tempName+"_Ratio");
  hRatio->SetTitle("Data/MC");
  //hRatio->Add(tMC,-1);
  hRatio->Divide(tMC);
  // Format the ratio plot and Draw it
  st->tdrDraw(hRatio,"");

  TH1* BkgOverBkg = (TH1*) hRatio->Clone("bkgOverbkg");
  BkgOverBkg->Divide(tMC, tMC);
  TGraphAsymmErrors* pullUncBandTot = new TGraphAsymmErrors((TH1*)BkgOverBkg->Clone("pulluncTot"));
  for(int ibin=1; ibin<=tDa->GetNbinsX(); ibin++) {
    if (tMC->GetBinContent(ibin)!=0) {
      pullUncBandTot->SetPointEYhigh(ibin-1,errorBand->GetErrorYhigh(ibin-1)/tMC->GetBinContent(ibin));
      pullUncBandTot->SetPointEYlow(ibin-1,errorBand->GetErrorYlow(ibin-1)/tMC->GetBinContent(ibin));   
    }
  }
  st->tdrDraw(pullUncBandTot,"2",kNone,kGray+2,kNone,kGray+2,3004,kGray+2);

  delete st;

  // Return the Canvas
  return can;
}

// ------------------------------------------------------------
// Find the first histo in groupedHistos that matches the title,
// return a null pointer if it can't find one
TH1 * FormattedPlot::findTitleInTH1Vector(TString title, vector<TH1*> groupedHistos) {
  
   // Loop over original histos
  for (unsigned int h=0; h < groupedHistos.size(); h ++)
    if (title.EqualTo(groupedHistos[h]->GetTitle())) 
      return groupedHistos[h];

    return 0;

}//findTitleInTH1Vector

// ------------------------------------------------------------
// Find the first histo in groupedHistos that contains the title,
// return a null pointer if it can't find one
TH1 * FormattedPlot::titleContainedInTH1Vector(TString title, vector<TH1*> groupedHistos) {
  
   // Loop over original histos
  for (unsigned int h=0; h < groupedHistos.size(); h ++)
    if (TString(groupedHistos[h]->GetTitle()).Contains(title,TString::kIgnoreCase))
      return groupedHistos[h];

    return 0;

}//findTitleInTH1Vector
 

// ------------------------------------------------------------
// This groups all the histograms with the same title together
// Since the titles are given by process->groupName it groups all processes
// with the same groupNames together
vector<TH1*> FormattedPlot::doGrouping(vector<PhysicsProcess*> procs){

  // The returning vector
  vector<TH1*> groupedHistos;
  
  // Loop over original histos
  for (unsigned int h=0; h < histos.size(); h ++){
    //cout << "\thistos[h]->GetName() = " << histos[h]->GetName() << endl;
    //cout << "\thistos[h]->GetTitle() = " << histos[h]->GetTitle() << endl;
    // Find the histo in groupedHistos that matches the title, return null pointer if it can't.
    TH1 * hgroup = findTitleInTH1Vector(histos[h]->GetTitle(),groupedHistos);
    
    // if it does not exist, add one to the groupedHistos
    if (hgroup==0){
      TString cloneName = histos[h]->GetName();
      cloneName += "_clone";
      groupedHistos.push_back((TH1*) histos[h]->Clone(cloneName));
    }
    // if it does exist just add this to the already existing one
    else 
      hgroup->Add(histos[h]);
  }

  return groupedHistos;

}//doGrouping

// ------------------------------------------------------------
void FormattedPlot::formatColors(vector<PhysicsProcess*> procs)
{
   for(unsigned int i = 0; i < histos.size(); i++) {
      int proc_index = 0;
      for(unsigned int p = 0; p < procs.size(); p++) {
        if(TString(histos[i]->GetName()).Contains(procs[p]->name,TString::kIgnoreCase))
          proc_index = p;
      }

      histos[i]->SetLineColor(kBlack);
      histos[i]->SetMarkerColor(((PlotterPhysicsProcess*)procs[proc_index])->color);
      histos[i]->SetFillColor(((PlotterPhysicsProcess*)procs[proc_index])->color);
      if(procs[proc_index]->name.Contains("data",TString::kIgnoreCase)) {
         histos[i]->SetMarkerStyle(20);
         histos[i]->SetMarkerSize(0.7);
      }
      else {
         histos[i]->SetFillStyle(1001); 
      }
   }
   /*
   for(unsigned int i = 0; i < procs.size(); i++)
   {
      //histos[i]->SetLineColor(((PlotterPhysicsProcess*)procs[i])->color);
      histos[i]->SetLineColor(kBlack);
      histos[i]->SetMarkerColor(((PlotterPhysicsProcess*)procs[i])->color);
      histos[i]->SetFillColor(((PlotterPhysicsProcess*)procs[i])->color);
      
      TString pname = procs[i]->name;
      pname.ToUpper();

      if(!pname.Contains("DATA"))
      {
         histos[i]->SetFillStyle(1001);   
      }
      if (pname.Contains("DATA"))
      { 
         histos[i]->SetMarkerStyle(20);
         histos[i]->SetMarkerSize(0.7);
      }
   }
   */
}

// ------------------------------------------------------------
void FormattedPlot::formatRatio(TH1* hRatio)
{
   hRatio->SetMinimum(0.0); //formerly -0.99
   hRatio->SetMaximum(2.0); //formerly -0.99
   hRatio->SetLineWidth(2);

   hRatio->GetXaxis()->SetLabelFont(42);
   hRatio->GetXaxis()->SetLabelOffset(0.007);
   hRatio->GetXaxis()->SetLabelSize(0.11);
   hRatio->GetXaxis()->SetTitleSize(0.12);
   hRatio->GetXaxis()->SetTitleFont(42);
   hRatio->GetXaxis()->SetRangeUser(range.first,range.second);

   hRatio->GetYaxis()->SetTitle("#frac{Data}{MC}");
   hRatio->GetYaxis()->SetTitleOffset(0.55);
   hRatio->GetYaxis()->SetNdivisions(105);
   hRatio->GetYaxis()->SetLabelFont(42);
   hRatio->GetYaxis()->SetLabelOffset(0.012);
   hRatio->GetYaxis()->SetLabelSize(0.11);
   hRatio->GetYaxis()->SetTitleSize(0.1);
   hRatio->GetYaxis()->SetTitleFont(42);

   hRatio->SetStats(0);
}

// ------------------------------------------------------------
void FormattedPlot::formatStack(THStack * stack, double maxi)
{
   for(unsigned int a=0; a<axisTitles.size(); a++)
   {
      if (a==0) // X axis
      {
         stack->GetXaxis()->SetTitleOffset(1.1);
         stack->GetXaxis()->SetTitle(axisTitles[a].c_str());
         stack->GetXaxis()->SetRangeUser(range.first,range.second);
         stack->GetXaxis()->SetLabelSize(0.06);
         stack->GetXaxis()->SetTitleSize(0.06);
         stack->GetXaxis()->SetLabelFont(42);
         stack->GetXaxis()->SetTitleFont(42);
      }
      else if (a==1) // Y axis
      {
         stack->GetYaxis()->SetTitleOffset(1.3);
         stack->GetYaxis()->SetTitle(axisTitles[a].c_str());
         stack->SetMaximum(maxi);
         stack->GetYaxis()->SetLabelSize(0.06);
         stack->GetYaxis()->SetTitleSize(0.06);
         stack->GetYaxis()->SetLabelFont(42);
         stack->GetYaxis()->SetTitleFont(42);
         stack->GetYaxis()->SetTitleOffset(0.85);
      }
   }
}

//________________________________________________________________________________
void FormattedPlot::drawKSandChi2Tests(TH1* totalData, TH1* all, pair<double, double> range, bool doTDR, Style* st){

    // Skip all this if either histo has no integral
    if (totalData->Integral() == 0){
      cout<<"\tWARNING in Plots::drawKSandChi2Tests() not drawing KS or chi2 because  "
	  <<totalData->GetName()<<" has zero integral"<<endl;
      return ;
    }
    // Skip all this if either histo has no integral
    if (all->Integral() == 0){
      cout<<"\tWARNING in Plots::drawKSandChi2Tests() not drawing KS or chi2 because  "
	  <<all->GetName()<<" has zero integral"<<endl;
      return ;
    }
  
   double x = (range.second - range.first)*0.45 + range.first;
   double y = max(totalData->GetMaximum(),all->GetMaximum())*1.1;

   double chi2 = 0;
   int NDF = 0;
   int igood;
   double chi2NDF = totalData->Chi2TestX(all,chi2,NDF,igood,"UW CHI2/NDF");

/*   
   TLatex * ks      = new TLatex(x, y     , Form("KSTest   = %5.4g", totalData->KolmogorovTest(all)));
   ks->SetName("ks");
   TLatex * chi2P   = new TLatex(x, y*0.92, Form("Chi2Prob = %5.4g", all->Chi2TestX(totalData,chi2,NDF,igood,"WW")));
   chi2P->SetName("chi2P");
   TLatex * chi2NDF = new TLatex(x, y*0.84, Form("Chi2/NDF = %5.4g", chi2/NDF));
   chi2NDF->SetName("chi2NDF");
   
   ks->Draw();
   chi2P->Draw();
   chi2NDF->Draw();
*/
   TLatex Tl;
   if(doTDR) {
      if(st) {
        TPaveText* text = st->tdrText(0.8,0.45,0.9,0.25);
        text->AddText(Form("KSTest   = %5.4g", totalData->KolmogorovTest(all)));
        text->AddText(Form("Chi2/NDF = %5.4g", chi2NDF));
        text->Draw("same");
      }
      else {
        Tl.SetNDC();
        Tl.SetTextFont(42);
        Tl.SetTextSize(0.10546875);
        Tl.DrawLatex(0.60, 0.49, Form("KS       = %5.3g", totalData->KolmogorovTest(all)));
        Tl.DrawLatex(0.60, 0.39, Form("#chi^{2}/NDF = %5.3g", chi2NDF));
      }
   }
   else{
     Tl.DrawLatex(x, y*1.00, Form("KSTest   = %5.4g", totalData->KolmogorovTest(all)));
     Tl.DrawLatex(x, y*0.92, Form("Chi2Prob = %5.4g", all->Chi2TestX(totalData,chi2,NDF,igood,"WW")));
     Tl.DrawLatex(x, y*0.84, Form("Chi2/NDF = %5.4g", chi2/NDF));
   }
}

//______________________________________________________________________________
void FormattedPlot::drawBinInfo() {
  TLatex latex; 
  latex.SetNDC();    
  latex.SetTextFont(42);                                     
  latex.SetTextSize(0.045);                                
  latex.DrawLatex(0.19, 0.72, Form("%s, %s",DEFS::getJetBinLabel(jetBin).c_str(),DEFS::getTagCatLabel(tagCat).c_str()));
  if(controlRegion == DEFS::HighKinBDT || controlRegion == DEFS::LowKinBDT) {
    latex.DrawLatex(0.19, 0.65, Form("BDT%s%.2f",(controlRegion==DEFS::HighKinBDT) ? ">" : "<",
                                     DefaultValues::getMedianPurity(jetBin,leptonCat,"KinBDT").first));
  }
}

//______________________________________________________________________________
void FormattedPlot::drawLumi(float intLum)
{                                
   TLatex latex; 
   latex.SetNDC();                                         
   latex.SetTextSize(0.063);                                
   latex.SetTextAlign(31); // align right                  
   latex.DrawLatex(0.97,0.936,"#sqrt{s} = 8 TeV");
   if (intLum > 0.)
   {
      latex.SetTextAlign(11); // align right
      latex.SetTextSize(0.043);
      latex.DrawLatex(0.5,0.936,"#int ");
      latex.SetTextSize(0.063);
      latex.DrawLatex(0.525,0.936,Form("L dt = %3.2f fb^{-1}",intLum));
   }
   latex.SetTextSize(0.063);
   latex.SetTextAlign(11); // align left
   latex.DrawLatex(0.113,0.936,"CMS preliminary 2012");
}

//______________________________________________________________________________
TGraphAsymmErrors* FormattedPlot::makeSystematicErrors(TH1* nominal, double totalNormalizationError) {
  TGraphAsymmErrors* ret = new TGraphAsymmErrors();
  for(int ibin=1; ibin<=nominal->GetNbinsX(); ibin++) {
    int ipoint = ibin-1;
    double errorYLow = sqrt(pow(nominal->GetBinError(ibin),2)+(pow(totalNormalizationError,2)*pow(nominal->GetBinContent(ibin),2)));
    double errorYHigh = sqrt(pow(nominal->GetBinError(ibin),2)+(pow(totalNormalizationError,2)*pow(nominal->GetBinContent(ibin),2)));
    ret->SetPoint(ipoint,nominal->GetBinCenter(ibin),nominal->GetBinContent(ibin));
    ret->SetPointError(ipoint,nominal->GetBinWidth(ibin)/2.0,nominal->GetBinWidth(ibin)/2.0,errorYLow,errorYHigh);
  }
  return ret;
}

//______________________________________________________________________________
void FormattedPlot::addSystematicErrors(TGraphAsymmErrors* g, TH1* nominal, TH1* sysUp, TH1* sysDown, string type)
{
  if(g->GetN()!=nominal->GetNbinsX()) {
    cout << "ERROR::Plots::addSystematicErrors The number of points in the error graph (" << g->GetN()
         << ") does not equal the number of bins in the nominal histogram (" << nominal->GetNbinsX() << ")." << endl;
    assert(g->GetN()==nominal->GetNbinsX());
  }

  bool printWarning = false;
  for(int ibin=1; ibin<=nominal->GetNbinsX(); ibin++) {
    int ipoint = ibin-1;
    double errorYLow = 0;
    double errorYHigh = 0;
    if(sysUp->GetEntries()>0 && sysDown->GetEntries()>0) {
      double errorPlus = sysUp->GetBinContent(ibin)-nominal->GetBinContent(ibin);
      double errorMinus = nominal->GetBinContent(ibin)-sysDown->GetBinContent(ibin);
      if(errorPlus > 0)  errorYHigh += pow(errorPlus,2);
      else               errorYLow += pow(errorPlus,2);
      if(errorMinus > 0) errorYLow += pow(errorMinus,2);
      else               errorYHigh += pow(errorMinus,2);
      errorYHigh = sqrt(pow(g->GetErrorYhigh(ipoint),2)+errorYHigh);
      errorYLow = sqrt(pow(g->GetErrorYlow(ipoint),2)+errorYLow);
      //errorYHigh = sqrt(pow(g->GetErrorYhigh(ipoint),2)+pow(TMath::Abs(sysDown->GetBinContent(ibin)-nominal->GetBinContent(ibin)),2));
      //errorYLow = sqrt(pow(g->GetErrorYlow(ipoint),2)+pow(TMath::Abs(sysUp->GetBinContent(ibin)-nominal->GetBinContent(ibin)),2));
      g->SetPointError(ipoint,nominal->GetBinWidth(ibin)/2.0,nominal->GetBinWidth(ibin)/2.0,errorYLow,errorYHigh);
    }
    else {
       printWarning = true;
    }
  }

  if(printWarning) {
     cout << "WARNING::Plots::addSystematicErrors Can't add systematic errors to band (" << type << ")." << endl;
  }
}

//______________________________________________________________________________
int FormattedPlot::RoundToNearest(int iNumberToRound, int iToNearest) {
  int iNearest = 0;
  bool bIsUpper = false;
 
  int iRest = iNumberToRound % iToNearest;
  if (iNumberToRound == 550) bIsUpper = true;
 
  if (bIsUpper == true) {
    iNearest = (iNumberToRound - iRest) + iToNearest;
    return iNearest;
  }
  else if (iRest >= (iToNearest/2)) {
    iNearest = (iNumberToRound - iRest) + iToNearest;
    return iNearest;
  }
  else if (iRest < (iToNearest/2)) {
    iNearest = (iNumberToRound - iRest);
    return iNearest;
  }
  
  return -1;
}
