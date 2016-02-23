#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"
#include "TAMUWW/Tools/interface/Plots.hh"
#include "TAMUWW/Tools/interface/PlotFiller.hh"
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"

// ROOT libraries
#include "TROOT.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TBenchmark.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF2.h"
#include "TH2D.h"

// C++ libraries
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <stdlib.h>

using namespace std;
using DEFS::LeptonCat;

typedef PlotFiller::MapOfPlots MapOfPlots;

////////////////////////////////////////////////////////////////////////////////
//  Local Functions
////////////////////////////////////////////////////////////////////////////////

/// returns a map containing all of the plots that will be made for each process and their specific attributes
MapOfPlots getPlots(DEFS::JetBin jetBin, DEFS::LeptonCat leptonCat);

/// prints a list of PhysicsProcesses
void printProcesses(vector <PhysicsProcess*> procs);

/// prints a list of histograms loaded into a plot

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
int main(int argc,char**argv) {
 
   // evaluate command-line / configuration file options
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;

   string 			ifilename			 = cl.getValue<string>	  ("ifilename");
   string 			idirname			 = cl.getValue<string>	  ("idirname",					"");
   string 			alternateDataFile	 = cl.getValue<string>	  ("alternateDataFile",			"");
   string 			ofilepath			 = cl.getValue<string> 	  ("ofilepath",				  "./");
   string 			suffix				 = cl.getValue<string>	  ("suffix",					"");
   string           lepCat               = cl.getValue<string>    ("lep",                   "both");
   DEFS::LeptonCat  leptonCat            = DEFS::getLeptonCat     (lepCat);
   string           ntype                = cl.getValue<string>    ("ntype",          "EventNtuple");
   DEFS::NtupleType ntupleType           = DEFS::getNtupleType    (ntype);
   string           jBin                 = cl.getValue<string>    ("jBin",                 "jets2");      
   DEFS::JetBin     jetBin               = DEFS::getJetBin        (jBin);
   string           tcat                 = cl.getValue<string>    ("tcat",                "pretag");
   DEFS::TagCat     tagCat               = DEFS::getTagCat(tcat);
   bool             include_data         = cl.getValue<bool>      ("include_data",            true);
   bool             include_systematics  = cl.getValue<bool>      ("include_systematics",     true);
   bool				debug				 = cl.getValue<bool> 	  ("debug",					 false);

   if (!cl.check()) return 0;
   cl.print();

   // Trying to speed up the code
   gEnv->SetValue("TFile.AsyncPrefetching", 1);

   // Check that the leptcat actually exists
   if (leptonCat == DEFS::none) {
      cout<<"plotter_x called with unknown lepton category "<<lepCat<<endl;
      return 1;
   }

   // The vector holding all processes.
   vector <PhysicsProcess*> procs = DefaultValues::getProcessesHiggs(jetBin,tagCat,include_data,include_systematics,
                                                                     true,ntupleType,leptonCat);

   if(debug)
   	 printProcesses(procs);

   // The vector containing all plots to be made
   MapOfPlots plots = getPlots(jetBin, leptonCat);

   // Open the input file
   TFile* ifile = TFile::Open(ifilename.c_str(),"READ");
   TDirectoryFile* idir = (TDirectoryFile*)ifile->GetDirectory(idirname.c_str());

   // Open the input file for data if necessary
   TFile* ifile_data(0);
   TDirectoryFile* data_dir(0);
   if(alternateDataFile!="") {
   	  ifile_data = TFile::Open(alternateDataFile.c_str(),"READ");
   	  data_dir = (TDirectoryFile*)ifile_data->GetDirectory("");
   }


   for (MapOfPlots::iterator p = plots.begin() ; p != plots.end() ; p++) {
      for (map<string,  Plot * >::iterator p2 = p->second.begin(); p2 != p->second.end(); p2++) {
      	 bool uniformBinWidth = true;
      	 double binWidth = p2->second->templateHisto->GetXaxis()->GetBinWidth(1);
      	 for(int ibin=1; ibin<=p2->second->templateHisto->GetXaxis()->GetNbins(); ibin++) {
      	 	if(binWidth!=p2->second->templateHisto->GetXaxis()->GetBinWidth(ibin))
      	 		uniformBinWidth = false;
      	 }

      	 double original_bin_width = p2->second->templateHisto->GetXaxis()->GetBinWidth(1);
         p2->second->loadHistogramsFromFile(idir, procs, leptonCat, idirname!="" ? true : false, data_dir);
         p2->second->setScaled(true);
         double bin_width_ratio = p2->second->templateHisto->GetXaxis()->GetBinWidth(1)/original_bin_width;

         if(uniformBinWidth) {
	         dynamic_cast<FormattedPlot*>(p2->second)->range.first *= bin_width_ratio;
		     dynamic_cast<FormattedPlot*>(p2->second)->range.second *= bin_width_ratio;
		 }
		 else {
		 	dynamic_cast<FormattedPlot*>(p2->second)->range.first = p2->second->templateHisto->GetXaxis()->GetXmin();
		 	dynamic_cast<FormattedPlot*>(p2->second)->range.second = p2->second->templateHisto->GetXaxis()->GetXmax();
		 }
      }
   }

    for(MapOfPlots::iterator p = plots.begin(); p != plots.end(); p++) {
    	for(map<string,  Plot * >::iterator p2 = p->second.begin(); p2 != p->second.end() ; p2++){
    		if(debug)
    			p2->second->printList();
    	    TCanvas* can = ((FormattedPlot*) p2->second)->getCanvasTDR(procs);
    	    TString canName = ofilepath+"/"+can->GetName();
    	    canName += "_"+DEFS::getLeptonCatString(leptonCat);
    	    cout << "\tSaving canvas " << canName << suffix << " ... ";
    	    can->SaveAs(canName+suffix.c_str()+".png");
    	    can->SaveAs(canName+suffix.c_str()+".root");
    	    //TString class_name = ((FormattedPlot*)p2->second)->templateHisto->ClassName();
    	    //if(!class_name.Contains("TH2")) {
    	    //   can->SaveAs(canName+".pdf");
    	    //}
    	    cout << "DONE" << endl << flush;
    	}
    }

    return 0;
}

////////////////////////////////////////////////////////////////////////////////
//  Implement Local Functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
MapOfPlots getPlots(DEFS::JetBin jetBin, DEFS::LeptonCat leptonCat){
   MapOfPlots plots;

   FormattedPlot* a = new FormattedPlot;

      //Goes in the label and tells us whether we are looking at electrons or muons
   TString lepStr = "_"+DEFS::getLeptonCatString(leptonCat);

   // The overlay of a scaled signal. For signalName pick the groupingName 
   // of one of the processes. Or just "" if you don't want a signal overlayed.
   //TString signalName = "ggH+WH+qqH(125)";
   TString signalName = "H(125)->WW";
   double signalFactor = 1.0;
   if(leptonCat == DEFS::electron) signalFactor = 3900;
   else if(leptonCat == DEFS::muon) signalFactor = 3900;
   else signalFactor = 3900;
   TString name = "";
   bool norm_data = false;

//   a = new FormattedPlot;
//   name = "MET";
//   a->templateHisto = new TH1D(name + lepStr, name,100,0,500);
//   a->axisTitles.push_back("Missing E_{T} [GeV]");
//   a->axisTitles.push_back("Number of Events / 5 GeV");
//   a->range = make_pair(30.,150.);
//   a->normToData = norm_data;
//   a->stacked = true; a->leptonCat = leptonCat;
//   a->overlaySignalName = signalName;
//   a->overlaySignalFactor = signalFactor;
//   plots[leptonCat][string(name)] = a;
//   
//   a = new FormattedPlot;
//   name = "KinBDT";
//   a->templateHisto = new TH1D(name + lepStr, name,
//                               //100,-1,1);
//                               DEFS::getNBinsX(jetBin,leptonCat,string(name)),
//                               &DEFS::getBinsX(jetBin,leptonCat,string(name))[0]);
//   a->axisTitles.push_back("KinBDT");
//   a->axisTitles.push_back("Number of Events");
//   a->range = make_pair(-1.0,1.0);
//   a->normToData = norm_data;
//   a->stacked = true; a->leptonCat = leptonCat;
//   a->overlaySignalName = signalName;
//   a->overlaySignalFactor = signalFactor;
//   plots[leptonCat][string(name)] = a;
//
//   a = new FormattedPlot;
//   name = "MEBDT";
//   a->templateHisto = new TH1D(name + lepStr, name,
//                               //100,-1,1);
//                               DEFS::getNBinsX(jetBin,leptonCat,string(name)),
//                               &DEFS::getBinsX(jetBin,leptonCat,string(name))[0]);
//   a->axisTitles.push_back("MEBDT");
//   a->axisTitles.push_back("Number of Events");
//   a->range = make_pair(-1.0,1.0);
//   a->normToData = norm_data;
//   a->stacked = true; a->leptonCat = leptonCat;
//   a->overlaySignalName = signalName;
//   a->overlaySignalFactor = signalFactor;
//   plots[leptonCat][string(name)] = a;
//
   a = new FormattedPlot;
   name = "KinMEBDT";
   a->templateHisto = new TH1D(name + lepStr, name,
                               //100,-1,1);
                               DEFS::getNBinsX(jetBin,leptonCat,string(name)),
                               &DEFS::getBinsX(jetBin,leptonCat,string(name))[0]);
   a->axisTitles.push_back("KinMEBDT");
   a->axisTitles.push_back("Number of Events");
   a->range = make_pair(-1.0,1.0);
   a->normToData = norm_data;
   a->stacked = true; a->leptonCat = leptonCat;
   a->overlaySignalName = signalName;
   a->overlaySignalFactor = signalFactor;
   plots[leptonCat][string(name)] = a;

   // return all the plots to be made
   return plots;
}

//______________________________________________________________________________
void printProcesses(vector <PhysicsProcess*> procs) {
   cout << "LIST OF PROCESSES:" << endl;
   for(unsigned int i=0; i<procs.size(); i++) {
   	cout << "\t" << procs[i]->name << endl;
   }
   cout << endl << endl << endl;
}