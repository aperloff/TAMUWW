// Our libraries
#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"
#include "TAMUWW/SpecialTools/interface/Table.hh"
#include "TAMUWW/SpecialTools/interface/TableRow.hh"
#include "TAMUWW/SpecialTools/interface/TableCellVal.hh"
#include "TAMUWW/SpecialTools/interface/Value.hh"
#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/Tools/interface/Plots.hh"
#include "TAMUWW/Tools/interface/PlotFiller.hh"
#include "TAMUWW/Tools/interface/PhysicsProcessNEW.hh"
#include "TAMUWW/Tools/interface/PUreweight.hh"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"

// ROOT libraries
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2D.h"
#include "TChain.h"

// C++ libraries
#include <iostream>
#include <map>
#include <string>
#include <vector>

using namespace std;
using DEFS::LeptonCat;

namespace UserFunctions
{
   TString cutRegion;
   DEFS::LeptonCat leptCat;
   PUreweight* puweight;
   bool doJetEnergyResolution;
   bool doPileUpReWeighting;

   // Is run once for each process before events are cut (initialize)
   void initEventFunc(EventNtuple* ntuple, const PhysicsProcessNEW* proc);
   // this function fills all of the plots for a given process
   void fillPlots(map<string, Plot*> &  plots, EventNtuple * ntuple, double weight = 1.0);
   // returns a boolean if the event passes the specified cuts
   bool eventPassCuts(EventNtuple * ntuple, const PhysicsProcessNEW*);
   // returns a double
   double weightFunc(EventNtuple* ntuple, const PhysicsProcessNEW* proc);
   // Is run once for each process before events (initializes PU Reweighting
   void processFunc(EventNtuple* ntuple, const PhysicsProcessNEW* proc);
}

//______________________________________________________________________________
void UserFunctions::fillPlots(map<string, Plot*> &  plots, EventNtuple * ntuple, double weight)
{
   plots["Mjj"]->Fill(ntuple->Mjj,weight);
   plots["LeptPt"]->Fill(ntuple->lLV[0].Pt(),weight);
   plots["LeptEta"]->Fill(ntuple->lLV[0].Eta(),weight);
   plots["LeptPhi"]->Fill(ntuple->lLV[0].Phi(),weight);
   plots["MET"]->Fill(ntuple->METLV[0].Et(),weight);
   plots["WmT"]->Fill(sqrt(pow(ntuple->lLV[0].Et()+ntuple->METLV[0].Et(), 2) - pow(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(), 2) - pow(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(), 2)),weight);
   plots["Jet1Pt"]->Fill(ntuple->jLV[0].Pt(),weight);
   plots["Jet1Eta"]->Fill(ntuple->jLV[0].Eta(),weight);
   plots["Jet1Phi"]->Fill(ntuple->jLV[0].Phi(),weight);
   plots["Jet2Pt"]->Fill(ntuple->jLV[1].Pt(),weight);
   plots["Jet2Eta"]->Fill(ntuple->jLV[1].Eta(),weight);
   plots["Jet2Phi"]->Fill(ntuple->jLV[1].Phi(),weight);
   plots["DeltaEtaJ1J2"]->Fill(TMath::Abs(ntuple->jLV[0].Eta() - ntuple->jLV[1].Eta()),weight);
   plots["Ptjj"]->Fill((ntuple->jLV[0] + ntuple->jLV[1]).Pt(),weight);
   plots["j1Pt_Mjj"]->Fill(ntuple->jLV[0].Pt() / ntuple->Mjj,weight);
   plots["j2Pt_Mjj"]->Fill(ntuple->jLV[1].Pt() / ntuple->Mjj,weight);
   plots["Mlvjj"]->Fill((ntuple->jLV[0] + ntuple->jLV[1] + ntuple->lLV[0] + ntuple->METLV[0]).M(),weight);
   plots["DeltaRLepMET"]->Fill(sqrt(pow(ntuple->lLV[0].Eta()-ntuple->METLV[0].Eta(),2)+pow(ntuple->lLV[0].Phi()-ntuple->METLV[0].Phi(),2)),weight);
   plots["EJ1EJ2"]->Fill(ntuple->jLV[0].E() * ntuple->jLV[1].E(),weight);
   plots["BetaJ1BetaJ2"]->Fill(ntuple->jLV[0].Beta() * ntuple->jLV[1].Beta(),weight);
   plots["DeltaRJ1J2"]->Fill(sqrt(pow(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta(),2)+pow(ntuple->jLV[0].Phi()-ntuple->jLV[1].Phi(),2)),weight);
   plots["AngleJ1J2"]->Fill(ntuple->jLV[0].Angle(ntuple->jLV[1].Vect()),weight);
   plots["DeltaPhi_LJ1"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]),weight);
   plots["DeltaPhi_J1J2"]->Fill(ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]),weight);   
   plots["DeltaPhi_LJ1_vs_J1J2"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]), ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]),weight);
   plots["npv"]->Fill(ntuple->npv,weight);
   plots["jjlvPhi"]->Fill((ntuple->jLV[0] + ntuple->jLV[1]).Phi() - (ntuple->lLV[0] + ntuple->METLV[0]).Phi(),weight);

}//fillPlots

//______________________________________________________________________________
// Return true if the event pass the cuts imposed to the given lepton category
bool UserFunctions::eventPassCuts(EventNtuple * ntuple, const PhysicsProcessNEW*){

  /*
  ntuple->leptonCat  |    req   | result
          muon       |   muon   | continue
          electron   |   muon   | false
          both       |   muon   | continue
          muon       | electron | false
          electron   | electron | continue
          both       | electron | continue
          muon       |   both   | continue
          electron   |   both   | continue
          both       |   both   | continue
  */
  // Kill onlye electron-muon combinations, nothing else
  if ( ntuple->leptonCat != leptCat && (leptCat != DEFS::both && ntuple->leptonCat != DEFS::both) )
    return false;
  
 /* if (ntuple->leptonCat == DEFS::muon)
    if ((ntuple->lLV[0].Pt()) <= 25.0 || (ntuple->METLV[0].Et()) <= 30.0)
      return false;
  
  if (ntuple->leptonCat == DEFS::electron) 
    if ((ntuple->lLV[0].Pt()) <= 30.0 || (ntuple->METLV[0].Et()) <= 35.0)
      return false;
 */

  if (cutRegion.Contains("signal")){
    if(ntuple->jLV[0].Et() < 30 || ntuple->jLV[1].Et() < 30)
      return false;
  }
  else if (cutRegion.Contains("control1")){
    if (ntuple->jLV[0].Et() > 30 || ntuple->jLV[1].Et() > 30)
      return false;
  }
  else if (cutRegion.Contains("control2")){
    if ((ntuple->jLV[0].Et() > 30 && ntuple->jLV[1].Et() > 30) || 
	(ntuple->jLV[0].Et() < 30 && ntuple->jLV[1].Et() < 30))
      return false;
  }
  else if (cutRegion.Contains("all")){
    return true;
  }

  //Implement FNAL cuts
  /*  if (ntuple->leptonCat == DEFS::muon)
    return ntuple->FNALcutsMuon();
  else if(ntuple->leptonCat == DEFS::electron)
    return ntuple->FNALcutsElectron();
  */

  // return true otherwise
  return true;  

}// eventPassCuts

double UserFunctions::weightFunc(EventNtuple* ntuple, const PhysicsProcessNEW* proc)
{
   TString auxName = proc->name;
   auxName.ToUpper();
   if(auxName.Contains("DATA") || auxName.Contains("QCD"))
   {
      return 1.0;
   }
   
   double weight = 1.0;

   // Pileup reweighting
   if (doPileUpReWeighting)
   weight *= puweight->getWeight(ntuple->tnpus[1]);
   return weight;
   
}

void UserFunctions::initEventFunc(EventNtuple* ntuple, const PhysicsProcessNEW* proc)
{
  TString auxName = proc->name;
  auxName.ToUpper();
  if(auxName.Contains("DATA") || auxName.Contains("QCD"))
    return;

   // Performs JER
  if (doJetEnergyResolution)
    ntuple->doJER();
}

void UserFunctions::processFunc(EventNtuple* ntuple, const PhysicsProcessNEW* proc)
{
 
  TString auxName = proc->name;
  auxName.ToUpper();
  if(auxName.Contains("DATA") || auxName.Contains("QCD"))
    return;
   
  // Initializes PU Reweighting
  string dataname = string("/uscms/home/amejia94/CMSSW_5_2_3_patch4/src/pileup12_noTrig.root");
  string MCname   = string(TString("/uscms/home/aperloff/nobackup/PS_outfiles_20120621_NTUPLES/"+proc->name+"/PS.root"));
  puweight = new PUreweight(dataname,MCname,"pileup_IsoMu24_eta2p1","PS/TPUDist");
}

///  fills the histograms and controls the output canvas and file for the rest of the program
void doPlotter(TString fileName, map<string, Plot*> & plots, vector<PhysicsProcessNEW*> procs, bool doJER, bool doPUrewt);

/// returns the cross section for the given process
double getCrossSection(TString channelName);

/// returns the number of Monte Carlo events for the given process
double getNumMCEvts(TString channelName);

/// returns a vector containing all of the processes that the program will use
vector<PhysicsProcessNEW*> getProcesses(DEFS::LeptonCat leptCat, double intLum);

/// returns a vector containing all of the processes that the program will use
vector<PhysicsProcessNEW*> getProcesses(vector<string>& processNames, double intLum);

/// returns a map containing all of the plots that will be made for each process and their specific attributes
map<string, Plot*> getPlots(string leptonCatString);

/// returns the Color_t for a specific process
Color_t getProcessColor(TString channelName);

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
int main(int argc,char**argv) {
 
   // evaluate command-line / configuration file options
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;

   vector<string> processNames      = cl.getVector<string> ("procs",           "");
   string         leptonCatString   = cl.getValue<string>  ("lep",         "both");
   double         intLum            = cl.getValue<double>  ("lum",           3600);
   bool           doJER             = cl.getValue<bool>    ("doJer",         true);
   bool           doPUrewt          = cl.getValue<bool>    ("doPUrewt",      true);
   TString        cutRegionCL       = cl.getValue<string>  ("cutRegion", "signal");

   if (!cl.check()) return 0;
   cl.print();
   
   UserFunctions::doJetEnergyResolution = doJER;
   UserFunctions::doPileUpReWeighting = doPUrewt;

   // The enum telling the program to run over electrons, muons, or both
   UserFunctions::leptCat = DEFS::getLeptonCat(leptonCatString);
   UserFunctions::cutRegion = cutRegionCL;
   
   // The vector containing all plots to be made
   map<string, Plot*> plots = getPlots(leptonCatString);

   cout << "Doing lepton category "<<DEFS::getLeptonCatString(UserFunctions::leptCat)<< endl;

   // The vector holding all processes.
   vector <PhysicsProcessNEW*> procs;
   
   if (!processNames.empty())
     procs = getProcesses(processNames,intLum);
   else
     procs = getProcesses(UserFunctions::leptCat,intLum);
      
   // Make all the plots and store to outputFile
   TString outFileName = "outputFile_";
   outFileName += DEFS::getLeptonCatString(UserFunctions::leptCat)+".root";
   doPlotter(outFileName, plots, procs, doJER, doPUrewt);
   
   plots["DeltaPhi_LJ1_vs_J1J2"]->saveHistogramsToFile("DeltaPhi_LJ1_vs_J1J2_single.root");

   return 0;

}//plotter()


////////////////////////////////////////////////////////////////////////////////
//  Implement Local Functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void doPlotter(TString fileName, map<string, Plot*> & plots, vector<PhysicsProcessNEW*> procs, bool doJER, bool doPUrewt) {

   PlotFiller pFill(plots, procs, &UserFunctions::fillPlots);
   pFill.setCutFunction(&UserFunctions::eventPassCuts);
   pFill.setWeightFunction(&UserFunctions::weightFunc);
   pFill.setProcessFunction(&UserFunctions::processFunc);
   pFill.setInitializeEventFunction(&UserFunctions::initEventFunc);
   pFill.setMaximumEventsDEBUG(2000); // TEST
   pFill.run();

   // Will all the info in the plots get the canvas and write it to file
   TFile * outFile = new TFile(fileName,"RECREATE");
   for ( map<string, Plot*>::iterator p = plots.begin(); p != plots.end() ; p++){
     
     TCanvas* can = ((FormattedPlot*) p->second)->getCanvas(procs);
     can->Write();			
     
     can->SaveAs(TString(can->GetName())+".png");
     can->SaveAs(TString(can->GetName())+".eps");
   }
  
   outFile->Close();

}//plotter(.., .., etc)

//______________________________________________________________________________
double getCrossSection(TString channelName)
{
   Table table;
   double xsec;

   string basePath = gSystem->pwd();
   basePath = basePath.substr(0,basePath.find("TAMUWW"));
   table.parseFromFile(basePath+string("/TAMUWW/ConfigFiles/Official/CrossSections_8TeV.txt"),"TableCellVal");
   TableCell * cell = table.getCellRowColumn(string(channelName),"CrossSection");
   if(cell){
      xsec = ((TableCellVal*)cell)->val.value;
      if (xsec==0)
         cout << "WARNING::getCrossSection::The cross section for " << channelName << " is 0.0 +/- 0.0" << endl;
      return xsec;
   } else{
      cout << "WARNING::getCrossSection::channelName not recognized. Returning -1 for the cross section." << endl << "The events will have the same scale as the MC sample, but on a negative scale." << endl << "Please check channel names." << endl;
      return -1.;
   }
}//getCrossSection


//______________________________________________________________________________
double getNumMCEvts(TString channelName)
{
   Table table;
   double value;

   string basePath = gSystem->pwd();
   basePath = basePath.substr(0,basePath.find("TAMUWW"));
   table.parseFromFile(basePath+string("/TAMUWW/ConfigFiles/Official/EventsFromMC_commonPATTuple_525.txt"),"TableCellVal");
   TableCell * cell =table.getCellRowColumn(string(channelName),"Events_PATtuple");
   if(cell){
      value = ((TableCellVal*)cell)->val.value;
      if (value==0)
         cout << "WARNING::getNumMCEvts::There are 0 events in the " << channelName << " MC." << endl;
      return value;
   } else{
      cout << "WARNING::getNumMCEvts::channelName not recognized. Returning -1 event from MC." << endl << "Please check channel names." << endl;
      return -1.;
   }
}//getNumMCEvts


//______________________________________________________________________________
// leptonCan = The lepton category for which the processes are needed.
// intLum    = The integrated luminosity in pb-1 for this given lepton category 
vector<PhysicsProcessNEW*> getProcesses(DEFS::LeptonCat leptonCat, double intLum){
   
  string basePath = "/uscms/home/aperloff/nobackup/PS_outfiles_20120621_NTUPLES/";
  vector <PhysicsProcessNEW*> procs;
  
  
  procs.push_back(new ColoredPhysicsProcessNEW("WJets",basePath+"WJets/PS.root", getCrossSection("WJets"),
		   intLum, getNumMCEvts("WJets"), getProcessColor("WJets")));
  procs.push_back(new ColoredPhysicsProcessNEW("WW",basePath+"WW/PS.root", getCrossSection("WW"),
			     intLum, getNumMCEvts("WW"), getProcessColor("WW")));
  procs.push_back(new ColoredPhysicsProcessNEW("WZ",basePath+"WZ/PS.root", getCrossSection("WZ"),
		   intLum, getNumMCEvts("WZ"), getProcessColor("WZ")));
  procs.push_back(new ColoredPhysicsProcessNEW("DYJets",basePath+"DYJets/PS.root", getCrossSection("DYJets"),
			   intLum, getNumMCEvts("DYJets"), getProcessColor("DYJets")));
  procs.push_back(new ColoredPhysicsProcessNEW("TTbar",basePath+"TTbar/PS.root", getCrossSection("TTbar"),
                           intLum, getNumMCEvts("TTbar"), getProcessColor("TTbar")));
  //here 
  ///////////////   QCD new code
  //  procs.push_back(new ColoredPhysicsProcessNEW("QCD (Mu)", basePath+"QCD_MuEnriched/PS.root", getCrossSection("QCD"),
  //		   intLum, getNumMCEvts("QCD_Pt20_MuEnriched"), getProcessColor("QCD_MuEnriched")));
  
  procs.push_back(new ColoredPhysicsProcessNEW("STopT_T",basePath+"STopT_T/PS.root", getCrossSection("STopT_T"),
					       intLum, getNumMCEvts("STopT_T"), getProcessColor("STopT_T")));
  //   procs.push_back(new ColoredPhysicsProcessNEW("STopT_Tbar",basePath+"STopT_Tbar/PS.root",
  //                            getCrossSection("STopT_Tbar"), intLum, getNumMCEvts("STopT_Tbar"),
  //                            getProcessColor("STopT_Tbar")));
  procs.push_back(new ColoredPhysicsProcessNEW("STopS_T",basePath+"STopS_T/PS.root", getCrossSection("STopS_T"),
			   intLum, getNumMCEvts("STopS_T"), getProcessColor("STopS_T")));
  procs.push_back(new ColoredPhysicsProcessNEW("STopS_Tbar",basePath+"STopS_Tbar/PS.root",
			   getCrossSection("STopS_Tbar"), intLum, getNumMCEvts("STopS_Tbar"),
			   getProcessColor("STopS_Tbar")));
  procs.push_back(new ColoredPhysicsProcessNEW("STopTW_T",basePath+"STopTW_T/PS.root",
			   getCrossSection("STopTW_T"), intLum, getNumMCEvts("STopTW_T"),
			   getProcessColor("STopTW_T")));
  procs.push_back(new ColoredPhysicsProcessNEW("STopTW_Tbar",basePath+"STopTW_Tbar/PS.root", 
			   getCrossSection("STopTW_Tbar"), intLum, getNumMCEvts("STopTW_Tbar"),
			   getProcessColor("STopTW_Tbar")));
  
  // The Data
  // The normalization of the data is done in the same way as the MC. That is 
  // we multiply the histograms by getScaleFactor in class proc in PlotterAux.hh
  //  double getScaleFactor(){return sigma * intLum / initial_events;}
  // FOR DATA, to guarrantee the proper normalization, we set
  // sigma  = 1 / collected luminosity
  // initial_events = 1;

  if (leptonCat == DEFS::electron || leptonCat == DEFS::both)
    procs.push_back(new ColoredPhysicsProcessNEW("Data (el)",basePath+"SingleEl_Data/PS.root",
			     1./1599, intLum, 1,
			     getProcessColor("SingleEl_Data"))); 

  if (leptonCat == DEFS::muon || leptonCat == DEFS::both)
    procs.push_back(new ColoredPhysicsProcessNEW("Data (mu)",basePath+"SingleMu_Data/PS.root",
			     1./1606, intLum, 1,
			     getProcessColor("SingleMu_Data")));
 
  return procs;
  
}//getProcesses


//______________________________________________________________________________
vector<PhysicsProcessNEW*> getProcesses(vector<string>& processNames, double intLum){
   
  //  string basePath = "/uscms_data/d3/ilyao/Summer11/428Full/Stage2SkimResults/"; 
   string basePath = "/uscms/home/aperloff/nobackup/PS_outfiles_20120621_NTUPLES/"; //602 to 613
   vector <PhysicsProcessNEW*> procs;

   for (unsigned int i=0; i<processNames.size(); i++){
      if (TString(processNames[i]).Contains("Data")==0)
         if (processNames[i].compare("TTbar")==0)
            procs.push_back(new ColoredPhysicsProcessNEW(processNames[i],basePath+processNames[i]+"/PS.root",
                                     getCrossSection(processNames[i]), intLum, getNumMCEvts(processNames[i]),
                                     getProcessColor(processNames[i])));
         else
            procs.push_back(new ColoredPhysicsProcessNEW(processNames[i],basePath+processNames[i]+"/PS.root",
                                     getCrossSection(processNames[i]), intLum, getNumMCEvts(processNames[i]),
                                     getProcessColor(processNames[i])));
      else
         procs.push_back(new ColoredPhysicsProcessNEW(processNames[i],basePath+processNames[i]+"/PS.root",1 ,1, 1,
                                  getProcessColor(processNames[i]))); // 1's so we don't normalize data.
   }

   return procs;

}//getProcesses


//______________________________________________________________________________
map<string, Plot*> getPlots(string leptonCatString){

   map<string, Plot*> plots;

   FormattedPlot* a = new FormattedPlot;
   bool norm_data = true;
   double pi = TMath::Pi();

   //Goes in the label and tells us whether we are looking at electrons or muons
   TString lep = TString(leptonCatString);

   //goes in the label and tells us what cuts we are applying
   string cut;
   if (lep.CompareTo("muon") == 0)
      cut = "_MET > 30, muPt > 25";
   else
      cut = "_MET > 35, elPt > 30";
   TString cuts = TString(cut);

   TString suffix = lep;
   
   a = new FormattedPlot;

   a->templateHisto = new TH1D("Mjj_" + suffix,"Mjj_" + lep +  cuts,500,0,500);
   a->axisTitles.push_back("M_{jj} [GeV]");
   a->axisTitles.push_back("Number of Events / 10 GeV");
   a->range = make_pair(0.,400.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["Mjj"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("LeptPt_" + suffix,"LeptPt_" + lep + cuts,1000,0,500);
   a->axisTitles.push_back("p_{T}^{lepton} [GeV]");
   a->axisTitles.push_back("Number of Events / 5 GeV");
   a->range = make_pair(20.,150.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["LeptPt"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;
 
   a->templateHisto = new TH1D("LeptEta_" + suffix,"LeptEta_" + lep + cuts,50,-5,5);
   a->axisTitles.push_back("#eta^{lepton} [Radians]");
   a->axisTitles.push_back("Number of Events / 0.2 Radians");
   a->range = make_pair(-3.,3.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["LeptEta"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("LeptPhi_" + suffix,"LeptPhi_" + lep + cuts,70,-3.5,3.5);
   a->axisTitles.push_back("#phi^{lepton} [Radians]");
   a->axisTitles.push_back("Number of Events / 0.1 Radians");
   a->range = make_pair(-3.5,3.5);
   a->normToData = norm_data;
   a->stacked = true;
   plots["LeptPhi"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("MET_" + suffix,"MET_" + lep + cuts,1000,0,500);
   a->axisTitles.push_back("Missing E_{T} [GeV]");
   a->axisTitles.push_back("Number of Events / 5 GeV");
   a->range = make_pair(30.,150.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["MET"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("WmT_" + suffix,"WmT_" + lep + cuts,1000,0,500);
   a->axisTitles.push_back("M_{T}^{W} [GeV]");
   a->axisTitles.push_back("Number of Events / 5 GeV");
   a->range = make_pair(0.,150.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["WmT"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("Jet1Pt_" + suffix,"Jet1Pt_" + lep + cuts,200,0,300);
   a->axisTitles.push_back("p_{T}^{jet_{1}} [GeV]");
   a->axisTitles.push_back("Number of Events / 5 GeV");
   a->range = make_pair(20.,200.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["Jet1Pt"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("Jet1Eta_" + suffix,"Jet1Eta_" + lep + cuts,50,-5,5);
   a->axisTitles.push_back("#eta^{jet_{1}} [Radians]");
   a->axisTitles.push_back("Number of Events / 0.2 Radians");
   a->range = make_pair(-3.,3.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["Jet1Eta"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("Jet1Phi_" + suffix,"Jet1Phi_" + lep + cuts,70,-3.5,3.5);
   a->axisTitles.push_back("#phi^{jet_{1}} [Radians]");
   a->axisTitles.push_back("Number of Events / 0.1 Radians");
   a->range = make_pair(-3.5,3.5);
   a->normToData = norm_data;
   a->stacked = true;
   plots["Jet1Phi"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("Jet2Pt_" + suffix,"Jet2Pt_" + lep + cuts,200,0,300);
   a->axisTitles.push_back("p_{T}^{jet_{2}} [GeV]");
   a->axisTitles.push_back("Number of Events / 5 GeV");
   a->range = make_pair(20.,100.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["Jet2Pt"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("Jet2Eta_" + suffix,"Jet2Eta_" + lep + cuts,50,-5,5);
   a->axisTitles.push_back("#eta^{jet_{1}} [Radians]");
   a->axisTitles.push_back("Number of Events / 0.2 Radians");
   a->range = make_pair(-3.,3.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["Jet2Eta"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("Jet2Phi_" + suffix,"Jet2Phi_" + lep + cuts,70,-3.5,3.5);
   a->axisTitles.push_back("#phi^{jet_{1}} [Radians]");
   a->axisTitles.push_back("Number of Events / 0.1 Radians");
   a->range = make_pair(-3.5,3.5);
   a->normToData = norm_data;
   a->stacked = true;
   plots["Jet2Phi"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("DeltaEtaJ1J2_" + suffix,"DeltaEtaJ1J2_" + lep + cuts,50,0,5);
   a->axisTitles.push_back("#eta^{jet_{1}} - #eta^{jet_{2}} [Radians]");
   a->axisTitles.push_back("Number of Events / 0.1 Radians");
   a->range = make_pair(0.,5.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["DeltaEtaJ1J2"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("Ptjj_" + suffix,"Ptjj_" + lep + cuts,100,0,300);
   a->axisTitles.push_back("p_{T}^{jj} [GeV]");
   a->axisTitles.push_back("Number of Events / 10 GeV");
   a->range = make_pair(0.,250.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["Ptjj"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("j1Pt_Mjj_" + suffix,"j1Pt_Mjj_" + lep + cuts,500,0,5);
   a->axisTitles.push_back("p_{T}^{jet_{1}}/M_{jj} [GeV]");
   a->axisTitles.push_back("Number of Events / 0.01 GeV");
   a->range = make_pair(0.,2.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["j1Pt_Mjj"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("j2Pt_Mjj_" + suffix,"j2Pt_Mjj_" + lep + cuts,500,0,5);
   a->axisTitles.push_back("p_{T}^{jet_{2}}/M_{jj} [GeV]");
   a->axisTitles.push_back("Number of Events / 0.01 GeV");
   a->range = make_pair(0.,1.5);
   a->normToData = norm_data;
   a->stacked = true;
   plots["j2Pt_Mjj"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("Mlvjj_" + suffix,"Mlvjj_" + lep + cuts,250,0,1000);
   a->axisTitles.push_back("M_{lvjj} [GeV]");
   a->axisTitles.push_back("Number of Events / 10 GeV");
   a->range = make_pair(50.,800.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["Mlvjj"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("DeltaRLepMET_" + suffix,"DeltaRLepMET_" + lep + cuts,50,0,10);
   a->axisTitles.push_back("#DeltaR of lep and MET [Radians]");
   a->axisTitles.push_back("Number of Events / 0.2 Radians");
   a->range = make_pair(0.,7.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["DeltaRLepMET"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("EJ1EJ2_" + suffix,"EJ1EJ2_" + lep + cuts,500,700,5000);
   a->axisTitles.push_back("E_{J1} * E_{J2}  [GeV]");
   a->axisTitles.push_back("Number of Events / 10 GeV");
   a->range = make_pair(0.,5000.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["EJ1EJ2"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("BetaJ1BetaJ2_" + suffix,"BetaJ1BetaJ2_" + lep + cuts,10000,0,10);
   a->axisTitles.push_back("#beta_{J1} * #beta_{J2} [GeV]");
   a->axisTitles.push_back("Number of Events / .01 GeV");
   a->range = make_pair(0.9,1.03);
   a->normToData = norm_data;
   a->stacked = true;
   plots["BetaJ1BetaJ2"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("DeltaRJ1J2_" + suffix,"DeltaRJ1J2_" + lep + cuts,50,0,10);
   a->axisTitles.push_back("#DeltaR of J1 and J2 [Radians]");
   a->axisTitles.push_back("Number of Events / 0.2 Radians");
   a->range = make_pair(0.,7.);
   a->normToData = norm_data;
   a->stacked = true;
   plots["DeltaRJ1J2"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("AngleJ1J2_" + suffix,"AngleJ1J2_" + lep + cuts,50,0,5);
   a->axisTitles.push_back("Angle between J1 and J2 [Radians]");
   a->axisTitles.push_back("Number of Events / 0.1 Radians");
   a->range = make_pair(-0.5,3.5);
   a->normToData = norm_data;
   a->stacked = true;
   plots["AngleJ1J2"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("jjlvPhi_" + suffix,"jjlvPhi__" + lep + cuts,70,-3.5,3.5);
   a->axisTitles.push_back("#phi(jj) - #phi(lv)");
   a->axisTitles.push_back("Number of Events");
   a->range = make_pair(-3.5,3.5);
   a->normToData = norm_data;
   a->stacked = true;
   plots["jjlvPhi"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("DeltaPhi_LJ1_" + suffix,"DeltaPhi_LJ1__" + lep + cuts,50,-10,10);
   a->axisTitles.push_back("#delta #phi(LJ1)");
   a->axisTitles.push_back("Number of Events / .2 Radians");
   a->range = make_pair(-pi,pi);
   a->normToData = norm_data;
   a->stacked = true;
   plots["DeltaPhi_LJ1"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("DeltaPhi_J1J2_" + suffix,"DeltaPhi_J1J2__" + lep + cuts,50,-10,10);
   a->axisTitles.push_back("#delta #phi(J1J2)");
   a->axisTitles.push_back("Number of Events / .2 Radians");
   a->range = make_pair(-pi,pi);
   a->normToData = norm_data;
   a->stacked = true;
   plots["DeltaPhi_J1J2"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH1D("npv_" + suffix,"npv__" + lep + cuts,100,0,100);
   a->axisTitles.push_back("npv");
   a->axisTitles.push_back("Number of Events");
   a->range = make_pair(0,100);
   a->normToData = norm_data;
   a->stacked = true;
   plots["npv"] = a;

   a->axisTitles.clear();
   a = new FormattedPlot;

   a->templateHisto = new TH2D(("DeltaPhi_LJ1_vs_J1J2_" + suffix),("DeltaPhi_LJ1_vs_J1J2__" + lep + cuts), 15,-pi,pi,15,-pi,pi);
   a->axisTitles.push_back("#delta #phi(J1J2) vs. #delta #phi(LJ1)");
   a->axisTitles.push_back("Number of Events / .2 Radians");
   a->range = make_pair(-pi,pi);
   a->normToData = norm_data;
   a->stacked = true;
   plots["DeltaPhi_LJ1_vs_J1J2"] = a;

   // return all the plots to be made
   return plots;

}//getPlots


//______________________________________________________________________________
Color_t getProcessColor(TString channelName){

  if (channelName.CompareTo("WW") == 0)
    return kRed;
   else if (channelName.CompareTo("WZ") == 0)
     return kBlue;
   else if (channelName.CompareTo("WJets") == 0)
     return kGreen+1;
   else if (channelName.CompareTo("DYJets") == 0)
     return kMagenta;
  //   else if (channelName.CompareTo("TTbar")==0)
  //      return kCyan;
     else if (channelName.CompareTo("QCD_MuEnriched") == 0)
    return kCyan;  
  //right now these aren't being used because we are combining all the STop's in PlotterAux.cc
   else if (channelName.CompareTo("STopT_T") == 0)
     return kGreen;
   else if (channelName.CompareTo("STopT_Tbar") == 0)
     return kCyan+3;
   else if (channelName.CompareTo("STopS_T") == 0)
     return kBlue;
   else if (channelName.CompareTo("STopS_Tbar") == 0)
     return kBlue+3;
   else if (channelName.CompareTo("STopTW_T") == 0)
     return kMagenta;
   else if (channelName.CompareTo("STopTW_Tbar") == 0)
     return kGreen+3;
   else if (channelName.CompareTo("STopTW_Tbar") == 0)
     return kGreen+3;
   else if (channelName.CompareTo("TTbar") == 0)
     return kCyan; 
  
  //these are still being used
  
   else if (channelName.CompareTo("SingleEl_Data") == 0)
     return kBlack;
   else if (channelName.CompareTo("SingleMu_Data") == 0)
     return kBlack;
   else{
     cout << "WARNING Plotter::GetProcessColor() Unknown process name=|"<<channelName
	  <<"|. Returning process color as kYellow." << endl;
   }

  return kYellow;

}//getProcessColor
