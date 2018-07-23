// Name: FillPlotsTester_x.cc
// Author: Travis Lamb

// This is a basic demonstration/test of PlotFiller's functionality

#include "TString.h"
#include "TChain.h"
#include "TFile.h"
#include "TROOT.h"
#include "TSystem.h"
#include "THStack.h"

#include "TAMUWW/Tools/interface/PUreweight.hh"
#include "TAMUWW/SpecialTools/interface/Table.hh"
#include "TAMUWW/SpecialTools/interface/TableRow.hh"
#include "TAMUWW/SpecialTools/interface/TableCellVal.hh"
#include "TAMUWW/SpecialTools/interface/Value.hh"
#include "TAMUWW/SpecialTools/interface/Defs.hh"

//#include "TAMUWW/Tools/interface/PlotMap.hh"

#include "TAMUWW/Tools/interface/TriggerEfficiency.hh"
#include "TAMUWW/Tools/interface/Plots.hh"
#include "TAMUWW/Tools/interface/PlotFiller.hh"
#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"

#include <vector>

using namespace std;

TriggerEfficiency* triggerEfficiency;
DEFS::LeptonCat leptonCategory;
PUreweight* puweight;

void fillPlots(map<string, Plot*> &  plots, EventNtuple * ntuple, METree * metree, MicroNtuple * mnt, vector<TString>, double weight = 1.0);
bool eventPassCuts(EventNtuple* ntuple, const PhysicsProcess* proc);
double weightFunc(EventNtuple* ntuple, const PhysicsProcess* proc)
{
   TString auxName = proc->name;
   auxName.ToUpper();
   if(auxName.Contains("DATA") || auxName.Contains("QCD"))
   {
      return 1.0;
   }
   
   double weight = 1.0;
   // Trigger Effeciency reweighting
   //weight *= triggerEfficiency->getWeight(ntuple->leptonCat, ntuple->lLV[0].Pt(), ntuple->lLV[0].Eta());
   // Pileup reweighting
   weight *= puweight->getWeight(ntuple->vLV[0].tnpus[1]);
   
   return weight;
}
void processFunc(EventNtuple* ntuple, const PhysicsProcess* proc)
{
   TString auxName = proc->name;
   auxName.ToUpper();
   if(auxName.Contains("DATA") || auxName.Contains("QCD"))
      return;
   
   // Initializes PU Reweighting
   string dataname = string("/uscms/home/ssinjini/nobackup/CMSSW_5_2_5/src/pileup12_noTrig.root");
   string MCname   = string(TString("/uscms/home/aperloff/nobackup/PS_outfiles_20120621_NTUPLES/"+proc->name+"/PS.root"));
   puweight = new PUreweight(dataname,MCname,"pileup_IsoMu24_eta2p1","PS/TPUDist");
}
void initEventFunc(EventNtuple* ntuple, const PhysicsProcess* proc)
{
   TString auxName = proc->name;
   auxName.ToUpper();
   if(auxName.Contains("DATA") || auxName.Contains("QCD"))
      return;
   
   // Performs JER
   ntuple->doJER();
}

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
      cout << "WARNING::getNumMCEvts::" << channelName << " not recognized. Returning -1 event from MC." << endl << "Please check channel names." << endl;
      return -1.;
   }
}//getNumMCEvts

//______________________________________________________________________________
Color_t getProcessColor(TString channelName){

  if (channelName.CompareTo("WW") == 0)
    return kPink;
  else if (channelName.CompareTo("WZ") == 0)
    return kBlue;
  else if (channelName.CompareTo("WJets") == 0)
    return kTeal+2;
  else if (channelName.CompareTo("DYJets") == 0)
    return kPink-8; //kViolet-5;
  //   else if (channelName.CompareTo("TTbar")==0)
  //      return kCyan;
  else if (channelName.CompareTo("QCD_MuEnriched") == 0)
    return kCyan;  
  //right now these aren't being used because we are combining all the STop's in PlotterAux.cc
  else if (channelName.CompareTo("STopT_T") == 0)
    return kOrange+1;
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
     return kAzure-2; 
   else if (channelName.CompareTo("H125") == 0)
     return kRed+2; 
  
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

vector<PhysicsProcess*> getProcesses(vector<string>& processNames, double intLum){
   
   //  string basePath = "/uscms_data/d3/ilyao/Summer11/428Full/Stage2SkimResults/"; 
   string basePath = "/uscms/home/aperloff/nobackup/PS_outfiles_20120621_NTUPLES/"; //602 to 613
   vector <PhysicsProcess*> procs;

   for (unsigned int i=0; i<processNames.size(); i++){
      if (TString(processNames[i]).Contains("Data")==0)
         if (processNames[i].compare("TTbar")==0)
            procs.push_back(new PhysicsProcess(processNames[i],basePath+processNames[i]+"/PS.root",
                                                  getCrossSection(processNames[i]), intLum, getNumMCEvts(processNames[i])));
         else
            procs.push_back(new PhysicsProcess(processNames[i],basePath+processNames[i]+"/PS.root",
                                                  getCrossSection(processNames[i]), intLum, getNumMCEvts(processNames[i])));
      else
         procs.push_back(new PhysicsProcess(processNames[i],basePath+processNames[i]+"/PS.root", 1, 1, 1)); // 1's so we don't normalize data.
   }

   return procs;

}//getProcesses

vector<PhysicsProcess*> getProcesses(DEFS::LeptonCat leptonCat, double intLum){
   
  string basePath = "/uscms/home/aperloff/nobackup/PS_outfiles_20120621_NTUPLES/";
  vector <PhysicsProcess*> procs;
  

  procs.push_back(new PlotterPhysicsProcess("WJets",basePath+"WJets/PS.root", getCrossSection("WJets"),
					       1.045*intLum, getNumMCEvts("WJets"), getProcessColor("WJets")));
  procs.push_back(new PlotterPhysicsProcess("WW",basePath+"WW/PS.root", getCrossSection("WW"),
					       intLum, getNumMCEvts("WW"), getProcessColor("WW")));
  procs.push_back(new PlotterPhysicsProcess("WZ",basePath+"WZ/PS.root", getCrossSection("WZ"),
					       intLum, getNumMCEvts("WZ"), getProcessColor("WZ")));
  procs.push_back(new PlotterPhysicsProcess("DYJets",basePath+"DYJets/PS.root", getCrossSection("DYJets"),
					       intLum, getNumMCEvts("DYJets"), getProcessColor("DYJets")));
  procs.push_back(new PlotterPhysicsProcess("TTbar",basePath+"TTbar/PS.root", getCrossSection("TTbar"),
					       intLum, getNumMCEvts("TTbar"), getProcessColor("TTbar")));
  //here 
  ///////////////   QCD new code
  procs.push_back(new PlotterPhysicsProcess("QCD", basePath+"QCD/PS.root", getCrossSection("QCD_Pt20_MuEnriched"),
					       0.754*intLum, getNumMCEvts("QCD"), getProcessColor("QCD"))); //QCD_Pt20_MuEnriched
  procs.push_back(new PlotterPhysicsProcess("STopT_T",basePath+"STopT_T/PS.root", getCrossSection("STopT_T"),
					       intLum, getNumMCEvts("STopT_T"), getProcessColor("STopT_T")));
  //   procs.push_back(new PlotterPhysicsProcess("STopT_Tbar",basePath+"STopT_Tbar/PS.root",
  //                            getCrossSection("STopT_Tbar"), intLum, getNumMCEvts("STopT_Tbar"),
  //                            getProcessColor("STopT_Tbar")));
  procs.push_back(new PlotterPhysicsProcess("STopS_T",basePath+"STopS_T/PS.root", getCrossSection("STopS_T"),
			   intLum, getNumMCEvts("STopS_T"), getProcessColor("STopS_T")));
  procs.push_back(new PlotterPhysicsProcess("STopS_Tbar",basePath+"STopS_Tbar/PS.root",
			   getCrossSection("STopS_Tbar"), intLum, getNumMCEvts("STopS_Tbar"),
			   getProcessColor("STopS_Tbar")));
  procs.push_back(new PlotterPhysicsProcess("STopTW_T",basePath+"STopTW_T/PS.root",
			   getCrossSection("STopTW_T"), intLum, getNumMCEvts("STopTW_T"),
			   getProcessColor("STopTW_T")));
  procs.push_back(new PlotterPhysicsProcess("STopTW_Tbar",basePath+"STopTW_Tbar/PS.root", 
			   getCrossSection("STopTW_Tbar"), intLum, getNumMCEvts("STopTW_Tbar"),
			   getProcessColor("STopTW_Tbar")));
  procs.push_back(new PlotterPhysicsProcess("H125",basePath+"H125/PS.root", 
                                               getCrossSection("WH125"), intLum, getNumMCEvts("WH125"),
                                               getProcessColor("H125"), "PS/jets2"));
  
  // The Data
  // The normalization of the data is done in the same way as the MC. That is 
  // we multiply the histograms by getScaleFactor in class proc in PlotterAux.hh
  //  double getScaleFactor(){return sigma * intLum / initial_events;}
  // FOR DATA, to guarrantee the proper normalization, we set
  // sigma  = 1 / collected luminosity
  // initial_events = 1;

  if (leptonCat == DEFS::electron || leptonCat == DEFS::both)
    procs.push_back(new PlotterPhysicsProcess("SingleEl_Data",basePath+"SingleEl_Data/PS.root",
			     1./intLum, intLum, 1,
			     getProcessColor("SingleEl_Data"))); 

  if (leptonCat == DEFS::muon || leptonCat == DEFS::both)
    procs.push_back(new PlotterPhysicsProcess("SingleMu_Data",basePath+"SingleMu_Data/PS.root",
			     1./intLum, intLum, 1,
			     getProcessColor("SingleMu_Data")));
 
  return procs;
}

bool eventPassCuts(EventNtuple* ntuple, const PhysicsProcess*){

   /*
     ntuple->leptonCat |   req    | result
     muon              |   muon   | continue
     electron          |   muon   | false
     both              |   muon   | continue
     muon              | electron | false
     electron          | electron | continue
     both              | electron | continue
     muon              |   both   | continue
     electron          |   both   | continue
     both              |   both   | continue
   */
   // Kill only electron-muon combinations, nothing else
   if ( ntuple->lLV[0].leptonCat != leptonCategory && (leptonCategory != DEFS::both && ntuple->lLV[0].leptonCat != DEFS::both) )
      return false;
  
   /* OLD CUTS
      if (ntuple->leptonCat == DEFS::muon)
      if ((ntuple->lLV[0].Pt()) <= 25.0 || (ntuple->METLV[0].Et()) <= 30.0)
      return false;
     
      if (ntuple->leptonCat == DEFS::electron) 
      if ((ntuple->lLV[0].Pt()) <= 30.0 || (ntuple->METLV[0].Et()) <= 35.0)
      return false;
   */
   
   /* // NEWEST CUTS
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
   */
   
   
   /* // FNAL CUTS
      if (ntuple->leptonCat == DEFS::muon)
         return ntuple->FNALcutsMuon();
      else if(ntuple->leptonCat == DEFS::electron)
         return ntuple->FNALcutsElectron();
   */
   

   // return true otherwise
   return true;  

}// eventPassCuts

void fillPlots(map<string, Plot*> &  plots, EventNtuple * ntuple, METree * metree, MicroNtuple * mnt, vector<TString>, double weight)
{
   plots["Mjj"]->Fill(ntuple->Mjj,weight);
   plots["MET"]->Fill(ntuple->METLV[0].Et(),weight);
   plots["WmT"]->Fill(sqrt(pow(ntuple->lLV[0].Et()+ntuple->METLV[0].Et(), 2) - pow(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(), 2) - pow(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(), 2)),weight);
   return;
}

int main()
{
   triggerEfficiency = new TriggerEfficiency();
   triggerEfficiency->createEfficiencies8TeV();
   
   leptonCategory = DEFS::electron;
   
   puweight = 0;
   
   map<string, Plot*> plots;
   ///////////////////////
   // From plotter_x.cc //
   ///////////////////////
   FormattedPlot* a = new FormattedPlot;

   //string leptonCatString = "muon";
   string leptonCatString = "electron";
   
   bool norm_data = false;

   //Goes in the label and tells us whether we are looking at electrons or muons
   TString lep = TString(leptonCatString);
   DEFS::LeptonCat leptonCat = DEFS::getLeptonCat(leptonCatString);

   //goes in the label and tells us what cuts we are applying
   string cut;
   if (lep.CompareTo("muon") == 0)
      cut = "_MET > 30, muPt > 25";
   else
      cut = "_MET > 35, elPt > 30";
   TString cuts = TString(cut);

   TString suffix = lep;
   
   a->templateHisto = new TH1D("MET_" + suffix,"MET_" + lep + cuts,500,0,500);
   a->axisTitles.push_back("Missing E_{T} [GeV]");
   a->axisTitles.push_back("Number of Events / 5 GeV");
   a->range = make_pair(30.,150.);
   a->normToData = norm_data;
   a->stacked = true;
   a->overlaySignalName = "H125";
   a->overlaySignalFactor = 500;
   plots["MET"] = a;
   
   a = new FormattedPlot;
   
   a->templateHisto = new TH1D("Mjj_" + suffix,"Mjj_" + lep +  cuts,500,0,500);
   a->axisTitles.push_back("M_{jj} [GeV]");
   a->axisTitles.push_back("Number of Events / 10 GeV");
   a->range = make_pair(0.,400.);
   a->normToData = norm_data;
   a->stacked = true;
   a->overlaySignalName = "H125";
   a->overlaySignalFactor = 500;
   plots["Mjj"] = a;
   
   a = new FormattedPlot;
   
   a->templateHisto = new TH1D("WmT_" + suffix,"WmT_" + lep + cuts,500,0,500);
   a->axisTitles.push_back("M_{T}^{W} [GeV]");
   a->axisTitles.push_back("Number of Events / 5 GeV");
   a->range = make_pair(0.,150.);
   a->normToData = norm_data;
   a->stacked = true;
   a->overlaySignalName = "H125";
   a->overlaySignalFactor = 500;
   plots["WmT"] = a;
   /////////
   // End //
   /////////
   
   vector<PhysicsProcess*> procs;
   ///////////////////////
   // From plotter_x.cc //
   ///////////////////////
   vector<string> processNames;
   processNames.push_back("WJets");
   processNames.push_back("QCD");
   processNames.push_back("WW");
   processNames.push_back("WZ");
   processNames.push_back("DYJets");
   processNames.push_back("TTbar");
   processNames.push_back("STopT_T");
   processNames.push_back("STopS_T");
   processNames.push_back("STopS_Tbar");
   processNames.push_back("STopTW_T");
   processNames.push_back("STopTW_Tbar");
   processNames.push_back("SingleEl_Data");
   
   procs = getProcesses(leptonCat, 3600);
   /////////
   // End //
   /////////
   
   PlotFiller procMan(plots, procs, &fillPlots);
   procMan.setCutFunction(&eventPassCuts);
   procMan.setWeightFunction(&weightFunc);
   procMan.setProcessFunction(&processFunc);
   procMan.setInitializeEventFunction(&initEventFunc);
   //procMan.setMaximumEventsDEBUG(2000); // TEST
   procMan.run();
   
   plots["MET"]->scaleToData(procs);
   plots["Mjj"]->scaleToData(procs);
   plots["WmT"]->scaleToData(procs);
   
   TString fileName = DEFS::getLeptonCatString(leptonCat)+".root";
   
   // Write everything to a file
   TFile * outFile = new TFile(fileName,"RECREATE");
   for ( map<string, Plot*>::iterator p = plots.begin(); p != plots.end() ; p++)
   {
      TCanvas* can = ((FormattedPlot*) p->second)->getCanvas(procs);
      can->Write();
      
      THStack* stackedMC = new THStack(((p->first) + "_StackedMC").c_str(),((p->first) + "_StackedMC").c_str());
      
      p->second->SetName((p->first).c_str());
      p->second->Write();
      
      for(unsigned int i = 0; i < p->second->histos.size(); i++)
      {
         p->second->histos[i]->Write();
         
         TString histName = p->second->histos[i]->GetName();
         histName.ToUpper();
         if(!histName.Contains("DATA"))
            stackedMC->Add(p->second->histos[i]);
      }
      stackedMC->Write();
   }
  
   outFile->Close();
   
   return 0;
}
