
// Title:     PlotFiller.cc
// Author:    Travis Lamb
// Started:   22nd of June, 2012
// Last Edit: 22nd of June, 2012

// Uses Plot and PhysicsProcess to return a plot for a given process and yield.

// Our libraries
#include "TAMUWW/Tools/interface/PlotFiller.hh"

using namespace std;

//#####################################################################
//#####################################################################
//#####################################################################

PlotFiller::PlotFiller(MapOfPlots &plotsTemp,
                       vector<PhysicsProcess*> &procsTemp,
                       void (*userFillFuncTemp) (MapOfPlots &, TString, EventNtuple*, METree*, MicroNtuple*, vector<TString>, map<TString,MVAVar>&, double)):
   plots(plotsTemp),
   processes(procsTemp),
   userFillFunc(userFillFuncTemp)
{
   userWeightFunc = &defaultWeightFunc;
   userCutFunc = &defaultCutFunc;
   userProcessFunc = &defaultProcessFunc;
   userInitEventFunc = &defaultInitEventFunc;
   
   numberOfEvents = 0;
   
   debug = false;

   event_benchmark = new TBenchmark();
   event_benchmark->Reset();
   func_benchmark = new TBenchmark();
   func_benchmark->Reset();
}

PlotFiller::~PlotFiller()
{
   delete event_benchmark;
   delete func_benchmark;
}

void PlotFiller::setWeightFunction(double (*userWeightFuncTemp) (EventNtuple*, MicroNtuple*, const PhysicsProcess*))
{
   userWeightFunc = userWeightFuncTemp;
}

void PlotFiller::setCutFunction(bool (*userCutFuncTemp) (EventNtuple*, MicroNtuple*, const PhysicsProcess*, map<TString,MVAVar>&, vector<TString>))
{
   userCutFunc = userCutFuncTemp;
}

void PlotFiller::setProcessFunction(void (*userProcessFuncTemp) (EventNtuple*, const PhysicsProcess*))
{
   userProcessFunc = userProcessFuncTemp;
}

void PlotFiller::setInitializeEventFunction(void (*userInitEventFuncTemp) (EventNtuple*, const PhysicsProcess*))
{
   userInitEventFunc = userInitEventFuncTemp;
}

void PlotFiller::setMaximumEventsDEBUG(unsigned int maxEvts)
{
   debugNumberOfEvents = maxEvts;
   debug = true;
}

void PlotFiller::setMVAWeightDir(TString MVAWD)
{
   MVAWeightDir = MVAWD;
}

void PlotFiller::setMVAMethods(vector<TString> MVAM)
{
   MVAMethods = MVAM;
}

void PlotFiller::setMVAVar(vector<TString> mvav)
{
   MVAV = mvav;
}

void PlotFiller::setMVASpec(vector<TString> mvas)
{
   MVAS = mvas;
}

void PlotFiller::setLimitBranches(int lb)
{
   limitBranches = lb;
}

void PlotFiller::run()
{
   cout << "\nPlotFiller::run Begin filling plots" << endl;
   for(unsigned int i = 0; i < processes.size(); i++)
   {
      cout << "\nDoing Process " << processes[i]->name << endl;
      cout << "\tFrom file " << processes[i]->fileName << endl;
      // Tell all plots to prepare for filling 
      for (MapOfPlots::iterator p = plots.begin() ; p != plots.end() ; p++) {
         TString suffix = "_"+processes[i]->name + "_" + DEFS::getLeptonCatString(p->first);
         for (map<string,  Plot * >::iterator p2 = p->second.begin(); p2 != p->second.end(); p2++)
            p2->second->prepareToFillProcess(suffix, processes[i]->groupName);
      }//for 
      
      // Create the eventNtuple, METree, and MicroNtuple
      EventNtuple * ntuple = 0;
      METree * metree = 0;
      MicroNtuple * mnt = 0;
      TChain * c = processes[i]->chain;
      // Try to speed up processing
      //c->SetCacheSize(10000000);
      //c->AddBranchToCache("*");
      
      if(limitBranches == 0) {
         if(!processes[i]->name.Contains("TTbar",TString::kIgnoreCase) && !processes[i]->name.Contains("WJets",TString::kIgnoreCase)) {
            cout << "\tPlotFiller::Turning off the genParticlecollection* branch(es) ... ";
            c->SetBranchStatus("genParticleCollection*",0);
            cout << "DONE" << endl;
         }
      }
      else if(limitBranches > 0) {
         cout << "\tPlotFiller::Turning off the epd* and m_tProbStat* branch(es) ... ";
         c->SetBranchStatus("epd*",0);
         c->SetBranchStatus("m_tProbStat*",0);
         cout << "DONE" << endl;
         cout << "\tPlotFiller::Turning off the triggerMap, m_tSize, m_run, m_event, and indexMap branch(es) ... ";
         c->SetBranchStatus("triggerMap",0);
         c->SetBranchStatus("m_tSize",0);
         c->SetBranchStatus("m_run",0);
         c->SetBranchStatus("m_event",0);
         c->SetBranchStatus("indexMap",0);
         cout << "DONE" << endl;
         cout << "\tPlotFiller::Turning off the genParticlecollection* branch(es) ... ";
         c->SetBranchStatus("genParticleCollection*",0);
         cout << "DONE" << endl;

         if(processes[i]->name.Contains("TTbar",TString::kIgnoreCase)) {
            cout << "\tPlotFiller::Turning on the genParticlecollection.{p4,pdgId}* branch(es) ... ";
            c->SetBranchStatus("genParticleCollection.p4*",1);
            c->SetBranchStatus("genParticleCollection.pdgId*",1);
            cout << "DONE" << endl;
         }
         if(limitBranches == 1 && processes[i]->name.Contains("WJets",TString::kIgnoreCase)) {
            cout << "\tPlotFiller::Turning on the genParticlecollection.{p4,pdgId}* branch(es) ... ";
            c->SetBranchStatus("genParticleCollection.p4*",1);
            c->SetBranchStatus("genParticleCollection.pdgId*",1);
            cout << "DONE" << endl;  
         }
      }

      // Set the branch address
      if (c->GetBranch("EvtTree")) {
         ntuple = new EventNtuple();
         c->SetBranchAddress("EvtTree",&ntuple);
      }
      else if (c->GetBranch("EvtNtuple")) {
         ntuple = new EventNtuple();
         c->SetBranchAddress("EvtNtuple",&ntuple);
      }
      else {
	      cout << "\tPlotFiller::run EvtTree and EvtNtuple branches not found." << endl
	           << "\tSetting EventNtuple pointer to 0x0." << endl;
      }
      if (c->GetBranch("METree")) {
         metree = new METree();
         c->SetBranchAddress("METree",&metree);
      }
      else {
	      cout << "\tPlotFiller::run METree branch not found." << endl
              << "\tSetting METree pointer to 0x0." << endl;
      }
      if (c->GetBranch("mnt")) {
         mnt = new MicroNtuple(2);
         c->SetBranchAddress("mnt",&mnt);
         if (!MVAWeightDir.IsNull() && !MVAMethods.empty() && !MVAV.empty()) {
            c->GetEntry(16);

            MVAVar::setVarMap(MVAVars);
            for(unsigned int iMVA=0; iMVA<MVAV.size(); iMVA++) {
               MVAVars[MVAV[iMVA]].use = true;
               MVAVars[MVAV[iMVA]].index = iMVA;
               MVAVars[MVAV[iMVA]].maxIndex = MVAV.size();
            }
            for(unsigned int iMVA=0; iMVA<MVAS.size(); iMVA++) {
               MVAVars[MVAS[iMVA]].use = true;
               MVAVars[MVAS[iMVA]].index = iMVA;
               MVAVars[MVAS[iMVA]].maxIndex = MVAS.size();
            }

            mnt->setMVAReader(MVAMethods,MVAWeightDir,MVAVars);
         }
      }
      else {
         cout << "\tPlotFiller::run mnt branch not found." << endl
              << "\tSetting MicroNtuple pointer to 0x0." << endl;
      }

      // The counter for how many events pass the cuts in each process
      unsigned int numProcEvts = 0;
      double sumW = 0;
   
      // Define the weight variable and set the default as 1.0
      double weight = 1.0;

      // Loop over events in the process
      if(debug && debugNumberOfEvents < c->GetEntries())
      {
         numberOfEvents = debugNumberOfEvents;
      }
      else
      {
         numberOfEvents = c->GetEntries();
      }
      
      cout << "\tProcessing " << numberOfEvents << " events (out of " << c->GetEntries() << ")." << endl;
      
      // This runs once for each process before the events are run.
      userProcessFunc(ntuple, processes[i]);

      for (unsigned int ev = 0 ; ev < numberOfEvents ; ev++)
      {
         //if debug, time the event
         if(debug && numberOfEvents<100) {
            event_benchmark->Start("event");
            func_benchmark->Start("GetEntry");
         }

         // Report every now and then
         ProgressBar::loadbar2(ev+1,numberOfEvents);
         //if ((ev % 10000) == 0)
         //   cout<<"\t\tevent "<<ev<<endl;
         
         // Get the given entry
         c->GetEntry(ev);

         if(debug && numberOfEvents<100)
            func_benchmark->Stop("GetEntry");

         // Runs before any cuts are made
         userInitEventFunc(ntuple, processes[i]);
         
         // Cut events here
         if(!userCutFunc(ntuple, mnt, processes[i], MVAVars, MVAMethods))
            continue;
         
         // reset the weight for each event, so the default will be 1.0
         weight = 1.0;
         weight *= userWeightFunc(ntuple, mnt, processes[i]);
         
         // Fill plots
         //cout << "Entry = " << ev << endl; 
         userFillFunc(plots, processes[i]->name, ntuple, metree, mnt, MVAMethods, MVAVars, weight);
         
         // Keep track of the total number of events & weights passing the cuts 
         sumW += weight;
         numProcEvts++;

         if(debug && numberOfEvents<100) {
            event_benchmark->Stop("event");
            float rt = 0, ct = 0;
            cout << endl << "PlotFiller::event_benchmark" << endl;
            vector<string> timers(1,"event");
            DefaultValues::printSummary(event_benchmark, 8, rt, ct, timers);
            event_benchmark->Reset();

            cout << endl << "PlotFiller::func_benchmark" << endl;
            vector<string> timers2(1,"GetEntry");
            DefaultValues::printSummary(func_benchmark, 8, rt, ct, timers2);
            func_benchmark->Reset();
         }
      }
      
      cout << endl << "Number of " << processes[i]->name << " events passing the cuts: " << numProcEvts << " with weights: " << sumW << endl;
   }
   cout << "\nPlotFiller::run DONE" << endl;
}
