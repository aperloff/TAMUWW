//Our libraries
#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"

//ROOT Libraries
#include <TStyle.h>
#include "TDirectory.h"
#include "TH1F.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TMath.h"
#include "TGraph.h"

//C++ libraries
#include <map>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <utility>
#include <exception>

using namespace std;

//A struct describing a given systematic
//Contains a boolean indicating if it is a shape systematic (0=no, 1=yes), the position in the combined vector of processes (empty if for all processes), and the values of the systematic
class systematic{
public:
  systematic() : name(""), shape(false) {}
  systematic(TString name_, bool shape_, bool symmetrize_, vector<TString> processes);// : name(name_), shape(shape_), symmetrize(symmetrize_), values(map<DEFS::PhysicsProcess::Type,pair<TString,TString> > (nvalues_,make_pair("-","-"))) {}
  systematic(TString name_, bool shape_, bool symmetrize_, vector<TString> processes, TString value_);// : name(name_), shape(shape_), symmetrize(symmetrize_), values(map<DEFS::PhysicsProcess::Type,pair<TString,TString> > (nvalues_,make_pair(value_,value_))) {}
  systematic(TString name_, bool shape_, bool symmetrize_, vector<TString> processes, DEFS::PhysicsProcess::Type process, pair<TString,TString> values_);
  systematic(TString name_, bool shape_, bool symmetrize_, map<DEFS::PhysicsProcess::Type,pair<TString,TString> > values_) : name(name_), shape(shape_), symmetrize(symmetrize_), values(values_) {}

  string getValueString(DEFS::PhysicsProcess::Type index);
  void updateValue(DEFS::PhysicsProcess::Type index, pair<TString,TString> value_, bool verbose);

  TString         name;
  bool            shape;
  bool            symmetrize;
  //vector<pair<TString,TString> > values;
  map<DEFS::PhysicsProcess::Type,pair<TString,TString> > values;
};

systematic::systematic(TString name_, bool shape_, bool symmetrize_, vector<TString> processes) : name(name_), shape(shape_), symmetrize(symmetrize_) {
  for(unsigned int i=0; i<processes.size(); i++) {
    values[DEFS::PhysicsProcess::getProcessType(string(processes[i]))] = make_pair("-","-");
  }
}

systematic::systematic(TString name_, bool shape_, bool symmetrize_, vector<TString> processes, TString value_) : name(name_), shape(shape_), symmetrize(symmetrize_) {
  for(unsigned int i=0; i<processes.size(); i++) {
    values[DEFS::PhysicsProcess::getProcessType(string(processes[i]))] = make_pair(value_,value_);
  }
}

systematic::systematic(TString name_, bool shape_, bool symmetrize_, vector<TString> processes, DEFS::PhysicsProcess::Type process, pair<TString,TString> values_) : name(name_), shape(shape_), symmetrize(symmetrize_) {
  for(unsigned int i=0; i<processes.size(); i++) {
    if(process == DEFS::PhysicsProcess::getProcessType(string(processes[i]))) {
      values[DEFS::PhysicsProcess::getProcessType(string(processes[i]))] = make_pair(values_.first,values_.second);
    }
    else {
      values[DEFS::PhysicsProcess::getProcessType(string(processes[i]))] = make_pair("-","-");
    }
  }
}

string systematic::getValueString(DEFS::PhysicsProcess::Type index) {
  std::stringstream ss;

  if(symmetrize && values[index].first != "-" && values[index].second != "-") {
    double diff1 = TMath::Abs(1.0-values[index].first.Atof());
    double diff2 = TMath::Abs(1.0-values[index].second.Atof());

    double bigger = max(diff1,diff2)==diff1 ? 1.0+diff1 : 1.0+diff2;//values[index].first.Atof() : values[index].second.Atof();

    values[index].first = Form("%.3f",bigger);
    values[index].second = Form("%.3f",bigger);
  }

  if(values[index].first == values[index].second) {
    ss << values[index].first;
  }
  else {
    ss << values[index].first << "/" << values[index].second;
  }
  return ss.str();
}

void systematic::updateValue(DEFS::PhysicsProcess::Type index, pair<TString,TString> value_, bool verbose) {
  if(values.count(index))
    values.at(index) = value_;
  else if(verbose)
    cout << "WARNING::systematic::updateValue There is no process " << DEFS::PhysicsProcess::getTypeString(index)
         << " in the values map for systematic with name " << name << endl;
}

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
int main(int argc,char**argv) {
 
  // evaluate command-line / configuration file options
  CommandLine cl;
  if (!cl.parse(argc,argv)) return 0;

  TString         outpath           = cl.getValue<TString> ("outpath");
  TString         basepath          = cl.getValue<TString> ("basepath",          "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_04_27_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/");
  TString         histogramFileBase = cl.getValue<TString> ("histogramFileBase", "histos");
  TString         suffix            = cl.getValue<TString> ("suffix",                  "");
  TString         variable          = cl.getValue<TString> ("variable",           "Mlvjj");
  TString         channelName       = cl.getValue<TString> ("channelName",       "TwoJ0B");
  TString         lepCat            = cl.getValue<TString> ("lepton",              "both");
  DEFS::LeptonCat leptonCat         = DEFS::getLeptonCat(string(lepCat));
  string          jetBinS           = cl.getValue<string>  ("jetBin",             "jets2");
  DEFS::JetBin    jetBin            = DEFS::getJetBin(jetBinS);
  string          tagcatS           = cl.getValue<string>  ("tagcat",            "pretag");
  DEFS::TagCat    tagcat            = DEFS::getTagCat(tagcatS);
  string          anaTypeS          = cl.getValue<string>  ("anaType",    "HiggsAnalysis");
  DEFS::Ana::Type anaType           = DEFS::Ana::getAnaType(anaTypeS);
  vector<TString> overrideProcesses = cl.getVector<TString>("overrideProcesses",       "");
  bool            batch             = cl.getValue<bool>    ("batch",                 true);
  bool            debug             = cl.getValue<bool>    ("debug",                false);
  bool            verbose           = cl.getValue<bool>    ("verbose",              false);

  if (!cl.check()) return 0;
  cl.print();

  //TString histogramFile = histogramFileBase+"_"+lepCat+"_withSumMC_SysNames.root";
  TString histogramFile = histogramFileBase+"_"+lepCat+"_"+jetBinS+"_SysNames"+suffix+".root";
  TFile* in_file = TFile::Open(basepath+histogramFile,"READ"); // Input the root file filled with histograms

  // Get the input processes
  vector<PhysicsProcess*> processes = DefaultValues::getProcessesHiggs(jetBin, tagcat, false, false, false, DEFS::MicroNtuple);
  if(debug) {
     processes.erase(processes.begin(),processes.begin()+11);
  }

  // Name the bin for the processes
  vector <TString> vecSig;
  vector <TString> vecBak;
  vector <TString> vecCombined;

  for(unsigned int iprocess=0; iprocess<processes.size(); iprocess++) {

    TString processName = processes[iprocess]->name;

    //Specifically list the processes to use
    if(overrideProcesses.size()>0) {
      bool skip = true;
      for(unsigned int ip=0; ip<overrideProcesses.size(); ip++) {
        if(processName.CompareTo(overrideProcesses[ip])==0) {
          skip=false;
          break;
        }
      }
      if(skip)
        continue;
    }

    //Skip Data
    if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::SingleEl_Data) {
       continue;
    }
    else if(DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::SingleMu_Data) {
       continue;
    }

    if(leptonCat==DEFS::electron && DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::QCD_MuFULL)
      continue;
    if(leptonCat==DEFS::muon && DEFS::PhysicsProcess::getProcessType(string(processName))==DEFS::PhysicsProcess::QCD_ElFULL)
      continue;


    if ((anaType == DEFS::Ana::WWAnalysis && (processName.EqualTo("WW") || processName.EqualTo("WZ"))) || 
        (anaType == DEFS::Ana::HiggsAnalysis && (processName.Contains("WH") || processName.Contains("ggH") || processName.Contains("qqH") || processName.Contains("ZH") || processName.Contains("TTH"))) ||
        (anaType == DEFS::Ana::UNKNOWN && (processName.Contains("WH") || processName.Contains("ggH") || processName.Contains("qqH")))){

      // Input each signal
      //vecSig.push_back(variable+"_"+processName+"_"+lepCat);
      vecSig.push_back(processName);
      
      cout<<"\tputting process "<<processName<<" into signal."<<endl;
    
    }
    else{ 
      // Input each background
      if(leptonCat==DEFS::both && processName.Contains("QCD_El")) {
         //vecBak.push_back(variable+"_QCD_"+lepCat);
        vecBak.push_back("QCD");
         cout<<"\tputting process QCD into background."<<endl;  
      }
      else if(leptonCat==DEFS::both && processName.Contains("QCD_Mu")) {
         continue;
      }
      else {
         //vecBak.push_back(variable+"_"+processName+"_"+lepCat);
        vecBak.push_back(processName);
         cout<<"\tputting process "<<processName<<" into background."<<endl;  
      }
    }
  }

  //Combine the siganl and background vectors
  vecCombined = vecSig;
  vector<TString>::iterator it = vecCombined.end();
  vecCombined.insert(it,vecBak.begin(),vecBak.end());

  // Creating the datacard
  ofstream outFile(Form("%s/DataCard_%s_%s_%s.txt",outpath.Data(),lepCat.Data(),channelName.Data(),variable.Data()));
  outFile << "#Simple datacard for ggH M125 Low mass H->WW->lnujj experiment" << endl
          << "#can list imax,jmax, kmax as specific values or * to make it auto calculate" << endl
          << "#expected limit command (with systematics) \"combine -M Asymptotic --significance *DATACARDNAME* -t -1 --expectSignal=1 -S 1 -m 125 --run expected\"" << endl
          << "#expected limit command (without systematics) \"combine -M Asymptotic --significance *DATACARDNAME* -t -1 --expectSignal=1 -S 0 -m 125 --run expected\"" << endl
          << "#observed limit command (with systematics) \"combine -M Asymptotic --significance *DATACARDNAME* -t -1 --expectSignal=1 -m 125\"" << endl
          << "#observed p-value/significance command \"combine -M ProfileLikelihood --significance *DATACARDNAME* --expectSignal=1 -S 1 -m 125\"" << endl
          << "#expected p-value/significance command \"combine -M ProfileLikelihood --significance *DATACARDNAME* -t -1 --expectSignal=1 -S 0 -m 125\"" << endl
          << "#To combine the datacards in each channel use the command \"combineCards.py Name1=card1.txt Name2=card2.txt .... > card.txt\"" << endl
          << "#MaxLikelidooh fit of datacard to get mlfit.root file \"combine -M MaxLikelihoodFit -t -1 --expectSignal=1 --minos=all DataCardName.txt\"" << endl
          << "#To get the post-fit nuisance parameters \"python diffNuisances.py -a --format=latex -g plots.root mlfit.root\"" << endl
          << "#Can use \"--setPhysicsModelParameters r=0.1 --setPhysicsModelParameterRanges r=0,5\" if you need to" << endl;

  // process counter line
  //int jmax = vecSig.size()+vecBak.size();

  // Input total number of bins
  outFile << "imax " << left << setw(2) << "1" << "  number of channels" << endl;
  
  // Input total number of processes in bin1 minus one
  //outFile << "jmax " << left << setw(2) << jmax - 1 << "  number of backgrounds minus 1" << endl;
  outFile << "jmax " << left << setw(2) << "*" << "  number of backgrounds minus 1 or *" << endl;
  
  // Input total number of systematics
  outFile << "kmax " << left << setw(2) << "*" << "  number of nuisance parameters (sources of systematic uncertainties)" << endl << endl;
  
  outFile << "#list shape uncertainties and the input files for them " << endl
          << "#NOTE: hists must be in root file with names exactly as listed in systematics" << endl
          << "#The channel name can be either the actual name or a * for all channels" << endl
          << "----------------------------------------------------------------------------------------------"  << endl << endl;

  // Input the root file name and whatever edits to the names you wish to make
  //outFile << "shapes " << left << setw(9) << "* "<< channelName << "\t" << basepath << "histos_" << lepCat << "_SysNames.root " << variable << "_$PROCESS_" << lepCat << " " << variable << "_$PROCESS_" << lepCat << "_$SYSTEMATIC" << endl;
  basepath = batch ? "" : basepath;
  outFile << "shapes " << left << setw(9) << "* "<< channelName << "\t" << basepath << histogramFile << " " << variable << "_$PROCESS_" << lepCat << " " << variable << "_$PROCESS_" << lepCat << "_$SYSTEMATIC" << endl;

  // Input the root file name and the data_obs file name
  TH1D* dataHist = 0;
  if(leptonCat==DEFS::electron) {
     outFile << "shapes " << left << setw(9) << "data_obs " << channelName << "\t" << basepath << histogramFile << " " << variable << "_SingleEl_Data_" << lepCat << endl << endl;
    dataHist = (TH1D*) in_file->Get(variable+"_SingleEl_Data_"+lepCat);
  }
  else if(leptonCat==DEFS::muon) {
     outFile << "shapes " << left << setw(9) << "data_obs " << channelName << "\t" << basepath << histogramFile << " " << variable << "_SingleMu_Data_" << lepCat << endl << endl;
    dataHist = (TH1D*) in_file->Get(variable+"_SingleMu_Data_"+lepCat);
  }
  else if(leptonCat==DEFS::both) {
     outFile << "shapes " << left << setw(9) << "data_obs " << channelName << "\t" << basepath << histogramFile << " " << variable << "_Data_" << lepCat << endl << endl;
    dataHist = (TH1D*) in_file->Get(variable+"_Data_"+lepCat);
  }

  outFile << "# bin = name of channel we are observing events in, i.e. 2+jets" << endl
          << "# observation = number of events seen in that bin" << endl << endl;  
  
  // Input the name of your bin
  outFile << left << setw(12) << "bin " << channelName << endl; 
  
  outFile << left << setw(12) << "observation " << std::setprecision(8) << dataHist->Integral() << endl << endl;
  
  outFile << "# bin = bin from above.  Each signal and background will have a column for every bin listed" << endl
          << "# process line 1: the name of the process.  This must match the syntax used in root files of the shape uncertainties" << endl
          << "# process line 2: assign a number to the process.  All signals get process = 0, all backgrounds get independent positive numbers" << endl
          << "# rate = expected number of events of the process in specified bin" << endl
          << "----------------------------------------------------------------------------------------------" << endl << endl;

  outFile << left << setw(8) << "bin "; 

  //Bin line
  for (unsigned int num_bins = 0 ; num_bins < (vecSig.size()+vecBak.size()); num_bins++){
    outFile << left << setw(max(11,max(vecCombined[num_bins].Length()+1,channelName.Length()+1))) << channelName;
  } //for bin line
  
  outFile << endl;
  
  outFile << left << setw(8) << "process ";

  // signal process name line
  for (unsigned int name_sig = 0 ; name_sig < vecSig.size(); name_sig++){
    outFile << left << setw(max(11,max(vecSig[name_sig].Length()+1,channelName.Length()+1))) << vecSig[name_sig];
  } // for sig process name line

  // background process name line
  for (unsigned int name_bak = 0 ; name_bak < vecBak.size(); name_bak++){
    outFile << left << setw(max(11,max(vecBak[name_bak].Length()+1,channelName.Length()+1))) << vecBak[name_bak];
  } // for bak process name line
  
  outFile << endl;
  
  outFile << left << setw(8) << "process ";

  // Making an array to output the signal process IDs
  for (int sig_num = vecSig.size()-1 ; sig_num >= 0; sig_num--){ 

    // Output a negative number (or zero) for each signal process
    //if(sig_num == int(vecSig.size()-1))
    //  outFile << left << setw(max(10,max(vecSig[sig_num].Length()+2,channelName.Length()+2))) << Form("%i",-int(sig_num));
    //if (sig_num == 0)
    //  outFile << left << setw(max(11,max(vecSig[sig_num].Length(),channelName.Length()))) << Form("%i",-int(sig_num));
    //else
      outFile << left << setw(max(11,max(vecSig[vecSig.size()-1-sig_num].Length()+1,channelName.Length()+1))) << Form("%i",-int(sig_num));

  } // for sig num ID

  // Making an array to output the background process IDs
  for (unsigned int bak_num = 0; bak_num < vecBak.size(); bak_num++){
    // Output a positive number for each background process
    outFile << left << setw(max(11,max(vecBak[bak_num].Length()+1,channelName.Length()+1))) << bak_num+1;
  } // for bak num ID

  outFile << endl;
  
  outFile << left << setw(8) << "rate ";

   //Printing process rates
  for (unsigned int iproc = 0 ; iproc < vecCombined.size(); iproc++){
    
    // Get the histograms out
    TH1D* proc = (TH1D*) in_file->Get(variable+"_"+vecCombined[iproc]+"_"+lepCat);
    
    // Integrate each histogram
    double data_store = proc->Integral(); 
    
    // Print signal integral data
    outFile << left << setw(max(11,max(vecCombined[iproc].Length()+1,channelName.Length()+1))) << std::setprecision(8) << data_store; 
    
  } // for signal integral loop

  outFile << endl << endl;
  
  outFile << "------------------------" << endl
          << "#List all systematic uncertainties" << endl
          << "#column 1 = sysName, column 2 = Dist type (options are lnN, shape), column 3+ = uncertainty value. 1.2 = 20% uncertainty. - means N/A " << endl << endl;


  // adding systematics
  //WH_HToZZ_M125 ZH_HToZZ_M125 TTH_HToZZ_M125 WH125_HToBB TTH_HToBB_M125 ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 WZ STopS_T STopS_Tbar STopT_T STopT_Tbar STopTW_T STopTW_Tbar TTbar WW QCD_ElFULL ZJets WJets     
  //0             1             2              3           4              5      6      7             8             9              10 11      12         13      14         15       16          17    18 19         20    21        
  bool make_symmetric = false;
  bool separate_shapes_and_rates = false;
  DEFS::PhysicsProcess::Type QCD_Enum;
  if(vecCombined.size()>=20)
    QCD_Enum = DEFS::PhysicsProcess::getProcessType(string(vecCombined[19]));
  else if(leptonCat==DEFS::electron) {
    QCD_Enum = DEFS::PhysicsProcess::getProcessType(string("QCD_ElFULL"));
  }
  else if(leptonCat==DEFS::muon) {
    QCD_Enum = DEFS::PhysicsProcess::getProcessType(string("QCD_MuFULL"));
  }
  else
    QCD_Enum = DEFS::PhysicsProcess::Type::UNKNOWN;
  vector<systematic>systematics;

  //Flat ratee uncernties across the board
  systematics.push_back(systematic("lumi_8TeV",false,false,vecCombined,"1.026"));
  systematics.back().updateValue(QCD_Enum, make_pair("-","-"), verbose); //QCD needs rate sys removed
  //systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("-","-"), verbose); //WJets needs rate sys removed

  if(leptonCat==DEFS::electron) {
    systematics.push_back(systematic("CMS_eff_e",false,false,vecCombined,"1.02"));
    systematics.back().updateValue(QCD_Enum, make_pair("-","-"), verbose); //QCD needs rate sys removed
    //systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("-","-"), verbose); //WJets needs rate sys removed
    systematics.push_back(systematic("CMS_eff_m",false,false,vecCombined));
  }
  else if(leptonCat==DEFS::muon) {
    systematics.push_back(systematic("CMS_eff_m",false,false,vecCombined,"1.02"));
    systematics.back().updateValue(QCD_Enum, make_pair("-","-"), verbose); //QCD needs rate sys removed
    //systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("-","-"), verbose); //WJets needs rate sys removed
    systematics.push_back(systematic("CMS_eff_e",false,false,vecCombined));
  }
  else if(leptonCat==DEFS::both) {
    systematics.push_back(systematic("CMS_eff_l",false,false,vecCombined,"1.02"));
    systematics.back().updateValue(QCD_Enum, make_pair("-","-"), verbose); //QCD needs rate sys removed
    //systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("-","-"), verbose); //WJets needs rate sys removed
  }

  systematics.push_back(systematic("CMS_scale_met",false,false,vecCombined,"1.002"));
  systematics.back().updateValue(QCD_Enum, make_pair("-","-"), verbose); //QCD needs rate sys removed
  //systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("-","-"), verbose); //WJets needs rate sys removed

  //rate uncertainties specific to a given process
  systematics.push_back(systematic("QCDscale_ggH",false,false,vecCombined)); //ggH scale uncertainty
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("0.922","1.072"), verbose); //ggH125
  systematics.push_back(systematic("QCDscale_qqH",false,false,vecCombined)); //qqH scale uncertainty
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.998","1.002"), verbose); //qqH125
  systematics.push_back(systematic("QCDscale_WH",false,false,vecCombined)); //WH scale uncertainty
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("0.99","1.01"), verbose); //WH125_HToBB
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.99","1.01"), verbose); //WH_HToWW_M125
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.99","1.01"), verbose); //WH_HToZZ_M125
  systematics.push_back(systematic("QCDscale_ZH",false,false,vecCombined)); //ZH scale uncertainty
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.969","1.031"), verbose); //ZH_HToWW_M125
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("0.969","1.031"), verbose); //ZH_HToZZ_M125
  systematics.push_back(systematic("QCDscale_TTH",false,false,vecCombined)); //TTH scale uncertainty
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("0.907","1.038"), verbose); //TTH_HToWW_M125
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("0.907","1.038"), verbose); //TTH_HToZZ_M125
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("0.907","1.038"), verbose); //TTH_HToBB_M125
  systematics.push_back(systematic("pdf_gg",false,false,vecCombined)); //gg->X PDF+alpha_s uncertainty
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("0.931","1.075"), verbose); //ggH125
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("0.919","1.081"), verbose); //TTH_HToWW_M125
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("0.919","1.081"), verbose); //TTH_HToZZ_M125
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("0.919","1.081"), verbose); //TTH_HToBB_M125
  systematics.push_back(systematic("pdf_qq",false,false,vecCombined)); //qq->X PDF+alpha_s uncertainty
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.972","1.026"), verbose); //qqH125
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("0.977","1.023"), verbose); //WH125_HToBB
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.977","1.023"), verbose); //WH_HToWW_M125
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.975","1.025"), verbose); //ZH_HToWW_M125
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.977","1.023"), verbose); //WH_HToWZZ_M125
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("0.975","1.025"), verbose); //ZH_HToWZZ_M125
  //systematics.push_back(systematic("CMS_QCDscale_pdf_wjets",false, vecCombined)); //Do not include when using QCD Eta Weight Uncertainties
  //systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets")) = "1.033";
  systematics.push_back(systematic("CMS_QCDscale_pdf_ttbar",false,false,vecCombined));
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("1.057","1.057"), verbose); //TTbar https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO#Top_quark_pair_cross_sections_at
  systematics.push_back(systematic("CMS_QCDscale_pdf_zjets",false,false,vecCombined));
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.034","1.034"), verbose); //ZJets
  systematics.push_back(systematic("CMS_QCDscale_pdf_singlet",false,false,vecCombined));
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.05","1.05"), verbose);
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.05","1.05"), verbose);
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.05","1.05"), verbose);
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.05","1.05"), verbose);
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("1.05","1.05"), verbose);
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("1.05","1.05"), verbose);
  systematics.push_back(systematic("CMS_QCDscale_pdf_diboson",false,false,vecCombined));
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("1.03","1.03"), verbose); //WZ
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("1.03","1.03"), verbose); //WW
  //systematics.push_back(systematic("CMS_hww_lnujj_qcdNorm",false,vecCombined)); //Do not include when using QCD Eta Weight Uncertainties
  //systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1.02","1.02"), verbose);
  systematics.push_back(systematic("CMS_hww_lnujj_PUWeight",false,make_symmetric,vecCombined)); 
  //systematics.push_back(systematic(TString("CMS_hww_lnujj_PUWeight")+TString(separate_shapes_and_rates ? "" : "_shape"),
  //                                 !separate_shapes_and_rates,separate_shapes_and_rates ? make_symmetric : true,vecCombined)); //Vary min bias xsec up and down by 7% (nominal 69400 up 74258 down 64542)
  if(jetBin==DEFS::jets2) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.984","1.019"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("0.959","1.035"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("1.012","1.014"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("0.995","1.014"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("0.981","1.049"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("0.995","1.011"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.984","1.024"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.983","1.015"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.980","1.027"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("1.094","0.897"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("0.992","1.014"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("0.950","1.059"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.006","1.010"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("0.994","1.011"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("0.992","1.017"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("0.993","1.009"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("0.998","1.009"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.982","1.023"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("0.990","1.015"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.039","0.971"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.951","0.963"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.993","1.011"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("0.997","1.020"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("1.007","1.023"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("0.987","1.016"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("0.987","0.987"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("0.997","1.010"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.993","1.012"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.978","1.026"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.959","1.042"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("0.943","1.078"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("0.990","1.016"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("0.984","1.011"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.021","0.990"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("0.991","1.015"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("0.989","1.014"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("0.992","1.014"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("0.991","1.015"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.989","1.017"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("0.990","1.015"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.045","0.966"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.999","1.005"), verbose);
    }
  }
  else if(jetBin==DEFS::jets3) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.982","1.024"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("1.065","0.931"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("1.005","1.020"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("0.982","1.026"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("0.934","1.051"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("0.993","1.013"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.998","1.010"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.983","1.022"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.976","1.031"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("0.897","1.139"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("1.004","1.001"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.037","0.981"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("0.978","1.019"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("0.993","1.009"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("0.994","1.005"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("0.988","1.017"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("0.968","1.035"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.989","1.016"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("0.998","1.007"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.036","0.975"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.941","0.954"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.984","1.020"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("0.974","1.038"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("0.926","1.078"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("0.992","1.016"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("1.000","1.007"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("1.001","1.005"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.993","1.011"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("1.004","1.003"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.965","1.026"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("0.890","1.107"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("0.999","1.006"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("0.993","1.013"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.004","1.011"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.000","1.011"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.009","1.001"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("0.996","1.008"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("0.993","1.012"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.991","1.014"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("1.000","1.006"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.040","0.970"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.998","1.005"), verbose);
    }
  }
  else if(jetBin==DEFS::jets4) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("1.024","0.991"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("1.008","1.005"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("0.976","1.010"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.021","0.983"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("1.003","1.003"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("1.010","0.997"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.985","1.019"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("1.001","1.011"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.948","1.071"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("0.934","1.084"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("1.003","1.005"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.211","0.857"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.128","0.927"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("0.995","1.005"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("0.995","0.990"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("0.985","1.017"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("1.000","1.001"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.991","1.014"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("0.993","1.013"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.050","0.962"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.928","0.920"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.988","1.025"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("1.037","0.980"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("0.973","1.028"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("0.978","1.024"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("0.999","1.010"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("0.993","1.011"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("1.006","1.005"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("1.021","0.994"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("1.042","0.988"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("0.947","1.058"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("1.004","1.003"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("0.953","1.069"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.089","0.960"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.019","0.991"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("0.995","1.011"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("0.991","1.011"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("0.987","1.018"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.991","1.016"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("1.000","1.006"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.057","0.959"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.995","1.017"), verbose);
    }
  }
  systematics.push_back(systematic("CMS_hww_lnujj_CSVWeight",false,make_symmetric,vecCombined));
  //systematics.push_back(systematic(TString("CMS_hww_lnujj_CSVWeight")+TString(separate_shapes_and_rates ? "" : "_shape"),
  //                                 !separate_shapes_and_rates,separate_shapes_and_rates ? make_symmetric : true,vecCombined));
  if(jetBin==DEFS::jets2) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.959","1.011"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("1.002","0.969"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("1.090","0.886"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.135","0.816"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("1.061","0.899"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("0.955","1.016"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.977","0.996"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.977","0.997"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.997","0.974"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("1.082","0.920"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("0.981","0.989"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.122","0.829"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.136","0.832"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.054","0.902"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.047","0.912"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("1.034","0.930"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("1.024","0.933"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("1.058","0.902"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("0.975","0.998"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.038","0.933"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.920","0.959"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.991","0.977"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("0.985","0.988"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("1.226","0.652"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.141","0.812"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("1.098","0.847"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("0.956","1.015"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.979","0.994"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.984","0.988"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.965","1.007"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("1.111","0.892"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("0.978","0.993"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.098","0.847"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.091","0.862"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.054","0.902"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.048","0.909"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("1.021","0.943"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("1.025","0.937"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("1.061","0.898"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("0.977","0.996"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.016","0.957"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.965","1.001"), verbose);
    }
  }
  else if(jetBin==DEFS::jets3) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.961","0.992"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("0.987","0.966"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("1.008","0.941"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.134","0.801"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("1.189","0.746"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("0.935","1.021"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.971","0.989"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.975","0.981"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.982","0.985"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("1.154","0.813"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("0.970","0.983"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.119","0.748"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.084","0.862"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.064","0.865"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.067","0.871"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("1.018","0.930"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("1.023","0.914"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("1.072","0.868"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("0.957","1.000"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.025","0.926"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.893","0.954"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.959","1.002"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("0.960","0.997"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("1.121","0.854"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.158","0.789"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("1.169","0.768"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("0.939","1.017"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.972","0.989"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.956","0.998"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.958","1.006"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("1.103","0.889"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("0.965","0.987"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.167","0.814"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.124","0.818"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.052","0.873"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.062","0.880"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("1.025","0.926"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("1.031","0.917"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("1.069","0.871"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("0.961","0.996"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.002","0.954"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.941","1.010"), verbose);
    }
  }
  else if(jetBin==DEFS::jets4) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.976","0.949"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("1.027","0.924"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("1.007","0.893"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.132","0.795"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("1.250","0.698"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("0.911","1.023"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.964","0.972"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.953","0.975"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.931","1.017"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("1.206","0.735"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("0.947","0.987"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.130","0.788"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("0.981","0.930"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.040","0.871"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.071","0.866"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("1.013","0.905"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("1.049","0.873"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("1.060","0.857"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("0.942","0.996"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("0.997","0.937"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.844","0.933"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.963","0.937"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("1.111","0.827"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("1.141","0.749"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.126","0.801"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("1.201","0.716"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("0.910","1.024"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("0.964","0.979"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.973","0.960"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.959","0.977"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("1.063","0.859"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("0.961","0.972"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.123","0.809"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.277","0.725"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.115","0.820"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.084","0.837"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("1.033","0.889"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("1.080","0.855"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("1.066","0.847"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("0.933","1.003"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("0.975","0.960"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.933","0.996"), verbose);
    }
  }
  systematics.push_back(systematic(TString("CMS_hww_lnujj_topPtWeight")+TString(separate_shapes_and_rates ? "" : "_shape"),
                                   !separate_shapes_and_rates,separate_shapes_and_rates ? make_symmetric : true,vecCombined)); //up as TopPtWeight^2 and down as no topPTWeight
  if(jetBin==DEFS::jets2) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.998","0.994"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.998","0.993"), verbose);
    }
  }
  else if(jetBin==DEFS::jets3) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.984","1.007"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.988","1.004"), verbose);
    }
  }
  else if(jetBin==DEFS::jets4) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.962","1.028"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.970","1.021"), verbose);
    }
  }
  systematics.push_back(systematic(TString("CMS_scale_j")+TString(separate_shapes_and_rates ? "" : "_shape"),
                                   !separate_shapes_and_rates,separate_shapes_and_rates ? make_symmetric : true,vecCombined));
  if(jetBin==DEFS::jets2) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("1.008","0.995"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("0.975","0.932"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("1.162","1.000"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.056","1.023"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("4.468","1.041"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("1.025","0.972"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("1.001","0.996"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.982","0.990"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("1.005","0.993"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("1.000","1.000"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("1.028","0.995"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.025","1.020"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.081","0.925"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.020","0.999"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.016","0.985"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("0.983","1.029"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("0.958","1.023"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.849","0.929"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("0.997","0.951"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.041","0.971"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1.002","0.960"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("0.997","0.995"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("0.969","0.979"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("0.940","0.994"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.058","1.045"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("6.666","1.002"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("1.029","0.973"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("1.008","1.004"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("0.989","1.003"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("0.997","0.999"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("1.000","0.888"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("1.023","0.987"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.018","1.068"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("0.970","1.006"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.008","0.995"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("0.994","0.994"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("1.005","1.027"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("0.981","1.020"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.848","0.936"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("0.997","0.951"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.034","0.967"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.938","0.895"), verbose);
    }
  }
  else if(jetBin==DEFS::jets3) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("1.080","0.947"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("1.042","1.049"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("1.011","1.060"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.098","1.042"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("4.517","1.112"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("1.056","0.954"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("1.044","0.972"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("1.036","0.980"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("1.031","0.991"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("0.742","1.131"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("1.062","0.967"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("0.989","0.992"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.043","1.014"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.021","0.935"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.051","0.963"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("1.048","1.000"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("1.012","0.982"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.947","0.986"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("1.031","0.939"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.035","0.938"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1.055","0.959"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("1.034","0.993"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("1.103","0.992"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("0.898","1.100"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.101","1.010"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("4.865","1.053"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("1.043","0.945"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("1.051","0.973"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("1.023","0.934"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("1.040","0.912"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("0.901","1.000"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("1.052","0.966"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.083","0.940"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.079","0.898"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.075","0.982"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.041","0.944"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("1.001","1.002"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("1.004","1.007"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("0.952","0.977"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("1.028","0.927"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.068","0.942"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.977","0.881"), verbose);
    }
  }
  else if(jetBin==DEFS::jets4) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("1.041","0.923"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("1.071","0.965"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("0.996","1.061"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.130","0.991"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("5.918","1.004"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("1.080","0.920"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("1.089","0.944"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("1.059","0.938"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("1.104","0.921"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("1.022","0.966"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("1.084","0.947"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.276","0.704"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("1.285","0.988"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.096","0.960"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.033","0.913"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("1.045","0.957"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("1.060","0.971"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("1.058","0.997"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("1.073","0.916"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.103","0.905"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1.055","0.931"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToZZ_M125"), make_pair("1.068","0.920"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToZZ_M125"), make_pair("1.008","0.998"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToZZ_M125"), make_pair("1.067","1.012"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH125_HToBB"), make_pair("1.206","1.001"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToBB_M125"), make_pair("4.355","1.007"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ggH125"), make_pair("1.089","0.919"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("qqH125"), make_pair("1.074","0.931"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WH_HToWW_M125"), make_pair("1.094","0.921"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZH_HToWW_M125"), make_pair("1.047","0.982"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTH_HToWW_M125"), make_pair("0.972","0.985"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WZ"), make_pair("1.076","0.933"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_T"), make_pair("1.061","0.949"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopS_Tbar"), make_pair("0.999","0.852"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_T"), make_pair("1.086","0.937"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopT_Tbar"), make_pair("1.167","0.922"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_T"), make_pair("1.060","0.962"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("STopTW_Tbar"), make_pair("1.063","0.942"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("1.055","0.998"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WW"), make_pair("1.061","0.903"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("ZJets"), make_pair("1.074","0.935"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.985","0.850"), verbose);
    }
  }
  systematics.push_back(systematic(TString("CMS_hww_lnujj_QCDEtaWeight")+TString(separate_shapes_and_rates ? "" : "_shape"),
                                   !separate_shapes_and_rates,separate_shapes_and_rates ? make_symmetric : true,vecCombined));
  //ratio of scale factors from fit
  if(leptonCat==DEFS::electron) {
    systematics.back().updateValue(QCD_Enum, make_pair("2.0431651785","1.1069364859"), verbose);
    systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.9992536956","1.0002392001"), verbose);
  }
  else if(leptonCat==DEFS::muon) {
    systematics.back().updateValue(QCD_Enum, make_pair("1.0990661404","1.0500213179"), verbose);
    systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("0.9999783397","1.0000103144"), verbose);
  }
  //ratio of histograms after fit
  /*
  if(jetBin==DEFS::jets2) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(QCD_Enum, make_pair("0.994","1.002"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1.001","1.000"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(QCD_Enum, make_pair("1.000","1.000"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1.000","1.000"), verbose);
    }
  }
  else if(jetBin==DEFS::jets3) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(QCD_Enum, make_pair("0.985","1.005"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1.001","1.000"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(QCD_Enum, make_pair("1.000","1.000"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1.000","1.000"), verbose);
    }
  }
  else if(jetBin==DEFS::jets4) {
    if(leptonCat==DEFS::electron) {
      systematics.back().updateValue(QCD_Enum, make_pair("0.970","1.010"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1.001","1.000"), verbose);
    }
    else if(leptonCat==DEFS::muon) {
      systematics.back().updateValue(QCD_Enum, make_pair("1.000","1.000"), verbose);
      systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1.000","1.000"), verbose);
    }
  }
  */

  //shape uncertainties
  systematics.push_back(systematic("CMS_hww_lnujj_matching_shape",true,true,vecCombined));
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1","1"), verbose);
  systematics.push_back(systematic("CMS_hww_lnujj_scale_shape",true,true,vecCombined));
  systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1","1"), verbose);
  if(separate_shapes_and_rates) {
     systematics.push_back(systematic("CMS_hww_lnujj_topPtWeight_shape",true,true,vecCombined));
     systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("TTbar"), make_pair("1","1"), verbose); 
     systematics.push_back(systematic("CMS_hww_lnujj_QCDEtaWeight_shape",true,true,vecCombined));
     systematics.back().updateValue(QCD_Enum, make_pair("1","1"), verbose);
     //systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("1","1"), verbose); //Only need QCD shape change, but still need WJets rate change
     systematics.push_back(systematic("CMS_scale_j_shape",true,true,vecCombined,"1"));
     systematics.back().updateValue(QCD_Enum, make_pair("-","-"), verbose); //QCD does not have a jes uncertainty
     //systematics.back().updateValue(DEFS::PhysicsProcess::getProcessType("WJets"), make_pair("-","-"), verbose); //Do not include WJets when using QCD Eta Weight Uncertainties
     systematics.push_back(systematic("CMS_hww_lnujj_PUWeight_shape",true,true,vecCombined,"1"));
     systematics.back().updateValue(QCD_Enum, make_pair("-","-"), verbose); //QCD does not have a PU uncertainty
     systematics.push_back(systematic("CMS_hww_lnujj_CSVWeight_shape",true,true,vecCombined,"1"));
     systematics.back().updateValue(QCD_Enum, make_pair("-","-"), verbose); //QCD does not have a CSV uncertainty
  }
  //statistical errors
  for(unsigned int iproc=0; iproc<vecCombined.size(); iproc++) {
    TString var_proc_lep = variable+"_"+vecCombined[iproc]+"_"+lepCat;
    TH1D* proc = (TH1D*) in_file->Get(var_proc_lep);
    int nbins = proc->GetNbinsX();
    for(int ibin=1; ibin<=nbins; ibin++) {
      TString sysName = "CMS_hww_lnujj_statError_"+channelName+"_"+var_proc_lep+"_bin_"+Form("%i",ibin);
      if(!in_file->Get(var_proc_lep+"_"+sysName+"Up") || !in_file->Get(var_proc_lep+"_"+sysName+"Up")) continue;
      //if(vecCombined[iproc].CompareTo("WJets")!=0) continue;
      systematics.push_back(systematic(sysName,true,true,vecCombined,DEFS::PhysicsProcess::getProcessType(string(vecCombined[iproc])), make_pair("1","1")));
    }
  }

  //Find the maximum width of the systematics' names
  int sysWidth = 0;
  for (unsigned int iSys=0; iSys!=systematics.size(); iSys++) {
    if(systematics[iSys].name.Length()>sysWidth) sysWidth = systematics[iSys].name.Length();
  }

  cout << "Printing the systematics ... ";
  for (unsigned int iSys=0; iSys!=systematics.size(); iSys++) {
    outFile << left << setw(sysWidth+1) << systematics[iSys].name;
    if(systematics[iSys].shape==false) outFile << left << setw(7) << "lnN";
    else                               outFile << left << setw(7) << "shape";

    for (unsigned int iproc = 0 ; iproc < (systematics[iSys].values.size()); iproc++){
      DEFS::PhysicsProcess::Type indexType = DEFS::PhysicsProcess::getProcessType(string(vecCombined[iproc]));
       if(systematics[iSys].shape==true && systematics[iSys].getValueString(indexType)!="-")
          outFile << left << setw(13) << "1";
       else if(systematics[iSys].name.Contains("QCDEtaWeight") && systematics[iSys].getValueString(indexType)!="-")
          outFile << left << setw(27) << systematics[iSys].getValueString(indexType);
       else
          outFile << left << setw(13) << systematics[iSys].getValueString(indexType);
    } //for luminosity loop

    outFile << endl << flush;
  }
  cout << "DONE" << endl;

  cout << "Closing outFile ... ";
  outFile.close();
  cout << "DONE" << endl;

  cout << "Closing the in_file ... ";
  in_file->Close();
  cout << "DONE" << endl;

  std::terminate();
  return 0;
  
} //LimitsDatacard
