//Our libraries
#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"
#include "TAMUWW/SpecialTools/interface/Table.hh"
#include "TAMUWW/SpecialTools/interface/TableRow.hh"
#include "TAMUWW/SpecialTools/interface/TableCellVal.hh"
#include "TAMUWW/SpecialTools/interface/Value.hh"
#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"
#include "TAMUWW/SpecialTools/interface/MVAVar.hh"
#include "TAMUWW/Tools/interface/mymath.hh"
#include "TAMUWW/Tools/interface/Plots.hh"
#include "TAMUWW/Tools/interface/PlotFiller.hh"
#include "TAMUWW/Tools/interface/PUreweight.hh"
#include "TAMUWW/Tools/interface/CSVreweight.hh"
#include "TAMUWW/Tools/interface/TTbarreweight.hh"
#include "SHyFT/TemplateMakers/bin/AngularVars.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"
#include "JetMETAnalysis/JetUtilities/interface/TProfileMDF.h"

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
#include "TChain.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TError.h"

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

namespace UserFunctions
{
   TString outDir;
   TString signalTitle;
   DEFS::ControlRegion controlRegion;
   DEFS::TagCat tagCat;
   DEFS::LeptonCat leptonCat;
   DEFS::JetBin jetBin;
   DEFS::NtupleType ntupleType;
   PUreweight* puweight;
   CSVreweight* csvweight;
   TTbarreweight* ttbarweight;
   double elIsoMin;
   double elIsoMax;
   double muIsoMin;
   double muIsoMax;
   bool doJER;
   bool doMETPhiCorrection;
   bool doPUreweight;
   TString pileupSystematic;
   bool doCSVreweight;
   TString doCSVsys;
   bool doTTbarreweight;
   TString doTTbarsys;
   bool doFNAL;
   bool doMetPhiWeight;
   TH1D* metPhiWeight= 0;
   bool fillTMDF;
   bool  QCDweight;
   TH1D* QCDWeightFunc = 0;
   TString QCDweightSys;
   bool CosThetaLweight;
   TH1D* CosThetaLWeightFunc = 0;
   TString CosThetaLweightSys;
   bool WPtWeight;
   TH1D* WPtWeightFunc = 0;
   bool  PtWeight;
   TH1D* PtWeightFunc = 0;
   bool  ZllPtWeight;
   TH1D* ZllPtWeightFunc = 0;
   bool WJweight;
   TH2D* WJetsWeightFunc= 0;
   TF2 * WJetsWeightTF2 = 0;
   bool  fill2D = 0;
   bool  fillNonStandard = 0;
   bool  fillEPD = 0;
   bool  fillBumpCrossCheck = 0;
   bool  fillWJetsFlavor = 0;
   bool  fillLimitTemplatesOnly = 0;
   bool  fillGen = 0;
   int   limitBranches = 0;
   bool  createMVADiscriminatorPlot = 0;
   bool  doBCJetWeight;
   bool verbose; // adds or takes away cout statements when running
   map<int,pair<double, double> > maxEventProbs;
   bool doBenchmarks = false;
   TBenchmark* func_benchmark  = new TBenchmark();
   TBenchmark* once_benchmark  = new TBenchmark();

   ////////////////////////////////////////////////////////////////////////////////
   //  User Functions
   ////////////////////////////////////////////////////////////////////////////////

   // Is run once for each process before events are cut (initialize)
   void initEventFunc(EventNtuple* ntuple, const PhysicsProcess* proc);
   // this function fills all of the plots for a given process
   void fillPlots(MapOfPlots &  plots, TString processName, EventNtuple * ntuple, 
                  METree * metree, MicroNtuple * mnt, vector<TString>, 
                  map<TString, MVAVar>& evMap, double weight = 1.0);
   // compute the combined BDT
   double computeCombinedBDT(EventNtuple * ntuple, MicroNtuple * mnt, map<TString,
                             MVAVar>& evMap, vector<TString>, TString valueType = "response");
   // returns a boolean if the event passes the specified cuts
   bool eventPassCuts(EventNtuple * ntuple, MicroNtuple * mnt, const PhysicsProcess*,
                      map<TString,MVAVar>& evMap, vector<TString> MVAMethods);
                                  
   // returns a double
   double weightFunc(EventNtuple* ntuple, MicroNtuple* mnt, const PhysicsProcess* proc);
   // Is run once for each process before events (initializes PU Reweighting
   void processFunc(EventNtuple* ntuple, const PhysicsProcess* proc);
   // concat any two streamable objects into a string
   template <class T, class U>
   std::string concatString(const T& obj1, const U& obj2);
}

////////////////////////////////////////////////////////////////////////////////
//  Implement User Functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void UserFunctions::fillPlots(MapOfPlots &  plots, TString processName, EventNtuple * ntuple,
                              METree * metree, MicroNtuple * mnt,vector<TString> MVAMethods,
                              map<TString,MVAVar>& evMap, double weight)
{
   if(UserFunctions::doBenchmarks)
      UserFunctions::func_benchmark->Start("fillPlots");
   //Weird events with  ntuple->lLV[0].leptonCat different than electrons
   // or muons are killed in eventCuts so they never reach this stage
   DEFS::LeptonCat leptonCat = ntuple->lLV[0].leptonCat;

   //WARNING!!! Make sure that ntuple, metree, and mnt are set before using!!! 
   // Example::metree and mnt will not be set if you are running on an MEInput file.
   if (ntuple && !UserFunctions::fillLimitTemplatesOnly) {
      double WmT = sqrt(pow(ntuple->lLV[0].Et()+ntuple->METLV[0].Pt(), 2) -
                        pow(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(), 2) -
                        pow(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(), 2));
      double Mjj = 0;
      double Mlv = 0;

      if(!processName.Contains("Data") && processName.Contains("QCD"))//&& processName.Contains("Enriched"))
         plots[leptonCat]["PUWeights"]->Fill(puweight->getWeight(ntuple->vLV[0].npv));
      else if(!processName.Contains("Data") && !processName.Contains("QCD"))
         plots[leptonCat]["PUWeights"]->Fill(puweight->getWeight(ntuple->vLV[0].tnpus[1]));

      if (ntuple->jLV.size()>1) {
         Mjj = (ntuple->jLV[0] + ntuple->jLV[1]).M();
         Mlv = (ntuple->lLV[0] + ntuple->METLV[0]).M();
         plots[leptonCat]["Mjj"]->Fill(Mjj,weight);
         plots[leptonCat]["Mlv"]->Fill(Mlv,weight);
         plots[leptonCat]["MjjmWmT"]->Fill(Mjj - WmT, weight);
         plots[leptonCat]["j1Pt_Mjj"]->Fill(ntuple->jLV[0].Pt() / ntuple->Mjj,weight);
         plots[leptonCat]["j2Pt_Mjj"]->Fill(ntuple->jLV[1].Pt() / ntuple->Mjj,weight);
         plots[leptonCat]["Jet2Pt"]->Fill(ntuple->jLV[1].Pt(),weight);
         plots[leptonCat]["Jet2Eta"]->Fill(ntuple->jLV[1].Eta(),weight);
         plots[leptonCat]["Jet2Phi"]->Fill(ntuple->jLV[1].Phi(),weight);
         plots[leptonCat]["DeltaEtaJ1J2"]->Fill(TMath::Abs(ntuple->jLV[0].Eta() - ntuple->jLV[1].Eta()),weight);
         plots[leptonCat]["Ptjj"]->Fill((ntuple->jLV[0] + ntuple->jLV[1]).Pt(),weight);
         plots[leptonCat]["Mlvjj"]->Fill((ntuple->jLV[0] + ntuple->jLV[1] + ntuple->lLV[0] 
                                          + ntuple->METLV[0]).M(),weight);
         plots[leptonCat]["DeltaPhi_J1J2"]->Fill(ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]),weight);
         plots[leptonCat]["DeltaRJ1J2"]->Fill(sqrt(pow(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta(),2)+
                                                   pow(ntuple->jLV[0].Phi()-ntuple->jLV[1].Phi(),2)),weight);
         if(UserFunctions::fillNonStandard) {
            plots[leptonCat]["EJ1EJ2"]->Fill(ntuple->jLV[0].E() * ntuple->jLV[1].E(),weight);
            plots[leptonCat]["BetaJ1BetaJ2"]->Fill(ntuple->jLV[0].Beta() * ntuple->jLV[1].Beta(),weight);
            plots[leptonCat]["AngleJ1J2"]->Fill(ntuple->jLV[0].Angle(ntuple->jLV[1].Vect()),weight);
            plots[leptonCat]["jjlvPhi"]->Fill((ntuple->jLV[0] + ntuple->jLV[1]).Phi() -
                                              (ntuple->lLV[0] + ntuple->METLV[0]).Phi(),weight);
         }
         if(UserFunctions::fill2D) {
            plots[leptonCat]["DeltaPhi_LJ1_vs_J1J2"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]), 
                                                           ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]),weight);
            plots[leptonCat]["DeltaEta_LJ1_vs_J1J2"]->Fill(ntuple->lLV[0].Eta()-ntuple->jLV[0].Eta(), 
                                                           ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta(),weight);
            plots[leptonCat]["DeltaPhi_LJ1_vs_J1J2_Subtracted"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]),
                                                                      ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]), 
                                                                      (ntuple->lLV[0].lQ)*weight);
         }
      }
      plots[leptonCat]["LeptPt"]->Fill(ntuple->lLV[0].Pt(),weight);
      plots[leptonCat]["LeptEta"]->Fill(ntuple->lLV[0].Eta(),weight);
      plots[leptonCat]["LeptPhi"]->Fill(ntuple->lLV[0].Phi(),weight);
      plots[leptonCat]["LeptPFIso"]->Fill(ntuple->lLV[0].lpfIso ,weight);
      plots[leptonCat]["MET"]->Fill(ntuple->METLV[0].Pt(),weight);
      plots[leptonCat]["METPhi"]->Fill(ntuple->METLV[0].Phi(),weight);
      plots[leptonCat]["Ptlv"]->Fill((ntuple->lLV[0]+ntuple->METLV[0]).Pt(),weight);
      plots[leptonCat]["WmT"]->Fill(WmT, weight);
      plots[leptonCat]["Jet1Pt"]->Fill(ntuple->jLV[0].Pt(),weight);
      plots[leptonCat]["Jet1Eta"]->Fill(ntuple->jLV[0].Eta(),weight);
      plots[leptonCat]["Jet1Phi"]->Fill(ntuple->jLV[0].Phi(),weight);
      plots[leptonCat]["DeltaPhi_LJ1"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]),weight);
      plots[leptonCat]["DeltaPhi_METJ1"]->Fill(ntuple->jLV[0].DeltaPhi(ntuple->METLV[0]),weight);
      plots[leptonCat]["DeltaPhi_LMET"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->METLV[0]),weight);
      plots[leptonCat]["npv"]->Fill(ntuple->vLV[0].npv,weight);

      if(UserFunctions::fillNonStandard) {
         plots[leptonCat]["SigmaIetaIeta"]->Fill(ntuple->lLV[0].eSigmaIetaIeta,weight);
         plots[leptonCat]["MVATrig"]->Fill(ntuple->lLV[0].emvaTrig,weight);
         plots[leptonCat]["MVANonTrig"]->Fill(ntuple->lLV[0].emvaNonTrig,weight);
         plots[leptonCat]["DeltaRLepMET"]->Fill(sqrt(pow(ntuple->lLV[0].Eta()-ntuple->METLV[0].Eta(),2)+
                                                     pow(ntuple->lLV[0].Phi()-ntuple->METLV[0].Phi(),2) ),weight);
         for (unsigned int j=0; j<ntuple->jLV.size() && j<32; j++) {
            plots[leptonCat][string(Form("JetEta_%luJets",ntuple->jLV.size()))]->Fill(ntuple->jLV[j].Eta(),weight);
            plots[leptonCat]["nJets_JetEta"]->Fill(ntuple->jLV[j].Eta(),ntuple->jLV.size(),weight);
         }
         plots[leptonCat]["WmT_Positive"]->Fill(WmT, weight);
         plots[leptonCat]["WmT_Negative"]->Fill(WmT, weight);
         plots[leptonCat]["WmT_Subtracted"]->Fill(WmT, (ntuple->lLV[0].lQ)*weight);
         if(processName.Contains("TTbar"))
            plots[leptonCat]["TTbarWeights"]->Fill(ttbarweight->getWeight(ntuple,0),1.0);
            //plots[leptonCat]["TTbarWeights"]->Fill(ttbarWeight(ntuple,0),1.0);
      }

      if(UserFunctions::fill2D) {
         plots[leptonCat]["MET_vs_LeptEta"]->Fill(ntuple->lLV[0].Eta(),
                                                  ntuple->METLV[0].Pt(),
                                                  weight);
         plots[leptonCat]["MET_vs_LeptPt"]->Fill(ntuple->lLV[0].Pt(),
                                                 ntuple->METLV[0].Pt(),
                                                 weight); 
         plots[leptonCat]["MET_vs_AbsLeptEta"]->Fill(fabs(ntuple->lLV[0].Eta()),
                                                     ntuple->METLV[0].Pt(),
                                                     weight);
         plots[leptonCat]["IsolationEnergyVsPt"]->Fill(ntuple->lLV[0].Pt(),ntuple->lLV[0].lpfIso*ntuple->lLV[0].Pt());
         plots[leptonCat]["IsolationEnergyVsPtCoshEta"]->Fill(ntuple->lLV[0].Pt()*TMath::CosH(ntuple->lLV[0].Eta()),ntuple->lLV[0].lpfIso*ntuple->lLV[0].Pt());
         pair<double,double> leptVsHadWMass = ntuple->onVsOffShellInclusive(verbose);
         plots[leptonCat]["MWjjVsMWlv"]->Fill(leptVsHadWMass.first,leptVsHadWMass.second,weight);
         if (ntuple->lLV[0].lQ == 1 && ntuple->jLV.size()>1){
            plots[leptonCat]["DeltaPhi_LJ1_vs_J1J2_Positive"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]), 
                                                                    ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]), weight);
            plots[leptonCat]["DeltaEta_LJ1_vs_J1J2_Positive"]->Fill(ntuple->lLV[0].Eta()-ntuple->jLV[0].Eta(), 
                                                                    ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta(),weight);
         }
         if (ntuple->lLV[0].lQ == -1 && ntuple->jLV.size()>1){
            plots[leptonCat]["DeltaPhi_LJ1_vs_J1J2_Negative"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]),
                                                                    ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]), weight);
            plots[leptonCat]["DeltaEta_LJ1_vs_J1J2_Negative"]->Fill(ntuple->lLV[0].Eta()-ntuple->jLV[0].Eta(), 
                                                                    ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta(),weight);
         }
      }

      if (UserFunctions::fillTMDF){
         vector<Double_t> coord;
         //coord.assign(((TProfileMDF*)plots[leptonCat]["lpt_leta_j1pt_j1eta_j2pt_j2eta"]->templateHisto)->GetNaxis(),0);
         //coord.assign(((TProfileMDF*)plots[leptonCat]["lpt_lpt_j1pt_j1pt_j2pt_j2pt"]->templateHisto)->GetNaxis(),0);
         coord.assign(((TProfileMDF*)plots[leptonCat]["Mjj_Mjj_Mt_MET_DeltaR_DeltaR"]->templateHisto)->GetNaxis(),0);
         //coord[0] = ntuple->lLV[0].Pt();
         coord[0] = Mjj;
         //coord[1] = TMath::Abs(ntuple->lLV[0].Eta());
         //coord[1] = ntuple->lLV[0].Pt();
         coord[1] = Mjj;
         //coord[2] = ntuple->jLV[0].Pt();
         coord[2] = WmT;
         //coord[3] = TMath::Abs(ntuple->jLV[0].Eta());
         //coord[3] = ntuple->jLV[0].Pt();
         coord[3] = ntuple->METLV[0].Pt();
         //coord[4] = ntuple->jLV[1].Pt();
         coord[4] = sqrt(pow(ntuple->lLV[0].Eta()-ntuple->jLV[0].Eta(),2)
                         +pow(ntuple->lLV[0].Phi()-ntuple->jLV[0].Phi(),2));
         //coord[5] = TMath::Abs(ntuple->jLV[1].Eta());
         //coord[5] = ntuple->jLV[1].Pt();
         coord[5] = sqrt(pow(ntuple->lLV[0].Eta()-ntuple->jLV[0].Eta(),2)
                         +pow(ntuple->lLV[0].Phi()-ntuple->jLV[0].Phi(),2));
         //plots[leptonCat]["lpt_leta_j1pt_j1eta_j2pt_j2eta"]->Fill(coord,1.0,weight);
         //plots[leptonCat]["lpt_lpt_j1pt_j1pt_j2pt_j2pt"]->Fill(coord,1.0,weight);
         plots[leptonCat]["Mjj_Mjj_Mt_MET_DeltaR_DeltaR"]->Fill(coord,1.0,weight);
      }

      if (UserFunctions::fillGen && !processName.Contains("Data") && !processName.Contains("QCD")) {
         TLorentzVector lepton;
         TLorentzVector boson;
         if(processName.Contains("ZJets") && !processName.Contains("ZJetsToLL")) {
            //lepton = ntuple->getGenVorDaughter(EventNtuple::LEPTON, 23, 0, false);
            //boson = ntuple->getGenVorDaughter(EventNtuple::BOSON, 23, 0, false);
         }
         else if(processName.Contains("WJets") || processName.Contains("ZJetsToLL") || processName.Contains("TTbar") ||
                 processName.Contains("STop") || processName.Contains("WW") || processName.Contains("WZ")) {
            lepton = ntuple->getGenVorDaughter(EventNtuple::LEPTON, 24, 0, false);
            boson = ntuple->getGenVorDaughter(EventNtuple::BOSON, 24, 0, false);
         }
         else {
            lepton = ntuple->getGenVorDaughter(EventNtuple::LEPTON, 25, 0, false);
            boson = ntuple->getGenVorDaughter(EventNtuple::BOSON, 25, 0, false);
         }
         plots[leptonCat]["GenLeptPt"]->Fill(lepton.Pt());
         plots[leptonCat]["GenLeptEta"]->Fill(lepton.Eta());
         plots[leptonCat]["GenLeptPhi"]->Fill(lepton.Phi());
         plots[leptonCat]["GenWmT"]->Fill(boson.Mt());
         plots[leptonCat]["GenVPt"]->Fill(boson.Pt());
         plots[leptonCat]["GenVEta"]->Fill(boson.Eta());
         plots[leptonCat]["GenVPhi"]->Fill(boson.Phi());
      }

   }
   else if (ntuple && UserFunctions::fillLimitTemplatesOnly) {
      plots[leptonCat]["MET"]->Fill(ntuple->METLV[0].Pt(),weight);
   }

   if (metree && !UserFunctions::fillLimitTemplatesOnly) {
      if(!UserFunctions::fillLimitTemplatesOnly) {
         //double norms[9] = {0.781996e-6,76.4008e-9,0.1939e-3,42.31735e-6,0.12688e-3,5.85109e-6,0.592615,0.66255e-9,0.208542e-9};
         //TString normNames[9] = {"WW","WZ","WLg","WLgsub","Wgg","WLL","QCD","ggH125","WH125"};
         for (int tep=0; tep<metree->getNProbStat(); tep++) {
            string name = UserFunctions::concatString("tEventProb",tep);
            if((int)maxEventProbs.size()!=metree->getNProbStat()) {
               maxEventProbs[tep] = DefaultValues::getMaxEventProbAndError(tep);
            }
   
            if(metree->getProbStat(tep)->tEventProb == 0) {
               cout << "ERROR::plotter_x tEventProb" << tep << " == 0" << endl
                    << "This will cause issues down the line" << endl
                    << "Please fix this and try again" << endl;
               //assert(metree->getProbStat(tep)->tEventProb!=0);
               continue;
            }
   
            if(maxEventProbs[tep].first>0) {
               plots[leptonCat][name]->Fill(TMath::Log10(metree->getProbStat(tep)->tEventProb)-
                                            TMath::Log10(maxEventProbs[tep].first),weight);
            }
            else {
               plots[leptonCat][name]->Fill(TMath::Log10(metree->getProbStat(tep)->tEventProb),weight);
            }
         }
      }
      if(UserFunctions::fillEPD) {
         if (metree->getNProbStat()>=3) {
            double tEventProb0 = metree->getProbStat(0)->tEventProb;
            double tEventProb1 = metree->getProbStat(1)->tEventProb;
            double tEventProb3 = metree->getProbStat(3)->tEventProb;
   
            double signal = (tEventProb0 / 0.8e-06) + (tEventProb1 / 0.1e-06);
            double back   = tEventProb3 / 0.75e-03;
            double epd    = signal /(signal+back);
            plots[leptonCat]["epdPretagWWandWZ_RE"]->Fill(epd,weight);
         }
      }
   }

   if (mnt) {
      if(UserFunctions::fillEPD && !UserFunctions::fillLimitTemplatesOnly) {
         plots[leptonCat]["epdPretagWWandWZ"]->Fill(-TMath::Log10(mnt->epdPretagWWandWZ),weight);
         plots[leptonCat]["epdPretagHiggs125"]->Fill(-TMath::Log10(mnt->epdPretagHiggs[6]),weight);
      }
      plots[leptonCat]["MEBDT"]->Fill(mnt->MEBDT,weight);
      plots[leptonCat]["KinBDT"]->Fill(mnt->KinBDT,weight);
      plots[leptonCat]["KinMEBDT"]->Fill(mnt->KinMEBDT,weight);
      if(processName.Contains("WJets") && UserFunctions::fillWJetsFlavor) {
         double minpt = 25.0;
         plots[leptonCat]["KinMEBDT_WpInclusive"]->Fill(mnt->KinMEBDT,weight);
         if(!ntuple->containsParticle(4,true,minpt) && !ntuple->containsParticle(5,true,minpt))
            plots[leptonCat]["KinMEBDT_WpLight"]->Fill(mnt->KinMEBDT,weight);
         else if(ntuple->containsParticle(4,true,minpt) && !ntuple->containsParticle(5,true,minpt))
            plots[leptonCat]["KinMEBDT_WpC"]->Fill(mnt->KinMEBDT,weight);
         else if(ntuple->containsParticle(5,true,minpt))
            plots[leptonCat]["KinMEBDT_WpB"]->Fill(mnt->KinMEBDT,weight);
         else {
            plots[leptonCat]["KinMEBDT_WpOther"]->Fill(mnt->KinMEBDT,weight);
         }
      }
      if(!UserFunctions::fillLimitTemplatesOnly && ntuple) {
         computeCombinedBDT(ntuple,mnt,evMap,MVAMethods,"response");
         plots[leptonCat]["jet1dRLep"]->Fill(evMap["jet1dRLep"].value,weight);
         plots[leptonCat]["jet2dRLep"]->Fill(evMap["jet2dRLep"].value,weight);
         plots[leptonCat]["jet3dRLep"]->Fill(evMap["jet3dRLep"].value,weight);
         plots[leptonCat]["leptonEtaCharge"]->Fill(evMap["leptonEtaCharge"].value,weight);
         plots[leptonCat]["ht"]->Fill(evMap["ht"].value,weight);
         plots[leptonCat]["Ptlnujj"]->Fill(evMap["Ptlnujj"].value,weight);
         plots[leptonCat]["dPhiMETJet"]->Fill(evMap["dPhiMETJet"].value,weight);
         plots[leptonCat]["dPhiMETLep"]->Fill(evMap["dPhiMETLep"].value,weight);
         plots[leptonCat]["dRlepjj"]->Fill(evMap["dRlepjj"].value,weight);
         plots[leptonCat]["minDPhiLepJet"]->Fill(evMap["minDPhiLepJet"].value,weight);
         plots[leptonCat]["CosTheta_l"]->Fill(evMap["CosTheta_l"].value,weight);
         plots[leptonCat]["CosTheta_j"]->Fill(evMap["CosTheta_j"].value,weight);
         plots[leptonCat]["CosTheta_WH"]->Fill(evMap["CosTheta_WH"].value,weight);
      }
      if(UserFunctions::fillBumpCrossCheck) {
         plots[leptonCat]["MEBDT_vs_KinBDT"]->Fill(mnt->KinBDT,mnt->MEBDT,weight);
      }

      if (!MVAMethods.empty() && ntuple) {
         plots[leptonCat]["MVADiscriminator"]->Fill(computeCombinedBDT(ntuple,mnt,evMap,MVAMethods,"response"),weight);
         if(!UserFunctions::fillLimitTemplatesOnly) {
            plots[leptonCat]["MVAProbability"]->Fill(computeCombinedBDT(ntuple,mnt,evMap,MVAMethods,"probability"),weight);
         }
         if(UserFunctions::fillBumpCrossCheck) {
            plots[leptonCat]["MVADiscriminator_fineBinning"]->Fill(computeCombinedBDT(ntuple,mnt,evMap,MVAMethods,"response"),weight);
         }
         if(UserFunctions::fillBumpCrossCheck && UserFunctions::limitBranches == 0) {
            //plots[leptonCat]["MVADiscriminator_fineBinning_0GenMatch"]->Fill(,weight);
            //plots[leptonCat]["MVADiscriminator_fineBinning_1GenMatch"]->Fill(,weight);
            //plots[leptonCat]["MVADiscriminator_fineBinning_2GenMatch"]->Fill(,weight);
         }
      }
   }
   if(UserFunctions::doBenchmarks) {
      UserFunctions::func_benchmark->Stop("fillPlots");
      float rt = 0, ct = 0;
      cout << endl << "PlotFiller::func_benchmark" << endl;
      vector<string> timers; 
      timers.push_back("weightFunc");
      timers.push_back("initEventFunc");
      timers.push_back("processFunc");
      timers.push_back("fillPlots");
      DefaultValues::printSummary(func_benchmark, 8, rt, ct, timers);
      func_benchmark->Reset();
   }
}//fillPlots

//______________________________________________________________________________
// compute the combined BDT
double UserFunctions::computeCombinedBDT(EventNtuple * ntuple, MicroNtuple * mnt, map<TString,
                                         MVAVar>& evMap, vector<TString> MVAMethods, TString valueType) {
   //Alexx's BDT Spectators
   evMap["event"].value = ntuple->event;
   evMap["lumi"].value = ntuple->lumi;
   evMap["run"].value = ntuple->run;

   //Joey's BDT Spectators
   //evMap["csvWeight_lep"].value = csvweight->getWeight(ntuple);
   //evMap["ttbarWeight"].value = ttbarweight->getWeight(ntuple,0);
   //evMap["puWeight_lep"].value = puweight->getWeight(ntuple->vLV[0].tnpus[1]);

   //Joey's BDT Variables
   evMap["leptonEtaCharge"].value = ntuple->lLV[0].Eta()*ntuple->lLV[0].lQ;
   evMap["leptonPt"].value = ntuple->lLV[0].Pt();
   evMap["lepMT"].value = TMath::Sign(TMath::Sqrt(TMath::Abs(TMath::Power(ntuple->lLV[0].Pt()+ntuple->METLV[0].Pt(),2)-
                                      (TMath::Power(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->lLV[0].Py()+
                                       ntuple->METLV[0].Py(),2)+TMath::Power(0,2)))),TMath::Power(ntuple->lLV[0].Pt()+ntuple->METLV[0].Pt(),2)-
                                       (TMath::Power(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+TMath::Power(ntuple->lLV[0].Py()+
                                        ntuple->METLV[0].Py(),2)+TMath::Power(0,2)));
   evMap["jet1dRLep"].value = ntuple->jLV[0].DRlj;
   evMap["jet2dRLep"].value = ntuple->jLV[1].DRlj;
   evMap["jet3dRLep"].value = ntuple->jLV[2].DRlj;
   double tmpHT = 0.0;
   for(unsigned int i=0; i<ntuple->jLV.size(); i++) {
      tmpHT+=ntuple->jLV[i].Et();
   }
   tmpHT+=ntuple->lLV[0].Pt();
   evMap["ht"].value = tmpHT;
   evMap["Ptlnujj"].value = TMath::Sqrt(TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+
                                        TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2));
   evMap["Mlnujj"].value = TMath::Sign(TMath::Sqrt(TMath::Abs(TMath::Power(ntuple->jLV[0].E()+ntuple->jLV[1].E()+ntuple->lLV[0].E()+ntuple->METLV[0].E(),2)-
                                       (TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+
                                        TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+
                                        TMath::Power(ntuple->jLV[0].Pz()+ntuple->jLV[1].Pz()+ntuple->lLV[0].Pz()+ntuple->METLV[0].Pz(),2)))),
                                       TMath::Power(ntuple->jLV[0].E()+ntuple->jLV[1].E()+ntuple->lLV[0].E()+ntuple->METLV[0].E(),2)-
                                       (TMath::Power(ntuple->jLV[0].Px()+ntuple->jLV[1].Px()+ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),2)+
                                        TMath::Power(ntuple->jLV[0].Py()+ntuple->jLV[1].Py()+ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),2)+
                                        TMath::Power(ntuple->jLV[0].Pz()+ntuple->jLV[1].Pz()+ntuple->lLV[0].Pz()+ntuple->METLV[0].Pz(),2)));
   evMap["dRlepjj"].value = mnt->dRlepjj;
   evMap["dPhiMETJet"].value = mnt->dPhiMETJet;
   evMap["dPhiJetJet"].value = mnt->dPhiJetJet;
   evMap["dPhiMETLep"].value = TMath::Abs(ntuple->METLV[0].Phi()-ntuple->lLV[0].Phi());
   evMap["minDPhiLepJet"].value = mnt->minDPhiLepJet;
   evMap["dEtaJetJet"].value = TMath::Abs(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta());
   evMap["CosTheta_l"].value = mnt->CosTheta_l;
   evMap["CosTheta_j"].value = mnt->CosTheta_j;
   evMap["CosTheta_WH"].value = mnt->CosTheta_WH;
   evMap["MEBDT"].value = mnt->MEBDT;
   //cout << "nJets => " << ntuple->jLV.size() << endl;
   //cout << "nBTag => " << ntuple->getNBTags() << endl;
   //cout << "jet1 => pt=" << ntuple->jLV[0].Pt() << " eta=" << ntuple->jLV[0].Eta() << " phi=" << ntuple->jLV[0].Phi() << " et=" << ntuple->jLV[0].Et() << endl;
   //cout << "jet2 => pt=" << ntuple->jLV[1].Pt() << " eta=" << ntuple->jLV[1].Eta() << " phi=" << ntuple->jLV[1].Phi() << " et=" << ntuple->jLV[1].Et() << endl;
   //for (std::map<TString,MVAVar>::iterator it=evMap.begin(); it!=evMap.end(); ++it)
   //   if(it->second.use)
   //      std::cout << it->first << " => " << it->second.value << '\n';
   //ME Extra BDT Variables

   //mnt->setExtraVarMVA(evMap);
   if(MVAMethods.empty()) {
      return -9999;
   }
   else {
      assert(!TMath::IsNaN(mnt->getMVAOutput(MVAMethods).front()[valueType]));
      return mnt->getMVAOutput(MVAMethods).front()[valueType];
   }
}

//______________________________________________________________________________
// Return true if the event pass the cuts imposed to the given lepton category
bool UserFunctions::eventPassCuts(EventNtuple * ntuple, MicroNtuple * mnt, const PhysicsProcess* proc,
                                  map<TString,MVAVar>& evMap, vector<TString> MVAMethods){  

   // Remove events categorized as "both" or "none", leaving only base categories of electron and muons
   if ( ntuple->lLV[0].leptonCat == DEFS::both ||  ntuple->lLV[0].leptonCat == DEFS::none)
      return false;

   // Remove events that have the wrong number of jets
   if ( (ntuple->jLV.size() != 1 && jetBin == DEFS::jet1) || (ntuple->jLV.size() != 2 && jetBin == DEFS::jets2) || (ntuple->jLV.size() != 3 && jetBin == DEFS::jets3) || (ntuple->jLV.size() < 4 && jetBin == DEFS::jets4))
      return false;

   // Remove events that should not be used from this particular process.
   //if (proc->leptonCat != DEFS::both && proc->leptonCat != ntuple->lLV[0].leptonCat)
   PhysParMap ppm = proc->sigma;
   if (ppm[ntuple->lLV[0].leptonCat] == 0.0)
      return false;

   // Remove events that are not requested from the command line  
   if ( ntuple->lLV[0].leptonCat != leptonCat && (leptonCat != DEFS::both) )
      return false;

   //MET Cut
   if ( ntuple->METLV[0].Pt() <= 25.0 )
      return false;

   // PFISO cut for QCD samples
   if (proc->name.Contains("QCD") && (proc->name.Contains("ElFULL") || proc->name.Contains("ElEnriched")) && (ntuple->lLV[0].lpfIso <= elIsoMin || ntuple->lLV[0].lpfIso >= elIsoMax) )
      return false;
   if (proc->name.Contains("QCD") && (proc->name.Contains("MuFULL") || proc->name.Contains("MuEnriched")) && (ntuple->lLV[0].lpfIso <= muIsoMin || ntuple->lLV[0].lpfIso >= muIsoMax) )
     return false;

   //Implement FNAL cuts
   if (doFNAL){
      if (ntuple->lLV[0].leptonCat == DEFS::muon){
         if(!ntuple->FNALcutsMuon())
            return false;
      }
      else if(ntuple->lLV[0].leptonCat == DEFS::electron)
         if (!ntuple->FNALcutsElectron())
            return false;
   }

   // regardless of cut region cut on minimum lepton Pt 
   //if(ntuple->lLV[0].leptonCat == DEFS::electron && ntuple->lLV[0].Pt() < 45)
   if(ntuple->lLV[0].leptonCat == DEFS::electron && ntuple->lLV[0].Pt() <= 30)
     return false;
   
   if(ntuple->lLV[0].leptonCat == DEFS::muon && ntuple->lLV[0].Pt() <= 25)
     return false;

   //X axis cuts
   if (controlRegion == DEFS::signal || controlRegion == DEFS::BDTBump || controlRegion == DEFS::BDTAntiBump ||
       controlRegion == DEFS::MVAEleID || controlRegion == DEFS::AntiMVAEleID || controlRegion == DEFS::FlatMVAEleID ||
       controlRegion == DEFS::HailMary || controlRegion == DEFS::HailMaryLoose) {

      //test if cutting much harder on eta of lepton will work
      //if(ntuple->lLV[0].Eta()>=1.3)
      //   return false;
      //if(ntuple->lLV[0].Eta()<1.3)
      //   return false;

     // leading jet with PT > 30, and second leading with at least 25 GeV
     if (ntuple->jLV[0].Pt() <= 30 || ntuple->jLV[1].Pt() <= 25)
        return false;

     if (controlRegion == DEFS::BDTBump || controlRegion == DEFS::BDTAntiBump) {
         double combinedBDT = computeCombinedBDT(ntuple,mnt,evMap,MVAMethods,"response");
         if(controlRegion == DEFS::BDTBump && combinedBDT<0.2)
            return false;
         else if(controlRegion == DEFS::BDTAntiBump && combinedBDT>0.0)
            return false;
     }

     if (controlRegion == DEFS::FlatMVAEleID) {
        if (ntuple->lLV[0].emvaTrig > 0.95 || ntuple->lLV[0].emvaTrig < 0.05) {
          return false;
        }
     }

     if (controlRegion == DEFS::MVAEleID) {
        double aeta = TMath::Abs(ntuple->lLV[0].Eta());
        double emvaTrig = ntuple->lLV[0].emvaTrig;
        double pfIso = ntuple->lLV[0].lpfIso;
        if ((aeta<0.8 && pfIso<0.093 && emvaTrig < 0.977) ||
            (aeta>0.8 && aeta<1.479 && pfIso<0.095 && emvaTrig<0.956) ||
            (aeta>1.479 && aeta<2.5 && pfIso<0.171 && emvaTrig<0.966) ) {
          return false;
        }
     }

     if (controlRegion == DEFS::AntiMVAEleID) {
        double aeta = TMath::Abs(ntuple->lLV[0].Eta());
        double emvaTrig = ntuple->lLV[0].emvaTrig;
        double pfIso = ntuple->lLV[0].lpfIso;
        if ((aeta<0.8 && pfIso<0.093 && emvaTrig > 0.977) ||
            (aeta>0.8 && aeta<1.479 && pfIso<0.095 && emvaTrig>0.956) ||
            (aeta>1.479 && aeta<2.5 && pfIso<0.171 && emvaTrig>0.966) ) {
          return false;
        }
     }

     int nBtag = ntuple->getNBTags();
     if(tagCat != DEFS::pretag){
 
        if(tagCat == DEFS::eq0tag && nBtag!=0) {
           return false;
        }
        if(tagCat == DEFS::eq1tag && nBtag!=1) {
           return false;
        }
        if(tagCat == DEFS::eq2tag && nBtag!=2) {
           return false;
        }
        if(tagCat == DEFS::ge1tag && nBtag<1) {
           return false;
        }
     }

     if (controlRegion == DEFS::HailMary) {
        double dRj1j2 = sqrt(pow(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta(),2)+pow(ntuple->jLV[0].Phi()-ntuple->jLV[1].Phi(),2));
        double ptjj = (ntuple->jLV[0] + ntuple->jLV[1]).Pt();
        double WmT = sqrt(pow(ntuple->lLV[0].Et()+ntuple->METLV[0].Pt(), 2) -
                        pow(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(), 2) -
                        pow(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(), 2));
        if (ntuple->lLV[0].leptonCat == DEFS::electron){
           if(ntuple->lLV[0].Eta() > 1.3                                   ||
              ntuple->lLV[0].Pt() < 37                                     ||
              ntuple->lLV[0].Pt() > 100                                    ||
              ntuple->jLV[0].Eta() > 2.0                                   ||
              WmT > 100                                                    ||
              ntuple->METLV[0].Pt() > 100                                  ||
              ntuple->vLV[0].npv > 25                                      ||
              dRj1j2 < 1.0                                                 ||
              ptjj < 10.0                                                  ||
              mnt->CosTheta_l > 0.6                                        ||
              mnt->CosTheta_j < 0.5                                        ||
              mnt->CosTheta_WH > 0.5                                       ||
              (ntuple->jLV.size()>=2 && ntuple->jLV[1].DRlj > TMath::Pi()) ||
              (ntuple->jLV.size()>=3 && ntuple->jLV[2].DRlj > TMath::Pi()) ||
              (ntuple->jLV.size()>=4 && ntuple->jLV[3].DRlj > TMath::Pi()) ||
              TMath::Log10(mnt->eventProb[19]) < -20) {
               return false;
           }
        }
        if (ntuple->lLV[0].leptonCat == DEFS::muon) {
            if(ntuple->lLV[0].Eta() > 1.3                                  ||
              ntuple->jLV[0].Eta() > 2.0                                   ||
              WmT > 100                                                    ||
              ntuple->vLV[0].npv > 25                                      ||
              dRj1j2 < 1.0                                                 ||
              ptjj > 80.0                                                  ||
              mnt->CosTheta_l > 0.6                                        ||
              mnt->CosTheta_j < 0.5                                        ||
              mnt->CosTheta_WH > 0.5                                       ||
              (ntuple->jLV.size()>=1 && ntuple->jLV[0].DRlj > TMath::Pi()) ||
              (ntuple->jLV.size()>=2 && ntuple->jLV[1].DRlj > TMath::Pi()) ||
              (ntuple->jLV.size()>=3 && ntuple->jLV[2].DRlj > TMath::Pi()) ||
              (ntuple->jLV.size()>=4 && ntuple->jLV[3].DRlj > TMath::Pi()) ||
              TMath::Log10(mnt->eventProb[19]) < -20) {
               return false;
           }
        }
     }
     if (controlRegion == DEFS::HailMaryLoose) {
        double dRj1j2 = sqrt(pow(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta(),2)+pow(ntuple->jLV[0].Phi()-ntuple->jLV[1].Phi(),2));
        double ptjj = (ntuple->jLV[0] + ntuple->jLV[1]).Pt();
        double WmT = sqrt(pow(ntuple->lLV[0].Et()+ntuple->METLV[0].Pt(), 2) -
                        pow(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(), 2) -
                        pow(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(), 2));
        if (ntuple->lLV[0].leptonCat == DEFS::electron){
           if(ntuple->lLV[0].Eta() > 1.3                                   ||
              ntuple->lLV[0].Pt() < 37                                     ||
              ntuple->lLV[0].Pt() > 100                                    ||
              ntuple->jLV[0].Eta() > 2.0                                   ||
              //WmT > 100                                                    ||
              //ntuple->METLV[0].Pt() > 100                                  ||
              ntuple->vLV[0].npv > 25                                      ||
              dRj1j2 < 1.0                                                 ||
              ptjj < 10.0                                                  ||
              //mnt->CosTheta_l > 0.6                                        ||
              mnt->CosTheta_l > 0.7                                        ||
              //mnt->CosTheta_j < 0.5                                        ||
              mnt->CosTheta_j < 0.0                                        ||
              //mnt->CosTheta_WH > 0.5                                       ||
              (ntuple->jLV.size()>=2 && ntuple->jLV[1].DRlj > TMath::Pi()) ||
              (ntuple->jLV.size()>=3 && ntuple->jLV[2].DRlj > TMath::Pi()) ||
              (ntuple->jLV.size()>=4 && ntuple->jLV[3].DRlj > TMath::Pi()) ||
              TMath::Log10(mnt->eventProb[19]) < -20) {
               return false;
           }
        }
        if (ntuple->lLV[0].leptonCat == DEFS::muon) {
            if(ntuple->lLV[0].Eta() > 1.3                                  ||
              ntuple->jLV[0].Eta() > 2.0                                   ||
              //WmT > 100                                                    ||
              ntuple->vLV[0].npv > 25                                      ||
              dRj1j2 < 1.0                                                 ||
              ptjj > 80.0                                                  ||
              //mnt->CosTheta_l > 0.6                                        ||
              mnt->CosTheta_l > 0.7                                        ||
              //mnt->CosTheta_j < 0.5                                        ||
              mnt->CosTheta_j < 0.0                                        ||
              //mnt->CosTheta_WH > 0.5                                       ||
              (ntuple->jLV.size()>=1 && ntuple->jLV[0].DRlj > TMath::Pi()) ||
              (ntuple->jLV.size()>=2 && ntuple->jLV[1].DRlj > TMath::Pi()) ||
              (ntuple->jLV.size()>=3 && ntuple->jLV[2].DRlj > TMath::Pi()) ||
              (ntuple->jLV.size()>=4 && ntuple->jLV[3].DRlj > TMath::Pi()) ||
              TMath::Log10(mnt->eventProb[19]) < -20) {
               return false;
           }
        }
     }
   }
   else if (controlRegion == DEFS::control1){
     
     if (ntuple->jLV[0].Pt() > 30 || ntuple->jLV[1].Pt() > 25)
       return false;
   }
   else if (controlRegion == DEFS::control2){
      if ((ntuple->jLV[0].Pt() > 30 && ntuple->jLV[1].Pt() > 25) || 
          (ntuple->jLV[0].Pt() < 30 && ntuple->jLV[1].Pt() < 25))
         return false;
   }
   else if (controlRegion == DEFS::control3){
     if(ntuple->jLV[0].Pt() > 30 && ntuple->jLV[1].Pt() > 25)
       return false;
   }
   else if (controlRegion == DEFS::control4){
      if(ntuple->jLV.size()<4)
         return false;
   }
   else if (controlRegion == DEFS::control5){
      int nBtag = 0;
      for(unsigned int j=0; j<ntuple->jLV.size();j++) {
         if(ntuple->jLV[j].jBtagCSV==1)
            nBtag++;
      }
      if(ntuple->jLV.size()<4 || nBtag>0) {
         return false;
      }
   }
   else if (controlRegion == DEFS::control6){
      if(ntuple->jLV.size()!=1)
         return false;
   }
   else if (controlRegion == DEFS::control7){
      if(ntuple->jLV.size()!=2)
         return false;
   }
   else if (controlRegion == DEFS::control8){
      if(ntuple->jLV.size()!=3)
         return false;
   }
   else if (controlRegion == DEFS::control9){
      if(ntuple->jLV.size()!=1 && ntuple->jLV.size()<4)
         return false;
   }
   else if (controlRegion == DEFS::UVaCuts){
      if(ntuple->lLV[0].leptonCat == DEFS::electron){
         if(ntuple->METLV[0].Pt() <= 25.0           ||
            ntuple->lLV[0].Pt() <= 30.0             ||
            TMath::Abs(ntuple->lLV[0].Eta()) >= 2.5 ||
            ntuple->jLV[0].Pt() <= 30.0             ||
            ntuple->jLV[1].Pt() <= 25.0             ||
            TMath::Abs(ntuple->jLV[0].Eta()) >= 2.4 ||
            TMath::Abs(ntuple->jLV[1].Eta()) >= 2.4)
            return false;
         if (proc->name.Contains("QCD") && (ntuple->lLV[0].lpfIso < 0.3 || ntuple->lLV[0].lpfIso > 0.7))
            return false;
      }
      if(ntuple->lLV[0].leptonCat == DEFS::muon){
         if(ntuple->METLV[0].Pt() <= 25.0           ||
            ntuple->lLV[0].Pt() <= 25.0             ||
            TMath::Abs(ntuple->lLV[0].Eta()) >= 2.1 ||
            ntuple->jLV[0].Pt() <= 30.0             ||
            ntuple->jLV[1].Pt() <= 25.0             ||
            TMath::Abs(ntuple->jLV[0].Eta()) >= 2.4 ||
            TMath::Abs(ntuple->jLV[1].Eta()) >= 2.4)
            return false;
         if (proc->name.Contains("QCD") && (ntuple->lLV[0].lpfIso < 0.3 || ntuple->lLV[0].lpfIso > 0.7))
            return false;
      }
   }
   else if (controlRegion == DEFS::event){
      // Keep only specific events
      if(ntuple->event!=2663602)
         return false;
   }
   else if (controlRegion == DEFS::Diboson){
      // Diboson analysis cuts
      TLorentzVector mt(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),0,ntuple->lLV[0].Et()+ntuple->METLV[0].Et());
      double wmt = sqrt(pow(ntuple->lLV[0].Et()+ntuple->METLV[0].Pt(), 2) - 
                        pow(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(), 2) - 
                        pow(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(), 2));

      if (ntuple->lLV[0].leptonCat == DEFS::muon) {
         if ((ntuple->lLV[0].Pt()) <= 25.0                                  ||
             (ntuple->METLV[0].Pt()) <= 25.0                                ||
             (ntuple->jLV[0].Pt()) <= 40.0                                  ||
             (ntuple->jLV[1].Pt()) <= 35.0                                  ||
             (ntuple->jLV[2].Pt()) > 30.0                                   ||
             (TMath::Abs(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta())) >= 1.5 || 
             ((ntuple->jLV[0]+ntuple->jLV[1]).Pt()) <= 70.0                 ||
             (TMath::Abs(ntuple->METLV[0].DeltaPhi(ntuple->jLV[0]))) <= 0.4 ||
             //(mt.M()) <= 30.0                                               ||
             (ntuple->lLV[0]+ntuple->METLV[0]).Pt() >= 200                  ||
             ntuple->lLV.size() > 1                                         ||
             wmt <= 30.0                                                     ) {
            return false;
         }
      }
      if (ntuple->lLV[0].leptonCat == DEFS::electron) {
         if ((ntuple->lLV[0].Pt()) <= 30.0                                  ||
             (ntuple->METLV[0].Pt()) <= 25.0                                ||
             (ntuple->jLV[0].Pt()) <= 40.0                                  ||
             (ntuple->jLV[1].Pt()) <= 35.0                                  ||
             (ntuple->jLV[2].Pt()) > 30.0                                   ||
             (TMath::Abs(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta())) >= 1.5 || 
             ((ntuple->jLV[0]+ntuple->jLV[1]).Pt()) <= 70.0                 ||
             (TMath::Abs(ntuple->METLV[0].DeltaPhi(ntuple->jLV[0]))) <= 0.4 ||
             //(mt.M()) <= 30.0                                               ||
             (ntuple->lLV[0]+ntuple->METLV[0]).Pt() >= 200                  ||
             ntuple->lLV.size() > 1                                         ||                       
             wmt <= 30.0                                                     ) {
            return false;
         }
      }
   }
   else if (controlRegion == DEFS::all)
      return false;

   return true;  

}// eventPassCuts

double UserFunctions::weightFunc(EventNtuple* ntuple, MicroNtuple* mnt, const PhysicsProcess* proc)
{
   if(UserFunctions::doBenchmarks)
      UserFunctions::func_benchmark->Start("weightFunc");

   TString auxName = proc->name;
   auxName.ToUpper();

   
   double weight = 1.0;

   if (doMetPhiWeight){
     
      if(auxName.Contains("DATA") || auxName.Contains("QCD")) {
        
         // TEST OF MET PHI WEIGHTING
         double metphi = ntuple->METLV[0].Phi();
         int bin     = metPhiWeight->FindBin(metphi);
         weight *= metPhiWeight->GetBinContent(bin);
	  	
      } //Data or QCD
	 
   } // doMetPhiWeight
  
   // Pileup reweighting
   if (doPUreweight){
      if(!auxName.Contains("DATA") && auxName.Contains("QCD")) //&& auxName.Contains("ENRICHED"))
         weight *= puweight->getWeight(ntuple->vLV[0].npv);
      else if(!auxName.Contains("DATA") && !auxName.Contains("QCD"))
         weight *= puweight->getWeight(ntuple->vLV[0].tnpus[1]);
   }

   //CSV reweighting
   if (doCSVreweight) {
      if (!auxName.Contains("DATA") && !auxName.Contains("QCD")) {
         if(doCSVsys.CompareTo("up",TString::kIgnoreCase)==0)
            weight *= TMath::Power(csvweight->getWeight(ntuple),2);
         else if(doCSVsys.CompareTo("down",TString::kIgnoreCase)==0)
            weight *= 1.0;
         else
            weight *= csvweight->getWeight(ntuple);
      }
   }
   
   //TTbar reweighting
   if (doTTbarreweight) {
      if(auxName.Contains("TTBAR")) {
         if(doTTbarsys.CompareTo("up",TString::kIgnoreCase)==0)
            weight *= TMath::Power(ttbarweight->getWeight(ntuple,proc),2);
         else if(doTTbarsys.CompareTo("down",TString::kIgnoreCase)==0)
            weight *= 1.0;
         else
            weight *= ttbarweight->getWeight(ntuple,proc);
         //weight *= ttbarWeight(ntuple,proc);
      }
   }

   // WJets weighting (specific to fix shape issues)
   if (WJweight && auxName.Contains("WJETS")){

      // find delta phi between lepton and leading jet and 
      // between two leading jets
      double dpLJ = ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]);
      double dpJJ = ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]);
      int bin     = WJetsWeightFunc->FindBin(dpLJ,dpJJ);
      weight     *= WJetsWeightFunc->GetBinContent(bin);
      
   } // wjets

   // QCD weighting (eta weighting)
   if (QCDweight && auxName.Contains("QCD")){
      
      // find lepton eta and apply weight
      double leptonEta = ntuple->lLV[0].Eta();
      //cout << "\tQCDWeight = " << QCDWeightFunc->GetBinContent(QCDWeightFunc->FindBin(fabs(leptonEta))) << " for eta = " << fabs(leptonEta) << endl;
      weight *= QCDWeightFunc->GetBinContent(QCDWeightFunc->FindBin(fabs(leptonEta)));
   }

   if (CosThetaLweight && auxName.Contains("WJETS")) {

      // find the CosTheta_l value and apply the weight
      if(mnt!=0) {
         if(CosThetaLweightSys.CompareTo("up",TString::kIgnoreCase)==0)
            weight *= TMath::Power(CosThetaLWeightFunc->GetBinContent(CosThetaLWeightFunc->FindBin(mnt->CosTheta_l)),2);
         else if(CosThetaLweightSys.CompareTo("down",TString::kIgnoreCase)==0)
            weight *= 1.0;
         else
            weight *= CosThetaLWeightFunc->GetBinContent(CosThetaLWeightFunc->FindBin(mnt->CosTheta_l));
      }
      else {
         Float_t Cos_dPhiWW = 0;
         Float_t Cos_dPhiWH = 0;
         Float_t CosTheta_l = 0;
         Float_t CosTheta_j = 0;
         Float_t CosTheta_WH = 0;
         Float_t JacksonAngle = 0;
         ntuple->getAngularVariables(Cos_dPhiWW,Cos_dPhiWH,CosTheta_l,
                                     CosTheta_j,CosTheta_WH,JacksonAngle);
         if(CosThetaLweightSys.CompareTo("up",TString::kIgnoreCase)==0)
            weight *= TMath::Power(CosThetaLWeightFunc->GetBinContent(CosThetaLWeightFunc->FindBin(CosTheta_l)),2);
         else if(CosThetaLweightSys.CompareTo("down",TString::kIgnoreCase)==0)
            weight*= 1.0;
         else
            weight *= CosThetaLWeightFunc->GetBinContent(CosThetaLWeightFunc->FindBin(CosTheta_l));
      }
   }

   if (WPtWeight && auxName.Contains("WJETS")) {
      weight *= WPtWeightFunc->GetBinContent(WPtWeightFunc->FindBin((ntuple->lLV[0]+ntuple->METLV[0]).Pt()));
   }

   if (PtWeight && !auxName.Contains("DATA") && auxName.Contains("QCD")) {
      weight *= PtWeightFunc->GetBinContent(PtWeightFunc->FindBin(ntuple->lLV[0].Pt()));
   }

   if (ZllPtWeight && auxName.Contains("ZJETSTOLL")) {
      weight *= ZllPtWeightFunc->GetBinContent(ZllPtWeightFunc->FindBin(ntuple->lLV[0].Pt()));
   }

   //if(doBCJetWeight && (auxName.Contains("WJETS") || auxName.Contains("ZJETS"))) {
   if(doBCJetWeight && auxName.Contains("WJETS")) {
      bool containsB = ntuple->containsParticle(5,true,25);
      bool containsC = ntuple->containsParticle(4,true,25);
      if(containsB || containsC) {
         weight *= 2.0;
      }
   }
   
   if(UserFunctions::doBenchmarks)
      UserFunctions::func_benchmark->Stop("weightFunc");

   return weight;

} // weightFunc

void UserFunctions::initEventFunc(EventNtuple* ntuple, const PhysicsProcess* proc)
{
   if(UserFunctions::doBenchmarks) {
      UserFunctions::func_benchmark->Reset();
      UserFunctions::func_benchmark->Start("initEventFunc");
   }

   TString auxName = proc->name;
   auxName.ToUpper();

   // Correct MET Phi
   if (doMETPhiCorrection) {
      if(auxName.Contains("DATA") || auxName.Contains("QCD")) {
         if(verbose) cout << "BEFORE\tPx=" << ntuple->METLV[0].Px() << "\tPy=" << ntuple->METLV[0].Py() << endl;
         //ntuple->doMETPhiCorrection("pfMEtSysShiftCorrParameters_2012runAvsNvtx_data");
         ntuple->doMETPhiCorrection("pfMEtSysShiftCorrParameters_2012runAvsNvtx_TAMUWW_data");
         if(verbose) cout << "AFTER  \tPx=" << ntuple->METLV[0].Px() << "\tPy=" << ntuple->METLV[0].Py() << endl;
      }
      else {
         if(verbose) cout << "BEFORE\tPx=" << ntuple->METLV[0].Px() << "\tPy=" << ntuple->METLV[0].Py() << endl;
         //ntuple->doMETPhiCorrection("pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc");
         ntuple->doMETPhiCorrection("pfMEtSysShiftCorrParameters_2012runAvsNvtx_TAMUWW_mc");
         if(verbose) cout << "AFTER  \tPx=" << ntuple->METLV[0].Px() << "\tPy=" << ntuple->METLV[0].Py() << endl;
      }
   }

   if(auxName.Contains("DATA") || auxName.Contains("QCD"))
      return;

   // Performs JER
   if (doJER)
      ntuple->doJER();

   if(UserFunctions::doBenchmarks)
      UserFunctions::func_benchmark->Stop("initEventFunc");
} // initEventFunc

//This is run once per process before looping over events
void UserFunctions::processFunc(EventNtuple* ntuple, const PhysicsProcess* proc)
{
   if(UserFunctions::doBenchmarks)
      UserFunctions::func_benchmark->Start("processFunc");

   TString auxName = proc->name;
   auxName.ToUpper();

   if (doMetPhiWeight){
     
      // get the filename
      TString filename = "data_weight.root";

      TString hname = "metphi";

      metPhiWeight = (TH1D*) DefaultValues::getConfigTH1(filename, hname, "metPhiWeightClone");
 	 
   }// if doMetPhiWeight
  
   if (QCDweight && auxName.Contains("QCD")){
      // Initializes QCD Reweighting
      cout << "Initializing QCD Weighting:" << endl;

      // get the filename
      TString filename = "QCDWeight_";
      if(QCDweightSys.CompareTo("up",TString::kIgnoreCase)==0) {
         filename += DEFS::getLeptonCatString(UserFunctions::leptonCat)+"_up.root";
      }
      else if(QCDweightSys.CompareTo("down",TString::kIgnoreCase)==0) {
         filename += DEFS::getLeptonCatString(UserFunctions::leptonCat)+"_down.root";
      }
      else {
         filename += DEFS::getLeptonCatString(UserFunctions::leptonCat)+".root";
      }
    
      // get the histogram from the inside
      TString hname = "QCDWeight_";
      hname += DEFS::getLeptonCatString(UserFunctions::leptonCat);

      QCDWeightFunc = (TH1D*) DefaultValues::getConfigTH1(filename,hname,"QCDweightClone");
      cout << "\t" << auxName << ": Using " << hname << " from " << filename << endl;
   }// QCD weight

   if (CosThetaLweight && (auxName.Contains("WJETS") || auxName.Contains("ZJETSTOLL"))) {
      // Initializes CosThetaL Reweighting
      cout << "Initializing CosTheta_l Weighting:" << endl;

      // get the filename
      TString filename = "CosThetaLWeight_eq1tag_TTbarSub_changedNames.root";

      // get the histogram from the inside
      TString hname = "CosThetaL_WeightHist_";
      hname += DEFS::getJetBinString(UserFunctions::jetBin)+"_"+DEFS::getLeptonCatString(UserFunctions::leptonCat);

      CosThetaLWeightFunc = (TH1D*) DefaultValues::getConfigTH1(filename,hname,"CosThetaLweightClone");
      cout << "\t" << auxName << ": Using " << hname << " from " << filename << endl;
   }

   if (WPtWeight && auxName.Contains("WJETS")) {
      cout << "Initializing WPt Weighting:" << endl;
      TString filename = "WPtWeights.root";
      string jb = DEFS::getJetBinString(UserFunctions::jetBin);
      jb[0] = toupper(jb[0]);
      TString hname = Form("ratioDataToWJets_%s_%s_Ptlv",jb.c_str(),DEFS::getLeptonCatString(UserFunctions::leptonCat).c_str());
      WPtWeightFunc = (TH1D*) DefaultValues::getConfigTH1(filename,hname,"WPtweightClone");
      cout << "\t" << auxName << ": Using " << hname << " from " << filename << endl;
   }

   if (PtWeight && !auxName.Contains("DATA") && auxName.Contains("QCD")) {
      TString filename = "PtWeight_";
      filename += DEFS::getLeptonCatString(UserFunctions::leptonCat)+".root";
      TString hname;
      if (auxName.Contains("QCD")) {
        //hname = "weights/QCDWeight_";
         hname = "QCDWeight_";
      } 
      else{
        hname = "MCWeight_";
      }
      hname += DEFS::getLeptonCatString(UserFunctions::leptonCat);
      PtWeightFunc = (TH1D*) DefaultValues::getConfigTH1(filename,hname,hname+"_Clone");
   }

   if (ZllPtWeight && auxName.Contains("ZJETSTOLL")) {
      cout << "Initializing ZllPt Weighting:" << endl;
      TString filename = "ZllPtWeights.root";
      TString hname = "ZllWeight_";
      hname += DEFS::getJetBinString(UserFunctions::jetBin)+"_"+DEFS::getLeptonCatString(UserFunctions::leptonCat);
      ZllPtWeightFunc = (TH1D*) DefaultValues::getConfigTH1(filename,hname,"ZllPtWeightClone");
      cout << "\t" << auxName << ": Using " << hname << " from " << filename << endl;
   }

   if (WJweight){
      if(auxName.Contains("WJETS")){

         // get the filename
         TString filename = "wjweight_";
         filename += DEFS::getLeptonCatString(UserFunctions::leptonCat)+".root";

         TString hname = "DeltaPhi_LJ1_vs_J1J2_";
         hname += DEFS::getLeptonCatString(UserFunctions::leptonCat);

         WJetsWeightFunc = (TH2D*) DefaultValues::getConfigTH2(filename, hname, "WJweightClone");
       
         /* In case we wanted to use a TF2
            wJetsWeightTF2 = new TF2("f2","[1]+( [2]*(x*cos([0])-y*sin([0])) + [3]*(x*cos([0])-y*sin([0]))**2 )" 
            "+( [4]*(x*sin([0])+y*cos([0])) + [5]*(x*sin([0])+y*cos([0]))**2 )",-10,10,-10,10);
            wJetsWeightTF2->SetParameters(-5.17471e-01, 7.79407e-01, 1.84199e-02, 1.10139e-01, 1.28591e-03, 2.18312e-02);  //redchi2 = 5.11
         */
      }
   }
  
   auxName.ToUpper();
   if(auxName.Contains("DATA"))// || auxName.Contains("QCD"))
      return;

   // Initializes PU Reweighting
   cout << "Initializing PU Reweighting:" << endl;
   //string dataname = DefaultValues::getConfigPath()+"pileup12_noTrig.root";
   string dataname;
   if(pileupSystematic.CompareTo("nominal")==0)
      dataname = DefaultValues::getConfigPath()+"pileup12_noTrig_minBiasXsec69400_coarseBinning_withAdditionalNPVHist.root";
   //PU Up
   else if(pileupSystematic.CompareTo("up")==0)
      dataname = DefaultValues::getConfigPath()+"pileup12_noTrig_minBiasXsec74258_coarseBinning_withAdditionalNPVHist.root";
   //PUDown
   else if(pileupSystematic.CompareTo("down")==0)
      dataname = DefaultValues::getConfigPath()+"pileup12_noTrig_minBiasXsec64542_coarseBinning_withAdditionalNPVHist.root";
   string MCname   = DefaultValues::getConfigPath()+"TPUDistributions.root";
   if(auxName.Contains("QCD") && (auxName.Contains("ElFULL") || auxName.Contains("ElEnriched"))) {
      cout << "\t" << auxName << ": Using pileup_SingleEl from " << dataname << " and " << TString("TPUDist_")+proc->name << " from " << MCname << endl;
      puweight = new PUreweight(dataname,MCname,"pileup_SingleEl",
                                string(TString("TPUDist_")+proc->name),
                                make_pair(1,10));
   }
   else if(auxName.Contains("QCD") && (auxName.Contains("MuFULL") || auxName.Contains("MuEnriched"))) {
      cout << "\t" << auxName << ": Using pileup_SingleMu from " << dataname << " and " << TString("TPUDist_")+proc->name << " from " << MCname << endl;
      puweight = new PUreweight(dataname,MCname,"pileup_SingleMu",
                                string(TString("TPUDist_")+proc->name),
                                make_pair(1,10));
   }
   else if(auxName.Contains("WH125_HTOBB")) {
      MCname = "root://cmseos.fnal.gov//store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/WH125_HToBB.root";
      cout << "\t" << auxName << ": Using pileup_noTrig from " << dataname << " and " << "TPUDist" << " from " << MCname << endl;
      puweight = new PUreweight(dataname,MCname,"pileup_noTrig",
                                string("PS/TPUDist"),
                                make_pair(1,10));      
   }
   else if(auxName.Contains("ZJETSTOLL_M10TO50")) {
      MCname = "root://cmseos.fnal.gov//store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/DYJetsToLL_M-10To50.root";
      cout << "\t" << auxName << ": Using pileup_noTrig from " << dataname << " and " << "TPUDist" << " from " << MCname << endl;
      puweight = new PUreweight(dataname,MCname,"pileup_noTrig",
                                string("PS/TPUDist"),
                                make_pair(1,10));      
   }
   else if(auxName.Contains("ZJETSTOLL_M50")) {
      MCname = "root://cmseos.fnal.gov//store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/DYJetsToLL_M-50.root";
      cout << "\t" << auxName << ": Using pileup_noTrig from " << dataname << " and " << "TPUDist" << " from " << MCname << endl;
      puweight = new PUreweight(dataname,MCname,"pileup_noTrig",
                                string("PS/TPUDist"),
                                make_pair(1,10));      
   }
   else {
      cout << "\t" << auxName << ": Using pileup_noTrig from " << dataname << " and " << TString("TPUDist_")+proc->name << " from " << MCname << endl;
      puweight = new PUreweight(dataname,MCname,"pileup_noTrig",
                                string(TString("TPUDist_")+proc->name),
                                make_pair(1,10));
   }

   if(auxName.Contains("QCD"))// && auxName.Contains("FULL"))
      return;

   if (doTTbarreweight && auxName.Contains("TTBAR")){
      ttbarweight = new TTbarreweight();
   }

   csvweight = new CSVreweight();

   if(UserFunctions::doBenchmarks)
      UserFunctions::func_benchmark->Stop("processFunc");

} // processFunc

template <class T, class U>
std::string UserFunctions::concatString(const T& obj1, const U& obj2)
{
   std::ostringstream output;
   output << obj1 << obj2;
   return output.str();
}//concatString

////////////////////////////////////////////////////////////////////////////////
//  Local Functions
////////////////////////////////////////////////////////////////////////////////

///  fills the histograms and controls the output canvas and file for the rest of the program
void doPlotter(MapOfPlots & plots, vector<PhysicsProcess*> procs, bool doJER, bool doPUrewt,
               bool doCSVrewt, bool doTTbarrewt, bool doFNAL, int maxEvts, bool WJweight,
               TString MVAWeightDir, vector<TString> MVAMethods, vector<TString> MVAVars,
               vector<TString> MVASpec, bool verbose);

// write the Canvases and plots to output files 
void writePlotsToFile(TString histoFileName, TString canvasFileName,
                      MapOfPlots & plots, vector<PhysicsProcess*> procs);

/// returns a map containing all of the plots that will be made for each process and their specific attributes
PlotFiller::MapOfPlots getPlots(DEFS::LeptonCat leptonCat, bool norm_data);

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
int main(int argc,char**argv) {
 
   // evaluate command-line / configuration file options
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;

   string           lepCat               = cl.getValue<string>    ("lep",                   "both");
   UserFunctions::leptonCat              = DEFS::getLeptonCat     (lepCat);
   string           ntype                = cl.getValue<string>    ("ntype",          "EventNtuple");
   UserFunctions::ntupleType             = DEFS::getNtupleType    (ntype);
   string           jBin                 = cl.getValue<string>    ("jBin",                 "jets2");      
   UserFunctions::jetBin                 = DEFS::getJetBin        (jBin);
   string           cutRegion            = cl.getValue<string>    ("cutRegion",              "all");
   UserFunctions::controlRegion          = DEFS::getControlRegion (cutRegion);
   string           tcat                 = cl.getValue<string>    ("tcat",                "pretag");
   UserFunctions::tagCat                 = DEFS::getTagCat(tcat);
   UserFunctions::elIsoMin               = cl.getValue<double>    ("elIsoMin",                 0.3);
   UserFunctions::elIsoMax               = cl.getValue<double>    ("elIsoMax",                 0.7);
   UserFunctions::muIsoMin               = cl.getValue<double>    ("muIsoMin",                 0.3);
   UserFunctions::muIsoMax               = cl.getValue<double>    ("muIsoMax",                 2.0);
   UserFunctions::doJER                  = cl.getValue<bool>      ("doJer",                  false);
   UserFunctions::doMETPhiCorrection     = cl.getValue<bool>      ("doMETPhi",               false);
   UserFunctions::doPUreweight           = cl.getValue<bool>      ("doPUrewt",                true);
   UserFunctions::pileupSystematic       = cl.getValue<TString>   ("pileupSystematic",   "nominal");
   UserFunctions::doCSVreweight          = cl.getValue<bool>      ("doCSVrewt",               true);
   UserFunctions::doCSVsys               = cl.getValue<TString>   ("doCSVsys",              "none");
   UserFunctions::doTTbarreweight        = cl.getValue<bool>      ("doTTbarrewt",             true);
   UserFunctions::doTTbarsys             = cl.getValue<TString>   ("doTTbarsys",            "none");
   UserFunctions::doFNAL                 = cl.getValue<bool>      ("doFNAL",                 false);
   UserFunctions::doMetPhiWeight         = cl.getValue<bool>      ("doMetPhiWeight",         false);
   UserFunctions::doBCJetWeight          = cl.getValue<bool>      ("doBCJetWeight",          false);
   UserFunctions::outDir                 = cl.getValue<string>    ("outDir",                   ".");
   UserFunctions::WJweight               = cl.getValue<bool>      ("WJweight",               false);
   UserFunctions::QCDweight              = cl.getValue<bool>      ("QCDweight",              false);
   UserFunctions::QCDweightSys           = cl.getValue<TString>   ("QCDweightSys",          "none");
   UserFunctions::CosThetaLweight        = cl.getValue<bool>      ("CosThetaLweight",        false);
   UserFunctions::CosThetaLweightSys     = cl.getValue<TString>   ("CosThetaLweightSys",    "none");
   UserFunctions::WPtWeight              = cl.getValue<bool>      ("WPtWeight",              false);
   UserFunctions::PtWeight               = cl.getValue<bool>      ("PtWeight",               false);
   UserFunctions::ZllPtWeight            = cl.getValue<bool>      ("ZllPtWeight",            false);
   //UserFunctions::PtShift
   UserFunctions::verbose                = cl.getValue<bool>      ("verbose",                false);
   UserFunctions::fillTMDF               = cl.getValue<bool>      ("fillTMDF",               false);
   UserFunctions::fill2D                 = cl.getValue<bool>      ("fill2D",                 false);
   UserFunctions::fillNonStandard        = cl.getValue<bool>      ("fillNonStandard",        false);
   UserFunctions::fillEPD                = cl.getValue<bool>      ("fillEPD",                false);
   UserFunctions::fillBumpCrossCheck     = cl.getValue<bool>      ("fillBumpCrossCheck",     false);
   UserFunctions::fillWJetsFlavor        = cl.getValue<bool>      ("fillWJetsFlavor",        false);
   UserFunctions::fillLimitTemplatesOnly = cl.getValue<bool>      ("fillLimitTemplatesOnly", false);
   UserFunctions::fillGen                = cl.getValue<bool>      ("fillGen",                false);
   bool             include_data         = cl.getValue<bool>      ("include_data",            true);
   bool             include_systematics  = cl.getValue<bool>      ("include_systematics",    false);
   bool             norm_data            = cl.getValue<bool>      ("norm_data",              false);
   int              maxEvts              = cl.getValue<int>       ("maxEvts",                    0);
   TString          MVAWeightDir         = cl.getValue<TString>   ("MVAWeightDir",              "");
   vector<TString>  MVAMethods           = cl.getVector<TString>  ("MVAMethods",                "");
   vector<TString>  MVAVars              = cl.getVector<TString>  ("MVAVars",                   "");
   vector<TString>  MVASpec              = cl.getVector<TString>  ("MVASpec",                   "");
   UserFunctions::limitBranches          = cl.getValue<int>       ("limitBranches",              0);
   bool             debug                = cl.getValue<bool>      ("debug",                  false);
   UserFunctions::doBenchmarks           = cl.getValue<bool>      ("diBenchmarks",           false);
   bool             bumpHunt             = cl.getValue<bool>      ("bumpHunt",               false);
   UserFunctions::signalTitle            = cl.getValue<TString>   ("signalTitle",     "H(125)->WW");
   bool             help                 = cl.getValue<bool>      ("help",                   false);

   if (help) {cl.print(); return 0;}
   if (!cl.check()) return 0;
   cl.print();

   // Trying to speed up the code
   gEnv->SetValue("TFile.AsyncPrefetching", 1);

   // Check that the leptcat actually exists
   if (UserFunctions::leptonCat == DEFS::none) {
      cout<<"plotter_x called with unknown lepton category "<<lepCat<<endl;
      return 1;
   }

   //Print error if limitBranches and fillGen are both on
   if(UserFunctions::limitBranches>=0 && UserFunctions::fillGen) {
      cout << "ERROR::plotter_x Can't fill the general gen level information when limitBranches>-1" << endl;
      return 2;
   }
  
   //Turn off EPD and bumpCrossCheck plots if the limitBranches flag is set
   if(UserFunctions::limitBranches > 0) {
      UserFunctions::fillEPD = false;
      UserFunctions::fillBumpCrossCheck = false;
      if(UserFunctions::limitBranches > 1)
         UserFunctions::fillWJetsFlavor = false;
   }

   //If no MVA methods, weightdir, or variables, don't create the MVADiscriminator plot
   if(MVAWeightDir.IsNull() || MVAMethods.size()==0 || MVAVars.size()==0) {
      UserFunctions::createMVADiscriminatorPlot = false;
   }
   else {
      UserFunctions::createMVADiscriminatorPlot = true;
   }

   // Check that the leptcat is not both as we are not ready for it. 
   // Problems if combining:
   // - we have to figure out how to include two lums.
   // - we have to figure out how to cut on one category but not the other.
   //if (UserFunctions::leptonCat == DEFS::both) {
   //   cout<<"plotter_x called with lepton category both. WE HAVE NOT FIXED THIS TO WORK YET."<<lepCat<<endl;
   //   return 1;
   //}

   // Tell ROOT to not print useless INFO
   gErrorIgnoreLevel = kWarning;

   TBenchmark* m_benchmark = new TBenchmark();
   m_benchmark->Reset();
   m_benchmark->Start("plotter_x");
   if(UserFunctions::doBenchmarks)
      UserFunctions::once_benchmark->Reset();

   cout << "Doing lepton category "<<lepCat << endl;

   // The vector containing all plots to be made
   MapOfPlots plots = getPlots(UserFunctions::leptonCat,norm_data);
   
   // The vector holding all processes.
   vector <PhysicsProcess*> procs;
   if(debug) {
      //procs.erase(procs.begin(),procs.begin()+17);
      //procs.erase(procs.begin()+1,procs.begin()+8);
      vector<DEFS::PhysicsProcess::Type> p;
      p.push_back(DEFS::PhysicsProcess::ZJetsToLL_M10To50);
      p.push_back(DEFS::PhysicsProcess::ZJetsToLL_M50);
      p.push_back(DEFS::PhysicsProcess::WJets);
      p.push_back(DEFS::PhysicsProcess::ZJets);
      procs = DefaultValues::getProcesses(p, UserFunctions::jetBin, UserFunctions::tagCat, true, UserFunctions::ntupleType);
   }
   else if(bumpHunt) {
      for(unsigned int i=0; i<procs.size(); i++) {
         if(!DEFS::PhysicsProcess::isHiggs(DEFS::PhysicsProcess::getProcessType(string(procs[i]->name)))) {
            cout << "Erasing process " << procs[i]->name << endl;
            procs.erase(procs.begin()+i);
            i--;
         }
      }
   }
   else if(UserFunctions::fillWJetsFlavor && !include_data) {
      vector<DEFS::PhysicsProcess::Type> p;
      p.push_back(DEFS::PhysicsProcess::WJets);
      procs = DefaultValues::getProcesses(p, UserFunctions::jetBin, UserFunctions::tagCat, true, UserFunctions::ntupleType);
   }
   else {
      procs = DefaultValues::getProcessesHiggs(UserFunctions::jetBin, UserFunctions::tagCat,
                                               include_data, include_systematics,true,
                                               UserFunctions::ntupleType, UserFunctions::leptonCat);
   }

   // Report Scale Factors
   cout << "Process " << setw(20) << "<name>" << " will be scaled by " << setw(15) << "<cross section>" << " * " << setw(12) << "<luminosity>"
        << " * " << setw(14) << "<scale factor>" << " * " << setw(17) << "<branching ratio>" << " / " << setw(15) << "<events in PAT>" << " = " << setw(13) << "<final value>" << endl; 
   for (unsigned p = 0; p< procs.size(); p++) {
      //Reset the scale factors for QCD/WJets if doing systematics
      if (UserFunctions::QCDweight && (procs[p]->name.Contains("QCD")||procs[p]->name.CompareTo("WJets")==0) && UserFunctions::QCDweightSys.CompareTo("none",TString::kIgnoreCase)!=0) {
         cout << "\tWARNING::Overridding the " << procs[p]->name << " scale factor." << endl
              << "\tMake sure this is what you want (usually okay if running systematic samples)." << endl;
         //procs[p]->scaleFactor[UserFunctions::leptonCat] = 1.0;
         if(UserFunctions::QCDweightSys.CompareTo("down",TString::kIgnoreCase)==0) {
            if(UserFunctions::leptonCat==DEFS::electron) {
               if(procs[p]->name.Contains("QCD")) {
                  procs[p]->scaleFactor[UserFunctions::leptonCat] = DefaultValues::getScaleFactor(procs[p]->name+"_QCDEtaWeightDown",UserFunctions::leptonCat);//0.508458;
               }
               else if(procs[p]->name.CompareTo("WJets")==0) {
                  procs[p]->scaleFactor[UserFunctions::leptonCat] = DefaultValues::getScaleFactor(procs[p]->name+"_QCDEtaWeightDown",UserFunctions::leptonCat);//1.04437;
               }
            }
            else if(UserFunctions::leptonCat==DEFS::muon) {
               if(procs[p]->name.Contains("QCD")) {
                  procs[p]->scaleFactor[UserFunctions::leptonCat] = DefaultValues::getScaleFactor(procs[p]->name+"_QCDEtaWeightDown",UserFunctions::leptonCat);//0.159824;
               }
               else if(procs[p]->name.CompareTo("WJets")==0) {
                  procs[p]->scaleFactor[UserFunctions::leptonCat] = DefaultValues::getScaleFactor(procs[p]->name+"_QCDEtaWeightDown",UserFunctions::leptonCat);//0.969496;
               }
            }
         }
         else if(UserFunctions::QCDweightSys.CompareTo("up",TString::kIgnoreCase)==0) {
            if(UserFunctions::leptonCat==DEFS::electron) {
               if(procs[p]->name.Contains("QCD")) {
                  procs[p]->scaleFactor[UserFunctions::leptonCat] = DefaultValues::getScaleFactor(procs[p]->name+"_QCDEtaWeightUp",UserFunctions::leptonCat);//0.27547;
               }
               else if(procs[p]->name.CompareTo("WJets")==0) {
                  procs[p]->scaleFactor[UserFunctions::leptonCat] = DefaultValues::getScaleFactor(procs[p]->name+"_QCDEtaWeightUp",UserFunctions::leptonCat);//1.0454;
               }
            }
            else if(UserFunctions::leptonCat==DEFS::muon) {
               if(procs[p]->name.Contains("QCD")) {
                  procs[p]->scaleFactor[UserFunctions::leptonCat] = DefaultValues::getScaleFactor(procs[p]->name+"_QCDEtaWeightUp",UserFunctions::leptonCat);//0.152692;
               }
               else if(procs[p]->name.CompareTo("WJets")==0) {
                  procs[p]->scaleFactor[UserFunctions::leptonCat] = DefaultValues::getScaleFactor(procs[p]->name+"_QCDEtaWeightUp",UserFunctions::leptonCat);//0.969527;
               }  
            }
         }
      }
      else if (!UserFunctions::QCDweight && (procs[p]->name.Contains("QCD")||procs[p]->name.CompareTo("WJets")==0) && UserFunctions::QCDweightSys.CompareTo("none",TString::kIgnoreCase)==0) {
         cout << "\tWARNING::Overridding the " << procs[p]->name << " scale factor." << endl
              << "\tThis is done because QCDWeight is set to false, not doing QCD weight systematics, and the process is either WJets or QCD." << endl;
         procs[p]->scaleFactor[UserFunctions::leptonCat] = 1.0;
      }
     cout<<"Process "<<setw(20)<<procs[p]->name<<" will be scaled by "<<setw(15)<<procs[p]->sigma[UserFunctions::leptonCat]<<" * "
         <<setw(12)<<procs[p]->intLum[UserFunctions::leptonCat]<<" * "<<setw(14)<<procs[p]->scaleFactor[UserFunctions::leptonCat]
         <<" * "<<setw(17)<<procs[p]->branching_ratio[UserFunctions::leptonCat]<<" / "<<setw(15)<<procs[p]->initial_events[UserFunctions::leptonCat]
         <<" = "<<setw(13)<<procs[p]->getScaleFactor(UserFunctions::leptonCat)<<endl;
   }

   // Fill all the plots 
   doPlotter(plots, procs, UserFunctions::doJER, UserFunctions::doPUreweight, UserFunctions::doCSVreweight, UserFunctions::doTTbarreweight,
             UserFunctions::doFNAL, maxEvts, UserFunctions::WJweight, MVAWeightDir,
             MVAMethods, MVAVars, MVASpec, UserFunctions::verbose);

   // Write output to file. The processs are passed to obtain the canvases
   if(!UserFunctions::outDir.EndsWith("/")) UserFunctions::outDir+="/";
   if(!gSystem->OpenDirectory(UserFunctions::outDir)) gSystem->mkdir(UserFunctions::outDir);
   TString histoFileName = UserFunctions::outDir + "histos_" + lepCat + "_" + jBin + ".root";
   TString canvasFileName = UserFunctions::outDir + "canvases_" + lepCat + "_" + jBin + ".root";
   writePlotsToFile(histoFileName, canvasFileName, plots, procs);

   if(UserFunctions::doBenchmarks) {
      cout << endl << "plotter_x::once_benchmark" << endl;
      float rt = 0, ct = 0;
      vector<string> timers;
      timers.push_back("getPlotsForLeptonCat");
      timers.push_back("getPlots");
      timers.push_back("doPlotterSetup");
      timers.push_back("writePlotsToFile");
      DefaultValues::printSummary(UserFunctions::once_benchmark, 8, rt, ct, timers);
   }

   m_benchmark->Stop("plotter_x");
   cout << "plotter_x" << endl << "\tCPU time = " << m_benchmark->GetCpuTime("plotter_x") << " s" << endl
        << "\tReal time = " << m_benchmark->GetRealTime("plotter_x") << " s" << endl;
   delete m_benchmark;

   return 0;

}//plotter()


////////////////////////////////////////////////////////////////////////////////
//  Implement Local Functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void writePlotsToFile(TString histoFileName, TString canvasFileName, 
                      MapOfPlots & plots,  vector<PhysicsProcess*> procs){
   if(UserFunctions::doBenchmarks)
      UserFunctions::once_benchmark->Start("writePlotsToFile");

   //Get the canvas and write them to file and as png and eps
   cout<<"Writing canvas to rootfile "<<canvasFileName<<endl;
   TFile * canOutFile = new TFile(canvasFileName,"RECREATE");
   for ( MapOfPlots::iterator p = plots.begin(); p != plots.end() ; p++)
      for ( map<string,  Plot * >::iterator p2 = p->second.begin(); 
            p2 != p->second.end() ; p2++){
            if(UserFunctions::fillLimitTemplatesOnly && 
               ((TString(p2->first).CompareTo("MVADiscriminator")!=0) && (TString(p2->first).CompareTo("MEBDT")!=0) &&
                (TString(p2->first).CompareTo("KinBDT")!=0) && !(TString(p2->first).Contains("KinMEBDT")) &&
                (TString(p2->first).CompareTo("MET")!=0)))
               continue;
         //try {
            TCanvas* can = ((FormattedPlot*) p2->second)->getCanvasTDR(procs);
            TString canName = UserFunctions::outDir+"/"+can->GetName();
            canName += "_"+DEFS::getLeptonCatString(UserFunctions::leptonCat);
            cout << "\tSaving canvas " << canName << " ... ";
            can->Write();
            can->SaveAs(canName+".png");
            TString class_name = ((FormattedPlot*)p2->second)->templateHisto->ClassName();
            if(!class_name.Contains("TH2")) {
               can->SaveAs(canName+".pdf");
            }
            //can->SaveAs(canName+".C");
            cout << "DONE" << endl << flush;
         //}
         //catch (const std::bad_alloc&) {
         //   cout << "\tERROR::writePlotsToFile::Caught a std::bad_alloc error for plot " << p2->first << endl;
         //}
      }
   canOutFile->Close();


   //Get the Histos and write them to file
   //NOTE: this needs to be done AFTER the writing of the canvases as the
   //      scaling of the histos is done at that stage. 
   cout<<"Writing histos to rootfile "<<histoFileName<<endl;
   bool firstH = true;
   for(MapOfPlots::iterator p = plots.begin(); p != plots.end() ; p++) {
      for(map<string,  Plot * >::iterator p2 = p->second.begin(); p2 != p->second.end() ; p2++) {
         if(firstH) {
            ((FormattedPlot*) p2->second)->saveHistogramsToFile(histoFileName,"RECREATE");
            firstH = false;
         }
         else {
            ((FormattedPlot*) p2->second)->saveHistogramsToFile(histoFileName,"UPDATE");
         }
      }// for histos inside plot
   }// for plots  

   if(UserFunctions::doBenchmarks)
      UserFunctions::once_benchmark->Stop("writePlotsToFile");
}//writePlotsToFile

//______________________________________________________________________________
void doPlotter(MapOfPlots & plots, vector<PhysicsProcess*> procs, bool doJER, bool doPUrewt,
               bool doCSVrewt, bool doTTbarrewt, bool doFNAL, int maxEvts, bool WJweight,
               TString MVAWeightDir, vector<TString> MVAMethods, vector<TString> MVAVars,
               vector<TString> MVASpec, bool verbose) {
   if(UserFunctions::doBenchmarks)
      UserFunctions::once_benchmark->Start("doPlotterSetup");

   PlotFiller pFill(plots, procs, &UserFunctions::fillPlots);
   pFill.setCutFunction(&UserFunctions::eventPassCuts);
   pFill.setWeightFunction(&UserFunctions::weightFunc);
   pFill.setProcessFunction(&UserFunctions::processFunc);
   pFill.setInitializeEventFunction(&UserFunctions::initEventFunc);
   pFill.setLimitBranches(UserFunctions::limitBranches);
   if (maxEvts>0)
      pFill.setMaximumEventsDEBUG(maxEvts); // TEST
   if (!MVAWeightDir.IsNull()) {
      pFill.setMVAWeightDir(MVAWeightDir);
      pFill.setMVAMethods(MVAMethods);
      pFill.setMVAVar(MVAVars);
      pFill.setMVASpec(MVASpec);
   }
   if(UserFunctions::doBenchmarks)
      UserFunctions::once_benchmark->Start("doPlotterSetup");

   pFill.run();

}//plotter(.., .., etc)

//______________________________________________________________________________
MapOfPlots getPlotsForLeptonCat(DEFS::LeptonCat leptonCat, bool norm_data){
   if(UserFunctions::doBenchmarks)
      UserFunctions::once_benchmark->Start("getPlotsForLeptonCat");

   MapOfPlots plots;

   FormattedPlot* a = new FormattedPlot;
   double pi = TMath::Pi();

   //Goes in the label and tells us whether we are looking at electrons or muons
   TString lepStr = "_"+DEFS::getLeptonCatString(leptonCat);

   // The overlay of a scaled signal. For signalName pick the groupingName 
   // of one of the processes. Or just "" if you don't want a signal overlayed.
   //TString signalName = "ggH+WH+qqH(125)";
   TString signalName = UserFunctions::signalTitle;
   double signalFactor = 1.0;
   if(leptonCat == DEFS::electron) signalFactor = 3900;
   else if(leptonCat == DEFS::muon) signalFactor = 3900;
   else signalFactor = 3900;

   //Double_t leptonptbinslow[9] = {20,25,30,35,40,50,70,100,1000};
   //Double_t leptonptbinshigh[10] = {20,50,55,60,65,70,80,100,120,1000};
   //Double_t jetptbinslow[9] = {20,25,30,35,40,50,70,100,1000};   
   //Double_t jetptbinshigh[10] = {20,50,80,100,110,120,130,140,160,1000};
   Double_t Mjjbinslow[11] = {20,40,50,60,70,80,90,100,200,300,1000};
   Double_t Mjjbinshigh[11] = {20,100,110,120,130,150,180,200,250,350,1000};
   Double_t Mtbinshigh[11] = {0,20,40,60,70,80,100,120,140,200,1000};
   Double_t METbinshigh[10] = {0,40,50,60,70,80,100,120,200,1000};
   Double_t DRlepjet1low[10] = {0.3,0.5,0.9,1.3,1.5,1.7,2.0,2.5,3.0,5.0};
   Double_t DRlepjet1high[11] = {0.3,2.0,2.25,2.5,2.75,3.0,3.25,3.5,4.0,4.5,5.0};
   Double_t etabins[83] = {
     -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664, 
     -3.489, -3.314, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172, 
     -2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305, 
     -1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522, 
     -0.435, -0.348, -0.261, -0.174, -0.087, 
     +0.000, 
     +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, 
     +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, 
     +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500, +2.650, 
     +2.853, +2.964, +3.139, +3.314, +3.489, +3.664, +3.839, +4.013, +4.191, 
     +4.363, +4.538, +4.716, +4.889, +5.191};
   
   Double_t absetabins[42] = {0.000, 
     +0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783, 
     +0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566, 
     +1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500, +2.650, 
     +2.853, +2.964, +3.139, +3.314, +3.489, +3.664, +3.839, +4.013, +4.191, 
     +4.363, +4.538, +4.716, +4.889, +5.191};
   
   vector<double> MET_binning;
   for(int i=0; i<=300; i+=5) {
      //if(i>300 && i%100!=0) continue;
      MET_binning.push_back(i);
   }
   MET_binning.push_back(1000);
     
   //goes in the label and tells us what cuts we are applying
   string cut = "_MET > 35, elPt > 30";
   TString cuts = TString(cut);
   TString name = "";

   a = new FormattedPlot;
   name = "MET";
   a->templateHisto = new TH1D(name + lepStr, name,MET_binning.size()-1,&MET_binning[0]);
   a->axisTitles.push_back("Missing E_{T} [GeV]");
   a->axisTitles.push_back("Number of Events / 5 GeV");
   a->range = make_pair(30.,150.);
   a->normToData = norm_data;
   a->stacked = true; a->leptonCat = leptonCat;
   a->overlaySignalName = signalName;
   a->overlaySignalFactor = signalFactor;
   plots[leptonCat][string(name)] = a;

   if(!UserFunctions::fillLimitTemplatesOnly) {
      if(UserFunctions::fillGen) {
         a = new FormattedPlot;
         name = "GenLeptPt";
         a->templateHisto = new TH1D(name + lepStr, name , MET_binning.size()-1,&MET_binning[0]);
         a->axisTitles.push_back("p_{T}^{gen lepton} [GeV]");
         a->axisTitles.push_back("Number of Events / 2 GeV");
         a->range = make_pair(20.,150.);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;
 
         a = new FormattedPlot;
         name = "GenLeptEta";
         a->templateHisto = new TH1D(name + lepStr, name, 82, etabins);
         a->axisTitles.push_back("#eta^{gen lepton}");
         a->axisTitles.push_back("Number of Events ");
         a->range = make_pair(-3.,3.);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;

         a = new FormattedPlot;
         name = "GenLeptPhi";
         a->templateHisto = new TH1D(name + lepStr, name,70,-3.5,3.5);
         a->axisTitles.push_back("#phi^{gen lepton} [Radians]");
         a->axisTitles.push_back("Number of Events / 0.1 Radians");
         a->range = make_pair(-3.5,3.5);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;

         a = new FormattedPlot;
         name = "GenWmT";
         a->templateHisto = new TH1D(name + lepStr, name,200,0,1000);
         a->axisTitles.push_back("M_{T}^{gen W} [GeV]");
         a->axisTitles.push_back("Number of Events");
         a->range = make_pair(0.0,150.0);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;

         a = new FormattedPlot;
         name = "GenVPt";
         a->templateHisto = new TH1D(name + lepStr, name,200,0,1000);
         a->axisTitles.push_back("p_{T}^{gen V} [GeV]");
         a->axisTitles.push_back("Number of Events");
         a->range = make_pair(0.0,150.0);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;

         a = new FormattedPlot;
         name = "GenVEta";
         a->templateHisto = new TH1D(name + lepStr, name,82,etabins);
         a->axisTitles.push_back("#eta^{gen V}");
         a->axisTitles.push_back("Number of Events");
         a->range = make_pair(-3.0,3.0);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;

         a = new FormattedPlot;
         name = "GenVPhi";
         a->templateHisto = new TH1D(name + lepStr, name,70,-3.5,3.5);
         a->axisTitles.push_back("#phi^{gen V} [Radians]");
         a->axisTitles.push_back("Number of Events / 0.1 Radians");
         a->range = make_pair(-3.5,3.5);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;
      }

      a = new FormattedPlot;
      name = "PUWeights";
      a->templateHisto = new TH1D(name + lepStr, name, 60,0,60);
      a->axisTitles.push_back("PU Weight");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,60);
      a->logxy = make_pair(false,false);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name  = "Mjj";
      a->templateHisto = new TH1D(name + lepStr, name,200,0.0,1000.0);
      a->axisTitles.push_back("M_{jj} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(40.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name  = "Mlv";
      a->templateHisto = new TH1D(name + lepStr, name,200,0.0,1000.0);
      a->axisTitles.push_back("M_{l#nu} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(40.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "MjjmWmT";
      a->templateHisto = new TH1D(name + lepStr, name,500,-250,250);
      a->axisTitles.push_back("M_{jj} - WmT [GeV]");
      a->axisTitles.push_back("Number of Events / GeV");
      a->range = make_pair(-100.,250.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   
      a = new FormattedPlot;
      name = "LeptPt";
      a->templateHisto = new TH1D(name + lepStr, name , MET_binning.size()-1,&MET_binning[0]);
      a->axisTitles.push_back("p_{T}^{lepton} [GeV]");
      a->axisTitles.push_back("Number of Events / 2 GeV");
      a->range = make_pair(20.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
 
      a = new FormattedPlot;
      name = "LeptEta";
      a->templateHisto = new TH1D(name + lepStr, name, 82, etabins);
      a->axisTitles.push_back("#eta^{lepton}");
      a->axisTitles.push_back("Number of Events ");
      a->range = make_pair(-3.,3.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "LeptPhi";
      a->templateHisto = new TH1D(name + lepStr, name,70,-3.5,3.5);
      a->axisTitles.push_back("#phi^{lepton} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(-3.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "LeptPFIso";
      a->templateHisto = new TH1D(name + lepStr, name,100,0,10);
      a->axisTitles.push_back("lepton PFIso");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,7);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "METPhi";
      a->templateHisto = new TH1D(name + lepStr, name, 62, -TMath::Pi(), TMath::Pi());
      a->axisTitles.push_back("Missing E_{T} Phi ");
      a->axisTitles.push_back("Number of Events" );
      a->range = make_pair(-3.5, 3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   
      a = new FormattedPlot;
      name = "Ptlv";
      a->templateHisto = new TH1D(name + lepStr, name,100,0,500);
      a->axisTitles.push_back("p_{T}^{W_{l#nu}} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(0.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "WmT";
      a->templateHisto = new TH1D(name + lepStr, name,100,0,500);
      a->axisTitles.push_back("M_{T}^{W} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(0.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "Jet1Pt";
      a->templateHisto = new TH1D(name + lepStr, name,100,0,300);
      a->axisTitles.push_back("p_{T}^{jet_{1}} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(20.,200.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "Jet1Eta";
      a->templateHisto = new TH1D(name + lepStr, name,50,-5,5);
      a->axisTitles.push_back("#eta^{jet_{1}} ");
      a->axisTitles.push_back("Number of Events ");
      a->range = make_pair(-3.,3.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "Jet1Phi";
      a->templateHisto = new TH1D(name + lepStr, name,70,-3.5,3.5);
      a->axisTitles.push_back("#phi^{jet_{1}} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(-3.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   
      a = new FormattedPlot;
      name = "Jet2Pt";
      a->templateHisto = new TH1D(name + lepStr, name,200,0,300);
      a->axisTitles.push_back("p_{T}^{jet_{2}} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(20.,100.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   
      a = new FormattedPlot;
      name = "Jet2Eta";
      a->templateHisto = new TH1D(name + lepStr, name,50,-5,5);
      a->axisTitles.push_back("#eta^{jet_{2}}");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-3.,3.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "Jet2Phi";
      a->templateHisto = new TH1D(name + lepStr, name,70,-3.5,3.5);
      a->axisTitles.push_back("#phi^{jet_{2}} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(-3.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   
      a = new FormattedPlot;
      name = "DeltaEtaJ1J2";
      a->templateHisto = new TH1D(name + lepStr, name,50,0,5);
      a->axisTitles.push_back("#eta^{jet_{1}} - #eta^{jet_{2}} ");
      a->axisTitles.push_back("Number of Events ");
      a->range = make_pair(0.,5.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   
      a = new FormattedPlot;
      name = "Ptjj";
      a->templateHisto = new TH1D(name + lepStr, name,100,0,300);
      a->axisTitles.push_back("p_{T}^{jj} [GeV]");
      a->axisTitles.push_back("Number of Events / 10 GeV");
      a->range = make_pair(0.,250.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "j1Pt_Mjj";
      a->templateHisto = new TH1D(name + lepStr, name,500,0,5);
      a->axisTitles.push_back("p_{T}^{jet_{1}}/M_{jj} [GeV]");
      a->axisTitles.push_back("Number of Events / 0.01 GeV");
      a->range = make_pair(0.,2.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   
      a = new FormattedPlot;
      name = "j2Pt_Mjj";
      a->templateHisto = new TH1D(name + lepStr, name,500,0,5);
      a->axisTitles.push_back("p_{T}^{jet_{2}}/M_{jj} [GeV]");
      a->axisTitles.push_back("Number of Events / 0.01 GeV");
      a->range = make_pair(0.,1.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "Mlvjj";
      a->templateHisto = new TH1D(name + lepStr, name,50,0,1000);
      a->axisTitles.push_back("M_{lvjj} [GeV]");
      a->axisTitles.push_back("Number of Events / 4 GeV");
      a->range = make_pair(50.,800.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "DeltaRJ1J2";
      a->templateHisto = new TH1D(name + lepStr, name,50,0,10);
      a->axisTitles.push_back("#DeltaR of J1 and J2 ");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0.,7.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   
      a = new FormattedPlot;
      name = "DeltaPhi_J1J2";
      a->templateHisto = new TH1D(name + lepStr, name,50,-10,10);
      a->axisTitles.push_back("#Delta #phi(J1J2)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   
      a = new FormattedPlot;
      name = "npv";
      a->templateHisto = new TH1D(name + lepStr, name,100,0,100);
      a->axisTitles.push_back("npv");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,40);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "jet1dRLep";
      a->templateHisto = new TH1D(name + lepStr, name,50,0,10);
      a->axisTitles.push_back("#DeltaR(Lepton,Jet1)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0.,7.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "jet2dRLep";
      a->templateHisto = new TH1D(name + lepStr, name,50,0,10);
      a->axisTitles.push_back("#DeltaR(Lepton,Jet2)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0.,7.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "jet3dRLep";
      a->templateHisto = new TH1D(name + lepStr, name,50,0,10);
      a->axisTitles.push_back("#DeltaR(Lepton,Jet3)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0.,7.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "leptonEtaCharge";
      a->templateHisto = new TH1D(name + lepStr, name,5,-2,2);
      a->axisTitles.push_back("Lepton Eta Charge");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-2,2);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;      

      a = new FormattedPlot;
      name = "ht";
      a->templateHisto = new TH1D(name + lepStr, name,40,0,800);
      a->axisTitles.push_back("ht");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,800);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a; 

      a = new FormattedPlot;
      name = "Ptlnujj";
      a->templateHisto = new TH1D(name + lepStr, name,25,0,150);
      a->axisTitles.push_back("p_{T}^{l#nujj}");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,150);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a; 

      a = new FormattedPlot;
      name = "dPhiMETJet";
      a->templateHisto = new TH1D(name + lepStr, name,62,-TMath::Pi(),TMath::Pi());
      a->axisTitles.push_back("#Delta#Phi(MET,Jet1)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-TMath::Pi(),TMath::Pi());
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a; 

      a = new FormattedPlot;
      name = "dPhiMETLep";
      a->templateHisto = new TH1D(name + lepStr, name,62,-TMath::Pi(),TMath::Pi());
      a->axisTitles.push_back("#Delta#Phi(MET,Lepton)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-TMath::Pi(),TMath::Pi());
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a; 

      a = new FormattedPlot;
      name = "dRlepjj";
      a->templateHisto = new TH1D(name + lepStr, name,50,0,10);
      a->axisTitles.push_back("#DeltaR(Lepton,Jet1Jet2)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0.,7.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "DeltaPhi_LJ1";
      a->templateHisto = new TH1D(name + lepStr, name,50,-10,10);
      a->axisTitles.push_back("#Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "DeltaPhi_METJ1";
      a->templateHisto = new TH1D(name + lepStr, name, 50,-10,10);
      a->axisTitles.push_back("#Delta #phi(MET,J1)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-pi,pi);
      a->logxy = make_pair(false,false);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "DeltaPhi_LMET";
      a->templateHisto = new TH1D(name + lepStr, name,50,-10,10);
      a->axisTitles.push_back("#Delta #phi(Lepton,MET)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "minDPhiLepJet";
      a->templateHisto = new TH1D(name + lepStr, name,62,-TMath::Pi(),TMath::Pi());
      a->axisTitles.push_back("min(#Delta#Phi(Lepton,Jet1))");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-TMath::Pi(),TMath::Pi());
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a; 

      a = new FormattedPlot;
      name = "CosTheta_l";
      a->templateHisto = new TH1D(name + lepStr, name,50,-1,1);
      a->axisTitles.push_back("Cos(#theta)_{l}");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-TMath::Pi(),TMath::Pi());
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a; 

      a = new FormattedPlot;
      name = "CosTheta_j";
      a->templateHisto = new TH1D(name + lepStr, name,50,-1,1);
      a->axisTitles.push_back("Cos(#theta)_{j}");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-TMath::Pi(),TMath::Pi());
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a; 

      a = new FormattedPlot;
      name = "CosTheta_WH";
      a->templateHisto = new TH1D(name + lepStr, name,50,-1,1);
      a->axisTitles.push_back("Cos(#theta)_{WH}");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-TMath::Pi(),TMath::Pi());
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      if(UserFunctions::ntupleType==DEFS::MicroNtuple || UserFunctions::ntupleType==DEFS::METree) {
         for (unsigned int tep=0; tep<MicroNtuple::nEventProb; tep++) {
            a = new FormattedPlot;
            name = TString(UserFunctions::concatString("tEventProb",tep));
            string xaxis = UserFunctions::concatString("log10(tEventProb[",tep) + "])";
            a->templateHisto = new TH1D(name + lepStr, name,180,-45,0);
            a->axisTitles.push_back(xaxis);
            a->axisTitles.push_back("Number of Events");
            if(UserFunctions::maxEventProbs[tep].first>0)
               a->range = make_pair(-13,0);
            else
               a->range = make_pair(-25,0);
            a->normToData = norm_data;
            a->stacked = true; a->leptonCat = leptonCat;
            a->overlaySignalName = signalName;
            a->overlaySignalFactor = signalFactor;
            plots[leptonCat][string(name)] = a;
         }

         a = new FormattedPlot;
         name = "MVAProbability";
         a->templateHisto = new TH1D(name + lepStr, name, 100,0,1);
         a->axisTitles.push_back("Probability");
         a->axisTitles.push_back("Number of Events");
         a->range = make_pair(0,1);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;
      }
   }

   if(UserFunctions::fillNonStandard) {
      a = new FormattedPlot;
      name = "SigmaIetaIeta";
      a->templateHisto = new TH1D(name + lepStr, name,200,0,0.1);
      a->axisTitles.push_back("lepton SigmaIetaIeta");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,0.05);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "MVATrig";
      a->templateHisto = new TH1D(name + lepStr, name,2000,0,1.0);
      a->axisTitles.push_back("lepton mvaTrig");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0.95,1.0);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "MVANonTrig";
      a->templateHisto = new TH1D(name + lepStr, name,2000,0,1.0);
      a->axisTitles.push_back("lepton mvaNonTrig");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0.8,1.0);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
      
      a = new FormattedPlot;
      name = "TTbarWeights";
      a->templateHisto = new TH1D(name + lepStr, name, 400,-2,2);
      a->axisTitles.push_back("TTbar Weight");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,1.3);
      a->logxy = make_pair(false,false);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "EJ1EJ2";
      a->templateHisto = new TH1D(name + lepStr, name,500,700,5000);
      a->axisTitles.push_back("E_{J1} * E_{J2}  [GeV]");
      a->axisTitles.push_back("Number of Events / 10 GeV");
      a->range = make_pair(0.,5000.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;


      a = new FormattedPlot;
      name  = "BetaJ1BetaJ2";
      a->templateHisto = new TH1D(name + lepStr, name,10000,0,10);
      a->axisTitles.push_back("#beta_{J1} * #beta_{J2} [GeV]");
      a->axisTitles.push_back("Number of Events / .01 GeV");
      a->range = make_pair(0.9,1.03);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "AngleJ1J2";
      a->templateHisto = new TH1D(name + lepStr, name,50,0,5);
      a->axisTitles.push_back("Angle between J1 and J2 [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(-0.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "jjlvPhi";
      a->templateHisto = new TH1D(name + lepStr, name,70,-3.5,3.5);
      a->axisTitles.push_back("#phi(jj) - #phi(lv)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-3.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "DeltaRLepMET";
      a->templateHisto = new TH1D(name + lepStr, name,50,0,10);
      a->axisTitles.push_back("#DeltaR of lep and MET ");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0.,7.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "WmT_Negative";
      a->templateHisto = new TH2D(name + lepStr, name ,15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Negative WmT");
      a->axisTitles.push_back("Number of Events ");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name  = "WmT_Positive";
      a->templateHisto = new TH2D(name + lepStr, name ,15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Positive WmT");
      a->axisTitles.push_back("Number of Events ");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "WmT_Subtracted";
      a->templateHisto = new TH2D(name + lepStr, name ,15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Subtracted WmT");
      a->axisTitles.push_back("Number of Events ");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "nJets_JetEta";
      a->templateHisto = new TH2D(name + lepStr, name,70,-3.5,3.5,20,0,20);
      a->axisTitles.push_back("#eta^{jet_{i}}");
      a->axisTitles.push_back("Number of Jets");
      a->range = make_pair(-TMath::Pi(),TMath::Pi());
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      for (unsigned int nj=0; nj<32; nj++) {
         a = new FormattedPlot;
         name = Form("JetEta_%uJets",nj);
         a->templateHisto = new TH1D(name + lepStr, name,70,-3.5,3.5);
         a->axisTitles.push_back("#eta^{jets} ");
         a->axisTitles.push_back("Number of Events ");
         a->range = make_pair(-pi,pi);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;
      }
   }
   
   if(UserFunctions::fill2D) {
      a = new FormattedPlot;
      name = "MET_vs_LeptEta";
      a->templateHisto = new TH2D(name + lepStr, name, 82, etabins, 1000,0,500);
      a->axisTitles.push_back("Lepton #eta");
      a->axisTitles.push_back("MET");
      a->axisTitles.push_back("Number of Events" );
      a->range = make_pair(-4,4);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "MET_vs_LeptPt";
      a->templateHisto = new TH2D(name + lepStr, name, 1000, 0, 500, 1000,0,500);
      a->axisTitles.push_back("p_{T}^{lepton} [GeV]");
      a->axisTitles.push_back("MET");
      a->axisTitles.push_back("Number of Events" );
      a->range = make_pair(20.0,150.0);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "MET_vs_AbsLeptEta";
      a->templateHisto = new TH2D(name + lepStr, name, 41, absetabins, 1000,0,500);
      a->axisTitles.push_back("Lepton #eta");
      a->axisTitles.push_back("MET");
      a->axisTitles.push_back("Number of Events" );
      a->range = make_pair(-4,4);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "MWjjVsMWlv";
      a->templateHisto = new TH2D(name + lepStr, name ,200,0,200,200,0,200);
      a->axisTitles.push_back("M_{W_{jj}}");
      a->axisTitles.push_back("M_{W_{l#nu}}");
      a->range = make_pair(0,200);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "DeltaPhi_LJ1_vs_J1J2_Positive";
      a->templateHisto = new TH2D(name + lepStr, name,15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Positive #Delta #phi(J1J2) vs. #Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "DeltaPhi_LJ1_vs_J1J2_Negative";
      a->templateHisto = new TH2D(name + lepStr, name ,15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Negative #Delta #phi(J1J2) vs. #Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "DeltaEta_LJ1_vs_J1J2_Positive";
      a->templateHisto = new TH2D(name + lepStr, name,40,-10,10,40,-10,10);
      a->axisTitles.push_back("Positive #Delta #eta(J1J2) vs. #Delta #eta(LJ1)");
      a->axisTitles.push_back("Number of Events / #eta unit");
      a->range = make_pair(-10,10);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "DeltaEta_LJ1_vs_J1J2_Negative";
      a->templateHisto = new TH2D(name + lepStr, name,40,-10,10,40,-10,10);
      a->axisTitles.push_back("Negative #Delta #eta(J1J2) vs. #Delta #eta(LJ1)");
      a->axisTitles.push_back("Number of Events / #eta unit");
      a->range = make_pair(-10,10);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "DeltaPhi_LJ1_vs_J1J2";
      a->templateHisto = new TH2D(name + lepStr, name, 15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("#Delta #phi(LJ1)");
      a->axisTitles.push_back("#Delta #phi(J1J2)");
      a->axisTitles.push_back("Number of Events" );
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name  = "DeltaPhi_LJ1_vs_J1J2_Subtracted";
      a->templateHisto = new TH2D(name + lepStr, name ,15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Subtracted #Delta #phi(J1J2) vs. #Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "DeltaEta_LJ1_vs_J1J2";
      a->templateHisto = new TH2D(name + lepStr, name,40,-10,10,40,-10,10);
      a->axisTitles.push_back("#Delta #eta(J1J2) vs. #Delta #eta(LJ1)");
      a->axisTitles.push_back("Number of Events / #eta unit");
      a->range = make_pair(-10,10);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "IsolationEnergyVsPt";
      a->templateHisto = new TH2D(name + lepStr, name,500,0,500,500,0,500);
      a->axisTitles.push_back("p_{T}^{Lepton} [GeV/c]");
      a->axisTitles.push_back("Isolation Energy [GeV]");
      a->range = make_pair(0,500);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "IsolationEnergyVsPtCoshEta";
      a->templateHisto = new TH2D(name + lepStr, name,1000,0,1000,500,0,500);
      a->axisTitles.push_back("p_{T} x cosh(#eta) [GeV]");
      a->axisTitles.push_back("Isolation Energy [GeV]");
      a->range = make_pair(0,1000);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   }

   if(UserFunctions::ntupleType==DEFS::MicroNtuple || UserFunctions::ntupleType==DEFS::METree) {
      if(UserFunctions::fillWJetsFlavor) {
         a = new FormattedPlot;
         name = "KinMEBDT_WpLight";
         a->templateHisto = new TH1D(name + lepStr, name,
                                     DEFS::getNBinsX(UserFunctions::jetBin,leptonCat,string("KinMEBDT")),
                                     &DEFS::getBinsX(UserFunctions::jetBin,leptonCat,string("KinMEBDT"))[0]);
         a->axisTitles.push_back("KinMEBDT");
         a->axisTitles.push_back("Number of Events");
         a->range = make_pair(-1.0,1.0);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;

         a = new FormattedPlot;
         name = "KinMEBDT_WpC";
         a->templateHisto = new TH1D(name + lepStr, name,
                                     DEFS::getNBinsX(UserFunctions::jetBin,leptonCat,string("KinMEBDT")),
                                     &DEFS::getBinsX(UserFunctions::jetBin,leptonCat,string("KinMEBDT"))[0]);
         a->axisTitles.push_back("KinMEBDT");
         a->axisTitles.push_back("Number of Events");
         a->range = make_pair(-1.0,1.0);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;

         a = new FormattedPlot;
         name = "KinMEBDT_WpB";
         a->templateHisto = new TH1D(name + lepStr, name,
                                     DEFS::getNBinsX(UserFunctions::jetBin,leptonCat,string("KinMEBDT")),
                                     &DEFS::getBinsX(UserFunctions::jetBin,leptonCat,string("KinMEBDT"))[0]);
         a->axisTitles.push_back("KinMEBDT");
         a->axisTitles.push_back("Number of Events");
         a->range = make_pair(-1.0,1.0);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;

         a = new FormattedPlot;
         name = "KinMEBDT_WpOther";
         a->templateHisto = new TH1D(name + lepStr, name,
                                     DEFS::getNBinsX(UserFunctions::jetBin,leptonCat,string("KinMEBDT")),
                                     &DEFS::getBinsX(UserFunctions::jetBin,leptonCat,string("KinMEBDT"))[0]);
         a->axisTitles.push_back("KinMEBDT");
         a->axisTitles.push_back("Number of Events");
         a->range = make_pair(-1.0,1.0);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;

         a = new FormattedPlot;
         name = "KinMEBDT_WpInclusive";
         a->templateHisto = new TH1D(name + lepStr, name,
                                     DEFS::getNBinsX(UserFunctions::jetBin,leptonCat,string("KinMEBDT")),
                                     &DEFS::getBinsX(UserFunctions::jetBin,leptonCat,string("KinMEBDT"))[0]);
         a->axisTitles.push_back("KinMEBDT");
         a->axisTitles.push_back("Number of Events");
         a->range = make_pair(-1.0,1.0);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = leptonCat;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[leptonCat][string(name)] = a;
      }

      a = new FormattedPlot;
      name = "KinBDT";
      a->templateHisto = new TH1D(name + lepStr, name,
                                  //100,-1,1);
                                  DEFS::getNBinsX(UserFunctions::jetBin,leptonCat,string(name)),
                                  &DEFS::getBinsX(UserFunctions::jetBin,leptonCat,string(name))[0]);
      a->axisTitles.push_back("KinBDT");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-1.0,1.0);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "MEBDT";
      a->templateHisto = new TH1D(name + lepStr, name,
                                  //100,-1,1);
                                  DEFS::getNBinsX(UserFunctions::jetBin,leptonCat,string(name)),
                                  &DEFS::getBinsX(UserFunctions::jetBin,leptonCat,string(name))[0]);
      a->axisTitles.push_back("MEBDT");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-1.0,1.0);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "KinMEBDT";
      a->templateHisto = new TH1D(name + lepStr, name,
                                  //100,-1,1);
                                  DEFS::getNBinsX(UserFunctions::jetBin,leptonCat,string(name)),
                                  &DEFS::getBinsX(UserFunctions::jetBin,leptonCat,string(name))[0]);
      a->axisTitles.push_back("KinMEBDT");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-1.0,1.0);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   }

   if(UserFunctions::createMVADiscriminatorPlot) {
      a = new FormattedPlot;
      name = "MVADiscriminator";
      a->templateHisto = new TH1D(name + lepStr, name,
                                  DEFS::getNBinsX(UserFunctions::jetBin,leptonCat,string(name)),
                                  &DEFS::getBinsX(UserFunctions::jetBin,leptonCat,string(name))[0]);
      a->axisTitles.push_back("MVADiscriminator");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-1.0,1.0);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   }

   if(UserFunctions::fillBumpCrossCheck) {
      a = new FormattedPlot;
      name = "MVADiscriminator_fineBinning";
      a->templateHisto = new TH1D(name + lepStr, name, 100, -1.0, 1.0);
      a->axisTitles.push_back("MVADiscriminator");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-1.0,1.0);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
/*
      a = new FormattedPlot;
      name = "MVADiscriminator_fineBinning_0GenMatch";
      a->templateHisto = new TH1D(name + lepStr, name, 100, -1.0, 1.0);
      a->axisTitles.push_back("MVADiscriminator");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-1.0,1.0);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "MVADiscriminator_fineBinning_1GenMatch";
      a->templateHisto = new TH1D(name + lepStr, name, 100, -1.0, 1.0);
      a->axisTitles.push_back("MVADiscriminator");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-1.0,1.0);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "MVADiscriminator_fineBinning_2GenMatch";
      a->templateHisto = new TH1D(name + lepStr, name, 100, -1.0, 1.0);
      a->axisTitles.push_back("MVADiscriminator");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-1.0,1.0);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
*/
      a = new FormattedPlot;
      name = "MEBDT_vs_KinBDT";
      a->templateHisto = new TH2D(name + lepStr, name,
                                  DEFS::getNBinsX(UserFunctions::jetBin,leptonCat,string("KinBDT")),
                                  &DEFS::getBinsX(UserFunctions::jetBin,leptonCat,string("KinBDT"))[0],
                                  DEFS::getNBinsX(UserFunctions::jetBin,leptonCat,string("MEBDT")),
                                  &DEFS::getBinsX(UserFunctions::jetBin,leptonCat,string("MEBDT"))[0]);
      a->axisTitles.push_back("KinBDT");
      a->axisTitles.push_back("MEBDT");
      a->axisTitles.push_back("Number of Events");
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;


   }

   if(UserFunctions::fillEPD) {
      a = new FormattedPlot;
      name = "epdPretagWWandWZ";
      a->templateHisto = new TH1D(name + lepStr, name, 400,-10.0,10.0);
      a->axisTitles.push_back("-log10(epdPretagWWandWZ)");
      a->axisTitles.push_back("Number of Events");
      //a->range = make_pair(-7.0,0.0);
      a->range = make_pair(0.0,7.0);
      a->logxy = make_pair(false,false);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "epdPretagHiggs125";
      a->templateHisto = new TH1D(name + lepStr, name, 400,-10.0,10.0);
      a->axisTitles.push_back("-log10(epdPretagHiggs125)");
      a->axisTitles.push_back("Number of Events");
      //a->range = make_pair(-7.0,0.0);
      a->range = make_pair(0.0,7.0);
      a->logxy = make_pair(false,false);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;

      a = new FormattedPlot;
      name = "epdPretagWWandWZ_RE";
      a->templateHisto = new TH1D(name + lepStr, name, 22,0,1);
      a->axisTitles.push_back("EPD(WW,WZ)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,1);
      a->logxy = make_pair(false,false);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   }

   if (UserFunctions::fillTMDF){
/*
  a = new FormattedPlot;
  name = "lpt_leta_j1pt_j1eta_j2pt_j2eta";
  a->templateHisto = new TProfileMDF(name + lepStr, name);
  ((TProfileMDF*)a->templateHisto)->AddAxis("lpt",8,leptonptbinslow);
  ((TProfileMDF*)a->templateHisto)->AddAxis("leta",10,0,2.5);
  ((TProfileMDF*)a->templateHisto)->AddAxis("j1pt",8,jetptbinslow);
  ((TProfileMDF*)a->templateHisto)->AddAxis("j1eta",10,0,2.5);
  ((TProfileMDF*)a->templateHisto)->AddAxis("j2pt",8,jetptbinslow);
  ((TProfileMDF*)a->templateHisto)->AddAxis("j2eta",10,0,2.5);
  ((TProfileMDF*)a->templateHisto)->Sumw2();
  a->normToData = norm_data;
  a->stacked = true; a->leptonCat = leptonCat;
  a->overlaySignalName = signalName;
  a->overlaySignalFactor = signalFactor;
  plots[leptonCat][string(name)] = a;
     

  name = "lpt_lpt_j1pt_j1pt_j2pt_j2pt";
  a->templateHisto = new TProfileMDF(name + lepStr, name);
  ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{lepton} [GeV]",8,leptonptbinslow);
  ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{lepton} [GeV]",9,leptonptbinshigh);
  ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{jet_{1}} [GeV]",8,jetptbinslow);
  ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{jet_{1}} [GeV]",9,jetptbinshigh);
  ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{jet_{2}} [GeV]",8,jetptbinslow);
  ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{jet_{2}} [GeV]",9,jetptbinshigh);
  ((TProfileMDF*)a->templateHisto)->Sumw2();
  a->normToData = norm_data;
  a->stacked = true; a->leptonCat = leptonCat;
  a->overlaySignalName = signalName;
  a->overlaySignalFactor = signalFactor;
  plots[leptonCat][string(name)] = a;
*/
      name = "Mjj_Mjj_Mt_MET_DeltaR_DeltaR_";
      a->templateHisto = new TProfileMDF(name + lepStr, name);
      ((TProfileMDF*)a->templateHisto)->AddAxis("M_{jj} [GeV]",10,Mjjbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("M_{jj} [GeV]",10,Mjjbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("M_{T}^{W} [GeV]",10,Mtbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("Missing E_{T} [GeV]",9,METbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("#DeltaR(#mu,jet1) [Radians]",9,DRlepjet1low);
      ((TProfileMDF*)a->templateHisto)->AddAxis("#DeltaR(#mu,jet1) [Radians]",10,DRlepjet1high);
      ((TProfileMDF*)a->templateHisto)->Sumw2();
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = leptonCat;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[leptonCat][string(name)] = a;
   }

   if(UserFunctions::doBenchmarks)
      UserFunctions::once_benchmark->Stop("getPlotsForLeptonCat");

   // return all the plots to be made
   return plots;

}//getPlotsForLeptonCat


//______________________________________________________________________________
MapOfPlots getPlots(DEFS::LeptonCat leptonCat, bool norm_data){
   if(UserFunctions::doBenchmarks)
      UserFunctions::once_benchmark->Start("getPlots");

   MapOfPlots plots;

   if (leptonCat == DEFS::electron || leptonCat == DEFS::both){
      MapOfPlots plots_ele = getPlotsForLeptonCat(DEFS::electron, norm_data);
      plots.insert(plots_ele.begin(),plots_ele.end());
   }
   if (leptonCat == DEFS::muon || leptonCat == DEFS::both){
      MapOfPlots plots_muo = getPlotsForLeptonCat(DEFS::muon, norm_data);
      plots.insert(plots_muo.begin(),plots_muo.end());
   }
   if(UserFunctions::doBenchmarks)
      UserFunctions::once_benchmark->Stop("getPlots");

   return plots;
}
