//Our libraries
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
#include "JetMETAnalysis/JetUtilities/interface/TProfileMDF.h"

// ROOT libraries
#include "TROOT.h"
#include "TSystem.h"
#include "TBenchmark.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2D.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TMath.h"

// C++ libraries
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <utility>

using namespace std;
using DEFS::LeptonCat;

namespace UserFunctions
{
   TString outDir;
   TString cutRegion;
   DEFS::LeptonCat leptonCat;
   PUreweight* puweight;
   bool doJER;
   bool doPUreweight;
   bool doFNAL;
   bool WJweight;
   TH2D* wJetsWeight= 0;
   bool verbose; // adds or takes away cout statements when running
   
   // Used if we run both leptons
   bool doElectron;
   bool doMuon;

////////////////////////////////////////////////////////////////////////////////
//  User Functions
////////////////////////////////////////////////////////////////////////////////

   // Is run once for each process before events are cut (initialize)
   void initEventFunc(EventNtuple* ntuple, const PhysicsProcessNEW* proc);
   // this function fills all of the plots for a given process
   void fillPlots(map<string, Plot*> &  plots, EventNtuple * ntuple, METree * metree, MicroNtuple * mnt, vector<TString>, double weight = 1.0);
   // returns a boolean if the event passes the specified cuts
   bool eventPassCuts(EventNtuple * ntuple, const PhysicsProcessNEW*);
   // returns a double
   double weightFunc(EventNtuple* ntuple, const PhysicsProcessNEW* proc);
   // Is run once for each process before events (initializes PU Reweighting
   void processFunc(EventNtuple* ntuple, const PhysicsProcessNEW* proc);
   // returns a number based on if the particle is a lepton or quark based on its pdgId
   int leptonOrQuark(int x);
   // returns the charge of a particle based on its pdgId
   double charge(int x);
   // returns the mass of the hadronic W and the leptonic W
   pair<double,double> onVsOffShell(EventNtuple * ntuple);
   pair<double,double> onVsOffShellInclusive(EventNtuple * ntuple);
   // returns a number rounded away from zero
   double round(double r);
   int round_int(double r);
   // concat any two streamable objects into a string
   template <class T, class U>
   std::string concatString(const T& obj1, const U& obj2);
}

////////////////////////////////////////////////////////////////////////////////
//  Implement User Functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void UserFunctions::fillPlots(map<string, Plot*> &  plots, EventNtuple * ntuple,  METree * metree, MicroNtuple * mnt, vector<TString> MVAMethods, double weight)
{
   //WARNING!!! Make sure that ntuple, metree, and mnt are set before using!!! Example::metree and mnt will not be set if you are running on an MEInput file.

   if (ntuple) {
      //-----------------ELECTRON--------------------//
      if (ntuple->lLV[0].leptonCat == DEFS::electron && (leptonCat == DEFS::electron || leptonCat == DEFS::both))
      {
      
         double Mjj = (ntuple->jLV[0] + ntuple->jLV[1]).M();
         double WmT = sqrt(pow(ntuple->lLV[0].Et()+ntuple->METLV[0].Et(), 2) - pow(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(), 2) - pow(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(), 2));
         pair<double,double> leptVsHadWMass = onVsOffShellInclusive(ntuple);

         if(!doElectron)
            return;

         plots["Mjj_electron"]->Fill(Mjj,weight);
         plots["LeptPt_electron"]->Fill(ntuple->lLV[0].Pt(),weight);
         plots["LeptEta_electron"]->Fill(ntuple->lLV[0].Eta(),weight);
         plots["LeptPhi_electron"]->Fill(ntuple->lLV[0].Phi(),weight);
         plots["MET_electron"]->Fill(ntuple->METLV[0].Et(),weight);
         plots["WmT_electron"]->Fill(WmT, weight);
         plots["MjjmWmT_electron"]->Fill(Mjj - WmT, weight);
         plots["Jet1Pt_electron"]->Fill(ntuple->jLV[0].Pt(),weight);
         plots["Jet1Eta_electron"]->Fill(ntuple->jLV[0].Eta(),weight);
         plots["Jet1Phi_electron"]->Fill(ntuple->jLV[0].Phi(),weight);
         plots["Jet2Pt_electron"]->Fill(ntuple->jLV[1].Pt(),weight);
         plots["Jet2Eta_electron"]->Fill(ntuple->jLV[1].Eta(),weight);
         plots["Jet2Phi_electron"]->Fill(ntuple->jLV[1].Phi(),weight);
         plots["DeltaEtaJ1J2_electron"]->Fill(TMath::Abs(ntuple->jLV[0].Eta() - ntuple->jLV[1].Eta()),weight);
         plots["Ptjj_electron"]->Fill((ntuple->jLV[0] + ntuple->jLV[1]).Pt(),weight);
         plots["j1Pt_Mjj_electron"]->Fill(ntuple->jLV[0].Pt() / ntuple->Mjj,weight);
         plots["j2Pt_Mjj_electron"]->Fill(ntuple->jLV[1].Pt() / ntuple->Mjj,weight);
         plots["Mlvjj_electron"]->Fill((ntuple->jLV[0] + ntuple->jLV[1] + ntuple->lLV[0] + ntuple->METLV[0]).M(),weight);
         plots["DeltaRLepMET_electron"]->Fill(sqrt(pow(ntuple->lLV[0].Eta()-ntuple->METLV[0].Eta(),2)+pow(ntuple->lLV[0].Phi()-ntuple->METLV[0].Phi(),2)),weight);
         plots["EJ1EJ2_electron"]->Fill(ntuple->jLV[0].E() * ntuple->jLV[1].E(),weight);
         plots["BetaJ1BetaJ2_electron"]->Fill(ntuple->jLV[0].Beta() * ntuple->jLV[1].Beta(),weight);
         plots["DeltaRJ1J2_electron"]->Fill(sqrt(pow(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta(),2)+pow(ntuple->jLV[0].Phi()-ntuple->jLV[1].Phi(),2)),weight);
         plots["AngleJ1J2_electron"]->Fill(ntuple->jLV[0].Angle(ntuple->jLV[1].Vect()),weight);
         plots["DeltaPhi_LJ1_electron"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]),weight);
         plots["DeltaPhi_METJ1_electron"]->Fill(ntuple->jLV[0].DeltaPhi(ntuple->METLV[0]),weight);
         plots["DeltaPhi_J1J2_electron"]->Fill(ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]),weight);
         plots["DeltaPhi_LJ1_vs_J1J2_electron"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]), ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]),weight);
         plots["npv_electron"]->Fill(ntuple->vLV[0].npv,weight);
         plots["jjlvPhi_electron"]->Fill((ntuple->jLV[0] + ntuple->jLV[1]).Phi() - (ntuple->lLV[0] + ntuple->METLV[0]).Phi(),weight);
         plots["MWjjVsMWlv_electron"]->Fill(leptVsHadWMass.first,leptVsHadWMass.second,weight);
         if (ntuple->lLV[0].lQ == 1){
            plots["DeltaPhi_LJ1_vs_J1J2_Positive_electron"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]), ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]), weight);
            plots["WmT_Positive_electron"]->Fill(WmT, weight);
         }
         if (ntuple->lLV[0].lQ == -1){
            plots["DeltaPhi_LJ1_vs_J1J2_Negative_electron"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]), ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]), weight);
            plots["WmT_Negative_electron"]->Fill(WmT, weight);
         }
         plots["WmT_Subtracted_electron"]->Fill(WmT, (ntuple->lLV[0].lQ)*weight);
         plots["DeltaPhi_LJ1_vs_J1J2_Subtracted_electron"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]), ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]), (ntuple->lLV[0].lQ)*weight);

         vector<Double_t> coord;
         //coord.assign(((TProfileMDF*)plots["lpt_leta_j1pt_j1eta_j2pt_j2eta_electron"]->templateHisto)->GetNaxis(),0);
         //coord.assign(((TProfileMDF*)plots["lpt_lpt_j1pt_j1pt_j2pt_j2pt_electron"]->templateHisto)->GetNaxis(),0);
         coord.assign(((TProfileMDF*)plots["Mjj_Mjj_Mt_MET_DeltaR_DeltaR_electron"]->templateHisto)->GetNaxis(),0);
         //coord[0] = ntuple->lLV[0].Pt();
         coord[0] = Mjj;
         //coord[1] = TMath::Abs(ntuple->lLV[0].Eta());
         //coord[1] = ntuple->lLV[0].Pt();
         coord[1] = Mjj;
         //coord[2] = ntuple->jLV[0].Pt();
         coord[2] = WmT;
         //coord[3] = TMath::Abs(ntuple->jLV[0].Eta());
         //coord[3] = ntuple->jLV[0].Pt();
         coord[3] = ntuple->METLV[0].Et();
         //coord[4] = ntuple->jLV[1].Pt();
         coord[4] = sqrt(pow(ntuple->lLV[0].Eta()-ntuple->jLV[0].Eta(),2)+pow(ntuple->lLV[0].Phi()-ntuple->jLV[0].Phi(),2));
         //coord[5] = TMath::Abs(ntuple->jLV[1].Eta());
         //coord[5] = ntuple->jLV[1].Pt();
         coord[5] = sqrt(pow(ntuple->lLV[0].Eta()-ntuple->jLV[0].Eta(),2)+pow(ntuple->lLV[0].Phi()-ntuple->jLV[0].Phi(),2));
         //plots["lpt_leta_j1pt_j1eta_j2pt_j2eta_electron"]->Fill(coord,1.0,weight);
         //plots["lpt_lpt_j1pt_j1pt_j2pt_j2pt_electron"]->Fill(coord,1.0,weight);
         plots["Mjj_Mjj_Mt_MET_DeltaR_DeltaR_electron"]->Fill(coord,1.0,weight);
      }
      // -------------MUON-------------------//
      else if (ntuple->lLV[0].leptonCat == DEFS::muon && (leptonCat == DEFS::muon || leptonCat == DEFS::both))
      {
         double Mjj = (ntuple->jLV[0] + ntuple->jLV[1]).M();
         double WmT = sqrt(pow(ntuple->lLV[0].Et()+ntuple->METLV[0].Et(), 2) - pow(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(), 2) - pow(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(), 2));
         pair<double,double> leptVsHadWMass = onVsOffShellInclusive(ntuple);
      
         if(!doMuon)
            return;

         plots["Mjj_muon"]->Fill(Mjj,weight);
         plots["LeptPt_muon"]->Fill(ntuple->lLV[0].Pt(),weight);
         plots["LeptEta_muon"]->Fill(ntuple->lLV[0].Eta(),weight);
         plots["LeptPhi_muon"]->Fill(ntuple->lLV[0].Phi(),weight);
         plots["MET_muon"]->Fill(ntuple->METLV[0].Et(),weight);
         plots["WmT_muon"]->Fill(WmT, weight);
         plots["MjjmWmT_muon"]->Fill(Mjj - WmT, weight);
         plots["Jet1Pt_muon"]->Fill(ntuple->jLV[0].Pt(),weight);
         plots["Jet1Eta_muon"]->Fill(ntuple->jLV[0].Eta(),weight);
         plots["Jet1Phi_muon"]->Fill(ntuple->jLV[0].Phi(),weight);
         plots["Jet2Pt_muon"]->Fill(ntuple->jLV[1].Pt(),weight);
         plots["Jet2Eta_muon"]->Fill(ntuple->jLV[1].Eta(),weight);
         plots["Jet2Phi_muon"]->Fill(ntuple->jLV[1].Phi(),weight);
         plots["DeltaEtaJ1J2_muon"]->Fill(TMath::Abs(ntuple->jLV[0].Eta() - ntuple->jLV[1].Eta()),weight);
         plots["Ptjj_muon"]->Fill((ntuple->jLV[0] + ntuple->jLV[1]).Pt(),weight);
         plots["j1Pt_Mjj_muon"]->Fill(ntuple->jLV[0].Pt() / ntuple->Mjj,weight);
         plots["j2Pt_Mjj_muon"]->Fill(ntuple->jLV[1].Pt() / ntuple->Mjj,weight);
         plots["Mlvjj_muon"]->Fill((ntuple->jLV[0] + ntuple->jLV[1] + ntuple->lLV[0] + ntuple->METLV[0]).M(),weight);
         plots["DeltaRLepMET_muon"]->Fill(sqrt(pow(ntuple->lLV[0].Eta()-ntuple->METLV[0].Eta(),2)+pow(ntuple->lLV[0].Phi()-ntuple->METLV[0].Phi(),2)),weight);
         plots["EJ1EJ2_muon"]->Fill(ntuple->jLV[0].E() * ntuple->jLV[1].E(),weight);
         plots["BetaJ1BetaJ2_muon"]->Fill(ntuple->jLV[0].Beta() * ntuple->jLV[1].Beta(),weight);
         plots["DeltaRJ1J2_muon"]->Fill(sqrt(pow(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta(),2)+pow(ntuple->jLV[0].Phi()-ntuple->jLV[1].Phi(),2)),weight);
         plots["AngleJ1J2_muon"]->Fill(ntuple->jLV[0].Angle(ntuple->jLV[1].Vect()),weight);
         plots["DeltaPhi_LJ1_muon"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]),weight);
         plots["DeltaPhi_METJ1_muon"]->Fill(ntuple->jLV[0].DeltaPhi(ntuple->METLV[0]),weight);
         plots["DeltaPhi_J1J2_muon"]->Fill(ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]),weight);   
         plots["DeltaPhi_LJ1_vs_J1J2_muon"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]), ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]),weight);
         plots["npv_muon"]->Fill(ntuple->vLV[0].npv,weight);
         plots["jjlvPhi_muon"]->Fill((ntuple->jLV[0] + ntuple->jLV[1]).Phi() - (ntuple->lLV[0] + ntuple->METLV[0]).Phi(),weight);
         plots["MWjjVsMWlv_muon"]->Fill(leptVsHadWMass.first,leptVsHadWMass.second,weight);
         if (ntuple->lLV[0].lQ == 1){
            plots["DeltaPhi_LJ1_vs_J1J2_Positive_muon"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]), ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]), weight);
            plots["WmT_Positive_muon"]->Fill(WmT, weight);
         }
         if (ntuple->lLV[0].lQ == -1){
            plots["DeltaPhi_LJ1_vs_J1J2_Negative_muon"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]), ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]), weight);
            plots["WmT_Negative_muon"]->Fill(WmT, weight);
         }
         plots["WmT_Subtracted_muon"]->Fill(WmT, (ntuple->lLV[0].lQ)*weight);
         plots["DeltaPhi_LJ1_vs_J1J2_Subtracted_muon"]->Fill(ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]), ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]), (ntuple->lLV[0].lQ)*weight);

         vector<Double_t> coord;
         //coord.assign(((TProfileMDF*)plots["lpt_leta_j1pt_j1eta_j2pt_j2eta_muon"]->templateHisto)->GetNaxis(),0);
         //coord.assign(((TProfileMDF*)plots["lpt_lpt_j1pt_j1pt_j2pt_j2pt_muon"]->templateHisto)->GetNaxis(),0);
         coord.assign(((TProfileMDF*)plots["Mjj_Mjj_Mt_MET_DeltaR_DeltaR_muon"]->templateHisto)->GetNaxis(),0);
         //coord[0] = ntuple->lLV[0].Pt();
         coord[0] = Mjj;
         //coord[1] = TMath::Abs(ntuple->lLV[0].Eta());
         //coord[1] = ntuple->lLV[0].Pt();
         coord[1] = Mjj;
         //coord[2] = ntuple->jLV[0].Pt();
         coord[2] = WmT;
         //coord[3] = TMath::Abs(ntuple->jLV[0].Eta());
         //coord[3] = ntuple->jLV[0].Pt();
         coord[3] = ntuple->METLV[0].Et();
         //coord[4] = ntuple->jLV[1].Pt();
         coord[4] = sqrt(pow(ntuple->lLV[0].Eta()-ntuple->jLV[0].Eta(),2)+pow(ntuple->lLV[0].Phi()-ntuple->jLV[0].Phi(),2));
         //coord[5] = TMath::Abs(ntuple->jLV[1].Eta());
         //coord[5] = ntuple->jLV[1].Pt();
         coord[5] = sqrt(pow(ntuple->lLV[0].Eta()-ntuple->jLV[0].Eta(),2)+pow(ntuple->lLV[0].Phi()-ntuple->jLV[0].Phi(),2));
         //plots["lpt_lpt_j1pt_j1pt_j2pt_j2pt_muon"]->Fill(coord,1.0,weight);
         plots["Mjj_Mjj_Mt_MET_DeltaR_DeltaR_muon"]->Fill(coord,1.0,weight);
      }
   }

   if (metree) { 
      //-----------------ELECTRON--------------------//
      if (ntuple->lLV[0].leptonCat == DEFS::electron && (leptonCat == DEFS::electron || leptonCat == DEFS::both)) {
         for (unsigned int tep=0; tep<15; tep++) {
            string name = UserFunctions::concatString("tEventProb",tep) + "_electron";
            //cout << name << " = " << metree->getProbStat(tep)->tEventProb << " (" << (Float_t*)(&metree->getProbStat(tep)->tEventProb)<< ")"<<endl;
            plots[name]->Fill(log(metree->getProbStat(tep)->tEventProb),weight);
         }

         double tEventProb0 = metree->getProbStat(0)->tEventProb;
         double tEventProb1 = metree->getProbStat(1)->tEventProb;
         double tEventProb3 = metree->getProbStat(3)->tEventProb;
         double signal = (tEventProb0 / 0.8e-06) + (tEventProb1 / 0.1e-06);
         double back   = tEventProb3 / 0.75e-03;
         double epd    = signal /(signal+back);
         plots["epdPretagWWandWZ_RE_electron"]->Fill(epd,weight);
      }
      // -------------MUON-------------------//
      else if (ntuple->lLV[0].leptonCat == DEFS::muon && (leptonCat == DEFS::muon || leptonCat == DEFS::both)) {
         for (unsigned int tep=0; tep<15; tep++) {
            string name = UserFunctions::concatString("tEventProb",tep) + "_muon";
            //cout << name << " = " << metree->getProbStat(tep)->tEventProb << " (" << (Float_t*)(&metree->getProbStat(tep)->tEventProb)<< ")"<<endl;
            plots[name]->Fill(log(metree->getProbStat(tep)->tEventProb),weight);
         }

         double tEventProb0 = metree->getProbStat(0)->tEventProb;
         double tEventProb1 = metree->getProbStat(1)->tEventProb;
         double tEventProb3 = metree->getProbStat(3)->tEventProb;
         double signal = (tEventProb0 / 0.8e-06) + (tEventProb1 / 0.1e-06);
         double back   = tEventProb3 / 0.75e-03;
         double epd    = signal /(signal+back);
         plots["epdPretagWWandWZ_RE_muon"]->Fill(epd,weight);
      }
   }

   if (mnt) {
      //-----------------ELECTRON--------------------//
      if (ntuple->lLV[0].leptonCat == DEFS::electron && (leptonCat == DEFS::electron || leptonCat == DEFS::both)) {
         plots["epdPretagWWandWZ_electron"]->Fill(log(mnt->epdPretagWWandWZ),weight);

         if (!MVAMethods.empty() && ntuple) {
            //cout << "MVADiscriminator_electron = " << mnt->getMVAOutput(MVAMethods).front()["response"] << endl;
            //cout << "MVAProbability_electron = " << mnt->getMVAOutput(MVAMethods).front()["probability"] << endl;
            //cout << "weight = " << weight << endl;
            mnt->setMjjMVA(ntuple->Mjj);
            plots["MVADiscriminator_electron"]->Fill(mnt->getMVAOutput(MVAMethods).front()["response"],weight);
            plots["MVAProbability_electron"]->Fill(mnt->getMVAOutput(MVAMethods).front()["probability"],weight);
         }
      }
      // -------------MUON-------------------//
      else if (ntuple->lLV[0].leptonCat == DEFS::muon && (leptonCat == DEFS::muon || leptonCat == DEFS::both)) {
         plots["epdPretagWWandWZ_muon"]->Fill(log(mnt->epdPretagWWandWZ),weight);

         if (!MVAMethods.empty() && ntuple) {
            //cout << "MVADiscriminator_muon = " << mnt->getMVAOutput(MVAMethods).front()["response"] << endl;
            //cout << "MVAProbability_muon = " << mnt->getMVAOutput(MVAMethods).front()["probability"] << endl;
            //cout << "weight = " << weight << endl;
            mnt->setMjjMVA(ntuple->Mjj);
            plots["MVADiscriminator_muon"]->Fill(mnt->getMVAOutput(MVAMethods).front()["response"],weight);
            plots["MVAProbability_muon"]->Fill(mnt->getMVAOutput(MVAMethods).front()["probability"],weight);
         }
      }
   }
    
}//fillPlots

//______________________________________________________________________________
// Return true if the event pass the cuts imposed to the given lepton category
bool UserFunctions::eventPassCuts(EventNtuple * ntuple, const PhysicsProcessNEW*){
   
  /*
  ntuple->lLV[0].leptonCat  |    req   | result
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
  // Kill only electron-muon combinations, nothing else
  if ( ntuple->lLV[0].leptonCat != leptonCat && (leptonCat != DEFS::both && ntuple->lLV[0].leptonCat != DEFS::both) )
    return false;

/*
  // Keep only specific events
  if (ntuple->event!=2663602)
     return false;
*/
/*
  unsigned int count = ntuple->jLV.size();
  if (ntuple->lLV[0].leptonCat == DEFS::muon)
     for (unsigned int i=0; i<ntuple->jLV.size(); i++) {
        if ((ntuple->jLV[i].Pt()) <= 25.0 )
           count--;
     }
  
  if (ntuple->lLV[0].leptonCat == DEFS::electron)
     for (unsigned int i=0; i<ntuple->jLV.size(); i++) {
        if ((ntuple->jLV[i].Pt()) <= 25.0 )
           count--;
     }
  //if(count<ntuple->jLV.size())
  if(count<2)
     return false;
*/
/*
  if (ntuple->leptonCat == DEFS::muon)
    if ((ntuple->lLV[0].Pt()) <= 25.0 || (ntuple->METLV[0].Et()) <= 30.0)
      return false;
  
  if (ntuple->leptonCat == DEFS::electron) 
    if ((ntuple->lLV[0].Pt()) <= 30.0 || (ntuple->METLV[0].Et()) <= 35.0)
      return false;
 */

/*
  // Diboson analysis cuts
  TLorentzVector mt(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(),ntuple->lLV[0].Py()+ntuple->METLV[0].Py(),0,ntuple->lLV[0].Et()+ntuple->METLV[0].Et());
  double wmt = sqrt(pow(ntuple->lLV[0].Et()+ntuple->METLV[0].Et(), 2) - pow(ntuple->lLV[0].Px()+ntuple->METLV[0].Px(), 2) - pow(ntuple->lLV[0].Py()+ntuple->METLV[0].Py(), 2));

  if (ntuple->leptonCat == DEFS::muon)
     if ((ntuple->lLV[0].Pt()) <= 25.0                                  ||
         (ntuple->METLV[0].Et()) <= 25.0                                ||
         (ntuple->jLV[0].Pt()) <= 35.0                                  ||
         (ntuple->jLV[1].Pt()) <= 35.0                                  ||
         (TMath::Abs(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta())) >= 1.5 || 
         ((ntuple->jLV[0]+ntuple->jLV[1]).Pt()) <= 20.0                 ||
         (ntuple->METLV[0].DeltaPhi(ntuple->jLV[0])) <= 0.4             ||
         (mt.M()) <= 30.0                                               ||
         wmt <= 50.0                                                     )
        return false;
  
  if (ntuple->leptonCat == DEFS::electron) 
     if ((ntuple->lLV[0].Pt()) <= 35.0                                  ||
         (ntuple->METLV[0].Et()) <= 30.0                                ||
         (ntuple->jLV[0].Pt()) <= 35.0                                  ||
         (ntuple->jLV[1].Pt()) <= 35.0                                  ||
         (TMath::Abs(ntuple->jLV[0].Eta()-ntuple->jLV[1].Eta())) >= 1.5 || 
         ((ntuple->jLV[0]+ntuple->jLV[1]).Pt()) <= 20.0                 ||
         (ntuple->METLV[0].DeltaPhi(ntuple->jLV[0])) <= 0.4             ||
         (mt.M()) <= 30.0                                               ||
         wmt <= 50.0                                                     )
        return false;
*/

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

  //X axis cuts
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
  else if (cutRegion.Contains("all"))
    return true;

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
   if (doPUreweight){
     weight *= puweight->getWeight(ntuple->vLV[0].tnpus[1]);
   }
   
   // WJets weighting (specific to fix shape issues)
   if (WJweight){
      if (auxName.Contains("WJETS")){
         double dpLJ = ntuple->lLV[0].DeltaPhi(ntuple->jLV[0]);
         double dpJJ = ntuple->jLV[0].DeltaPhi(ntuple->jLV[1]);
         int bin     = wJetsWeight->FindBin(dpLJ,dpJJ);
         weight *= wJetsWeight->GetBinContent(bin);
      }
   }

   return weight;
}

void UserFunctions::initEventFunc(EventNtuple* ntuple, const PhysicsProcessNEW* proc)
{
  TString auxName = proc->name;
  auxName.ToUpper();
  if(auxName.Contains("DATA") || auxName.Contains("QCD"))
    return;

   // Performs JER
  if (doJER)
    ntuple->doJER();
}

void UserFunctions::processFunc(EventNtuple* ntuple, const PhysicsProcessNEW* proc)
{
   doMuon = true;
   doElectron = true;
   if(((PlotterPhysicsProcessNEW*)proc)->leptonCat == DEFS::electron)
      doMuon = false;
   if(((PlotterPhysicsProcessNEW*)proc)->leptonCat == DEFS::muon)
      doElectron = false;
   
  TString auxName = proc->name;
  auxName.ToUpper();

  if (WJweight){
     TFile* Subs;
     if(auxName.Contains("WJETS")){
        Subs = TFile::Open(TString("/uscms/home/amejia94/CMSSW_5_2_5/src/subtracted_"+DEFS::getLeptonCatString(UserFunctions::leptonCat)+".root"));
        if (!Subs->IsOpen()) {
           cout << "\tERROR::Weight file (subtracted_"+DEFS::getLeptonCatString(UserFunctions::leptonCat)+".root) was not opened" << endl
                << "\tWJets histograms will not be filled" << endl;
           return;
        }
        TH2D* auxh  = (TH2D*) Subs->Get("hdSubtract");
        if (!auxh) {
           cout << "\tERROR::Weight hist hdSubtract could not be found int subtracted_"+DEFS::getLeptonCatString(UserFunctions::leptonCat)+".root" << endl
                << "\tWJets histograms will not be filled" << endl;
           return;
        }
        wJetsWeight = (TH2D*) auxh->Clone();
        if (!wJetsWeight) {
           cout << "\tERROR::Weight hist wJetsWeight was not cloned properly" << endl
                << "\tWJets histograms will not be filled" << endl;
           return;
        }
        wJetsWeight->SetDirectory(0);
      Subs->Close();
     }
  }
  
  auxName.ToUpper();
  if(auxName.Contains("DATA") || auxName.Contains("QCD"))
    return;
   
  // Initializes PU Reweighting
  string dataname = string("/uscms/home/aperloff/MatrixElement/CMSSW_5_3_2_patch4/src/TAMUWW/Tools/bin/pileup12_noTrig.root");
  string MCname   = string(TString(proc->fileName));
  puweight = new PUreweight(dataname,MCname,"pileup_noTrig","PS/TPUDist");
}

double UserFunctions::charge(int x) {
   int absx = abs(x);
   if (absx==2 || absx==4 || absx==6)
      return (x > 0) ? (2./3.) : ((x < 0) ? -(2./3.) : 0);
   else if (absx==1 || absx==3 || absx==5)
      return (x > 0) ? (-1./3.) : ((x < 0) ? (1./3.) : 0);
   else if (absx==11 || absx==13 || absx==15)
      return (x > 0) ? -1 : ((x < 0) ? 1 : 0);
   else if (absx==24)
      return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
   else
      return 0;
}//charge

int UserFunctions::leptonOrQuark(int x)
{
   int absx = abs(x);
   if (absx>10 && absx < 17)
      return 0;
   else if (absx>0 && absx < 7)
      return 1;
   else
      return -1;
}//leptonOrQuark


pair<double,double> UserFunctions::onVsOffShell(EventNtuple * ntuple) {
   if (ntuple->genParticleCollection.size()==0) {
      //cout << "WARNING::No genParticleCollection present." << endl;
      return make_pair(0.0,0.0);
   }

   pair<int,int> fW = make_pair(0,0); //first = sign, second = position                                     
   pair<int,int> sW = make_pair(0,0);
   int sl  = 0;
   double cj1 = 0;
   double cj2 = 0;
   int sjj = 0;

   for (unsigned int i=0; i<ntuple->genParticleCollection.size(); i++) {
      if (abs(ntuple->genParticleCollection[i].pdgId)==24) {
         if (fW.first==0) {
            fW.first = charge(ntuple->genParticleCollection[i].pdgId);
            fW.second = i;
         }
         else {
            sW.first = charge(ntuple->genParticleCollection[i].pdgId);
            sW.second = i;
         }
      }
      if (abs(ntuple->genParticleCollection[i].pdgId)==11 || abs(ntuple->genParticleCollection[i].pdgId)==13 ||
          abs(ntuple->genParticleCollection[i].pdgId)==15) {
         sl = charge(ntuple->genParticleCollection[i].pdgId);
         if (sl==0)
            cout << "sl==0 and lept pdgId = " << (ntuple->genParticleCollection[i].pdgId) << endl;
      }
      if (abs(ntuple->genParticleCollection[i].pdgId)>0 && abs(ntuple->genParticleCollection[i].pdgId)<9) {
         if (cj1!=0.0 && cj2!=0.0)
            cout << "WARNING::There were more than two jets found. Both cj1 and cj2 are non-zero." << endl;
         /*
         cout << "pdgId" << i << " = " << ntuple->genParticleCollection[i].pdgId << "\tmother_pdgIds = ";
         for (int j=0; j<ntuple->genParticleCollection[i].numberOfMothers; j++) {
            cout << ntuple->genParticleCollection[ntuple->genParticleCollection[i].motherPositions[j]].pdgId << ", ";
         }
         cout << endl;
         */
         bool Wjet = true;
         for (int j=0; j<ntuple->genParticleCollection[i].numberOfMothers; j++) {
            if (abs(ntuple->genParticleCollection[ntuple->genParticleCollection[i].motherPositions[j]].pdgId)!=24)
               Wjet = false;
         }

         if (Wjet) {
            if (cj1==0)
               cj1 = charge(ntuple->genParticleCollection[i].pdgId);
            else {
               cj2 = charge(ntuple->genParticleCollection[i].pdgId);
            }
         }
      }
   }

   if (cj1!=0.0 && cj2!=0.0) {
      sjj = round_int(cj1 + cj2);
      //cout << "sjj = " << sjj << "\tcj1 = " << cj1 << "\tcj2 = " << cj2 << "\tcj1+cj2 = " << cj1+cj2 
      //     << "\tround_int(cj1+cj2) = " << round_int(cj1+cj2) << endl;
   }

   if (fW.first==sl && sW.first==sjj) {
      //h->Fill(ntuple->genParticleCollection[sW.second].p4.M(),
      //        ntuple->genParticleCollection[fW.second].p4.M());
      return make_pair(ntuple->genParticleCollection[sW.second].p4.M(),
                       ntuple->genParticleCollection[fW.second].p4.M());
   }
   else if (fW.first==sjj && sW.first==sl) {
      //h->Fill(ntuple->genParticleCollection[fW.second].p4.M(),
      //        ntuple->genParticleCollection[sW.second].p4.M());
      return make_pair(ntuple->genParticleCollection[fW.second].p4.M(),
                       ntuple->genParticleCollection[sW.second].p4.M());
   }
   else {
     //cout << "WARNING::Unable to determine which W is hadronic and which W is leptonic." << endl
     //      << "fW.first = " << fW.first << "\tsW.first = " << sW.first << "\tsl = " << sl << "\tcj1 = "
     //      << cj1 << "\tcj2 = " << cj2 << "\tsjj = " << sjj << endl;
      return make_pair(0.0,0.0);
   }

}//onVsOffShell

pair<double,double> UserFunctions::onVsOffShellInclusive(EventNtuple * ntuple)
{
   if (ntuple->genParticleCollection.size()==0) {
      //cout << "WARNING::No genParticleCollection present." << endl;
      return make_pair(0.0,0.0);
   }
   
   //0=leptonic, 1=hadronic, 2=both, 3=neither/not set/unknown
   vector<pair<int,int> > W;  //first = category, second = position
   int lc = 0; //lepton counter
   int qc = 0; //quark counter
   for (unsigned int i=0; i<ntuple->genParticleCollection.size(); i++) {
      if (abs(ntuple->genParticleCollection[i].pdgId)!=25)
	continue;
      else {
	if (verbose) 
	  cout << "H->";
	for (unsigned int j=0; j<ntuple->genParticleCollection[i].daughterPositions.size(); j++) {
	  if (ntuple->genParticleCollection[i].daughterPositions[j]<=500 && abs(ntuple->genParticleCollection[ntuple->genParticleCollection[i].daughterPositions[j]].pdgId)==24) {
	    if (verbose) 
	      cout << "W";
	    W.push_back(make_pair(3,ntuple->genParticleCollection[i].daughterPositions[j]));
	  }
	} //loop through the daughters of the Higgs
	if (verbose)
	  cout << "->";
    if(W.size()==2) {
       for (unsigned int j=0; j<W.size(); j++) {
          lc = 0;
          qc = 0;
          for (unsigned int k=0; k<ntuple->genParticleCollection[W[j].second].daughterPositions.size(); k++) {
             if (W[j].first==2)
                continue;
             else if (ntuple->genParticleCollection[W[j].second].daughterPositions[k]<=500 && leptonOrQuark(ntuple->genParticleCollection[ntuple->genParticleCollection[W[j].second].daughterPositions[k]].pdgId)==0) {
                if (abs(ntuple->genParticleCollection[ntuple->genParticleCollection[W[j].second].daughterPositions[k]].pdgId)==11 && ntuple->genParticleCollection[ntuple->genParticleCollection[W[j].second].daughterPositions[k]].p4.Pt()<=27) {
                   if (verbose)
                      cout << "CUT" << endl;
                   return make_pair(-1.0,-1.0);
                }
                if (abs(ntuple->genParticleCollection[ntuple->genParticleCollection[W[j].second].daughterPositions[k]].pdgId)==13 && ntuple->genParticleCollection[ntuple->genParticleCollection[W[j].second].daughterPositions[k]].p4.Pt()<=24) {
                   if (verbose)
                      cout << "CUT" << endl;
                   return make_pair(-1.0,-1.0);
                }
                lc++;
             }
             else if (ntuple->genParticleCollection[W[j].second].daughterPositions[k]<=500 && leptonOrQuark(ntuple->genParticleCollection[ntuple->genParticleCollection[W[j].second].daughterPositions[k]].pdgId)==1) {
                qc++;
             }
             if (ntuple->genParticleCollection[W[j].second].daughterPositions[k]<=500 && ntuple->genParticleCollection[ntuple->genParticleCollection[W[j].second].daughterPositions[k]].pdgId!=0)
                if (verbose)
                   cout << ntuple->genParticleCollection[ntuple->genParticleCollection[W[j].second].daughterPositions[k]].pdgId
                        << ",";
          } //loop through the daughters of the W
          if (lc==2 && qc==0)
             W[j].first = 0;
          else if (lc==0 && qc==2)
             W[j].first = 1;
          else if (lc>0 && qc>0)
             W[j].first = 2;
          else {
             cout << "Something is wonky!!!" << endl;
             W[j].first = 3;
          }
       } //loop through all of the W
    } //if there are only 2 W
	if (verbose)
	  cout << endl;
      } //find the Higgs
   } //loop through the gen particles
   
   if (W.size()==0) {
     //cout << "WARNING::No W decaying from a Higgs was found." << endl;
     return make_pair(-1.0,-1.0);
   }
   else if (W[0].first==0 && W[1].first==1)
     return make_pair(ntuple->genParticleCollection[W[1].second].p4.M(),
		      ntuple->genParticleCollection[W[0].second].p4.M());
   else if (W[0].first==1 && W[1].first==0)
     return make_pair(ntuple->genParticleCollection[W[0].second].p4.M(),
		      ntuple->genParticleCollection[W[1].second].p4.M());
   else {
     cout << "WARNING::Unable to determine which W is hadronic and which W is leptonic." << endl
	  << "W[0].first = " << W[0].first << "\tW1[].first = " << W[1].first << endl
	  << "W[0].second = " << W[0].second << "\tW1[].second = " << W[1].second << endl;
     return make_pair(-1.0,-1.0);
   }
}


double UserFunctions::round(double r) {
   return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}//round

int UserFunctions::round_int(double r) {
   return (r > 0.0) ? (r + 0.5) : (r - 0.5); 
}//round

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
void doPlotter(TString fileName, map<string, Plot*> & plots, vector<PhysicsProcessNEW*> procs, bool doJER, bool doPUrewt, bool doFNAL, int maxEvts, bool WJweight, TString MVAWeightDir, vector<TString> MVAMethods, bool verbose);

/// returns the cross section for the given process
double getCrossSection(TString channelName);

/// returns the branching ratio for the given process
double getBranchingRatio(TString channelName);

/// returns the number of Monte Carlo events for the given process
double getNumMCEvts(TString channelName);

/// returns a vector containing all of the processes that the program will use
vector<PhysicsProcessNEW*> getProcesses(DEFS::LeptonCat leptonCat, double intLum);

/// returns a vector containing all of the processes that the program will use
vector<PhysicsProcessNEW*> getProcesses(vector<string>& processNames, DEFS::LeptonCat leptonCat, double intLum);

/// returns a map containing all of the plots that will be made for each process and their specific attributes
map<string, Plot*> getPlots(string leptonCatString, bool norm_data);

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
   double         intLum            = cl.getValue<double>  ("lum",          10000);  
   string         lepCat            = cl.getValue<string>  ("lep",         "both");
   UserFunctions::doJER             = cl.getValue<bool>    ("doJer",         true);
   UserFunctions::doPUreweight      = cl.getValue<bool>    ("doPUrewt",      true);
   UserFunctions::doFNAL            = cl.getValue<bool>    ("doFNAL",       false);
   UserFunctions::cutRegion         = cl.getValue<string>  ("cutRegion",    "all");
   UserFunctions::outDir            = cl.getValue<string>  ("outDir",         ".");
   UserFunctions::leptonCat         = DEFS::getLeptonCat   (lepCat);
   UserFunctions::WJweight          = cl.getValue<bool>    ("WJweight",     false);
   UserFunctions::verbose           = cl.getValue<bool>    ("verbose",      false);
   bool           norm_data         = cl.getValue<bool>    ("norm_data",    false);
   int            maxEvts           = cl.getValue<int>     ("maxEvts",          0);
   TString        MVAWeightDir      = cl.getValue<TString> ("MVAWeightDir",    "");
   vector<TString> MVAMethods       = cl.getVector<TString>("MVAMethods",      "");

   if (!cl.check()) return 0;
   cl.print();
   
   TBenchmark* m_benchmark = new TBenchmark();
   m_benchmark->Reset();
   m_benchmark->Start("event");

   // The vector containing all plots to be made
   map<string, Plot*> plots = getPlots(lepCat,norm_data);
   
   cout << "Doing lepton category "<<DEFS::getLeptonCatString(UserFunctions::leptonCat)<< endl;

   // The vector holding all processes.
   vector <PhysicsProcessNEW*> procs;
   
   if (!processNames.empty())
      procs = getProcesses(processNames,UserFunctions::leptonCat,intLum);
   else
      procs = getProcesses(UserFunctions::leptonCat,intLum);
      
   // Make all the plots and store to outputFile
   if(!UserFunctions::outDir.EndsWith("/"))
      UserFunctions::outDir+="/";
   if(!gSystem->OpenDirectory(UserFunctions::outDir))
      gSystem->mkdir(UserFunctions::outDir);
   TString outFileName = UserFunctions::outDir + "outputFile_";
   outFileName += DEFS::getLeptonCatString(UserFunctions::leptonCat)+".root";
   doPlotter(outFileName, plots, procs, UserFunctions::doJER, UserFunctions::doPUreweight, UserFunctions::doFNAL, maxEvts, UserFunctions::WJweight, MVAWeightDir, MVAMethods, UserFunctions::verbose);
   
   if (UserFunctions::leptonCat == DEFS::muon || UserFunctions::leptonCat == DEFS::both)
   {
      plots["DeltaPhi_LJ1_vs_J1J2_muon"]->saveHistogramsToFile(UserFunctions::outDir + "DeltaPhi_LJ1_vs_J1J2_mu.root");
      plots["DeltaPhi_LJ1_vs_J1J2_Subtracted_muon"]->saveHistogramsToFile(UserFunctions::outDir + "DeltaPhi_LJ1_vs_J1J2_Subtracted_muon.root");
      //plots["lpt_leta_j1pt_j1eta_j2pt_j2eta_muon"]->saveHistogramsToFile(UserFunctions::outDir + "lpt_leta_j1pt_j1eta_j2pt_j2eta_muon.root");
      //plots["lpt_lpt_j1pt_j1pt_j2pt_j2pt_muon"]->saveHistogramsToFile(UserFunctions::outDir + "lpt_lpt_j1pt_j1pt_j2pt_j2pt_muon.root";)
      plots["Mjj_Mjj_Mt_MET_DeltaR_DeltaR_muon"]->saveHistogramsToFile(UserFunctions::outDir + "Mjj_Mjj_Mt_MET_DeltaR_DeltaR_muon.root");
   }
   if (UserFunctions::leptonCat == DEFS::electron || UserFunctions::leptonCat == DEFS::both)
   {
      plots["DeltaPhi_LJ1_vs_J1J2_electron"]->saveHistogramsToFile(UserFunctions::outDir + "DeltaPhi_LJ1_vs_J1J2_el.root");
      plots["DeltaPhi_LJ1_vs_J1J2_Subtracted_electron"]->saveHistogramsToFile(UserFunctions::outDir + "DeltaPhi_LJ1_vs_J1J2_Subtracted_electron.root");
      //plots["lpt_leta_j1pt_j1eta_j2pt_j2eta_electron"]->saveHistogramsToFile(UserFunctions::outDir + "lpt_leta_j1pt_j1eta_j2pt_j2eta_electron.root");
      //plots["lpt_lpt_j1pt_j1pt_j2pt_j2pt_electron"]->saveHistogramsToFile(UserFunctions::outDir + "lpt_lpt_j1pt_j1pt_j2pt_j2pt_electron.root");
      plots["Mjj_Mjj_Mt_MET_DeltaR_DeltaR_electron"]->saveHistogramsToFile(UserFunctions::outDir + "Mjj_Mjj_Mt_MET_DeltaR_DeltaR_electron.root");
   }
   
   m_benchmark->Stop("event");
   cout << "plotter_x" << endl << "\tCPU time = " << m_benchmark->GetCpuTime("event") << " s" << endl
        << "\tReal time = " << m_benchmark->GetRealTime("event") << " s" << endl;
   delete m_benchmark;

   return 0;

}//plotter()


////////////////////////////////////////////////////////////////////////////////
//  Implement Local Functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void doPlotter(TString fileName, map<string, Plot*> & plots, vector<PhysicsProcessNEW*> procs, bool doJER, bool doPUrewt, bool doFNAL, int maxEvts, bool WJweight, TString MVAWeightDir, vector<TString> MVAMethods, bool verbose) {

   PlotFiller pFill(plots, procs, &UserFunctions::fillPlots);
   pFill.setCutFunction(&UserFunctions::eventPassCuts);
   pFill.setWeightFunction(&UserFunctions::weightFunc);
   pFill.setProcessFunction(&UserFunctions::processFunc);
   pFill.setInitializeEventFunction(&UserFunctions::initEventFunc);
   if (maxEvts>0)
      pFill.setMaximumEventsDEBUG(maxEvts); // TEST
   if (!MVAWeightDir.IsNull()) {
      pFill.setMVAWeightDir(MVAWeightDir);
      pFill.setMVAMethods(MVAMethods);
   }
   pFill.run();

   // Will all the info in the plots get the canvas and write it to file
   TFile * outFile = new TFile(fileName,"RECREATE");
   for ( map<string, Plot*>::iterator p = plots.begin(); p != plots.end() ; p++)
   {
      TCanvas* can = ((FormattedPlot*) p->second)->getCanvas(procs);
      
      can->Write();
      
      can->SaveAs(UserFunctions::outDir+"/"+can->GetName()+".png");
      can->SaveAs(UserFunctions::outDir+"/"+can->GetName()+".eps");
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
      cout << "WARNING::getCrossSection::channelName " << channelName << " not recognized. Returning -1 for the cross section." << endl << "The events will have the same scale as the MC sample, but on a negative scale." << endl << "Please check channel names." << endl;
      return -1.;
   }
}//getCrossSection

//______________________________________________________________________________
double getBranchingRatio(TString channelName)
{
   Table table;
   double xsec;

   string basePath = gSystem->pwd();
   basePath = basePath.substr(0,basePath.find("TAMUWW"));
   table.parseFromFile(basePath+string("/TAMUWW/ConfigFiles/Official/BranchingRatios_8TeV.txt"),"TableCellVal");
   TableCell * cell = table.getCellRowColumn(string(channelName),"BranchingRatio");
   if(cell){
      xsec = ((TableCellVal*)cell)->val.value;
      if (xsec==0)
         cout << "WARNING::getBranchingRatio::The branching ratio for " << channelName << " is 0.0 +/- 0.0" << endl;
      return xsec;
   } else{
      cout << "WARNING::getBranchingRatio::channelName " << channelName << " not recognized. Returning -1 for the branching ratio." << endl << "The events will have the same scale as the MC sample, but on a negative scale." << endl << "Please check channel names." << endl;
      return -1.;
   }
}//getBranchingRatio

//______________________________________________________________________________
double getNumMCEvts(TString channelName)
{
   Table table;
   double value;

   string basePath = gSystem->pwd();
   basePath = basePath.substr(0,basePath.find("TAMUWW"));
   table.parseFromFile(basePath+string("/TAMUWW/ConfigFiles/Official/EventsFromMC_commonPATTuple_532.txt"),"TableCellVal");
   TableCell * cell =table.getCellRowColumn(string(channelName),"Events_PATtuple");
   if(cell){
      value = ((TableCellVal*)cell)->val.value;
      if (value==0)
         cout << "WARNING::getNumMCEvts::There are 0 events in the " << channelName << " MC." << endl;
      return value;
   } else{
      cout << "WARNING::getNumMCEvts::channelName " << channelName << " not recognized. Returning -1 event from MC." << endl << "Please check channel names." << endl;
      return -1.;
   }
}//getNumMCEvts


//______________________________________________________________________________
// leptonCan = The lepton category for which the processes are needed.
// intLum    = The integrated luminosity in pb-1 for this given lepton category 
vector<PhysicsProcessNEW*> getProcesses(DEFS::LeptonCat leptonCat, double intLum){
   
  string basePath = "/uscms/home/aperloff/nobackup/PS_outfiles_20121102_NTUPLES/";
  vector <PhysicsProcessNEW*> procs;
  
  procs.push_back(new PlotterPhysicsProcessNEW("WW",basePath+"WW.root", getCrossSection("WW"),
                                               intLum, getNumMCEvts("WW"), getProcessColor("WW"), "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("WZ",basePath+"WZ.root", getCrossSection("WZ"),
                                               intLum, getNumMCEvts("WZ"), getProcessColor("WZ"), "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("ZZ",basePath+"ZZ.root", getCrossSection("ZZ"),
                                               intLum, getNumMCEvts("ZZ"), getProcessColor("ZZ"), "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("DYJets",basePath+"DYJets.root", getCrossSection("DYJets"),
                                               intLum, getNumMCEvts("DYJets"), getProcessColor("DYJets"),
                                               "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("TTbar",basePath+"TTbar.root", getCrossSection("TTbar"),
                                               intLum, getNumMCEvts("TTbar"), getProcessColor("TTbar"),
                                               "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("STopT_T",basePath+"STopT_T.root", getCrossSection("STopT_T"),
                                               intLum, getNumMCEvts("STopT_T"), getProcessColor("STopT_T"),
                                               "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("STopT_Tbar",basePath+"STopT_Tbar.root",
                                               getCrossSection("STopT_Tbar"), intLum, getNumMCEvts("STopT_Tbar"),
                                               getProcessColor("STopT_Tbar"), "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("STopS_T",basePath+"STopS_T.root", getCrossSection("STopS_T"),
                                               intLum, getNumMCEvts("STopS_T"), getProcessColor("STopS_T"),
                                               "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("STopS_Tbar",basePath+"STopS_Tbar.root",
                                               getCrossSection("STopS_Tbar"), intLum, getNumMCEvts("STopS_Tbar"),
                                               getProcessColor("STopS_Tbar"), "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("STopTW_T",basePath+"STopTW_T.root",
                                               getCrossSection("STopTW_T"), intLum, getNumMCEvts("STopTW_T"),
                                               getProcessColor("STopTW_T"), "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("STopTW_Tbar",basePath+"STopTW_Tbar.root", 
                                               getCrossSection("STopTW_Tbar"), intLum, getNumMCEvts("STopTW_Tbar"),
                                               getProcessColor("STopTW_Tbar"), "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("ggH125",basePath+"ggH125.root", 
                                               getCrossSection("ggH125")*getBranchingRatio("ggH125"),
                                               intLum, getNumMCEvts("ggH125"), getProcessColor("ggH125"), 
                                               "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("qqH125",basePath+"qqH125.root", 
                                               getCrossSection("qqH125")*getBranchingRatio("qqH125"), 
                                               intLum, getNumMCEvts("ggH125"), getProcessColor("qqH125"), 
                                               "PS/jets2p"));
  procs.push_back(new PlotterPhysicsProcessNEW("WH125",basePath+"WH125.root", 
                                               getCrossSection("WH125")*getBranchingRatio("WH125"), 
                                               intLum, getNumMCEvts("WH125"), getProcessColor("WH125"), 
                                               "PS/jets2p"));

  // The Data
  // The normalization of the data is done in the same way as the MC. That is 
  // we multiply the histograms by getScaleFactor in class proc in PlotterAux.hh
  //  double getScaleFactor(){return sigma * intLum / initial_events;}
  // FOR DATA, to guarrantee the proper normalization, we set
  // sigma  = 1 / collected luminosity
  // initial_events = 1;
  
  if (leptonCat == DEFS::electron || leptonCat == DEFS::both){
     procs.push_back(new PlotterPhysicsProcessNEW("WJets_electron",basePath+"WJets.root",1.0*getCrossSection("WJets"),
     //procs.push_back(new PlotterPhysicsProcessNEW("WJets_electron","/uscms_data/d2/aperloff/PS_outfiles_WJets_20121115/WJets.root",1.0*getCrossSection("WJets"),
                                                 intLum, getNumMCEvts("WJets"), getProcessColor("WJets"),
                                                 "PS/jets2p", DEFS::electron));
    procs.push_back(new PlotterPhysicsProcessNEW("Data_electron",basePath+"SingleEl_Data_19p148fb.root",
                                                 1./intLum, intLum, 1,
                                                 getProcessColor("SingleEl_Data"),
                                                 "PS/jets2p", DEFS::electron)); 
    //procs.push_back(new PlotterPhysicsProcessNEW("QCD_electronEnriched", basePath+"QCD_ElEnriched.root",
    //                                             1.0*getCrossSection("QCD_ElEnriched"), intLum,
    //                                             getNumMCEvts("QCD_ElEnriched"),
    //                                             getProcessColor("QCD_ElEnriched"),
    //                                             "PS/jets2p", DEFS::electron));
  }  
  
  if (leptonCat == DEFS::muon || leptonCat == DEFS::both){
     procs.push_back(new PlotterPhysicsProcessNEW("WJets_muon",basePath+"WJets.root",1.0*getCrossSection("WJets"),
     //procs.push_back(new PlotterPhysicsProcessNEW("WJets_muon","/uscms_data/d2/aperloff/PS_outfiles_WJets_20121115/WJets.root",1.0*getCrossSection("WJets"),
                                                 intLum, getNumMCEvts("WJets"), getProcessColor("WJets"),
                                                 "PS/jets2p", DEFS::muon));
    procs.push_back(new PlotterPhysicsProcessNEW("Data_muon",basePath+"SingleMu_Data_19p279fb.root",
                                                 1./intLum, intLum, 1,
                                                 getProcessColor("SingleMu_Data"),
                                                 "PS/jets2p", DEFS::muon));
    //procs.push_back(new PlotterPhysicsProcessNEW("QCD_muonEnriched", basePath+"QCD_MuEnriched.root",
    //                                             1.0*getCrossSection("QCD_MuEnriched"), intLum,
    //                                             getNumMCEvts("QCD_MuEnriched"),
    //                                             getProcessColor("QCD_MuEnriched"),
    //                                             "PS/jets2p", DEFS::muon));
  }
  return procs;
  
}//getProcesses


//______________________________________________________________________________
vector<PhysicsProcessNEW*> getProcesses(vector<string>& processNames, DEFS::LeptonCat leptonCat, double intLum){
   
   string basePath = "/uscms/home/aperloff/nobackup/PS_outfiles_20120910_NTUPLES/";
   vector <PhysicsProcessNEW*> procs;

   for (unsigned int i=0; i<processNames.size(); i++){
      if (TString(processNames[i]).Contains("Data")==1) {
         if (leptonCat == DEFS::electron || leptonCat == DEFS::both)
            procs.push_back(new PlotterPhysicsProcessNEW("SingleEl_Data",basePath+processNames[i]+".root",
                                                         1.0/intLum, intLum, 1.0,
                                                         getProcessColor(processNames[i]),
                                                         "PS/jets2p",DEFS::electron));
         if (leptonCat == DEFS::muon || leptonCat == DEFS::both)
            procs.push_back(new PlotterPhysicsProcessNEW("SingleMu_Data",basePath+processNames[i]+".root",
                                                         1.0/intLum, intLum, 1.0,
                                                         getProcessColor(processNames[i]),
                                                         "PS/jets2p",DEFS::muon));
      }
      else if (TString(processNames[i]).Contains("WJets")==1) {
         if (leptonCat == DEFS::electron || leptonCat == DEFS::both)
            procs.push_back(new PlotterPhysicsProcessNEW("WJets_electron",basePath+processNames[i]+".root",
                                                         getCrossSection(processNames[i]), intLum, 
                                                         getNumMCEvts(processNames[i]),
                                                         getProcessColor(processNames[i]),
                                                         "PS/jets2p",DEFS::electron));
         if (leptonCat == DEFS::muon || leptonCat == DEFS::both)
            procs.push_back(new PlotterPhysicsProcessNEW("WJets_muon",basePath+processNames[i]+".root",
                                                         getCrossSection(processNames[i]), intLum,
                                                         getNumMCEvts(processNames[i]),
                                                         getProcessColor(processNames[i]),
                                                         "PS/jets2p",DEFS::muon));
      }
      else if (TString(processNames[i]).Contains("QCD")==1) {
         if (leptonCat == DEFS::electron || leptonCat == DEFS::both)
            procs.push_back(new PlotterPhysicsProcessNEW("QCD_ElEnriched", basePath+"QCD_ElEnriched.root",
                                                         getCrossSection("QCD_ElEnriched"), intLum,
                                                         getNumMCEvts("QCD_ElEnriched"),
                                                         getProcessColor("QCD_ElEnriched"),
                                                         "PS/jets2p", DEFS::electron));
         if (leptonCat == DEFS::muon || leptonCat == DEFS::both)
            procs.push_back(new PlotterPhysicsProcessNEW("QCD_MuEnriched", basePath+"QCD_MuEnriched.root",
                                                         getCrossSection("QCD_MuEnriched"), intLum,
                                                         getNumMCEvts("QCD_MuEnriched"),
                                                         getProcessColor("QCD_MuEnriched"),
                                                         "PS/jets2p", DEFS::muon));
      }
      else
         procs.push_back(new PlotterPhysicsProcessNEW(processNames[i],basePath+processNames[i]+".root",
                                                      getCrossSection(processNames[i]), intLum, 
                                                      getNumMCEvts(processNames[i]),
                                                      getProcessColor(processNames[i])));
   }

   return procs;

}//getProcesses


//______________________________________________________________________________
map<string, Plot*> getPlots(string leptonCatString, bool norm_data){

   map<string, Plot*> plots;

   FormattedPlot* a = new FormattedPlot;
   double pi = TMath::Pi();

   //Goes in the label and tells us whether we are looking at electrons or muons
   TString lep = TString(leptonCatString);

   double signalFactor = 5000;
   TString signalName = "ggH+WH+VBF";

   Double_t leptonptbinslow[9] = {20,25,30,35,40,50,70,100,1000};
   Double_t leptonptbinshigh[10] = {20,50,55,60,65,70,80,100,120,1000};
   Double_t jetptbinslow[9] = {20,25,30,35,40,50,70,100,1000};   
   Double_t jetptbinshigh[10] = {20,50,80,100,110,120,130,140,160,1000};
   Double_t Mjjbinslow[11] = {20,40,50,60,70,80,90,100,200,300,1000};
   Double_t Mjjbinshigh[11] = {20,100,110,120,130,150,180,200,250,350,1000};
   Double_t Mtbinshigh[11] = {0,20,40,60,70,80,100,120,140,200,1000};
   Double_t METbinshigh[10] = {0,40,50,60,70,80,100,120,200,1000};
   Double_t DRlepjet1low[10] = {0.3,0.5,0.9,1.3,1.5,1.7,2.0,2.5,3.0,5.0};
   Double_t DRlepjet1high[11] = {0.3,2.0,2.25,2.5,2.75,3.0,3.25,3.5,4.0,4.5,5.0};

   //---------------------ELECTRON-------------------//
   if (lep.CompareTo("electron") == 0 || lep.CompareTo("both") == 0)
   {
      TString electron = "electron";
      
      //goes in the label and tells us what cuts we are applying
      string cut = "_MET > 35, elPt > 30";
      TString cuts = TString(cut);
      
      a = new FormattedPlot;

      //a->templateHisto = new TH1D("Mjj_" + electron,"Mjj_" + electron +  cuts,500,0,500);
      //a->templateHisto = new TH1D("Mjj_" + electron,"Mjj_" + electron +  cuts,270,30,300);
      a->templateHisto = new TH1D("Mjj_" + electron,"Mjj_" + electron +  cuts,22,40,150);
      a->axisTitles.push_back("M_{jj} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      //a->range = make_pair(0.,400.);
      //a->range = make_pair(30.,300.);
      a->range = make_pair(40.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Mjj_electron"] = a;

      a = new FormattedPlot;

      a->templateHisto = new TH1D("MjjmWmT_" + electron,"MjjmWmT_" + electron +  cuts,500,-250,250);
      a->axisTitles.push_back("M_{jj} - WmT [GeV]");
      a->axisTitles.push_back("Number of Events / GeV");
      a->range = make_pair(-100.,250.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["MjjmWmT_electron"] = a;
   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("LeptPt_" + electron,"LeptPt_" + electron + cuts,1000,0,500);
      a->axisTitles.push_back("p_{T}^{lepton} [GeV]");
      a->axisTitles.push_back("Number of Events / 0.5 GeV");
      a->range = make_pair(20.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["LeptPt_electron"] = a;

   
      a = new FormattedPlot;
 
      a->templateHisto = new TH1D("LeptEta_" + electron,"LeptEta_" + electron + cuts,50,-5,5);
      a->axisTitles.push_back("#eta^{lepton} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.2 Radians");
      a->range = make_pair(-3.,3.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["LeptEta_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("LeptPhi_" + electron,"LeptPhi_" + electron + cuts,70,-3.5,3.5);
      a->axisTitles.push_back("#phi^{lepton} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(-3.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["LeptPhi_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("MET_" + electron,"MET_" + electron + cuts,1000,0,500);
      a->axisTitles.push_back("Missing E_{T} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(30.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["MET_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("WmT_" + electron,"WmT_" + electron + cuts,1000,0,500);
      a->axisTitles.push_back("M_{T}^{W} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(0.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["WmT_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Jet1Pt_" + electron,"Jet1Pt_" + electron + cuts,200,0,300);
      a->axisTitles.push_back("p_{T}^{jet_{1}} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(20.,200.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Jet1Pt_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Jet1Eta_" + electron,"Jet1Eta_" + electron + cuts,50,-5,5);
      a->axisTitles.push_back("#eta^{jet_{1}} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.2 Radians");
      a->range = make_pair(-3.,3.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Jet1Eta_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Jet1Phi_" + electron,"Jet1Phi_" + electron + cuts,70,-3.5,3.5);
      a->axisTitles.push_back("#phi^{jet_{1}} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(-3.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Jet1Phi_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Jet2Pt_" + electron,"Jet2Pt_" + electron + cuts,200,0,300);
      a->axisTitles.push_back("p_{T}^{jet_{2}} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(20.,100.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Jet2Pt_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Jet2Eta_" + electron,"Jet2Eta_" + electron + cuts,50,-5,5);
      a->axisTitles.push_back("#eta^{jet_{1}} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.2 Radians");
      a->range = make_pair(-3.,3.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Jet2Eta_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Jet2Phi_" + electron,"Jet2Phi_" + electron + cuts,70,-3.5,3.5);
      a->axisTitles.push_back("#phi^{jet_{1}} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(-3.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Jet2Phi_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("DeltaEtaJ1J2_" + electron,"DeltaEtaJ1J2_" + electron + cuts,50,0,5);
      a->axisTitles.push_back("#eta^{jet_{1}} - #eta^{jet_{2}} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(0.,5.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaEtaJ1J2_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Ptjj_" + electron,"Ptjj_" + electron + cuts,100,0,300);
      a->axisTitles.push_back("p_{T}^{jj} [GeV]");
      a->axisTitles.push_back("Number of Events / 10 GeV");
      a->range = make_pair(0.,250.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Ptjj_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("j1Pt_Mjj_" + electron,"j1Pt_Mjj_" + electron + cuts,500,0,5);
      a->axisTitles.push_back("p_{T}^{jet_{1}}/M_{jj} [GeV]");
      a->axisTitles.push_back("Number of Events / 0.01 GeV");
      a->range = make_pair(0.,2.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["j1Pt_Mjj_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("j2Pt_Mjj_" + electron,"j2Pt_Mjj_" + electron + cuts,500,0,5);
      a->axisTitles.push_back("p_{T}^{jet_{2}}/M_{jj} [GeV]");
      a->axisTitles.push_back("Number of Events / 0.01 GeV");
      a->range = make_pair(0.,1.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["j2Pt_Mjj_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Mlvjj_" + electron,"Mlvjj_" + electron + cuts,250,0,1000);
      a->axisTitles.push_back("M_{lvjj} [GeV]");
      a->axisTitles.push_back("Number of Events / 10 GeV");
      a->range = make_pair(50.,800.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Mlvjj_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("DeltaRLepMET_" + electron,"DeltaRLepMET_" + electron + cuts,50,0,10);
      a->axisTitles.push_back("#DeltaR of lep and MET [Radians]");
      a->axisTitles.push_back("Number of Events / 0.2 Radians");
      a->range = make_pair(0.,7.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaRLepMET_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("EJ1EJ2_" + electron,"EJ1EJ2_" + electron + cuts,500,700,5000);
      a->axisTitles.push_back("E_{J1} * E_{J2}  [GeV]");
      a->axisTitles.push_back("Number of Events / 10 GeV");
      a->range = make_pair(0.,5000.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["EJ1EJ2_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("BetaJ1BetaJ2_" + electron,"BetaJ1BetaJ2_" + electron + cuts,10000,0,10);
      a->axisTitles.push_back("#beta_{J1} * #beta_{J2} [GeV]");
      a->axisTitles.push_back("Number of Events / .01 GeV");
      a->range = make_pair(0.9,1.03);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["BetaJ1BetaJ2_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("DeltaRJ1J2_" + electron,"DeltaRJ1J2_" + electron + cuts,50,0,10);
      a->axisTitles.push_back("#DeltaR of J1 and J2 [Radians]");
      a->axisTitles.push_back("Number of Events / 0.2 Radians");
      a->range = make_pair(0.,7.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaRJ1J2_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("AngleJ1J2_" + electron,"AngleJ1J2_" + electron + cuts,50,0,5);
      a->axisTitles.push_back("Angle between J1 and J2 [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(-0.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["AngleJ1J2_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("jjlvPhi_" + electron,"jjlvPhi__" + electron + cuts,70,-3.5,3.5);
      a->axisTitles.push_back("#phi(jj) - #phi(lv)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-3.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["jjlvPhi_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("DeltaPhi_LJ1_" + electron,"DeltaPhi_LJ1__" + electron + cuts,50,-10,10);
      a->axisTitles.push_back("#Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaPhi_LJ1_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("DeltaPhi_J1J2_" + electron,"DeltaPhi_J1J2__" + electron + cuts,50,-10,10);
      a->axisTitles.push_back("#Delta #phi(J1J2)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaPhi_J1J2_electron"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("npv_" + electron,"npv__" + electron + cuts,100,0,100);
      a->axisTitles.push_back("npv");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,40);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["npv_electron"] = a;


      a = new FormattedPlot;

      a->templateHisto = new TH2D(("DeltaPhi_LJ1_vs_J1J2_" + electron),("DeltaPhi_LJ1_vs_J1J2__" + electron + cuts), 15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("#Delta #phi(J1J2) vs. #Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaPhi_LJ1_vs_J1J2_electron"] = a;

      a = new FormattedPlot;
   
      a->templateHisto = new TH2D(("MWjjVsMWlv_" + electron),("MWjjVsMWlv_" + electron + cuts),200,0,200,200,0,200);
      a->axisTitles.push_back("M_{W_{jj}}");
      a->axisTitles.push_back("M_{W_{l#nu}}");
      a->range = make_pair(0,200);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["MWjjVsMWlv_electron"] = a;


      a = new FormattedPlot;

      a->templateHisto = new TH2D(("DeltaPhi_LJ1_vs_J1J2_Positive_" + electron),("DeltaPhi_LJ1_vs_J1J2_Positive__" + electron + cuts),15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Positive #Delta #phi(J1J2) vs. #Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaPhi_LJ1_vs_J1J2_Positive_electron"] = a;


      a = new FormattedPlot;

      a->templateHisto = new TH2D(("DeltaPhi_LJ1_vs_J1J2_Negative_" + electron),("DeltaPhi_LJ1_vs_J1J2_Negative__" + electron + cuts),15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Negative #Delta #phi(J1J2) vs. #Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaPhi_LJ1_vs_J1J2_Negative_electron"] = a;


      a = new FormattedPlot;

      a->templateHisto = new TH2D(("DeltaPhi_LJ1_vs_J1J2_Subtracted_" + electron),("DeltaPhi_LJ1_vs_J1J2_Subtracted__" + electron + cuts),15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Subtracted #Delta #phi(J1J2) vs. #Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaPhi_LJ1_vs_J1J2_Subtracted_electron"] = a;

      //////////////test

      a = new FormattedPlot;

      a->templateHisto = new TH2D(("WmT_Negative_" + electron),("WmT_Negative__" + electron + cuts),15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Negative WmT");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["WmT_Negative_electron"] = a;

      a = new FormattedPlot;

      a->templateHisto = new TH2D(("WmT_Positive_" + electron),("WmT_Positive__" + electron + cuts),15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Positive WmT");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["WmT_Positive_electron"] = a;

      a = new FormattedPlot;

      a->templateHisto = new TH2D(("WmT_Subtracted_" + electron),("WmT_Subtracted__" + electron + cuts),15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Subtracted WmT");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::electron;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["WmT_Subtracted_electron"] = a;

      for (unsigned int tep=0; tep<15; tep++) {
         a = new FormattedPlot;
         TString name = TString(UserFunctions::concatString("tEventProb",tep)) + "_" + electron;
         TString title = name + cuts;
         string xaxis = UserFunctions::concatString("log(tEventProb[",tep) + "])";
         a->templateHisto = new TH1D(name,title,70,-45,0);
         a->axisTitles.push_back(xaxis);
         a->axisTitles.push_back("Number of Events");
         if (tep==14)
            a->range = make_pair(-45,0);
         else
            a->range = make_pair(-35,0);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = DEFS::both;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[string(name)] = a;
      }

      a = new FormattedPlot;
      
      TString name = "MVADiscriminator_" + electron;
      TString title = name + cuts;
      //a->templateHisto = new TH1D(name,title,2000,-10,10);
      a->templateHisto = new TH1D(name,title,70,-0.05,0.65);
      a->axisTitles.push_back("MVADiscriminator");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-0.05,0.65);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;
      
      a = new FormattedPlot;
      
      name = "MVAProbability_" + electron;
      title = name + cuts;
      a->templateHisto = new TH1D(name,title,100,0,1);
      a->axisTitles.push_back("Probability");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,1);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;

      a = new FormattedPlot;
      
      name = "epdPretagWWandWZ_" + electron;
      title = name + cuts;
      //a->templateHisto = new TH1D(name,title,400,-20,20);
      //a->templateHisto = new TH1D(name,title,110,-11.0,0.0);
      a->templateHisto = new TH1D(name,title,22,-11.0,0.0);
      a->axisTitles.push_back("log(epdPretagWWandWZ)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-11,0);
      a->logxy = make_pair(false,false);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;

      a = new FormattedPlot;
      
      name = "epdPretagWWandWZ_RE_" + electron;
      title = name + cuts;
      a->templateHisto = new TH1D(name,title,22,0,1);
      a->axisTitles.push_back("EPD(WW,WZ)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,1);
      a->logxy = make_pair(false,false);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;

      a = new FormattedPlot;

      name = "DeltaPhi_METJ1_" + electron;
      title = name + cuts;
      a->templateHisto = new TH1D(name,title,50,-10,10);
      a->axisTitles.push_back("#Delta #phi(MET,J1)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-pi,pi);
      a->logxy = make_pair(false,false);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;

      a = new FormattedPlot;
      /*
      name = "lpt_leta_j1pt_j1eta_j2pt_j2eta_" + electron;
      title = "lpt_leta_j1pt_j1eta_j2pt_j2eta_" + electron + cuts;
      a->templateHisto = new TProfileMDF("lpt_leta_j1pt_j1eta_j2pt_j2eta","lpt_leta_j1pt_j1eta_j2pt_j2eta_");
      ((TProfileMDF*)a->templateHisto)->AddAxis("lpt",8,leptonptbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("leta",10,0,2.5);
      ((TProfileMDF*)a->templateHisto)->AddAxis("j1pt",8,jetptbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("j1eta",10,0,2.5);
      ((TProfileMDF*)a->templateHisto)->AddAxis("j2pt",8,jetptbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("j2eta",10,0,2.5);
      ((TProfileMDF*)a->templateHisto)->Sumw2();
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;
      */
      /*
      name = "lpt_lpt_j1pt_j1pt_j2pt_j2pt_" + electron;
      title = "lpt_lpt_j1pt_j1pt_j2pt_j2pt_" + electron + cuts;
      a->templateHisto = new TProfileMDF("lpt_lpt_j1pt_j1pt_j2pt_j2pt","lpt_lpt_j1pt_j1pt_j2pt_j2pt_");
      ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{lepton} [GeV]",8,leptonptbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{lepton} [GeV]",9,leptonptbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{jet_{1}} [GeV]",8,jetptbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{jet_{1}} [GeV]",9,jetptbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{jet_{2}} [GeV]",8,jetptbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{jet_{2}} [GeV]",9,jetptbinshigh);
      ((TProfileMDF*)a->templateHisto)->Sumw2();
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;
      */
      name = "Mjj_Mjj_Mt_MET_DeltaR_DeltaR_" + electron;
      title = "Mjj_Mjj_Mt_MET_DeltaR_DeltaR_" + electron + cuts;
      a->templateHisto = new TProfileMDF("Mjj_Mjj_Mt_MET_DeltaR_DeltaR","Mjj_Mjj_Mt_MET_DeltaR_DeltaR_");
      ((TProfileMDF*)a->templateHisto)->AddAxis("M_{jj} [GeV]",10,Mjjbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("M_{jj} [GeV]",10,Mjjbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("M_{T}^{W} [GeV]",10,Mtbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("Missing E_{T} [GeV]",9,METbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("#DeltaR(#mu,jet1) [Radians]",9,DRlepjet1low);
      ((TProfileMDF*)a->templateHisto)->AddAxis("#DeltaR(#mu,jet1) [Radians]",10,DRlepjet1high);
      ((TProfileMDF*)a->templateHisto)->Sumw2();
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;
   }
   //----------------------------MUON------------------------------//
   if (lep.CompareTo("muon") == 0 || lep.CompareTo("both") == 0)
   {
      TString muon = "muon";
      
      //goes in the label and tells us what cuts we are applying
      string cut = "_MET > 30, muPt > 25";
      TString cuts = TString(cut);
      
      a = new FormattedPlot;
      
      //a->templateHisto = new TH1D("Mjj_" + muon,"Mjj_" + muon +  cuts,500,0,500);
      //a->templateHisto = new TH1D("Mjj_" + muon,"Mjj_" + muon +  cuts,270,30,300);
      a->templateHisto = new TH1D("Mjj_" + muon,"Mjj_" + muon +  cuts,22,40,150);
      a->axisTitles.push_back("M_{jj} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      //a->range = make_pair(0.,400.);
      //a->range = make_pair(30.,300.);
      a->range = make_pair(40.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Mjj_muon"] = a;


      a = new FormattedPlot;

      a->templateHisto = new TH1D("MjjmWmT_" + muon,"MjjmWmT_" + muon +  cuts,500,-250,250);
      a->axisTitles.push_back("M_{jj} - WmT [GeV]");
      a->axisTitles.push_back("Number of Events / GeV");
      a->range = make_pair(-100.,250.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["MjjmWmT_muon"] = a;
   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("LeptPt_" + muon,"LeptPt_" + muon + cuts,1000,0,500);
      a->axisTitles.push_back("p_{T}^{lepton} [GeV]");
      a->axisTitles.push_back("Number of Events / 0.5 GeV");
      a->range = make_pair(20.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["LeptPt_muon"] = a;

   
      a = new FormattedPlot;
 
      a->templateHisto = new TH1D("LeptEta_" + muon,"LeptEta_" + muon + cuts,50,-5,5);
      a->axisTitles.push_back("#eta^{lepton} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.2 Radians");
      a->range = make_pair(-3.,3.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["LeptEta_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("LeptPhi_" + muon,"LeptPhi_" + muon + cuts,70,-3.5,3.5);
      a->axisTitles.push_back("#phi^{lepton} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(-3.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["LeptPhi_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("MET_" + muon,"MET_" + muon + cuts,1000,0,500);
      a->axisTitles.push_back("Missing E_{T} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(30.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["MET_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("WmT_" + muon,"WmT_" + muon + cuts,1000,0,500);
      a->axisTitles.push_back("M_{T}^{W} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(0.,150.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["WmT_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Jet1Pt_" + muon,"Jet1Pt_" + muon + cuts,200,0,300);
      a->axisTitles.push_back("p_{T}^{jet_{1}} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(20.,200.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Jet1Pt_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Jet1Eta_" + muon,"Jet1Eta_" + muon + cuts,50,-5,5);
      a->axisTitles.push_back("#eta^{jet_{1}} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.2 Radians");
      a->range = make_pair(-3.,3.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Jet1Eta_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Jet1Phi_" + muon,"Jet1Phi_" + muon + cuts,70,-3.5,3.5);
      a->axisTitles.push_back("#phi^{jet_{1}} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(-3.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Jet1Phi_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Jet2Pt_" + muon,"Jet2Pt_" + muon + cuts,200,0,300);
      a->axisTitles.push_back("p_{T}^{jet_{2}} [GeV]");
      a->axisTitles.push_back("Number of Events / 5 GeV");
      a->range = make_pair(20.,100.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Jet2Pt_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Jet2Eta_" + muon,"Jet2Eta_" + muon + cuts,50,-5,5);
      a->axisTitles.push_back("#eta^{jet_{1}} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.2 Radians");
      a->range = make_pair(-3.,3.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Jet2Eta_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Jet2Phi_" + muon,"Jet2Phi_" + muon + cuts,70,-3.5,3.5);
      a->axisTitles.push_back("#phi^{jet_{1}} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(-3.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Jet2Phi_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("DeltaEtaJ1J2_" + muon,"DeltaEtaJ1J2_" + muon + cuts,50,0,5);
      a->axisTitles.push_back("#eta^{jet_{1}} - #eta^{jet_{2}} [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(0.,5.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaEtaJ1J2_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Ptjj_" + muon,"Ptjj_" + muon + cuts,100,0,300);
      a->axisTitles.push_back("p_{T}^{jj} [GeV]");
      a->axisTitles.push_back("Number of Events / 10 GeV");
      a->range = make_pair(0.,250.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Ptjj_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("j1Pt_Mjj_" + muon,"j1Pt_Mjj_" + muon + cuts,500,0,5);
      a->axisTitles.push_back("p_{T}^{jet_{1}}/M_{jj} [GeV]");
      a->axisTitles.push_back("Number of Events / 0.01 GeV");
      a->range = make_pair(0.,2.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["j1Pt_Mjj_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("j2Pt_Mjj_" + muon,"j2Pt_Mjj_" + muon + cuts,500,0,5);
      a->axisTitles.push_back("p_{T}^{jet_{2}}/M_{jj} [GeV]");
      a->axisTitles.push_back("Number of Events / 0.01 GeV");
      a->range = make_pair(0.,1.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["j2Pt_Mjj_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("Mlvjj_" + muon,"Mlvjj_" + muon + cuts,250,0,1000);
      a->axisTitles.push_back("M_{lvjj} [GeV]");
      a->axisTitles.push_back("Number of Events / 10 GeV");
      a->range = make_pair(50.,800.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["Mlvjj_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("DeltaRLepMET_" + muon,"DeltaRLepMET_" + muon + cuts,50,0,10);
      a->axisTitles.push_back("#DeltaR of lep and MET [Radians]");
      a->axisTitles.push_back("Number of Events / 0.2 Radians");
      a->range = make_pair(0.,7.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaRLepMET_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("EJ1EJ2_" + muon,"EJ1EJ2_" + muon + cuts,500,700,5000);
      a->axisTitles.push_back("E_{J1} * E_{J2}  [GeV]");
      a->axisTitles.push_back("Number of Events / 10 GeV");
      a->range = make_pair(0.,5000.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["EJ1EJ2_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("BetaJ1BetaJ2_" + muon,"BetaJ1BetaJ2_" + muon + cuts,10000,0,10);
      a->axisTitles.push_back("#beta_{J1} * #beta_{J2} [GeV]");
      a->axisTitles.push_back("Number of Events / .01 GeV");
      a->range = make_pair(0.9,1.03);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["BetaJ1BetaJ2_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("DeltaRJ1J2_" + muon,"DeltaRJ1J2_" + muon + cuts,50,0,10);
      a->axisTitles.push_back("#DeltaR of J1 and J2 [Radians]");
      a->axisTitles.push_back("Number of Events / 0.2 Radians");
      a->range = make_pair(0.,7.);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaRJ1J2_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("AngleJ1J2_" + muon,"AngleJ1J2_" + muon + cuts,50,0,5);
      a->axisTitles.push_back("Angle between J1 and J2 [Radians]");
      a->axisTitles.push_back("Number of Events / 0.1 Radians");
      a->range = make_pair(-0.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["AngleJ1J2_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("jjlvPhi_" + muon,"jjlvPhi__" + muon + cuts,70,-3.5,3.5);
      a->axisTitles.push_back("#phi(jj) - #phi(lv)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-3.5,3.5);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["jjlvPhi_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("DeltaPhi_LJ1_" + muon,"DeltaPhi_LJ1__" + muon + cuts,50,-10,10);
      a->axisTitles.push_back("#Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaPhi_LJ1_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("DeltaPhi_J1J2_" + muon,"DeltaPhi_J1J2__" + muon + cuts,50,-10,10);
      a->axisTitles.push_back("#Delta #phi(J1J2)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaPhi_J1J2_muon"] = a;

   
      a = new FormattedPlot;

      a->templateHisto = new TH1D("npv_" + muon,"npv__" + muon + cuts,100,0,100);
      a->axisTitles.push_back("npv");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,40);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["npv_muon"] = a;


      a = new FormattedPlot;

      a->templateHisto = new TH2D(("DeltaPhi_LJ1_vs_J1J2_" + muon),("DeltaPhi_LJ1_vs_J1J2__" + muon + cuts), 15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("#Delta #phi(J1J2) vs. #Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaPhi_LJ1_vs_J1J2_muon"] = a;

      a = new FormattedPlot;
   
      a->templateHisto = new TH2D(("MWjjVsMWlv_" + muon),("MWjjVsMWlv_" + muon + cuts),200,0,200,200,0,200);
      a->axisTitles.push_back("M_{W_{jj}}");
      a->axisTitles.push_back("M_{W_{l#nu}}");
      a->range = make_pair(0,200);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["MWjjVsMWlv_muon"] = a;


      a = new FormattedPlot;

      a->templateHisto = new TH2D(("DeltaPhi_LJ1_vs_J1J2_Positive_" + muon),("DeltaPhi_LJ1_vs_J1J2_Positive__" + muon + cuts),15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Positive #Delta #phi(J1J2) vs. #Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaPhi_LJ1_vs_J1J2_Positive_muon"] = a;


      a = new FormattedPlot;

      a->templateHisto = new TH2D(("DeltaPhi_LJ1_vs_J1J2_Negative_" + muon),("DeltaPhi_LJ1_vs_J1J2_Negative__" + muon + cuts),15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Negative #Delta #phi(J1J2) vs. #Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaPhi_LJ1_vs_J1J2_Negative_muon"] = a;


      a = new FormattedPlot;

      a->templateHisto = new TH2D(("DeltaPhi_LJ1_vs_J1J2_Subtracted_" + muon),("DeltaPhi_LJ1_vs_J1J2_Subtracted__" + muon + cuts),15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Subtracted #Delta #phi(J1J2) vs. #Delta #phi(LJ1)");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["DeltaPhi_LJ1_vs_J1J2_Subtracted_muon"] = a;

      //////////////test

      a = new FormattedPlot;

      a->templateHisto = new TH2D(("WmT_Negative_" + muon),("WmT_Negative__" + muon + cuts),15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Negative WmT");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["WmT_Negative_muon"] = a;

      a = new FormattedPlot;

      a->templateHisto = new TH2D(("WmT_Positive_" + muon),("WmT_Positive__" + muon + cuts),15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Positive WmT");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["WmT_Positive_muon"] = a;

      a = new FormattedPlot;

      a->templateHisto = new TH2D(("WmT_Subtracted_" + muon),("WmT_Subtracted__" + muon + cuts),15,-pi,pi,15,-pi,pi);
      a->axisTitles.push_back("Subtracted WmT");
      a->axisTitles.push_back("Number of Events / .2 Radians");
      a->range = make_pair(-pi,pi);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::muon;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots["WmT_Subtracted_muon"] = a;

      for (unsigned int tep=0; tep<15; tep++) {
         a = new FormattedPlot;
         TString name = TString(UserFunctions::concatString("tEventProb",tep)) + "_" + muon;
         TString title = name + cuts;
         string xaxis = UserFunctions::concatString("log(tEventProb[",tep) + "])";
         a->templateHisto = new TH1D(name,title,70,-45,0);
         a->axisTitles.push_back(xaxis);
         a->axisTitles.push_back("Number of Events");
         if (tep==14)
            a->range = make_pair(-45,0);
         else
            a->range = make_pair(-35,0);
         a->normToData = norm_data;
         a->stacked = true; a->leptonCat = DEFS::both;
         a->overlaySignalName = signalName;
         a->overlaySignalFactor = signalFactor;
         plots[string(name)] = a;
      }

      a = new FormattedPlot;
      
      TString name = "MVADiscriminator_" + muon;
      TString title = name + cuts;
      //a->templateHisto = new TH1D(name,title,2000,-10,10);
      a->templateHisto = new TH1D(name,title,70,-0.05,0.65);
      a->axisTitles.push_back("MVADiscriminator");
      a->axisTitles.push_back("Number of Events");
      //a->range = make_pair(-1,1);
      a->range = make_pair(-0.05,0.65);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;
      
      a = new FormattedPlot;
      
      name = "MVAProbability_" + muon;
      title = name + cuts;
      a->templateHisto = new TH1D(name,title,100,0,1);
      a->axisTitles.push_back("Probability");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,1);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;

      a = new FormattedPlot;
      
      name = "epdPretagWWandWZ_" + muon;
      title = name + cuts;
      //a->templateHisto = new TH1D(name,title,400,-20,20);
      //a->templateHisto = new TH1D(name,title,110,-11,0);
      a->templateHisto = new TH1D(name,title,22,-11,0);
      a->axisTitles.push_back("log(epdPretagWWandWZ)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-11,0);
      a->logxy = make_pair(false,false);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;

      a = new FormattedPlot;
      
      name = "epdPretagWWandWZ_RE_" + muon;
      title = name + cuts;
      a->templateHisto = new TH1D(name,title,22,0,1);
      a->axisTitles.push_back("EPD(WW,WZ)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(0,1);
      a->logxy = make_pair(false,false);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;

      a = new FormattedPlot;

      name = "DeltaPhi_METJ1_" + muon;
      title = name + cuts;
      a->templateHisto = new TH1D(name,title,50,-10,10);
      a->axisTitles.push_back("#Delta #phi(MET,J1)");
      a->axisTitles.push_back("Number of Events");
      a->range = make_pair(-pi,pi);
      a->logxy = make_pair(false,false);
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;

      a = new FormattedPlot;
      /*
      name = "lpt_leta_j1pt_j1eta_j2pt_j2eta_" + muon;
      title = "lpt_leta_j1pt_j1eta_j2pt_j2eta_" + muon + cuts;
      a->templateHisto = new TProfileMDF("lpt_leta_j1pt_j1eta_j2pt_j2eta","lpt_leta_j1pt_j1eta_j2pt_j2eta_");
      ((TProfileMDF*)a->templateHisto)->AddAxis("lpt",8,leptonptbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("leta",10,0,2.5);
      ((TProfileMDF*)a->templateHisto)->AddAxis("j1pt",8,jetptbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("j1eta",10,0,2.5);
      ((TProfileMDF*)a->templateHisto)->AddAxis("j2pt",8,jetptbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("j2eta",10,0,2.5);
      ((TProfileMDF*)a->templateHisto)->Sumw2();
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;
      */
      /*
      name = "lpt_lpt_j1pt_j1pt_j2pt_j2pt_" + muon;
      title = "lpt_lpt_j1pt_j1pt_j2pt_j2pt_" + muon + cuts;
      a->templateHisto = new TProfileMDF("lpt_lpt_j1pt_j1pt_j2pt_j2pt","lpt_lpt_j1pt_j1pt_j2pt_j2pt_");
      ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{lepton} [GeV]",8,leptonptbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{lepton} [GeV]",9,leptonptbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{jet_{1}} [GeV]",8,jetptbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{jet_{1}} [GeV]",9,jetptbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{jet_{2}} [GeV]",8,jetptbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("p_{T}^{jet_{2}} [GeV]",9,jetptbinshigh);
      ((TProfileMDF*)a->templateHisto)->Sumw2();
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;
      */
      name = "Mjj_Mjj_Mt_MET_DeltaR_DeltaR_" + muon;
      title = "Mjj_Mjj_Mt_MET_DeltaR_DeltaR_" + muon + cuts;
      a->templateHisto = new TProfileMDF("Mjj_Mjj_Mt_MET_DeltaR_DeltaR","Mjj_Mjj_Mt_MET_DeltaR_DeltaR_");
      ((TProfileMDF*)a->templateHisto)->AddAxis("M_{jj} [GeV]",10,Mjjbinslow);
      ((TProfileMDF*)a->templateHisto)->AddAxis("M_{jj} [GeV]",10,Mjjbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("M_{T}^{W} [GeV]",10,Mtbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("Missing E_{T} [GeV]",9,METbinshigh);
      ((TProfileMDF*)a->templateHisto)->AddAxis("#DeltaR(#mu,jet1) [Radians]",9,DRlepjet1low);
      ((TProfileMDF*)a->templateHisto)->AddAxis("#DeltaR(#mu,jet1) [Radians]",10,DRlepjet1high);
      ((TProfileMDF*)a->templateHisto)->Sumw2();
      a->normToData = norm_data;
      a->stacked = true; a->leptonCat = DEFS::both;
      a->overlaySignalName = signalName;
      a->overlaySignalFactor = signalFactor;
      plots[string(name)] = a;
   }

   // return all the plots to be made
   return plots;

}//getPlots


//______________________________________________________________________________
Color_t getProcessColor(TString channelName){

  if (channelName.CompareTo("WW") == 0)
    return kPink;
  else if (channelName.CompareTo("WZ") == 0)
    return kBlue;
  else if (channelName.CompareTo("ZZ") == 0)
    return kBlue-1;
  else if (channelName.CompareTo("WJets") == 0)
    return kTeal+2;
  else if (channelName.CompareTo("DYJets") == 0)
    return kPink-8;
  else if (channelName.CompareTo("QCD_ElEnriched") == 0)
    return kYellow+1;
  else if (channelName.CompareTo("QCD_MuEnriched") == 0)
    return kYellow+1;
  //right now these aren't being used because we are combining all the STop's in Plots.cc
  else if (channelName.CompareTo("STopT_T") == 0)
    return kOrange+1;
  else if (channelName.CompareTo("STopT_Tbar") == 0)
    return kOrange+1; //kCyan+3;
  else if (channelName.CompareTo("STopS_T") == 0)
    return kOrange+1; //kBlue;
  else if (channelName.CompareTo("STopS_Tbar") == 0)
    return kOrange+1; //kBlue+3;
  else if (channelName.CompareTo("STopTW_T") == 0)
    return kOrange+1; //kMagenta;
  else if (channelName.CompareTo("STopTW_Tbar") == 0)
    return kOrange+1; //kGreen+3;
  else if (channelName.CompareTo("STopTW_Tbar") == 0)
    return kOrange+1; //kGreen+3;
   else if (channelName.CompareTo("TTbar") == 0)
     return kAzure-2; 
   else if (channelName.CompareTo("ggH125") == 0)
     return kRed+2; 
   else if (channelName.CompareTo("qqH125") == 0)
     return kRed+2; 
   else if (channelName.CompareTo("WH125") == 0)
     return kRed+2; 
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
