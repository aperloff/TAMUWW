////////////////////////////////////////////////////////////////////////////////
//
// PerformSelectionCrab
// --------------------
//
//                         03/30/2012 Ilya Osipenkov <ilyao@fnal.gov> 
//                                    Alexx Perloff  <aperloff@physics.tamu.edu>
////////////////////////////////////////////////////////////////////////////////
// CMS 
// Stage 2 Skimming of the PAT-tuples to be used in Matrix Element calculations.
// WW/H->WW CrossSection Measurement using Matrix Element 
////////////////////////////////////////////////////////////////////////////////
// Designed to work with 4_2_8 PAT-tuples and tested on V3 HWW sync exercises.
// Use to skim the PATNtuples and create a custom made .root file containing the needed information.
////////////////////////////////////////////////////////////////////////////////

//
// CMS includes
//
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/LorentzVectorFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerPath.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Provenance/interface/RunAuxiliary.h"
#include "DataFormats/Provenance/interface/LuminosityBlockAuxiliary.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include <Math/VectorUtil.h>
#include "Muon/MuonAnalysisTools/interface/MuonEffectiveArea.h"
#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexSorter.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//
// User Defined Includes
//
#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"
#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/SpecialTools/interface/Table.hh"
#include "TAMUWW/SpecialTools/interface/TableRow.hh"
#include "TAMUWW/SpecialTools/interface/TableCellInt.hh"
#include "TAMUWW/SpecialTools/interface/TableCellVal.hh"
#include "TAMUWW/SpecialTools/interface/TableCellText.hh"
#include "TAMUWW/SpecialTools/interface/Value.hh"
#include "TAMUWW/Selection/bin/RunEventSet.h"

//
// ROOT includes
//
#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TRegexp.h"

//
// Standard Library Includes
//
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <memory>

//
// Namespace
//
using namespace std;
//using namespace edm;
//using namespace reco;

////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////
class PerformSelection : public edm::EDAnalyzer {
public:
   // construction/destruction
   explicit PerformSelection(const edm::ParameterSet& iConfig);
   virtual ~PerformSelection();
   
private:
   // member functions
   void beginJob();
   void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup);
   void endJob();
   // trigger functions
   void triggerSelection();
   // event functions
   void setFlags();
   // vertex functions
   void vertexSelection();
   // jet functions
   void jetClear();
   void jetSelection();
   void computeInvariantMass();
   double getPtD(vector<reco::PFCandidatePtr> pfCandidates);
   // b-tag functions
   // muon functions
   void muonClear();
   void muonSelection();
   // electron functions
   void eleClear();
   void eleSelection();
   bool CutBasedEID(vector< pat::Electron >::const_iterator elIter, bool tight = true);
   bool MVABasedEID(vector< pat::Electron >::const_iterator elIter, bool tight = true);
   bool eleIso(vector< pat::Electron >::const_iterator elIter, bool tight = true);
   bool eleVtx(vector< pat::Electron >::const_iterator elIter, bool tight = true);
   bool eleEP(vector< pat::Electron >::const_iterator elIter);
   // lepton functions
   // MET functions
   void metSelection();
   pair<double,double> getMETPhiCorrection_JetPhi(TString eraType);
   pair<double,double> getMETPhiCorrection_NPV(TString eraType);
   void doMETPhiCorrection(TString eraType_npv, TString eraType_jp);
   // additional (local) functions
   void getCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup);
   void makeTPUDist();
   void setDRlj();
   void setThetalj1pj2();
   void saveGenPart();
   bool matchGenParticles(const reco::Candidate* p1, const reco::Candidate* p2);
   /// matches jets to generator level objects
   pair<int, TLorentzVector> matchToGen(double eta, double phi);
   void fillJetMap(map<Int_t, Int_t> & jetMap, bool jets);
   vector<int> matchToGen_particleCollection(bool jets);
   double getJERfactor(double pt, double eta, double ptgen);
   void ptSort(vector<Jet> & vec);
   void swap(TLorentzVector & x, TLorentzVector & y);
   int  max_position(vector<Jet> & vec, int from, int to);
   int getBin(double x, const vector<double> boundaries, int length);
   /// increments the specified tables
   void incrementCounter(int nCut, unsigned int nJets, Table* t1, Table* t2, Table* t3 = 0);
   void printEventInformation(bool print, int cLevel, bool muon);
   void printJetInformation(); vector<string> jstreams; stringstream jstream;
   void printLeptonInformation(); vector<string> lstreams; stringstream lstream;

private:
   //
   // member data
   //
   // program variables
   bool Data;
   bool saveGenParticles;
   int  particleStatus;
   bool saveMETPhiPlots;
   bool noMETCut;
   bool invertEID;
   bool noMVAIsoCut;
   bool PFlowLoose;
   bool elONLY;
   bool muONLY;
   bool OnVsOffShell;
   bool StoreJets0;
   bool StoreJet1;
   int SQWaT_Version;
   bool doTrackerIso;
   bool doDetectorIso;
   bool doPFIso;
   bool doMVAeleSel;
   bool doJER;
   bool doMETPhi;
   bool doJESUncertainty;
   string JESUncertainty;
   string JESUncertaintyType;
   string JESUncertaintyFile;
   JetCorrectionUncertainty* jecUnc;
   double uncert;
   bool printEventInfo;
   bool printJetInfo;
   bool printLeptonInfo;
   // file variables
   TString outtablefilename;
   ofstream outtablefile;
   // table variables
   Table* tableEl;
   Table* tableMu;
   Table* tableLp;
   TString outtablenameEl;
   TString outtablenameMu;
   TString outtablenameLp;
   // tree variables
   TTree* EvtTree_0Jets;
   TTree* EvtTree_1Jets;
   TTree* EvtTree_2pJets;
   EventNtuple* EvtNtuple;
   // histogram variables
   TH1D* TPUDist;
   TH1D* PFIsoDist;
   //TH2D* PhotonIsoVsEta;
   //TH2D* NeutralHadronIsoVsEta;
   //TH2D* ChargedHadronIsoVsEta;
   //TH2D* EffectiveAreaVsEta;
   //TH2D* rhoPrimeVsEta;
   //TH2D* ptVsEta;
   TH1D* METPhi_BeforeCut;
   TH1D* METPhi_AfterCut;
   TH2D* METMagVsMETPhi_BeforeCut;
   TH2D* METMagVsMETPhi_AfterCut;
   TGraph* METxVsMETy_BeforeCut;
   TGraph* METxVsMETy_AfterCut;
   TH2D* METxVsNPV;
   TH2D* METyVsNPV;
   TH2D* METxVsNPV_0J;
   TH2D* METyVsNPV_0J;
   TH2D* METxVsNPV_1J;
   TH2D* METyVsNPV_1J;
   TH2D* METxVsJetPt;
   TH2D* METyVsJetPt;
   TH2D* METparaVsJPhi;
   TH2D* METperpVsJPhi;
   TH2D* METparaVsJPhi_1J;
   TH2D* METperpVsJPhi_1J;
   TH2D* METparaVsJPhi_1J_Pt[100];
   TH2D* METperpVsJPhi_1J_Pt[100];
   // handle InputTags/sources
   edm::InputTag triggerSource;
   edm::InputTag vtxSource;
   edm::InputTag genParticleSource;
   edm::InputTag pfJetSource;
   edm::InputTag electronSource;
   edm::InputTag muonSource;
   edm::InputTag METSource;
   edm::InputTag rhoSource;
   edm::InputTag pileupSource;
   // handles
   edm::Handle<pat::TriggerEvent> triggerHandle;
   edm::Handle<reco::VertexCollection> vtxHandle;
   edm::Handle<reco::GenParticleCollection> genParticleHandle;
   edm::Handle<vector<pat::Jet> > pfJetHandle;
   edm::Handle<vector<pat::Electron> > electronHandle;
   edm::Handle<vector<pat::Muon> > muonHandle;
   edm::Handle<vector<pat::MET> > METHandle;
   edm::Handle<double> rhoHandle;
   edm::Handle<vector<PileupSummaryInfo> > pileupHandle;
   // trigger variables
   bool mu_passTrigger;
   bool el_passTrigger;
   vector<string> muTrigger;
   vector<string> eleTrigger;
   bool MCpTrigger;
   bool passTrigger;
   map<string,bool> triggerMap;
   // event variables
   long runNumber;
   long eventNumber;
   long lumiNumber;
   long bxNumber;
   long orbitNumber;    
   long storeNumber; 
   long timeNumber;  
   long lumiSegmentNumber;
   int elcnt_Prim;
   int elcnt_Loose;
   int mucnt_Prim;
   int mucnt_Loose;
   int leptonCnt;
   int jcnt_tot;
   int jNEntries;
   int EvtTotCount;

   bool mu_passAll;
   bool mu_passStandard;
   bool mu_passFlag;
   bool el_passAll;
   bool el_passStandard;
   bool el_passFlag;
   // vertex variables
   reco::Vertex pv;
   bool PVfound;
   int vtxcnt;
   vector<Vertex> vp4;
   // jet variables
   vector <Jet> jp4;
   TVector3 jjp3;
   double j_pt;
   double j_ptMin;
   double j_eta;
   double j_aetaMax;
   double j_phi;
   double j_DRlepton;
   double j_DRelMin;
   double JERCor;
   double jesUncScale;
   double muPrim_DRjMin;
   double adphi;
   double CHEFMin;
   double CEMEFMax;
   double NHEFMax;    
   double NEMEFMax;
   unsigned int NDaughtersMin;
   int CMultiplicityMin;

   // b-tag variables
   int nBtagSSV, nBtagTC, nBtagCSV, nBtagCSVMVA;
   double bDiscriminatorSSVMin, bDiscriminatorTCMin, bDiscriminatorCSVMin, bDiscriminatorCSVMVAMin;

   // muon variables
   // Primary Muons (used in muon selection)
   double mu_pt;
   double mu_eta;
   double mu_phi;
   double mu_TrkIso;   
   double mu_DetIso;
   double mu_PFIso;

   double muPrim_ptMin;
   double muPrim_aetaMax;
   double muPrim_itNHits;
   double muPrim_dBMax;
   double muPrim_dzMax;
   double muPrim_TrkIsoMax;
   double muPrim_DetIsoMax;
   double muPrim_PFIsoMax;
   double muPrim_normChi2Max;
   int muPrim_nValidHitsMin;
   int muPrim_nValidMuonHitsMin;
   int muPrim_nValidPixelHitsMin;
   int muPrim_nMatchedStationsMin;

   bool muisGmTm;
   bool muisPFm;
   bool muisTightPrompt;
   bool mu_vtxPass;
   bool isProperlyIsolatedMu;
   double muDetIsoConeSize;
   double muPFIsoConeSize;

   // Loose Muons (used in both electron and muon selection)
   double muLoose_ptMin;
   double muLoose_aetaMax;
   double muLoose_TrkIsoMax;
   double muLoose_DetIsoMax;
   double muLoose_PFIsoMax;

   // electron variables
   // Primary Electrons (used in electron selection):
   double el_pt;
   double el_eta;
   double el_phi;
   double el_aetasc;
   double el_TrkIso;
   double el_DetIso;
   double el_PFIso;
   double el_mvaTrig;
   double el_mvaNonTrig;
   double elDetIsoConeSize;
   double elPFIsoConeSize;

   double elPrim_ptMin;
   double elPrim_aetaMax;
   double elPrim_aetascExcludeMax;
   double elPrim_aetascExcludeMin;
   double el_sigmaIetaIeta;
   double el_aDeltaPhi;
   double el_aDeltaEta;
   double el_HoE;
   // Equivalent to the EID cut (https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID)
   double elPrim_sigmaIetaIetaMaxEB;
   double elPrim_aDeltaPhiMaxEB;
   double elPrim_aDeltaEtaMaxEB;
   double elPrim_HoEMaxEB;
   double elPrim_sigmaIetaIetaMaxEE;
   double elPrim_aDeltaPhiMaxEE;
   double elPrim_aDeltaEtaMaxEE;
   double elPrim_HoEMaxEE;
   double elPrim_dBMax;
   double elPrim_dzMax;
   double elPrim_EinvMinusPinvMax;
   int elPrim_nMissingHits;
   double elPrim_TrkIsoMax;
   double elPrim_DetIsoMax;
   double elPrim_PFIsoMax;
   double elPrim_PFIsoMin_invertEID;

   // Loose Electrons
   double elLoose_ptMin;
   double elLoose_dBMax;
   double elLoose_dzMax;

   double elLoose_TrkIsoMax;
   double elLoose_DetIsoMax;
   double elLoose_PFIsoMax;
   double elLoose_sigmaIetaIetaMaxEB;
   double elLoose_aDeltaPhiMaxEB;
   double elLoose_aDeltaEtaMaxEB;
   double elLoose_HoEMaxEB;
   double elLoose_sigmaIetaIetaMaxEE;
   double elLoose_aDeltaPhiMaxEE;
   double elLoose_aDeltaEtaMaxEE;
   double elLoose_HoEMaxEE;

   // Electron Flags
   bool el_aetaPass;
   bool el_IsoPass;
   bool el_vtxPass;
   bool el_EinvMinusPinvPass;
   bool el_convPass;
   bool el_MVAEIDPass;
   bool el_CutEIDPass;
   bool el_PassTight;
   bool el_PassLoose;

   // lepton variables
   vector<Lepton> lp4;
   //vector<TLorentzVector> lp4;
   TVector3 lp3;
   // MET variables
   vector<MET> METp4;
   double MET_PtMin;
   bool MET_Pass;
   // additional variables
   double Mjj;
   int lQ;
   double lEta;
   //constants
   double etaBarrelMax;
   int NPtBins;
   vector<double> vpt;
};
