////////////////////////////////////////////////////////////////////////////////
//
// PerformSelectionCrab
// --------------------
//
//                         01/27/2012 Ilya Osipenkov <ilyao@fnal.gov> 
//                                    Alexx Perloff  <aperloff@physics.tamu.edu>
////////////////////////////////////////////////////////////////////////////////
// CMS 
// Stage 2 Skimming of the PAT-tuples to be used in Matrix Element calculations.
// WW/H->WW CrossSection Measurement using Matrix Element 
////////////////////////////////////////////////////////////////////////////////
// Designed to work with 4_2_8 PAT-tuples and tested on V3 HWW sync exercises.
// Use to skim the PATNtuples and create a custom made .root file containing the needed information.
////////////////////////////////////////////////////////////////////////////////

#include "TAMUWW/Selection/interface/PerformSelectionCrab.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
PerformSelection::PerformSelection(const edm::ParameterSet& iConfig)
   // Sample initialization
   //: s1(0) //An int defined in the header
   //, s2(false) //A bool defined in the header
   //, s3(0.1) //A double defined in the header
{
   //-----Input Tags For Handles
   triggerSource              =         iConfig.getParameter<edm::InputTag>   ("triggerSource");
   vtxSource                  =         iConfig.getParameter<edm::InputTag>   ("vtxSource");
   genParticleSource          =         iConfig.getParameter<edm::InputTag>   ("genParticleSource");
   pfJetSource                =         iConfig.getParameter<edm::InputTag>   ("pfJetSource");
   electronSource             =         iConfig.getParameter<edm::InputTag>   ("electronSource");
   muonSource                 =         iConfig.getParameter<edm::InputTag>   ("muonSource");
   METSource                  =         iConfig.getParameter<edm::InputTag>   ("METSource");
   rhoSource                  =         iConfig.getParameter<edm::InputTag>   ("rhoSource");
   pileupSource               =         iConfig.getParameter<edm::InputTag>   ("pileupSource");

   //-----Trigger Information
   muTrigger                  =         iConfig.getParameter<vector<string> > ("muTrigger");
   eleTrigger                 =         iConfig.getParameter<vector<string> > ("eleTrigger");

   //-----Program Level Inputs
   outtablefilename           = TString(iConfig.getParameter<string>          ("outtablefilename"));
   outtablenameEl             = TString(iConfig.getParameter<string>          ("outtablenameEl"));
   outtablenameMu             = TString(iConfig.getParameter<string>          ("outtablenameMu"));
   outtablenameLp             = TString(iConfig.getParameter<string>          ("outtablenameLp"));

   printEventInfo             =         iConfig.getParameter<bool>            ("printEventInfo");
   printJetInfo               =         iConfig.getParameter<bool>            ("printJetInfo");
   printLeptonInfo            =         iConfig.getParameter<bool>            ("printLeptonInfo");
   Data                       =         iConfig.getParameter<bool>            ("Data");
   saveGenParticles           =         iConfig.getParameter<bool>            ("saveGenParticles");
   particleStatus             =         iConfig.getParameter<int>             ("particleStatus");
   saveMETPhiPlots            =         iConfig.getParameter<bool>            ("saveMETPhiPlots");
   noMETCut                   =         iConfig.getParameter<bool>            ("noMETCut");
   invertEID                  =         iConfig.getParameter<bool>            ("invertEID");
   PFlowLoose                 =         iConfig.getParameter<bool>            ("PFlowLoose");
   elONLY                     =         iConfig.getParameter<bool>            ("elONLY");
   muONLY                     =         iConfig.getParameter<bool>            ("muONLY");
   OnVsOffShell               =         iConfig.getParameter<bool>            ("OnVsOffShell");
   StoreJets0                 =         iConfig.getParameter<bool>            ("StoreJets0");
   StoreJet1                  =         iConfig.getParameter<bool>            ("StoreJet1");

   MCpTrigger                 =         iConfig.getParameter<bool>            ("MCpTrigger");
   SQWaT_Version              =         iConfig.getParameter<int>             ("SQWaT_Version");
   doTrackerIso               =         iConfig.getParameter<bool>            ("doTrackerIso");
   doDetectorIso              =         iConfig.getParameter<bool>            ("doDetectorIso");
   doPFIso                    =         iConfig.getParameter<bool>            ("doPFIso");
   doMVAeleSel                =         iConfig.getParameter<bool>            ("doMVAeleSel");
   noMVAIsoCut                =         iConfig.getParameter<bool>            ("noMVAIsoCut");
   doJER                      =         iConfig.getParameter<bool>            ("doJER");
   doMETPhi                   =         iConfig.getParameter<bool>            ("doMETPhi");
   doJESUncertainty           =         iConfig.getParameter<bool>            ("doJESUncertainty");
   JESUncertainty             =         iConfig.getParameter<string>          ("JESUncertainty");
   JESUncertaintyType         =         iConfig.getParameter<string>          ("JESUncertaintyType");
   JESUncertaintyFile         =         iConfig.getParameter<string>          ("JESUncertaintyFile");

   //-----Event Variable Inputs
   elcnt_Prim                 =         iConfig.getParameter<int>             ("elcnt_Prim");
   elcnt_Loose                =         iConfig.getParameter<int>             ("elcnt_Loose");
   mucnt_Prim                 =         iConfig.getParameter<int>             ("mucnt_Prim");
   mucnt_Loose                =         iConfig.getParameter<int>             ("mucnt_Loose");
   leptonCnt                  =         iConfig.getParameter<int>             ("leptonCnt");
   jcnt_tot                   =         iConfig.getParameter<int>             ("jcnt_tot");
   EvtTotCount                =         iConfig.getParameter<int>             ("EvtTotCount");
   mu_passAll                 =         iConfig.getParameter<bool>            ("mu_passAll");
   mu_passStandard            =         iConfig.getParameter<bool>            ("mu_passStandard");
   mu_passFlag                =         iConfig.getParameter<bool>            ("mu_passFlag");
   el_passAll                 =         iConfig.getParameter<bool>            ("el_passAll");
   el_passStandard            =         iConfig.getParameter<bool>            ("el_passStandard");
   el_passFlag                =         iConfig.getParameter<bool>            ("el_passFlag");

   //-----Trigger Variable Inputs
   mu_passTrigger             =         iConfig.getParameter<bool>            ("mu_passTrigger");
   el_passTrigger             =         iConfig.getParameter<bool>            ("el_passTrigger");

   //-----Vertex Variable Inputs
   PVfound                    =         iConfig.getParameter<bool>            ("PVfound");
   vtxcnt                     =         iConfig.getParameter<int>             ("vtxcnt");

   //-----Jet Variable Inputs
   j_ptMin                    =         iConfig.getParameter<double>          ("j_ptMin");
   j_aetaMax                  =         iConfig.getParameter<double>          ("j_aetaMax");
   j_DRelMin                  =         iConfig.getParameter<double>          ("j_DRelMin");
   muPrim_DRjMin              =         iConfig.getParameter<double>          ("muPrim_DRjMin");
   CHEFMin                    =         iConfig.getParameter<double>          ("CHEFMin");
   CEMEFMax                   =         iConfig.getParameter<double>          ("CEMEFMax");
   NHEFMax                    =         iConfig.getParameter<double>          ("NHEFMax");
   NEMEFMax                   =         iConfig.getParameter<double>          ("NEMEFMax");
   NDaughtersMin              =         iConfig.getParameter<unsigned int>    ("NDaughtersMin");
   CMultiplicityMin           =         iConfig.getParameter<int>             ("CMultiplicityMin");

   //-----B-Tag Variable Inputs
   nBtagSSV                   =         iConfig.getParameter<int>             ("nBtagSSV");
   nBtagTC                    =         iConfig.getParameter<int>             ("nBtagTC");
   bDiscriminatorSSVMin       =         iConfig.getParameter<double>          ("bDiscriminatorSSVMin");
   bDiscriminatorTCMin        =         iConfig.getParameter<double>          ("bDiscriminatorTCMin");
   bDiscriminatorCSVMin       =         iConfig.getParameter<double>          ("bDiscriminatorCSVMin");
   bDiscriminatorCSVMVAMin    =         iConfig.getParameter<double>          ("bDiscriminatorCSVMVAMin");

   //-----Muon Variable Inputs
   muPrim_ptMin               =         iConfig.getParameter<double>          ("muPrim_ptMin");
   muPrim_aetaMax             =         iConfig.getParameter<double>          ("muPrim_aetaMax");
   muPrim_dBMax               =         iConfig.getParameter<double>          ("muPrim_dBMax");
   muPrim_dzMax               =         iConfig.getParameter<double>          ("muPrim_dzMax");
   muPrim_TrkIsoMax           =         iConfig.getParameter<double>          ("muPrim_TrkIsoMax");
   muPrim_DetIsoMax           =         iConfig.getParameter<double>          ("muPrim_DetIsoMax");
   muPrim_PFIsoMax            =         iConfig.getParameter<double>          ("muPrim_PFIsoMax");
   muPrim_normChi2Max         =         iConfig.getParameter<double>          ("muPrim_normChi2Max");
   muPrim_nValidHitsMin       =         iConfig.getParameter<int>             ("muPrim_nValidHitsMin");
   muPrim_nValidMuonHitsMin   =         iConfig.getParameter<int>             ("muPrim_nValidMuonHitsMin");
   muPrim_nValidPixelHitsMin  =         iConfig.getParameter<int>             ("muPrim_nValidPixelHitsMin");
   muPrim_nMatchedStationsMin =         iConfig.getParameter<int>             ("muPrim_nMatchedStationsMin");
   isProperlyIsolatedMu       =         iConfig.getParameter<bool>            ("isProperlyIsolatedMu");
   muDetIsoConeSize           =         iConfig.getParameter<double>          ("muDetIsoConeSize");
   muPFIsoConeSize            =         iConfig.getParameter<double>          ("muPFIsoConeSize");
   muLoose_ptMin              =         iConfig.getParameter<double>          ("muLoose_ptMin");
   muLoose_aetaMax            =         iConfig.getParameter<double>          ("muLoose_aetaMax");
   muLoose_TrkIsoMax          =         iConfig.getParameter<double>          ("muLoose_TrkIsoMax");
   muLoose_DetIsoMax          =         iConfig.getParameter<double>          ("muLoose_DetIsoMax");
   muLoose_PFIsoMax           =         iConfig.getParameter<double>          ("muLoose_PFIsoMax");

   //-----Electron Variable Inputs
   elPrim_ptMin               =         iConfig.getParameter<double>          ("elPrim_ptMin");
   elPrim_aetaMax             =         iConfig.getParameter<double>          ("elPrim_aetaMax");
   elPrim_aetascExcludeMax    =         iConfig.getParameter<double>          ("elPrim_aetascExcludeMax");
   elPrim_aetascExcludeMin    =         iConfig.getParameter<double>          ("elPrim_aetascExcludeMin");
   elPrim_sigmaIetaIetaMaxEB  =         iConfig.getParameter<double>          ("elPrim_sigmaIetaIetaMaxEB");
   elPrim_aDeltaPhiMaxEB      =         iConfig.getParameter<double>          ("elPrim_aDeltaPhiMaxEB");
   elPrim_aDeltaEtaMaxEB      =         iConfig.getParameter<double>          ("elPrim_aDeltaEtaMaxEB");
   elPrim_HoEMaxEB            =         iConfig.getParameter<double>          ("elPrim_HoEMaxEB"); 
   elPrim_sigmaIetaIetaMaxEE  =         iConfig.getParameter<double>          ("elPrim_sigmaIetaIetaMaxEE");
   elPrim_aDeltaPhiMaxEE      =         iConfig.getParameter<double>          ("elPrim_aDeltaPhiMaxEE");
   elPrim_aDeltaEtaMaxEE      =         iConfig.getParameter<double>          ("elPrim_aDeltaEtaMaxEE");
   elPrim_HoEMaxEE            =         iConfig.getParameter<double>          ("elPrim_HoEMaxEE");
   elPrim_dBMax               =         iConfig.getParameter<double>          ("elPrim_dBMax");
   elPrim_dzMax               =         iConfig.getParameter<double>          ("elPrim_dzMax");
   elPrim_EinvMinusPinvMax    =         iConfig.getParameter<double>          ("elPrim_EinvMinusPinvMax");
   elPrim_TrkIsoMax           =         iConfig.getParameter<double>          ("elPrim_TrkIsoMax");
   elPrim_DetIsoMax           =         iConfig.getParameter<double>          ("elPrim_DetIsoMax");
   elPrim_PFIsoMax            =         iConfig.getParameter<double>          ("elPrim_PFIsoMax");
   elPrim_PFIsoMin_invertEID  =         iConfig.getParameter<double>          ("elPrim_PFIsoMin_invertEID");

   elLoose_ptMin              =         iConfig.getParameter<double>          ("elLoose_ptMin");
   elLoose_sigmaIetaIetaMaxEB =         iConfig.getParameter<double>          ("elLoose_sigmaIetaIetaMaxEB");
   elLoose_aDeltaPhiMaxEB     =         iConfig.getParameter<double>          ("elLoose_aDeltaPhiMaxEB");
   elLoose_aDeltaEtaMaxEB     =         iConfig.getParameter<double>          ("elLoose_aDeltaEtaMaxEB");
   elLoose_HoEMaxEB           =         iConfig.getParameter<double>          ("elLoose_HoEMaxEB");
   elLoose_sigmaIetaIetaMaxEE =         iConfig.getParameter<double>          ("elLoose_sigmaIetaIetaMaxEE");
   elLoose_aDeltaPhiMaxEE     =         iConfig.getParameter<double>          ("elLoose_aDeltaPhiMaxEE");
   elLoose_aDeltaEtaMaxEE     =         iConfig.getParameter<double>          ("elLoose_aDeltaEtaMaxEE");
   elLoose_HoEMaxEE           =         iConfig.getParameter<double>          ("elLoose_HoEMaxEE"); 
   elLoose_dBMax              =         iConfig.getParameter<double>          ("elLoose_dBMax");
   elLoose_dzMax              =         iConfig.getParameter<double>          ("elLoose_dzMax");
   elLoose_TrkIsoMax          =         iConfig.getParameter<double>          ("elLoose_TrkIsoMax");
   elLoose_DetIsoMax          =         iConfig.getParameter<double>          ("elLoose_DetIsoMax");
   elLoose_PFIsoMax           =         iConfig.getParameter<double>          ("elLoose_PFIsoMax");
   
   elDetIsoConeSize           =         iConfig.getParameter<double>          ("elDetIsoConeSize");
   elPFIsoConeSize            =         iConfig.getParameter<double>          ("elPFIsoConeSize");

   el_aetaPass                =         iConfig.getParameter<bool>            ("el_aetaPass");

   //-----Lepton Variable Inputs

   //-----MET Variable Inputs
   if (invertEID)
      MET_PtMin               =         iConfig.getParameter<double>          ("MET_PtMin_invertEID");
   else
      MET_PtMin               =         iConfig.getParameter<double>          ("MET_PtMin");

   //-----Additional Variable Inputs
   Mjj                        =         iConfig.getParameter<double>          ("Mjj");
   lQ                         =         iConfig.getParameter<int>             ("lQ");
   lEta                       =         iConfig.getParameter<double>          ("lEta");

   //-----Constant Inputs
   etaBarrelMax               =         iConfig.getParameter<double>          ("etaBarrelMax");
   NPtBins                    =         iConfig.getParameter<int>             ("NPtBins");
   vpt                        =         iConfig.getParameter<vector<double> > ("vpt");

   if ( iConfig.exists("JESUncertainty") ){
      if (JESUncertainty == "up" || JESUncertainty == "down") {
         jecUnc =  new JetCorrectionUncertainty(*(new JetCorrectorParameters(JESUncertaintyFile, JESUncertaintyType)));
      }
   }
}


//______________________________________________________________________________
PerformSelection::~PerformSelection() {}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void PerformSelection::beginJob() {
   edm::Service<TFileService> fs;
   if (!fs) throw edm::Exception(edm::errors::Configuration,
                                 "TFileService missing from configuration!");
   
   vector<string> cuts;
   for (unsigned int c = 0; c<DEFS::nCutLevel; c++) {
      cuts.push_back(DEFS::getCutLevelString((DEFS::CutLevel)c));
   }
   vector<string> jets;
   for (unsigned int j = 0; j<DEFS::NJETS; j++) {
      jets.push_back(DEFS::getJetBinString((DEFS::JetBin)j));
   }

   // create the tree
   tableEl = fs->make<Table>(string(outtablenameEl),cuts,jets,"TableCellVal");
   tableMu = fs->make<Table>(string(outtablenameMu),cuts,jets,"TableCellVal");
   tableLp = fs->make<Table>(string(outtablenameLp),cuts,jets,"TableCellVal");
   EvtNtuple = new EventNtuple();
   if (StoreJets0) {
      EvtTree_0Jets = fs->make<TTree>(getJetBinString(DEFS::jets0).c_str(), "Output tree");
      EvtTree_0Jets->Branch("EvtNtuple", "EventNtuple", &EvtNtuple);
   }
   if (StoreJet1) {
      EvtTree_1Jets = fs->make<TTree>(getJetBinString(DEFS::jet1).c_str(), "Output tree");
      EvtTree_1Jets->Branch("EvtNtuple", "EventNtuple", &EvtNtuple);
   }
   EvtTree_2pJets = fs->make<TTree>((getJetBinString(DEFS::jets2)+"p").c_str(), "Output tree for matrix element");
   EvtTree_2pJets->Branch("EvtNtuple", "EventNtuple", &EvtNtuple);
   if (!Data)
      TPUDist = fs->make<TH1D>("TPUDist","TPUDist",600,0,60);
   PFIsoDist = fs->make<TH1D>("PFIsoDist","PFIsoDist",1000,0,2);
   if (saveMETPhiPlots) {
      METPhi_BeforeCut = fs->make<TH1D>("METPhi_BeforeCut","METPhi_BeforeCut",63,-TMath::Pi(),TMath::Pi());
      METPhi_AfterCut = fs->make<TH1D>("METPhi_AfterCut","METPhi_AfterCut",63,-TMath::Pi(),TMath::Pi());
      METMagVsMETPhi_BeforeCut = fs->make<TH2D>("METMagVsMETPhi_BeforeCut","METMagVsMETPhi_BeforeCut",63,-TMath::Pi(),TMath::Pi(),2000,-1000,1000);
      METMagVsMETPhi_AfterCut = fs->make<TH2D>("METMagVsMETPhi_AfterCut","METMagVsMETPhi_AfterCut",63,-TMath::Pi(),TMath::Pi(),2000,-1000,1000);
      METxVsMETy_BeforeCut = fs->make<TGraph>();
      METxVsMETy_AfterCut = fs->make<TGraph>();
      METxVsNPV = fs->make<TH2D>("METxVsNPV","METxVsNPV",60,0,60,2000,-1000,1000);
      METyVsNPV = fs->make<TH2D>("METyVsNPV","METyVsNPV",60,0,60,2000,-1000,1000);
      METxVsNPV->Sumw2();
      METyVsNPV->Sumw2();
      METxVsNPV_0J = fs->make<TH2D>("METxVsNPV_0J","METxVsNPV_0J",60,0,60,2000,-1000,1000);
      METyVsNPV_0J = fs->make<TH2D>("METyVsNPV_0J","METyVsNPV_0J",60,0,60,2000,-1000,1000);
      METxVsNPV_0J->Sumw2();
      METyVsNPV_0J->Sumw2();
      METxVsNPV_1J = fs->make<TH2D>("METxVsNPV_1J","METxVsNPV_1J",60,0,60,2000,-1000,1000);
      METyVsNPV_1J = fs->make<TH2D>("METyVsNPV_1J","METyVsNPV_1J",60,0,60,2000,-1000,1000);
      METxVsNPV_1J->Sumw2();
      METyVsNPV_1J->Sumw2();
      METxVsJetPt = fs->make<TH2D>("METxVsJetPt","METxVsJetPt",200,0,1000,2000,-1000,1000);
      METyVsJetPt = fs->make<TH2D>("METyVsJetPt","METyVsJetPt",200,0,1000,2000,-1000,1000);
      METxVsJetPt->Sumw2();
      METyVsJetPt->Sumw2();
      METparaVsJPhi = fs->make<TH2D>("METparaVsJPhi","METparaVsJPhi",80,-4,4,2000,-1000,1000);
      METperpVsJPhi = fs->make<TH2D>("METperpVsJPhi","METperpVsJPhi",80,-4,4,2000,-1000,1000);
      METparaVsJPhi->Sumw2();
      METperpVsJPhi->Sumw2();
      METparaVsJPhi_1J = fs->make<TH2D>("METparaVsJPhi_1J","METparaVsJPhi_1J",80,-4,4,2000,-1000,1000);
      METperpVsJPhi_1J = fs->make<TH2D>("METperpVsJPhi_1J","METperpVsJPhi_1J",80,-4,4,2000,-1000,1000);
      METparaVsJPhi_1J->Sumw2();
      METperpVsJPhi_1J->Sumw2();
      for(int i=0; i<NPtBins; i++) {
         METparaVsJPhi_1J_Pt[i] = fs->make<TH2D>(Form("METparaVsJPhi_1J_%.0fto%.0fGeV",vpt[i],vpt[i+1]),Form("METparaVsJPhi_1J_%.0fto%.0fGeV",vpt[i],vpt[i+1]),80,-4,4,2000,-1000,1000);
         METparaVsJPhi_1J_Pt[i]->Sumw2();
         METperpVsJPhi_1J_Pt[i] = fs->make<TH2D>(Form("METperpVsJPhi_1J_%.0fto%.0fGeV",vpt[i],vpt[i+1]),Form("METperpVsJPhi_1J_%.0fto%.0fGeV",vpt[i],vpt[i+1]),80,-4,4,2000,-1000,1000);
         METperpVsJPhi_1J_Pt[i]->Sumw2();
      }
   }

/*
   PhotonIsoVsEta = fs->make<TH2D>("PhotonIsoVsEta","PhotonIsoVsEta",1000,0,2,200,-5,5);
   NeutralHadronIsoVsEta = fs->make<TH2D>("NeutralHadronIsoVsEta","NeutralHadronIsoVsEta",1000,0,2,200,-5,5);
   ChargedHadronIsoVsEta = fs->make<TH2D>("ChargedHadronIsoVsEta","ChargedHadronIsoVsEta",1000,0,2,200,-5,5);
   EffectiveAreaVsEta = fs->make<TH2D>("EffectiveAreaVsEta","EffectiveAreaVsEta",1000,0,2,200,-5,5);
   rhoPrimeVsEta = fs->make<TH2D>("rhoPrimeVsEta","rhoPrimeVsEta",1000,0,2,200,-5,5);
   ptVsEta = fs->make<TH2D>("ptVsEta","ptVsEta",400,0,2000,200,-5,5);
*/
}


//______________________________________________________________________________
void PerformSelection::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
   
   //
   //find the runNumber and eventNumber here
   //
   eventNumber       = iEvent.id().event();
   runNumber         = iEvent.id().run();
   lumiNumber        = iEvent.id().luminosityBlock();
   bxNumber          = iEvent.eventAuxiliary().bunchCrossing();
   orbitNumber       = iEvent.eventAuxiliary().orbitNumber();
   storeNumber       = iEvent.eventAuxiliary().storeNumber();
   timeNumber        = iEvent.eventAuxiliary().time().value();
   lumiSegmentNumber = iEvent.eventAuxiliary().luminosityBlock();
      
   //if (eventNumber!=2535884 && eventNumber!=3093466 && eventNumber!=9993322 && eventNumber!=3145918 && eventNumber!=3146025 && eventNumber!=4938770 && eventNumber!=4464867)
   //if (eventNumber!=513862 && eventNumber!=515700 && eventNumber!=516645 && eventNumber!=527553 && eventNumber!=629459)
   //return;

   //
   // Create a set containing all the run and event numbers
   // so we can keep track of duplicate events.
   //
   RunEventSet runEventSet;
   //
   //Check and remove duplicates here
   //
   if (runEventSet.alreadySeen(runNumber, eventNumber)) {
      cout << "WARNING, event Run = " << runNumber << ", Event = " << eventNumber
           <<" is duplicated" << endl << "We will skip this event." << endl;
      return;
   }
      
   //
   // Load the objects of interest
   //
   getCollections(iEvent, iSetup);

   //
   // Make the distribution of the true pileup interactions
   //
   if (!Data)
      makeTPUDist();

   //
   //Perform The Selection
   //

   //
   // Trigger Selection
   //
   triggerSelection();

   //
   // Vertex Selection
   //
   vertexSelection();

   //
   // Electron Selection
   //
   eleClear();
   eleSelection();

   //
   // Muon Selection (No Jet Cleaning)
   //
   muonClear();
   muonSelection();

   //
   // MET Selection
   //
   metSelection();

   //
   // Jet Selection (without muon cleaning)
   // Electron Cleaning (i.e. remove the jets which are likely to be electorns using the DR cut).
   // Using PF Jets
   //
   jetClear();
   jetSelection();
   // Make sure the jets are sorted by pt
   if (jp4.size()>1)
      ptSort(jp4);

   if(doMETPhi) {
      if(Data) {
         //doMETPhiCorrection("pfMEtSysShiftCorrParameters_2012runAvsNvtx_TAMUWW_data","pfMEtSysShiftCorrParameters_2012runAvsJPhi_TAMUWW_data");
         doMETPhiCorrection("pfMEtSysShiftCorrParameters_2012runAvsNvtx_TAMUWW_data","");
      }
      else {
         //doMETPhiCorrection("pfMEtSysShiftCorrParameters_2012runAvsNvtx_TAMUWW_mc","pfMEtSysShiftCorrParameters_2012runAvsJPhi_TAMUWW_mc");
         doMETPhiCorrection("pfMEtSysShiftCorrParameters_2012runAvsNvtx_TAMUWW_mc","");
      }
   }

   //
   // Set any further quantities that require all of the physics objects to be selected first
   //
   // Compute the invariant mass, when only two jets are present
   computeInvariantMass();
   setDRlj();
   setThetalj1pj2();

   //
   // Record The Results
   //
      
   //
   // Record the number of events passing the cuts for each jet count
   //
   el_passStandard=false;
   el_passAll=false;
   mu_passStandard=false;
   mu_passAll=false;

   printEventInformation(printEventInfo, -1, true);
   for(unsigned int Nj=0; Nj<DEFS::NJETS; Nj++) {
      if((unsigned int)jcnt_tot==Nj || (Nj==(DEFS::NJETS-1) && (unsigned int)jcnt_tot>(DEFS::NJETS-1))) {
         incrementCounter(1,Nj,tableEl,tableMu,tableLp);
         printEventInformation(printEventInfo, 0, true);
         
         if (el_passTrigger) {
            incrementCounter(2,Nj,tableEl,tableLp);
            printEventInformation(printEventInfo, 1, true);

            if (PVfound) {
               incrementCounter(3,Nj,tableEl,tableLp);
               printEventInformation(printEventInfo, 2, true);
               //
               // Record the electron info
               //
               if (elcnt_Prim>=leptonCnt) {
                  incrementCounter(4,Nj,tableEl,tableLp);
                  if (!muONLY) printEventInformation(printEventInfo, 3, false);
                  if (mucnt_Loose==0) {
                     incrementCounter(5,Nj,tableEl,tableLp);
                     if (!muONLY) printEventInformation(printEventInfo, 4, false);
                     if (elcnt_Loose==leptonCnt) {
                        incrementCounter(6,Nj,tableEl,tableLp);
                        if (!muONLY) printEventInformation(printEventInfo, 5, false);
                        el_passStandard=true;
                        if(saveMETPhiPlots) {
                           METPhi_BeforeCut->Fill(METp4[0].Phi());
                           METMagVsMETPhi_BeforeCut->Fill(METp4[0].Phi(),METp4[0].Mag());
                           METxVsNPV->Fill(vp4[0].npv,METp4[0].X());
                           METyVsNPV->Fill(vp4[0].npv,METp4[0].Y());
                           METxVsMETy_BeforeCut->SetPoint(METxVsMETy_BeforeCut->GetN(),METp4[0].X(),METp4[0].Y());
                           if(jp4.size()==0) {
                              METxVsNPV_0J->Fill(vp4[0].npv,METp4[0].X());
                              METyVsNPV_0J->Fill(vp4[0].npv,METp4[0].Y());
                           }
                           if(jp4.size()==1) {
                              METxVsNPV_1J->Fill(vp4[0].npv,METp4[0].X());
                              METyVsNPV_1J->Fill(vp4[0].npv,METp4[0].Y());
                           }
                           if(jp4.size()>0 && METp4.size()>0) {
                              METxVsJetPt->Fill(jp4[0].Pt(),METp4[0].X());
                              METyVsJetPt->Fill(jp4[0].Pt(),METp4[0].Y());
                              METparaVsJPhi->Fill(jp4[0].Phi(),((METp4[0].X()*jp4[0].X())+(METp4[0].Y()*jp4[0].Y()))/TMath::Sqrt(TMath::Power(jp4[0].X(),2)+TMath::Power(jp4[0].Y(),2)));
                              METperpVsJPhi->Fill(jp4[0].Phi(),((METp4[0].X()*jp4[0].Y())-(METp4[0].Y()*jp4[0].X()))/TMath::Sqrt(TMath::Power(jp4[0].X(),2)+TMath::Power(jp4[0].Y(),2)));
                              if(jp4.size()==1) {
                                 METparaVsJPhi_1J->Fill(jp4[0].Phi(),((METp4[0].X()*jp4[0].X())+(METp4[0].Y()*jp4[0].Y()))/TMath::Sqrt(TMath::Power(jp4[0].X(),2)+TMath::Power(jp4[0].Y(),2)));
                                 METperpVsJPhi_1J->Fill(jp4[0].Phi(),((METp4[0].X()*jp4[0].Y())-(METp4[0].Y()*jp4[0].X()))/TMath::Sqrt(TMath::Power(jp4[0].X(),2)+TMath::Power(jp4[0].Y(),2)));
                                 METparaVsJPhi_1J_Pt[getBin(jp4[0].Pt(),vpt,NPtBins)]->Fill(jp4[0].Phi(),((METp4[0].X()*jp4[0].X())+(METp4[0].Y()*jp4[0].Y()))/TMath::Sqrt(TMath::Power(jp4[0].X(),2)+TMath::Power(jp4[0].Y(),2)));
                                 METperpVsJPhi_1J_Pt[getBin(jp4[0].Pt(),vpt,NPtBins)]->Fill(jp4[0].Phi(),((METp4[0].X()*jp4[0].Y())-(METp4[0].Y()*jp4[0].X()))/TMath::Sqrt(TMath::Power(jp4[0].X(),2)+TMath::Power(jp4[0].Y(),2)));
                              }
                           }
                        }
                        if (METp4[0].Pt() > MET_PtMin) {
                           incrementCounter(7,Nj,tableEl,tableLp);
                           if (!muONLY) printEventInformation(printEventInfo, 6, false);
                           el_passAll=true;
                           //METp4[0].SetX(METp4[0].Px() + 0.179189);
                           //METp4[0].SetY(METp4[0].Py() + 3.32732);
                           if(saveMETPhiPlots) {
                              METxVsMETy_AfterCut->SetPoint(METxVsMETy_BeforeCut->GetN(),METp4[0].X(),METp4[0].Y());
                              METPhi_AfterCut->Fill(METp4[0].Phi());
                              METMagVsMETPhi_AfterCut->Fill(METp4[0].Phi(),METp4[0].Mag());
                           }
                           //
                           // Record the Btag info
                           //
                           if (nBtagSSV==0)
                              incrementCounter(14,Nj,tableEl,tableLp);
                           else {
                              if (nBtagSSV==1)
                                 incrementCounter(15,Nj,tableEl,tableLp);
                              else {
                                 if (nBtagSSV==2)
                                    incrementCounter(16,Nj,tableEl,tableLp);
                                 else
                                    incrementCounter(17,Nj,tableEl,tableLp);
                              }
                           }//BTags

                           /////////////Analysis Level Cuts
                           if (lp4[0].leptonCat==DEFS::electron) {
                              incrementCounter(8,Nj,tableEl,tableLp);
                              if (lp4[0].Pt()>30) {
                                 incrementCounter(9,Nj,tableEl,tableLp);
                                 if (METp4[0].Pt()>25) {
                                    incrementCounter(10,Nj,tableEl,tableLp);
                                    if (jp4.size()>0 && jp4[0].Pt()>30) {
                                       incrementCounter(11,Nj,tableEl,tableLp);
                                       if (jp4.size()>1 && jp4[1].Pt()>25) {
                                          incrementCounter(12,Nj,tableEl,tableLp);
                                          if (lp4[0].lpfIso>0.3 && lp4[0].lpfIso<0.7) {
                                             incrementCounter(13,Nj,tableEl,tableLp);
                                          }
                                       }
                                    }
                                 }
                              }
                           }

                           /*
                           /////////////FNAL CUTs
                           if (METp4[0].Et() > 30) {
                              incrementCounter(8,Nj,tableEl,tableLp);
                              TLorentzVector mt(lp4[0].Px()+METp4[0].Px(),lp4[0].Py()+METp4[0].Py(),0,lp4[0].Et()+METp4[0].Et());
                              if (mt.M() > 30) {
                                 incrementCounter(9,Nj,tableEl,tableLp);
                                 if (lp4[0].Pt() > 30) {
                                    incrementCounter(10,Nj,tableEl,tableLp);
                                    if (abs(lp4[0].Eta()) < 2.5) {
                                       incrementCounter(11,Nj,tableEl,tableLp);
                                       if (jp4.size()!=0 && jp4[0].Pt() > 30) {
                                          incrementCounter(12,Nj,tableEl,tableLp);
                                          if (jp4.size() > 1 && jp4[1].Pt() > 30) {
                                             incrementCounter(13,Nj,tableEl,tableLp);
                                             if (abs(jp4[0].Eta()) < 2.4) {
                                                incrementCounter(14,Nj,tableEl,tableLp);
                                                if (abs(jp4[1].Eta()) < 2.4) {
                                                   incrementCounter(15,Nj,tableEl,tableLp);
                                                   bool dr = true;
                                                   for (unsigned int i = 0; i<Nj; i++){
                                                      TLorentzVector j(jp4[i].Px(),jp4[i].Py(),jp4[i].Pz(),jp4[i].E());
                                                      if (j.DeltaR(lp4[0]) <= 0.3)
                                                         dr = false;
                                                      else
                                                         continue;
                                                   }
                                                   if (dr) {
                                                      incrementCounter(16,Nj,tableEl,tableLp);
                                                      TLorentzVector jj(jp4[0].Px()+jp4[1].Px(),jp4[0].Py()+jp4[1].Py(),jp4[0].Pz()+jp4[1].Pz(),jp4[0].E()+jp4[1].E());
                                                      TLorentzVector j4;
                                                      if (jp4.size()>3) {
                                                         j4.SetPx(jp4[3].Px());
                                                         j4.SetPy(jp4[3].Py());
                                                         j4.SetPz(jp4[3].Pz());
                                                         j4.SetE(jp4[3].E());
                                                      }
                                                      else {
                                                         j4.SetPx(0);
                                                         j4.SetPy(0);
                                                         j4.SetPz(0);
                                                         j4.SetE(0);
                                                      }
                                                      if (j4.Pt()<30 && jj.M()>65 && jj.M()<95) {
                                                         incrementCounter(17,Nj,tableEl,tableLp);
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }//First FNAL Cut
                           */
                        }
                     }
                  }
               }//electrons
            } //if (PVfound)
         } //if (el_passTrigger)

         if (mu_passTrigger) {
            incrementCounter(2,Nj,tableMu,tableLp);
            printEventInformation(printEventInfo, 1, true);
            if (PVfound) {
               //printEventInformation(true, -1, true);
               //printEventInformation(true, 7, true);
               incrementCounter(3,Nj,tableMu,tableLp);
               printEventInformation(printEventInfo, 2, true);
               //
               // Record the muon info
               //
               if (mucnt_Prim>=leptonCnt) {
                  incrementCounter(4,Nj,tableMu,tableLp);
                  if (!elONLY) printEventInformation(printEventInfo, 3, true);
                  if (mucnt_Loose==leptonCnt) {
                     incrementCounter(5,Nj,tableMu,tableLp);
                     if (!elONLY)printEventInformation(printEventInfo, 4, true);
                     if (elcnt_Loose==0) {
                        incrementCounter(6,Nj,tableMu,tableLp);
                        if (!elONLY) printEventInformation(printEventInfo, 5, true);
                        mu_passStandard=true;
                        if (saveMETPhiPlots) {
                           METPhi_BeforeCut->Fill(METp4[0].Phi());
                           METMagVsMETPhi_BeforeCut->Fill(METp4[0].Phi(),METp4[0].Mag());
                           METxVsNPV->Fill(vp4[0].npv,METp4[0].X());
                           METyVsNPV->Fill(vp4[0].npv,METp4[0].Y());
                           METxVsMETy_BeforeCut->SetPoint(METxVsMETy_BeforeCut->GetN(),METp4[0].X(),METp4[0].Y());
                           if(jp4.size()==0) {
                              METxVsNPV_0J->Fill(vp4[0].npv,METp4[0].X());
                              METyVsNPV_0J->Fill(vp4[0].npv,METp4[0].Y());
                           }
                           if(jp4.size()==1) {
                              METxVsNPV_1J->Fill(vp4[0].npv,METp4[0].X());
                              METyVsNPV_1J->Fill(vp4[0].npv,METp4[0].Y());
                           }
                           if(jp4.size()>0 && METp4.size()>0) {
                              METxVsJetPt->Fill(jp4[0].Pt(),METp4[0].X());
                              METyVsJetPt->Fill(jp4[0].Pt(),METp4[0].Y());
                              METparaVsJPhi->Fill(jp4[0].Phi(),((METp4[0].X()*jp4[0].X())+(METp4[0].Y()*jp4[0].Y()))/TMath::Sqrt(TMath::Power(jp4[0].X(),2)+TMath::Power(jp4[0].Y(),2)));
                              METperpVsJPhi->Fill(jp4[0].Phi(),((METp4[0].X()*jp4[0].Y())-(METp4[0].Y()*jp4[0].X()))/TMath::Sqrt(TMath::Power(jp4[0].X(),2)+TMath::Power(jp4[0].Y(),2)));
                              if(jp4.size()==1) {
                                 METparaVsJPhi_1J->Fill(jp4[0].Phi(),((METp4[0].X()*jp4[0].X())+(METp4[0].Y()*jp4[0].Y()))/TMath::Sqrt(TMath::Power(jp4[0].X(),2)+TMath::Power(jp4[0].Y(),2)));
                                 METperpVsJPhi_1J->Fill(jp4[0].Phi(),((METp4[0].X()*jp4[0].Y())-(METp4[0].Y()*jp4[0].X()))/TMath::Sqrt(TMath::Power(jp4[0].X(),2)+TMath::Power(jp4[0].Y(),2)));
                                 METparaVsJPhi_1J_Pt[getBin(jp4[0].Pt(),vpt,NPtBins)]->Fill(jp4[0].Phi(),((METp4[0].X()*jp4[0].X())+(METp4[0].Y()*jp4[0].Y()))/TMath::Sqrt(TMath::Power(jp4[0].X(),2)+TMath::Power(jp4[0].Y(),2)));
                                 METperpVsJPhi_1J_Pt[getBin(jp4[0].Pt(),vpt,NPtBins)]->Fill(jp4[0].Phi(),((METp4[0].X()*jp4[0].Y())-(METp4[0].Y()*jp4[0].X()))/TMath::Sqrt(TMath::Power(jp4[0].X(),2)+TMath::Power(jp4[0].Y(),2)));
                              }
                           }
                        }
                        if (METp4[0].Pt() > MET_PtMin) {
                           incrementCounter(7,Nj,tableMu,tableLp);
                           if (!elONLY) printEventInformation(printEventInfo, 6, true);
                           mu_passAll=true;
                           //METp4[0].SetX(METp4[0].Px() + 0.179189);
                           //METp4[0].SetY(METp4[0].Py() + 3.32732);
                           if(saveMETPhiPlots) {
                              METxVsMETy_AfterCut->SetPoint(METxVsMETy_BeforeCut->GetN(),METp4[0].X(),METp4[0].Y());
                              METPhi_AfterCut->Fill(METp4[0].Phi());
                              METMagVsMETPhi_AfterCut->Fill(METp4[0].Phi(),METp4[0].Mag());
                           }
                           //
                           // Record the Btag info
                           //
                           if (nBtagSSV==0 )
                              incrementCounter(14,Nj,tableMu,tableLp);
                           else {
                              if (nBtagSSV==1) 
                                 incrementCounter(15,Nj,tableMu,tableLp);
                              else {
                                 if (nBtagSSV==2)
                                    incrementCounter(16,Nj,tableMu,tableLp);
                                 else
                                    incrementCounter(17,Nj,tableMu,tableLp);
                              }
                           }//BTags
                           
                           /////////////Analysis Level Cuts
                           if (lp4[0].leptonCat==DEFS::muon) {
                              incrementCounter(8,Nj,tableMu,tableLp);
                              if (lp4[0].Pt()>25) {
                                 incrementCounter(9,Nj,tableMu,tableLp);
                                 if (METp4[0].Pt()>25) {
                                    incrementCounter(10,Nj,tableMu,tableLp);
                                    if (jp4.size()>0 && jp4[0].Pt()>30) {
                                       incrementCounter(11,Nj,tableMu,tableLp);
                                       if (jp4.size()>1 && jp4[1].Pt()>25) {
                                          incrementCounter(12,Nj,tableMu,tableLp);
                                          if (lp4[0].lpfIso>0.3 && lp4[0].lpfIso<2.0) {
                                             incrementCounter(13,Nj,tableMu,tableLp);
                                          }
                                       }
                                    }
                                 }
                              }
                           }

                           /*
                           /////////////FNAL CUTs
                           if (METp4[0].Et() > 25) {
                              incrementCounter(8,Nj,tableMu,tableLp);
                              TLorentzVector mt(lp4[0].Px()+METp4[0].Px(),lp4[0].Py()+METp4[0].Py(),0,lp4[0].Et()+METp4[0].Et());
                              if (mt.M() > 30) {
                                 incrementCounter(9,Nj,tableMu,tableLp);
                                 if (lp4[0].Pt() > 25) {
                                    incrementCounter(10,Nj,tableMu,tableLp);
                                    if (abs(lp4[0].Eta()) < 2.1) {
                                       incrementCounter(11,Nj,tableMu,tableLp);
                                       if (jp4.size()!=0 && jp4[0].Pt() > 30) {
                                          incrementCounter(12,Nj,tableMu,tableLp);
                                          if (jp4.size() > 1 && jp4[1].Pt() > 30) {
                                             incrementCounter(13,Nj,tableMu,tableLp);
                                             if (abs(jp4[0].Eta()) < 2.4) {
                                                incrementCounter(14,Nj,tableMu,tableLp);
                                                if (abs(jp4[1].Eta()) < 2.4) {
                                                   incrementCounter(15,Nj,tableMu,tableLp);
                                                   bool dr = true;
                                                   for (unsigned int i = 0; i<Nj; i++){
                                                      TLorentzVector j(jp4[i].Px(),jp4[i].Py(),jp4[i].Pz(),jp4[i].E());
                                                      if (j.DeltaR(lp4[0]) <= 0.3)
                                                         dr = false;
                                                      else
                                                         continue;
                                                   }
                                                   if (dr) {
                                                      incrementCounter(16,Nj,tableMu,tableLp);
                                                      TLorentzVector jj(jp4[0].Px()+jp4[1].Px(),jp4[0].Py()+jp4[1].Py(),jp4[0].Pz()+jp4[1].Pz(),jp4[0].E()+jp4[1].E());
                                                      TLorentzVector j4;
                                                      if (jp4.size()>3) {
                                                         j4.SetPx(jp4[3].Px());
                                                         j4.SetPy(jp4[3].Py());
                                                         j4.SetPz(jp4[3].Pz());
                                                         j4.SetE(jp4[3].E());
                                                      }
                                                      else {
                                                         j4.SetPx(0);
                                                         j4.SetPy(0);
                                                         j4.SetPz(0);
                                                         j4.SetE(0);
                                                      }
                                                      if (j4.Pt()<30 && jj.M()>65 && jj.M()<95) {
                                                         incrementCounter(17,Nj,tableMu,tableLp);
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }//First FNAL Cut
                           */
                        }
                     }
                  }
               } //if(mucnt_Prim>=1 && PVfound)
            } //if (PVfound)
         } //if (mu_passTrigger||el_passTrigger)
      } //if(jcnt_tot==Nj || (Nj==(DEFS::NJETS-1) && jcnt_tot>(DEFS::NJETS-1))) 
   } //for(Int_t Nj=0; Nj<DEFS::NJETS;Nj++)
   printEventInformation(printEventInfo, 7, true);
   if (printJetInfo)
      printJetInformation();
   if (printLeptonInfo)
      printLeptonInformation();

   //
   // Fill the Ntuple to be used in Matrix Element computations
   //    and the other ntuples for future analysis
   //    (requiring 2 jets and 1 tight lepton passing all of the cuts)
   // Usual option: fill in the elements which pass all of the cuts, *except* the MET cut
   //
   setFlags();
   if (lp4.size()>0) {
      for(unsigned int i=0; i<lp4.size(); i++){
         PFIsoDist->Fill(lp4[i].lpfIso);
      }
   }
   if ((el_passFlag || mu_passFlag) || OnVsOffShell) {
      //
      // Event Level Information
      //
      EvtNtuple->run=runNumber;
      EvtNtuple->event=eventNumber;
      EvtNtuple->lumi=lumiNumber;
      EvtNtuple->rho=*rhoHandle;       

      if (!Data && saveGenParticles)
         saveGenPart();

      //
      // Trigger Infromation
      //
      EvtNtuple->triggerMap = triggerMap;

      //
      // Vertex Information
      //
      EvtNtuple->vLV.clear();
      EvtNtuple->vLV = vp4;

      //
      // Jet Information
      //
      EvtNtuple->jLV.clear();
      EvtNtuple->jLV = jp4;

      // Transfer Function Variables
      // Do the matching between the pat::Jets and the reco::genParticles
      pair<int, TLorentzVector> matchedPair;
      if (!Data) {
         vector<int> pdgId = matchToGen_particleCollection(true);
         if(pdgId.size()!=jp4.size()) cout <<"PerformSelection::Analyze::WARNING:: sizes don't match!!!"  << endl
                                           << "recJets.size()=" << jp4.size() << "\t matches.size()=" << pdgId.size() << endl;
         for (unsigned int j = 0; j<jp4.size(); j++) {
            matchedPair = matchToGen(jp4[j].Eta(),jp4[j].Phi());
            EvtNtuple->jLV[j].matchedGenParticles = matchedPair.second;
            if(pdgId.size()==jp4.size()) {
               EvtNtuple->jLV[j].matchedpdgId = pdgId[j];//matchedPair.first;
            }
            EvtNtuple->jLV[j].matchedDeltaR = reco::deltaR(matchedPair.second.Eta(), matchedPair.second.Phi(), jp4[j].Eta(), jp4[j].Phi());
         }
      }
      EvtNtuple->Mjj=Mjj;

      //
      // Lepton Information
      //
      EvtNtuple->lLV.clear();
      EvtNtuple->lLV = lp4;
      if (abs(lEta)<etaBarrelMax) {
         //
         // barrel
         //
         EvtNtuple->lLV[0].ldetComp = 0;
      } 
      else {
         //
         // endcap
         //
         EvtNtuple->lLV[0].ldetComp = 1;
      }
      if (el_passAll)
         EvtNtuple->lLV[0].leptonCat_passAll=DEFS::electron;
      if (mu_passAll)
         EvtNtuple->lLV[0].leptonCat_passAll=DEFS::muon;
      /*
      if (el_passFlag)
         EvtNtuple->leptonCat=DEFS::electron;
      if (mu_passFlag)
         EvtNtuple->leptonCat=DEFS::muon;
      */

      //
      // MET Information
      //
      EvtNtuple->METLV.clear();
      EvtNtuple->METLV = METp4;

      //
      // Fill The Tree
      //
      int treeSwitch = jcnt_tot;
      if (treeSwitch > 2)
         treeSwitch = 2;
      if (OnVsOffShell)
         treeSwitch = 2;
      switch (treeSwitch) {
         case 0:
            if (StoreJets0)
               EvtTree_0Jets->Fill();
            break;
         case 1:
            if (StoreJet1)
               EvtTree_1Jets->Fill();
            break;
         case 2:
            EvtTree_2pJets->Fill();
            break;
         default:
            cout << "WARNING::Value of jcnt_tot does not correspond to a known tree." << endl
                 << "\tNo tree will be filled with event " << eventNumber << "." << endl;
      }

   } //if((el_passFlag || mu_passFlag)) //&& jcnt_tot==2)

   EvtTotCount++;

   return;
}


//______________________________________________________________________________
void PerformSelection::endJob()
{
   for (unsigned int Nj=0; Nj<DEFS::NJETS;Nj++) {
      (*(*tableLp)(DEFS::getCutLevelString((DEFS::CutLevel)0),DEFS::getJetBinString((DEFS::JetBin)Nj))) = Value(EvtTotCount,0.0);
      (*(*tableMu)(DEFS::getCutLevelString((DEFS::CutLevel)0),DEFS::getJetBinString((DEFS::JetBin)Nj))) = Value(EvtTotCount,0.0);
      (*(*tableEl)(DEFS::getCutLevelString((DEFS::CutLevel)0),DEFS::getJetBinString((DEFS::JetBin)Nj))) = Value(EvtTotCount,0.0);
   }
}


////////////////////////////////////////////////////////////////////////////////
// implement trigger functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void PerformSelection::triggerSelection() {
   ///Do not apply trigger selection to MC
   mu_passTrigger = true;
   el_passTrigger = true;
   if (Data || MCpTrigger) {
      triggerMap.clear();
      mu_passTrigger = false;
      el_passTrigger = false;
      pat::TriggerEvent const * trig = &*triggerHandle;
      
      const pat::TriggerPathCollection* triggerPaths = trig->paths();
      for (unsigned int i=0; i<triggerPaths->size(); i++) {
         if ((*triggerPaths)[i].wasRun() && !(*triggerPaths)[i].wasError())
            triggerMap[(*triggerPaths)[i].name()] = (*triggerPaths)[i].wasAccept();
      }
      
      if (trig->wasRun() && trig->wasAccept()) {
         //
         // match muon trigger names to our wild card expression                                                         
         //
         for (unsigned int iTrig=0; iTrig<muTrigger.size(); iTrig++) {
            TRegexp  muTrigRegexp(muTrigger[iTrig]);
            bool matchMuTrigName = false;
            for(pat::TriggerPathCollection::const_iterator iPath = trig->paths()->begin(); iPath != trig->paths()->end(); ++iPath){
               TString thisTrigPath = iPath->name();
               matchMuTrigName =  thisTrigPath.Contains(muTrigRegexp);
               if(matchMuTrigName == true){
                  pat::TriggerPath const * muPath = trig->path(iPath->name());
                  if (muPath != 0 && muPath->wasAccept()) {
                     mu_passTrigger = true;
                  }
               }
            }
         }
               
         //
         // match electron trigger names to our wild card expression
         //
         for (unsigned int iTrig=0; iTrig<eleTrigger.size(); iTrig++) {
            TRegexp  eleTrigRegexp(eleTrigger[iTrig]);
            bool matchElTrigName = false;
            for(pat::TriggerPathCollection::const_iterator iPath = trig->paths()->begin(); iPath != trig->paths()->end(); ++iPath){
               TString thisTrigPath = iPath->name();
               matchElTrigName =  thisTrigPath.Contains(eleTrigRegexp);
               if(matchElTrigName == true){
                  pat::TriggerPath const * elePath = trig->path(iPath->name());
                  if (elePath != 0 && elePath->wasAccept()) {
                     el_passTrigger = true;
                  }
               }
            }
         }
      }
   }
}


////////////////////////////////////////////////////////////////////////////////
// implement event functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void PerformSelection::setFlags()
{
   el_passFlag=el_passAll;
   mu_passFlag=mu_passAll;
   if(noMETCut) {
      el_passFlag=el_passStandard;
      mu_passFlag=mu_passStandard;
   }
   if(elONLY)
      mu_passFlag=false;
   if(muONLY)
      el_passFlag=false;
}


////////////////////////////////////////////////////////////////////////////////
// implement vertex functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void PerformSelection::vertexSelection()
{
   vp4.clear();
   PVfound = (vtxHandle->begin() != vtxHandle->end());
   if (PVfound) {
      pv = (*vtxHandle->begin());
      if (pv.ndof()>=4 && fabs(pv.z())<=24 && fabs(pv.position().Rho())<=2.0) {
         vtxcnt++;
         vp4.push_back(Vertex(pv.x(),pv.y(),pv.z(),0));
         vp4.back().vNDOF = pv.ndof();

         //
         // Pileup information
         // 
         const reco::VertexCollection::const_iterator vtxEnd = vtxHandle->end();
         for (reco::VertexCollection::const_iterator vtxIter = vtxHandle->begin(); vtxEnd != vtxIter; ++vtxIter) {
            if (!vtxIter->isFake() && vtxIter->ndof()>=4 && fabs(vtxIter->z())<=24 && fabs(vtxIter->position().Rho())<=2.0)
               ++vp4.back().npv;
         }
         vp4.back().tnpus.clear();
         vp4.back().npus.clear();
         if (!Data) {
            for(unsigned int i=0; i<pileupHandle->size(); i++) {
               vp4.back().tnpus.push_back((*pileupHandle)[i].getTrueNumInteractions());
               vp4.back().npus.push_back((*pileupHandle)[i].getPU_NumInteractions());
            }
         }
      }
      else {
         PVfound=false; 
         vp4.push_back(Vertex(-999,-999,-999,-999));
      }
   }
   else {
      throw cms::Exception("InvalidInput") << " There needs to be at least one primary vertex in the event."
                                           << std::endl;
   }
}


////////////////////////////////////////////////////////////////////////////////
// implement jet functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void PerformSelection::jetClear() {
   jcnt_tot=0;
   j_DRlepton=999;
   Mjj=-1;
     
   jNEntries=0;
   jp4.clear();
   nBtagSSV=0;
   nBtagTC=0;
   jstreams.clear();
}

//______________________________________________________________________________
void PerformSelection::jetSelection() {
   //
   // Define the constants
   //
   double pi_=TMath::Pi();

   const vector< pat::Jet >::const_iterator pfJetEnd = pfJetHandle->end();
   for (vector<pat::Jet>::const_iterator jetIter = pfJetHandle->begin(); pfJetEnd != jetIter; ++jetIter) { 

      math::XYZTLorentzVector vvv = jetIter->p4();

      jesUncScale = 1.0;
      if (doJESUncertainty) {
         if (JESUncertainty == "none") {
            jesUncScale *= 1.0;
         } 
         else {
            jecUnc->setJetEta(jetIter->eta());
            jecUnc->setJetPt(jetIter->pt());
            
            uncert = jecUnc->getUncertainty(true);
            
            if (JESUncertainty == "up") {
               jesUncScale *= (1 + uncert);
            } 
            else if (JESUncertainty == "down") {
               jesUncScale *= (1 - uncert);
            }
         }
      }

      JERCor = 1.0;
      if (doJER && !Data) {
         // get the JER correction factor for this jet
         (jetIter->genJet() && vvv.pt()>10)? JERCor = getJERfactor(jetIter->pt(), jetIter->eta(), jetIter->genJet()->pt()) : JERCor = 1.0;

         // recompute the met for this change
         double newMetX = (1 - JERCor)*jetIter->p4().X() + METp4[0].X();
         double newMetY = (1 - JERCor)*jetIter->p4().Y() + METp4[0].Y();

         METp4[0].SetX(newMetX);
         METp4[0].SetY(newMetY);
      }

      vvv *= jesUncScale;      
      vvv *= JERCor;

      j_pt=vvv.pt();
      j_eta=vvv.eta();
      j_phi=vvv.phi();

      //
      // loop through all of the leptons to check the DR with the current jet
      //
      bool jDR = true;
      for (unsigned int i=0; i<lp4.size(); i++) {
         adphi=abs(lp4[i].Phi()-j_phi);
         //
         // we're on a circle
         //
         if (adphi>pi_)
            adphi=2*pi_-adphi;
         
         j_DRlepton=pow((j_eta-lp4[i].Eta())*(j_eta-lp4[i].Eta())+adphi*adphi,0.5);
         
         if ((lp4[i].leptonCat==DEFS::electron && j_DRlepton<j_DRelMin) || (lp4[i].leptonCat==DEFS::muon && j_DRlepton<muPrim_DRjMin))
            jDR = false;
         else
            continue;
      }

      double genJetPt = -999;
      if(jetIter->genJet()) genJetPt = jetIter->genJet()->pt();
      jstream << "Jet number = " << jp4.size()-1 << endl
              << "fancyJES_ = " << JESUncertainty << ", uncert = " << uncert << endl
              << "Original PT = " << jetIter->pt() << ", Eta = " << jetIter->eta() << ", ptScale = " << jesUncScale*JERCor << ", GenJet Pt = " << genJetPt << endl
              << "Scaled Jet PT = " << vvv.Pt() << ", Eta = " << vvv.Eta() << endl;
      jstreams.push_back(jstream.str());
      jstream.str("");

      //
      // Apply The Cut
      //
      if (j_pt>j_ptMin && abs(j_eta)<j_aetaMax && jDR &&
          jetIter->neutralHadronEnergyFraction()<NHEFMax && jetIter->neutralEmEnergyFraction()<NEMEFMax && 
          jetIter->numberOfDaughters()>NDaughtersMin && jetIter->chargedHadronEnergyFraction()>CHEFMin && 
          jetIter->chargedMultiplicity()>CMultiplicityMin && jetIter->chargedEmEnergyFraction()<CEMEFMax) {
         jcnt_tot++;
         jNEntries++;
         //
         // store the parameters
         //

         //
         // Corrected (reco) jet
         //
         //math::XYZTLorentzVector vvv = jetIter->p4();
         jp4.push_back(Jet(vvv.Px(),vvv.Py(),vvv.Pz(),vvv.E()));
         jp4.back().jChargedMultiplicity = jetIter->chargedMultiplicity();
         jp4.back().jNeutralMultiplicity = jetIter->neutralMultiplicity();
         jp4.back().jNeutralHadronEnergyFraction = jetIter->neutralHadronEnergyFraction();
         jp4.back().jNeutralEmEnergyFraction = jetIter->neutralEmEnergyFraction();
         jp4.back().jChargedHadronEnergyFraction = jetIter->chargedHadronEnergyFraction();
         jp4.back().jChargedEmEnergyFraction = jetIter->chargedEmEnergyFraction();
         jp4.back().jNumberOfDaughters = jetIter->numberOfDaughters();
         jp4.back().partonFlavour = jetIter->partonFlavour();
         //PFConstituents not in SQWaT PATtuples created by V00-00-04 or earlier
         if (SQWaT_Version >= 5)
            jp4.back().jPtD = getPtD(jetIter->getPFConstituents());
         else
            jp4.back().jPtD = 0;

         //
         // Uncorrected (raw) jet
         //
         math::XYZTLorentzVector vv = jetIter->correctedJet("Uncorrected").p4();
         jp4.back().rawLV = TLorentzVector(vv.Px(),vv.Py(),vv.Pz(),vv.E());
         
         //
         // Gen (ref) jet
         //
         if (!Data) {
            if (jetIter->genJet()) {
               math::XYZTLorentzVector v = jetIter->genJet()->p4();
               jp4.back().refLV = TLorentzVector(v.Px(),v.Py(),v.Pz(),v.E());
               jp4.back().refpdgId = jetIter->genJet()->pdgId();
            }
            else {
               TLorentzVector lv(0,0,0,0);
               jp4.back().refLV = lv;
               jp4.back().refpdgId = 0;
            }
         }

         /*
         jstream << "pt " << j_pt << " eta " << j_eta << " j_DRlepton " << j_DRlepton
                 << " NHEF " << jetIter->neutralHadronEnergyFraction()
                 << " NEMEF " << jetIter->neutralEmEnergyFraction()
                 << " CHEF " << jetIter->chargedHadronEnergyFraction()
                 << " CEMEF " << jetIter->chargedEmEnergyFraction()
                 << " num daughters " << jetIter->numberOfDaughters()
                 << " charged multiplicity " << jetIter->chargedMultiplicity()
                 << endl;
         jstreams.push_back(jstream.str());
         jstream.str("");
         */

         //
         // b-tagging
         //
         jp4.back().jBtagDiscriminatorSSV = jetIter->bDiscriminator("simpleSecondaryVertexHighEffBJetTags");
         jp4.back().jBtagDiscriminatorTC = jetIter->bDiscriminator("trackCountingHighEffBJetTags");
         jp4.back().jBtagDiscriminatorCSV = jetIter->bDiscriminator("combinedSecondaryVertexBJetTags");
         jp4.back().jBtagDiscriminatorCSVMVA = jetIter->bDiscriminator("combinedSecondaryVertexMVABJetTags");
         if (jetIter->bDiscriminator("simpleSecondaryVertexHighEffBJetTags")>bDiscriminatorSSVMin) {
            jp4.back().jBtagSSV = 1;
            nBtagSSV++;
         } 
         else
            jp4.back().jBtagSSV = 0;
         if (jetIter->bDiscriminator("trackCountingHighEffBJetTags")>bDiscriminatorTCMin) {
            jp4.back().jBtagTC = 1;
            nBtagTC++;
         } 
         else
            jp4.back().jBtagTC = 0;
         if (jetIter->bDiscriminator("combinedSecondaryVertexBJetTags")>bDiscriminatorCSVMin) {
            jp4.back().jBtagCSV = 1;
            nBtagCSV++;
         } 
         else
            jp4.back().jBtagCSV = 0;
         if (jetIter->bDiscriminator("combinedSecondaryVertexMVABJetTags")>bDiscriminatorCSVMVAMin) {
            jp4.back().jBtagCSVMVA = 1;
            nBtagCSVMVA++;
         } 
         else
            jp4.back().jBtagCSVMVA = 0;
      } //for(Selection Cuts)
   } //for(jetIter)
}

//______________________________________________________________________________
void PerformSelection::computeInvariantMass() {
   if (jcnt_tot==2)
      Mjj=(jp4[0]+jp4[1]).M();
   else
      Mjj=-2;
}

//______________________________________________________________________________
double PerformSelection::getPtD(vector<reco::PFCandidatePtr> pfCandidates) {
   double sumPt_cands=0.;
   double sumPt2_cands=0.;
   
   for (std::vector<reco::PFCandidatePtr>::const_iterator jt = pfCandidates.begin();
        jt != pfCandidates.end(); ++jt) {
      
      // Convert particle momentum to normal TLorentzVector, wrong type :(
      math::XYZTLorentzVectorD const& p4t = (*jt)->p4();
      TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
      if(p4.Pt()!=0){
	     sumPt_cands += p4.Pt();
	     sumPt2_cands += (p4.Pt()*p4.Pt());
      }
   }

   if (sumPt_cands != 0)
      return sqrt( sumPt2_cands )/sumPt_cands;
   else
      return 0.0;
}


////////////////////////////////////////////////////////////////////////////////
// implement b-tag functions
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// implement muon functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void PerformSelection::muonClear() {
   mucnt_Prim=0;
   mucnt_Loose=0;
}

//______________________________________________________________________________
void PerformSelection::muonSelection() {
   const vector< pat::Muon >::const_iterator muEnd = muonHandle->end();
   for(vector<pat::Muon>::const_iterator muIter = muonHandle->begin(); muEnd != muIter; ++muIter) { 
      mu_pt=muIter->pt();
      mu_eta=muIter->eta();
      mu_phi=muIter->phi();

      //
      // calculate the relative/tracker isolation of the muon
      //
      //Double_t rhoPrime = std::max((*rhoHandle),0.0);
      //Double_t AEff = MuonEffectiveArea::GetMuonEffectiveArea(MuonEffectiveArea::kMuGammaAndNeutralHadronIso04, mu_eta, MuonEffectiveArea::kMuEAData2011);
      mu_PFIso = (muIter->pfIsolationR04().sumChargedHadronPt + max(0.0, muIter->pfIsolationR04().sumNeutralHadronEt + muIter->pfIsolationR04().sumPhotonEt - 0.5*(muIter->pfIsolationR04().sumPUPt))) / mu_pt;
      mu_DetIso = (muIter->ecalIso() + muIter->hcalIso() + muIter->trackIso() -
                   (*rhoHandle * TMath::Pi() * pow(muDetIsoConeSize,2)));
      mu_DetIso /= mu_pt;
      mu_TrkIso=(muIter->trackIso())/mu_pt;
      
      muisGmTm=false;
      if (muIter->isGlobalMuon()/* && muIter->isTrackerMuon()*/)
         muisGmTm=true;

      muisPFm=false;
      if (muIter->isPFMuon())
         muisPFm=true;

      muisTightPrompt=false;
      if (muIter->globalTrack().isNonnull()) {
         if (muIter->globalTrack()->normalizedChi2()<muPrim_normChi2Max &&
            muIter->globalTrack()->hitPattern().numberOfValidMuonHits()>muPrim_nValidMuonHitsMin)
            muisTightPrompt=true;
      }

      //if (muIter->innerTrack().isNonnull())
         //muPrim_itNHits=muIter->innerTrack()->numberOfValidHits();
      //else
      //muPrim_itNHits=-1;
      if (muIter->innerTrack().isNonnull())
         muPrim_itNHits=muIter->globalTrack()->hitPattern().trackerLayersWithMeasurement();
      else
         muPrim_itNHits=-1;
      
       
      mu_vtxPass=false;
      if (abs(muIter->innerTrack()->dz(vtxHandle->at(0).position()))<muPrim_dzMax && muIter->dB()<muPrim_dBMax)
         mu_vtxPass=true;
       
      //
      // Apply The Primary Cuts
      //
      if (mu_pt>muPrim_ptMin && abs(mu_eta)<muPrim_aetaMax &&
          muisGmTm && muisPFm && muisTightPrompt && mu_vtxPass &&
          muPrim_itNHits>muPrim_nValidHitsMin &&
          muIter->innerTrack()->hitPattern().numberOfValidPixelHits()>muPrim_nValidPixelHitsMin &&
          muIter->numberOfMatchedStations()>muPrim_nMatchedStationsMin) {
         isProperlyIsolatedMu=false;
         
         if ( mu_PFIso<muPrim_PFIsoMax || noMVAIsoCut) //((doDetectorIso && mu_DetIso<muPrim_DetIsoMax)||(doTrackerIso && mu_TrkIso<muPrim_TrkIsoMax)) )
            isProperlyIsolatedMu=true;
            
         //
         // Record the muon quantities:
         //
         if (isProperlyIsolatedMu) {
            mucnt_Prim++;
            math::XYZTLorentzVector lv = muIter->p4();
            lp4.push_back(Lepton(lv.Px(),lv.Py(),lv.Pz(),lv.E())); 
            lp4.back().lQ = muIter->charge();
            lp4.back().leptonCat = DEFS::muon;
            lp4.back().ltrkIso = mu_TrkIso;
            lp4.back().ldetIso = mu_DetIso;
            lp4.back().lpfIso = mu_PFIso;
            if (doDetectorIso)
               lp4.back().lTotIso = mu_DetIso;
            else if (doTrackerIso)
               lp4.back().lTotIso = mu_TrkIso;
            else
               lp4.back().lTotIso = mu_PFIso;
            lp4.back().mIsGlobal = muisGmTm;
            lp4.back().mIsPF = muisPFm;
            lp4.back().mIsTracker = muIter->isTrackerMuon();
            lp4.back().mNormalizedChi2 = muIter->globalTrack()->normalizedChi2();
            lp4.back().mNumberOfValidMuonHits = muIter->globalTrack()->hitPattern().numberOfValidMuonHits();
            lp4.back().mNumberValidHits = muPrim_itNHits;
            lp4.back().mNumberOfValidPixelHits = muIter->innerTrack()->hitPattern().numberOfValidPixelHits();
            lp4.back().mNumberOfMatchedStations = muIter->numberOfMatchedStations();
            lp4.back().ldz = abs(muIter->innerTrack()->dz(vtxHandle->at(0).position()));
            lp4.back().ldB = muIter->dB();

            lstream << "mu_tight pt " << mu_pt << " eta " << abs(mu_eta) << " mu_vtxPass " << mu_vtxPass << " dz " << abs(muIter->innerTrack()->dz(vtxHandle->at(0).position())) << " d0 " << abs(muIter->dB()) << " PFIso " << mu_PFIso << " global " << muIter->isGlobalMuon() << " pf " << muIter->isPFMuon() << " nMuonHits " << muIter->globalTrack()->hitPattern().numberOfValidMuonHits() << " nHits " << muIter->innerTrack()->numberOfValidHits() << " nPixelHits " << muIter->innerTrack()->hitPattern().numberOfValidPixelHits() << " Chi2 " << muIter->globalTrack()->normalizedChi2() << " nMatchedStations " << muIter->numberOfMatchedStations() << endl;
            lstreams.push_back(lstream.str());
            lstream.str("");
         }
      }

      //
      // Apply The Loose Cuts
      //
      if (mu_pt>muLoose_ptMin && abs(mu_eta)<muLoose_aetaMax &&
          (muIter->isGlobalMuon() || muIter->isTrackerMuon()) && muisPFm &&
          //((doDetectorIso && mu_DetIso<muLoose_DetIsoMax)||(doTrackerIso && mu_TrkIso<muLoose_TrkIsoMax))
          (mu_PFIso<muLoose_PFIsoMax || noMVAIsoCut) ) {
         mucnt_Loose++;

         lstream << "mu_loose pt " << mu_pt << " eta " << abs(mu_eta) << " mu_vtxPass " << mu_vtxPass << " dz " << abs(muIter->innerTrack()->dz(vtxHandle->at(0).position())) << " d0 " << abs(muIter->dB()) << " PFIso " << mu_PFIso << " global " << muIter->isGlobalMuon() << " pf " << muIter->isPFMuon() << " nMuonHits " << muIter->globalTrack()->hitPattern().numberOfValidMuonHits() << " nHits " << muIter->innerTrack()->numberOfValidHits() << " nPixelHits " << muIter->innerTrack()->hitPattern().numberOfValidPixelHits() << " Chi2 " << muIter->globalTrack()->normalizedChi2() << " nMatchedStations " << muIter->numberOfMatchedStations() << endl;
         lstreams.push_back(lstream.str());
         lstream.str("");
      }
   }//for muons
}


////////////////////////////////////////////////////////////////////////////////
// implement electron functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void PerformSelection::eleClear() {
   elcnt_Prim=0;
   elcnt_Loose=0;
   lstreams.clear();
   lp4.clear();
}

//______________________________________________________________________________
void PerformSelection::eleSelection() {
   const vector< pat::Electron >::const_iterator elEnd = electronHandle->end();
   for (vector< pat::Electron >::const_iterator elIter = electronHandle->begin(); elEnd != elIter; ++elIter) { 
      el_pt=elIter->p4().Pt();
      el_eta=elIter->eta();
      el_phi=elIter->phi();
      el_aetasc=abs(elIter->superCluster()->eta());

      el_aetaPass=false;
      if (abs(el_eta)<elPrim_aetaMax)
         el_aetaPass=true; //noscv
      if (el_aetasc<elPrim_aetascExcludeMax && el_aetasc>elPrim_aetascExcludeMin)
         el_aetaPass=false;

      //
      // Tight EID Cuts
      //
      el_CutEIDPass = CutBasedEID(elIter,true);
      el_EinvMinusPinvPass = eleEP(elIter);

      //
      // Tight Isolation Cuts
      //
      el_IsoPass = eleIso(elIter,true);

      //
      // Tight Vertex Cuts
      //
      el_vtxPass = eleVtx(elIter,true);

      //
      // Tight conversion Cuts
      //
      el_convPass = elIter->passConversionVeto();

      //
      // Tight MVA based EID
      //
      el_MVAEIDPass = MVABasedEID(elIter,true);
         
      //
      // Combine the Tight Cuts
      //
      if (doMVAeleSel)
         el_PassTight = (el_pt>elPrim_ptMin) && el_aetaPass && el_MVAEIDPass && el_vtxPass && el_convPass;
      else
         el_PassTight = (el_pt>elPrim_ptMin) && el_aetaPass && el_CutEIDPass && el_vtxPass && el_EinvMinusPinvPass && el_convPass && el_IsoPass;

      //
      // Apply The Primary Cuts
      //
      if (el_PassTight) {
         //
         // Record the electron quantities:
         //
         elcnt_Prim++;
         math::XYZTLorentzVector lv = elIter->p4();
         lp4.push_back(Lepton(lv.Px(),lv.Py(),lv.Pz(),lv.E()));
         lp4.back().lQ = elIter->charge();
         lp4.back().leptonCat = DEFS::electron;
         lp4.back().ltrkIso = el_TrkIso;
         lp4.back().ldetIso = el_DetIso;
         lp4.back().lpfIso = el_PFIso;
         if (doDetectorIso)
            lp4.back().lTotIso = el_DetIso;
         else if(doTrackerIso)
            lp4.back().lTotIso = el_TrkIso;
         else
            lp4.back().lTotIso = el_PFIso;
         lp4.back().lphotonIso = elIter->photonIso();
         lp4.back().lchargedHadronIso = elIter->chargedHadronIso();
         lp4.back().lneutralHadronIso = elIter->neutralHadronIso();
         lp4.back().lAEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, elIter->eta(), ElectronEffectiveArea::kEleEAData2012);
         lp4.back().ldz = abs(elIter->gsfTrack()->dz(pv.position()));
         lp4.back().ldB = abs(elIter->dB());
         lp4.back().emvaTrig = el_mvaTrig;
         lp4.back().emvaNonTrig = el_mvaNonTrig;
         lp4.back().ePassConversionVeto = el_convPass;
         lp4.back().eIsEB = elIter->isEB();
         lp4.back().eIsEE = elIter->isEE();
         lp4.back().eSuperClusterEta = elIter->superCluster()->eta();
         lp4.back().eEcalEnergy = elIter->ecalEnergy();
         lp4.back().eESuperClusterOverP = elIter->eSuperClusterOverP();
         lp4.back().eSigmaIetaIeta = elIter->sigmaIetaIeta();
         lp4.back().eDeltaPhi = elIter->deltaPhiSuperClusterTrackAtVtx();
         lp4.back().eDeltaEta = elIter->deltaEtaSuperClusterTrackAtVtx();
         lp4.back().eHadronicOverEm = elIter->hadronicOverEm();

         lstream << "el_tight pt " << el_pt << " eta " << abs(el_eta) << " dz " << abs(elIter->gsfTrack()->dz(pv.position())) << " d0 " << abs(elIter->dB()) << " PFIso " << el_PFIso << " EB " << elIter->isEB() << " Passed conversion veto " << el_convPass << " sigmaIetaIeta " << elIter->sigmaIetaIeta() << " DeltaPhi " << abs(elIter->deltaPhiSuperClusterTrackAtVtx()) << " DeltaEta " << abs(elIter->deltaEtaSuperClusterTrackAtVtx()) << " HoE " << elIter->hadronicOverEm() << " EinvMinusPinv " << el_EinvMinusPinvPass << " nHits " << elIter->gsfTrack()->trackerExpectedHitsInner().numberOfHits() << endl; 
         lstreams.push_back(lstream.str());
         lstream.str("");
      }

      //
      // Loose EID Cuts
      //
      el_CutEIDPass = CutBasedEID(elIter,false);

      //
      // Loose Isolation Cuts
      //
      el_IsoPass = eleIso(elIter,false);

      //
      // Loose Vertex Cuts
      //
      el_vtxPass = eleVtx(elIter,false);

      //
      // Loose MVA based EID
      //
      el_MVAEIDPass = MVABasedEID(elIter,false);

      //
      // Combine the Loose Cuts
      //
      if (doMVAeleSel)
         el_PassLoose = (el_pt>elLoose_ptMin) && el_aetaPass && el_vtxPass && el_MVAEIDPass && el_convPass;
      else
         el_PassLoose = (el_pt>elLoose_ptMin) && el_aetaPass && el_vtxPass && el_CutEIDPass && el_IsoPass;

      //
      // Apply The Loose Cuts
      //
      if (el_PassLoose) {

         elcnt_Loose++;

         lstream << "el_loose pt " << el_pt << " eta " << abs(el_eta) << " dz " << abs(elIter->gsfTrack()->dz(pv.position())) << " d0 " << abs(elIter->dB()) << " PFIso " << el_PFIso << " EB " << elIter->isEB() << " Passed conversion veto " << el_convPass << " sigmaIetaIeta " << elIter->sigmaIetaIeta() << " DeltaPhi " << abs(elIter->deltaPhiSuperClusterTrackAtVtx()) << " DeltaEta " << abs(elIter->deltaEtaSuperClusterTrackAtVtx()) << " HoE " << elIter->hadronicOverEm() << " EinvMinusPinv " << el_EinvMinusPinvPass << " nHits " << elIter->gsfTrack()->trackerExpectedHitsInner().numberOfHits() << endl; 
         lstreams.push_back(lstream.str());
         lstream.str("");
      }
   }
}

//______________________________________________________________________________
bool PerformSelection::CutBasedEID(vector< pat::Electron >::const_iterator elIter, bool tight){
   el_sigmaIetaIeta=elIter->sigmaIetaIeta();
   el_aDeltaPhi=abs(elIter->deltaPhiSuperClusterTrackAtVtx());
   el_aDeltaEta=abs(elIter->deltaEtaSuperClusterTrackAtVtx());
   el_HoE=elIter->hadronicOverEm();
   if (tight) {
      if (invertEID) {
         if (elIter->isEB()) {
            if ( !(el_sigmaIetaIeta<elPrim_sigmaIetaIetaMaxEB) || 
                 !(el_aDeltaPhi<elPrim_aDeltaPhiMaxEB) || 
                 !(el_aDeltaEta<elPrim_aDeltaEtaMaxEB) || 
                 !(el_HoE<elPrim_HoEMaxEB) )
               return true;
         }
         if (elIter->isEE()) {
            if ( !(el_sigmaIetaIeta<elPrim_sigmaIetaIetaMaxEE) || 
                 !(el_aDeltaPhi<elPrim_aDeltaPhiMaxEE) || 
                 !(el_aDeltaEta<elPrim_aDeltaEtaMaxEE) || 
                 !(el_HoE<elPrim_HoEMaxEE) )
               return true;
         }
      }
      else {
         if (elIter->isEB()) {
            if ( (el_sigmaIetaIeta<elPrim_sigmaIetaIetaMaxEB) && 
                 (el_aDeltaPhi<elPrim_aDeltaPhiMaxEB) && 
                 (el_aDeltaEta<elPrim_aDeltaEtaMaxEB) && 
                 (el_HoE<elPrim_HoEMaxEB) )
               return true;
         }
         if (elIter->isEE()) {
            if ( (el_sigmaIetaIeta<elPrim_sigmaIetaIetaMaxEE) && 
                 (el_aDeltaPhi<elPrim_aDeltaPhiMaxEE) && 
                 (el_aDeltaEta<elPrim_aDeltaEtaMaxEE) && 
                 (el_HoE<elPrim_HoEMaxEE) )
               return true;
         }
      }
   }
   else {
      if (elIter->isEB()) {
         if ( (el_sigmaIetaIeta<elLoose_sigmaIetaIetaMaxEB) && 
              (el_aDeltaPhi<elLoose_aDeltaPhiMaxEB) && 
              (el_aDeltaEta<elLoose_aDeltaEtaMaxEB) && 
              (el_HoE<elLoose_HoEMaxEB) )
            return true;
      }
      if (elIter->isEE()) {
         if ( (el_sigmaIetaIeta<elLoose_sigmaIetaIetaMaxEE) && 
              (el_aDeltaPhi<elLoose_aDeltaPhiMaxEE) && 
              (el_aDeltaEta<elLoose_aDeltaEtaMaxEE) )
            return true;
      }
   }
   return false;
}

//______________________________________________________________________________
bool PerformSelection::MVABasedEID(vector< pat::Electron >::const_iterator elIter, bool tight){
   el_mvaTrig = elIter->electronID("mvaTrigV0");
   el_mvaNonTrig = elIter->electronID("mvaNonTrigV0");
   if(noMVAIsoCut) {
      return true;
   }
   if (tight) {
      if (invertEID) {
         if (el_PFIso>elPrim_PFIsoMin_invertEID && elIter->electronID("mvaTrigV0") > -1.0)
            return true;
      }
      else {
         if ( (el_aetasc < 0.8 && elIter->electronID("mvaTrigV0") > 0.977 && el_PFIso < 0.093)
              || (el_aetasc > 0.8 && el_aetasc < 1.479 && elIter->electronID("mvaTrigV0") > 0.956 && el_PFIso < 0.095)
              || (el_aetasc > 1.479 && el_aetasc < 2.5 && elIter->electronID("mvaTrigV0") > 0.966 && el_PFIso < 0.171) )
            return true;
      }
   }
   else {
      if (invertEID) {
         if (el_PFIso>elPrim_PFIsoMin_invertEID && elIter->electronID("mvaNonTrigV0") > -1.0)
            return true;
      }
      else {
         if ( (el_aetasc < 0.8 && elIter->electronID("mvaNonTrigV0") > 0.877 && el_PFIso < 0.426)
             || (el_aetasc > 0.8 && el_aetasc < 1.479 && elIter->electronID("mvaNonTrigV0") > 0.811 && el_PFIso < 0.481)
              || (el_aetasc > 1.479 && el_aetasc < 2.5 && elIter->electronID("mvaNonTrigV0") > 0.707 && el_PFIso < 0.390) )
            return true;
      }
   }
   return false;
}

//______________________________________________________________________________
bool PerformSelection::eleIso(vector< pat::Electron >::const_iterator elIter, bool tight) {
   //
   // calculate the pf/detector/tracker isolation of the electron
   //
   //PF Isolation
   Double_t rhoPrime = std::max((*rhoHandle),0.0);
   Double_t AEff = ElectronEffectiveArea::GetElectronEffectiveArea(ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03, elIter->eta(), ElectronEffectiveArea::kEleEAData2012);
   //if (invertEID)
   //   el_PFIso = (elIter->pfIsolationVariables().chargedHadronIso+ max(0.,elIter->pfIsolationVariables().neutralHadronIso+elIter->pfIsolationVariables().photonIso - AEff*rhoPrime))/elIter->p4().Pt();
   //else
   el_PFIso = (elIter->chargedHadronIso() + max( 0.0, elIter->neutralHadronIso() + elIter->photonIso() - AEff*rhoPrime )) / elIter->pt();
   //Detector Isolation
   if (elIter->isEB()) {
      el_DetIso = ( elIter->dr03TkSumPt() + std::max(0.,elIter->dr03EcalRecHitSumEt() -1.) 
                    + elIter->dr03HcalTowerSumEt() - (*rhoHandle)*TMath::Pi()*pow(elDetIsoConeSize,2)) / elIter->p4().Pt();
   }
   else {
      el_DetIso = ( elIter->dr03TkSumPt() + elIter->dr03EcalRecHitSumEt()
                    + elIter->dr03HcalTowerSumEt() - (*rhoHandle)*TMath::Pi()*pow(elDetIsoConeSize,2)) / elIter->p4().Pt();
   }
   //Tracker Isolation
   el_TrkIso=(elIter->dr03TkSumPt())/elIter->p4().Pt();


   if (tight) {
      if (invertEID) {
         if (el_PFIso>elPrim_PFIsoMin_invertEID)
            return true;
      }
      else {
         if ( (doDetectorIso && el_DetIso<elPrim_DetIsoMax)||(doTrackerIso && el_TrkIso<elPrim_TrkIsoMax)||(doPFIso && el_PFIso<elPrim_PFIsoMax) )
            return true;
      }
   }
   else {
      if ( (doDetectorIso && el_DetIso<elLoose_DetIsoMax)||(doTrackerIso && el_TrkIso<elLoose_TrkIsoMax)||(doPFIso && el_PFIso<elLoose_PFIsoMax) )
         return true;
      }
   return false;
}

//______________________________________________________________________________
bool PerformSelection::eleVtx(vector< pat::Electron >::const_iterator elIter, bool tight) {
   if (tight) {
      if (abs(elIter->gsfTrack()->dz(pv.position()))<elPrim_dzMax && abs(elIter->dB())<elPrim_dBMax)
         return true;
   }
   else {
      if (abs(elIter->gsfTrack()->dz(pv.position()))<elLoose_dzMax && abs(elIter->dB())<elLoose_dBMax)
         return true;
   }
   return false;
}

//______________________________________________________________________________
bool PerformSelection::eleEP(vector< pat::Electron >::const_iterator elIter) {
   if (abs((1/elIter->ecalEnergy())-(1/(elIter->ecalEnergy()/elIter->eSuperClusterOverP())))<elPrim_EinvMinusPinvMax)
      return true;
   else
      return false;
}


////////////////////////////////////////////////////////////////////////////////
// implement lepton functions
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// implement MET functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void PerformSelection::metSelection() {
   //MET_Pass = false;
   METp4.clear();

   math::XYZTLorentzVector v = METHandle->front().p4();
   METp4.push_back(TLorentzVector(v.Px(),v.Py(),v.Pz(),v.E()));

   //MET_Pass = METp4[0].Et() > MET_EtMin;
}

//______________________________________________________________________________
pair<double,double> PerformSelection::getMETPhiCorrection_NPV(TString eraType) {
   pair<double,double> METCor = std::make_pair(0,0);

   //From SingleElectron RunD v1 p8
   /*
   METxVsNPV_0J
   ****************************************
   Minimizer is Linear
   Chi2                      =      37.4033
   NDf                       =           40
   p0                        =     0.201051   +/-   0.10382
   p1                        =     0.426627   +/-   0.0062184

   METyVsNPV_0J
   ****************************************
   Minimizer is Linear
   Chi2                      =      38.4546
   NDf                       =           40
   p0                        =    -0.913497   +/-   0.104427
   p1                        =    -0.231199   +/-   0.0062539
   */
   if (eraType.CompareTo("pfMEtSysShiftCorrParameters_2012runAvsNvtx_TAMUWW_data")==0) {
      //METCor.first  = +1.62143 + (0.489774*vp4[0].npv);
      //METCor.second = -1.41580 - (0.314019*vp4[0].npv);
      METCor.first  = +2.01051137751020431e-01 + (4.26626581124211535e-01*vp4[0].npv);
      METCor.second = -9.13496683888803562e-01 - (2.31198930335691499e-01*vp4[0].npv);
   }

   //From WJets_part1 sample
   /*
   METxVsNPV_0J
   ****************************************
   Minimizer is Linear
   Chi2                      =      43.6772
   NDf                       =           48
   p0                        =     0.290594   +/-   0.051253
   p1                        =  -0.00352926   +/-   0.00320029
   
   METyVsNPV_0J
   ****************************************
   Minimizer is Linear
   Chi2                      =      87.9294
   NDf                       =           48
   p0                        =     0.301828   +/-   0.0515022
   p1                        =    -0.199738   +/-   0.00321429
   */
   else if (eraType.CompareTo("pfMEtSysShiftCorrParameters_2012runAvsNvtx_TAMUWW_mc")==0) {
      METCor.first  = +2.90593612787009126e-01 - (3.52925900447996705e-03*vp4[0].npv);
      METCor.second = +3.01827976672935872e-01 - (1.99738437329797819e-01*vp4[0].npv);
   }
   else {
      METCor.first  = 0;
      METCor.second = 0;
   }

   METCor.first*=-1;
   METCor.second*=-1;

   return  METCor;
}

//______________________________________________________________________________
pair<double,double> PerformSelection::getMETPhiCorrection_JetPhi(TString eraType) {
   pair<double,double> METCor = std::make_pair(0,0);

   if (eraType.CompareTo("pfMEtSysShiftCorrParameters_2012runAvsJPhi_TAMUWW_data")==0) {
      METCor.first  = +0.0 + (0.0*vp4[0].npv);
      METCor.second = +0.0 - (0.0*vp4[0].npv);
   }

   //From WJets_part1 sample
   /*
     TProfile* p = new TProfile("p","p",80,-4,4)
     jet1->Draw("((METLV[0].fP.fX*jLV[0].fP.fX)+(METLV[0].fP.fY*jLV[0].fP.fY))/TMath::Sqrt(TMath::Power(jLV[0].fP.fX,2)+TMath::Power(jLV[0].fP.fY,2)):jLV[0].Phi()>>p","","")

     Parallel Component
     root [11]  FCN=51.5277 FROM MIGRAD    STATUS=CONVERGED      61 CALLS          62 TOTAL
     EDM=5.05296e-07    STRATEGY= 1      ERROR MATRIX ACCURATE
     EXT PARAMETER                                   STEP         FIRST
     NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
     1  p0          -3.75112e+00   7.75361e-02   2.64725e-04  -7.32284e-04
     2  p1           9.88865e-01   1.08493e-02   3.70402e-05  -8.51341e-02
     3  p2          -1.74725e+01   5.32091e-02   1.88268e-04   8.01995e-03

     TProfile* p2 = new TProfile("p2","p2",80,-4,4)
     jet1->Draw("((METLV[0].fP.fX*jLV[0].fP.fY)-(METLV[0].fP.fY*jLV[0].fP.fX))/TMath::Sqrt(TMath::Power(jLV[0].fP.fX,2)+TMath::Power(jLV[0].fP.fY,2)):jLV[0].Phi()>>p2","","")
     
     Perpendicular Component
     EDM=9.02122e-09    STRATEGY= 1      ERROR MATRIX ACCURATE
     EXT PARAMETER                                   STEP         FIRST
     NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
     1  p0           3.79709e+00   7.89306e-02   3.08818e-04   1.33880e-03
     2  p1           1.02242e+00   2.34708e-02   5.06172e-05  -1.50613e-03
     3  p2           2.90239e-02   9.78729e-02   2.20350e-04   1.61557e-03


    x'=x \cos \theta - y \sin \theta
    y'=x \sin \theta + y \cos \theta 
    */
   /*
     Attempt 2: No MET cut
     METparaVsJPhi_pfx
     [0]*sin([1]*x)+[2]
     FCN=61.0499 FROM MIGRAD    STATUS=CONVERGED      60 CALLS          61 TOTAL
     EDM=3.58299e-07    STRATEGY= 1      ERROR MATRIX ACCURATE
     EXT PARAMETER                                   STEP         FIRST
     NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
     1  p0          -2.98529e+00   5.65000e-02   2.09765e-04  -1.32546e-02
     2  p1           9.90985e-01   9.89128e-03   3.67579e-05  -6.38539e-02
     3  p2          -1.44748e+01   3.87877e-02   1.49187e-04  -1.50330e-03

     METperpVsJPhi_pfx
     [0]*cos([1]*x)
     FCN=82.6718 FROM MIGRAD    STATUS=CONVERGED      83 CALLS          84 TOTAL
     EDM=1.7879e-07    STRATEGY= 1      ERROR MATRIX ACCURATE
     EXT PARAMETER                                   STEP         FIRST
     NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
     1  p0           2.97787e+00   5.01198e-02   2.13845e-04   1.11023e-02
     2  p1           1.01031e+00   1.02443e-02   4.36243e-05  -1.06512e-02
    */
   
   else if (eraType.CompareTo("pfMEtSysShiftCorrParameters_2012runAvsJPhi_TAMUWW_mc")==0) {
      if(jp4.size()>0) {
         //METCor.first  = (-3.75112*TMath::Sin(0.988865*jp4[0].Phi()) * TMath::Cos(jp4[0].Phi())) + (3.79709*TMath::Cos(1.02242*jp4[0].Phi()) * TMath::Sin(jp4[0].Phi()));
         //METCor.second = (-3.75112*TMath::Sin(0.988865*jp4[0].Phi()) * TMath::Sin(jp4[0].Phi())) - (3.79709*TMath::Cos(1.02242*jp4[0].Phi()) * TMath::Cos(jp4[0].Phi()));
         METCor.first  = (-2.98529*TMath::Sin(0.990985*jp4[0].Phi()) * TMath::Cos(jp4[0].Phi())) + (2.97787*TMath::Cos(1.01031*jp4[0].Phi()) * TMath::Sin(jp4[0].Phi()));
         METCor.second = (-2.98529*TMath::Sin(0.990985*jp4[0].Phi()) * TMath::Sin(jp4[0].Phi())) - (2.97787*TMath::Cos(1.01031*jp4[0].Phi()) * TMath::Cos(jp4[0].Phi()));
      }
      else {
         METCor.first  = 0;
         METCor.second = 0;
      }
   }
   else {
      METCor.first  = 0;
      METCor.second = 0;
   }
   
   METCor.first*=-1;
   METCor.second*=-1;

   return  METCor;
}

//______________________________________________________________________________
void PerformSelection::doMETPhiCorrection(TString eraType_npv, TString eraType_jp) {
   //Loop over the MET lorentz vectors
   for (unsigned int m=0; m < METp4.size(); m++){
      
      // get the correction factor for this MET
      pair<double,double> cor_npv = getMETPhiCorrection_NPV(eraType_npv);
      pair<double,double> cor_jp = getMETPhiCorrection_JetPhi(eraType_jp);
      //cout << "\tPx before correction = " << METp4[0].Px() << endl
      //     << "\tPy before correction = " << METp4[0].Py() << endl
      //     << "\tPx_jp correction = " << cor.first << endl
      //     << "\tPy_jp correction = " << cor.second << endl
      //     << "\tPx_npv correction = " << cor.first << endl
      //     << "\tPy_npv correction = " << cor.second << endl;
      METp4[0].SetX(METp4[0].Px()+cor_jp.first+cor_npv.first);
      METp4[0].SetY(METp4[0].Py()+cor_jp.second+cor_npv.second);
      //cout << "\tPx after correction = " << METp4[0].Px() << endl
      //     << "\tPy after correction = " << METp4[0].Py() << endl;
   }// for METs
}


////////////////////////////////////////////////////////////////////////////////
// implement additional functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
void PerformSelection::getCollections(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
   iEvent.getByLabel(triggerSource, triggerHandle);
   assert ( triggerHandle.isValid() );

   iEvent.getByLabel(vtxSource, vtxHandle);
   assert ( vtxHandle.isValid() ); 

   if ( !Data ) {
      iEvent.getByLabel(genParticleSource, genParticleHandle);
      assert ( genParticleHandle.isValid() );

      iEvent.getByLabel(pileupSource, pileupHandle);
      assert ( pileupHandle.isValid() );
   }

   iEvent.getByLabel(pfJetSource, pfJetHandle);
   assert ( pfJetHandle.isValid() );

   iEvent.getByLabel(electronSource, electronHandle);
   assert ( electronHandle.isValid() );
 
   iEvent.getByLabel(muonSource, muonHandle);
   assert ( muonHandle.isValid() );

   iEvent.getByLabel(METSource, METHandle);
   assert ( METHandle.isValid() );

   iEvent.getByLabel(rhoSource, rhoHandle);
   assert ( rhoHandle.isValid() );
}

//______________________________________________________________________________
void PerformSelection::makeTPUDist() {
   TPUDist->Fill((*pileupHandle)[1].getTrueNumInteractions());
}

//______________________________________________________________________________
void PerformSelection::setDRlj() {
   for (int i=0; i<jcnt_tot; i++)
      if (lp4.size()>0)
         jp4[i].DRlj = lp4[0].DeltaR(jp4[i]);
}

//______________________________________________________________________________
void PerformSelection::setThetalj1pj2() {
   if (jp4.size()>=2 && lp4.size()>0) {
      lp3.SetXYZ(lp4[0].Px(),lp4[0].Py(),lp4[0].Pz());
      jjp3.SetXYZ(jp4[0].Px()+jp4[1].Px(),jp4[0].Py()+jp4[1].Py(),jp4[0].Pz()+jp4[1].Pz());
      lp4[0].Thetalj1pj2 = lp3.Dot(jjp3)/(lp3.Mag()*jjp3.Mag());
   }
   else if (lp4.size()>0)
      lp4[0].Thetalj1pj2 = 0;
   else
      return;
}

//______________________________________________________________________________
void PerformSelection::saveGenPart() {
   EvtNtuple->genParticleCollection.clear();
   map<size_t,pair<int,int> > memMap;
   int counter = 0;
   for (reco::GenParticleCollection::const_iterator genParticle = genParticleHandle->begin();genParticle != genParticleHandle->end();++genParticle, ++counter) {
      memMap[size_t((const reco::Candidate*)(& *genParticle))] = make_pair(counter, 9999);
      if (genParticle->status() == particleStatus) {
         GenParticle g;
         g.charge = genParticle->charge();
         math::XYZTLorentzVectorD const& p4t = genParticle->p4();
         TLorentzVector p4(p4t.px(), p4t.py(), p4t.pz(), p4t.energy());
         g.p4 = p4;
         g.vtx.SetXYZ(genParticle->vertex().X(),genParticle->vertex().Y(),genParticle->vertex().Z());
         g.pdgId = genParticle->pdgId();
         g.status = genParticle->status();
         g.particlePosition = size_t((const reco::Candidate*)(& *genParticle));
         g.numberOfMothers = genParticle->numberOfMothers();
         g.numberOfDaughters = genParticle->numberOfDaughters();
         for (unsigned int i = 0; i<genParticle->numberOfMothers(); i++) {
            g.motherPositions.push_back(size_t(genParticle->mother(i)));
         }
         for (unsigned int i = 0; i<genParticle->numberOfDaughters(); i++) {
            g.daughterPositions.push_back(size_t(genParticle->daughter(i)));
         }
         EvtNtuple->genParticleCollection.push_back(g);
         memMap[size_t((const reco::Candidate*)(& *genParticle))].second = EvtNtuple->genParticleCollection.size() - 1;
         //cout << "Particle Position = " << memMap[size_t((const reco::Candidate*)(& *genParticle))].second << endl;
         //cout << "Particle pointer = " << size_t((const reco::Candidate*)(& *genParticle)) << end;
         //cout << "Particle key = " << EvtNtuple->genParticleCollection.particlePosition << endl;
         //cout << "EvtNtuple->genParticleCollection.size() = " << EvtNtuple->genParticleCollection.size() << endl;
/*
          for (reco::GenParticleCollection::const_iterator genParticle2 = genParticleHandle->begin();genParticle2 != genParticleHandle->end();++genParticle2) {
             if (genParticle->mother(0)==(const reco::Candidate*)(& *genParticle2))
             //if (matchGenParticles(genParticle->mother(0),(const reco::Candidate*)(& *genParticle2)))
                cout << "Found mother match at: " << genParticle->mother(0) << endl;
          }
*/
      }
   }
   for(unsigned int i=0; i< EvtNtuple->genParticleCollection.size(); i++) {
      //cout << "Particle Position (" << EvtNtuple->genParticleCollection[i].particlePosition << "," << memMap[EvtNtuple->genParticleCollection[i].particlePosition].second << ")" << endl;
      EvtNtuple->genParticleCollection[i].particlePosition = memMap[EvtNtuple->genParticleCollection[i].particlePosition].second;
      for(unsigned int j=0; j<EvtNtuple->genParticleCollection[i].motherPositions.size(); j++) {
         //cout << "Expected Mother Position  = " << memMap[EvtNtuple->genParticleCollection[i].motherPositions[j]].second << endl;
         EvtNtuple->genParticleCollection[i].motherPositions[j] = memMap[EvtNtuple->genParticleCollection[i].motherPositions[j]].second; 
         //cout << "Actual Mother Position  = " << EvtNtuple->genParticleCollection[i].motherPositions[j] << endl;
      }
      for(unsigned int j=0; j<EvtNtuple->genParticleCollection[i].daughterPositions.size(); j++) {
         //cout << "Expected Daughter Position  = " << memMap[EvtNtuple->genParticleCollection[i].daughterPositions[j]].second << endl;
         EvtNtuple->genParticleCollection[i].daughterPositions[j] = memMap[EvtNtuple->genParticleCollection[i].daughterPositions[j]].second;
         //cout << "Actual Daughter Position  = " << EvtNtuple->genParticleCollection[i].daughterPositions[j] << endl;
      }
   }
}

//______________________________________________________________________________
bool PerformSelection::matchGenParticles(const reco::Candidate* p1, const reco::Candidate* p2) {
   if (p1->status()==p2->status() && p1->pdgId()==p2->pdgId() && p1->charge()==p2->charge()) {
      math::XYZTLorentzVectorD const& tl1p4 = p1->p4();
      math::XYZTLorentzVectorD const& tl2p4 = p2->p4();
      if (tl1p4.px()==tl2p4.px() && tl1p4.py()==tl2p4.py() && tl1p4.pz()==tl2p4.pz() && tl1p4.energy()==tl2p4.energy())
         return true;
   }
   return false;
}

//______________________________________________________________________________
pair<int, TLorentzVector> PerformSelection::matchToGen(double eta, double phi) {
   double bestDeltaR = 9999;
   int bestPDGID = 0;
   TLorentzVector theGenObject(0,0,0,0);
   TLorentzVector theBestGenObject(0,0,0,0);
   for (reco::GenParticleCollection::const_iterator genParticle = genParticleHandle->begin();genParticle != genParticleHandle->end();++genParticle) {
      if (genParticle->pt()==0)
         continue;
      if (genParticle->status() == particleStatus) {
         theGenObject.SetPxPyPzE(genParticle->px(),genParticle->py(),genParticle->pz(),genParticle->energy());
         if (theGenObject.Pt() > 1.0e-6 && reco::deltaR(theGenObject.Eta(), theGenObject.Phi(), eta, phi) < bestDeltaR) {
            bestDeltaR = reco::deltaR(theGenObject.Eta(), theGenObject.Phi(), eta, phi);
            theBestGenObject = theGenObject;
            bestPDGID = genParticle->pdgId();
         }
      }
   }
   pair<int, TLorentzVector> GenMatchedInformation(bestPDGID,theBestGenObject);
   return GenMatchedInformation;
}

//______________________________________________________________________________
void PerformSelection::fillJetMap(map<Int_t, Int_t> & jetMap, bool jets){
   typedef map<double, pair<Int_t, Int_t> > ITJ;
   TLorentzVector theGenObject(0,0,0,0);

   // Create an aux map with the dR of all the possible
   // combinations of jets in both events
   map<double, pair<Int_t, Int_t> >  auxMap;
   for (unsigned int recIndex = 0; recIndex<jp4.size(); recIndex++) {
      bool matchMade = false;
      int genIndex = -1;
      for (reco::GenParticleCollection::const_iterator genParticle = genParticleHandle->begin();genParticle != genParticleHandle->end();++genParticle) {
         genIndex++;
         if (genParticle->pt()==0)
            continue;
         if (jets && (fabs(genParticle->pdgId())>6 && fabs(genParticle->pdgId())!=21))
            continue;
         /*
         bool skip = false;
         for (unsigned int i = 0; i<genParticle->numberOfMothers(); i++) {
            //g.motherPositions.push_back(size_t(genParticle->mother(i)));
            //cout << "Part pdgId = " << genParticle->pdgId() << "\tMother pdgId=" << genParticle->mother(i)->pdgId() << endl;

            //if (genParticle->mother(i)->pdgId()==2212)
            //   continue;
            if (genParticle->mother(i)->pdgId()==genParticle->pdgId()){
               //cout << "WARNING::PerformSelection::fillJetMap Both genParticle and mother have same pdgId (" << genParticle->pdgId() << ")" << endl;
               skip = true;
            }
         }
         if(skip)
            continue;
         */
         if (genParticle->status() == particleStatus) {
            theGenObject.SetPxPyPzE(genParticle->px(),genParticle->py(),genParticle->pz(),genParticle->energy());
            if(theGenObject.Pt()>1.0e-6) {
               double dR = reco::deltaR(theGenObject.Eta(), theGenObject.Phi(), jp4[recIndex].Eta(), jp4[recIndex].Phi());
               //if(auxMap.find(dR) != auxMap.end())
               //   cout << "PerformSelection::fillJetMap::WARNING Duplicate DeltaR when doing matching." << endl;
               auxMap[dR] = std::make_pair(recIndex, genIndex);
               matchMade = true;
               //if (recIndex == 1)
               //   cout << "RecIndex == 1 and a match made and entered into auxMap[" << dR << "] = std::make_pair(" << recIndex << "," << genIndex << ")" << endl;
            }
         }
      }
      if (!matchMade) {
         cout << "No match" << endl;
         auxMap[-1-recIndex] = std::make_pair(recIndex,-1-recIndex);
      }
   }

   // First clear the map for this new set of events
   jetMap.clear();

   // 1-Find the pair of jets with the smallest dr. Add them to resulting map
   // 2-Remove the pair from the map, rinse and repeat.
   while (auxMap.size() > 0){

      // 1- The first element which is the one with the smallest dR. Get the jet indexes
      int j1 = auxMap.begin()->second.first;
      int j2 = auxMap.begin()->second.second;

      // Add to the results
      //if (auxMap.begin()->first < 0.25 && fabs(jp4[j1].refLV.Pt()-(genParticleHandle->begin()+j2)->pt())<1)
      jetMap[j1] = j2;
      //cout << "RecJet = " << j1 << "\tGenJet = " << j2 << endl;
      
      // 2- Now remove all elements from the auxMap that contain either the first or second jet
      ITJ::iterator itr = auxMap.begin();
      while(itr != auxMap.end()){
         if (itr->second.first == j1 || itr->second.second == j2)
            auxMap.erase(itr++);
         else
            ++itr;

      }//while removing

   }//while

   //Double check that each rec jet has a gen jet match. If not, put a placeholder in.
   for (unsigned int iRec = 0; iRec<jp4.size(); iRec++) {
      if (jetMap.find(iRec)!=jetMap.end())
         continue;
      else {
         cout << "PerformSelection::fillJetMap::WARNING Did not find gen particle match for reco jet index " << iRec << endl;
         jetMap[iRec] = -1-iRec;
      }
   }
}//fillJetMap

//______________________________________________________________________________
vector<int> PerformSelection::matchToGen_particleCollection(bool jets) {
   typedef map<Int_t, Int_t> ITJ;
   map<Int_t, Int_t> jetMap;
   fillJetMap(jetMap, jets);
   vector<int> ret(jetMap.size(),0);

   for (ITJ::const_iterator itj = jetMap.begin(); itj != jetMap.end(); itj++) {
      int jrec = itj->first;
      int jgen = itj->second;

      if (jgen>=0)
         ret[jrec] = (genParticleHandle->begin()+jgen)->pdgId();
      else {
         ret[jrec] = 0;
      }
   }

   //if(jetMap.size()!=jp4.size()) cout <<"PerformSelection::matchToGen_particleCollection::WARNING:: sizes don't match!!!"  << endl
   //                                << "recJets.size()=" << jp4.size() << "\t matches.size()=" << jetMap.size() << endl;
   //if(ret.size()!=jp4.size()) cout <<"PerformSelection::matchToGen_particleCollection::WARNING:: sizes don't match!!!"  << endl
   //                                << "recJets.size()=" << jp4.size() << "\t matches.size()=" << ret.size() << endl;
   return ret;   
}

//______________________________________________________________________________
double PerformSelection::getJERfactor(double pt, double eta, double ptgen){

  double jer = 1;

  if (fabs(eta) < 0.5)
    jer = 1.052;
  else if (fabs(eta) < 1.1)
    jer = 1.057;
  else if (fabs(eta) < 1.7)
    jer = 1.096;
  else if (fabs(eta) < 2.3)
    jer = 1.134;
  else if (fabs(eta) < 5)
    jer = 1.288;

  double corr = ptgen / pt;

  return  max(0.0,corr + jer * (1 - corr));

}

//______________________________________________________________________________
void PerformSelection::ptSort(vector<Jet> & vec) {
   for (unsigned int next = 0; next < vec.size()-1; next++) {
      unsigned int max_pos = (unsigned int)max_position(vec, next, vec.size()-1);
      if (max_pos != next)
         swap(vec[max_pos],vec[next]);
   }
}

//______________________________________________________________________________
void PerformSelection::swap(TLorentzVector & x, TLorentzVector & y) {
   TLorentzVector temp = x;
   x = y;
   y = temp;
}

//______________________________________________________________________________
int PerformSelection::max_position(vector<Jet> & vec, int from, int to) {
   int max_pos = from;
   int i;
   for (i = from + 1; i <= to; i++)
      if (vec[i].Pt() > vec[max_pos].Pt()) max_pos = i;
   return max_pos;
}

//______________________________________________________________________________
int PerformSelection::getBin(double x, const vector<double> boundaries, int length) {
   int i;
   int n = length;
   if (n<=0) return -1;
   if (x<boundaries[0] || x>=boundaries[n])
      return -1;
   for(i=0;i<n;i++)
   {
      if (x>=boundaries[i] && x<boundaries[i+1])
         return i;
   }
   return 0;
}

//______________________________________________________________________________
void PerformSelection::incrementCounter(int nCut, unsigned int nJets, Table* t1, Table* t2, Table* t3) {
   (*(*t1)(DEFS::getCutLevelString((DEFS::CutLevel)nCut),DEFS::getJetBinString((DEFS::JetBin)nJets)))++;
   (*(*t2)(DEFS::getCutLevelString((DEFS::CutLevel)nCut),DEFS::getJetBinString((DEFS::JetBin)nJets)))++;   
   if (t3 != 0)
      (*(*t3)(DEFS::getCutLevelString((DEFS::CutLevel)nCut),DEFS::getJetBinString((DEFS::JetBin)nJets)))++;   
}

//______________________________________________________________________________
void PerformSelection::printEventInformation(bool print, int cLevel, bool muon) {
   if (!print) return;
   if (cLevel == -1)
      cout << "Event " << eventNumber << " NJets " << jcnt_tot << " Primary Muons " << mucnt_Prim
           << " Primary Electrons " << elcnt_Prim << " Loose Muons " << mucnt_Loose-mucnt_Prim << " Loose Electrons " 
           << elcnt_Loose-elcnt_Prim;
   else if (cLevel == 0)
      cout << " c0";
   else if (cLevel == 1)
      cout << " c1";
   else if (cLevel == 2)
      cout << " c2";
   else if (cLevel == 3) {
      if (muon)
         cout << " c3";  
      else if (!muon)
         cout << " c3";  
      else
         cout << " WARNING::Unknown primary lepton type!";
   }
   else if (cLevel == 4)
      cout << " c4";
   else if (cLevel == 5)
      cout << " c5";
   else if (cLevel == 6)
      cout << " c6";
   else if (cLevel == 7)
      cout << endl;
   else
      cout << "\tWARNING::Unknown cut level!" << endl;
}

//______________________________________________________________________________
void PerformSelection::printJetInformation() {
   for(unsigned int i=0; i<jstreams.size(); i++)
      cout << "Jet " << i << " " << jstreams[i];
}

//______________________________________________________________________________
void PerformSelection::printLeptonInformation() {
   for(unsigned int i=0; i<lstreams.size(); i++)
      cout << lstreams[i];
}


////////////////////////////////////////////////////////////////////////////////
// define PerformSelection as a plugin
////////////////////////////////////////////////////////////////////////////////

DEFINE_FWK_MODULE(PerformSelection);
