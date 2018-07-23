#ifndef DEFS_DEF
#define DEFS_DEF

//ROOT libraries
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

//C++ libraries
#include <iostream>
#include <string>
#include <map>
#include <vector>

// This namespace holds the definitions of all physics processes 
// author: Ricardo Eusebi, Feb 12, 2009
// modified: Osipenkov, Ilya

namespace DEFS{

  // ---------------------------------------------------------------
  //            ALL ABOUT THE ANALYSIS WE WANT TO DO 
  // ---------------------------------------------------------------  
  namespace Ana {
    enum Type{
      SingleTopAnalysis,
      HiggsAnalysis,
      WWAnalysis,
      UNKNOWN
    };
  
    // A routine that returns the string given the type 
    std::string getTypeString(Type);
 
     // A routine that returns the type given a string
    Type getAnaType(std::string str);

  }// Ana namespace
  
  // ---------------------------------------------------------------
  //            ALL ABOUT THE PHYSICSPROCESS TYPES
  // ---------------------------------------------------------------
  namespace PhysicsProcess {
    enum Type {WH100 , WH105 , WH110 , WH115 , WH120 , WH125 , WH130 , WH135 , WH140 , WH145 , WH150 ,
               WH160 , WH170 , WH180 , WH190 , WH200 , WH250 , WH300 , WH350 , WH400 , WH450 , WH500 ,
               WH550 , WH600 , WH700 , WH800 , WH900 , WH1000 , 
               ggH100 , ggH105 , ggH110 , ggH115 , ggH120 , ggH125 , ggH130 , ggH135 , ggH140 , ggH145 , ggH150 ,
               ggH160 , ggH170 , ggH180 , ggH190 , ggH200 , ggH250 , ggH300 , ggH350 , ggH400 , ggH450 , ggH500 ,
               ggH550 , ggH600 , ggH700 , ggH800 , ggH900 , ggH1000 , 
               qqH100 , qqH105 , qqH110 , qqH115 , qqH120 , qqH125 , qqH130 , qqH135 , qqH140 , qqH145 , qqH150 ,
               qqH160 , qqH170 , qqH180 , qqH190 , qqH200 , qqH250 , qqH300 , qqH350 , qqH400 , qqH450 , qqH500 ,
               qqH550 , qqH600 , qqH700 , qqH800 , qqH900 , qqH1000 ,
               WH_ZH_TTH_HToZZ_M125, WH125_HToBB, WH125_HToZG, WH_ZH_TTH_HToWW_M125, TTH_Inclusive_M125, TTH_HToBB_M125,
               WH_HToZZ_M125, ZH_HToZZ_M125, TTH_HToZZ_M125, WH_HToWW_M125, ZH_HToWW_M125, TTH_HToWW_M125, 
               STopS_T , STopS_Tbar , STopT_T , STopT_Tbar , STopTW_T , STopTW_Tbar , TTbar , TTbarLJ, TTbarDil ,
               Wcc , WJets , WJets_part2 , W1Jets, W2Jets, W3Jets, W4Jets, WLg , Wgg , WLL , WLb , Wbb ,
               WW , WZbb , WZ , ZZ , ZJets , Ztautau , ZJetsToLL_M50 , ZJetsToLL_M10To50,
               QCD_ElEnriched , QCD_ElFULL, QCD_MuEnriched , QCD_MuFULL ,
	             QCD_Pt20to30_EMEnriched, QCD_Pt30to80_EMEnriched, QCD_Pt80to170_EMEnriched, 
	             QCD_Pt170to250_EMEnriched, QCD_Pt250to350_EMEnriched, QCD_Pt350_EMEnriched,
               SingleEl_Data , SingleMu_Data,
               SingleEl_ZJetsToLL , SingleMu_ZJetsToLL,
               SingleEl_WlnuJets , SingleMu_WlnuJets,
               WZ_JESUp, WZ_JESDown, WW_JESUp, WW_JESDown,
               STopS_T_JESUp , STopS_Tbar_JESUp , STopT_T_JESUp , STopT_Tbar_JESUp , STopTW_T_JESUp , STopTW_Tbar_JESUp ,
               STopS_T_JESDown , STopS_Tbar_JESDown , STopT_T_JESDown , STopT_Tbar_JESDown , STopTW_T_JESDown , STopTW_Tbar_JESDown ,
               TTbar_JESUp, TTbar_JESDown,
               ZJets_JESUp, ZJets_JESDown,
               WJets_JESUp, WJets_JESDown, WJets_matchingup, WJets_matchingdown, WJets_scaleup, WJets_scaledown,
               WH_ZH_TTH_HToWW_M125_JESUp, WH_ZH_TTH_HToWW_M125_JESDown, WH_ZH_TTH_HToZZ_M125_JESUp, WH_ZH_TTH_HToZZ_M125_JESDown, 
               WH_HToWW_M125_JESUp, WH_HToWW_M125_JESDown, WH_HToZZ_M125_JESUp, WH_HToZZ_M125_JESDown,
               ZH_HToWW_M125_JESUp, ZH_HToWW_M125_JESDown, ZH_HToZZ_M125_JESUp, ZH_HToZZ_M125_JESDown,
               TTH_HToWW_M125_JESUp, TTH_HToWW_M125_JESDown, TTH_HToZZ_M125_JESUp, TTH_HToZZ_M125_JESDown,
               WH_HToBB_M125_JESUp, WH_HToBB_M125_JESDown, TTH_HToBB_M125_JESUp, TTH_HToBB_M125_JESDown,
               ggH125_JESUp, ggH125_JESDown, qqH125_JESUp, qqH125_JESDown, 
               UNKNOWN};

    // A routine that returns the type given a string 
    Type getProcessType(std::string str);///
    
    // A routine that returns the type given a string 
    std::string getTypeString(Type );
    
    // A routine that tells whether this process is Higgs or not
    bool isHiggs(Type type);

    // A routine that returns the title for a given process type
    std::string getTypeTitle(Type type);

    // A  routine that returns the Color_t for a given process type
    Color_t getProcessColor(Type type);
    
  }// PhysicsProcess namespace
  typedef DEFS::PhysicsProcess::Type PhysicsProcessType ;

  // ---------------------------------------------------------------
  //            ALL ABOUT THE DIFFERENT EVENT CATEGORIES
  // ---------------------------------------------------------------
  // Under no condition set the value of the category, otherwise nLeptonCat 
  // will not be working on some classes like PhysicsProcessForOpt
  // In addition, make sure that all the categories that data/MC has are here
  enum LeptonCat {none, muon, electron, both};
  static const unsigned int nLeptonCat = 4; 

  // A routine that returns the string given the LeptonCat 
  std::string getLeptonCatString(LeptonCat );

  // for historical reasons
  std::string getEventCatString(LeptonCat a);

  // A routine that returns the LeptonCat given the string
  LeptonCat getLeptonCat(std::string type);

  // ---------------------------------------------------------------
  //            ALL ABOUT THE JET TYPES
  // ---------------------------------------------------------------
  enum JetBin {jets0, jet1, jets2, jets3, jets4, jets5};
  
  // for historical reasons
  static const unsigned int NJETS = 5;

  // A routine that returns a string given the type
  std::string getJetBinString(JetBin type);

  // A routine that returns a label given the type
  std::string getJetBinLabel(JetBin type);

  // A routine that returns a JetBin given
  JetBin getJetBin(std::string str);

  // A routine that returns a jetBin given an integer number of jets
  JetBin getJetBin(int nJets, JetBin inclusive_bin = jets4);

  //With Signal cuts
    enum NBinsX {nbinsx_default=100, nbinsx_jets2_electron_KinBDT=45, nbinsx_jets2_electron_MEBDT=27, nbinsx_jets2_electron_KinMEBDT=37,
               nbinsx_jets2_muon_KinBDT=47, nbinsx_jets2_muon_MEBDT=23, nbinsx_jets2_muon_KinMEBDT=37,//nbinsx_jets2_muon_MEBDT=30
               nbinsx_jets3_electron_KinBDT=34, nbinsx_jets3_electron_MEBDT=25, nbinsx_jets3_electron_KinMEBDT=29,
               nbinsx_jets3_muon_KinBDT=34, nbinsx_jets3_muon_MEBDT=25, nbinsx_jets3_muon_KinMEBDT=29,
               nbinsx_jets4_electron_KinBDT=25, nbinsx_jets4_electron_MEBDT=20, nbinsx_jets4_electron_KinMEBDT=20,
               nbinsx_jets4_muon_KinBDT=25, nbinsx_jets4_muon_MEBDT=19, nbinsx_jets4_muon_KinMEBDT=20};
  //With HailMary cuts
  /*
  enum NBinsX {nbinsx_default=100, nbinsx_jets2_electron_KinBDT=30, nbinsx_jets2_electron_MEBDT=20, nbinsx_jets2_electron_KinMEBDT=27,
               nbinsx_jets2_muon_KinBDT=30, nbinsx_jets2_muon_MEBDT=18, nbinsx_jets2_muon_KinMEBDT=25,
               nbinsx_jets3_electron_KinBDT=21, nbinsx_jets3_electron_MEBDT=14, nbinsx_jets3_electron_KinMEBDT=17,
               nbinsx_jets3_muon_KinBDT=18, nbinsx_jets3_muon_MEBDT=13, nbinsx_jets3_muon_KinMEBDT=16,
               nbinsx_jets4_electron_KinBDT=8, nbinsx_jets4_electron_MEBDT=7, nbinsx_jets4_electron_KinMEBDT=8,
               nbinsx_jets4_muon_KinBDT=6, nbinsx_jets4_muon_MEBDT=5, nbinsx_jets4_muon_KinMEBDT=5};
  */
  // A routing which returns the number of x bins based on the JetBin
  int getNBinsX(JetBin type, LeptonCat lcat, std::string variable);

  // A routine that returns an array of the bin bondaries for the 1-D BDT plots
  std::vector<Double_t> getBinsX(JetBin type, LeptonCat lcat, std::string variable);

  // ---------------------------------------------------------------
  //            ALL ABOUT THE TAGGING CATEGORIES
  // ---------------------------------------------------------------

  // The tagging categories per event
  // TSV = tight SecVtx
  // JP5 = JP < 5%
  enum TagCat {

    pretag,   // The pretag category
    eq0tag,   // Exactly zero SVX tags, i.e. the untag category.
    eq1tag,   // Exactly one SVX tag 
    eq2tag,   // Exactly two SVX tags 
    ge0tag,   // Greater equal than 0 SVX tags 
    ge1tag,   // Greater equal than 1 SVX tag 
    ge2tag,   // Greater equal than 2 SVX tags 

    // For the tagging with JP in descendant exclusive orthogonal category
    TSVTSV,   // for some reason this is not exactly equal to eq2TSV
    TSVJP5,   // one jet with a tight SVX tag and the other one with a loose JP tag which is not a tight SVX
    TSVNOJP5, // one jet with a tight SVX tag and the other one without a loose JP or tight SVX tag.
    JP5NOTSV // one jet with a Loose SVT tag and neither jet is a tight SVX tag.
   

    // For the tagging with Loose in descendant exclusive orthogonal category
    // TSVTSV  //  for some reason this is not exaclty equal to eq2TSV
    ///TSVLSV,    // one jet with a tight SVX tag and the other one with a loose SecVTx tag which is not a tight SVX
    //TSVNOLSV,  // one jet with a tight SVX tag and the other one without a loose SVX tag.   
    //LSVNOTSV   // one jet with a Loose SVT tag and neither jet is a tight SVX tag.

    // In principle we can have four choices, each with 4 exclusive categories in which in a 
    // given column a given line is always orthogonal to all lines above:
    //       | CHOICE A    | CHOICE B  |  CHOICE C | CHOICE D
    // cat 0 | TSVTSV      | TSVTSV    | TSVTSV    | TSVTSV   
    // cat 1 | TSVJP5      | TSVLSV    | TSVJP5    | TSVLSV   
    // cat 2 | TSVNOJP5    | TSVNOLSV  | TSVNOJP5  | TSVNOLSV 
    // cat 3 | JP5NOTSV    | LSVNOTSV  | LSVNOTSV  | JP5NOTSV    
    
  };

  // A routine that returns a string given the type
  std::string getTagCatString(TagCat type);

  // A routine that returns a label given the type
  std::string getTagCatLabel(TagCat type);

  // A routine that returns a jetBin given a string
  TagCat getTagCat(std::string str);

  // ---------------------------------------------------------------
  //            ALL ABOUT THE CUT LEVELS
  // ---------------------------------------------------------------

   //enum CutLevel {NPATtupleEvts, c0, c1, c2, c3, c4, c5, c6, FNAL1, FNAL2, FNAL3, FNAL4, FNAL5, FNAL6, FNAL7, FNAL8, FNAL9, FNAL10, BTag0, BTag1, BTag2, BTag3p};
   //static const unsigned int nCutLevel = 22;
   enum CutLevel {NPATtupleEvts, c0, c1, c2, c3, c4, c5, c6, ac1, ac2, ac3, ac4, ac5, ac6, BTag0, BTag1, BTag2, BTag3p};
   static const unsigned int nCutLevel = 18;
   static const unsigned int nFinalCutLevel = 7;

  // A routine that returns a string given the type
  std::string getCutLevelString(CutLevel type);

  // A routine that returns a CutLevel given a string
  CutLevel getCutLevel(std::string str);

  // ---------------------------------------------------------------
  //            ALL ABOUT THE CONTROL REGIONS
  // ---------------------------------------------------------------

   enum ControlRegion {all, signal, control1, control2, control3, control4, control5, control6, control7, control8, control9, UVaCuts, event, Diboson, MVAEleID, AntiMVAEleID, FlatMVAEleID, BDTBump, BDTAntiBump, HailMary, HailMaryLoose, LowKinBDT, HighKinBDT, Andrea, None};
  static const unsigned int nControlRegion = 23;

  //A routine that returns a string given the type
  std::string getControlRegionString(ControlRegion type);

  // A routine that returns a CutLevel given a string
  ControlRegion getControlRegion(std::string str);

  // ---------------------------------------------------------------
  //            ALL ABOUT THE NTUPLE TYPES
  // ---------------------------------------------------------------
   enum NtupleType {EventNtuple, METree, MicroNtuple, Other};
   static const unsigned int nNtupleType = 4;
   
   // A routine that returns a string given the type
   std::string getNtupleTypeString(NtupleType type);
   
   // A routine that returns a NtupleType given a string
   NtupleType getNtupleType(std::string str);
   
   std::string getTreeName(NtupleType type, JetBin nJets = jets2);

  // ---------------------------------------------------------------
  //            ALL ABOUT THE UNIVERSITY
  // ---------------------------------------------------------------
   enum University {TAMU, UVa, Combined, OtherUniversity};
   static const unsigned int nUniversity = 4;

   // A routine that returns a string given the University
   std::string getUniversityString(University univ);

   // A routine that returns a University given a string
   University getUniversity(std::string str);

  // ---------------------------------------------------------------
  //            ALL ABOUT PDG ID Codes
  // ---------------------------------------------------------------
   enum PdgId {
     /// Special wildcard particle name
     ANY = 10000, 
     /// @name Charged leptons
     ELECTRON = 11, POSITRON = -ELECTRON, EMINUS = ELECTRON, EPLUS = POSITRON, MUON = 13, ANTIMUON = -MUON, TAU = 15, ANTITAU = -TAU,
     /// @name Neutrinos
     NU_E = 12, NU_EBAR = -NU_E, NU_MU = 14, NU_MUBAR = -NU_MU, NU_TAU = 16, NU_TAUBAR = -NU_TAU, 
     /// @name Bosons
     PHOTON = 22, GAMMA = PHOTON, GLUON = 21, WPLUSBOSON = 24, WMINUSBOSON = -WPLUSBOSON, WPLUS = WPLUSBOSON, WMINUS = WMINUSBOSON, Z0BOSON = 23, ZBOSON = Z0BOSON, Z0 = Z0BOSON, HIGGSBOSON = 25, HIGGS = HIGGSBOSON, 
     /// @name Quarks
     DQUARK = 1, UQUARK = 2, SQUARK = 3, CQUARK = 4, BQUARK = 5, TQUARK = 6, ANTIDQUARK = -DQUARK, ANTIUQUARK = -UQUARK, ANTISQUARK = -SQUARK, ANTICQUARK = -CQUARK, ANTIBQUARK = -BQUARK, ANTITQUARK = -TQUARK, d = DQUARK, u = UQUARK, s = SQUARK, c = CQUARK, b = BQUARK, t = TQUARK, dbar = -d, ubar = -u, sbar = -s, cbar = -c, bbar = -b, tbar = -t, 
     /// @name Nucleons
     PROTON = 2212, ANTIPROTON = -PROTON, PBAR = ANTIPROTON, NEUTRON = 2112, ANTINEUTRON = -NEUTRON, 
     /// @name Light mesons
     PI0 = 111, PIPLUS = 211, PIMINUS = -PIPLUS, K0L = 130, K0S = 310, KPLUS = 321, KMINUS = -KPLUS, ETA = 221, ETAPRIME = 331, PHI = 333, OMEGA = 223, 
     /// @name Charmonia
     ETAC = 441, JPSI = 443, PSI2S = 100443, 
     /// @name Charm mesons
     D0 = 421, DPLUS = 411, DMINUS = -DPLUS, DSPLUS = 431, DSMINUS = -DSPLUS, 
     /// @name Bottomonia
     ETAB = 551, UPSILON1S = 553, UPSILON2S = 100553, UPSILON3S = 200553, UPSILON4S = 300553, 
     /// @name b mesons
     B0 = 511, BPLUS = 521, BMINUS = -BPLUS, B0S = 531, BCPLUS = 541, BCMINUS = -BCPLUS, 
     /// @name Baryons
     LAMBDA = 3122, SIGMA0 = 3212, SIGMAPLUS = 3222, SIGMAMINUS = 3112, LAMBDACPLUS = 4122, LAMBDACMINUS = 4122, LAMBDAB = 5122, XI0 = 3322, XIMINUS = 3312, XIPLUS = -XIMINUS, OMEGAMINUS = 3334, OMEGAPLUS = -OMEGAMINUS, 
     /// @name Exotic/weird stuff
     REGGEON = 110, POMERON = 990, ODDERON = 9990, GRAVITON = 39, NEUTRALINO1 = 1000022, GRAVITINO = 1000039, GLUINO = 1000021
     /// @todo Add axion, black hole remnant, etc. on demand
   };
   static const unsigned int nPdgId = 105;

   /// Return a PdgId as a named string
   std::string toParticleName(PdgId p);

   /// Return a name as a PdgId
   PdgId toParticleId(const std::string pname);

   /// Return a PdgId from an int
   PdgId getParticleIdFromInt(int p);

}

#endif

