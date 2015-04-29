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
               qqH550 , qqH600 , qqH700 , qqH800 , qqH900 , qqH1000 , qqH125_JESUp, qqH125_JESDown,
               WH_ZH_TTH_HToZZ_M125, WH125_HToBB, WH125_HToZG, WH_ZH_TTH_HToWW_M125, TTH_Inclusive_M125, TTH_HToBB_M125, WH_HToBB_M125_JESUp, WH_HToBB_M125_JESDown, WH_ZH_TTH_HToWW_M125_JESUp, WH_ZH_TTH_HToWW_M125_JESDown, WH_ZH_TTH_HToZZ_M125_JESUp, WH_ZH_TTH_HToZZ_M125_JESDown, 
               WH_HToZZ_M125, ZH_HToZZ_M125, TTH_HToZZ_M125, WH_HToWW_M125, ZH_HToWW_M125, TTH_HToWW_M125, 
               STopS_T , STopS_Tbar , STopT_T , STopT_Tbar , STopTW_T , STopTW_Tbar , TTbar , TTbarLJ, TTbarDil , TTbar_JESUp, TTbar_JESDown,
               Wcc , WJets , WJets_part2 , W1Jets, W2Jets, W3Jets, W4Jets, WLg , Wgg , WLL , WLb , Wbb , WJets_JESUp, WJets_JESDown, WJets_matchingup, WJets_matchingdown, WJets_scaleup, WJets_scaledown, 
               WW , WZbb , WZ , ZZ , ZJets , Ztautau , 
               QCD_ElEnriched , QCD_ElFULL, QCD_MuEnriched , QCD_MuFULL ,
	             QCD_Pt20to30_EMEnriched, QCD_Pt30to80_EMEnriched, QCD_Pt80to170_EMEnriched, 
	             QCD_Pt170to250_EMEnriched, QCD_Pt250to350_EMEnriched, QCD_Pt350_EMEnriched,
               SingleEl_Data , SingleMu_Data, 
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

  // A routine that returns a JetBin given
  JetBin getJetBin(std::string str);

  enum NBinsX {nbinsx_default=100, nbinsx_jets2=59, nbinsx_jets3=45, nbinsx_jets4=41};

  // A routing which returns the number of x bins based on the JetBin
  int getNBinsX(JetBin type);

  // A routine that returns an array of the bin bondaries for the 1-D BDT plots
  Double_t* getBinsX(JetBin type);

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

   enum ControlRegion {all, signal, control1, control2, control3, control4, control5, control6, control7, control8, control9, UVaCuts, event, Diboson, MVAEleID, AntiMVAEleID, FlatMVAEleID, None};
  static const unsigned int nControlRegion = 18;

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
     enum University {TAMU, UVa, OtherUniversity};
     static const unsigned int nUniversity = 3;

     // A routine that returns a string given the University
     std::string getUniversityString(University univ);

     // A routine that returns a University given a string
     University getUniversity(std::string str);

}

#endif

