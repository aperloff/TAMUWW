#ifndef MICRONTUPLE_HH
#define MICRONTUPLE_HH

#include "TObject.h"
#include "TArrayD.h"
#include "TMath.h"
#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

#include <string>
#include <map>
#include <vector>
#include <sstream>

#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/MEPATNtuple/interface/ProbsForEPD.hh"
#include "TAMUWW/MatrixElement/interface/EventProbDefs.hh"
#include "TAMUWW/SpecialTools/interface/MVAVar.hh"

typedef std::map<double, int> indexMap1;
typedef std::map<DEFS::EP::Type , indexMap1> indexMap2;

class MicroNtuple : public TObject
{
public:


  MicroNtuple(); 
  MicroNtuple(int nJets);
  MicroNtuple(const MicroNtuple&);
  MicroNtuple& operator=(const MicroNtuple&);
  ~MicroNtuple();
      
  enum {nEventProb = 62};

  // ===================================
  //  Member Variables
  // ===================================

  Int_t nJets;
  Double_t eventProb[nEventProb];
  Double_t eventMaxProb[nEventProb];

  //Double_t* bProb; //[nJets]
  std::vector<Double_t> bProb;
  Double_t weight;

  Double_t epd1tag;
  Double_t epd2tag;
  Double_t epd1tagSchan;
  Double_t epd2tagSchan;
  Double_t epd1tagTchan;
  Double_t epd2tagTchan;
  Double_t epdPretagWWandWZ;
  Double_t epd0tagWWandWZ;
  Double_t epd1tagWWandWZ;
  Double_t epd2tagWWandWZ;


  enum {nWHmasses = 13};
  Double_t epd1tagWH[nWHmasses];
  Double_t epd2tagWH[nWHmasses];
 
  enum {nHWWmasses = 34};
  Double_t epd1tagHWW[nHWWmasses];
  Double_t epd2tagHWW[nHWWmasses];

  enum {nHiggsMasses = 49};
  Double_t absorbError[nHiggsMasses];
   //TArrayD epdPretagHiggs;
  Double_t epdPretagHiggs[nHiggsMasses];
  Double_t epd1tagHiggs[nHiggsMasses];
  Double_t epd2tagHiggs[nHiggsMasses];

  enum EPDType {kSchan = 0, kTchan = 1, kCombined = 2};

  static void printIndexMap(indexMap2 & im);
  static void checkIndexMap(indexMap2 & im);
  static void setIndexMap(const indexMap2 & im){indexMap = im;}
  static indexMap2 getIndexMap(){return indexMap;}
  static indexMap2 indexMap;

  // ===================================
  //  Members Methods
  // ===================================

  // This method functions as a map which maps from the eventProbs to the ProbsForEPD.
  // If you don't care about the higgs mass set it to zero. It returns the raw meProbs.
  static ProbsForEPD getEventProbs(double mhiggs, const double evtProbs[nEventProb]) ;

  // This method just calls the above method with evtProbs as a parameter
  ProbsForEPD getEventProbs(double mhiggs) const;

  // Get the probabilities for each process
  static ProbsForEPD getProbsForEPD(const ProbsForEPD & meProbs, 
				    const ProbsForEPD & coeffs,
				    const std::vector<double> bProb,  
				    DEFS::TagCat tagcat);
  
  

  // ===================================
  //  Methods for the particular analyses
  // ===================================

  // Methods for the WH 
  static const ProbsForEPD getWHEPDCoefficients(double mhiggs, DEFS::TagCat tagcat, int nJets); //  WH coefficients 
  Double_t calcWHEPD(DEFS::TagCat tagcat, double mass) const;
  static Double_t calcWHEPD(DEFS::TagCat tagcat,  double mass,
			    const ProbsForEPD & meProbs, 
			    const std::vector<double> bProb, Int_t arraySize);

  // Methods for the H->WW 
  static const ProbsForEPD getHWWEPDCoefficients(double mhiggs, DEFS::TagCat tagcat, int nJets); //  H->WW coefficients 
  Double_t calcHWWEPD(DEFS::TagCat tagcat, double mass) const;
  static Double_t calcHWWEPD(DEFS::TagCat tagcat,  double mass,
			     const ProbsForEPD & meProbs, 
			     const std::vector<double> bProb, Int_t arraySize);

  // Methods for the ggH+qqH+WH
  static const ProbsForEPD getHiggsEPDCoefficients(double mhiggs, DEFS::TagCat tagcat, int nJets); //  ggH+qqH+WH coefficients 
  Double_t calcHiggsEPD(DEFS::TagCat tagcat, double mass) const;
  static Double_t calcHiggsEPD(DEFS::TagCat tagcat,  double mass,
                               const ProbsForEPD & meProbs, 
                               const std::vector<double> bProb, Int_t arraySize);

  // Methods for the Diboson
  static const ProbsForEPD getWZEPDCoefficients(DEFS::TagCat tagcat, int nJets); // Diboson coefficients
  Double_t calcWZEPD(DEFS::TagCat tagcat) const;
  static Double_t calcWZEPD(DEFS::TagCat tagcat, 
			    const ProbsForEPD & meProbs, 
			    const std::vector<double> bProb, Int_t arraySize);

  // Methods for the SingleTop
  static const ProbsForEPD getSingleTopEPDCoefficients (DEFS::TagCat tagcat, EPDType type);
  Double_t calcSingleTopEPD(DEFS::TagCat tagcat, EPDType type = kCombined) const;
  static Double_t calcSingleTopEPD(DEFS::TagCat tagcat, EPDType type,
				   const ProbsForEPD & meProbs,
				   const std::vector<double> bProb, Int_t arraySize);


  // ===================================
  //  Other methods
  // ===================================
  static Double_t triggerTO(double detector, double etalep, double metraw,
			    double etlep);


  void clear();

  // ===================================
  //  TMVA Methods 
  // ===================================
   Int_t size;
   Int_t run;
   Int_t lumi;
   Int_t event;
   Float_t dPhiJetJet;
   Float_t dPhiMETJet;
   Float_t minDPhiMETJet;
   Float_t minDPhiLepJet;
   Float_t Cos_dPhiWW;
   Float_t Cos_dPhiWH;
   Float_t CosTheta_l;
   Float_t CosTheta_j;
   Float_t CosTheta_WH;
   Float_t JacksonAngle;
   Float_t JacobePeak;
   Float_t dRlepjj;
   Float_t sumJetEt;
   Float_t KinBDT;
   Float_t MEBDT;
   Float_t KinMEBDT;
   TMVA::Reader *reader; //!

   // This method sets up the MVA reader.
   // The weights read in by this method must be previously trained.
   // dir is the directory to the weightfiles.
   void setMVAReader(std::vector<TString> MVAMethods, TString dir,
                     std::map<TString,MVAVar>& MVAVars, bool debug = false);

   // This method returns the response, error, probability, and rarity of 
   // the MVA classifier specified by the user.
   // effS is used for the efficiency calculator for cut method.
   // The outputs can be accessed using the map keys "response", "error", "probability", and "rarity".
   std::vector<std::map<TString,Double_t> > getMVAOutput(std::vector<TString> MVAMethods, Double_t effS = 0.7);
   std::vector<Double_t> getMVAOutput(std::vector<TString> MVAMethods, TString specificOutput, Double_t effS = 0.7);
   Double_t getMVAOutput(TString MVAMethod, TString specificOutput, Double_t effS = 0.7);

   // This method returns the efficiency for the CutsGA MVA method.
   // nSelCutsGA is the number of cuts that passed.
   Double_t getCutsGAEfficiency(Int_t nSelCutsGA, Int_t nTreeEntries);

   // Retrieve cuts for particular signal efficiency
   // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer
   void getMVACuts(std::vector<TString> &inputVars, std::vector<TString> &inputLabels, 
                   std::vector<TString> &inputTitles, std::vector<Double_t> &cutsMin, 
                   std::vector<Double_t> &cutsMax, Double_t effS = 0.7);

  ClassDef(MicroNtuple, 19)

};

#endif
