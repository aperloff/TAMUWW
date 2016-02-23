/////////////////////////////////////////////////////////////////////////////////////////////////
//// CMS 
//// WW CrossSection Measurement using Matrix Element 
//// Created by Osipenkov, Ilya : ilyao@fnal.gov
/////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EVENTNTUPLE_HH
#define EVENTNTUPLE_HH

//
// CMSSW Includes
//
#include "DataFormats/Math/interface/deltaPhi.h"

//
// User Defined Includes
//
#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "ElectroWeakAnalysis/VPlusJets/interface/QGLikelihoodCalculator.h"
#include "TAMUWW/MEPATNtuple/interface/PhysicsObject.hh"
#include "TAMUWW/MEPATNtuple/interface/Vertex.hh"
#include "TAMUWW/MEPATNtuple/interface/Jet.hh"
#include "TAMUWW/MEPATNtuple/interface/Lepton.hh"
#include "TAMUWW/MEPATNtuple/interface/MET.hh"

//
// ROOT includes
//
#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TRegexp.h"

//
// Standard Library Includes
//
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cassert>
#include <utility>

//
// Namespace
//
using std::vector;
using std::map;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::setw;
using std::left;

////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////
class GenParticle
{
public:
   //
   // Construction/Destruction
   //
   GenParticle();
   ~GenParticle();

   double charge;
   TLorentzVector p4;
   TVector3 vtx;
   int pdgId;
   int status;
   size_t particlePosition;
   int numberOfMothers;
   int numberOfDaughters;
   vector<size_t> motherPositions;
   vector<size_t> daughterPositions;

   ClassDef(GenParticle,1)
};

class EventNtuple : public TObject
{
public:
   //
   // Construction/Destruction
   //
   EventNtuple();
   ~EventNtuple();

   //
   // Grouped Cuts
   //
   bool baseCuts();
   bool FNALcutsElectron();
   bool FNALcutsMuon();

   //
   // Individual Cuts
   //
   bool METEtMin(double EtMin = 25.0);

   //
   // QGLikelihood
   //
   double getQGLikelihood(unsigned int index, QGLikelihoodCalculator* qglikeli = 0);
   vector < double > getQGLikelihoods(QGLikelihoodCalculator* qglikeli = 0);

   //
   // Jet Energy Resolution
   //
   void doJER();
   double getJERfactor(double pt, double eta, double ptgen);

   //
   // MET-Phi Corrections
   //
   void doMETPhiCorrection(TString eraType = "pfMEtSysShiftCorrParameters_2012runAvsNvtx_data");
   pair<double,double> getMETPhiCorrection(TString eraType = "pfMEtSysShiftCorrParameters_2012runAvsNvtx_data");

   //
   // Generator-Level Particle Information
   //
   // Print Event Information
   void printDecayInformation(int decayParticle, Int_t instance, Int_t depth, TString option = "");
   void printHiggsDecayInformation();
   // Check if a specific final state particle was contained in the event
   bool containsParticle(int pdgId, bool doAbs = true, double ptmin = 0.0);
   string particleNameFromInt(int p);
   // returns a number based on if the particle is a lepton or quark based on its pdgId
   enum particleType{LEPTON,NEUTRINO,QUARK,BOSON,UNKNOWN};
   int leptonNeutrinoOrQuark(int x);
   // returns the charge of a particle based on its pdgId
   double charge(int x);
   // returns a number rounded away from zero
   double round(double r);
   int round_int(double r);
   // returns the mass of the hadronic W and the leptonic W
   pair<double,double> onVsOffShell();
   pair<double,double> onVsOffShellInclusive(bool verbose = false);
   // form a string which contains the Lorentz vector information in the form (Pt, Eta, Phi, E)
   string getLVString(const TLorentzVector &lv);
   // return the TLorentzVector associated with the lepton from a W decay
   TLorentzVector getGenVorDaughter(particleType type, int expectedVPDGID=24, int *pdgid = 0, bool verbose = false);

   //
   // Trigger Information
   //
   // Get wasAccept information for a specific trigger
   bool triggered(char * triggerName, bool andor = true);
   // Get wasAccept information for a specific trigger
   bool findSpecificTrigger(string triggerName);
   // Print wasAccept information for a specific trigger
   void printSpecificTrigger(string triggerName);
   // Get wasAccept information for multiple triggers. Return the && if true and || if false.
   bool findTriggers(TString triggerName, bool andor = true);
   // Print wasAccept information for multiple triggers
   void printTriggers(TString triggerName);

   //
   // B-Tag Information
   //
   int getNBTags();

   //
   // Calculate some more complicated quantities
   //
   double getDeltaPhiJetJet();
   double getDeltaPhiMETJet();
   double getMinDeltaPhiMETJet();
   void getAngularVariables(Float_t &cosdPhiWW, Float_t &cosdPhiWH, Float_t &costhetal, 
                            Float_t &costhetaj, Float_t &costhetaWH, Float_t &jacksonAngle,
                            bool verbose = false);
   double getJacobePeak();
   double getDeltaRlepjj();
   double getMinDPhiLepJet();
   double getSumJetEt();

   //Needed for ME
   int run;
   int event;
   int lumi;
   double rho;
   vector < Vertex > vLV;
   vector < Jet > jLV;
   double Mjj;
   vector < Lepton > lLV;
   vector < MET > METLV;
   vector < GenParticle > genParticleCollection;

   map<string,bool> triggerMap;

   ClassDef(EventNtuple,11)
};

#endif
