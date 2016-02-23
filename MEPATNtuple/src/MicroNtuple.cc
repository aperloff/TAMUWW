#include "TAMUWW/MEPATNtuple/interface/MicroNtuple.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "TF1.h"

using std::cout;
using std::endl;
using std::make_pair;
using std::vector;
using std::map;

indexMap2 MicroNtuple::indexMap;

MicroNtuple::MicroNtuple()
{
  MicroNtuple(2);
}

MicroNtuple::MicroNtuple(int ijets) : nJets(ijets)
{
   //bProb = new Double_t[ijets];
   bProb = vector<Double_t>(ijets,0.0);
   clear();
}

MicroNtuple::MicroNtuple(const MicroNtuple& rhs)
{
  reader = rhs.reader;

  nJets = rhs.nJets;

  //delete [] bProb;
  //bProb = new Double_t[nJets]; //[ntag][nJets]
  bProb = vector<Double_t>(nJets,0.0);

  for (unsigned i = 0; i < nEventProb; ++i) {
    eventProb[i] = rhs.eventProb[i];
    eventMaxProb[i] = rhs.eventMaxProb[i];
  }
    
  
  for (int i = 0; i < nJets; ++i)
    bProb[i] = rhs.bProb[i];

  weight = rhs.weight;
  epd1tag = rhs.epd1tag;
  epd2tag = rhs.epd2tag;
  epd1tagSchan = rhs.epd1tagSchan;
  epd2tagSchan = rhs.epd2tagSchan;
  epd1tagTchan = rhs.epd1tagTchan;
  epd2tagTchan = rhs.epd2tagTchan;
  epd0tagWWandWZ = rhs.epd0tagWWandWZ;
  epd1tagWWandWZ = rhs.epd1tagWWandWZ;
  epd2tagWWandWZ = rhs.epd2tagWWandWZ;

  for (unsigned i=0; i<nWHmasses; i++) {
    epd1tagWH[i] = rhs.epd1tagWH[i];
    epd2tagWH[i] = rhs.epd2tagWH[i];
  }
 
  for (unsigned i=0; i<nHWWmasses; i++) {
    epd1tagHWW[i] = rhs.epd1tagHWW[i];
    epd2tagHWW[i] = rhs.epd2tagHWW[i];
  }

  for (unsigned i=0; i<nHiggsMasses; i++) {
    absorbError[i] = rhs.absorbError[i];
    epdPretagHiggs[i] = rhs.epdPretagHiggs[i];
    epd1tagHiggs[i] = rhs.epd1tagHiggs[i];
    epd2tagHiggs[i] = rhs.epd2tagHiggs[i];
  }

//  MjjMVA = rhs.MjjMVA;
  size = rhs.size;
  run = rhs.run;
  event = rhs.event;
}

MicroNtuple& MicroNtuple::operator=(const MicroNtuple& rhs)
{
   reader = rhs.reader;

   nJets = rhs.nJets;

   //delete [] bProb;
   bProb.clear();
   //bProb = new Double_t[nJets];
   bProb = vector<Double_t>(nJets,0.0);

   for (unsigned i = 0; i < nEventProb; ++i) {
     eventProb[i] = rhs.eventProb[i];
     eventMaxProb[i] = rhs.eventMaxProb[i];
   }

   for (int i = 0; i < nJets; ++i)
     bProb[i] = rhs.bProb[i];

   weight = rhs.weight;
   epd1tag = rhs.epd1tag;
   epd2tag = rhs.epd2tag;
   epd1tagSchan = rhs.epd1tagSchan;
   epd2tagSchan = rhs.epd2tagSchan;
   epd1tagTchan = rhs.epd1tagTchan;
   epd2tagTchan = rhs.epd2tagTchan;
   epd0tagWWandWZ = rhs.epd0tagWWandWZ;
   epd1tagWWandWZ = rhs.epd1tagWWandWZ;
   epd2tagWWandWZ = rhs.epd2tagWWandWZ;

   for (unsigned i=0; i<nWHmasses; i++) {
     epd1tagWH[i] = rhs.epd1tagWH[i];
     epd2tagWH[i] = rhs.epd2tagWH[i];
   }
 
   for (unsigned i=0; i<nHWWmasses; i++) {
     epd1tagHWW[i] = rhs.epd1tagHWW[i];
     epd2tagHWW[i] = rhs.epd2tagHWW[i];
   }

   for (unsigned i=0; i<nHiggsMasses; i++) {
     absorbError[i] = rhs.absorbError[i];
     epdPretagHiggs[i] = rhs.epdPretagHiggs[i];
     epd1tagHiggs[i] = rhs.epd1tagHiggs[i];
     epd2tagHiggs[i] = rhs.epd2tagHiggs[i];
   }

//   MjjMVA = rhs.MjjMVA;
   size = rhs.size;
   run = rhs.run;
   event = rhs.event;

   return *this;
}

MicroNtuple::~MicroNtuple()
{
  // I don't understand ROOT. Technically I should be able to 
  // delete bProb right here, however the destructor of MicroNtuple 
  // seems to be called by ROOT when using TTreeFormula even though
  // the constructor was never called. Another ROOT mistery...
  //if (bProb)
    //delete [] bProb;
}



void MicroNtuple::clear()
{
  for (unsigned i = 0; i < nEventProb; ++i) {
    eventProb[i] = 0;
    eventMaxProb[i] = 0;
  }
  

  //delete [] bProb;
  bProb.clear();
  //bProb = new Double_t[nJets];
  bProb = vector<Double_t>(nJets,0.0);
  //for (int i = 0; i < nJets; ++i)
    //bProb[i] = 0.;
  
  weight  = 0;
  epd1tag = 0;
  epd2tag = 0;
  epd1tagSchan = 0;
  epd2tagSchan = 0;
  epd1tagTchan = 0;
  epd2tagTchan = 0;
  epdPretagWWandWZ = 0;
  epd0tagWWandWZ = 0;
  epd1tagWWandWZ = 0;
  epd2tagWWandWZ = 0;

  for (unsigned i=0; i<nWHmasses; i++) {
    epd1tagWH[i] = 0.;
    epd2tagWH[i] = 0.;
  }
 
  for (unsigned i=0; i<nHWWmasses; i++) {
    epd1tagHWW[i] = 0.;
    epd2tagHWW[i] = 0.;
  }

  for (unsigned i=0; i<nHiggsMasses; i++) {
    absorbError[i] = 0.;
    epdPretagHiggs[i] = 0.;
    epd1tagHiggs[i] = 0.;
    epd2tagHiggs[i] = 0.;
  }

  reader = 0;
//  MjjMVA = 0;
  size = 0;
  run = 0;
  event = 0;

  KinBDT = -9999;
  MEBDT = -9999;
  KinMEBDT = -9999;
}

void MicroNtuple::printIndexMap(indexMap2 & im) {
  cout << "MicroNtuple::printIndexMap" << endl;
  for(indexMap2::iterator it=im.begin(); it!=im.end(); ++it) {
    for(indexMap1::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2) {
      cout << "tmeType = " << it->first << " tmeParam = " << it2->first << " index = " << it2->second << endl;
    }
  }
}

void MicroNtuple::checkIndexMap(indexMap2 & im) {
   bool printInfo = false;
   std::stringstream sstream;
   for(indexMap2::iterator it=im.begin(); it!=im.end(); ++it) {
      for(indexMap1::iterator it2=it->second.begin(); it2!=it->second.end(); ++it2) {
         printInfo = false;
         if(it2->first < 0) {
            sstream << "\tWARNING MicroNtuple::checkIndexMap a mass in the indexMap is less than zero" << endl;
            printInfo = true;
         }
         if(it2->second < 0) {
            sstream << "\tWARNING MicroNtuple::checkIndexMap an index in the indexMap is less than zero" << endl;
            printInfo = true;
         }
         if(printInfo)
            sstream << "\t\ttmeType = " << it->first << " tmeParam = " << it2->first << " index = " << it2->second << endl;
      }
   }

   if(!sstream.str().empty()) {
      cout << "MicroNtuple::checkIndexMap" << endl
           << sstream.str();
   }
}

// This method functions as a map which maps from the eventProbs to the ProbsForEPD.
// If you don't care about the higgs mass set it to zero.
ProbsForEPD  MicroNtuple::getEventProbs(double mhiggs, const double evtProb[nEventProb]) {


  ProbsForEPD res;

  //printIndexMap(indexMap);

  // The HWW which depends on the higgs mass
  if(indexMap.find(DEFS::EP::HWW)!=indexMap.end() && indexMap.find(DEFS::EP::HWW)->second.find(mhiggs)!=indexMap.find(DEFS::EP::HWW)->second.end())
     res.hww   = evtProb[indexMap[DEFS::EP::HWW][mhiggs]];
  else
     res.hww   = 0;

  // The WH which depends on the higgs mass
  if(indexMap.find(DEFS::EP::WH)!=indexMap.end() && indexMap.find(DEFS::EP::WH)->second.find(mhiggs)!=indexMap.find(DEFS::EP::WH)->second.end())
     res.wh    = evtProb[indexMap[DEFS::EP::WH][mhiggs]];
  else
     res.wh    = 0;

  // The rest which don't depend on the higgs mass
  if(indexMap.find(DEFS::EP::STopS)!=indexMap.end() && indexMap.find(DEFS::EP::STopS)->second.find(172.5)!=indexMap.find(DEFS::EP::STopS)->second.end())
     res.schan  = evtProb[indexMap[DEFS::EP::STopS][172.5]];
  else
     res.schan  = 0;
  if(indexMap.find(DEFS::EP::STopT)!=indexMap.end() && indexMap.find(DEFS::EP::STopT)->second.find(172.5)!=indexMap.find(DEFS::EP::STopT)->second.end())
     res.tchan  = evtProb[indexMap[DEFS::EP::STopT][172.5]];
  else
     res.tchan  = 0;
  if(indexMap.find(DEFS::EP::TTbar)!=indexMap.end() && indexMap.find(DEFS::EP::TTbar)->second.find(172.5)!=indexMap.find(DEFS::EP::TTbar)->second.end())
     res.tt     = evtProb[indexMap[DEFS::EP::TTbar][172.5]];
  else
     res.tt     = 0;
  if(indexMap.find(DEFS::EP::WLight)!=indexMap.end() && indexMap.find(DEFS::EP::WLight)->second.find(0)!=indexMap.find(DEFS::EP::WLight)->second.end())
     res.wlight = evtProb[indexMap[DEFS::EP::WLight][0]];
  else
     res.wlight = 0;
  if(indexMap.find(DEFS::EP::Wbb)!=indexMap.end() && indexMap.find(DEFS::EP::Wbb)->second.find(0)!=indexMap.find(DEFS::EP::Wbb)->second.end())
     res.wbb = evtProb[indexMap[DEFS::EP::Wbb][0]];
  else
     res.wbb = 0;
  if(indexMap.find(DEFS::EP::Wgg)!=indexMap.end() && indexMap.find(DEFS::EP::Wgg)->second.find(0)!=indexMap.find(DEFS::EP::Wgg)->second.end())
     res.wgg = evtProb[indexMap[DEFS::EP::Wgg][0]];
  else
     res.wgg = 0;
  if(indexMap.find(DEFS::EP::ZLight)!=indexMap.end() && indexMap.find(DEFS::EP::ZLight)->second.find(0)!=indexMap.find(DEFS::EP::ZLight)->second.end())
     res.zlight = evtProb[indexMap[DEFS::EP::ZLight][0]];
  else
     res.zlight = 0;
  if(indexMap.find(DEFS::EP::WW)!=indexMap.end() && indexMap.find(DEFS::EP::WW)->second.find(0)!=indexMap.find(DEFS::EP::WW)->second.end())
     res.ww     = evtProb[indexMap[DEFS::EP::WW][0]];
  else
     res.ww     = 0;
  if(indexMap.find(DEFS::EP::WZ)!=indexMap.end() && indexMap.find(DEFS::EP::WZ)->second.find(0)!=indexMap.find(DEFS::EP::WZ)->second.end())
     res.wz     = evtProb[indexMap[DEFS::EP::WZ][0]];
  else
     res.wz     = 0;
  if(indexMap.find(DEFS::EP::QCD)!=indexMap.end() && indexMap.find(DEFS::EP::QCD)->second.find(0)!=indexMap.find(DEFS::EP::QCD)->second.end())
     res.qcd    = evtProb[indexMap[DEFS::EP::QCD][0]];
  else
     res.qcd    = 0;

  //res.show("MicroNtuple::getEventProbs");

  return res;

}//getEventProbs


// This method just calls the overloaded getEventProbs(mhiggs, evtProbs[nEventProb])
// with the MicroNtuple's eventProb
ProbsForEPD MicroNtuple::getEventProbs(double mhiggs) const {
  return getEventProbs(mhiggs,eventProb);
}



//------------------------------------------------------------------------------
// This is the core of the code where all the parameters are mixed to obtain the EPD
ProbsForEPD MicroNtuple::getProbsForEPD(const ProbsForEPD & meProbs, 
					const ProbsForEPD & coeffs,
					const vector<double> bProb, DEFS::TagCat tagcat){


  // The returning probabilities
  ProbsForEPD probs;

  // Start with whatever we have for the ME's 
  probs = meProbs;
  //probs.show(std::string("Original meProbs"));
  //ProbsForEPD ncCoeffs = coeffs;
  //ncCoeffs.show(std::string("Coeffs"));
  // Multiply each by their correspondent parameters
  probs *= coeffs;
  //probs.show(std::string("final Probs"));
  // for the pretag category, in which we don't care about tags at all
  // just return as we have it 
  if (tagcat == DEFS::pretag) 
    return probs;

  // Otherwise for the tagged categories do this: 
  // If bProb[i] represents the probability that jet i is a b-jet 
  // Then for processes in which the number of actual b-quarks at 
  // tree level is:
  //   0 => multiply by (1 - bProb[0]) * (1 - bProb[1]);
  //   1 => multiply by bProb[0] * (1 - bProb[1]) + bProb[1] * (1 - bProb[0]);
  //   2 => multiply by bProb[0] * bProb[1];
  probs.wh     *=      bProb[0]  *  bProb[1];
  probs.hww    *= (1 - bProb[0]) * (1 - bProb[1]);
  probs.schan  *=      bProb[0]  *  bProb[1];
  probs.tchan  *=  0.5*(bProb[0] * (1 - bProb[1]) + bProb[1] * (1 - bProb[0]));
  probs.tchan2 *=  0.5*(bProb[0] * (1 - bProb[1]) + bProb[1] * (1 - bProb[0]));
  probs.tt     *=       bProb[0] * bProb[1];
  probs.wlight *= (1 - bProb[0]) * (1 - bProb[1]); 
  probs.zlight *= (1 - bProb[0]) * (1 - bProb[1]);
  probs.wbb    *=       bProb[0] * bProb[1];
  probs.wc     *= (1 - bProb[0]) * (1 - bProb[1]);
  probs.wgg    *= (1 - bProb[0]) * (1 - bProb[1]);
  probs.ww     *= (1 - bProb[0]) * (1 - bProb[1]);
  probs.wz     *= (1 - bProb[0]) * (1 - bProb[1]);
  probs.zz     *= (1 - bProb[0]) * (1 - bProb[1]);
  probs.qcd    *= (1 - bProb[0]) * (1 - bProb[1]);

  return probs;

}//getProbsForEPD



//------------------------------------------------------------------------------
//-----------------------------  THE WH METHODS      ---------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
Double_t MicroNtuple::calcWHEPD(DEFS::TagCat tagcat, double mhiggs) const
{

  // Get the meProbs
  ProbsForEPD meProbs = getEventProbs(mhiggs);
  return calcWHEPD(tagcat, mhiggs, meProbs, bProb, nJets);
  
}//calcWHEPD

//------------------------------------------------------------------------------
Double_t MicroNtuple::calcWHEPD(DEFS::TagCat tagcat, double mass,
				const ProbsForEPD & meProbs,
				const vector<double> bProb, Int_t nJets){

  // Get the Higgs Coefficients
  ProbsForEPD coeffs = getWHEPDCoefficients(mass, tagcat, nJets);  
  
  // Get the probabilities for the EPD
  ProbsForEPD probs = getProbsForEPD(meProbs, coeffs, bProb, tagcat);
  
  // There is less than three events per MC sample in which all the 
  // probabilities are zero.  To avoid returning a NaN just 
  // return zero for those
  if (probs.wh == 0) return 0;
  
  // Return the wh probability
  return probs.wh / (probs.wh + probs.schan + probs.tchan + probs.tchan2 + probs.wbb + 
		     probs.wc + probs.qcd + probs.tt +
		     probs.ww + probs.wz);      

}// calcWHEPD


//------------------------------------------------------------------------------
const ProbsForEPD MicroNtuple::getWHEPDCoefficients(double mhiggs, DEFS::TagCat tagcat, int nJets){

  // The returning object
  ProbsForEPD coeffs;

  if (tagcat == DEFS::pretag){
  
    if      (mhiggs <= 115) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 120) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 125) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 130) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else                    coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    
  }else if(tagcat == DEFS::eq0tag){

    if      (mhiggs <= 115) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 120) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 125) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 130) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else                    coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);

  }else if (tagcat ==  DEFS::eq1tag){

    if      (mhiggs <= 115) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 120) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 125) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 130) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else                    coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);

  } else if (tagcat == DEFS::eq2tag){

    if      (mhiggs <= 115) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 120) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 125) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 130) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else                    coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);

  } else
    throw std::runtime_error("ERROR MicroNtuple::getWHEPDCoefficients. Invalid tag requested");

  return coeffs;

}//getWHEPDCoefficients



//------------------------------------------------------------------------------
//-----------------------------  THE H->WW METHODS      ------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
Double_t MicroNtuple::calcHWWEPD(DEFS::TagCat tagcat, double mhiggs) const
{

  // Get the meProbs
  ProbsForEPD meProbs = getEventProbs(mhiggs);
  return calcWHEPD(tagcat, mhiggs, meProbs, bProb, nJets);
  
}//calcWHEPD

//------------------------------------------------------------------------------
Double_t MicroNtuple::calcHWWEPD(DEFS::TagCat tagcat, double mass,
				const ProbsForEPD & meProbs,
				const vector<double> bProb, Int_t nJets){

   
   // Translate mass to WH index
   //unsigned indexWH  = getHWWMassIndex(mass);
  
  cout<<" \t\t SFSG 1.1"<<endl;
  
   // Get the Higgs Coefficients
   ProbsForEPD coeffs = getHWWEPDCoefficients(mass, tagcat, nJets);  

   cout<<" \t\t SFSG 1.2"<<endl;

   // Get the probabilities for the EPD
   ProbsForEPD probs = getProbsForEPD(meProbs, coeffs, bProb, tagcat);

  cout<<" \t\t SFSG 1.3"<<endl;

   // There is less than three events per MC sample in which all the 
   // probabilities are zero.  To avoid returning a NaN just 
   // return zero for those
   if (probs.wh == 0) return 0;

   // Return the wh probability
   return probs.hww / (probs.hww + probs.schan + probs.tchan + probs.tchan2 + probs.wbb + 
		      probs.wc + probs.qcd + probs.tt +
		      probs.ww + probs.wz);      

}// calcHWWEPD

//------------------------------------------------------------------------------
const ProbsForEPD MicroNtuple::getHWWEPDCoefficients(double mhiggs, DEFS::TagCat tagcat, int nJets){

  // The returning object
  ProbsForEPD coeffs;

  if (tagcat == DEFS::pretag){
  
    if      (mhiggs <= 120) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 130) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 140) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 150) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else                    coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
 
  }else if (tagcat == DEFS::eq0tag){

    if      (mhiggs <= 120) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 130) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 140) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 150) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else                    coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);

  }else if (tagcat ==  DEFS::eq1tag){   

    if      (mhiggs <= 120) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 130) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 140) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 150) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else                    coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);

  } else if (tagcat == DEFS::eq2tag){

    if      (mhiggs <= 120) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 130) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 140) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 150) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else                    coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);

  } else
    throw std::runtime_error("ERROR MicroNtuple::getWHEPDCoefficients. Invalid tag requested");
    
  return coeffs;

}//getHWWEPDCoefficients



//------------------------------------------------------------------------------
//-----------------------------  THE Higgs METHODS      ---------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
Double_t MicroNtuple::calcHiggsEPD(DEFS::TagCat tagcat, double mhiggs) const
{

  // Get the meProbs
  ProbsForEPD meProbs = getEventProbs(mhiggs);
  return calcHiggsEPD(tagcat, mhiggs, meProbs, bProb, nJets);
  
}//calcHiggsEPD

//------------------------------------------------------------------------------
Double_t MicroNtuple::calcHiggsEPD(DEFS::TagCat tagcat, double mass,
				const ProbsForEPD & meProbs,
				const vector<double> bProb, Int_t nJets){

  // Get the Higgs Coefficients
  ProbsForEPD coeffs = getHiggsEPDCoefficients(mass, tagcat, nJets);  
  
  // Get the probabilities for the EPD
  ProbsForEPD probs = getProbsForEPD(meProbs, coeffs, bProb, tagcat);
  
  // There is less than three events per MC sample in which all the 
  // probabilities are zero.  To avoid returning a NaN just 
  // return zero for those
  if (probs.wh + probs.hww == 0) return 0;
  
  // Return the higgs probability
  //cout << "Mass Higgs " << mass << endl;
  //cout << "Prob wh " << probs.wh << endl;
  //cout << "Prob hww " << probs.hww << endl;
  //cout << "Prob wz " << probs.wz << endl;
  //cout << "Prob ww "<< probs.ww <<  endl;
  //cout << "Prob schan " << probs.schan << endl;
  //cout << "Prob tchan " << probs.tchan << endl;
  //cout << "Prob tchan2 " << probs.tchan2 << endl;
  //cout << "Prob wbb " << probs.wbb << endl;
  //cout << "Prob wc " << probs.wc << endl;
  //cout << "Prob qcd " << probs.qcd << endl;
  //cout << "Prob tt "<< probs.tt << endl;
  //cout << "EPD higgs " << (probs.wh + probs.hww) / (probs.wh + probs.hww + probs.schan + probs.tchan + probs.tchan2 + probs.wbb + 
  //                                                  probs.wc + probs.qcd + probs.tt +
  //                                                  probs.ww + probs.wz) << endl << endl; 
  return (probs.wh + probs.hww) / (probs.wh + probs.hww + probs.schan + probs.tchan + probs.tchan2 + probs.wbb + 
                                   probs.wc + probs.qcd + probs.tt +
                                   probs.ww + probs.wz);      

}// calcHiggsEPD


//------------------------------------------------------------------------------
const ProbsForEPD MicroNtuple::getHiggsEPDCoefficients(double mhiggs, DEFS::TagCat tagcat, int nJets){

  // The returning object
  ProbsForEPD coeffs;

  if (tagcat == DEFS::pretag){

    if      (mhiggs <= 115) coeffs = ProbsForEPD(1.75919e+09,2.95395e+09,3.28345e+07,121402,0,0,5158.93,147810,1.89051e+07,0,7882.16,1.33396e+06,6.3964e+06,0,0.408211);
    else if (mhiggs <= 120) coeffs = ProbsForEPD(1.75919e+09,2.95395e+09,3.28345e+07,121402,0,0,5158.93,147810,1.89051e+07,0,7882.16,1.33396e+06,6.3964e+06,0,0.408211);
    else if (mhiggs <= 125) coeffs = ProbsForEPD(1.75919e+09,2.95395e+09,3.28345e+07,121402,0,0,5158.93,147810,1.89051e+07,0,7882.16,1.33396e+06,6.3964e+06,0,0.408211);
    else if (mhiggs <= 130) coeffs = ProbsForEPD(1.75919e+09,2.95395e+09,3.28345e+07,121402,0,0,5158.93,147810,1.89051e+07,0,7882.16,1.33396e+06,6.3964e+06,0,0.408211);
    else                    coeffs = ProbsForEPD(1.75919e+09,2.95395e+09,3.28345e+07,121402,0,0,5158.93,147810,1.89051e+07,0,7882.16,1.33396e+06,6.3964e+06,0,0.408211);

    //if      (mhiggs <= 115) coeffs = ProbsForEPD(1.71676e+07,1.26194e+10,4.43136e+07,5.92018e+06,0,333.333,4000,458825,1.66778e+06,0,8000,1.25323e+07,3.70829e+07,0,1.41181);
    //else if (mhiggs <= 120) coeffs = ProbsForEPD(1.71676e+07,1.26194e+10,4.43136e+07,5.92018e+06,0,333.333,4000,458825,1.66778e+06,0,8000,1.25323e+07,3.70829e+07,0,1.41181);
    //else if (mhiggs <= 125) coeffs = ProbsForEPD(1.71676e+07,1.26194e+10,4.43136e+07,5.92018e+06,0,333.333,4000,458825,1.66778e+06,0,8000,1.25323e+07,3.70829e+07,0,1.41181);
    //else if (mhiggs <= 130) coeffs = ProbsForEPD(1.71676e+07,1.26194e+10,4.43136e+07,5.92018e+06,0,333.333,4000,458825,1.66778e+06,0,8000,1.25323e+07,3.70829e+07,0,1.41181);
    //else                    coeffs = ProbsForEPD(1.71676e+07,1.26194e+10,4.43136e+07,5.92018e+06,0,333.333,4000,458825,1.66778e+06,0,8000,1.25323e+07,3.70829e+07,0,1.41181);
    
  }else if (tagcat == DEFS::eq0tag){

    if      (mhiggs <= 115) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 120) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 125) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 130) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else                    coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);

  }else if (tagcat ==  DEFS::eq1tag){

    if      (mhiggs <= 115) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 120) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 125) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 130) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else                    coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);

  } else if (tagcat == DEFS::eq2tag){

    if      (mhiggs <= 115) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 120) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 125) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else if (mhiggs <= 130) coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
    else                    coeffs = ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);

  } else
    throw std::runtime_error("ERROR MicroNtuple::getHiggsEPDCoefficients. Invalid tag requested");

  return coeffs;

}//getHiggsEPDCoefficients



//------------------------------------------------------------------------------
//-----------------------------  THE WZ METHODS      ---------------------------
//------------------------------------------------------------------------------
Double_t MicroNtuple::calcWZEPD(DEFS::TagCat tagcat) const
{

  // Get the meProbs
  ProbsForEPD meProbs = getEventProbs(0);
  return calcWZEPD(tagcat, meProbs, bProb, nJets);
  
}

//------------------------------------------------------------------------------
// Without Higgs
Double_t MicroNtuple::calcWZEPD(DEFS::TagCat tagcat, const ProbsForEPD & meProbs,
				const vector<double> bProb, Int_t nJets){

  
  // Get the WZ Coefficients.
  ProbsForEPD coeffs = getWZEPDCoefficients(tagcat, nJets);  
  
  // Get the probabilities for the EPD
  ProbsForEPD probs = getProbsForEPD(meProbs, coeffs, bProb, tagcat);

  // Return the WZ probability
  //cout << "Prob wz " << probs.wz << endl;
  //cout << "Prob ww "<< probs.ww <<  endl;
  //cout << "Prob schan " << probs.schan << endl;
  //cout << "Prob tchan " << probs.tchan << endl;
  //cout << "Prob tchan2 " << probs.tchan2 << endl;
  //cout << "Prob wbb " << probs.wbb << endl;
  //cout << "Prob wc " << probs.wc << endl;
  //cout << "Prob qcd " << probs.qcd << endl;
  //cout << "Prob tt "<< probs.tt << endl;
  //cout << "EPD WZ " << (probs.wz+probs.ww) / (probs.schan + probs.tchan + probs.tchan2 + probs.wbb + 
  //                                            probs.wc + probs.qcd + probs.tt +
  //                                            probs.ww + probs.wz) << endl << endl;
  return (probs.wz+probs.ww) / (probs.schan + probs.tchan + probs.tchan2 + probs.wbb + 
				probs.wc + probs.qcd + probs.tt +
				probs.ww + probs.wz);    
        
}//calcWZEPD 

//------------------------------------------------------------------------------
const ProbsForEPD MicroNtuple::getWZEPDCoefficients(DEFS::TagCat tagcat, int nJets){

  // The returning object
  ProbsForEPD coeffs;

  if(tagcat == DEFS::pretag)
     return ProbsForEPD(4.79519e+09,1.50931e+09,1.13096e+08,2.00503e+07,0,0,5158.93,147810,81575.4,0,7882.16,3.73303e+07,1.68125e+08,0,0.00163767); // DEFS::TagCat=0 figOfMerit=2.05129
     //return ProbsForEPD(2.85714e+09,2e+10,5.37021e+07,1.3351e+07,0,333.333,4000,357143,1.6571e+07,0,8000,2.86161e+06,5.03532e+06,0,0.259127); // DEFS::TagCat=0 figOfMerit=11.0944
  else if (tagcat == DEFS::eq0tag)
    return ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
  else if (tagcat == DEFS::eq1tag)
    return ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);
  else if (tagcat != DEFS::eq2tag)
    throw std::runtime_error("ERROR MicroNtuple::getWZEPDCoefficients. Invalid tag requested");


  // the two tag case
  return ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0);

}//getWZEPDCoefficients



//------------------------------------------------------------------------------
//-----------------------------  THE SingleTop METHODS      --------------------
//------------------------------------------------------------------------------

Double_t MicroNtuple::calcSingleTopEPD(DEFS::TagCat tagcat, EPDType type) const {

  // Get the meProbs
  ProbsForEPD meProbs = getEventProbs(0);
  return calcSingleTopEPD(tagcat, type, meProbs, bProb, nJets);
  
} // calcSingleTopEPD



//------------------------------------------------------------------------------
const ProbsForEPD MicroNtuple::getSingleTopEPDCoefficients (DEFS::TagCat tagcat, EPDType type){
  
  static const ProbsForEPD params1tag[3] = {ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0),
					    ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0),
					    ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0)};
  
  static const ProbsForEPD params2tag[3] = {ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0),
					    ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0),
					    ProbsForEPD(2.5, 0, 0, .5, .01, .165, 0, 0, 0, 0,0,0,0,0,0)};

  if (tagcat == DEFS::eq1tag)
    return params1tag[type];
  else if (tagcat != DEFS::eq2tag)
    throw std::runtime_error("ERROR MicroNtuple::getSingleTopEPDCoefficients. Invalid tag requested");

  return  params2tag[type];

}//getSingleTopEPDCoefficients

//------------------------------------------------------------------------------
//EPD for SingleTop.
Double_t MicroNtuple::calcSingleTopEPD(DEFS::TagCat tagcat, EPDType type,
				       const ProbsForEPD & meProbs,
				       const vector<double> bProb, Int_t nJets){
 
  if (nJets == 2){

    // Get the Higgs Coefficients
    const ProbsForEPD coeffs = getSingleTopEPDCoefficients(tagcat,type);  
    
    // Get the probabilities for the EPD. the indexWH is set to zero and it's irrelevant here 
    ProbsForEPD probs = getProbsForEPD(meProbs, coeffs, bProb, tagcat);
    
    // There is less than three events per MC sample in which all the probabilities are zero. 
    // To avoid returning a NaN just return zero for those
    if (probs.schan + probs.tchan + probs.tchan2 == 0) return 0;

    // Return the wh probability. Interestingly there is not diboson here
    return (probs.schan+probs.tchan+probs.tchan2)/ (probs.schan +  probs.tchan + probs.tchan2 + probs.wbb + 
						    probs.wc + probs.qcd + probs.tt );       
   
   
  }//njets==2
   
  return 0;

}//calcEPD


//------------------------------------------------------------------------------
//-----------------------------  OTHER METHODS      ---------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
Double_t MicroNtuple::triggerTO(double detector, double etalep, double metraw,
                                double etlep)
{
  return 1;
}


//------------------------------------------------------------------------------
//-----------------------------  TMVA METHODS  ---------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
void MicroNtuple::setMVAReader(vector<TString> MVAMethods, TString dir,
                               map<TString,MVAVar>& MVAVars, bool debug) {
   // This loads the library
   TMVA::Tools::Instance();

   // --- Create the Reader object
   reader = new TMVA::Reader( "!Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   
   TString var;
   //Add non-higgs MEs
   size = 0;
   map<TString,MVAVar>::iterator firstUsed;
   map<TString,MVAVar>::iterator firstUsedSpectator;
   for(firstUsed = MVAVars.begin(); firstUsed!=MVAVars.end(); ++firstUsed) {
      if(firstUsed->second.use && !firstUsed->second.isSpectator) break;
   }
   for(firstUsedSpectator = MVAVars.begin(); firstUsedSpectator!=MVAVars.end(); ++firstUsedSpectator) {
      if(firstUsedSpectator->second.use && firstUsedSpectator->second.isSpectator) break;
   }

   while(size!=firstUsed->second.maxIndex) {
      for(map<TString,MVAVar>::iterator it=MVAVars.begin(); it!=MVAVars.end(); ++it) {
         if(it->second.use && !it->second.isSpectator && it->second.index == size) {
            if(debug) cout << "Using variable " << it->second.name << " in the TMVA::Reader" << endl; 
            reader->AddVariable(it->second.definition,&it->second.value);
            size++;
         }
      }
   }
   size = 0;
   while(size!=firstUsedSpectator->second.maxIndex) {
      for(map<TString,MVAVar>::iterator it=MVAVars.begin(); it!=MVAVars.end(); ++it) {
         if(it->second.use && it->second.isSpectator && it->second.index == size) {
            if(debug) cout << "Using spectator " << it->second.name << " in the TMVA::Reader" << endl;
            reader->AddSpectator(it->second.name,&it->second.value);
            size++;
         }
      }
   }

   // --- Book the MVA methods
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (unsigned int i=0; i<MVAMethods.size(); i++) {
      TString methodName = MVAMethods[i] + TString(" method");
      TString weightfile;
      if(dir.Contains(".xml"))
        weightfile = dir;
      else
        weightfile = dir + prefix + TString("_") + MVAMethods[i] + TString(".weights.xml");
      reader->BookMVA( methodName, weightfile );
   }
}

//------------------------------------------------------------------------------
vector<map<TString,Double_t> > MicroNtuple::getMVAOutput(vector<TString> MVAMethods, Double_t effS) {
   // This loads the library
   TMVA::Tools::Instance();
   
   // Retrieve MVA output, error, probability, and rarity.
   vector<map<TString,Double_t> > outputs;

   for (unsigned int i=0; i<MVAMethods.size(); i++) {
      map<TString,Double_t> output;
      TString methodName = MVAMethods[i] + TString(" method");  
      
      if (MVAMethods[i].CompareTo("CutsGA")==0)
         output["response"] = reader->EvaluateMVA( methodName, effS );
      else
         output["response"] = reader->EvaluateMVA( methodName );
      output["error"] = reader->GetMVAError();
      output["probability"] = reader->GetProba ( methodName );
      output["rarity"] = reader->GetRarity( methodName );

      outputs.push_back(output);
   }

   return outputs;
}

//------------------------------------------------------------------------------
vector<Double_t> MicroNtuple::getMVAOutput(std::vector<TString> MVAMethods, TString specificOutput, Double_t effS) {
   vector<Double_t> res;
   vector<map<TString,Double_t> > ret = getMVAOutput(MVAMethods, effS);
   for (unsigned int i=0; i<ret.size(); i++) {
      res.push_back(ret[i][specificOutput]);
   }
   return res;
}

//------------------------------------------------------------------------------
Double_t MicroNtuple::getMVAOutput(TString MVAMethod, TString specificOutput, Double_t effS) {
   vector<TString> mvam(1,MVAMethod);
   return getMVAOutput(mvam,specificOutput,effS).front();
}

//------------------------------------------------------------------------------
Double_t MicroNtuple::getCutsGAEfficiency(Int_t nSelCutsGA, Int_t nTreeEntries){
   return Double_t(nSelCutsGA)/nTreeEntries;
}

//------------------------------------------------------------------------------
void MicroNtuple::getMVACuts(std::vector<TString> &inputVars, std::vector<TString> &inputLabels, 
                             std::vector<TString> &inputTitles, std::vector<Double_t> &cutsMin, 
                             std::vector<Double_t> &cutsMax, Double_t effS) {
   // This loads the library
   TMVA::Tools::Instance();

   TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;
   
   if (mcuts) {      
      mcuts->GetCuts( effS, cutsMin, cutsMax );
      for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
         inputVars.push_back(mcuts->GetInputVar(ivar));
         inputLabels.push_back(mcuts->GetInputLabel(ivar));
         inputTitles.push_back(mcuts->GetInputTitle(ivar));
      }
   }
}

ClassImp(MicroNtuple)
