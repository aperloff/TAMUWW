/////////////////////////////////////////////////////////////////////////////////////////////////
//// CMS 
//// WW CrossSection Measurement using Matrix Element 
//// Created by Osipenkov, Ilya : ilyao@fnal.gov
/////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef EVENTNTUPLE_HH
#define EVENTNTUPLE_HH

#include "TLorentzVector.h"
#include "TObject.h"
//#include "TROOT.h"
//#include "TTree.h"
#include <string>
#include <vector>
//#include <math>
using std::vector;

///Headers for Transfer Function generation
#include "DataFormats/Math/interface/deltaR.h"
#include <utility>

class EventNtuple : public TObject
{
public:

  EventNtuple();
  ~EventNtuple();


  //Needed for ME
  //TLorentzVector lLV;
  vector < TLorentzVector > jLV, METLV, lLV;
  vector < TLorentzVector > matchedGenParticles;
  vector < int > matchedpdgId;
  vector < double > matchedDeltaR;
  vector < int > jBtag;
  int lQ;
  int ldetComp;
  int run;
  int event;

  //Additional
  double        Mjj;
  int           passStd;
  int           passAll;

  double        DRlj1;
  double        DRlj2;
  double        Thetalj1pj2;

  double        lTotIso;
  double        lecalIso;
  double        lhcalIso;
  double        ltrkIso;

  double        METEt;
  double        lPhi;
  ClassDef(EventNtuple, 2)

};

#endif
