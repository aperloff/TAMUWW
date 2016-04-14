#ifndef STYLE_HH
#define STYLE_HH

#include "TStyle.h"
#include "TROOT.h"
#include "TPad.h"
#include "TFrame.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "THStack.h"
#include "TPaveText.h"
#include "TString.h"

#include <iostream>
#include <string>
#include <assert.h>

using std::cout;
using std::endl;
using std::string;

class Style {
public:
  Style() {reset_globals();}
  ~Style() {;}
 
  void tdrDraw(TH1* h, string opt,
               int marker=kFullCircle, int mcolor = kBlack,
               int lstyle=kSolid, int lcolor=-1,
               int fstyle=1001, int fcolor=kYellow+1);
  void tdrDraw(TProfile* h, string opt,
               int marker=kFullCircle, int mcolor = kBlack,
               int lstyle=kSolid, int lcolor=-1,
               int fstyle=1001, int fcolor=kYellow+1);
  void tdrDraw(THStack* s, string opt);
  void tdrDraw(TGraph* g, string opt,
               int marker=kFullCircle, int mcolor = kBlack,
               int lstyle=kSolid, int lcolor=-1,
               int fstyle=1001, int fcolor=kYellow+1);
  TLegend *tdrLeg(double x1, double y1, double x2, double y2);
  TPaveText *tdrText(double x1, double y1, double x2, double y2, int alignment = 21);

  void tdrGrid(bool gridOn);

  void fixOverlay();

  void setTDRStyle();

  void CMS_lumi( TPad* pad, int iPeriod=3, int iPosX=10 );
  void reset_globals();

  TCanvas* tdrCanvas(const char* canvName, TH1D *h,
                            int iPeriod = 2, int iPos = 11,
                            bool square = kRectangular);
  TCanvas* tdrDiCanvas(const char* canvName, TH1D *hup, TH1D *hdw,
                              int iPeriod = 2, int iPos = 11, double lumiOverride = -1);

  void cmsPrel(int energy = 8, double intLumi=-1, bool wide = false);

  void cmsFinal(double intLumi=-1, bool wide = false);

private:
  static const bool kSquare = true;
  static const bool kRectangular = false;

  // colors to use
  //EColor tdrColors[13] = {kBlack, kBlue, kRed, kGreen, kYellow, kMagenta, kCyan,
  //                        kOrange, kSpring, kTeal, kAzure, kViolet, kPink};

  TString cmsText;
  float cmsTextFont;  // default is helvetic-bold
  
  bool writeExtraText;//false;
  TString extraText;
  TString extraText2; // For Simulation Preliminary on two lines
  float extraTextFont;  // default is helvetica-italics
  
  // text sizes and text offsets with respect to the top frame
  // in unit of the top margin size
  float lumiTextSize;
  float lumiTextOffset;
  float cmsTextSize;
  float cmsTextOffset;  // only used in outOfFrame version
  
  float relPosX;
  float relPosY;
  float relExtraDY;

  // ratio of "CMS" and extra text size
  float extraOverCmsTextSize;
  
  TString lumi_13TeV;
  TString lumi_8TeV;
  TString lumi_7TeV;
  
  bool drawLogo;
};

#endif
