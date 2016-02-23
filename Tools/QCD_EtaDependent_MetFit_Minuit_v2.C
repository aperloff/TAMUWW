#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TFractionFitter.h"
#include "TFitResult.h"
#include "TMinuit.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TGraphErrors.h"

#include <iostream>
#include <vector>


using namespace std;

double veta[] = { -5.191, -4.889, -4.716, -4.538, -4.363, -4.191, -4.013, -3.839, -3.664
		  ,-3.489, -3.314, -3.139, -2.964, -2.853, -2.650, -2.500, -2.322, -2.172 
		  ,-2.043, -1.930, -1.830, -1.740, -1.653, -1.566, -1.479, -1.392, -1.305 
		  ,-1.218, -1.131, -1.044, -0.957, -0.879, -0.783, -0.696, -0.609, -0.522 
		  ,-0.435, -0.348, -0.261, -0.174, -0.087 
		  ,+0.000 
		  ,+0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783
		  ,+0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566 
		  ,+1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500, +2.650 
		  ,+2.853, +2.964, +3.139, +3.314, +3.489, +3.664, +3.839, +4.013, +4.191 
		  ,+4.363, +4.538, +4.716, +4.889, +5.191};

double vabseta[] = { +0.000 
		     ,+0.087, +0.174, +0.261, +0.348, +0.435, +0.522, +0.609, +0.696, +0.783
		     ,+0.879, +0.957, +1.044, +1.131, +1.218, +1.305, +1.392, +1.479, +1.566 
		     ,+1.653, +1.740, +1.830, +1.930, +2.043, +2.172, +2.322, +2.500, +2.650 
		     ,+2.853, +2.964, +3.139, +3.314, +3.489, +3.664, +3.839, +4.013, +4.191 
		     ,+4.363, +4.538, +4.716, +4.889, +5.191};


//global so i can use it in the minimization function
TH1D * met_s, * met_q, * met_w, * met_t, * met_r;
double ndf;

// ----------------------------------------------------------------------------
// The minimization function 
void minimizationFunction(Int_t &npar, Double_t *gin, Double_t &f,
			  Double_t *par, Int_t iflag){


  double met_sint = met_s->Integral();
  double met_qint = met_q->Integral();
  double met_wint = met_w->Integral();
  double met_tint = met_t->Integral();


  //calculate chisquare
  Double_t chi2 = 0;
  ndf = 0;
  for (int ibin = 0 ; ibin <= met_s->GetNbinsX()+1 ; ibin++) {

    //skip events without entries
    if (met_s->GetBinContent(ibin) == 0  ||
	met_q->GetBinContent(ibin) == 0  ||
	met_w->GetBinContent(ibin) == 0  ||
    met_t->GetBinContent(ibin) == 0)
      continue;

    double bin_num =  met_s->GetBinContent(ibin) 
      - par[0] * met_q->GetBinContent(ibin) * met_sint / met_qint 
      - par[1] * met_w->GetBinContent(ibin) * met_sint / met_wint 
      - par[2] * met_t->GetBinContent(ibin) * met_sint / met_tint
      - met_r->GetBinContent(ibin);
  
    // In the denominator use the error of the bin, which is Poisson scaled
    double bin_den2 = 0 
      //+  pow(met_s->GetBinError(ibin),2) 
      +  pow(par[0] * met_q->GetBinError(ibin) * met_sint / met_qint,2)
      +  pow(par[1] * met_w->GetBinError(ibin) * met_sint / met_wint,2);
      +  pow(par[2] * met_t->GetBinError(ibin) * met_sint / met_tint,2);
    if (met_r->GetBinContent(ibin)>0)
      bin_den2 += pow(met_r->GetBinError(ibin),2);   

    if (bin_den2>0){
      chi2 += pow(bin_num,2)/bin_den2;
      ndf ++;
    }

  }//for

  ndf = ndf-2;
  f = chi2;

}// minimizationFunction



void getHistosFromDrawOutput(TH2D *& hd, TH2D *& hw, TH2D *& hq,  int rebinEta, int rebinMet){

  // Data: open the input file and get all histos
  TFile *_file0 = TFile::Open("output_signal.root");
  hd =  (TH2D*) _file0->Get("h2");
  hd->Rebin2D(rebinEta, rebinMet);

  // WJets: open the input file and get all histos
  TFile *_file1 = TFile::Open("output_wjets.root");
  hw =  (TH2D*) _file1->Get("h2");
  hw->Rebin2D(rebinEta, rebinMet);

  // QCD: open the input file and get all histos
  TFile *_file2 = TFile::Open("output_nomvacutnew.root");
  //TFile *_file2 = TFile::Open("output_mvacut.root");
  hq =  (TH2D*) _file2->Get("h2");
  hq->Rebin2D(rebinEta, rebinMet);

}//getHistosFromDrawOutput

void getHistosFromPlotterOutput(TH2D *& hd, TH2D *& hw, TH2D *& hq, TH2D *& ht, TH2D  *& hr, int rebinEta, int rebinMet, bool full){

  bool removesw2 = false;// true;

  // list of background processes. WJets and QCD are handled separately
  TString hn = "MET_vs_AbsLeptEta_";

  // Name of the data histo
  TString data = hn+"SingleEl_Data_electron"; // 19148 pb-1

  // Name of the WJets histo
  TString wj = hn+"WJets_electron";

  // Name of the QCD histo
  TString qcd;
  if(full)
     qcd = hn+"QCD_ElFULL_electron"; // 19148 pb-1
  else
     qcd = hn+"QCD_ElEnriched_electron"; // 1169 pb-1

  // Name of the TTbar histo
  TString ttbar = hn+"TTbar_electron";

  // Vector containing all other histos. Proper normalization 
  // to the luminosity of the signal data is assumed
  vector<TString> backs;
  backs.push_back(hn+"STopS_T_electron");
  backs.push_back(hn+"STopS_Tbar_electron");
  backs.push_back(hn+"STopT_T_electron");
  backs.push_back(hn+"STopT_Tbar_electron");
  backs.push_back(hn+"STopTW_T_electron");
  backs.push_back(hn+"STopTW_Tbar_electron");
  backs.push_back(hn+"ZJets_electron");
  backs.push_back(hn+"WW_electron");
  backs.push_back(hn+"WZ_electron");
  backs.push_back(hn+"ggH125_electron");
  backs.push_back(hn+"qqH125_electron");
  backs.push_back(hn+"WH125_electron");
   
  // Open the input file and get all histos
  //  TFile *_file0 = TFile::Open("/uscms_data/d2/eusebi/validationForQCD/histos_electron.root");
  TFile *_file0 = TFile::Open("./histos_electron.root");

  // Data
  hd =  (TH2D*) _file0->Get(data);
  if (removesw2) hd->GetSumw2()->Set(0); 
  hd->Rebin2D(rebinEta, rebinMet);
 
  //WJets
  hw = (TH2D*) _file0->Get(wj);
  if (removesw2) hw->GetSumw2()->Set(0); 
  hw->Rebin2D(rebinEta, rebinMet);

  //QCD 
  hq = (TH2D*) _file0->Get(qcd);
  if (removesw2) hq->GetSumw2()->Set(0); 
  hq->Rebin2D(rebinEta, rebinMet);

  //TTbar
  ht = (TH2D*) _file0->Get(ttbar);
  if (removesw2) ht->GetSumw2()->Set(0);
  ht->Rebin2D(rebinEta, rebinMet);

  // Remainding Processes
  // Substract all other backgrounds from the data
  hr = 0;
  for (unsigned int b = 0 ; b < backs.size();b++){
    TH2D * aux = (TH2D*) _file0->Get(backs[b]);
    aux->Rebin2D(rebinEta,rebinMet);
    if (hr==0) hr = aux;
    else hr->Add(aux);
  }//for

}// getHistosFromPlotterOutput


void QCD_EtaDependent_MetFit_Minuit(bool full = true){

  // ============ GET THE INPUT HISTOS ==============
  int rebinMet = 20; // best value(s): 5 or 20
  int rebinEta = 2; // best value(s): 2 or 41
  bool fixQCDFraction = false;

  double lum_data = 19148; // units of pb-1
  double lum_qcd;
  if(full)
     lum_qcd = 18937.6012264;// units of pb-1
     //lum_qcd = 19148;// units of pb-1
  else
     lum_qcd = 1143.95; // units of pb-1
     //lum_qcd = 1169; // units of pb-1

  // ============ CREATE THE OUTPUT HISTOS ==============
  TFile *fout1 = new TFile("QCD_EtaDependent_MetFit_output.root","RECREATE");

  // Data, WJets, QCD
  TH2D * hd, *hw, *hq, *ht, *hr;
  getHistosFromPlotterOutput(hd,hw,hq,ht,hr,rebinEta,rebinMet,full);
  //getHistosFromDrawOutput(hd,hw,hq, rebinEta,rebinMet);

  fout1->cd();

  // set the styles
  hd->UseCurrentStyle();
  hw->UseCurrentStyle();
  hq->UseCurrentStyle();
  ht->UseCurrentStyle();
  hr->UseCurrentStyle();

  // The integral in data
  TH1D * da_int = hd->ProjectionX("da_int");
  da_int->SetTitle("da_int");
  da_int->Reset();

  // The integral in qcd 
  TH1D * qcd_int = hd->ProjectionX("qcd_int");
  qcd_int->SetTitle("qcd_int");
  qcd_int->Reset();

  // The integral in wj
  TH1D * wj_int = hd->ProjectionX("wj_int");
  wj_int->SetTitle("wj_int");
  wj_int->Reset();

  // The integral in ttbar
  TH1D * ttbar_int = hd->ProjectionX("ttbar_int");
  ttbar_int->SetTitle("ttbar_int");
  ttbar_int->Reset();  

  // The number of entries in data
  TH1D * da_nev = hd->ProjectionX("da_nev");
  da_nev->SetTitle("da_nev");
  da_nev->Reset();

  // The number of entries in qcd
  TH1D * qcd_nev = hd->ProjectionX("qcd_nev");
  qcd_nev->SetTitle("qcd_nev");
  qcd_nev->Reset();

  // The number of entries in wj
  TH1D * wj_nev = hd->ProjectionX("wj_nev");
  wj_nev->SetTitle("wj_nev");
  wj_nev->Reset();

  // The number of entries in ttbar
  TH1D * ttbar_nev = hd->ProjectionX("ttbar_nev");
  ttbar_nev->SetTitle("ttbar_nev");
  ttbar_nev->Reset();

  // The relative qcd fraction according to the fit 
  TH1D * qcd_fr = hd->ProjectionX("qcd_fr");
  qcd_fr->SetTitle("qcd_fr");
  qcd_fr->Reset();

  // The relative wjets fraction according to the fit 
  TH1D * wj_fr = hd->ProjectionX("wj_fr");
  wj_fr->SetTitle("wj_fr");
  wj_fr->Reset();

  // The relative ttbar fraction according to the fit 
  TH1D * ttbar_fr = hd->ProjectionX("ttbar_fr");
  ttbar_fr->SetTitle("ttbar_fr");
  ttbar_fr->Reset();

  // The chi2 of the fit 
  TH1D * chi2 = hd->ProjectionX("chi2");
  chi2->SetTitle("chi2");
  chi2->Reset();

  // The qcd yield according to the fit 
  TH1D * qcd_yield = hd->ProjectionX("qcd_yield");
  qcd_yield->SetTitle("qcd_yield;#eta;number of events / 19148 fb^{-1}");
  qcd_yield->Reset();

  // The wjets yield according to the fit 
  TH1D * wj_yield = hd->ProjectionX("wj_yield");
  wj_yield->SetTitle("wj_yield;#eta;number of events / 19148 fb^{-1}");
  wj_yield->Reset();

  // The wjets yield according to the fit 
  TH1D * ttbar_yield = hd->ProjectionX("ttbar_yield");
  ttbar_yield->SetTitle("ttbar_yield;#eta;number of events / 19148 fb^{-1}");
  ttbar_yield->Reset();

  // The qcd production xs according to the fit 
  TH1D * qcd_xs = hd->ProjectionX("qcd_xs");
  qcd_xs->SetTitle("qcd_xs");
  qcd_xs->Reset();

  // The wjets production xs according to the fit 
  TH1D * wj_xs = hd->ProjectionX("wj_xs");
  wj_xs->SetTitle("wj_xs");
  wj_xs->Reset();

  // The ttbar production xs according to the fit 
  TH1D * ttbar_xs = hd->ProjectionX("ttbar_xs");
  ttbar_xs->SetTitle("ttbar_xs");
  ttbar_xs->Reset();

  // The qcd production scale factor according to the fit 
  TH1D * qcd_sf = hd->ProjectionX("qcd_sf");
  qcd_sf->SetTitle(";|#eta|;s_{QCD}");
  qcd_sf->Reset();

  // The wjets production scale factor according to the fit 
  TH1D * wj_sf = hd->ProjectionX("wj_sf");
  wj_sf->SetTitle(";|#eta|;N^{WJets}_{measured} /N^{WJets}_{expected} ");
  wj_sf->Reset();

  // The wjets production scale factor according to the fit 
  TH1D * ttbar_sf = hd->ProjectionX("ttbar_sf");
  ttbar_sf->SetTitle(";|#eta|;N^{TTbar}_{measured} /N^{TTbar}_{expected} ");
  ttbar_sf->Reset();

  TGraphErrors *wjsm_sf = new TGraphErrors();
  wjsm_sf->SetPoint(0,0.05,1.0);
  wjsm_sf->SetPointError(0,0,0.0256);
  wjsm_sf->SetPoint(1,5.0,1.0);
  wjsm_sf->SetPointError(1,0,0.0256);

  wjsm_sf->SetFillColor(kGreen);
  wjsm_sf->SetMarkerColor(0);


  TGraphErrors *ttbarsm_sf = new TGraphErrors();
  ttbarsm_sf->SetPoint(0,0.05,1.0);
  ttbarsm_sf->SetPointError(0,0,0.0256);
  ttbarsm_sf->SetPoint(1,5.0,1.0);
  ttbarsm_sf->SetPointError(1,0,0.0256);
  ttbarsm_sf->SetFillColor(kGreen);
  ttbarsm_sf->SetMarkerColor(0);
  

  // The qcd production xs according to the fit 
  TH1D * qcd_close = hd->ProjectionX("qcd_close");
  qcd_close->SetTitle("qcd_close;#eta");
  qcd_close->Reset();

  // The wjets production xs according to the fit 
  TH1D * wj_close = hd->ProjectionX("wj_close");
  wj_close->SetTitle("wj_close;#eta");
  wj_close->Reset();

  // The wjets production xs according to the fit 
  TH1D * ttbar_close = hd->ProjectionX("ttbar_close");
  ttbar_close->SetTitle("ttbar_close;#eta");
  ttbar_close->Reset();


  // ============ AUXILIARY VARIABLES ==============

  // Show stack every <skip_stack> fits in range |eta|<3
  int skip_stack = 1;
  int min_ieta = qcd_fr->FindBin(-3);
  int max_ieta = qcd_fr->FindBin(3);

  // the max number of canvases to be drawn is then 
  int ncan = max_ieta-min_ieta+1;

  //Create the Canvas and divide it
  TCanvas * cfits = new TCanvas("fits","fits",1200,1200);
  int ncols = ceil(sqrt(ncan));
  int nrows = ceil(ncan/ncols);
  cfits->Divide(ncols,nrows);
  int icfits=1; // counter for canvas pad


  // ============ LOOP OVER ETAS AND DO THE MET FITS ==============
  for (int ieta = min_ieta; ieta <= max_ieta ; ieta++){
      
    double eta = hd->GetXaxis()->GetBinLowEdge(ieta);
    double etap1 = hd->GetXaxis()->GetBinLowEdge(ieta+1);

    cout<<"******* Fitting ieta="<<ieta<<": "<<eta<<" << eta <<"<< etap1 << endl;

    // get the met distributions for that eta
    met_s = hd->ProjectionY(Form("s_ieta%i",ieta),ieta,ieta);
    met_q = hq->ProjectionY(Form("q_ieta%i",ieta),ieta,ieta);
    met_w = hw->ProjectionY(Form("w_ieta%i",ieta),ieta,ieta);
    met_t = ht->ProjectionY(Form("t_ieta%i",ieta),ieta,ieta);
    met_r = hr->ProjectionY(Form("r_ieta%i",ieta),ieta,ieta);

    if ( met_s->GetEntries()==0 || met_q->GetEntries()==0 || met_t->GetEntries()==0 || met_w->GetEntries()==0 ){
      cout<<"\t skip this eta due to empty entries"<<endl; 
      continue;
    }


    // ============ SETUP MINUIT
    TMinuit minuit(3);
  
    // Set the minimization function
    minuit.SetFCN(minimizationFunction);
  
    Double_t arglist[10];
    Int_t ierflg = 0;
  
    // Make Minuit run in a quiet (i.e. non-verbose) mode
    minuit.SetPrintLevel(-1);
    arglist[0] = 0;
    minuit.mnexcm("SET NOWarnings",&arglist[0],1,ierflg);
    arglist[0] = -1;
    minuit.mnexcm("SET PRint",&arglist[0],1,ierflg);

  
    // Set the error as half a point around the NLL
    arglist[0] = 1;
    minuit.mnexcm("SET ERR", arglist ,1,ierflg);
  
    // Set starting values and step sizes for parameters
    // (par# , name, start val , step , min val, max val, ierflg);
    if (fixQCDFraction) {
      double qcdfix = met_q->Integral()/met_s->Integral();
      minuit.mnparm( 0, "qcd_sf", qcdfix, 0.01, 0, 1, ierflg); // start at 3%
      minuit.FixParameter(0);
    }else
      minuit.mnparm( 0, "qcd_sf", 0.10 , 0.01 , 0 , 1 , ierflg); // start at 3%
    
    minuit.mnparm( 1, "wj_sf", 0.50 , 0.01 , 0 , 1 , ierflg); // start at 60%
    minuit.mnparm( 2, "ttbar_sf", 0.40 , 0.01 , 0 , 1 , ierflg); // start at 40%
  
    // Do the minimization
    arglist[0] = 1; //1 for 1 sigma.
    //minuit.mnexcm("SET NOG",arglist,0,ierflg); //
    arglist[0] = 5000;
    arglist[1] = 0.01;
    minuit.mnexcm("MIGRAD", arglist ,2,ierflg); // -
    minuit.mnexcm("MINOS",arglist,2,ierflg);
  
    Int_t status =0;
    if (status==0) {                       
      double val0=0, err0=0;
      double val1=0, err1=0;
      double val2=0, err2=0;
    
      // get the results of the fit here
    
      // Retrieve the central value (central) and errors (up,down)
      double central, err_down, err_up, junk, junk1, junk2;
      minuit.GetParameter(0, val0, junk);
      minuit.mnerrs(0, err_up, err_down, junk1, junk2);
      err0 =  0.5*(fabs(err_up)+fabs(err_down));

      minuit.GetParameter(1, val1, junk);
      minuit.mnerrs(0, err_up, err_down, junk1, junk2);
      err1 =  0.5*(fabs(err_up)+fabs(err_down));

      minuit.GetParameter(2, val2, junk);
      minuit.mnerrs(0, err_up, err_down, junk1, junk2);
      err2 =  0.5*(fabs(err_up)+fabs(err_down));         

      cout<<"\t qcd_sf ="<<val0<<" +/- "<<err0<<endl;
      cout<<"\t wj_sf  ="<<val1<<" +/- "<<err1<<endl;
      cout<<"\t ttbar_sf  ="<<val2<<" +/- "<<err2<<endl;
      double met_sint = met_s->Integral();
      double met_qint = met_q->Integral();
      double met_wint = met_w->Integral();
      double met_tint = met_t->Integral();

      Double_t amin,edm,errdef;
      Int_t nvpar,nparx,icstat;
      minuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
 

      double redchi2 = amin/ndf; 
      chi2->SetBinContent(ieta,redchi2);
      chi2->SetBinError(ieta,redchi2/sqrt(ndf));   
      cout<<"\t chi2   = "<<amin<<" ndf = "<<ndf<<endl;

      // fill QCD info
      qcd_int->SetBinContent(ieta, met_q->Integral());
      qcd_nev->SetBinContent(ieta, met_q->GetEntries());
      qcd_fr->SetBinContent(ieta,val0);
      qcd_fr->SetBinError(ieta,err0);

      qcd_yield->SetBinContent(ieta,val0*met_s->Integral());
      double rel2 = pow(err0/val0,2) + 1.0/met_s->Integral();
      qcd_yield->SetBinError(ieta,sqrt(rel2) * qcd_yield->GetBinContent(ieta));

      qcd_xs->SetBinContent(ieta, val0*met_s->Integral()/lum_data);
      qcd_xs->SetBinError(ieta, err0*met_s->Integral()/lum_data);

      qcd_sf->SetBinContent(ieta, val0*met_s->Integral()*lum_qcd / (lum_data*met_q->GetEntries()));
      rel2 = pow(err0/val0,2) 
	+ 1./met_s->GetEntries()
	+ 1./met_q->GetEntries();
      qcd_sf->SetBinError  (ieta, sqrt(rel2) * qcd_sf->GetBinContent(ieta));

      qcd_close->SetBinContent (ieta, val0*met_s->Integral()*lum_qcd / (lum_data*met_q->Integral()));
      qcd_close->SetBinError(ieta, sqrt(rel2)*qcd_close->GetBinContent(ieta));

      //double rel2 = pow(err0/val0,2) + pow(sqrt(met_q->Integral())/met_q->Integral(),2);
      //qcd_xs->SetBinError(ieta, sqrt(rel2) * qcd_xs->GetBinContent(ieta));
      met_q->Scale(val0*met_s->Integral()/met_q->Integral()); // scale to plot on stack

      // fill WJets info
      wj_int->SetBinContent(ieta, met_w->Integral());
      wj_nev->SetBinContent(ieta, met_w->GetEntries());
      wj_fr->SetBinContent(ieta,val1);
      wj_fr->SetBinError(ieta,err1);
      wj_yield->SetBinContent(ieta,val1*met_s->Integral());
      wj_yield->SetBinError(ieta,err1*met_s->Integral());
      wj_xs->SetBinContent(ieta, val1*met_s->Integral()/lum_data);
      wj_xs->SetBinError(ieta, err1*met_s->Integral()/lum_data);
      wj_sf->SetBinContent(ieta, val1*met_s->Integral() /met_w->Integral());
      rel2 = pow(err1/val1,2) 
	+ 1./met_s->GetEntries()
	+ 1./met_w->GetEntries()
	+ pow(960./37509.,2); // added error on expected WJets XS
      wj_sf->SetBinError(ieta, sqrt(rel2)*wj_sf->GetBinContent(ieta));
      wj_close->SetBinContent(ieta, val1*met_s->Integral() /met_w->Integral());
      rel2 = pow(err1/val1,2)
	+ 1.0/met_s->GetEntries()
	+ 1.0/met_w->GetEntries();
      wj_close->SetBinError(ieta, sqrt(rel2)*wj_close->GetBinContent(ieta));
      cout<<"\t ieta="<<ieta<<" wj closure ="<< wj_close->GetBinContent(ieta)
	  <<" +/- "<<wj_close->GetBinError(ieta)<<endl;
      
      // rel2 = pow(err1/val1,2) + pow(sqrt(met_w->Integral())/met_w->Integral(),2);
      // wj_xs->SetBinError(ieta, sqrt(rel2) * wj_xs->GetBinContent(ieta));
      met_w->Scale(val1*met_s->Integral()/met_w->Integral()); // scale to plot on stack

      // fill TTbar info
      ttbar_int->SetBinContent(ieta, met_t->Integral());
      ttbar_nev->SetBinContent(ieta, met_t->GetEntries());
      ttbar_fr->SetBinContent(ieta,val2);
      ttbar_fr->SetBinError(ieta,err2);
      ttbar_yield->SetBinContent(ieta,val2*met_s->Integral());
      ttbar_yield->SetBinError(ieta,err2*met_s->Integral());
      ttbar_xs->SetBinContent(ieta, val2*met_s->Integral()/lum_data);
      ttbar_xs->SetBinError(ieta, err2*met_s->Integral()/lum_data);
      ttbar_sf->SetBinContent(ieta, val2*met_s->Integral() /met_t->Integral());
      rel2 = pow(err2/val2,2) 
	+ 1./met_s->GetEntries()
	+ 1./met_t->GetEntries()
	+ pow(14./225.197,2); // added error on expected TTbar XS
      ttbar_sf->SetBinError(ieta, sqrt(rel2)*ttbar_sf->GetBinContent(ieta));
      ttbar_close->SetBinContent(ieta, val2*met_s->Integral() /met_t->Integral());
      rel2 = pow(err2/val2,2)
	+ 1.0/met_s->GetEntries()
	+ 1.0/met_t->GetEntries();
      ttbar_close->SetBinError(ieta, sqrt(rel2)*ttbar_close->GetBinContent(ieta));
      cout<<"\t ieta="<<ieta<<" ttbar closure ="<< ttbar_close->GetBinContent(ieta)
	  <<" +/- "<<ttbar_close->GetBinError(ieta)<<endl;
      
      // rel2 = pow(err1/val1,2) + pow(sqrt(met_w->Integral())/met_w->Integral(),2);
      // wj_xs->SetBinError(ieta, sqrt(rel2) * wj_xs->GetBinContent(ieta));
      met_t->Scale(val2*met_s->Integral()/met_t->Integral()); // scale to plot on stack

      // fill Data info
      da_int->SetBinContent(ieta, met_s->Integral());
      da_nev->SetBinContent(ieta, met_s->GetEntries());

  
      // Draw the stack when in range     
      if (min_ieta <= ieta && ieta <= max_ieta && ((ieta-min_ieta)%skip_stack)==0){

	// create the stack
	THStack * st = new THStack(Form("stack_%i",ieta),Form("%4.3f<#eta<%4.3f",eta,etap1));
	met_r->SetFillColor(kBlue);
	met_r->SetLineColor(kBlue);
	met_r->SetFillStyle(1001);
	st->Add(met_r);
	met_t->SetFillColor(kMagenta);
	met_t->SetLineColor(kMagenta);
	met_t->SetFillStyle(1001);
	st->Add(met_t);
	met_w->SetFillColor(kGreen);
	met_w->SetLineColor(kGreen);
	met_w->SetFillStyle(1001);
	st->Add(met_w);
	met_q->SetFillColor(kRed);
	met_q->SetLineColor(kRed);
	met_q->SetFillStyle(1001);
	st->Add(met_q);
	cfits->cd(icfits);
	st->Draw("HIST");
	st->GetXaxis()->SetTitle("MET");
	st->GetXaxis()->SetRangeUser(0,200);
	met_s->Draw("Epsame");
    
	// Draw the legend too
	TLegend *legend=new TLegend(0.66,0.6,0.99,0.9);
	legend->SetFillStyle(0);
	legend->AddEntry(met_s,"Data","lpe");
	legend->AddEntry(met_q,"QCD","f");
	legend->AddEntry(met_w,"WJets","f");
	legend->AddEntry(met_t,"TTbar","f");
	legend->AddEntry(met_r,"Others","f");
	legend->Draw();

	// The eta-region text
	TLatex * lat = new TLatex(0.4,0.96,st->GetTitle());
	lat->SetNDC();
	lat->SetTextSize(0.054);
	lat->Draw();

	icfits++;

      }  

    } else
      cout<<"*** fit failed for ieta="<<ieta<<" with status="<<status<<endl;
   
    //delete fit;

  }// for eta bins

  // also save the chi2 in cfits
  cfits->cd(icfits);
  chi2->GetYaxis()->SetTitle("#chi^{2}/NDF");
  chi2->GetXaxis()->SetTitle("|#eta|");

  chi2->Draw("E");

  // ============ PLOT ALL RESULTS TO CANVASES  ==============
  TCanvas * cint = new TCanvas("Stats","Stats",1100,800);
  cint->Divide(4,2);
  cint->cd(1); da_int->Draw("E");
  cint->cd(2); qcd_int->Draw("E");
  cint->cd(3); wj_int->Draw("E");
  cint->cd(4); ttbar_int->Draw("E");
  cint->cd(5); da_nev->Draw("E");
  cint->cd(6); qcd_nev->Draw("E");
  cint->cd(7); wj_nev->Draw("E");
  cint->cd(8); ttbar_nev->Draw("E");
  cint->SaveAs("Stats.pdf");
  cint->Write();
  da_int->Write();
  qcd_int->Write();
  wj_int->Write();
  ttbar_int->Write();
  da_nev->Write();
  qcd_nev->Write();
  wj_nev->Write();
  ttbar_nev->Write();

  TCanvas * cfr = new TCanvas("fraction","fraction",1200,400);
  cfr->Divide(3,1);
  cfr->cd(1); qcd_fr->Draw("E");
  cfr->cd(2); wj_fr->Draw("E");
  cfr->cd(3); ttbar_fr->Draw("E");
  cfr->SaveAs("FractionsVsEta.pdf");
  cfr->Write();
  qcd_fr->Write();
  wj_fr->Write();
  ttbar_fr->Write();

  TCanvas * cyi = new TCanvas("yields","yields",1200,400);
  cyi->Divide(3,1); 
  cyi->cd(1); qcd_yield->Draw("E");
  cyi->cd(2); wj_yield->Draw("E");
  cyi->cd(3); ttbar_yield->Draw("E");
  cyi->SaveAs("EventYieldsVsEta.pdf");
  cyi->Write();
  qcd_yield->Write();
  wj_yield->Write();
  ttbar_yield->Write();

  TCanvas * cxs = new TCanvas("prod xs","prod xs",1200,400);
  cxs->Divide(3,1); 
  cxs->cd(1); qcd_xs->Draw("E");
  cxs->cd(2); wj_xs->Draw("E");
  cxs->cd(3); ttbar_xs->Draw("E");
  cxs->SaveAs("EventXSVsEta.pdf");
  cxs->Write();
  qcd_xs->Write();
  wj_xs->Write();
  ttbar_xs->Write();

  TCanvas * cwe = new TCanvas("sfs","sfs",1200,400);
  cwe->Divide(3,1); 
  cwe->cd(1); qcd_sf->Draw("E");
  cwe->cd(2); wj_sf->Draw("E"); wjsm_sf->Draw("sameE3"); wj_sf->Draw("sameE");
  cwe->cd(3); ttbar_sf->Draw("E"); ttbarsm_sf->Draw("sameE3"); ttbar_sf->Draw("sameE");
  cwe->SaveAs("SfVsEta.pdf");
  cwe->Write();
  qcd_sf->Write();
  wj_sf->Write();
  ttbar_sf->Write();

  TCanvas * ccl = new TCanvas("Closure","closure" ,1200,400);
  ccl->Divide(3,1); 
  ccl->cd(1); qcd_close->Draw("E");
  ccl->cd(2); wj_close->Draw("E");
  ccl->cd(3); ttbar_close->Draw("E");
  ccl->SaveAs("ClosureVsEta.pdf");
  ccl->Write();
  qcd_close->Write();
  wj_close->Write();
  ttbar_close->Write();

  TCanvas * cchi2 = new TCanvas("chi2","chi2",400,400);
  chi2->Draw("E");
  cchi2->SaveAs("chi2.pdf");
  chi2->Write();
  cchi2->Write();

  cfits->SaveAs("AllMetFits.pdf");
  cfits->SaveAs("AllMetFits.png");
  cfits->Write();
  
  // ============ CREATE THE OUTPUT HISTOS ==============
  TFile *fout = new TFile("QCDWeight_electron.root","RECREATE");
  
  TH1 * qcdweight = (TH1*) qcd_sf->Clone("QCDWeight_electron");
  qcdweight->SetDirectory(0);
  qcdweight->Write();
  fout->Close();

  fout1->Close();

}//QCD_EtaDependent_MetFit
