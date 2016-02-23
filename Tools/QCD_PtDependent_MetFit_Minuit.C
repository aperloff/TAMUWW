//There is a known bug where the Rebin2D function does not correctly recomput the underflow and overflow bins.

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
#include "TMath.h"

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

double vpt[] = {0,30,31,32,33,34,35,36,37,38,39,40,50,60,70,80,90,100,110,120,130,150,170,200,250,500};
//double npt = sizeof(vpt)/sizeof(double);
double npt = 25;

//global so i can use it in the minimization function
TH1D * met_s, * met_q, * met_w, * met_r;
double ndf;

// ----------------------------------------------------------------------------
// The minimization function 
void minimizationFunction(Int_t &npar, Double_t *gin, Double_t &f,
			  Double_t *par, Int_t iflag){


  double met_sint = met_s->Integral();
  double met_qint = met_q->Integral();
  double met_wint = met_w->Integral();


  //calculate chisquare
  Double_t chi2 = 0;
  ndf = 0;
  for (int ibin = 0 ; ibin <= met_s->GetNbinsX()+1 ; ibin++) {

    //skip events without entries
    if (met_s->GetBinContent(ibin) == 0  ||
	met_q->GetBinContent(ibin) == 0  ||
	met_w->GetBinContent(ibin) == 0 )
      continue;

    double bin_num =  met_s->GetBinContent(ibin) 
      - par[0] * met_q->GetBinContent(ibin) * met_sint / met_qint 
      - par[1] * met_w->GetBinContent(ibin) * met_sint / met_wint 
      - met_r->GetBinContent(ibin);
  
    // In the denominator use the error of the bin, which is Poisson scaled
    double bin_den2 = 0 
      //+  pow(met_s->GetBinError(ibin),2) 
      +  pow(par[0] * met_q->GetBinError(ibin) * met_sint / met_qint,2)
      +  pow(par[1] * met_w->GetBinError(ibin) * met_sint / met_wint,2);
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

TH2D* Rebin2D_Pt_simple(TH2* old, int rebinMet){
   //create a new TH2 with your bin arrays spec
   Int_t nbinsy    = old->GetNbinsY();
   Double_t ymin  = old->GetYaxis()->GetXmin();
   Double_t ymax  = old->GetYaxis()->GetXmax();
   if ((rebinMet <= 0) || (rebinMet > nbinsy)) {
      Error("Rebin", "Illegal value of ngroup=%d",rebinMet);
      return 0;
   }
   Int_t newbinsy = nbinsy/rebinMet;
   TH2D *h = new TH2D(old->GetName()+TString("oldrebin"),old->GetTitle(),npt,vpt,newbinsy,ymin,ymax);
   TAxis *xaxis = old->GetXaxis();
   TAxis *yaxis = old->GetYaxis();
   for (int j=1;j<=yaxis->GetNbins();j++) {
      for (int i=1;i<=xaxis->GetNbins();i++) {
         h->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),old->GetBinContent(i,j));
      }
   }

   return h;
}

TH2* Rebin2D(TH2* old, Int_t nxgroup, Int_t nygroup, const char*newname, const Double_t *xbins, const Double_t *ybins, TString Options = "") {

  Options.ToUpper();
  bool verbose = false;
  if(Options.Contains("VERBOSE")) {
    verbose = true;
  }

   Int_t nxbins   = old->GetXaxis()->GetNbins();
   Double_t xmin  = old->GetXaxis()->GetXmin();
   Double_t xmax  = old->GetXaxis()->GetXmax();
   Int_t nybins   = old->GetYaxis()->GetNbins();
   Double_t ymin  = old->GetYaxis()->GetXmin();
   Double_t ymax  = old->GetYaxis()->GetXmax();
   if ((nxgroup <= 0) || (nxgroup > nxbins)) {
      Error("Rebin", "Illegal value of nxgroup=%d",nxgroup);
      return 0;
   }
   if ((nygroup <= 0) || (nygroup > nybins)) {
      Error("Rebin", "Illegal value of nygroup=%d",nygroup);
      return 0;
   }

   if (!newname && (xbins || ybins)) {
      Error("Rebin","if xbins or ybins are specified, newname must be given");
      return 0;
   }

   Int_t newbinsx = nxbins/nxgroup;
   Int_t newbinsy = nybins/nygroup;
   if (!xbins) {
      Int_t nbg = nxbins/nxgroup;
      if (nbg*nxgroup != nxbins) {
         Warning("Rebin", "nxgroup=%d is not an exact divider of nxbins=%d.",nxgroup,nxbins);
      }
   }
   else {
   // in the case that xbins is given (rebinning in variable bins), ngroup is
   // the new number of bins and number of grouped bins is not constant.
   // when looping for setting the contents for the new histogram we
   // need to loop on all bins of original histogram.  Then set ngroup=nbins
      newbinsx = nxgroup;
      nxgroup = nxbins;
   }
   if (!ybins) {
      Int_t nbg = nybins/nygroup;
      if (nbg*nygroup != nybins) {
         Warning("Rebin", "nygroup=%d is not an exact divider of nybins=%d.",nygroup,nybins);
      }
   }
   else {
   // in the case that xbins is given (rebinning in variable bins), ngroup is
   // the new number of bins and number of grouped bins is not constant.
   // when looping for setting the contents for the new histogram we
   // need to loop on all bins of original histogram.  Then set ngroup=nbins
      newbinsy = nygroup;
      nygroup = nybins;
   }

   // Save old bin contents into a new array
   Double_t entries = old->GetEntries();
   Double_t* oldBins = new Double_t[(nxbins+2)*(nybins+2)];
   Int_t i, j, binx, biny;
   for (binx=0;binx<nxbins+2;binx++) {
     for (biny=0;biny<nybins+2;biny++) { 
         oldBins[binx*(nybins+2)+biny] = old->GetBinContent(binx,biny);
     }
   }
   Double_t* oldErrors = 0;
   if (old->GetSumw2N() != 0) {
      oldErrors = new Double_t[(nxbins+2)*(nybins+2)];
      for (binx=0;binx<nxbins+2;binx++) {
          for (biny=0;biny<nybins+2;biny++) {
            oldErrors[binx*(nybins+2)+biny] = old->GetBinError(binx,biny);
          }
      }
   }

   // create a clone of the old histogram if newname is specified
   TH2 *hnew = old;
   if ((newname && strlen(newname) > 0) || xbins || ybins) {
      hnew = (TH2*)old->Clone(newname);
   }

   //reset kCanRebin bit to avoid a rebinning in SetBinContent
   Int_t bitRebin = hnew->TestBit(old->kCanRebin);
   hnew->SetBit(old->kCanRebin,0);

   // save original statistics
   const Int_t kNstat2D = 7;
   Double_t stat[kNstat2D];
   old->GetStats(stat);
   bool resetStat = false;

   // change axis specs and rebuild bin contents array::RebinAx
   if(!xbins && (newbinsx*nxgroup != nxbins)) {
      xmax = old->GetXaxis()->GetBinUpEdge(newbinsx*nxgroup);
      resetStat = true; //stats must be reset because top bins will be moved to overflow bin
   }
   if(!ybins && (newbinsy*nygroup != nybins)) {
      ymax = old->GetYaxis()->GetBinUpEdge(newbinsy*nygroup);
      resetStat = true; //stats must be reset because top bins will be moved to overflow bin
   }
   // save the TAttAxis members (reset by SetBins)
   Int_t    xnDivisions  = old->GetXaxis()->GetNdivisions();
   Color_t  xaxisColor   = old->GetXaxis()->GetAxisColor();
   Color_t  xlabelColor  = old->GetXaxis()->GetLabelColor();
   Style_t  xlabelFont   = old->GetXaxis()->GetLabelFont();
   Float_t  xlabelOffset = old->GetXaxis()->GetLabelOffset();
   Float_t  xlabelSize   = old->GetXaxis()->GetLabelSize();
   Float_t  xtickLength  = old->GetXaxis()->GetTickLength();
   Float_t  xtitleOffset = old->GetXaxis()->GetTitleOffset();
   Float_t  xtitleSize   = old->GetXaxis()->GetTitleSize();
   Color_t  xtitleColor  = old->GetXaxis()->GetTitleColor();
   Style_t  xtitleFont   = old->GetXaxis()->GetTitleFont();
   Int_t    ynDivisions  = old->GetYaxis()->GetNdivisions();
   Color_t  yaxisColor   = old->GetYaxis()->GetAxisColor();
   Color_t  ylabelColor  = old->GetYaxis()->GetLabelColor();
   Style_t  ylabelFont   = old->GetYaxis()->GetLabelFont();
   Float_t  ylabelOffset = old->GetYaxis()->GetLabelOffset();
   Float_t  ylabelSize   = old->GetYaxis()->GetLabelSize();
   Float_t  ytickLength  = old->GetYaxis()->GetTickLength();
   Float_t  ytitleOffset = old->GetYaxis()->GetTitleOffset();
   Float_t  ytitleSize   = old->GetYaxis()->GetTitleSize();
   Color_t  ytitleColor  = old->GetYaxis()->GetTitleColor();
   Style_t  ytitleFont   = old->GetYaxis()->GetTitleFont();

   if(!xbins && !ybins) { //&& (old->GetXaxis()->GetXbins()->GetSize() > 0) && (old->GetYaxis()->GetXbins()->GetSize() > 0)){ // variable bin sizes
      if (verbose) cout << "Case 1" << endl;
      Double_t *binsx = new Double_t[newbinsx+1];
      Double_t *binsy = new Double_t[newbinsy+1];
      for(i = 0; i <= newbinsx; ++i) {
        binsx[i] = old->GetXaxis()->GetBinLowEdge(1+i*nxgroup);
      }
      for(i = 0; i <= newbinsy; ++i) {
        binsy[i] = old->GetYaxis()->GetBinLowEdge(1+i*nygroup);
      }
      hnew->SetBins(newbinsx,binsx,newbinsy,binsy); //this also changes errors array (if any)
      delete [] binsx;
      delete [] binsy;
   } else if (!xbins && ybins) { //&& (old->GetXaxis()->GetXbins()->GetSize() > 0)) {
      if (verbose) cout << "Case 2" << endl;
      Double_t *binsx = new Double_t[newbinsx+1];
      for(i = 0; i <= newbinsx; ++i) {
        binsx[i] = old->GetXaxis()->GetBinLowEdge(1+i*nxgroup);
      }
      hnew->SetBins(newbinsx,binsx,newbinsy,ybins); //this also changes errors array (if any)
      delete [] binsx;
   } else if (!ybins && xbins) {//&& (old->GetYaxis()->GetXbins()->GetSize() > 0)) {
      if (verbose) cout << "Case 3" << endl;
      Double_t *binsy = new Double_t[newbinsy+1];
      for(i = 0; i <= newbinsy; ++i) {
        binsy[i] = old->GetYaxis()->GetBinLowEdge(1+i*nygroup);
      }
      hnew->SetBins(newbinsx,xbins,newbinsy,binsy); //this also changes errors array (if any)
      delete [] binsy;
   } else if (xbins && ybins) {
      if (verbose) cout << "Case 4" << endl;
      hnew->SetBins(newbinsx,xbins,newbinsy,ybins);
   } else {
      if (verbose) cout << "Case 5" << endl;
      hnew->SetBins(newbinsx,xmin,xmax,newbinsy,ymin,ymax);
   }

   // Restore axis attributes
   old->GetXaxis()->SetNdivisions(xnDivisions);
   old->GetXaxis()->SetAxisColor(xaxisColor);
   old->GetXaxis()->SetLabelColor(xlabelColor);
   old->GetXaxis()->SetLabelFont(xlabelFont);
   old->GetXaxis()->SetLabelOffset(xlabelOffset);
   old->GetXaxis()->SetLabelSize(xlabelSize);
   old->GetXaxis()->SetTickLength(xtickLength);
   old->GetXaxis()->SetTitleOffset(xtitleOffset);
   old->GetXaxis()->SetTitleSize(xtitleSize);
   old->GetXaxis()->SetTitleColor(xtitleColor);
   old->GetXaxis()->SetTitleFont(xtitleFont);
   old->GetYaxis()->SetNdivisions(ynDivisions);
   old->GetYaxis()->SetAxisColor(yaxisColor);
   old->GetYaxis()->SetLabelColor(ylabelColor);
   old->GetYaxis()->SetLabelFont(ylabelFont);
   old->GetYaxis()->SetLabelOffset(ylabelOffset);
   old->GetYaxis()->SetLabelSize(ylabelSize);
   old->GetYaxis()->SetTickLength(ytickLength);
   old->GetYaxis()->SetTitleOffset(ytitleOffset);
   old->GetYaxis()->SetTitleSize(ytitleSize);
   old->GetYaxis()->SetTitleColor(ytitleColor);
   old->GetYaxis()->SetTitleFont(ytitleFont);

   // copy merged bin contents (ignore under/overflows)
   // Start merging only once the new lowest edge is reached
   if (nxgroup != 1 || nygroup != 1) {
     Int_t startbinx = 1;
     Int_t startbiny = 1;
     const Double_t newxmin = hnew->GetXaxis()->GetBinLowEdge(1);
     const Double_t newymin = hnew->GetYaxis()->GetBinLowEdge(1);
     while( old->GetXaxis()->GetBinCenter(startbinx) < newxmin && startbinx <= nxbins ) {
      startbinx++;
    }
    while( old->GetYaxis()->GetBinCenter(startbiny) < newymin && startbiny <= nybins ) {
      startbiny++;
    }
    Int_t oldbinx = startbinx;
    Int_t oldbiny = startbiny;
    if(verbose) {
      cout << "startbinx = " << oldbinx << endl;
      cout << "startbiny = " << oldbiny << endl;
    }

    Double_t binContent, binError;

    for (binx = 1;binx<=newbinsx;binx++) {
      oldbiny = startbiny;
      Int_t ixmax = nxgroup;
      Double_t xbinmax = hnew->GetXaxis()->GetBinUpEdge(binx);
      for (biny = 1;biny<=newbinsy;biny++) {
        binContent = 0;
        binError   = 0;
        Int_t iymax = nygroup;
        Double_t ybinmax = hnew->GetYaxis()->GetBinUpEdge(biny);
        for (i=0;i<nxgroup;i++) {
          if (verbose) cout << "i = " << i << endl;
          if( ((hnew == old) && (oldbinx+i > nxbins)) || ((hnew != old) && (old->GetXaxis()->GetBinCenter(oldbinx+i) > xbinmax))) {
            ixmax = i;
            if(verbose) {
              cout << "WARNING::Before X Break!!!!" << endl;
              //cout << "old->GetXaxis()->GetBinCenter(oldbinx+i) > xbinmax\t" << (old->GetXaxis()->GetBinCenter(oldbinx+i) > xbinmax) << endl;
              //cout << "old->GetXaxis()->GetBinCenter(oldbinx+i) = " << old->GetXaxis()->GetBinCenter(oldbinx+i) << endl;
              //cout << "xbinmax = " << xbinmax << endl;
              //cout << "oldbinx = " << oldbinx << endl;
              //cout << "i = " << i << endl;
              //cout << "nxgroup = " << nxgroup << endl;
            }
            break;
          }
          for (j=0;j<nygroup;j++) {
            if (verbose) cout << "j = " << j << endl;
            if( ((hnew == old) && (oldbiny+j > nybins)) || ((hnew != old) && (old->GetYaxis()->GetBinCenter(oldbiny+j) > ybinmax))) {
              iymax = j;
              if(verbose) {
                cout << "WARNING::Before Y Break!!!!" << endl;
              //cout << "hnew==old = " << (hnew==old) << endl;
              //cout << "oldbinx+i > nxbins || oldbiny+j > nybins\t" << (oldbinx+i > nxbins || oldbiny+j > nybins) << endl;
              //cout << "old->GetYaxis()->GetBinCenter(oldbiny+j) > ybinmax\t" << (old->GetYaxis()->GetBinCenter(oldbiny+j) > ybinmax) << endl;
              //cout << "old->GetYaxis()->GetBinCenter(oldbiny+j) = " << old->GetYaxis()->GetBinCenter(oldbiny+j) << endl;
              //cout << "ybinmax = " << ybinmax << endl;
              //cout << "oldbiny = " << oldbiny << endl;
              //cout << "j = " << j << endl;
              //cout << "nygroup = " << nygroup << endl;
              }
              break;
            }
            binContent += oldBins[(oldbiny+j) + (oldbinx+i)*(nybins+2)];
            if (oldErrors) binError += oldErrors[(oldbiny+j)+(oldbinx+i)*(nybins+2)]*oldErrors[(oldbiny+j)+(oldbinx+i)*(nybins+2)];
          }
        }
        if (verbose) cout << "binx = " << binx << "\tbiny = " << biny << endl;
        hnew->SetBinContent(binx,biny,binContent);
        if (oldErrors) hnew->SetBinError(binx,biny,TMath::Sqrt(binError));
        oldbiny += iymax;
      }
      oldbinx += ixmax;
    }

    //  recompute under/overflow contents in y for the new  x bins 
    Double_t binContent0, binContent2;
    Double_t binError0, binError2;
    oldbinx = 1;
    for (binx = 1; binx<=newbinsx; binx++) {
     binContent0 = binContent2 = 0;
     binError0 = binError2 = 0;
     for (i=0; i<nxgroup; i++) {
      if (oldbinx+i > nxbins) break;
            //N.B  convention used for index is opposite than TH1::GetBin(ix,iy)
            Int_t ufbin = (oldbinx+i)*(nybins+2);   // index for y underflow bins 
            Int_t ofbin = (oldbinx+i)*(nybins+2) + (nybins+1);   // index for y overflow bins 
            binContent0 += oldBins[ufbin];
            binContent2 += oldBins[ofbin];
            if (oldErrors)  { 
             binError0 += oldErrors[ufbin] * oldErrors[ufbin];
             binError2 += oldErrors[ofbin] * oldErrors[ofbin];
           }
         }
         hnew->SetBinContent(binx,0,binContent0);
         hnew->SetBinContent(binx,newbinsy+1,binContent2);
         if (oldErrors) { 
          hnew->SetBinError(binx,0,TMath::Sqrt(binError0));
          hnew->SetBinError(binx,newbinsy+1,TMath::Sqrt(binError2) );
        }
        oldbinx += nxgroup;
      }

      //  recompute under/overflow contents in x for the new y bins
      Int_t oldybin = 1;
      for (biny = 1; biny<=newbinsy; biny++) {
       binContent0 = binContent2 = 0;
       binError0 = binError2 = 0;
       for (i=0; i<nygroup; i++) {
        if (oldbiny+i > nybins) break;
            Int_t ufbin = (oldbiny+i);   // global index for x underflow bins 
            Int_t ofbin = (nxbins+1)*(nybins+2) + (oldbiny+i);   // global index for x overflow bins 
            binContent0 += oldBins[ufbin];
            binContent2 += oldBins[ofbin];
            if (oldErrors)  { 
             binError0 += oldErrors[ufbin] * oldErrors[ufbin];
             binError2 += oldErrors[ofbin] * oldErrors[ofbin];
           }
         }
         hnew->SetBinContent(0,biny,binContent0);
         hnew->SetBinContent(newbinsx+1,biny,binContent2);
         if (oldErrors) { 
          hnew->SetBinError(0,biny, TMath::Sqrt(binError0));
          hnew->SetBinError(newbinsx+1, biny, TMath::Sqrt(binError2));
        }
        oldbiny += nygroup;
      }
    }

   // restore statistics and entries modified by SetBinContent
   hnew->SetEntries(entries);
   if (!resetStat) hnew->PutStats(stat);

   delete [] oldBins;
   if (oldErrors) delete [] oldErrors;
   return hnew;
}

void Rebin2DTest(TString Options = "") {
  /*
  Options:
    Verbose
  */
  TH2D* h = new TH2D("h","h",4,0,40,4,0,40);
  h->Sumw2();
  h->Fill(5,5);
  h->Fill(15,15);
  h->Fill(25,25);
  h->Fill(35,25);
  h->Fill(35,10);
  h->Fill(35,10);
  h->Fill(35,10);
  h->Fill(-1,-1);
  h->Fill(50,50);
  h->Fill(-1,50);
  h->Fill(50,-1);
  Double_t x[3] = {0,30,40};
  TH2D* hh = (TH2D*)Rebin2D(h,2,2,"test",x,0,Options);
  //Print options: range (no uf/of), all (with uf/of)
  h->Print("all");
  hh->Print("all");
}


void getHistosFromDrawOutput(TH2D *& hd, TH2D *& hw, TH2D *& hq,  int rebinPt, int rebinMet){

  // Data: open the input file and get all histos
  TFile *_file0 = TFile::Open("output_signal.root");
  hd =  (TH2D*) _file0->Get("h2");
  hd->Rebin2D(rebinPt, rebinMet);

  // WJets: open the input file and get all histos
  TFile *_file1 = TFile::Open("output_wjets.root");
  hw =  (TH2D*) _file1->Get("h2");
  hw->Rebin2D(rebinPt, rebinMet);

  // QCD: open the input file and get all histos
  TFile *_file2 = TFile::Open("output_nomvacutnew.root");
  //TFile *_file2 = TFile::Open("output_mvacut.root");
  hq =  (TH2D*) _file2->Get("h2");
  hq->Rebin2D(rebinPt, rebinMet);

}//getHistosFromDrawOutput

void getHistosFromPlotterOutput(TH2D *& hd, TH2D *& hw, TH2D *& hq, TH2D  *& hr, int rebinPt, int rebinMet, bool full){

  bool removesw2 = false;// true;

  // list of background processes. WJets and QCD are handled separately
  TString hn = "MET_vs_LeptPt_";

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

  // Vector containing all other histos. Proper normalization 
  // to the luminosity of the signal data is assumed
  vector<TString> backs;
  backs.push_back(hn+"STopS_T_electron");
  backs.push_back(hn+"STopS_Tbar_electron");
  backs.push_back(hn+"STopT_T_electron");
  backs.push_back(hn+"STopT_Tbar_electron");
  backs.push_back(hn+"STopTW_T_electron");
  backs.push_back(hn+"STopTW_Tbar_electron");
  backs.push_back(hn+"TTbar_electron");
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
  //hd->Rebin2D(rebinPt, rebinMet);
  hd = (TH2D*)Rebin2D(hd,npt,rebinMet,hd->GetName(),vpt,0,"");

  //WJets
  hw = (TH2D*) _file0->Get(wj);
  if (removesw2) hw->GetSumw2()->Set(0); 
  //hw->Rebin2D(rebinPt, rebinMet);
  hw = (TH2D*)Rebin2D(hw,npt,rebinMet,hw->GetName(),vpt,0,"");

  //QCD 
  hq = (TH2D*) _file0->Get(qcd);
  if (removesw2) hq->GetSumw2()->Set(0); 
  //hq->Rebin2D(rebinPt, rebinMet);
  hq = (TH2D*)Rebin2D(hq,npt,rebinMet,hq->GetName(),vpt,0,"");

  // Remainding Processes
  // Substract all other backgrounds from the data
  hr = 0;
  for (unsigned int b = 0 ; b < backs.size();b++){
    TH2D * aux = (TH2D*) _file0->Get(backs[b]);
    //aux->Rebin2D(rebinPt,rebinMet);
    aux = (TH2D*)Rebin2D(aux,npt,rebinMet,aux->GetName(),vpt,0,"");
    if (hr==0) hr = aux;
    else hr->Add(aux);
  }//for

}// getHistosFromPlotterOutput


void QCD_PtDependent_MetFit_Minuit(bool full = true){

  // ============ GET THE INPUT HISTOS ==============
  int rebinMet = 20; // best value(s): 5 or 20
  int rebinPt = 40; // best value(s): 40
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
  TFile *fout1 = new TFile("QCD_PtDependent_MetFit_output.root","RECREATE");

  // Data, WJets, QCD
  TH2D * hd, *hw, *hq, *hr;
  getHistosFromPlotterOutput(hd,hw,hq,hr,rebinPt,rebinMet,full);
  //getHistosFromDrawOutput(hd,hw,hq, rebinPt,rebinMet);

  fout1->cd();

  // set the styles
  hd->UseCurrentStyle();
  hw->UseCurrentStyle();
  hq->UseCurrentStyle();
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

  // The relative qcd fraction according to the fit 
  TH1D * qcd_fr = hd->ProjectionX("qcd_fr");
  qcd_fr->SetTitle("qcd_fr");
  qcd_fr->Reset();

  // The relative wjets fraction according to the fit 
  TH1D * wj_fr = hd->ProjectionX("wj_fr");
  wj_fr->SetTitle("wj_fr");
  wj_fr->Reset();

  // The chi2 of the fit 
  TH1D * chi2 = hd->ProjectionX("chi2");
  chi2->SetTitle("chi2");
  chi2->Reset();


  // The qcd yield according to the fit 
  TH1D * qcd_yield = hd->ProjectionX("qcd_yield");
  qcd_yield->SetTitle("qcd_yield;p_{T} [GeV];number of events / 19148 fb^{-1}");
  qcd_yield->Reset();

  // The wjets yield according to the fit 
  TH1D * wj_yield = hd->ProjectionX("wj_yield");
  wj_yield->SetTitle("wj_yield;p_{T} [GeV];number of events / 19148 fb^{-1}");
  wj_yield->Reset();

  // The qcd production xs according to the fit 
  TH1D * qcd_xs = hd->ProjectionX("qcd_xs");
  qcd_xs->SetTitle("qcd_xs");
  qcd_xs->Reset();

  // The wjets production xs according to the fit 
  TH1D * wj_xs = hd->ProjectionX("wj_xs");
  wj_xs->SetTitle("wj_xs");
  wj_xs->Reset();

  // The qcd production xs according to the fit 
  TH1D * qcd_sf = hd->ProjectionX("qcd_sf");
  qcd_sf->SetTitle(";p_{T} [GeV];s_{QCD}");
  qcd_sf->Reset();

  // The wjets production xs according to the fit 
  TH1D * wj_sf = hd->ProjectionX("wj_sf");
  wj_sf->SetTitle(";p_{T} [GeV];N^{WJets}_{measured} /N^{WJets}_{expected} ");
  wj_sf->Reset();

  TGraphErrors *wjsm_sf = new TGraphErrors();
  wjsm_sf->SetPoint(0,0.05,1.0);
  wjsm_sf->SetPointError(0,0,0.0256);
  wjsm_sf->SetPoint(1,500.0,1.0);
  wjsm_sf->SetPointError(1,0,0.0256);

  wjsm_sf->SetFillColor(kGreen);
  wjsm_sf->SetMarkerColor(0);
  

  // The qcd production xs according to the fit 
  TH1D * qcd_close = hd->ProjectionX("qcd_close");
  qcd_close->SetTitle("qcd_close;p_{T} [GeV]");
  qcd_close->Reset();

  // The wjets production xs according to the fit 
  TH1D * wj_close = hd->ProjectionX("wj_close");
  wj_close->SetTitle("wj_close;p_{T} [GeV]");
  wj_close->Reset();


  // ============ AUXILIARY VARIABLES ==============

  // Show stack every <skip_stack> fits in range pt<500
  int skip_stack = 1;
  int min_ipt = qcd_fr->FindBin(0);
  int max_ipt = qcd_fr->FindBin(500);

  // the max number of canvases to be drawn is then 
  int ncan = max_ipt-min_ipt+1;

  //Create the Canvas and divide it
  TCanvas * cfits = new TCanvas("fits","fits",1200,1200);
  int ncols = ceil(sqrt(ncan));
  int nrows = ceil(double(ncan)/ncols);
  cfits->Divide(ncols,nrows);
  int icfits=1; // counter for canvas pad


  // ============ LOOP OVER PtS AND DO THE MET FITS ==============
  for (int ipt = min_ipt; ipt <= max_ipt ; ipt++){
      
    double pt = hd->GetXaxis()->GetBinLowEdge(ipt);
    double ptp1 = hd->GetXaxis()->GetBinLowEdge(ipt+1);

    cout<<"******* Fitting ipt="<<ipt<<": "<<pt<<" << pt <<"<< ptp1 << endl;

    // get the met distributions for that pt
    met_s = hd->ProjectionY(Form("s_ipt%i",ipt),ipt,ipt);
    met_q = hq->ProjectionY(Form("q_ipt%i",ipt),ipt,ipt);
    met_w = hw->ProjectionY(Form("w_ipt%i",ipt),ipt,ipt);
    met_r = hr->ProjectionY(Form("r_ipt%i",ipt),ipt,ipt);

    if ( met_s->GetEntries()==0 || met_q->GetEntries()==0 || met_w->GetEntries()==0 ){
      cout<<"\t skip this pt due to empty entries"<<endl; 
      continue;
    }


    // ============ SETUP MINUIT
    TMinuit minuit(2);
  
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
      minuit.mnparm( 0, "qcd_sf", 0.03 , 0.01 , 0 , 1 , ierflg); // start at 3%
    
    minuit.mnparm( 1, "wj_sf", 0.75 , 0.01 , 0 , 1 , ierflg); // start at 75%
  
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
    
      // get the results of the fit here
    
      // Retrieve the central value (central) and errors (up,down)
      double central, err_down, err_up, junk, junk1, junk2;
      minuit.GetParameter(0, val0, junk);
      minuit.mnerrs(0, err_up, err_down, junk1, junk2);
      err0 =  0.5*(fabs(err_up)+fabs(err_down));

      minuit.GetParameter(1, val1, junk);
      minuit.mnerrs(0, err_up, err_down, junk1, junk2);
      err1 =  0.5*(fabs(err_up)+fabs(err_down));

      cout<<"\t qcd_sf ="<<val0<<" +/- "<<err0<<endl;
      cout<<"\t wj_sf  ="<<val1<<" +/- "<<err1<<endl;

      double met_sint = met_s->Integral();
      double met_qint = met_q->Integral();
      double met_wint = met_w->Integral();

      Double_t amin,edm,errdef;
      Int_t nvpar,nparx,icstat;
      minuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
 

      double redchi2 = amin/ndf; 
      chi2->SetBinContent(ipt,redchi2);
      chi2->SetBinError(ipt,redchi2/sqrt(ndf));   
      cout<<"\t chi2   = "<<amin<<" ndf = "<<ndf<<endl;

      // fill QCD info
      qcd_int->SetBinContent(ipt, met_q->Integral());
      qcd_nev->SetBinContent(ipt, met_q->GetEntries());
      qcd_fr->SetBinContent(ipt,val0);
      qcd_fr->SetBinError(ipt,err0);

      qcd_yield->SetBinContent(ipt,val0*met_s->Integral());
      double rel2 = pow(err0/val0,2) + 1.0/met_s->Integral();
      qcd_yield->SetBinError(ipt,sqrt(rel2) * qcd_yield->GetBinContent(ipt));

      qcd_xs->SetBinContent(ipt, val0*met_s->Integral()/lum_data);
      qcd_xs->SetBinError(ipt, err0*met_s->Integral()/lum_data);

      qcd_sf->SetBinContent(ipt, val0*met_s->Integral()*lum_qcd / (lum_data*met_q->GetEntries()));
      rel2 = pow(err0/val0,2) 
	+ 1./met_s->GetEntries()
	+ 1./met_q->GetEntries();
      qcd_sf->SetBinError  (ipt, sqrt(rel2) * qcd_sf->GetBinContent(ipt));

      qcd_close->SetBinContent (ipt, val0*met_s->Integral()*lum_qcd / (lum_data*met_q->Integral()));
      qcd_close->SetBinError(ipt, sqrt(rel2)*qcd_close->GetBinContent(ipt));

      //double rel2 = pow(err0/val0,2) + pow(sqrt(met_q->Integral())/met_q->Integral(),2);
      //qcd_xs->SetBinError(ipt, sqrt(rel2) * qcd_xs->GetBinContent(ipt));
      met_q->Scale(val0*met_s->Integral()/met_q->Integral()); // scale to plot on stack

      // fill WJets info
      wj_int->SetBinContent(ipt, met_w->Integral());
      wj_nev->SetBinContent(ipt, met_w->GetEntries());
      wj_fr->SetBinContent(ipt,val1);
      wj_fr->SetBinError(ipt,err1);
      wj_yield->SetBinContent(ipt,val1*met_s->Integral());
      wj_yield->SetBinError(ipt,err1*met_s->Integral());
      wj_xs->SetBinContent(ipt, val1*met_s->Integral()/lum_data);
      wj_xs->SetBinError(ipt, err1*met_s->Integral()/lum_data);
      wj_sf->SetBinContent(ipt, val1*met_s->Integral() /met_w->Integral());
      rel2 = pow(err1/val1,2) 
	+ 1./met_s->GetEntries()
	+ 1./met_w->GetEntries()
	+ pow(960./37509.,2); // added error on expected WJets XS
      wj_sf->SetBinError(ipt, sqrt(rel2)*wj_sf->GetBinContent(ipt));
      wj_close->SetBinContent(ipt, val1*met_s->Integral() /met_w->Integral());
      rel2 = pow(err1/val1,2)
	+ 1.0/met_s->GetEntries()
	+ 1.0/met_w->GetEntries();
      wj_close->SetBinError(ipt, sqrt(rel2)*wj_close->GetBinContent(ipt));
      cout<<"\t ipt="<<ipt<<" wj closure ="<< wj_close->GetBinContent(ipt)
	  <<" +/- "<<wj_close->GetBinError(ipt)<<endl;
      
      // rel2 = pow(err1/val1,2) + pow(sqrt(met_w->Integral())/met_w->Integral(),2);
      // wj_xs->SetBinError(ipt, sqrt(rel2) * wj_xs->GetBinContent(ipt));
      met_w->Scale(val1*met_s->Integral()/met_w->Integral()); // scale to plot on stack


      da_int->SetBinContent(ipt, met_s->Integral());
      da_nev->SetBinContent(ipt, met_s->GetEntries());

  
      // Draw the stack when in range     
      if (min_ipt <= ipt && ipt <= max_ipt && ((ipt-min_ipt)%skip_stack)==0){

	// create the stack
	THStack * st = new THStack(Form("stack_%i",ipt),Form("%4.3f<p_{T}<%4.3f",pt,ptp1));
	met_r->SetFillColor(kBlue);
	met_r->SetLineColor(kBlue);
	met_r->SetFillStyle(1001);
	st->Add(met_r);
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
	legend->AddEntry(met_r,"Others","f");
	legend->Draw();

	// The pt-region text
	TLatex * lat = new TLatex(0.4,0.96,st->GetTitle());
	lat->SetNDC();
	lat->SetTextSize(0.054);
	lat->Draw();

	icfits++;

      }  

    } else
      cout<<"*** fit failed for ipt="<<ipt<<" with status="<<status<<endl;
   
    //delete fit;

  }// for pt bins

  // also save the chi2 in cfits
  cfits->cd(icfits);
  chi2->GetYaxis()->SetTitle("#chi^{2}/NDF");
  chi2->GetXaxis()->SetTitle("p_{T} [GeV]");

  chi2->Draw("E");

  // ============ PLOT ALL RESULTS TO CANVASES  ==============
  TCanvas * cint = new TCanvas("Stats","Stats",800,800);
  cint->Divide(3,2);
  cint->cd(1); da_int->Draw("E");
  cint->cd(2); qcd_int->Draw("E");
  cint->cd(3); wj_int->Draw("E");
  cint->cd(4); da_nev->Draw("E");
  cint->cd(5); qcd_nev->Draw("E");
  cint->cd(6); wj_nev->Draw("E");
  cint->SaveAs("Stats.pdf");
  cint->SaveAs("Stats.eps");
  cint->Write();
  da_int->Write();
  qcd_int->Write();
  wj_int->Write();
  da_nev->Write();
  qcd_nev->Write();
  wj_nev->Write();

  TCanvas * cfr = new TCanvas("fraction","fraction",800,400);
  cfr->Divide(2,1);
  cfr->cd(1); qcd_fr->Draw("E");
  cfr->cd(2); wj_fr->Draw("E");
  cfr->SaveAs("FractionsVsPt.pdf");
  cfr->SaveAs("FractionsVsPt.eps");
  cfr->Write();
  qcd_fr->Write();
  wj_fr->Write();

  TCanvas * cyi = new TCanvas("yields","yields",800,400);
  cyi->Divide(2,1); 
  cyi->cd(1); qcd_yield->Draw("E");
  cyi->cd(2); wj_yield->Draw("E");
  cyi->SaveAs("EventYieldsVsPt.pdf");
  cyi->SaveAs("EventYieldsVsPt.eps");
  cyi->Write();
  qcd_yield->Write();
  wj_yield->Write();

  TCanvas * cxs = new TCanvas("prod xs","prod xs",800,400);
  cxs->Divide(2,1); 
  cxs->cd(1); qcd_xs->Draw("E");
  cxs->cd(2); wj_xs->Draw("E");
  cxs->SaveAs("EventXSVsPt.pdf");
  cxs->SaveAs("EventXSVsPt.eps");
  cxs->Write();
  qcd_xs->Write();
  wj_xs->Write();

  TCanvas * cwe = new TCanvas("sfs","sfs",800,400);
  cwe->Divide(2,1); 
  cwe->cd(1); qcd_sf->Draw("E");
  cwe->cd(2); wj_sf->Draw("E"); wjsm_sf->Draw("sameE3"); wj_sf->Draw("sameE");
  cwe->SaveAs("SfVsPt.pdf");
  cwe->SaveAs("SfVsPt.eps");
  cwe->Write();
  qcd_sf->Write();
  wj_sf->Write();

  TCanvas * ccl = new TCanvas("Closure","closure" ,800,400);
  ccl->Divide(2,1); 
  ccl->cd(1); qcd_close->Draw("E");
  ccl->cd(2); wj_close->Draw("E");
  ccl->SaveAs("ClosureVsPt.pdf");
  ccl->SaveAs("ClosureVsPt.eps");
  ccl->Write();
  qcd_close->Write();
  wj_close->Write();

  TCanvas * cchi2 = new TCanvas("chi2","chi2",400,400);
  chi2->Draw("E");
  cchi2->SaveAs("chi2.pdf");
  cchi2->SaveAs("chi2.eps");
  chi2->Write();
  cchi2->Write();

  cfits->SaveAs("AllMetFits.pdf");
  cfits->SaveAs("AllMetFits.png");
  cfits->SaveAs("AllMetFits.eps");
  cfits->Write();
  
  // ============ CREATE THE OUTPUT HISTOS ==============
  TFile *fout = new TFile("QCDWeight_electron.root","RECREATE");
  
  TH1 * qcdweight = (TH1*) qcd_sf->Clone("QCDWeight_electron");
  qcdweight->SetDirectory(0);
  qcdweight->Write();
  fout->Close();

  fout1->Close();

}//QCD_PtDependent_MetFit
