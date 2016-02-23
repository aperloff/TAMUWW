////////////////////////////////////////////////////////////////////////////////
//
// MakeWPtWeights_x
// ----------------
//
//                          08/04/2015 Alexx Perloff <aperloff@Physics.tamu.edu>
////////////////////////////////////////////////////////////////////////////////

// ROOT libraries
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TString.h"
#include "TLegend.h"
#include <TVirtualFitter.h>
#include <TMath.h>
#include "TSpectrum.h"
#include "TLine.h"

// C++ libraries
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

// My libraries
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"
#include "TAMUWW/SpecialTools/interface/Defs.hh"

#define nCballPar 8 //#pars for cball fn

using namespace std;
using namespace TMath;

////////////////////////////////////////////////////////////////////////////////
// define local functions
////////////////////////////////////////////////////////////////////////////////

/// default fit with gaussian in niter iteration of mean
void fit_gaussian(TH1F*& hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter);

/// optional double sided crystal ball fit to response distributions
int fit_dscb(TH1F*& hrsp,
             const double nsigma,
             const double jtptmin,
             const int niter,
             const string alg);

/// double sided crystal ball function definition
double fnc_dscb(double*xx,double*pp);

/// test this...
void guesstimate_fitrange(TH1* h,double& min,double& max,const string alg);


/// check if a vector of strings contains a certain element
bool contains(const vector<string>& collection,const string& element);
void adjust_fitrange(TH1* h,double& min,double& max);
template <class T>
bool from_string(T& t, 
                 const std::string& s, 
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}

//My implementation of a one sided Crystal Ball function
double CrystalBall (double* x, double* par);

//Kevin Pedro's implementation of a double-sided Crystal Ball function
Double_t cball(Double_t *x, Double_t *par);



////////////////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////////////////

int main(int argc,char**argv)
{
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;

   TString basepath = cl.getValue<TString> ("basepath", "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_04_27_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal_AllPlots/");
   TString rebin    = cl.getValue<TString> ("rebin",   "rebin");
   TString fittype  = cl.getValue<TString> ("fittype", "Kevin");

   if (!cl.check()) 
      return 0;
   cl.print();
   
   vector<TString> jetBins;
   jetBins.push_back("Jets2");
   jetBins.push_back("Jets3");
   jetBins.push_back("Jets4");

   vector<TString> leptons;
   leptons.push_back("electron");
   leptons.push_back("muon");

   vector<TString> shapes;
   shapes.push_back("Ptlv");

   TString processNames[] = {"ggH125","qqH125","WH125_HToBB","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125",
                             "WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","TTH_HToBB_M125",
                             "WJets","ZJets","TTbar","QCD_ElFULL","QCD_MuFULL",
                             "WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","WW"};
   vector<TString> processes(processNames, processNames + sizeof(processNames) / sizeof(TString));

   TFile* ofile = TFile::Open("./WPtWeights.root","RECREATE");

   for(unsigned int j = 0; j<jetBins.size(); j++) {
      for(unsigned int l = 0; l<leptons.size(); l++) {
         TFile* ifile = TFile::Open(basepath+"/"+jetBins[j]+"/"+leptons[l]+"/histos_"+leptons[l]+".root","READ");
         for(unsigned int s = 0; s<shapes.size(); s++) {
            TString canName = Form("WPtRatio_%s_%s_%s",jetBins[j].Data(),leptons[l].Data(),shapes[s].Data());
            cout << "Making canvas " << canName << " ..." << endl;
            TCanvas * c = new TCanvas(canName,canName,600,600);
            TPad* pad = (TPad*)c->GetPad(0);
            TLegend * leg = new TLegend(0.6,0.6,0.9,0.9);

            TH1D* data;
            if(leptons[l].CompareTo("electron")==0)
               data = (TH1D*)ifile->Get(Form("%s_%s_%s",shapes[s].Data(),"SingleEl_Data",leptons[l].Data()));
            else if(leptons[l].CompareTo("muon")==0)
               data = (TH1D*)ifile->Get(Form("%s_%s_%s",shapes[s].Data(),"SingleMu_Data",leptons[l].Data()));

            TH1D* WJets;

            for(unsigned int p = 0; p<processes.size(); p++) {
               if(leptons[l].CompareTo("electron")==0 && processes[p].CompareTo("QCD_MuFULL")==0)
                  continue;
               if(leptons[l].CompareTo("muon")==0 && processes[p].CompareTo("QCD_ElFULL")==0)                  
                  continue;

               //Get input histograms
               TString ihname = Form("%s_%s_%s",shapes[s].Data(),processes[p].Data(),leptons[l].Data());
               cout << "\tGetting histogram " << ihname << " ..." << endl;
               ifile->cd();
               assert(ifile->Get(ihname));
               TH1D* tmp = (TH1D*)ifile->Get(ihname);
               if(tmp) {
                  if(processes[p].CompareTo("WJets")!=0) {
                     data->Add(tmp,-1.0);
                  }
                  else {
                     WJets = tmp;
                  }
               }
            }

            
            TH1F* ratio;
            TString ratioName = Form("ratioDataToWJets_%s_%s_%s",jetBins[j].Data(),leptons[l].Data(),shapes[s].Data());
            if(rebin.CompareTo("rebin")==0) {
               data->Rebin(4);
               WJets->Rebin(4);
               ratio = (TH1F*)data->Clone(ratioName);
               ratio->Divide(WJets);
            }
            else if (rebin.CompareTo("variableBinning")==0){
               //Set new bins for histograms to smooth tails
               //Removes the errors right now.
               vector<double> bin_boundaries;
               for(int ibin=0; ibin<data->GetNbinsX()+1; ibin++) {
                  if(data->GetBinLowEdge(ibin)<150.0)
                     bin_boundaries.push_back(data->GetBinLowEdge(ibin));
                  else
                     break;
               }
               bin_boundaries.push_back(200);
               bin_boundaries.push_back(250);
               bin_boundaries.push_back(300);
               bin_boundaries.push_back(350);
               bin_boundaries.push_back(400);
               bin_boundaries.push_back(450);
               bin_boundaries.push_back(500);

               TH1F* data_rebinned = new TH1F("data_rebinned","data_rebinned",bin_boundaries.size()-1,&bin_boundaries[0]);
               TH1F* wjets_rebinned = new TH1F("wjets_rebinned","wjets_rebinned",bin_boundaries.size()-1,&bin_boundaries[0]);
               for(int ibin=0; ibin<data->GetNbinsX()+1; ibin++) {
                  double data_bin_center = data->GetBinCenter(ibin);
                  double data_bin_content = data->GetBinContent(ibin);
                  double data_old_rebinned_value = data_rebinned->GetBinContent(data_rebinned->FindBin(data_bin_center));
                  data_rebinned->SetBinContent(ibin,data_old_rebinned_value+data_bin_content);

                  double wjets_bin_content = WJets->GetBinContent(ibin);
                  double wjets_old_rebinned_value = wjets_rebinned->GetBinContent(wjets_rebinned->FindBin(data_bin_center));
                  wjets_rebinned->SetBinContent(ibin,wjets_old_rebinned_value+wjets_bin_content);
               }
               ratio = (TH1F*)data_rebinned->Clone(ratioName);
               ratio->Divide(wjets_rebinned);
            }
            else {
               ratio = (TH1F*)data->Clone(ratioName);
               ratio->Divide(WJets);
            }

            ratio->Draw();
            
            if(fittype.CompareTo("JRA")==0) {
               //JRA implementation
               int fitstatus(0);
               //fit_gaussian(ratio,1.5,1.0,10);
               fitstatus = fit_dscb(ratio,1.5,1.0,10,"PF");
               TF1* fitfnc = (TF1*) ratio->GetListOfFunctions()->Last();
               if (0!=fitfnc && 0==fitstatus) fitfnc->ResetBit(TF1::kNotDraw);
            }
            else if(fittype.CompareTo("CBF")==0) {
               //My implementation
               TF1* fit = new TF1("fit",CrystalBall,0,500,5);
               fit->SetParameters(1.0,1.0,75.0,25.0,1.0);
               ratio->Fit(fit);
               fit->Draw("same");
            }
            else if(fittype.CompareTo("Kevin")==0) {
               //Kevin Pedro's implementation
            
               //get values from histo
               Double_t m = ratio->GetMean();
               Double_t me = ratio->GetMeanError();
               //Double_t m = ratio->GetBinCenter(ratio->GetMaximumBin()); //peak
               Double_t s = ratio->GetRMS();
               Double_t se = ratio->GetRMSError();
               Int_t N = ratio->GetEntries();
   
               std::vector<Double_t> stats(3,0);
               std::vector<Double_t> stat_err(3,0);
               stats[0] = N;
               stat_err[0] = 0;
               stats[1] = m;
               stat_err[1] = me;
               stats[2] = s;
               stat_err[2] = se;

               //find peak
               int nbins = 50;//100;
               double Emin = 5;
               double Emax = 200;
               TSpectrum *spec = new TSpectrum(5);
               if(nbins < 100) spec->Search(ratio,6,"nobackground nodraw goff"); //turn off background removal when nbins too small
               else spec->Search(ratio,6,"nodraw goff");
               Float_t* xpos = spec->GetPositionX();
               Float_t* ypos = spec->GetPositionY();
               Double_t p = xpos[0];
               Double_t ph = ypos[0];
               std::cout << "peak: " << p << std::endl;
               std::cout << "peak height: " << ph << std::endl;
   
               //setup fitting function & do fit
               TF1* gfit = new TF1("DSCBF",cball,Emin,Emax,nCballPar);
               gfit->SetParameters(ph,p,s,1,1.1,1,1.1);
            
               //limits on parameters: 0 < a < 10, 1 < n < 200
               gfit->SetParLimits(3,0,10);
               gfit->SetParLimits(5,0,10);
               gfit->SetParLimits(4,1.01,200);
               gfit->SetParLimits(6,1.01,200);
            
               //formatting
               gfit->SetLineColor(kRed);
               gfit->SetMarkerColor(kRed);
               gfit->SetLineWidth(2);
               //fit
               //ratio->Fit(gfit,"LNQRB");
               ratio->Fit(gfit,"LQRB");

               //TLegend* leg = new TLegend(xmin,0.78,xmin+0.2,0.88);
               leg->AddEntry(ratio,"Ratio");
               leg->AddEntry(gfit,"Fit");
               leg->SetFillColor(0);
               leg->SetBorderSize(0);
               leg->SetTextSize(0.05);
               leg->SetTextFont(42);
               leg->Draw("same");
         
               c->Update();
         
               //left line
               Double_t bndL = gfit->GetParameter(1) - gfit->GetParameter(2)*gfit->GetParameter(3);
               TLine* aLline = new TLine(bndL,pad->GetUymin(),bndL,pad->GetUymax());
               aLline->SetLineStyle(2);
               aLline->SetLineWidth(3);
               aLline->SetLineColor(kBlue);
               aLline->Draw("same");
         
               //left gaussian
               TF1* gsnL = new TF1("gsn","gaus",Emin,bndL);
               gsnL->SetParameters(gfit->GetParameter(0),gfit->GetParameter(1),gfit->GetParameter(2));
               gsnL->SetLineColor(kRed);
               gsnL->SetMarkerColor(kRed);
               gsnL->SetLineWidth(2);
               gsnL->SetLineStyle(2);
               gsnL->Draw("same");

               //line
               Double_t bndR = gfit->GetParameter(1) + gfit->GetParameter(2)*gfit->GetParameter(5);
               TLine* aRline = new TLine(bndR,pad->GetUymin(),bndR,pad->GetUymax());
               aRline->SetLineStyle(2);
               aRline->SetLineWidth(3);
               aRline->SetLineColor(kBlue);
               aRline->Draw("same");
         
               //right gaussian
               TF1* gsnR = new TF1("gsn","gaus",bndR,Emax);
               gsnR->SetParameters(gfit->GetParameter(0),gfit->GetParameter(1),gfit->GetParameter(2));
               gsnR->SetLineColor(kRed);
               gsnR->SetMarkerColor(kRed);
               gsnR->SetLineWidth(2);
               gsnR->SetLineStyle(2);
               gsnR->Draw("same");
            }

            ofile->cd();
            c->Write();
            ratio->Write();
         }
         ifile->Close();
      }
   }

   ofile->Close();

   return 0;
}

////////////////////////////////////////////////////////////////////////////////
// implement local functions
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
int fit_dscb(TH1F*& hrsp,
             const double nsigma,
             const double jtptmin,
             const int niter,
             const string alg)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_dscb()"<<endl;return -1;
  }

  // first use a gaussian to constrain crystal ball gaussian core

  fit_gaussian(hrsp,nsigma,jtptmin,niter);

  TF1* fgaus = hrsp->GetFunction("fgaus");
  if (0==fgaus) {
    hrsp->GetListOfFunctions()->Delete();
    return -1;
  }

  // implementation of the low pt bias threshold

  string histname = hrsp->GetName();
  double ptRefMax(1.0),rspMax(0.0);

  int pos1     = histname.find("RefPt");
  int pos2     = histname.find("to",pos1);
  string ss    = histname.substr(pos1+5,pos2);
  if (from_string(ptRefMax,ss,std::dec)) {
    if (histname.find("RelRsp")==0)
      rspMax = jtptmin/ptRefMax;
    if (histname.find("AbsRsp")==0)
      rspMax = jtptmin-ptRefMax;
  }

  double fitrange_min(0.0);
  if (alg.find("pf")!=string::npos) fitrange_min = std::max(rspMax,0.3);
  else if (alg.find("PF")!=string::npos) fitrange_min = std::max(rspMax,0.3);
  else fitrange_min = std::max(rspMax,0.2);
  double fitrange_max = 1.7;

  adjust_fitrange(hrsp,fitrange_min,fitrange_max);
  //guesstimate_fitrange(hrsp,fitrange_min,fitrange_max,alg);

  
  TF1* fdscb = new TF1("fdscb",fnc_dscb,fitrange_min,fitrange_max,7);

  // set the std values

  double norm = fgaus->GetParameter(0);
  double mean = fgaus->GetParameter(1);
  double sigma= fgaus->GetParameter(2);

  double aone(2.0),atwo(2.0),pone(5.0),ptwo(5.0);
  TVirtualFitter::SetDefaultFitter("Minuit2");

  int fitstatus(0); 
  for (unsigned i=0;i<4;i++) {

    fdscb->SetParameter(0,norm); // N
    fdscb->SetParameter(1,mean); // mean
    fdscb->SetParameter(2,sigma);// sigma
    fdscb->SetParameter(3,aone); // a1
    fdscb->SetParameter(4,pone); // p1
    fdscb->SetParameter(5,atwo); // a2
    fdscb->SetParameter(6,ptwo); // p2

    fdscb->FixParameter(1,mean);
    fdscb->FixParameter(2,sigma);
    
    if (i>0) fdscb->FixParameter(3,aone);
    else fdscb->SetParLimits(3,1.,5.);

    if (i>1) fdscb->FixParameter(5,atwo);
    else fdscb->SetParLimits(5,1.,5.);

    fdscb->SetParLimits(4,0.,25.);
    fdscb->SetParLimits(6,0.,25.);

    fitstatus = hrsp->Fit(fdscb,"RQB0+");
    
    if (0==fitstatus) i=999;

    delete fdscb;
    fdscb = hrsp->GetFunction("fdscb");

    if (0==fdscb) return -1;
      
    norm  = fdscb->GetParameter(0);
    aone  = fdscb->GetParameter(3);
    pone  = fdscb->GetParameter(4);
    atwo  = fdscb->GetParameter(5);
    ptwo  = fdscb->GetParameter(6);

  }

  // reset sigma and mean to gauss values...
  fdscb->SetParameter(1,fgaus->GetParameter(1));
  fdscb->SetParError(1,fgaus->GetParError(1));
  fdscb->SetParameter(2,fgaus->GetParameter(2));
  fdscb->SetParError(2,fgaus->GetParError(2));


  if (0!=fitstatus){
    cout<<"fit_fdscb() to "<<hrsp->GetName()
      <<" failed. Fitstatus: "<<fitstatus<<endl;
    hrsp->GetFunction("fdscb")->Delete();
  }

  return fitstatus;
}

//______________________________________________________________________________
double fnc_dscb(double*xx,double*pp)
{
  double x   = xx[0];
  // gaussian core
  double N   = pp[0];//norm
  double mu  = pp[1];//mean
  double sig = pp[2];//variance
  // transition parameters
  double a1  = pp[3];
  double p1  = pp[4];
  double a2  = pp[5];
  double p2  = pp[6];
  
  double u   = (x-mu)/sig;
  double A1  = TMath::Power(p1/TMath::Abs(a1),p1)*TMath::Exp(-a1*a1/2);
  double A2  = TMath::Power(p2/TMath::Abs(a2),p2)*TMath::Exp(-a2*a2/2);
  double B1  = p1/TMath::Abs(a1) - TMath::Abs(a1);
  double B2  = p2/TMath::Abs(a2) - TMath::Abs(a2);

  double result(N);
  if      (u<-a1) result *= A1*TMath::Power(B1-u,-p1);
  else if (u<a2)  result *= TMath::Exp(-u*u/2);
  else            result *= A2*TMath::Power(B2+u,-p2);
  return result;
}


//______________________________________________________________________________
void fit_gaussian(TH1F*& hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_gaussian()"<<endl;return;
  }
  
  string histname = hrsp->GetName();
  double mean     = hrsp->GetMean();
  double rms      = hrsp->GetRMS();
  double ptRefMax(1.0),rspMax(0.0);     

  double norm  = hrsp->GetMaximumStored();
  double peak  = mean;
  double sigma = rms;
  int pos1     = histname.find("RefPt");
  int pos2     = histname.find("to",pos1);
  string ss    = histname.substr(pos1+5,pos2);
  if (from_string(ptRefMax,ss,std::dec)) {
    if (histname.find("RelRsp")==0)
      rspMax = jtptmin/ptRefMax;
    if (histname.find("AbsRsp")==0)
      rspMax = jtptmin-ptRefMax;
  }
  double xmin  = hrsp->GetXaxis()->GetXmin();
  double xmax  = hrsp->GetXaxis()->GetXmax();
  TF1* fitfnc(0); int fitstatus(-1);
  for (int iiter=0;iiter<niter;iiter++) {
    vector<double> vv;
    vv.push_back(rspMax);
    vv.push_back(xmin);
    vv.push_back(peak-nsigma*sigma);   
    double fitrange_min = *std::max_element(vv.begin(),vv.end());
    double fitrange_max = std::min(xmax,peak+nsigma*sigma);
    adjust_fitrange(hrsp,fitrange_min,fitrange_max);
    fitfnc = new TF1("fgaus","gaus",fitrange_min,fitrange_max);
    fitfnc->SetParNames("N","#mu","#sigma");
    fitfnc->SetParameter(0,norm);
    fitfnc->SetParameter(1,peak);
    fitfnc->SetParameter(2,sigma);
    fitstatus = hrsp->Fit(fitfnc,"RQ0");
    delete fitfnc;
    fitfnc = hrsp->GetFunction("fgaus");
    //fitfnc->ResetBit(TF1::kNotDraw);
    if (fitfnc) {
       norm  = fitfnc->GetParameter(0);
       peak  = fitfnc->GetParameter(1);
       sigma = fitfnc->GetParameter(2);
    }
  }
  if(hrsp->GetFunction("fgaus")==0)
    {
      cout << "No function recorded in histogram " << hrsp->GetName() << endl;
    }
  if (0!=fitstatus){
    cout<<"fit_gaussian() to "<<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus
        <<" - FNC deleted."<<endl;
    hrsp->GetListOfFunctions()->Delete();
  }
}


//______________________________________________________________________________
bool contains(const vector<string>& collection,const string& element)
{
  vector<string>::const_iterator it;
  for (it=collection.begin();it!=collection.end();++it)
    if ((*it)==element) return true;
  return false;
}


//______________________________________________________________________________
void adjust_fitrange(TH1* h,double& min,double& max)
{
  int imin=1; while (h->GetBinLowEdge(imin)<min) imin++;
  int imax=1; while (h->GetBinLowEdge(imax)<max) imax++;
  while ((imax-imin)<8) {
    if (imin>1) {imin--; min = h->GetBinCenter(imin); }
    if (imax<h->GetNbinsX()-1) { imax++; max=h->GetBinCenter(imax); }
  }
}

//______________________________________________________________________________
void guesstimate_fitrange(TH1* h,double& min,double& max, const string alg)
{

  //hh: tried a variety of different possibilities here

  const double norm = h->GetEffectiveEntries();
  if (norm<=0.0) return;
  double lsum(0.0); double hsum(0.0);
  int nmax = h->GetNbinsX();
  if (nmax<3) return;
  int imin=1; while (imin<nmax) {
    lsum+=h->GetBinContent(imin);
    if (alg.find("pf")!=string::npos) {
      if (lsum/norm>.0005) break;
    }
    else {
      if (lsum/norm>.00005) break;
    }
    imin++;
  } 
  int imax=nmax-1; while (imax>0) {
    hsum+=h->GetBinContent(imax);
    if (lsum/norm>.000001) break;
    imin--;
  }
  min = h->GetBinLowEdge(imin);
  max = h->GetBinLowEdge(imax);
}

//______________________________________________________________________________
double CrystalBall (double* x, double* par){ 
  //http://en.wikipedia.org/wiki/Crystal_Ball_function 
  double xcur = x[0]; 
  double alpha = par[0]; 
  double n = par[1]; 
  double mu = par[2]; 
  double sigma = par[3]; 
  double N = par[4]; 
  TF1* exp = new TF1("exp","exp(x)",1e-20,1e20); 
  double A; double B; 
  if (alpha < 0){ 
    A = pow((n/(-1*alpha)),n)*exp->Eval((-1)*alpha*alpha/2); 
    B = n/(-1*alpha) + alpha;} 
  else { 
    A = pow((n/alpha),n)*exp->Eval((-1)*alpha*alpha/2); 
    B = n/alpha - alpha;}

  double f;
  if ((xcur-mu)/sigma > (-1)*alpha)  //original
    f = N*exp->Eval((-1)*(xcur-mu)*(xcur-mu)/(2*sigma*sigma)); 
  else 
    f = N*A*pow((B- (xcur-mu)/sigma),(-1*n)); //original 
  delete exp;

  return f; 
}  


//------------------------------------------
//Double-sided Crystal Ball function
//parameters:
//N, mu, sigma, aL, nL, aR, nR, C
//0,  1,     2,  3,  4,  5,  6, 7
Double_t cball(Double_t *x, Double_t *par){
   //ensure sigma > 0 and a > 0
   //Double_t N = 1/(sigma*(n/a*1/(n-1)*Exp(-a*a/2) + Sqrt(Pi()/2)*(1+Erf(a/Sqrt(2))))); //normalized N
   Double_t N = par[0]; //let N float
   Double_t mu = par[1];
   par[2] = Abs(par[2]);
   Double_t sigma = par[2];
   par[3] = Abs(par[3]);
   Double_t aL = par[3];
   par[4] = (par[4]>1) ? par[4] : 1.01; //n>1 required
   Double_t nL = par[4];
   par[5] = Abs(par[5]);
   Double_t aR = par[5];
   par[6] = (par[6]>1) ? par[6] : 1.01; //n>1 required
   Double_t nR = par[6];   
   Double_t arg = (x[0]-mu)/sigma;
   Double_t C = par[7]; //let C float
   
   //left tail
   if(arg <= -aL){
      return N*Power(nL/aL,nL)*Exp(-Power(aL,2)/2)*Power(nL/aL-aL-arg,-nL)+C;
   }
   //right tail
   else if(arg >= aR){
      return N*Power(nR/aR,nR)*Exp(-Power(aR,2)/2)*Power(nR/aR-aR+arg,-nR)+C;
   }
   //core
   else{
      return N*Exp(-Power(arg,2)/2)+C;
   }

}


