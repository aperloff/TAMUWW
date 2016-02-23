#include "tmvaglob.C"
#include "Math/Polynomial.h"
#include "Math/Interpolator.h"

void plot_ROC_FOM( TFile* file, Int_t type = 2,
                   TDirectory* BinDir, TString method )
{
  // input:   - Input file (result from TMVA),
  //          - type = 1 --> plot efficiency(B) versus eff(S)
  //                 = 2 --> plot rejection (B) versus efficiency (S)

  Bool_t __PLOT_LOGO__  = kTRUE;
  Bool_t __SAVE_IMAGE__ = kTRUE;

  // the coordinates
  Float_t x1 = 0;
  Float_t x2 = 1;
  Float_t y1 = 0;
  Float_t y2 = 0.8;

  // reverse order if "rejection"
  if (type == 2) {
    Float_t z = y1;
    y1 = 1 - y2;
    y2 = 1 - z;
    //      cout << "--- type==2: plot background rejection versus signal efficiency" << endl;
  }
  else {
    //  cout << "--- type==1: plot background efficiency versus signal efficiency" << endl;
  }
  // create canvas
  TCanvas* c = new TCanvas( "c", "the canvas", 200, 0, 650, 500 );

  // global style settings
  c->SetGrid();
  c->SetTicks();

  // legend
  Float_t x0L = 0.107,     y0H = 0.899;
  //Float_t dxL = 0.457-x0L, dyH = 0.22;
  Float_t dxL = 0.6-x0L, dyH = 0.66;
  if (type == 2) {
    x0L = 0.15;
    y0H = 1 - y0H + dyH + 0.07;
  }
  TLegend *legend = new TLegend( x0L, y0H-dyH, x0L+dxL, y0H );
  //legend->SetTextSize( 0.05 );
  legend->SetHeader( "MVA Method:" );
  legend->SetMargin( 0.4 );

  TString xtit = "Signal efficiency";
  TString ytit = "Background efficiency";
  if (type == 2) ytit = "Background rejection";
  TString ftit = ytit + " versus " + xtit;

  if (TString(BinDir->GetName()).Contains("multicut")){
    ftit += "  Bin: ";
    ftit += (BinDir->GetTitle());
  }

  // draw empty frame
  if(gROOT->FindObject("frame")!=0) gROOT->FindObject("frame")->Delete();
  TH2F* frame = new TH2F( "frame", ftit, 500, x1, x2, 500, y1, y2 );
  frame->GetXaxis()->SetTitle( xtit );
  frame->GetYaxis()->SetTitle( ytit );
  TMVAGlob::SetFrameStyle( frame, 1.0 );

  frame->Draw();

  Int_t color = 1;
  Int_t nmva  = 0;
  TKey *key, *hkey;

  TString hNameRef = Form("effBvsS_%s_simple",method.Data());
  if (type == 2) hNameRef = Form("rejBvsS_%s_simple",method.Data());

  // begin section from Joey
  Int_t nbins = 1000;
  Float_t xmin = -1.0, xmax = 1.0;
  TTree* tmvaTestTree = (TTree*)gDirectory->Get("TestTree");
  TTree* tmvaTrainTree = (TTree*)gDirectory->Get("TrainTree");
  TH1F * TMVA_BDT_bkg = new TH1F(Form("TMVA_%s_bkg",method.Data()),
                                 Form("Bkg %s ",method.Data()),nbins, xmin, xmax);
  TH1F * TMVA_BDT_sig = new TH1F(Form("TMVA_%s_sig",method.Data()),
                                 Form("Sig %s ",method.Data()),nbins, xmin, xmax);
  Float_t BDT_test, BDT_train;
  Int_t classID_test, classID_train;
  tmvaTestTree->SetBranchAddress(method,&BDT_test);
  tmvaTrainTree->SetBranchAddress(method,&BDT_train);
  tmvaTestTree->SetBranchAddress("classID",&classID_test);
  tmvaTrainTree->SetBranchAddress("classID",&classID_train);

  for(int a = 0; a < tmvaTestTree->GetEntries();a++){
    tmvaTestTree->GetEntry(a);
    if(classID_test == 1) TMVA_BDT_sig->Fill(BDT_test);
    if(classID_test == 0) TMVA_BDT_bkg->Fill(BDT_test);
    if(classID_test != 0 && classID_test != 1) cout << "We have a problem in ID Classifier (classID_test = " << classID_test << ")" << endl;
  }
  for(int a = 0; a < tmvaTrainTree->GetEntries();a++){
    tmvaTrainTree->GetEntry(a);
    if(classID_train == 1) TMVA_BDT_sig->Fill(BDT_train);
    if(classID_train == 0) TMVA_BDT_bkg->Fill(BDT_train);
    if(classID_train != 0 && classID_train != 1) cout << "We have a problem in ID Classifier (classID_train = " << classID_train << ")" << endl;
  }

  TGraph * rocBDT = new TGraph();
  rocBDT->SetNameTitle(Form("roc%s",method.Data()),ftit);
  TGraph * rocBDT_int = new TGraph();
  rocBDT_int->SetNameTitle(Form("roc%s_int",method.Data()),ftit+" (Interpolated Values)");
  TGraph * bestLine = new TGraph();
  bestLine->SetNameTitle("bestLineGraph","bestLineGraph");

  vector<double> revX, revY;
  Float_t backRej_bdt = 0, sigEff_bdt = 0, minDist = 9999, tempDist=9999;
  Int_t minPoint = 0;
  double minPointX = 0.0, minPointY = 0.0;
  for(int x=0; x < nbins; x++ ){
  //Calculate and fill ROC curve for BDT                                                                                                            
    backRej_bdt = (TMVA_BDT_bkg->Integral(0,x))/TMVA_BDT_bkg->Integral();
    sigEff_bdt = (TMVA_BDT_sig->Integral() - TMVA_BDT_sig->Integral(0,x))/TMVA_BDT_sig->Integral();
    rocBDT->SetPoint(x,sigEff_bdt,backRej_bdt);

    if(revX.size()==0 || revX.front()!=sigEff_bdt) { //|| revY.front()!=backRej_bdt) {
      revX.insert(revX.begin(),sigEff_bdt);
      revY.insert(revY.begin(),backRej_bdt);
    }

    tempDist = TMath::Sqrt(TMath::Power(1-sigEff_bdt,2)+TMath::Power(1-backRej_bdt,2));
    if (tempDist<minDist) {
      minDist = tempDist;
      minPoint = x;
      minPointX = sigEff_bdt;
      minPointY = backRej_bdt;
    }
  }
  // end section from Joey

  for(unsigned int j=0; j<revX.size(); j++) {
    cout << "j = " << j << "\t(" << revX[j] << "," << revY[j] << ")" << endl;
    if(j!=revX.size()-1 && revX[j]>revX[j+1])
      cout << "WARNING::ROC_FOM x values not monotonically increasing!" << endl;
    //if(revX[j] == revX[j+1] && revY[j] == revY[j+1] && j!=revX.size()-1) {
    //  revX.erase(revX.begin()+j+1);
    //  revY.erase(revY.begin()+j+1);
    //  j--;
    //} 
  }

  ROOT::Math::Interpolator inter(revX,revY, ROOT::Math::Interpolation::kCSPLINE); //kCSPLINE,kPOLYNOMIAL
  inter.SetData(revX, revY);
  int nIntPoints = 100;
  double intStep = 1.0/nIntPoints;
  for(int i=0; i<nIntPoints; i++) {
    rocBDT_int->SetPoint(i,intStep*i,inter.Eval(intStep*i));
  }

  //Calculate the integral under the curve using the trapezoid rule
  double area_under_curve = 0.0;
  for(int i=0; i<rocBDT->GetN()-1; i++) {
    double diffX = rocBDT->GetX()[i]-rocBDT->GetX()[i+1];
    double sumY = rocBDT->GetY()[i]+rocBDT->GetY()[i+1];
    area_under_curve += (diffX*sumY);
  }
  area_under_curve/=2.0;

  bestLine->SetPoint(0,minPointX,minPointY);
  bestLine->SetPoint(1,1,1);

  rocBDT->Draw("APC");
  rocBDT->SetLineWidth(3);
  rocBDT->GetXaxis()->SetTitle(xtit);
  rocBDT->GetYaxis()->SetTitle(ytit);
  rocBDT_int->Draw("SAME P");
  rocBDT_int->SetMarkerColor(kRed);
  rocBDT_int->SetMarkerStyle(6);
  bestLine->Draw("SAME LP");
  bestLine->SetLineColor(kBlue);
  bestLine->SetMarkerStyle(20);
  bestLine->SetMarkerSize(1);
  legend->AddEntry(rocBDT,method,"lp");
  legend->AddEntry(rocBDT_int,Form("%s (Interpolated Points)",method.Data()),"p");
  legend->AddEntry(bestLine,Form("Distance to (1,1) = %.4f",minDist),"l");
  legend->AddEntry((TObject*)0,Form("Best Point = (%.4f,%.4f)",minPointX,minPointY),"");
  legend->AddEntry((TObject*)0,Form("Area Under Curve = %.4f",area_under_curve),"");

  // rescale legend box size
  // current box size has been tuned for 3 MVAs + 1 title
  if (type == 1) {
    dyH *= (1.0 + Float_t(nmva - 3.0)/4.0);
    legend->SetY1( y0H - dyH );
  }
  else {
    dyH *= (Float_t(TMath::Min(10,nmva) - 3.0)/4.0);
    legend->SetY2( y0H + dyH);
  }

  // redraw axes
  frame->Draw("sameaxis");
  legend->Draw("same");

  // ============================================================

  if (__PLOT_LOGO__) TMVAGlob::plot_logo();

  // ============================================================

  c->Update();

  TString fname = "plots/" + hNameRef;
  if (__SAVE_IMAGE__) TMVAGlob::imgconv( c, fname );

  return;
}

void ROC_FOM( TString fin = "TMVA.root", TString method = "BDT",
              Int_t type = 2, Bool_t useTMVAStyle = kTRUE )
{
  // argument: type = 1 --> plot efficiency(B) versus eff(S)
  //           type = 2 --> plot rejection (B) versus efficiency (S)

  // set style and remove existing canvas'
  TMVAGlob::Initialize( useTMVAStyle );

  // checks if file with name "fin" is already open, and if not opens one
  TFile* file = TMVAGlob::OpenFile( fin );

  plot_ROC_FOM( file, type, gDirectory, method );

  return;
}
