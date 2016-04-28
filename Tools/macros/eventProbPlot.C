void simpleEPDPlot() {
   TCanvas* c = new TCanvas("c","c",800,800);
   c->cd()->SetGrid();
   TLegend* l = new TLegend(0.50,0.70,0.90,0.90);
   TFile* fWJets = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWJets_BaseCuts.root","READ");
   METree->SetLineColor(kBlue);
   drawCustomEPD(l,"WJets","");
   TFile* fWH125 = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWH125_BaseCuts.root","READ");
   METree->SetLineColor(kRed);
   drawCustomEPD(l,"WH125","same");
   TFile* fggH125 = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microggH125_BIG_BaseCuts.root","READ");
   METree->SetLineColor(kGreen);
   drawCustomEPD(l,"ggH125","same");

   l->Draw("same");

   c->SaveAs("simpleEPDPlot.png");
   c->SaveAs("simpleEPDPlot.eps");
}
/*
 EventProbs | MatrixElementType | PhysicsProcessType |    tmeParam |                  MaxEventProb |                 MeanEventProb |    MedianEventProb
          0 |                WW |                 WW |     0 +/- 0 |  7.819953e-07 +/- 0.006416942 | 4.478285e-09 +/- 9.720642e-06 | 2.391546e-10 +/- 0
          1 |                WZ |                 WZ |     0 +/- 0 |  7.640095e-08 +/- 0.006219349 | 2.967587e-10 +/- 8.629397e-06 | 4.906673e-12 +/- 0
          2 |              WZbb |                 WZ |     0 +/- 0 |  7.743198e-08 +/- 0.009021043 | 3.394647e-10 +/- 1.247041e-05 | 6.567607e-12 +/- 0
          3 |               WLg |              WJets |     0 +/- 0 | 0.0001938386 +/- 0.0005559995 | 1.124578e-06 +/- 9.116343e-06 | 1.300328e-07 +/- 0
          4 |               WLg |              WJets |     0 +/- 0 |  4.231735e-05 +/- 0.001675434 | 3.522228e-07 +/- 8.794365e-06 | 4.422508e-08 +/- 0
          5 |               Wgg |              WJets |     0 +/- 0 |  0.0001268688 +/- 0.001890234 | 5.372542e-07 +/- 7.983288e-06 | 3.309296e-08 +/- 0
          6 |               WLL |              WJets |     0 +/- 0 | 5.851085e-06 +/- 0.0005302307 | 4.278642e-08 +/- 7.532599e-06 | 6.100224e-09 +/- 0
          7 |               WLb |              WJets | 172.5 +/- 0 |  9.279982e-07 +/- 0.009765282 | 4.570092e-09 +/- 1.006993e-05 | 4.925576e-10 +/- 0
          8 |               Wbb |              WJets |     0 +/- 0 |  3.737306e-07 +/- 0.005721496 | 3.412428e-09 +/- 7.677563e-06 | 3.196485e-10 +/- 0
          9 |            ZLight |              ZJets |     0 +/- 0 |   6.765438e-06 +/- 0.01249524 | 7.571615e-08 +/- 5.002138e-05 | 2.014442e-08 +/- 0
         10 |             STopT |            STopT_T | 172.5 +/- 0 |  1.080755e-06 +/- 0.007142552 | 5.561144e-09 +/- 2.112516e-05 | 4.289709e-10 +/- 0
         10 |             STopT |         STopT_Tbar | 172.5 +/- 0 |  3.863697e-07 +/- 0.009917325 | 2.947556e-09 +/- 2.887322e-05 | 2.239826e-10 +/- 0
         11 |             STopS |            STopS_T | 172.5 +/- 0 |    1.5464e-08 +/- 0.009877172 | 3.346865e-10 +/- 7.010831e-05 | 4.970496e-11 +/- 0
         11 |             STopS |         STopS_Tbar | 172.5 +/- 0 |  1.340438e-08 +/- 0.008848001 | 1.753939e-10 +/- 9.402285e-05 | 3.076242e-11 +/- 0
         12 |               QCD |         QCD_ElFULL |     0 +/- 0 |     0.3190735 +/- 0.006974173 |  0.001470852 +/- 8.683943e-06 |  1.04846e-05 +/- 0
         12 |               QCD |         QCD_MuFULL |     0 +/- 0 |      1.054128 +/- 0.007607975 |  0.003589481 +/- 2.084871e-05 | 1.244808e-05 +/- 0
         19 |               HWW |             ggH125 |   125 +/- 0 |   6.625523e-10 +/- 0.00913204 | 5.296594e-12 +/- 9.963569e-05 | 3.016789e-13 +/- 0
         54 |                WH |              WH125 |   125 +/- 0 |  2.085422e-10 +/- 0.006971693 | 3.676811e-13 +/- 2.859221e-05 | 1.481892e-15 +/- 0

*/
void drawCustomEPD(TLegend* l, TString title, TString option) {
   TH1F *htemp;
   //Normalized to the mean
   //Zeros for tchan2, Wc, TTbar
   //TString epdString = "((eventProb[19]/5.296594e-12) + (eventProb[54]/3.676811e-13)) / ((eventProb[19]/5.296594e-12) + (eventProb[54]/3.676811e-13) + (eventProb[11]/3.346865e-10) + (eventProb[10]/5.561144e-09) + 0 + (eventProb[8]/3.412428e-09) + 0 + (eventProb[12]/0.001470852) + 0 + (eventProb[0]/4.478285e-09) + (eventProb[1]/2.967587e-10))";
   //double meNormMean[6] = {1.029e-6, 3.235e-7, 5.084e-7, 3.929e-8,5.043e-12,3.64e-13};
   //TString epdString = "((eventProb[19]/5.043e-12) + (eventProb[54]/3.64e-13)) / ((eventProb[19]/5.043e-12) + (eventProb[54]/3.64e-13) + (eventProb[3]/1.029e-6) + (eventProb[4]/3.235e-7) + (eventProb[5]/5.084e-7) + (eventProb[6]/3.929e-8))";
   //Normalized to the median
   //TString epdString = "((eventProb[19]/3.001251e-13) + (eventProb[54]/1.481892e-15)) / ((eventProb[19]/3.001251e-13) + (eventProb[54]/1.481892e-15) + (eventProb[3]/1.300328e-07) + (eventProb[4]/4.422508e-08) + (eventProb[5]/3.309222e-08) + (eventProb[6]/6.101245e-09))";
   //Same MEs as Saideth and Andrea. HWW125, WH125, WLg, WLg sub, Wgg, WLL
   ////TString epdString = "((eventProb[19]/5.296594e-12) + (eventProb[54]/3.676811e-13)) / ((eventProb[19]/5.296594e-12) + (eventProb[54]/3.676811e-13) + (eventProb[3]/1.124578e-06) + (eventProb[4]/3.522228e-07) + (eventProb[5]/5.372542e-07) + (eventProb[6]/4.278642e-08))";
   
   //Ricardos EPD test as of 05/29/2014
   //(num=WH, ggH, QCD, WLL, Wgg, WLg_sub, WLg den=same+WW+WZ) (normalize epd by mean) (do on WJets sample an on Higgs sample)
   TString num = "(eventProb[19]/5.296594e-12) + (eventProb[54]/3.676811e-13) + (eventProb[12]/0.001470852) + (eventProb[6]/4.278642e-08) + (eventProb[5]/5.372542e-07) + (eventProb[3]/1.124578e-06) + (eventProb[4]/3.522228e-07)";
   TString den = num + " + (eventProb[0]/4.478285e-09) + (eventProb[1]/2.967587e-10)";
   TString epdString = Form("(%s)/(%s)",num.Data(),den.Data());

   METree->Draw(epdString+">>"+title,"Entry$<40000",option);
   htemp = (TH1F*)gPad->GetPrimitive(title);
   //htemp->GetYaxis()->SetRangeUser(0,0.2);
   htemp->GetXaxis()->SetTitleSize(0.04);
   htemp->GetXaxis()->SetTitleOffset(1.5);
   htemp->GetXaxis()->SetTitle("EPD");
   htemp->GetYaxis()->SetTitle("# of Evts");
   l->AddEntry(htemp,title,"l");
}
/*************************************************************************************************************/
void epdPlot() {
   TCanvas* c = new TCanvas("c","c",800,800);
   c->cd()->SetLogx();
   TLegend* l = new TLegend(0.19,0.47,0.43,0.92);
   TPaveText *pt = new TPaveText(.05,.96,.95,0.99,"brNDC");
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetBorderSize(0);
   pt->AddText("epdPretagHiggs[6] = Optimized EPD for ggH(125)+qqH(125)+WH(125)");

   TFile* fQCD         = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microQCD_Electron_Dp6p7_BaseCuts.root","READ");
   METree->SetLineColor(kBlack);
   drawEPDHiggs(l,"QCD","norm");
   TFile* fSTopS_Tbar = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopS_Tbar_BaseCuts.root","READ");
   METree->SetLineColor(kRed);
   drawEPDHiggs(l,"STopS_Tbar","same norm");
   TFile* fSTopTW_T = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopTW_T_BaseCuts.root","READ");
   METree->SetLineColor(kOrange);
   drawEPDHiggs(l,"STopTW_T","same norm");
   TFile* fSTopTW_Tbar = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopTW_Tbar_BaseCuts.root","READ");
   METree->SetLineColor(kGreen);
   drawEPDHiggs(l,"STopTW_Tbar","same norm");
   TFile* fSTopT_T = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopT_T_BaseCuts.root","READ");
   METree->SetLineColor(kBlue);
   drawEPDHiggs(l,"STopT_T","same norm");
   TFile* fSTopT_Tbar = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopT_Tbar_BaseCuts.root","READ");
   METree->SetLineColor(kViolet);
   drawEPDHiggs(l,"STopT_Tbar","same norm");
   TFile* fTTbar = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microTTbar_BaseCuts.root","READ");
   METree->SetLineColor(kMagenta);
   drawEPDHiggs(l,"TTbar","same norm");
   TFile* fWH125 = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWH125_BaseCuts.root","READ");
   METree->SetLineColor(kCyan);
   drawEPDHiggs(l,"WH125","same norm");
   TFile* fWJets = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWJets_BaseCuts.root","READ");
   METree->SetLineColor(36);
   drawEPDHiggs(l,"WJets","same norm");
   TFile* fWW = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWW_BaseCuts.root","READ");
   METree->SetLineColor(30);
   drawEPDHiggs(l,"WW","same norm");
   TFile* fWZ = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWZ_BaseCuts.root","READ");
   METree->SetLineColor(46);
   drawEPDHiggs(l,"WZ","same norm");
   TFile* fZJets = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microZJets_BaseCuts.root","READ");
   METree->SetLineColor(42);
   drawEPDHiggs(l,"ZJets","same norm");
   TFile* fggH125 = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microggH125_BaseCuts.root","READ");
   METree->SetLineColor(40);
   drawEPDHiggs(l,"ggH125","same norm");
   TFile* fqqH125 = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microqqH125_BaseCuts.root","READ");
   METree->SetLineColor(41);
   drawEPDHiggs(l,"qqH125","same norm");

   l->Draw("same");
   pt->Draw("same");

   c->SaveAs("OptimizedEPDs.eps");
   c->SaveAs("OptimizedEPDs.png");
}
void drawEPDHiggs(TLegend* l, TString title, TString option) {
   TH1F *htemp;
   METree->Draw("-TMath::Log10(epdPretagHiggs[6])>>epdHiggs_"+title,"",option);
   htemp = (TH1F*)gPad->GetPrimitive("epdHiggs_"+title);
   htemp->GetYaxis()->SetRangeUser(0,0.2);
   htemp->GetXaxis()->SetTitleSize(0.04);
   htemp->GetXaxis()->SetTitleOffset(1.5);
   htemp->GetXaxis()->SetTitle("-log_{10}(epdPretagHiggs(125))");
   l->AddEntry(htemp,title,"l");
}
/*************************************************************************************************************/
void eventProbPlot() {
   TCanvas* c = new TCanvas("c","c",1500,900);
   c->Divide(5,3);

   TFile* fQCD         = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microQCD_Electron_Dp6p7_BaseCuts.root","READ");
   c->cd(1);
   doThis("eventProbs for QCD Sample");
   c->Update();
   TFile* fSTopS_Tbar  = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopS_Tbar_BaseCuts.root","READ");
   c->cd(2);
   doThis("eventProbs for STopS_Tbar Sample");
   c->Update();
   TFile* fSTopTW_T    = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopTW_T_BaseCuts.root","READ");
   c->cd(3);
   doThis("eventProbs for STopTW_T Sample");
   c->Update();
   TFile* fSTopTW_Tbar = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopTW_Tbar_BaseCuts.root","READ");
   c->cd(4);
   doThis("eventProbs for STopTW_Tbar Sample");
   c->Update();
   TFile* fSTopT_T     = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopT_T_BaseCuts.root","READ");
   c->cd(5);
   doThis("eventProbs for STopT_T Sample");
   c->Update();
   TFile* fSTopT_Tbar  = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopT_Tbar_BaseCuts.root","READ");
   c->cd(6);
   doThis("eventProbs for STopT_Tbar Sample");
   c->Update();
   TFile* fTTbar       = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microTTbar_BaseCuts.root","READ");
   c->cd(7);
   doThis("eventProbs for TTbar Sample");
   c->Update();
   TFile* fWH125       = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWH125_BaseCuts.root","READ");
   c->cd(8);
   doThis("eventProbs for WH125 Sample");
   c->Update();
   TFile* fWJets       = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWJets_BaseCuts.root","READ");
   c->cd(9);
   doThis("eventProbs for WJets Sample");
   c->Update();
   TFile* fWW          = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWW_BaseCuts.root","READ");
   c->cd(10);
   doThis("eventProbs for WW Sample");
   c->Update();
   TFile* fWZ          = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWZ_BaseCuts.root","READ");
   c->cd(11);
   doThis("eventProbs for WZ Sample");
   c->Update();
   TFile* fZJets       = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microZJets_BaseCuts.root","READ");
   c->cd(12);
   doThis("eventProbs for ZJets Sample");
   c->Update();
   TFile* fggH125      = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microggH125_BaseCuts.root","READ");
   c->cd(13);
   doThis("eventProbs for ggH125 Sample");
   c->Update();
   TFile* fqqH125      = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microqqH125_BaseCuts.root","READ");
   c->cd(14);
   doThis("eventProbs for qqH125 Sample");
   c->Update();

   c->SaveAs("eventProbs.eps");
   c->SaveAs("eventProbs.png");
}
/*************************************************************************************************************/
void doThis(TString title){
   cout << "\tPlotting " << title << " ... " << std::flush;
   gPad->SetTitle(title);
   TPaveText *pt = new TPaveText(.05,.94,.95,0.99,"brNDC");
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetBorderSize(0);
   pt->AddText(title);

   TH1F *htemp;
   TLegend* l = new TLegend(0.2,0.4,0.45,0.85);
   METree->SetLineColor(kBlack);
   METree->Draw("TMath::Log10(eventProb[0])-TMath::Log10(0.781996e-6)>>WW","","norm");
   htemp = (TH1F*)gPad->GetPrimitive("WW");
   htemp->SetTitle(title);
   if(title.Contains("Jets"))
      htemp->GetYaxis()->SetRangeUser(0,0.1);
   else
      htemp->GetYaxis()->SetRangeUser(0,0.07);
   htemp->GetXaxis()->SetTitleSize(0.04);
   htemp->GetXaxis()->SetTitleOffset(1.2);
   htemp->GetXaxis()->SetTitle("log_{10}(eventProb_{i})-log_{10}(eventProb_{i}^{MAX})");
   l->AddEntry(htemp,"WW","l");
   METree->SetLineColor(kRed);
   METree->Draw("TMath::Log10(eventProb[1])-TMath::Log10(76.4008e-9)>>WZ","","same norm");
   htemp = (TH1F*)gPad->GetPrimitive("WZ");
   l->AddEntry(htemp,"WZ","l");
   METree->SetLineColor(kOrange);
   METree->Draw("TMath::Log10(eventProb[3])-TMath::Log10(0.1939e-3)>>WLg","","same norm");
   htemp = (TH1F*)gPad->GetPrimitive("WLg");
   l->AddEntry(htemp,"WLg","l");
   METree->SetLineColor(kGreen);
   METree->Draw("TMath::Log10(eventProb[4])-TMath::Log10(42.31735e-6)>>WLgsub","","same norm");
   htemp = (TH1F*)gPad->GetPrimitive("WLgsub");
   l->AddEntry(htemp,"WLg sub","l");
   METree->SetLineColor(kBlue);
   METree->Draw("TMath::Log10(eventProb[5])-TMath::Log10(0.12688e-3)>>Wgg","","same norm");
   htemp = (TH1F*)gPad->GetPrimitive("Wgg");
   l->AddEntry(htemp,"Wgg","l");
   METree->SetLineColor(kViolet);
   METree->Draw("TMath::Log10(eventProb[6])-TMath::Log10(5.85109e-6)>>WLL","","same norm");
   htemp = (TH1F*)gPad->GetPrimitive("WLL");
   l->AddEntry(htemp,"WLL","l");
   METree->SetLineColor(kMagenta);
   METree->Draw("TMath::Log10(eventProb[12])-TMath::Log10(0.592615)>>QCD","","same norm");
   htemp = (TH1F*)gPad->GetPrimitive("QCD");
   l->AddEntry(htemp,"QCD","l");
   METree->SetLineColor(kCyan);
   METree->Draw("TMath::Log10(eventProb[19])-TMath::Log10(0.66255e-9)>>ggH125","","same norm");
   htemp = (TH1F*)gPad->GetPrimitive("ggH125");
   l->AddEntry(htemp,"ggH125","l");
   METree->SetLineColor(kAzure);
   METree->Draw("TMath::Log10(eventProb[54])-TMath::Log10(0.208542e-9)>>WH125","","same norm");
   htemp = (TH1F*)gPad->GetPrimitive("WH125");
   l->AddEntry(htemp,"WH125","l");
   l->Draw("same");
   pt->Draw("same");

   cout << "DONE" << endl;
}
/*************************************************************************************************************/
TString samples[14] = {"QCD","STopS_Tbar","STopTW_T","STopTW_Tbar","STopT_T","STopT_Tbar","TTbar","WH125","WJets","WW","WZ","ZJets","ggH125","qqH125"};
Color_t colors[14] = {kBlack,kRed,kOrange,kGreen,kBlue,kViolet,kMagenta,kCyan,36,30,46,42,40,41};
TString ep[9] = {"WW","WZ","WLg","WLg sub","Wgg","WLL","QCD","ggH125","WH125"};
TString eps[9] = {"0","1","3","4","5","6","12","19","54"}

void eventProbPlot_2() {
   TFile* ofile = new TFile("eventProbs_2.root","RECREATE");
   TCanvas* c = new TCanvas("c","c",1200,950);
   c->Divide(3,3);

   TFile* fQCD         = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microQCD_Electron_Dp6p7_BaseCuts.root","READ");
   TFile* fSTopS_Tbar  = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopS_Tbar_BaseCuts.root","READ");
   TFile* fSTopTW_T    = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopTW_T_BaseCuts.root","READ");
   TFile* fSTopTW_Tbar = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopTW_Tbar_BaseCuts.root","READ");
   TFile* fSTopT_T     = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopT_T_BaseCuts.root","READ");
   TFile* fSTopT_Tbar  = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microSTopT_Tbar_BaseCuts.root","READ");
   TFile* fTTbar       = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microTTbar_BaseCuts.root","READ");
   TFile* fWH125       = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWH125_BaseCuts.root","READ");
   TFile* fWJets       = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWJets_BaseCuts.root","READ");
   TFile* fWW          = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWW_BaseCuts.root","READ");
   TFile* fWZ          = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microWZ_BaseCuts.root","READ");
   TFile* fZJets       = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microZJets_BaseCuts.root","READ");
   TFile* fggH125      = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microggH125_BaseCuts.root","READ");
   TFile* fqqH125      = new TFile("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/microNtuples_optimized/microqqH125_BaseCuts.root","READ");

   TFile* fArray[14];
   fArray[0] = fQCD;
   fArray[1] = fSTopS_Tbar;
   fArray[2] = fSTopTW_T;
   fArray[3] = fSTopTW_Tbar;
   fArray[4] = fSTopT_T;
   fArray[5] = fSTopT_Tbar;
   fArray[6] = fTTbar;
   fArray[7] = fWH125;
   fArray[8] = fWJets;
   fArray[9] = fWW;
   fArray[10] = fWZ;
   fArray[11] = fZJets;
   fArray[12] = fggH125;
   fArray[13] = fqqH125;

   for(unsigned int iep=0; iep<9; iep++) {
      c->cd(iep+1);
      doThis_2(iep,ep[iep]+" eventProbs",fArray);
      c->Update();
   }

   c->SaveAs("eventProbs_2.eps");
   c->SaveAs("eventProbs_2.pdf");
   c->SaveAs("eventProbs_2.png");

   ofile->cd();
   c->Write();
   ofile->Close();
}

void doThis_2(int iep, TString title, TFile* fArray[14]){
   cout << "\tPlotting " << title << " ... " << std::flush;
   gPad->SetTitle(title);
   TPaveText *pt = new TPaveText(.05,.94,.95,0.99,"brNDC");
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetBorderSize(0);
   pt->AddText(title);

   TH1F *htemp;
   TLegend* l = new TLegend(0.2,0.4,0.45,0.85);
   double min = -30.0;
   double max = 0.0;
   bool minSet = false;
   for(unsigned int isample=0; isample<14; isample++) {
      fArray[isample]->cd();
      METree->SetLineColor(colors[isample]);
      if(isample==0) {
         METree->Draw("TMath::Log10(eventProb["+eps[iep]+"])>>"+samples[isample]+"_"+ep[iep],"","norm");
      }
      else {
         METree->Draw("TMath::Log10(eventProb["+eps[iep]+"])>>"+samples[isample]+"_"+ep[iep],"","same norm");
      }
      htemp = (TH1F*)gPad->GetPrimitive(samples[isample]+"_"+ep[iep]);
      htemp->SetTitle(title);
      htemp->GetYaxis()->SetRangeUser(0,0.14);
      htemp->GetXaxis()->SetTitle("log_{10}(eventProb["+ep[iep]+"])");
      if(isample == 0) {
         for(int i=1; i<=htemp->GetNbinsX(); i++) {
            if(htemp->GetBinContent(i)!=0) {
               if(htemp->GetBinContent(i+1)==0) {
                  max = ceil(htemp->GetBinCenter(i))+1.0;
               }
               if(minSet == false) {
                  min = floor(htemp->GetBinCenter(i))-1.0;
                  minSet = true;
               }
            }
         }
      }
      htemp->GetXaxis()->SetRangeUser(min,max);
      l->AddEntry(htemp,samples[isample],"l");
   }
   l->Draw("same");
   pt->Draw("same");

   cout << "DONE" << endl;
}

/*
***********************************************
*    Row   * Instance * m_tProbSt * m_tProbSt *
***********************************************
*        0 *        0 *       400 *         0 *WW
*        0 *        1 *       401 *         0 *WZ
*        0 *        2 *       403 *         0 *WZbb
*        0 *        3 *       301 *         0 *WLg
*        0 *        4 *       301 *         0 *WLg sub
*        0 *        5 *       302 *         0 *Wgg
*        0 *        6 *       303 *         0 *WLL
*        0 *        7 *       304 *     172.5 *WLb
*        0 *        8 *       305 *         0 *Wbb
*        0 *        9 *       350 *         0 *ZLight
*        0 *       10 *       200 *     172.5 *STopT
*        0 *       11 *       201 *     172.5 *STopS
*        0 *       12 *       500 *         0 *QCD
*        0 *       13 *       150 *       110 *HWW
*        0 *       14 *       150 *       115 *HWW
*        0 *       15 *       150 *       120 *HWW
*        0 *       16 *       150 *       122 *HWW
*        0 *       17 *       150 *       123 *HWW
*        0 *       18 *       150 *       124 *HWW
*        0 *       19 *       150 *       125 *HWW
*        0 *       20 *       150 *       126 *HWW
*        0 *       21 *       150 *       127 *HWW
*        0 *       22 *       150 *       128 *HWW
*        0 *       23 *       150 *       130 *HWW
*        0 *       24 *       150 *       135 *HWW
*        0 *       25 *       150 *       140 *HWW
*        0 *       26 *       150 *       150 *HWW
*        0 *       27 *       150 *       160 *HWW
*        0 *       28 *       150 *       170 *HWW
*        0 *       29 *       150 *       180 *HWW
*        0 *       30 *       150 *       190 *HWW
*        0 *       31 *       150 *       200 *HWW
*        0 *       32 *       150 *       250 *HWW
*        0 *       33 *       150 *       300 *HWW
*        0 *       34 *       150 *       350 *HWW
*        0 *       35 *       150 *       400 *HWW
*        0 *       36 *       150 *       450 *HWW
*        0 *       37 *       150 *       500 *HWW
*        0 *       38 *       150 *       550 *HWW
*        0 *       39 *       150 *       600 *HWW
*        0 *       40 *       150 *       650 *HWW
*        0 *       41 *       150 *       700 *HWW
*        0 *       42 *       150 *       750 *HWW
*        0 *       43 *       150 *       800 *HWW
*        0 *       44 *       150 *       850 *HWW
*        0 *       45 *       150 *       900 *HWW
*        0 *       46 *       150 *       950 *HWW
*        0 *       47 *       150 *      1000 *HWW
*        0 *       48 *       100 *       110 *WH
*        0 *       49 *       100 *       115 *WH
*        0 *       50 *       100 *       120 *WH
*        0 *       51 *       100 *       122 *WH
*        0 *       52 *       100 *       123 *WH
*        0 *       53 *       100 *       124 *WH
*        0 *       54 *       100 *       125 *WH
*        0 *       55 *       100 *       126 *WH
*        0 *       56 *       100 *       127 *WH
*        0 *       57 *       100 *       128 *WH
*        0 *       58 *       100 *       130 *WH
*        0 *       59 *       100 *       135 *WH
*        0 *       60 *       100 *       140 *WH
*        0 *       61 *       100 *       150 *WH
*/
