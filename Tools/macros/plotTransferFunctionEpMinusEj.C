//Old TF Files
root -n TFAnalysis_outputHistograms1D.root
TFHistograms->cd()
gStyle->SetOptStat(0)
GGenJetEnergyMinusJetEnergy_0->SetTitle("")
GGenJetEnergyMinusJetEnergy_0->SetLineColor(kBlack)
UDSGenJetEnergyMinusJetEnergy_0->SetLineColor(kRed)
BGenJetEnergyMinusJetEnergy_0->SetLineColor(kBlue)
GGenJetEnergyMinusJetEnergy_0->SetLineWidth(2)
UDSGenJetEnergyMinusJetEnergy_0->SetLineWidth(2)
BGenJetEnergyMinusJetEnergy_0->SetLineWidth(2)
GGenJetEnergyMinusJetEnergy_0->GetXaxis()->SetTitle("E_{p}-E_{j} (GeV)")
GGenJetEnergyMinusJetEnergy_0->GetXaxis()->SetRangeUser(-200,100)
GGenJetEnergyMinusJetEnergy_0->GetYaxis()->SetRangeUser(0,0.14)
GGenJetEnergyMinusJetEnergy_0->Draw("hist")
UDSGenJetEnergyMinusJetEnergy_0->Draw("histsame")
BGenJetEnergyMinusJetEnergy_0->Draw("histsame")
TLegend* leg = new TLegend(0.7,0.7,0.89,0.89)
leg->SetLineWidth(0)
leg->SetLineStyle(0)
leg->AddEntry(GGenJetEnergyMinusJetEnergy_0,"Gluons","l")
leg->AddEntry(UDSGenJetEnergyMinusJetEnergy_0,"Light Jets","l")
leg->AddEntry(BGenJetEnergyMinusJetEnergy_0,"b Jets","l")
leg->Draw("same")
c1->SaveAs("Ep-Ej.eps")
c1->SaveAs("Ep-Ej.png")
c1->SaveAs("Ep-Ej.pdf")

//New Official TF Files
TFile *_file0 = TFile::Open("TF_TTbarMG.root")
gStyle->SetOptStat(0)
TF_Histo1D_G_00_24->SetTitle("")
TF_Histo1D_G_00_24->SetLineColor(kBlack)
TF_Histo1D_UDS_00_24->SetLineColor(kRed)
TF_Histo1D_B_00_24->SetLineColor(kBlue)
TF_Histo1D_G_00_24->SetLineWidth(2)
TF_Histo1D_UDS_00_24->SetLineWidth(2)
TF_Histo1D_B_00_24->SetLineWidth(2)
TF_Histo1D_G_00_24->GetXaxis()->SetTitle("E_{j}-E_{p} (GeV)")
TF_Histo1D_G_00_24->GetXaxis()->SetRangeUser(-100,100)
TF_Histo1D_G_00_24->GetYaxis()->SetRangeUser(0,0.14)
double GInt = TF_Histo1D_G_00_24->Integral()
double UDSInt = TF_Histo1D_UDS_00_24->Integral()
double BInt = TF_Histo1D_B_00_24->Integral()
TF_Histo1D_G_00_24->Scale(1.0/GInt)
TF_Histo1D_UDS_00_24->Scale(1.0/UDSInt)
TF_Histo1D_B_00_24->Scale(1.0/BInt)
TF_Histo1D_G_00_24->Draw("hist")
TF_Histo1D_UDS_00_24->Draw("histsame")
TF_Histo1D_B_00_24->Draw("histsame")
TLegend* leg = new TLegend(0.7,0.7,0.89,0.89)
leg->SetLineWidth(0)
leg->SetLineStyle(0)
leg->AddEntry(TF_Histo1D_G_00_24,"Gluons","l")
leg->AddEntry(TF_Histo1D_UDS_00_24,"Light Jets","l")
leg->AddEntry(TF_Histo1D_B_00_24,"b Jets","l")
leg->Draw("same")
c1->SaveAs("Ej-Ep.eps")
c1->SaveAs("Ej-Ep.png")
c1->SaveAs("Ej-Ep.pdf")
