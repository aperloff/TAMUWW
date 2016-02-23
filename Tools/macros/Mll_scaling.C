/*
Plot Mll for M50 and M10To50
Fit line to M50
scale M10To50 to that the distribution matches that line

Fit panel with range [50,75] ==> Set upper range to 70 and fit panel will change it by 5

 FCN=19.587 FROM MIGRAD    STATUS=CONVERGED      48 CALLS          49 TOTAL
                     EDM=7.21937e-07    STRATEGY= 1      ERROR MATRIX ACCURATE
  EXT PARAMETER                                   STEP         FIRST
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
   1  Constant     3.98610e+00   8.70060e-02   1.88004e-05   1.41603e-01
   2  Slope        6.16084e-02   1.32273e-03   2.85833e-07   9.28709e+00

Line at 47.5 = 1004.74942448678223
10To50: 440.440
50: 323.471

line = a*10To50+50
a=(line-50)/10To50
a = 1.5468132424
*/

void {
   TFile* _file0 = TFile::Open("WlnuJets_M-10To50.root","READ");
   TFile* _file1 = TFile::Open("WlnuJets_M-50.root","READ");
   _file0->cd("PS");
   TH1D* Mll10To50 = Mll;
   _file1->cd("PS");
   TH1D* Mll50 = Mll;
   Mll10To50->Scale(1.5468132424);
   Mll50->Add(Mll10To50);
   Mll50->Draw();
}
