// TAMUWW libraries
#include "TAMUWW/Tools/interface/CSVreweight.hh"

using std::cout;
using std::endl;
using std::string;



//------------------------------------------------------------------------
CSVreweight::CSVreweight(){

   openFiles();
   
   getHistosFromFiles();

   closeFiles();

}

//------------------------------------------------------------------------
double CSVreweight::getWeight(EventNtuple * ntuple){

   double csvWgthf = 1.;
   double csvWgtC  = 1.;
   double csvWgtlf = 1.;
   double csvWgtTotal = 1.;
   int iSysHF = 0;
   int iSysC = 0;
   int iSysLF = 0;

   //Now, calculate the weights themselves.
   for (int iJet = 0; iJet < (int)ntuple->jLV.size(); ++iJet) {

      double csv = ntuple->jLV[iJet].jBtagDiscriminatorCSV;
      double jetPt = ntuple->jLV[iJet].Pt();
      double jetAbsEta = fabs( ntuple->jLV[iJet].Eta() );
      int flavor = abs( ntuple->jLV[iJet].partonFlavour );

      int iPt = -1; int iEta = -1;
      if (jetPt >=24.99 && jetPt<40) iPt = 0;    //changed jetPt >=29.99 to jetPt >=24.99 because our jet cut is lower
      else if (jetPt >=40 && jetPt<60) iPt = 1;
      else if (jetPt >=60 && jetPt<100) iPt = 2;
      else if (jetPt >=100 && jetPt<160) iPt = 3;
      else if (jetPt >=160 && jetPt<10000) iPt = 4;

      if (jetAbsEta >=0 &&  jetAbsEta<0.8) iEta = 0;
      else if ( jetAbsEta>=0.8 && jetAbsEta<1.6) iEta = 1;
      else if ( jetAbsEta>=1.6 && jetAbsEta<2.41) iEta = 2;

      if (iPt < 0 || iEta < 0){
         std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;
         //cout << "flavor is " << flavor << " csv is " << csv << " jetPt is " << iPt << " eta is " << iEta <<endl;
         continue;
      }
      if (abs(flavor) == 5 ){
         int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
         double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
         if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;
      }
      else if( abs(flavor) == 4 ){
         int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
         double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
 
         if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;
      }
      else {
         if (iPt >=2) iPt=2;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
         int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
         double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
         if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;
      }

      csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;
      //cout << "Total weight is " << csvWgtTotal << endl;
   }

   return csvWgtTotal;
    
}//getWeight

//-----------------------------------------------------------------------
bool CSVreweight::checkFiles(){

   if (!f_CSVwgt_HF->IsOpen()){
      cout << "Error in CSVreweight::checkFiles Could not open file " << DefaultValues::getConfigPath()+"csv_rwt_hf_IT.root" << endl;
      return false;
   }
   if (!f_CSVwgt_LF->IsOpen()){
      cout << "Error in CSVreweight::checkFiles Could not open file " << DefaultValues::getConfigPath()+"csv_rwt_lf_IT.root" << endl;
      return false;
   }

   return true;

}

//-----------------------------------------------------------------------
void CSVreweight::openFiles(){

   f_CSVwgt_HF = TFile::Open(TString(DefaultValues::getConfigPath()+"csv_rwt_hf_IT.root"));
   f_CSVwgt_LF = TFile::Open(TString(DefaultValues::getConfigPath()+"csv_rwt_lf_IT.root"));

   assert(checkFiles());

   return;
}

//-----------------------------------------------------------------------
void CSVreweight::closeFiles(){
   f_CSVwgt_HF->Close();
   f_CSVwgt_LF->Close();
}

//-----------------------------------------------------------------------
void CSVreweight::getHistosFromFiles(){

   //loop through uncertainty options.  For now use case 0, i.e. nominal.
   for( int iSys=0; iSys<1; iSys++ ){
      TString syst_csv_suffix_hf = "final";
      TString syst_csv_suffix_c = "final";
      TString syst_csv_suffix_lf = "final";

      TH1D * tempH;
      TString tempS;

      for( int iPt=0; iPt<5; iPt++ ) {
         tempS = Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data());
         tempH = (TH1D*)f_CSVwgt_HF->Get(tempS);
         if (!tempH){
            cout << "Error in CSVreweight::getHistosFromFiles Could not get histo " << tempS << " from " << f_CSVwgt_HF->GetName() << endl;
            return; 
         }
         h_csv_wgt_hf[iSys][iPt] = (TH1D*)tempH->Clone();
         h_csv_wgt_hf[iSys][iPt]->SetDirectory(0);
      }
      //if( iSys<5 ){
      for( int iPt=0; iPt<5; iPt++ ) {
         tempS = Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data());
         tempH = (TH1D*)f_CSVwgt_HF->Get(tempS);
         if (!tempH){
            cout << "Error in CSVreweight::getHistosFromFiles Could not get histo " << tempS << " from " << f_CSVwgt_HF->GetName() << endl;
            return; 
         }
         c_csv_wgt_hf[iSys][iPt] = (TH1D*)tempH->Clone();
         c_csv_wgt_hf[iSys][iPt]->SetDirectory(0);
      }
      //}
      for( int iPt=0; iPt<3; iPt++ ){
         for( int iEta=0; iEta<3; iEta++ ) {
            tempS = Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data());
            tempH = (TH1D*)f_CSVwgt_LF->Get(tempS);
            if (!tempH){
               cout << "Error in CSVreweight::getHistosFromFiles Could not get histo " << tempS << " from " << f_CSVwgt_LF->GetName() << endl;
               return; 
            }
            h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)tempH->Clone();
            h_csv_wgt_lf[iSys][iPt][iEta]->SetDirectory(0);
         }
      }
   }

}

