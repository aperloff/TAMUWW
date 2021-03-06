#include "TAMUWW/Tools/interface/BackgroundEstimator.hh"
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"

using namespace std;

int BackgroundEstimator_x(string lep, string cat, string inFileLoc, string outFileLoc)
{
   BackgroundEstimator histFit(lep, cat, inFileLoc, outFileLoc);
   
   histFit.readHistograms();
   
   // DEBUG
   double sum = 0;
   cout << "Data Integral: ";
   cout << histFit.dataHistogram->Integral() << endl;
   for(map<string, TH1D*>::iterator mapit = histFit.monteCarloHistograms.begin(); mapit != histFit.monteCarloHistograms.end(); mapit++)
   {
      cout << mapit->first << " Integral: ";
      cout << mapit->second->Integral() << endl;
      
      sum += mapit->second->Integral();
   }
   cout << "SUM: " << sum << endl;
   
   histFit.fitHistograms();
   
   // DEBUG
   sum = 0;
   cout << "Data Integral: ";
   cout << histFit.dataHistogram->Integral() << endl;
   for(map<string, TH1D*>::iterator mapit = histFit.monteCarloHistograms.begin(); mapit != histFit.monteCarloHistograms.end(); mapit++)
   {
      cout << mapit->first << " Integral: ";
      cout << mapit->second->Integral() << endl;
      
      sum += mapit->second->Integral();
   }
   cout << "SUM: " << sum << endl;
   
   histFit.writeHistograms();
   
   return 0;
}

//##################################################
//###################### MAIN ######################
//##################################################

int main(int argc,char**argv)
{
   CommandLine cl;
   if (!cl.parse(argc,argv)) return 0;

   string leptonCL     = cl.getValue<string> ("lepton",        "");
   string objectCL     = cl.getValue<string> ("object",        "");
   string inFileLocCL  = cl.getValue<string> ("readLocation",  "");
   string outFileLocCL = cl.getValue<string> ("writeLocation", "");
   
   if (!cl.check()) 
      return 0;
   cl.print();
   
   // Check that lepton is set to either muon or electron
   if(leptonCL != "muon" && leptonCL != "electron")
   {
      cout << "ERROR: lepton was not properly set. Options are electron and muon." << endl;
      return 1;
   }
   // Check that object is set to either MET or WmT
   if(objectCL != "MET" && objectCL != "WmT")
   {
      cout << "ERROR: object was not properly set. Options are MET and WmT." << endl;
      return 1;
   }
   
   if(inFileLocCL == "")
      inFileLocCL = "outputFile_" + leptonCL + ".root";
   
   if(outFileLocCL == "")
      outFileLocCL = "BackgroundEstimation_" + leptonCL + "_" + objectCL + ".root";
   
   BackgroundEstimator_x(leptonCL, objectCL, inFileLocCL, outFileLocCL);
   return 0;
}
