#include "TAMUWW/Tools/interface/Fitter.hh"
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"

using namespace std;

int Fitter_x(string lep, string cat, string inFileLoc, string outFileLoc, vector<string> fproc,
             vector<string> signals, double FOM=3.0, string minimization = "chi2")
{
   Fitter histFit(lep, cat, fproc, inFileLoc, outFileLoc);
   histFit.setMinimization(minimization);
   histFit.readHistograms();
   
   // DEBUG
   histFit.printDataIntegral();
   histFit.printMCIntegrals();
   cout << "SUM: " << histFit.getMCSum() << endl;

   histFit.addSigBkgHistograms(signals);
   histFit.fitHistograms();
   
   // DEBUG
   histFit.printDataIntegral();
   histFit.printMCIntegrals();
   cout << "SUM: " << histFit.getMCSum() << endl;

   if (signals.size()>0) {
      //histFit.addSigBkgHistograms(signals);
      histFit.printSignalIntegral();
      histFit.printBackgroundIntegral();
      cout << "FOM: " << histFit.getFOM(FOM) << endl;
   }

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

   string         leptonCL     = cl.getValue<string>  ("lepton",                         "");
   string         objectCL     = cl.getValue<string>  ("object",                         "");
   string         inFileLocCL  = cl.getValue<string>  ("readLocation",                   "");
   string         outFileLocCL = cl.getValue<string>  ("writeLocation",                  "");
   vector<string> fproc        = cl.getVector<string> ("fproc",               "WJets:::QCD");
   vector<string> signals      = cl.getVector<string> ("signals", "ggH125:::qqH125:::WH125");
   double         fom          = cl.getValue<double>  ("fom",                           3.0);
   string         minimization = cl.getValue<string>  ("minimization",               "chi2");
   
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
/*
   if(objectCL != "MET" && objectCL != "WmT")
   {
      cout << "ERROR: object was not properly set. Options are MET and WmT." << endl;
      return 1;
   }
*/   
   if(inFileLocCL == "")
      inFileLocCL = "histos_" + leptonCL + ".root";
   
   if(outFileLocCL == "")
      outFileLocCL = "Fitter_" + leptonCL + "_" + objectCL + ".root";

   Fitter_x(leptonCL, objectCL, inFileLocCL, outFileLocCL, fproc, signals, fom, minimization);

   return 0;
}
