// This Program Takes input Monte Carlo histograms of several processes
// and fits them to a data histogram. It only scales two processes.
// 
// Written by Alexx Perloff
// Started on the 16th of October, 2012

// TAMU Libraries
#include "TAMUWW/Tools/interface/Fitter.hh"

using namespace std;
using namespace ROOT::Minuit2;

//##################################################
//################ PUBLIC FUNCTIONS ################
//##################################################

Fitter::Fitter(string lepton, string object, vector<string> fproc, string inFileLoc, string outFileLoc)
{
   // Set default debug options
   rebinSizeDEBUG = 1;
   debug = false;

   leptonName = lepton;
   objectName = object;
   initializeFileLocations(inFileLoc, outFileLoc);
   initializeHistNames();

   // The default values (before fitting) of Chi2 and parameters
   scaleParameters = vector<double> (2,1.0);
   reducedChiSquared = 0.0;

   //resultStack = new THStack((objectName + "_MonteCarlo").c_str(), (objectName + "_MonteCarlo").c_str());

   HistogramsFitter::fProcessNames = fproc;
   for(unsigned int nproc=0; nproc<HistogramsFitter::fProcessNames.size(); nproc++) {
      HistogramsFitter::fProcessXsec.push_back(DV.getCrossSectionAndError(HistogramsFitter::fProcessNames[nproc]).first);
      HistogramsFitter::fProcessXsecError.push_back(DV.getCrossSectionAndError(HistogramsFitter::fProcessNames[nproc]).second);
   }
}

Fitter::~Fitter() {
   HistogramsFitter::fProcessNames.clear();
   HistogramsFitter::fProcessXsec.clear();
   HistogramsFitter::fProcessXsecError.clear();
   
   // Destroy the MC histograms
   for(map<string,TH1D*>::iterator it=HistogramsFitter::monteCarloHistograms.begin(); it!=HistogramsFitter::monteCarloHistograms.end(); it++) {
      (*it).second->~TH1D();
      if((*it).second)
         delete (*it).second;
   }
   HistogramsFitter::monteCarloHistograms.clear();

   // Destroy data histogram
   if(HistogramsFitter::dataHistogram) {
      HistogramsFitter::dataHistogram->~TH1D();
      delete HistogramsFitter::dataHistogram;
   }

   // Destroy other histograms
   if(HistogramsFitter::signalHistogram) {
      HistogramsFitter::signalHistogram->~TH1D();
      delete HistogramsFitter::signalHistogram;
   }

   if(HistogramsFitter::backgroundHistogram) {
      HistogramsFitter::backgroundHistogram->~TH1D();
      delete HistogramsFitter::backgroundHistogram;
   }
}

void Fitter::readHistograms()
{
   TH1::AddDirectory(kFALSE);
   TFile* rootInFile = TFile::Open(inRootFileLocation.c_str(), "READ");
   if (!rootInFile->IsOpen()) {  cout<<"Can't open "<<inRootFileLocation<<endl; return; }

   for(unsigned int i = 0; i < histNames.size(); i++)
   {
      string name = histNames[i];
      TString hist = objectName + "_" + name + "_" + leptonName;
      cout << "\tFitter::readHistograms Getting histogram named " << hist << endl;
      
      if(gDirectory->Get(hist))
         HistogramsFitter::monteCarloHistograms[name] = (TH1D*) (gDirectory->Get(hist)->Clone());
      else
         cout << "\t\tFitter::readHistograms Cannot find histogram " << hist << endl
              << "\t\t\tSkipping this histogram." << endl;

      if (debug)
         HistogramsFitter::monteCarloHistograms[name]->Rebin(rebinSizeDEBUG);
   }

   cout << "\tFitter::readHistograms Getting histograms named " << TString(objectName + "_SingleEl_Data_" + leptonName) 
        << " and " << TString(objectName + "_SingleMu_Data_" + leptonName)  << endl;
   if(gDirectory->Get(TString(objectName + "_SingleEl_Data_" + leptonName)) && gDirectory->Get(TString(objectName + "_SingleMu_Data_" + leptonName))) {
      HistogramsFitter::dataHistogram = (TH1D*) (gDirectory->Get(TString(objectName + "_SingleEl_Data_" + leptonName))->Clone());
      HistogramsFitter::dataHistogram->Add((TH1D*) (gDirectory->Get(TString(objectName + "_SingleMu_Data_" + leptonName))->Clone()));
   }
   else if(gDirectory->Get(TString(objectName + "_SingleEl_Data_" + leptonName)) && !gDirectory->Get(TString(objectName + "_SingleMu_Data_" + leptonName)))
      HistogramsFitter::dataHistogram = (TH1D*) (gDirectory->Get(TString(objectName + "_SingleEl_Data_" + leptonName))->Clone());
   else if(!gDirectory->Get(TString(objectName + "_SingleEl_Data_" + leptonName)) && gDirectory->Get(TString(objectName + "_SingleMu_Data_" + leptonName)))
      HistogramsFitter::dataHistogram = (TH1D*) (gDirectory->Get(TString(objectName + "_SingleMu_Data_" + leptonName))->Clone());
   else
      cout << "\t\tFitter::readHistograms Cannot find histograms " << TString(objectName + "_SingleEl_Data_" + leptonName) 
           << " and " << TString(objectName + "_SingleMu_Data_" + leptonName) << endl
           << "\t\t\tSkipping these histograms." << endl;

   if (debug)
      HistogramsFitter::dataHistogram->Rebin(rebinSizeDEBUG);
   
   rootInFile->Close();
}

void Fitter::fitHistograms()
{
   vector<double> parameters = fitAndReturnParameters();
   
   for (unsigned int i=0; i<HistogramsFitter::fProcessNames.size(); i++) {
      HistogramsFitter::monteCarloHistograms[HistogramsFitter::fProcessNames[i]]->Scale(parameters[i]);
   }
   
   //mcStack->Modified();
}

void Fitter::writeHistograms()
{
   if (outRootFileLocation == "")
   {
      cout << "ERROR: writeLocation has not been set. Cannot write Histograms." << endl;
      return;
   }

   TFile* rootInFile = TFile::Open(inRootFileLocation.c_str());
   if (!rootInFile->IsOpen()) {  cout<<"Can't open "<<inRootFileLocation<<endl; return; }
   TIter next(rootInFile->GetListOfKeys());
   TKey* key(0);
   //vector<string> plotPrefixes = getPlotNames();
   
   vector<TCanvas*> inputCanvases;
   
   TCanvas* inputCanvas;
   TPad* inputPad;
   THStack* inputMCStack;
   TList* inputMCList;
   //TH1D* histQCD;
   //TH1D* histWJets;
   vector<TH1D*> histsToScale;
   
   cout << "\tFitter::writeHistograms Scaling input histograms...";
   //for(unsigned int i = 0; i < plotPrefixes.size(); i++) {
   while ((key=(TKey*)next())) {
      if (strcmp(key->GetClassName(),"TCanvas")!=0) continue;

      //string prefix = plotPrefixes[i] + "_" + leptonName;
      string prefix = key->GetName();

      //inputCanvas = (TCanvas*)(key->ReadObj())->Clone();
      inputCanvas = (TCanvas*)(gDirectory->Get((prefix).c_str()))->Clone();
      inputPad = (TPad *) inputCanvas->GetPrimitive((prefix + "_1").c_str());
      inputMCStack = (THStack*) inputPad->GetPrimitive((prefix + "_" + leptonName + "_stackMC").c_str());
      if(!inputMCStack) {
         cout << endl << "ERROR Fitter::writeHistograms inputMCStack not found!" << endl;
         return;
      }
      inputMCList = inputMCStack->GetHists();

      for (unsigned int nproc=0; nproc<HistogramsFitter::fProcessNames.size(); nproc++) {
         /*
         if(HistogramsFitter::fProcessNames[nproc] == "QCD") {
            if(leptonName == "electron")
               //histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_QCD_ElEnriched_electron_clone").c_str()));
               histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_QCD_ElFULL_electron_clone").c_str()));
            else if (leptonName == "muon")
               histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_QCD_MuEnriched_clone").c_str()));
            else
               histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_QCD").c_str()));
         }
         else if(HistogramsFitter::fProcessNames[nproc] == "WJets") {
           if(leptonName == "electron")
               histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_WJets_electron_clone").c_str()));
            else if (leptonName == "muon")
               histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_WJets_muon_clone").c_str()));
            else
               histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_WJets").c_str()));
         }
         */
         //else histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_" + HistogramsFitter::fProcessNames[nproc]+ "_" + leptonName + "_clone").c_str()));
         histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_" + HistogramsFitter::fProcessNames[nproc]+ "_" + leptonName + "_clone").c_str()));

         histsToScale[nproc]->Scale(scaleParameters[nproc]);
      }

      inputMCStack->Modified();
      
      inputCanvases.push_back(inputCanvas);
   }
   cout << "DONE" << endl; 

   cout << "\tFitter::writeHistograms Writing output histograms...";
   TFile rootOutFile;
   rootOutFile.Open(outRootFileLocation.c_str(), "RECREATE");
   
   //canvas->Write();
   
   for(unsigned int i = 0; i < inputCanvases.size(); i++)
   {
      inputCanvases[i]->Write();
   }
   HistogramsFitter::signalHistogram->Write();
   HistogramsFitter::backgroundHistogram->Write();
   
   cout << "DONE" << endl;

   rootOutFile.Close();
   rootInFile->Close();
}

void Fitter::addSigBkgHistograms(vector<string> sig) {
   HistogramsFitter::signalHistogram = new TH1D("signalHistogram","signalHistograms",
                                                HistogramsFitter::monteCarloHistograms.begin()->second->GetNbinsX(),
                                                HistogramsFitter::monteCarloHistograms.begin()->second->GetXaxis()->GetXmin(),
                                                HistogramsFitter::monteCarloHistograms.begin()->second->GetXaxis()->GetXmax());
   HistogramsFitter::signalHistogram->Sumw2();
   HistogramsFitter::backgroundHistogram = new TH1D("backgroundHistogram","backgroundHistogram",
                                                    HistogramsFitter::monteCarloHistograms.begin()->second->GetNbinsX(),
                                                    HistogramsFitter::monteCarloHistograms.begin()->second->GetXaxis()->GetXmin(),
                                                    HistogramsFitter::monteCarloHistograms.begin()->second->GetXaxis()->GetXmax());
   HistogramsFitter::backgroundHistogram->Sumw2();

   for(map<string,TH1D*>::iterator it=HistogramsFitter::monteCarloHistograms.begin(); it!=HistogramsFitter::monteCarloHistograms.end(); it++) {
      if (DV.vfind(sig,(*it).first)!=-1) {
         HistogramsFitter::signalHistogram->Add((*it).second);
         cout << "Fitter::Added the " << (*it).first << " histogram to the signalHistogram." << endl;
         cout << "\t" << setw(26) << (*it).first << " Entries: " << (*it).second->GetEntries() << endl;
         cout << "\t" << setw(26) << "signalHistogram Entries: " << HistogramsFitter::signalHistogram->GetEntries() << endl;
      }
      else {
         HistogramsFitter::backgroundHistogram->Add((*it).second);
         cout << "Fitter::Added the " << (*it).first << " histogram to the backgroundHistogram." << endl;
         cout << "\t" << setw(26) << (*it).first << " Entries: " << (*it).second->GetEntries() << endl;
         cout << "\t" << setw(26) << "backgroundHistogram Entries: " << HistogramsFitter::backgroundHistogram->GetEntries() << endl; 
      }
   }
}

void Fitter::setWriteLocation(string outFileLoc)
{
   outRootFileLocation = outFileLoc;
}

void Fitter::setReadLocation(string inFileLoc)
{
   inRootFileLocation = inFileLoc;
}

void Fitter::setRebinSizeDEBUG(unsigned int rebinSize)
{
   debug = true;
   rebinSizeDEBUG = rebinSize;
}

vector<double> Fitter::getParameters()
{
   return scaleParameters;
}

double Fitter::getReducedChiSquared()
{
   return reducedChiSquared;
}

double Fitter::getFOM(double FOM) {
   if(HistogramsFitter::signalHistogram->GetEntries()==0 || HistogramsFitter::backgroundHistogram->GetEntries()==0) {
      cout << "Fitter::WARNING Either the signal or background histogram has zero entries!" << endl;
      cout << "\t" << setw(22) << "Signal Entries: " << HistogramsFitter::signalHistogram->GetEntries() << endl;    
      cout << "\t" << setw(22) << "Background Entries: " << HistogramsFitter::backgroundHistogram->GetEntries() << endl;
   }

   // set signal bins to zero if background bin is zero
   // set background bins to zero if signal bin is zero
   for(int i=1; i<=HistogramsFitter::backgroundHistogram->GetNbinsX(); i++) {
      if (HistogramsFitter::signalHistogram->GetBinContent(i)>0 && HistogramsFitter::backgroundHistogram->GetBinContent(i)==0) {
         HistogramsFitter::signalHistogram->SetBinContent(i,0.0);
         HistogramsFitter::signalHistogram->SetBinError(i,0.0);
      }
   }

   if (FOM==1) {
      double fom1 = 0.0;
      for(int i=1; i<=HistogramsFitter::backgroundHistogram->GetNbinsX(); i++)
         fom1 += TMath::Power((HistogramsFitter::signalHistogram->GetBinContent(i)/TMath::Sqrt(HistogramsFitter::backgroundHistogram->GetBinContent(i))),2);
      return TMath::Sqrt(fom1);
   }
   else if (FOM==2) {
      double fom2 = 0.0;
      for(int i=1; i<=HistogramsFitter::backgroundHistogram->GetNbinsX(); i++)
         fom2 += TMath::Power((HistogramsFitter::signalHistogram->GetBinContent(i)/TMath::Sqrt(HistogramsFitter::signalHistogram->GetBinContent(i)+HistogramsFitter::backgroundHistogram->GetBinContent(i))),2);
      return TMath::Sqrt(fom2);
   }
   else if (FOM==2.5)
      return FigureOfMerit::usingChi2(HistogramsFitter::signalHistogram,
                                      HistogramsFitter::backgroundHistogram,0.0001);
   else if (FOM==3)
      return FigureOfMerit::usingShapeFromTemplates(HistogramsFitter::signalHistogram,
                                                    HistogramsFitter::backgroundHistogram,0.0001);
   else {
      cout << "Fitter::WARNING The FOM type was not specified. Using defaul (FOM3)." << endl;
      return FigureOfMerit::usingShapeFromTemplates(HistogramsFitter::signalHistogram,
                                                    HistogramsFitter::backgroundHistogram,0.0001);
   }
}

//##################################################
//############### PRIVATE FUNCTIONS ################
//##################################################

void Fitter::initializeFileLocations(string inFileLoc, string outFileLoc)
{
   inRootFileLocation = inFileLoc;
   outRootFileLocation = outFileLoc;
}

void Fitter::initializeHistNames()
{
   //histNames.push_back("Diboson");
   histNames.push_back("WW");
   histNames.push_back("WZ");
   //histNames.push_back("ZZ");
   //histNames.push_back("STop");
   histNames.push_back("STopS_T");
   histNames.push_back("STopS_Tbar");
   histNames.push_back("STopT_T");
   histNames.push_back("STopT_Tbar");
   histNames.push_back("STopTW_T");
   histNames.push_back("STopTW_Tbar");
   histNames.push_back("WJets");
   histNames.push_back("ZJets");
   histNames.push_back("TTbar");
   //histNames.push_back("QCD_ElEnriched");
   histNames.push_back("QCD_ElFULL");
   //histNames.push_back("QCD_MuFULL");
   histNames.push_back("ggH125");
   histNames.push_back("qqH125");
   histNames.push_back("WH125");
}

vector<double> Fitter::fitAndReturnParameters()
{
   ROOT::Math::Functor funcFit(&fitFunc,2);
   Minuit2Minimizer* minFit = new Minuit2Minimizer(kMigrad);
   
   // Tolerance and printouts
   minFit->SetPrintLevel(3);
   minFit->SetStrategy(2);
   minFit->SetMaxFunctionCalls(100000);
   minFit->SetMaxIterations(100000);
   //EDM must be less than 0.001*tolerance*up for minimization to be successful
   minFit->SetTolerance(0.001);
   //1.0 for chi-square fit and 1 stdev errors
   //4.0 for chi-square fit and 2 stdev errors
   //0.5 for negative-log-likelihood function
   minFit->SetErrorDef(1.0);
   minFit->SetValidError(true);
   
   // Fitting
   minFit->SetFunction(funcFit);
   
   double parameter[2];
   if(DV.vfind(HistogramsFitter::fProcessNames,"QCD")<0) {
      parameter[0] = 1.0;
      parameter[1] = 1.0;
   }
   else {
      parameter[0] = 1.0;
      parameter[1] =  0.0001;
   }
   double step[2] = {0.0001, 0.0001};
   double lower[2] = {0.0, 0.0};
   
   //minFit->SetVariable(0, "Par0", parameter[0], step[0]);
   //minFit->SetVariable(1, "Par1", parameter[1], step[1]);
   minFit->SetLowerLimitedVariable(0, "Par0", parameter[0], step[0], lower[0]);
   minFit->SetLowerLimitedVariable(1, "Par1", parameter[1], step[1], lower[1]);
   
   cout << endl << "Beginning Fit..." << endl;
   if(minFit->Minimize())
      cout << endl << "Fit Successful!" << endl;
   else
      cout << endl << "Fit Failed." << endl;
   
   cout << endl << "##### FIT RESULTS #####" << endl;
   minFit->PrintResults();
   
   reducedChiSquared = minFit->MinValue() / (HistogramsFitter::dataHistogram->GetNbinsX() - 2);
   
   cout << "Chi2\t  = " << minFit->MinValue() << endl;
   cout << "NDF\t  = " << HistogramsFitter::dataHistogram->GetNbinsX() - 2 << endl;
   cout << "Chi2/NDF\t  = " << reducedChiSquared << endl;
   
   scaleParameters[0] = minFit->X()[0];
   scaleParameters[1] = minFit->X()[1];
   
   vector<double> ret;
   for(unsigned int i=0; i<minFit->NDim(); i++) {
      ret.push_back(minFit->X()[i]);
   }

   return ret;
}

vector<string> Fitter::getPlotNames()
{
   vector<string> names;
   
   names.push_back("AngleJ1J2");
   names.push_back("BetaJ1BetaJ2");
   names.push_back("DeltaEtaJ1J2");
   names.push_back("DeltaPhi_J1J2");
   names.push_back("DeltaPhi_LJ1");
   names.push_back("DeltaPhi_LJ1_vs_J1J2");
   //names.push_back("DeltaPhi_LJ1_vs_J1J2_Positive");
   //names.push_back("DeltaPhi_LJ1_vs_J1J2_Negative");
   //names.push_back("DeltaPhi_LJ1_vs_J1J2_Subtracted");
   names.push_back("DeltaRJ1J2");
   names.push_back("DeltaRLepMET");
   names.push_back("EJ1EJ2");
   names.push_back("Jet1Eta");
   names.push_back("Jet1Phi");
   names.push_back("Jet1Pt");
   names.push_back("Jet2Eta");
   names.push_back("Jet2Phi");
   names.push_back("Jet2Pt");
   names.push_back("LeptEta");
   names.push_back("LeptPhi");
   names.push_back("LeptPt");
   names.push_back("MET");
   names.push_back("Mjj");
   names.push_back("MjjmWmT");
   names.push_back("Mlvjj");
   names.push_back("MWjjVsMWlv");
   names.push_back("Ptjj");
   names.push_back("WmT");
   //names.push_back("WmT_Positive");
   //names.push_back("WmT_Negative");
   //names.push_back("WmT_Subtracted");
   names.push_back("j1Pt_Mjj");
   names.push_back("j2Pt_Mjj");
   names.push_back("jjlvPhi");
   //names.push_back("npv");
   
   return names;
}

void Fitter::printDataIntegral(bool includeOverflowUnderflow) {
   double i = 0;
   if(includeOverflowUnderflow)
      i = HistogramsFitter::dataHistogram->Integral(0, HistogramsFitter::dataHistogram->GetNbinsX()+1);
   else
      i = HistogramsFitter::dataHistogram->Integral(1, HistogramsFitter::dataHistogram->GetNbinsX());

   cout << "Data Integral: "
        << i << endl;
}

void Fitter::printMCIntegrals(bool includeOverflowUnderflow) {
   for(map<string, TH1D*>::iterator mapit = HistogramsFitter::monteCarloHistograms.begin(); mapit != HistogramsFitter::monteCarloHistograms.end(); mapit++)
   {
      cout << mapit->first << " Integral: ";
      if(mapit->second) {
         if(includeOverflowUnderflow)
            cout << mapit->second->Integral(0,mapit->second->GetNbinsX()+1) << endl;
         else
            cout << mapit->second->Integral(1,mapit->second->GetNbinsX()) << endl;
      }
      else
         cout << "ERROR Fitter::printMCIntegrals" << mapit->first << " histogram not found!" << endl;
   }
}

void Fitter::printSignalIntegral(bool includeOverflowUnderflow) {
   double i = 0;
   if(includeOverflowUnderflow)
      i = HistogramsFitter::signalHistogram->Integral(0, HistogramsFitter::signalHistogram->GetNbinsX()+1);
   else
      i = HistogramsFitter::signalHistogram->Integral(1, HistogramsFitter::signalHistogram->GetNbinsX());

   cout << endl <<"Signal Integral: " 
        << i << endl;
}

void Fitter::printBackgroundIntegral(bool includeOverflowUnderflow) {
   double i = 0;
   if(includeOverflowUnderflow)
      i =  HistogramsFitter::backgroundHistogram->Integral(0, HistogramsFitter::backgroundHistogram->GetNbinsX()+1);
   else
      i =  HistogramsFitter::backgroundHistogram->Integral(1, HistogramsFitter::backgroundHistogram->GetNbinsX());

   cout << endl <<"Background Integral: " 
        << i << endl;
}

double Fitter::getMCSum(bool includeOverflowUnderflow) {
   double sum = 0;
   for(map<string, TH1D*>::iterator mapit = HistogramsFitter::monteCarloHistograms.begin(); mapit != HistogramsFitter::monteCarloHistograms.end(); mapit++)
   {
      if(mapit->second) {
         if(includeOverflowUnderflow)
            sum += mapit->second->Integral(0,mapit->second->GetNbinsX()+1);
         else
            sum += mapit->second->Integral(1,mapit->second->GetNbinsX());
      }
      else {
         sum += 0;
         cout << "ERROR Fitter::GetMCSum " << mapit->first << " histogram not found!" << endl;
      }
   }
   return sum;
}

//##################################################
//################ STATIC FUNCTIONS ################
//##################################################

double Fitter::fitFunc(const double *par)
{
   HistogramsFitter::dataHistogram->GetNbinsX();
   double chiSquare = 0;

   for(int bin = 1; bin <= HistogramsFitter::dataHistogram->GetNbinsX(); bin++)
   {
      double mc = 0;
      double mcError2 = 0;
      double data = 0;
      double dataError2 = 0;

      for (unsigned int nproc=0; nproc<HistogramsFitter::fProcessNames.size(); nproc++) {
         mc += par[nproc]*HistogramsFitter::monteCarloHistograms[HistogramsFitter::fProcessNames[nproc]]->GetBinContent(bin);
         mcError2 += pow(par[nproc]*HistogramsFitter::monteCarloHistograms[HistogramsFitter::fProcessNames[nproc]]->GetBinError(bin), 2);
      }

      data += HistogramsFitter::dataHistogram->GetBinContent(bin);
      dataError2 += pow(HistogramsFitter::dataHistogram->GetBinError(bin), 2);

      for(map<string, TH1D*>::iterator it = HistogramsFitter::monteCarloHistograms.begin(); it != HistogramsFitter::monteCarloHistograms.end(); it++)
      {
         bool notFProc = true;
         for (unsigned int nproc=0; nproc<HistogramsFitter::fProcessNames.size(); nproc++) {
            if ((*it).first == HistogramsFitter::fProcessNames[nproc])
               notFProc = false;
         }
         if(notFProc)
         {
            mc += (*it).second->GetBinContent(bin);
            mcError2 += pow((*it).second->GetBinError(bin), 2);
         }
      }

      if(mcError2 == 0 || mc == 0 || data == 0)
         continue;
      else
      {
         double chi = (data - mc) / sqrt(mcError2 + dataError2);
         chiSquare += chi*chi;

         // DEBUG
         // if (fabs(par[0] - par[1]) > 0.1)
         //    cout << "Bin:" << bin
         //         << "   Data:" << data 
         //         << "   dataError:" << sqrt(dataError2) 
         //         << "   MC:" << mc 
         //         << "   MCError:" << sqrt(mcError2) 
         //         << "   Par0:" << par[0] 
         //         << "   Par1:" << par[1] 
         //         << "   Chi2:" << chiSquare << endl;
      }
   }

   // cout << endl << "Par0:" << par[0] << " Par1:" << par[1] << " Chi2:" << chiSquare << endl << endl;
   //Add Gaussian Constraints
   for(unsigned int nproc=0; nproc<HistogramsFitter::fProcessNames.size(); nproc++) {
      if(HistogramsFitter::fProcessXsec[nproc]!=0 && HistogramsFitter::fProcessXsecError[nproc]!=0) {
         //cout << "\tchiSquare=" << chiSquare << "\tHistogramsFitter::fProcessXsec[nproc]=" << HistogramsFitter::fProcessXsec[nproc] << "\tHistogramsFitter::fProcessXsecError[nproc]=" << HistogramsFitter::fProcessXsecError[nproc] << endl; 
         //chiSquare*=TMath::Gaus(par[nproc],1,HistogramsFitter::fProcessXsecError[nproc]/HistogramsFitter::fProcessXsec[nproc]);
         chiSquare+=TMath::Power((par[nproc]-1.0)/(HistogramsFitter::fProcessXsecError[nproc]/HistogramsFitter::fProcessXsec[nproc]),2);
      }
   }

   return chiSquare;
}
