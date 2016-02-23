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

Fitter::Fitter(string lepton, string object, vector<string> fproc, string inFileLoc, string outFileLoc) :
minFit(0), funcFit(0)
{
   // Set default debug options
   rebinSizeDEBUG = 1;
   debug = false;
   minimizationTechnique = "chi2";

   leptonName = lepton;
   objectName = object;
   initializeFileLocations(inFileLoc, outFileLoc);
   initializeHistNames();

   // The default values (before fitting) of Chi2 and parameters
   scaleParameters = vector<double> (fproc.size(),1.0);
   reducedChiSquared = 0.0;

   //resultStack = new THStack((objectName + "_MonteCarlo").c_str(), (objectName + "_MonteCarlo").c_str());

   fProcessNames = fproc;
   for(unsigned int nproc=0; nproc<fProcessNames.size(); nproc++) {
      fProcessXsec.push_back(DV.getCrossSectionAndError(fProcessNames[nproc]).first);
      fProcessXsecError.push_back(DV.getCrossSectionAndError(fProcessNames[nproc]).second);
   }
}

Fitter::~Fitter() {
   fProcessNames.clear();
   fProcessXsec.clear();
   fProcessXsecError.clear();
   
   // Destroy the MC histograms
   for(map<string,TH1D*>::iterator it=monteCarloHistograms.begin(); it!=monteCarloHistograms.end(); it++) {
      (*it).second->~TH1D();
      if((*it).second)
         delete (*it).second;
   }
   monteCarloHistograms.clear();

   // Destroy data histogram
   if(dataHistogram) {
      dataHistogram->~TH1D();
      delete dataHistogram;
   }

   // Destroy other histograms
   if(signalHistogram) {
      signalHistogram->~TH1D();
      delete signalHistogram;
   }

   if(backgroundHistogram) {
      backgroundHistogram->~TH1D();
      delete backgroundHistogram;
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
         monteCarloHistograms[name] = (TH1D*) (gDirectory->Get(hist)->Clone());
      else
         cout << "\t\tFitter::readHistograms Cannot find histogram " << hist << endl
              << "\t\t\tSkipping this histogram." << endl;

      if (debug)
         monteCarloHistograms[name]->Rebin(rebinSizeDEBUG);
   }

   cout << "\tFitter::readHistograms Getting histograms named " << TString(objectName + "_SingleEl_Data_" + leptonName) 
        << " and " << TString(objectName + "_SingleMu_Data_" + leptonName)  << endl;
   if(gDirectory->Get(TString(objectName + "_SingleEl_Data_" + leptonName)) && gDirectory->Get(TString(objectName + "_SingleMu_Data_" + leptonName))) {
      dataHistogram = (TH1D*) (gDirectory->Get(TString(objectName + "_SingleEl_Data_" + leptonName))->Clone());
      dataHistogram->Add((TH1D*) (gDirectory->Get(TString(objectName + "_SingleMu_Data_" + leptonName))->Clone()));
   }
   else if(gDirectory->Get(TString(objectName + "_SingleEl_Data_" + leptonName)) && !gDirectory->Get(TString(objectName + "_SingleMu_Data_" + leptonName)))
      dataHistogram = (TH1D*) (gDirectory->Get(TString(objectName + "_SingleEl_Data_" + leptonName))->Clone());
   else if(!gDirectory->Get(TString(objectName + "_SingleEl_Data_" + leptonName)) && gDirectory->Get(TString(objectName + "_SingleMu_Data_" + leptonName)))
      dataHistogram = (TH1D*) (gDirectory->Get(TString(objectName + "_SingleMu_Data_" + leptonName))->Clone());
   else
      cout << "\t\tFitter::readHistograms Cannot find histograms " << TString(objectName + "_SingleEl_Data_" + leptonName) 
           << " and " << TString(objectName + "_SingleMu_Data_" + leptonName) << endl
           << "\t\t\tSkipping these histograms." << endl;

   if (debug)
      dataHistogram->Rebin(rebinSizeDEBUG);
   
   rootInFile->Close();
}

void Fitter::fitHistograms()
{
   vector<double> parameters;
   if(minimizationTechnique == "TFractionFitter")
      parameters = fitTFractionFitter();
   else
      parameters = fitMinimization();
   
   for (unsigned int i=0; i<fProcessNames.size(); i++) {
      monteCarloHistograms[fProcessNames[i]]->Scale(parameters[i]);
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

      for (unsigned int nproc=0; nproc<fProcessNames.size(); nproc++) {
         /*
         if(fProcessNames[nproc] == "QCD") {
            if(leptonName == "electron")
               //histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_QCD_ElEnriched_electron_clone").c_str()));
               histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_QCD_ElFULL_electron_clone").c_str()));
            else if (leptonName == "muon")
               histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_QCD_MuEnriched_clone").c_str()));
            else
               histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_QCD").c_str()));
         }
         else if(fProcessNames[nproc] == "WJets") {
           if(leptonName == "electron")
               histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_WJets_electron_clone").c_str()));
            else if (leptonName == "muon")
               histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_WJets_muon_clone").c_str()));
            else
               histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_WJets").c_str()));
         }
         */
         //else histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_" + fProcessNames[nproc]+ "_" + leptonName + "_clone").c_str()));
         histsToScale.push_back((TH1D*)inputMCList->FindObject((prefix + "_" + fProcessNames[nproc]+ "_" + leptonName + "_clone").c_str()));

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
   signalHistogram->Write();
   backgroundHistogram->Write();
   
   cout << "DONE" << endl;

   rootOutFile.Close();
   rootInFile->Close();
}

void Fitter::addSigBkgHistograms(vector<string> sig) {
   signalHistogram = new TH1D("signalHistogram","signalHistograms",
                                                monteCarloHistograms.begin()->second->GetNbinsX(),
                                                monteCarloHistograms.begin()->second->GetXaxis()->GetXmin(),
                                                monteCarloHistograms.begin()->second->GetXaxis()->GetXmax());
   signalHistogram->Sumw2();
   backgroundHistogram = new TH1D("backgroundHistogram","backgroundHistogram",
                                                    monteCarloHistograms.begin()->second->GetNbinsX(),
                                                    monteCarloHistograms.begin()->second->GetXaxis()->GetXmin(),
                                                    monteCarloHistograms.begin()->second->GetXaxis()->GetXmax());
   backgroundHistogram->Sumw2();

   for(map<string,TH1D*>::iterator it=monteCarloHistograms.begin(); it!=monteCarloHistograms.end(); it++) {
      if (DV.vfind(sig,(*it).first)!=-1) {
         signalHistogram->Add((*it).second);
         cout << "Fitter::Added the " << (*it).first << " histogram to the signalHistogram." << endl;
         cout << "\t" << setw(26) << (*it).first << " Entries: " << (*it).second->GetEntries() << endl;
         cout << "\t" << setw(26) << "signalHistogram Entries: " << signalHistogram->GetEntries() << endl;
      }
      else {
         backgroundHistogram->Add((*it).second);
         cout << "Fitter::Added the " << (*it).first << " histogram to the backgroundHistogram." << endl;
         cout << "\t" << setw(26) << (*it).first << " Entries: " << (*it).second->GetEntries() << endl;
         cout << "\t" << setw(26) << "backgroundHistogram Entries: " << backgroundHistogram->GetEntries() << endl; 
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
   if(signalHistogram->GetEntries()==0 || backgroundHistogram->GetEntries()==0) {
      cout << "Fitter::WARNING Either the signal or background histogram has zero entries!" << endl;
      cout << "\t" << setw(22) << "Signal Entries: " << signalHistogram->GetEntries() << endl;    
      cout << "\t" << setw(22) << "Background Entries: " << backgroundHistogram->GetEntries() << endl;
   }

   // set signal bins to zero if background bin is zero
   // set background bins to zero if signal bin is zero
   for(int i=1; i<=backgroundHistogram->GetNbinsX(); i++) {
      if (signalHistogram->GetBinContent(i)>0 && backgroundHistogram->GetBinContent(i)==0) {
         signalHistogram->SetBinContent(i,0.0);
         signalHistogram->SetBinError(i,0.0);
      }
   }

   if (FOM==1) {
      double fom1 = 0.0;
      for(int i=1; i<=backgroundHistogram->GetNbinsX(); i++)
         fom1 += TMath::Power((signalHistogram->GetBinContent(i)/TMath::Sqrt(backgroundHistogram->GetBinContent(i))),2);
      return TMath::Sqrt(fom1);
   }
   else if (FOM==2) {
      double fom2 = 0.0;
      for(int i=1; i<=backgroundHistogram->GetNbinsX(); i++)
         fom2 += TMath::Power((signalHistogram->GetBinContent(i)/TMath::Sqrt(signalHistogram->GetBinContent(i)+backgroundHistogram->GetBinContent(i))),2);
      return TMath::Sqrt(fom2);
   }
   else if (FOM==2.5)
      return FigureOfMerit::usingChi2(signalHistogram,
                                      backgroundHistogram,0.0001);
   else if (FOM==3)
      return FigureOfMerit::usingShapeFromTemplates(signalHistogram,
                                                    backgroundHistogram,0.0001);
   else {
      cout << "Fitter::WARNING The FOM type was not specified. Using defaul (FOM3)." << endl;
      return FigureOfMerit::usingShapeFromTemplates(signalHistogram,
                                                    backgroundHistogram,0.0001);
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
   if(TString(leptonName).CompareTo("electron",TString::kIgnoreCase)==0 ||
      TString(leptonName).CompareTo("both",TString::kIgnoreCase)==0) {
      //histNames.push_back("QCD_ElEnriched");
      histNames.push_back("QCD_ElFULL");
   }
   else if(TString(leptonName).CompareTo("muon",TString::kIgnoreCase)==0 ||
           TString(leptonName).CompareTo("both",TString::kIgnoreCase)==0) {
      histNames.push_back("QCD_MuFULL");
   }
   histNames.push_back("ggH125");
   histNames.push_back("qqH125");
   histNames.push_back("WH125_HToBB");
}

vector<double> Fitter::fitMinimization()
{
   if(minimizationTechnique == "chi2") {
      funcFit = new ROOT::Math::Functor(this,&Fitter::fitFunc,2);
   }
   else if(minimizationTechnique == "KS") {
      funcFit = new ROOT::Math::Functor(this,&Fitter::fitKSFunc,2);
   }
   minFit = new ROOT::Minuit2::Minuit2Minimizer( ROOT::Minuit2::kMigrad );
   
   // Tolerance and printouts
   minFit->SetPrintLevel(3);
   //strategy = (0,1,2)
   //0 - Don't wast any function calls to determine where you are. When there are a lot of variables and precision is not as important.
   //1 - Intermediate (default)
   //2 - When the function is evaluated quickly and precision is important.
   minFit->SetStrategy(2);
   minFit->SetMaxFunctionCalls(100000);
   minFit->SetMaxIterations(100000);
   //EDM must be less than 0.001*tolerance*up for minimization to be successful
   minFit->SetTolerance(0.001);
   //1.0 for chi-square fit and 1 stdev errors
   //4.0 for chi-square fit and 2 stdev errors
   //0.5 for negative-log-likelihood function
   //Value in the minimization function used to compute the parameter errors. The default is to get the uncertainties
   //at the 68% CL is a value of 1 for a chi-squared function minimization and 0.5 for a log-likelihood function.
   if(minimizationTechnique == "chi2") {
      minFit->SetErrorDef(1.0);
   }
   else if(minimizationTechnique == "LL") {
      minFit->SetErrorDef(0.5);
   }
   minFit->SetValidError(true);
   
   // Fitting
   minFit->SetFunction(*funcFit);
   
   int iproc = 0;
   for(map<string, TH1D*>::iterator it = monteCarloHistograms.begin(); it != monteCarloHistograms.end(); it++) {
      if(DefaultValues::vfind(fProcessNames,it->first)>-1 && (it->first=="QCD_ElFULL" || it->first=="QCD_MuFULL")) {
         minFit->SetLowerLimitedVariable(iproc, it->first.c_str(), 0.364, 0.000001, 0.0);
         iproc++;
      }
      else if(DefaultValues::vfind(fProcessNames,it->first)>-1){
         minFit->SetLowerLimitedVariable(iproc, it->first.c_str(), 1.0, 0.000001, 0.0);
         iproc++;
      }
      else {
         continue;
      }
   }
   //double parameter[2];
   //if(DefaultValues::vfind(fProcessNames,"QCD")<0) {
   //   parameter[0] = 1.0;
   //   parameter[1] = 1.0;
   //}
   //else {
   //   parameter[0] = 1.0;
   //   parameter[1] =  0.0001;
   //}
   //double step[2] = {0.0001, 0.0001};
   //double lower[2] = {0.0, 0.0};
   
   //minFit->SetVariable(0, "Par0", parameter[0], step[0]);
   //minFit->SetVariable(1, "Par1", parameter[1], step[1]);
   //minFit->SetLowerLimitedVariable(0, "Par0", parameter[0], step[0], lower[0]);
   //minFit->SetLowerLimitedVariable(1, "Par1", parameter[1], step[1], lower[1]);
   
   cout << endl << "Beginning Fit..." << endl;
   if(minFit->Minimize()) {
      cout << endl << "Fit Successful!" << endl;
      //minFit->PrintResults();
      reducedChiSquared = minFit->MinValue() / (dataHistogram->GetNbinsX() - 2);
      cout << "Chi2\t  = " << minFit->MinValue() << endl;
      cout << "NDF\t  = " << dataHistogram->GetNbinsX() - 2 << endl;
      cout << "Chi2/NDF\t  = " << reducedChiSquared << endl;

      for(unsigned int i=0; i<minFit->NDim(); i++) {
         scaleParameters[i] = minFit->X()[i];
      }
   }
   else {
      cout << endl << "Fit Failed." << endl;
      for(unsigned int i=0; i<minFit->NDim(); i++) {
         scaleParameters[i] = 1.0;
      }
   }

   return scaleParameters;
}

vector<double> Fitter::fitTFractionFitter()
{
   TObjArray *mc = new TObjArray(monteCarloHistograms.size());
   for(map<string, TH1D*>::iterator it = monteCarloHistograms.begin(); it != monteCarloHistograms.end(); it++) {
      mc->Add(it->second);
   }
   TFractionFitter* fit = new TFractionFitter(dataHistogram, mc,"V");
   TVirtualFitter* vFit = fit->GetFitter();
   Foption_t fopt;
   fopt.StoreResult = 1;
   fopt.Like = 1;
   vFit->SetFitOption(fopt);

   int iproc = 0;
   for(map<string, TH1D*>::iterator it = monteCarloHistograms.begin(); it != monteCarloHistograms.end(); it++) {
      bool notFProc = true;
      for (unsigned int nproc=0; nproc<fProcessNames.size(); nproc++) {
         if (it->first == fProcessNames[nproc])
            notFProc = false;
      }
      if(notFProc) {
         vFit->SetParameter(iproc,it->first.c_str(),it->second->Integral()/dataHistogram->Integral(),0.0,0.0,1.0);
         vFit->FixParameter(iproc);
      }
      else if(!notFProc && (it->first=="QCD_ElFULL" || it->first=="QCD_MuFULL")) {
         vFit->SetParameter(iproc,it->first.c_str(),0.0001,0.0001,0.0,2.0);
      }
      else if(!notFProc && it->first=="WJets") {
         vFit->SetParameter(iproc,it->first.c_str(),1.0,0.0001,0.0,2.0);
      }
      iproc++;
   }

   cout << endl << "Beginning Fit..." << endl;
   TFitResultPtr status = fit->Fit();
   if((int)status == 0) {
      cout << endl << "Fit Successful!" << endl;
      TH1F* result = (TH1F*) fit->GetPlot();
      //dataHistogram->Draw("Ep");
      //result->Draw("same");
      //cin.get();
      result->Write();
   }
   else
      cout << endl << "Fit Failed." << endl;
   
   cout << endl << "##### FIT RESULTS #####" << endl;
   //vFit->PrintResults();
   status->Print();
   
   for (unsigned int nproc=0; nproc<fProcessNames.size(); nproc++) {
      if(fProcessNames[nproc]=="QCD_ElFULL" || fProcessNames[nproc]=="QCD_MuFULL")
         cout << "Scale factor (" << fProcessNames[nproc] << "): " << fit->GetMCPrediction(1)->Integral()/monteCarloHistograms[fProcessNames[nproc]]->Integral() << endl;
      else
         cout << "Scale factor (" << fProcessNames[nproc] << "): " << fit->GetMCPrediction(0)->Integral()/monteCarloHistograms[fProcessNames[nproc]]->Integral() << endl;
   }

   reducedChiSquared = fit->GetChisquare() / fit->GetNDF();
   cout << "Chi2\t  = " << fit->GetChisquare() << endl;
   cout << "NDF\t  = " << fit->GetNDF() << endl;
   cout << "Chi2/NDF\t  = " << reducedChiSquared << endl << endl;

   double error = 0;
   for(unsigned int i = 0; i<scaleParameters.size(); i++) {
      fit->GetResult(i, scaleParameters[i], error);
   }
   
   return scaleParameters;
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
      i = dataHistogram->Integral(0, dataHistogram->GetNbinsX()+1);
   else
      i = dataHistogram->Integral(1, dataHistogram->GetNbinsX());

   cout << "Data Integral: "
        << i << endl;
}

void Fitter::printMCIntegrals(bool includeOverflowUnderflow) {
   for(map<string, TH1D*>::iterator mapit = monteCarloHistograms.begin(); mapit != monteCarloHistograms.end(); mapit++)
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
      i = signalHistogram->Integral(0, signalHistogram->GetNbinsX()+1);
   else
      i = signalHistogram->Integral(1, signalHistogram->GetNbinsX());

   cout << endl <<"Signal Integral: " 
        << i << endl;
}

void Fitter::printBackgroundIntegral(bool includeOverflowUnderflow) {
   double i = 0;
   if(includeOverflowUnderflow)
      i =  backgroundHistogram->Integral(0, backgroundHistogram->GetNbinsX()+1);
   else
      i =  backgroundHistogram->Integral(1, backgroundHistogram->GetNbinsX());

   cout << endl <<"Background Integral: " 
        << i << endl;
}

double Fitter::getMCSum(bool includeOverflowUnderflow) {
   double sum = 0;
   for(map<string, TH1D*>::iterator mapit = monteCarloHistograms.begin(); mapit != monteCarloHistograms.end(); mapit++)
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
   dataHistogram->GetNbinsX();
   double chiSquare = 0;

   for(int bin = 1; bin <= dataHistogram->GetNbinsX(); bin++)
   {
      double mc = 0;
      double mcError2 = 0;
      double data = 0;
      double dataError2 = 0;

      for (unsigned int nproc=0; nproc<fProcessNames.size(); nproc++) {
         mc += par[nproc]*monteCarloHistograms[fProcessNames[nproc]]->GetBinContent(bin);
         mcError2 += pow(par[nproc]*monteCarloHistograms[fProcessNames[nproc]]->GetBinError(bin), 2);
      }

      data += dataHistogram->GetBinContent(bin);
      dataError2 += pow(dataHistogram->GetBinError(bin), 2);

      for(map<string, TH1D*>::iterator it = monteCarloHistograms.begin(); it != monteCarloHistograms.end(); it++)
      {
         bool notFProc = true;
         for (unsigned int nproc=0; nproc<fProcessNames.size(); nproc++) {
            if ((*it).first == fProcessNames[nproc])
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

   //Add Gaussian Constraints
   for(unsigned int nproc=0; nproc<fProcessNames.size(); nproc++) {
      if(fProcessXsec[nproc]!=0 && fProcessXsecError[nproc]!=0) {
         chiSquare+=TMath::Power((par[nproc]-1.0)/(fProcessXsecError[nproc]/fProcessXsec[nproc]),2);
      }
   }

   return chiSquare;
}

double Fitter::fitKSFunc(const double *par) {
   TH1D* mc = (TH1D*)dataHistogram->Clone("mc");
   mc->Reset();
   int iproc = 0;
   for(map<string, TH1D*>::iterator it = monteCarloHistograms.begin(); it != monteCarloHistograms.end(); it++) {
      if(DefaultValues::vfind(fProcessNames,it->first)>-1) {
         //cout << it->first << "\tint_before=";
         TH1D* tmp = (TH1D*)it->second->Clone((it->first+"_clone").c_str());
         //cout << tmp->Integral() << "\tscale_factor=" << par[iproc] << "\tint_after=";
         tmp->Scale(par[iproc]);
         //cout << tmp->Integral() << endl;
         mc->Add(tmp);
         iproc++;
      }
      else {
         mc->Add(it->second);
      }
   }
   //cout << "Data=" << dataHistogram->Integral() << "\tmc=" << mc->Integral() << endl;
   return -dataHistogram->KolmogorovTest(mc);
}
