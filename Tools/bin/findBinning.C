////////////////////////////////////////////////////////////////////////////////
//
// findBinning
// -----------
//
//                          05/16/2015 Alexx Perloff <aperloff@Physics.tamu.edu>
////////////////////////////////////////////////////////////////////////////////

//ROOT Libraries
#include "TROOT.h"
#include "TSystem.h"
#include "TObject.h"
#include "TClass.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "THStack.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TArrayD.h"
#include "TString.h"

//Standard Template Libraries
#include <iostream>
#include <iomanip>
#include <vector>

//TAMUWW Includes
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// declare global variables
////////////////////////////////////////////////////////////////////////////////
TString signalNames[] = {"ggH125","qqH125","WH125_HToBB","WH125_HToBB_M125","WH_ZH_TTH_HToWW","WH_ZH_TTH_HToZZ","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125"/*,"WH_ZH_TTH_HToWW_M125"*/,"WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125"/*,"WH_ZH_TTH_HToZZ_M125"*/};
TString backgroundNames[] = {"TTH_HToBB","TTH_HToBB_M125","WJets","ZJets","TTbar","Diboson","singleTop","QCD","QCD_ElFULL","QCD_MuFULL","WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","WW"};

////////////////////////////////////////////////////////////////////////////////
// declare local functions
////////////////////////////////////////////////////////////////////////////////
bool checkKey(TClass* cl, TString name, TString prefix, TString suffix, TString jetBin, bool fromMC = false, bool verbose = false);

bool checkSums(TH1F* signal, TH1F* background, TH1F* data, bool verbose = false);

bool checkMergeConditions(TH1F* signal, TH1F* background, TH1F* data, int ibin);

void setBinEdges(TH1F* hist, TArrayD* binEdges, int binMerge = -1);

TH1F* rebinHistogram(TH1F* hist, TArrayD* binEdges);

void printBinning(TArrayD* binEdges);

////////////////////////////////////////////////////////////////////////////////
// main
////////////////////////////////////////////////////////////////////////////////

//Example commands:
// .L findBinning.C+
// findBinning("/uscms/home/jdg6y/Files_forAlexx/FinelyBinnedBDTHists/VarJetBins_2Jets_Scaled.root","","","hist_","",false)
// findBinning("/uscms/home/aperloff/nobackup/Summer12ME8TeV/2015_03_06_ValidationPlots_PU_CSV_TTbar_QCDEta_Scale/Jets2/electron/histos_electron.root","","","KinBDT_","_electron",false)

//______________________________________________________________________________
void findBinning(TString filename, TString idir = "", TString jetBin = "",
                 TString prefix = "", TString suffix = "", bool sumMC = false,
                 bool fromMC = false, bool verbose = false) {

	//Open the file containing the histograms to be analyzed
	TFile* ifile = TFile::Open(filename,"READ");

	//Move the gDirectory to the input directory, assuming is need to be moved
	if (!idir.IsNull()) {
		ifile->cd(idir);
	}
	
	//Get a list of the histograms, in a given directory, which will be analyzed
	//Create an iterator for that list
	//Create a key pointer to access each key in the list
	//Create a pointer to the object which is being accessed
	TIter next(gDirectory->GetListOfKeys());
	TKey *key;
	TH1F* h;

	//Create the signal and background stacks as well as vectors containing the names of the processes for each
	TH1F* signalSum;
	TH1F* backgroundSum;
	TH1F* dataSum;
	vector<TString> sNames(signalNames, signalNames + sizeof(signalNames) / sizeof(TString));
	vector<TString> bNames(backgroundNames, backgroundNames + sizeof(backgroundNames) / sizeof(TString));

	//Add in the prefixes and suffixes for the histogram names
	for(unsigned int i=0; i<sNames.size(); i++) {
		sNames[i] = prefix+sNames[i]+suffix;
	}
	for(unsigned int i=0; i<bNames.size(); i++) {
		bNames[i] = prefix+bNames[i]+suffix;
	}

	//Loop through the list of histograms
	TString name;
	bool isFirst = true;
	while((key=(TKey*)next())) {

		//Store some information about the next key
        TClass *cl = gROOT->GetClass(key->GetClassName());
		name = key->GetName();

		//Do some sanity checks
		if(!checkKey(cl,name,prefix,suffix,jetBin,fromMC,verbose))
			continue;

		//Access the histogram object
		h = (TH1F*)key->ReadObj();

		//If this is the first histogram, use it to clone into the sum histograms
		if(isFirst) {
			if(verbose) cout << "Cloning into the signalSum and backgroundSum histograms ... ";
			signalSum = (TH1F*)h->Clone("signalSum");
			signalSum->Reset();
			signalSum->SetNameTitle("signalSum","signalSum");
			signalSum->Sumw2();
			backgroundSum = (TH1F*)h->Clone("backgroundSum");
			backgroundSum->Reset();
			backgroundSum->SetNameTitle("backgroundSum","backgroundSum");
			backgroundSum->Sumw2();
			dataSum = (TH1F*)h->Clone("dataSum");
			dataSum->Reset();
			dataSum->SetNameTitle("dataSum","dataSum");
			dataSum->Sumw2();
			isFirst = false;
			if(verbose) cout << "DONE" << endl;
		}

		if(DefaultValues::vfind(sNames,name)>=0) {
			cout << "Adding " << name << " to the signal sum." << endl;
			signalSum->Add(h);
		}
		else if(DefaultValues::vfind(bNames,name)>=0) {
			cout << "Adding " << name << " to the background sum." << endl;
			backgroundSum->Add(h);

			//Basically checks if the files are coming from Joey or Alexx
            //Only use sum of MC if sumMC flag is set. Otherwise use a different histogram
			if(prefix.CompareTo("hist_")==0 && sumMC) {
               cout << "Adding " << name << " to the data sum." << endl;
               dataSum->Add(h);
            }
		}
		//Don't double count the data if the files are coming from Joey
		//else if(name.Contains("data",TString::kIgnoreCase) && !prefix.CompareTo("hist_")==0) {
        else if(name.Contains("data",TString::kIgnoreCase)) {
			cout << "Adding " << name << " to the data sum." << endl;
			dataSum->Add(h);
		}
		else {
			cout << "ERROR::Unknown category (signal,background) for histogram " << name << endl
				 << "Please fix and try again" << endl;
		    return;
		}
	}

	//Do some sanity checks on the sums
	cout << endl << "Doing some sanity checks on the histogram sums ... ";
	if(!checkSums(signalSum,backgroundSum,dataSum,verbose))
		return;
	cout << "DONE" << endl;

	//Format the histograms
	signalSum->SetFillStyle(1001);
	signalSum->SetFillColor(kGreen);
	signalSum->SetLineColor(kGreen);
	backgroundSum->SetFillStyle(1001);
	backgroundSum->SetFillColor(kRed);
	backgroundSum->SetLineColor(kRed);
	dataSum->SetFillStyle(0);
	dataSum->SetFillColor(kBlack);
	dataSum->SetLineColor(kBlack);
	dataSum->SetMarkerStyle(8);

	//Create the canvas and draw the original sums
	TCanvas* c = new TCanvas("c","c",800,400);
	c->Divide(2,1);
	c->cd(1)->SetLogy();
	THStack* originalStack = new THStack("originalStack","originalStack");
	originalStack->Add((TH1F*)backgroundSum->Clone("backgroundSumOriginal"),"HIST");
	originalStack->Add((TH1F*)signalSum->Clone("signalSumOriginal"),"HIST");
	originalStack->Draw();
	((TH1F*)dataSum->Clone("dataSumOriginal"))->Draw("same");
	TLegend* originalLegend = new TLegend(0.6,0.7,0.9,0.9);
	originalLegend->AddEntry(backgroundSum,"Sum of backgrounds","f");
	originalLegend->AddEntry(signalSum,"Sum of signals","f");
	originalLegend->AddEntry(dataSum,"Sum of the data","p");
	originalLegend->Draw("same");

    int mergeCounter = 0;
	TArrayD* currentBinning = new TArrayD();
	TArrayD* resetBinning = new TArrayD();
	cout << endl << "Starting to loop through the background sum bins ... " << endl; 
	for(int ibin=1; ibin<=backgroundSum->GetXaxis()->GetNbins(); ibin++) {
		cout << "\tDoing bin " << ibin << " ... "; 

		//Set the current lower/upper bin edges
		setBinEdges(backgroundSum,currentBinning);
		if(verbose) printBinning(currentBinning);
		
		//Conditions for when to merge bins
		if(ibin!=backgroundSum->GetXaxis()->GetNbins() && checkMergeConditions(signalSum,backgroundSum,dataSum,ibin)) {

			cout << "Merging bin " << ibin << " with bin " << ibin+1 << " which makes a range of "
				 << currentBinning->At(ibin-1) << " to ";
            mergeCounter++;
			if(ibin+1>=currentBinning->GetSize())
				cout << currentBinning->At(currentBinning->GetSize()-1) << " ... ";
			else
				cout << currentBinning->At(ibin+1) << " ... ";

			//Merge the correct bins (i.e. ibin with ibin+1)
			setBinEdges(backgroundSum,resetBinning,ibin);
			if(verbose) printBinning(resetBinning);

			//Now reset the background, signal, and data histograms to have the new binning
			backgroundSum = rebinHistogram(backgroundSum,resetBinning);
			signalSum     = rebinHistogram(signalSum,resetBinning);
			dataSum       = rebinHistogram(dataSum,resetBinning);

			//Make sure when you go through the loop that the bin number (ibin) doesn't increment
			//But don't do this on the last bin
			ibin--;
		}
		//If you are on the last bin, then check if you need to be merged with the previous bin
		//You won't need to continue to loop because by merging the last bin with a bin that didn't
		//    need to be merged you are guaranteed to satisfy the merge conditions.			else {
		else if(ibin==backgroundSum->GetXaxis()->GetNbins() && checkMergeConditions(signalSum,backgroundSum,dataSum,ibin)) {
			cout << "SPECIAL::Merging bin " << ibin-1 << " with bin " << ibin 
				 << " which makes a range of " << currentBinning->At(ibin-2) << " to ";
            mergeCounter++;
			if(ibin+1>=currentBinning->GetSize())
				cout << currentBinning->At(currentBinning->GetSize()-1) << " ... ";
			else
				cout << currentBinning->At(ibin+1) << " ... ";

			//Merge the final two bins
			setBinEdges(backgroundSum,resetBinning,ibin-1);

			//Rebin the histograms
			backgroundSum = rebinHistogram(backgroundSum,resetBinning);
			signalSum     = rebinHistogram(signalSum,resetBinning);
			dataSum       = rebinHistogram(dataSum,resetBinning);
		}

		cout << "DONE" << endl;
	}

	//Find the final binning
	setBinEdges(backgroundSum,currentBinning);

	//Print the final binning
	cout << "////////////////////////////////////////////////////////////////////////////////" << endl
		 << "// the final binning (" << mergeCounter << " merges)" << endl
		 << "////////////////////////////////////////////////////////////////////////////////" << endl;
	printBinning(currentBinning);

	//Draw the final histograms
	c->cd(2)->SetLogy();
	THStack* finalStack = new THStack("finalStack","finalStack");
	finalStack->Add((TH1F*)backgroundSum->Clone("backgroundSumFinal"),"HIST");
	finalStack->Add((TH1F*)signalSum->Clone("signalSumFinal"),"HIST");
	finalStack->Draw();
	((TH1F*)dataSum->Clone("dataSumFinal"))->Draw("same");
	TLegend* finalLegend = new TLegend(0.6,0.7,0.9,0.9);
	finalLegend->AddEntry(backgroundSum,"Sum of backgrounds","f");
	finalLegend->AddEntry(signalSum,"Sum of signals","f");
	finalLegend->AddEntry(dataSum,"Sum of data","p");
	finalLegend->Draw("same");

	return;
}

////////////////////////////////////////////////////////////////////////////////
// implement local functions
////////////////////////////////////////////////////////////////////////////////
//______________________________________________________________________________
bool checkKey(TClass* cl, TString name, TString prefix, TString suffix, TString jetBin, bool fromMC, bool verbose) {
	TString className = cl->ClassName();

	//Check that the key is a histogram
	if(!cl->InheritsFrom("TH1")) {
		if(verbose) cout << "WARNING::The class name (" << className << ") of " << name 
		   				 << " does not inherit from a TH1 and will be skipped." << endl;
		return false;
	}

	//Check that the histogram is 1D
	if(!className.Contains(1) && (className.Contains(2)||className.Contains(3)) ) {
		if(verbose) cout << "WARNING::The histogram " << name << " is not of dimension 1 and will be skipped." << endl;
		return false;
	}

	//Check that the histogram name does not contain any of the systematics histograms
	if(name.Contains("JESDown",TString::kIgnoreCase) || name.Contains("JESUp",TString::kIgnoreCase)
       || name.Contains("JECDown",TString::kIgnoreCase) || name.Contains("JECUp",TString::kIgnoreCase)
	   || name.Contains("EtaWUp",TString::kIgnoreCase) || name.Contains("EtaWDown",TString::kIgnoreCase)
	   || name.Contains("matchingup",TString::kIgnoreCase) || name.Contains("matchingdown",TString::kIgnoreCase)
	   || name.Contains("scaleup",TString::kIgnoreCase) || name.Contains("scaledown",TString::kIgnoreCase) 
	   || name.Contains("PtUp",TString::kIgnoreCase) || name.Contains("PtDown",TString::kIgnoreCase)       || name.Contains("CMS_hww_lnujj",TString::kIgnoreCase) || name.Contains("CMS_scale_j",TString::kIgnoreCase)) {
		if(verbose) cout << "WARNING::The histogram " << name << " contains data or is a systematic and will be skipped." << endl;
		return false;
	}

	//Check that the histogram name contains the prefix and suffix
	if(!name.BeginsWith(prefix,TString::kIgnoreCase) || !name.Contains(suffix,TString::kIgnoreCase)) {
       //if((!fromMC && !name.Contains("data_Electron",TString::kIgnoreCase) && !name.Contains("data_Muon",TString::kIgnoreCase)) || (fromMC && !name.Contains("data_fromMC",TString::kIgnoreCase))) {
          if(verbose) cout << "WARNING::The histogram " << name << " does not contain the prefix (" << prefix
                           << ") or suffix (" << suffix << ") and will be skipped." << endl;
          return false;
          //}
	}

	//Check that the jet bin names are all the same
	//If the user inputs a file with more than one jet bin, then they must specify
	//    a jetBin name as an input parameter.
	if(!jetBin.IsNull() && !name.Contains(jetBin)) {
		if(verbose) cout << "WARNING::The histogram " << name << " is not of the correct jet bin and will be skipped." << endl;
		return false;
	}

	return true;
}

//______________________________________________________________________________
bool checkSums(TH1F* signal, TH1F* background, TH1F* data, bool verbose) {
	//Check that the histograms exist
    if(!signal) {
		if(verbose) cout << endl << "ERROR::The signal histogram does not exist." << endl;
		return false;
	}
	if(!background) {
		if(verbose) cout << endl << "ERROR::The background histogram does not exist." << endl;
		return false;
	}
	if(!data) {
		if(verbose) cout << endl << "ERROR::The data histogram does not exist." << endl;
		return false;
	}

	//Check that the histograms have at least some bins
	if(signal->GetXaxis()->GetNbins()==0 || background->GetXaxis()->GetNbins()==0 || data->GetXaxis()->GetNbins()==0) {
		if(verbose) cout << endl << "ERROR::The signal, background and/or data histograms do not have any bins." << endl;
		return false;
	}

	//Check if the histograms have the same number of bins
	if(signal->GetXaxis()->GetNbins()!=background->GetXaxis()->GetNbins()
	   || signal->GetXaxis()->GetNbins()!=data->GetXaxis()->GetNbins()
	   || background->GetXaxis()->GetNbins()!=data->GetXaxis()->GetNbins()) {
		if(verbose) cout << endl << "ERROR::The signal, background and/or data histograms do not have the same number of bins." << endl;
		return false;
	}

	//Check if the histograms have the same upper and lower limits
	if(signal->GetXaxis()->GetXmin()!=background->GetXaxis()->GetXmin() ||
	   signal->GetXaxis()->GetXmax()!=background->GetXaxis()->GetXmax() ||
	   signal->GetXaxis()->GetXmin()!=data->GetXaxis()->GetXmin() ||
	   signal->GetXaxis()->GetXmax()!=data->GetXaxis()->GetXmax()) {
		if(verbose) cout << endl << "ERROR:: The signal, background and/or data histograms do not have the same upper and lower bounds." << endl;
		return false;
	}

    //Check that the integrals of all of the histograms are greater than zero
    if(signal->Integral()==0 || background->Integral()==0 || data->Integral()==0) {
       if(verbose) cout << endl << "ERROR:: The signal, background, and/or data histograms have an integral equal to zero." << endl; 
       return false;
    }

	return true;
}

//______________________________________________________________________________
bool checkMergeConditions(TH1F* signal, TH1F* background, TH1F* data, int ibin) {
	if((background->GetBinContent(ibin)>0 && background->GetBinError(ibin)/background->GetBinContent(ibin)>0.1)
		|| (signal->GetBinContent(ibin)>0 && background->GetBinContent(ibin)==0)
		|| (background->GetBinContent(ibin)==0 && signal->GetBinContent(ibin)==0)
        || (data->GetBinContent(ibin)==0) 
        || (data->GetBinError(ibin)/data->GetBinContent(ibin)>0.1) ) {
       //cout << endl << background->GetBinContent(ibin) << " +\- "<< background->GetBinError(ibin) << " = " << background->GetBinError(ibin)/background->GetBinContent(ibin) << endl;
		return true;
	}
    //cout << endl << background->GetBinContent(ibin) << " +\- "<< background->GetBinError(ibin) << " = " << background->GetBinError(ibin)/background->GetBinContent(ibin) << endl;
	return false;
}

//______________________________________________________________________________
void setBinEdges(TH1F* hist, TArrayD* binEdges, int binMerge) {
	binEdges->Reset();
	if(binMerge>-1)
		binEdges->Set(hist->GetNbinsX());
	else
		binEdges->Set(hist->GetNbinsX()+1);

	for(int ibin=1; ibin<=hist->GetNbinsX(); ibin++) {
		if(binMerge>-1 && ibin==binMerge+1) {
			continue;
		}
		else if(binMerge>-1 && ibin>binMerge+1)
			binEdges->SetAt(hist->GetBinLowEdge(ibin),ibin-2);
		else
			binEdges->SetAt(hist->GetBinLowEdge(ibin),ibin-1);
	}
	//Set the upper edge of the last bin
	if(binMerge>-1)
		binEdges->SetAt(hist->GetBinLowEdge(hist->GetNbinsX())+hist->GetBinWidth(hist->GetNbinsX()),hist->GetNbinsX()-1);
	else
		binEdges->SetAt(hist->GetBinLowEdge(hist->GetNbinsX())+hist->GetBinWidth(hist->GetNbinsX()),hist->GetNbinsX());
}

//______________________________________________________________________________
TH1F* rebinHistogram(TH1F* hist, TArrayD* binEdges) {
	return (TH1F*)hist->Rebin(binEdges->GetSize()-1,hist->GetName(),binEdges->GetArray());	
}

//______________________________________________________________________________
void printBinning(TArrayD* binEdges) {
	//Double_t xbins[25] = {...}
	cout << "Double_t xbins[" << binEdges->GetSize() << "] = {";
	for(int i=0; i<binEdges->GetSize(); i++) {
		cout << binEdges->At(i);
		if(i!=binEdges->GetSize()-1)
			cout << ",";
	}
	cout << "}" << endl;
}



