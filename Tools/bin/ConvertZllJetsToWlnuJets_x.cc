//
// User Defined Includes
//
#include "TAMUWW/MEPATNtuple/interface/EventNtuple.hh"
#include "JetMETAnalysis/JetUtilities/interface/CommandLine.h"
#include "TAMUWW/SpecialTools/interface/Defs.hh"
#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"
#include "TAMUWW/SpecialTools/interface/PhysicsProcess.hh"
#include "TAMUWW/SpecialTools/interface/ProgressBar.hh"
// Might want to use this to modify the lepton/MET to convert from Z (91.2 GeV) to W (80.4 GeV)
//http://cmslxr.fnal.gov/source/PhysicsTools/RecoUtils/src/CandMassKinFitter.cc?v=CMSSW_6_1_1

//
// CMSSW Includes
//
#include "PhysicsTools/RecoUtils/interface/CandMassKinFitter.h"
#include "PhysicsTools/KinFitter/interface/TKinFitter.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintM.h"
#include "PhysicsTools/KinFitter/interface/TFitConstraintEp.h"
#include "PhysicsTools/KinFitter/interface/TAbsFitParticle.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEMomDev.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleMCCart.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtEtaPhi.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEtThetaPhi.h"
#include "PhysicsTools/KinFitter/interface/TFitParticleEScaledMomDev.h"

//
// ROOT includes
//
#include "TROOT.h"
#include "TSystem.h"
#include "TEnv.h"
#include "TObject.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile3D.h"
#include "TRandom3.h"
#include "TString.h"
#include "TF1.h"
#include "TMath.h"
#include "TBenchmark.h"
#include "TMatrixD.h"

//
// Standard Library Includes
//
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <cassert>
#include <utility>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <exception>
#include <chrono>
#include <random> // for std::random_device and std::mt19937
#include <stdlib.h> // system, NULL, EXIT_FAILURE
#include <memory> // for std::shared_pointer

using namespace std;
using namespace reco;

////////////////////////////////////////////////////////////////////////////////
//  Typedefs
////////////////////////////////////////////////////////////////////////////////
typedef pair<map<int,int>,map<int,int> > pair_map;

////////////////////////////////////////////////////////////////////////////////
//  Declare Global Variables
////////////////////////////////////////////////////////////////////////////////

const double ZMASS = 91.1876; //GeV
const double WMASS = 80.3980; //GeV
enum Constraint { kWMass = 1, kNeutrinoMass, kSumPt };

////////////////////////////////////////////////////////////////////////////////
//  Declare Local Functions
////////////////////////////////////////////////////////////////////////////////

// Copy the selected lepton to the refLV of the MET and the original METLV to the rawLV position
// This way we have access to this information later on when doing reweighting, studies, etc.
void saveOriginalInformation(EventNtuple *ntuple, int randIndex);

// Return the value of the MET resolution in the x and y directions
pair<double,double> getMETResolution(EventNtuple *ntuple, int randIndex, double sumETMult);

// Move a randomly selected lepton to the MET (add the Lorentz vectors) and remove that lepton 
//  from the lepton vector
// The higher that probIndex1 is the more likely it is that the lower pT lepton will be moved
//  to the MET Lorentz vector (i.e. the percentage of time that the higher pT lepton will be kept).
//  The first probIndex1 is for muons and the second is for electrons.
// Returns true is the lower pT lepton was moved to the MET (i.e. randIndex==1)
// For the MET resolution smearing you have three methods to choose from real, simple, and coupled.
//  The real one evaluates a function from an 8 TeV analysis to get a realistic estimation of the MET
//  resolution. The simple one multiplies the MET.Px(y) by a random number from a Gaussian distribution
//  centered around 1 with a sigma of METResolution. The coupled method uses a similar method from real,
//  but only randomly chooses Px and then calculates Py assuming the pT will stay the same.
bool moveLeptonToMET(EventNtuple *ntuple, pair<double,double> probIndex1, pair_map &leptonSelectionCounter,
                     string METResMethod, double sumETMult, double METResolution, double METPhiResolution,
                     bool resetMETMass, bool genOnly, bool debug = false);

// Change one of the gen lepton's pdgId to be that of a neutrino
// Change either the higher or lower pT gen particle based on what happened with the reconstructed leptons
void setGenLeptonToNeutrino(EventNtuple *ntuple, bool lowerPtLepton, bool debug = false);

// Rebalance the event using the TKinFitter class
// The event needs rebalancing because of the shift from Z (91.2 GeV) to W (80.4 GeV)
// Will need a modification of lepton/MET due to the change in momther particle mass
// Will need a modigication of the jets because they will be recoiling against something lighter
// Currently not working!!!
void rebalanceEvent(EventNtuple *ntuple, bool debug = false);

// Function to print the KinFitter quantities
void print(TKinFitter *fitter);

// Simply scale the lepton and jet Lorentz vectors to account for Z->W (not the MET which comes from PU)
void scaleFourVectors(EventNtuple *ntuple, bool debug = false);

// Scale the V and daugter 4-vectors to account for Z->W, just as we did for the reconstructed quantities
void scaleGenVFourVectorAndDaugters(EventNtuple *ntuple, bool debug = false);

// If the oFilePath didn't point to an EOS directory, then move the output file to EOS
void moveOutputFile(string oFilePath, string oProcessName, string suffix);

// Counts the number of events which have been processed by using the leptonSelectionCounter
// Use the leptonCat from DEFS to decide if the sum should be for muons, electrons, or both
int leptonSelectionCounterSum(pair_map &leptonSelectionCounter, DEFS::LeptonCat leptonCat = DEFS::both);

////////////////////////////////////////////////////////////////////////////////
//  main
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
int main(int argc,char**argv)
{
    //
    // Setup the command line options and check that they make sense
    //
    CommandLine cl;
    if (!cl.parse(argc,argv)) return 0;

    string         iProcessName     = cl.getValue<string>  ("iProcessName", "ZJetsToLL_M50");
    string         oProcessName     = cl.getValue<string>  ("oProcessName",      "WlnuJets");
    string         suffix           = cl.getValue<string>  ("suffix",               "_M-50");
    string         oFilePath        = cl.getValue<string>  ("oFilePath",                 "");
    int            nEntriesMax      = cl.getValue<int>     ("nEntriesMax",               -1);
    vector<double> probIndex1       = cl.getVector<double> ("probIndex1",       "0.5:::0.5");
    double         sumETMult        = cl.getValue<double>  ("sumETMult",                1.0);
    double         METResolution    = cl.getValue<double>  ("METResolution",           -1.0);
    double         METPhiResolution = cl.getValue<double>  ("METPhiResolution",         0.0);
    string         METResMethod     = cl.getValue<string>  ("METResMethod",              "");
    bool           resetMETMass     = cl.getValue<bool>    ("resetMETMass",            true);
    bool           rebalanceEvt     = cl.getValue<bool>    ("rebalanceEvt",           false);
    bool           genOnly          = cl.getValue<bool>    ("genOnly",                false);
    bool           doMove           = cl.getValue<bool>    ("doMove",                  true);
    bool           debug            = cl.getValue<bool>    ("debug",                  false);
    bool           help             = cl.getValue<bool>    ("help",                   false);

    if (help) {cl.print(); return 0;}
    if (!cl.check())       return 0;
    cl.print();

    //
    // Setup the program benchmarking
    //
    TBenchmark* m_benchmark = new TBenchmark();
    m_benchmark->Reset();
    m_benchmark->Start("sample");

    //
    // Trying to speed up the code by allowing it to prefetch data
    //
    gEnv->SetValue("TFile.AsyncPrefetching", 1);

    //
    // Some basic checks on the settings, specifically that there are only two values in the probIndex1 vector,
    //  one for muons and one for electrons
    //
    if(probIndex1.size()!=2) {
        cout << "ERROR::ConvertZllJetsToWlnuJets::main There must be two and only two values in the probIndex1 vector"
             << ", one for muons and one for electrons." << endl;
        std::terminate();
    }

    //
    // Open the input physics process and the input file
    //
    map<DEFS::LeptonCat, Table> normTable;
    // get the table with the files location
    Table fileTable = DefaultValues::getFileLocationTable(DEFS::pretag);
    PhysicsProcess* iProcess = DefaultValues::getSingleProcess(DEFS::PhysicsProcess::getProcessType(iProcessName),DEFS::jets2,
                                                               normTable, fileTable, false, DEFS::EventNtuple);
    
    //
    // Get the chain from the input physics process and set the input ntuple
    //
    TChain * c = iProcess->chain;
    // Turn this on if you want to limit the number of branches saved
    //c->SetBranchStatus();
    EventNtuple * ntuple = 0;
    if (c->GetBranch("EvtTree")) {
       c->SetBranchAddress("EvtTree", &ntuple);
       ntuple = new EventNtuple();
    }
    else if (c->GetBranch("EvtNtuple")) {
       c->SetBranchAddress("EvtNtuple", &ntuple);
       ntuple = new EventNtuple();
    }
    else {
       cout << "\tConvertZllJetsToWlnuJets::main EvtTree and EvtNtuple branches not found." << endl
            << "\tThe program cannot continue." << endl;
       std::terminate();
    }

    //
    // Find the input/output file path+name
    //
    if(oFilePath.empty()) {
        oFilePath = string(c->GetCurrentFile()->GetName());
        oFilePath = oFilePath.substr(0, oFilePath.rfind("/")+1);
    }
    if(oFilePath.back()!='/') oFilePath += "/";
    oFilePath += oProcessName + suffix + ".root";

    //
    // Open the output file and setup the output tree/ntuple
    //
    TFile* ofile = TFile::Open(oFilePath.c_str(),"RECREATE");
    ofile->mkdir("PS");
    ofile->cd("PS");
    TTree* otree = (TChain*)c->CloneTree(0);
    TH1D* Mll = new TH1D("Mll","Mll (L*#sigma*BR*SF/init. evt. scaled)",200,0,1000);
    Mll->Sumw2();

    //
    // Set the number of events to run over
    //
    int nEntries = (nEntriesMax<0 || c->GetEntries()<nEntriesMax) ? c->GetEntries() : nEntriesMax;

    //
    // Keep track of which lepton is kept
    //
    pair_map leptonSelectionCounter;
    leptonSelectionCounter.first.insert(pair<int,int>(0,0));
    leptonSelectionCounter.first.insert(pair<int,int>(1,0));
    leptonSelectionCounter.second.insert(pair<int,int>(0,0));
    leptonSelectionCounter.second.insert(pair<int,int>(1,0));

    //
    // Loop over the events
    //
    for(int ientry=0; ientry<nEntries; ientry++) {
       if(!debug)
          ProgressBar::loadbar2(ientry+1,nEntries);
        //
        // Load the input ntuple
        //
        c->GetEntry(ientry);

        //
        // Some basic checks of the starting vectors
        //
        if(!genOnly && ntuple->jLV.size()==0 && ntuple->lLV.size()==0 && ntuple->METLV.size()==0) {
            cout << "ERROR::ConvertZllJetsToWlnuJets::main Something is wrong. If the jet, lepton, and MET Lorentz vectors "
                 << "all have size zero and the \"genOnly\" flag is not set then the ntuple must be faulty. Otherwise you may "
                 << "have forgotten to set the \"genOnly\" flag." << endl;
            std::terminate();
        }
        if(!genOnly && ntuple->METLV.size()>1) {
           cout << "ERROR::ConvertZllJetsToWlnuJets::main There is more than one entry in the MET vector." << endl;
           std::terminate();
        }
        if(!genOnly && ntuple->lLV.size()!=2) {
           cout << "ERROR::ConvertZllJetsToWlnuJets::main There are not two entries to the lepton vector." << endl
                << "\tSince the idea is to convert from Z->ll to W->lnu we need two leptons for this process to work." << endl;
           std::terminate();
        }
        if(!genOnly && ntuple->jLV.size()==0) {
            cout << "WARNING::ConvertZllJetsToWlnuJets::main There are no jets in this event." << endl
                 << "\tThis might signal that there is something wrong with this ntuple." << endl;
            std::terminate();
        }
        if(genOnly && ntuple->genParticleCollection.size()==0) {
            cout << "ERROR::ConvertZllJetsToWlnuJets::main There is no genParticleCollection." << endl;
            std::terminate();
        }

        //
        // Fill the Mll histogram which will be used to match the M-10To50 and M-50 samples together
        //
        if(ntuple->lLV.size()>1) {
            TLorentzVector ll = ntuple->lLV[0]+ntuple->lLV[1];
            Mll->Fill(ll.M());
        }

        //
        // Do the modifications to the ntuple
        //

        //
        // Scale the lepton and jet Lorentz vectors
        //
        if(!rebalanceEvt && !genOnly) {
            scaleFourVectors(ntuple, debug);
        }
        scaleGenVFourVectorAndDaugters(ntuple, debug);

        //
        // Move one of the leptons to the MET Lorentz vector
        //
        bool lowerPtLepton = moveLeptonToMET(ntuple, make_pair(probIndex1[0],probIndex1[1]), leptonSelectionCounter, METResMethod, 
                                             sumETMult, METResolution, METPhiResolution, resetMETMass,
                                             genOnly, debug);
        setGenLeptonToNeutrino(ntuple, lowerPtLepton, debug);

        //
        // Zero out the MET z component (really coming from the lepton)
        // This is just a double check as the same thing should have occured in moveLeptonToMET()
        //
        if(ntuple->METLV.size()>0) 
            ntuple->METLV[0].SetPz(0);

        //
        // Rebalance the event
        //
        if(rebalanceEvt)
            rebalanceEvent(ntuple, debug);

        //
        // Fill the tree
        //
        otree->Fill();
    }
    cout << endl;

    //
    // Scale the final histograms to the correct lumi*cross section*branching ratio*scale factor/num initial events
    //
    Mll->Scale(iProcess->getScaleFactor(DEFS::both));

    //
    // Write the modified tree and any histograms to the output file
    //
    otree->Write();
    Mll->Write();

    //
    // Calculate the output file's size and clean up
    //
    double fileSize = ofile->GetSize()/(1024.0*1024.0*1024.0);
    ofile->Close();

    //
    // Copy output file to EOS and delete the original if the copy is successful
    //
    if(doMove && oFilePath.find("root://")==string::npos)
        moveOutputFile(oFilePath, oProcessName, suffix);

    //
    // Print the program benchmarks and report the persentage of time we kept the higher/lower pT lepton
    //
    m_benchmark->Stop("sample");
    cout << "ConvertZllJetsToWlnuJetsNtuple_x" << endl
         << "\tMuon" << endl
         << "\t\tHigher pT lepton kept " << (double)leptonSelectionCounter.first[1]/(leptonSelectionCounterSum(leptonSelectionCounter,DEFS::muon))*100 << "% of the time" << endl
         << "\t\tLower pT lepton kept " << (double)leptonSelectionCounter.first[0]/(leptonSelectionCounterSum(leptonSelectionCounter,DEFS::muon))*100 << "% of the time" << endl
         << "\tElectron" << endl
         << "\t\tHigher pT lepton kept " << (double)leptonSelectionCounter.second[1]/(leptonSelectionCounterSum(leptonSelectionCounter,DEFS::electron))*100 << "% of the time" << endl
         << "\t\tLower pT lepton kept " << (double)leptonSelectionCounter.second[0]/(leptonSelectionCounterSum(leptonSelectionCounter,DEFS::electron))*100 << "% of the time" << endl
         << "\tFile size = " << fileSize << " GB" << endl
         << "\tCPU time = " << m_benchmark->GetCpuTime("sample") << " s" << endl
         << "\tReal time = " << m_benchmark->GetRealTime("sample") << " s" << endl;
    delete m_benchmark;

    return 0;
}

//______________________________________________________________________________
void saveOriginalInformation(EventNtuple *ntuple, int randIndex) {
    ntuple->METLV[0].rawLV = ntuple->METLV[0];
    ntuple->METLV[0].refLV = ntuple->lLV[randIndex];
}

//______________________________________________________________________________
pair<double,double> getMETResolution(EventNtuple *ntuple, int randIndex, double sumETMult) {
    //
    // The function for the MET resolution versus sumEt
    //
    static TF1 sigma_MET_func("sigma_MET","[0]+([1]*TMath::Sqrt(x))",0,2500);

    //
    // Calculate the sumEt
    // This is an approximation using (the leptons and) the jets as we don't have all of the PF particles
    //  available in our ntuple
    // We don't usually (ever) include the leptons from the Z as they should not be included per
    //  CHEF2013_JMEPerformance_Chayanit_Final.pdf
    //
    double sumEt = 0;
    //for(unsigned int ilepton=0; ilepton<ntuple->lLV.size(); ilepton++){
    //    sumEt += ntuple->lLV[ilepton].Et();
    //}
    for(unsigned int ijet=0; ijet<ntuple->jLV.size(); ijet++){
        sumEt += ntuple->jLV[ijet].Et();
    }
    sumEt*=sumETMult;

    //
    // Get sigma for MET in x
    //
    double sigma_MET_x = 0;
    if(ntuple->lLV[randIndex].leptonCat == DEFS::electron) {
        // If event comes from Z-->e+e-
        sigma_MET_func.SetParameters(0.05,0.63);
    }
    else if(ntuple->lLV[randIndex].leptonCat == DEFS::muon) {
        // If event comes from Z-->mu+mu-
        sigma_MET_func.SetParameters(0.87,0.62);
    }
    sigma_MET_x = sigma_MET_func.Eval(sumEt);

    //
    // Get sigma for MET in y
    //
    double sigma_MET_y = 0;
    if(ntuple->lLV[randIndex].leptonCat == DEFS::electron) {
        // If event comes from Z-->e+e-
        sigma_MET_func.SetParameters(0.90,0.59);
    }
    else if(ntuple->lLV[randIndex].leptonCat == DEFS::muon) {
        // If event comes from Z-->mu+mu-
        sigma_MET_func.SetParameters(1.42,0.60);
    }
    sigma_MET_y = sigma_MET_func.Eval(sumEt);

    return make_pair(sigma_MET_x,sigma_MET_y);
}

//______________________________________________________________________________
bool moveLeptonToMET(EventNtuple *ntuple, pair<double,double> probIndex1, pair_map &leptonSelectionCounter,
                     string METResMethod, double sumETMult, double METResolution, double METPhiResolution,
                     bool resetMETMass, bool genOnly, bool debug) {
    // Randomly choose 0 or 1 for the index of the first or second lepton
    // Generate a random number between 0 and 1
    // http://www.cplusplus.com/reference/random/mt19937/
    // http://www.cplusplus.com/reference/random/random_device/
    // Uncomment the next two lines if we want a random seed
    //std::random_device rd;
    //std::mt19937 mersenne(rd());
    // initialize our mersenne twister with a deterministic seed
    static unsigned seed = 1234567890;
    static std::mt19937 mersenne(seed);
    // Notice though that this modulo operation does not generate uniformly distributed random
    // numbers in the span (since in most cases this operation makes lower numbers slightly more likely).
    // Having %2 is okay since the choices are 0 and 1 (odd and even), but other numbers for the modulus
    // can cause problems
    //int randIndex = mersenne()%2;
    //bernoulli_distribution::param_type pp(0.9)
    //distribution.param(pp)
    // Using a single seed/generator and two different distributions for muons and electrons
    // It would be okay to use a single distribution except that the probIndex1 is not the same and using a
    //  single distribution would necessitate us resetting the param_type every time we go through the loop.
    // http://stackoverflow.com/questions/21327249/how-to-generate-uncorrelated-random-sequences-using-c
    static std::bernoulli_distribution distributionMuon    (probIndex1.first);
    static std::bernoulli_distribution distributionElectron(probIndex1.second);
    int randIndex = (ntuple->lLV[0].leptonCat==1) ? (int)distributionMuon(mersenne) : (int)distributionElectron(mersenne);
    (ntuple->lLV[0].leptonCat==1) ? leptonSelectionCounter.first[randIndex]++ : leptonSelectionCounter.second[randIndex]++;
    if(debug) {
        cout << "leptonCat: " << ntuple->lLV[0].leptonCat << endl;
        cout << "probIndex1: " << ((ntuple->lLV[0].leptonCat==1) ? probIndex1.first : probIndex1.second) << endl;
        cout << "seed: " << seed << " randIndex: " << randIndex << endl;
    }

    //
    // If we are running over only generator level information than we should return the randIndex and
    //  go no further.
    //
    if(genOnly) return randIndex;

    //
    // Copy the selected lepton to the refLV of the MET and the original METLV to the rawLV position 
    //
    saveOriginalInformation(ntuple, randIndex);

    // If METResMethod!="" then smear the selected lepton before it is added to the MET
    // We need to smear the MET because the lepton has a better resolution than the MET
    // Don't want to do it after because then we'd be smearing any MET that was already in the event
    // Based on https://indico.in2p3.fr/event/7691/session/2/contribution/25/material/0/0.pdf
    // sigma(MET_x,MET_y)=sigma_0+sigma_s*sqrt(sumEt)
    static std::mt19937 generator(seed-7890);
    if(METResMethod!="") {
        pair<double,double> sigma_MET;
        if(METResolution<=0.0) {
            sigma_MET = getMETResolution(ntuple, randIndex, sumETMult);
        }
        else {
            sigma_MET = make_pair(METResolution,METResolution);
        }

        double new_px=0.0, new_py=0.0, new_phi=0.0;
        if(METResMethod.find("coupled")!=string::npos) {
            double sigma_MET_xy = (sigma_MET.first+sigma_MET.second)/2.0;
            std::normal_distribution<double> n_distribution_px (ntuple->lLV[randIndex].Px(),sigma_MET_xy);
            new_px = n_distribution_px(generator);
            new_py = sqrt(abs(pow(ntuple->lLV[randIndex].Pt(),2)-pow(new_px,2)));
            //new_py*=(isless(pow(ntuple->lLV[randIndex].Pt(),2)-pow(new_px,2),0.0)) ? -1.0 : 1.0;
            new_py*=(isless(ntuple->lLV[randIndex].Py(),0.0)) ? -1 : 1;
            if(TMath::IsNaN(new_py)) {
                cout << "ERROR::moveLeptonToMET The new METLV.Py() value is a NaN" << endl;
                std::terminate();
            }
        }
        else if(METResMethod.find("simple")!=string::npos) {
            std::normal_distribution<double> n_distribution (1.0,sigma_MET.first);
            new_px = ntuple->lLV[randIndex].Px()*n_distribution(generator);
            new_py = ntuple->lLV[randIndex].Py()*n_distribution(generator);
        }
        else if(METResMethod.find("real")!=string::npos){
            std::normal_distribution<double> n_distribution_px (ntuple->lLV[randIndex].Px(),sigma_MET.first);
            std::normal_distribution<double> n_distribution_py (ntuple->lLV[randIndex].Py(),sigma_MET.second);
            new_px = n_distribution_px(generator);
            new_py = n_distribution_py(generator);
        }
        else if(METResMethod.find("simulated")!=string::npos) {
            
        }
        else {
            cout << "ERROR::moveLeptonToMET You must choose a known method for METResMethod. The options are real, simple, and coupled." << endl;
            std::terminate();
        }
        
        ntuple->lLV[randIndex].SetPxPyPzE(new_px, new_py, ntuple->lLV[randIndex].Pz(),ntuple->lLV[randIndex].E());
        //ntuple->lLV[randIndex].SetPxPyPzE(new_px, new_py, ntuple->lLV[randIndex].Pz(),
        //                                  sqrt(new_M2+pow(new_px,2)+pow(new_py,2)+pow(ntuple->lLV[randIndex].Pz(),2)));
        //ntuple->lLV[randIndex].SetXYZM(new_px, new_py, ntuple->lLV[randIndex].Pz(),new_M);
        //ntuple->lLV[randIndex].SetXYZM(new_px, new_py, 0.0, new_M);

        if(METResMethod.find("phi")!=string::npos) {
            std::normal_distribution<double> n_distribution_phi (1.0,METPhiResolution);
            new_phi = ntuple->lLV[randIndex].Phi()*n_distribution_phi(generator);
            ntuple->lLV[randIndex].SetPtEtaPhiE(ntuple->lLV[randIndex].Pt(), ntuple->lLV[randIndex].Eta(), new_phi, ntuple->lLV[randIndex].E());
        }
    }

    //Set Pz of the chosen lepton to be zero
    ntuple->lLV[randIndex].SetPz(0);

    // Add the lepton Lorentz vector to the MET Lorentz vector
    if(debug) {
        cout << "MET (before): ";
        ntuple->METLV[0].Print();
    }
    ntuple->METLV[0] += ntuple->lLV[randIndex];
    if(debug) {
        cout << "MET (after) ";
        ntuple->METLV[0].Print();
        cout << "Adding: ";
        ntuple->lLV[randIndex].Print();
    }

    // If requested, reset the MET mass distribution to be more like that of W+Jets
    if(resetMETMass) {
        double new_M=0.0, new_M2=0.0;
        static TF1 f("f","[3]*ROOT::Math::exponential_cdf_c(x+[0],[1],[2])+(x>1000)*[7]*ROOT::Math::exponential_pdf(x+[4],[5],[6])",0,40000);
        if(leptonSelectionCounterSum(leptonSelectionCounter,DEFS::both)==1) {
            cout << "Setting the M2 PRNG parameters for the first time ... " << flush;
            f.SetParameters(-5.56021e+01, 4.13745e-03, 5.56021e+01, 6.61098e+03,-9.99958e+02, 1.27623e-03, 4.99997e+01, 2.90200e+04);
            f.SetNpx(500);
            cout << "DONE" << endl;
        }
        new_M2 = f.GetRandom();
        if(new_M2<0) {
            cout << "ERROR::moveLeptonToMET The initial M2 value was less than 0." << endl;
            std::terminate();
        }
        static std::bernoulli_distribution met_m_pm(0.44);
        new_M = met_m_pm(generator) ? sqrt(new_M2) : -1.0*sqrt(new_M2);
        ntuple->METLV[0].SetXYZM(ntuple->METLV[0].Px(), ntuple->METLV[0].Py(), ntuple->METLV[0].Pz(), new_M);
    }

    // Remove the selected lepton from lLV
    if(debug) {
        cout << "Erasing: ";
        (ntuple->lLV.begin()+randIndex)->Print();
    }
    ntuple->lLV.erase(ntuple->lLV.begin()+randIndex);
    if(debug) {
        cout << "Remaining: ";
        ntuple->lLV.begin()->Print();
    }

    // Check that the size of the lepton vector is 1 after the erasure
    if(ntuple->lLV.size()!=1) {
       cout << "ERROR::ConvertZllJetsToWlnuJets::main The number of entries in the lepton vector is not 1." << endl
            << "\tThere should only be 1 entry after erasing the lepton that we added to the MET." << endl;
       std::terminate();
    }

    return randIndex;
}

//______________________________________________________________________________
/*
Z -> l+ l-
W -> l+ l-
if(choose l+ -> neutrino)
    W- -> nubar  l-
else if(choose l- -> neutrino)
    W+ -> nu l+

Test cases:
    If l==11 and choose -11 to move to neutrino
    daughterToChange_initialPDGID == -11
    signOfWAndLeptonToKeep = 1; //isless(-11,0.0) ? 1 : -1;
    24*=1; //==24
    daughterToChange_initialPDGID += daughterToChange_sign; // -11+-1 == -12
    24 —> 11 , -12

    if l==11 and choose 11 to move to neutrino
    daughterToChange_initialPDGID == 11
    daughterToChange_sign = 1; // isless(11,0.0) ? -1 : 1;
    signOfWAndLeptoToKeep = -1; // isless(11,0.0) ? 1 : -1;
    24*=-1; // ==-24
    daughterToChange_initialPDGID += daughterToChange_sign; // 11+1 == 12
    24 —> -11, 12
*/
void setGenLeptonToNeutrino(EventNtuple *ntuple, bool lowerPtLepton, bool debug) {
    //
    // Find the W boson
    // Then find the pT and position of all daughters
    //
    int boson_position = 0;
    map<double, int> daughters; //key==pT, value==position in genParticleCollection
    unsigned int size = min((int)ntuple->genParticleCollection.size(),50);
    for (unsigned int i=0; i<size; i++) {
        if (abs(ntuple->genParticleCollection[i].pdgId)!=24)
            continue;
        else {
            boson_position = i;
            if(debug) cout << "W" << ntuple->getLVString(ntuple->genParticleCollection[i].p4) << "->";
            for(unsigned int j=0; j<ntuple->genParticleCollection[i].daughterPositions.size(); j++) {
                int currentDaughterPosition = ntuple->genParticleCollection[i].daughterPositions[j];
                if (currentDaughterPosition<=500 &&
                    ntuple->leptonNeutrinoOrQuark(ntuple->genParticleCollection[currentDaughterPosition].pdgId)==EventNtuple::LEPTON) {
                    daughters[ntuple->genParticleCollection[currentDaughterPosition].p4.Pt()] = currentDaughterPosition;
                    if (debug)
                        cout << ntuple->genParticleCollection[currentDaughterPosition].pdgId
                             << ntuple->getLVString(ntuple->genParticleCollection[currentDaughterPosition].p4) << ",";
                }
            }
        }
    }
    if(daughters.size()!=2) {
        cout << "ERROR::setGenLeptonToNeutrino::The number of daughters from the W boson is " << daughters.size() << " and not 2!" << endl; 
        std::terminate();
    }
    int daughterToChange_position = (lowerPtLepton) ? daughters.begin()->second : daughters.rbegin()->second;
    int daughterToChange_initialPDGID = ntuple->genParticleCollection[daughterToChange_position].pdgId;
    int daughterToChange_sign = (isless(daughterToChange_initialPDGID,0.0)) ? -1 : 1;
    int signOfWAndLeptonToKeep = (isless(daughterToChange_initialPDGID,0.0)) ? 1 : -1;
    if(debug) {
        cout << "daughterToChange_position: " << daughterToChange_position << endl;
        cout << "daughterToChange_initialPDGID: " << daughterToChange_initialPDGID << endl;
        cout << "daughterToChange_sign: " << daughterToChange_sign << endl;
        cout << "signOfWAndLeptonToKeep: " << signOfWAndLeptonToKeep << endl; 
    }

    //
    // Change the W boson sign to match that of the remaining lepton
    //
    if(debug) cout << "V PDGID before: " << ntuple->genParticleCollection[boson_position].pdgId << endl;
    ntuple->genParticleCollection[boson_position].pdgId*=signOfWAndLeptonToKeep;
    if(debug) cout << "V PDGID after: " << ntuple->genParticleCollection[boson_position].pdgId << endl;

    //
    // Change the pdgId of the lepton begin converted to a neutrino
    //
    if(debug) cout << "lepton PDGID before: " << ntuple->genParticleCollection[daughterToChange_position].pdgId << endl;
    ntuple->genParticleCollection[daughterToChange_position].pdgId += daughterToChange_sign;
    if(debug) cout << "lepton PDGID after: " << ntuple->genParticleCollection[daughterToChange_position].pdgId << endl;
}

//______________________________________________________________________________
void print(TKinFitter *fitter) {
    std::cout << "=============================================" << std ::endl;
    std::cout << "-> Number of measured Particles  : " << fitter->nbMeasParticles() << std::endl;
    std::cout << "-> Number of unmeasured particles: " << fitter->nbUnmeasParticles() << std::endl;
    std::cout << "-> Number of constraints         : " << fitter->nbConstraints() << std::endl;
    std::cout << "-> Number of degrees of freedom  : " << fitter->getNDF() << std::endl;
    std::cout << "-> Number of parameters A        : " << fitter->getNParA() << std::endl;
    std::cout << "-> Number of parameters B        : " << fitter->getNParB() << std::endl;
    std::cout << "-> Maximum number of iterations  : " << fitter->getMaxNumberIter() << std::endl;
    std::cout << "-> Maximum deltaS                : " << fitter->getMaxDeltaS() << std::endl;
    std::cout << "-> Maximum F                     : " << fitter->getMaxF() << std::endl;
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++" << std ::endl;
    std::cout << "-> Status                        : " << fitter->getStatus() << std::endl;
    std::cout << "-> Number of iterations          : " << fitter->getNbIter() << std::endl;
    std::cout << "-> S                             : " << fitter->getS() << std::endl;
    std::cout << "-> F                             : " << fitter->getF() << std::endl;
    std::cout << "=============================================" << std ::endl;
}

//______________________________________________________________________________
void rebalanceEvent(EventNtuple *ntuple, bool debug) {
    /// internally use simple boolean for this constraint to reduce the per-event computing time
    bool constrainSumPt_ = true;

    TKinFitter* fitter_ = new TKinFitter("TopKinFitter", "TopKinFitter");
    fitter_->setMaxNbIter(500);
    fitter_->setMaxDeltaS(5e-05);
    fitter_->setMaxF(0.0001);
    fitter_->setVerbosity(0);
    //fitter_->setMaxNbIter(30);
    //fitter_->setMaxDeltaS(5e-02);
    //fitter_->setMaxF(0.1);
    //fitter_->setVerbosity(3);

    TMatrixD empty3x3(3,3);
    TMatrixD empty4x4(4,4);
    empty3x3.Zero();
    empty3x3(0,0) = 100.0;
    empty3x3(1,1) = 0.002;
    empty3x3(2,2) = 0.0004;
    empty4x4.Zero();
    empty4x4(0,0) = 100.1;
    empty4x4(1,1) = 100.0;
    empty4x4(2,2) = 100.0;
    empty4x4(3,3) = 0.1;
    TAbsFitParticle* lepton_   = new TFitParticleEScaledMomDev("Lepton",   "Lepton",   0, &empty3x3);
    TAbsFitParticle* neutrino_ = new TFitParticleEScaledMomDev("Neutrino", "Neutrino", 0, &empty3x3);
    TAbsFitParticle* hadQ1_    = new TFitParticleEMomDev      ("Jet1",     "Jet1",     0, &empty4x4);
    TAbsFitParticle* hadQ2_    = new TFitParticleEMomDev      ("Jet2",     "Jet2",     0, &empty4x4);

    /// supported constraints
    std::map<Constraint, TFitConstraintM*> massConstr_;
    TFitConstraintEp* sumPxConstr_;
    TFitConstraintEp* sumPyConstr_;
    massConstr_[kWMass         ] = new TFitConstraintM("WMass",         "WMass",         0, 0, WMASS);
    massConstr_[kNeutrinoMass  ] = new TFitConstraintM("NeutrinoMass",  "NeutrinoMass",  0, 0,    0.);
    sumPxConstr_                 = new TFitConstraintEp("SumPx",        "SumPx", 0, TFitConstraintEp::pX, 0.);
    sumPyConstr_                 = new TFitConstraintEp("SumPy",        "SumPy", 0, TFitConstraintEp::pY, 0.);

    massConstr_[kWMass       ]->addParticles1(lepton_, neutrino_);
    massConstr_[kNeutrinoMass]->addParticle1 (neutrino_);
    sumPxConstr_->addParticles(lepton_, neutrino_, hadQ1_, hadQ2_);
    sumPyConstr_->addParticles(lepton_, neutrino_, hadQ1_, hadQ2_);

    // add measured particles
    fitter_->addMeasParticle(hadQ1_);
    fitter_->addMeasParticle(hadQ2_);
    fitter_->addMeasParticle(lepton_);
    fitter_->addMeasParticle(neutrino_);

    // add constraints
    fitter_->addConstraint(massConstr_[kWMass]);
    fitter_->addConstraint(massConstr_[kNeutrinoMass]);
    if(constrainSumPt_) {
      fitter_->addConstraint(sumPxConstr_);
      fitter_->addConstraint(sumPyConstr_);
    }

    /*
    /// object used to construct the covariance matrices for the individual particles
    CovarianceMatrix* covM_;

    // initialize helper class used to bring the resolutions into covariance matrices
    if(udscResolutions_->size() &&  bResolutions_->size() && lepResolutions_->size() && metResolutions_->size())
      covM_ = new CovarianceMatrix(*udscResolutions_, *bResolutions_, *lepResolutions_, *metResolutions_,
                   *jetEnergyResolutionScaleFactors_, *jetEnergyResolutionEtaBinning_);
    else
      covM_ = new CovarianceMatrix();
    */

    // initialize particles
    const TLorentzVector p4HadQ1( ntuple->jLV[0].Px(), ntuple->jLV[0].Py(), ntuple->jLV[0].Pz(), ntuple->jLV[0].E() );
    const TLorentzVector p4HadQ2( ntuple->jLV[1].Px(), ntuple->jLV[1].Py(), ntuple->jLV[1].Pz(), ntuple->jLV[1].E() );
    const TLorentzVector p4Lepton  ( ntuple->lLV[0].Px(), ntuple->lLV[0].Py(), ntuple->lLV[0].Pz(), ntuple->lLV[0].E() );
    const TLorentzVector p4Neutrino( ntuple->METLV[0].Px(), ntuple->METLV[0].Py(), 0, ntuple->METLV[0].Et() );

    // initialize covariance matrices
    //TMatrixD covHadQ1 = covM_->setupMatrix(hadP, jetParam_);
    //TMatrixD covHadQ2 = covM_->setupMatrix(hadQ, jetParam_);
    //TMatrixD covLepton   = covM_->setupMatrix(lepton  , lepParam_);
    //TMatrixD covNeutrino = covM_->setupMatrix(neutrino, metParam_);

    // set the kinematics of the objects to be fitted
    hadQ1_->setIni4Vec( &p4HadQ1 );
    hadQ2_->setIni4Vec( &p4HadQ2 );
    lepton_->setIni4Vec( &p4Lepton );
    neutrino_->setIni4Vec( &p4Neutrino );

    /*
    hadQ1_->setCovMatrix( &covHadP );
    hadQ2_->setCovMatrix( &covHadQ );
    lepton_  ->setCovMatrix( &covLepton   );
    neutrino_->setCovMatrix( &covNeutrino );
    */

    if(constrainSumPt_){
      // setup Px and Py constraint for curent event configuration so that sum Pt will be conserved
      sumPxConstr_->setConstraint( p4HadQ1.Px() + p4HadQ2.Px() + p4Lepton.Px() + p4Neutrino.Px() );
      sumPyConstr_->setConstraint( p4HadQ1.Py() + p4HadQ2.Py() + p4Lepton.Py() + p4Neutrino.Py() );
    }

    // now do the fit
    fitter_->fit();

    // read back the resulting particles if the fit converged
    if(fitter_->getStatus()==0){
        // read back jet kinematics
        ntuple->jLV[0] = *(hadQ1_->getCurr4Vec());
        ntuple->jLV[1] = *(hadQ2_->getCurr4Vec());

        // read back lepton kinematics
        ntuple->lLV[0] = *(lepton_->getCurr4Vec());

        // read back the MET kinematics
        ntuple->METLV[0] = *(neutrino_->getCurr4Vec());
    }
    if(debug) {
        cout << "Fitter Status: " << fitter_->getStatus() << endl;
        print(fitter_);
    }

}

//______________________________________________________________________________
void scaleFourVectors(EventNtuple *ntuple, bool debug) {
    if(ntuple->lLV.size()>0) {
        for(unsigned int ilepton=0; ilepton<ntuple->lLV.size(); ilepton++) {
            ntuple->lLV[ilepton]*=(WMASS/ZMASS);
        }
    }
    if(ntuple->jLV.size()>0) {
        for(unsigned int ijet=0; ijet<ntuple->jLV.size(); ijet++) {
            ntuple->jLV[ijet]*=(WMASS/ZMASS);
        }
    }
}

//______________________________________________________________________________
void scaleGenVFourVectorAndDaugters(EventNtuple *ntuple, bool debug) {
    //
    // Change the boson type and scale it's 4-vector
    // Then change it's daughters' 4-vectors
    //
    unsigned int size = min((int)ntuple->genParticleCollection.size(),50);
    for (unsigned int i=0; i<size; i++) {
        if (abs(ntuple->genParticleCollection[i].pdgId)!=23)
            continue;
        else {            
            ntuple->genParticleCollection[i].pdgId = 24;
            ntuple->genParticleCollection[i].p4*=(WMASS/ZMASS);
            for(unsigned int j=0; i<ntuple->genParticleCollection[i].daughterPositions.size(); j++) {
                ntuple->genParticleCollection[ntuple->genParticleCollection[i].daughterPositions[j]].p4*=(WMASS/ZMASS);
            }
        }
    }
}

//______________________________________________________________________________
void moveOutputFile(string oFilePath, string oProcessName, string suffix) {
    cout << "Checking if system processor is available ... ";
    if (system(NULL)) cout << "Ok" << endl;
    else exit (EXIT_FAILURE);

    cout << "Check for the evironment variable MEInput ... ";
    char* MEInput;
    MEInput = getenv ("SMEInput");
    if (MEInput!=NULL) cout << "Ok" << endl;
    else exit (EXIT_FAILURE);

    string command = "xrdcp " + oFilePath + " " + MEInput + "WlnuJetsTest/" + oProcessName + suffix + ".root";
    cout << "Will use the following command to copy the output file to EOS:" << endl
         << "\t" << command << endl;
    int copyStatus = system(command.c_str());
    if(copyStatus == 0) {
        cout << "File " << oFilePath << " successfully copied to " << MEInput << "WlnuJetsTest/" << oProcessName << "!" << endl;
        string rmCommand = "rm " + oFilePath;
        cout << "Will use the following command to rm the original output file (keeping the newly copied one):" << endl
             << "\t" << rmCommand << endl;
        int rmStatus = system(rmCommand.c_str());
        if(rmStatus == 0) {
            cout << "File " << oFilePath << " successfully deleted!" << endl;
            return;
        }
        else {
            cout << "ERROR::moveOutputFile Delete of file " << oFilePath << " was unsuccesful!" << endl
                 << "\tThe file must be manually deleted!" << endl;
            exit (EXIT_FAILURE);
        }
    }
    else {
        cout << "ERROR::moveOutputFile Copy unsuccesful and returned the error code " << copyStatus << endl;
        exit (EXIT_FAILURE);
    }
}

//______________________________________________________________________________
int leptonSelectionCounterSum(pair_map &leptonSelectionCounter, DEFS::LeptonCat leptonCat) {
    if(leptonCat == DEFS::muon)
        return leptonSelectionCounter.first[0]+leptonSelectionCounter.first[1];
    else if(leptonCat == DEFS::electron)
        return leptonSelectionCounter.second[0]+leptonSelectionCounter.second[1];
    else if(leptonCat == DEFS::both)
        return leptonSelectionCounter.first[0]+leptonSelectionCounter.first[1]+leptonSelectionCounter.second[0]+leptonSelectionCounter.second[1];
    else {
        cout << "WARNING::ConvertZllJetsToWlnuJetsNtuple_x::leptonSelectionCounterSum Unable to determine the leptonCat and leptonSelectionCounterSum" << endl;
        return -1;
    }
}
