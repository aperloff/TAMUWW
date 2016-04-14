#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include <vector>
#include <map>

#include "TAMUWW/MEPATNtuple/interface/METree.hh"

using namespace std;

void mvatest() {
   TFile* f = TFile::Open("/uscms_data/d2/aperloff/Spring12ME7TeV/MEResults/microNtuples_oldStructure/microWW_EPDv01.root");
   TTree* t = (TTree*)f->Get("METree");
   METree* metree = new METree();
   t->SetBranchAddress("METree",&metree);
   vector<TString> v;
   v.push_back("BDT");
   t->GetEntry(10);
   metree->setMVAReader(v,"/uscms_data/d2/aperloff/2012_10_12_TMVA_plots_cuts_muons_reducedVariables/weights/");
   t->GetEntry(10);
   cout << "Entry 10" << endl;
   cout << "\tvector return (resp): " << metree->getMVAOutput(v,"response")[0] << endl;
   cout << "\tdouble return (resp): " << metree->getMVAOutput("BDT","response") << endl;
   cout << "\tvector return (prob): " << metree->getMVAOutput(v,"probability")[0] << endl;
   cout << "\tdouble return (prob): " << metree->getMVAOutput("BDT","probability") << endl;
   t->GetEntry(20);
   cout << "Entry 20" << endl;
   cout << "\tvector return (resp): " << metree->getMVAOutput(v,"response")[0] << endl;
   cout << "\tdouble return (resp): " << metree->getMVAOutput("BDT","response") << endl;
   cout << "\tvector return (prob): " << metree->getMVAOutput(v,"probability")[0] << endl;
   cout << "\tdouble return (prob): " << metree->getMVAOutput("BDT","probability") << endl;
   t->GetEntry(30);
   cout << "Entry 30" << endl;
   cout << "\tvector return (resp): " << metree->getMVAOutput(v,"response")[0] << endl;
   cout << "\tdouble return (resp): " << metree->getMVAOutput("BDT","response") << endl;
   cout << "\tvector return (prob): " << metree->getMVAOutput(v,"probability")[0] << endl;
   cout << "\tdouble return (prob): " << metree->getMVAOutput("BDT","probability") << endl;
}
