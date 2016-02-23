#include "TROOT.h"
#include "TSystem.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"

#include <iostream>
#include <string>

using namespace std;

void renameTTree(string filename, string directory,
                 string treename_old, string treename_new) {
   TFile* f = TFile::Open(filename.c_str(),"UPDATE");
   gDirectory->ls();
   f->cd(directory.c_str());
   gDirectory->ls();
   TTree *t;
   f->GetObject((directory+"/"+treename_old).c_str(),t);
   cout << "old tree memory address: " << t << endl;
   cout << "old tree entries: " << t->GetEntries() << endl;
   if(t){
      t->SetName(treename_new.c_str());
      t->Write();
   }
   gDirectory->ls();
   gDirectory->Delete((treename_old+";*").c_str());
   gDirectory->ls();
   f->Close();
   f = TFile::Open(filename.c_str(),"READ");
   gDirectory->ls();
   gDirectory->cd(directory.c_str());
   gDirectory->ls();
   f->GetObject((directory+"/"+treename_new).c_str(),t);
   cout << "new tree memory address: " << t << endl;
   cout << "new tree entries: " << t->GetEntries() << endl;
}
