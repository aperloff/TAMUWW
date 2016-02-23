// Name: PhysicsProcess.cc
// Author: Travis Lamb

// This class will eventually replace an older Physics class. It is a base class that simply holds data members.

#include "TAMUWW/SpecialTools/interface/PhysicsProcess.hh"

// ##################################################
// ################ PHYSICS PROCESS #################
// ##################################################

PhysicsProcess::PhysicsProcess (string procName,
                                string groupingName,
                                string fileNameTEMP,
                                string treeName):
   name(procName),
   groupName(groupingName),
   fileName(fileNameTEMP),
   chain (0)
{
   TString currentDir = gDirectory->GetPathStatic();
   TFile * file = TFile::Open(fileName,"READ");//TFile::Open(fileName);
   if (!file)
   {
      cout << "ERROR PhysicsProcess::PhysicsProcess() could not open file " << fileName << endl;
      return;
   }
   gDirectory->cd(currentDir);
   /*
   if (!file->cd("PS"))
   {
      cout << "ERROR PhysicsProcess::PhysicsProcess() could not CD into directory PS in file " << fileName << endl;
      return;
   }
   */

   chain = (TChain*) file->Get(treeName.c_str());
   if (chain == 0)
   {
      cout << "ERROR PhysicsProcess::PhysicsProcess() could not find tree named " << treeName << " in file " << fileName << endl;
      return;
   }

   PhysParMap blankD;
   PhysParMapUI blankUI;
   for(unsigned int i = 0; i<DEFS::nLeptonCat; i++) {
      blankD[(DEFS::LeptonCat)i] = 0.0;
      blankUI[(DEFS::LeptonCat)i] = 1;
   }
   sigma = blankD;
   branching_ratio = blankD;
   intLum = blankD;
   initial_events = blankUI;
   scaleFactor = blankD;
}

void PhysicsProcess::setPhysicsParameters(PhysParMap cross_section, PhysParMap lum, PhysParMap br, PhysParMapUI in_ev, PhysParMap sf){
   for(PhysParMap::iterator it = cross_section.begin(); it!=cross_section.end(); it++) {
      sigma[it->first] = it->second;
   }
   for(PhysParMap::iterator it = br.begin(); it!=br.end(); it++) {
      branching_ratio[it->first] = it->second;
   }
   for(PhysParMap::iterator it = lum.begin(); it!=lum.end(); it++) {
      intLum[it->first] = it->second;
   }
   for(PhysParMapUI::iterator it = in_ev.begin(); it!=in_ev.end(); it++) {
      initial_events[it->first] = it->second;
   }
   for(PhysParMap::iterator it = sf.begin(); it!=sf.end(); it++) {
      scaleFactor[it->first] = it->second;
   }
}

void PhysicsProcess::fillMETreeIndexMap(){
   METree* metree = new METree();
   TBranch * branch = chain->GetBranch("METree");
   if (branch) {
      chain->SetBranchAddress("METree",&metree);
      chain->GetEntry(0);
      indexMapOfME = metree->fillIndexMap();
      chain->ResetBranchAddress(branch);
   }
   else {
      cout << "ERROR PhysicsProcess::fillMETreeIndexMap() could not find \"METree\" in the chain to fill the indexMap" << endl;
   }

  delete metree;
  return;
}

// ##################################################
// ############ COLORED PHYSICS PROCESS #############
// ##################################################

PlotterPhysicsProcess::PlotterPhysicsProcess (string procName,
                                              string groupingName,
                                              string fileNameTEMP,
                                              int col,
                                              string treeName):
   PhysicsProcess(procName, groupingName, fileNameTEMP, treeName),
   color(col){
   
}
