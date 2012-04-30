#include "TAMUWW/SpecialTools/interface/DefaultValues.hh"

#include <iostream>
#include "TAMUWW/SpecialTools/interface/Table.hh"
#include "TAMUWW/SpecialTools/interface/FileLocationTable.hh"
#include "TAMUWW/SpecialTools/interface/TableCellVal.hh"
#include "TAMUWW/SpecialTools/interface/TableCellText.hh"
#include "TAMUWW/SpecialTools/interface/PhysicsProcess.hh"

using std::string;
using std::vector;
using std::map;
using std::cout;
using std::endl;



// ----------------------------------------------------------------------------
// This method returns the table with the event expectation for the evt/tag category
Table DefaultValues::getNormTable(DEFS::LeptonCat evtcat, DEFS::TagCat tagcat){

  Table table("NormTable");
  
  string eventEstimatesFile = "TAMUWW/ConfigFiles/Official/EventEstimates_";

  // add the tag name and the ".txt" at the end
  eventEstimatesFile += DEFS::getEventCatString(evtcat)+"_";
  eventEstimatesFile += DEFS::getTagCatString(tagcat);
  eventEstimatesFile += ".txt";

  if(!table.parseFromFile(eventEstimatesFile))
    cout<<"ERROR  DefaultValues::getNormTable() cannot parse config file "
	<<eventEstimatesFile<<endl;

  return table;
  
}//getNormTable

// ----------------------------------------------------------------------------
// This method returns the appropriate table of file location for each jet bin
// and tag category
Table DefaultValues::getFileLocationTable(DEFS::TagCat tagcat){ 

  // The location of the table with the file locations
  string fileLocationFile = "TAMUWW/ConfigFiles/Official/FileLocation_";
  
  // add the tag name and the ".txt" at the end
  fileLocationFile += DEFS::getTagCatString(tagcat);
  fileLocationFile += ".txt";

  // Create the table and parse the contents of the file
  FileLocationTable table("FileLocationTable");  
  if(!table.parseFromFile(fileLocationFile))
    cout<<"ERROR  DefaultValues::getFileLocationTable() cannot parse config file "
	<<fileLocationFile<<endl;


  // make sure you add the basepath to the table
  table.addBasePath();

  return table;
  
}//getFileLocationTable



// ----------------------------------------------------------------------------
vector < PhysicsProcess * > DefaultValues::getProcesses(vector<DEFS::PhysicsProcessType> processName,
							DEFS::JetBin jetBin, 
							DEFS::TagCat tagcat){

  // The returning vector of processes
  vector<PhysicsProcess*>  proc;

  // get the table with the expected number of 
  map<DEFS::LeptonCat, Table> normTable;
  normTable[DEFS::muon    ] = getNormTable(DEFS::muon    ,tagcat);
  normTable[DEFS::electron] = getNormTable(DEFS::electron,tagcat);

  // get the table with the files location
  Table fileTable = getFileLocationTable(tagcat);

  // Loop over all the process names
  for (unsigned int prn = 0; prn < processName.size(); prn++){

    PhysicsProcess * pr = getSingleProcess(processName[prn], jetBin, normTable, fileTable);
    if (pr == 0) {
      cout<<"ERROR DefaultValues::getProcesses couldnot add process"<<endl;
      continue;
    }

    proc.push_back(pr);
    
  }//for

  return proc;

}// getProcesses


// ----------------------------------------------------------------------------
// Return a PhysicsProcess object for process "process" from the given inputs.
// For some reson it does not compile when I use the signature
// (..., const Table & normTable, const Table & fileTable, ...) 
PhysicsProcess * DefaultValues::getSingleProcess(DEFS::PhysicsProcessType process,
						 DEFS::JetBin jetBin,
						 map<DEFS::LeptonCat, Table> normTable,
						 Table fileTable){

    // get the process name
  string prName = DEFS::PhysicsProcess::getTypeString(process);
  
  // get the name of the jetBin
  string jetBinName = DEFS::getJetBinString(jetBin);
    
  // find the file location for that process
  TableCellText * cellFile = (TableCellText *) fileTable.getCellRowColumn(prName,"FilePath");
  
  // make sure we found the cell
  if (cellFile == 0){
    cout<<"ERROR DefaultValues::getSingleProcess Table "<<fileTable.getTableOrigin()
	<<" does not have row "<<prName
	<<" and column FilePath"<<endl;
    cout<<" SKIPPING PROCESS "<<prName<<endl;
    return 0;
  }
  
  // Create the chain to that file    
  TChain * chain ;
  chain  = new TChain("METree");
  chain->Add(cellFile->text.c_str());
  
  // Create the PhysicsProcess 
  PhysicsProcess *  proc =  new PhysicsProcess(prName, prName, chain);

  // Tell it the formula to get the categories from its own data
  proc->setCategory("EvtTree.leptonCat");//"h.det");
  
  // Set the expected number of events for each category
  // iterating over the map of normTables.
  for ( map<DEFS::LeptonCat, Table>::iterator it = normTable.begin();
	it != normTable.end(); it++){

    // Get the cell from the table
    TableCellVal * cellNorm = (TableCellVal *) it->second.getCellRowColumn(prName, jetBinName);
    
    // make sure we found the cell
    if (cellNorm == 0){
      cout<<"ERROR DefaultValues::getSingleProcess normalization Table " <<it->second.getTableOrigin()
	  <<" for LeptonCat="<<DEFS::getEventCatString(it->first)
	  <<" does not have row "<<prName
	  <<" and column "<<jetBinName<<endl;
      cout<<" SKIPPING PROCESS "<<prName<<endl;
      return 0;
    }
 
    proc->setCategoryNorm(it->first, cellNorm->val);

  }//for map of categories

  // and return it.
  return proc;

}//getSingleProcess


// ----------------------------------------------------------------------------
string  DefaultValues::getWeightForCategory(DEFS::TagCat tagcat, DEFS::PhysicsProcessType type, int iDet){

  // The returning weights, for the moment just "1"
  return "1";

  /*
    string wei ;
  // Set the detector type
  if (iDet ==  TopLepType::TCEM)
    wei += "(h.det==1)";
  else if (iDet ==  TopLepType::TPHX)
    wei += "(h.det==2)";
  else if (iDet ==  TopLepType::TCMUP)
    wei += "(h.det==3)";
  else if (iDet ==  TopLepType::TCMX)
    wei += "(h.det==4)";
  else if (iDet == 5)
    wei += "(h.det>=5)";
  else 
    wei += "(h.det>0)";

  // Set the QCD veto. For PHX add a h.met>25 cut
  wei += "*((h.det==2)*(h.met>25)+(h.det != 2))";


  // For data require 
  if (type == DEFS::Data ){
    if (tagcat == DEFS::pretag )    
      wei += "*1";
    else if (tagcat == DEFS::eq0TSV)
      wei += "*(h.ntag==0)";
    else if (tagcat == DEFS::eq1TSV)
      wei += "*(h.ntag==1)";
    else if (tagcat == DEFS::eq2TSV)
      wei += "*(h.ntag==2)";      
    else
      cout<<"ERROR DefaultValues::getWeightForCategory (data) called with tagcat="<<tagcat<<endl;     
  }
  // For MC require 
  else {

    // I have to apply this for the 1 and the 2 tags. Do I need it for the 0 tag ??
    if (type != DEFS::WLight && type != DEFS::NonW)
      if (tagcat != DEFS::pretag &&
	  tagcat != DEFS::eq0TSV )
	wei += "*(h.ntag>0)";

    // Pretag
    //if (tag==-1) 
    //wei += "*h.wgt*(h.tagProb0*h.passQCD + h.tagProb1*h.passQCD + h.tagProb2*( (h.det==2)*h.passQCD+(h.det!=2) ) )";

    //switch on ntags
    if (tagcat == DEFS::pretag )      // pretag sample 
      wei += "*(h.wgt*h.passQCD)";
    else if (tagcat == DEFS::eq0TSV )      // untag sample 
      wei += "*(h.tagProb0*h.wgt*h.passQCD)";
    else if (tagcat == DEFS::eq1TSV) // single tag sample
      wei += "*(h.tagProb1*h.wgt*h.passQCD)";
    else if (tagcat == DEFS::eq2TSV) // double tags sample
      wei += "*(h.tagProb2*h.wgt*( (h.det==2)*h.passQCD+(h.det!=2) ))";      
    else
      cout<<"ERROR DefaultValues::getWeightForCategory (MC) called with tagcat="<<tagcat<<endl;     

  }

  // For MC only, use triggerTO for tight, but getMet2JetTurnOn_v2 for loose leptons
  if (type != DEFS::Data)
    wei += "*(((h.det<5) * h.triggerTO)+ (h.det>=5) * 0.982 /(1.0 + exp(-(h.corrVtxMet-42.78)/3.681)))";


  return wei;
  */

}//getWeightForCategory

// ----------------------------------------------------------------------------
vector < PhysicsProcess * > DefaultValues::getProcessesWW(DEFS::JetBin jetBin,
							 DEFS::TagCat tagcat, 
							 bool include_data ){

  vector<DEFS::PhysicsProcess::Type> procs;


  procs.push_back(DEFS::PhysicsProcess::STopS   );
  procs.push_back(DEFS::PhysicsProcess::STopT   );
  procs.push_back(DEFS::PhysicsProcess::STopTW  );
  procs.push_back(DEFS::PhysicsProcess::TTbar   );
  //procs.push_back(DEFS::PhysicsProcess::TTbarLJ );
  //procs.push_back(DEFS::PhysicsProcess::TTbarDil); 
  //procs.push_back(DEFS::PhysicsProcess::Wbb     );
  //procs.push_back(DEFS::PhysicsProcess::Wcc     );
  //procs.push_back(DEFS::PhysicsProcess::WLight  );
  procs.push_back(DEFS::PhysicsProcess::WJets   ); 
  procs.push_back(DEFS::PhysicsProcess::Zjets   );
  procs.push_back(DEFS::PhysicsProcess::Ztautau );
  procs.push_back(DEFS::PhysicsProcess::QCDMu               );
  procs.push_back(DEFS::PhysicsProcess::QCDEl_Pt30to80      );
  procs.push_back(DEFS::PhysicsProcess::QCDEl_Pt80to170     );
  procs.push_back(DEFS::PhysicsProcess::QCDEl_BCtoE30to80   );
  procs.push_back(DEFS::PhysicsProcess::QCDEl_BCtoE80to170  );
  //procs.push_back(DEFS::PhysicsProcess::QCD250  );
  procs.push_back(DEFS::PhysicsProcess::WW      );
  procs.push_back(DEFS::PhysicsProcess::WZ      );
  //procs.push_back(DEFS::PhysicsProcess::ZZ      );

  if (include_data)
    procs.push_back(DEFS::PhysicsProcess::Data    );

  return getProcesses(procs, jetBin, tagcat);

}//getProcessesWW
