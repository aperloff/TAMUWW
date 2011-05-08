/////////////////////////////////////////////////////////////////////////////////////////////////
//// CMS 
//// WW CrossSection Measurement using Matrix Element 
//// Modified by Osipenkov, Ilya : ilyao@fnal.gov
/////////////////////////////////////////////////////////////////////////////////////////////////
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "TAMUWW/Integrator/interface/CubaIntegrator.hh"
#include "TAMUWW/Integrator/interface/RootIntegrator.hh"
#include "TAMUWW/MatrixElement/interface/EventFile.hh"
#include "TAMUWW/MatrixElement/interface/MEConstants.hh"
#include "TAMUWW/MatrixElement/interface/MEJob.hh"

#include "TAMUWW/MatrixElement/interface/sChannelEventProb2Jet.hh"
#include "TAMUWW/MatrixElement/interface/tChannelEventProb2Jet.hh"
#include "TAMUWW/MatrixElement/interface/ttEventProb2Jet.hh"
#include "TAMUWW/MatrixElement/interface/WbbEventProb2Jet.hh"
#include "TAMUWW/MatrixElement/interface/WcEventProb2Jet.hh"
#include "TAMUWW/MatrixElement/interface/WggEventProb2Jet.hh"
#include "TAMUWW/MatrixElement/interface/WHEventProb2Jet.hh"
#include "TAMUWW/MatrixElement/interface/WWEventProb2Jet.hh"
#include "TAMUWW/MatrixElement/interface/WZEventProb2Jet.hh"
#include "TAMUWW/MatrixElement/interface/HWWEventProb2Jet.hh"
#include "TAMUWW/MatrixElement/interface/WLightEventProb2Jet.hh"
#include "TAMUWW/MatrixElement/interface/ZLightEventProb2Jet.hh"
#include "TAMUWW/MatrixElement/interface/QCDEventProb2Jet.hh"

#include "TAMUWW/MatrixElement/interface/WHEventProb3Jet.hh"
#include "TAMUWW/MatrixElement/interface/sChannelEventProb3Jet.hh"
#include "TAMUWW/MatrixElement/interface/tChannelEventProb3Jet.hh"
#include "TAMUWW/MatrixElement/interface/ttEventProb3Jet.hh"
#include "TAMUWW/MatrixElement/interface/WbbEventProb3Jet.hh"

using std::atoi;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

int main(int argc, char* argv[]){
  if (argc != 7)
   {
     cerr << "Usage: ./run_MatrixElement <inputfile> <outputfile> <treename> [nEvents] [nSkip] [nPrescale]\n";
     return 1;
   }
  
  string inputFilename(argv[1]);
  string outputFilename(argv[2]);
  string treeName(argv[3]);
  int nEvents   = atoi(argv[4]);
  int nSkip     = atoi(argv[5]);
  int nPrescale = atoi(argv[6]);
  
  cout << " ++++ Reading from " << inputFilename << " ++++\n"
       << " ++++ Writing to " << outputFilename << " ++++\n"
       << " ++++ Reading tree " << treeName << " ++++\n"
       << " ++++ Running over " << nEvents << " events ++++\n"
       << " ++++ Skipping the first " << nSkip << " events ++++\n"
       << " ++++ Prescale by " << nPrescale << " ++++" << endl;

// ****************************
// **** Transfer Functions ****
// ****************************

  // NEW TF's
  edm::FileInPath tf_ttbar_uds_file("TAMUWW/ConfigFiles/Official/TransferFunctions/TF_TTbarMG_UDS.txt");
  edm::FileInPath tf_ttbar_g_file("TAMUWW/ConfigFiles/Official/TransferFunctions/TF_TTbarMG_G.txt");
  edm::FileInPath tf_ttbar_b_file("TAMUWW/ConfigFiles/Official/TransferFunctions/TF_TTbarMG_B.txt");
 
  // The double-gaussian transfer functions
  DGTransferFunction bTF    ( tf_ttbar_b_file.fullPath());
  DGTransferFunction lightTF( tf_ttbar_uds_file.fullPath());
  DGTransferFunction gluonTF( tf_ttbar_g_file.fullPath());
  

  
// ------------------------------------------------------------------
  
// *************************
// **** Input File Type ****
// *************************
  //RecoRootEventFile inputFile(inputFilename, "EvtTree", nEvents, nSkip, nPrescale, false);

  EventNtupleEventFile inputFile(inputFilename, treeName, nEvents, nSkip, nPrescale, false);

  // MadEventFile inputFile(inputFilename, nEvents, nSkip, nPrescale, false);
  // SmearedMadEventFile inputFile(inputFilename, nEvents, nSkip, nPrescale,false);
  //   SingleTopNtupleEventFile inputFile(inputFilename, "SingleTopNtuple",
  //                                      nEvents, nSkip, nPrescale, false);
  //  inputFile.setJetCorrLevel(5);
  
// ------------------------------------------------------------------

// ********************
// **** Integrator ****
// ********************
// The type of integration
//   NullIntegrator integrator;
  
  RootIntegrator rootInt;
  rootInt.setNcomp(1);
  rootInt.setEpsilon(.01, 0);
  rootInt.setVerbose(0);
  rootInt.setSampleSet(true);
  rootInt.setPseudoRandom(true);
  rootInt.setNeval(0, 1000000);
  
//   CubaInt::Cuhre integrator;
//   integrator.setCubatureRule(9);

//   CubaInt::Vegas integrator;
//   integrator.setParams(1000, 500);

//   CubaInt::Suave integrator;
//   integrator.setNew(1000);
//   integrator.setFlatness(25);

  CubaInt::Divonne divonneInt;
  divonneInt.setPartitioningRule(47);
  divonneInt.setIntegrationRule(1);
  divonneInt.setRefinementRule(1);
  divonneInt.setMaxPass(5);
  divonneInt.setBorder(0);
  divonneInt.setChi2(10);
  divonneInt.setMinDev(.25);
//   divonneInt.setMaxPeaks(4);
//   divonneInt.setPeakFinder(MEJob::sm_peakFinder);
  divonneInt.setNcomp(1);
  divonneInt.setEpsilon(.01, 0);
  divonneInt.setVerbose(0);
  divonneInt.setSampleSet(true);
  divonneInt.setPseudoRandom(true);
  divonneInt.setNeval(0, static_cast<int>(1e7));


// ------------------------------------------------------------------

  MEJob job(inputFile, outputFilename);

// ************************
// **** Dynamic Bounds ****
// ************************
  bool dynBounds;
  dynBounds = true;
  job.setDynamicBounds(dynBounds);
  
  cout << "Using b-jet " << bTF.getName() << endl;
  cout << "Using light-jet " << lightTF.getName() << endl;
  cout << "Using gluon-jet " << gluonTF.getName() << endl;
  if (dynBounds)
    cout << "Using dynamic bounds" << endl;
  
// ------------------------------------------------------------------

// *************************
// **** Matrix Elements ****
// *************************
  vector<EventProb*> eventProbs2jet;
  vector<EventProb*> eventProbs3jet;


  eventProbs2jet.push_back(new WWEventProb2Jet(rootInt, lightTF));
  eventProbs2jet.push_back(new WZEventProb2Jet(rootInt, lightTF));
  eventProbs2jet.push_back(new WLightEventProb2Jet(rootInt, lightTF));
  eventProbs2jet.push_back(new ZLightEventProb2Jet(rootInt, lightTF));
  eventProbs2jet.push_back(new ttEventProb2Jet(divonneInt, bTF));
  eventProbs2jet.back()->setBounds(3, 0, MEConstants::beamEnergy);
  eventProbs2jet.back()->setBounds(4, 0, TMath::TwoPi());
  eventProbs2jet.back()->setBounds(5, 0, TMath::Pi());
  eventProbs2jet.push_back(new tChannelEventProb2Jet(rootInt, bTF, lightTF));
  eventProbs2jet.push_back(new sChannelEventProb2Jet(rootInt, bTF));
  eventProbs2jet.push_back(new QCDEventProb2Jet(rootInt, gluonTF));


  ///HWW:
  eventProbs2jet.push_back(new HWWEventProb2Jet(rootInt, lightTF, 120));
  eventProbs2jet.push_back(new HWWEventProb2Jet(rootInt, lightTF, 130));
  eventProbs2jet.push_back(new HWWEventProb2Jet(rootInt, lightTF, 140));
  eventProbs2jet.push_back(new HWWEventProb2Jet(rootInt, lightTF, 150));
  eventProbs2jet.push_back(new HWWEventProb2Jet(rootInt, lightTF, 160));

  ///WH:
  eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF, 115));
  eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF, 120));
  eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF, 125));
  eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF, 130));
  eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF, 135));


  // // Simplify the computation for now
//   eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF, 100)); //0
// eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF, 105)); //1
// //   eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF_WH110GeV, 110)); //2
// //   eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF_WH115GeV, 115)); //3
// //   eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF_WH120GeV, 120)); //4
// //   eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF_WH125GeV, 125)); //5
// //   eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF_WH130GeV, 130)); //6
// //   eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF_WH135GeV, 135)); //7
// //   eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF_WH140GeV, 140)); //8
// //   eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF_WH145GeV, 145)); //9
// //   eventProbs2jet.push_back(new WHEventProb2Jet(rootInt, bTF_WH150GeV, 150)); //10

// //   //s-channel
// //   eventProbs2jet.push_back(new sChannelEventProb2Jet(rootInt, bTF_schan)); //11 
// //   eventProbs3jet.push_back(new sChannelEventProb3Jet(divonneInt, bTF, gluonTF)); //11

// //   //t-channel
// //   eventProbs2jet.push_back(new tChannelEventProb2Jet(rootInt, bTF_tchan, lightTF_tchan)); //12
// //   eventProbs3jet.push_back(new tChannelEventProb3Jet(divonneInt, bTF, lightTF)); //12

// // //   //t-channel alternate
// // //   eventProbs2jet.push_back(new tChannelEventProb2JetAlt(rootInt, bTF_tchan, lightTF_tchan)); //13
// // //   eventProbs3jet.push_back(new tChannelEventProb3JetAlt(divonneInt, bTF, lightTF)); //13

// //   //Wbbbar
// //   eventProbs2jet.push_back(new WbbEventProb2Jet(rootInt,  bTF_WBB)); //14
// //   eventProbs3jet.push_back(new WbbEventProb3Jet(divonneInt, bTF, gluonTF)); //14

// //   //Wcc ME, use same as Wbb ME but with a different TF
// //   eventProbs2jet.push_back(new WbbEventProb2Jet(rootInt,  cTF_WCC)); //15
// //   eventProbs3jet.push_back(new WbbEventProb3Jet(divonneInt, cTF_WCC, gluonTF)); //15
  
// //   //Wc, use a Dummy class that is never evaluated for the 3 jets events
// //   eventProbs2jet.push_back(new WcEventProb2Jet(divonneInt, cTF_WC, gluonTF_WC));//16
// //   eventProbs3jet.push_back(new DummyEventProb("Dummy WC",rootInt)); //16

// //   //Wjg, use a Dummy class that is never evaluated for the 3 jets events
// //   eventProbs2jet.push_back(new WjgEventProb2Jet(divonneInt, lightTF_Wjg,gluonTF_Wjg)); //17
// //   eventProbs3jet.push_back(new DummyEventProb("Dummy Wjg",rootInt)); //17

// //   //Wgg, use a Dummy class that is never evaluated for the 3 jets events
// //   eventProbs2jet.push_back(new WggEventProb2Jet(divonneInt, gluonTF_Wgg)); //18
// //   eventProbs3jet.push_back(new DummyEventProb("Dummy Wgg",rootInt)); //18

// //   // tt-bar matrix element 2-jet
// //   eventProbs2jet.push_back(new ttEventProb2Jet(divonneInt, bTF_ttbar)); //19
// //   eventProbs2jet.back()->setBounds(3, 0, MEConstants::beamEnergy);
// //   eventProbs2jet.back()->setBounds(4, 0, TMath::TwoPi());
// //   eventProbs2jet.back()->setBounds(5, 0, TMath::Pi());
  
// //   // tt-bar matrix element 3-jet
// //   eventProbs3jet.push_back(new ttEventProb3Jet(divonneInt, bTF, lightTF)); //19
// //   eventProbs3jet.back()->setBounds(3, 0, MEConstants::beamEnergy);
// //   eventProbs3jet.back()->setBounds(4, 0, TMath::TwoPi());
// //   eventProbs3jet.back()->setBounds(5, 0, TMath::Pi());
  
// //   // Diboson WW, use a Dummy class that is never evaluated for the 3 jets events
// //   eventProbs2jet.push_back(new WWEventProb2Jet(rootInt, lightTF_diboson)); //20
// //   eventProbs3jet.push_back(new DummyEventProb("Dummy WW", rootInt)); //20

// //   // Diboson WZ, use a Dummy class that is never evaluated for the 3 jets events
// //   eventProbs2jet.push_back(new WZEventProb2Jet(rootInt, lightTF_diboson));//21
// //   eventProbs3jet.push_back(new DummyEventProb("Dummy WZ", rootInt));//21
 
  
//   // ********************************************************
//   // The following are the ME's for Ricardo's correlated EPD
//   // ********************************************************
//   vector<EventProb*> evtPrb2jetCEPD;

//   if (1){

//     //WH matrix elements
//     evtPrb2jetCEPD.push_back(new WHEventProb2Jet(rootInt, bTF_WH100GeV, 100)); //0
     
//     //s-channel
//     evtPrb2jetCEPD.push_back(new sChannelEventProb2Jet(rootInt, bTF_schan)); //1 
    
//     //t-channel
//     evtPrb2jetCEPD.push_back(new tChannelEventProb2Jet(rootInt, bTF_tchan, lightTF_tchan)); //2
    
//     //t-channel alternate
//     evtPrb2jetCEPD.push_back(new tChannelEventProb2JetAlt(rootInt, bTF_tchan, lightTF_tchan)); //3
    
//     //Wbbbar
//     evtPrb2jetCEPD.push_back(new WbbEventProb2Jet(rootInt,  bTF_WBB)); //4
    
//     //Wcc ME, use same as Wbb ME but with a different TF
//     evtPrb2jetCEPD.push_back(new WbbEventProb2Jet(rootInt,  cTF_WCC)); //5
    
//     //Wc, use a Dummy class that is never evaluated for the 3 jets events
//     evtPrb2jetCEPD.push_back(new WcEventProb2Jet(divonneInt, cTF_WC, gluonTF_WC));//6
    
//     //Wjg, use a Dummy class that is never evaluated for the 3 jets events
//     evtPrb2jetCEPD.push_back(new WjgEventProb2Jet(divonneInt, lightTF_Wjg,gluonTF_Wjg)); //7
    
//     //Wgg, use a Dummy class that is never evaluated for the 3 jets events
//     evtPrb2jetCEPD.push_back(new WggEventProb2Jet(divonneInt, gluonTF_Wgg)); //8
    
//     // tt-bar matrix element 2-jet
//     evtPrb2jetCEPD.push_back(new ttEventProb2Jet(divonneInt, bTF_ttbar)); //9
//     evtPrb2jetCEPD.back()->setBounds(3, 0, MEConstants::beamEnergy);
//     evtPrb2jetCEPD.back()->setBounds(4, 0, TMath::TwoPi());
//     evtPrb2jetCEPD.back()->setBounds(5, 0, TMath::Pi());
    
//     // Diboson WW, use a Dummy class that is never evaluated for the 3 jets events
//     evtPrb2jetCEPD.push_back(new WWEventProb2Jet(rootInt, lightTF_diboson)); //10
    
//     // Diboson WZ, use a Dummy class that is never evaluated for the 3 jets events
//     evtPrb2jetCEPD.push_back(new WZEventProb2Jet(rootInt, lightTF_diboson));//11
    
//   }else { // DEBUG
//     //evtPrb2jetCEPD.push_back(new sChannelEventProb2Jet(divonneInt, bTF));
//     // In principle there is no need to set the bounds here as 
//     // long as it is done in the EPDEventProb2Jet  object below
//     evtPrb2jetCEPD.push_back(new ttEventProb2Jet(rootInt, bTF)); //18
//     evtPrb2jetCEPD.back()->setBounds(3, 0, MEConstants::beamEnergy);
//     evtPrb2jetCEPD.back()->setBounds(4, 0, TMath::TwoPi());
//     evtPrb2jetCEPD.back()->setBounds(5, 0, TMath::Pi());

//     if (0){
//       //eventProbs2jet.push_back(new sChannelEventProb2Jet(rootInt, bTF)); //11

//       eventProbs2jet.push_back(new ttEventProb2Jet(divonneInt, bTF)); //18
//       eventProbs2jet.back()->setBounds(3, 0, MEConstants::beamEnergy);
//       eventProbs2jet.back()->setBounds(4, 0, TMath::TwoPi());
//       eventProbs2jet.back()->setBounds(5, 0, TMath::Pi());
//     }
    
//  }

//   divonneInt.setEpsilon(.01, 0);
//   divonneInt.setNeval(10000, 1000000);

//   rootInt.setEpsilon(.10, 0);
//   rootInt.setNeval(10000, 100000);

//   //For each of the 11 Higgs masses
//   for (int mh=100; mh <101; mh += 5){

//     //Finnaly put the vector into the new EPD class, and the new EPD into the 
//     // MEJob vector 
//     eventProbs2jet.push_back(new EPDEventProb2Jet(divonneInt, evtPrb2jetCEPD,
// 						  bTF, 175, mh)); 
//     eventProbs2jet.back()->setBounds(3, 0, MEConstants::beamEnergy);
//     eventProbs2jet.back()->setBounds(4, 0, TMath::TwoPi());
//     eventProbs2jet.back()->setBounds(5, 0, TMath::Pi());
//   }

//   // ********************************************************
//   // END of ME's for Ricardo's correlated EPD
//   // ********************************************************



  // ***************************
  // **** Loading into MEJob
  // ***************************

  // Load the matrix elements for 2-jets into the job
  for (vector<EventProb*>::const_iterator it = eventProbs2jet.begin();
       it != eventProbs2jet.end(); ++it){

      (*it)->setBounds(0, -MEConstants::beamEnergy, MEConstants::beamEnergy);
      (*it)->setBounds(1, 0, MEConstants::beamEnergy);
      (*it)->setBounds(2, 0, MEConstants::beamEnergy);
      
      job.addEventProb(**it, 2);
  }
  
  // Load the matrix elements for 3-jets into the job
  for (vector<EventProb*>::const_iterator it = eventProbs3jet.begin();
       it != eventProbs3jet.end(); ++it)   {
    
    (*it)->setBounds(0, -MEConstants::beamEnergy, MEConstants::beamEnergy);
    (*it)->setBounds(1, 0, MEConstants::beamEnergy);
    (*it)->setBounds(2, 0, MEConstants::beamEnergy);
    
    job.addEventProb(**it, 3); // 3 means 3 jets
    
  }//for


  // ***************************
  // **** Command execution ****
  // ***************************
  job.loopOverEvents();
  // job.massScan(100, 150, 10, "Top");
  // job.massScan(100, 150, 50, "Higgs");
  // job.massScan(100, 250, 50,"Top");
  // job.massScan(50, 200, 50,"Top");


  // ***************************
  // ****   Cleaning up     ****
  // ***************************
 
   // Clean up for the 2-jet ME's
   for (vector<EventProb*>::const_iterator it = eventProbs2jet.begin();
        it != eventProbs2jet.end(); ++it)
     delete *it;

   // Clean up for the 3-jet ME's
   for (vector<EventProb*>::const_iterator it = eventProbs3jet.begin();
        it != eventProbs3jet.end(); ++it)
     delete *it;


   return 0;


}


