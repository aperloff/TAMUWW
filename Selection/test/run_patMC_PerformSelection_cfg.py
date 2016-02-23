#
# This file contains the Top PAG reference selection work-flow for mu + jets analysis.
# as defined in
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopLeptonPlusJetsRefSel_mu#Selection_Version_SelV4_valid_fr
#

import sys

import FWCore.ParameterSet.Config as cms

process = cms.Process( 'PAT' )


### ======================================================================== ###
###                                                                          ###
###                                 Constants                                ###
###                            (user job steering)                           ###
###                                                                          ###
### ======================================================================== ###


### Data or MC?
runOnMC = True

### Standard and PF work flow

# Standard
runStandardPAT = False
usePFJets      = True
useCaloJets    = False

# PF2PAT
runPF2PAT = True

### Switch on/off selection steps

# Step 0a
useTrigger      = False
# Step 0b
useGoodVertex   = True
# Step 1a
useLooseMuon    = False
# Step 1b
useTightMuon    = False
# Step 2
useMuonVeto     = False
# Step 3
useElectronVeto = False
# Step 4a
use1Jet         = False
# Step 4b
use2Jets        = False
# Step 4c
use3Jets        = False
# Step 5
use4Jets        = False
# Step 6
useBTag         = False

addTriggerMatching = False

### Reference selection

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets import *
# Muons general
muonsUsePV             = True
#muonEmbedTrack         = True
#muonJetsDR             = 0.3
# Standard mouns
muonCut                = 'isGlobalMuon && pt > 10. && abs(eta) < 2.5'
#looseMuonCut           = ''
#tightMuonCut           = ''
# PF muons
muonCutPF              = 'isGlobalMuon && pt > 10. && abs(eta) < 2.5'
#looseMuonCutPF         = ''
#tightMuonCutPF         = ''
# Standard electrons
electronCut            = 'et > 15. && abs(eta) < 2.5'
# PF electrons
electronCutPF          = 'et > 15. && abs(eta) < 2.5'
# Standard photons -> PFones are missing!
photonCut              = 'et > 20. && abs(eta) < 2.5'
# Calo jets
#jetCut                 = ''
jetCutCA8              = 'pt > 80.'
jetCutCA12             = 'pt > 80.'
jetCutAK7              = 'pt > 80.'
jetCutAK8              = 'pt > 80.'
# PF jets
#jetCutPF               = ''
#jetMuonsDRPF           = 0.1

# Fat Jets
isFatJetVal = False

# Trigger and trigger object
#triggerSelectionData       = ''
#triggerObjectSelectionData = ''
#triggerSelectionMC       = ''
#triggerObjectSelectionMC = ''

### Particle flow
### takes effect only, if 'runPF2PAT' = True

postfix = 'PFlow' # needs to be a non-empty string, if 'runStandardPAT' = True

# subtract charged hadronic pile-up particles (from wrong PVs)
# effects also JECs
usePFnoPU       = True # before any top projection
usePfIsoLessCHS = True # switch to new PF isolation with L1Fastjet CHS

# other switches for PF top projections (default: all 'True')
useNoMuon     = True # before electron top projection
useNoElectron = True # before jet top projection
useNoJet      = True # before tau top projection
useNoTau      = False # before MET top projection

# cuts used in top projections
# vertices
#pfVertices  = 'goodOfflinePrimaryVertices'
#pfD0Cut     = 0.2
#pfDzCut     = 0.5
# muons
#pfMuonSelectionCut = 'pt > 5.'
useMuonCutBasePF = False # use minimal (veto) muon selection cut on top of 'pfMuonSelectionCut'
pfMuonIsoConeR03 = False
pfElectronIsoConeR03 = True
#pfMuonCombIsoCut = 0.2
pfMuonCombIsoCut = 9999 #added on 2013/05/22
# electrons
#pfElectronSelectionCut  = 'pt > 5. && gsfTrackRef.isNonnull && gsfTrackRef.trackerExpectedHitsInner.numberOfLostHits < 2'
useElectronCutBasePF  = False # use minimal (veto) electron selection cut on top of 'pfElectronSelectionCut'
#pfElectronnIsoConeR03 = False
#pfElectronCombIsoCut  = 0.2
pfElectronCombIsoCut  = 9999 # added on 2013/03/19

### JEC levels

# levels to be accessible from the jets
# jets are corrected to L3Absolute (MC), L2L3Residual (data) automatically, if enabled here
# and remain uncorrected, if none of these levels is enabled here
useL1FastJet    = True  # needs useL1Offset being off, error otherwise
useL1Offset     = False # needs useL1FastJet being off, error otherwise
useL2Relative   = True
useL3Absolute   = True
useL2L3Residual = True  # takes effect only on data
useL5Flavor     = False
useL7Parton     = False
if runOnMC:
  inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'])
else:
  inputJetCorrLabel = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute' , 'L2L3Residual'])

### Input

# list of input files
useRelVals = False # if 'False', "inputFiles" is used
inputFiles = [] # overwritten, if "useRelVals" is 'True'

# maximum number of events
maxInputEvents = 1000 # reduce for testing

### Conditions

# GlobalTags (w/o suffix '::All')
globalTagData = 'GR_R_53_V10'
globalTagMC   = 'START53_V7E'

### Output

# output file
# outputFile = 'patRefSel_muJets.root'
outputFile = 'pat_53x_test_v03_data.root'

# event frequency of Fwk report
fwkReportEvery = 100

# switch for 'TrigReport'/'TimeReport' at job end
wantSummary = True


###                              End of constants                            ###
###                                                                          ###
### ======================================================================== ###


###
### Basic configuration
###

process.load( "TopQuarkAnalysis.Configuration.patRefSel_basics_cff" )
process.MessageLogger.cerr.FwkReport.reportEvery = fwkReportEvery
process.options.wantSummary = wantSummary
if runOnMC:
  process.GlobalTag.globaltag = globalTagMC   + '::All'
else:
  process.GlobalTag.globaltag = globalTagData + '::All'


###
### Input configuration
###

process.load( "TopQuarkAnalysis.Configuration.patRefSel_inputModule_cfi" )
if useRelVals:
  from PhysicsTools.PatAlgos.tools.cmsswVersionTools import pickRelValInputFiles
  if runOnMC:
    inputFiles = pickRelValInputFiles( cmsswVersion  = 'CMSSW_5_2_5_cand1'
                                     , relVal        = 'RelValTTbar'
                                     , globalTag     = 'START52_V9'
                                     , maxVersions   = 1
                                     )
  else:
    inputFiles = pickRelValInputFiles( cmsswVersion  = 'CMSSW_5_2_5_cand1'
                                     , relVal        = 'SingleMu'
                                     , dataTier      = 'RECO'
                                     , globalTag     = 'GR_R_52_V7_RelVal_mu2011B'
                                     , maxVersions   = 1
                                     )

if runOnMC:
  #inputFiles = cms.untracked.vstring('file:/afs/cern.ch/work/d/dimatteo/public/TestSamples/RelValTTbar.root')
  #inputFiles = cms.untracked.vstring('/store/mc/Summer12_DR53X/WH_ZH_TTH_HToWW_M-125_8TeV-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/0CE9D2F9-3CF5-E111-A26E-00215E21DC78.root')
  inputFiles = cms.untracked.vstring('/store/user/lnujj/PatTuples_8TeV_53X/custodio/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/SQWaT_PAT_53X_Summer12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_100_1_qoN.root')
else:
  inputFiles = cms.untracked.vstring('file:/afs/cern.ch/work/d/dimatteo/public/TestSamples/singleEle_532.root')         

process.source.fileNames = inputFiles
process.maxEvents.input  = maxInputEvents


###
### Output configuration
###

process.load( "TopQuarkAnalysis.Configuration.patRefSel_outputModule_cff" )
# output file name
process.out.fileName = outputFile
# event content
#from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
#process.out.outputCommands += patEventContent
# clear event selection
process.out.SelectEvents.SelectEvents = []


###
### Cleaning and trigger selection configuration
###

### Trigger selection
if runOnMC:
  triggerSelection = triggerSelectionMC
else:
  if useRelVals:
    triggerSelection = triggerSelectionDataRelVals
  else:
    triggerSelection = triggerSelectionData
from TopQuarkAnalysis.Configuration.patRefSel_triggerSelection_cff import triggerResults
process.step0a = triggerResults.clone(
  triggerConditions = [ triggerSelection ]
)

### Good vertex selection
process.load( "TopQuarkAnalysis.Configuration.patRefSel_goodVertex_cfi" )
process.step0b = process.goodOfflinePrimaryVertices.clone( filter = True )

### Event cleaning
process.load( 'TopQuarkAnalysis.Configuration.patRefSel_eventCleaning_cff' )
process.step0c = cms.Sequence()
if not runOnMC:
  process.step0c += process.eventCleaningData


###
### PAT/PF2PAT configuration
###

pfSuffix = 'PF'
if runStandardPAT and runPF2PAT:
  if postfix == '':
    sys.exit( 'ERROR: running standard PAT and PF2PAT in parallel requires a defined "postfix" for PF2PAT' )
  if usePFJets:
    if postfix == 'Add' + pfSuffix or postfix == jetAlgo + pfSuffix:
      sys.exit( 'ERROR: running standard PAT with additional PF jets  and PF2PAT in parallel does not allow for the "postfix" %s'%( postfix ) )
if not runStandardPAT and not runPF2PAT:
  sys.exit( 'ERROR: standard PAT and PF2PAT are both switched off' )

process.load( "PhysicsTools.PatAlgos.patSequences_cff" )
from PhysicsTools.PatAlgos.tools.coreTools import *

# Add PAT trigger information to the configuration
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, hltProcess = '*' )

### Check JECs

# JEC set
jecSet      = jecSetBase + 'Calo'
jecSetAddPF = jecSetBase + pfSuffix
jecSetPF    = jecSetAddPF
if usePFnoPU:
  jecSetPF += 'chs'

# JEC levels
if useL1FastJet and useL1Offset:
  sys.exit( 'ERROR: switch off either "L1FastJet" or "L1Offset"' )
jecLevels = []
if useL1FastJet:
  jecLevels.append( 'L1FastJet' )
if useL1Offset:
  jecLevels.append( 'L1Offset' )
if useL2Relative:
  jecLevels.append( 'L2Relative' )
if useL3Absolute:
  jecLevels.append( 'L3Absolute' )
if useL2L3Residual and not runOnMC:
  jecLevels.append( 'L2L3Residual' )
if useL5Flavor:
  jecLevels.append( 'L5Flavor' )
if useL7Parton:
  jecLevels.append( 'L7Parton' )

### Switch configuration

if runPF2PAT:
  if useMuonCutBasePF:
    pfMuonSelectionCut += ' && %s'%( muonCutBase )
  if useElectronCutBasePF:
    pfElectronSelectionCut += ' && %s'%( electronCutBase )
  from PhysicsTools.PatAlgos.tools.pfTools import usePF2PAT
  usePF2PAT( process
           , runPF2PAT      = runPF2PAT
           , runOnMC        = runOnMC
           , jetAlgo        = jetAlgo
           , postfix        = postfix
           , jetCorrections = ( jecSetPF
                              , jecLevels
                              )
           , pvCollection   = cms.InputTag( pfVertices )
           , typeIMetCorrections = True
           )
  applyPostfix( process, 'pfNoPileUp'  , postfix ).enable = usePFnoPU
  applyPostfix( process, 'pfNoMuon'    , postfix ).enable = useNoMuon
  applyPostfix( process, 'pfNoElectron', postfix ).enable = useNoElectron
  applyPostfix( process, 'pfNoJet'     , postfix ).enable = useNoJet
  applyPostfix( process, 'pfNoTau'     , postfix ).enable = useNoTau
  if useL1FastJet:
    applyPostfix( process, 'pfPileUp'   , postfix ).checkClosestZVertex = False
    applyPostfix( process, 'pfPileUpIso', postfix ).checkClosestZVertex = usePfIsoLessCHS
    applyPostfix( process, 'pfJets', postfix ).doAreaFastjet = True
    applyPostfix( process, 'pfJets', postfix ).doRhoFastjet  = False
  applyPostfix( process, 'pfMuonsFromVertex'    , postfix ).d0Cut    = pfD0Cut
  applyPostfix( process, 'pfMuonsFromVertex'    , postfix ).dzCut    = pfDzCut
  applyPostfix( process, 'pfSelectedMuons'      , postfix ).cut = pfMuonSelectionCut
  applyPostfix( process, 'pfIsolatedMuons'      , postfix ).isolationCut = pfMuonCombIsoCut
  if pfMuonIsoConeR03:
    applyPostfix( process, 'pfIsolatedMuons', postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'muPFIsoValueCharged03' + postfix )
                                                                                                  )
    applyPostfix( process, 'pfIsolatedMuons', postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'muPFIsoValuePU03' + postfix )
    applyPostfix( process, 'pfIsolatedMuons', postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
                                                                                                  , cms.InputTag( 'muPFIsoValueGamma03' + postfix )
                                                                                                  )
    applyPostfix( process, 'patMuons', postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'muPFIsoValueNeutral03' + postfix )
    applyPostfix( process, 'patMuons', postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'muPFIsoValuePU03' + postfix )
    applyPostfix( process, 'patMuons', postfix ).isolationValues.pfPhotons          = cms.InputTag( 'muPFIsoValueGamma03' + postfix )
    applyPostfix( process, 'patMuons', postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'muPFIsoValueCharged03' + postfix )
  applyPostfix( process, 'pfElectronsFromVertex'    , postfix ).d0Cut    = pfD0Cut
  applyPostfix( process, 'pfElectronsFromVertex'    , postfix ).dzCut    = pfDzCut
  applyPostfix( process, 'pfSelectedElectrons'      , postfix ).cut = pfElectronSelectionCut
  applyPostfix( process, 'pfIsolatedElectrons'      , postfix ).isolationCut = pfElectronCombIsoCut
  if pfElectronIsoConeR03:
    applyPostfix( process, 'pfIsolatedElectrons', postfix ).isolationValueMapsCharged  = cms.VInputTag( cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )
                                                                                                       )
    applyPostfix( process, 'pfIsolatedElectrons', postfix ).deltaBetaIsolationValueMap = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
    applyPostfix( process, 'pfIsolatedElectrons', postfix ).isolationValueMapsNeutral  = cms.VInputTag( cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
                                                                                                      , cms.InputTag( 'elPFIsoValueGamma03PFId'   + postfix )
                                                                                                      )
    applyPostfix( process, 'patElectrons', postfix ).isolationValues.pfNeutralHadrons   = cms.InputTag( 'elPFIsoValueNeutral03PFId' + postfix )
    applyPostfix( process, 'patElectrons', postfix ).isolationValues.pfPUChargedHadrons = cms.InputTag( 'elPFIsoValuePU03PFId' + postfix )
    applyPostfix( process, 'patElectrons', postfix ).isolationValues.pfPhotons          = cms.InputTag( 'elPFIsoValueGamma03PFId' + postfix )
    applyPostfix( process, 'patElectrons', postfix ).isolationValues.pfChargedHadrons   = cms.InputTag( 'elPFIsoValueCharged03PFId' + postfix )

from TopQuarkAnalysis.Configuration.patRefSel_refMuJets_cfi import *

# remove MC matching, object cleaning, objects etc.
jecLevelsCalo = copy.copy( jecLevels )
if runStandardPAT:
  if not runOnMC:
    runOnData( process )
  # subsequent jet area calculations needed for L1FastJet on RECO jets
  if useCaloJets and useL1FastJet:
    if useRelVals:
      process.ak5CaloJets = ak5CaloJets.clone( doAreaFastjet = True )
      process.ak5JetID    = ak5JetID.clone()
      process.ak5CaloJetSequence = cms.Sequence(
        process.ak5CaloJets
      * process.ak5JetID
      )
      from PhysicsTools.PatAlgos.tools.jetTools import switchJetCollection
      switchJetCollection( process
                         , cms.InputTag( jetAlgo.lower() + 'CaloJets' )
                         , doJTA            = True
                         , doBTagging       = True
                         , jetCorrLabel     = ( jecSet, jecLevels )
                         , doType1MET       = False
                         , genJetCollection = cms.InputTag( jetAlgo.lower() + 'GenJets' )
                         , doJetID          = True
                         )
    else:
      print 'WARNING patRefSel_muJets_test_cfg.py:'
      print '        L1FastJet JECs are not available for AK5Calo jets in this data due to missing jet area computation;'
      print '        switching to   L1Offset   !!!'
      jecLevelsCalo.insert( 0, 'L1Offset' )
      jecLevelsCalo.remove( 'L1FastJet' )
      process.patJetCorrFactors.levels = jecLevelsCalo
      #process.patJetCorrFactors.useRho = False # FIXME: does not apply
  if usePFJets:
    if useL1FastJet:
      process.ak5PFJets = ak5PFJets.clone( doAreaFastjet = True )
      

    from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
    from PhysicsTools.PatAlgos.tools.metTools import addPfMET
    addJetCollection( process
                    , cms.InputTag( jetAlgo.lower() + pfSuffix + 'Jets' )
                    , jetAlgo
                    , pfSuffix
                    , doJTA            = True
                    , doBTagging       = True
                    , jetCorrLabel     = ( jecSetAddPF, jecLevels )
                    , doType1MET       = False
                    , doL1Cleaning     = False
                    , doL1Counters     = True
                    , genJetCollection = cms.InputTag( jetAlgo.lower() + 'GenJets' )
                    , doJetID          = True
                    )
    
    addPfMET( process
            , jetAlgo + pfSuffix
            )
    removeSpecificPATObjects( process
                            , names = [ 'Taus' ]
                            ) # includes 'removeCleaning'
if runPF2PAT:
  if not runOnMC:
    runOnData( process
             , names = [ 'PFAll' ]
             , postfix = postfix
             )
    removeMCMatching( process, ['All'] )

  ## Photon PF isolation Producer - phoPFIso       
  process.load('EGamma.EGammaAnalysisTools.photonIsoProducer_cfi')

  if runOnMC and isFatJetVal:
    ####################################
    ###### gen jets for boosted ana #####
    ####################################
    
    process.load("RecoJets.Configuration.GenJetParticles_cff")
    from RecoJets.JetProducers.ca4GenJets_cfi import ca4GenJets
    from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
    process.ca8GenJetsNoNu = ca4GenJets.clone( rParam = cms.double(0.8),
                                             src = cms.InputTag("genParticlesForJetsNoNu"))
    process.ca12GenJetsNoNu = ca4GenJets.clone( rParam = cms.double(1.2),
                                              src = cms.InputTag("genParticlesForJetsNoNu"))
    
    process.ak7GenJetsNoNu = ak5GenJets.clone( rParam = cms.double(0.7),
                                           src = cms.InputTag("genParticlesForJetsNoNu"))      

    process.ak8GenJetsNoNu = ak5GenJets.clone( rParam = cms.double(0.8),
                                           src = cms.InputTag("genParticlesForJetsNoNu"))      
    
    process.additionalGenJets = cms.Sequence( process.genParticlesForJetsNoNu * process.ca8GenJetsNoNu * process.ca12GenJetsNoNu * process.ak7GenJetsNoNu * process.ak8GenJetsNoNu)
    

  ####################################
  ###### CA jets for boosted ana #####
  ####################################

  from RecoJets.JetProducers.ca4PFJets_cfi import ca4PFJets
    
  process.ca8PFJetsPFlow = ca4PFJets.clone(
    rParam = cms.double(0.8),
    src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax = cms.double(6.0),
    Ghost_EtaMax = cms.double(7.0)
  )
     
  process.ca12PFJetsPFlow = ca4PFJets.clone(
    rParam = cms.double(1.2),
    src = cms.InputTag('pfNoElectron'+postfix),
    doAreaFastjet = cms.bool(True),
    doRhoFastjet = cms.bool(True),
    Rho_EtaMax = cms.double(6.0),
    Ghost_EtaMax = cms.double(7.0)
  )

  ####################################
  ###### AKt jets for boosted ana #####
  ####################################
    
  process.ak7PFJetsPFlow = process.pfJetsPFlow.clone( 
    rParam = cms.double(0.7)
  )

  process.ak8PFJetsPFlow = process.pfJetsPFlow.clone( 
    rParam = cms.double(0.7)
  )

  if isFatJetVal:

    ###############################
    ###### AK 0.5 jets groomed ####
    ###############################

    from RecoJets.JetProducers.ak5PFJetsTrimmed_cfi import ak5PFJetsTrimmed
    process.ak5TrimmedPFlow = ak5PFJetsTrimmed.clone(
      src = process.pfJetsPFlow.src,
      doAreaFastjet = cms.bool(True)
      )

    from RecoJets.JetProducers.ak5PFJetsFiltered_cfi import ak5PFJetsFiltered
    process.ak5FilteredPFlow = ak5PFJetsFiltered.clone(
      src = process.pfJetsPFlow.src,
      doAreaFastjet = cms.bool(True)
      )

    from RecoJets.JetProducers.ak5PFJetsPruned_cfi import ak5PFJetsPruned
    process.ak5PrunedPFlow = ak5PFJetsPruned.clone(
      src = process.pfJetsPFlow.src,
      doAreaFastjet = cms.bool(True)
      )

    from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters

    ###############################
    ###### AK 0.7 jets groomed ####
    ###############################

    process.ak7TrimmedPFlow = process.ak5TrimmedPFlow.clone(
  	src = process.pfJetsPFlow.src, 
	rParam = cms.double(0.7)
    )

    process.ak7FilteredPFlow = process.ak5FilteredPFlow.clone(
	src = process.pfJetsPFlow.src,
	rParam = cms.double(0.7)
	)

    process.ak7PrunedPFlow = process.ak5PrunedPFlow.clone(
	src = process.pfJetsPFlow.src,
	rParam = cms.double(0.7)
    )


    ###############################
    ###### AK 0.8 jets groomed ####
    ###############################

    process.ak8TrimmedPFlow = process.ak5TrimmedPFlow.clone(
	src = process.pfJetsPFlow.src,
	rParam = cms.double(0.8)
    )

    process.ak8FilteredPFlow = process.ak5FilteredPFlow.clone(
	src = process.pfJetsPFlow.src,
	rParam = cms.double(0.8)
	)

    process.ak8PrunedPFlow = process.ak5PrunedPFlow.clone(
	src = process.pfJetsPFlow.src,
	rParam = cms.double(0.8)
    )

    ###############################
    ###### CA8 Pruning Setup ######
    ###############################


    # Pruned PF Jets
    process.caPrunedPFlow = process.ak5PrunedPFlow.clone(
	jetAlgorithm = cms.string("CambridgeAachen"),
	rParam       = cms.double(0.8)
    )


    ###############################
    ###### CA8 Filtered Setup #####
    ###############################


    # Filtered PF Jets
    process.caFilteredPFlow = ak5PFJetsFiltered.clone(
	src = cms.InputTag('pfNoElectron'+postfix),
	jetAlgorithm = cms.string("CambridgeAachen"),
	rParam       = cms.double(1.2),
	writeCompound = cms.bool(True),
	doAreaFastjet = cms.bool(True),
	jetPtMin = cms.double(100.0)
    )


  removeSpecificPATObjects( process
                          , names = [ 'Photons' ]
                          , postfix = postfix
                          ) # includes 'removeCleaning'

# JetCorrFactorsProducer configuration has to be fixed in standard work flow after a call to 'runOnData()':
if runStandardPAT:
  process.patJetCorrFactors.payload = jecSet
  process.patJetCorrFactors.levels  = jecLevelsCalo

# additional event content has to be (re-)added _after_ the call to 'removeCleaning()':
process.out.outputCommands += ['drop *'
#                              , 'keep *'
                              , 'keep *_ak5GenJets_*_SIM'
#                              , 'keep *_selectedPatPhotons__*'
                              , 'keep *_phoPFIso_*_*'
                              , 'keep *_photon*_*_*'
                              , 'keep *_selectedPatElectrons__*'
                              , 'keep *_selectedPatElectronsPFlow__*'
                              , 'keep *_selectedPatMuons__*'
                              , 'keep *_selectedPatMuonsPFlow__*'
                              , 'keep *_selectedPatJetsAK5PF_*_*'
                              , 'keep *_selectedPatJetsPFlow_*_*'
                              , 'keep *_patMETsAK5PF__*'
                              , 'keep *_patMETsPFlow__*'
                              , 'keep *_selectedPatPFParticlesPFlow__*'
                              , 'keep edmTriggerResults_*_*_*'
                              , 'keep *_hltTriggerSummaryAOD_*_*'
                              , 'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*'
                              # vertices and beam spot
                              , 'keep *_offlineBeamSpot_*_*'
                              , 'keep *_offlinePrimaryVertices*_*_*'
                              , 'keep *_goodOfflinePrimaryVertices*_*_*'
                              , 'keep patTriggerObjects_patTrigger_*_*'
                              , 'keep patTriggerFilters_patTrigger_*_*'
                              , 'keep patTriggerPaths_patTrigger_*_*'
                              , 'keep patTriggerEvent_patTriggerEvent_*_*'
                              , 'keep *_l1extraParticles_*_*'
                              , 'keep *_kt6PFJetsCentralChargedPileUp_rho_*'
                              , 'keep *_kt6PFJetsCentralNeutral_rho_*'
                              , 'keep *_kt6PFJetsCentralNeutralTight_rho_*'
                              , 'keep double_kt6PFJetsCentral_rho_RECO'
                              , 'keep *_patConversions*_*_*'
                              , 'keep *_allConversions__*'
                              , 'keep *_gsfElectron*__*'
                              , 'keep recoGsfTracks_electronGsfTracks__*'
                              , 'keep patMETs_patType1CorrectedPFMetPFlow__*'
                              , 'keep patMETs_patType1p2CorrectedPFMetPFlow__*'
                              , 'keep patMETs_patType1CorrectedPFMet__*' # misleading name: with smearing and shift correction
                              , 'keep patMETs_patType1p2CorrectedPFMet__*' # misleading name: with smearing and shift correction
                              , 'keep *_selectedSmearedPatJetsPFlow_*_*'
                              , 'keep *_pfMEtMVA_*_*'
                              , 'keep *_pu*JetId*_*_*'
                              , 'keep *_pu*JetMva*_*_*'
                              , 'keep *_selectedPatJetsCA8PF__*'
                              , 'keep *_selectedPatJetsCA12PF__*'
                              , 'keep *_selectedPatJetsAK7PF__*'
                              , 'keep *_selectedPatJetsAK8PF__*'
                              , 'keep *_pfInputsCA8_*_*'
                              , 'keep *_pfInputsCA12_*_*'
                              , 'keep *_pfInputsAK7_*_*'
                              , 'keep *_pfInputsAK8_*_*'
                              , 'keep *_patJetCorrFactorsPFlow_*_*'
                              , 'keep *_patJetCorrFactorsAK7PF_*_*'
                              , 'keep *_patJetCorrFactorsAK8PF_*_*'
                              , 'keep *_patJetCorrFactorsCA8PF_*_*'
                              , 'keep *_patJetCorrFactorsCA12PF_*_*'
                              , 'drop *_*PF_caloTowers_*'
                                ]
if isFatJetVal:
  process.out.outputCommands += [ 'keep *_selectedPatJetsCA8PrunedPF__*'
                              , 'keep *_selectedPatJetsCA12FilteredPF__*'
                              , 'keep *_selectedPatJetsAK7PrunedPF__*'
                              , 'keep *_selectedPatJetsAK7TrimmedPF__*'
                              , 'keep *_selectedPatJetsAK7FilteredPF__*'
                              , 'keep *_selectedPatJetsAK8PrunedPF__*'
                              , 'keep *_selectedPatJetsAK8TrimmedPF__*'
                              , 'keep *_selectedPatJetsAK8FilteredPF__*'
                              , 'keep recoGenJets_ca8GenJetsNoNu_*_*'
                              , 'keep recoGenJets_ca12GenJetsNoNu_*_*'
                              , 'keep recoGenJets_ak7GenJetsNoNu_*_*'
                              , 'keep recoGenJets_ak8GenJetsNoNu_*_*'
                                ]


keepSC = True
keepTK = True

if keepSC:
  process.out.outputCommands += [ 'keep *_correctedHybridSuperClusters_*_*'
                                , 'keep *_correctedMulti5x5SuperClustersWithPreshower_*_*'
                                ]
if keepTK:
  process.highPtTracks = cms.EDFilter (
    "TrackSelector",
    src = cms.InputTag ('generalTracks'),
    cut = cms.string ('pt > 15')
  )
#  process.out.outputCommands += [ 'keep *_generalTracks_*_*'
  process.out.outputCommands += [ 'keep *_highPtTracks_*_*'
                                ]
if runOnMC:
  process.out.outputCommands += [ 'keep GenEventInfoProduct_*_*_*'
                                , 'keep recoGenParticles_*_*_*'
                                , 'keep *_addPileupInfo_*_*'
                                , 'keep LHEEventProduct_*_*_*'
                                ]


###
### Additional configuration
###

goodPatJetsAddPFLabel = 'goodPatJets' + jetAlgo + pfSuffix

if runStandardPAT:

  ### Muons

  #process.intermediatePatMuons = intermediatePatMuons.clone()
  #process.loosePatMuons        = loosePatMuons.clone()
  process.step1a               = step1a.clone()
  #process.tightPatMuons        = tightPatMuons.clone()
  process.step1b               = step1b.clone()
  process.step2                = step2.clone( src = cms.InputTag( 'selectedPatMuons' ) )

  if usePFJets:
    loosePatMuonsAddPF = loosePatMuons.clone()
    loosePatMuonsAddPF.checkOverlaps.jets.src = cms.InputTag( goodPatJetsAddPFLabel )
    setattr( process, 'loosePatMuons' + jetAlgo + pfSuffix, loosePatMuonsAddPF )
    step1aAddPF = step1a.clone( src = cms.InputTag( 'loosePatMuons' + jetAlgo + pfSuffix ) )
    setattr( process, 'step1a' + jetAlgo + pfSuffix, step1aAddPF )
    tightPatMuonsAddPF = tightPatMuons.clone( src = cms.InputTag( 'loosePatMuons' + jetAlgo + pfSuffix ) )
    setattr( process, 'tightPatMuons' + jetAlgo + pfSuffix, tightPatMuonsAddPF )
    step1bAddPF = step1b.clone( src = cms.InputTag( 'tightPatMuons' + jetAlgo + pfSuffix ) )
    setattr( process, 'step1b' + jetAlgo + pfSuffix, step1bAddPF )

  ### Jets

  process.kt6PFJets = kt6PFJets.clone( src          = cms.InputTag( 'particleFlow' )
                                     , doRhoFastjet = True
                                     )
                                     
                                   
  # compute FastJet rho to correct isolation
  #process.kt6PFJetsForIsolation = kt6PFJets.clone()
  #process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)
                
  #process.patDefaultSequence.replace( process.patJetCorrFactors
  #                                  , process.kt6PFJets * process.kt6PFJetsForIsolation * process.patJetCorrFactors
  #                                  )
  process.patDefaultSequence.replace( process.patJetCorrFactors
                                     , process.kt6PFJets * process.patJetCorrFactors
                                     )
  #process.out.outputCommands.append( 'keep double_kt6PFJets_*_' + process.name_() )
  process.out.outputCommands.append( 'keep double_*_*_' + process.name_() )
  if useL1FastJet:
    process.patJetCorrFactors.useRho = True
    if usePFJets:
      getattr( process, 'patJetCorrFactors' + jetAlgo + pfSuffix ).useRho = True

  process.goodPatJets = goodPatJets.clone()
  process.step4a      = step4a.clone()
  process.step4b      = step4b.clone()
  process.step4c      = step4c.clone()
  process.step5       = step5.clone()

  if usePFJets:
    goodPatJetsAddPF = goodPatJets.clone( src = cms.InputTag( 'selectedPatJets' + jetAlgo + pfSuffix ) )
    setattr( process, goodPatJetsAddPFLabel, goodPatJetsAddPF )
    step4aAddPF = step4a.clone( src = cms.InputTag( goodPatJetsAddPFLabel ) )
    setattr( process, 'step4a' + jetAlgo + pfSuffix, step4aAddPF )
    step4bAddPF = step4b.clone( src = cms.InputTag( goodPatJetsAddPFLabel ) )
    setattr( process, 'step4b' + jetAlgo + pfSuffix, step4bAddPF )
    step4cAddPF = step4c.clone( src = cms.InputTag( goodPatJetsAddPFLabel ) )
    setattr( process, 'step4c' + jetAlgo + pfSuffix, step4cAddPF )
    step5AddPF = step5.clone( src = cms.InputTag( goodPatJetsAddPFLabel ) )
    setattr( process, 'step5' + jetAlgo + pfSuffix, step5AddPF )

  ### Electrons

  process.step3 = step3.clone( src = cms.InputTag( 'selectedPatElectrons' ) )

if runPF2PAT:

  ### Muons

  #intermediatePatMuonsPF = intermediatePatMuons.clone( src = cms.InputTag( 'selectedPatMuons' + postfix ) )
  #setattr( process, 'intermediatePatMuons' + postfix, intermediatePatMuonsPF )

  #loosePatMuonsPF = loosePatMuons.clone( src = cms.InputTag( 'intermediatePatMuons' + postfix ) )
  #setattr( process, 'loosePatMuons' + postfix, loosePatMuonsPF )
  #getattr( process, 'loosePatMuons' + postfix ).checkOverlaps.jets.src = cms.InputTag( 'goodPatJets' + postfix )

  #step1aPF = step1a.clone( src = cms.InputTag( 'loosePatMuons' + postfix ) )
  #setattr( process, 'step1a' + postfix, step1aPF )

  #tightPatMuonsPF = tightPatMuons.clone( src = cms.InputTag( 'loosePatMuons' + postfix ) )
  #setattr( process, 'tightPatMuons' + postfix, tightPatMuonsPF )

  #step1bPF = step1b.clone( src = cms.InputTag( 'tightPatMuons' + postfix ) )
  #setattr( process, 'step1b' + postfix, step1bPF )

  #step2PF = step2.clone( src = cms.InputTag( 'selectedPatMuons' + postfix ) )
  #setattr( process, 'step2' + postfix, step2PF )

  ### Jets

  kt6PFJetsPF = kt6PFJets.clone( doRhoFastjet = True )
  setattr( process, 'kt6PFJets' + postfix, kt6PFJetsPF )
  getattr( process, 'patPF2PATSequence' + postfix).replace( getattr( process, 'pfNoElectron' + postfix )
                                                          , getattr( process, 'pfNoElectron' + postfix ) * getattr( process, 'kt6PFJets' + postfix )
                                                          )
  for module in (
        getattr(process,"ca8PFJets" + postfix),
        getattr(process,"ca12PFJets" + postfix),
        getattr(process,"ak7PFJets" + postfix),
        getattr(process,"ak8PFJets" + postfix),
        ) :
        getattr(process,"patPF2PATSequence"+postfix).replace( getattr(process,"pfNoElectron"+postfix), getattr(process,"pfNoElectron"+postfix)*module )

  if isFatJetVal:
    for module in (
        getattr(process,"ak5Trimmed" + postfix),
        getattr(process,"ak5Filtered" + postfix),
        getattr(process,"ak5Pruned" + postfix),
        getattr(process,"ak7Trimmed" + postfix),
        getattr(process,"ak7Filtered" + postfix),
        getattr(process,"ak7Pruned" + postfix),
        getattr(process,"ak8Trimmed" + postfix),
        getattr(process,"ak8Filtered" + postfix),
        getattr(process,"ak8Pruned" + postfix),
        getattr(process,"caPruned" + postfix),
        getattr(process,"caFiltered" + postfix),
        ) :
        getattr(process,"patPF2PATSequence"+postfix).replace( getattr(process,"pfNoElectron"+postfix), getattr(process,"pfNoElectron"+postfix)*module )

  from PhysicsTools.PatAlgos.tools.pfTools import *

  addJetCollection(process,
              cms.InputTag('ca8PFJetsPFlow'), # Jet collection; must be already in the event when patLayer0 sequence is executed
              'CA8', 'PF',
              doJTA=True, # Run Jet-Track association & JetCharge
              doBTagging=True, # Run b-tagging
              jetCorrLabel= ( jecSetAddPF, jecLevels ),
              doType1MET=True,
              doL1Cleaning=False,
              doL1Counters=False,
              #genJetCollection = cms.InputTag("ca8GenJetsNoNu"),
              doJetID = False
              )

  addJetCollection(process,
              cms.InputTag('ca12PFJetsPFlow'), # Jet collection; must be already in the event when patLayer0 sequence is executed
              'CA12', 'PF',
              doJTA=True, # Run Jet-Track association & JetCharge
              doBTagging=True, # Run b-tagging
              jetCorrLabel= ( jecSetAddPF, jecLevels ),
              doType1MET=True,
              doL1Cleaning=False,
              doL1Counters=False,
              #genJetCollection = cms.InputTag("ca12GenJetsNoNu"),
              doJetID = False
              )

  addJetCollection(process,
              cms.InputTag('ak7PFJetsPFlow'), # Jet collection; must be already in the event when patLayer0 sequence is executed
              'AK7', 'PF',
              doJTA=True, # Run Jet-Track association & JetCharge
              doBTagging=True, # Run b-tagging
              jetCorrLabel= ( jecSetAddPF, jecLevels ),
              doType1MET=True,
              doL1Cleaning=False,
              doL1Counters=False,
              #genJetCollection = cms.InputTag("ak7GenJetsNoNu"),
              doJetID = False
              )

  addJetCollection(process, 
			 cms.InputTag('ak8PFJetsPFlow'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
			 'AK8', 'PF',
			 doJTA=False,            # Run Jet-Track association & JetCharge
			 doBTagging=False,       # Run b-tagging
			 jetCorrLabel=inputJetCorrLabel,
			 doType1MET=False,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 #genJetCollection = cms.InputTag("ak8GenJetsNoNu"),
			 doJetID = False
			 )

  if isFatJetVal:

    addJetCollection(process, 
                 cms.InputTag('caPrunedPFlow'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
                 'CA8Pruned', 'PF',
                 doJTA=False,            # Run Jet-Track association & JetCharge
                 doBTagging=False,       # Run b-tagging
                 jetCorrLabel=inputJetCorrLabel,
                 doType1MET=False,
                 doL1Cleaning=False,
                 doL1Counters=False,
                 #genJetCollection = cms.InputTag("ca8GenJetsNoNu"),
                 doJetID = False
                 )

    addJetCollection(process, 
			 cms.InputTag('caFilteredPFlow'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
			 'CA12Filtered', 'PF',
			 doJTA=False,            # Run Jet-Track association & JetCharge
			 doBTagging=False,       # Run b-tagging
			 jetCorrLabel=inputJetCorrLabel,
			 doType1MET=False,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 #genJetCollection = cms.InputTag("ca8GenJetsNoNu"),
			 doJetID = False
			 )


    addJetCollection(process, 
			 cms.InputTag('ak5PrunedPFlow'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
			 'AK5Pruned', 'PF',
			 doJTA=False,            # Run Jet-Track association & JetCharge
			 doBTagging=False,       # Run b-tagging
			 jetCorrLabel=inputJetCorrLabel,
			 doType1MET=False,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 #genJetCollection = cms.InputTag("ak5GenJetsNoNu"),
			 doJetID = False
			 )

    addJetCollection(process, 
			 cms.InputTag('ak5FilteredPFlow'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
			 'AK5Filtered', 'PF',
			 doJTA=False,            # Run Jet-Track association & JetCharge
			 doBTagging=False,       # Run b-tagging
			 jetCorrLabel=inputJetCorrLabel,
			 doType1MET=False,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 #genJetCollection = cms.InputTag("ak5GenJetsNoNu"),
			 doJetID = False
			 )
    addJetCollection(process, 
			 cms.InputTag('ak5TrimmedPFlow'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
			 'AK5Trimmed', 'PF',
			 doJTA=False,            # Run Jet-Track association & JetCharge
			 doBTagging=False,       # Run b-tagging
			 jetCorrLabel=inputJetCorrLabel,
			 doType1MET=False,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 #genJetCollection = cms.InputTag("ak5GenJetsNoNu"),
			 doJetID = False
			 )

    addJetCollection(process, 
			 cms.InputTag('ak7PrunedPFlow'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
			 'AK7Pruned', 'PF',
			 doJTA=False,            # Run Jet-Track association & JetCharge
			 doBTagging=False,       # Run b-tagging
			 jetCorrLabel=inputJetCorrLabel,
			 doType1MET=False,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 #genJetCollection = cms.InputTag("ak7GenJetsNoNu"),
			 doJetID = False
			 )
    addJetCollection(process, 
			 cms.InputTag('ak7FilteredPFlow'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
			 'AK7Filtered', 'PF',
			 doJTA=False,            # Run Jet-Track association & JetCharge
			 doBTagging=False,       # Run b-tagging
			 jetCorrLabel=inputJetCorrLabel,
			 doType1MET=False,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 #genJetCollection = cms.InputTag("ak7GenJetsNoNu"),
			 doJetID = False
			 )
    addJetCollection(process, 
			 cms.InputTag('ak7TrimmedPFlow'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
			 'AK7Trimmed', 'PF',
			 doJTA=False,            # Run Jet-Track association & JetCharge
			 doBTagging=False,       # Run b-tagging
			 jetCorrLabel=inputJetCorrLabel,
			 doType1MET=False,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 #genJetCollection = cms.InputTag("ak7GenJetsNoNu"),
			 doJetID = False
			 )


    addJetCollection(process, 
			 cms.InputTag('ak8PrunedPFlow'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
			 'AK8Pruned', 'PF',
			 doJTA=False,            # Run Jet-Track association & JetCharge
			 doBTagging=False,       # Run b-tagging
			 jetCorrLabel=inputJetCorrLabel,
			 doType1MET=False,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 #genJetCollection = cms.InputTag("ak8GenJetsNoNu"),
			 doJetID = False
			 )

    addJetCollection(process, 
			 cms.InputTag('ak8FilteredPFlow'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
			 'AK8Filtered', 'PF',
			 doJTA=False,            # Run Jet-Track association & JetCharge
			 doBTagging=False,       # Run b-tagging
			 jetCorrLabel=inputJetCorrLabel,
			 doType1MET=False,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 #genJetCollection = cms.InputTag("ak8GenJetsNoNu"),
			 doJetID = False
			 )
    addJetCollection(process, 
			 cms.InputTag('ak8TrimmedPFlow'),         # Jet collection; must be already in the event when patLayer0 sequence is executed
			 'AK8Trimmed', 'PF',
			 doJTA=False,            # Run Jet-Track association & JetCharge
			 doBTagging=False,       # Run b-tagging
			 jetCorrLabel=inputJetCorrLabel,
			 doType1MET=False,
			 doL1Cleaning=False,
			 doL1Counters=False,
			 #genJetCollection = cms.InputTag("ak8GenJetsNoNu"),
			 doJetID = False
			 )


##  from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector
##  process.goodPatJetsPFlow = cms.EDFilter("PFJetIDSelectionFunctorFilter",
##                                        filterParams = pfJetIDSelector.clone(),
##                                        src = cms.InputTag("selectedPatJetsPFlow")
##                                        )


  process.pfInputsCA8 = cms.EDProducer(
      "CandViewNtpProducer", 
      src = cms.InputTag('selectedPatJetsCA8PF', 'pfCandidates'),
      lazyParser = cms.untracked.bool(True),
      eventInfo = cms.untracked.bool(False),
      variables = cms.VPSet(
          cms.PSet(
              tag = cms.untracked.string("px"),
              quantity = cms.untracked.string("px")
              ),
          cms.PSet(
              tag = cms.untracked.string("py"),
              quantity = cms.untracked.string("py")
              ),
          cms.PSet(
              tag = cms.untracked.string("pz"),
              quantity = cms.untracked.string("pz")
              ),
          cms.PSet(
              tag = cms.untracked.string("energy"),
              quantity = cms.untracked.string("energy")
              ),
          cms.PSet(
              tag = cms.untracked.string("pdgId"),
              quantity = cms.untracked.string("pdgId")
              )
          )
  )
  
  process.pfInputsCA12 = cms.EDProducer(
      "CandViewNtpProducer", 
      src = cms.InputTag('selectedPatJetsCA12PF', 'pfCandidates'),
      lazyParser = cms.untracked.bool(True),
      eventInfo = cms.untracked.bool(False),
      variables = cms.VPSet(
          cms.PSet(
              tag = cms.untracked.string("px"),
              quantity = cms.untracked.string("px")
              ),
          cms.PSet(
              tag = cms.untracked.string("py"),
              quantity = cms.untracked.string("py")
              ),
          cms.PSet(
              tag = cms.untracked.string("pz"),
              quantity = cms.untracked.string("pz")
              ),
          cms.PSet(
              tag = cms.untracked.string("energy"),
              quantity = cms.untracked.string("energy")
              ),
          cms.PSet(
              tag = cms.untracked.string("pdgId"),
              quantity = cms.untracked.string("pdgId")
              )
          )
  )

  process.pfInputsAK7 = cms.EDProducer(
      "CandViewNtpProducer",
      src = cms.InputTag('selectedPatJetsAK7PF', 'pfCandidates'),
      lazyParser = cms.untracked.bool(True),
      eventInfo = cms.untracked.bool(False),
      variables = cms.VPSet(
          cms.PSet(
              tag = cms.untracked.string("px"),
              quantity = cms.untracked.string("px")
              ),
          cms.PSet(
              tag = cms.untracked.string("py"),
              quantity = cms.untracked.string("py")
              ),
          cms.PSet(
              tag = cms.untracked.string("pz"),
              quantity = cms.untracked.string("pz")
              ),
          cms.PSet(
              tag = cms.untracked.string("energy"),
              quantity = cms.untracked.string("energy")
              ),
          cms.PSet(
              tag = cms.untracked.string("pdgId"),
              quantity = cms.untracked.string("pdgId")
              )
          )
  )

  process.pfInputsAK8 = cms.EDProducer(
      "CandViewNtpProducer",
      src = cms.InputTag('selectedPatJetsAK8PF', 'pfCandidates'),
      lazyParser = cms.untracked.bool(True),
      eventInfo = cms.untracked.bool(False),
      variables = cms.VPSet(
          cms.PSet(
              tag = cms.untracked.string("px"),
              quantity = cms.untracked.string("px")
              ),
          cms.PSet(
              tag = cms.untracked.string("py"),
              quantity = cms.untracked.string("py")
              ),
          cms.PSet(
              tag = cms.untracked.string("pz"),
              quantity = cms.untracked.string("pz")
              ),
          cms.PSet(
              tag = cms.untracked.string("energy"),
              quantity = cms.untracked.string("energy")
              ),
          cms.PSet(
              tag = cms.untracked.string("pdgId"),
              quantity = cms.untracked.string("pdgId")
              )
          )
  )


  process.selectedPatJetsCA8PF.cut = jetCutCA8
  process.selectedPatJetsCA12PF.cut = jetCutCA12
  process.selectedPatJetsAK7PF.cut = jetCutAK7
  process.selectedPatJetsAK8PF.cut = jetCutAK8
  if isFatJetVal:
    process.selectedPatJetsCA8PrunedPF.cut = jetCutAK7
    process.selectedPatJetsCA12FilteredPF.cut = jetCutAK7
    process.selectedPatJetsAK5PrunedPF.cut = jetCutAK7
    process.selectedPatJetsAK5FilteredPF.cut = jetCutAK7
    process.selectedPatJetsAK5TrimmedPF.cut = jetCutAK7
    process.selectedPatJetsAK7PrunedPF.cut = jetCutAK7
    process.selectedPatJetsAK7FilteredPF.cut = jetCutAK7
    process.selectedPatJetsAK7TrimmedPF.cut = jetCutAK7
    process.selectedPatJetsAK8PrunedPF.cut = jetCutAK7
    process.selectedPatJetsAK8FilteredPF.cut = jetCutAK7
    process.selectedPatJetsAK8TrimmedPF.cut = jetCutAK7

  for icorr in [process.patJetCorrFactorsCA8PF,
	      process.patJetCorrFactorsAK7PF,
	      process.patJetCorrFactorsAK8PF,
              process.patJetCorrFactorsCA12PF] :
    icorr.rho = cms.InputTag('kt6PFJets' + postfix, 'rho')

  for module in [process.patJetCorrFactorsCA8PF,
              process.patJetCorrFactorsAK7PF,
              process.patJetCorrFactorsAK8PF,
              process.patJetCorrFactorsCA12PF]:
    module.primaryVertices = "goodOfflinePrimaryVertices"

  if isFatJetVal:
    for icorr in [process.patJetCorrFactorsCA8PrunedPF,
              process.patJetCorrFactorsCA12FilteredPF,
              process.patJetCorrFactorsAK5PrunedPF,
	      process.patJetCorrFactorsAK5FilteredPF,
	      process.patJetCorrFactorsAK5TrimmedPF,
	      process.patJetCorrFactorsAK7PrunedPF,
	      process.patJetCorrFactorsAK7FilteredPF,
	      process.patJetCorrFactorsAK7TrimmedPF,
	      process.patJetCorrFactorsAK8PrunedPF,
	      process.patJetCorrFactorsAK8FilteredPF,
	      process.patJetCorrFactorsAK8TrimmedPF] :
              icorr.rho = cms.InputTag('kt6PFJets' + postfix, 'rho')

    for module in  [process.patJetCorrFactorsCA8PrunedPF,
              process.patJetCorrFactorsCA12FilteredPF,
              process.patJetCorrFactorsAK5PrunedPF,
              process.patJetCorrFactorsAK5FilteredPF,
              process.patJetCorrFactorsAK5TrimmedPF,
              process.patJetCorrFactorsAK7PrunedPF,
              process.patJetCorrFactorsAK7FilteredPF,
              process.patJetCorrFactorsAK7TrimmedPF,
              process.patJetCorrFactorsAK8PrunedPF,
              process.patJetCorrFactorsAK8FilteredPF,
              process.patJetCorrFactorsAK8TrimmedPF] :
      module.primaryVertices = "goodOfflinePrimaryVertices"

  if useL1FastJet:
    applyPostfix( process, 'patJetCorrFactors', postfix ).rho = cms.InputTag( 'kt6PFJets' + postfix, 'rho' )
  process.out.outputCommands.append( 'keep double_kt6PFJets' + postfix + '_*_' + process.name_() )

  #goodPatJetsPF = goodPatJets.clone( src = cms.InputTag( 'selectedPatJets' + postfix ) )
  #setattr( process, 'goodPatJets' + postfix, goodPatJetsPF )
  #getattr( process, 'goodPatJets' + postfix ).checkOverlaps.muons.src = cms.InputTag( 'intermediatePatMuons' + postfix )

  #step4aPF = step4a.clone( src = cms.InputTag( 'goodPatJets' + postfix ) )
  #setattr( process, 'step4a' + postfix, step4aPF )
  #step4bPF = step4b.clone( src = cms.InputTag( 'goodPatJets' + postfix ) )
  #setattr( process, 'step4b' + postfix, step4bPF )
  #step4cPF = step4c.clone( src = cms.InputTag( 'goodPatJets' + postfix ) )
  #setattr( process, 'step4c' + postfix, step4cPF )
  #step5PF = step5.clone( src = cms.InputTag( 'goodPatJets' + postfix ) )
  #setattr( process, 'step5'  + postfix, step5PF  )

  ### Electrons

  step3PF = step3.clone( src = cms.InputTag( 'selectedPatElectrons' + postfix ) )
  setattr( process, 'step3' + postfix, step3PF )

process.out.outputCommands.append( 'keep *_intermediatePatMuons*_*_*' )
process.out.outputCommands.append( 'keep *_loosePatMuons*_*_*' )
process.out.outputCommands.append( 'keep *_tightPatMuons*_*_*' )
process.out.outputCommands.append( 'keep *_goodPatJets*_*_*' )


###
### Selection configuration
###

if runStandardPAT:

  ### Muons

  process.patMuons.usePV      = muonsUsePV
  process.patMuons.embedTrack = muonEmbedTrack

  process.selectedPatMuons.cut = muonCut

  #process.intermediatePatMuons.preselection = looseMuonCut

  #process.loosePatMuons.checkOverlaps.jets.deltaR = muonJetsDR
  if usePFJets:
    print "not using loose muons"
    #getattr( process, 'loosePatMuons' + jetAlgo + pfSuffix ).checkOverlaps.jets.deltaR = muonJetsDR
  
  #process.tightPatMuons.preselection = tightMuonCut
  if usePFJets:
    print "not using tight muons"
    #getattr( process, 'tightPatMuons' + jetAlgo + pfSuffix ).preselection = tightMuonCut

  ### Jets

  process.goodPatJets.preselection = jetCut
  if usePFJets:
    getattr( process, goodPatJetsAddPFLabel ).preselection               = jetCutPF
    getattr( process, goodPatJetsAddPFLabel ).checkOverlaps.muons.deltaR = jetMuonsDRPF

  ### Electrons

  process.patElectrons.electronIDSources = electronIDSources

  process.selectedPatElectrons.cut = electronCut

  ### Photons

  process.selectedPatPhotons.cut = photonCut

if runPF2PAT:

  applyPostfix( process, 'patMuons', postfix ).usePV      = muonsUsePV
  applyPostfix( process, 'patMuons', postfix ).embedTrack = muonEmbedTrack

  applyPostfix( process, 'selectedPatMuons', postfix ).cut = muonCutPF

  #getattr( process, 'intermediatePatMuons' + postfix ).preselection = looseMuonCutPF

  #getattr( process, 'loosePatMuons' + postfix ).preselection              = looseMuonCutPF
  #getattr( process, 'loosePatMuons' + postfix ).checkOverlaps.jets.deltaR = muonJetsDR

  #getattr( process, 'tightPatMuons' + postfix ).preselection = tightMuonCutPF

  ### Jets

  #getattr( process, 'goodPatJets' + postfix ).preselection               = jetCutPF
  #getattr( process, 'goodPatJets' + postfix ).checkOverlaps.muons.deltaR = jetMuonsDRPF

  ### Electrons

  applyPostfix( process, 'patElectrons', postfix ).electronIDSources = electronIDSources

  applyPostfix( process, 'selectedPatElectrons', postfix ).cut = electronCutPF


###
### Trigger matching
###

if addTriggerMatching:

  if runOnMC:
    triggerObjectSelection = triggerObjectSelectionMC
  else:
    if useRelVals:
      triggerObjectSelection = triggerObjectSelectionDataRelVals
    else:
      triggerObjectSelection = triggerObjectSelectionData
  ### Trigger matching configuration
  from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger
  from TopQuarkAnalysis.Configuration.patRefSel_triggerMatching_cfi import patMuonTriggerMatch
  from PhysicsTools.PatAlgos.tools.trigTools import *
  if runStandardPAT:
    triggerProducer = patTrigger.clone()
    setattr( process, 'patTrigger', triggerProducer )
    process.triggerMatch = patMuonTriggerMatch.clone( matchedCuts = triggerObjectSelection )
    switchOnTriggerMatchEmbedding( process
                                 , triggerMatchers = [ 'triggerMatch' ]
                                 )
    removeCleaningFromTriggerMatching( process )
    #process.intermediatePatMuons.src = cms.InputTag( 'selectedPatMuonsTriggerMatch' )
  if runPF2PAT:
    triggerProducerPF = patTrigger.clone()
    setattr( process, 'patTrigger' + postfix, triggerProducerPF )
    triggerMatchPF = patMuonTriggerMatch.clone( matchedCuts = triggerObjectSelection )
    setattr( process, 'triggerMatch' + postfix, triggerMatchPF )
    switchOnTriggerMatchEmbedding( process
                                 , triggerProducer = 'patTrigger' + postfix
                                 , triggerMatchers = [ 'triggerMatch' + postfix ]
                                 , sequence        = 'patPF2PATSequence' + postfix
                                 , postfix         = postfix
                                 )
    removeCleaningFromTriggerMatching( process
                                     , sequence = 'patPF2PATSequence' + postfix
                                     )
    #getattr( process, 'intermediatePatMuons' + postfix ).src = cms.InputTag( 'selectedPatMuons' + postfix + 'TriggerMatch' )


###
### Scheduling
###

# MVA electron ID

process.load( "EGamma.EGammaAnalysisTools.electronIdMVAProducer_cfi" )
process.eidMVASequence = cms.Sequence(
  process.mvaTrigV0
+ process.mvaNonTrigV0
)

#MVA met sequence
from JetMETCorrections.Configuration.JetCorrectionServices_cff import *
from JetMETCorrections.Configuration.DefaultJEC_cff import *
if runOnMC:
  jetCorrectorsForMVAMET = 'ak5PFL1FastL2L3'
else:
  jetCorrectorsForMVAMET = 'ak5PFL1FastL2L3Residual'
process.load('RecoMET.METProducers.mvaPFMET_cff')
#process.load('SQWaT.TopPAT.metProducerSequence_cff')
process.calibratedAK5PFJetsForPFMEtMVA.correctors = cms.vstring(jetCorrectorsForMVAMET) # NOTE: use "ak5PFL1FastL2L3" for MC / "ak5PFL1FastL2L3Residual" for Data

#search for the hardest lepton
process.goodLeptons = cms.EDProducer("CandViewMerger",
  src = cms.VInputTag(cms.InputTag("selectedPatElectronsPFlow"),cms.InputTag("selectedPatMuonsPFlow"))
  )
process.largestGoodLeptonCand = cms.EDFilter("LargestPtCandSelector",
   src = cms.InputTag("goodLeptons"),
   maxNumber = cms.uint32( 1 )
 )

process.pfMEtMVA.srcLeptons = cms.VInputTag( 'largestGoodLeptonCand' )
process.pfMEtMVA.verbosity = cms.int32(0)

#Jet Met smearing and shift correction
process.load("JetMETCorrections.Type1MET.pfMETsysShiftCorrections_cfi")
if runOnMC:
  process.pfMEtSysShiftCorrParameters = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_mc
  doJetSmearing = True
  jetCorrectorsForMetUnc = 'L3Absolute'
else:
  process.pfMEtSysShiftCorrParameters = process.pfMEtSysShiftCorrParameters_2012runAvsNvtx_data
  doJetSmearing = False
  jetCorrectorsForMetUnc = 'L2L3Residual'

from SQWaT.TopPAT.metUncertaintyTools import runMEtUncertainties
runMEtUncertainties(process, 
                    'patElectronsPFlow',    #Electron input collection
                    None,                   #Photon input collection
                    'patMuonsPFlow',        #Muon input collection
                    'patTausPFlow',         #Tau input collection
                    'patJetsPFlow',         #Jet input collection 
                    0.3,                    #DR cone for Jet-lepton cleaning
                    jetCorrectorsForMetUnc, #Jet correction label
                    doJetSmearing,          #Smear jets 
                    True,                   #Apply type 1 correction 
                    False,                  #Apply type 1 + type 2 correction
                    False,                  #Apply the sequence to MVA MET
                    False,                  #Apply type 0 correction
                    process.pfMEtSysShiftCorrParameters, #MET shift params 
                    True,                   #Apply MET shift correction
                    'PhysicsTools/PatUtils/data/pfJetResolutionMCtoDataCorrLUT.root', #Jet smearing file name
                    'pfJetResolutionMCtoDataCorrLUT', #Jet smearing histo
                    'particleFlow',      #PF input candidates
                    'AK5PF',             #Just a name? verify
                    1.0,                 #Sigma variations
                    False,               #Modify automatically the pat sequence 
                    None )               #out process name

#Select the smeared jets with pt > 15
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
process.selectedSmearedPatJetsPFlow = selectedPatJets.clone(src='smearedPatJetsPFlow', cut='pt>15')

#Pileup jet id
process.load("SQWaT.TopPAT.pujetidsequence_cff")
  
# The additional sequence

if runStandardPAT:
  process.patAddOnSequence = cms.Sequence(
  #process.intermediatePatMuons
  #* process.goodPatJets
  #* process.loosePatMuons
  #* process.tightPatMuons
  process.pfInputsCA8 *
  process.pfInputsCA12 *
  process.pfInputsAK7
  )
if runPF2PAT:
  patAddOnSequence = cms.Sequence(
    #getattr( process, 'intermediatePatMuons' + postfix )
  #* getattr( process, 'goodPatJets' + postfix )
  #* getattr( process, 'loosePatMuons' + postfix )
  #* getattr( process, 'tightPatMuons' + postfix )
  )
  setattr( process, 'patAddOnSequence' + postfix, patAddOnSequence )

  process.pfInputsSequence = cms.Sequence(
   process.pfInputsCA8 *
   process.pfInputsCA12 *
   process.pfInputsAK7 *
   process.pfInputsAK8
  )

if runPF2PAT:
  process.patConversions = cms.EDProducer("PATConversionProducer",
      electronSource = cms.InputTag("selectedPatElectronsPFlow")  
      )


# The paths
if runStandardPAT:

  if useCaloJets:

    process.p = cms.Path()
    if useTrigger:
      process.p += process.step0a
    process.p += process.goodOfflinePrimaryVertices
    if useGoodVertex:
      process.p += process.step0b
    process.p += process.step0c
    process.p += process.eidMVASequence
    if useL1FastJet and useRelVals:
      process.p += process.ak5CaloJetSequence
    process.p += process.patDefaultSequence
    if usePFJets:
      process.p.remove( getattr( process, 'patJetCorrFactors'                    + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'jetTracksAssociatorAtVertex'          + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'impactParameterTagInfos'              + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'secondaryVertexTagInfos'              + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'softMuonTagInfos'                     + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'secondaryVertexNegativeTagInfos'      + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'jetBProbabilityBJetTags'              + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'jetProbabilityBJetTags'               + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'trackCountingHighPurBJetTags'         + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'trackCountingHighEffBJetTags'         + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'simpleSecondaryVertexHighEffBJetTags' + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'simpleSecondaryVertexHighPurBJetTags' + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'combinedSecondaryVertexBJetTags'      + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'combinedSecondaryVertexMVABJetTags'   + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'softMuonBJetTags'                     + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'softMuonByPtBJetTags'                 + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'softMuonByIP3dBJetTags'               + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'simpleSecondaryVertexNegativeHighEffBJetTags' + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'simpleSecondaryVertexNegativeHighPurBJetTags' + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'negativeTrackCountingHighEffJetTags'          + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'negativeTrackCountingHighPurJetTags'          + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patJetCharge'                         + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patJetPartonMatch'                    + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patJetGenJetMatch'                    + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patJetPartonAssociation'              + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patJetFlavourAssociation'             + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patJets'                              + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'patMETs'                              + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'selectedPatJets'                      + jetAlgo + pfSuffix ) )
      process.p.remove( getattr( process, 'countPatJets'                         + jetAlgo + pfSuffix ) )
    process.p += process.patAddOnSequence
    if useLooseMuon:
      process.p += process.step1a
    if useTightMuon:
      process.p += process.step1b
    if useMuonVeto:
      process.p += process.step2
    if useElectronVeto:
      process.p += process.step3
    if use1Jet:
      process.p += process.step4a
    if use2Jets:
      process.p += process.step4b
    if use3Jets:
      process.p += process.step4c
    if use4Jets:
      process.p += process.step5
    process.out.SelectEvents.SelectEvents.append( 'p' )

  if usePFJets:

    pAddPF = cms.Path()
    if useTrigger:
      pAddPF += process.step0a
    pAddPF += process.goodOfflinePrimaryVertices
    if useGoodVertex:
      pAddPF += process.step0b
    pAddPF += process.step0c
    pAddPF += process.eidMVASequence
    pAddPF += process.phoPFIso
#    pAddPF += process.additionalPatJets
    if useL1FastJet:
      pAddPF += process.ak5PFJets
    pAddPF += process.patDefaultSequence
    pAddPF.remove( process.patJetCorrFactors )
    if useCaloJets and useL1FastJet and useRelVals:
      pAddPF.remove( process.jetTracksAssociatorAtVertex )
      pAddPF.remove( process.impactParameterTagInfosAOD )
      pAddPF.remove( process.secondaryVertexTagInfosAOD )
      pAddPF.remove( process.softMuonTagInfosAOD )
      pAddPF.remove( process.secondaryVertexNegativeTagInfosAOD )
      pAddPF.remove( process.jetBProbabilityBJetTagsAOD )
      pAddPF.remove( process.jetProbabilityBJetTagsAOD )
      pAddPF.remove( process.trackCountingHighPurBJetTagsAOD )
      pAddPF.remove( process.trackCountingHighEffBJetTagsAOD )
      pAddPF.remove( process.simpleSecondaryVertexHighEffBJetTagsAOD )
      pAddPF.remove( process.simpleSecondaryVertexHighPurBJetTagsAOD )
      pAddPF.remove( process.combinedSecondaryVertexBJetTagsAOD )
      pAddPF.remove( process.combinedSecondaryVertexMVABJetTagsAOD )
      pAddPF.remove( process.softMuonBJetTagsAOD )
      pAddPF.remove( process.softMuonByPtBJetTagsAOD )
      pAddPF.remove( process.softMuonByIP3dBJetTagsAOD )
      pAddPF.remove( process.simpleSecondaryVertexNegativeHighEffBJetTagsAOD )
      pAddPF.remove( process.simpleSecondaryVertexNegativeHighPurBJetTagsAOD )
      pAddPF.remove( process.negativeTrackCountingHighEffJetTagsAOD )
      pAddPF.remove( process.negativeTrackCountingHighPurJetTagsAOD )                        
    pAddPF.remove( process.patJetCharge )
    pAddPF.remove( process.patJetPartonMatch )
    pAddPF.remove( process.patJetGenJetMatch )
    pAddPF.remove( process.patJetPartonAssociation )
    pAddPF.remove( process.patJetFlavourAssociation )
    pAddPF.remove( process.patJets )
    pAddPF.remove( process.patMETs )
    pAddPF.remove( process.selectedPatJets )
    pAddPF.remove( process.countPatJets )
    pAddPF += process.patAddOnSequence
    # pAddPF.replace( process.loosePatMuons
#                   , getattr( process, 'loosePatMuons' + jetAlgo + pfSuffix )
#                   )
#     pAddPF.replace( process.tightPatMuons
#                   , getattr( process, 'tightPatMuons' + jetAlgo + pfSuffix )
#                   )
#     pAddPF.replace( process.goodPatJets
#                   , getattr( process, 'goodPatJets' + jetAlgo + pfSuffix )
#                   )
    if useLooseMuon:
      pAddPF += getattr( process, 'step1a' + jetAlgo + pfSuffix )
    if useTightMuon:
      pAddPF += getattr( process, 'step1b' + jetAlgo + pfSuffix )
    if useMuonVeto:
      pAddPF += process.step2
    if useElectronVeto:
      pAddPF += process.step3
    if use1Jet:
      pAddPF += getattr( process, 'step4a' + jetAlgo + pfSuffix )
    if use2Jets:
      pAddPF += getattr( process, 'step4b' + jetAlgo + pfSuffix )
    if use3Jets:
      pAddPF += getattr( process, 'step4c' + jetAlgo + pfSuffix )
    if use4Jets:
      pAddPF += getattr( process, 'step5' + jetAlgo + pfSuffix )
    setattr( process, 'p' + jetAlgo + pfSuffix, pAddPF )
    process.out.SelectEvents.SelectEvents.append( 'p' + jetAlgo + pfSuffix )

if runPF2PAT:
      
  pPF = cms.Path()
  if useTrigger:
    pPF += process.step0a
  pPF += process.goodOfflinePrimaryVertices
  if useGoodVertex:
    pPF += process.step0b
  pPF += process.step0c
  pPF += process.eidMVASequence
  pPF += process.phoPFIso
  if runOnMC and isFatJetVal: 
    pPF += process.additionalGenJets 
  pPF += getattr( process, 'patPF2PATSequence' + postfix )
  pPF += process.patDefaultSequence
  pPF += process.pfInputsSequence
  pPF += process.metUncertaintySequence
  if runOnMC:
    pPF += process.selectedSmearedPatJetsPFlow
  else:
    process.puJetIdSequenceChs.remove( process.puSmearedJetIdChs )
    process.puJetIdSequenceChs.remove( process.puSmearedJetMvaChs )
  pPF += process.puJetIdSequenceChs
  pPF += process.goodLeptons
  pPF += process.largestGoodLeptonCand
  pPF += process.pfMEtMVAsequence
  pPF += getattr( process, 'patAddOnSequence' + postfix )
  pPF += process.patConversions
  if useLooseMuon:
    pPF += getattr( process, 'step1a' + postfix )
  if useTightMuon:
    pPF += getattr( process, 'step1b' + postfix )
  if useMuonVeto:
    pPF += getattr( process, 'step2' + postfix )
  if useElectronVeto:
    pPF += getattr( process, 'step3' + postfix )
  if use1Jet:
    pPF += getattr( process, 'step4a' + postfix )
  if use2Jets:
    pPF += getattr( process, 'step4b' + postfix )
  if use3Jets:
    pPF += getattr( process, 'step4c' + postfix )
  if use4Jets:
    pPF += getattr( process, 'step5' + postfix )


#PerformSelectionCrab part of code

#!
#! PARAMETERS
#!

  process.load('CommonTools.UtilAlgos.TFileService_cfi')
  process.TFileService.fileName=cms.string('PS.root')

  process.PS = cms.EDAnalyzer('PerformSelection')

    #-----Input Tags For Handles
    #
    # From PF2PAT, standard postfix=PFlow
    #
  process.PS.vtxSource         = cms.InputTag('offlinePrimaryVertices')               # InputTag for vertex collection
  process.PS.genParticleSource = cms.InputTag('genParticles')                         # InputTag for genParticle collection
  process.PS.pfJetSource       = cms.InputTag('selectedPatJetsPFlow')                 # InputTag for jet collection
  process.PS.electronSource    = cms.InputTag('selectedPatElectronsPFlow')            # InputTag for electron collection
  process.PS.muonSource        = cms.InputTag('selectedPatMuonsPFlow')                # InputTag for muon collection
  process.PS.METSource         = cms.InputTag('patMETsPFlow')                         # InputTag for MET collection
  process.PS.rhoSource         = cms.InputTag("kt6PFJetsPFlow", "rho")                # InputTag for event rho source
  process.PS.triggerSource     = cms.InputTag('patTriggerEvent')                      # InputTag for trigger collection
  process.PS.pileupSource      = cms.InputTag('addPileupInfo')                        # InputTag for pileupSummaryInfos collection
  
      #-----Trigger Information
  process.PS.muTrigger         = cms.vstring('HLT_IsoMu24_eta2p1_v*')                 # Muon trigger name
  #process.PS.muTrigger         = cms.vstring('HLT_Mu24_eta2p1_v*')                    # Muon trigger name
  #process.PS.muTrigger         = cms.vstring('HLT_IsoMu24_*',"HLT_IsoMu30_*")                 # Muon trigger name
  process.PS.eleTrigger        = cms.vstring('HLT_Ele27_WP80_v*')                     # Electron trigger name
  #process.PS.eleTrigger        = cms.vstring('HLT_Ele27_*','HLT_Ele32_*')                     # Electron trigger name 
  
      #-----Program Level Inputs
  process.PS.outtablefilename  = cms.string("outputTable.txt") # name of the file containing the output/efficiency tables
  process.PS.outtablenameEl    = cms.string("ElectronTable")   # name of the event selection table for electrons
  process.PS.outtablenameMu    = cms.string("MuonTable")       # name of the event selection table for muons
  process.PS.outtablenameLp    = cms.string("LeptonTable")     # name of the event selection table for leptons
  
  process.PS.printEventInfo    = cms.bool(False)               # tells the program to cout which cuts the events passed and what objects they contained
  process.PS.printJetInfo      = cms.bool(False)               # tells the program to cout the relevant information for the jets (for debugging purposes)
  process.PS.printLeptonInfo   = cms.bool(False)               # tells the program to cout the relevant information for the leptons (for debugging purposes)
  process.PS.SQWaT_Version     = cms.int32(14)                 # version number used to create the common PATtuples
  process.PS.doTrackerIso      = cms.bool(False)               # if true, use a tracker based isolation cut
  process.PS.doDetectorIso     = cms.bool(False)               # if true, use a rho corrected, detector based isolation cut
  process.PS.doPFIso           = cms.bool(True)                # if true, use a rho corrected, PF based isolation cut
  process.PS.doMVAeleSel       = cms.bool(True)                # if true, use the MVA based electron selection
  process.PS.noMVAIsoCut       = cms.bool(False)               # if true, turn off the MVA and PF isolation requirements to make a "Full" sample
  process.PS.doJER             = cms.bool(True)                # if true, JER corrections will be applied to all of the jets before passing selection
  process.PS.doMETPhi          = cms.bool(True)                # if true, the inherent METx and METy shifts will be corrected for using hard coded functions only needed before CMSSW_6_2_X
  process.PS.doJESUncertainty  = cms.bool(False)               # if true, turns on the jet uncertainty scaling
  process.PS.JESUncertainty    = cms.string("up")            # if none, the jet energy uncertainty is stored
                                                               # if up, the jet energy is scaled up by the uncertainty
                                                               # if down, the jet energy is scaled down by the uncertainty
  process.PS.JESUncertaintyType= cms.string("TotalNoTime")     # JES uncertainty collection
  process.PS.JESUncertaintyFile= cms.string("Winter14_V5_DATA_UncertaintySources_AK5PFchs.txt")
  
  process.PS.Data              = cms.bool(False)               # is the dataset from real data or Monte Carlo
  process.PS.MCpTrigger        = cms.bool(True)                # if true, only event that pass the trigger requirements will be saved
  process.PS.saveGenParticles  = cms.bool(True)                # save the generated particle information for hard scatter decays
  process.PS.particleStatus    = cms.int32(3)                  # The hard scatter particle status (eliminates the shower particles and protons). Pythia6 = 3, herwigpp = 2
  process.PS.saveMETPhiPlots   = cms.bool(False)               # save the TH1D and TH2D plots that have to do with MET Phi Corrections.
  process.PS.noMETCut          = cms.bool(True)               # disregard the MET cut when storing the events
  process.PS.invertEID         = cms.bool(False)               # electrons which *fail* at least two of the EID requirements will be kept instead
  process.PS.PFlowLoose        = cms.bool(False)               # use the collections with the PF2PAT postfix=PFlowLoose (instead of PFlow)
  process.PS.elONLY            = cms.bool(False)               # only save the output if the lepton is an electron
  process.PS.muONLY            = cms.bool(False)               # only save the output if the lepton is a muon
  process.PS.OnVsOffShell      = cms.bool(False)               # save all the trees in the 2 jets tree and don't use any cuts
  process.PS.StoreJets0        = cms.bool(False)               # save the jets0 tree
  process.PS.StoreJet1         = cms.bool(True)                # save the jet1 tree
  
      #-----Event Variable Inputs
  process.PS.elcnt_Prim                 = cms.int32(0)
  process.PS.elcnt_Loose                = cms.int32(0)
  process.PS.mucnt_Prim                 = cms.int32(0)
  process.PS.mucnt_Loose                = cms.int32(0)
  process.PS.leptonCnt                  = cms.int32(2)
  process.PS.jcnt_tot                   = cms.int32(0)
  process.PS.EvtTotCount                = cms.int32(0)
  process.PS.mu_passAll                 = cms.bool(False)
  process.PS.mu_passStandard            = cms.bool(False)
  process.PS.mu_passFlag                = cms.bool(False)
  process.PS.el_passAll                 = cms.bool(False)
  process.PS.el_passStandard            = cms.bool(False)
  process.PS.el_passFlag                = cms.bool(False)
  
      #-----Trigger Variable Inputs
  process.PS.mu_passTrigger             = cms.bool(True)
  process.PS.el_passTrigger             = cms.bool(True)
  
      #-----Vertex Variable Inputs
  process.PS.PVfound                    = cms.bool(False)
  process.PS.vtxcnt                     = cms.int32(0)
  
      #-----Jet Variable Inputs
  process.PS.j_ptMin                    = cms.double(25.0) #25.0 #raised from 20 to reduce size
  process.PS.j_aetaMax                  = cms.double(2.4)
  process.PS.j_DRelMin                  = cms.double(0.3)
  process.PS.muPrim_DRjMin              = cms.double(0.3)
  process.PS.CHEFMin                    = cms.double(0.0)
  process.PS.CEMEFMax                   = cms.double(0.99)
  process.PS.NHEFMax                    = cms.double(0.99)
  process.PS.NEMEFMax                   = cms.double(0.99)
  process.PS.NDaughtersMin              = cms.uint32(1)
  process.PS.CMultiplicityMin           = cms.int32(0)
  
      #-----B-Tag Variable Inputs
  process.PS.nBtagSSV                   = cms.int32(0)
  process.PS.nBtagTC                    = cms.int32(0)
  process.PS.bDiscriminatorSSVMin       = cms.double(1.74)     # https://twiki.cern.ch/twiki/bin/view/CMS/BTagPerformanceOP
  process.PS.bDiscriminatorTCMin        = cms.double(3.3)      # https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagging?redirectedfrom=CMS.SWGuideBTagging#Algorithms
  process.PS.bDiscriminatorCSVMin      = cms.double(0.679)     # corresponds to the 'medium' cuts
  process.PS.bDiscriminatorCSVMVAMin   = cms.double(0.679)
  
      #-----Muon Variable Inpust
  process.PS.muPrim_ptMin               = cms.double(24.0) # raised from 20 to reduce size
  process.PS.muPrim_aetaMax             = cms.double(2.1)
  process.PS.muPrim_dBMax               = cms.double(0.02) #Changed for FNAL numbers
  process.PS.muPrim_dzMax               = cms.double(0.5)
  process.PS.muPrim_TrkIsoMax           = cms.double(0.05)
  process.PS.muPrim_DetIsoMax           = cms.double(0.10)
  process.PS.muPrim_PFIsoMax            = cms.double(0.12)
  process.PS.muPrim_normChi2Max         = cms.double(10.0)
  process.PS.muPrim_nValidHitsMin       = cms.int32(5)
  process.PS.muPrim_nValidMuonHitsMin   = cms.int32(0)
  process.PS.muPrim_nValidPixelHitsMin  = cms.int32(0)
  process.PS.muPrim_nMatchedStationsMin = cms.int32(1)
  process.PS.isProperlyIsolatedMu       = cms.bool(False)
  process.PS.muDetIsoConeSize           = cms.double(0.3)
  process.PS.muPFIsoConeSize            = cms.double(0.4)
  process.PS.muLoose_ptMin              = cms.double(10.0)
  process.PS.muLoose_aetaMax            = cms.double(2.5)
  process.PS.muLoose_TrkIsoMax          = cms.double(0.1)
  process.PS.muLoose_DetIsoMax          = cms.double(0.2)
  process.PS.muLoose_PFIsoMax           = cms.double(0.2)
  
  
      #-----Electron Variable Inputs
  process.PS.elPrim_ptMin               = cms.double(27.0) # raised from 20 to reduce size
  process.PS.elPrim_aetaMax             = cms.double(2.5)
  process.PS.elPrim_aetascExcludeMax    = cms.double(1.5660)
  process.PS.elPrim_aetascExcludeMin    = cms.double(1.4442)
  process.PS.elPrim_sigmaIetaIetaMaxEB  = cms.double(0.01)
  process.PS.elPrim_aDeltaPhiMaxEB      = cms.double(0.03)
  process.PS.elPrim_aDeltaEtaMaxEB      = cms.double(0.004)
  process.PS.elPrim_HoEMaxEB            = cms.double(0.12)
  process.PS.elPrim_sigmaIetaIetaMaxEE  = cms.double(0.03)
  process.PS.elPrim_aDeltaPhiMaxEE      = cms.double(0.02)
  process.PS.elPrim_aDeltaEtaMaxEE      = cms.double(0.005)
  process.PS.elPrim_HoEMaxEE            = cms.double(0.10)
  process.PS.elPrim_dBMax               = cms.double(0.02)
  process.PS.elPrim_dzMax               = cms.double(0.1)
  process.PS.elPrim_EinvMinusPinvMax    = cms.double(0.05)
  process.PS.elPrim_nMissingHits        = cms.int32(0)
  process.PS.elPrim_TrkIsoMax           = cms.double(0.1)
  process.PS.elPrim_DetIsoMax           = cms.double(0.05)
  process.PS.elPrim_PFIsoMax            = cms.double(0.1)
  process.PS.elPrim_PFIsoMin_invertEID  = cms.double(0.3)
  
  process.PS.elLoose_ptMin              = cms.double(15.0)
  process.PS.elLoose_sigmaIetaIetaMaxEB = cms.double(0.01)
  process.PS.elLoose_aDeltaPhiMaxEB     = cms.double(0.8)
  process.PS.elLoose_aDeltaEtaMaxEB     = cms.double(0.007)
  process.PS.elLoose_HoEMaxEB           = cms.double(0.15)
  process.PS.elLoose_sigmaIetaIetaMaxEE = cms.double(0.03)
  process.PS.elLoose_aDeltaPhiMaxEE     = cms.double(0.7)
  process.PS.elLoose_aDeltaEtaMaxEE     = cms.double(0.01)
  process.PS.elLoose_HoEMaxEE           = cms.double(0.07)
  process.PS.elLoose_dBMax              = cms.double(0.04)
  process.PS.elLoose_dzMax              = cms.double(0.2)
  process.PS.elLoose_TrkIsoMax          = cms.double(0.2)
  process.PS.elLoose_DetIsoMax          = cms.double(0.2)
  process.PS.elLoose_PFIsoMax           = cms.double(0.15)
  
  process.PS.el_aetaPass                = cms.bool(False)
  
  process.PS.elDetIsoConeSize           = cms.double(0.3)
  process.PS.elPFIsoConeSize            = cms.double(0.3)
  
      #-----Lepton Variable Inputs
  
      #-----MET Variable Inputs
  process.PS.MET_PtMin                  = cms.double(25)
  process.PS.MET_PtMin_invertEID        = cms.double(20)
  
      #-----Additional Variable Inputs
  process.PS.lTotIso                    = cms.double(-1)
  process.PS.lecalIso                   = cms.double(-10)
  process.PS.lhcalIso                   = cms.double(-10)
  process.PS.ltrkIso                    = cms.double(-1)
  process.PS.ldetIso                    = cms.double(-1)
  process.PS.lpfIso                     = cms.double(-1)
  process.PS.Mjj                        = cms.double(-1)
  process.PS.lQ                         = cms.int32(-10)
  process.PS.lEta                       = cms.double(-1)
  
      #-----Constant Inputs
  process.PS.etaBarrelMax               = cms.double(1.442) # the Maximum value of eta for the Barrel
  process.PS.NPtBins                    = cms.int32(13)
  process.PS.vpt                        = cms.vdouble(25,35,45,55,70,85,100,150,200,400,600,800,1000,1500)
  
  #!
  #! ADD TO PATH
  #!
  pPF += process.PS

  del process.outpath

  setattr( process, 'p' + postfix, pPF )
  process.out.SelectEvents.SelectEvents.append( 'p' + postfix )
