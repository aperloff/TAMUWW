import FWCore.ParameterSet.Config as cms

#!
#! DEFAULT PARAMETERS (CHANGE WHEN NEEDED)
#!
#import TAMUWW.Selection.Defaults_cff as Defaults;


#!
#! PROCESS
#!
process = cms.Process("PerformSelectionProcess")


#!
#! SERVICES
#!
#process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('PS.root')


#!
#! INPUT
#!
inputFiles = cms.untracked.vstring(
##################
# WZ Monte Carlo #
##################
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_9_1_Lph.root',
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_99_1_pzW.root',
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_98_1_77z.root',
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_97_1_GnJ.root',
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_96_1_wKp.root',
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_95_1_2Ke.root',
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_94_1_Aj2.root',
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_93_1_LYz.root'
#####################
# WJets Monte Carlo #
#####################
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_518_1_QXj.root',
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_517_1_N8m.root',
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_510_1_87t.root'
######################
# ggH125 Monte Carlo #
######################
#    '/store/user/aperloff/LQ-ggh125_SIM/SQWaT_PAT_52X_Summer12_v1/290326670ba15ca0752d90668da7d2ec/pat_52x_test_1_1_5x8.root'
#####################
# WH120 Monte Carlo #
#####################
#    'dcap::///pnfs/cms/WAX/11/store/user/ssengupt/H_M120/output/ssengupt/WH_ZH_TTH_HToWW_M-120_8TeV-pythia6/SQWaT_PAT_52X_wH_zH_ttH_120_Summer12_v2/290326670ba15ca0752d90668da7d2ec/pat_52x_test_66_1_Abv.root'
#############
# Muon Data #
#############
#    '/store/user/lnujj/PATtuples/sboutle/SingleMu/SQWaT_PAT_52X_2012A-PromptReco-v1_v5//79dd4263edc1b886eee0f23fb76cda69/pat_52x_test_9_1_jB0.root',
#    '/store/user/lnujj/PATtuples/sboutle/SingleMu/SQWaT_PAT_52X_2012A-PromptReco-v1_v5//79dd4263edc1b886eee0f23fb76cda69/pat_52x_test_9_1_RdL.root',
#    '/store/user/lnujj/PATtuples/sboutle/SingleMu/SQWaT_PAT_52X_2012A-PromptReco-v1_v5//79dd4263edc1b886eee0f23fb76cda69/pat_52x_test_9_1_RD8.root'
#################
# Electron Data #
#################
#    '/store/user/lnujj/PATtuples/sboutle/SingleElectron/SQWaT_PAT_52X_2012A-PromptReco-v1_v5//79dd4263edc1b886eee0f23fb76cda69/pat_52x_test_9_1_qLJ.root',
#    '/store/user/lnujj/PATtuples/sboutle/SingleElectron/SQWaT_PAT_52X_2012A-PromptReco-v1_v5//79dd4263edc1b886eee0f23fb76cda69/pat_52x_test_9_1_fwF.root'
    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = inputFiles )


#!
#! PARAMETERS
#!
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
process.PS.eleTrigger        = cms.vstring('HLT_Ele27_WP80_v*')                     # Electron trigger name 

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

process.PS.Data              = cms.bool(False)               # is the dataset from real data or Monte Carlo
process.PS.MCpTrigger        = cms.bool(True)                # if true, only event that pass the trigger requirements will be saved
process.PS.saveGenParticles  = cms.bool(True)                # save the generated particle information for hard scatter decays
process.PS.noMETCut          = cms.bool(False)               # disregard the MET cut when storing the events
process.PS.invertEID         = cms.bool(False)               # electrons which *fail* at least two of the EID requirements will be kept instead
process.PS.PFlowLoose        = cms.bool(False)               # use the collections with the PF2PAT postfix=PFlowLoose (instead of PFlow)
process.PS.elONLY            = cms.bool(False)               # only save the output if the lepton is an electron
process.PS.muONLY            = cms.bool(False)               # only save the output if the lepton is a muon
process.PS.OnVsOffShell      = cms.bool(False)               # save all the trees in the 2 jets tree and don't use any cuts
process.PS.StoreJets01       = cms.bool(False)               # save the jets0 and jet1 trees

    #-----Event Variable Inputs
process.PS.elcnt_Prim                 = cms.int32(0)
process.PS.elcnt_Loose                = cms.int32(0)
process.PS.mucnt_Prim                 = cms.int32(0)
process.PS.mucnt_Loose                = cms.int32(0)
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
process.PS.elPrim_PFIsoMin_invertEID  = cms.double(0.15)

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
process.PS.MET_EtMin                  = cms.double(25)
process.PS.MET_EtMin_invertEID        = cms.double(20)

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

#!
#! PATH
#!
process.p = cms.Path(process.PS)

#process.options.wantSummary = True
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
