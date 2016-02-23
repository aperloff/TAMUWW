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
######################
# DYJets Monte Carlo #
######################
'/store/user/lnujj/PatTuples_8TeV_53X/custodio/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/SQWaT_PAT_53X_Summer12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_100_1_qoN.root'
##################
# WZ Monte Carlo #
##################
#'/store/user/lnujj/PatTuples_8TeV_53X/jdg6y/goodell/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_43X_Fall12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_94_1_jEA.root',
#'/store/user/lnujj/PatTuples_8TeV_53X/jdg6y/goodell/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_43X_Fall12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_95_1_4be.root',
#'/store/user/lnujj/PatTuples_8TeV_53X/jdg6y/goodell/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_43X_Fall12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_96_1_hKm.root',
#'/store/user/lnujj/PatTuples_8TeV_53X/jdg6y/goodell/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_43X_Fall12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_97_1_KXp.root',
#'/store/user/lnujj/PatTuples_8TeV_53X/jdg6y/goodell/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_43X_Fall12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_98_1_0Yo.root',
#'/store/user/lnujj/PatTuples_8TeV_53X/jdg6y/goodell/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_43X_Fall12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_99_1_gd7.root',
#'/store/user/lnujj/PatTuples_8TeV_53X/jdg6y/goodell/WZ_TuneZ2star_8TeV_pythia6_tauola/SQWaT_PAT_43X_Fall12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_9_1_cBY.root'
#####################
# WJets Monte Carlo #
#####################
#	'/store/user/lnujj/PatTuples_8TeV_53X-v1/jdamgov/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/SQWaT_PAT_53X_Summer12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_217_0_c0L.root',
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_518_1_QXj.root',
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_517_1_N8m.root',
#    '/store/user/lnujj/PatTuples_8TeV_v2/ajkumar/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/SQWaT_PAT_52X_Summer12_v2//290326670ba15ca0752d90668da7d2ec/pat_52x_test_510_1_87t.root'
######################
# ggH125 Monte Carlo #
######################
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_10_1_3lN.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_11_1_luM.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_12_1_N49.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_13_1_uI3.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_14_1_7qv.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_15_1_Gay.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_16_1_uPL.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_17_1_cT6.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_18_1_Y7m.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_19_1_HpJ.root',
#	   '/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_1_1_wvm.root',
#	   '/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_20_1_7FF.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_21_1_2uu.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_22_1_74F.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_23_1_0lW.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_24_1_DCU.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_25_1_gcN.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_26_1_6iZ.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_27_1_Wdb.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_28_1_obe.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_29_1_mt9.root',
#	   '/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_2_1_K4H.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_30_1_l1g.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_31_1_1NF.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_32_1_1p9.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_33_1_V62.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_34_1_L1Z.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_35_1_16h.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_36_1_OrL.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_37_1_RnK.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_3_1_Yh1.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_4_1_ZP4.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_5_1_EU6.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_6_1_KPS.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_7_1_RDd.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_8_1_rMd.root',
#	'/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_9_1_KPu.root'
#####################
# WH125 Monte Carlo #
#####################
#	'/store/user/lnujj/PatTuples_8TeV_53X/jdamgov/WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp/SQWaT_PAT_53X_Summer12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_445_1_x22.root',
#	'/store/user/lnujj/PatTuples_8TeV_53X/jdamgov/WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp/SQWaT_PAT_53X_Summer12_v1/829f288d768dd564418efaaf3a8ab9aa/pat_53x_test_v03_828_1_5Rd.root'
#########################
# WH_ZH_TTH Monte Carlo #
#########################
#	'/store/mc/Summer12_DR53X/WH_ZH_TTH_HToWW_M-125_8TeV-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/0CE9D2F9-3CF5-E111-A26E-00215E21DC78.root',
#	'/store/mc/Summer12_DR53X/WH_ZH_TTH_HToWW_M-125_8TeV-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/186F9B2D-45F5-E111-81B3-E41F1318148C.root'
#############
# Muon Data #
#############
#	'/store/user/lnujj/PatTuples_8TeV_53X/dimatteo/SingleMu/SQWaT_PAT_53X_2012C-24Aug2012-v1/3e4086321697e2c39c90dad08848274b/pat_53x_test_v03_data_90_1_elI.root',
#	'/store/user/lnujj/PatTuples_8TeV_53X/dimatteo/SingleMu/SQWaT_PAT_53X_2012C-24Aug2012-v1/3e4086321697e2c39c90dad08848274b/pat_53x_test_v03_data_91_1_OmO.root',
#	'/store/user/lnujj/PatTuples_8TeV_53X/dimatteo/SingleMu/SQWaT_PAT_53X_2012C-24Aug2012-v1/3e4086321697e2c39c90dad08848274b/pat_53x_test_v03_data_97_1_WN4.root'
#################
# Electron Data #
#################
#	'/store/user/eusebi/PatTuples_8TeV_53X/goodell/SingleElectron/SQWaT_PAT_53X_Fall12_Electron2012D_p7_loose/40abd2242537ef01a4e068cea4bb511b/pat_53x_test_v03_9_1_Odn.root',
#	'/store/user/eusebi/PatTuples_8TeV_53X/goodell/SingleElectron/SQWaT_PAT_53X_Fall12_Electron2012D_p7_loose/40abd2242537ef01a4e068cea4bb511b/pat_53x_test_v03_86_1_9Fu.root',
#	'/store/user/eusebi/PatTuples_8TeV_53X/goodell/SingleElectron/SQWaT_PAT_53X_Fall12_Electron2012D_p7_loose/40abd2242537ef01a4e068cea4bb511b/pat_53x_test_v03_88_1_D3l.root',
#	'/store/user/eusebi/PatTuples_8TeV_53X/goodell/SingleElectron/SQWaT_PAT_53X_Fall12_Electron2012D_p7_loose/40abd2242537ef01a4e068cea4bb511b/pat_53x_test_v03_94_1_3Y2.root',
#	'dcap:///pnfs/cms/WAX/11/store/user/lnujj/PatTuples_8TeV_53X/jdg6y/goodell/SingleElectron/SQWaT_PAT_53X_Fall12_v1_Electron2012B_Test/40abd2242537ef01a4e068cea4bb511b/pat_53x_test_v03_4_1_eDN.root',
#	'dcap:///pnfs/cms/WAX/11/store/user/lnujj/PatTuples_8TeV_53X/jdg6y/goodell/SingleElectron/SQWaT_PAT_53X_Fall12_v1_Electron2012B_Test/40abd2242537ef01a4e068cea4bb511b/pat_53x_test_v03_1_1_nnt.root'
#	'/store/user/lnujj/PatTuples_8TeV_53X/dimatteo/SingleElectron/SQWaT_PAT_53X_2012B-13Jul2012-v1_part1/3e4086321697e2c39c90dad08848274b/pat_53x_test_v03_data_107_1_FAG.root',
#	'/store/user/lnujj/PatTuples_8TeV_53X/dimatteo/SingleElectron/SQWaT_PAT_53X_2012B-13Jul2012-v1_part1/3e4086321697e2c39c90dad08848274b/pat_53x_test_v03_data_99_1_Xl0.root',
#	'/store/user/lnujj/PatTuples_8TeV_53X/dimatteo/SingleElectron/SQWaT_PAT_53X_2012B-13Jul2012-v1_part1/3e4086321697e2c39c90dad08848274b/pat_53x_test_v03_data_98_1_bH1.root',
#	'/store/user/lnujj/PatTuples_8TeV_53X/dimatteo/SingleElectron/SQWaT_PAT_53X_2012B-13Jul2012-v1_part1/3e4086321697e2c39c90dad08848274b/pat_53x_test_v03_data_97_1_2js.root'
#############
# Joey Skim #
#############
#	'/store/user/jdg6y/HWW125_Part8_JECup_muskim/goodell/LQ-ggh125_Part8_GENSIM/SHYFT_pat532_HWW125_Part8_JECup_muskim/9fd1d441a8eb7f2c4f1541582aa0dd17/ljmet_skim_100_1_58H.root'
#	'file:/eos/uscms/store/user/aperloff/MatrixElement/LQ-ggh125_BIG_PAT_part8/aperloff/LQ-ggh125_Part8_GENSIM/SQWaT_PAT_53X_Spring14_HWW125_Part8/33519632dfb71bf616e1fd5cd9d10cef/pat_53x_test_v03_11_1_luM.root'
    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5000))
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
process.PS.JESUncertainty    = cms.string("up")              # if none, the jet energy uncertainty is stored
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
process.PS.NPtBins                       = cms.int32(13)
process.PS.vpt                           = cms.vdouble(25,35,45,55,70,85,100,150,200,400,600,800,1000,1500)

#!
#! PATH
#!
process.p = cms.Path(process.PS)

#process.options.wantSummary = True
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
