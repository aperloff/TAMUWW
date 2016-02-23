import FWCore.ParameterSet.Config as cms

#!
#! DEFAULT PARAMETERS (CHANGE WHEN NEEDED)
#!
#import TAMUWW.Selection.Defaults_cff as Defaults;


#!
#! PROCESS
#!
process = cms.Process("SaveGenParticlesProcess")


#!
#! SERVICES
#!
#process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('SGP_WJets.root')


#!
#! INPUT
#!
#inputFiles = cms.untracked.vstring(
######################
# ggH125 Monte Carlo #
######################
#    '/store/user/aperloff/LQ-ggh125_SIM/LQ-ggh125_AODSIM/51bfef3f471c3bed6331eb644e144b83/gghlvjj-125_AODSIM_24_1_aNW.root'
#)
######################
# DYJets Monte Carlo #
######################
#process.load("TAMUWW.Selection.DYJetsToLL_M50_8TeV_53X_pat_cff")
#####################
# WJets Monte Carlo #
#####################
process.load("TAMUWW.Selection.WJetsToLNu_8TeV_53X_pat_cff")


process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000000))
#process.source = cms.Source("PoolSource",
#                            skipEvents = cms.untracked.uint32(0),
#                            fileNames = inputFiles )


#!
#! PARAMETERS
#!
process.PS = cms.EDAnalyzer('SaveGenParticles')

process.PS.genParticleSource = cms.InputTag('genParticles') # InputTag for genParticle collection

#!
#! PATH
#!
process.p = cms.Path(process.PS)

#process.options.wantSummary = True
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
