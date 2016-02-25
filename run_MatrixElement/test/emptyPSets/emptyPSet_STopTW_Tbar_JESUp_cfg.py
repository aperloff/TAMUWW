import FWCore.ParameterSet.Config as cms
import os

#!
#! PROCESS
#!
process = cms.Process("MatrixElementProcess")


#!
#! SERVICES
#!
#process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.load('CommonTools.UtilAlgos.TFileService_cfi')
process.TFileService.fileName=cms.string('STopTW_Tbar_JESUp.root')


#!
#! INPUT
#!
inputFiles = cms.untracked.vstring(
    'root://cmsxrootd.fnal.gov//store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/STopTW_Tbar_JESUp.root'
    )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.source = cms.Source("PoolSource",
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = inputFiles )

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
