from run_PerformSelection_DataPreTriggerTemplate import *

#-----Trigger Information
process.PS.muTrigger                  = cms.vstring('HLT_IsoMu17_v*','HLT_Mu30_v*')          # Muon trigger name(s)
process.PS.eleTrigger                 = cms.vstring('HLT_Ele*')          # Electron trigger name(s)


#!
#! PATH
#!
process.p = cms.Path(process.PS)

#process.options.wantSummary = True
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
