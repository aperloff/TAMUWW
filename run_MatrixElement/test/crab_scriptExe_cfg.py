from WMCore.Configuration import Configuration
config = Configuration()
#from CRABClient.UserUtilities import config
#config = config()

config.section_("General")
config.General.requestName = 'run_MatrixElement_qqH125_v2'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
#config.JobType.pluginName = 'Analysis'
config.JobType.pluginName = 'PrivateMC'
config.JobType.scriptExe = 'CondorScript_WW.sh'
config.JobType.scriptArgs = ["numberOfEvents=3"]
config.JobType.psetName = 'emptyPSet_cfg.py'
config.JobType.allowNonProductionCMSSW = True
config.JobType.inputFiles = ['FrameworkJobReport.xml','../../../../bin/slc6_amd64_gcc472/run_MatrixElement','../data/cteq5l.tbl', '../data/cteq6l.tbl', '../data/TF_TTbarMG_B_00_24.txt', '../data/TF_TTbarMG_G_00_24.txt', '../data/TF_TTbarMG_UDS_00_24.txt']
config.JobType.outputFiles = ['WW.root']
#config.JobType.pyCfgParams = ['jobNumber','inputFilename=root://cmsxrootd.fnal.gov//store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/ggH125_BIG.root','totalNumberOfEvents=2','eventsPerJob=2']

config.section_("Data")
config.Data.userInputFiles = ['root://cmsxrootd.fnal.gov//store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/WW.root']
config.Data.primaryDataset = ['WW_METestJob']
#config.Data.inputDataset = '/QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6/Fall14DR73-PFr1_NoPU_RAW_GSFromFall14_MCRUN2_73_V11-v1/GEN-SIM-RECO'
#config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 1
NJOBS = 2  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = False
#config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
#config.Data.publishDataName = 'test_run_MatrixElement'
config.Data.outLFN = '/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEResults/rootOutput/test/'
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
