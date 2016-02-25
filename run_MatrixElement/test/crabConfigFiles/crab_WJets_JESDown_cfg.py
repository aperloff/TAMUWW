from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'run_MatrixElement_WJets_JESDown'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'PrivateMC'
config.JobType.scriptExe = 'JobScript.sh'
config.JobType.scriptArgs = ["numberOfEvents=100","SampleName=WJets_JESDown"]
config.JobType.psetName = 'emptyPSet_WJets_JESDown_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.inputFiles = ['FrameworkJobReport.xml','../../../../bin/slc6_amd64_gcc472/run_MatrixElement','../data/cteq5l.tbl', '../data/cteq6l.tbl', '../data/TF_TTbarMG_B_00_24.txt', '../data/TF_TTbarMG_G_00_24.txt', '../data/TF_TTbarMG_UDS_00_24.txt']
config.JobType.outputFiles = ['WJets_JESDown.root']

config.section_("Data")
config.Data.userInputFiles = ['root://cmsxrootd.fnal.gov//store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/WJets_JESDown.root']
config.Data.primaryDataset = 'WJets_JESDown'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 1
NJOBS = 6116  # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEResults/rootOutput/'
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
