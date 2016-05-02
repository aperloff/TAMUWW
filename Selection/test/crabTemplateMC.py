from WMCore.Configuration import Configuration
config = Configuration()
config.section_("General")
config.General.requestName = 'REQUESTNAME'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'run_patMC_PerformSelection_cfg.py'
config.JobType.allowNonProductionCMSSW = True
config.JobType.inputFiles = ['Winter14_V5_DATA_UncertaintySources_AK5PFchs.txt']

config.section_("Data")
config.Data.inputDataset = 'INPUTDATASET'
config.Data.inputDBS = 'https://cmsweb.cern.ch/dbs/prod/DBSURL/DBSReader/'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = FILESPERJOB
config.Data.publication = False
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter/'
config.Data.publishDataName = 'PUBLISHDATANAME'
config.Data.outLFN = '/store/user/aperloff/MatrixElement/'
config.Data.ignoreLocality = True

config.section_("Site")
config.Site.storageSite = 'T3_US_FNALLPC'
