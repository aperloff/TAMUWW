from os import system

#jobs=[#'REQUESTNAME','INPUTDATASET','PUBLISHDATANAME','LBSPERJOB','DBSURL'],
jobs=[#'REQUESTNAME','INPUTDATASET','PUBLISHDATANAME','FILESPERJOB','DBSURL'],
#	['WJets_part1_JESUp','/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM','WJets_part1_JESUp',1,'global'],
#	['WJets_part2_JESUp','/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball/Summer12_DR53X-PU_S10_START53_V7A-v2/AODSIM','WJets_part2_JESUp',1,'global'],
#	['TTbar_JESUp','/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM','TTbar_JESUp',1,'global'],
#	['qqH125_JESUp','/LQ-vbf125_GENSIM/ajkumar-LQ-qqh125_AODSIM_Summer12_START53_V7A-c8f8ed334db8a7d6f56c62266b1dfa5b/USER','qqH125_JESUp',4,'phys02'],
#        ['SingleMu_A_ReReco',  '/SingleMu/Run2012A-recover-06Aug2012-v1/AOD',      'SingleMu_A_ReReco',  5000,'global','JSON/SingleMuon/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt'],
#        ['SingleMu_A',         '/SingleMu/Run2012A-13Jul2012-v1/AOD',              'SingleMu_A',         5000,'global','JSON/SingleMuon/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt'],
#        ['SingleMu_B_p1',      '/SingleMu/Run2012B-13Jul2012-v1/AOD',              'SingleMu_B_p1',      7500,'global','JSON/SingleMuon/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt'],
#        ['SingleMu_B_p2',      '/SingleMu/Run2012B-13Jul2012-v1/AOD',              'SingleMu_B_p2',      7500,'global','JSON/SingleMuon/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt'],
#        ['SingleMu_C_v1',      '/SingleMu/Run2012C-24Aug2012-v1/AOD',              'SingleMu_C_v1',      3,'global','JSON/SingleMuon/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt'],
        ['SingleMu_C_v2_p1_p1','/SingleMu/Run2012C-PromptReco-v2/AOD',             'SingleMu_C_v2_p1_p1',4,'global','JSON/SingleMuon/Cert_190456-200601_8TeV_PromptReco_Collisions12_JSON_v2.txt'],
        ['SingleMu_C_v2_p1_p2','/SingleMu/Run2012C-PromptReco-v2/AOD',             'SingleMu_C_v2_p1_p2',4,'global','JSON/SingleMuon/Cert_190456-200601_8TeV_PromptReco_Collisions12_JSON_v2.txt'],
        ['SingleMu_C_v2_p2',   '/SingleMu/Run2012C-PromptReco-v2/AOD',             'SingleMu_C_v2_p2',   3,'global','JSON/SingleMuon/Cert_200519-202016_8TeV_PromptReco_Collisions12_JSON.txt'],
        ['SingleMu_C_v2_p3',   '/SingleMu/Run2012C-PromptReco-v2/AOD',             'SingleMu_C_v2_p3',   3,'global','JSON/SingleMuon/Cert_194631-203002_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleMu_D_v1_p1',   '/SingleMu/Run2012D-PromptReco-v1/AOD',             'SingleMu_D_v1_p1',   3,'global','JSON/SingleMuon/Cert_203894-205618_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleMu_D_v1_p2',   '/SingleMu/Run2012D-PromptReco-v1/AOD',             'SingleMu_D_v1_p2',   3,'global','JSON/SingleMuon/Cert_205620-206539_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleMu_D_v1_p3',   '/SingleMu/Run2012D-PromptReco-v1/AOD',             'SingleMu_D_v1_p3',   2,'global','JSON/SingleMuon/Cert_205339-206088_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleMu_D_v1_p4',   '/SingleMu/Run2012D-PromptReco-v1/AOD',             'SingleMu_D_v1_p4',   3,'global','JSON/SingleMuon/Cert_200961-206940_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleMu_D_v1_p5',   '/SingleMu/Run2012D-PromptReco-v1/AOD',             'SingleMu_D_v1_p5',   2,'global','JSON/SingleMuon/Cert_206744-207469_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleMu_D_v1_p6',   '/SingleMu/Run2012D-PromptReco-v1/AOD',             'SingleMu_D_v1_p6',   2,'global','JSON/SingleMuon/Cert_mu_207477-207898_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleMu_D_v1_p7',   '/SingleMu/Run2012D-PromptReco-v1/AOD',             'SingleMu_D_v1_p7',   3,'global','JSON/SingleMuon/Cert_207889-208357_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleMu_D_v1_p8',   '/SingleMu/Run2012D-PromptReco-v1/AOD',             'SingleMu_D_v1_p8',   1,'global','JSON/SingleMuon/Cert_194480-208686_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleEl_A_ReReco',  '/SingleElectron/Run2012A-recover-06Aug2012-v1/AOD','SingleEl_A_ReReco',  5,'global','JSON/SingleElectron/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt'],
#        ['SingleEl_A',         '/SingleElectron/Run2012A-13Jul2012-v1/AOD',        'SingleEl_A',         5000,'global','JSON/SingleElectron/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt'],
#        ['SingleEl_B_p1',      '/SingleElectron/Run2012B-13Jul2012-v1/AOD',        'SingleEl_B_p1',      8000,'global','JSON/SingleElectron/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt'],
#        ['SingleEl_B_p2',      '/SingleElectron/Run2012B-13Jul2012-v1/AOD',        'SingleEl_B_p2',      8500,'global','JSON/SingleElectron/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON.txt'],
#        ['SingleEl_C_v1',      '/SingleElectron/Run2012C-24Aug2012-v1/AOD',        'SingleEl_C_v1',      5000,'global','JSON/SingleElectron/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt'],
#        ['SingleEl_C_v2_p1',   '/SingleElectron/Run2012C-PromptReco-v2/AOD',       'SingleEl_C_v2_p1',   7000,'global','JSON/SingleElectron/Cert_198934-200532_8TeV_PromptReco_Collisions12_JSON_v2.txt'],
#        ['SingleEl_C_v2_p2',   '/SingleElectron/Run2012C-PromptReco-v2/AOD',       'SingleEl_C_v2_p2',   8500,'global','JSON/SingleElectron/Cert_200519-202016_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleEl_C_v2_p3',   '/SingleElectron/Run2012C-PromptReco-v2/AOD',       'SingleEl_C_v2_p3',   8500,'global','JSON/SingleElectron/Cert_194631-203002_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleEl_D_v1_p1',   '/SingleElectron/Run2012D-PromptReco-v1/AOD',       'SingleEl_D_v1_p1',   5000,'global','JSON/SingleElectron/Cert_203894-205618_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleEl_D_v1_p2',   '/SingleElectron/Run2012D-PromptReco-v1/AOD',       'SingleEl_D_v1_p2',   5000,'global','JSON/SingleElectron/Cert_205620-206539_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleEl_D_v1_p3',   '/SingleElectron/Run2012D-PromptReco-v1/AOD',       'SingleEl_D_v1_p3',   5000,'global','JSON/SingleElectron/Cert_205339-206088_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleEl_D_v1_p4',   '/SingleElectron/Run2012D-PromptReco-v1/AOD',       'SingleEl_D_v1_p4',   2300,'global','JSON/SingleElectron/Cert_200961-206940_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleEl_D_v1_p5',   '/SingleElectron/Run2012D-PromptReco-v1/AOD',       'SingleEl_D_v1_p5',   5000,'global','JSON/SingleElectron/Cert_206744-207469_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleEl_D_v1_p6',   '/SingleElectron/Run2012D-PromptReco-v1/AOD',       'SingleEl_D_v1_p6',   5000,'global','JSON/SingleElectron/Cert_207477-207898_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleEl_D_v1_p7',   '/SingleElectron/Run2012D-PromptReco-v1/AOD',       'SingleEl_D_v1_p7',   5000,'global','JSON/SingleElectron/Cert_207889-208357_8TeV_PromptReco_Collisions12_JSON.txt'],
#        ['SingleEl_D_v1_p8',   '/SingleElectron/Run2012D-PromptReco-v1/AOD',       'SingleEl_D_v1_p8',   1300,'global','JSON/SingleElectron/Cert_194480-208686_8TeV_PromptReco_Collisions12_JSON.txt'],
    ]

MC = False
if MC:
    template='crabTemplateMC.py'
else:
    template='crabTemplateData.py'

#---------------------------------------------------------------------------------------------------------------------------------------------------

for job in jobs:
    requestName=job[0]
    inputDS=job[1]
    publishDN=job[2]
    Units=str(job[3])
    dbsUrl=job[4]
    Mask=job[5]
    
    config=requestName+'.py'
    system('cp '+template+' '+config)
    system('sed \"s%REQUESTNAME%'+requestName+'%g\" '+config+' --in-place')
    system('sed \"s%INPUTDATASET%'+inputDS+'%g\" '+config+' --in-place')
    system('sed \"s%PUBLISHDATANAME%'+publishDN+'%g\" '+config+' --in-place')
    if MC:
        system('sed \"s%FILESPERJOB%'+Units+'%g\" '+config+' --in-place')
    else:
        system('sed \"s%UNITSPERJOB%'+Units+'%g\" '+config+' --in-place')
        system('sed \"s%LUMIMASK%'+Mask+'%g\" '+config+' --in-place')
    system('sed \"s%DBSURL%'+dbsUrl+'%g\" '+config+' --in-place')

    #system('crab submit '+config)
