from os import system
from ROOT import *
from math import ceil

#basepath = "/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/"
#basepath = "root://cmseos.fnal.gov//store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/WlnuJetsTest/"
basepath = "root://cmseos.fnal.gov//store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/"

processes = [
#"qqH125", #submitted 05/08/2015 #needed completion job (qqH125_0.root) because skipped first 100 events #finished on 
#"qqH125_JESUp", #submitted 05/08/2015 
#"qqH125_JESDown", #submitted 05/08/2015
#"ggH125_BIG_JESUp", #submitted 05/08/2015
#"ggH125_BIG_JESDown", #submitted 05/08/2015
#"WH_ZH_TTH_HToWW_M125_JESUp", #submitted 05/09/2015
#"WH_ZH_TTH_HToWW_M125_JESDown", #submitted 05/09/2015
#"WH_ZH_TTH_HToZZ_M125_JESUp", #submitted 05/09/2015
#"WH_ZH_TTH_HToZZ_M125_JESDown", #submitted 05/09/2015
#"WH_HToBB_M125_JESUp", #submitted 05/09/2015
#"WH_HToBB_M125_JESDown", #submitted 05/09/2015
#"STopS_T_JESUp", #submitted 05/13/2015
#"STopS_T_JESDown", #submitted 05/13/2015
#"STopS_Tbar_JESUp", #submitted 05/13/2015
#"STopS_Tbar_JESDown", #submitted 05/13/2015
#"STopT_T_JESUp", #submitted 05/14/2015
#"STopT_T_JESDown", #submitted 05/14/2015
#"STopT_Tbar_JESUp", #submitted 05/14/2015
#"STopT_Tbar_JESDown", #submitted 05/14/2015
#"STopTW_T_JESUp", #submitted 05/13/2015
#"STopTW_T_JESDown", #submitted 05/13/2015
#"STopTW_Tbar_JESUp", #submitted 05/13/2015
#"STopTW_Tbar_JESDown", #submitted 05/13/2015
#"WZ_JESUp", #submitted 05/10/2015 #needed to run job 2956 manually as the job became stale and unable to resubmit
#"WZ_JESDown", #submitted 05/10/2015
#"DYJets_JESUp", #submitted 05/11/2015
#"DYJets_JESDown", #submitted 05/11/2015
#"WJets_JESUp", #submitted 05/10/2015
#"WJets_JESDown", #submitted 05/10/2015
#"TTH_HToBB_M125_JESUp", #submitted 05/11/2015
#"TTH_HToBB_M125_JESDown", #submitted 05/11/2015
#"TTbar_JESUp", #submitted 05/14/2015
#"TTbar_JESDown", #submitted 05/14/2015
#"WW_JESUp", #submitted 05/14/2015
#"WW_JESDown", #submitted 05/14/2015
#"WJets_matchingdown",
#"WJets_matchingup",
#"WJets_scaledown",
#"WJets_scaleup"
#"WlnuJets_M-10To50_HighPtLeptonKept062",
#"WlnuJets_M-50_HighPtLeptonKept062",
#"WlnuJets_M-10To50_HighPtLeptonKept070",
#"WlnuJets_M-50_HighPtLeptonKept070",
#"WlnuJets_M-10To50_HighPtLeptonKept075",
#"WlnuJets_M-50_HighPtLeptonKept075",
"WlnuJets_M-10To50_lowerSubLeadingJet_HighPtLeptonKept050mu061el",
"WlnuJets_M-50_lowerSubLeadingJet_HighPtLeptonKept050mu061el"
]


template='crab_scriptExe_cfg.py'
template2='emptyPSet_cfg.py'

nevents = 100.0

for iprocess in processes:
	print "Doing Process " + iprocess + " ... "

	ifile=TFile.Open(basepath+iprocess+".root","READ")
	t=ifile.Get("PS/jets2p")
	n=t.GetEntries()
	njobs = ceil(n/nevents)
	ifile.Close()

	config='crab_'+iprocess+'_cfg.py'
	system('cp '+template+' '+config)
	system('sed s%PROCESSNAME%'+iprocess+'%g '+config+' --in-place')
	system('sed s%NUMBEROFJOBS%'+str(int(njobs))+'%g '+config+' --in-place')
	system('sed s%NUMBEROFEVENTS%'+str(int(nevents))+'%g '+config+' --in-place')

	config2='emptyPSet_'+iprocess+'_cfg.py'
	system('cp '+template2+' '+config2)
	system('sed s%PROCESSNAME%'+iprocess+'%g '+config2+' --in-place')

    #system('crab submit '+config)

