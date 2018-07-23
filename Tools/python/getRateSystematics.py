import os, sys, getopt, argparse, fnmatch, errno, subprocess, tempfile, time, getopt, argparse, glob
from ROOT import *
from OrderedDict26 import *

def getRatesFromSeparateFiles():
	leptons = ["electron","muon"]
	#variables = ["KinBDT","MEBDT","KinMEBDT"]
	variables = ["KinMEBDT"]
	jetBinTmp = {
		"jets2":"TwoJ0B",
		"jets3":"ThreeJ0B",
		"jets4":"FourJ0B"
	}
	jetBin = OrderedDict(sorted(jetBinTmp.items(), key=lambda t: t[0]))
	BDTCat = ["HighKinBDT","LowKinBDT"]

	processes = ["WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","WH125_HToBB","TTH_HToBB_M125","ggH125","qqH125","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125","WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","WW","QCD_ElFULL","ZJets","WJets"]
	#processes = ["TTbar"]
	
	#nominalFilePath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/"
	nominalFilePath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/"
	#UpPath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/PUWeightSysUp/"
	#UpPath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/TopPtWeightSysUp/"
	#UpPath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/CSVWeightSysUp/"
	#UpPath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/QCDEtaWeightUp/"
	#UpPath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/TopPtWeightSysUp/"
	UpPath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/PUWeightSysUp/"
	#DownPath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/PUWeightSysDown/"
	#DownPath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/TopPtWeightSysDown/"
	#DownPath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/CSVWeightSysDown/"
	#DownPath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/QCDEtaWeightDown/"
	#DownPath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/TopPtWeightSysDown/"
	DownPath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/PUWeightSysDown/"
	
	for iv, v in enumerate(variables):
		if iv==0:
			print "if(variable.CompareTo(\""+v+"\")==0) {"
		else:
			print "else if(variable.CompareTo(\""+v+"\")==0) {"
		for ij, ijet in enumerate(jetBin):
			if ij==0:
				print "\tif(jetBin==DEFS::"+ijet.lower()+") {"
			else:
				print "\telse if(jetBin==DEFS::"+ijet.lower()+") {"	
			for il, l in enumerate(leptons):	
				if il==0:
					print "\t\tif(leptonCat==DEFS::"+l+") {"
				else:
					print "\t\telse if(leptonCat==DEFS::"+l+") {"
				#print "Doing variable="+v+" lepton="+l+" ijet="+ijet
				nominalFile = TFile(nominalFilePath+ijet+"/"+l+"/"+BDTCat[0]+"/histos_"+l+"_"+ijet+".root","READ")
				UpFile = TFile(UpPath+ijet+"/"+l+"/"+BDTCat[0]+"/histos_"+l+"_"+ijet+".root","READ")
				DownFile = TFile(DownPath+ijet+"/"+l+"/"+BDTCat[0]+"/histos_"+l+"_"+ijet+".root","READ")
	
				for iprocess, proc in enumerate(processes):
					if l=="muon" and proc=="QCD_ElFULL":
						proc = "QCD_MuFULL"
					if proc=="TTbar":
						iprocess = 17
					#if (proc=="WJets" or proc=="QCD_ElFULL" or proc=="QCD_MuFULL") and ("QCDEtaWeight" not in UpPath) and ("QCDEtaWeight" not in DownPath):
					if (proc=="QCD_ElFULL" or proc=="QCD_MuFULL") and ("QCDEtaWeight" not in UpPath) and ("QCDEtaWeight" not in DownPath):
						continue
					if "TopPtWeightSys" in UpPath and proc!="TTbar":
						continue
					nominalHistogram = nominalFile.Get(v+"_"+proc+"_"+l)
					UpHistogram = UpFile.Get(v+"_"+proc+"_"+l)
					DownHistogram = DownFile.Get(v+"_"+proc+"_"+l)
	
					nominalInt = nominalHistogram.Integral(0,DownHistogram.GetNbinsX()+1)
					UpInt = UpHistogram.Integral(0,DownHistogram.GetNbinsX()+1)
					DownInt = DownHistogram.Integral(0,DownHistogram.GetNbinsX()+1)
					
					if proc=="WH_HToBB_M125":
						proc="WH125_HToBB"

					if UpInt==0 or DownInt==0:
						print "\t\t\tsystematics.back().updateValue(DEFS::PhysicsProcess::getProcessType(\""+proc+"\"), make_pair(\"-\",\"-\"), verbose);"
					else:
						print "\t\t\tsystematics.back().updateValue(DEFS::PhysicsProcess::getProcessType(\""+proc+"\"), make_pair(\""+str("%.3f" % (DownInt/nominalInt))+"\",\""+str("%.3f" % (UpInt/nominalInt))+"\"), verbose);"
				print "\t\t}"
			print "\t}"
		print "}"

def getJESRates():
	leptons = ["electron","muon"]
	variables = ["KinBDT","MEBDT","KinMEBDT"]
	jetBinTmp = {
		"Jets2":"TwoJ0B",
		"Jets3":"ThreeJ0B",
		"Jets4":"FourJ0B"
	}
	jetBin = OrderedDict(sorted(jetBinTmp.items(), key=lambda t: t[0]))	
	processes = ["WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","WH125_HToBB","TTH_HToBB_M125","ggH125","qqH125","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125","WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","WW","QCD_ElFULL","ZJets","WJets"]
	nominalFilePath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/"

	for iv, v in enumerate(variables):
		if iv==0:
			print "if(variable.CompareTo(\""+v+"\")==0) {"
		else:
			print "else if(variable.CompareTo(\""+v+"\")==0) {"
		for ij, ijet in enumerate(jetBin):
			if ij==0:
				print "\tif(jetBin==DEFS::"+ijet.lower()+") {"
			else:
				print "\telse if(jetBin==DEFS::"+ijet.lower()+") {"	
			for il, l in enumerate(leptons):	
				if il==0:
					print "\t\tif(leptonCat==DEFS::"+l+") {"
				else:
					print "\t\telse if(leptonCat==DEFS::"+l+") {"

				nominalFile = TFile(nominalFilePath+ijet+"/"+l+"/histos_"+l+".root","READ")
				for iprocess, proc in enumerate(processes):
					if l=="muon" and proc=="QCD_ElFULL":
						proc = "QCD_MuFULL"
					if proc=="TTbar":
						iprocess = 17
					nominalHistogram = nominalFile.Get(v+"_"+proc+"_"+l)
					if proc=="WH125_HToBB":
						proc="WH_HToBB_M125"
					if "QCD" in proc:
						continue
					UpHistogram = nominalFile.Get(v+"_"+proc+"_JESUp_"+l)
					DownHistogram = nominalFile.Get(v+"_"+proc+"_JESDown_"+l)
	
					nominalInt = nominalHistogram.Integral(0,DownHistogram.GetNbinsX()+1)
					UpInt = UpHistogram.Integral(0,DownHistogram.GetNbinsX()+1)
					DownInt = DownHistogram.Integral(0,DownHistogram.GetNbinsX()+1)
					
					if proc=="WH_HToBB_M125":
						proc="WH125_HToBB"

					if UpInt==0 or DownInt==0:
						print "\t\t\tsystematics.back().updateValue(DEFS::PhysicsProcess::getProcessType(\""+proc+"\"), make_pair(\"-\",\"-\"), verbose);"
					else:
						print "\t\t\tsystematics.back().updateValue(DEFS::PhysicsProcess::getProcessType(\""+proc+"\"), make_pair(\""+str("%.3f" % (nominalInt/DownInt))+"\",\""+str("%.3f" % (nominalInt/UpInt))+"\"), verbose);"
				print "\t\t}"
			print "\t}"
		print "}"



getRatesFromSeparateFiles()
#getJESRates()





