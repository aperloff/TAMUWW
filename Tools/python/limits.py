import os, sys, getopt, argparse, fnmatch, errno, subprocess, tempfile, time, getopt, argparse, glob
from copy import copy
from ROOT import *
from subprocess import call
from glob import glob
from OrderedDict26 import *

class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'
        BOKGREEN = '\033[1m\033[92m'
        BWARNING = '\033[1m\033[93m'
        BFAIL = '\033[1m\033[91m'

def addSumMCToFile(background_only):
	#leptons = ["electron","muon"]
	leptons = ["electron"]
	#variables = ["KinBDT","MEBDT","KinMEBDT"]
	variables = ["MET"]
	jetBinTmp = {
		"jets2":"TwoJ0B",
		#"jets3":"ThreeJ0B",
		#"jets4":"FourJ0B"
	}
	jetBin = OrderedDict(sorted(jetBinTmp.items(), key=lambda t: t[0]))	
	backgrounds = ["WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","WW","QCD_ElFULL","ZJets","WJets"]	
	signals = ["WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","WH125_HToBB","TTH_HToBB_M125","ggH125","qqH125","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125"]	
	processes = []
	nominalFilePath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/"

	if background_only==True:
		processes = backgrounds
	else:
		processes = signals+backgrounds

	for ij, ijet in enumerate(jetBin):
		for il, l in enumerate(leptons):	
			print "Doing lepton="+l+" ijet="+ijet
			nominalFile = TFile(nominalFilePath+ijet+"/"+l+"/histos_"+l+"_"+ijet+".root","READ")
			os.system("cp "+nominalFilePath+ijet+"/"+l+"/histos_"+l+"_"+ijet+"_SysNames.root"+" "+nominalFilePath+ijet+"/"+l+"/histos_"+l+"_"+ijet+"_SysNames_withSumMC.root")
			updatedFile = TFile(nominalFilePath+ijet+"/"+l+"/histos_"+l+"_"+ijet+"_SysNames_withSumMC.root","UPDATE")

			for iv, v in enumerate(variables):
				print "\tDoing variable="+v
				sumMC = TH1D()
				validationString = ""

				for iprocess, proc in enumerate(processes):
					if l=="muon" and proc=="QCD_ElFULL":
						proc = "QCD_MuFULL"
					nominalHistogram = nominalFile.Get(v+"_"+proc+"_"+l)

					if iprocess==0:
						sumMC = nominalHistogram.Clone(v+"_Single"+l[:2].title()+"_Data_"+l)
						sumMC.Reset()

					validationString+=str(nominalHistogram.Integral())+" + "
					sumMC.Add(nominalHistogram)

					updatedFile.cd()
					#nominalHistogram.Write()

				validationString+=" = "+str(sumMC.Integral())
				updatedFile.cd()
				sumMC.Write()
				print "\t\t"+validationString

def combineProcesses(nominalFilename,newFilename,l,v,ijet):
	processesTmp = {
		"WH_ZH_TTH_HToZZ_M125" : ["WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125"],
		"WH125_HToBB"		   : ["WH125_HToBB"],
		"TTH_HToBB_M125"	   : ["TTH_HToBB_M125"],
		"ggH125" 			   : ["ggH125"],
		"qqH125" 			   : ["qqH125"],
		"WH_ZH_TTH_HToWW_M125" : ["WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125"],
		"STop"		 		   : ["STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar"],
		"TTbar"      		   : ["TTbar"],
		"QCD_ElFULL" 		   : ["QCD_ElFULL"],
		"ZJets" 	 		   : ["ZJets"],
		"WJets" 	 		   : ["WJets"],
		"SingleEl_Data"        : ["SingleEl_Data"],
		"Diboson"	 		   : ["WZ","WW"]
	}
	if l=="muon":
		processesTmp["QCD_ElFULL"] = ["QCD_MuFULL"]
		processesTmp["SingleEl_Data"] = ["SingleMu_Data"]
	processes = OrderedDict(sorted(processesTmp.items(), key=lambda t: t[0]))	

	nominalFile = TFile(nominalFilename,"READ")
	newFile = TFile(newFilename,"RECREATE")

	sumMC = 0
	data = 0
	for key,value in processes.iteritems():
		print "\t\tDoing group "+key
		print "\t\t\tSumming:"
		sum = TH1D()

		for counter, ivalue in enumerate(value):
			print "\t\t\t\t"+ivalue
			originalHistogram = nominalFile.Get(v+"_"+ivalue+"_"+l)
			if counter==0:
				sum = originalHistogram.Clone(v+"_"+key+"_"+l)
				sum.Reset()
			sum.Add(originalHistogram)

		newFile.cd()
		sum.Write()
		if key == "SingleEl_Data":
			data = originalHistogram.Integral()
		else:
			sumMC+=sum.Integral()

	if abs(sumMC-data)>0.1:
		print "Something went WRONG!!!! (sumMC="+str(sumMC)+",data="+str(data)+")"

def normalizeHistogramsToNominal(nominal,up,down):
	up.Scale(nominal.Integral()/up.Integral())
	down.Scale(nominal.Integral()/down.Integral())

def getRateAndAddMergeInHistogram(sourceHistogram,destinationHistogram):
	sourceRate = sourceHistogram.Integral()
	destinationRate = destinationHistogram.Integral()
	destinationHistogram.Scale((destinationRate+sourceRate)/destinationRate)

def checkForZeroBin(sourceHistogram):
	hasZeroBin = False
	for ibin in range(1,sourceHistogram.GetNbinsX()+1):
		if sourceHistogram.GetBinContent(ibin)<=0:
			print "WARNING::checkForZeroBin The bin content for bin "+str(ibin)+" in "+str(sourceHistogram.GetName())+" is less than or equal to zero ("+str(sourceHistogram.GetBinContent(ibin))+")"
			print "\tThis will cause problems for combine"
			hasZeroBin = True
	return hasZeroBin

def checkForNegativeBins(sourceHistogram):
	for ibin in range(1,sourceHistogram.GetNbinsX()+1):
		if sourceHistogram.GetBinContent(ibin)<0:
			print "WARNING::checkForNegativeBins The bin content for bin "+str(ibin)+" in "+str(sourceHistogram.GetName())+" is negative ("+str(sourceHistogram.GetBinContent(ibin))+")"
			print "\tThis will cause problems for combine"
			exit(2);

def maxDifferenceShapes(nominal,up,down):
	for ibin in range(1,nominal.GetNbinsX()+1):
		nominalUpDiff = abs(nominal.GetBinContent(ibin)-up.GetBinContent(ibin))
		nominalDownDiff = abs(nominal.GetBinContent(ibin)-up.GetBinContent(ibin))
		up.SetBinContent(ibin,nominal.GetBinContent(ibin)+max(nominalUpDiff,nominalDownDiff))
		down.SetBinContent(ibin,max(0.0,nominal.GetBinContent(ibin)-max(nominalUpDiff,nominalDownDiff)))

def setupShapeSystematicsInNominalFile(use_sum_MC, background_only, simplify, combineProc, addStatUncert, scaleSignal=1.0, printSum=False):
	leptons = ["electron","muon"]
	#leptons = ["electron"]
	#variables = ["KinBDT","MEBDT","KinMEBDT","MET"]
	variables = ["KinMEBDT"]
	jetBinTmp = {
		"jets2":"TwoJ0B",
		"jets3":"ThreeJ0B",
		"jets4":"FourJ0B"
	}
	jetBin = OrderedDict(sorted(jetBinTmp.items(), key=lambda t: t[0]))	
	datas = ["SingleEl_Data","SingleMu_Data"]
	signals = ["WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","WH125_HToBB","TTH_HToBB_M125","ggH125","qqH125","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125"]
	backgrounds = ["WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","WW","ZJets","WJets","QCD_ElFULL","QCD_MuFULL"]
	processes = []
	if simplify==1:
		processes = ["ggH125","WJets","SingleEl_Data"]
	elif simplify==2:
		processes = ["ggH125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","WJets","SingleEl_Data"]
	elif simplify==3:
		processes = ["ggH125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","WW","WJets","SingleEl_Data"]
	elif simplify==4: #This one has problems
		processes = ["ggH125","WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","WJets","SingleEl_Data"]
	elif simplify==5:
		processes = ["ggH125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","WJets","SingleEl_Data"]
	elif simplify==6:
		processes = ["ggH125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","ZJets","WJets","SingleEl_Data"]
	elif simplify==7: #This one has problems
		processes = ["ggH125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","WW","WJets","SingleEl_Data"]
	elif simplify==8:
		processes = ["ggH125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","ZJets","WJets","SingleEl_Data"]
	elif simplify==99:
		processes = ["WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","WH125_HToBB","TTH_HToBB_M125","ggH125","qqH125","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","ZJets","WJets","QCD_ElFULL","SingleEl_Data"]
	else:
		processes = ["WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","WH125_HToBB","TTH_HToBB_M125","ggH125","qqH125","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125","WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","WW","ZJets","WJets","QCD_ElFULL","SingleEl_Data"]
	shapeSystematics = {
		"JES"       : "CMS_scale_j_shape",
		"matching"  : "CMS_hww_lnujj_matching_shape",
		"scale"     : "CMS_hww_lnujj_scale_shape",
		"TopPt"     : "CMS_hww_lnujj_topPtWeight_shape",
		"QCDEtaW"   : "CMS_hww_lnujj_QCDEtaWeight_shape",
		"CSVWeight" : "CMS_hww_lnujj_CSVWeight_shape",
		"PUWeight"  : "CMS_hww_lnujj_PUWeight_shape"
	}
	suffix = ""
	if use_sum_MC and not background_only:
		suffix = "_SumMC"
	elif use_sum_MC and background_only:
		suffix = "_SumBackgrounds"
	if simplify:
		suffix += "_Simplify"
	if scaleSignal!=1.0:
		suffix += "_SignalsScaled"

	#nominalFilePath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal"
	nominalFilePath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_12_11_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal"

	for ij, ijet in enumerate(jetBin):
		for il, l in enumerate(leptons):	
			print "Doing lepton="+l+" ijet="+ijet
			nominalFilename = nominalFilePath+"/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+".root"
			nominalFile = TFile(nominalFilename,"READ")
			newFilename = nominalFilePath+"/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+"_SysNames"+suffix+".root"
			newFile = TFile(newFilename,"RECREATE")

			for iv, v in enumerate(variables):
				print "\tDoing variable="+v
				sumBackgrounds = None
				sumSignals = None
				sumData = None
				sumMC = TH1D()
				validationString = ""
				formatString = "\t\t%-20s%-20s"
				print formatString % ("Process","Integral",)
				print formatString % ("-"*7,"-"*8,)

				for processCount, process in enumerate(processes):
					if l=="muon" and process=="QCD_ElFULL":
						process = "QCD_MuFULL"
					if l=="muon" and process=="SingleEl_Data":
						process = "SingleMu_Data"

					originalHistogram = nominalFile.Get(v+"_"+process+"_"+l)
					if "Data" in process:
						print bcolors.OKGREEN+formatString % (process,originalHistogram.Integral(),)+bcolors.ENDC
					elif process in signals:
						print bcolors.BFAIL+formatString % (process,originalHistogram.Integral(),)+bcolors.ENDC
					else:
						print bcolors.BWARNING+formatString % (process,originalHistogram.Integral(),)+bcolors.ENDC

					if process=="WH125_HToBB":
						upHistogram = nominalFile.Get(v+"_WH_HToBB_M125_JESUp_"+l)
						downHistogram = nominalFile.Get(v+"_WH_HToBB_M125_JESDown_"+l)
					else:	
						upHistogram = nominalFile.Get(v+"_"+process+"_JESUp_"+l)
						downHistogram = nominalFile.Get(v+"_"+process+"_JESDown_"+l)

					if upHistogram and downHistogram:
						upHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["JES"]+"Up")
						downHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["JES"]+"Down")
						#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)

					if simplify==99 and not sumMC and processes=="WJets":
						wzHistogram = nominalFile.Get(v+"_WZ_"+l)
						wwHistogram = nominalFile.Get(v+"_WW_"+l)
						wzRate = wzHistogram.Integral()
						wwRate = wwHistogram.Integral()
						originalRate = originalHistogram.Integral()
						originalHistogram.Scale((originalRate+wzRate+wwRate)/originalRate)
					'''
					if process == "WJets":
						qcdHistogram = nominalFile.Get(v+"_QCD_ElFULL_"+l)
						#getRateAndAddMergeInHistogram(qcdHistogram,originalHistogram)
						qcdRate = qcdHistogram.Integral()
						originalRate = originalHistogram.Integral()
						originalHistogram.Scale((originalRate+qcdRate)/originalRate)
					if process == "QCD_ElFULL":
						originalHistogram.Scale(0)
					'''
					
					newFile.cd()

					'''
					Check for negative bins which might cause problems for combine
					'''
					if originalHistogram:
						checkForNegativeBins(originalHistogram)
					else:
						print "ERROR::setupShapeSystematicsInNominalFile can't find the nominal histogram"
						print "\tEither the process was not included when creating the templates or there were no events after all of the cuts"
						continue

					if process in signals and scaleSignal!=1.0:
						originalHistogram.Scale(scaleSignal)
						if upHistogram and downHistogram:
							upHistogram.Scale(scaleSignal)
							downHistogram.Scale(scaleSignal)

					if not use_sum_MC or (process != "SingleEl_Data" and process != "SingleMu_Data"):
						originalHistogram.Write()
					if upHistogram and downHistogram:
						checkForNegativeBins(upHistogram)
						checkForNegativeBins(downHistogram)
						upHistogram.Write()
						downHistogram.Write()

					'''
					Do CSV weight uncertainty
					'''
					if "QCD" not in process and "Data" not in process:
						CSVUpFile = TFile(nominalFilePath[:-8]+"/CSVWeightSysUp/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+".root","READ")
						upHistogram = CSVUpFile.Get(v+"_"+process+"_"+l)
						CSVDownFile = TFile(nominalFilePath[:-8]+"/CSVWeightSysDown/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+".root","READ")
						downHistogram = CSVDownFile.Get(v+"_"+process+"_"+l)
						#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
						upHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["CSVWeight"]+"Up")
						downHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["CSVWeight"]+"Down")
						newFile.cd()
						checkForNegativeBins(upHistogram)
						checkForNegativeBins(downHistogram)
						upHistogram.Write()
						downHistogram.Write()

					'''
					Do PU weight uncertainty
					'''
					if "QCD" not in process and "Data" not in process:
						CSVUpFile = TFile(nominalFilePath[:-8]+"/PUWeightSysUp/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+".root","READ")
						upHistogram = CSVUpFile.Get(v+"_"+process+"_"+l)
						CSVDownFile = TFile(nominalFilePath[:-8]+"/PUWeightSysDown/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+".root","READ")
						downHistogram = CSVDownFile.Get(v+"_"+process+"_"+l)
						#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
						upHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["PUWeight"]+"Up")
						downHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["PUWeight"]+"Down")
						newFile.cd()
						checkForNegativeBins(upHistogram)
						checkForNegativeBins(downHistogram)
						upHistogram.Write()
						downHistogram.Write()

					if process == "WJets":
						upHistogram = nominalFile.Get(v+"_"+process+"_matchingup_"+l)
						downHistogram = nominalFile.Get(v+"_"+process+"_matchingdown_"+l)
						if upHistogram and downHistogram:
							#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
							#maxDifferenceShapes(originalHistogram,upHistogram,downHistogram)
							upHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["matching"]+"Up")
							downHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["matching"]+"Down")
							newFile.cd()
							checkForNegativeBins(upHistogram)
							checkForNegativeBins(downHistogram)
							upHistogram.Write()
							downHistogram.Write()

						upHistogram = nominalFile.Get(v+"_"+process+"_scaleup_"+l)
						downHistogram = nominalFile.Get(v+"_"+process+"_scaledown_"+l)
						if upHistogram and downHistogram:
							#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
							upHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["scale"]+"Up")
							downHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["scale"]+"Down")
							newFile.cd()
							checkForNegativeBins(upHistogram)
							checkForNegativeBins(downHistogram)
							upHistogram.Write()
							downHistogram.Write()

						QCDEtaWUpFile = TFile(nominalFilePath[:-8]+"/QCDEtaWeightUp/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+".root","READ")
						upHistogram = QCDEtaWUpFile.Get(v+"_"+process+"_"+l)
						QCDEtaWDownFile = TFile(nominalFilePath[:-8]+"/QCDEtaWeightDown/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+".root","READ")
						downHistogram = QCDEtaWDownFile.Get(v+"_"+process+"_"+l)
						#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
						upHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["QCDEtaW"]+"Up")
						downHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["QCDEtaW"]+"Down")
						newFile.cd()
						checkForNegativeBins(upHistogram)
						checkForNegativeBins(downHistogram)
						upHistogram.Write()
						downHistogram.Write()

					if process == "QCD_ElFULL" or process == "QCD_MuFULL":
						QCDEtaWUpFile = TFile(nominalFilePath[:-8]+"/QCDEtaWeightUp/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+".root","READ")
						upHistogram = QCDEtaWUpFile.Get(v+"_"+process+"_"+l)
						QCDEtaWDownFile = TFile(nominalFilePath[:-8]+"/QCDEtaWeightDown/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+".root","READ")
						downHistogram = QCDEtaWDownFile.Get(v+"_"+process+"_"+l)
						#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
						upHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["QCDEtaW"]+"Up")
						downHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["QCDEtaW"]+"Down")
						newFile.cd()
						checkForNegativeBins(upHistogram)
						checkForNegativeBins(downHistogram)
						upHistogram.Write()
						downHistogram.Write()

					if process == "TTbar":
						TTbarUpFile = TFile(nominalFilePath[:-8]+"/TopPtWeightSysUp/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+".root","READ")
						upHistogram = TTbarUpFile.Get(v+"_"+process+"_"+l)
						TTbarDownFile = TFile(nominalFilePath[:-8]+"/TopPtWeightSysDown/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+".root","READ")
						downHistogram = TTbarDownFile.Get(v+"_"+process+"_"+l)
						#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
						upHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["TopPt"]+"Up")
						downHistogram.SetName(v+"_"+process+"_"+l+"_"+shapeSystematics["TopPt"]+"Down")
						newFile.cd()
						checkForNegativeBins(upHistogram)
						checkForNegativeBins(downHistogram)
						upHistogram.Write()
						downHistogram.Write()

					if processCount!=0 and "Data" not in process:
						validationString += " + "
					if any(signal in originalHistogram.GetName() for signal in signals):
						if not sumSignals:
							sumSignals = originalHistogram.Clone(originalHistogram.GetName()+"_clone")
						else:
							sumSignals.Add(originalHistogram)
						if not background_only:
							validationString+=bcolors.BFAIL+str(originalHistogram.Integral())+bcolors.ENDC#+" + "
					elif any(background in originalHistogram.GetName() for background in backgrounds):
						if not sumBackgrounds:
							sumBackgrounds = originalHistogram.Clone(originalHistogram.GetName()+"_clone")
						else:
							sumBackgrounds.Add(originalHistogram)
						validationString+=bcolors.BWARNING+str(originalHistogram.Integral())+bcolors.ENDC#+" + "
					elif any(data in originalHistogram.GetName() for data in datas):
						if not sumData:
							sumData = originalHistogram.Clone(originalHistogram.GetName()+"_clone")
						else:
							sumData.Add(originalHistogram)

					tmpProcesses = copy(processes)
					if l=="muon" and ("QCD_ElFULL" in tmpProcesses or "SingleEl_Data" in tmpProcesses):
						tmpProcesses[tmpProcesses.index("QCD_ElFULL")] = "QCD_MuFULL"
						tmpProcesses[tmpProcesses.index("SingleEl_Data")] = "SingleMu_Data"

					if tmpProcesses.index(process)==0:
						sumMC = originalHistogram.Clone(v+"_Single"+l[:2].title()+"_Data_"+l)
						sumMC.Reset()

				checkForZeroBin(sumBackgrounds)
				sumMC.Add(sumBackgrounds)
				if not background_only:
					sumMC.Add(sumSignals)

				if printSum:
					validationString+=" = "+bcolors.OKBLUE+str(sumMC.Integral())+bcolors.ENDC
					print "\t\t"+validationString
				else:
					print bcolors.OKBLUE+formatString % ("sumMC",sumMC.Integral(),)+bcolors.ENDC

				if use_sum_MC:
					sumMC.Write()

				newFile.Close()
				nominalFile.Close()
				if combineProc:
					combineProcesses(newFilename,nominalFilePath+"/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+"_SysNames"+suffix+"_CombinedProc.root",l,v,ijet)

				if addStatUncert:
					newFilenameTwo = nominalFilePath+"/"+ijet+"/"+l+"/histos_"+l+"_"+ijet+"_SysNames"+suffix+"_StatShapes"+".root"
					matcher=v+"_*_"+l
					command = "python $CMSSW_BASE/src/TAMUWW/Tools/python/add_stat_shapes.py --filter "+matcher+" --prefix CMS_hww_lnujj_statError --suffix "+jetBin[ijet]+" --threshold 0.05 "+newFilename+" "+newFilenameTwo
					os.system(command)

def produceLimitsDatacard():
	debug = False
	leptons = ["electron","muon"]
	#variables = ["KinBDT","MEBDT","KinMEBDT","MET"]
	variables = ["KinMEBDT"]
	jetBin = {
		"jets2":"TwoJ0B",
		"jets3":"ThreeJ0B",
		"jets4":"FourJ0B"
	}
	batch = " -batch true"
	suffix = ""
	suffix+=" -suffix _StatShapes"
	#simplify = 90
	simplify = 0
	processes = []
	overrideProcesses = ""
	if simplify>0:
		if simplify==1:
			processes = ["ggH125","WJets","SingleEl_Data"]
		elif simplify==2:
			processes = ["ggH125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","WJets","SingleEl_Data"]
		elif simplify==3:
			processes = ["ggH125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","WW","WJets","SingleEl_Data"]
		elif simplify==4: #This one has problems
			processes = ["ggH125","WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","WJets","SingleEl_Data"]
		elif simplify==5:
			processes = ["ggH125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","WJets","SingleEl_Data"]
		elif simplify==6:
			processes = ["ggH125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","ZJets","WJets","SingleEl_Data"]
		elif simplify==7: #This one has problems
			processes = ["ggH125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","WW","WJets","SingleEl_Data"]
		elif simplify==8:
			processes = ["ggH125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","ZJets","WJets","SingleEl_Data"]
		elif simplify==90:
			processes = ["ggH125","qqH125","WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","WW","ZJets","WJets","QCD_ElFULL","SingleEl_Data"]
		elif simplify==99:
			processes = ["WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","WH125_HToBB","TTH_HToBB_M125","ggH125","qqH125","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","ZJets","WJets","QCD_ElFULL","SingleEl_Data"]
		overrideProcesses =" -overrideProcesses "+' '.join(processes)
		suffix+=""
		#suffix+=" -suffix _signalsScaled"
		#suffix+=" -suffix _SumBackgrounds"
		#suffix+=" -suffix _Simplify"
		#suffix+=" -suffix _SumMC_Simplify"
		#suffix+=" -suffix _SumBackgrounds_Simplify"

	for l in leptons:
		for v in variables:
			for ijet in jetBin:
	
				basepath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_12_11_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/"+ijet+"/"+l+"/"
				if simplify>0:
					tmpProcesses = copy(processes)
					if l=="muon" and ("QCD_ElFULL" in tmpProcesses or "SingleEl_Data" in tmpProcesses):
						tmpProcesses[tmpProcesses.index("QCD_ElFULL")] = "QCD_MuFULL"
						tmpProcesses[tmpProcesses.index("SingleEl_Data")] = "SingleMu_Data"
					overrideProcesses =" -overrideProcesses "+' '.join(tmpProcesses)
	
				command = "ProduceLimitsDatacard_x -outpath ./ -basepath "+basepath+" -variable "+v+" -channelName "+jetBin[ijet]+" -lepton "+l+" -jetBin "+ijet.lower()+" -tagcat eq0tag"+overrideProcesses+suffix+batch
				if debug:
					print "***********"
					print "* TESTING *"
					print "***********"
					print command
				else:
					print command
					os.system(command)

#setupShapeSystematicsInNominalFile(use_sum_MC=False,background_only=False,simplify=0,combineProc=False,addStatUncert=True,scaleSignal=1.0, printSum=False)
#addSumMCToFile(False)
produceLimitsDatacard()
