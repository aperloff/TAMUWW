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
	leptons = ["electron","muon"]
	variables = ["KinBDT","MEBDT","KinMEBDT","MET"]
	jetBinTmp = {
		"jets2":"TwoJ0B",
		"jets3":"ThreeJ0B",
		"jets4":"FourJ0B"
	}
	BDTCat = ["HighKinBDT","LowKinBDT"]

	jetBin = OrderedDict(sorted(jetBinTmp.items(), key=lambda t: t[0]))	
	backgrounds = ["WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","WW","QCD_ElFULL","ZJets","WJets"]	
	signals = ["WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","WH125_HToBB","TTH_HToBB_M125","ggH125","qqH125","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125"]	
	processes = []
	#nominalFilePath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_10_29_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/"
	nominalFilePath = "/uscms_data/d2/aperloff//Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/"

	if background_only==True:
		processes = backgrounds
	else:
		processes = signals+backgrounds

	for ij, ijet in enumerate(jetBin):
		for il, l in enumerate(leptons):	
			for ibc, bc in enumerate(BDTCat): 
				print "Doing lepton="+l+" ijet="+ijet
				nominalFile = TFile(nominalFilePath+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
				os.system("cp "+nominalFilePath+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+"_"+bc+"_SysNames.root"+" "+nominalFilePath+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+"_"+bc+"_SysNames_withSumMC.root")
				updatedFile = TFile(nominalFilePath+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+"_"+bc+"_SysNames_withSumMC.root","UPDATE")

				for iv, v in enumerate(variables):
					print "\tDoing variable="+v
					sumMC = TH1D()
					validationString = ""

					for iprocess, proc in enumerate(processes):
						if l=="muon" and proc=="QCD_ElFULL":
							proc = "QCD_MuFULL"
						nominalHistogram = nominalFile.Get(v+"_"+proc+"_"+l)

						if iprocess==0:
							sumMC = nominalHistogram.Clone(v+"_Single"+l[:2].title()+"_Data_"+l+"_asimov")
							sumMC.Reset()

						validationString+=str(nominalHistogram.Integral())+" + "
						sumMC.Add(nominalHistogram)

						updatedFile.cd()
						#nominalHistogram.Write()

					validationString+=" = "+str(sumMC.Integral())
					updatedFile.cd()
					sumMC.Write()
					print "\t\t"+validationString

def combineLeptons():
	'''
	Example:
	from limits import combineLeptons
	combineLeptons("/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal_AllPlots/jets2/electron/HighKinBDT/histos_electron_jets2.root","/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal_AllPlots/jets2/muon/HighKinBDT/histos_muon_jets2.root","histos.root")
	'''
	jetBinTmp = {
		"jets2":"TwoJ0B",
		"jets3":"ThreeJ0B",
		"jets4":"FourJ0B"
	}
	jetBin = OrderedDict(sorted(jetBinTmp.items(), key=lambda t: t[0]))	
	basename = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal_AllPlots/"

	for ij, ijet in enumerate(jetBin):
		print "Doing ijet="+ijet
		electronFilename = basename+"/"+ijet+"/electron/histos_electron_"+ijet+"_SysNames.root"
		electronFile = TFile(electronFilename,"READ")
		muonFilename = basename+"/"+ijet+"/muon/histos_muon_"+ijet+"_SysNames.root"
		muonFile = TFile(muonFilename,"READ")
		histDict = OrderedDict()

		for key in electronFile.GetListOfKeys():
			h = TH1D(electronFile.Get(key.GetName()))
			if "BDT" in h.GetName():
				continue
			histDict[h.GetName()] = h
		#print "Adding histograms ... "
		for key in muonFile.GetListOfKeys():
			h = TH1D(muonFile.Get(key.GetName()))
			nameToFind = h.GetName()
			if "BDT" in nameToFind:
				continue
			if "muon" in nameToFind:
				nameToFind = nameToFind.replace("muon","electron")
			if "SingleMu_Data" in nameToFind:
				nameToFind = nameToFind.replace("Mu","El")
			if "QCD_MuFULL" in nameToFind:
				nameToFind = nameToFind.replace("Mu","El")
			#print "\t",nameToFind
			histDict[nameToFind].Add(h)
	
		newFilename = basename+"/"+ijet+"/both/histos_both_"+ijet+".root"
		newFile = TFile(newFilename,"RECREATE")
		#print "Saving histograms ..."
		for key, value in histDict.iteritems():
			#print "\t",key
			value.Write()
	
		newFile.Close()
		electronFile.Close()
		muonFile.Close()

def combineProcesses(nominalFilename,newFilename,l,v,ijet):
	'''
	from limits import combineProcesses
	for j in ["jets2","jets3","jets4"]:
	    combineProcesses("/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/"+j+"/histos_"+j+"_SysNames.root","/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/combinedByCategory/histos_"+j+"_SysNames_combinedProcesses.root",["electron","muon"],"KinMEBDT",j)
	hadd histos_jets2_SysNames_combinedProcesses.root histos_electron_jets2_SysNames_combinedProcesses.root histos_muon_jets2_SysNames_combinedProcesses.root
	'''
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
		"Diboson"	 		   : ["WZ","WW"],
		"Total Background"     : ["WZ","WW","WJets","ZJets","QCD_ElFULL","QCD_MuFULL","TTbar","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar"],
		"Total \HWW"           : ["ggH125","qqH125","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125"],
		"Total Volunteer"      : ["WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","WH125_HToBB","TTH_HToBB_M125"],
	}
	if len(l)>1:
		processesTmp["QCD_ElFULL"] = ["QCD_ElFULL","QCD_MuFULL"]
		processesTmp["SingleEl_Data"] = ["SingleEl_Data","SingleMu_Data"]
	elif l[0]=="muon":
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
			for ilep, lep in enumerate(l):
				originalHistogram = nominalFile.Get(v+"_"+ivalue+"_"+lep)
				if originalHistogram==None:
					continue
				print "\t\t\t\t"+ivalue
				if counter==0 and ilep==0:
					sum = originalHistogram.Clone(v+"_"+key)
					sum.Reset()
				sum.Add(originalHistogram)

		err = Double(0.0)
		integral = sum.IntegralAndError(1,sum.GetNbinsX(),err) 
		print "\t\t\tIntegralAndError:",'$'+'{0:.2f}'.format(integral)+"\pm"+'{0:.2f}'.format(err)+'$'

		newFile.cd()
		sum.Write()
		if key == "SingleEl_Data":
			data = sum.Integral()
		elif kin not in ["Total Background","Total \HWW","Total Volunteer"]:
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
	#variables = ["KinBDT","MEBDT","KinMEBDT"]
	variables = ["KinMEBDT"]
	jetBinTmp = {
		"jets2":"TwoJ0B",
		"jets3":"ThreeJ0B",
		"jets4":"FourJ0B"
	}
	jetBin = OrderedDict(sorted(jetBinTmp.items(), key=lambda t: t[0]))	
	#BDTCat = ["HighKinBDT","LowKinBDT"]
	BDTCat = ["HailMaryLoose"]
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
		processes = ["WH_HToZZ_M125","ZH_HToZZ_M125","TTH_HToZZ_M125","WH125_HToBB","TTH_HToBB_M125","ggH125","qqH125","WH_HToWW_M125","ZH_HToWW_M125","TTH_HToWW_M125","WZ","STopS_T","STopS_Tbar","STopT_T","STopT_Tbar","STopTW_T","STopTW_Tbar","TTbar","WW","QCD_ElFULL","ZJets","WJets","SingleEl_Data"]
	shapeSystematics = {
		"JES"       : "CMS_scale_j_shape",
		"matching"  : "CMS_hww_lnujj_matching_shape",
		"scale"     : "CMS_hww_lnujj_scale_shape",
		"TopPt"     : "CMS_hww_lnujj_topPtWeight_shape",
		"QCDEtaW"   : "CMS_hww_lnujj_QCDEtaWeight_shape",
		"CSVWeight" : "CMS_hww_lnujj_CSVWeight_shape",
		"PUWeight"  : "CMS_hww_lnujj_PUWeight_shape",
		"CosThetaL" : "CMS_hww_lnujj_CosThetaLWeight_shape"
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
	#nominalFilePath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_12_11_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal"
	#nominalFilePath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal"
	nominalFilePath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_11_14_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale_HailMaryLoose/nominal"

	for ij, ijet in enumerate(jetBin):
		for il, l in enumerate(leptons):
			for ibc, bc in enumerate(BDTCat): 
				print "Doing lepton="+l+" ijet="+ijet+" BDTCat="+bc
				iopath = nominalFilePath+"/"+ijet+"/"+l+"/"+bc+"/"
				nominalFilename = iopath+"histos_"+l+"_"+ijet+".root"
				nominalFile = TFile(nominalFilename,"READ")
				newFilename = nominalFilePath+"/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+"_"+bc+"_SysNames"+suffix+".root"
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
							upHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["JES"]+"Up",shapeSystematics["JES"]+"Up")
							downHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["JES"]+"Down",shapeSystematics["JES"]+"Down")
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
							CSVUpFile = TFile(nominalFilePath[:-8]+"/CSVWeightSysUp/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
							upHistogram = CSVUpFile.Get(v+"_"+process+"_"+l)
							CSVDownFile = TFile(nominalFilePath[:-8]+"/CSVWeightSysDown/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
							downHistogram = CSVDownFile.Get(v+"_"+process+"_"+l)
							#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
							upHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["CSVWeight"]+"Up",shapeSystematics["CSVWeight"]+"Up")
							downHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["CSVWeight"]+"Down",shapeSystematics["CSVWeight"]+"Down")
							newFile.cd()
							checkForNegativeBins(upHistogram)
							checkForNegativeBins(downHistogram)
							upHistogram.Write()
							downHistogram.Write()
	
						'''
						Do PU weight uncertainty
						'''
						if "QCD" not in process and "Data" not in process:
							CSVUpFile = TFile(nominalFilePath[:-8]+"/PUWeightSysUp/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
							upHistogram = CSVUpFile.Get(v+"_"+process+"_"+l)
							CSVDownFile = TFile(nominalFilePath[:-8]+"/PUWeightSysDown/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
							downHistogram = CSVDownFile.Get(v+"_"+process+"_"+l)
							#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
							upHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["PUWeight"]+"Up",shapeSystematics["PUWeight"]+"Up")
							downHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["PUWeight"]+"Down",shapeSystematics["PUWeight"]+"Down")
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
								upHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["matching"]+"Up",shapeSystematics["matching"]+"Up")
								downHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["matching"]+"Down",shapeSystematics["matching"]+"Down")
								newFile.cd()
								checkForNegativeBins(upHistogram)
								checkForNegativeBins(downHistogram)
								upHistogram.Write()
								downHistogram.Write()
	
							upHistogram = nominalFile.Get(v+"_"+process+"_scaleup_"+l)
							downHistogram = nominalFile.Get(v+"_"+process+"_scaledown_"+l)
							if upHistogram and downHistogram:
								#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
								upHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["scale"]+"Up",shapeSystematics["scale"]+"Up")
								downHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["scale"]+"Down",shapeSystematics["scale"]+"Down")
								newFile.cd()
								checkForNegativeBins(upHistogram)
								checkForNegativeBins(downHistogram)
								upHistogram.Write()
								downHistogram.Write()
	
							QCDEtaWUpFile = TFile(nominalFilePath[:-8]+"/QCDEtaWeightUp/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
							upHistogram = QCDEtaWUpFile.Get(v+"_"+process+"_"+l)
							QCDEtaWDownFile = TFile(nominalFilePath[:-8]+"/QCDEtaWeightDown/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
							downHistogram = QCDEtaWDownFile.Get(v+"_"+process+"_"+l)
							#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
							upHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["QCDEtaW"]+"Up",shapeSystematics["QCDEtaW"]+"Up")
							downHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["QCDEtaW"]+"Down",shapeSystematics["QCDEtaW"]+"Down")
							newFile.cd()
							checkForNegativeBins(upHistogram)
							checkForNegativeBins(downHistogram)
							upHistogram.Write()
							downHistogram.Write()

							CosThetaLUpFile = TFile(nominalFilePath[:-8]+"/CosThetaLWeightSysUp/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
							upHistogram = CosThetaLUpFile.Get(v+"_"+process+"_"+l)
							CosThetaLDownFile = TFile(nominalFilePath[:-8]+"/CosThetaLWeightSysDown/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
							downHistogram = CosThetaLDownFile.Get(v+"_"+process+"_"+l)
							#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
							upHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["CosThetaL"]+"Up",shapeSystematics["CosThetaL"]+"Up")
							downHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["CosThetaL"]+"Down",shapeSystematics["CosThetaL"]+"Down")
							newFile.cd()
							checkForNegativeBins(upHistogram)
							checkForNegativeBins(downHistogram)
							upHistogram.Write()
							downHistogram.Write()
	
						if process == "QCD_ElFULL" or process == "QCD_MuFULL":
							QCDEtaWUpFile = TFile(nominalFilePath[:-8]+"/QCDEtaWeightUp/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
							upHistogram = QCDEtaWUpFile.Get(v+"_"+process+"_"+l)
							QCDEtaWDownFile = TFile(nominalFilePath[:-8]+"/QCDEtaWeightDown/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
							downHistogram = QCDEtaWDownFile.Get(v+"_"+process+"_"+l)
							#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
							upHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["QCDEtaW"]+"Up",shapeSystematics["QCDEtaW"]+"Up")
							downHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["QCDEtaW"]+"Down",shapeSystematics["QCDEtaW"]+"Down")
							newFile.cd()
							checkForNegativeBins(upHistogram)
							checkForNegativeBins(downHistogram)
							upHistogram.Write()
							downHistogram.Write()
	
						if process == "TTbar":
							TTbarUpFile = TFile(nominalFilePath[:-8]+"/TopPtWeightSysUp/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
							upHistogram = TTbarUpFile.Get(v+"_"+process+"_"+l)
							TTbarDownFile = TFile(nominalFilePath[:-8]+"/TopPtWeightSysDown/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+".root","READ")
							downHistogram = TTbarDownFile.Get(v+"_"+process+"_"+l)
							#normalizeHistogramsToNominal(originalHistogram,upHistogram,downHistogram)
							upHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["TopPt"]+"Up",shapeSystematics["TopPt"]+"Up")
							downHistogram.SetNameTitle(v+"_"+process+"_"+l+"_"+shapeSystematics["TopPt"]+"Down",shapeSystematics["TopPt"]+"Down")
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
	
					if combineProc:
						combineProcesses(newFilename,nominalFilePath+"/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+"_"+bc+"_SysNames"+suffix+"_CombinedProc.root",l,v,ijet)

				newFile.Close()
				nominalFile.Close()
	
				if addStatUncert:
					minimums = {
						"jets2" : {"electron" : 37, "muon": 37},
						"jets3" : {"electron" : 29, "muon": 29},
						"jets4" : {"electron" : 30, "muon": 30},
					}
					newFilenameTwo = nominalFilePath+"/"+ijet+"/"+l+"/"+bc+"/histos_"+l+"_"+ijet+"_"+bc+"_SysNames"+suffix+"_StatShapes"+".root"
					matcher="*_*_"+l
					command = "python $CMSSW_BASE/src/TAMUWW/Tools/python/add_stat_shapes.py --filter "+matcher+" --prefix CMS_hww_lnujj_statError --suffix "+jetBin[ijet]+"_"+l+"_"+bc+" --minimum-bin "+minimums[ijet][l]+" --scale 1.0 --threshold 0.05 "+newFilename+" "+newFilenameTwo
					os.system(command)

def combineBDTBins(debug=False):
	leptons = ["electron","muon"]
	jetBin = {
		"jets2":"TwoJ0B",
		"jets3":"ThreeJ0B",
		"jets4":"FourJ0B"
	}
	BDTCat = ["HighKinBDT","LowKinBDT"]

	basepath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/"
	#basepath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal_AllPlots/"

	for ijet in jetBin:
		for l in leptons:
			outputPath = basepath+ijet+"/"+l+"/"
			outputFilename = "histos_"+l+"_"+ijet+"_SysNames.root"
			inputFiles = ""
			for ibc, bc in enumerate(BDTCat): 
				inputFiles+=outputPath+bc+"/histos_"+l+"_"+ijet+"_"+bc+"_SysNames.root "
				#For use with nominal_AllPlots because systematics not run with all plots
				#inputFiles+=outputPath+bc+"/histos_"+l+"_"+ijet+".root "

			command = "hadd -f " + outputPath+outputFilename + " " + inputFiles
			if debug:
				print "***********"
				print "* TESTING *"
				print "***********"
				print command
			else:
				print command
				os.system(command)

def produceLimitsDatacard(doStatErrors=False, debug=False):
	leptons = ["electron","muon"]
	#variables = ["KinBDT","MEBDT","KinMEBDT","MET"]
	#variables = ["KinBDT","MEBDT","KinMEBDT"]
	variables = ["KinMEBDT"]
	jetBin = {
		"jets2":"TwoJ0B",
		"jets3":"ThreeJ0B",
		"jets4":"FourJ0B"
	}
	#BDTCat = ["HighKinBDT","LowKinBDT",""]
	BDTCat = ["HailMaryLoose"]
	batch = " -batch true"
	suffix = ""
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
				for ibc, bc in enumerate(BDTCat): 

					#basepath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2015_12_11_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/"+ijet+"/"+l+"/"
					#basepath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/nominal/"+ijet+"/"+l+"/"+bc+"/"
					basepath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_11_14_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale_HailMaryLoose/nominal/"+ijet+"/"+l+"/"+bc+"/"
					if simplify>0:
						tmpProcesses = copy(processes)
						if l=="muon" and ("QCD_ElFULL" in tmpProcesses or "SingleEl_Data" in tmpProcesses):
							tmpProcesses[tmpProcesses.index("QCD_ElFULL")] = "QCD_MuFULL"
							tmpProcesses[tmpProcesses.index("SingleEl_Data")] = "SingleMu_Data"
						overrideProcesses =" -overrideProcesses "+' '.join(tmpProcesses)
	
					channelName = jetBin[ijet]+"_"+l
					if bc != "":
						channelName += "_"+bc
						suffix = " -suffix _"+bc

					command = "ProduceLimitsDatacard_x -outpath "+basepath+" -basepath "+basepath+" -variable "+v+" -channelName "+channelName+" -lepton "+l+" -jetBin "+ijet.lower()+" -tagcat eq0tag -doStatErrors "+('true' if doStatErrors else 'false')+overrideProcesses+suffix+batch
					if debug:
						print "***********"
						print "* TESTING *"
						print "***********"
						print command
					else:
						print command
						os.system(command)

def combineCardsByJetBin(debug=False):
	leptons = ["electron","muon"]
	#variables = ["KinBDT","MEBDT","KinMEBDT","MET"]
	variables = ["KinBDT","MEBDT","KinMEBDT"]
	#variables = ["KinMEBDT"]
	jetBin = {
		"jets2":"TwoJ0B",
		"jets3":"ThreeJ0B",
		"jets4":"FourJ0B"
	}
	#BDTCat = ["HighKinBDT","LowKinBDT"]
	BDTCat = [""]

	basepath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/"
	outputFolder = "datacards/"

	for iv, v in enumerate(variables):
		for ijet in jetBin:
			outputPath = basepath+outputFolder
			outputFilename = "DataCard_"+jetBin[ijet]+"_CombinedLepton_"+v+".txt"
			bins = ""
			for l in leptons:
				for ibc, bc in enumerate(BDTCat): 
					pathToCurrentCard=basepath+"nominal/"+ijet+"/"+l+"/"+bc+"/"
					if bc=="":
						currentCardName="DataCard_"+jetBin[ijet]+"_"+l+"_"+v+".txt"
						channelName = jetBin[ijet]+"_"+l
					else:
						currentCardName="DataCard_"+jetBin[ijet]+"_"+l+"_"+bc+"_"+v+".txt"
						channelName = jetBin[ijet]+"_"+l+"_"+bc
					bins+=channelName+"="+pathToCurrentCard+currentCardName+" "

			command = "combineCards.py "+bins+" > "+outputPath+outputFilename
			if debug:
				print "***********"
				print "* TESTING *"
				print "***********"
				print command
			else:
				print command
				os.system(command)

def combineCardsByLeptonBin(debug=False):
	leptons = ["electron","muon"]
	#variables = ["KinBDT","MEBDT","KinMEBDT","MET"]
	variables = ["KinBDT","MEBDT","KinMEBDT"]
	#variables = ["KinMEBDT"]
	jetBin = {
		"jets2":"TwoJ0B",
		"jets3":"ThreeJ0B",
		"jets4":"FourJ0B"
	}
	#BDTCat = ["HighKinBDT","LowKinBDT"]
	BDTCat = [""]

	basepath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/"
	outputFolder = "datacards/"

	for v in variables:
		for l in leptons:
			outputPath = basepath+outputFolder
			outputFilename = "DataCard_CombinedJet_"+l+"_"+v+".txt"
			bins = ""
			for ijet in jetBin:
				for ibc, bc in enumerate(BDTCat): 
					pathToCurrentCard=basepath+"nominal/"+ijet+"/"+l+"/"+bc+"/"
					if bc=="":
						currentCardName="DataCard_"+jetBin[ijet]+"_"+l+"_"+v+".txt"
						channelName = jetBin[ijet]+"_"+l
					else:
						currentCardName="DataCard_"+jetBin[ijet]+"_"+l+"_"+bc+"_"+v+".txt"
						channelName = jetBin[ijet]+"_"+l+"_"+bc
					bins+=channelName+"="+pathToCurrentCard+currentCardName+" "

			command = "combineCards.py "+bins+" > "+outputPath+outputFilename
			if debug:
				print "***********"
				print "* TESTING *"
				print "***********"
				print command
			else:
				print command
				os.system(command)

def combineAllCards(debug=False):
	leptons = ["electron","muon"]
	#variables = ["KinBDT","MEBDT","KinMEBDT","MET"]
	#variables = ["KinBDT","MEBDT","KinMEBDT"]
	variables = ["KinMEBDT"]
	jetBin = {
		"jets2":"TwoJ0B",
		"jets3":"ThreeJ0B",
		"jets4":"FourJ0B"
	}
	#BDTCat = ["HighKinBDT","LowKinBDT"]
	#BDTCat = ["HailMaryLoose"]
	BDTCat = [""]

	basepath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/"
	#basepath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_11_14_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale_HailMaryLoose/"
	outputFolder = "datacards/"

	for v in variables:
		outputPath = basepath+outputFolder
		outputFilename = "DataCard_CombinedJet_CombinedLepton_"+v+".txt"
		bins = ""
		for ijet in jetBin:
			for l in leptons:
				if ijet=="jets4" and l=="electron":
					continue
				#if ijet=="jets2" and l=="electron":
				#	continue
				for ibc, bc in enumerate(BDTCat): 
					pathToCurrentCard=basepath+"nominal/"+ijet+"/"+l+"/"+bc+"/"
					if bc=="":
						currentCardName="DataCard_"+jetBin[ijet]+"_"+l+"_"+v+".txt"
						channelName = jetBin[ijet]+"_"+l
					else:
						currentCardName="DataCard_"+jetBin[ijet]+"_"+l+"_"+bc+"_"+v+".txt"
						channelName = jetBin[ijet]+"_"+l+"_"+bc
					bins+=channelName+"="+pathToCurrentCard+currentCardName+" "

		command = "combineCards.py "+bins+" > "+outputPath+outputFilename
		if debug:
			print "***********"
			print "* TESTING *"
			print "***********"
			print command
		else:
			print command
			os.system(command)

def makeLinks(debug=False):
	leptons = ["electron","muon"]
	#variables = ["KinBDT","MEBDT","KinMEBDT","MET"]
	#variables = ["KinBDT","MEBDT","KinMEBDT"]
	variables = ["KinMEBDT"]
	jetBin = {
		"jets2":"TwoJ0B",
		"jets3":"ThreeJ0B",
		"jets4":"FourJ0B"
	}
	#BDTCat = ["HighKinBDT","LowKinBDT"]
	BDTCat = ["HailMaryLoose"]

	#basepath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_09_18_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale/"
	basepath = "/uscms_data/d2/aperloff/Summer12ME8TeV/2017_11_14_LimitHistograms_PU_CSV_TTbar_QCDEta_Scale_HailMaryLoose/"
	outputFolder = "datacards/"
	outputPath = basepath+outputFolder

	for ijet in jetBin:
		for l in leptons:
			for v in variables:
				for ibc, bc in enumerate(BDTCat): 
					pathToCurrentInfo=basepath+"nominal/"+ijet+"/"+l+"/"+bc+"/"
					if bc=="":
						pathToCurrentDatacard = pathToCurrentInfo+"DataCard_"+jetBin[ijet]+"_"+l+"_"+v+".txt"
						pathToCurrentFile = pathToCurrentInfo + "histos_"+l+"_"+ijet+"_SysNames.root"
					else:
						pathToCurrentDatacard=pathToCurrentInfo+"DataCard_"+jetBin[ijet]+"_"+l+"_"+bc+"_"+v+".txt"
						pathToCurrentFile = pathToCurrentInfo + "histos_"+l+"_"+ijet+"_"+bc+"_SysNames.root"
					command1 = "ln -s "+pathToCurrentDatacard +" "+ outputPath+os.path.basename(pathToCurrentDatacard)
					command2 = "ln -s "+pathToCurrentFile     +" "+ outputPath+os.path.basename(pathToCurrentFile)
					if debug:
						print "***********"
						print "* TESTING *"
						print "***********"
						print command1
						print command2
					else:
						print command1
						print command2
						os.system(command1)
						os.system(command2)

#setupShapeSystematicsInNominalFile(use_sum_MC=False,background_only=False,simplify=0,combineProc=False,addStatUncert=False,scaleSignal=1.0, printSum=False)
#addSumMCToFile(True)
#combineBDTBins(False)
#produceLimitsDatacard(doStatErrors=False, debug=False)
#combineCardsByJetBin(debug=False)
#combineCardsByLeptonBin(debug=False)
combineAllCards(debug=False)
#makeLinks(debug=False)
#combineLeptons()