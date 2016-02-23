import os, fnmatch, math
from ROOT import *

basepath = "/eos/uscms/store/user/"
eusebi = "eusebi/Winter12to13ME8TeV/rootOutput/"
eusebi2 = "eusebi/MatrixElement/Summer12ME8TeV/MEResults/rootOutput/"
aperloff = "aperloff/MatrixElement/Summer12ME8TeV/MEResults/rootOutput/"

template_path = "/uscms_data/d2/aperloff/YOURWORKINGAREA/MatrixElement/gitty/CMSSW_5_3_2_patch5/src/TAMUWW/Tools/python/"
template=template_path+'microTemplate_BaseCuts.jdl'
template2=template_path+'microTemplate_BaseCuts.sh'

processes = {
#"WW"							: (eusebi,   ["WW"]),						#Submitted on 05/18/2015
#"WZ"							: (eusebi,   ["WZ"]),						#Submitted on 05/18/2015
#"WJets"							: (eusebi,   ["WJets"]),				#Submitted on 05/18/2015
#"ZJets"							: (eusebi,   ["ZJets"]),				#Submitted on 05/18/2015
#"TTbar"							: (eusebi,   ["TTbar"]),					#Submitted on 05/18/2015
#"STopS_T"						: (eusebi,   ["STopS_T"]),					#Submitted on 05/18/2015
#"STopS_Tbar"					: (eusebi,   ["STopS_Tbar"]),				#Submitted on 05/18/2015
#"STopT_T"						: (eusebi,   ["STopT_T"]),					#Submitted on 05/18/2015
#"STopT_Tbar"					: (eusebi,   ["STopT_Tbar"]),				#Submitted on 05/18/2015
#"STopTW_T"						: (eusebi,   ["STopTW_T"]),					#Submitted on 05/18/2015
#"STopTW_Tbar"					: (eusebi,   ["STopTW_Tbar"]),				#Submitted on 05/18/2015
#"SingleEl_Full_Subset"			: (eusebi,   ["SingleEl_Full_Subset"]),		#Submitted on 05/18/2015
#"SingleMu_Full"					: (eusebi,   ["SingleMu_Full"]),		#Submitted on 05/18/2015
#"SingleEl_Data_19p148fb"		: (eusebi,   ["SingleEl_Data_19p148fb"]),	#Submitted on 05/18/2015
#"SingleMu_Data_19p279fb"		: (eusebi,   ["SingleMu_Data_19p279fb"]),	#Submitted on 05/18/2015
#"ggH125_BIG"					: (eusebi,   ["ggH125_BIG"]),				#Submitted on 05/18/2015
#"qqH125"						: (aperloff, ["qqH125"]),					#Submitted on 05/15/2015
#"WH_ZH_TTH_HToWW_M125"			: (eusebi,   ["WH_ZH_TTH_HToWW_M125"]),		#Submitted on 05/18/2015
#"WH_ZH_TTH_HToZZ_M125"			: (eusebi,   ["WH_ZH_TTH_HToZZ_M125"]),		#Submitted on 05/18/2015
#"WH125_HToBB"					: (eusebi,   ["WH125_HToBB"]),				#Submitted on 05/18/2015
#"TTH_HToBB_M125"				: (eusebi,   ["TTH_HToBB_M125"]),			#Submitted on 05/18/2015
"WW_JESUp"						: (eusebi, ["WW_JESUp_partialHadd"]),
"WW_JESDown"					: (eusebi, ["WW_JESDown_partialHadd"]),
#"WZ_JESUp"						: (aperloff, ["WZ_JESUp"]),					#Submitted on 05/20/2015
#"WZ_JESDown"					: (aperloff, ["WZ_JESDown"]),				#Submitted on 05/20/2015
"WJets_JESUp"					: (eusebi, ["WJets_JESUp_partialHadd"]),
"WJets_JESDown"					: (eusebi, ["WJets_JESDown_partialHadd"]), 
#"WJets_matchingup"              : (aperloff,   ["WJets_matchingup"]),        #Submitted on 06/04/2015
#"WJets_matchingdown"            : (aperloff,   ["WJets_matchingdown"]),      #Submitted on 06/04/2015
#"WJets_scaleup"                 : (aperloff,   ["WJets_scaleup"]),           #Submitted on 06/04/2015
#"WJets_scaledown"               : (aperloff,   ["WJets_scaledown"]),         #Submitted on 06/04/2015
#"DYJets_JESUp"					: (aperloff, ["DYJets_JESUp"]),             #Submitted on 06/04/2015
#"DYJets_JESDown"				: (aperloff, ["DYJets_JESDown"]),           #Submitted on 06/04/2015
"TTbar_JESUp"					: (eusebi, ["TTbar_JESUp_partialHadd"]),
"TTbar_JESDown"					: (eusebi, ["TTbar_JESDown_partialHadd"]),
#"STopS_T_JESUp"					: (aperloff, ["STopS_T_JESUp"]),				#Submitted on 05/20/2015
#"STopS_T_JESDown"				: (aperloff, ["STopS_T_JESDown"]),				#Submitted on 05/20/2015
#"STopS_Tbar_JESUp"				: (aperloff, ["STopS_Tbar_JESUp"]),				#Submitted on 05/20/2015
#"STopS_Tbar_JESDown"			: (aperloff, ["STopS_Tbar_JESDown"]),			#Submitted on 05/20/2015
#"STopT_T_JESUp"					: (aperloff, ["STopT_T_JESUp"]),				#Submitted on 05/20/2015
#"STopT_T_JESDown"				: (aperloff, ["STopT_T_JESDown"]),				#Submitted on 05/20/2015
#"STopT_Tbar_JESUp"				: (aperloff, ["STopT_Tbar_JESUp"]),				#Submitted on 05/20/2015
#"STopT_Tbar_JESDown"			: (aperloff, ["STopT_Tbar_JESDown"]),			#Submitted on 05/20/2015
#"STopTW_T_JESUp"					: (aperloff, ["STopTW_T_JESUp"]),				#Submitted on 05/20/2015
#"STopTW_T_JESDown"				: (aperloff, ["STopTW_T_JESDown"]),				#Submitted on 05/20/2015
#"STopTW_Tbar_JESUp"				: (aperloff, ["STopTW_Tbar_JESUp"]),			#Submitted on 05/20/2015
#"STopTW_Tbar_JESDown"			: (aperloff, ["STopTW_Tbar_JESDown"]),			#Submitted on 05/20/2015
#"ggH125_BIG_JESUp"				: (aperloff, ["ggH125_BIG_JESUp"]),				#Submitted on 05/20/2015
#"ggH125_BIG_JESDown"			: (aperloff, ["ggH125_BIG_JESDown"]),			#Submitted on 05/20/2015
#"qqH125_JESUp"					: (aperloff, ["qqH125_JESUp"]),					#Submitted on 05/20/2015
#"qqH125_JESDown"				: (aperloff, ["qqH125_JESDown"]),				#Submitted on 05/20/2015
#"WH_HToBB_M125_JESUp"			: (aperloff, ["WH_HToBB_M125_JESUp"]),			#Submitted on 05/20/2015
#"WH_HToBB_M125_JESDown"			: (aperloff, ["WH_HToBB_M125_JESDown"]),		#Submitted on 05/20/2015
#"WH_ZH_TTH_HToWW_M125_JESUp"	: (aperloff, ["WH_ZH_TTH_HToWW_M125_JESUp"]),	#Submitted on 05/20/2015
#"WH_ZH_TTH_HToWW_M125_JESDown"	: (aperloff, ["WH_ZH_TTH_HToWW_M125_JESDown"]),	#Submitted on 05/20/2015
#"WH_ZH_TTH_HToZZ_M125_JESUp"	: (aperloff, ["WH_ZH_TTH_HToZZ_M125_JESUp"]),	#Submitted on 05/20/2015
#"WH_ZH_TTH_HToZZ_M125_JESDown"	: (aperloff, ["WH_ZH_TTH_HToZZ_M125_JESDown"]),	#Submitted on 05/20/2015
#"TTH_HToBB_M125_JESUp"			: (aperloff, ["TTH_HToBB_M125_JESUp"]),			#Submitted on 05/20/2015
#"TTH_HToBB_M125_JESDown"		: (aperloff, ["TTH_HToBB_M125_JESDown"])		#Submitted on 05/20/2015
}

def long_substr(data):
    substr = ''
    if len(data) > 1 and len(data[0]) > 0:
        for i in range(len(data[0])):
            for j in range(len(data[0])-i+1):
                if j > len(substr) and is_substr(data[0][i:i+j], data):
                    substr = data[0][i:i+j]
    return substr

def is_substr(find, data):
    if len(data) < 1 and len(find) < 1:
        return False
    for i in range(len(data)):
        if find not in data[i]:
            return False
    return True

def make_scripts(splitJob):
	for key,value in processes.iteritems():
		processnames = ""
		inputpaths = ""
		adddir = ""
		startentry = 0
		endentry = -1
		outputsuffix = "\" \""
		queueNumber = 1
		entriesPerJob = 50000
	
		print "Doing process group " + key + " ... "
	
		tmpStr = "\tCombining processes "
		first = True
		for jprocess in value[1]:
			#Find the part of the path after basepath+value[0]+jprocess
			matches = []
			for root, dirnames, filenames in os.walk(basepath+value[0]+jprocess):
				for filename in fnmatch.filter(filenames, '*.root'):
					#matches.append(os.path.join(root, filename).replace(basepath+value[0]+jprocess+"/",""))
					matches.append(os.path.join(root, filename))

			common_dir = os.path.dirname(long_substr(matches))
			#print "common_dir = " + common_dir
			matches = next(os.walk(common_dir))[1]
			if len(matches)<2:
				matches = []
			else:
				matches = [str(x) + "/" for x in next(os.walk(common_dir))[1]]
				#matches = [common_dir+"/"+str(x)+"/" for x in next(os.walk(common_dir))[1]]
	
			if first:
				tmpStr+=jprocess
				processnames+=jprocess
				inputpaths = common_dir+"/"
				adddir = " ".join(matches)
				first = False
			else:
				tmpStr+=", "+jprocess
				processnames+=" "+jprocess
				inputpaths+= common_dir+"/"
				adddir = " " + " ".join(matches)
		if adddir.isspace() or adddir == "":
			adddir = "/"
		print tmpStr
		print "\t"+processnames
		print "\tGetting files from "+basepath+value[0]
		print "\tinputpaths = "+inputpaths
		print "\tadddir = "+adddir
	
		config='micro'+key+'_BaseCuts.jdl'
		print '\tcp '+template+' '+config
		os.system('cp '+template+' '+config)
		os.system('sed s%PROCESSGROUP%'+key+'%g '+config+' --in-place')
		if splitJob:
			mergeFile = TFile(basepath+"aperloff/MatrixElement/Summer12ME8TeV/MEInput/"+processnames+".root","READ")
			mergeTree = mergeFile.Get("PS/jets2p")
			totalEntries = mergeTree.GetEntries()
			queueNumber = math.ceil(totalEntries/float(entriesPerJob))
		else:
			mergeFile = TFile(basepath+"aperloff/MatrixElement/Summer12ME8TeV/MEInput/"+processnames+".root","READ")
			mergeTree = mergeFile.Get("PS/jets2p")
			entriesPerJob = mergeTree.GetEntries()+1
		os.system('sed s%EVENTSPERJOB%'+str(entriesPerJob)+'%g '+config+' --in-place')
		os.system('sed s%QUEUE%'+str(int(queueNumber))+'%g '+config+' --in-place')

		config2='micro'+key+'_BaseCuts.sh'
		print '\tcp '+template2+' '+config2
		os.system('cp '+template2+' '+config2)
		os.system('sed \"s%PROCESSNAMES%'+processnames+'%g\" '+config2+' --in-place')
		os.system('sed \"s%ADDDIR%'+adddir+'%g\" '+config2+' --in-place')
		os.system('sed \"s%INPUTPATH%'+inputpaths+'%g\" '+config2+' --in-place')
#		os.system('sed \"s%UPDATEMICRONTUPLE%/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEResults/2015_02_27_microNtuples_optimized_BDT/micro'+processnames+'_BaseCuts.root%g\" '+config2+' --in-place')

		'''
		if splitJob:
			mergeFile = TFile(basepath+"aperloff/MatrixElement/Summer12ME8TeV/MEInput/"+processnames+".root","READ")
			mergeTree = mergeFile.Get("PS/jets2p")
			totalEntries = mergeTree.GetEntries()
			entriesPerJob = 10000
			entry_split = math.ceil(totalEntries/float(entriesPerJob))
			print "\tMergeTree has "+str(totalEntries)+" entries"
			print "\tSplitting mergeTree entries into prieces of size "+str(entriesPerJob)
			print "\tThere are now "+str(entry_split)+" jobs"
			for i in range(entry_split):
				startentry=i*entriesPerJob
				endentry=min((i+1)*entriesPerJob,totalEntries+1)
				outputsuffix = "_"+str(i)
				print "\tJob "+str(i)+": STARTENTRY="+str(startentry)+"\tENDENTRY="+str(endentry)+"\tOUTPUTSUFFIX="+str(outputsuffix)

				os.system('sed s%QUEUE%'+str(queueNumber)+'%g '+config+' --in-place')

				os.system('sed \"s%STARTENTRY%'+str(startentry)+'%g\" '+config+' --in-place')
				os.system('sed \"s%ENDENTRY%'+str(endentry)+'%g\" '+config+' --in-place')
				os.system('sed \"s%OUTPUTSIFFIX%'+str(outputsuffix)+'%g\" '+config+' --in-place')
		else:

			os.system('sed \"s%STARTENTRY%'+str(startentry)+'%g\" '+config+' --in-place')
			os.system('sed \"s%ENDENTRY%'+str(endentry)+'%g\" '+config+' --in-place')
			os.system('sed \"s%OUTPUTSIFFIX%'+str(outputsuffix)+'%g\" '+config+' --in-place')
		'''
	
		#os.system('condor_submit '+config)

def make_filelist():
	for key,value in processes.iteritems():
		processnames = ""
		inputpaths = ""
		adddir = ""
	
		print "Doing process group " + key + " ... "
	
		tmpStr = "\tCombining processes "
		first = True

		for jprocess in value[1]:
			outputFilename = jprocess+"_fileList.txt"
			f = open(outputFilename, 'w')
			#Find the part of the path after basepath+value[0]+jprocess
			matches = []
			for root, dirnames, filenames in os.walk(basepath+value[0]+jprocess):
				for filename in fnmatch.filter(filenames, '*.root'):
					matches.append(os.path.join(root, filename))
					f.write("root://cmseos.fnal.gov//"+os.path.join(root, filename).replace("/eos/uscms","")+"\n")
			common_dir = os.path.dirname(long_substr(matches))
			f.close()
			if "eusebi" not in value[0]:
				os.system('mv '+outputFilename+' '+common_dir)
			else:
				print "WARNING::multiMakeMicroNtuple Must move file lists by hand as input files are in Dr. Ricardo Eusebi's eos area."

make_filelist()
#make_scripts(False)
