#! /usr/bin/env python
import os, sys, getopt, argparse, fnmatch, errno, subprocess, tempfile
from subprocess import call
from glob import glob
from time import sleep

#output_file = "/eos/uscms/store/user/eusebi/Winter12to13ME8TeV/rootOutput/WJets_partialHadd/WJets_part"
#output_file = "/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/2015_06_05_microNtuples_optimized/SingleMu_Data_19p279fb_partialHadd/SingleMu_Data_19p279fb_part"
#output_file = "/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/2015_06_05_microNtuples_optimized/tmp2/TTbar_JESUp_partialHadd/TTbar_JESUp_part"
#output_file = "/uscms_data/d2/aperloff/DYJetsToLL_M-10to50_script"
#output_file = "/uscms_data/d2/aperloff/WlnuJets_M-50_HighPtLeptonKept062"
#output_file = "/uscms_data/d2/aperloff/WlnuJets_M-50_HighPtLeptonKept075"
output_file = "/uscms_data/d2/aperloff/DYJetsToLL_M-10to50_lowerSubLeadingJet"

USEGLOB = False
#files = glob("/eos/uscms/store/user/aperloff/MatrixElement/PS_outfiles_20130916_MC13/WJets_part1_nominal/PS_*.root")
#files = glob("/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEResults/rootOutput/WJets_JESUp/WJets_JESUp*.root")
#files = glob("root://cmseos.fnal.gov//store/user/eusebi/Winter12to13ME8TeV/rootOutput/SingleMu_Data_19p279fb/SingleMu_Data_19p279fb*.root")
#files = glob("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/2015_06_05_microNtuples_optimized/tmp2/TTbar_JESUp/TTbar_JESUp*.root")
#files = os.system("xrdfs root://cmseos.fnal.gov/ ls -u $STORE/MatrixElement/PS_outfiles_20150202_MC18/")
#files = subprocess.check_output(["xrdfs","root://cmseos.fnal.gov/","ls -u","$STORE/MatrixElement/PS_outfiles_20150202_MC18/"], stderr=subprocess.STDOUT)
#lsCmd = "xrdfs root://cmseos.fnal.gov/ ls -u /store/user/aperloff/MatrixElement/PS_outfiles_20150202_MC18/"
#lsCmd = "eos root://cmseos.fnal.gov/ find /store/user/aperloff/MatrixElement/Summer12ME8TeV/MEResults/rootOutput/WlnuJets_M-50_HighPtLeptonKept062/crab_run_MatrixElement_WlnuJets_M-50_HighPtLeptonKept062/160309_233526/ | grep -v failed | grep -v log | grep -F .root"
#lsCmd = "eos root://cmseos.fnal.gov/ find /store/user/aperloff/MatrixElement/Summer12ME8TeV/MEResults/rootOutput/WlnuJets_M-50_HighPtLeptonKept070/crab_run_MatrixElement_WlnuJets_M-50_HighPtLeptonKept070_v2/160320_062201/ | grep -v failed | grep -v log | grep -F .root"
#lsCmd = "eos root://cmseos.fnal.gov/ find /store/user/aperloff/MatrixElement/Summer12ME8TeV/MEResults/rootOutput/WlnuJets_M-50_HighPtLeptonKept075/crab_run_MatrixElement_WlnuJets_M-50_HighPtLeptonKept075_v2/160320_062308/ | grep -v failed | grep -v log | grep -F .root"
lsCmd = "eos root://cmseos.fnal.gov/ find /store/user/aperloff/MatrixElement/PS_outfiles_20150202_MC20/ | grep -F .root"
p = subprocess.Popen(lsCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
files, err = p.communicate()

DEBUG = 0
if DEBUG:
	print "Files direct from ls:\n"
	print files
	print "lines:",str(files.count("\n"))

if USEGLOB:
	file_split = len(files)/900
else:
	files = files.split("\n")
	files = filter(bool, files)
	file_split = len(files)/900
	files = [f.replace('/eos/uscms', 'root://cmseos.fnal.gov/') for f in files]
	if DEBUG:
		print "Files after splitting and filtering:\n"
		print files
		print "length:",str(len(files))
		print "Number times to split:",file_split

for i in range(1+file_split):
	current_set_of_files = files[i*900:min(len(files),(i+1)*900)]
	if DEBUG:
		print len(current_set_of_files)
	command = "nohup hadd " + output_file+str(i)+".root "
	for f in current_set_of_files:
		command += f+" "
	if DEBUG:
		print command
	script=open(output_file+str(i)+".sh","w")
	script.write("#! /usr/bin/sh\n\n")
	script.write(command)
	script.write("\n\necho \"HADD COMPLETE\"\n")
	script.close()
	out=open(output_file+str(i)+".stdout","w")
	subprocess.Popen(["bash",output_file+str(i)+".sh"],stdout=out,stderr=subprocess.STDOUT)
	#subprocess.Popen(["nohup","hadd",output_file+str(i)+".root",input_files],stdout=out,stderr=subprocess.STDOUT)
