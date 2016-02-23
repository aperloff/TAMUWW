#! /usr/bin/env python
import os, sys, getopt, argparse, fnmatch, errno, subprocess, tempfile, math, pprint
from subprocess import call
from glob import glob
from time import sleep
from ROOT import *

basepath = ""
output_dir = ""
output_root = ""
output_logs = ""
output_updateScripts = ""

def set_basepath():
	global basepath
	basepath = "/eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEResults/2015_07_17_microNtuples_optimized_BDT/rootFiles/"

def set_output_directories():
	global output_dir
	global output_root
	global output_logs
	global output_updateScripts
	output_dir = "/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/2015_07_23_microNtuples_optimized/"
	output_root = output_dir+"rootFiles/"
	output_logs = output_dir+"logs/"
	output_updateScripts = output_dir+"updateScripts/"

def update_files(debug):
	files = glob(basepath+"/*.root")
	
	print str(len(files))+" files to update\n"
	
	print "***************"
	print "* INPUT FILES *"
	print "***************"
	pprint.pprint(files)
	
	print "\n"
	print "************"
	print "* COMMANDS *"
	print "************"
	
	for i in files:
		ifile = i.replace(basepath,"")
		iprocess = ifile.replace("micro","")
		iprocess = iprocess.replace("_BaseCuts.root","")
		ilog = output_logs+ifile.replace(".root",".stdout")
		iscript = output_updateScripts+ifile.replace(".root",".sh")
	
		command = "nohup makeMicroNtuple_x -inputPaths "+basepath+" -addDir / -outputPath "+output_root+" -processes "+iprocess+" -fillBDT true -debug false -updateMicroNtuple "+i
	
		print command+"\n"
	
		if not debug:
			script=open(iscript,"w")
			script.write("#! /usr/bin/sh\n\n")
			script.write(command)
			script.write("\necho \"UPDATE COMPLETE\"\n")
			script.close()
			out=open(ilog,"w")
			subprocess.Popen(["bash",iscript],stdout=out,stderr=subprocess.STDOUT)

def check_logs_for_completeness():
	processes_complete = []
	processes_incomplete = []

	files = glob(output_logs+"/*.stdout")

	for i in files:
		logfile = open(i, 'r')
		loglist = logfile.readlines()
		logfile.close()
		found = False

		iprocess = i.replace(output_logs,"")
		iprocess = iprocess.replace("micro","")
		iprocess = iprocess.replace("_BaseCuts.stdout","")

		for line in loglist:
			if str("UPDATE COMPLETE") in line:
				found = True

		if found:
			processes_complete.append(iprocess)
		else:
			processes_incomplete.append(iprocess)

	print "***************************"
	print "* COMPLETE UPDATES ("+str(len(processes_complete))+"/"+str(len(files))+") *"
	print "***************************"
	pprint.pprint(processes_complete)

	print "*****************************"
	print "* INCOMPLETE UPDATES ("+str(len(processes_incomplete))+"/"+str(len(files))+")*"
	print "*****************************"
	pprint.pprint(processes_incomplete)

def validate_files():
	new_files = glob(output_root+"/*.root")
	old_files = []
	num_invalid_files = 0
	for counter, ifile in enumerate(new_files):
		print "("+str(counter+1)+"/"+str(len(new_files))+") Validating file "+ifile+" ... "
		new_file = TFile(ifile,"READ")
		new_tree = new_file.Get("METree")
		new_entry_count = new_tree.GetEntries()
		old_files.append(ifile.replace(output_root,basepath))
		old_file = TFile(ifile.replace(output_root,basepath),"READ")
		old_tree = old_file.Get("METree")
		old_entry_count = old_tree.GetEntries()

		if new_entry_count != old_entry_count:
			print "\tWARNING::validate_files The new_entry_count ("+str(new_entry_count)+") does not equal the old_entry_count ("+str(old_entry_count)+")"
			num_invalid_files+=1

	print "WARNING::There are "+str(num_invalid_files)+" invalid files"

set_basepath()
set_output_directories()
#update_files(False)
check_logs_for_completeness()
#validate_files()
