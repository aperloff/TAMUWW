#! /usr/bin/env python
import os, sys, getopt, argparse, fnmatch, errno, subprocess, tempfile
from subprocess import call
from glob import glob

output_folder = "/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEResults/rootOutput/WJets_scaleup/"

for root,dirnames,filenames in os.walk("/uscms_data/d2/aperloff/Summer12ME8TeV/MEResults/2015_06_05_microNtuples_optimized/tmp2/WJets_scaleup/"):
	for filename in fnmatch.filter(filenames,'*.root'):
		ifile = os.path.join(root,filename)
		#print ifile
		command = "xrdcp "+ifile+" root://cmseos.fnal.gov/"+output_folder+filename
		print command
		#os.system(command)
