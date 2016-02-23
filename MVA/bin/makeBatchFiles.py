#! /usr/bin/env python
import os, re, sys, subprocess
from optparse import OptionParser

def vararg_callback(option, opt_str, value, parser):
    assert value is None
    value = []

    def floatable(str):
    	try:
            float(str)
            return True
        except ValueError:
            return False

    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break
        # stop on -a, but not on -3 or -3.0
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break
        value.append(arg)

    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)


def are_eq(a, b):
    return set(a) == set(b) and len(a) == len(b)

def makeBatchFiles():
  #******************   template file  **********************************
  jdlFile = open('MVABatchJob_Template.jdl', 'r')
  shFile  = open('MVABatchJob_Template.sh', 'r')
  
  jdlFileContents = jdlFile.read(-1)
  shFileContents  = shFile.read(-1)
  print '--------------- TEMPLATE (jdl) :  -----------------'
  print jdlFileContents
  print '--------------- TEMPLATE (sh) :  -----------------'
  print shFileContents
  p0  = re.compile(r'METHOD')
  p1  = re.compile(r'JETBIN')
  p2  = re.compile(r'TAGCAT')
  p3  = re.compile(r'LEPTONCAT')
  p4  = re.compile(r'SIGNALS')
  p5  = re.compile(r'SIGNALSSHORT')
  p6  = re.compile(r'BACKGROUNDS')
  p7  = re.compile(r'EVENTPROBS')
  p8  = re.compile(r'EVENTPROBSSHORT')
  p9  = re.compile(r'KINVAR')
  p10 = re.compile(r'NKINVAR')
  p11 = re.compile(r'WORKINGDIRECTORY')

  #******************   definitions  **********************************
  script_prefix 		= 'MVABatchJob'
  JETBIN        		= options.jetbin
  EVENTPROB_FULL_LIST = ['0','1','2','3','4','5','6','7','8','9','10','11','12','19','54']
  SIGNALS_FULL_LIST   = ['ggH125','qqH125','WH_HToWW_M125','ZH_HToWW_M125','TTH_HToWW_M125']
  #SIGNALS_FULL_LIST   = ['ggH125','qqH125','WH_ZH_TTH_HToWW_M125']
  
  #********************  parsing  *************************************
  
  SIGNALSSHORT      = ''
  SIGNALSSTRING     = ''
  BACKGROUNDS       = ''
  BACKGROUNDSSTRING = ''
  EVENTPROBSSHORT   = ''
  EVENTPROBSSTRING  = ''
  KINVARSTRING      = ''
  METHODSSTRING     = ''
  
  last = len(options.signals) - 1
  for i, s in enumerate(options.signals) :
  	SIGNALSSHORT = SIGNALSSHORT + s
  	if i == last :
  		SIGNALSSTRING = SIGNALSSTRING + s
  	else :
  		SIGNALSSTRING = SIGNALSSTRING + s + ' '
  if (are_eq(SIGNALS_FULL_LIST,options.signals)) :
  	SIGNALSSHORT = 'HToWW'
  
  last = len(options.backgrounds) - 1
  for i, b in enumerate(options.backgrounds) :
  	BACKGROUNDS = BACKGROUNDS + b
  	if i == last :
  		BACKGROUNDSSTRING = BACKGROUNDSSTRING + b
  	else :
  		BACKGROUNDSSTRING = BACKGROUNDSSTRING + b + ' '
  
  last = len(options.eventprobs) - 1
  for i, e in enumerate(options.eventprobs) :
  	EVENTPROBSSHORT = EVENTPROBSSHORT + e
  	if i == last :
  		EVENTPROBSSTRING = EVENTPROBSSTRING + e
  	else :
  		EVENTPROBSSTRING = EVENTPROBSSTRING + e + ' '
  if(are_eq(EVENTPROB_FULL_LIST,options.eventprobs) or options.eventprobs[0] == "all") :
  	EVENTPROBSSHORT  = 'all'
  	EVENTPROBSSTRING = '\" \"'
  if(len(options.eventprobs)==0 or options.eventprobs[0] == "no" or options.eventprobs[0] == "-1") :
  	EVENTPROBSSHORT = 'no'
  	EVENTPROBSSTRING = "-1"
  
  last = len(options.kinvar) - 1
  for i, k in enumerate(options.kinvar) :
  	if i == last :
  		KINVARSTRING = KINVARSTRING + k
  	else :
  		KINVARSTRING = KINVARSTRING + k + ' '
  if(len(options.kinvar)==0) :
  	KINVARSTRING = '\" \"'
  NKINVAR = len(options.kinvar)
  
  last = len(options.methods) - 1
  for i, m in enumerate(options.methods) :
  	if i == last :
  		METHODSSTRING = METHODSSTRING + m
  	else :
  		METHODSSTRING = METHODSSTRING + m + ' '
  
  #*******************************************************************
  
  WORKINGDIRECTORY = '2015_07_17_TMVA_output_' + options.jetbin + '_' + options.tagcat + '_' + options.leptoncat
  WORKINGDIRECTORY = WORKINGDIRECTORY + '_' + SIGNALSSHORT + '_' + BACKGROUNDS + '_' + EVENTPROBSSHORT + 'EvtProbs_' + str(NKINVAR) + 'KinVar'
  
  if not os.path.isdir(WORKINGDIRECTORY):
  	os.makedirs(WORKINGDIRECTORY)
  	os.makedirs(WORKINGDIRECTORY+'/plots')
  
  k1 = p11.sub(WORKINGDIRECTORY, jdlFileContents)
  k1outfile = script_prefix + '.jdl'
  print '--------------------------------------'
  print 'Condor file for TMVA job : ' + WORKINGDIRECTORY
  print 'Written to ' + WORKINGDIRECTORY + "/" + k1outfile
  FILE = open(WORKINGDIRECTORY+"/"+k1outfile,"w")
  FILE.write(k1)
  
  k1 = p1.sub(options.jetbin,    shFileContents)
  k2 = p2.sub(options.tagcat,    k1)
  k3 = p3.sub(options.leptoncat, k2)
  k4 = p4.sub(SIGNALSSTRING,     k3)
  k5 = p6.sub(BACKGROUNDSSTRING, k4)
  k6 = p7.sub(EVENTPROBSSTRING,  k5)
  k7 = p9.sub(KINVARSTRING,      k6)
  k8 = p0.sub(METHODSSTRING,     k7)
  k8outfile = script_prefix + '.sh'
  print '--------------------------------------'
  print 'Script file for TMVA job : ' + WORKINGDIRECTORY
  print 'Written to ' + WORKINGDIRECTORY + "/" + k8outfile
  FILE = open(WORKINGDIRECTORY+"/"+k8outfile,"w")
  FILE.write(k8)

parser = OptionParser(version="%prog 1.0")

parser.add_option('-m','--methods', metavar='METHOD', action='callback',
                  dest='methods', callback=vararg_callback,
                  help='The methods to use in the TMVA training')
parser.add_option('-j','--jetbin', metavar='JETBIN', action='store',
                  dest='jetbin', type="choice", choices=["jets2","jets3","jets4"],
                  help='Input DEFS::JetBin')
parser.add_option('-t','--tagcat', metavar='TAGCAT', action='store',
                  dest='tagcat', type="choice", choices=["pretag","eq0tag","eq1tag","eq2tag","ge0tag","ge1tag","ge2tag"],
                  help='Input DEFS::TagCat')
parser.add_option('-l','--leptoncat', metavar='LEPTONCAT', action='store',
                  dest='leptoncat', type="choice", choices=["none","electron","muon","both"],
                  help='Input DEFS::LeptonCat')
parser.add_option('-s','--signals', metavar='signals', action='callback',
                  dest='signals', callback=vararg_callback,
                  help='The signal samples to use')
parser.add_option('-b','--backgrounds', metavar='backgrounds', action='callback',
                  dest='backgrounds', callback=vararg_callback,
                  help='The background samples to use')
parser.add_option('-e','--eventprobs', metavar='eventprobs', action='callback',
                  dest='eventprobs', callback=vararg_callback,
                  help='The event probabilities to train over. Use \'all\' to train over all event probabilities and \'-1\' to train over none of them')
parser.add_option('-k','--kinvar', metavar='kinvar', action='callback',
                  dest='kinvar', callback=vararg_callback,
                  help='The kinematic variables to train over.')

(options, args) = parser.parse_args()

makeBatchFiles()
