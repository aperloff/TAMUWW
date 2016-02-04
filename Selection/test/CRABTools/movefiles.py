#! /usr/bin/env python
import os, sys, getopt, argparse, fnmatch, errno, subprocess, tempfile
from subprocess import call

class Error(EnvironmentError):
    pass

siteDBDict = {
#   alias: ('host', local_file_system, 'server', 'path_to_user_area','local_home','xrootd_path')    
    'T1_US_FNAL_EOS'   : (False, 'cmseos.fnal.gov'    , '/srm/v2/server?SFN=' ,  '/eos/uscms/store/user/',                 '', 'root://cmseos.fnal.gov//store/user/'),
	'T1_US_FNAL_LPC'   : (True, ''                   , ''                    ,                        '', os.environ['HOME'], ''),
    'T3_US_TAMU'       : (False, 'srm.brazos.tamu.edu', '/srm/v2/server?SFN=' , '/fdata/hepx/store/user/',                 '', '')
}

#program name available through the %(prog)s command
#can use prog="" in the ArgumentParser constructor
#can use the type=int option to make the parameters integers
#can use the action='append' option to make a list of options
parser = argparse.ArgumentParser(description="Transfer files from one location on the OSG to another",
                                 epilog="And those are the options available. Deal with it.")
group = parser.add_mutually_exclusive_group()
parser.add_argument("STARTserver", help="The name of the server where the files are initially located",
                    choices=siteDBDict.keys())
parser.add_argument("STARTpath", help="The location of the files within the users store or resilient areas")
parser.add_argument("ENDserver", help="The name of the server to which the files should be copied",
                    choices=siteDBDict.keys())
parser.add_argument("ENDpath", help="The location of the files within the users fdata area")
parser.add_argument("-d","--debug", help="Shows some extra information in order to debug this program",
					action="store_true")
parser.add_argument("-diff","--diff",
                    help="Tells the program to do a diff between the two directories and only copy the missing files. Only works for two local directories.",
                    action="store_true")
parser.add_argument("-i","--ignore", help="Patterns of files/folders to ignore",
                    nargs='+', type=str, default=())
protocol = parser.add_mutually_exclusive_group()
protocol.add_argument("-lcg", "--lcg", help="Tells the program to use the lcg-cp command rather than srmcp",
                    action="store_true")
protocol.add_argument("-x", "--xrdcp", help="Tells the program to use the xrdcp command rather than the lcg or srm protocols",
                      action="store_true")
parser.add_argument("-pnfs", "--pnfs", help="dCache location of original files", choices=["store","resilient"],
                    default="store")
group.add_argument("-q", "--quiet", help="decrease output verbosity to minimal amount",
                   action="store_true")
parser.add_argument("-r", "--recursive", help="Recursively copies directories and files",
                    action="store_true")
parser.add_argument("-s", "--sample", help="Shared portion of the name of the files to be copied", nargs='+',
                    default=["*"])
parser.add_argument("-srt","--sendreceive_timeout", help="Sets the send/recieve timeout for the lcg command",
					default="20")
parser.add_argument("-str","--streams", help="The number of transfer streams", default="15")
parser.add_argument("-su", "--start_user", help="username of the person transfering the files",
                    default=os.environ['USER'])
parser.add_argument("-eu", "--end_user", help="username of the person transfering the files",
                    default=os.environ['USER'])
group.add_argument("-v", "--verbose", help="Increase output verbosity of lcg-cp (-v) or srm (-debug) commands",
                    action="store_true")
parser.add_argument('--version', action='version', version='%(prog)s 1.7')
args = parser.parse_args()

if(args.debug):
     print 'Number of arguments:', len(sys.argv), 'arguments.'
     print 'Argument List:', str(sys.argv)
     print "Argument ", args

SAMPLE = args.sample
START_USER = args.start_user
END_USER = args.end_user
PNFS = args.pnfs
LCG = args.lcg
XRDCP = args.xrdcp
SRT = args.sendreceive_timeout
START = args.STARTserver
STARTpath = args.STARTpath
END = args.ENDserver
ENDpath = args.ENDpath
DIFF = args.diff
VERBOSE = args.verbose
QUIET = args.quiet
RECURSIVE = args.recursive
IGNORE = tuple(args.ignore)
DEPTH = 1
STREAMS = args.streams
scommand = ""
ecommand = ""

def run_checks():
    global RECURSIVE
    global END

    if RECURSIVE :
        global DEPTH
        DEPTH = 9999

	if os.system("voms-proxy-info")!=0 :
		print "WARNING::You must have a valid proxy for this script to work.\nRunning \"voms-proxy-init -voms cms\"...\n"
		call("voms-proxy-init -voms cms -valid 192:00", shell=True)
		if os.system("voms-proxy-info")!=0 :
			print "ERROR::Sorry, but I still could not find your proxy.\nWithout a valid proxy, this program will fail spectacularly.\nThe program will now exit." 
			sys.exit()

def init_commands():
    if(LCG):
        if(VERBOSE):
            scommand = "lcg-cp -v -b -n "+STREAMS+" --sendreceive-timeout "+SRT+" --srm-timeout 60 -D srmv2"
        else:
            scommand = "lcg-cp -b -n "+STREAMS+" --sendreceive-timeout "+SRT+" --srm-timeout 60 -D srmv2"
    elif (XRDCP):
        if(VERBOSE):
            scommand = "xrdcp -v"
        elif(QUIET):
            scommand = "xrdcp -s"
        else:
            scommand = "xrdcp"
    else:
        if(VERBOSE):
            scommand = "srmcp -2 -pushmode=true -debug=true"
        else:
            scommand = "srmcp -2 -pushmode=true"

    if XRDCP and siteDBDict[START][0]:
        scommand += " \""+STARTpath+"/"
    elif siteDBDict[START][0]:
        scommand += " \"file:////"+siteDBDict[START][2]+siteDBDict[START][3]+"/"+STARTpath+"/"
    elif(XRDCP):
        scommand += " \""+siteDBDict[START][5]+"/"+START_USER+"/"+STARTpath+"/"
    else:
        scommand += " \"srm://"+siteDBDict[START][1]+":8443"+siteDBDict[START][2]+siteDBDict[START][3]+"/"+START_USER+"/"+STARTpath+"/"

    if(siteDBDict[END][0]):
        ecommand = "\"file:////"+siteDBDict[END][1]+siteDBDict[END][2]+siteDBDict[END][3]+"/"+ENDpath+"/"
    elif(XRDCP):
        ecommand = "\""+siteDBDict[END][5]+"/"+END_USER+"/"+ENDpath+"/"
    else:
        ecommand = "\"srm://"+siteDBDict[END][1]+":8443"+siteDBDict[END][2]+siteDBDict[END][3]+"/"+END_USER+"/"+ENDpath+"/"

    return (scommand,ecommand)


def make_directory(END, path):
    if END.find("FNAL") > 0 :
        if not os.path.exists(path) :
            print "making directory "+path
            os.makedirs(path)
    else:
        if LCG :
            cmd = "lcg-ls -v -b -D srmv2 \"srm://"+siteDBDict[END][1]+":8443"+siteDBDict[END][2]+path+"\""
        else:
            cmd = "srmls -2 -pushmode=true \"srm://"+siteDBDict[END][1]+":8443"+siteDBDict[END][2]+path+"\""
        if os.system(cmd)!=0 :
            print "making directory "+path
            cmd = "srmmkdir -2 -pushmode=true \"srm://"+siteDBDict[END][1]+":8443"+siteDBDict[END][2]+path+"\""
            os.system(cmd)

def get_list_of_files(SAMPLE, path):
    FILES_UNFILTERED = []
    if START.find("FNAL") > 0 and os.environ.get('HOSTNAME',"not found").find("fnal.gov") > 0 :
        FILES_UNFILTERED = os.listdir(path)
    else:
        #options = '-2 -pushmode=true "srm://'+siteDBDict[START][0]+':8443'+siteDBDict[START][1]+path+'"'
		#proc = subprocess.Popen(['srmls',options], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
        options = 'lcg-ls "srm://'+siteDBDict[START][1]+':8443'+siteDBDict[START][2]+path+'"'
        proc = subprocess.Popen(options, shell=True, stdout=subprocess.PIPE).communicate()[0]
		# if output is too long
        #proc = subprocess.Popen(['srmls',options], stdout=tempfile.TemporaryFile()).communicate()[0]
        FILES_UNFILTERED = proc.splitlines()
        FILES_UNFILTERED = [x.strip(' ') for x in FILES_UNFILTERED]
        FILES_UNFILTERED = [x.strip('\n') for x in FILES_UNFILTERED]
	FILES_UNFILTERED = [x.replace("//", "/") for x in FILES_UNFILTERED]
	#FILES_UNFILTERED.remove('Picked up _JAVA_OPTIONS: -Xmx1024m')
	FILES_UNFILTERED = [x.split(' ',1)[1] if x.split(' ',1)[0].isdigit() else x for x in FILES_UNFILTERED]
	FILES_UNFILTERED = [x.replace((siteDBDict[START][3]+"/"+START_USER+"/"+STARTpath+"/").replace("//","/"), "") for x in FILES_UNFILTERED]
	FILES_UNFILTERED = [x for x in FILES_UNFILTERED if x != '']
    FILES = []
    #if(len(SAMPLE)==1):
    #    SAMPLE = ["*"]
    for ss in SAMPLE:
        FILES += fnmatch.filter(FILES_UNFILTERED, '*'+ss+'*')
    return FILES

def filter_list_of_files(SAMPLE, FILES_UNFILTERED):
    FILES = []
    if(len(SAMPLE)==1):
        SAMPLE = ["*"]
    for ss in SAMPLE:
        FILES += fnmatch.filter(FILES_UNFILTERED, '*'+ss+'*')
    return FILES

def ignore_patterns(patterns):
    """Function that can be used as copytree() ignore parameter.

    Patterns is a sequence of glob-style patterns
    that are used to exclude files"""
    def _ignore_patterns(path, names):
        ignored_names = []
        for pattern in patterns:
            if len(fnmatch.filter([path],'*'+pattern+'*')) > 0:
                ignored_names.extend(names)
            #names_extended = [os.path.join(path,n) for n in names]
            #ignored_names.extend(fnmatch.filter(names_extended, '*'+pattern+'*'))
            else:
                ignored_names.extend(fnmatch.filter(names, '*'+pattern+'*'))
        return set(ignored_names)
    return _ignore_patterns

def copytree(src, dst, DEPTH, CURRENTdepth, symlinks=False, ignore=None):
    global args

    if args.debug :
        print "CURRENTdepth="+str(CURRENTdepth)
    if CURRENTdepth >= DEPTH :
        return
    
    make_directory(END,dst)

    global SAMPLE
    FILES = get_list_of_files(SAMPLE, src)
    if DIFF :
        FILES1 = FILES
        if END.find("FNAL") > 0 and os.environ.get('HOSTNAME',"not found").find("fnal.gov") > 0 :
            FILES2 = os.listdir(dst)
        else:
            #options = '-2 -pushmode=true "srm://'+siteDBDict[END][0]+':8443'+siteDBDict[END][1]+dst+'"'
            #proc = subprocess.Popen(['srmls',options], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            options = '"srm://'+siteDBDict[END][1]+':8443'+siteDBDict[END][2]+dst+'"'
            proc = subprocess.Popen(['lcg-ls',options], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0]
            proc_splitlines = proc.splitlines()
            FILES2 = []
            for f in proc_splitlines :
                if len(f.split()) > 0 and f.split()[0].isdigit() :
                    FILES2.append(os.path.basename(f.split()[1]))
        FILES_DIFF = [f for f in FILES1 if f not in FILES2 or os.path.isdir(dst+"/"+f)]
        print "Adding an additional "+str(len(FILES_DIFF))+" files/folders"
        FILES = FILES_DIFF
    
    if ignore is not None:
        ignored_names = ignore(src, FILES)
    else:
        ignored_names = set()

    errors = []
    for FILE in FILES:
        if FILE in ignored_names:
            continue
        srcname = os.path.join(src, FILE)
        dstname = os.path.join(dst, FILE)
        try:
            if symlinks and os.path.islink(srcname):
                linkto = os.readlink(srcname)
                os.symlink(linkto,dstname)
            elif os.path.isdir(srcname):
                copytree(srcname, dstname, DEPTH, CURRENTdepth+1, symlinks, ignore)
            else:
                global TOPsrc
                srel = os.path.relpath(src,TOPsrc)+"/"
                scommand, ecommand = init_commands()
                print "Copying file: "+FILE
                command = scommand+srel+FILE+"\" "+ecommand+srel+FILE+"\""
                print command
                os.system(command)
                
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except Error, err:
            errors.extend(err.args[0])
        except EnvironmentError, why:
            errors.append((srcname, dstname, str(why)))
    if errors:
        raise Error, errors
    
run_checks()

TOPsrc = ""
if(siteDBDict[START][0]):
	TOPsrc = STARTpath+"/"
else:
	TOPsrc = siteDBDict[START][4]+"/"+siteDBDict[START][3]+"/"+START_USER+"/"+STARTpath+"/"

TOPdst = siteDBDict[END][4]+"/"+siteDBDict[END][3]+"/"+END_USER+"/"+ENDpath+"/"

if args.debug:
	print "DEPTH="+str(DEPTH)
#copytree(TOPsrc,TOPdst, DEPTH, 0, False, ignore_patterns('DYJets','WJets','TTJets','STopT_T','ggH1?0'))
copytree(TOPsrc,TOPdst, DEPTH, 0, False, ignore_patterns(IGNORE))
