#! /usr/bin/env python
import os, sys, getopt, argparse, fnmatch, errno, re, subprocess, tempfile, shlex, checkfiles
from glob import glob
from time import sleep
from string import translate,maketrans,punctuation 
from textwrap import dedent

'''
From http://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text
and https://bitbucket.org/ruamel/std.argparse
'''
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        # this is the RawTextHelpFormatter._split_lines
        if text.startswith('R|'):
            return text[2:].splitlines()  
        return argparse.HelpFormatter._split_lines(self, text, width)

def run_checks(input_path):
    # That the path to the output files exists
    stores = ("/store/user/","/store/group/")
    if any(x in input_path for x in stores) and ("/eos/uscms/" not in input_path):
        print "Using XROOTD\n"
    else:
        print "Using local files or the FUSE mount\n"
        if not os.path.exists(input_path):
            sys.exit()

def do_hadd(input_path, pattern, file_types):
    files = checkfiles.get_list_of_files(pattern=pattern, file_types=file_types, path=input_path)

    file_split = len(files)/900

    stores = ("/store/user/","/store/group/")
    if any(x in input_path for x in stores) and ("/eos/uscms/" not in input_path):
        files = [f.replace('/eos/uscms', 'root://cmseos.fnal.gov/') for f in files]
    else:
        files = [input_path+f for f in files]

    if DEBUG:
        print "Files after splitting and filtering:\n"
        print files
        print "length:",str(len(files))
        print "Number times to split:",file_split

    output_file = args.output_path+"/"+args.output_prefix
    
#    for i in range(1+file_split):
#        current_set_of_files = files[i*900:min(len(files),(i+1)*900)]
#        if DEBUG:
#            print len(current_set_of_files)
#        command = "nohup hadd " + output_file+str(i)+".root "
#        for f in current_set_of_files:
#            command += f+" "
#        if DEBUG:
#            print command
#        script=open(output_file+str(i)+".sh","w")
#        script.write("#! /usr/bin/sh\n\n")
#        script.write(command)
#        script.write("\n\necho \"HADD COMPLETE\"\n")
#        script.close()
#        out=open(output_file+str(i)+".stdout","w")
#        subprocess.Popen(["bash",output_file+str(i)+".sh"],stdout=out,stderr=subprocess.STDOUT)
#        #subprocess.Popen(["nohup","hadd",output_file+str(i)+".root",input_files],stdout=out,stderr=subprocess.STDOUT)

    # Create the nohup hadd processes broken into groups of 900 files
    child_filenames = []
    procs = []
    for i in range(1+file_split):
        current_set_of_files = files[i*900:min(len(files),(i+1)*900)]
        if DEBUG:
            print "The number of files currently being hadded is",len(current_set_of_files)
        child_filenames.append(output_file+str(i))
        command = "nohup hadd " + child_filenames[-1] + ".root "
        for f in current_set_of_files:
            command += f+" "
        if DEBUG:
            print "The current command is",command
        out=open(output_file+str(i)+".stdout","w")
        procs.append(subprocess.Popen(shlex.split(command), shell=False, stdout=out, stderr=subprocess.STDOUT))

    return_code_sum = 0
    for iproc in procs:
        return_code_sum+=iproc.wait()
    out.close()

    # This section of code deals with making the final ntuple from the child ROOT files
    out=open(output_file+".stdout","a")
    if return_code_sum==0:
        command = "nohup hadd " + output_file + ".root "
        for f in child_filenames:
            command += f+".root "
        if DEBUG:
            print "The current command is",command
        final_proc = subprocess.Popen(shlex.split(command), shell=False, stdout=out, stderr=subprocess.STDOUT)
        return_code_ntuple = final_proc.wait()
    else:
        out.write("ERROR::hadd.py One of the child hadd processes did not finish successfully. Cannnot hadd the child ROOT files into a single ntuple.\n")
        exit(-1)

    #This section of code deals with reporting the final result and cleanup
    if return_code_ntuple==0:
        out.write("\nHADD COMPLETE\n")
        out.write("Removing the child files ... \n")
        for ifile in child_filenames:
            out.write("\t"+ifile+".[root,stdout] ... ")
            os.remove(ifile+".root")
            os.remove(ifile+".stdout")
            out.write("DONE\n")
    else:
        out.write("ERROR::hadd.py There was a problem hadding the child ROOT files into a single ntuple.\n\tThe child jobs probably finished successfully (see the logs to be sure), but something when wrong with the final process.\n")
        exit(-2)
    out.close()

def main(input_path, pattern, file_types, ignore, output_path, output_prefix, total_number_of_jobs, debug):
    global DEBUG
    DEBUG = debug

    output_file = args.output_path+"/"+args.output_prefix
    with open(output_file+".stdout", 'w') as sys.stdout:
        run_checks(input_path)

        if total_number_of_jobs > 0:
            return_code = checkfiles.main(input_path=input_path, pattern=pattern, file_types=file_types, symlinks=False,
                                              ignore=ignore, total_number_of_jobs=total_number_of_jobs)
            if return_code!=0:
                print "ERROR checkfiles.py returns a non-zero exit code.\n\tThis means that there were either missing or duplicate files.\n\tFor safety reasons this code will discontinue without doing an hadd."
                exit(-1)

    do_hadd(input_path=input_path, pattern=pattern, file_types=file_types)

if __name__ == '__main__':
    #program name available through the %(prog)s command
    #can use prog="" in the ArgumentParser constructor
    #can use the type=int option to make the parameters integers
    #can use the action='append' option to make a list of options
    parser = argparse.ArgumentParser(description="""\
                                     Hadd the output of a bunch of CRAB files. Should work for both CRAB2 and CRAB3
                                     Example: [nohup] python hadd.py /store/user/aperloff/MatrixElement/SingleElectron/SingleEl_A_ReReco/ ./ test -j 362 -p PS_ [&]
                                     """,
                                     epilog="And those are the options available. Deal with it.\nNOTE: This program doesn't yet check the validity of the files (i.e. are the ROOT files corrupt).",
                                     formatter_class=SmartFormatter)
    group = parser.add_mutually_exclusive_group()
    parser.add_argument("input_path", help="The location of the files within the user/group store area or a local area.")
    parser.add_argument("output_path", help="The path to put the output ROOT files and script files (i.e. /uscms_data/d2/<username>/).")
    parser.add_argument("output_prefix", help="The prefix to append to the output files. This is usually the process name (i.e. DYJetsToLL_M-10to50_lowerSubLeadingJet or WlnuJets_M-50_lowerSubLeadingJet_HighPtLeptonKept050mu061el).")
    parser.add_argument("-d", "--debug", help="Shows some extra information in order to debug this program.",
                        action="store_true")
    parser.add_argument("-i", "--ignore", help="Patterns of files/folders to ignore. Implemented for checkfiles.py, but not fully implemented for hadd.py.",
                        nargs='+', type=str, default=())
    parser.add_argument("-j", "--total_number_of_jobs", help="This is the number of files that CRAB should have created and is corresponds to the highest file number.",
                        type=int, default=0)
    parser.add_argument("-p", "--pattern", help="Shared portion of the name of the files.",
                        nargs='*', default=["*"])
    parser.add_argument("-t", "--file_types", help="Specify the type (extension) of the files to check. By default this program only checks the .root output of a CRAB job.",
                        nargs='+', default=[".root"])
    group.add_argument("-v", "--verbose", help="Increase output verbosity of lcg-cp (-v) or srm (-debug) commands.",
                        action="store_true")
    parser.add_argument('--version', action='version', version='%(prog)s 1.1', help=dedent("""\
                        R|Prints the current version of the code.
                        V1.0 Manually edit the input and output locations
                        V1.1 More options, no need to manually edit this file,
                             checks the CRAB job first to make sure all files are present,
                             HADD the split jobs into a single file.
                        """))
    args = parser.parse_args()
    
    if(args.debug):
         print 'Number of arguments:', len(sys.argv), 'arguments.'
         print 'Argument List:', str(sys.argv)
         print "Argument ", args

    main(input_path=args.input_path, pattern=args.pattern, file_types=args.file_types, ignore=args.ignore, 
         output_path=args.output_path, output_prefix=args.output_prefix, total_number_of_jobs=args.total_number_of_jobs,
         debug=args.debug)
