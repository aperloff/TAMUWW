#! /usr/bin/env python
import os, sys, getopt, argparse, fnmatch, errno
from subprocess import call
from string import translate,maketrans,punctuation 

class Error(EnvironmentError):
    pass

#program name available through the %(prog)s command
#can use prog="" in the ArgumentParser constructor
#can use the type=int option to make the parameters integers
#can use the action='append' option to make a list of options
parser = argparse.ArgumentParser(description="Check the output or res folder of a CRAB job an make sure that all of the files, and only one version of a file, exist.",
                                 epilog="And those are the options available. Deal with it.\nNOTE: This program doesn't yet check the validity of the files (i.e. are the ROOT files corrupt).")
group = parser.add_mutually_exclusive_group()
parser.add_argument("input_path", help="The location of the files within the users store or resilient areas")
parser.add_argument("total_number_of_jobs", help="This is the number of files that CRAB should have created and is corresponds to the highest file number",
                    type=int)
parser.add_argument("-d","--debug", help="Shows some extra information in order to debug this program",
					action="store_true")
parser.add_argument("-i","--ignore", help="Patterns of files/folders to ignore",
                    nargs='+', type=str, default=())
parser.add_argument("-p", "--pattern", help="Shared portion of the name of the files",
                    nargs='*', default=["*"])
parser.add_argument("-t", "--file_types", help="Specify the type (extension) of the files to check. By default this program only checks the .root output of a CRAB job.",
                    nargs='+', default=[".root"])
group.add_argument("-v", "--verbose", help="Increase output verbosity of lcg-cp (-v) or srm (-debug) commands",
                    action="store_true")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
args = parser.parse_args()

if(args.debug):
     print 'Number of arguments:', len(sys.argv), 'arguments.'
     print 'Argument List:', str(sys.argv)
     print "Argument ", args

input_path = args.input_path
total_number_of_jobs = args.total_number_of_jobs
pattern = args.pattern
file_types = args.file_types
VERBOSE = args.verbose
IGNORE = tuple(args.ignore)

def run_checks():
    # That the path to the output files exists
    global input_path
    if not os.path.exists(input_path):
        sys.exit()

def get_list_of_files(pattern, file_types, path):
    FILES_UNFILTERED = os.listdir(path)
    FILES = []
    if(len(pattern)==1):
        pattern = ["*"]
    for pp in pattern:
        for ft in file_types:
            FILES += fnmatch.filter(FILES_UNFILTERED, '*'+pp+'*'+ft)
    return FILES

def ignore_patterns(patterns):
    """Function that can be used as checkfiles() ignore parameter.

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

def get_file_number(filename):
    T = maketrans(punctuation, ' '*len(punctuation))
    pieces = translate(filename, T).split()
    if (pieces[0] == "crab"):
        return int(pieces[2])
    elif (pieces[0] == "Watchdog"):
        return int(pieces[1])
    elif (pieces[0] == "CMSSW"):
        return int(pieces[1])
    else:
        return int(filename.split('_')[len(filename.split('_'))-3])

def get_file_submission(filename):
    T = maketrans(punctuation, ' '*len(punctuation))
    pieces = translate(filename, T).split()
    if (pieces[-1:]=='.root'):
        return int(filename.split('_')[len(filename.split('_'))-2])
    else:
        return -9999

def print_header(header):
    print header + "\n" + "=" * len(header)
    
def print_list(list_to_print, header):
    print_header(header)
    print str(list_to_print) + "\n"
    
def print_dictionary(dict_to_print, header):
    print_header(header)
    print str(dict_to_print) + "\n"

def checkfiles(src, pattern, file_types, symlinks=False, ignore=None):

    print "\nChecking files in " + str(src) + " with the pattern(s) " + '[%s]' % ', '.join(map(str, pattern)) + " and the file type(s) " + '[%s]' % ', '.join(map(str, file_types)) + ", but ignoring the files with the pattern(s) " + '[%s]' % ', '.join(map(str, IGNORE)) + "\n"
    
    FILES = get_list_of_files(pattern, file_types, src)
    FILES_FOUND = {}
    FILES_WITH_DUPLICATES = []
    FILES_MISSING = []
    
    if ignore is not None:
        ignored_names = ignore(src, FILES)
    else:
        ignored_names = set()
        
    errors = []
    for FILE in FILES:
        if FILE in ignored_names:
            continue
        srcname = os.path.join(src, FILE)
        try:
            if not symlinks and os.path.islink(srcname):
                continue;
            else:
                if get_file_number(FILE) in FILES_FOUND:
                    FILES_FOUND[get_file_number(FILE)].append(FILE)
                    FILES_WITH_DUPLICATES.append(get_file_number(FILE))
                else:
                    FILES_FOUND[get_file_number(FILE)] = [FILE]
                
        # catch the Error from checkfiles so that we can continue with other files
        except Error, err:
            errors.extend(err.args[0])
        except EnvironmentError, why:
            errors.append((srcname, str(why)))
    if errors:
        raise Error, errors

    global total_number_of_jobs
    for i in range(1,total_number_of_jobs+1):
        if i not in FILES_FOUND:
            FILES_MISSING.append(i)
        else:
            continue;    

    print_dictionary(FILES_FOUND,"FILES FOUND ("+str(len(FILES_FOUND))+")")
    print_list(FILES_WITH_DUPLICATES,"FILES WITH DUPLICATES ("+str(len(FILES_WITH_DUPLICATES))+")")
    print_list(FILES_MISSING,"MISSING FILES ("+str(len(FILES_MISSING))+")")
    
run_checks()

checkfiles(input_path, pattern, file_types, False, ignore_patterns(IGNORE))
