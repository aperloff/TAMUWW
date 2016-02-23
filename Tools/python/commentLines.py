#! /usr/bin/env python
import os, sys, getopt, argparse, fnmatch, errno, subprocess, tempfile, logging, shutil
from subprocess import call

log = logging.getLogger('commentLines')

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

def string_contains(string_to_check, array_of_options):
    for option in array_of_options:
        if option in string_to_check:
            return (True, option)
    return (False, "")

def line_is_commented(line):
    if line[0]=='#':
        return True
    else:
        return False

def line_is_systematic(line):
    if "lnN" in line or " shape " in line:
        return True
    else:
        return False

def un_comment(input, output, comment, uncomment):
    line_counter = 0
    sys_counter = 0
    commented_sys = 0

    for line in input:
        if string_contains(line, comment)[0] and not line_is_commented(line):
            log.info(bcolors.OKGREEN+"Commenting the line containing "+string_contains(line, comment)[1]+bcolors.ENDC)
            output.write('#'+line)
            if line_is_systematic(line):
                commented_sys+=1
        elif string_contains(line, comment)[0] and line_is_commented(line):
            log.info(bcolors.BWARNING+"The line containing \""+string_contains(line, comment)[1]+"\" is already commented!"+bcolors.ENDC)
            output.write(line)
            if line_is_systematic(line):
                commented_sys+=1
        elif string_contains(line, uncomment)[0] and not line_is_commented(line):
            log.info(bcolors.BWARNING+"The line containing \""+string_contains(line, uncomment)[1]+"\" is already uncommented!"+bcolors.ENDC)
            output.write(line)
        elif string_contains(line, uncomment)[0] and line_is_commented(line):
            log.info(bcolors.OKBLUE+"Uncommenting the line containing "+string_contains(line, uncomment)[1]+bcolors.ENDC)
            output.write(line[1:])
        else:
            output.write(line)

        line_counter+=1
        if line_is_systematic(line):
            sys_counter+=1

    log.info(bcolors.UNDERLINE+"Total lines: "+str(line_counter)+bcolors.ENDC)
    log.info(bcolors.UNDERLINE+"Systematic lines (Total): "+str(sys_counter)+bcolors.ENDC)
    log.info(bcolors.UNDERLINE+"Systematic lines (Uncommented): "+str(sys_counter-commented_sys)+bcolors.ENDC)

def main(inputfilename, outputfilename, comment, uncomment):
    try:
        input = open(args.input, 'r')
    except IOError:
        print bcolors.BFAIL + "Can't open input file: %s" % inputfilename + bcolors.ENDC
        exit(1)
    try:
        output = open(args.output, 'w')
    except IOError:
        print bcolors.BFAIL + "Can't open output file: %s" % outputfilename + bcolors.ENDC
        exit(2)
    un_comment(input, output, comment, uncomment)

if __name__ == "__main__":
    #program name available through the %(prog)s command
    #can use prog="" in the ArgumentParser constructor
    #can use the type=int option to make the parameters integers
    #can use the action='append' option to make a list of options
    parser = argparse.ArgumentParser(description="This program will comment or uncomment lines which are specified by a portion of the line. It was originally intended to (un)comment datacards to turn on and off certain systematics.",
                                     epilog="And those are the options available. Deal with it.")
    parser.add_argument("input", help="Input datacard")
    parser.add_argument("output", help="Output datacard")
    parser.add_argument("-c", "--comment", help="The names of the systematics to turn off",
                        nargs='*', default=[])
    parser.add_argument("-u", "--uncomment", help="The names of the systematics to turn on",
                        nargs='*', default=[])
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('--verbose', action='store_true',
                        help='Print debug output.')
    args = parser.parse_args()

    logging.basicConfig(
        stream=sys.stderr,
        level=logging.INFO if not args.verbose else logging.DEBUG)

    if args.input == args.output:
        log.info("Modifying in place!  Backing up input file...")
        shutil.copy(args.input, args.input.replace('.txt', '.txt.bak'))
        args.output = args.output.replace('.txt', '.tmp.txt')

    log.info("(Un)commenting datacard. input: %s output: %s",
             args.input, args.output)
    main(args.input, args.output, args.comment, args.uncomment)
    log.info("Moving temprorary output to final destination")
    shutil.move(args.output, args.output.replace('.tmp.txt', '.txt'))
