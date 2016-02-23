#!/bin/tcsh
echo "Starting job on " `date` #Only to display the starting of production date
echo "Running on " `uname -a` #Only to display the machine where the job is running
echo "System release " `cat /etc/redhat-release` #And the system release
echo "CMSSW on Condor"

setenv SCRAM_ARCH slc5_amd64_gcc462
setenv X509_USER_PROXY /uscms/home/aperloff/x509up_u43841
source /cvmfs/cms.cern.ch/cmsset_default.csh
set CONDORDIR=$PWD

scram p CMSSW_5_3_2_patch5
cd CMSSW_5_3_2_patch5
mkdir -p src/TAMUWW
mkdir -p lib/${SCRAM_ARCH}

mv ../ConfigFiles src/TAMUWW
mv ../lib*.so lib/slc5_amd64_gcc462

scram b
eval `scram runtime -csh`
setenv BASEPATH ${PWD}
echo "BASEPATH="$BASEPATH

@ JobNumber = $argv[1]
echo "JobNumber=$JobNumber"
@ NEvtsPerJob = $argv[2]
echo "NEvtsPerJob=$NEvtsPerJob"
@ TheJob = 0 + $JobNumber
echo "TheJob=$TheJob"
@ StartEvt = $TheJob * $NEvtsPerJob
echo "StartEvt=$StartEvt"
@ EndEvt = ( $TheJob + 1 ) * $NEvtsPerJob
echo "EndEvt=$EndEvt"

../makeMicroNtuple_x -inputPaths INPUTPATH -addDir ADDDIR -outputPath $CONDORDIR -outputSuffix _$JobNumber -processes PROCESSNAMES -mergeNewEventNtuple /eos/uscms/store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/ -fillBDT true -debug false -useXROOTD false -startEntry $StartEvt -endEntry $EndEvt 
#-updateMicroNtuple UPDATEMICRONTUPLE

date
