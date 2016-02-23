#!/bin/tcsh
setenv SCRAM_ARCH slc5_amd64_gcc462
setenv X509_USER_PROXY /uscms/home/aperloff/x509up_u43841
source /cvmfs/cms.cern.ch/cmsset_default.sh
set CONDORDIR=$PWD
set CMSSW_VERSION=CMSSW_5_3_2_patch4

scram p $CMSSW_VERSION

mkdir -p $CMSSW_VERSION/src/TAMUWW
mkdir -p $CMSSW_VERSION/src/TAMUWW/MVA
mkdir -p $CMSSW_VERSION/lib/${SCRAM_ARCH}
mkdir -p plots
mkdir -p weights

mv ConfigFiles $CMSSW_VERSION/src/TAMUWW
mv macros $CMSSW_VERSION/src/TAMUWW/MVA
mv lib*.so $CMSSW_VERSION/lib/slc5_amd64_gcc462
cp $CMSSW_VERSION/src/TAMUWW/MVA/macros/tmva_logo.gif .

cd CMSSW_5_3_2_patch4
scram b
eval `scram runtime -csh`
cd ..

setenv BASEPATH ${PWD}
echo "BASEPATH="$BASEPATH

set START_TIME=`/bin/date +%m_%d_%y_%H_%M_%S`

#./TMVAClassification_x -myMethodList BDT -lep both -custom false -signals ggH125 qqH125 WH125 WH_HToZZ_M125 ZH_HToZZ_M125 TTH_HToZZ_M125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 TTH_HToBB_M125 -backgrounds " " -train true -plot true -batch true
#./TMVAClassification_x -myMethodList BDT -lep electron -custom false -signals ggH125 -train false -plot true -ofile TMVA1 -batch true
#./TMVAClassification_x -myMethodList BDT -lep both -jetBin jets2 -tagcat pretag -custom false -signals ggH125 -train true -plot true -batch true



date

