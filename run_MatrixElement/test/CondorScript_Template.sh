#!/bin/sh
echo "================= Dumping Input files ===================="
python -c "import PSet; print '\n'.join(list(PSet.process.source.fileNames))"

echo "========================= pwd ============================"
pwd
echo "======================== ls ./ ==========================="
ls
echo "======== ls CMSSW_5_3_27/lib/slc6_amd64_gcc472 ==========="
ls CMSSW_5_3_27/lib/slc6_amd64_gcc472
echo "=== mv lib/lib*.so CMSSW_5_3_27/lib/slc6_amd64_gcc472 ===="
echo "mv lib/slc6_amd64_gcc472/lib*.so CMSSW_5_3_27/lib/slc6_amd64_gcc472"
mv lib/slc6_amd64_gcc472/lib*.so CMSSW_5_3_27/lib/slc6_amd64_gcc472
echo "======= cp TF_TTbarMG_*_00_24.txt CMSSW_5_3_27/src ======="
echo "cp TF_TTbarMG_*_00_24.txt CMSSW_5_3_27/src"
cp TF_TTbarMG_*_00_24.txt CMSSW_5_3_27/src
echo "======== ls CMSSW_5_3_27/lib/slc6_amd64_gcc472 ==========="
ls CMSSW_5_3_27/lib/slc6_amd64_gcc472
echo "=================== cd CMSSW_5_3_27 ======================"
echo "cd CMSSW_5_3_27"
cd CMSSW_5_3_27
echo "========================= scram =========================="
echo "eval \`scram runtime -sh\`"
eval `scram runtime -sh`
echo "========================= cd - ==========================="
echo "cd -"
cd -
echo "==================== Add Search Path ====================="
echo "export CMSSW_SEARCH_PATH \`pwd\`:$CMSSW_SEARCH_PATH"
export CMSSW_SEARCH_PATH=`pwd`:$CMSSW_SEARCH_PATH

echo "=========================================================="
echo "**********************************************************"
echo "Starting" 
JobNumber=$1
echo "JobNumber=$JobNumber"
NEvtsPerJob=${2#*=}
echo "NEvtsPerJob=$NEvtsPerJob"
TheJob=$JobNumber
echo "TheJob=$TheJob"
StartEvt=`expr $TheJob \* $NEvtsPerJob`
echo "StartEvt=$StartEvt"
./run_MatrixElement root://cmsxrootd.fnal.gov//store/user/aperloff/MatrixElement/Summer12ME8TeV/MEInput/WW.root WW.root jets2p $NEvtsPerJob $StartEvt 1 1 PS 0 MissingEventsFiles/microWW_missingEvents.txt
echo "Finished"
echo "**********************************************************"
echo "=========================================================="

echo "========================= pwd ============================"
pwd
echo "======================== ls ./ ==========================="
ls