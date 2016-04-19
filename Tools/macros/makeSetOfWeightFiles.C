#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"

#include <iostream>
#include <vector>
#include <string>

#include "TAMUWW/Tools/macros/MakeZJetToLLGenNeutrinoWeightFile.C"

using namespace std;

void makeSetOfWeightFiles() {
   for(int i=40; i<=70; i++) {
      MakeZJetsToLLGenNeutrinoWeightFile("/uscms_data/d2/aperloff/Summer12ME8TeV/CompareConvertedZllJetsToWJets/MakeGenNeutrinoPtWeights/",Form("0%i",i));
   }
}
