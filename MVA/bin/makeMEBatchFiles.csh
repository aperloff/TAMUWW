#! /bin/tcsh

python makeBatchFiles.py --methods BDT --jetbin jets2 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs all --kinvar
python makeBatchFiles.py --methods BDT --jetbin jets3 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs all --kinvar
python makeBatchFiles.py --methods BDT --jetbin jets4 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs all --kinvar
