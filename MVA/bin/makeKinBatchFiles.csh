#! /bin/tcsh

python makeBatchFiles.py --methods BDT --jetbin jets2 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs no --kinvar leptonPt lepMT jet1dRLep jet2dRLep ht Ptlnujj dRlepjj dPhiMETJet dPhiJetJet CosTheta_l CosTheta_WH
python makeBatchFiles.py --methods BDT --jetbin jets3 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs no --kinvar leptonPt leptonEtaCharge jet2dRLep jet3dRLep ht Mlnujj dRlepjj dPhiMETJet dEtaJetJet minDPhiLepJet CosTheta_l CosTheta_j CosTheta_WH
python makeBatchFiles.py --methods BDT --jetbin jets4 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs no --kinvar leptonEtaCharge jet2dRLep jet3dRLep ht Mlnujj dPhiMETJet dPhiMETLep
