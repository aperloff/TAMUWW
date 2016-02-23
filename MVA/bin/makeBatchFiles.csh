#***2jets***
python makeBatchFiles.py --methods BDT --jetbin jets2 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs all --kinvar
python makeBatchFiles.py --methods BDT --jetbin jets2 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs all --kinvar leptonPt lepMT jet1dRLep jet2dRLep ht Ptlnujj dRlepjj dPhiMETJet dPhiJetJet CosTheta_l CosTheta_WH
python makeBatchFiles.py --methods BDT --jetbin jets2 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs no --kinvar leptonPt lepMT jet1dRLep jet2dRLep ht Ptlnujj dRlepjj dPhiMETJet dPhiJetJet CosTheta_l CosTheta_WH
python makeBatchFiles.py --methods BDT --jetbin jets2 --leptoncat both --tagcat eq1tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets TTbar --eventprobs all --kinvar
python makeBatchFiles.py --methods BDT --jetbin jets2 --leptoncat both --tagcat eq1tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets TTbar --eventprobs all --kinvar ht leptonPt jet2dRLep jet1dRLep CosTheta_l dEtaJetJet lepMT dPhiMETJet dPhiJetJet
python makeBatchFiles.py --methods BDT --jetbin jets2 --leptoncat both --tagcat eq1tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets TTbar --eventprobs no --kinvar ht leptonPt jet2dRLep jet1dRLep CosTheta_l dEtaJetJet lepMT dPhiMETJet dPhiJetJet

python makeBatchFiles.py --methods BDT --jetbin jets2 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs all --kinvar KinBDT
python makeBatchFiles.py --methods BDT --jetbin jets2 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs no --kinvar leptonPt lepMT jet1dRLep jet2dRLep ht Ptlnujj dRlepjj dPhiMETJet dPhiJetJet CosTheta_l CosTheta_WH MEBDT
python makeBatchFiles.py --methods BDT --jetbin jets2 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs no --kinvar MEBDT KinBDT

#***3jets***
python makeBatchFiles.py --methods BDT --jetbin jets3 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs all --kinvar
python makeBatchFiles.py --methods BDT --jetbin jets3 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs all --kinvar leptonPt leptonEtaCharge jet2dRLep jet3dRLep ht Mlnujj dRlepjj dPhiMETJet dEtaJetJet minDPhiLepJet CosTheta_l CosTheta_j CosTheta_WH
python makeBatchFiles.py --methods BDT --jetbin jets3 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs no --kinvar leptonPt leptonEtaCharge jet2dRLep jet3dRLep ht Mlnujj dRlepjj dPhiMETJet dEtaJetJet minDPhiLepJet CosTheta_l CosTheta_j CosTheta_WH
python makeBatchFiles.py --methods BDT --jetbin jets3 --leptoncat both --tagcat eq1tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets TTbar --eventprobs all --kinvar
python makeBatchFiles.py --methods BDT --jetbin jets3 --leptoncat both --tagcat eq1tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets TTbar --eventprobs all --kinvar dPhiMETJet CosTheta_j lepMT CosTheta_l jet2Pt ht jet3dRLep leptonPt dRlepjj CosTheta_WH JacobePeak
python makeBatchFiles.py --methods BDT --jetbin jets3 --leptoncat both --tagcat eq1tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets TTbar --eventprobs no --kinvar dPhiMETJet CosTheta_j lepMT CosTheta_l jet2Pt ht jet3dRLep leptonPt dRlepjj CosTheta_WH JacobePeak

python makeBatchFiles.py --methods BDT --jetbin jets3 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs all --kinvar KinBDT
python makeBatchFiles.py --methods BDT --jetbin jets3 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs no --kinvar leptonPt leptonEtaCharge jet2dRLep jet3dRLep ht Mlnujj dRlepjj dPhiMETJet dEtaJetJet minDPhiLepJet CosTheta_l CosTheta_j CosTheta_WH MEBDT
python makeBatchFiles.py --methods BDT --jetbin jets3 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs no --kinvar MEBDT KinBDT

#***4jets***
python makeBatchFiles.py --methods BDT --jetbin jets4 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs all --kinvar
python makeBatchFiles.py --methods BDT --jetbin jets4 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs all --kinvar leptonEtaCharge jet2dRLep jet3dRLep ht Mlnujj dPhiMETJet dPhiMETLep
python makeBatchFiles.py --methods BDT --jetbin jets4 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs no --kinvar leptonEtaCharge jet2dRLep jet3dRLep ht Mlnujj dPhiMETJet dPhiMETLep
python makeBatchFiles.py --methods BDT --jetbin jets4 --leptoncat both --tagcat eq1tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets TTbar --eventprobs all --kinvar
python makeBatchFiles.py --methods BDT --jetbin jets4 --leptoncat both --tagcat eq1tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets TTbar --eventprobs all --kinvar dPhiMETJet dRlepjj jet2dRLep lepMT sumJetEt leptonPt ht jet3dRLep JacobePeak jet4dRLep CosTheta_l nJets CosTheta_j
python makeBatchFiles.py --methods BDT --jetbin jets4 --leptoncat both --tagcat eq1tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets TTbar --eventprobs no --kinvar dPhiMETJet dRlepjj jet2dRLep lepMT sumJetEt leptonPt ht jet3dRLep JacobePeak jet4dRLep CosTheta_l nJets CosTheta_j

python makeBatchFiles.py --methods BDT --jetbin jets4 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs all --kinvar KinBDT
python makeBatchFiles.py --methods BDT --jetbin jets4 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs no --kinvar leptonEtaCharge jet2dRLep jet3dRLep ht Mlnujj dPhiMETJet dPhiMETLep MEBDT
python makeBatchFiles.py --methods BDT --jetbin jets4 --leptoncat both --tagcat eq0tag --signals ggH125 qqH125 WH_HToWW_M125 ZH_HToWW_M125 TTH_HToWW_M125 --backgrounds WJets --eventprobs no --kinvar MEBDT KinBDT
