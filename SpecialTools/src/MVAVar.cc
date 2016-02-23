#include "TAMUWW/SpecialTools/interface/MVAVar.hh"

using std::cout;
using std::endl;
using std::map;

// ----------------------------------------------------------------------------
// Default Constructor
MVAVar::MVAVar() {
   name = "";
   definition = "";
   unit = "";
   type = 'I';
   value = -9999;
   use = false;
   isSpectator = false;
   index = -9999;
   maxIndex = -9999;
}

// ----------------------------------------------------------------------------
// Constructor
MVAVar::MVAVar(TString n, TString def, TString un, char t,
               float v, bool u, bool s) {
   name = n;
   definition = def;
   unit = un;
   type = t;
   value = v;
   use = u;
   isSpectator = s;
   index = -9999;
   maxIndex = -9999;
}

// ----------------------------------------------------------------------------
MVAVar::~MVAVar() {

}

// ----------------------------------------------------------------------------
// A routine that returns the TSize
int MVAVar::getTSize() {
   //int size(0);
   //TFile* ftemp = TFile::Open(ifilePath + "micro" + ifilesSignal[0] + "_EPDv01.root");
   //TTree* ttemp = (TTree*)ftemp->Get("METree");
   //m_tSize = ttemp->GetBranch("m_tSize")->GetLeaf("m_tSize")->GetValue(0);
   //ttemp->SetBranchAddress("m_tSize",&m_tSize);
   //ttemp->GetEntry(0);
   //return m_tSize;
   return TSize;
   //return Micro::nEventProb;
}

// ----------------------------------------------------------------------------
// A routine that sets a map of MVA variables
void MVAVar::setVarMap(map<TString,MVAVar> & varMap) {
   for(int i=0; i<62; i++) {
      varMap[Form("EventProb%i",i)]       = MVAVar(Form("EventProb%i",i),Form("eventProb[%i]",i));
      varMap[Form("EventMaxProb%i",i)]    = MVAVar(Form("EventMaxProb%i",i),Form("eventMaxProb[%i]",i));
      varMap[Form("logEventProb%i",i)]    = MVAVar(Form("logEventProb%i",i),Form("TMath::Log(eventProb[%i])",i));
      varMap[Form("logEventMaxProb%i",i)] = MVAVar(Form("logEventMaxProb%i",i),Form("TMath::Log(eventMaxProb[%i])",i));

      //varMap[Form("eventProb%i",i)]       = Form("eventMaxProb%i := eventProb[%i]",i,i);
      //varMap[Form("eventMaxProb%i",i)]    = Form("eventMaxProb%i := eventMaxProb[%i]",i,i);
      //varMap[Form("logEventProb%i",i)]    = Form("logEventMaxProb%i := TMath::Log(eventProb[%i])",i,i);
      //varMap[Form("logEventMaxProb%i",i)] = Form("logEventMaxProb%i := TMath::Log(eventMaxProb[%i])",i,i);
   }

   varMap["leptonPt"]        = MVAVar(Form("leptonPt"),Form("lLV[0].Pt()"),"GeV");
   varMap["leptonEtaCharge"] = MVAVar(Form("leptonEtaCharge"),Form("lLV[0].Eta()*lLV[0].lQ"),"#etaCharge");
   varMap["lepMT"]           = MVAVar(Form("lepMT"),Form("TMath::Sign(TMath::Sqrt(TMath::Abs(TMath::Power(lLV[0].Pt()+METLV[0].Pt(),2)-(TMath::Power(lLV[0].Px()+METLV[0].Px(),2)+TMath::Power(lLV[0].Py()+METLV[0].Py(),2)+TMath::Power(0,2)))),TMath::Power(lLV[0].Pt()+METLV[0].Pt(),2)-(TMath::Power(lLV[0].Px()+METLV[0].Px(),2)+TMath::Power(lLV[0].Py()+METLV[0].Py(),2)+TMath::Power(0,2)))"),"GeV");
   varMap["lepMet"]          = MVAVar(Form("lepMet"),Form("lepMet"),"GeV");
   varMap["jet1Pt"]          = MVAVar(Form("jet1Pt"),Form("jLV[0].Pt()"),"GeV");
   varMap["jet2Pt"]          = MVAVar(Form("jet2Pt"),Form("jLV[1].Pt()"),"GeV");
   varMap["jet1dRLep"]       = MVAVar(Form("jet1dRLep"),Form("jLV[0].DRlj"),"dR");
   varMap["jet2dRLep"]       = MVAVar(Form("jet2dRLep"),Form("jLV[1].DRlj"),"dR");
   varMap["jet3dRLep"]       = MVAVar(Form("jet3dRLep"),Form("jLV[2].DRlj"),"dR");
   varMap["jet4dRLep"]       = MVAVar(Form("jet4dRLep"),Form("jLV[3].DRlj"),"dR");
   varMap["nBTagsCSVLep"]    = MVAVar(Form("nBTagsCSVLep"),Form("Sum$(jLV.jBtagCSV)"),"nBTags");
   varMap["jet1discrCSV"]    = MVAVar(Form("jet1discrCSV"),Form("jLV[0].jBtagDiscriminatorCSV"),"CSV Discr");
   varMap["ht"]              = MVAVar(Form("ht_lep"),Form("lLV[0].Pt()+Sum$(jLV.Et())"),"GeV");
   varMap["htLep"]           = MVAVar(Form("ht_lep"),Form("lLV[0].Pt()+METLV[0].Pt()"),"GeV");
   varMap["htAll"]           = MVAVar(Form("ht_all"),Form("lLV[0].Pt()+METLV[0].Pt()+Sum$(jLV.Et())"),"GeV");
   varMap["Mlnujj"]          = MVAVar(Form("Mlnujj_lep"),Form("TMath::Sign(TMath::Sqrt(TMath::Abs(TMath::Power(jLV[0].E()+jLV[1].E()+lLV[0].E()+METLV[0].E(),2)-(TMath::Power(jLV[0].Px()+jLV[1].Px()+lLV[0].Px()+METLV[0].Px(),2)+TMath::Power(jLV[0].Py()+jLV[1].Py()+lLV[0].Py()+METLV[0].Py(),2)+TMath::Power(jLV[0].Pz()+jLV[1].Pz()+lLV[0].Pz()+METLV[0].Pz(),2)))),TMath::Power(jLV[0].E()+jLV[1].E()+lLV[0].E()+METLV[0].E(),2)-(TMath::Power(jLV[0].Px()+jLV[1].Px()+lLV[0].Px()+METLV[0].Px(),2)+TMath::Power(jLV[0].Py()+jLV[1].Py()+lLV[0].Py()+METLV[0].Py(),2)+TMath::Power(jLV[0].Pz()+jLV[1].Pz()+lLV[0].Pz()+METLV[0].Pz(),2)))"),"GeV");
   varMap["Mjj"]             = MVAVar(Form("Mjj"),Form("Mjj"),"GeV");
   varMap["Ptlnujj"]         = MVAVar(Form("PtlnujjLep"),Form("TMath::Sqrt(TMath::Power(jLV[0].Px()+jLV[1].Px()+lLV[0].Px()+METLV[0].Px(),2)+TMath::Power(jLV[0].Py()+jLV[1].Py()+lLV[0].Py()+METLV[0].Py(),2))"),"GeV");
   varMap["dRlepjj"]         = MVAVar(Form("dRLepJet"),Form("dRlepjj"),"dR");
   varMap["dPhiMETJet"]      = MVAVar(Form("dPhiMetJet_lep"),Form("dPhiMETJet"),"d#Phi(met,j)");
   varMap["dPhiMETLep"]      = MVAVar(Form("dPhiMetLep_lep"),Form("TMath::Abs(METLV[0].Phi()-lLV[0].Phi())"),"d#Phi(met,l)");
   varMap["dEtaJetJet"]      = MVAVar(Form("dEtaJetJet_lep"),Form("TMath::Abs(jLV[0].Eta()-jLV[1].Eta())"),"d#eta(j,j)");
   varMap["dPhiJetJet"]      = MVAVar(Form("dPhiJJ_lep"),Form("dPhiJetJet"),"d#Phi(j,j)");
   varMap["minDPhiLepJet"]   = MVAVar(Form("minDPhiLep"),Form("minDPhiLepJet"),"min d#Phi");
   varMap["sumJetEt"]        = MVAVar(Form("sumJetEt"),Form("Sum$(jLV.Et())"),"#Sigma(Et)");
   varMap["CosTheta_l"]      = MVAVar(Form("CosThetaL_lep"),Form("CosTheta_l"),"cos(#theta)");
   varMap["CosTheta_j"]      = MVAVar(Form("CosThetaJ_lep"),Form("CosTheta_j"),"cos(#theta)");
   varMap["CosTheta_WH"]     = MVAVar(Form("CosThetaWH_lep"),Form("CosTheta_WH"),"cos(#theta)");
   varMap["JacobePeak"]      = MVAVar(Form("JacobePeak"),Form("JacobePeak"),"JPeak");
   varMap["CosdPhiWW"]       = MVAVar(Form("CosdPhiWW"),Form("CosdPhiWW_lep"),"cos(d#Phi)");
   varMap["nJets"]           = MVAVar(Form("nJetsLep"),Form("TMath::Min(@jLV.size(),10)"),"nJets");
   varMap["MEBDT"]           = MVAVar(Form("MEBDT"),Form("MEBDT"),"MEBDT",'F');
   varMap["KinBDT"]          = MVAVar(Form("KinBDT"),Form("KinBDT"),"KinBDT",'F');
   varMap["run"]             = MVAVar(Form("run"),Form("run"),"run",'F',-9999,false,true);
   varMap["lumi"]            = MVAVar(Form("lumi"),Form("lumi"),"lumi",'F',-9999,false,true);
   varMap["event"]           = MVAVar(Form("event"),Form("event"),"event",'F',-9999,false,true);
   varMap["csvWeight_lep"]   = MVAVar(Form("csvWeight_lep"),Form("csvWeight_lep"),"Weight",'F',-9999,false,true);
   varMap["ttbarWeight"]     = MVAVar(Form("ttbarWeight"),Form("ttbarWeight"),"Weight",'F',-9999,false,true);
   varMap["puWeight_lep"]    = MVAVar(Form("puWeight_lep"),Form("puWeight_lep"),"Weight",'F',-9999,false,true);
}

// ----------------------------------------------------------------------------
// A routine that returns an iterator to the first variable set to be used
// There is an option to check for first spectator or regular variable
map<TString,MVAVar>::iterator MVAVar::getFirstUsed(std::map<TString,MVAVar> & varMap, bool spectator) {
   map<TString,MVAVar>::iterator it;
   for(it=varMap.begin(); it!=varMap.end(); ++it) {
      if(it->second.use) {
         if(spectator && it->second.isSpectator) return it;
         if(!spectator && !it->second.isSpectator) return it;
      }
   }
   if(!spectator) {
      cout << "ERROR::MVAVar::getFirstUsed There are no MVA variables that are set to be used in this spectator category (" << spectator << ")" << endl;
      assert(it!=varMap.end());
   }
   return varMap.end();
}

// ----------------------------------------------------------------------------
// A routine that sets the variables to use based on the input from a table
void MVAVar::setUseFromTable(std::map<TString,MVAVar> & varMap, TString tableRow) {
   vector<TString> MVAV;
   vector<TString> MVAS;
   //cout << "\tvarMap.size() (1)= " << varMap.size() << endl; 
   DefaultValues::getMVAVar(tableRow,MVAV,MVAS);

   if (varMap.size()==0) {
      cout << "WARNING::MCAVar::setUseFromTable you should always set your varMap before calling this function." << endl;
      MVAVar::setVarMap(varMap);
   }
   for(unsigned int iMVA=0; iMVA<MVAV.size(); iMVA++) {
      varMap[MVAV[iMVA]].use = true;
      varMap[MVAV[iMVA]].index = iMVA;
      varMap[MVAV[iMVA]].maxIndex = MVAV.size();
      //cout << "\tiMVA = " << iMVA << "\tMVAV[iMVA] = " << MVAV[iMVA] << endl;
   }
   for(unsigned int iMVA=0; iMVA<MVAS.size(); iMVA++) {
      varMap[MVAS[iMVA]].use = true;
      varMap[MVAS[iMVA]].index = iMVA;
      varMap[MVAS[iMVA]].maxIndex = MVAS.size();
      //cout << "\tiMVA = " << iMVA << "\tMVAS[iMVA] = " << MVAV[iMVA] << endl;
   }
   //cout << "\tMVAV.size() = " << MVAV.size() << endl;
   //cout << "\tMVAS.size() = " << MVAS.size() << endl;
   //cout << "\tvarMap.size() (2)= " << varMap.size() << endl; 
}
