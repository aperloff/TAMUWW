#include "TAMUWW/SpecialTools/interface/Defs.hh"

namespace DEFS {

  using std::cout;
  using std::endl;
  using std::string;
  using std::map;
  using std::vector;
  
  namespace Ana{

    //---------------------------------------------------------------------------
    // A routine that returns the type given a string 
    string getTypeString(Type type){
      
      //the returning string
      if (type == SingleTopAnalysis)  return "Single Top Analysis";
      else if (type == HiggsAnalysis) return "Higgs Analysis";
      else if (type == WWAnalysis) 	  return "WW Analysis";
      else if (type == UNKNOWN)       return "UNKNOWN Analysis";

      cout<<"ERROR DEFS::getProcessTypeString AnalysisType "<<type<<" is unknown"<<endl;   
      
      return "";
      
    }//getAnalysisTypeString

    Type getAnaType(std::string str) {
       if(str == "SingleTopAnalysis")  return DEFS::Ana::SingleTopAnalysis;
       else if(str == "HiggsAnalysis") return DEFS::Ana::HiggsAnalysis;
       else if(str == "WWAnalysis")    return DEFS::Ana::WWAnalysis;
       else if(str == "UNKNOWN")       return DEFS::Ana::UNKNOWN;

       cout<<" ERROR  DEFS::getAnaType str ="<<str<<" not found"<<endl;

       return  DEFS::Ana::UNKNOWN;
    }//getAnaType
    
  }// namespace

  namespace PhysicsProcess {

    //---------------------------------------------------------------------------
    DEFS::PhysicsProcessType getProcessType(std::string str){
  
      if(str == "WH100")    	return  DEFS::PhysicsProcess::WH100;
      else if(str == "WH105")	return  DEFS::PhysicsProcess::WH105;
      else if(str == "WH110")	return  DEFS::PhysicsProcess::WH110;
      else if(str == "WH115")	return  DEFS::PhysicsProcess::WH115;
      else if(str == "WH120")	return  DEFS::PhysicsProcess::WH120;
      else if(str == "WH125")	return  DEFS::PhysicsProcess::WH125;
      else if(str == "WH130")	return  DEFS::PhysicsProcess::WH130;
      else if(str == "WH135")	return  DEFS::PhysicsProcess::WH135;
      else if(str == "WH140")	return  DEFS::PhysicsProcess::WH140;
      else if(str == "WH145")	return  DEFS::PhysicsProcess::WH145;
      else if(str == "WH150")	return  DEFS::PhysicsProcess::WH150;
      else if(str == "WH160")	return  DEFS::PhysicsProcess::WH160;
      else if(str == "WH170")	return  DEFS::PhysicsProcess::WH170;
      else if(str == "WH180")	return  DEFS::PhysicsProcess::WH180;
      else if(str == "WH190")	return  DEFS::PhysicsProcess::WH190;
      else if(str == "WH200")	return  DEFS::PhysicsProcess::WH200;
      else if(str == "WH250")	return  DEFS::PhysicsProcess::WH250;
      else if(str == "WH300")	return  DEFS::PhysicsProcess::WH300;
      else if(str == "WH350")	return  DEFS::PhysicsProcess::WH350;
      else if(str == "WH400")	return  DEFS::PhysicsProcess::WH400;
      else if(str == "WH450")	return  DEFS::PhysicsProcess::WH450;
      else if(str == "WH500")	return  DEFS::PhysicsProcess::WH500;
      else if(str == "WH550")	return  DEFS::PhysicsProcess::WH550;
      else if(str == "WH600")	return  DEFS::PhysicsProcess::WH600;
      else if(str == "WH700")	return  DEFS::PhysicsProcess::WH700;
      else if(str == "WH800")	return  DEFS::PhysicsProcess::WH800;
      else if(str == "WH900")	return  DEFS::PhysicsProcess::WH900;
      else if(str == "WH1000")	return  DEFS::PhysicsProcess::WH1000;

      else if(str == "ggH110")	return  DEFS::PhysicsProcess::ggH110;
      else if(str == "ggH115")	return  DEFS::PhysicsProcess::ggH115;
      else if(str == "ggH120")	return  DEFS::PhysicsProcess::ggH120;
      else if(str == "ggH125")	return  DEFS::PhysicsProcess::ggH125;
      else if(str == "ggH130")	return  DEFS::PhysicsProcess::ggH130;
      else if(str == "ggH135")	return  DEFS::PhysicsProcess::ggH135;
      else if(str == "ggH140")	return  DEFS::PhysicsProcess::ggH140;
      else if(str == "ggH145")	return  DEFS::PhysicsProcess::ggH145;
      else if(str == "ggH150")	return  DEFS::PhysicsProcess::ggH150;
      else if(str == "ggH160")	return  DEFS::PhysicsProcess::ggH160;
      else if(str == "ggH170")	return  DEFS::PhysicsProcess::ggH170;
      else if(str == "ggH180")	return  DEFS::PhysicsProcess::ggH180;
      else if(str == "ggH190")	return  DEFS::PhysicsProcess::ggH190;
      else if(str == "ggH200")	return  DEFS::PhysicsProcess::ggH200;
      else if(str == "ggH250")	return  DEFS::PhysicsProcess::ggH250;
      else if(str == "ggH300")	return  DEFS::PhysicsProcess::ggH300;
      else if(str == "ggH350")	return  DEFS::PhysicsProcess::ggH350;
      else if(str == "ggH400")	return  DEFS::PhysicsProcess::ggH400;
      else if(str == "ggH450")	return  DEFS::PhysicsProcess::ggH450;
      else if(str == "ggH500")	return  DEFS::PhysicsProcess::ggH500;
      else if(str == "ggH550")	return  DEFS::PhysicsProcess::ggH550;
      else if(str == "ggH600")	return  DEFS::PhysicsProcess::ggH600;
      else if(str == "ggH700")	return  DEFS::PhysicsProcess::ggH700;
      else if(str == "ggH800")	return  DEFS::PhysicsProcess::ggH800;
      else if(str == "ggH900")	return  DEFS::PhysicsProcess::ggH900;
      else if(str == "ggH1000")	return  DEFS::PhysicsProcess::ggH1000;

      else if(str == "qqH110")	return  DEFS::PhysicsProcess::qqH110;
      else if(str == "qqH115")	return  DEFS::PhysicsProcess::qqH115;
      else if(str == "qqH120")	return  DEFS::PhysicsProcess::qqH120;
      else if(str == "qqH125")	return  DEFS::PhysicsProcess::qqH125;
      else if(str == "qqH130")	return  DEFS::PhysicsProcess::qqH130;
      else if(str == "qqH135")	return  DEFS::PhysicsProcess::qqH135;
      else if(str == "qqH140")	return  DEFS::PhysicsProcess::qqH140;
      else if(str == "qqH145")	return  DEFS::PhysicsProcess::qqH145;
      else if(str == "qqH150")	return  DEFS::PhysicsProcess::qqH150;
      else if(str == "qqH160")	return  DEFS::PhysicsProcess::qqH160;
      else if(str == "qqH170")	return  DEFS::PhysicsProcess::qqH170;
      else if(str == "qqH180")	return  DEFS::PhysicsProcess::qqH180;
      else if(str == "qqH190")	return  DEFS::PhysicsProcess::qqH190;
      else if(str == "qqH200")	return  DEFS::PhysicsProcess::qqH200;
      else if(str == "qqH250")	return  DEFS::PhysicsProcess::qqH250;
      else if(str == "qqH300")	return  DEFS::PhysicsProcess::qqH300;
      else if(str == "qqH350")	return  DEFS::PhysicsProcess::qqH350;
      else if(str == "qqH400")	return  DEFS::PhysicsProcess::qqH400;
      else if(str == "qqH450")	return  DEFS::PhysicsProcess::qqH450;
      else if(str == "qqH500")	return  DEFS::PhysicsProcess::qqH500;
      else if(str == "qqH550")	return  DEFS::PhysicsProcess::qqH550;
      else if(str == "qqH600")	return  DEFS::PhysicsProcess::qqH600;
      else if(str == "qqH700")	return  DEFS::PhysicsProcess::qqH700;
      else if(str == "qqH800")	return  DEFS::PhysicsProcess::qqH800;
      else if(str == "qqH900")	return  DEFS::PhysicsProcess::qqH900;
      else if(str == "qqH1000")	return  DEFS::PhysicsProcess::qqH1000;

      else if(str == "WH_ZH_TTH_HToZZ_M125") return DEFS::PhysicsProcess::WH_ZH_TTH_HToZZ_M125;
      else if(str == "WH_ZH_TTH_HToWW_M125") return DEFS::PhysicsProcess::WH_ZH_TTH_HToWW_M125;
      else if(str == "WH_HToZZ_M125")        return DEFS::PhysicsProcess::WH_HToZZ_M125;
      else if(str == "ZH_HToZZ_M125")        return DEFS::PhysicsProcess::ZH_HToZZ_M125;
      else if(str == "TTH_HToZZ_M125")       return DEFS::PhysicsProcess::TTH_HToZZ_M125;
      else if(str == "WH_HToWW_M125")        return DEFS::PhysicsProcess::WH_HToWW_M125;
      else if(str == "ZH_HToWW_M125")        return DEFS::PhysicsProcess::ZH_HToWW_M125;
      else if(str == "TTH_HToWW_M125")       return DEFS::PhysicsProcess::TTH_HToWW_M125;
      else if(str == "WH125_HToBB")          return DEFS::PhysicsProcess::WH125_HToBB;
      else if(str == "WH125_HToZG")          return DEFS::PhysicsProcess::WH125_HToZG;
      else if(str == "TTH_Inclusive_M125")   return DEFS::PhysicsProcess::TTH_Inclusive_M125;
      else if(str == "TTH_HToBB_M125")       return DEFS::PhysicsProcess::TTH_HToBB_M125;

      else if(str == "ggH125_JESUp")                 return  DEFS::PhysicsProcess::ggH125_JESUp;
      else if(str == "ggH125_JESDown")               return  DEFS::PhysicsProcess::ggH125_JESDown;
      else if(str == "qqH125_JESUp")                 return  DEFS::PhysicsProcess::qqH125_JESUp;
      else if(str == "qqH125_JESDown")               return  DEFS::PhysicsProcess::qqH125_JESDown;
      else if(str == "WH_ZH_TTH_HToZZ_M125_JESUp")   return DEFS::PhysicsProcess::WH_ZH_TTH_HToZZ_M125_JESUp;
      else if(str == "WH_ZH_TTH_HToZZ_M125_JESDown") return DEFS::PhysicsProcess::WH_ZH_TTH_HToZZ_M125_JESDown;
      else if(str == "WH_ZH_TTH_HToWW_M125_JESUp")   return DEFS::PhysicsProcess::WH_ZH_TTH_HToWW_M125_JESUp;
      else if(str == "WH_ZH_TTH_HToWW_M125_JESDown") return DEFS::PhysicsProcess::WH_ZH_TTH_HToWW_M125_JESDown;
      else if(str == "WH_HToZZ_M125_JESUp")          return DEFS::PhysicsProcess::WH_HToZZ_M125_JESUp;
      else if(str == "WH_HToZZ_M125_JESDown")        return DEFS::PhysicsProcess::WH_HToZZ_M125_JESDown;
      else if(str == "WH_HToWW_M125_JESUp")          return DEFS::PhysicsProcess::WH_HToWW_M125_JESUp;
      else if(str == "WH_HToWW_M125_JESDown")        return DEFS::PhysicsProcess::WH_HToWW_M125_JESDown;
      else if(str == "ZH_HToZZ_M125_JESUp")          return DEFS::PhysicsProcess::ZH_HToZZ_M125_JESUp;
      else if(str == "ZH_HToZZ_M125_JESDown")        return DEFS::PhysicsProcess::ZH_HToZZ_M125_JESDown;
      else if(str == "ZH_HToWW_M125_JESUp")          return DEFS::PhysicsProcess::ZH_HToWW_M125_JESUp;
      else if(str == "ZH_HToWW_M125_JESDown")        return DEFS::PhysicsProcess::ZH_HToWW_M125_JESDown;
      else if(str == "TTH_HToZZ_M125_JESUp")         return DEFS::PhysicsProcess::TTH_HToZZ_M125_JESUp;
      else if(str == "TTH_HToZZ_M125_JESDown")       return DEFS::PhysicsProcess::TTH_HToZZ_M125_JESDown;
      else if(str == "TTH_HToWW_M125_JESUp")         return DEFS::PhysicsProcess::TTH_HToWW_M125_JESUp;
      else if(str == "TTH_HToWW_M125_JESDown")       return DEFS::PhysicsProcess::TTH_HToWW_M125_JESDown;
      else if(str == "WH_HToBB_M125_JESUp")          return DEFS::PhysicsProcess::WH_HToBB_M125_JESUp;
      else if(str == "WH_HToBB_M125_JESDown")        return DEFS::PhysicsProcess::WH_HToBB_M125_JESDown;
      else if(str == "TTH_HToBB_M125_JESUp")         return DEFS::PhysicsProcess::TTH_HToBB_M125_JESUp;
      else if(str == "TTH_HToBB_M125_JESDown")       return DEFS::PhysicsProcess::TTH_HToBB_M125_JESDown;

      else if(str == "STopS_T")	                  return  DEFS::PhysicsProcess::STopS_T;
      else if(str == "STopS_Tbar")	              return  DEFS::PhysicsProcess::STopS_Tbar;
      else if(str == "STopT_T")	                  return  DEFS::PhysicsProcess::STopT_T;
      else if(str == "STopT_Tbar")	              return  DEFS::PhysicsProcess::STopT_Tbar;
      else if(str == "STopTW_T")	                return  DEFS::PhysicsProcess::STopTW_T;
      else if(str == "STopTW_Tbar")	              return  DEFS::PhysicsProcess::STopTW_Tbar;
      else if(str == "Wcc")	                      return  DEFS::PhysicsProcess::Wcc;
      else if(str == "WJets")                     return  DEFS::PhysicsProcess::WJets;
      else if(str == "W1Jets")                    return  DEFS::PhysicsProcess::W1Jets;
      else if(str == "W2Jets")                    return  DEFS::PhysicsProcess::W2Jets;
      else if(str == "W3Jets")                    return  DEFS::PhysicsProcess::W3Jets;
      else if(str == "W4Jets")                    return  DEFS::PhysicsProcess::W4Jets;
      else if(str == "WLg")	                      return  DEFS::PhysicsProcess::WLg;
      else if(str == "Wgg")	                      return  DEFS::PhysicsProcess::Wgg;
      else if(str == "WLL")	                      return  DEFS::PhysicsProcess::WLL;
      else if(str == "Wbb")	                      return  DEFS::PhysicsProcess::WLb;
      else if(str == "WLb")	                      return  DEFS::PhysicsProcess::Wbb;
      else if(str == "QCD_ElEnriched")	          return  DEFS::PhysicsProcess::QCD_ElEnriched;
      else if(str == "QCD_ElFULL")                return  DEFS::PhysicsProcess::QCD_ElFULL;
      else if(str == "QCD_MuEnriched")	          return  DEFS::PhysicsProcess::QCD_MuEnriched;
      else if(str == "QCD_MuFULL")                return  DEFS::PhysicsProcess::QCD_MuFULL;
      else if(str == "QCD_Pt20to30_EMEnriched")   return  DEFS::PhysicsProcess::QCD_Pt20to30_EMEnriched;
      else if(str == "QCD_Pt30to80_EMEnriched")   return  DEFS::PhysicsProcess::QCD_Pt30to80_EMEnriched;
      else if(str == "QCD_Pt80to170_EMEnriched")  return  DEFS::PhysicsProcess::QCD_Pt80to170_EMEnriched;
      else if(str == "QCD_Pt170to250_EMEnriched") return  DEFS::PhysicsProcess::QCD_Pt170to250_EMEnriched;
      else if(str == "QCD_Pt250to350_EMEnriched") return  DEFS::PhysicsProcess::QCD_Pt250to350_EMEnriched;
      else if(str == "QCD_Pt350_EMEnriched")      return  DEFS::PhysicsProcess::QCD_Pt350_EMEnriched;
      else if(str == "WW")	                      return  DEFS::PhysicsProcess::WW;
      else if(str == "WZ")	                      return  DEFS::PhysicsProcess::WZ;
      else if(str == "WZbb")	                  return  DEFS::PhysicsProcess::WZbb;
      else if(str == "ZZ")	                      return  DEFS::PhysicsProcess::ZZ;
      else if(str == "TTbar")	                  return  DEFS::PhysicsProcess::TTbar;
      else if(str == "TTbar_JESUp")               return  DEFS::PhysicsProcess::TTbar_JESUp;
      else if(str == "TTbar_JESDown")             return  DEFS::PhysicsProcess::TTbar_JESDown;
      else if(str == "TTbarLJ")	                  return  DEFS::PhysicsProcess::TTbarLJ;
      else if(str == "TTbarDil")                  return  DEFS::PhysicsProcess::TTbarLJ;
      else if(str == "ZJets")                     return  DEFS::PhysicsProcess::ZJets;
      else if(str == "Ztautau")	                  return  DEFS::PhysicsProcess::Ztautau;
      else if(str == "ZJetsToLL_M50")             return  DEFS::PhysicsProcess::ZJetsToLL_M50;
      else if(str == "ZJetsToLL_M10To50")         return  DEFS::PhysicsProcess::ZJetsToLL_M10To50;
      else if(str == "SingleEl_Data")	          return  DEFS::PhysicsProcess::SingleEl_Data;
      else if(str == "SingleMu_Data")	          return  DEFS::PhysicsProcess::SingleMu_Data;

      else if(str == "STopS_T_JESUp")       return  DEFS::PhysicsProcess::STopS_T_JESUp;
      else if(str == "STopS_T_JESDown")     return  DEFS::PhysicsProcess::STopS_T_JESDown;
      else if(str == "STopS_Tbar_JESUp")    return  DEFS::PhysicsProcess::STopS_Tbar_JESUp;
      else if(str == "STopS_Tbar_JESDown")  return  DEFS::PhysicsProcess::STopS_Tbar_JESDown;
      else if(str == "STopT_T_JESUp")       return  DEFS::PhysicsProcess::STopT_T_JESUp;
      else if(str == "STopT_T_JESDown")     return  DEFS::PhysicsProcess::STopT_T_JESDown;
      else if(str == "STopT_Tbar_JESUp")    return  DEFS::PhysicsProcess::STopT_Tbar_JESUp;
      else if(str == "STopT_Tbar_JESDown")  return  DEFS::PhysicsProcess::STopT_Tbar_JESDown;
      else if(str == "STopTW_T_JESUp")      return  DEFS::PhysicsProcess::STopTW_T_JESUp;
      else if(str == "STopTW_T_JESDown")    return  DEFS::PhysicsProcess::STopTW_T_JESDown;
      else if(str == "STopTW_Tbar_JESUp")   return  DEFS::PhysicsProcess::STopTW_Tbar_JESUp;
      else if(str == "STopTW_Tbar_JESDown") return  DEFS::PhysicsProcess::STopTW_Tbar_JESDown;
      else if(str == "WJets_matchingup")    return  DEFS::PhysicsProcess::WJets_matchingup;
      else if(str == "WJets_matchingdown")  return  DEFS::PhysicsProcess::WJets_matchingdown;
      else if(str == "WJets_scaleup")       return  DEFS::PhysicsProcess::WJets_scaleup;
      else if(str == "WJets_scaledown")     return  DEFS::PhysicsProcess::WJets_scaledown;
      else if(str == "WJets_JESUp")         return  DEFS::PhysicsProcess::WJets_JESUp;
      else if(str == "WJets_JESDown")       return  DEFS::PhysicsProcess::WJets_JESDown;
      else if(str == "WW_JESUp")            return  DEFS::PhysicsProcess::WW_JESUp;
      else if(str == "WW_JESDown")          return  DEFS::PhysicsProcess::WW_JESDown;
      else if(str == "WZ_JESUp")            return  DEFS::PhysicsProcess::WZ_JESUp;
      else if(str == "WZ_JESDown")          return  DEFS::PhysicsProcess::WZ_JESDown;
      else if(str == "TTbar_JESUp")         return  DEFS::PhysicsProcess::TTbar_JESUp;
      else if(str == "TTbar_JESDown")       return  DEFS::PhysicsProcess::TTbar_JESDown;
      else if(str == "ZJets_JESUp")         return  DEFS::PhysicsProcess::ZJets_JESUp;
      else if(str == "ZJets_JESDown")       return  DEFS::PhysicsProcess::ZJets_JESDown;

      else if(str == "UNKNOWN")              return  DEFS::PhysicsProcess::UNKNOWN;
 
      cout<<" ERROR  DEFS::getProcessType str ="<<str<<" not found"<<endl;
  
      return  DEFS::PhysicsProcess::UNKNOWN;

    }//getProcessType

    //---------------------------------------------------------------------------
    string getTypeString(Type type){

      //the returning string
      if (type == WH100)	return string("WH100");
      else if (type == WH105)	return string("WH105");
      else if (type == WH110)	return string("WH110");
      else if (type == WH115)	return string("WH115");
      else if (type == WH120)	return string("WH120");
      else if (type == WH125)	return string("WH125");
      else if (type == WH130)	return string("WH130");
      else if (type == WH135)	return string("WH135");
      else if (type == WH140)	return string("WH140");
      else if (type == WH145)	return string("WH145");
      else if (type == WH150)	return string("WH150");
      else if (type == WH160)	return string("WH160");
      else if (type == WH170)	return string("WH170");
      else if (type == WH180)	return string("WH180");
      else if (type == WH190)	return string("WH190");
      else if (type == WH200)	return string("WH200");
      else if (type == WH250)	return string("WH250");
      else if (type == WH300)	return string("WH300");
      else if (type == WH350)	return string("WH350");
      else if (type == WH400)	return string("WH400");
      else if (type == WH450)	return string("WH450");
      else if (type == WH500)	return string("WH500");
      else if (type == WH550)	return string("WH550");
      else if (type == WH600)	return string("WH600");
      else if (type == WH700)	return string("WH700");
      else if (type == WH800)	return string("WH800");
      else if (type == WH900)	return string("WH900");
      else if (type == WH1000)return string("WH1000");

      else if (type == ggH110)	 return string("ggH110");
      else if (type == ggH115)	 return string("ggH115");
      else if (type == ggH120)	 return string("ggH120");
      else if (type == ggH125)	 return string("ggH125");
      else if (type == ggH130)	 return string("ggH130");
      else if (type == ggH135)	 return string("ggH135");
      else if (type == ggH140)	 return string("ggH140");
      else if (type == ggH145)	 return string("ggH145");
      else if (type == ggH150)	 return string("ggH150");
      else if (type == ggH160)	 return string("ggH160");
      else if (type == ggH170)	 return string("ggH170");
      else if (type == ggH180)	 return string("ggH180");
      else if (type == ggH190)	 return string("ggH190");
      else if (type == ggH200)	 return string("ggH200");
      else if (type == ggH250)	 return string("ggH250");
      else if (type == ggH300)	 return string("ggH300");
      else if (type == ggH350)	 return string("ggH350");
      else if (type == ggH400)	 return string("ggH400");
      else if (type == ggH450)	 return string("ggH450");
      else if (type == ggH500)	 return string("ggH500");
      else if (type == ggH550)	 return string("ggH550");
      else if (type == ggH600)	 return string("ggH600");
      else if (type == ggH700)	 return string("ggH700");
      else if (type == ggH800)	 return string("ggH800");
      else if (type == ggH900)	 return string("ggH900");
      else if (type == ggH1000)  return string("ggH1000");

      else if (type == qqH110)	 return string("qqH110");
      else if (type == qqH115)	 return string("qqH115");
      else if (type == qqH120)	 return string("qqH120");
      else if (type == qqH125)	 return string("qqH125");
      else if (type == qqH130)	 return string("qqH130");
      else if (type == qqH135)	 return string("qqH135");
      else if (type == qqH140)	 return string("qqH140");
      else if (type == qqH145)	 return string("qqH145");
      else if (type == qqH150)	 return string("qqH150");
      else if (type == qqH160)	 return string("qqH160");
      else if (type == qqH170)	 return string("qqH170");
      else if (type == qqH180)	 return string("qqH180");
      else if (type == qqH190)	 return string("qqH190");
      else if (type == qqH200)	 return string("qqH200");
      else if (type == qqH250)	 return string("qqH250");
      else if (type == qqH300)	 return string("qqH300");
      else if (type == qqH350)	 return string("qqH350");
      else if (type == qqH400)	 return string("qqH400");
      else if (type == qqH450)	 return string("qqH450");
      else if (type == qqH500)	 return string("qqH500");
      else if (type == qqH550)	 return string("qqH550");
      else if (type == qqH600)	 return string("qqH600");
      else if (type == qqH700)	 return string("qqH700");
      else if (type == qqH800)	 return string("qqH800");
      else if (type == qqH900)	 return string("qqH900");
      else if (type == qqH1000)  return string("qqH1000");

      else if (type == WH_ZH_TTH_HToZZ_M125) return string("WH_ZH_TTH_HToZZ_M125");
      else if (type == WH_ZH_TTH_HToWW_M125) return string("WH_ZH_TTH_HToWW_M125");
      else if (type == WH_HToZZ_M125)        return string("WH_HToZZ_M125");
      else if (type == ZH_HToZZ_M125)        return string("ZH_HToZZ_M125");
      else if (type == TTH_HToZZ_M125)       return string("TTH_HToZZ_M125");
      else if (type == WH_HToWW_M125)        return string("WH_HToWW_M125");
      else if (type == ZH_HToWW_M125)        return string("ZH_HToWW_M125");
      else if (type == TTH_HToWW_M125)       return string("TTH_HToWW_M125");
      else if (type == WH125_HToBB)          return string("WH125_HToBB");
      else if (type == WH125_HToZG)          return string("WH125_HToZG");
      else if (type == TTH_Inclusive_M125)   return string("TTH_Inclusive_M125");
      else if (type == TTH_HToBB_M125)       return string("TTH_HToBB_M125");

      else if (type == ggH125_JESUp)                 return string("ggH125_JESUp");
      else if (type == ggH125_JESDown)               return string("ggH125_JESDown");
      else if (type == qqH125_JESUp)                 return string("qqH125_JESUp");
      else if (type == qqH125_JESDown)               return string("qqH125_JESDown");
      else if (type == WH_ZH_TTH_HToZZ_M125_JESUp)   return string("WH_ZH_TTH_HToZZ_M125_JESUp");
      else if (type == WH_ZH_TTH_HToZZ_M125_JESDown) return string("WH_ZH_TTH_HToZZ_M125_JESDown");
      else if (type == WH_ZH_TTH_HToWW_M125_JESUp)   return string("WH_ZH_TTH_HToWW_M125_JESUp");
      else if (type == WH_ZH_TTH_HToWW_M125_JESDown) return string("WH_ZH_TTH_HToWW_M125_JESDown");
      else if (type == WH_HToZZ_M125_JESUp)          return string("WH_HToZZ_M125_JESUp");
      else if (type == WH_HToZZ_M125_JESDown)        return string("WH_HToZZ_M125_JESDown");
      else if (type == WH_HToWW_M125_JESUp)          return string("WH_HToWW_M125_JESUp");
      else if (type == WH_HToWW_M125_JESDown)        return string("WH_HToWW_M125_JESDown");
      else if (type == ZH_HToZZ_M125_JESUp)          return string("ZH_HToZZ_M125_JESUp");
      else if (type == ZH_HToZZ_M125_JESDown)        return string("ZH_HToZZ_M125_JESDown");
      else if (type == ZH_HToWW_M125_JESUp)          return string("ZH_HToWW_M125_JESUp");
      else if (type == ZH_HToWW_M125_JESDown)        return string("ZH_HToWW_M125_JESDown");
      else if (type == TTH_HToZZ_M125_JESUp)         return string("TTH_HToZZ_M125_JESUp");
      else if (type == TTH_HToZZ_M125_JESDown)       return string("TTH_HToZZ_M125_JESDown");
      else if (type == TTH_HToWW_M125_JESUp)         return string("TTH_HToWW_M125_JESUp");
      else if (type == TTH_HToWW_M125_JESDown)       return string("TTH_HToWW_M125_JESDown");
      else if (type == WH_HToBB_M125_JESDown)        return string("WH_HToBB_M125_JESDown");
      else if (type == WH_HToBB_M125_JESUp)          return string("WH_HToBB_M125_JESUp");
      else if (type == TTH_HToBB_M125_JESUp)         return string("TTH_HToBB_M125_JESUp");
      else if (type == TTH_HToBB_M125_JESDown)       return string("TTH_HToBB_M125_JESDown");

      else if (type == STopS_T)	                   return string("STopS_T");
      else if (type == STopS_Tbar)	               return string("STopS_Tbar");
      else if (type == STopT_T)	                   return string("STopT_T");
      else if (type == STopT_Tbar)	               return string("STopT_Tbar");
      else if (type == STopTW_T)	                 return string("STopTW_T");
      else if (type == STopTW_Tbar)	               return string("STopTW_Tbar");
      else if (type == Wcc) 	                     return string("Wcc");
      else if (type == WJets)	                     return string("WJets");
      else if (type == W1Jets)                     return string("W1Jets");
      else if (type == W2Jets)                     return string("W2Jets");
      else if (type == W3Jets)                     return string("W3Jets");
      else if (type == W4Jets)                     return string("W4Jets");
      else if (type == WLg)	                       return string("WLg");
      else if (type == Wgg)	                       return string("Wgg");
      else if (type == WLL)	                       return string("WLL");
      else if (type == Wbb)	                       return string("Wbb");
      else if (type == QCD_ElEnriched)	           return string("QCD_ElEnriched");
      else if (type == QCD_ElFULL)                 return string("QCD_ElFULL");
      else if (type == QCD_MuEnriched)	           return string("QCD_MuEnriched");
      else if (type == QCD_MuFULL)                 return string("QCD_MuFULL");
      else if (type == QCD_Pt20to30_EMEnriched )   return string ("QCD_Pt20to30_EMEnriched");
      else if (type == QCD_Pt30to80_EMEnriched )   return string ("QCD_Pt30to80_EMEnriched");
      else if (type == QCD_Pt80to170_EMEnriched )  return string ("QCD_Pt80to170_EMEnriched");
      else if (type == QCD_Pt170to250_EMEnriched ) return string ("QCD_Pt170to250_EMEnriched");
      else if (type == QCD_Pt250to350_EMEnriched ) return string ("QCD_Pt250to350_EMEnriched");
      else if (type == QCD_Pt350_EMEnriched )      return string ("QCD_Pt350_EMEnriched");
      else if (type == WW)	                       return string("WW");
      else if (type == WZ)	                       return string("WZ");
      else if (type == ZZ)	                       return string("ZZ");
      else if (type == TTbar)	                   return string("TTbar");
      else if (type == TTbarLJ)	                   return string("TTbarLJ");
      else if (type == TTbarDil)                   return string("TTbarDil");
      else if (type == ZJets)	                   return string("ZJets");
      else if (type == Ztautau)	                   return string("Ztautau");
      else if (type == ZJetsToLL_M50)              return string("ZJetsToLL_M50");
      else if (type == ZJetsToLL_M10To50)          return string("ZJetsToLL_M10To50");
      else if (type == SingleEl_Data)              return string("SingleEl_Data");
      else if (type == SingleMu_Data)	             return string("SingleMu_Data");

      else if (type == STopTW_Tbar_JESUp)   return string("STopTW_Tbar_JESUp");
      else if (type == STopTW_Tbar_JESDown) return string("STopTW_Tbar_JESDown");
      else if (type == STopTW_T_JESUp)      return string("STopTW_T_JESUp");
      else if (type == STopTW_T_JESDown)    return string("STopTW_T_JESDown");
      else if (type == STopT_Tbar_JESUp)    return string("STopT_Tbar_JESUp");
      else if (type == STopT_Tbar_JESDown)  return string("STopT_Tbar_JESDown");
      else if (type == STopT_T_JESUp)       return string("STopT_T_JESUp");
      else if (type == STopT_T_JESDown)     return string("STopT_T_JESDown");
      else if (type == STopS_Tbar_JESUp)    return string("STopS_Tbar_JESUp");
      else if (type == STopS_Tbar_JESDown)  return string("STopS_Tbar_JESDown");
      else if (type == STopS_T_JESUp)       return string("STopS_T_JESUp");
      else if (type == STopS_T_JESDown)     return string("STopS_T_JESDown");
      else if (type == WJets_matchingup)    return string("WJets_matchingup");
      else if (type == WJets_matchingdown)  return string("WJets_matchingdown");
      else if (type == WJets_scaleup)       return string("WJets_scaleup");
      else if (type == WJets_scaledown)     return string("WJets_scaledown");
      else if (type == WJets_JESUp)         return string("WJets_JESUp");
      else if (type == WJets_JESDown)       return string("WJets_JESDown");
      else if (type == WZ_JESUp)            return string("WZ_JESUp");
      else if (type == WZ_JESDown)          return string("WZ_JESDown");
      else if (type == WW_JESUp)            return string("WW_JESUp");
      else if (type == WW_JESDown)          return string("WW_JESDown");
      else if (type == TTbar_JESUp)         return string("TTbar_JESUp");
      else if (type == TTbar_JESDown)       return string("TTbar_JESDown");
      else if (type == ZJets_JESUp)         return string("ZJets_JESUp");
      else if (type == ZJets_JESDown)       return string("ZJets_JESDown");

      else if (type == UNKNOWN)              return string("UNKNOWN");

      else cout<<" ERROR  DEFS::getProcessTypeString type ="<<type<<" not found"<<endl;

      return "ERROR in DEFS::getProcessTypeString ";

    }//getTypeString


    //---------------------------------------------------------------------------
    bool isHiggs(Type type){

      //the returning string
      if (type == WH100 || type == WH105 || type == WH110 ||
          type == WH115 || type == WH120 || type == WH125 ||
          type == WH130 || type == WH135 || type == WH140 ||
          type == WH145 || type == WH150)
         return true;

      if ( type == ggH100 || type == ggH105 || type == ggH110 || type == ggH115 || type == ggH120 ||
	   type == ggH125 || type == ggH130 || type == ggH135 || type == ggH140 || type == ggH145 || 
	   type == ggH150 || type == ggH160 || type == ggH170 || type == ggH180 || type == ggH190 || 
	   type == ggH200 || type == ggH250 || type == ggH300 || type == ggH350 || type == ggH400 || 
	   type == ggH450 || type == ggH500 || type == ggH550 || type == ggH600 || type == ggH700 || 
	   type == ggH800 || type == ggH900 || type == ggH1000)
         return true;

      if ( type == qqH100 || type == qqH105 || type == qqH110 || type == qqH115 || type == qqH120 ||
	   type == qqH125 || type == qqH130 || type == qqH135 || type == qqH140 || type == qqH145 ||
	   type == qqH150 || type == qqH160 || type == qqH170 || type == qqH180 || type == qqH190 ||
	   type == qqH200 || type == qqH250 || type == qqH300 || type == qqH350 || type == qqH400 ||
	   type == qqH450 || type == qqH500 || type == qqH550 || type == qqH600 || type == qqH700 ||
	   type == qqH800 || type == qqH900 || type == qqH1000)
         return true;

      if ( type == WH_ZH_TTH_HToWW_M125 || type == WH_ZH_TTH_HToZZ_M125 || type == TTH_HToBB_M125 || 
          type == TTH_Inclusive_M125 || type == WH_HToWW_M125 || type == ZH_HToWW_M125 || type == TTH_HToWW_M125 || 
          type == WH_HToZZ_M125 || type == ZH_HToZZ_M125 || type == TTH_HToZZ_M125)
         return true;

      if ( type == ggH125_JESUp || type == ggH125_JESDown || type == qqH125_JESUp || type == qqH125_JESDown ||
          type == WH125_HToBB || type == WH_HToBB_M125_JESUp || type == WH_HToBB_M125_JESDown ||
          type == TTH_HToBB_M125_JESUp || type == TTH_HToBB_M125_JESDown ||
          type == WH_ZH_TTH_HToWW_M125_JESUp || type == WH_ZH_TTH_HToWW_M125_JESDown || 
          type == WH_ZH_TTH_HToZZ_M125_JESUp || type == WH_ZH_TTH_HToZZ_M125_JESDown ||
          type == WH_HToWW_M125_JESUp || type == WH_HToWW_M125_JESDown || 
          type == WH_HToZZ_M125_JESUp || type == WH_HToZZ_M125_JESDown ||
          type == ZH_HToWW_M125_JESUp || type == ZH_HToWW_M125_JESDown || 
          type == ZH_HToZZ_M125_JESUp || type == ZH_HToZZ_M125_JESDown ||
          type == TTH_HToWW_M125_JESUp || type == TTH_HToWW_M125_JESDown || 
          type == TTH_HToZZ_M125_JESUp || type == TTH_HToZZ_M125_JESDown)
         return true;

      return false;

    }//isHiggs
     
    //---------------------------------------------------------------------------
    // A routine that returns the title for a given process type
    std::string getTypeTitle(Type type) {
       if (type == WW || type == WZ || type == ZZ)
          return string("WW+WZ+ZZ");
       else if (type == STopS_T || type == STopS_Tbar || type == STopT_T || type == STopT_Tbar ||
                type == STopTW_T || type == STopTW_Tbar)
          return string("Single Top");
       else if (type == ggH125 || type == qqH125 || type == WH_HToWW_M125 || type == ZH_HToWW_M125 ||
                type == TTH_HToWW_M125 || type == WH_ZH_TTH_HToWW_M125)
            return string("H(125)->WW");
       else if (type == WH125 || type == WH125_HToBB || type == TTH_HToBB_M125)
            return string("H(125)->bb");
       else if (type == WH_HToZZ_M125 || type == ZH_HToZZ_M125 || type == TTH_HToZZ_M125 ||
                type == WH_ZH_TTH_HToZZ_M125)
            return string("H(125)->ZZ");
       else if (type == TTH_Inclusive_M125)
            return string("TTH(125)");
       //else if (type == ggH125 || type == qqH125 || type == WH125)
       //   return string("ggH+WH+qqH(125)");
       //else if (type == WH_ZH_TTH_HToWW_M125 || type == WH_ZH_TTH_HToZZ_M125 || type == TTH_HToBB_M125 ||
       //         type == TTH_Inclusive_M125)
       //   return string("WH_ZH_TTH(125)");
       else if (type == WJets || type == W1Jets || type == W2Jets || type == W3Jets || type == W4Jets)
          return string("W+Jets");
       else if (type == ZJets || type == ZJetsToLL_M50 || type == ZJetsToLL_M10To50)
          return string("Z+Jets");
       else if (type == SingleEl_Data || type == SingleMu_Data)
          return string("DATA");
       else if (type == QCD_ElEnriched || type == QCD_MuEnriched || type == QCD_ElFULL || type == QCD_MuFULL ||
		type == QCD_Pt20to30_EMEnriched   || type == QCD_Pt30to80_EMEnriched   ||
		type == QCD_Pt80to170_EMEnriched  || type == QCD_Pt170to250_EMEnriched ||
		type == QCD_Pt250to350_EMEnriched || type == QCD_Pt350_EMEnriched )
          return string("QCD");
       else if (type == WJets_matchingup)
          return string("Systematic_matchingup");
       else if (type == WJets_matchingdown)
          return string("Systematic_matchingdown");
       else if (type == WJets_scaleup)
          return string("Systematic_scaleup");
       else if (type == WJets_scaledown)
          return string("Systematic_scaledown");
       else if (type == STopS_T_JESUp || type == STopS_Tbar_JESUp || type == STopT_T_JESUp || type == STopT_Tbar_JESUp || type == STopTW_T_JESUp || type == STopTW_Tbar_JESUp || 
                type == WJets_JESUp || type == WW_JESUp || type == WZ_JESUp || type == TTbar_JESUp || type == ZJets_JESUp ||
                type == ggH125_JESUp || type == qqH125_JESUp || 
                type == WH_ZH_TTH_HToWW_M125_JESUp || 
                type == WH_HToWW_M125_JESUp || type == ZH_HToWW_M125_JESUp || type == TTH_HToWW_M125_JESUp || 
                type == WH_HToBB_M125_JESUp || type == TTH_HToBB_M125_JESUp || 
                type == WH_ZH_TTH_HToZZ_M125_JESUp || 
                type == WH_HToZZ_M125_JESUp || type == ZH_HToZZ_M125_JESUp || type == TTH_HToZZ_M125_JESUp)
          return string("Systematic_JESUp");
       else if (type == STopS_T_JESDown || type == STopS_Tbar_JESDown || type == STopT_T_JESDown || type == STopT_Tbar_JESDown || type == STopTW_T_JESDown || type == STopTW_Tbar_JESDown || 
                type == WJets_JESDown || type == WW_JESDown || type == WZ_JESDown || type == TTbar_JESDown || type == ZJets_JESDown ||
                type == ggH125_JESDown || type == qqH125_JESDown ||
                type == WH_ZH_TTH_HToWW_M125_JESDown ||
                type == WH_HToWW_M125_JESDown || type == ZH_HToWW_M125_JESDown || type == TTH_HToWW_M125_JESDown ||
                type == WH_HToBB_M125_JESDown || type == TTH_HToBB_M125_JESDown ||
                type == WH_ZH_TTH_HToZZ_M125_JESDown ||
                type == WH_HToZZ_M125_JESDown || type == ZH_HToZZ_M125_JESDown || type == TTH_HToZZ_M125_JESDown)
          return string("Systematic_JESDown");
       else
          return getTypeString(type);
    }//getTypeTitle

    //---------------------------------------------------------------------------
    // A  routine that returns the Color_t for a given process type
    Color_t getProcessColor(Type type){
       
       if (type == WW)                              return kMagenta+1;//kPink;
       else if (type == WW_JESUp)                   return kMagenta+1;//kPink;
       else if (type == WW_JESDown)                 return kMagenta+1;//kPink;
       else if (type == WZ)                         return kCyan+1;//kBlue;
       else if (type == WZ_JESUp)                   return kCyan+1;//kBlue;
       else if (type == WZ_JESDown)                 return kCyan+1;//kBlue;
       else if (type == ZZ)                         return kYellow;//kBlue-1;
       else if (type == WJets)                      return kGreen-3;//kTeal+2;
       else if (type == WJets_matchingup)           return kGreen-3;//kTeal+2;
       else if (type == WJets_matchingdown)         return kGreen-3;//kTeal+2;
       else if (type == WJets_scaleup)              return kGreen-3;//kTeal+2;
       else if (type == WJets_scaledown)            return kGreen-3;//kTeal+2;
       else if (type == WJets_JESUp)                return kGreen-3;//kTeal+2;
       else if (type == WJets_JESDown)              return kGreen-3;//kTeal+2;
       else if (type == W1Jets)                     return kGreen-3;//kTeal+2;
       else if (type == W2Jets)                     return kGreen-3;//kTeal+2;
       else if (type == W3Jets)                     return kGreen-3;//kTeal+2;
       else if (type == W4Jets)                     return kGreen-3;//kTeal+2;
       else if (type == ZJets)                      return kBlue-3;//kPink-8;
       else if (type == ZJets_JESUp)                return kBlue-3;//kPink-8;
       else if (type == ZJets_JESDown)              return kBlue-3;//kPink-8;
       else if (type == ZJetsToLL_M50)              return kBlue-3;//kPink-8;
       else if (type == ZJetsToLL_M10To50)          return kBlue-3;//kPink-8;
       else if (type == QCD_ElEnriched)             return kRed-6;//kYellow+1;
       else if (type == QCD_ElFULL)                 return kRed-6;//kYellow+1;
       else if (type == QCD_MuEnriched)             return kRed-6;//kYellow+1;
       else if (type == QCD_MuFULL)                 return kRed-6;//kYellow+1;
       else if (type == QCD_Pt20to30_EMEnriched )   return kRed-6;//kYellow+1;
       else if (type == QCD_Pt30to80_EMEnriched )   return kRed-6;//kYellow+1;
       else if (type == QCD_Pt80to170_EMEnriched )  return kRed-6;//kYellow+1;
       else if (type == QCD_Pt170to250_EMEnriched ) return kRed-6;//kYellow+1;
       else if (type == QCD_Pt250to350_EMEnriched ) return kRed-6;//kYellow+1;
       else if (type == QCD_Pt350_EMEnriched )      return kRed-6;//kYellow+1;
       else if (type == STopTW_Tbar)                return kOrange+1; //kGreen+3;
       else if (type == STopTW_Tbar_JESUp)          return kOrange+1; //kGreen+3;
       else if (type == STopTW_Tbar_JESDown)        return kOrange+1; //kGreen+3;
       else if (type == STopTW_T)                   return kOrange+1; //kMagenta;
       else if (type == STopTW_T_JESUp)             return kOrange+1; //kMagenta;
       else if (type == STopTW_T_JESDown)           return kOrange+1; //kMagenta;
       else if (type == STopS_Tbar)                 return kOrange+1; //kBlue+3;
       else if (type == STopS_Tbar_JESUp)           return kOrange+1; //kBlue+3;
       else if (type == STopS_Tbar_JESDown)         return kOrange+1; //kBlue+3;
       else if (type == STopS_T)                    return kOrange+1; //kBlue;
       else if (type == STopS_T_JESUp)              return kOrange+1; //kBlue;
       else if (type == STopS_T_JESDown)            return kOrange+1; //kBlue;
       else if (type == STopT_Tbar)                 return kOrange+1; //kCyan+3;
       else if (type == STopT_Tbar_JESUp)           return kOrange+1; //kCyan+3;
       else if (type == STopT_Tbar_JESDown)         return kOrange+1; //kCyan+3;
       else if (type == STopT_T)                    return kOrange+1;
       else if (type == STopT_T_JESUp)              return kOrange+1;
       else if (type == STopT_T_JESDown)            return kOrange+1;
       else if (type == TTbar)                      return kRed-3;//kAzure-2; 
       else if (type == TTbar_JESUp)                return kRed-3;//kAzure-2; 
       else if (type == TTbar_JESDown)              return kRed-3;//kAzure-2; 
       else if (type == ggH125)                     return kMagenta-3;//kRed-4;//was kRed+2; 
       else if (type == ggH125_JESUp)               return kMagenta-3;//kRed-4;//was kRed+2; 
       else if (type == ggH125_JESDown)             return kMagenta-3;//kRed-4;//was kRed+2; 
       else if (type == qqH125)                     return kMagenta-3;//kRed-4;
       else if (type == qqH125_JESUp)               return kMagenta-3;//kRed-4;
       else if (type == qqH125_JESDown)             return kMagenta-3;//kRed-4;
       else if (type == WH125)                      return kMagenta-3;//kRed-4;
       else if (type == WH125_HToBB)                return kMagenta-3;
       else if (type == WH_HToBB_M125_JESUp)        return kMagenta-3;
       else if (type == WH_HToBB_M125_JESDown)      return kMagenta-3;
       else if (type == WH_ZH_TTH_HToWW_M125)       return kMagenta-3;
       else if (type == WH_ZH_TTH_HToWW_M125_JESUp)   return kMagenta-3;
       else if (type == WH_ZH_TTH_HToWW_M125_JESDown) return kMagenta-3;
       else if (type == WH_ZH_TTH_HToZZ_M125)       return kMagenta-3;
       else if (type == WH_ZH_TTH_HToZZ_M125_JESUp)   return kMagenta-3;
       else if (type == WH_ZH_TTH_HToZZ_M125_JESDown) return kMagenta-3;
       else if (type == TTH_HToWW_M125)             return kMagenta-3;
       else if (type == TTH_HToWW_M125_JESUp)       return kMagenta-3;
       else if (type == TTH_HToWW_M125_JESDown)     return kMagenta-3;
       else if (type == ZH_HToWW_M125)              return kMagenta-3;
       else if (type == ZH_HToWW_M125_JESUp)        return kMagenta-3;
       else if (type == ZH_HToWW_M125_JESDown)      return kMagenta-3;
       else if (type == WH_HToWW_M125)              return kMagenta-3;
       else if (type == WH_HToWW_M125_JESUp)        return kMagenta-3;
       else if (type == WH_HToWW_M125_JESDown)      return kMagenta-3;
       else if (type == TTH_HToZZ_M125)             return kMagenta-3;
       else if (type == TTH_HToZZ_M125_JESUp)       return kMagenta-3;
       else if (type == TTH_HToZZ_M125_JESDown)     return kMagenta-3;
       else if (type == ZH_HToZZ_M125)              return kMagenta-3;
       else if (type == ZH_HToZZ_M125_JESUp)        return kMagenta-3;
       else if (type == ZH_HToZZ_M125_JESDown)      return kMagenta-3;
       else if (type == WH_HToZZ_M125)              return kMagenta-3;
       else if (type == WH_HToZZ_M125_JESUp)        return kMagenta-3;
       else if (type == WH_HToZZ_M125_JESDown)      return kMagenta-3;
       else if (type == TTH_HToBB_M125)             return kMagenta-3;
       else if (type == TTH_HToBB_M125_JESUp)       return kMagenta-3;
       else if (type == TTH_HToBB_M125_JESDown)     return kMagenta-3;
       else if (type == TTH_Inclusive_M125)         return kMagenta-3;
       else if (type == SingleEl_Data)              return kBlack;
       else if (type == SingleMu_Data)              return kBlack;
       else{
          cout << "WARNING Plotter::GetProcessColor() Unknown process name=|"<<getTypeString(type)
               <<"|. Returning process color as kBlack." << endl;
       }
       
       return kBlack;
       
    }//getProcessColor

  } //namespace   
   
  //---------------------------------------------------------------------------
  // A routine that returns the string given the Lepton Category
  string getLeptonCatString(LeptonCat type){
    if (type == electron)  return "electron";
    else if (type == muon) return "muon";
    else if (type == both) return "both";
    else if (type == none) return "none";

    cout<<"ERROR in DEFS::getLeptonCatString type="<<type<<" not defined"<<endl;

    return "";

  }//getLeptonCatString

  //---------------------------------------------------------------------------
  // for historical reasons
  std::string getEventCatString(LeptonCat a){
    return  getLeptonCatString(a);
  }


  // A routine that returns the string given the Lepton Category
  LeptonCat getLeptonCat(string type){
    TString s = type;
    s.ToUpper();

    if (s.Contains("NONE"))          return none;
    else if (s.Contains("MUON"))     return muon;
    else if (s.Contains("ELECTRON")) return electron;
    else if (s.Contains("BOTH"))     return both;
    else cout<<"ERROR in DEFS::getLeptonCatString type="<<type<<" not defined!. Returning none."<<endl;

    return none;

  }//getLeptonCat


  //---------------------------------------------------------------------------
  string getJetBinString(JetBin jBin){

    //the returning string
    if (jBin == jets0)      return string("jets0");
    else if (jBin == jet1)  return string("jet1");
    else if (jBin == jets2) return string("jets2");
    else if (jBin == jets3) return string("jets3");
    else if (jBin == jets4) return string("jets4");
    else if (jBin == jets5) return string("jets5");
    else cout<<" ERROR  DEFS::geJetBinString jBin ="<<jBin<<" not found"<<endl;

    return "ERROR in DEFS::getJetBinString ";

  }//getJetBinString

  //---------------------------------------------------------------------------
  JetBin getJetBin(std::string str){

    if (str == "jets0")      return jets0;
    else if (str == "jet1")  return jet1;
    else if (str == "jets2") return jets2;
    else if (str == "jets3") return jets3;
    else if (str == "jets4") return jets4;
    else if (str == "jets5") return jets5;

    cout<<" ERROR  DEFS::getJetBin str ="<<str<<" not found"<<endl;
    return  jets0;

  }//getJetBin

  //---------------------------------------------------------------------------
  JetBin getJetBin(int nJets, JetBin inclusive_bin){

    if(nJets >= int(inclusive_bin)) return inclusive_bin;

    if (nJets == 0)      return jets0;
    else if (nJets == 1) return jet1;
    else if (nJets == 2) return jets2;
    else if (nJets == 3) return jets3;
    else if (nJets == 4) return jets4;
    else if (nJets == 5) return jets5;

    else cout<<" ERROR  DEFS::geJetBin jet(s)"<<nJets<<" not found"<<endl;
    return jets0;

  }//getJetBin

  //---------------------------------------------------------------------------
  int getNBinsX(JetBin jBin, LeptonCat lcat, std::string variable){

    if (jBin == jets2 && lcat==electron && variable=="KinBDT")        return int(nbinsx_jets2_electron_KinBDT);
    else if (jBin == jets2 && lcat==electron && variable=="MEBDT")    return int(nbinsx_jets2_electron_MEBDT);
    else if (jBin == jets2 && lcat==electron && variable=="KinMEBDT") return int(nbinsx_jets2_electron_KinMEBDT);
    else if (jBin == jets2 && lcat==muon && variable=="KinBDT")       return int(nbinsx_jets2_muon_KinBDT);
    else if (jBin == jets2 && lcat==muon && variable=="MEBDT")        return int(nbinsx_jets2_muon_MEBDT);
    else if (jBin == jets2 && lcat==muon && variable=="KinMEBDT")     return int(nbinsx_jets2_muon_KinMEBDT);
    else if (jBin == jets3 && lcat==electron && variable=="KinBDT")   return int(nbinsx_jets3_electron_KinBDT);
    else if (jBin == jets3 && lcat==electron && variable=="MEBDT")    return int(nbinsx_jets3_electron_MEBDT);
    else if (jBin == jets3 && lcat==electron && variable=="KinMEBDT") return int(nbinsx_jets3_electron_KinMEBDT);
    else if (jBin == jets3 && lcat==muon && variable=="KinBDT")       return int(nbinsx_jets3_muon_KinBDT);
    else if (jBin == jets3 && lcat==muon && variable=="MEBDT")        return int(nbinsx_jets3_muon_MEBDT);
    else if (jBin == jets3 && lcat==muon && variable=="KinMEBDT")     return int(nbinsx_jets3_muon_KinMEBDT);
    else if (jBin == jets4 && lcat==electron && variable=="KinBDT")   return int(nbinsx_jets4_electron_KinBDT);
    else if (jBin == jets4 && lcat==electron && variable=="MEBDT")    return int(nbinsx_jets4_electron_MEBDT);
    else if (jBin == jets4 && lcat==electron && variable=="KinMEBDT") return int(nbinsx_jets4_electron_KinMEBDT);
    else if (jBin == jets4 && lcat==muon && variable=="KinBDT")       return int(nbinsx_jets4_muon_KinBDT);
    else if (jBin == jets4 && lcat==muon && variable=="MEBDT")        return int(nbinsx_jets4_muon_MEBDT);
    else if (jBin == jets4 && lcat==muon && variable=="KinMEBDT")     return int(nbinsx_jets4_muon_KinMEBDT);
    else cout << "ERROR in DEFS::getNBinsX jBin="<<jBin<<" not found!. Returning default number of bins."<<endl;

    return int(nbinsx_default);
  }//getNBinsX

  //---------------------------------------------------------------------------
  vector<Double_t> getBinsX(JetBin jBin, LeptonCat lcat, std::string variable) {
    //With Signal cuts
    if (jBin==jets2 && lcat==electron && variable=="KinBDT") {
      vector<Double_t> ret = {-1,-0.64,-0.6,-0.58,-0.56,-0.54,-0.52,-0.5,-0.48,-0.46,-0.44,-0.42,-0.4,-0.38,-0.36,
                              -0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,
                              -0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,1};
      return ret;
    }
    else if (jBin==jets2 && lcat==electron && variable=="MEBDT") {
      vector<Double_t> ret = {-1,-0.4,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,
                              -0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,1};
      return ret;
    }
    else if (jBin==jets2 && lcat==electron && variable=="KinMEBDT") {
      vector<Double_t> ret = {-1,-0.46,-0.44,-0.42,-0.4,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,
                              -0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,
                              0.24,1};
      return ret;
    }
    else if (jBin==jets2 && lcat==muon && variable=="KinBDT") {
      vector<Double_t> ret = {-1,-0.66,-0.62,-0.6,-0.58,-0.56,-0.54,-0.52,-0.5,-0.48,-0.46,-0.44,-0.42,-0.4,-0.38,-0.36,
                              -0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,
                              -0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,1};
      return ret;
    }
    else if (jBin==jets2 && lcat==muon && variable=="MEBDT") {
      vector<Double_t> ret = {-1,-0.42,-0.4,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,
                              -0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.18,1};
      return ret;
    }
    else if (jBin==jets2 && lcat==muon && variable=="KinMEBDT") {
      vector<Double_t> ret = {-1,-0.46,-0.44,-0.42,-0.4,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,
                              -0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,
                              0.24,1};
      return ret;
    }
    else if (jBin==jets3 && lcat==electron && variable=="KinBDT") {
      vector<Double_t> ret = {-1,-0.46,-0.42,-0.4,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,
                              -0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,1};
      return ret;
    }
    else if (jBin==jets3 && lcat==electron && variable=="MEBDT") {
      vector<Double_t> ret = {-1,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,
                              0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,1};
      return ret;
    }
    else if (jBin==jets3 && lcat==electron && variable=="KinMEBDT") {
      vector<Double_t> ret = {-1,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,
                              -0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,1};
      return ret;
    }
    else if (jBin==jets3 && lcat==muon && variable=="KinBDT") {
      vector<Double_t> ret = {-1,-0.46,-0.42,-0.4,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,
                              -0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,1};
      return ret;
    }
    else if (jBin==jets3 && lcat==muon && variable=="MEBDT") {
      vector<Double_t> ret = {-1,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,
                              0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,1};
      return ret;
    }
    else if (jBin==jets3 && lcat==muon && variable=="KinMEBDT") {
      vector<Double_t> ret = {-1,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,
                              -0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,1};
      return ret;
    }
    else if (jBin==jets4 && lcat==electron && variable=="KinBDT") {
      vector<Double_t> ret = {-1,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,
                              0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,1};
      return ret;
    }
    else if (jBin==jets4 && lcat==electron && variable=="MEBDT") {
      vector<Double_t> ret = {-1,-0.24,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,
                              0.08,0.1,0.12,0.14,1};
      return ret;
    }
    else if (jBin==jets4 && lcat==electron && variable=="KinMEBDT") {
      vector<Double_t> ret = {-1,-0.22,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,1};
      return ret;
    }
    else if (jBin==jets4 && lcat==muon && variable=="KinBDT") {
      vector<Double_t> ret = {-1,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,
                              -0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,1};
      return ret;
    }
    else if (jBin==jets4 && lcat==muon && variable=="MEBDT") {
      vector<Double_t> ret = {-1,-0.22,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,1};
      return ret;
    }
    else if (jBin==jets4 && lcat==muon && variable=="KinMEBDT") {
      vector<Double_t> ret = {-1,-0.22,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,1};
      return ret;
    }
    //With HailMary cuts
    /*
    if (jBin==jets2 && lcat==electron && variable=="KinBDT") {
      vector<Double_t> ret = {-1,-0.42,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,
                              -0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,1};
      return ret;
    }
    else if (jBin==jets2 && lcat==electron && variable=="MEBDT") {
      vector<Double_t> ret = {-1,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,
                              0,0.02,0.04,1};
      return ret;
    }
    else if (jBin==jets2 && lcat==electron && variable=="KinMEBDT") {
      vector<Double_t> ret = {-1,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,
                              -0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,1};
      return ret;
    }
    else if (jBin==jets2 && lcat==muon && variable=="KinBDT") {
      vector<Double_t> ret = {-1,-0.42,-0.38,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,
                              -0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.2,1};
      return ret;
    }
    else if (jBin==jets2 && lcat==muon && variable=="MEBDT") {
      vector<Double_t> ret = {-1,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,1};
      return ret;
    }
    else if (jBin==jets2 && lcat==muon && variable=="KinMEBDT") {
      vector<Double_t> ret = {-1,-0.32,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,
                              0.04,0.06,0.08,0.1,0.12,0.14,0.16,1};
      return ret;
    }
    else if (jBin==jets3 && lcat==electron && variable=="KinBDT") {
      vector<Double_t> ret = {-1,-0.26,-0.22,-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,
                              0.12,0.16,1};
      return ret;
    }
    else if (jBin==jets3 && lcat==electron && variable=="MEBDT") {
      vector<Double_t> ret = {-1,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,1};
      return ret;
    }
    else if (jBin==jets3 && lcat==electron && variable=="KinMEBDT") {
      vector<Double_t> ret = {-1,-0.18,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,1};
      return ret;
    }
    else if (jBin==jets3 && lcat==muon && variable=="KinBDT") {
      vector<Double_t> ret = {-1,-0.22,-0.18,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,1};
      return ret;
    }
    else if (jBin==jets3 && lcat==muon && variable=="MEBDT") {
      vector<Double_t> ret = {-1,-0.14,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,1};
      return ret;
    }
    else if (jBin==jets3 && lcat==muon && variable=="KinMEBDT") {
      vector<Double_t> ret = {-1,-0.16,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,1};
      return ret;
    }
    else if (jBin==jets4 && lcat==electron && variable=="KinBDT") {
      vector<Double_t> ret = {-1,-0.12,-0.08,-0.04,0,0.04,0.08,0.14,1};
      return ret;
    }
    else if (jBin==jets4 && lcat==electron && variable=="MEBDT") {
      vector<Double_t> ret = {-1,-0.04,0,0.02,0.04,0.06,0.08,1};
      return ret;
    }
    else if (jBin==jets4 && lcat==electron && variable=="KinMEBDT") {
      vector<Double_t> ret = {-1,-0.08,-0.04,-0.02,0.02,0.04,0.08,0.12,1};
      return ret;
    }
    else if (jBin==jets4 && lcat==muon && variable=="KinBDT") {
      vector<Double_t> ret = {-1,-0.1,-0.06,-0.02,0.02,0.06,1};
      return ret;
    }
    else if (jBin==jets4 && lcat==muon && variable=="MEBDT") {
      vector<Double_t> ret = {-1,-0.04,0,0.04,0.08,1};
      return ret;
    }
    else if (jBin==jets4 && lcat==muon && variable=="KinMEBDT") {
      vector<Double_t> ret = {-1,-0.04,0,0.04,0.08,1};
      return ret;
    }
    */
    else cout << "ERROR in DEFS::getBinsX jBin="<<jBin<<" not found! Returning the default bin boundaries."<<endl;

    vector<Double_t> ret = {-1,-0.98,-0.96,-0.94,-0.92,-0.9,-0.88,-0.86,-0.84,-0.82,-0.8,-0.78,-0.76,-0.74,
                            -0.72,-0.7,-0.68,-0.66,-0.64,-0.62,-0.6,-0.58,-0.56,-0.54,-0.52,-0.5,-0.48,-0.46,
                            -0.44,-0.42,-0.4,-0.38,-0.36,-0.34,-0.32,-0.3,-0.28,-0.26,-0.24,-0.22,-0.2,-0.18,
                            -0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,
                            0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,
                            0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,
                            0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1};
    return ret;
  }//getBinsX

  //---------------------------------------------------------------------------
  // A routine that returns a string given the type
  string getTagCatString(TagCat type){

    if (type == pretag)  return "pretag";
    if (type == eq0tag)  return "eq0tag";
    if (type == eq1tag)  return "eq1tag";
    if (type == eq2tag)  return "eq2tag";
    if (type == ge0tag)  return "ge0tag";
    if (type == ge1tag)  return "ge1tag";
    if (type == ge2tag)  return "ge2tag";

    cout <<"ERROR  PhysicsProcess::getTagCatString cannot find the given type"<<endl;
    return "ERROR";

  }
  
  // A routine that returns a jetBin given a string
  DEFS::TagCat getTagCat(std::string str){
    
    if (str == "pretag")  return pretag;
    if (str == "eq0tag")  return eq0tag;
    if (str == "eq1tag")  return eq1tag;
    if (str == "eq2tag")  return eq2tag;
    if (str == "ge0tag")  return ge0tag;
    if (str == "ge1tag")  return ge1tag;
    if (str == "ge2tag")  return ge2tag;

    cout<<"ERROR  PhysicsProcess::getTagCat cannot find the given string"<<endl;
    return pretag;
    
  }

  //---------------------------------------------------------------------------
  string getCutLevelString(CutLevel type) {
     
    if (type == NPATtupleEvts) return "NPATtupleEvts";
    else if (type == c0)       return "c0:NJets";
    else if (type == c1)       return "c1:HLT+Kin";
    else if (type == c2)       return "c2:VtxCut";
    else if (type == c3)       return "c3:PrimaryEl/Mu";
    else if (type == c4)       return "c4:NotLooseMu";
    else if (type == c5)       return "c5:NotLooseEl";
    else if (type == c6)       return "c6:METE";
    else if (type == ac1)      return "ac1:LepFlavor";
    else if (type == ac2)      return "ac2:LepPt";
    else if (type == ac3)      return "ac3:MET";
    else if (type == ac4)      return "ac4:JetPt1";
    else if (type == ac5)      return "ac5:JetPt2";
    else if (type == ac6)      return "ac6:QCDPFIso";
/*
    else if (type == FNAL1)    return "FNAL:MET";
    else if (type == FNAL2)    return "FNAL:Mt";
    else if (type == FNAL3)    return "FNAL:LepPt";
    else if (type == FNAL4)    return "FNAL:LepEta";
    else if (type == FNAL5)    return "FNAL:JetPt1";
    else if (type == FNAL6)    return "FNAL:JetPt2";
    else if (type == FNAL7)    return "FNAL:JetEta1";
    else if (type == FNAL8)    return "FNAL:JetEta2";
    else if (type == FNAL9)    return "FNAL:drJetLep";
    else if (type == FNAL10)   return "FNAL:Mjj";
*/
    else if (type == BTag0)    return "BTag0";
    else if (type == BTag1)    return "BTag1";
    else if (type == BTag2)    return "BTag2";
    else if (type == BTag3p)   return "BTag3+";

    cout <<"ERROR  DEFS::getCutLevelString cannot find the given type"<<endl;
    return "ERROR";

   }

   DEFS::CutLevel getCutLevel(std::string str) {

    if (str == "NPATtupleEvts")        return NPATtupleEvts;
    else if (str == "c0:NJets")        return c0;
    else if (str == "c1:HLT+Kin")      return c1;
    else if (str == "c2:VtxCut")       return c2;
    else if (str == "c3:PrimaryEl/Mu") return c3;
    else if (str == "c4:NotLooseMu")   return c4;
    else if (str == "c5:NotLooseEl")   return c5;
    else if (str == "c6:METE")         return c6;
    else if (str == "ac1:LepFlavor")   return ac1;
    else if (str == "ac2:LepPt")       return ac2;
    else if (str == "ac3:MET")         return ac3;
    else if (str == "ac4:JetPt1")      return ac4;
    else if (str == "ac5:JetPt2")      return ac5;
    else if (str == "ac6:QCDPFIso")    return ac6;
/*    
    else if (str == "FNAL:MET")        return FNAL1;
    else if (str == "FNAL:Mt")         return FNAL2;
    else if (str == "FNAL:LepPt")      return FNAL3;
    else if (str == "FNAL:LepEta")     return FNAL4;
    else if (str == "FNAL:JetPt1")     return FNAL5;
    else if (str == "FNAL:JetPt2")     return FNAL6;
    else if (str == "FNAL:JetEta1")    return FNAL7;
    else if (str == "FNAL:JetEta2")    return FNAL8;
    else if (str == "FNAL:drJetLep")   return FNAL9;
    else if (str == "FNAL:Mjj")        return FNAL10;
*/
    else if (str == "BTag0")           return BTag0;
    else if (str == "BTag1")           return BTag1;
    else if (str == "BTag2")           return BTag2;
    else if (str == "BTag3+")          return BTag3p;

    cout<<"ERROR  DEFS::getCutLevel cannot find the given string"<<endl;
    return NPATtupleEvts;

   }
   
   //---------------------------------------------------------------------------
   // A routine that returns a string given the type
   std::string getControlRegionString(ControlRegion type){
      
      if (type == all)                return "all";
      else if (type == signal)        return "signal";
      else if (type == control1)      return "control1";
      else if (type == control2)      return "control2";
      else if (type == control3)      return "control3";
      else if (type == control4)      return "control4";
      else if (type == control5)      return "control5";
      else if (type == control6)      return "control6";
      else if (type == control7)      return "control7";
      else if (type == control8)      return "control8";
      else if (type == control9)      return "control9";
      else if (type == UVaCuts)       return "UVaCuts";
      else if (type == event)         return "event";
      else if (type == Diboson)       return "Diboson";
      else if (type == MVAEleID)      return "MVAEleID";
      else if (type == AntiMVAEleID)  return "AntiMVAEleID";
      else if (type == FlatMVAEleID)  return "FlatMVAEleID";
      else if (type == BDTBump)       return "BDTBump";
      else if (type == BDTAntiBump)   return "BDTAntiBump";
      else if (type == HailMary)      return "HailMary";
      else if (type == HailMaryLoose) return "HailMaryLoose";
      else if (type == None)          return "None";

      cout<<"ERROR DEFS::getControlregionString cannot find the given type"<<endl;
      return "ERROR";

   }

   // A routine that returns a ControlRegion given a string
   DEFS::ControlRegion getControlRegion(std::string str){
      
      if (str == "all")                return all;
      else if (str == "signal")        return signal;
      else if (str == "control1")      return control1;
      else if (str == "control2")      return control2;
      else if (str == "control3")      return control3;
      else if (str == "control4")      return control4;
      else if (str == "control5")      return control5;
      else if (str == "control6")      return control6;
      else if (str == "control7")      return control7;
      else if (str == "control8")      return control8;
      else if (str == "control9")      return control9;
      else if (str == "UVaCuts")       return UVaCuts;
      else if (str == "event")         return event;
      else if (str == "Diboson")       return Diboson;
      else if (str == "MVAEleID")      return MVAEleID;
      else if (str == "AntiMVAEleID")  return AntiMVAEleID;
      else if (str == "FlatMVAEleID")  return FlatMVAEleID;
      else if (str == "BDTBump")       return BDTBump;
      else if (str == "BDTAntiBump")   return BDTAntiBump;
      else if (str == "HailMary")      return HailMary;
      else if (str == "HailMaryLoose") return HailMaryLoose;
      else if (str == "None")          return None;

      cout<<"ERROR  DEFS::getControlRegion cannot find the given string"<<endl;
      return None;

   }

   //---------------------------------------------------------------------------
   // A routine that returns a string given the type
   std::string getNtupleTypeString(NtupleType type){
      
      if (type == EventNtuple)      return "EventNtuple";
      else if (type == METree)      return "METree";
      else if (type == MicroNtuple) return "MicroNtuple";
      else if (type == Other)       return "Other";

      cout<<"ERROR DEFS::getNtupleTypeString cannot find the given type"<<endl;
      return "ERROR";

   }

   // A routine that returns a NtupleType given a string
   DEFS::NtupleType getNtupleType(std::string str){
      
      if (str == "EventNtuple")      return EventNtuple;
      else if (str == "MEtree")      return METree;
      else if (str == "MicroNtuple") return MicroNtuple;
      else if (str == "Other")       return Other;

      cout<<"ERROR  DEFS::getNtupleType cannot find the given string"<<endl;
      return Other;

   }

  // A routine that returns the treename for a given NtupleType
  std::string getTreeName(NtupleType type, JetBin nJets){
     
     if (type == EventNtuple) {
        if(nJets == jets0)         return "PS/jets0";
        else if (nJets == jet1)    return "PS/jet1";
        else                       return "PS/jets2p";
     }
     else if (type == METree)      return "METree";
     else if (type == MicroNtuple) return "METree";
     else if (type == Other)       return "Other";
     
     cout<<"ERROR DEFS::getTreeName cannot find the given type"<<endl;
     return "ERROR";
     
  }

   //---------------------------------------------------------------------------
   // A routine that returns a string given the University
   std::string getUniversityString(University univ) {
      if (univ == TAMU)                 return "TAMU";
      else if (univ == UVa)             return "UVa";
      else if (univ == Combined)        return "Combined";
      else if (univ == OtherUniversity) return "OtherUniversity";

      cout<<"ERROR DEFS::getUniversityString cannot find the given university"<<endl;
      return "ERROR";

   }

   // A routine that returns a University given a string
   DEFS::University getUniversity(std::string str) {

      if (str == "TAMU")                 return TAMU;
      else if (str == "UVa")             return UVa;
      else if (str == "Combined")        return Combined;
      else if (str == "OtherUniversity") return OtherUniversity;

      cout<<"ERROR  DEFS::getUniversity cannot find the given string"<<endl;
      return OtherUniversity;

   }

   //---------------------------------------------------------------------------
   /// Print a PdgId as a named string
   std::string toParticleName(PdgId p) {
      if (p == ELECTRON)         return "ELECTRON";
      else if (p == POSITRON)    return "POSITRON";
      else if (p == MUON)        return "MUON";
      else if (p == ANTIMUON)    return "ANTIMUON";
      else if (p == TAU)         return "TAU";
      else if (p == ANTITAU)     return "ANTITAU";
      else if (p == NU_E)        return "NU_E";
      else if (p == NU_EBAR)     return "NU_EBAR";
      else if (p == NU_MU)       return "NU_MU";
      else if (p == NU_MUBAR)    return "NU_MUBAR";
      else if (p == NU_TAU)      return "NU_TAU";
      else if (p == NU_TAUBAR)   return "NU_TAUBAR";

      else if (p == d)           return "d";
      else if (p == u)           return "u";
      else if (p == s)           return "s";
      else if (p == c)           return "c";
      else if (p == b)           return "b";
      else if (p == t)           return "t";
      else if (p == dbar)        return "dbar";
      else if (p == ubar)        return "ubar";
      else if (p == sbar)        return "sbar";
      else if (p == cbar)        return "cbar";
      else if (p == bbar)        return "bbar";
      else if (p == tbar)        return "tbar";

      else if (p == DQUARK)      return "DQUARK";
      else if (p == UQUARK)      return "UQUARK";
      else if (p == SQUARK)      return "SQUARK";
      else if (p == CQUARK)      return "CQUARK";
      else if (p == BQUARK)      return "BQUARK";
      else if (p == TQUARK)      return "TQUARK";
      else if (p == ANTIDQUARK)  return "ANTIDQUARK";
      else if (p == ANTIUQUARK)  return "ANTIUQUARK";
      else if (p == ANTISQUARK)  return "ANTISQUARK";
      else if (p == ANTICQUARK)  return "ANTICQUARK";
      else if (p == ANTIBQUARK)  return "ANTIBQUARK";
      else if (p == ANTITQUARK)  return "ANTITQUARK";

      else if (p == PROTON)      return "PROTON";
      else if (p == ANTIPROTON)  return "ANTIPROTON";
      else if (p == NEUTRON)     return "NEUTRON";
      else if (p == ANTINEUTRON) return "ANTINEUTRON";
      else if (p == PIPLUS)      return "PIPLUS";
      else if (p == PIMINUS)     return "PIMINUS";

      else if (p == PHOTON)      return "PHOTON";
      else if (p == WPLUSBOSON)  return "WPLUSBOSON";
      else if (p == WMINUSBOSON) return "WMINUSBOSON";
      else if (p == ZBOSON)      return "ZBOSON";
      else if (p == GLUON)       return "GLUON";
      else if (p == HIGGS)       return "HIGGS";
      else if (p == ANY)         return "*";

      cout << "ERROR DEFS::toParticleName cannot find the given PdgId" << endl;
      return "*";
   }



   /// Print a PdgId as a named string
   DEFS::PdgId toParticleId(const std::string& pname) {
      if (pname == "ELECTRON")         return ELECTRON;
      else if (pname == "POSITRON")    return POSITRON;
      else if (pname == "MUON")        return MUON;
      else if (pname == "ANTIMUON")    return ANTIMUON;     
      else if (pname == "TAU")         return TAU;
      else if (pname == "ANTITAU")     return ANTITAU;
      else if (pname == "NU_E")        return NU_E;
      else if (pname == "NU_EBAR")     return NU_EBAR;
      else if (pname == "NU_MU")       return NU_MU;
      else if (pname == "NU_MUBAR")    return NU_MUBAR;
      else if (pname == "NU_TAU")      return NU_TAU;
      else if (pname == "NU_TAUBAR")   return NU_TAUBAR;

      else if (pname == "d")           return d;
      else if (pname == "u")           return u;
      else if (pname == "s")           return s;
      else if (pname == "c")           return c;
      else if (pname == "b")           return b;
      else if (pname == "t")           return t;
      else if (pname == "dbar")        return dbar;
      else if (pname == "ubar")        return ubar;
      else if (pname == "sbar")        return sbar;
      else if (pname == "cbar")        return cbar;
      else if (pname == "bbar")        return bbar;
      else if (pname == "tbar")        return tbar;

      else if (pname == "DQUARK")      return DQUARK;
      else if (pname == "UQUARK")      return UQUARK;
      else if (pname == "SQUARK")      return SQUARK;
      else if (pname == "CQUARK")      return CQUARK;
      else if (pname == "BQUARK")      return BQUARK;
      else if (pname == "TQUARK")      return TQUARK;
      else if (pname == "ANTIDQUARK")  return ANTIDQUARK;
      else if (pname == "ANTIUQUARK")  return ANTIUQUARK;
      else if (pname == "ANTISQUARK")  return ANTISQUARK;
      else if (pname == "ANTICQUARK")  return ANTICQUARK;
      else if (pname == "ANTIBQUARK")  return ANTIBQUARK;
      else if (pname == "ANTITQUARK")  return ANTITQUARK;

      else if (pname == "PROTON")      return PROTON;
      else if (pname == "ANTIPROTON")  return ANTIPROTON;
      else if (pname == "NEUTRON")     return NEUTRON;
      else if (pname == "ANTINEUTRON") return ANTINEUTRON;
      else if (pname == "PIPLUS")      return PIPLUS;
      else if (pname == "PIMINUS")     return PIMINUS;

      else if (pname == "PHOTON")      return PHOTON;
      else if (pname == "WPLUSBOSON")  return WPLUSBOSON;
      else if (pname == "WMINUSBOSON") return WMINUSBOSON;
      else if (pname == "ZBOSON")      return ZBOSON;
      else if (pname == "GLUON")       return GLUON;
      else if (pname == "HIGGS")       return HIGGS;
      else if (pname == "*")           return ANY;

    cout << "ERROR DEFS::toParticleId cannot find the given string" << endl;
    return ANY;

   }

   /// Return a PdgId from an int
   DEFS::PdgId getParticleIdFromInt(int p) {
      if (p == 11)         return ELECTRON;
      else if (p == -11)   return POSITRON;
      else if (p == 13)    return MUON;
      else if (p == -13)   return ANTIMUON;     
      else if (p == 15)    return TAU;
      else if (p == -15)   return ANTITAU;
      else if (p == 12)    return NU_E;
      else if (p == -12)   return NU_EBAR;
      else if (p == 14)    return NU_MU;
      else if (p == -14)   return NU_MUBAR;
      else if (p == 16)    return NU_TAU;
      else if (p == -16)   return NU_TAUBAR;

      else if (p == 1)     return d;
      else if (p == 2)     return u;
      else if (p == 3)     return s;
      else if (p == 4)     return c;
      else if (p == 5)     return b;
      else if (p == 6)     return t;
      else if (p == -1)    return dbar;
      else if (p == -2)    return ubar;
      else if (p == -3)    return sbar;
      else if (p == -4)    return cbar;
      else if (p == -5)    return bbar;
      else if (p == -6)    return tbar;

      else if (p == 1)     return DQUARK;
      else if (p == 2)     return UQUARK;
      else if (p == 3)     return SQUARK;
      else if (p == 4)     return CQUARK;
      else if (p == 5)     return BQUARK;
      else if (p == 6)     return TQUARK;
      else if (p == -1)    return ANTIDQUARK;
      else if (p == -2)    return ANTIUQUARK;
      else if (p == -3)    return ANTISQUARK;
      else if (p == -4)    return ANTICQUARK;
      else if (p == -5)    return ANTIBQUARK;
      else if (p == -6)    return ANTITQUARK;

      else if (p == 2212)  return PROTON;
      else if (p == -2212) return ANTIPROTON;
      else if (p == 2112)  return NEUTRON;
      else if (p == -2112) return ANTINEUTRON;
      else if (p == 211)   return PIPLUS;
      else if (p == -211)  return PIMINUS;

      else if (p == 22)    return PHOTON;
      else if (p == -22)   return PHOTON;
      else if (p == 24)    return WPLUSBOSON;
      else if (p == -24)   return WMINUSBOSON;
      else if (p == 23)    return ZBOSON;
      else if (p == -23)   return ZBOSON;
      else if (p == 21)    return GLUON;
      else if (p == -21)   return GLUON;
      else if (p == 25)    return HIGGS;
      else if (p == -25)   return HIGGS;
      else if (p == 10000) return ANY;

      cout << "ERROR DEFS::getParticleIdFromInt cannot find the given PdgId" << endl;
      return ANY;
   }

}// end of namespace DEFS
