///////////////////////////////////////////////////////////////////////
//// Author: Ilya Osipenkov Texas A&M University.
//// The diagrams can be compared with madgraph using the processes:
//// Two b quark final state: u d~ -> e+ ve b b~ (d~ u -> e+ ve b b~) for W+, d u~ -> e- ve~ b b~ (u~ d -> e- ve~ b b~) for W-
//////////////////////////////////////////////////////////////////////
#include "TAMUWW/MatrixElement/interface/WZtobbEventProb2Jet.hh"

#include <iostream>

#include <cmath>
#include <vector>

#include "TAMUWW/MatrixElement/interface/DHELASWrapper.hh"
#include "TAMUWW/MatrixElement/interface/MEConstants.hh"
#include "TAMUWW/MatrixElement/interface/PartonColl.hh"
#include "TAMUWW/MatrixElement/interface/TransferFunction.hh"

//#define MADGRAPH_TEST
using std::vector;
// using std::cout;
// using std::endl;

#ifdef MADGRAPH_TEST
extern "C"
{
  void* wpzbbm_(double[][4], const double*, double*);// lepQ>0, u d~ -> e+ ve b b~ 
  void* wpzbbaltm_(double[][4], const double*, double*);// lepQ>0, d~ u -> e+ ve b b~ instead
  void* wmzbbm_(double[][4], const double*, double*);// lepQ<0, d u~ -> e- ve~ b b~
  void* wmzbbaltm_(double[][4], const double*, double*);// lepQ<0, u~ d -> e- ve~ b b~ instead
}
#endif

// ------------------------------------------------------------------
WZtobbEventProb2Jet::WZtobbEventProb2Jet(Integrator& integrator,
                                   const TransferFunction& tf) :
  EventProb2Jet(DEFS::EP::WZbb, integrator, 3, 4, tf), 
  swapPartonMom(false), alphas_process(0.13) //Take the alphas_process value from MadGraph or use MEConstants::alphas
{}

// ------------------------------------------------------------------


// ------------------------------------------------------------------
void WZtobbEventProb2Jet::changeVars(const vector<double>& parameters)
{
   TLorentzVector& jet1 = getPartonColl()->getJet(0);
   TLorentzVector& jet2 = getPartonColl()->getJet(1);

   jet1.SetRho(parameters[1]);
   jet1.SetE(parameters[1]);
   jet2.SetRho(parameters[2]);
   jet2.SetE(parameters[2]);

   getPartonColl()->getNeutrino().SetPz(parameters[0]);
}

// ------------------------------------------------------------------


double WZtobbEventProb2Jet::matrixElement() const
{
  double answer = 0;
   typedef SimpleArray<DHELAS::HelArray, 1> Array1;
   typedef SimpleArray<DHELAS::HelArray, 2> Array2;
   typedef SimpleArray<DHELAS::HelArray, 4> Array4;

   using MEConstants::zMass;
   using MEConstants::wMass;
   using MEConstants::wWidth;
   using MEConstants::zWidth;
   using MEConstants::bMass;

   //   MEConstants::PrintAllConstants(alphas_process);

   const PartonColl* partons = getPartonColl();


   doublecomplex factorGWF[2]   = {doublecomplex(MEConstants::gwf, 0),
				   doublecomplex(0, 0)};
   doublecomplex zFactorU[2] = {doublecomplex(-MEConstants::gzu1, 0),
                                doublecomplex(-MEConstants::gzu2, 0)};
   doublecomplex zFactorD[2] = {doublecomplex(-MEConstants::gzd1, 0),
                                doublecomplex(-MEConstants::gzd2, 0)};
   doublecomplex factorZWW[2] = {doublecomplex(-MEConstants::gwwz, 0),
                               doublecomplex(0, 0)};

   enum {vecSize = 4};
   typedef SimpleArray<doublecomplex, vecSize> OutputType;

   if (partons->getLepCharge() > 0)
   {
      // Calculate the lepton only once per integration
      static Array1 vec3;
      static double lepE = 0;
      if (lepE != partons->getLepton().E())
      {
         vec3 = DHELAS::ixxxxx<1>(partons->getLepton(), 0, -1);
         lepE = partons->getLepton().E();
      }
      Array1 vec1;
      Array1 vec2;
      if ( !swapPartonMom ) {
	vec1 = DHELAS::ixxxxx<1>(partons->getParton1(), 0, +1);
        vec2 = DHELAS::oxxxxx<1>(partons->getParton2(), 0, -1);
      } else {
	vec1 = DHELAS::ixxxxx<1>(partons->getParton2(), 0, +1);
        vec2 = DHELAS::oxxxxx<1>(partons->getParton1(), 0, -1);
      }
      Array1 vec4 = DHELAS::oxxxxx<1>(partons->getNeutrino(), 0, 1);
      Array2 vec5 = DHELAS::oxxxxx<2>(partons->getJet(0), bMass, +1);
      Array2 vec6 = DHELAS::ixxxxx<2>(partons->getJet(1), bMass, -1);

      Array1 vec7 = DHELAS::jioxxx(vec3, vec4, factorGWF, wMass, wWidth);
      Array1 vec8 = DHELAS::fvoxxx(vec2, vec7, factorGWF, 0, 0);
      Array1 vec9 = DHELAS::jioxxx(vec1, vec8, zFactorU, zMass, zWidth);
      OutputType output1 = DHELAS::iovxxx(vec6, vec5, vec9, zFactorD);

      Array1 vec10 = DHELAS::fvixxx(vec1, vec7, factorGWF, 0, 0);
      Array1 vec11 = DHELAS::jioxxx(vec10, vec2, zFactorD, zMass, zWidth);
      OutputType output2 = DHELAS::iovxxx(vec6, vec5, vec11, zFactorD);

      Array1 vec12 = DHELAS::jioxxx(vec1, vec2, factorGWF, wMass, wWidth);
      Array1 vec13 = DHELAS::jvvxxx(vec12, vec7, factorZWW, zMass, zWidth);
      OutputType output3 = DHELAS::iovxxx(vec6, vec5, vec13, zFactorD);

      for (unsigned i = 0; i < vecSize; ++i)
      {
	doublecomplex temp1 = -output1[i]-output2[i]-output3[i];
	double m1 = 9.0*std::norm(temp1);

	answer+= m1;
	//cout << "current helicity 'amplitude'" << m1 << endl;
      }

   }
   else
   {
     // Calculate the lepton only once per integration
     static Array1 vec3;
     static double lepE = 0;
     if (lepE != partons->getLepton().E())
       {
         vec3 = DHELAS::oxxxxx<1>(partons->getLepton(), 0, +1);
         lepE = partons->getLepton().E();
       }
     Array1 vec1;
     Array1 vec2;
     if ( !swapPartonMom ) {
       vec1 = DHELAS::ixxxxx<1>(partons->getParton1(), 0, +1);
       vec2 = DHELAS::oxxxxx<1>(partons->getParton2(), 0, -1);
     } else {
       vec1 = DHELAS::ixxxxx<1>(partons->getParton2(), 0, +1);
       vec2 = DHELAS::oxxxxx<1>(partons->getParton1(), 0, -1);
     }
     Array1 vec4 = DHELAS::ixxxxx<1>(partons->getNeutrino(), 0, -1);
     Array2 vec5 = DHELAS::oxxxxx<2>(partons->getJet(0), bMass, +1);
     Array2 vec6 = DHELAS::ixxxxx<2>(partons->getJet(1), bMass, -1);
     
     Array1 vec7 = DHELAS::jioxxx(vec4, vec3, factorGWF, wMass, wWidth);
     Array1 vec8 = DHELAS::fvoxxx(vec2, vec7, factorGWF, 0, 0);
     Array1 vec9 = DHELAS::jioxxx(vec1, vec8, zFactorD, zMass, zWidth);
     OutputType output1 = DHELAS::iovxxx(vec6, vec5, vec9, zFactorD);

     Array1 vec10 = DHELAS::fvixxx(vec1, vec7, factorGWF, 0, 0);
     Array1 vec11 = DHELAS::jioxxx(vec10, vec2, zFactorU, zMass, zWidth);
     OutputType output2 = DHELAS::iovxxx(vec6, vec5, vec11, zFactorD);

     Array1 vec12 = DHELAS::jioxxx(vec1, vec2, factorGWF, wMass, wWidth);
     Array1 vec13 = DHELAS::jvvxxx(vec7, vec12, factorZWW, zMass, zWidth);
     OutputType output3 = DHELAS::iovxxx(vec6, vec5, vec13, zFactorD);

      for (unsigned i = 0; i < vecSize; ++i)
      {
	doublecomplex temp1 = output1[i]+output2[i]+output3[i];
	double m1 = 9.0*std::norm(temp1);

	answer+= m1;
	//cout << "current helicity 'amplitude'" << m1 << endl;
      }

   }

   answer /= 36; //relative weight
//   std::cerr << "New Answer: " << answer << std::endl;
#ifdef MADGRAPH_TEST
  // -----------------------------------------------------------
  // This code reports our answer as well as the madgraph answer
  // -----------------------------------------------------------
  // Report our answers
  cout<<" My answer= "<<answer<<endl;

  // Make a fortran array, format which is needed for the fortran calls
  double fortranArray[6][4];
  makeFortranArray(fortranArray);
   
  // Evalute the matrix element according to madgraph
  //double mhiggs = m_massHiggs; // to get rid of the const identifier
  double mw = wMass; // to get rid of the const identifier
  //  double mz = zMass;
  double an = 0;
  if (partons->getLepCharge() > 0) {
    if ( !swapPartonMom ) {
      wpzbbm_(fortranArray, &mw, &an);
    } else {
      wpzbbaltm_(fortranArray, &mw, &an);
    }
  } else {
    if ( !swapPartonMom ) {
      wmzbbm_(fortranArray, &mw, &an);
    } else {
      wmzbbaltm_(fortranArray, &mw, &an);
    }
  }

  cout << "Madgraph answer= " << an << endl;
   
  // Exit right away
  exit(1);  

#endif

  return answer;
  
}

void WZtobbEventProb2Jet::setPartonTypes() const
{
   if (getMeasuredColl()->getLepCharge() > 0) {
     if ( !swapPartonMom ) {
       getMeasuredColl()->setParton1Type(kUp);
       getMeasuredColl()->setParton2Type(kAntiDown);
     } else {
       getMeasuredColl()->setParton1Type(kAntiDown);
       getMeasuredColl()->setParton2Type(kUp);
     }
   }
   else
   {
     if ( !swapPartonMom ) {
       getMeasuredColl()->setParton1Type(kDown);
       getMeasuredColl()->setParton2Type(kAntiUp);
     } else {
       getMeasuredColl()->setParton1Type(kAntiUp);
       getMeasuredColl()->setParton2Type(kDown);
     }
   }

}

void WZtobbEventProb2Jet::setJetTypes()
{
  m_JetType[0]=kBottom;
  m_JetType[1]=kAntiBottom;

  if ( getSwappedJet0Jet1Status() ) {
    int tempType=m_JetType[0];
    m_JetType[0]=m_JetType[1];
    m_JetType[1]=tempType;
  }
}

void WZtobbEventProb2Jet::getScale(double& scale1, double& scale2) const
{
   double scale = getPartonColl()->sHat();
   if (scale < 0)
      scale1 = scale2 = 0;
   else
      scale1 = scale2 = std::sqrt(scale);
}

bool WZtobbEventProb2Jet::onSwitch()
{
  
  switch (getLoop()) {
  case 0:
    swapPartonMom=false;
    setSwapJet0Jet1Status(false);
    //swapPartonMom=true; //when testing alternate functions
    break;
  case 1:
    swapJets(0, 1);
    setSwapJet0Jet1Status(true);
    break;
  case 2:
    swapPartonMom=true;
    break;
  case 3:
    swapJets(0, 1);
    setSwapJet0Jet1Status(false);
    break;
  default:
    return false;
  }

  return true;

}
