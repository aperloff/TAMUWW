#include "TAMUWW/MatrixElement/interface/WHEventProb2Jet.hh"

#include <iostream>
#include <vector>

#include <cmath>

#include "TAMUWW/MatrixElement/interface/DHELASWrapper.hh"
#include "TAMUWW/MatrixElement/interface/MEConstants.hh"
#include "TAMUWW/MatrixElement/interface/PartonColl.hh"
#include "TAMUWW/MatrixElement/interface/TransferFunction.hh"

using std::vector;

extern "C"
{
   void* mywh_(double[][4], double*, const double*, double*);
   void* mywh2_(double[][4], double*, const double*, double*);
}

WHEventProb2Jet::WHEventProb2Jet(Integrator& integrator,
                                 const TransferFunction& tf,
                                 double higgsMass) :
  EventProb2Jet("WH", integrator, 3, 1, tf){ 
  setHiggsMassAndWidth(higgsMass);
}

// ------------------------------------------------------------------
// This method sets the Higgs mass and the Higgs width
void  WHEventProb2Jet::setHiggsMassAndWidth(double mHiggs) {

  // Set the mass
  m_massHiggs = mHiggs;

  // Use the theoretical Higgs width for the given mass 
  // multiplied by a factor of 100
  m_widthHiggs = 100.0 * calcHiggsWidth(mHiggs);

}//setHiggsMassAndWidth


// ------------------------------------------------------------------
void WHEventProb2Jet::getPeaks(VecVecDouble& answer, const double bounds[]) const
{
   using PeterFunctions::Math::square;

   PartonColl temp(*getMeasuredColl());
   temp.setMet();

   const TLorentzVector& jet1 = temp.getJet(0);
   const TLorentzVector& jet2 = temp.getJet(1);

   vector<double> lower, upper;
   for (unsigned i = 0; i < getDimension(); ++i)
   {
      lower.push_back(bounds[2 * i]);
      upper.push_back(bounds[2 * i + 1]);
   }
   adjustBounds(lower);
   adjustBounds(upper);

   if (jet1.E() < lower[1] || jet1.E() > upper[1])
      return;
   if (jet2.E() < lower[2] || jet2.E() > upper[2])
      return;

   
   TLorentzVector& lepton = temp.getLepton();
   TLorentzVector neutrino = temp.getNeutrino();

   const double mW = MEConstants::wMass;


   double a1 = lepton.Px() * neutrino.Px() + lepton.Py() * neutrino.Py();
   double b1 = square(mW) / 2;
   double c1 = lepton.Pz();
   double d1 = lepton.E();
   double e1 = square(lepton.Px()) + square(lepton.Py());

   double a = square(d1) - square(c1);
   double b = -2 * c1 * (a1 + b1);
   double c = e1 * square(d1) - square(a1) - square(b1) - 2 * a1 * b1;

   double discriminant = square(b) - 4 * a * c;

   vector<double> pz;
   if (discriminant < 0)
   {
      pz.push_back(-b / (2 * a));
   }
   else
   {
      pz.push_back((-b + sqrt(discriminant)) / (2 * a));
      pz.push_back((-b - sqrt(discriminant)) / (2 * a));
   }

   TLorentzVector totalJets = jet1 + jet2;
   double totalE = std::sqrt(square(m_massHiggs) + square(totalJets.P()));

//   std::cerr << "totalE " << totalE << std::endl;

   for (vector<double>::const_iterator it = pz.begin(); it != pz.end(); ++it)
   {
      if (*it > lower[0] && *it < upper[0])
      {
         double jet1Epeak[2] = {jet1.E(), totalE - jet2.E()};
         double jet2Epeak[2] = {totalE - jet1.E(), jet2.E()};
         for (unsigned i = 0; i < 2; ++i)
         {
//            std::cerr << "jets: " << jet1Epeak[0] << " " << jet2Epeak[1] << "   " << lower[1] << " " << upper[1] << " " << lower[2] << " " << upper[2] << std::endl;

            if (jet1Epeak[i] > lower[1] && jet1Epeak[i] < upper[1]
                && jet2Epeak[i] > lower[2] && jet2Epeak[i] < upper[2])
            {
               double phaseSpaceCoord[3] = {*it, jet1Epeak[i], jet2Epeak[i]};
               vector<double> phaseSpaceVec(phaseSpaceCoord,
                                            phaseSpaceCoord + 3);
               answer.push_back(phaseSpaceVec);
               std::cerr << "Returning peak: " << phaseSpaceCoord[0] << " "
                         << phaseSpaceCoord[1] << " " << phaseSpaceCoord[2]
                         << std::endl;
            }
         }
      }
   }

//   exit(1);
   return;
}

void WHEventProb2Jet::setupIntegral()
{
   std::cout << "\tHiggs mass: " << m_massHiggs
	     <<"  width: "       << m_widthHiggs <<std::endl;
}

double WHEventProb2Jet::matrixElement() const
{
   typedef SimpleArray<DHELAS::HelArray, 1> Array1;
   typedef SimpleArray<DHELAS::HelArray, 2> Array2;

   using MEConstants::bMass;
   using MEConstants::wMass;
   using MEConstants::wWidth;

   const PartonColl* partons = getPartonColl();

   double answer = 0;

   enum {vecSize = 4};
   typedef SimpleArray<doublecomplex, vecSize> OutputType;

   doublecomplex factor[2] = {doublecomplex(MEConstants::gwf, 0),
                              doublecomplex(0, 0)};

   doublecomplex factorGWWH[1] = {doublecomplex(MEConstants::gwwh, 0)};

   doublecomplex factorGHBOT[2] = {doublecomplex(MEConstants::ghbot, 0),
                                   doublecomplex(MEConstants::ghbot, 0)};

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
      
      Array1 vec1 = DHELAS::ixxxxx<1>(partons->getParton1(), 0, 1);
      Array1 vec2 = DHELAS::oxxxxx<1>(partons->getParton2(), 0, -1);
//      Array1 vec3 = DHELAS::ixxxxx<1>(partons->getLepton(), 0, -1);
      Array1 vec4 = DHELAS::oxxxxx<1>(partons->getNeutrino(), 0, 1);
      Array2 vec5 = DHELAS::oxxxxx<2>(partons->getJet(0), bMass, 1);
      Array2 vec6 = DHELAS::ixxxxx<2>(partons->getJet(1), bMass, -1);

      Array1 vec7 = DHELAS::jioxxx(vec1, vec2, factor, wMass, wWidth);
      Array1 vec8 = DHELAS::jioxxx(vec3, vec4, factor, wMass, wWidth);
      Array1 vec9 = DHELAS::hvvxxx(vec7, vec8, factorGWWH, m_massHiggs,
                                   m_widthHiggs);

//      std::cout << "vec6: " << vec6[0] << std::endl;

      OutputType output = DHELAS::iosxxx(vec6, vec5, vec9, factorGHBOT);

      for (unsigned i = 0; i < vecSize; ++i)
      {
         answer += std::norm(-output[i]) * 9;
      }
   }
   else
   {
      // Calculate the lepton only once per integration
      static Array1 vec3;
      static double lepE = 0;
      if (lepE != partons->getLepton().E())
      {
         vec3 = DHELAS::oxxxxx<1>(partons->getLepton(), 0, 1);
         lepE = partons->getLepton().E();
      }

      Array1 vec1 = DHELAS::ixxxxx<1>(partons->getParton1(), 0, 1);
      Array1 vec2 = DHELAS::oxxxxx<1>(partons->getParton2(), 0, -1);
//      Array1 vec3 = DHELAS::oxxxxx<1>(partons->getLepton(), 0, 1);
      Array1 vec4 = DHELAS::ixxxxx<1>(partons->getNeutrino(), 0, -1);
      Array2 vec5 = DHELAS::ixxxxx<2>(partons->getJet(0), bMass, -1);
      Array2 vec6 = DHELAS::oxxxxx<2>(partons->getJet(1), bMass, 1);

      Array1 vec7 = DHELAS::jioxxx(vec1, vec2, factor, wMass, wWidth);
      Array1 vec8 = DHELAS::jioxxx(vec4, vec3, factor, wMass, wWidth);
      Array1 vec9 = DHELAS::hvvxxx(vec8, vec7, factorGWWH, m_massHiggs,
                                   m_widthHiggs);

      OutputType output = DHELAS::iosxxx(vec5, vec6, vec9, factorGHBOT);

      for (unsigned i = 0; i < vecSize; ++i)
      {
         answer += std::norm(-output[i]) * 9;
      }
   }

   answer /= 36;

   //For debugging
   //Fortran array for slow method
   
   //std::cerr << "New answer: " << answer << std::endl;
   
   //double fortranArray[6][4];
   //makeFortranArray(fortranArray);

   //const double wMass = MEConstants::wMass;

   //if (partons->getLepCharge() > 0)
   //   mywh_(fortranArray, &higgsMass, &wMass, &answer);
   //else
   //   mywh2_(fortranArray, &higgsMass, &wMass, &answer);
   //std::cerr << "slow matrix element: " << answer << std::endl;
   //exit(1);
   

   return answer;

}


void WHEventProb2Jet::setPartonTypes() const
{
   if (getMeasuredColl()->getLepCharge() > 0)
   {
      getMeasuredColl()->setParton1Type(kUp);
      getMeasuredColl()->setParton2Type(kDown);
   }
   else
   {
      getMeasuredColl()->setParton1Type(kDown);
      getMeasuredColl()->setParton2Type(kUp);
   }
}

void WHEventProb2Jet::getScale(double& scale1, double&scale2) const
{
   double scale = getPartonColl()->sHat();
   if (scale < 0)
      scale1 = scale2 = 0;
   else
      scale1 = scale2 = std::sqrt(scale);
}
