/**
 * \class CorrThreeBodyCondition
 *
 * Description: evaluation of a three-body correlation condition (= three-muon invariant mass). 
 *                                                                                                                                                                 
 * Implementation:                                                                                                                                                           
 *    <TODO: enter implementation details>                                                                                                                 
 *                                                                                                                                     
 * \author: Elisa Fontanesi - Boston University                                                                                                                                                            
 * Starting from CorrelationTemplate.h written by Vasile Mihai Ghete - HEPHY Vienna      
 *
 */

// this class header
#include "L1Trigger/L1TGlobal/interface/CorrCondition.h"
#include "L1Trigger/L1TGlobal/interface/CorrThreeBodyCondition.h"

// system include files
#include <iostream>
#include <iomanip>

#include <string>
#include <vector>
#include <algorithm>

// user include files
//   base classes
#include "L1Trigger/L1TGlobal/interface/CorrelationTemplate.h"
#include "L1Trigger/L1TGlobal/interface/CorrelationThreeBodyTemplate.h"
#include "L1Trigger/L1TGlobal/interface/ConditionEvaluation.h"

#include "L1Trigger/L1TGlobal/interface/MuCondition.h"
#include "L1Trigger/L1TGlobal/interface/MuonTemplate.h"
#include "L1Trigger/L1TGlobal/interface/GlobalScales.h"
#include "L1Trigger/L1TGlobal/interface/GlobalBoard.h"

#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

// constructors
//     default
l1t::CorrThreeBodyCondition::CorrThreeBodyCondition() : ConditionEvaluation() {}

//     from base template condition (from event setup usually)
l1t::CorrThreeBodyCondition::CorrThreeBodyCondition(const GlobalCondition* corrThreeBodyTemplate,
                                  const GlobalCondition* cond0Condition,
                                  const GlobalCondition* cond1Condition,
                                  const GlobalCondition* cond2Condition,
                                  const GlobalBoard* ptrGTB)
    : ConditionEvaluation(),
      m_gtCorrelationThreeBodyTemplate(static_cast<const CorrelationThreeBodyTemplate*>(corrThreeBodyTemplate)),
      m_gtCond0(cond0Condition),
      m_gtCond1(cond1Condition),
      m_gtCond2(cond2Condition),
      m_uGtB(ptrGTB) {}

// copy constructor
void l1t::CorrThreeBodyCondition::copy(const l1t::CorrThreeBodyCondition& cp) {
  m_gtCorrelationThreeBodyTemplate = cp.gtCorrelationThreeBodyTemplate();
  m_uGtB = cp.getuGtB();

  m_condMaxNumberObjects = cp.condMaxNumberObjects();
  m_condLastResult = cp.condLastResult();
  m_combinationsInCond = cp.getCombinationsInCond();

  m_verbosity = cp.m_verbosity;
}

l1t::CorrThreeBodyCondition::CorrThreeBodyCondition(const l1t::CorrThreeBodyCondition& cp) : ConditionEvaluation() { copy(cp); }

// destructor
l1t::CorrThreeBodyCondition::~CorrThreeBodyCondition() {
  // empty
}

// equal operator
l1t::CorrThreeBodyCondition& l1t::CorrThreeBodyCondition::operator=(const l1t::CorrThreeBodyCondition& cp) {
  copy(cp);
  return *this;
}

///   set the pointer to uGT GlobalBoard
void l1t::CorrThreeBodyCondition::setuGtB(const GlobalBoard* ptrGTB) { m_uGtB = ptrGTB; }

void l1t::CorrThreeBodyCondition::setScales(const GlobalScales* sc) { m_gtScales = sc; }

// try all object permutations and check their correlations, if required
const bool l1t::CorrThreeBodyCondition::evaluateCondition(const int bxEval) const {
  std::ostringstream myCout; //EF
  m_gtCorrelationThreeBodyTemplate->print(myCout);  //EF
  LogDebug("L1TGlobal") << "Three-body Correlation Condition Evaluation..." << std::endl; 

  bool condResult = false;
  bool reqObjResult = false;

  // number of objects in condition (it is 3, no need to retrieve from condition template) and their type
  int nObjInCond = 3;
  std::vector<GlobalObject> cndObjTypeVec(nObjInCond);

  // evaluate first the three subconditions (Type1s)
  const GtConditionCategory cond0Categ = m_gtCorrelationThreeBodyTemplate->cond0Category();
  const GtConditionCategory cond1Categ = m_gtCorrelationThreeBodyTemplate->cond1Category();
  const GtConditionCategory cond2Categ = m_gtCorrelationThreeBodyTemplate->cond2Category();

  const MuonTemplate* corrMuon = nullptr;

  // FIXME copying is slow...
  CombinationsInCond cond0Comb;
  CombinationsInCond cond1Comb;
  CombinationsInCond cond2Comb;

  int cond0bx(0);
  int cond1bx(0);
  int cond2bx(0);

  // FIRST OBJECT
  if (cond0Categ == CondMuon) {
    std::cout << "\n --------------------- First muon checks ---------------------" << std::endl; 
    corrMuon = static_cast<const MuonTemplate*>(m_gtCond0);
    MuCondition muCondition(
			    corrMuon, m_uGtB, 0, 0);
 
    muCondition.evaluateConditionStoreResult(bxEval);
    reqObjResult = muCondition.condLastResult();

    cond0Comb = (muCondition.getCombinationsInCond());
    cond0bx = bxEval + (corrMuon->condRelativeBx());
    cndObjTypeVec[0] = (corrMuon->objectType())[0];
    
    //EF if (m_verbosity) {
    std::ostringstream myCout;
    muCondition.print(myCout);
    std::cout << myCout.str() << std::endl;
      //LogDebug("L1TGlobal") << myCout.str() << std::endl;
      //EF }
  }
  
  else {
    // Interested only in three-muon correlations
    LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 0" << std::endl;
    return false;
  } 
  
  // return if first subcondition is false
  if (!reqObjResult) {
    LogDebug("L1TGlobal") << "\n  First subcondition false, second subcondition not evaluated and not printed."
                          << std::endl;
    return false;
  }
  
  // SECOND OBJECT
  reqObjResult = false;

  //switch (cond1Categ) {
  if (cond1Categ == CondMuon) {
    std::cout << "\n --------------------- Second muon checks ---------------------" << std::endl; 
    corrMuon = static_cast<const MuonTemplate*>(m_gtCond1);
    MuCondition muCondition(
			    corrMuon, m_uGtB, 0, 0);
    
    muCondition.evaluateConditionStoreResult(bxEval);
    reqObjResult = muCondition.condLastResult();

    cond1Comb = (muCondition.getCombinationsInCond());
    cond1bx = bxEval + (corrMuon->condRelativeBx());
    cndObjTypeVec[1] = (corrMuon->objectType())[0];
    
    //EF if (m_verbosity) {
    std::ostringstream myCout;
    muCondition.print(myCout);
    std::cout << myCout.str() << std::endl;
    //LogDebug("L1TGlobal") << myCout.str() << std::endl;
    //EF }
  }
  
  else {
    // Interested only in three-muon correlations
    LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 1" << std::endl;
    return false;
  }
  
  // return if second subcondition is false
  if (!reqObjResult) {
    LogDebug("L1TGlobal") << "\n  Second subcondition false, third subcondition not evaluated and not printed."
                          << std::endl;
    return false;
  }

  // THIRD OBJECT
  reqObjResult = false;

  if (cond2Categ == CondMuon) {
    std::cout << "\n --------------------- Third muon checks ---------------------" << std::endl; 
    corrMuon = static_cast<const MuonTemplate*>(m_gtCond2);
    MuCondition muCondition(
			    corrMuon, m_uGtB, 0, 0);
 
    muCondition.evaluateConditionStoreResult(bxEval);
    reqObjResult = muCondition.condLastResult();

    cond2Comb = (muCondition.getCombinationsInCond());
    cond2bx = bxEval + (corrMuon->condRelativeBx()); 
    cndObjTypeVec[2] = (corrMuon->objectType())[0];    

    //EF if (m_verbosity) {
    std::ostringstream myCout;
    muCondition.print(myCout);
    std::cout << myCout.str() << std::endl;
    //LogDebug("L1TGlobal") << myCout.str() << std::endl;
    //EF }
    
  }
  
  else {
    // Interested only in three-muon correlations
    LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 2" << std::endl;
    return false;
  } 
  
  // return if third subcondition is false
  if (!reqObjResult) {
    return false;
  } else {
    std::cout << "\n"
      //LogDebug("L1TGlobal") << "\n"
                          << "Found three objects satisfying subconditions: evaluate three-body correlation requirements.\n"
                          << std::endl;
  }

  // since we have three good legs get the correlation parameters
  CorrelationThreeBodyTemplate::CorrelationThreeBodyParameter corrPar = *(m_gtCorrelationThreeBodyTemplate->correlationThreeBodyParameter());

  //EF: Vector to store muon indexes                                                                                                              
  vector<int> muIndexes_01;                                                                                                                                
  vector<int> muIndexes_02;                                                                                                                                
  vector<int> muIndexes_12;                                                                                                                                
  //EF: Vector to store dimuon invariant masses 
  vector<long long> dimuInvMass_01;                                                                                                         
  vector<long long> dimuInvMass_02;                                                                                                         
  vector<long long> dimuInvMass_12;                                                                                                         
  muIndexes_01.clear();                                                                                                      
  muIndexes_02.clear();                                                                                                      
  muIndexes_12.clear();                                          
  dimuInvMass_01.clear();                                                                                                                
  dimuInvMass_02.clear();                                                                                                                
  dimuInvMass_12.clear();                                                                                                                
 
  // vector to store the indices of the objects from the combination evaluated in the condition
  SingleCombInCond objectsInComb;
  objectsInComb.reserve(nObjInCond);
  LogDebug("L1TGlobal") << "\n"
			<< "Number of objects considered in the condition is nObjInCond=" << nObjInCond << std::endl; 

  // clear the m_combinationsInCond vector
  (combinationsInCond()).clear();

  // pointers to objects
  const BXVector<const l1t::Muon*>* candMuVec = nullptr;

  // make the conversions of the indices, depending on the combination of objects involved
  // (via pair index)

  int phiIndex0 = 0;
  double phi0Phy = 0.;
  int phiIndex1 = 0;
  double phi1Phy = 0.;
  int phiIndex2 = 0;
  double phi2Phy = 0.;

  int etaIndex0 = 0;
  double eta0Phy = 0.;
  int etaBin0 = 0;
  int etaIndex1 = 0;
  double eta1Phy = 0.;
  int etaBin1 = 0;
  int etaIndex2 = 0;
  double eta2Phy = 0.;
  int etaBin2 = 0;

  int etIndex0 = 0;
  int etBin0 = 0;
  double et0Phy = 0.;
  int etIndex1 = 0;
  int etBin1 = 0;
  double et1Phy = 0.;
  int etIndex2 = 0;
  int etBin2 = 0;
  double et2Phy = 0.;

  // Determine the number of phi bins to get cutoff at pi
  int phiBound = 0;
  const GlobalScales::ScaleParameters& par = m_gtScales->getMUScales();
  phiBound = (int)((par.phiMax - par.phiMin) / par.phiStep) / 2;
  LogDebug("L1TGlobal") << "Phi Bound = " << phiBound << std::endl;

  // Keep track of objects for LUTS
  std::string lutObj0 = "NULL";
  std::string lutObj1 = "NULL";
  std::string lutObj2 = "NULL";

  // EF
  std::cout << "  Subcondition 0: std::vector<SingleCombInCond> size: " << (cond0Comb.size()) << std::endl;
  std::cout << "  Subcondition 1: std::vector<SingleCombInCond> size: " << (cond1Comb.size()) << std::endl;
  std::cout << "  Subcondition 2: std::vector<SingleCombInCond> size: " << (cond2Comb.size()) << std::endl;
  //LogTrace("L1TGlobal") << "  Subcondition 0: std::vector<SingleCombInCond> size: " << (cond0Comb.size()) << std::endl;
  //LogTrace("L1TGlobal") << "  Subcondition 1: std::vector<SingleCombInCond> size: " << (cond1Comb.size()) << std::endl;
  //LogTrace("L1TGlobal") << "  Subcondition 2: std::vector<SingleCombInCond> size: " << (cond2Comb.size()) << std::endl;


  ////////////////////////////////
  // LOOP OVER ALL COMBINATIONS //
  ////////////////////////////////
  // BLW: Optimization issue: potentially making the same comparison twice
  //                          if both legs are the same object type.
  unsigned int preShift = 0;
  // *** FIRST PAIR: 0-1
  for (std::vector<SingleCombInCond>::const_iterator it0Comb = cond0Comb.begin(); it0Comb != cond0Comb.end();
       it0Comb++) {
    LogDebug("L1TGlobal") << "Looking at the first subcondition of the pair 0-1" << std::endl;
    // Type1s: there is 1 object only, no need for a loop, index 0 should be OK in (*it0Comb)[0]
    // ... but add protection to not crash
    int obj0Index = -1;

    if (!(*it0Comb).empty()) {
      obj0Index = (*it0Comb)[0];
    } else {
      LogTrace("L1TGlobal") << "\n  SingleCombInCond (*it0Comb).size() " << ((*it0Comb).size()) << std::endl;
      return false;
    }

    // FIRST OBJECT: Collect the information on the first leg of the correlation
    if (cond0Categ == CondMuon) {
        lutObj0 = "MU";
        candMuVec = m_uGtB->getCandL1Mu();
        phiIndex0 = (candMuVec->at(cond0bx, obj0Index))->hwPhiAtVtx();  //(*candMuVec)[obj0Index]->phiIndex();
        etaIndex0 = (candMuVec->at(cond0bx, obj0Index))->hwEtaAtVtx();
        etIndex0 = (candMuVec->at(cond0bx, obj0Index))->hwPt();
        int etaBin0 = etaIndex0;
        if (etaBin0 < 0) etaBin0 = m_gtScales->getMUScales().etaBins.size() + etaBin0;  //twos complement
	
        etBin0 = etIndex0;
        int ssize = m_gtScales->getMUScales().etBins.size();
        if (etBin0 >= ssize) {
          etBin0 = ssize - 1;
          LogTrace("L1TGlobal") << "muon0 hw et" << etBin0 << " out of scale range. Setting to maximum.";
        }
	
        // Determine Floating Pt numbers for floating point calculation
        std::pair<double, double> binEdges = m_gtScales->getMUScales().phiBins.at(phiIndex0);
        phi0Phy = 0.5 * (binEdges.second + binEdges.first);
        binEdges = m_gtScales->getMUScales().etaBins.at(etaBin0);
        eta0Phy = 0.5 * (binEdges.second + binEdges.first);
        binEdges = m_gtScales->getMUScales().etBins.at(etBin0);
        et0Phy = 0.5 * (binEdges.second + binEdges.first);
	
        LogDebug("L1TGlobal") << "Found all quantities for the muon 0 in pair 0-1" << std::endl;
      }
    
    else {
        // Interested only in three-muon correlations
        LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 0" << std::endl;
        return false;
      }

    // SECOND OBJECT: Now loop over the second leg to get its information
    for (std::vector<SingleCombInCond>::const_iterator it1Comb = cond1Comb.begin(); it1Comb != cond1Comb.end();
         it1Comb++) {
      LogDebug("L1TGlobal") << "Looking at the second subcondition" << std::endl;
      int obj1Index = -1;

      if (!(*it1Comb).empty()) {
        obj1Index = (*it1Comb)[0];
      } else {
        LogTrace("L1TGlobal") << "\n  SingleCombInCond (*it1Comb).size() " << ((*it1Comb).size()) << std::endl;
        return false;
      }

      // Check to avoid the two legs either being the same object given that the type is the same (muon)
      if (cndObjTypeVec[0] == cndObjTypeVec[1] && obj0Index == obj1Index && cond0bx == cond1bx) {
        LogDebug("L1TGlobal") << "Corr Condition looking at same leg...skip" << std::endl;
        continue;
      }

      if (cond1Categ == CondMuon) {
          lutObj1 = "MU";
          candMuVec = m_uGtB->getCandL1Mu();
          phiIndex1 = (candMuVec->at(cond1bx, obj1Index))->hwPhiAtVtx();  //(*candMuVec)[obj0Index]->phiIndex();
          etaIndex1 = (candMuVec->at(cond1bx, obj1Index))->hwEtaAtVtx();
          etIndex1 = (candMuVec->at(cond1bx, obj1Index))->hwPt();
          etaBin1 = etaIndex1;
          if (etaBin1 < 0)
            etaBin1 = m_gtScales->getMUScales().etaBins.size() + etaBin1;
	  
          etBin1 = etIndex1;
          int ssize = m_gtScales->getMUScales().etBins.size();
          if (etBin1 >= ssize) {
            LogTrace("L1TGlobal") << "muon2 hw et" << etBin1 << " out of scale range.  Setting to maximum.";
            etBin1 = ssize - 1;
	  }
	  
          // Determine Floating Pt numbers for floating point calculation
          std::pair<double, double> binEdges = m_gtScales->getMUScales().phiBins.at(phiIndex1);
          phi1Phy = 0.5 * (binEdges.second + binEdges.first);
          binEdges = m_gtScales->getMUScales().etaBins.at(etaBin1);
          eta1Phy = 0.5 * (binEdges.second + binEdges.first);
          binEdges = m_gtScales->getMUScales().etBins.at(etBin1);
          et1Phy = 0.5 * (binEdges.second + binEdges.first); 

	  LogDebug("L1TGlobal") << "Found all quantities for the muon 1 in pair 0-1" << std::endl;
      } 
      else {
	  // Interested only in three-muon correlations
	  LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 1" << std::endl;
          return false;
        }
      
      //if (m_verbosity) {
      //  LogDebug("L1TGlobal") << "    First correlation pair [" << l1TGtObjectEnumToString(cndObjTypeVec[0]) << ", "
      std::cout << "\n ### EF First correlation pair [" << l1TGtObjectEnumToString(cndObjTypeVec[0]) << ", "
		<< l1TGtObjectEnumToString(cndObjTypeVec[1]) << "] with collection indices [" << obj0Index
		<< ", " << obj1Index << "] "
		<< " has: \n"
		<< "     Et  value   = [" << etIndex0 << ", " << etIndex1 << "]\n"
		<< "     phi indices = [" << phiIndex0 << ", " << phiIndex1 << "]\n"
		<< "     eta indices = [" << etaIndex0 << ", " << etaIndex1 << "]\n"
		<< std::endl;
      //}

      // Now perform the desired correlation on these two objects. Assume true until we find a contradition                                  
      
      // clear the indices in the combination                                                                                   
      objectsInComb.clear();
      objectsInComb.push_back(obj0Index);
      objectsInComb.push_back(obj1Index);
      
      // if we get here all checks were successful for this combination                                                     
      // set the general result for evaluateCondition to "true"                                                                                

      // Delta eta and phi calculations.                                                                                                         
      double deltaPhiPhy = fabs(phi1Phy - phi0Phy);
      if (deltaPhiPhy > M_PI)
	deltaPhiPhy = 2. * M_PI - deltaPhiPhy;
      double deltaEtaPhy = fabs(eta1Phy - eta0Phy);
      
      // Determine the integer based delta eta and delta phi                                                          
      int deltaPhiFW = abs(phiIndex0 - phiIndex1);
      if (deltaPhiFW >= phiBound)
	deltaPhiFW = 2 * phiBound - deltaPhiFW;
      std::string lutName = lutObj0;
      lutName += "-";
      lutName += lutObj1;
      long long deltaPhiLUT = m_gtScales->getLUT_DeltaPhi(lutName, deltaPhiFW);
      unsigned int precDeltaPhiLUT = m_gtScales->getPrec_DeltaPhi(lutName);
      
      int deltaEtaFW = abs(etaIndex0 - etaIndex1);
      long long deltaEtaLUT = 0;
      unsigned int precDeltaEtaLUT = 0;
      deltaEtaLUT = m_gtScales->getLUT_DeltaEta(lutName, deltaEtaFW);
      precDeltaEtaLUT = m_gtScales->getPrec_DeltaEta(lutName);
      
      //LogDebug("L1TGlobal") << "Obj0 phiFW = " << phiIndex0 << " Obj1 phiFW = " << phiIndex1 << "\n"
      std::cout << "Obj0 phiFW = " << phiIndex0 << " Obj1 phiFW = " << phiIndex1 << "\n"
		<< "    DeltaPhiFW = " << deltaPhiFW << "\n"
		<< "    LUT Name = " << lutName << " Prec = " << precDeltaPhiLUT
		<< "    DeltaPhiLUT = " << deltaPhiLUT << "\n"
		<< "Obj0 etaFW = " << etaIndex0 << " Obj1 etaFW = " << etaIndex1 << "\n"
		<< "    DeltaEtaFW = " << deltaEtaFW << "\n"
		<< "    LUT Name = " << lutName << " Prec = " << precDeltaEtaLUT
		<< "    DeltaEtaLUT = " << deltaEtaLUT << std::endl;
      
      if (corrPar.corrCutType & 0x8 || corrPar.corrCutType & 0x10) {
	//invariant mass calculation based on                                                                                                         
	// M = sqrt(2*p1*p2(cosh(eta1-eta2) - cos(phi1 - phi2)))                                                                                   
	// but we calculate (1/2)M^2                                                                                                                         
	//                                                                                                                                                            
	double cosDeltaPhiPhy = cos(deltaPhiPhy);
	double coshDeltaEtaPhy = cosh(deltaEtaPhy);
	if (corrPar.corrCutType & 0x10) coshDeltaEtaPhy = 1.;
	double massSqPhy = et0Phy * et1Phy * (coshDeltaEtaPhy - cosDeltaPhiPhy);
	
	long long cosDeltaPhiLUT = m_gtScales->getLUT_DeltaPhi_Cos(lutName, deltaPhiFW);
	unsigned int precCosLUT = m_gtScales->getPrec_DeltaPhi_Cos(lutName);
	
	long long coshDeltaEtaLUT;
	if (corrPar.corrCutType & 0x10) {
	  coshDeltaEtaLUT = 1 * pow(10, precCosLUT);
	} else {
	  coshDeltaEtaLUT = m_gtScales->getLUT_DeltaEta_Cosh(lutName, deltaEtaFW);
	  unsigned int precCoshLUT = m_gtScales->getPrec_DeltaEta_Cosh(lutName);
	  if (precCoshLUT - precCosLUT != 0)
	    LogDebug("L1TGlobal") << "Warning: Cos and Cosh LUTs on different Precision" << std::endl;
	}
	
	std::string lutName = lutObj0;
	lutName += "-ET";
	long long ptObj0 = m_gtScales->getLUT_Pt("Mass_" + lutName, etIndex0);
	unsigned int precPtLUTObj0 = m_gtScales->getPrec_Pt("Mass_" + lutName);
	
	lutName = lutObj1;
	lutName += "-ET";
	long long ptObj1 = m_gtScales->getLUT_Pt("Mass_" + lutName, etIndex1);
	unsigned int precPtLUTObj1 = m_gtScales->getPrec_Pt("Mass_" + lutName);
	
	// Pt and angles are at different precision.                                                                                                 
	long long massSq = ptObj0 * ptObj1 * (coshDeltaEtaLUT - cosDeltaPhiLUT);
	
	//Note: There is an assumption here that Cos and Cosh have the same precision                                                                         
	preShift = precPtLUTObj0 + precPtLUTObj1 + precCosLUT - corrPar.precMassCut;
	
	std::cout << "    Testing Invariant Mass for the first pair 0-1 (" << lutObj0 << "," << lutObj1 << ") ["
	  //LogDebug("L1TGlobal") << "    Testing Invariant Mass for the first pair 0-1 (" << lutObj0 << "," << lutObj1 << ") ["
			      << (long long)(corrPar.minMassCutValue * pow(10, preShift)) << ","
			      << (long long)(corrPar.maxMassCutValue * pow(10, preShift))
			      << "] with precision = " << corrPar.precMassCut << "\n"
			      << "    deltaPhiLUT  = " << deltaPhiLUT << "  cosLUT  = " << cosDeltaPhiLUT << "\n"
			      << "    deltaEtaLUT  = " << deltaEtaLUT << "  coshLUT = " << coshDeltaEtaLUT << "\n"
			      << "    etIndex0     = " << etIndex0 << "     pt0LUT  = " << ptObj0
			      << "    PhyEt0       = " << et0Phy << "\n"
			      << "    etIndex1     = " << etIndex1 << "     pt1LUT  = " << ptObj1
			      << "    PhyEt1       = " << et1Phy << "\n"
			      << "    massSq/2     = " << massSq << "\n"
			      << "    Precision Shift = " << preShift << "\n"
			      << "    massSq   (shift)= " << (massSq / pow(10, preShift + corrPar.precMassCut)) << "\n"
			      << "    deltaPhiPhy  = " << deltaPhiPhy << "  cos() = " << cosDeltaPhiPhy << "\n"
			      << "    deltaEtaPhy  = " << deltaEtaPhy << "  cosh()= " << coshDeltaEtaPhy << "\n"
			      << "    massSqPhy/2  = " << massSqPhy
			      << "    sqrt(|massSq|) = " << sqrt(fabs(2. * massSqPhy)) << std::endl;
	
	//EF:
	//std::cout << "Object pair has a pT of " << ptObj0 << " and " << ptObj1 << std::endl;                                                                                         
	muIndexes_01.push_back(etIndex0);                                                                                                    
	muIndexes_01.push_back(etIndex1);                                                                                                                                   
	dimuInvMass_01.push_back(massSq);                     
      }
      
      //EF
      /*std::cout << muIndexes_01.size() << std::endl;                                                                                                           
      std::cout << dimuInvMass_01.size() << std::endl;
      std::sort(dimuInvMass_01.begin(), dimuInvMass_01.end());
      if (muIndexes_01.size() != 0)                                                                                                                                    
	{                                                                                                                                                                  
	  for (long unsigned int i=0; i < muIndexes_01.size(); i++)                                                                                                       
	    { std::cout << "Position " << i << "  with index " << muIndexes_01.at(i) <<  std::endl; }                                       
	  for (long unsigned int m=0; m < dimuInvMass_01.size(); m++)                                                                                                         
	    { std::cout << "Position " << m << ": Inv mass calculated to be " << dimuInvMass_01.at(m) <<  " = " << dimuInvMass_01_GeV.at(m) << " GeV" << std::endl; }           
	}
      */                                                                                                                                               
    }  //end loop over second leg                                                                                                                                    
  }  //end loop over first leg                                                                                                                                                                             


  // *** SECOND PAIR: 0-2
  for (std::vector<SingleCombInCond>::const_iterator it0Comb = cond0Comb.begin(); it0Comb != cond0Comb.end();
       it0Comb++) {
    LogDebug("L1TGlobal") << "Looking at the first subcondition of the pair 0-2" << std::endl;
    // Type1s: there is 1 object only, no need for a loop, index 0 should be OK in (*it0Comb)[0]
    // ... but add protection to not crash
    int obj0Index = -1;

    if (!(*it0Comb).empty()) {
      obj0Index = (*it0Comb)[0];
    } else {
      LogTrace("L1TGlobal") << "\n  SingleCombInCond (*it0Comb).size() " << ((*it0Comb).size()) << std::endl;
      return false;
    }

    // FIRST OBJECT: Collect the information on the first leg of the correlation
    if (cond0Categ == CondMuon) {
        lutObj0 = "MU";
        candMuVec = m_uGtB->getCandL1Mu();
        phiIndex0 = (candMuVec->at(cond0bx, obj0Index))->hwPhiAtVtx();  //(*candMuVec)[obj0Index]->phiIndex();
        etaIndex0 = (candMuVec->at(cond0bx, obj0Index))->hwEtaAtVtx();
        etIndex0 = (candMuVec->at(cond0bx, obj0Index))->hwPt();
        int etaBin0 = etaIndex0;
        if (etaBin0 < 0) 
	  etaBin0 = m_gtScales->getMUScales().etaBins.size() + etaBin0;  //twos complement
	
        etBin0 = etIndex0;
        int ssize = m_gtScales->getMUScales().etBins.size();
        if (etBin0 >= ssize) {
          etBin0 = ssize - 1;
          LogTrace("L1TGlobal") << "muon0 hw et" << etBin0 << " out of scale range.  Setting to maximum.";
        }
	
        // Determine Floating Pt numbers for floating point calculation
        std::pair<double, double> binEdges = m_gtScales->getMUScales().phiBins.at(phiIndex0);
        phi0Phy = 0.5 * (binEdges.second + binEdges.first);
        binEdges = m_gtScales->getMUScales().etaBins.at(etaBin0);
        eta0Phy = 0.5 * (binEdges.second + binEdges.first);
        binEdges = m_gtScales->getMUScales().etBins.at(etBin0);
        et0Phy = 0.5 * (binEdges.second + binEdges.first);
	
        LogDebug("L1TGlobal") << "Found all quantities for the muon 0 in pair 0-2" << std::endl;
      }
    
    else {
        // Interested only in three-muon correlations
        LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 0" << std::endl;
        return false;
      }

    // SECOND OBJECT: Now loop over the second leg to get its information
    for (std::vector<SingleCombInCond>::const_iterator it2Comb = cond2Comb.begin(); it2Comb != cond2Comb.end();
         it2Comb++) {
      LogDebug("L1TGlobal") << "Looking at the third subcondition" << std::endl;
      int obj2Index = -1;

      if (!(*it2Comb).empty()) {
        obj2Index = (*it2Comb)[0];
      } else {
        LogTrace("L1TGlobal") << "\n  SingleCombInCond (*it2Comb).size() " << ((*it2Comb).size()) << std::endl;
        return false;
      }

      // Check to avoid the two legs either being the same object given that the type is the same (muon)
      if (cndObjTypeVec[0] == cndObjTypeVec[2] && obj0Index == obj2Index && cond0bx == cond2bx) {
        LogDebug("L1TGlobal") << "Corr Condition looking at same leg...skip" << std::endl;
        continue;
      }

      if (cond2Categ == CondMuon) {
          lutObj2 = "MU";
          candMuVec = m_uGtB->getCandL1Mu();
          phiIndex2 = (candMuVec->at(cond2bx, obj2Index))->hwPhiAtVtx();  //(*candMuVec)[obj0Index]->phiIndex();
          etaIndex2 = (candMuVec->at(cond2bx, obj2Index))->hwEtaAtVtx();
          etIndex2 = (candMuVec->at(cond2bx, obj2Index))->hwPt();
          etaBin2 = etaIndex2;
          if (etaBin2 < 0)
            etaBin2 = m_gtScales->getMUScales().etaBins.size() + etaBin2;
	  
          etBin2 = etIndex2;
          int ssize = m_gtScales->getMUScales().etBins.size();
          if (etBin2 >= ssize) {
            LogTrace("L1TGlobal") << "muon2 hw et" << etBin2 << " out of scale range.  Setting to maximum.";
            etBin2 = ssize - 2;
	  }
	  
          // Determine Floating Pt numbers for floating point calculation
          std::pair<double, double> binEdges = m_gtScales->getMUScales().phiBins.at(phiIndex2);
          phi2Phy = 0.5 * (binEdges.second + binEdges.first);
          binEdges = m_gtScales->getMUScales().etaBins.at(etaBin2);
          eta2Phy = 0.5 * (binEdges.second + binEdges.first);
          binEdges = m_gtScales->getMUScales().etBins.at(etBin2);
          et2Phy = 0.5 * (binEdges.second + binEdges.first); 

	  LogDebug("L1TGlobal") << "Found all quantities for the muon 2 in pair 0-2" << std::endl;
      } 
      else {
	  // Interested only in three-muon correlations
	  LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 2" << std::endl;
          return false;
        }
      
      //if (m_verbosity) {
      //  LogDebug("L1TGlobal") << "    Correlation pair [" << l1TGtObjectEnumToString(cndObjTypeVec[0]) << ", "
      std::cout << "\n ### EF Second correlation pair [" << l1TGtObjectEnumToString(cndObjTypeVec[0]) << ", "
		<< l1TGtObjectEnumToString(cndObjTypeVec[2]) << "] with collection indices [" << obj0Index
		<< ", " << obj2Index << "] "
		<< " has: \n"
		<< "     Et  value   = [" << etIndex0 << ", " << etIndex2 << "]\n"
		<< "     phi indices = [" << phiIndex0 << ", " << phiIndex2 << "]\n"
		<< "     eta indices = [" << etaIndex0 << ", " << etaIndex2 << "]\n"
		<< std::endl;
      //}

      // Now perform the desired correlation on these two objects. Assume true until we find a contradition                                  
      
      // clear the indices in the combination                                                                                   
      objectsInComb.clear();
      objectsInComb.push_back(obj0Index);
      objectsInComb.push_back(obj2Index);
      
      // if we get here all checks were successful for this combination                                                     
      // set the general result for evaluateCondition to "true"                                                                                

      // Delta eta and phi calculations.                                                                                                         
      double deltaPhiPhy = fabs(phi2Phy - phi0Phy);
      if (deltaPhiPhy > M_PI)
	deltaPhiPhy = 2. * M_PI - deltaPhiPhy;
      double deltaEtaPhy = fabs(eta2Phy - eta0Phy);
      
      // Determine the integer based delta eta and delta phi                                                          
      int deltaPhiFW = abs(phiIndex0 - phiIndex2);
      if (deltaPhiFW >= phiBound)
	deltaPhiFW = 2 * phiBound - deltaPhiFW;
      std::string lutName = lutObj0;
      lutName += "-";
      lutName += lutObj2;
      long long deltaPhiLUT = m_gtScales->getLUT_DeltaPhi(lutName, deltaPhiFW);
      unsigned int precDeltaPhiLUT = m_gtScales->getPrec_DeltaPhi(lutName);
      
      int deltaEtaFW = abs(etaIndex0 - etaIndex2);
      long long deltaEtaLUT = 0;
      unsigned int precDeltaEtaLUT = 0;
      deltaEtaLUT = m_gtScales->getLUT_DeltaEta(lutName, deltaEtaFW);
      precDeltaEtaLUT = m_gtScales->getPrec_DeltaEta(lutName);
      
      //LogDebug("L1TGlobal") << "Obj0 phiFW = " << phiIndex0 << " Obj2 phiFW = " << phiIndex2 << "\n"
      std::cout << "-----------------------------------" 
		<< "Obj0 phiFW = " << phiIndex0 << " Obj2 phiFW = " << phiIndex2 << "\n"
		<< "    DeltaPhiFW = " << deltaPhiFW << "\n"
		<< "    LUT Name = " << lutName << " Prec = " << precDeltaPhiLUT
		<< "    DeltaPhiLUT = " << deltaPhiLUT << "\n"
		<< "-----------------------------------" 
		<< "Obj0 etaFW = " << etaIndex0 << " Obj2 etaFW = " << etaIndex2 << "\n"
		<< "    DeltaEtaFW = " << deltaEtaFW << "\n"
		<< "    LUT Name = " << lutName << " Prec = " << precDeltaEtaLUT
		<< "    DeltaEtaLUT = " << deltaEtaLUT 
		<< "-----------------------------------" 
		<< std::endl;
      
      if (corrPar.corrCutType & 0x8 || corrPar.corrCutType & 0x10) {
	//invariant mass calculation based on                                                                                                         
	// M = sqrt(2*p1*p2(cosh(eta1-eta2) - cos(phi1 - phi2)))                                                                                   
	// but we calculate (1/2)M^2                                                                                                                         
	//                                                                                                                                                            
	double cosDeltaPhiPhy = cos(deltaPhiPhy);
	double coshDeltaEtaPhy = cosh(deltaEtaPhy);
	if (corrPar.corrCutType & 0x10)
	  coshDeltaEtaPhy = 1.;
	double massSqPhy = et0Phy * et2Phy * (coshDeltaEtaPhy - cosDeltaPhiPhy);
	
	long long cosDeltaPhiLUT = m_gtScales->getLUT_DeltaPhi_Cos(lutName, deltaPhiFW);
	unsigned int precCosLUT = m_gtScales->getPrec_DeltaPhi_Cos(lutName);
	
	long long coshDeltaEtaLUT;
	if (corrPar.corrCutType & 0x10) {
	  coshDeltaEtaLUT = 1 * pow(10, precCosLUT);
	} else {
	  coshDeltaEtaLUT = m_gtScales->getLUT_DeltaEta_Cosh(lutName, deltaEtaFW);
	  unsigned int precCoshLUT = m_gtScales->getPrec_DeltaEta_Cosh(lutName);
	  if (precCoshLUT - precCosLUT != 0)
	    LogDebug("L1TGlobal") << "Warning: Cos and Cosh LUTs on different Precision" << std::endl;
	}
	// if (corrPar.corrCutType & 0x10) coshDeltaEtaLUT=1*pow(10,precCosLUT);                                                                   
	
	std::string lutName = lutObj0;
	lutName += "-ET";
	long long ptObj0 = m_gtScales->getLUT_Pt("Mass_" + lutName, etIndex0);
	unsigned int precPtLUTObj0 = m_gtScales->getPrec_Pt("Mass_" + lutName);
	
	lutName = lutObj2;
	lutName += "-ET";
	long long ptObj2 = m_gtScales->getLUT_Pt("Mass_" + lutName, etIndex2);
	unsigned int precPtLUTObj2 = m_gtScales->getPrec_Pt("Mass_" + lutName);
	
	// Pt and Angles are at different precision.                                                                                                 
	long long massSq = ptObj0 * ptObj2 * (coshDeltaEtaLUT - cosDeltaPhiLUT);
	
	//Note: There is an assumption here that Cos and Cosh have the same precision                                                                         
	preShift = precPtLUTObj0 + precPtLUTObj2 + precCosLUT - corrPar.precMassCut;
	
	std::cout << "    Testing Invariant Mass for the second pair 0-2 (" << lutObj0 << "," << lutObj2 << ") ["
	  //LogDebug("L1TGlobal") << "    Testing Invariant Mass for the second pair 0-2 (" << lutObj0 << "," << lutObj2 << ") ["
			      << (long long)(corrPar.minMassCutValue * pow(10, preShift)) << ","
			      << (long long)(corrPar.maxMassCutValue * pow(10, preShift))
			      << "] with precision = " << corrPar.precMassCut << "\n"
			      << "    deltaPhiLUT  = " << deltaPhiLUT << "  cosLUT  = " << cosDeltaPhiLUT << "\n"
			      << "    deltaEtaLUT  = " << deltaEtaLUT << "  coshLUT = " << coshDeltaEtaLUT << "\n"
			      << "    etIndex0     = " << etIndex0 << "    pt0LUT      = " << ptObj0
			      << " PhyEt0 = " << et0Phy << "\n"
			      << "    etIndex2     = " << etIndex2 << "    pt2LUT      = " << ptObj2
			      << " PhyEt2 = " << et2Phy << "\n"
			      << "    massSq/2     = " << massSq << "\n"
			      << "    Precision Shift = " << preShift << "\n"
			      << "    massSq   (shift)= " << (massSq / pow(10, preShift + corrPar.precMassCut)) << "\n"
			      << "    deltaPhiPhy  = " << deltaPhiPhy << "  cos() = " << cosDeltaPhiPhy << "\n"
			      << "    deltaEtaPhy  = " << deltaEtaPhy << "  cosh()= " << coshDeltaEtaPhy << "\n"
			      << "    massSqPhy/2  = " << massSqPhy
			      << "  sqrt(|massSq|) = " << sqrt(fabs(2. * massSqPhy)) << std::endl;
	
	//EF                                                                                        
	//std::cout << "Object pair has a pT of " << ptObj0 << " and " << ptObj2 << std::endl;                                                                                         
	muIndexes_02.push_back(etIndex0);                                                                                                    
	muIndexes_02.push_back(etIndex2);                                                                                                                                   
	dimuInvMass_02.push_back(massSq);                                                                                                       
      }
      //EF                                                                                                                                       
      /*std::cout << muIndexes_02.size() << std::endl;                                                                                                           
      std::cout << dimuInvMass_02.size() << std::endl;
      if (muIndexes_02.size() != 0)                                                                                                                                    
	{                                                                                                                                                                  
	  for (long unsigned int i=0; i < muIndexes_02.size(); i++)                                                                                                       
	    { std::cout << "Position " << i << "  with index " << muIndexes_02.at(i) <<  std::endl; }                                                              
	  for (long unsigned int m=0; m < dimuInvMass_02.size(); m++)                                                                                                         
	    { std::cout << "Position " << m << ": Inv mass calculated to be " << dimuInvMass_02.at(m) <<  " = " << dimuInvMass_02_GeV.at(m) << " GeV" << std::endl; }                   
	}                                                                                                                                              
      */
    }  //end loop over second leg                                                                                                                                    
  }  //end loop over first leg                                                                                                                                                                             
  

  // *** THIRD PAIR: 1-2
  for (std::vector<SingleCombInCond>::const_iterator it1Comb = cond1Comb.begin(); it1Comb != cond1Comb.end();
       it1Comb++) {
    LogDebug("L1TGlobal") << "Looking at the first subcondition of the pair 1-2" << std::endl;
    // Type1s: there is 1 object only, no need for a loop, index 1 should be OK in (*it1Comb)[0]
    // ... but add protection to not crash
    int obj1Index = -1;

    if (!(*it1Comb).empty()) {
      obj1Index = (*it1Comb)[0];
    } else {
      LogTrace("L1TGlobal") << "\n  SingleCombInCond (*it1Comb).size() " << ((*it1Comb).size()) << std::endl;
      return false;
    }

    // FIRST OBJECT: Collect the information on the first leg of the correlation
    if (cond1Categ == CondMuon) {
        lutObj1 = "MU";
        candMuVec = m_uGtB->getCandL1Mu();
        phiIndex1 = (candMuVec->at(cond1bx, obj1Index))->hwPhiAtVtx();  //(*candMuVec)[obj1Index]->phiIndex();
        etaIndex1 = (candMuVec->at(cond1bx, obj1Index))->hwEtaAtVtx();
        etIndex1 = (candMuVec->at(cond1bx, obj1Index))->hwPt();
        int etaBin1 = etaIndex1;
        if (etaBin1 < 0) 
	  etaBin1 = m_gtScales->getMUScales().etaBins.size() + etaBin1;  //twos complement
        // LogDebug("L1TGlobal") << "Muon phi" << phiIndex1 << " eta " << etaIndex1 << " etaBin1 = " << etaBin1  << " et " << etIndex1 << std::endl;
	
        etBin1 = etIndex1;
        int ssize = m_gtScales->getMUScales().etBins.size();
        if (etBin1 >= ssize) {
          etBin1 = ssize - 1;
          LogTrace("L1TGlobal") << "muon1 hw et" << etBin1 << " out of scale range.  Setting to maximum.";
        }
	
        // Determine Floating Pt numbers for floating point calculation
        std::pair<double, double> binEdges = m_gtScales->getMUScales().phiBins.at(phiIndex1);
        phi1Phy = 0.5 * (binEdges.second + binEdges.first);
        binEdges = m_gtScales->getMUScales().etaBins.at(etaBin1);
        eta1Phy = 0.5 * (binEdges.second + binEdges.first);
        binEdges = m_gtScales->getMUScales().etBins.at(etBin1);
        et1Phy = 0.5 * (binEdges.second + binEdges.first);
	
        LogDebug("L1TGlobal") << "Found all quantities for the muon 1 in pair 1-2" << std::endl;
      } //break;
    
    else {
        // Interested only in three-muon correlations
        LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 0" << std::endl;
        return false;
      }

    // SECOND OBJECT: Now loop over the second leg to get its information
    for (std::vector<SingleCombInCond>::const_iterator it2Comb = cond2Comb.begin(); it2Comb != cond2Comb.end();
         it2Comb++) {
      LogDebug("L1TGlobal") << "Looking at the third subcondition" << std::endl;
      int obj2Index = -1;

      if (!(*it2Comb).empty()) {
        obj2Index = (*it2Comb)[0];
      } else {
        LogTrace("L1TGlobal") << "\n  SingleCombInCond (*it2Comb).size() " << ((*it2Comb).size()) << std::endl;
        return false;
      }

      // Check to avoid the two legs either being the same object given that the type is the same (muon)
      if (cndObjTypeVec[1] == cndObjTypeVec[2] && obj1Index == obj2Index && cond1bx == cond2bx) {
        LogDebug("L1TGlobal") << "Corr Condition looking at same leg...skip" << std::endl;
        continue;
      }

      if (cond2Categ == CondMuon) {
          lutObj2 = "MU";
          candMuVec = m_uGtB->getCandL1Mu();
          phiIndex2 = (candMuVec->at(cond2bx, obj2Index))->hwPhiAtVtx();  //(*candMuVec)[obj1Index]->phiIndex();
          etaIndex2 = (candMuVec->at(cond2bx, obj2Index))->hwEtaAtVtx();
          etIndex2 = (candMuVec->at(cond2bx, obj2Index))->hwPt();
          etaBin2 = etaIndex2;
          if (etaBin2 < 0)
            etaBin2 = m_gtScales->getMUScales().etaBins.size() + etaBin2;
          // LogDebug("L1TGlobal") << "Muon phi" << phiIndex2 << " eta " << etaIndex2 << " etaBin2 = " << etaBin2  << " et " << etIndex2 << std::endl;
	  
          etBin2 = etIndex2;
          int ssize = m_gtScales->getMUScales().etBins.size();
          if (etBin2 >= ssize) {
            LogTrace("L1TGlobal") << "muon2 hw et" << etBin2 << " out of scale range.  Setting to maximum.";
            etBin2 = ssize - 2;
	  }
	  
          // Determine Floating Pt numbers for floating point calculation
          std::pair<double, double> binEdges = m_gtScales->getMUScales().phiBins.at(phiIndex2);
          phi2Phy = 0.5 * (binEdges.second + binEdges.first);
          binEdges = m_gtScales->getMUScales().etaBins.at(etaBin2);
          eta2Phy = 0.5 * (binEdges.second + binEdges.first);
          binEdges = m_gtScales->getMUScales().etBins.at(etBin2);
          et2Phy = 0.5 * (binEdges.second + binEdges.first); 

	  LogDebug("L1TGlobal") << "Found all quantities for the muon 2 in pair 1-2" << std::endl;
      } 
      else {
	  // Interested only in three-muon correlations
	  LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 2" << std::endl;
          return false;
        }
      
      //if (m_verbosity) {
      //  LogDebug("L1TGlobal") << "    Correlation pair [" << l1TGtObjectEnumToString(cndObjTypeVec[1]) << ", "
      std::cout << "\n ### EF Third correlation pair [" << l1TGtObjectEnumToString(cndObjTypeVec[1]) << ", "
		<< l1TGtObjectEnumToString(cndObjTypeVec[2]) << "] with collection indices [" << obj1Index
		<< ", " << obj2Index << "] "
		<< " has: \n"
		<< "     Et  value   = [" << etIndex1 << ", " << etIndex2 << "]\n"
		<< "     phi indices = [" << phiIndex1 << ", " << phiIndex2 << "]\n"
		<< "     eta indices = [" << etaIndex1 << ", " << etaIndex2 << "]\n"
		<< std::endl;
      //}

      // Now perform the desired correlation on these two objects. Assume true until we find a contradition                                  
      
      // clear the indices in the combination                                                                                   
      objectsInComb.clear();
      objectsInComb.push_back(obj1Index);
      objectsInComb.push_back(obj2Index);
      
      // if we get here all checks were successful for this combination                                                     
      // set the general result for evaluateCondition to "true"                                                                                

      // Delta eta and phi calculations.                                                                                                         
      double deltaPhiPhy = fabs(phi2Phy - phi1Phy);
      if (deltaPhiPhy > M_PI)
	deltaPhiPhy = 2. * M_PI - deltaPhiPhy;
      double deltaEtaPhy = fabs(eta2Phy - eta1Phy);
      
      // Determine the integer based delta eta and delta phi                                                          
      int deltaPhiFW = abs(phiIndex1 - phiIndex2);
      if (deltaPhiFW >= phiBound)
	deltaPhiFW = 2 * phiBound - deltaPhiFW;
      std::string lutName = lutObj1;
      lutName += "-";
      lutName += lutObj2;
      long long deltaPhiLUT = m_gtScales->getLUT_DeltaPhi(lutName, deltaPhiFW);
      unsigned int precDeltaPhiLUT = m_gtScales->getPrec_DeltaPhi(lutName);
      
      int deltaEtaFW = abs(etaIndex1 - etaIndex2);
      long long deltaEtaLUT = 0;
      unsigned int precDeltaEtaLUT = 0;
      deltaEtaLUT = m_gtScales->getLUT_DeltaEta(lutName, deltaEtaFW);
      precDeltaEtaLUT = m_gtScales->getPrec_DeltaEta(lutName);
      
      //LogDebug("L1TGlobal") << "Obj1 phiFW = " << phiIndex1 << " Obj2 phiFW = " << phiIndex2 << "\n"
      std::cout << "Obj1 phiFW = " << phiIndex1 << " Obj2 phiFW = " << phiIndex2 << "\n"
		<< "    DeltaPhiFW = " << deltaPhiFW << "\n"
		<< "    LUT Name = " << lutName << " Prec = " << precDeltaPhiLUT
		<< "    DeltaPhiLUT = " << deltaPhiLUT << "\n"
		<< "Obj1 etaFW = " << etaIndex1 << " Obj2 etaFW = " << etaIndex2 << "\n"
		<< "    DeltaEtaFW = " << deltaEtaFW << "\n"
		<< "    LUT Name = " << lutName << " Prec = " << precDeltaEtaLUT
		<< "    DeltaEtaLUT = " << deltaEtaLUT << std::endl;
      
      if (corrPar.corrCutType & 0x8 || corrPar.corrCutType & 0x10) {
	//invariant mass calculation based on                                                                                                         
	// M = sqrt(2*p1*p2(cosh(eta1-eta2) - cos(phi1 - phi2)))                                                                                   
	// but we calculate (1/2)M^2                                                                                                                         
	//                                                                                                                                                            
	double cosDeltaPhiPhy = cos(deltaPhiPhy);
	double coshDeltaEtaPhy = cosh(deltaEtaPhy);
	if (corrPar.corrCutType & 0x10)
	  coshDeltaEtaPhy = 1.;
	double massSqPhy = et1Phy * et2Phy * (coshDeltaEtaPhy - cosDeltaPhiPhy);
	
	long long cosDeltaPhiLUT = m_gtScales->getLUT_DeltaPhi_Cos(lutName, deltaPhiFW);
	unsigned int precCosLUT = m_gtScales->getPrec_DeltaPhi_Cos(lutName);
	
	long long coshDeltaEtaLUT;
	if (corrPar.corrCutType & 0x10) {
	  coshDeltaEtaLUT = 1 * pow(10, precCosLUT);
	} else {
	  coshDeltaEtaLUT = m_gtScales->getLUT_DeltaEta_Cosh(lutName, deltaEtaFW);
	  unsigned int precCoshLUT = m_gtScales->getPrec_DeltaEta_Cosh(lutName);
	  if (precCoshLUT - precCosLUT != 0)
	    LogDebug("L1TGlobal") << "Warning: Cos and Cosh LUTs on different Precision" << std::endl;
	}
	
	std::string lutName = lutObj1;
	lutName += "-ET";
	long long ptObj1 = m_gtScales->getLUT_Pt("Mass_" + lutName, etIndex1);
	unsigned int precPtLUTObj1 = m_gtScales->getPrec_Pt("Mass_" + lutName);
	
	lutName = lutObj2;
	lutName += "-ET";
	long long ptObj2 = m_gtScales->getLUT_Pt("Mass_" + lutName, etIndex2);
	unsigned int precPtLUTObj2 = m_gtScales->getPrec_Pt("Mass_" + lutName);
	
	// Pt and Angles are at different precision.                                                                                                 
	long long massSq = ptObj1 * ptObj2 * (coshDeltaEtaLUT - cosDeltaPhiLUT);
	
	//Note: There is an assumption here that Cos and Cosh have the same precision                                                                         
	preShift = precPtLUTObj1 + precPtLUTObj2 + precCosLUT - corrPar.precMassCut;
	
	std::cout << "    Testing Invariant Mass for the third pair 1-2 (" << lutObj1 << "," << lutObj2 << ") ["
	  //LogDebug("L1TGlobal") << "    Testing Invariant Mass for the third pair 1-2 (" << lutObj1 << "," << lutObj2 << ") ["
			      << (long long)(corrPar.minMassCutValue * pow(10, preShift)) << ","
			      << (long long)(corrPar.maxMassCutValue * pow(10, preShift))
			      << "] with precision = " << corrPar.precMassCut << "\n"
			      << "    deltaPhiLUT  = " << deltaPhiLUT << "  cosLUT  = " << cosDeltaPhiLUT << "\n"
			      << "    deltaEtaLUT  = " << deltaEtaLUT << "  coshLUT = " << coshDeltaEtaLUT << "\n"
			      << "    etIndex1     = " << etIndex1 << "    pt1LUT      = " << ptObj1
			      << " PhyEt1 = " << et1Phy << "\n"
			      << "    etIndex2     = " << etIndex2 << "    pt2LUT      = " << ptObj2
			      << " PhyEt2 = " << et2Phy << "\n"
			      << "    massSq/2     = " << massSq << "\n"
			      << "    Precision Shift = " << preShift << "\n"
			      << "    massSq   (shift)= " << (massSq / pow(10, preShift + corrPar.precMassCut)) << "\n"
			      << "    deltaPhiPhy  = " << deltaPhiPhy << "  cos() = " << cosDeltaPhiPhy << "\n"
			      << "    deltaEtaPhy  = " << deltaEtaPhy << "  cosh()= " << coshDeltaEtaPhy << "\n"
			      << "    massSqPhy/2  = " << massSqPhy
			      << "  sqrt(|massSq|) = " << sqrt(fabs(2. * massSqPhy)) << std::endl;

	//EF                                                                                        
	std::cout << "Object pair has a pT of " << ptObj1 << " and " << ptObj2 << std::endl;                                                                                         
	muIndexes_12.push_back(etIndex1);                                                                                                    
	muIndexes_12.push_back(etIndex2);                                                                                                                                   
	dimuInvMass_12.push_back(massSq);                                                                                                                     
      }
      //EF                                                                                                                                       
      std::cout << muIndexes_12.size() << std::endl;                                                                                                           
      std::cout << dimuInvMass_12.size() << std::endl;
      if (muIndexes_12.size() != 0)                                                                                                                                    
	{                                                                                                                                                                  
	  for (long unsigned int i=0; i < muIndexes_12.size(); i++)                                                                                                       
	    { std::cout << "Position " << i << "  with index " << muIndexes_12.at(i) <<  std::endl; }                                                    
	  for (long unsigned int m=0; m < dimuInvMass_12.size(); m++)                             
	    { std::cout << "Position " << m << ": Inv mass calculated to be " << dimuInvMass_12.at(m) <<  " = " << dimuInvMass_12_GeV.at(m) << " GeV" << std::endl; }                   
	}
                                                                                                   
    }  //end loop over second leg                                                                                                                                    
  }  //end loop over first leg                                                                                                                    
  
  // EF 
  long long massSq_3body = 0;
  bool reqResult = false;
  vector<long long> trimuInvMass;
  for (auto it_01 = dimuInvMass_01.cbegin(); it_01 != dimuInvMass_01.cend(); it_01++){
    for (auto it_02 = dimuInvMass_02.cbegin(); it_02 != dimuInvMass_02.cend(); it_02++){
      for (auto it_12 = dimuInvMass_12.cbegin(); it_12 != dimuInvMass_12.cend(); it_12++){
	if((*it_01 != *it_02) && (*it_01 != *it_12) && (*it_02 != *it_12))
	  { 
	    LogDebug("L1TGlobal") << "Iterators are: " << *it_01 << ", " << *it_02 << ", " << *it_12 << std::endl;
	    //std::cout << "Iterators are: " << *it_01 << ", " << *it_02 << ", " << *it_12 << std::endl;
	    massSq_3body = *it_01 + *it_02 + *it_12;
	    reqResult = true;
	    LogDebug("L1TGlobal") << "\n -----> EF: massSq_3body = " << massSq_3body << std::endl;
	    //std::cout << "\n -----> EF: massSq_3body = " << massSq_3body << std::endl;
	    
	    if (massSq_3body >= 0 && massSq_3body >= (long long)(corrPar.minMassCutValue * pow(10, preShift)) &&
		massSq_3body <= (long long)(corrPar.maxMassCutValue * pow(10, preShift))) {
	      LogDebug("L1TGlobal") << "    Passed Invariant Mass Cut ["
				    << (long long)(corrPar.minMassCutValue * pow(10, preShift)) << ","
				    << (long long)(corrPar.maxMassCutValue * pow(10, preShift)) << "]" << std::endl;
	      trimuInvMass.push_back(massSq_3body);
	    } 
	    else {
	      LogDebug("L1TGlobal") << "    Failed Invariant Mass Cut ["
				    << (long long)(corrPar.minMassCutValue * pow(10, preShift)) << ","
				    << (long long)(corrPar.maxMassCutValue * pow(10, preShift)) << "]" << std::endl;
	      reqResult = false;
	    }
	    
	  }
	else
	  {
	    LogDebug("L1TGlobal") << " Same objects considered, skip!" << std::endl;	
	  }
      }
    }
  }

  /*if (trimuInvMass.size() != 0)                                                                                                                                    
    {                                                                                                                                                                  
      for (long unsigned int i=0; i < trimuInvMass.size(); i++)                                                                                                       
	{ std::cout << "trimuInvMass position " << i << "  with value " << trimuInvMass.at(i) <<  std::endl; }                                       
    }                                                                                                                                               
  */
  if (reqResult) {
    condResult = true;
    (combinationsInCond()).push_back(objectsInComb);  
  }
  
  if (m_verbosity && condResult) {
    LogDebug("L1TGlobal") << " pass(es) the three-body correlation condition.\n" << std::endl;
  }
  return condResult;
}



/**
 * checkObjectParameter - Compare a single particle with a numbered condition.
 *
 * @param iCondition The number of the condition.
 * @param cand The candidate to compare.
 *
 * @return The result of the comparison (false if a condition does not exist).
 */

const bool l1t::CorrThreeBodyCondition::checkObjectParameter(const int iCondition, const l1t::L1Candidate& cand) const {
  return true;
}

void l1t::CorrThreeBodyCondition::print(std::ostream& myCout) const {
  myCout << "Dummy Print for CorrThreeBodyCondition" << std::endl;
  m_gtCorrelationThreeBodyTemplate->print(myCout);

  ConditionEvaluation::print(myCout);
}
