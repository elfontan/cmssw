/**
 * \class CorrThreeBodyCondition
 *
 * Description: evaluation of a three-body correlation condition (= three-muon invariant mass). 
 *
 */

// this class header
#include "L1Trigger/L1TGlobal/interface/CorrCondition.h" // FIXME EF => Kept for initial tests: probably TO BE REMOVED LATER
#include "L1Trigger/L1TGlobal/interface/CorrThreeBodyCondition.h"

// system include files
#include <iostream>
#include <iomanip>

#include <string>
#include <vector>
#include <algorithm>

// user include files
//   base classes
#include "L1Trigger/L1TGlobal/interface/CorrelationTemplate.h" // FIXME EF => Kept for initial tests: probably TO BE REMOVED LATER
#include "L1Trigger/L1TGlobal/interface/CorrelationThreeBodyTemplate.h"
#include "L1Trigger/L1TGlobal/interface/ConditionEvaluation.h"

#include "L1Trigger/L1TGlobal/interface/MuCondition.h"
#include "L1Trigger/L1TGlobal/interface/CaloCondition.h"
#include "L1Trigger/L1TGlobal/interface/EnergySumCondition.h"
#include "L1Trigger/L1TGlobal/interface/MuonTemplate.h"
#include "L1Trigger/L1TGlobal/interface/CaloTemplate.h"
#include "L1Trigger/L1TGlobal/interface/EnergySumTemplate.h"
#include "L1Trigger/L1TGlobal/interface/GlobalScales.h"

#include "DataFormats/L1Trigger/interface/L1Candidate.h"

#include "L1Trigger/L1TGlobal/interface/GlobalBoard.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

// constructors
//     default
l1t::CorrThreeBodyCondition::CorrThreeBodyCondition() : ConditionEvaluation() {}

//     from base template condition (from event setup usually)
l1t::CorrThreeBodyCondition::CorrThreeBodyCondition(const GlobalCondition* corrTemplate,
                                  const GlobalCondition* cond0Condition,
                                  const GlobalCondition* cond1Condition,
                                  const GlobalCondition* cond2Condition,
                                  const GlobalBoard* ptrGTB)
    : ConditionEvaluation(),
      m_gtCorrelationThreeBodyTemplate(static_cast<const CorrelationThreeBodyTemplate*>(corrTemplate)),
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

// try all object permutations and check spatial correlations, if required
const bool l1t::CorrThreeBodyCondition::evaluateCondition(const int bxEval) const {
  // std::cout << "m_isDebugEnabled = " << m_isDebugEnabled << std::endl;
  // std::cout << "m_verbosity = " << m_verbosity << std::endl;

  //std::ostringstream myCout;
  //m_gtCorrelationThreeBodyTemplate->print(myCout); 
  //LogDebug("L1TGlobal")
  //std::cout << "Correlation Condition Evaluation \n" << myCout.str() << std::endl; 

  bool condResult = false;
  bool reqObjResult = false;

  // number of objects in condition (it is 2, no need to retrieve from
  // condition template) and their type
  int nObjInCond = 3;
  std::vector<GlobalObject> cndObjTypeVec(nObjInCond);

  // evaluate first the two subconditions (Type1s)

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
  //switch (cond0Categ) {
  if (cond0Categ == CondMuon) {
      corrMuon = static_cast<const MuonTemplate*>(m_gtCond0);
      MuCondition muCondition(
			      corrMuon, m_uGtB, 0, 0);  //BLW these are counts that don't seem to be used...perhaps remove
      
      muCondition.evaluateConditionStoreResult(bxEval);
      reqObjResult = muCondition.condLastResult();
      
      cond0Comb = (muCondition.getCombinationsInCond());
      cond0bx = bxEval + (corrMuon->condRelativeBx());
      cndObjTypeVec[0] = (corrMuon->objectType())[0];
      
      if (m_verbosity) {
        std::ostringstream myCout;
        muCondition.print(myCout);
	
        LogDebug("L1TGlobal") << myCout.str() << std::endl;
      }
    } //break;
  
  else {
      // Interested only in three-muon correlations
      LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 0" << std::endl;
      return false;
    } //break;
  //} //Switch removed
  
  // return if first subcondition is false
  if (!reqObjResult) {
    LogDebug("L1TGlobal") << "\n  First sub-condition false, second sub-condition not evaluated and not printed."
                          << std::endl;
    return false;
  }

  // SECOND OBJECT
  reqObjResult = false;

  //switch (cond1Categ) {
  if (cond1Categ == CondMuon) {
      corrMuon = static_cast<const MuonTemplate*>(m_gtCond1);
      MuCondition muCondition(
			      corrMuon, m_uGtB, 0, 0);  //BLW these are counts that don't seem to be used...perhaps remove
      
      muCondition.evaluateConditionStoreResult(bxEval);
      reqObjResult = muCondition.condLastResult();
      
      cond1Comb = (muCondition.getCombinationsInCond());
      cond1bx = bxEval + (corrMuon->condRelativeBx());
      
      cndObjTypeVec[1] = (corrMuon->objectType())[0];
      
      if (m_verbosity) {
        std::ostringstream myCout;
        muCondition.print(myCout);
	
        LogDebug("L1TGlobal") << myCout.str() << std::endl;
      }
    } //break;
  
  else {
    // Interested only in three-muon correlations
    LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 1" << std::endl;
    return false;
  } //break;
  //} //Switch removed
  
  // return if second subcondition is false
  if (!reqObjResult) {
    LogDebug("L1TGlobal") << "\n  Second sub-condition false, third sub-condition not evaluated and not printed."
                          << std::endl;
    return false;
  }

  // THIRD OBJECT
  reqObjResult = false;

  if (cond2Categ == CondMuon) {
      corrMuon = static_cast<const MuonTemplate*>(m_gtCond2);
      MuCondition muCondition(
			      corrMuon, m_uGtB, 0, 0);  //BLW these are counts that don't seem to be used...perhaps remove
      
      muCondition.evaluateConditionStoreResult(bxEval);
      reqObjResult = muCondition.condLastResult();
      
      cond2Comb = (muCondition.getCombinationsInCond());
      cond2bx = bxEval + (corrMuon->condRelativeBx()); 
      
      cndObjTypeVec[2] = (corrMuon->objectType())[0];
      
      if (m_verbosity) {
        std::ostringstream myCout;
        muCondition.print(myCout);
        LogDebug("L1TGlobal") << myCout.str() << std::endl;
      }
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
    LogDebug("L1TGlobal") << "\n"
                          << "    All subconditions true for object requirements."
                          << "    Evaluate three-body correlation requirements.\n"
                          << std::endl;
  }

  // since we have two good legs get the correlation parameters
  CorrelationThreeBodyTemplate::CorrelationThreeBodyParameter corrPar = *(m_gtCorrelationThreeBodyTemplate->correlationParameter());

  // vector to store the indices of the calorimeter objects
  // from the combination evaluated in the condition
  SingleCombInCond objectsInComb;
  objectsInComb.reserve(nObjInCond);

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
  //int phiIndex2 = 0;
  //double phi2Phy = 0.;

  int etaIndex0 = 0;
  double eta0Phy = 0.;
  int etaBin0 = 0;
  int etaIndex1 = 0;
  double eta1Phy = 0.;
  int etaBin1 = 0;
  int etaIndex2 = 0;
  //double eta2Phy = 0.;
  int etaBin2 = 0;

  int etIndex0 = 0;
  int etBin0 = 0;
  double et0Phy = 0.;
  int etIndex1 = 0;
  int etBin1 = 0;
  double et1Phy = 0.;
  int etIndex2 = 0;
  int etBin2 = 0;
  //double et2Phy = 0.;

  // Determine the number of phi bins to get cutoff at pi
  int phiBound = 0;
  if (cond0Categ == CondMuon || cond1Categ == CondMuon || cond2Categ == CondMuon) {
    const GlobalScales::ScaleParameters& par = m_gtScales->getMUScales();
    //phiBound = par.phiBins.size()/2;
    phiBound = (int)((par.phiMax - par.phiMin) / par.phiStep) / 2;
  } else {
    //Assumes all calorimeter objects are on same phi scale
    const GlobalScales::ScaleParameters& par = m_gtScales->getEGScales();
    //phiBound = par.phiBins.size()/2;
    phiBound = (int)((par.phiMax - par.phiMin) / par.phiStep) / 2;
  }
  LogDebug("L1TGlobal") << "Phi Bound = " << phiBound << std::endl;

  // Keep track of objects for LUTS
  std::string lutObj0 = "NULL";
  std::string lutObj1 = "NULL";
  std::string lutObj2 = "NULL";

  LogTrace("L1TGlobal") << "  Subcondition 0: std::vector<SingleCombInCond> size: " << (cond0Comb.size()) << std::endl;
  LogTrace("L1TGlobal") << "  Subcondition 1: std::vector<SingleCombInCond> size: " << (cond1Comb.size()) << std::endl;
  LogTrace("L1TGlobal") << "  Subcondition 2: std::vector<SingleCombInCond> size: " << (cond2Comb.size()) << std::endl;

  // LOOP OVER ALL COMBINATIONS
  //
  // BLW: Optimization issue: potentially making the same comparison twice
  //                          if both legs are the same object type.
  for (std::vector<SingleCombInCond>::const_iterator it0Comb = cond0Comb.begin(); it0Comb != cond0Comb.end();
       it0Comb++) {
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
    // switch (cond0Categ) {
    if (cond0Categ == CondMuon) {
        lutObj0 = "MU";
        candMuVec = m_uGtB->getCandL1Mu();
        phiIndex0 = (candMuVec->at(cond0bx, obj0Index))->hwPhiAtVtx();  //(*candMuVec)[obj0Index]->phiIndex();
        etaIndex0 = (candMuVec->at(cond0bx, obj0Index))->hwEtaAtVtx();
        etIndex0 = (candMuVec->at(cond0bx, obj0Index))->hwPt();
        int etaBin0 = etaIndex0;
        if (etaBin0 < 0) 
	  etaBin0 = m_gtScales->getMUScales().etaBins.size() + etaBin0;  //twos complement
        // LogDebug("L1TGlobal") << "Muon phi" << phiIndex0 << " eta " << etaIndex0 << " etaBin0 = " << etaBin0  << " et " << etIndex0 << std::endl;
	
        etBin0 = etIndex0;
        int ssize = m_gtScales->getMUScales().etBins.size();
        if (etBin0 >= ssize) {
          etBin0 = ssize - 1;
          LogTrace("L1TGlobal") << "muon0 hw et" << etBin0 << " out of scale range.  Setting to maximum.";
        }
	
        // Determine Floating Pt numbers for floating point caluclation
        std::pair<double, double> binEdges = m_gtScales->getMUScales().phiBins.at(phiIndex0);
        phi0Phy = 0.5 * (binEdges.second + binEdges.first);
        binEdges = m_gtScales->getMUScales().etaBins.at(etaBin0);
        eta0Phy = 0.5 * (binEdges.second + binEdges.first);
        binEdges = m_gtScales->getMUScales().etBins.at(etBin0);
        et0Phy = 0.5 * (binEdges.second + binEdges.first);
	
        LogDebug("L1TGlobal") << "Found all quantities for the muon 0" << std::endl;
      } //break;
    
    else {
        // Interested only in three-muon correlations
        LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 0" << std::endl;
        return false;
      } //break;
    //}  //end check on first leg type => Switch removed

    // SECOND OBJECT: Now loop over the second leg to get its information
    for (std::vector<SingleCombInCond>::const_iterator it1Comb = cond1Comb.begin(); it1Comb != cond1Comb.end();
         it1Comb++) {
      LogDebug("L1TGlobal") << "Looking at second Condition" << std::endl;
      int obj1Index = -1;

      if (!(*it1Comb).empty()) {
        obj1Index = (*it1Comb)[0];
      } else {
        LogTrace("L1TGlobal") << "\n  SingleCombInCond (*it1Comb).size() " << ((*it1Comb).size()) << std::endl;
        return false;
      }

      //If we are dealing with the same object type avoid the two legs
      // either being the same object
      if (cndObjTypeVec[0] == cndObjTypeVec[1] && obj0Index == obj1Index && cond0bx == cond1bx) {
        LogDebug("L1TGlobal") << "Corr Condition looking at same leg...skip" << std::endl;
        continue;
      }

      //switch (cond1Categ) {
      if (cond1Categ == CondMuon) {
          lutObj1 = "MU";
          candMuVec = m_uGtB->getCandL1Mu();
          phiIndex1 = (candMuVec->at(cond1bx, obj1Index))->hwPhiAtVtx();  //(*candMuVec)[obj0Index]->phiIndex();
          etaIndex1 = (candMuVec->at(cond1bx, obj1Index))->hwEtaAtVtx();
          etIndex1 = (candMuVec->at(cond1bx, obj1Index))->hwPt();
          etaBin1 = etaIndex1;
          if (etaBin1 < 0)
            etaBin1 = m_gtScales->getMUScales().etaBins.size() + etaBin1;
          //		   LogDebug("L1TGlobal") << "Muon phi" << phiIndex1 << " eta " << etaIndex1 << " etaBin1 = " << etaBin1  << " et " << etIndex1 << std::endl;
	  
          etBin1 = etIndex1;
          int ssize = m_gtScales->getMUScales().etBins.size();
          if (etBin1 >= ssize) {
            LogTrace("L1TGlobal") << "muon2 hw et" << etBin1 << " out of scale range.  Setting to maximum.";
            etBin1 = ssize - 1;
	  }
	  
          // Determine Floating Pt numbers for floating point caluclation
          std::pair<double, double> binEdges = m_gtScales->getMUScales().phiBins.at(phiIndex1);
          phi1Phy = 0.5 * (binEdges.second + binEdges.first);
          binEdges = m_gtScales->getMUScales().etaBins.at(etaBin1);
          eta1Phy = 0.5 * (binEdges.second + binEdges.first);
          binEdges = m_gtScales->getMUScales().etBins.at(etBin1);
          et1Phy = 0.5 * (binEdges.second + binEdges.first); 
        } //break;
      
      else {
	  // Interested only in three-muon correlations
	  LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 1" << std::endl;
          return false;
        } break;
      //}  //end switch on second leg => Switch removed

      // THIRD OBJECT: Finally loop over the third leg to get its information
      for (std::vector<SingleCombInCond>::const_iterator it2Comb = cond2Comb.begin(); it2Comb != cond2Comb.end();
	   it2Comb++) {
	LogDebug("L1TGlobal") << "Looking at the third object for the three-body condition" << std::endl;
	int obj2Index = -1;
	
	if (!(*it2Comb).empty()) {
	  obj2Index = (*it2Comb)[0];
	} else {
	  LogTrace("L1TGlobal") << "\n  SingleCombInCond (*it2Comb).size() " << ((*it2Comb).size()) << std::endl;
	  return false;
	}

	//If we are dealing with the same object type avoid the two legs
	// either being the same object
	if ((cndObjTypeVec[0] == cndObjTypeVec[2] && obj0Index == obj2Index && cond0bx == cond2bx) || 
	    (cndObjTypeVec[1] == cndObjTypeVec[2] && obj1Index == obj2Index && cond1bx == cond2bx)) {
	  LogDebug("L1TGlobal") << "Corr Condition looking at same leg...skip" << std::endl;
	  continue;
	}
	
	if (cond2Categ == CondMuon) {
          lutObj2 = "MU";
          candMuVec = m_uGtB->getCandL1Mu();
          //phiIndex2 = (candMuVec->at(cond2bx, obj2Index))->hwPhiAtVtx();  //(*candMuVec)[obj0Index]->phiIndex();
          etaIndex2 = (candMuVec->at(cond2bx, obj2Index))->hwEtaAtVtx();
          etIndex2 = (candMuVec->at(cond2bx, obj2Index))->hwPt();
          etaBin2 = etaIndex2;
          if (etaBin2 < 0)
            etaBin2 = m_gtScales->getMUScales().etaBins.size() + etaBin2;
	  
          etBin2 = etIndex2;
          int ssize = m_gtScales->getMUScales().etBins.size();
          if (etBin2 >= ssize) {
            LogTrace("L1TGlobal") << "muon3 hw et" << etBin2 << " out of scale range.  Setting to maximum.";
            etBin2 = ssize - 1;
	  }
	  
	  // EF FIXME
          // Determine Floating Pt numbers for floating point calculation
          //std::pair<double, double> binEdges = m_gtScales->getMUScales().phiBins.at(phiIndex2);
          //phi2Phy = 0.5 * (binEdges.second + binEdges.first);
          //binEdges = m_gtScales->getMUScales().etaBins.at(etaBin2);
          //eta2Phy = 0.5 * (binEdges.second + binEdges.first);
          //binEdges = m_gtScales->getMUScales().etBins.at(etBin2);
          //et2Phy = 0.5 * (binEdges.second + binEdges.first); 
        } 
	
	else {
	  // Interested only in three-muon correlations
	  LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 2" << std::endl;
          return false;
        };
	
	//if (m_verbosity) {
	//  LogDebug("L1TGlobal") << "    Correlation pair [" << l1TGtObjectEnumToString(cndObjTypeVec[0]) << ", "
	std::cout << "    Correlation pair [" << l1TGtObjectEnumToString(cndObjTypeVec[0]) << ", "
		  << l1TGtObjectEnumToString(cndObjTypeVec[1]) << "] with collection indices [" << obj0Index
		  << ", " << obj1Index << "] "
		  << " has: \n"
		  << "     Et  value   = [" << etIndex0 << ", " << etIndex1 << "]\n"
		  << "     phi indices = [" << phiIndex0 << ", " << phiIndex1 << "]\n"
		  << "     eta indices = [" << etaIndex0 << ", " << etaIndex1 << "]\n"
		  << std::endl;
	//}

	std::cout << ">>>>>> ELISA: 3MUON EVENT" << std::endl;
	// Now perform the desired correlation on these two objects. Assume true until we find a contradition
	bool reqResult = true;
	
	// clear the indices in the combination
	objectsInComb.clear();
	objectsInComb.push_back(obj0Index);
	objectsInComb.push_back(obj1Index);
	objectsInComb.push_back(obj2Index);
	
	// if we get here all checks were successful for this combination of three muons
	// set the general result for evaluateCondition to "true"
	
	// Delta eta and phi calculations needed to evaluate the three-body invariant mass.
	// ...need to revise this to line up with firmware calculations.
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
	//if (!etSumCond) {
	deltaEtaLUT = m_gtScales->getLUT_DeltaEta(lutName, deltaEtaFW);
	precDeltaEtaLUT = m_gtScales->getPrec_DeltaEta(lutName); //}
	
	//
	LogDebug("L1TGlobal") << "Obj0 phiFW = " << phiIndex0 << " Obj1 phiFW = " << phiIndex1 << "\n"
			      << "    DeltaPhiFW = " << deltaPhiFW << "\n"
			      << "    LUT Name = " << lutName << " Prec = " << precDeltaPhiLUT
			      << "  DeltaPhiLUT = " << deltaPhiLUT << "\n"
			      << "Obj0 etaFW = " << etaIndex0 << " Obj1 etaFW = " << etaIndex1 << "\n"
			      << "    DeltaEtaFW = " << deltaEtaFW << "\n"
			      << "    LUT Name = " << lutName << " Prec = " << precDeltaEtaLUT
			      << "  DeltaEtaLUT = " << deltaEtaLUT << std::endl;
		
	if (corrPar.corrCutType & 0x8 || corrPar.corrCutType & 0x10) {
	  //invariant mass calculation based on
	  // M = sqrt(2*p1*p2(cosh(eta1-eta2) - cos(phi1 - phi2)))
	  // but we calculate (1/2)M^2
	  //
	  double cosDeltaPhiPhy = cos(deltaPhiPhy);
	  double coshDeltaEtaPhy = cosh(deltaEtaPhy);
	  if (corrPar.corrCutType & 0x10)
	    coshDeltaEtaPhy = 1.;
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
	  // if (corrPar.corrCutType & 0x10) coshDeltaEtaLUT=1*pow(10,precCosLUT);
	  
	  std::string lutName = lutObj0;
	  lutName += "-ET";
	  long long ptObj0 = m_gtScales->getLUT_Pt("Mass_" + lutName, etIndex0);
	  unsigned int precPtLUTObj0 = m_gtScales->getPrec_Pt("Mass_" + lutName);
	  
	  lutName = lutObj1;
	  lutName += "-ET";
	  long long ptObj1 = m_gtScales->getLUT_Pt("Mass_" + lutName, etIndex1);
	  unsigned int precPtLUTObj1 = m_gtScales->getPrec_Pt("Mass_" + lutName);

	  // EF FIXME 
	  //lutName = lutObj2;
	  //lutName += "-ET";
	  //long long ptObj2 = m_gtScales->getLUT_Pt("Mass_" + lutName, etIndex2);
	  //unsigned int precPtLUTObj2 = m_gtScales->getPrec_Pt("Mass_" + lutName);
	  
	  // Pt and Angles are at different precission.
	  long long massSq = ptObj0 * ptObj1 * (coshDeltaEtaLUT - cosDeltaPhiLUT);
	  
	  //Note: There is an assumption here that Cos and Cosh have the same precission
	  unsigned int preShift = precPtLUTObj0 + precPtLUTObj1 + precCosLUT - corrPar.precMassCut;
	  
	  std::cout << "####################################\n";
	  //LogDebug("L1TGlobal") << "    Testing Invariant Mass (" << lutObj0 << "," << lutObj1 << ") ["
	  std::cout << "    Testing Invariant Mass (" << lutObj0 << "," << lutObj1 << ") ["
		    << (long long)(corrPar.minMassCutValue * pow(10, preShift)) << ","
		    << (long long)(corrPar.maxMassCutValue * pow(10, preShift))
		    << "] with precision = " << corrPar.precMassCut << "\n"
		    << "    deltaPhiLUT  = " << deltaPhiLUT << "  cosLUT  = " << cosDeltaPhiLUT << "\n"
		    << "    deltaEtaLUT  = " << deltaEtaLUT << "  coshLUT = " << coshDeltaEtaLUT << "\n"
		    << "    etIndex0     = " << etIndex0 << "    pt0LUT      = " << ptObj0
		    << " PhyEt0 = " << et0Phy << "\n"
		    << "    etIndex1     = " << etIndex1 << "    pt1LUT      = " << ptObj1
		    << " PhyEt1 = " << et1Phy << "\n"
		    << "    massSq/2     = " << massSq << "\n"
		    << "    Precision Shift = " << preShift << "\n"
		    << "    massSq   (shift)= " << (massSq / pow(10, preShift + corrPar.precMassCut)) << "\n"
		    << "    deltaPhiPhy  = " << deltaPhiPhy << "  cos() = " << cosDeltaPhiPhy << "\n"
		    << "    deltaEtaPhy  = " << deltaEtaPhy << "  cosh()= " << coshDeltaEtaPhy << "\n"
		    << "    massSqPhy/2  = " << massSqPhy
		    << "  sqrt(|massSq|) = " << sqrt(fabs(2. * massSqPhy)) << std::endl;
	  
	  
	  //if(preShift>0) massSq /= pow(10,preShift);
	  
	  if (massSq >= 0 && massSq >= (long long)(corrPar.minMassCutValue * pow(10, preShift)) &&
	      massSq <= (long long)(corrPar.maxMassCutValue * pow(10, preShift))) {
	    // LogDebug("L1TGlobal") << "    Passed Invariant Mass Cut [" // EF 
	    cout << "    Passed Invariant Mass Cut ["
		 << (long long)(corrPar.minMassCutValue * pow(10, preShift)) << ","
		 << (long long)(corrPar.maxMassCutValue * pow(10, preShift)) << "]" << std::endl;
	    
	  } else {
	    LogDebug("L1TGlobal") << "    Failed Invariant Mass Cut ["
				  << (long long)(corrPar.minMassCutValue * pow(10, preShift)) << ","
				  << (long long)(corrPar.maxMassCutValue * pow(10, preShift)) << "]" << std::endl;
	    reqResult = false;
	  }
	}
	
	if (reqResult) {
	  condResult = true;
	  (combinationsInCond()).push_back(objectsInComb);
	}
	
      }  //end loop over third leg
    }  //end loop over second leg
  }  //end loop over first leg
  
  if (m_verbosity && condResult) {
    LogDebug("L1TGlobal") << " pass(es) the correlation condition.\n" << std::endl;
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
