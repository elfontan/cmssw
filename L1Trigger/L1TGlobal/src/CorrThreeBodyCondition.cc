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

  // evaluate first the two sub-conditions (Type1s)

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

  // second object
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

  // third object
  reqObjResult = false;

  if (cond2Categ == CondMuon) {
      corrMuon = static_cast<const MuonTemplate*>(m_gtCond2);
      MuCondition muCondition(
			      corrMuon, m_uGtB, 0, 0);  //BLW these are counts that don't seem to be used...perhaps remove
      
      muCondition.evaluateConditionStoreResult(bxEval);
      reqObjResult = muCondition.condLastResult();
      
      cond2Comb = (muCondition.getCombinationsInCond());
      cond1bx = bxEval + (corrMuon->condRelativeBx()); //FIXME
      
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
  // /*
  //Vector to store muon indexes:
  vector<int> muIndexes;
  //Vector to store dimuon invariant masses:
  vector<long long> muInvMass;  
  muIndexes.clear();
  muInvMass.clear();
  //*/

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

  int etaIndex0 = 0;
  double eta0Phy = 0.;
  int etaIndex1 = 0;
  double eta1Phy = 0.;
  int etaBin0 = 0;
  int etaBin1 = 0;

  int etIndex0 = 0;
  int etBin0 = 0;
  double et0Phy = 0.;
  int etIndex1 = 0;
  int etBin1 = 0;
  double et1Phy = 0.;

  int chrg0 = -1;
  int chrg1 = -1;

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

  LogTrace("L1TGlobal") << "  Sub-condition 0: std::vector<SingleCombInCond> size: " << (cond0Comb.size()) << std::endl;
  LogTrace("L1TGlobal") << "  Sub-condition 1: std::vector<SingleCombInCond> size: " << (cond1Comb.size()) << std::endl;
  LogTrace("L1TGlobal") << "  Sub-condition 2: std::vector<SingleCombInCond> size: " << (cond2Comb.size()) << std::endl;

  // loop over all combinations which produced individually "true" as Type1s
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

    // Collect the information on the first leg of the correlation
    //switch (cond0Categ) {
    if (cond0Categ == CondMuon) {
        lutObj0 = "MU";
        candMuVec = m_uGtB->getCandL1Mu();
        phiIndex0 = (candMuVec->at(cond0bx, obj0Index))->hwPhiAtVtx();  //(*candMuVec)[obj0Index]->phiIndex();
        etaIndex0 = (candMuVec->at(cond0bx, obj0Index))->hwEtaAtVtx();
        etIndex0 = (candMuVec->at(cond0bx, obj0Index))->hwPt();
        chrg0 = (candMuVec->at(cond0bx, obj0Index))->hwCharge();
        int etaBin0 = etaIndex0;
        if (etaBin0 < 0)
          etaBin0 = m_gtScales->getMUScales().etaBins.size() + etaBin0;  //twos complement
        //		LogDebug("L1TGlobal") << "Muon phi" << phiIndex0 << " eta " << etaIndex0 << " etaBin0 = " << etaBin0  << " et " << etIndex0 << std::endl;
	
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


    // Now loop over the second leg to get its information
    for (std::vector<SingleCombInCond>::const_iterator it1Comb = cond1Comb.begin(); it1Comb != cond1Comb.end();
         it1Comb++) {
      LogDebug("L1TGlobal") << "Looking at second Condition" << std::endl;
      // Type1s: there is 1 object only, no need for a loop (*it1Comb)[0]
      // ... but add protection to not crash
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
          chrg1 = (candMuVec->at(cond1bx, obj1Index))->hwCharge();
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
	  LogDebug("L1TGlobal") << "CondMuon not satisfied for Leg 0" << std::endl;
          return false;
        } break;
      //}  //end switch on second leg => Switch removed


      //if (m_verbosity) {
      //  LogDebug("L1TGlobal") << "    Correlation pair [" << l1TGtObjectEnumToString(cndObjTypeVec[0]) << ", "
      std::cout << "    Correlation pair [" << l1TGtObjectEnumToString(cndObjTypeVec[0]) << ", "
                              << l1TGtObjectEnumToString(cndObjTypeVec[1]) << "] with collection indices [" << obj0Index
                              << ", " << obj1Index << "] "
                              << " has: \n"
                              << "     Et  value   = [" << etIndex0 << ", " << etIndex1 << "]\n"
                              << "     phi indices = [" << phiIndex0 << ", " << phiIndex1 << "]\n"
                              << "     eta indices = [" << etaIndex0 << ", " << etaIndex1 << "]\n"
                              << "     chrg        = [" << chrg0 << ", " << chrg1 << "]\n"
                              << std::endl;
      //}

      // Now perform the desired correlation on these two objects. Assume true until we find a contradition
      bool reqResult = true;

      // clear the indices in the combination
      objectsInComb.clear();

      objectsInComb.push_back(obj0Index);
      objectsInComb.push_back(obj1Index);

      // if we get here all checks were successful for this combination
      // set the general result for evaluateCondition to "true"

      // These all require some delta eta and phi calculations.  Do them first...for now real calculation but need to
      // revise this to line up with firmware calculations.
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

      // If there is a delta eta, check it.
      if (corrPar.corrCutType & 0x1) {
        unsigned int preShift = precDeltaEtaLUT - corrPar.precEtaCut;
        LogDebug("L1TGlobal") << "    Testing Delta Eta Cut (" << lutObj0 << "," << lutObj1 << ") ["
                              << (long long)(corrPar.minEtaCutValue * pow(10, preShift)) << ","
                              << (long long)(corrPar.maxEtaCutValue * pow(10, preShift))
                              << "] with precision = " << corrPar.precEtaCut << "\n"
                              << "    deltaEtaLUT = " << deltaEtaLUT << "\n"
                              << "    Precision Shift = " << preShift << "\n"
                              << "    deltaEta (shift)= " << (deltaEtaLUT / pow(10, preShift + corrPar.precEtaCut))
                              << "\n"
                              << "    deltaEtaPhy = " << deltaEtaPhy << std::endl;

        //if(preShift>0) deltaEtaLUT /= pow(10,preShift);
        if (deltaEtaLUT >= (long long)(corrPar.minEtaCutValue * pow(10, preShift)) &&
            deltaEtaLUT <= (long long)(corrPar.maxEtaCutValue * pow(10, preShift))) {
          LogDebug("L1TGlobal") << "    Passed Delta Eta Cut ["
                                << (long long)(corrPar.minEtaCutValue * pow(10, preShift)) << ","
                                << (long long)(corrPar.maxEtaCutValue * pow(10, preShift)) << "]" << std::endl;

        } else {
          LogDebug("L1TGlobal") << "    Failed Delta Eta Cut ["
                                << (long long)(corrPar.minEtaCutValue * pow(10, preShift)) << ","
                                << (long long)(corrPar.maxEtaCutValue * pow(10, preShift)) << "]" << std::endl;
          reqResult = false;
        }
      }

      //if there is a delta phi check it.
      if (corrPar.corrCutType & 0x2) {
        unsigned int preShift = precDeltaPhiLUT - corrPar.precPhiCut;
        LogDebug("L1TGlobal") << "    Testing Delta Phi Cut (" << lutObj0 << "," << lutObj1 << ") ["
                              << (long long)(corrPar.minPhiCutValue * pow(10, preShift)) << ","
                              << (long long)(corrPar.maxPhiCutValue * pow(10, preShift))
                              << "] with precision = " << corrPar.precPhiCut << "\n"
                              << "    deltaPhiLUT = " << deltaPhiLUT << "\n"
                              << "    Precision Shift = " << preShift << "\n"
                              << "    deltaPhi (shift)= " << (deltaPhiLUT / pow(10, preShift + corrPar.precPhiCut))
                              << "\n"
                              << "    deltaPhiPhy = " << deltaPhiPhy << std::endl;

        //if(preShift>0) deltaPhiLUT /= pow(10,preShift);
        if (deltaPhiLUT >= (long long)(corrPar.minPhiCutValue * pow(10, preShift)) &&
            deltaPhiLUT <= (long long)(corrPar.maxPhiCutValue * pow(10, preShift))) {
          LogDebug("L1TGlobal") << "    Passed Delta Phi Cut ["
                                << (long long)(corrPar.minPhiCutValue * pow(10, preShift)) << ","
                                << (long long)(corrPar.maxPhiCutValue * pow(10, preShift)) << "]" << std::endl;

        } else {
          LogDebug("L1TGlobal") << "    Failed Delta Phi Cut ["
                                << (long long)(corrPar.minPhiCutValue * pow(10, preShift)) << ","
                                << (long long)(corrPar.maxPhiCutValue * pow(10, preShift)) << "]" << std::endl;
          reqResult = false;
        }
      }

      if (corrPar.corrCutType & 0x4) {
        //Assumes Delta Eta and Delta Phi LUTs have the same precision
        unsigned int preShift = 2 * precDeltaPhiLUT - corrPar.precDRCut;
        double deltaRSqPhy = deltaPhiPhy * deltaPhiPhy + deltaEtaPhy * deltaEtaPhy;
        long long deltaRSq = deltaEtaLUT * deltaEtaLUT + deltaPhiLUT * deltaPhiLUT;

        LogDebug("L1TGlobal") << "    Testing Delta R Cut (" << lutObj0 << "," << lutObj1 << ") ["
                              << (long long)(corrPar.minDRCutValue * pow(10, preShift)) << ","
                              << (long long)(corrPar.maxDRCutValue * pow(10, preShift))
                              << "] with precision = " << corrPar.precDRCut << "\n"
                              << "    deltaPhiLUT = " << deltaPhiLUT << "\n"
                              << "    deltaEtaLUT = " << deltaEtaLUT << "\n"
                              << "    deltaRSqLUT = " << deltaRSq << "\n"
                              << "    Precision Shift = " << preShift << "\n"
                              << "    deltaRSqLUT (shift)= " << (deltaRSq / pow(10, preShift + corrPar.precDRCut))
                              << "\n"
                              << "    deltaRSqPhy = " << deltaRSqPhy << std::endl;

        //if(preShift>0) deltaRSq /= pow(10,preShift);
        if (deltaRSq >= (long long)(corrPar.minDRCutValue * pow(10, preShift)) &&
            deltaRSq <= (long long)(corrPar.maxDRCutValue * pow(10, preShift))) {
          LogDebug("L1TGlobal") << "    Passed Delta R Cut [" << (long long)(corrPar.minDRCutValue * pow(10, preShift))
                                << "," << (long long)(corrPar.maxDRCutValue * pow(10, preShift)) << "]" << deltaRSq
                                << std::endl;

        } else {
          LogDebug("L1TGlobal") << "    Failed Delta R Cut [" << (int)(corrPar.minDRCutValue * pow(10, preShift)) << ","
                                << (long long)(corrPar.maxDRCutValue * pow(10, preShift)) << "]" << deltaRSq
                                << std::endl;
          reqResult = false;
        }
      }

      if (corrPar.corrCutType & 0x20) {
        // Two body pt: pt^2 = pt1^2+pt2^2+2*pt1*pt2*(cos(phi1)*cos(phi2)+sin(phi1)*sin(phi2)).

        LogDebug("L1TGlobal") << " corrPar.corrCutType: " << corrPar.corrCutType << "\n";

        //calculate math sins and cosines for debugging
        double cosPhi1Phy = cos(phi0Phy);
        double sinPhi1Phy = sin(phi0Phy);
        double cosPhi2Phy = cos(phi1Phy);
        double sinPhi2Phy = sin(phi1Phy);

        double tbptSqPhy = et0Phy * et0Phy + et1Phy * et1Phy +
                           2 * et0Phy * et1Phy * (cosPhi1Phy * cosPhi2Phy + sinPhi1Phy * sinPhi2Phy);
        // get values from LUT's

        const std::string& lutName0 = lutObj0;
        unsigned int precCosLUT0 = m_gtScales->getPrec_Cos(lutName0);
        unsigned int precSinLUT0 = m_gtScales->getPrec_Sin(lutName0);

        const std::string& lutName1 = lutObj1;
        unsigned int precCosLUT1 = m_gtScales->getPrec_Cos(lutName1);
        unsigned int precSinLUT1 = m_gtScales->getPrec_Sin(lutName1);

        if (precCosLUT0 - precCosLUT1 != 0)
          LogDebug("L1TGlobal") << "Warning: Cos LUTs for TwoBodyPt on different Precision" << std::endl;
        if (precSinLUT0 - precSinLUT1 != 0)
          LogDebug("L1TGlobal") << "Warning: Sin LUTs for TwoBodyPt on different Precision" << std::endl;
        if (precSinLUT0 - precCosLUT1 != 0)
          LogDebug("L1TGlobal") << "Warning: Sin and Cos LUTs for TwoBodyPt on different Precision" << std::endl;
        if (precSinLUT1 - precCosLUT0 != 0)
          LogDebug("L1TGlobal") << "Warning: Sin and Cos LUTs for TwoBodyPt on different Precision" << std::endl;

        long long cosPhi1LUT = m_gtScales->getLUT_Cos(lutName0, phiIndex0);
        long long sinPhi1LUT = m_gtScales->getLUT_Sin(lutName0, phiIndex0);

        long long cosPhi2LUT = m_gtScales->getLUT_Cos(lutName1, phiIndex1);
        long long sinPhi2LUT = m_gtScales->getLUT_Sin(lutName1, phiIndex1);

        // now get pt LUTs
        std::string lutName = lutObj0;
        lutName += "-ET";
        long long ptObj0 = m_gtScales->getLUT_Pt("TwoBody_" + lutName, etIndex0);
        unsigned int precPtLUTObj0 = m_gtScales->getPrec_Pt("TwoBody_" + lutName);

        lutName = lutObj1;
        lutName += "-ET";
        long long ptObj1 = m_gtScales->getLUT_Pt("TwoBody_" + lutName, etIndex1);
        unsigned int precPtLUTObj1 = m_gtScales->getPrec_Pt("TwoBody_" + lutName);

        LogTrace("L1TGlobal") << " TBPT Trig precisions:\t " << precCosLUT0 << "\t" << precCosLUT1 << "\t"
                              << precSinLUT0 << "\t" << precSinLUT1;
        LogTrace("L1TGlobal") << " TBPT Pt precisions:\t " << precPtLUTObj0 << "\t" << precPtLUTObj1;
        LogTrace("L1TGlobal") << " TBPT Pt cut:\t " << corrPar.minTBPTCutValue << "\tPrecTBPTCut\t"
                              << corrPar.precTBPTCut;
        LogTrace("L1TGlobal") << " TBPT Pt1*Pt1 -- Phys:\t " << et0Phy * et0Phy << "\tHW:\t"
                              << ptObj0 * ptObj0 * (pow(10, 6));
        LogTrace("L1TGlobal") << " TBPT Pt2*Pt2 -- Phys:\t " << et1Phy * et1Phy << "\tHW:\t"
                              << ptObj1 * ptObj1 * (pow(10, 6));
        LogTrace("L1TGlobal") << " TBPT 2Pt1*Pt2 -- Phys:\t " << 2 * et0Phy * et1Phy << "\tHW:\t"
                              << 2 * (ptObj0 * pow(10, 0)) * (ptObj1 * pow(10, 0));
        LogTrace("L1TGlobal") << " TBPT Trig -- Phys:\t " << cosPhi1Phy * cosPhi2Phy + sinPhi1Phy * sinPhi2Phy
                              << "\tHW:\t" << cosPhi1LUT * cosPhi2LUT + sinPhi1LUT * sinPhi2LUT;

        //double tbptSqPhy =   et0Phy*et0Phy             + et1Phy*et1Phy + 2*et0Phy*et1Phy*(cosPhi1Phy*cosPhi2Phy + sinPhi1Phy*sinPhi2Phy);
        long long tbptSqHW = ptObj0 * ptObj0 * (pow(10, 2 * precCosLUT0)) +
                             ptObj1 * ptObj1 * (pow(10, 2 * precCosLUT0)) +
                             2 * ptObj0 * ptObj1 * (cosPhi1LUT * cosPhi2LUT + sinPhi1LUT * sinPhi2LUT);

        unsigned int preShift = precPtLUTObj0 + precPtLUTObj1 + 2 * precCosLUT0;

        LogTrace("L1TGlobal") << "TBPT Result -- Phys: " << tbptSqPhy << "\tHW: " << tbptSqHW << "\tShifted\t"
                              << tbptSqHW / pow(10, preShift) << std::endl;

        preShift = preShift - corrPar.precTBPTCut;

        LogDebug("L1TGlobal")
            << "    Testing Two Body Pt Cut (" << lutObj0 << "," << lutObj1 << ") ["
            << (long long)(corrPar.minTBPTCutValue * pow(10, preShift)) << ","
            << (long long)(corrPar.maxTBPTCutValue * pow(10, preShift)) << "] with precision = " << corrPar.precTBPTCut
            << "\n"
            << "    etIndex0     = " << etIndex0 << "    pt0LUT      = " << ptObj0 << " PhyEt0 = " << et0Phy << "\n"
            << "    etIndex1     = " << etIndex1 << "    pt1LUT      = " << ptObj1 << " PhyEt1 = " << et1Phy << "\n"
            << "    Precision Shift = " << preShift << "\n"
            << "    Sin(phi1): LUT/Phys\t " << sinPhi1LUT << " / " << sinPhi1Phy << "\n"
            << "    Sin(phi2): LUT/Phys\t " << sinPhi2LUT << " / " << sinPhi2Phy << "\n"
            << "    Cos(phi1): LUT/Phys\t " << cosPhi1LUT << " / " << cosPhi1Phy << "\n"
            << "    Cos(phi2): LUT/Phys\t " << cosPhi2LUT << " / " << cosPhi2Phy
            << "\n"

            //    << "    deltaPhiLUT = " << deltaPhiLUT << "\n"
            //    << "    deltaEtaLUT = " << deltaEtaLUT << "\n"
            //    << "    deltaRSqLUT = " << deltaRSq <<  "\n"
            //    << "    Precision Shift = " << preShift << "\n"
            //    << "    deltaRSqLUT (shift)= " << (deltaRSq/pow(10,preShift+corrPar.precDRCut))	<< "\n"
            //    << "    deltaRSqPhy = " << deltaRSqPhy
            << std::endl;

        if (tbptSqHW > 0. && tbptSqHW >= (long long)(corrPar.minTBPTCutValue * pow(10, preShift))) {
          LogDebug("L1TGlobal") << "    Passed Two Body pT Cut ["
                                << (long long)(corrPar.minTBPTCutValue * pow(10, preShift)) << "]"
                                << "\twith value: " << tbptSqHW << "\n"
                                << "\tPhysics Cut[" << corrPar.minTBPTCutValue / pow(10, corrPar.precTBPTCut)
                                << "]\tPhysics Value: " << tbptSqPhy << std::endl;

        } else {
          LogDebug("L1TGlobal") << "    Failed Two Body pT Cut ["
                                << (long long)(corrPar.minTBPTCutValue * pow(10, preShift)) << "]"
                                << "\t with value: " << tbptSqHW << "\n"
                                << "\tPhysics Cut[" << corrPar.minTBPTCutValue / pow(10, corrPar.precTBPTCut)
                                << "]\tPhysics Value: " << tbptSqPhy << std::endl;
          reqResult = false;
        }
      }

      std::cout << ">>>>>>>>>>>>> ELISA: " << corrPar.corrCutType << endl; //EF
      if (corrPar.corrCutType & 0x8 || corrPar.corrCutType & 0x10) {
	std::cout << ">>>>>>>>>>>>> ELISA: PASSED!" << endl; //EF
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
          //EF LogDebug("L1TGlobal") << "    Passed Invariant Mass Cut [" 
          cout << "    Passed Invariant Mass Cut ["
                                << (long long)(corrPar.minMassCutValue * pow(10, preShift)) << ","
                                << (long long)(corrPar.maxMassCutValue * pow(10, preShift)) << "]" << std::endl;

        } else {
          LogDebug("L1TGlobal") << "    Failed Invariant Mass Cut ["
                                << (long long)(corrPar.minMassCutValue * pow(10, preShift)) << ","
                                << (long long)(corrPar.maxMassCutValue * pow(10, preShift)) << "]" << std::endl;
          reqResult = false;
        }
	
	//EF /*
	if (cond0Categ == CondMuon && cond1Categ == CondMuon) {
	  std::cout << "Object pair is: " << lutObj0 << ", " << lutObj1 << std::endl;
	  muIndexes.push_back(etIndex0);
	  muIndexes.push_back(etIndex1);
	  muInvMass.push_back(massSq);
	}
	//EF */
      }
      //EF /*
      std::cout << muIndexes.size() << std::endl;
      std::cout << muInvMass.size() << std::endl;

      if (muIndexes.size() != 0)
	{
	  for (long unsigned int i=0; i < muIndexes.size(); i++)
	    { std::cout << "Position " << i << "  with index " << muIndexes.at(i) <<  std::endl; }
	  for (long unsigned int m=0; m < muInvMass.size(); m++)
	    { std::cout << "Position " << m << "  with index " << muInvMass.at(m) <<  std::endl; }
	}
      //EF */

      // For Muon-Muon Correlation Check the Charge Correlation if requested
      bool chrgCorrel = true;
      if (cond0Categ == CondMuon && cond1Categ == CondMuon) {
        // Check for like-sign
        if (corrPar.chargeCorrelation == 2 && chrg0 != chrg1)
          chrgCorrel = false;
        // Check for opp-sign
        if (corrPar.chargeCorrelation == 4 && chrg0 == chrg1)
          chrgCorrel = false;
      }

      if (reqResult & chrgCorrel) {
        condResult = true;
        (combinationsInCond()).push_back(objectsInComb);
      }

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
