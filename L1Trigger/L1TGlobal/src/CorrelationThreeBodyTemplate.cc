/**
 * \class CorrelationThreeBodyTemplate
 *
 *
 * Description: L1 Global Trigger correlation template.
 * Includes spatial correlation for two objects of different type.
 *
 * Implementation:
 *    <TODO: enter implementation details>
 *
 * \author: Vasile Mihai Ghete - HEPHY Vienna
 *
 * $Date$
 * $Revision$
 *
 */

// this class header
#include "L1Trigger/L1TGlobal/interface/CorrelationThreeBodyTemplate.h"

// system include files

#include <iostream>
#include <iomanip>

// user include files
//   base class

// forward declarations

// constructors
//   default

CorrelationThreeBodyTemplate::CorrelationThreeBodyTemplate() : GlobalCondition() {
  m_condCategory = l1t::CondCorrelationThreeBody; 
  m_condType = l1t::Type3s;
  m_condChipNr = -1;

  // there are in fact two objects
  int nObjects = nrObjects();

  if (nObjects > 0) {
    m_objectType.reserve(nObjects);
  }

  m_cond0Category = l1t::CondNull;
  m_cond1Category = l1t::CondNull;
  m_cond2Category = l1t::CondNull;
  m_cond0Index = -1;
  m_cond1Index = -1;
  m_cond2Index = -1;
}

//   from condition name
CorrelationThreeBodyTemplate::CorrelationThreeBodyTemplate(const std::string& cName) : GlobalCondition(cName) {
  m_condCategory = l1t::CondCorrelation; // FIXME EF: CondCorrelationThreeBody
  m_condType = l1t::Type3s; // FIXME EF: l1t::Type2cor 
  m_condChipNr = -1;

  // there are in fact two objects
  int nObjects = nrObjects();

  if (nObjects > 0) {
    m_objectType.reserve(nObjects);
  }

  m_cond0Category = l1t::CondNull;
  m_cond1Category = l1t::CondNull;
  m_cond2Category = l1t::CondNull;
  m_cond0Index = -1;
  m_cond1Index = -1;
  m_cond2Index = -1;
}

//   from condition name, the category of first sub-condition, the category of the
//   second sub-condition, the index of first sub-condition in the cor* vector,
//   the index of second sub-condition in the cor* vector
CorrelationThreeBodyTemplate::CorrelationThreeBodyTemplate(const std::string& cName,
							   const l1t::GtConditionCategory& cond0Cat,
							   const l1t::GtConditionCategory& cond1Cat,
							   const l1t::GtConditionCategory& cond2Cat,
							   const int cond0Index,
							   const int cond1index,
							   const int cond2index)
  : GlobalCondition(cName),
    m_cond0Category(cond0Cat),
    m_cond1Category(cond1Cat),
    m_cond2Category(cond2Cat),
    m_cond0Index(cond0Index),
    m_cond1Index(cond1index),
    m_cond2Index(cond2index)
    
{
  m_condCategory = l1t::CondCorrelation;
  m_condType = l1t::Type3s; // FIXME EF: l1t::Type2cor 
  m_condChipNr = -1;
  
  // there are in fact two objects
  int nObjects = nrObjects();
  
  if (nObjects > 0) {
    m_objectType.resize(nObjects);
  }
}

// copy constructor
CorrelationThreeBodyTemplate::CorrelationThreeBodyTemplate(const CorrelationThreeBodyTemplate& cp) : GlobalCondition(cp.m_condName) { copy(cp); }

// destructor
CorrelationThreeBodyTemplate::~CorrelationThreeBodyTemplate() {
  // empty now
}

// assign operator
CorrelationThreeBodyTemplate& CorrelationThreeBodyTemplate::operator=(const CorrelationThreeBodyTemplate& cp) {
  copy(cp);
  return *this;
}

// set the category of the two sub-conditions
void CorrelationThreeBodyTemplate::setCond0Category(const l1t::GtConditionCategory& condCateg) { m_cond0Category = condCateg; }
void CorrelationThreeBodyTemplate::setCond1Category(const l1t::GtConditionCategory& condCateg) { m_cond1Category = condCateg; }
void CorrelationThreeBodyTemplate::setCond2Category(const l1t::GtConditionCategory& condCateg) { m_cond2Category = condCateg; }

// set the index of the two sub-conditions in the cor* vector from menu
void CorrelationThreeBodyTemplate::setCond0Index(const int& condIndex) { m_cond0Index = condIndex; }
void CorrelationThreeBodyTemplate::setCond1Index(const int& condIndex) { m_cond1Index = condIndex; }
void CorrelationThreeBodyTemplate::setCond2Index(const int& condIndex) { m_cond2Index = condIndex; }

// set the correlation parameters of the condition
void CorrelationThreeBodyTemplate::setCorrelationThreeBodyParameter(const CorrelationThreeBodyParameter& corrParameter) {
  m_correlationParameter = corrParameter;
}

void CorrelationThreeBodyTemplate::print(std::ostream& myCout) const {
  myCout << "\n  CorrelationThreeBodyTemplate print..." << std::endl;

  GlobalCondition::print(myCout);

  myCout << "\n  First subcondition category:  " << m_cond0Category << std::endl;
  myCout << "  Second subcondition category: " << m_cond1Category << std::endl;
  myCout << "  Third subcondition category: " << m_cond2Category << std::endl;

  myCout << "\n  First subcondition index:  " << m_cond0Index << std::endl;
  myCout << "  Second subcondition index: " << m_cond1Index << std::endl;
  myCout << "  Third subcondition index: " << m_cond2Index << std::endl;

  myCout << "\n  Correlation parameters "
         << "[ hex ]" << std::endl;

  myCout << "    Cut Type:  " << m_correlationParameter.corrCutType << std::endl;
  myCout << "    minEtaCutValue        = " << std::dec << m_correlationParameter.minEtaCutValue << std::endl;
  myCout << "    maxEtaCutValue        = " << std::dec << m_correlationParameter.maxEtaCutValue << std::endl;
  myCout << "    precEtaCut            = " << std::dec << m_correlationParameter.precEtaCut << std::endl;
  myCout << "    minPhiCutValue        = " << std::dec << m_correlationParameter.minPhiCutValue << std::endl;
  myCout << "    maxPhiCutValue        = " << std::dec << m_correlationParameter.maxPhiCutValue << std::endl;
  myCout << "    precPhiCut            = " << std::dec << m_correlationParameter.precPhiCut << std::endl;
  myCout << "    minDRCutValue         = " << std::dec << m_correlationParameter.minDRCutValue << std::endl;
  myCout << "    maxDRCutValue         = " << std::dec << m_correlationParameter.maxDRCutValue << std::endl;
  myCout << "    precDRCut             = " << std::dec << m_correlationParameter.precDRCut << std::endl;
  myCout << "    minMassCutValue       = " << std::dec << m_correlationParameter.minMassCutValue << std::endl;
  myCout << "    maxMassCutValue       = " << std::dec << m_correlationParameter.maxMassCutValue << std::endl;
  myCout << "    precMassCut           = " << std::dec << m_correlationParameter.precMassCut << std::endl;
  myCout << "    minTBPTCutValue       = " << std::dec << m_correlationParameter.minTBPTCutValue << std::endl;
  myCout << "    maxTBPTCutValue       = " << std::dec << m_correlationParameter.maxTBPTCutValue << std::endl;
  myCout << "    precTBPTCut           = " << std::dec << m_correlationParameter.precTBPTCut << std::endl;
  myCout << "    chargeCorrelation  = " << std::dec << m_correlationParameter.chargeCorrelation << std::endl;

  // reset to decimal output
  myCout << std::dec << std::endl;
}

void CorrelationThreeBodyTemplate::copy(const CorrelationThreeBodyTemplate& cp) {
  m_condName = cp.condName();
  m_condCategory = cp.condCategory();
  m_condType = cp.condType();
  m_objectType = cp.objectType();
  m_condGEq = cp.condGEq();
  m_condChipNr = cp.condChipNr();

  m_cond0Category = cp.cond0Category();
  m_cond1Category = cp.cond1Category();
  m_cond2Category = cp.cond2Category();
  m_cond0Index = cp.cond0Index();
  m_cond1Index = cp.cond1Index();
  m_cond2Index = cp.cond2Index();

  m_correlationParameter = *(cp.correlationParameter());
}

// output stream operator
std::ostream& operator<<(std::ostream& os, const CorrelationThreeBodyTemplate& result) {
  result.print(os);
  return os;
}
