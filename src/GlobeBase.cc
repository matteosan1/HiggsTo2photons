#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeBase.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

GlobeBase::GlobeBase(const edm::ParameterSet& pset) {

  type_ = pset.getParameter<std::string>("Type");
  prefix_ = pset.getParameter<std::string>("Prefix");

  if (pset.exists("DebugLevel"))
    debug_level = pset.getParameter<int>("DebugLevel");
  
  gCUT = new GlobeCuts(pset);
}

void GlobeBase::defineBranch(GlobeAnalyzer* ana) {
  
  if (prefix_ != "")
    prefix_ = "_"+prefix_;
}

GlobeBase::~GlobeBase() {
  delete gCUT;
}
