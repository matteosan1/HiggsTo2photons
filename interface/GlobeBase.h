#ifndef GLOBEBASE_H
#define GLOBEBASE_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <string>

class GlobeAnalyzer;

class GlobeBase : public edm::EDAnalyzer {
 public:
  GlobeBase() {};
  GlobeBase(const edm::ParameterSet& pSet, std::string name="") { prefix = name;};
  virtual ~GlobeBase() {};

  virtual void defineBranch(GlobeAnalyzer* ana) {};
  virtual void analyze(const edm::Event&, const edm::EventSetup&) {};

 protected:
  std::string prefix; 
};

#endif
