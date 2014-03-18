#ifndef GLOBEBASE_H
#define GLOBEBASE_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TTree.h"
#include <string>

class GlobeAnalyzer;
class GlobeCuts;

class GlobeBase : public edm:EDAnalyzer {
 public:

  GlobeBase(const edm::ParameterSet&);
  ~GlobeAnalyzer() {};
  
  virtual std::string type() { return type_; };
  virtual std::string prefix() { return prefix_; };
  
  virtual void defineBranch(GlobeAbalyzer* ana) {};
  virtual void defineLumiBranch(TTree* lumitree) {};
  virtual void definePathBranch(TTree* tree2) {};
  
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {};
  
 private:
  std::string type_;
  std::string prefix_;
  GlobeCuts* gCUT;
  int debug_level;
  int order;
};

#endif
