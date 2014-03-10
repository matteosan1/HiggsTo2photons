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

class GlobeBase : public edm:EDAnalyzer {
public:

GlobeBase(const edm::ParameterSet&) { type_ = "Base"; prefix_ = ""; };
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
};

#endif
