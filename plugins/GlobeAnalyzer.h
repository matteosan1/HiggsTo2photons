// -*- C++ -*-
//
// Package:    GlobeAnalyzer
// Class:      GlobeAnalyzer
// 
/**\class GlobeAnalyzer GlobeAnalyzer.cc HiggsAnalysis/HiggsTo2photons/src/GlobeAnalyzer.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Matteosan SANI
//         Created:  Thu Feb  7 10:14:43 CET 2008
// $Id: GlobeAnalyzer.h,v 1.2 2012/04/23 21:04:36 sani Exp $
//
//

#ifndef GLOBEANALYZER_H
#define GLOBEANALYZER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCommon.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeElectrons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeMuons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePhotons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeConversions.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalClusters.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCaloTowers.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeHcal.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeL1.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeVertex.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeMET.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePFCandidates.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeSimHits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeSimTracks.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeTracks.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGsfTracks.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeTrackingParticles.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenerator.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenParticles.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenVertices.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenJets.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeJets.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalHits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeSelector.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeHLT.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeLeptons.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeHT.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeReducedGen.h"
//#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePAT.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePdfWeights.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeRho.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePileup.h"

#include "TFile.h"
#include "TTree.h"

// system include files
#include <memory>
#include <vector>
#include <string>
#include <algorithm>

typedef std::map<std::string, GlobeBase*(*)()> map_type;

class GlobeAnalyzer : public edm::EDAnalyzer {
public:
  explicit GlobeAnalyzer(const edm::ParameterSet&);
  ~GlobeAnalyzer();
  void readConfiguration(const edm::ParameterSet& iConfig);
  void defineBranch();
  void fillTree();

  int findModule(std::string type, std::string prefix="");

  template <class T> void Branch(const char* name, T* address, const char* leaflist, Int_t bufsize = 32000) {
    for (unsigned int i=0; i<branchesToSkim.size(); i++) {
      if (strcmp(name, branchesToSkim[i].c_str()) == 0) {
	std::cout << "Switching off " << branchesToSkim[i] << " branch." << std::endl;
	return;
      }
    }
    
    tree->Branch(name, address, leaflist, bufsize);
  }

  template <class T> void Branch(const char* name, const char* classname, T** obj, Int_t bufsize = 32000, Int_t splitlevel = 99) {
  
     for (unsigned int i=0; i<branchesToSkim.size(); i++) {
       if (strcmp(name, branchesToSkim[i].c_str()) == 0) {
	 std::cout << "Switching off " << branchesToSkim[i] << " branch." << std::endl;
	 return;
       }
     }
     tree->Branch(name, classname, obj, bufsize, splitlevel);
  }

  std::vector<GlobeBase*> container;

private:
  template<typename T> GlobeBase* createInstance(const edm::ParameterSet& pset) { return new T(pset); }
  void defineMap();
  void beginJob();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endLuminosityBlock(const edm::LuminosityBlock &, const edm::EventSetup &);
  void endJob();

  std::string fileName;
      
  TFile *file;
  TTree *tree, *tree2, *lumitree;

  std::vector<std::string> *parameters;
  std::vector<std::string> *hlt_path_names, *reduced_path;
  std::vector<int>* reduced_index;
  std::string jobmaker;
  std::vector<std::string> globalCountersNames;
  std::vector<std::string> branchesToSkim;

  std::vector<int> globalCounters, globalCountersPerLumi;
  
  int version, type, sel_events, tot_events; 
  int selector_bits;

  int debug_level;
  Int_t nProcessedEvents;

  map_type map;
};

#endif
