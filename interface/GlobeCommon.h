#ifndef GLOBECOMMON_H
#define GLOBECOMMON_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeBase.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

class GlobeCommon : public GlobeBase {
 public:
  
  GlobeCommon(const edm::ParameterSet& iConfig);
  ~GlobeCommon() {};

  void defineBranch(GlobeAnalyzer* ana);
  void defineLumiBranch(TTree* tree);
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endLumiBlock(const edm::LuminosityBlock &, const edm::EventSetup &);

  Int_t lumis;
  Int_t bx;
  Int_t event;
  Int_t run;
  Int_t process_id;
  Float_t weight;
  Float_t pthat;
  bool doParticleGun;
  
 private:
  edm::InputTag generatorColl;

};

#endif
