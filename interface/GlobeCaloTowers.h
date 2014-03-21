#ifndef GLOBECALOTOWERS_H
#define GLOBECALOTOWERS_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeBase.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include <iostream>

class GlobeCaloTowers : public GlobeBase {
 public:
  
  GlobeCaloTowers(const edm::ParameterSet&);
  virtual ~GlobeCaloTowers() {};

  void defineBranch(GlobeAnalyzer* ana);
  void analyze(const edm::Event&, const edm::EventSetup&);

  // variables

  TClonesArray *ct_p4;
  Int_t ct_n;
  Float_t ct_emEnergy[MAX_CALOTOWERS];
  Float_t ct_hadEnergy[MAX_CALOTOWERS];
  Float_t ct_outerEnergy[MAX_CALOTOWERS];
  Int_t ct_emL1[MAX_CALOTOWERS];
  Int_t ct_hadL1[MAX_CALOTOWERS];
  Int_t ct_size[MAX_CALOTOWERS];

 private:
  edm::InputTag calotowerColl; 
};


#endif
