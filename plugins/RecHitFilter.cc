// Package:    InvariantMassFilter
// Class:      InvariantMassFilter
//
// Original Author:  Matteo Sani,40 3-A02,+41227671577,
//         Created:  Thu Jun 22 13:56:58 CEST 2010

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Reconstruction Classes
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

// Geometry
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"

class RecHitFilter : public edm::EDFilter {
public:
  explicit RecHitFilter(const edm::ParameterSet&);
  ~RecHitFilter();

private:
  virtual void beginJob() ;
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  edm::InputTag ebColl_;
  edm::InputTag eeColl_;
  
  Float_t timeCut_;
};


RecHitFilter::RecHitFilter(const edm::ParameterSet& iConfig) {
  
  ebColl_  = iConfig.getParameter<edm::InputTag>("EBCollection");
  eeColl_  = iConfig.getParameter<edm::InputTag>("EECollection");
  
  timeCut_ =  iConfig.getParameter<double>("timeCut");
  
  produces<EcalRecHitCollection> ("timeCleanedRecHits");
}

RecHitFilter::~RecHitFilter()
{}

bool RecHitFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  //edm::ESHandle<CaloGeometry> geoHandle;
  //iSetup.get<CaloGeometryRecord>().get(geoHandle);
  //const CaloGeometry& geometry = *geoHandle;
  //const CaloSubdetectorGeometry *geometry_p;
  //std::auto_ptr<const CaloSubdetectorTopology> topology;
  edm::Handle<EcalRecHitCollection> rhcH[2];
    
  std::auto_ptr<EcalRecHitCollection> hits(new EcalRecHitCollection);
    
  
  for (unsigned int i=0; i<2; i++) {
    if (i==0)
      iEvent.getByLabel(ebColl_, rhcH[i]); 
    else
      iEvent.getByLabel(eeColl_, rhcH[i]);  

    
    const EcalRecHitCollection* recHits = rhcH[i].product();
      
    if (recHits->size() == 0)
      continue;
    
    //if ((*recHits)[0].id().subdetId() == EcalBarrel) {
    //  geometry_p = geometry.getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
    //  topology.reset(new EcalBarrelTopology(geoHandle));
    //} else if ((*recHits)[0].id().subdetId() == EcalEndcap) {
    //  geometry_p = geometry.getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
    //  topology.reset(new EcalEndcapTopology(geoHandle));
    //} else if ((*recHits)[0].id().subdetId() == EcalPreshower) {
    //  geometry_p = geometry.getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
    //  topology.reset(new EcalPreshowerTopology (geoHandle));
    //} else throw(std::runtime_error("\n\nProducer encountered invalied ecalhitcollection type.\n\n"));
      
    
    EcalRecHitCollection::const_iterator it;
    for (it = recHits->begin(); it != recHits->end(); it++){
      if (fabs((*it).time()) > timeCut_)
	continue;
      
      hits->push_back(*it);
    }
  }
  
  iEvent.put(hits, "timeCleanedRecHits");
  
  return true;
}

void RecHitFilter::beginJob()
{}

void RecHitFilter::endJob() 
{}

//define this as a plug-in
DEFINE_FWK_MODULE(RecHitFilter);
