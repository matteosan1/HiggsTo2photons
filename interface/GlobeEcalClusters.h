#ifndef GLOBEECALCLUSTERS_H
#define GLOBEECALCLUSTERS_H

#include "FWCore/Framework/interface/Event.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"

#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/ClusterShapeFwd.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionBaseClass.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterFunctionFactory.h"

#include "Math/VectorUtil.h"
#include <iostream>

class CaloSubdetectorTopology;
class GlobeAnalyzer;

class GlobeEcalClusters {
 public:
  
  GlobeEcalClusters(const edm::ParameterSet&, const char* n="unused");
  virtual ~GlobeEcalClusters() {};

  void defineBranch(GlobeAnalyzer* ana);
  bool analyze(const edm::Event&, const edm::EventSetup&);

protected:
  void analyzeSuperClusters(edm::Handle<reco::SuperClusterCollection> scH, edm::Handle<reco::BasicClusterCollection> bcH,EcalClusterLazyTools& ecalTools,  bool isBarrel, int type);
  void analyzeBasicClusters(edm::Handle<reco::BasicClusterCollection> bcH, EcalClusterLazyTools& ecalTools, int type);
  void fillBasicCluster(reco::BasicClusterRef bc, EcalClusterLazyTools& ecalTools, int type);

public:

  std::vector<float> getESHits(double X, double Y, double Z, std::map<DetId, EcalRecHit> rechits_map, const CaloGeometry& geometry, CaloSubdetectorTopology *topology_p, int row=0);
  std::vector<float> getESShape(std::vector<float> ESHits0);

  /** fills the association of the basic cluster to rechits */
  void fillBasicClusterRecHits(const std::vector<std::pair<DetId,float > > &hits);

  //----------------------------------------
protected:
  /** handles to superclusters */
  edm::Handle<reco::SuperClusterCollection> superClustersHybridH; 
  edm::Handle<reco::SuperClusterCollection> superClustersEndcapH; 

  /** handles to basic clusters */
  edm::Handle<reco::BasicClusterCollection> hybridClustersBarrelH; 
  edm::Handle<reco::BasicClusterCollection> basicClustersEndcapH; 
  edm::Handle<reco::BasicClusterCollection> basicClustersBarrelH;


  /** rechits of the current event (initialized in analyze(..)) */
  const EcalRecHitCollection *barrelRecHits;
  const EcalRecHitCollection *endcapRecHits;

  /** calorimeter topology (for the current event) */
  const CaloTopology *topology;

  const CaloGeometry *geometry;

public:
  //----------------------------------------
  // variables for the output ROOT tree 
  //----------------------------------------

  // SUPER CLUSTERS
  TClonesArray *sc_p4;
  TClonesArray *sc_xyz;
  TClonesArray *bc_p4;
  TClonesArray *bc_xyz;

  Int_t sc_n;
  Float_t sc_pre[MAX_SUPERCLUSTERS];
  Float_t sc_raw[MAX_SUPERCLUSTERS];
  Int_t sc_nbc[MAX_SUPERCLUSTERS];
  Int_t sc_bcseedind[MAX_SUPERCLUSTERS];
  Int_t sc_bcind[MAX_SUPERCLUSTERS][MAX_SUPERCLUSTER_BASICCLUSTERS];
  Int_t sc_barrel[MAX_SUPERCLUSTERS];
  Float_t sc_sphi[MAX_SUPERCLUSTERS];
  Float_t sc_seta[MAX_SUPERCLUSTERS];
  Float_t sc_brem[MAX_SUPERCLUSTERS];
  Float_t sc_r9[MAX_SUPERCLUSTERS];
  Float_t sc_bccrackcorr[MAX_SUPERCLUSTERS][MAX_SUPERCLUSTER_BASICCLUSTERS];
  Float_t sc_bclocalcorr[MAX_SUPERCLUSTERS][MAX_SUPERCLUSTER_BASICCLUSTERS];

  // BASIC CLUSTERS
  Int_t bc_n;
  Int_t bc_nhits[MAX_BASICCLUSTERS];
  Int_t bc_type[MAX_BASICCLUSTERS];
  
  /** first index is the basic cluster index, second index
      is the index of the rechit within this basic cluster,
      value is the detid of the rechit belonging to this
      basic cluster */
  // Int_t bc_hitdetid[bc_n][MAX_ECALRECHITS];
  std::vector<std::vector<Int_t> > *bc_hitdetid;

  Float_t bc_s1[MAX_BASICCLUSTERS];
  Float_t bc_s4[MAX_BASICCLUSTERS];
  Float_t bc_s9[MAX_BASICCLUSTERS];
  Float_t bc_s25[MAX_BASICCLUSTERS];
  Float_t bc_sipip[MAX_BASICCLUSTERS];
  Float_t bc_sieie[MAX_BASICCLUSTERS];
  Float_t bc_sieip[MAX_BASICCLUSTERS];
  Float_t bc_chx[MAX_BASICCLUSTERS];
  
  Float_t bc_2x5_max[MAX_BASICCLUSTERS];
  Float_t bc_5x1_sam[MAX_BASICCLUSTERS];
  Int_t bc_seed[MAX_BASICCLUSTERS];
  
 private:
  const char* nome;
  GlobeCuts *gCUT;
  
  // SUPER CLUSTERS
  edm::InputTag hybridSuperClusterColl; 
  edm::InputTag barrelSuperClusterColl; 
  edm::InputTag endcapSuperClusterColl; 

  // BASIC CLUSTERS
  edm::InputTag barrelHybridClusterColl; 
  edm::InputTag barrelBasicClusterColl; 
  edm::InputTag endcapBasicClusterColl; 
  edm::InputTag barrelHybridClusterShapeColl; 
  edm::InputTag barrelBasicClusterShapeColl; 
  edm::InputTag endcapBasicClusterShapeColl; 
  edm::InputTag ecalHitEBColl;
  edm::InputTag ecalHitEEColl;

  EcalClusterFunctionBaseClass *CrackCorrFunc;
  EcalClusterFunctionBaseClass *LocalCorrFunc;

  int debug_level;
};


#endif
