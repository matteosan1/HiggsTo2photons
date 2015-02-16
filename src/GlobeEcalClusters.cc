#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalClusters.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeEcalHits.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

//#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include <vector>
#include <algorithm>
#include <numeric>

//----------------------------------------------------------------------

GlobeEcalClusters::GlobeEcalClusters(const edm::ParameterSet& iConfig, const char* n): 
  bc_hitdetid(new vector<vector<Int_t> >), 
  nome(n) {
  
  debug_level = iConfig.getParameter<int>("Debug_Level");
  gCUT = new GlobeCuts(iConfig);
  
  edm::ParameterSet psetSC = iConfig.getParameter<edm::ParameterSet>("SuperClusterCuts");
  edm::ParameterSet psetBC = iConfig.getParameter<edm::ParameterSet>("BasicClusterCuts");
  
  //super clusters
  hybridSuperClusterColl = iConfig.getParameter<edm::InputTag>("HybridSuperClusterColl");
  endcapSuperClusterColl = iConfig.getParameter<edm::InputTag>("EndcapSuperClusterColl");
  
  //basic clusters
  barrelHybridClusterColl = iConfig.getParameter<edm::InputTag>("BarrelHybridClusterColl");
  barrelBasicClusterColl = iConfig.getParameter<edm::InputTag>("BarrelBasicClusterColl");
  endcapBasicClusterColl = iConfig.getParameter<edm::InputTag>("EndcapBasicClusterColl");
  ecalHitEBColl = iConfig.getParameter<edm::InputTag>("EcalHitEBColl");
  ecalHitEEColl = iConfig.getParameter<edm::InputTag>("EcalHitEEColl"); 

  CrackCorrFunc = EcalClusterFunctionFactory::get()->create("EcalClusterCrackCorrection", iConfig);
  LocalCorrFunc = EcalClusterFunctionFactory::get()->create("EcalClusterLocalContCorrection",iConfig);
}

//----------------------------------------------------------------------

void GlobeEcalClusters::defineBranch(GlobeAnalyzer* ana) {

  char a2[100];
  sc_p4 = new TClonesArray("TLorentzVector", MAX_SUPERCLUSTERS);
  sc_xyz = new TClonesArray("TVector3", MAX_SUPERCLUSTERS);
  bc_p4 = new TClonesArray("TLorentzVector", MAX_BASICCLUSTERS);
  bc_xyz = new TClonesArray("TVector3", MAX_SUPERCLUSTERS);
  
  //SC hybrid in barrel and island in endcap
  ana->Branch("sc_n", &sc_n, "sc_n/I");
  ana->Branch("sc_p4", "TClonesArray", &sc_p4, 32000, 0);
  ana->Branch("sc_xyz", "TClonesArray", &sc_xyz, 32000, 0);
  ana->Branch("sc_pre", &sc_pre, "sc_pre[sc_n]/F");
  ana->Branch("sc_raw", &sc_raw, "sc_raw[sc_n]/F");
  ana->Branch("sc_barrel", &sc_barrel, "sc_barrel[sc_n]/I");
  ana->Branch("sc_sphi", &sc_sphi, "sc_sphi[sc_n]/F");
  ana->Branch("sc_seta", &sc_seta, "sc_seta[sc_n]/F");
  ana->Branch("sc_brem", &sc_brem, "sc_brem[sc_n]/F");
  ana->Branch("sc_r9", &sc_r9, "sc_r9[sc_n]/F");
  ana->Branch("sc_nbc", &sc_nbc, "sc_nbc[sc_n]/I");
  ana->Branch("sc_bcseedind", &sc_bcseedind, "sc_bcseedind[sc_n]/I");
  sprintf (a2, "sc_bcind[sc_n][%d]/I", MAX_SUPERCLUSTER_BASICCLUSTERS);
  ana->Branch("sc_bcind", &sc_bcind, a2);
  sprintf (a2, "sc_bccrackcorr[sc_n][%d]/F", MAX_SUPERCLUSTER_BASICCLUSTERS);
  ana->Branch("sc_bccrackcorr",&sc_bccrackcorr, a2);
  sprintf (a2, "sc_bclocalcorr[sc_n][%d]/F", MAX_SUPERCLUSTER_BASICCLUSTERS);
  ana->Branch("sc_bclocalcorr",&sc_bclocalcorr, a2);

  //basic clusters
  ana->Branch("bc_n", &bc_n, "bc_n/I");
  ana->Branch("bc_p4", "TClonesArray", &bc_p4, 32000, 0);
  ana->Branch("bc_xyz", "TClonesArray", &bc_xyz, 32000, 0);
  ana->Branch("bc_nhits", &bc_nhits,"bc_nhits[bc_n]/I");

  // see http://root.cern.ch/phpBB3/viewtopic.php?t=8185
  ana->Branch("bc_hitdetid","vector<vector<Int_t> >",&bc_hitdetid);

  ana->Branch("bc_s1", &bc_s1, "bc_s1[bc_n]/F");
  ana->Branch("bc_chx", &bc_chx, "bc_chx[bc_n]/F");
  ana->Branch("bc_s4", &bc_s4, "bc_s4[bc_n]/F");
  ana->Branch("bc_s9", &bc_s9, "bc_s9[bc_n]/F");
  ana->Branch("bc_s25", &bc_s25, "bc_s25[bc_n]/F");
  ana->Branch("bc_sipip", &bc_sipip, "bc_sipip[bc_n]/F");
  ana->Branch("bc_sieie", &bc_sieie, "bc_sieie[bc_n]/F");
  ana->Branch("bc_sieip", &bc_sieip, "bc_sieip[bc_n]/F");
  ana->Branch("bc_type", &bc_type, "bc_type[bc_n]/I");//type 1 = hybrid, 2 = island endcap, 3 = island barrel.
}

//----------------------------------------------------------------------

bool GlobeEcalClusters::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // get the collection geometry:
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  geometry = geoHandle.product();

  // Get EcalRecHits
  edm::Handle<EBRecHitCollection> pEBRecHitH;
  edm::Handle<EERecHitCollection> pEERecHitH;
  //edm::Handle<ESRecHitCollection> pESRecHitH; 
  
  iEvent.getByLabel(ecalHitEBColl, pEBRecHitH);
  iEvent.getByLabel(ecalHitEEColl, pEERecHitH);
  
  barrelRecHits = pEBRecHitH.product();
  endcapRecHits = pEERecHitH.product();

  EcalClusterLazyTools lazyTool(iEvent, iSetup, ecalHitEBColl, ecalHitEEColl);   

  //edm::ESHandle<CaloTopology> theCaloTopo;
  //edm::ESHandle<CaloTopology> pTopology;
  //iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  //topology = theCaloTopo.product();

  iEvent.getByLabel(hybridSuperClusterColl,superClustersHybridH);
  iEvent.getByLabel(endcapSuperClusterColl, superClustersEndcapH);

  //const std::string instan="hybridBarrelBasicClusters";
  iEvent.getByLabel(barrelHybridClusterColl, hybridClustersBarrelH);
  if (barrelBasicClusterColl.encode() != "") 
    iEvent.getByLabel(barrelBasicClusterColl, basicClustersBarrelH);

  iEvent.getByLabel(endcapBasicClusterColl, basicClustersEndcapH);

  CrackCorrFunc->init(iSetup);
  LocalCorrFunc->init(iSetup);

  if (debug_level > 9) {
    std::cout << "GlobeEcalClusters: superClustersEndcapH->size() "<< superClustersEndcapH->size() << std::endl;
    std::cout << "GlobeEcalClusters: hybridClustersBarrelH->size() "<< hybridClustersBarrelH->size() << std::endl;
    std::cout << "GlobeEcalClusters: basicClustersEndcapH->size() "<< basicClustersEndcapH->size() << std::endl;
  }
  
  //----------------------------------------    
  // analyze super clusters
  //----------------------------------------  
  bc_p4->Clear();
  bc_xyz->Clear();
  bc_hitdetid->clear();

  sc_p4->Clear();
  sc_xyz->Clear();

  bc_n = 0;  
  sc_n = 0;
  
  analyzeSuperClusters(superClustersHybridH, hybridClustersBarrelH, lazyTool, true, 1);
  analyzeSuperClusters(superClustersEndcapH, basicClustersEndcapH, lazyTool, false, 2);

  //----------------------------------------  
  // analyze basic clusters 
  //----------------------------------------  

  //analyzeBasicClusters(hybridClustersBarrelH, 1);
  //analyzeBasicClusters(basicClustersEndcapH, 2);

  return true;
}

void GlobeEcalClusters::analyzeSuperClusters(edm::Handle<reco::SuperClusterCollection> scH, edm::Handle<reco::BasicClusterCollection> bcH, EcalClusterLazyTools& ecalTools, bool isBarrel, int type) {

  for(unsigned int i=0; i<scH->size(); i++) {
    
    if(sc_n >= MAX_SUPERCLUSTERS) {
      std::cout << "GlobeEcalCluster: WARNING too many superclusters. " << MAX_SUPERCLUSTERS << " allowed." << std::endl;
      return;
    }
    
    reco::SuperClusterRef sc(scH, i);

    sc_bcseedind[sc_n] = -1;
    // get index to seed basic cluster
    for(unsigned int j=0; j<bcH->size(); ++j) {
      reco::BasicClusterRef basic(bcH, j);
      if (&(*sc->seed()) == &(*basic)) {
	sc_bcseedind[sc_n] = bc_n;
	fillBasicCluster(basic, ecalTools, type);
        break;
      }
    }
    
    // apply the cuts
    if(gCUT->cut(*sc))
      continue;
    
    //Old Style:
    float phi = sc->phi();
    float theta = (2*atan(exp(-sc->eta())));
    float en = sc->energy();
    float px = en*sin(theta)*cos(phi);
    float py = en*sin(theta)*sin(phi);
    float pz = en*cos(theta);
    
    new ((*sc_p4)[sc_n]) TLorentzVector();
    ((TLorentzVector *)sc_p4->At(sc_n))->SetXYZT(px, py, pz, en);
    
    new ((*sc_xyz)[sc_n]) TVector3();
    ((TVector3 *)sc_xyz->At(sc_n))->SetXYZ(sc->position().x(), sc->position().y(), sc->position().z());
    
    sc_raw[sc_n] = sc->rawEnergy();
    sc_pre[sc_n] = sc->preshowerEnergy();

    // sigma_phi and sigma_eta as defined in reco::SuperCluster
    sc_sphi[sc_n] = sc->phiWidth();
    sc_seta[sc_n] = sc->etaWidth();
    if (sc_seta[sc_n]>0) 
      sc_brem[sc_n] = sc_sphi[sc_n] / sc_seta[sc_n];
    else 
      sc_brem[sc_n]=-1; // something is not ok with sigma_eta, in this case

    // SC r9
    if (sc->rawEnergy()>0) 
      sc_r9[sc_n] = ecalTools.e3x3(*(sc->seed()));//, &(*barrelRecHits), &(*topology)) / sc->rawEnergy();
    else 
      sc_r9[sc_n] = -1;

    //SEED BC 
    sc_nbc[sc_n]= sc->clustersSize();
    
    // get indices to basic clusters
    if (sc->clustersSize() > 0) { 
      int limit = 0;
      
      for(reco::CaloCluster_iterator itClus = sc->clustersBegin(); itClus != sc->clustersEnd(); ++itClus) {
	
        if (limit >= MAX_SUPERCLUSTER_BASICCLUSTERS) {
          std::cout << "GlobeEcalCluster: WARNING too many basiclusters. in basicClustersBarrelH (" << 
            MAX_SUPERCLUSTER_BASICCLUSTERS << " allowed). Event has " << sc->clustersSize() << std::endl;
          break;
        } 
	
        for(unsigned int j=0; j<bcH->size(); ++j) {
          reco::BasicClusterRef basic(bcH, j);
          if (&(**itClus) == &(*basic)) {
	    sc_bcind[sc_n][limit] = bc_n;
	    fillBasicCluster(basic, ecalTools, isBarrel);	    
            break;
          }
        }
	
	const reco::CaloClusterPtr cc = *itClus;
	sc_bccrackcorr[sc_n][limit] = CrackCorrFunc->getValue(*cc);
	sc_bclocalcorr[sc_n][limit] = LocalCorrFunc->getValue(*cc);

        limit++;
      }
    }
    
    if (isBarrel)
      sc_barrel[sc_n] = 1;
    else
      sc_barrel[sc_n] = 0;
    sc_n++;
  }
}

void GlobeEcalClusters::fillBasicCluster(reco::BasicClusterRef bc, EcalClusterLazyTools& ecalTools, int type) {
  
  if(bc_n >= MAX_BASICCLUSTERS) {
    std::cout << "GlobeEcalCluster: WARNING too many basicclusters. " << MAX_BASICCLUSTERS << " allowed. (Missing SEED Cluster)" << std::endl;
    return;
  }
      
  float phi = bc->phi();
  float theta = (2*atan(exp(-bc->eta())));
  float en = bc->energy();
  float px = en*sin(theta)*cos(phi);
  float py = en*sin(theta)*sin(phi);
  float pz = en*cos(theta);
  
  new ((*bc_p4)[bc_n]) TLorentzVector();
  ((TLorentzVector *)bc_p4->At(bc_n))->SetXYZT(px, py, pz, en);
  
  new ((*bc_xyz)[bc_n]) TVector3();
  ((TVector3 *)bc_xyz->At(bc_n))->SetXYZ(bc->position().x(), bc->position().y(), bc->position().z());
  
  std::vector<std::pair<DetId,float > > hits = bc->hitsAndFractions();
  bc_nhits[bc_n] = hits.size(); 

  // fill the rechits associated to the basic cluster
  fillBasicClusterRecHits(hits);
  
  bc_s1[bc_n]  = ecalTools.eMax(*(bc));//, &(*endcapRecHits)); 
  bc_s4[bc_n]  = ecalTools.e2x2(*(bc));//, &(*endcapRecHits), &(*topology)); 
  bc_s9[bc_n]  = ecalTools.e3x3(*(bc));//, &(*endcapRecHits), &(*topology)); 
  bc_s25[bc_n] = ecalTools.e5x5(*(bc));//, &(*endcapRecHits), &(*topology)); 
  
  //std::vector<float> rook_vect;
  //rook_vect.push_back(ecalTools.eLeft(*(bc)));//, &(*endcapRecHits), &(*topology)));
  //rook_vect.push_back(ecalTools.eTop(*(bc)));//, &(*endcapRecHits), &(*topology)));
  //rook_vect.push_back(ecalTools.eBottom(*(bc)));//, &(*endcapRecHits), &(*topology)));
  //rook_vect.push_back(ecalTools.eRight(*(bc)));//, &(*endcapRecHits), &(*topology)));
  bc_chx[bc_n] = 0;//std::accumulate(rook_vect.begin(), rook_vect.end(), 0.);
  
  std::vector<float> vCov = ecalTools.localCovariances(*(bc));//, &(*endcapRecHits), &(*topology));
  bc_sieie[bc_n] = sqrt(vCov[0]);
  bc_sieip[bc_n] = vCov[1];
  bc_sipip[bc_n] = sqrt(vCov[2]);

  bc_type[bc_n] = type;
  bc_n++;
}

void GlobeEcalClusters::analyzeBasicClusters(edm::Handle<reco::BasicClusterCollection> bcH, EcalClusterLazyTools& ecalTools, int type) { 

  if (barrelBasicClusterColl.encode() != "") {
    
    for(unsigned int i=0; i< bcH->size(); i++) {
       
      if(bc_n >= MAX_BASICCLUSTERS) {
	std::cout << "GlobeEcalCluster: WARNING too many basicclusters. " << MAX_BASICCLUSTERS << " allowed." << std::endl;
	break;
      }
      
      reco::BasicClusterRef bc(bcH, i);
      fillBasicCluster(bc, ecalTools, type);
    }
  }
}

std::vector<float> GlobeEcalClusters::getESHits(double X, double Y, double Z, std::map<DetId, EcalRecHit> rechits_map, const CaloGeometry& geometry, CaloSubdetectorTopology *topology_p, int row) {
  std::vector<float> esHits;

  const GlobalPoint point(X,Y,Z);

  const CaloSubdetectorGeometry *geometry_p ;
  geometry_p = geometry.getSubdetectorGeometry (DetId::Ecal,EcalPreshower) ;

  DetId esId1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 1);
  DetId esId2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 2);
  ESDetId esDetId1 = (esId1 == DetId(0)) ? ESDetId(0) : ESDetId(esId1);
  ESDetId esDetId2 = (esId2 == DetId(0)) ? ESDetId(0) : ESDetId(esId2);

  std::map<DetId, EcalRecHit>::iterator it;
  ESDetId next;
  ESDetId strip1;
  ESDetId strip2;

  strip1 = esDetId1;
  strip2 = esDetId2;
    
  EcalPreshowerNavigator theESNav1(strip1, topology_p);
  theESNav1.setHome(strip1);
    
  EcalPreshowerNavigator theESNav2(strip2, topology_p);
  theESNav2.setHome(strip2);

  if (row == 1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.north();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.east();
  } else if (row == -1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.south();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.west();
  }

  // Plane 1 
  if (strip1 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip1);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
    else esHits.push_back(0);
    //cout<<"center : "<<strip1<<" "<<it->second.energy()<<endl;      

    // east road 
    for (int i=0; i<15; ++i) {
      next = theESNav1.east();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"east "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"east "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // west road 
    theESNav1.setHome(strip1);
    theESNav1.home();
    for (int i=0; i<15; ++i) {
      next = theESNav1.west();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"west "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"west "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  if (strip2 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip2);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
    else esHits.push_back(0);
    //cout<<"center : "<<strip2<<" "<<it->second.energy()<<endl;      

    // north road 
    for (int i=0; i<15; ++i) {
      next = theESNav2.north();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"north "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;  
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"north "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // south road 
    theESNav2.setHome(strip2);
    theESNav2.home();
    for (int i=0; i<15; ++i) {
      next = theESNav2.south();
      if (next != ESDetId(0)) {
        it = rechits_map.find(next);
        if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
        else esHits.push_back(0);
        //cout<<"south "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
        //cout<<"south "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  return esHits;
}

std::vector<float> GlobeEcalClusters::getESShape(std::vector<float> ESHits0)
{
  std::vector<float> esShape;
    
  const int nBIN = 21;
  float esRH_F[nBIN];
  float esRH_R[nBIN];
  for (int idx=0; idx<nBIN; idx++) {
    esRH_F[idx] = 0.;
    esRH_R[idx] = 0.;
  }

  for(int ibin=0; ibin<((nBIN+1)/2); ibin++) {
    if (ibin==0) {
      esRH_F[(nBIN-1)/2] = ESHits0[ibin];
      esRH_R[(nBIN-1)/2] = ESHits0[ibin+31];
    } else {
      esRH_F[(nBIN-1)/2+ibin] = ESHits0[ibin];
      esRH_F[(nBIN-1)/2-ibin] = ESHits0[ibin+15];
      esRH_R[(nBIN-1)/2+ibin] = ESHits0[ibin+31];
      esRH_R[(nBIN-1)/2-ibin] = ESHits0[ibin+31+15];
    }
  } 

  // ---- Effective Energy Deposit Width ---- //
  double EffWidthSigmaXX = 0.;
  double EffWidthSigmaYY = 0.;
  double totalEnergyXX   = 0.;
  double totalEnergyYY   = 0.;
  double EffStatsXX      = 0.;
  double EffStatsYY      = 0.;
  for (int id_X=0; id_X<21; id_X++) {
    totalEnergyXX  += esRH_F[id_X];
    EffStatsXX     += esRH_F[id_X]*(id_X-10)*(id_X-10);
    totalEnergyYY  += esRH_R[id_X];
    EffStatsYY     += esRH_R[id_X]*(id_X-10)*(id_X-10);
  }
  EffWidthSigmaXX  = (totalEnergyXX>0.)  ? sqrt(fabs(EffStatsXX  / totalEnergyXX))   : 0.;
  EffWidthSigmaYY  = (totalEnergyYY>0.)  ? sqrt(fabs(EffStatsYY  / totalEnergyYY))   : 0.;

  esShape.push_back(EffWidthSigmaXX);
  esShape.push_back(EffWidthSigmaYY);
    
  return esShape;
}

//----------------------------------------------------------------------

void 
GlobeEcalClusters::fillBasicClusterRecHits(const std::vector<std::pair<DetId,float > > &hits)
{
  bc_hitdetid->push_back(vector<Int_t>());
  for (unsigned bcindex = 0; bcindex < hits.size(); ++bcindex)
    bc_hitdetid->back().push_back(hits[bcindex].first.rawId());
}

//----------------------------------------------------------------------
