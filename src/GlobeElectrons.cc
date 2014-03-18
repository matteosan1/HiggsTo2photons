#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeElectrons.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Tools.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/EcalBarrelTopology.h"
#include "Geometry/CaloTopology/interface/EcalEndcapTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"

#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollection.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"x

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//#include "HiggsAnalysis/HiggsTo2photons/interface/PFIsolation.h"
//#include "HiggsAnalysis/HiggsTo2photons/interface/Mustache.h"

#include "Utilities/General/interface/FileInPath.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include <cstdlib>
#include <iostream>


GlobeElectrons::GlobeElectrons(const edm::ParameterSet& iConfig) {
  
  GlobeBase::GlobeBase(iConfig);

  electronColl = iConfig.getParameter<edm::InputTag>("ElectronColl");
  doAodSim = iConfig.getParameter<bool>("doAodSim");

  conversionColl = iConfig.getParameter<edm::InputTag>("ConvertedPhotonColl");
  beamSpotColl = iConfig.getParameter<edm::InputTag>("BeamSpot");
  caloTowerColl = iConfig.getParameter<edm::InputTag>("CaloTowerColl");
  trackColl = iConfig.getParameter<edm::InputTag>("TrackColl");
  trackColl2 = iConfig.getParameter<edm::InputTag>("TrackColl3");
  vertexColl = iConfig.getParameter<edm::InputTag>("VertexColl_std");
  rhoCollection = iConfig.getParameter<edm::InputTag>("rhoCollection_algo1");

  hybridSuperClusterColl = iConfig.getParameter<edm::InputTag>("HybridSuperClusterColl");
  endcapSuperClusterColl = iConfig.getParameter<edm::InputTag>("EndcapSuperClusterColl");
  ecalHitEBColl = iConfig.getParameter<edm::InputTag>("EcalHitEBColl");
  ecalHitEEColl = iConfig.getParameter<edm::InputTag>("EcalHitEEColl");
  ecalHitESColl = iConfig.getParameter<edm::InputTag>("EcalHitESColl");
  hcalHitColl = iConfig.getParameter<edm::InputTag>("HcalHitsBEColl");

  eleRegressionFilename = iConfig.getParameter<std::string> ("eleRegressionFileName");  
  eleRegressionType     = iConfig.getParameter<int> ("eleRegressionType");  

  eleRegression = new ElectronEnergyRegressionEvaluate();
  char filename[500], filename2[500];
  char* descr = getenv("CMSSW_BASE");
  const char* version[] = {"V1", "V2"};
  sprintf(filename, "%s/src/EgammaAnalysis/ElectronTools/data/%s_%s.root", descr, eleRegressionFilename.c_str(), version[eleRegressionType]);
  sprintf(filename2, "EgammaAnalysis/ElectronTools/data/%s_%s.root", eleRegressionFilename.c_str(), version[eleRegressionType]);
  if (fexist(filename)) {
    sprintf(filename, "http://cern.ch/sani/%s_%s.root", eleRegressionFilename.c_str(), version[eleRegressionType]);
    eleRegression->initialize(filename, ElectronEnergyRegressionEvaluate::ElectronEnergyRegressionType::kNoTrkVar);
  } else {
    eleRegression->initialize(filename2, ElectronEnergyRegressionEvaluate::ElectronEnergyRegressionType::kNoTrkVar);
  }

  mvaNonTrigWeightFiles = iConfig.getParameter<std::vector<std::string> >("electronNonTrigMVAWeightFileNames");
  mvaTrigWeightFiles    = iConfig.getParameter<std::vector<std::string> >("electronTrigMVAWeightFileNames");
  
  for(unsigned int j=0; j<mvaTrigWeightFiles.size(); j++)
    myManualCatWeightsTrig.push_back(edm::FileInPath(mvaTrigWeightFiles[j]).fullPath());
  
  for(unsigned int j=0; j<mvaNonTrigWeightFiles.size(); j++)
    myManualCatWeightsNonTrig.push_back(edm::FileInPath(mvaNonTrigWeightFiles[j]).fullPath());
  
  myMVANonTrig = new EGammaMvaEleEstimator();
  myMVANonTrig->initialize("BDT",
           EGammaMvaEleEstimator::kNonTrig,
           true, // use manual cat
           myManualCatWeightsNonTrig);
  
  myMVATrig = new EGammaMvaEleEstimator();
  myMVATrig->initialize("BDT",
           EGammaMvaEleEstimator::kTrig,
           true, // use manual cat
           myManualCatWeightsTrig);

  energyCorrectionsFromDB = iConfig.getParameter<bool> ("energyCorrectionsFromDB"); 
  energyRegFilename       = iConfig.getParameter<std::string> ("energyCorrectionsFileNameEle");  
  regressionVersion       = iConfig.getParameter<std::string> ("energyCorrectionsVersion");
}

GlobeElectrons::~GlobeElectrons() {

  delete myMVANonTrig;
  delete myMVATrig;
  GlobeBase::~GlobeBase();
}

void GlobeElectrons::defineBranch(GlobeAnalyzer* ana) {
  
  GlobeBase::defineBranch(ana);
  
  el_sc = new TClonesArray("TLorentzVector", MAX_ELECTRONS);
  el_p4 = new TClonesArray("TLorentzVector", MAX_ELECTRONS);
  el_p4_corr = new TClonesArray("TLorentzVector", MAX_ELECTRONS);
  el_momvtx = new TClonesArray("TVector3", MAX_ELECTRONS);  
  el_momvtxconst = new TClonesArray("TVector3", MAX_ELECTRONS);
  el_momcalo = new TClonesArray("TVector3", MAX_ELECTRONS);
  el_momout = new TClonesArray("TVector3", MAX_ELECTRONS);
  el_posvtx = new TClonesArray("TVector3", MAX_ELECTRONS);
  el_poscalo = new TClonesArray("TVector3", MAX_ELECTRONS);

  el_catbased = new std::vector<std::vector<int> >;
  el_schits = new std::vector<std::vector<UInt_t> >;
  el_bchits = new std::vector<std::vector<UInt_t> >;

  char a1[50], a2[50];
    sprintf(a1, "el%s_n", prefix_.c_str());
  sprintf(a2, "el%s_n/I", prefix_.c_str());
  ana->Branch(a1, &el_n, a2);
  
  sprintf(a1, "el%s_sc", prefix_.c_str());
  ana->Branch(a1, "TClonesArray", &el_sc, 32000, 0);
  
  sprintf(a1, "el%s_p4", prefix_.c_str());
  ana->Branch(a1, "TClonesArray", &el_p4, 32000, 0);

  sprintf(a1, "el%s_p4_corr", prefix_.c_str());
  ana->Branch(a1, "TClonesArray", &el_p4_corr, 32000, 0);

  sprintf(a1, "el%s_momvtx", prefix_.c_str());
  ana->Branch(a1, "TClonesArray", &el_momvtx, 32000, 0);

  sprintf(a1, "el%s_momvtxconst", prefix_.c_str());
  ana->Branch(a1, "TClonesArray", &el_momvtxconst, 32000, 0);

  sprintf(a1, "el%s_momcalo", prefix_.c_str());
  ana->Branch(a1, "TClonesArray", &el_momcalo, 32000, 0);

  sprintf(a1, "el%s_momout", prefix_.c_str());
  ana->Branch(a1, "TClonesArray", &el_momout, 32000, 0);
  
  sprintf(a1, "el%s_posvtx", prefix_.c_str());
  ana->Branch(a1, "TClonesArray", &el_posvtx, 32000, 0);

  sprintf(a1, "el%s_poscalo", prefix_.c_str());
  ana->Branch(a1, "TClonesArray", &el_poscalo, 32000, 0);
  
  sprintf(a1, "el%s_eopin", prefix_.c_str());
  sprintf(a2, "el%s_eopin[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_eopin, a2);
  
  sprintf(a1, "el%s_eseedopout", prefix_.c_str());
  sprintf(a2, "el%s_eseedopout[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_eseedopout, a2);
  
  sprintf(a1, "el%s_pout", prefix_.c_str());
  sprintf(a2, "el%s_pout[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_pout, a2);
  
  sprintf(a1, "el%s_pin", prefix_.c_str());
  sprintf(a2, "el%s_pin[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_pin, a2);
  
  sprintf(a1, "el%s_e1x5", prefix_.c_str());
  sprintf(a2, "el%s_e1x5[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_e1x5, a2);
  
  sprintf(a1, "el%s_e5x5", prefix_.c_str());
  sprintf(a2, "el%s_e5x5[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_e5x5, a2);

  sprintf(a1, "el%s_e2x5", prefix_.c_str());
  sprintf(a2, "el%s_e2x5[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_e2x5, a2);
 
  sprintf(a1, "el%s_sipip", prefix_.c_str());
  sprintf(a2, "el%s_sipip[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_sipip, a2);
  
  sprintf(a1, "el%s_sieie", prefix_.c_str());
  sprintf(a2, "el%s_sieie[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_sieie, a2);

  sprintf(a1, "el%s_sieiesc", prefix_.c_str());
  sprintf(a2, "el%s_sieiesc[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_sieiesc, a2);

  sprintf(a1, "el%s_eseffsixix", prefix_.c_str());
  sprintf(a2, "el%s_eseffsixix[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_eseffsixix, a2);

  sprintf(a1, "el%s_eseffsiyiy", prefix_.c_str());
  sprintf(a2, "el%s_eseffsiyiy[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_eseffsiyiy, a2);

  sprintf(a1, "el%s_eseedopin", prefix_.c_str());
  sprintf(a2, "el%s_eseedopin[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_eseedopin, a2);
  
  sprintf(a1, "el%s_fbrem", prefix_.c_str());
  sprintf(a2, "el%s_fbrem[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_fbrem, a2);

  sprintf(a1, "el%s_nbrem", prefix_.c_str());
  sprintf(a2, "el%s_nbrem[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_nbrem, a2);

  sprintf(a1, "el%s_hoe", prefix_.c_str());
  sprintf(a2, "el%s_hoe[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_hoe, a2);

  sprintf(a1, "el%s_hoed1", prefix_.c_str());
  sprintf(a2, "el%s_hoed1[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_hoed1, a2);

  sprintf(a1, "el%s_hoed2", prefix_.c_str());
  sprintf(a2, "el%s_hoed2[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_hoed2, a2);

  sprintf(a1, "el%s_hoebc", prefix_.c_str());
  sprintf(a2, "el%s_hoebc[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_hoebc, a2);

  sprintf(a1, "el%s_hoebcd1", prefix_.c_str());
  sprintf(a2, "el%s_hoebcd1[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_hoebcd1, a2);

  sprintf(a1, "el%s_hoebcd2", prefix_.c_str());
  sprintf(a2, "el%s_hoebcd2[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_hoebcd2, a2);
  
  sprintf(a1, "el%s_detain", prefix_.c_str());
  sprintf(a2, "el%s_detain[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_detain, a2);
  
  sprintf(a1, "el%s_dphiin", prefix_.c_str());
  sprintf(a2, "el%s_dphiin[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_dphiin, a2);
  
  sprintf(a1, "el%s_detaout", prefix_.c_str());
  sprintf(a2, "el%s_detaout[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_detaout, a2);
  
  sprintf(a1, "el%s_dphiout", prefix_.c_str());
  sprintf(a2, "el%s_dphiout[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_dphiout, a2);
   
  sprintf(a1, "el%s_class", prefix_.c_str());
  sprintf(a2, "el%s_class[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_class, a2);
 
  sprintf(a1, "el%s_crack", prefix_.c_str());
  sprintf(a2, "el%s_crack[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_crack, a2);
   
  sprintf(a1, "el%s_nambtk", prefix_.c_str());
  sprintf(a2, "el%s_nambtk[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_nambtk, a2);

  sprintf(a1, "el%s_scind", prefix_.c_str());
  sprintf(a2, "el%s_scind[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_scind, a2);
  
  sprintf(a1, "el%s_z0", prefix_.c_str());
  sprintf(a2, "el%s_z0[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_z0, a2);
  
  sprintf(a1, "el%s_d0", prefix_.c_str());
  sprintf(a2, "el%s_d0[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_d0, a2);

  sprintf(a1, "el%s_chi2", prefix_.c_str());
  sprintf(a2, "el%s_chi2[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_chi2, a2);
  
  sprintf(a1, "el%s_mva_nontrig", prefix_.c_str());
  sprintf(a2, "el%s_mva_nontrig[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_mva_nontrig, a2);
  
  sprintf(a1, "el%s_mva_trig", prefix_.c_str());
  sprintf(a2, "el%s_mva_trig[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_mva_trig, a2);
  
  sprintf(a1, "el%s_ch_gsf", prefix_.c_str());
  sprintf(a2, "el%s_ch_gsf[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_ch_gsf, a2);

  sprintf(a1, "el%s_ch_scpix", prefix_.c_str());
  sprintf(a2, "el%s_ch_scpix[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_ch_scpix, a2);

  sprintf(a1, "el%s_charge", prefix_.c_str());
  sprintf(a2, "el%s_charge[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_charge, a2);
  
  sprintf(a1, "el%s_losthits", prefix_.c_str());
  sprintf(a2, "el%s_losthits[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_losthits, a2);
  
  sprintf(a1, "el%s_validhits", prefix_.c_str());
  sprintf(a2, "el%s_validhits[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_validhits, a2);

  sprintf(a1, "el%s_hp_expin", prefix_.c_str());
  sprintf(a2, "el%s_hp_expin[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_hp_expin, a2);

  sprintf(a1, "el%s_hp_expout", prefix_.c_str());
  sprintf(a2, "el%s_hp_expout[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_hp_expout, a2);

  sprintf(a1, "el%s_catbased", prefix_.c_str());
  ana->Branch(a1, "std::vector<std::vector<int> >", &el_catbased);

  sprintf(a1, "el%s_schits", prefix_.c_str());
  ana->Branch(a1, "std::vector<std::vector<UInt_t> >", &el_schits);

  sprintf(a1, "el%s_bchits", prefix_.c_str());
  ana->Branch(a1, "std::vector<std::vector<UInt_t> >", &el_bchits);

  sprintf(a1, "el%s_tkind", prefix_.c_str());
  sprintf(a2, "el%s_tkind[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_tkind, a2); 

  sprintf(a1, "el%s_pfiso_neutral", prefix_.c_str());
  sprintf(a2, "el%s_pfiso_neutral[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_pfiso_neutral, a2);

  sprintf(a1, "el%s_pfiso_charged", prefix_.c_str());
  sprintf(a2, "el%s_pfiso_charged[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_pfiso_charged, a2);
  
  sprintf(a1, "el%s_pfiso_photon", prefix_.c_str());
  sprintf(a2, "el%s_pfiso_photon[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_pfiso_photon, a2);

  sprintf(a1, "el%s_hcaliso03", prefix_.c_str());
  sprintf(a2, "el%s_hcaliso03[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_hcaliso03, a2);

  //sprintf(a1, "el%s_hcalsolidiso03", prefix_.c_str());
  //sprintf(a2, "el%s_hcalsolidiso03[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  //ana->Branch(a1, &el_hcalsolidiso03, a2);

  sprintf(a1, "el%s_ecaliso03", prefix_.c_str());
  sprintf(a2, "el%s_ecaliso03[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_ecaliso03, a2);
  
  sprintf(a1, "el%s_tkiso03", prefix_.c_str());
  sprintf(a2, "el%s_tkiso03[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_tkiso03, a2);

  sprintf(a1, "el%s_hcaliso04", prefix_.c_str());
  sprintf(a2, "el%s_hcaliso04[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_hcaliso04, a2);

  //sprintf(a1, "el%s_hcalsolidiso04", prefix_.c_str());
  //sprintf(a2, "el%s_hcalsolidiso04[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  //ana->Branch(a1, &el_hcalsolidiso04, a2);

  sprintf(a1, "el%s_hcalbciso03", prefix_.c_str());
  sprintf(a2, "el%s_hcalbciso03[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_hcalbciso03, a2);

  sprintf(a1, "el%s_hcalbciso04", prefix_.c_str());
  sprintf(a2, "el%s_hcalbciso04[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_hcalbciso04, a2);

  sprintf(a1, "el%s_ecaliso04", prefix_.c_str());
  sprintf(a2, "el%s_ecaliso04[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_ecaliso04, a2);
  
  sprintf(a1, "el%s_tkiso04", prefix_.c_str());
  sprintf(a2, "el%s_tkiso04[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_tkiso04, a2);
  
  sprintf(a1, "el%s_tkdrv", prefix_.c_str());
  sprintf(a2, "el%s_tkdrv[el%s_n]/O", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_tkdrv, a2);
 
  sprintf(a1, "el%s_ecaldrv", prefix_.c_str());
  sprintf(a2, "el%s_ecaldrv[el%s_n]/O", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_ecaldrv, a2);

  sprintf(a1, "el%s_ip_ctf", prefix_.c_str());
  sprintf(a2, "el%s_ip_ctf[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_ip_ctf, a2);

  sprintf(a1, "el%s_ip_gsf", prefix_.c_str());
  sprintf(a2, "el%s_ip_gsf[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_ip_gsf, a2);

  sprintf(a1, "el%s_dist", prefix_.c_str());
  sprintf(a2, "el%s_dist[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_dist, a2);

  sprintf(a1, "el%s_dcot", prefix_.c_str());
  sprintf(a2, "el%s_dcot[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_dcot, a2);

  sprintf(a1, "el%s_hp_1pxb", prefix_.c_str());
  sprintf(a2, "el%s_hp_1pxb[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_1pxb, a2);

  sprintf(a1, "el%s_hp_1pxf", prefix_.c_str());
  sprintf(a2, "el%s_hp_1pxf[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_1pxf, a2);

  sprintf(a1, "el%s_conv", prefix_.c_str());
  sprintf(a2, "el%s_conv[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_conv, a2);

  sprintf(a1, "el%s_corr_energy", prefix_.c_str());
  sprintf(a2, "el%s_corr_energy[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_corr_energy, a2);

  sprintf(a1, "el%s_corr_energyerr", prefix_.c_str());
  sprintf(a2, "el%s_corr_energyerr[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_corr_energyerr, a2);

  sprintf(a1, "el%s_calib_energy", prefix_.c_str());
  sprintf(a2, "el%s_calib_energy[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_calib_energy, a2);

  sprintf(a1, "el%s_calib_energyerr", prefix_.c_str());
  sprintf(a2, "el%s_calib_energyerr[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_calib_energyerr, a2);

  sprintf(a1, "el%s_regr_energy", prefix_.c_str());
  sprintf(a2, "el%s_regr_energy[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_regr_energy, a2);

  sprintf(a1, "el%s_regr_energyerr", prefix_.c_str());
  sprintf(a2, "el%s_regr_energyerr[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_regr_energyerr, a2);

  sprintf(a1, "el%s_eleopout", prefix_.c_str());
  sprintf(a2, "el%s_eleopout[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_eleopout, a2);

  sprintf(a1, "el%s_detaeleout", prefix_.c_str());
  sprintf(a2, "el%s_detaeleout[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_detaeleout, a2);

  sprintf(a1, "el%s_kfhits", prefix_.c_str());
  sprintf(a2, "el%s_kfhits[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_kfhits, a2);

  sprintf(a1, "el%s_kfchi2", prefix_.c_str());
  sprintf(a2, "el%s_kfchi2[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_kfchi2, a2);

  sprintf(a1, "el%s_psenergy", prefix_.c_str());
  sprintf(a2, "el%s_psenergy[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_psenergy, a2);  

  sprintf(a1, "el%s_passmvapresel", prefix_.c_str());
  sprintf(a2, "el%s_passmvapresel[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_passmvapresel, a2);

  sprintf(a1, "el%s_passcutpresel", prefix_.c_str());
  sprintf(a2, "el%s_passcutpresel[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_passcutpresel, a2);

  sprintf(a1, "el%s_psenergypf", prefix_.c_str());
  sprintf(a2, "el%s_psenergypf[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_psenergypf, a2);  
  
  sprintf(a1, "el%s_nbrempf", prefix_.c_str());
  sprintf(a2, "el%s_nbrempf[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_nbrempf, a2);

  sprintf(a1, "el%s_eseedpf", prefix_.c_str());
  sprintf(a2, "el%s_eseedpf[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_eseedpf, a2);

  sprintf(a1, "el%s_epf", prefix_.c_str());
  sprintf(a2, "el%s_epf[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_epf, a2);

  sprintf(a1, "el%s_psly1", prefix_.c_str());
  sprintf(a2, "el%s_psly1[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_psly1, a2);

  sprintf(a1, "el%s_psly2", prefix_.c_str());
  sprintf(a2, "el%s_psly2[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_psly2, a2);

  sprintf(a1, "el%s_psnstriply1", prefix_.c_str());
  sprintf(a2, "el%s_psnstriply1[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_psnstriply1, a2);

  sprintf(a1, "el%s_psnstriply2", prefix_.c_str());
  sprintf(a2, "el%s_psnstriply2[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_psnstriply2, a2);

  sprintf(a1, "el%s_D0Vtx", prefix_.c_str());
  sprintf(a2, "el%s_D0Vtx[el%s_n][100]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_D0Vtx, a2);

  sprintf(a1, "el%s_DZVtx", prefix_.c_str());
  sprintf(a2, "el%s_DZVtx[el%s_n][100]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_DZVtx, a2);

  sprintf(a1, "el%s_must", prefix_.c_str());
  sprintf(a2, "el%s_must[el%s_n]/F", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_must, a2);

  sprintf(a1, "el%s_mustnc", prefix_.c_str());
  sprintf(a2, "el%s_mustnc[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_mustnc, a2);

  sprintf(a1, "el%s_r9", prefix_.c_str());
  sprintf(a2, "el%s_r9[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_r9, a2);

  sprintf(a1, "el%s_gsfchi2", prefix_.c_str());
  sprintf(a2, "el%s_gsfchi2[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_gsfchi2, a2);
  
  sprintf(a1, "el%s_ip3d", prefix_.c_str());
  sprintf(a2, "el%s_ip3d[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_ip3d, a2);

  sprintf(a1, "el%s_ip3d_err", prefix_.c_str());
  sprintf(a2, "el%s_ip3d_err[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_ip3d_err, a2);

  sprintf(a1, "el%s_ip3d_sig", prefix_.c_str());
  sprintf(a2, "el%s_ip3d_sig[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_ip3d_sig, a2);

  sprintf(a1, "el%s_sc_time", prefix_.c_str());
  sprintf(a2, "el%s_sc_time[el%s_n]/I", prefix_.c_str(), prefix_.c_str());
  ana->Branch(a1, &el_sc_time, a2);

  sprintf(a1, "el%s_conv_vtxProb", prefix_.c_str());
  sprintf(a2, "el%s_conv_vtxProb[el%s_n]/F", prefix_.c_str(), prefix_.c_str());

  ana->Branch(a1, &el_conv_vtxProb, a2);
}

bool GlobeElectrons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // take collections
  edm::Handle<reco::GsfElectronCollection> elH;
  iEvent.getByLabel(electronColl, elH);

  edm::Handle<reco::GsfElectronCollection> calibEleH;
  iEvent.getByLabel("calibratedElectrons", "calibratedGsfElectrons", calibEleH);

  edm::Handle<edm::ValueMap<double>> calibEnergyH;
  iEvent.getByLabel("calibratedElectrons" ,"eneRegForGsfEle", calibEnergyH);
  const edm::ValueMap<double>* calibEnergy = calibEnergyH.product();
  
  edm::Handle<edm::ValueMap<double>> calibEnergyErrH;
  iEvent.getByLabel("calibratedElectrons" ,"eneErrorRegForGsfEle", calibEnergyErrH);
  const edm::ValueMap<double>* calibEnergyErr = calibEnergyErrH.product();
  
  edm::Handle<edm::ValueMap<double>> corrEnergyH;
  iEvent.getByLabel("eleRegressionEnergy" ,"eneRegForGsfEle", corrEnergyH);
  const edm::ValueMap<double>* corrEnergy = corrEnergyH.product();
  
  edm::Handle<edm::ValueMap<double>> corrEnergyErrH;
  iEvent.getByLabel("eleRegressionEnergy" ,"eneErrorRegForGsfEle", corrEnergyErrH);
  const edm::ValueMap<double>* corrEnergyErr = corrEnergyErrH.product();

  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel(conversionColl, hConversions);

  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel(beamSpotColl, bsHandle);
  const reco::BeamSpot &thebs = *bsHandle.product();

  edm::Handle<edm::SortedCollection<CaloTower> > ctH;
  iEvent.getByLabel(caloTowerColl, ctH);

  edm::Handle<reco::TrackCollection> tkH;
  iEvent.getByLabel(trackColl, tkH);

  edm::Handle<reco::GsfTrackCollection> tkH2;
  iEvent.getByLabel(trackColl2, tkH2);
  
  edm::Handle<reco::SuperClusterCollection> superClustersBarrelH; 
  iEvent.getByLabel(hybridSuperClusterColl,superClustersBarrelH);
  
  edm::Handle<reco::SuperClusterCollection> superClustersEndcapH; 
  iEvent.getByLabel(endcapSuperClusterColl, superClustersEndcapH);

  edm::Handle<reco::VertexCollection> vtxH;
  iEvent.getByLabel(vertexColl, vtxH);

  //edm::Handle<double> rhoHandle;
  //iEvent.getByLabel(rhoCollection, rhoHandle);
  //double rho = *(rhoHandle.product());
  
  // transient track builder needed for ele ID MVA
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", trackBuilder_);  
  const TransientTrackBuilder thebuilder = *(trackBuilder_.product());
  
  /*
  edm::Handle<reco::PFCandidateCollection> pfHandle;
  iEvent.getByLabel(pfColl, pfHandle);

  edm::Handle<reco::PFCandidateCollection> pfHandlePu;
  iEvent.getByLabel("pfPileUp", pfHandlePu);
  */

  edm::ESHandle<TransientTrackBuilder> hTransientTrackBuilder;
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", hTransientTrackBuilder);
  transientTrackBuilder = hTransientTrackBuilder.product();

  //IsoDepositVals electronIsoVals(3);
  //for (size_t j = 0; j<inputTagIsoValElectronsPFId_.size(); ++j) {
  //  iEvent.getByLabel(inputTagIsoValElectronsPFId_[j], electronIsoVals[j]);
  //}

  edm::ESHandle<CaloTopology> theCaloTopo;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  const CaloTopology *topology = theCaloTopo.product();

  edm::Handle<EBRecHitCollection> pEBRecHitH;
  edm::Handle<EERecHitCollection> pEERecHitH;
  iEvent.getByLabel(ecalHitEBColl, pEBRecHitH);
  iEvent.getByLabel(ecalHitEEColl, pEERecHitH);
  const EcalRecHitCollection *barrelRecHits = pEBRecHitH.product();
  const EcalRecHitCollection *endcapRecHits = pEERecHitH.product();
  
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  const CaloGeometry& geometry = *geoHandle;
  const CaloSubdetectorGeometry *geometryES = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
  CaloSubdetectorTopology *topology_p = 0;
  if (geometryES) 
    topology_p = new EcalPreshowerTopology(geoHandle);

  EcalClusterLazyTools ecalLazyTool(iEvent, iSetup, ecalHitEBColl, ecalHitEEColl);
  edm::Handle<EcalRecHitCollection> ESRecHits;
  iEvent.getByLabel(ecalHitESColl , ESRecHits);

  rechits_map_.clear();
  if (ESRecHits.isValid()) {
    EcalRecHitCollection::const_iterator it;
    for (it = ESRecHits->begin(); it != ESRecHits->end(); ++it) {
      // remove bad ES rechits
      if (it->recoFlag()==1 || it->recoFlag()==14 || (it->recoFlag()<=10 && it->recoFlag()>=5)) continue;
      //Make the map of DetID, EcalRecHit pairs
      rechits_map_.insert(std::make_pair(it->id(), *it));
    }
  }
  
  el_p4_corr->Clear();
  el_sc->Clear();
  el_p4->Clear(); 
  el_momvtx->Clear();
  el_momvtxconst->Clear();
  el_momcalo->Clear();
  el_posvtx->Clear();
  el_poscalo->Clear(); 
  el_catbased ->clear();
  el_schits->clear();
  el_bchits->clear();
  el_n = 0;

  if (debug_level > 9)
    std::cout << "GlobeElectrons: Electron collection size: "<< elH->size() << std::endl;

  for(reco::GsfElectronCollection::const_iterator igsf = elH->begin(); igsf != elH->end(); igsf++) {
    
    if (el_n >= MAX_ELECTRONS) {
      std::cout << "GlobeElectrons: WARNING TOO MANY ELECTRONS: " << elH->size() << " (allowed " << MAX_ELECTRONS << ")" << std::endl;
      break;
    }
    
    reco::GsfElectron egsf = reco::GsfElectron(*igsf);
    
    if(gCUT->cut(egsf))
      continue;
    
    float phi = egsf.superCluster()->phi();
    float theta = (2*atan(exp(-egsf.superCluster()->eta())));
    float en = egsf.ecalEnergy();
    float px = en*sin(theta)*cos(phi);
    float py = en*sin(theta)*sin(phi);
    float pz = en*cos(theta);

    new ((*el_sc)[el_n]) TLorentzVector();
    ((TLorentzVector *)el_sc->At(el_n))->SetXYZT(px, py, pz, en);
    
    new ((*el_p4)[el_n]) TLorentzVector();
    ((TLorentzVector *)el_p4->At(el_n))->SetXYZT(egsf.px(), egsf.py(), egsf.pz(), egsf.energy());

    new ((*el_momvtxconst)[el_n]) TVector3();
    ((TVector3 *)el_momvtxconst->At(el_n))->SetXYZ(egsf.trackMomentumAtVtxWithConstraint().x(), 
                                              egsf.trackMomentumAtVtxWithConstraint().y(), egsf.trackMomentumAtVtxWithConstraint().z());
    
    new ((*el_momvtx)[el_n]) TVector3();
    ((TVector3 *)el_momvtx->At(el_n))->SetXYZ(egsf.trackMomentumAtVtx().x(), 
                                              egsf.trackMomentumAtVtx().y(), egsf.trackMomentumAtVtx().z());
    
    new ((*el_momcalo)[el_n]) TVector3();
    ((TVector3 *)el_momcalo->At(el_n))->SetXYZ(egsf.trackMomentumAtCalo().x(), 
                                               egsf.trackMomentumAtCalo().y(), egsf.trackMomentumAtCalo().z());

    new ((*el_momout)[el_n]) TVector3();
    ((TVector3 *)el_momout->At(el_n))->SetXYZ(egsf.trackMomentumOut().x(), 
                                               egsf.trackMomentumOut().y(), egsf.trackMomentumOut().z());
    
    new ((*el_posvtx)[el_n]) TVector3();
    ((TVector3 *)el_posvtx->At(el_n))->SetXYZ(egsf.trackPositionAtVtx().x(), 
                                              egsf.trackPositionAtVtx().y(), egsf.trackPositionAtVtx().z());

    new ((*el_poscalo)[el_n]) TVector3();
    ((TVector3 *)el_poscalo->At(el_n))->SetXYZ(egsf.trackPositionAtCalo().x(), 
                                               egsf.trackPositionAtCalo().y(), egsf.trackPositionAtCalo().z());
    
    el_pout[el_n] = egsf.trackMomentumOut().R();
    el_pin[el_n] = egsf.trackMomentumAtVtx().R();
    
    el_1pxb[el_n] = 0;
    el_1pxf[el_n] = 0;

    if (egsf.gsfTrack()->hitPattern().hasValidHitInFirstPixelBarrel())
      el_1pxb[el_n] = 1;

    if (egsf.gsfTrack()->hitPattern().hasValidHitInFirstPixelEndcap())
      el_1pxf[el_n] = 1;

    el_fbrem[el_n] = egsf.fbrem();
          
    el_conv_vtxProb[el_n] = 0;
    for (reco::ConversionCollection::const_iterator conv = hConversions->begin(); conv!= hConversions->end(); ++conv) {
      reco::Vertex vtx = conv->conversionVertex();
      if (vtx.isValid()) {
	if (ConversionTools::matchesConversion(egsf, *conv)) {
	  el_conv_vtxProb[el_n] = TMath::Prob( vtx.chi2(), vtx.ndof() );
	  break;
	}
      }
    }
    
    bool passconversionveto = !ConversionTools::hasMatchedConversion(egsf, hConversions, thebs.position());
    el_conv[el_n] = int(passconversionveto);

    el_eseedopout[el_n] = egsf.eSeedClusterOverPout();
    el_eseedopin[el_n] = egsf.eSeedClusterOverP();
    el_eopin[el_n] = egsf.eSuperClusterOverP();
    el_eleopout[el_n] = egsf.eEleClusterOverPout();
    el_detain[el_n] = egsf.deltaEtaSuperClusterTrackAtVtx();
    el_dphiin[el_n] = egsf.deltaPhiSuperClusterTrackAtVtx();
    el_detaout[el_n] = egsf.deltaEtaSeedClusterTrackAtCalo();
    el_dphiout[el_n] = egsf.deltaPhiSeedClusterTrackAtCalo();
    el_detaeleout[el_n] = egsf.deltaEtaEleClusterTrackAtCalo();
    el_nambtk[el_n] = egsf.ambiguousGsfTracksSize();
    el_class[el_n] = egsf.classification();
    el_nbrem[el_n] = egsf.numberOfBrems();
    el_e5x5[el_n] = egsf.e5x5();
    el_e2x5[el_n] = egsf.e2x5Max();
    el_e1x5[el_n] = egsf.e1x5();
    el_sieie[el_n] = egsf.sigmaIetaIeta();
    //el_1oe_1op[el_n] = 1./egsf.ecalEnergy() - 1./egsf.trackMomentumAtVtx().R();
    el_1oe_1op[el_n] = 1./egsf.ecalEnergy() - 1./egsf.p();
    el_r9[el_n] = ecalLazyTool.e3x3(*(egsf.superCluster()->seed())) / egsf.superCluster()->rawEnergy();

    el_sc_time[el_n] = ecalLazyTool.SuperClusterTime(*(egsf.superCluster()), iEvent);

    el_must[el_n] = -9999.;
    el_mustnc[el_n] = -1;
    reco::Mustache m;
    m.MustacheID(*(egsf.superCluster()), el_mustnc[el_n], el_must[el_n]);

    // Regression Correction
    if (!ecorr_.IsInitialized()) {
      if (energyCorrectionsFromDB and regressionVersion == "V3") {
	std::cout << "DB version available only for V2" << std::endl;
	energyCorrectionsFromDB = false;
      }
      
      if (!energyCorrectionsFromDB) {
	char filename[500];
	char* descr = getenv("CMSSW_BASE");
	sprintf(filename, "%s/src/HiggsAnalysis/HiggsTo2photons/data/%s", descr, energyRegFilename.c_str());
	if (fexist(filename)) {
	  sprintf(filename, "http://cern.ch/sani/%s", energyRegFilename.c_str());
	  ecorr_.Initialize(iSetup, filename);
	} else {
	  ecorr_.Initialize(iSetup, filename);
	} 
      } else {
	ecorr_.Initialize(iSetup, "wgbrph", true); //FIXME no wgbrele in DB
      }
    }

    if (regressionVersion == "V3") {
      //std::pair<double,double> cor = ecorr_.CorrectedEnergyWithErrorV3(egsf, *vtxH, rho, ecalLazyTool, iSetup);
      el_regr_energy[el_n]    = 0;//cor.first;
      el_regr_energyerr[el_n] = 0;//cor.second;
    } else {
      std::pair<double,double> cor = ecorr_.CorrectedEnergyWithError(egsf, *vtxH, ecalLazyTool, iSetup);
      el_regr_energy[el_n]    = cor.first;
      el_regr_energyerr[el_n] = cor.second;
    }
    
    reco::GsfElectronRef calibEleRef(calibEleH, std::distance(elH->begin(), igsf));
    reco::GsfElectronRef myElectronRef(elH, std::distance(elH->begin(), igsf));
    el_corr_energy[el_n]     = (*corrEnergy)[myElectronRef];
    el_corr_energyerr[el_n]  = (*corrEnergyErr)[myElectronRef];
    el_calib_energy[el_n]    = (*calibEnergy)[calibEleRef];
    el_calib_energyerr[el_n] = (*calibEnergyErr)[calibEleRef];

    new ((*el_p4_corr)[el_n]) TLorentzVector();
    ((TLorentzVector *)el_p4_corr->At(el_n))->SetXYZT(calibEleRef->px(), calibEleRef->py(), calibEleRef->pz(), calibEleRef->energy());

    // ES variables
    el_eseffsixix[el_n] = 0.;
    el_eseffsiyiy[el_n] = 0.;
    if (ESRecHits.isValid() && (fabs(egsf.superCluster()->eta()) > 1.6 && fabs(egsf.superCluster()->eta()) < 3)) {
      std::vector<float> elESHits0 = gES->getESHits(egsf.superCluster()->x(), egsf.superCluster()->y(), egsf.superCluster()->z(), rechits_map_, geometry, topology_p, 0);
      std::vector<float> elESShape = gES->getESShape(elESHits0);
      el_eseffsixix[el_n] = elESShape[0];
      el_eseffsiyiy[el_n] = elESShape[1];
    }

    el_hoe[el_n] = egsf.hcalOverEcal();
    el_hoed1[el_n] = egsf.hcalDepth1OverEcal();
    el_hoed2[el_n] = egsf.hcalDepth2OverEcal();

    el_hoebc[el_n] = egsf.hcalOverEcalBc();
    el_hoebcd1[el_n] = egsf.hcalDepth1OverEcalBc();
    el_hoebcd2[el_n] = egsf.hcalDepth2OverEcalBc();

    el_d0[el_n] = egsf.gsfTrack()->d0();
    el_z0[el_n] = egsf.gsfTrack()->dz(); 
    el_chi2[el_n] = egsf.gsfTrack()->chi2();
    el_dof[el_n] = egsf.gsfTrack()->ndof();
    el_validhits[el_n] = egsf.gsfTrack()->numberOfValidHits();
    el_losthits[el_n] = egsf.gsfTrack()->numberOfLostHits();

    math::XYZPoint vtxPoint(0.0,0.0,0.0);
    if (vtxH->size() != 0) {
      reco::VertexRef vtx(vtxH, 0);
      vtxPoint = math::XYZPoint(vtx->x(),vtx->y(),vtx->z());
    }
    el_ip_gsf[el_n] = (-1.)*egsf.gsfTrack()->dxy(vtxPoint);

    if (!doAodSim && (trackColl2.encode() != "electronGsfTracks")) {
      for(unsigned int j=0; j<tkH2->size(); j++) { 
        reco::GsfTrackRef tk2(tkH2, j);
        std::pair<unsigned int, float> result = sharedHits(*tk2, *(egsf.gsfTrack()));
        
        if (result.second > .999) {
          el_hp_expin[el_n] = tk2->trackerExpectedHitsInner().numberOfHits();
          el_hp_expout[el_n] = tk2->trackerExpectedHitsOuter().numberOfHits();
          break;
        }
      }
    } else {
      el_hp_expin[el_n] = egsf.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
      el_hp_expout[el_n] = egsf.gsfTrack()->trackerExpectedHitsOuter().numberOfHits();
    }
 
    el_charge[el_n] = egsf.charge();
    
    el_ch_scpix[el_n] = egsf.scPixCharge();
    el_ch_gsf[el_n] = egsf.gsfTrack()->charge();

    if (egsf.isEB())
      el_crack[el_n] = 0;

    if (egsf.isEE())
      el_crack[el_n] = 1;

    if (egsf.isEBEEGap())
      el_crack[el_n] = 3;

    if (egsf.isEBGap())
      el_crack[el_n] = 4;

    if (egsf.isEEGap())
      el_crack[el_n] = 5;

    //int cmsTkind = -1;
    el_scind[el_n] = -1;
    
    el_gsfchi2[el_n] = egsf.gsfTrack()->normalizedChi2();

    if(egsf.closestCtfTrackRef().isNonnull()) {
      const double gsfsign = ((-egsf.gsfTrack()->dxy(vtxPoint)) >=0 ) ? 1. : -1.;
      const reco::TransientTrack& tt = transientTrackBuilder->build(egsf.gsfTrack()); 
      reco::VertexRef vtx(vtxH, 0);
      const std::pair<bool, Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt, *vtx);

      if (ip3dpv.first) {
	double ip3d = gsfsign*ip3dpv.second.value();
	el_ip3d_err[el_n] = ip3dpv.second.error();  
	el_ip3d[el_n] = ip3d; 
	el_ip3d_sig[el_n] = ip3d/el_ip3d_err[el_n];
      }
      
      el_kfhits[el_n] = egsf.closestCtfTrackRef()->hitPattern().trackerLayersWithMeasurement();
      el_kfchi2[el_n] = egsf.closestCtfTrackRef()->normalizedChi2();
      for(unsigned int j=0; j<tkH->size(); ++j) {
        reco::TrackRef tk(tkH, j);
        if(gCUT->cut(*tk))
          continue; 
        if (tk == egsf.closestCtfTrackRef()) {
          el_tkind[el_n] = j;
          el_ip_ctf[el_n] = (-1.)*egsf.closestCtfTrackRef()->dxy(vtxPoint);
        }
      }
    } else {
      el_kfhits[el_n] = -1;
      el_kfchi2[el_n] = 0.;
      el_tkind[el_n] = -1;
      el_ip_ctf[el_n] = -9999.;
      el_ip3d[el_n] = -999.;
      el_ip3d_err[el_n] = -999.;
      el_ip3d_sig[el_n] = 0.0;
    }
    
    int index = 0;
    // loop over the two SC collections
    for(int z = 0; z<2; ++z) {
      if (z == 0) {
        for(reco::SuperClusterCollection::size_type j = 0; j<superClustersBarrelH->size(); ++j) {
          
          reco::SuperClusterRef cluster(superClustersBarrelH, j);
          // apply the cuts
          if(gCUT->cut(*cluster))
	    continue;
          // passed cuts

          if (&(*egsf.superCluster()) == &(*cluster)) {
            el_scind[el_n] = index; 
	    el_sieiesc[el_n] = sqrt(EcalClusterTools::scLocalCovariances(*(cluster), &(*barrelRecHits), &(*topology))[0]);
	    std::vector<float> vCov = ecalLazyTool.localCovariances(*(cluster->seed()));
	    //std::vector<float> vCov = EcalClusterTools::localCovariances( *(cluster->seed()), &(*barrelRecHits), &(*topology));
	    if (!isnan(vCov[2]))
	      el_sipip[el_n] = sqrt(vCov[2]);
	    el_sieip[el_n] = vCov[1];
            break;
          }
          index++;
        }
      }
      
      if (z == 1) {
        for(reco::SuperClusterCollection::size_type j = 0; j<superClustersEndcapH->size(); ++j) {
          
          reco::SuperClusterRef cluster(superClustersEndcapH, j);
          // apply the cuts
          if(gCUT->cut(*cluster))continue;
          // passed cuts 

          if (&(*(egsf.superCluster())) == &(*cluster)) {
            el_scind[el_n] = index;
            el_sieiesc[el_n] = sqrt(EcalClusterTools::scLocalCovariances(*(cluster), &(*endcapRecHits), &(*topology))[0]);
	    std::vector<float> vCov = ecalLazyTool.localCovariances(*(cluster->seed()));
	    if (!isnan(vCov[2]))
	      el_sipip[el_n] = sqrt(vCov[2]);
	    el_sieip[el_n] = vCov[1];
            break;
          }
          index++;
        }
      }
    }

    if (el_scind[el_n] == -1) {
      if(fabs(egsf.superCluster()->eta()) < 1.479) {
        el_sieiesc[el_n] = sqrt(EcalClusterTools::scLocalCovariances(*(egsf.superCluster()), &(*barrelRecHits), &(*topology))[0]);
      } else {
        el_sieiesc[el_n] = sqrt(EcalClusterTools::scLocalCovariances(*(egsf.superCluster()), &(*endcapRecHits), &(*topology))[0]);
      }
    }

    std::vector<UInt_t> schits, bchits;    
    for (reco::CaloCluster_iterator clusterIt=egsf.superCluster()->clustersBegin(); clusterIt!=egsf.superCluster()->clustersEnd(); clusterIt++) {
      for (unsigned int j=0; j<(*clusterIt)->hitsAndFractions().size(); j++)
	schits.push_back(((*clusterIt)->hitsAndFractions())[j].first);
      el_schits->push_back(schits);
    }

    for (unsigned int j=0; j<egsf.superCluster()->seed()->hitsAndFractions().size(); j++)
      bchits.push_back((egsf.superCluster()->seed()->hitsAndFractions())[j].first);
    el_bchits->push_back(bchits);

    el_mva_nontrig[el_n] = myMVANonTrig->mvaValue(egsf, 
						  vtxH->front(), 
						  thebuilder,
						  ecalLazyTool,
						  false);
    
    el_mva_trig[el_n] = myMVATrig->mvaValue(egsf,
					    vtxH->front(),
					    thebuilder,
					    ecalLazyTool,
					    false);
    
    //std::cout<<"el_n el_mva_trig el_mva_nontrig "<<el_n<<" "<<el_mva_trig[el_n]<<" "<<el_mva_nontrig[el_n]<<std::endl;

    /* NEW variables */
    
    if (!egsf.pflowSuperCluster().isNull()) {
      el_nbrempf[el_n] = egsf.pflowSuperCluster()->clustersSize();
      el_eseedpf[el_n] = egsf.pflowSuperCluster()->seed()->energy();
      el_epf[el_n] = egsf.pflowSuperCluster()->energy();
      el_psenergypf[el_n] = egsf.pflowSuperCluster()->preshowerEnergy();
    } else {
      el_nbrempf[el_n] = -1;
      el_eseedpf[el_n] = -1;
      el_epf[el_n] = -1;
      el_psenergypf[el_n] = -1;
    }
    
    el_passcutpresel[el_n] = (Int_t)egsf.passingCutBasedPreselection();
    el_passmvapresel[el_n] = (Int_t)egsf.passingMvaPreselection();
    el_psenergy[el_n] = egsf.superCluster()->preshowerEnergy();  

    el_psly1[el_n] = 0;
    el_psly2[el_n] = 0;
    el_psnstriply1[el_n] = 0;
    el_psnstriply2[el_n] = 0;
    for (reco::CaloCluster_iterator it=egsf.superCluster()->preshowerClustersBegin();
	 it != egsf.superCluster()->preshowerClustersEnd(); it++) {

      if (ESDetId(((*it)->hitsAndFractions())[0].first).plane() == 1) {
	el_psly1[el_n] += (*it)->energy();
	el_psnstriply1[el_n] += (*it)->hitsAndFractions().size();
      } else {
	el_psly2[el_n] += (*it)->energy();
	el_psnstriply2[el_n] += (*it)->hitsAndFractions().size();
      }
    }

    el_ecaldrv[el_n] = egsf.ecalDrivenSeed();
    el_tkdrv[el_n] = egsf.trackerDrivenSeed();
    
    el_pfiso_charged[el_n] = egsf.pfIsolationVariables().chargedHadronIso;//(*(electronIsoVals[0].product()))[myElectronRef]; //
    el_pfiso_photon[el_n]  = egsf.pfIsolationVariables().photonIso	 ;//(*(electronIsoVals[1].product()))[myElectronRef]; //
    el_pfiso_neutral[el_n] = egsf.pfIsolationVariables().neutralHadronIso;//(*(electronIsoVals[2].product()))[myElectronRef]; //

    el_tkiso04[el_n]   = egsf.dr04TkSumPt();
    el_ecaliso04[el_n] = egsf.dr04EcalRecHitSumEt();
    el_hcaliso04[el_n] = egsf.dr04HcalTowerSumEt();    

    el_tkiso03[el_n]   = egsf.dr03TkSumPt();
    el_ecaliso03[el_n] = egsf.dr03EcalRecHitSumEt();
    el_hcaliso03[el_n] = egsf.dr03HcalTowerSumEt();

    //EgammaTowerIsolation hcaliso03(0.3, 0., 0, -1, ctH.product());
    //EgammaTowerIsolation hcaliso04(0.4, 0., 0, -1, ctH.product());
    //
    //el_hcalsolidiso03[el_n] = hcaliso03.getTowerEtSum(&(*(egsf.superCluster())) );
    //el_hcalsolidiso04[el_n] = hcaliso04.getTowerEtSum(&(*(egsf.superCluster())) );
    //
    //el_hcalbciso03[el_n] = egsf.dr03HcalTowerSumEtBc();
    //el_hcalbciso04[el_n] = egsf.dr04HcalTowerSumEtBc();

    // Fill out electron identification
    std::vector<edm::Handle<edm::ValueMap<float> > > eIDVM(9); 
    std::vector<int> results;
    
    for(unsigned int j=0; j<eIDLabels.size(); j++) {
      if (iEvent.getByLabel(eIDLabels[j], eIDVM[j])) {  
        const edm::ValueMap<float>& eIDmapTemp = *eIDVM[j];
        reco::GsfElectronRef electronRef(elH, std::distance(elH->begin(), igsf));
        results.push_back((Int_t)eIDmapTemp[electronRef]);
      } else {
        results.push_back(-1);
      }
    }
    
    el_catbased->push_back(results);

    el_dist[el_n] = (egsf.convDist() == -9999.? 9999:egsf.convDist());
    el_dcot[el_n] = (egsf.convDcot() == -9999.? 9999:egsf.convDcot());
  
    // loop through vertices for d0 and dZ w.r.t. each vertex
    // need number of vertices and vertices' positions
    int maxV = std::min(100, (int)vtxH->size());
    //std::cout << maxV << std::endl;
    for(int iv=0; iv<maxV; iv++){
      reco::VertexRef v(vtxH, iv);
      //std::cout << v->x() << std::endl;
      //std::cout << v->y() << std::endl;
      //std::cout << v->z() << std::endl;
      math::XYZPoint vtxPoint = math::XYZPoint(v->x(), v->y(), v->z());
      
      el_D0Vtx[el_n][iv] = egsf.gsfTrack()->dxy(vtxPoint);
      el_DZVtx[el_n][iv] = egsf.gsfTrack()->dz(vtxPoint);
    }
    
    el_n++;
  } // end of loop over electrons (reco::GsfElectronCollection)


  delete topology_p;
  return true;
}

std::pair<unsigned int, float> GlobeElectrons::sharedHits(const reco::Track& trackA, const reco::Track& trackB) {
  
  unsigned int shared = 0;
  for(trackingRecHit_iterator tkHitA = trackA.recHitsBegin(); tkHitA !=trackA.recHitsEnd(); ++tkHitA){
    for(trackingRecHit_iterator tkHitB = trackB.recHitsBegin();
        tkHitB !=trackB.recHitsEnd(); ++tkHitB){
      if( (**tkHitA).isValid() && (**tkHitB).isValid() &&(**tkHitA).sharesInput( &(**tkHitB),TrackingRecHit::all)) {
        shared++;
        break;
      }
    }
  }

  float fraction = (float) shared/std::min(trackA.found(),trackB.found());
  return std::make_pair(shared,fraction);
}



