#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePhotons.h"

#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/Tools.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/EgammaCandidates/interface/PhotonPi0DiscriminatorAssociation.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "DataFormats/EgammaTrackReco/interface/TrackCaloClusterAssociation.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFClusterWidthAlgo.h"
//#include "HiggsAnalysis/HiggsTo2photons/interface/pfFrixioneIso.h"
//#include "HiggsAnalysis/HiggsToGammaGamma/interface/PhotonFix.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/PFIsolation.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Mustache.h"
//#include "RecoEgamma/EgammaTools/interface/ggPFPhotons.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TFile.h"
#include <cstdlib>

#include "RecoEcal/EgammaCoreTools/interface/MyClusterTools.h"

void GlobePhotons::checkSetup(const edm::EventSetup& iSetup) {

  // Initialise the Correction Scheme
  fEtaCorr->init(iSetup);
  CrackCorr->init(iSetup);
  LocalCorr->init(iSetup);

  // Transform Track into TransientTrack (needed by the Vertex fitter)
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theTTkBuilder);

  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  geometry = *geoHandle;
    
  //edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
  //iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);
  //sevLevel = sevlv.product(); 

  if (haveTowers) {
    hcalHelper->checkSetup(iSetup);
    hcalHelperPflow->checkSetup(iSetup);
  }
}

GlobePhotons::~GlobePhotons() { 
  if (haveTowers) {
    delete hcalHelper;
    delete hcalHelperPflow;
  }
}

GlobePhotons::GlobePhotons(const edm::ParameterSet& iConfig, const char* n): nome(n) {

  isInitialized = false;

  debug_level = iConfig.getParameter<int>("Debug_Level");
  doFastSim = iConfig.getParameter<bool>("doFastSim");
  doAodSim = iConfig.getParameter<bool>("doAodSim");

  // PHOTONS 
  photonCollStd =  iConfig.getParameter<edm::InputTag>("PhotonCollStd");
  photonCollPf =  iConfig.getParameter<edm::InputTag>("PhotonCollPf");
  // SUPER CLUSTERS
  hybridSuperClusterColl = iConfig.getParameter<edm::InputTag>("HybridSuperClusterColl");
  endcapSuperClusterColl = iConfig.getParameter<edm::InputTag>("EndcapSuperClusterColl");

  ecalHitEBColl = iConfig.getParameter<edm::InputTag>("EcalHitEBColl");
  ecalHitEEColl = iConfig.getParameter<edm::InputTag>("EcalHitEEColl");
  ecalHitESColl = iConfig.getParameter<edm::InputTag>("EcalHitESColl");
  ecalHitEBCleanedColl = iConfig.getParameter<edm::InputTag>("EcalHitEBCleanedColl");
  ecalHitEECleanedColl = iConfig.getParameter<edm::InputTag>("EcalHitEECleanedColl");

  hcalBEColl =  iConfig.getParameter<edm::InputTag>("HcalHitsBEColl");
  //hcalFColl =  iConfig.getParameter<edm::InputTag>("HcalHitsFColl");
  //hcalHoColl =  iConfig.getParameter<edm::InputTag>("HcalHitsHoColl");

  convertedPhotonColl =  iConfig.getParameter<edm::InputTag>("ConvertedPhotonColl");
  beamSpotColl =  iConfig.getParameter<edm::InputTag>("BeamSpot");
  electronColl =  iConfig.getParameter<edm::InputTag>("ElectronColl_std");

  rhoCollection = iConfig.getParameter<edm::InputTag>("rhoCollection_algo1");
  vtxCollection = iConfig.getParameter<edm::InputTag>("VertexColl_std");
  tkCollection  = iConfig.getParameter<edm::InputTag>("tkColl");

  //hcalHitColl = iConfig.getParameter<edm::InputTag>("HcalHitsBEColl");
  pfColl = iConfig.getParameter<edm::InputTag>("PFCandidateColl");

  edm::ParameterSet isoVals03  = iConfig.getParameter<edm::ParameterSet> ("isolationValues03");
  inputTagIsoVals03_.push_back(isoVals03.getParameter<edm::InputTag>("pfChargedHadrons"));
  inputTagIsoVals03_.push_back(isoVals03.getParameter<edm::InputTag>("pfPhotons"));
  inputTagIsoVals03_.push_back(isoVals03.getParameter<edm::InputTag>("pfNeutralHadrons"));
  
  edm::ParameterSet isoVals04  = iConfig.getParameter<edm::ParameterSet> ("isolationValues04");
  inputTagIsoVals04_.push_back(isoVals04.getParameter<edm::InputTag>("pfChargedHadrons"));
  inputTagIsoVals04_.push_back(isoVals04.getParameter<edm::InputTag>("pfPhotons"));
  inputTagIsoVals04_.push_back(isoVals04.getParameter<edm::InputTag>("pfNeutralHadrons"));

  energyCorrectionsFromDB = iConfig.getParameter<bool> ("energyCorrectionsFromDB"); 
  energyRegFilename       = iConfig.getParameter<std::string> ("energyCorrectionsFileNamePho"); 
  regressionVersion       = iConfig.getParameter<std::string> ("energyCorrectionsVersion");

  // get the Correction Functions
  fEtaCorr  = EcalClusterFunctionFactory::get()->create("EcalClusterEnergyCorrection",iConfig);
  CrackCorr = EcalClusterFunctionFactory::get()->create("EcalClusterCrackCorrection",iConfig);
  LocalCorr = EcalClusterFunctionFactory::get()->create("EcalClusterLocalContCorrection",iConfig);

  // Initialise PhotonFix
  //PhotonFix::initialiseParameters(iConfig);

  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);
  gES  = new GlobeEcalClusters(iConfig);

  pho_pfiso_mycharged01 = new std::vector<std::vector<float> >();
  pho_pfiso_mycharged02 = new std::vector<std::vector<float> >();
  pho_pfiso_mycharged03 = new std::vector<std::vector<float> >();
  pho_pfiso_mycharged04 = new std::vector<std::vector<float> >();
  pho_pfiso_mycharged05 = new std::vector<std::vector<float> >();
  pho_pfiso_mycharged06 = new std::vector<std::vector<float> >();

  pho_schits = new std::vector<std::vector<UInt_t> >;
  pho_bchits = new std::vector<std::vector<UInt_t> >;

  //pho_frixiso = new std::vector<std::vector<float> >();
 
  cicPhotonId = new CiCPhotonID(iConfig);
  //pfFrixIso = new pfFrixioneIso(); 

  haveTowers = iConfig.getParameter<bool>("doCaloTower");
  
  if (haveTowers) {
    hcalCfg.hOverEConeSize = 0.15;
    hcalCfg.useTowers = true;
    hcalCfg.hcalTowers = iConfig.getParameter<edm::InputTag>("CaloTowerColl");
    hcalCfg.hOverEPtMin = 0;
    hcalCfgPflow.hOverEConeSize = 0.15;
    hcalCfgPflow.useTowers = true ;
    hcalCfgPflow.hcalTowers = iConfig.getParameter<edm::InputTag>("CaloTowerColl");
    hcalCfgPflow.hOverEPtMin = 0;
    
    hcalHelper = new ElectronHcalHelper(hcalCfg);
    hcalHelperPflow = new ElectronHcalHelper(hcalCfgPflow);
  }
}


void GlobePhotons::defineBranch(GlobeAnalyzer* ana) {

  pho_p4 = new TClonesArray("TLorentzVector", MAX_PHOTONS);
  pho_calopos = new TClonesArray("TVector3", MAX_PHOTONS);

  ana->Branch("pho_n", &pho_n, "pho_n/I");

  //Correction Schemes
  ana->Branch("pho_feta",&pho_feta,"pho_feta[pho_n][5]/F");
  ana->Branch("pho_crackcorr",&pho_crackcorr,"pho_crackcorr[pho_n]/F");
  ana->Branch("pho_localcorr",&pho_localcorr,"pho_localcorr[pho_n]/F");

  //pf
  ana->Branch("pho_isPFPhoton", &pho_isPFPhoton,"pho_isPFPhoton[pho_n]/I");
  ana->Branch("pho_isPFElectron", &pho_isPFElectron,"pho_isPFElectron[pho_n]/I");

  //fiducial flags
  ana->Branch("pho_isEB",&pho_isEB,"pho_isEB[pho_n]/I");
  ana->Branch("pho_isEE",&pho_isEE,"pho_isEE[pho_n]/I");
  ana->Branch("pho_isEBGap",&pho_isEBGap,"pho_isEBGap[pho_n]/I");
  ana->Branch("pho_isEEGap",&pho_isEEGap,"pho_isEEGap[pho_n]/I");
  ana->Branch("pho_isEBEEGap",&pho_isEBEEGap,"pho_isEBEEGap[pho_n]/I");
  ana->Branch("pho_isEBEtaGap",&pho_isEBEtaGap,"pho_isEBEtaGap[pho_n]/I");    
  ana->Branch("pho_isEBPhiGap",&pho_isEBPhiGap,"pho_isEBPhiGap[pho_n]/I");    
  ana->Branch("pho_isEEDeeGap",&pho_isEEDeeGap,"pho_isEEDeeGap[pho_n]/I");    
  ana->Branch("pho_isEERingGap",&pho_isEERingGap,"pho_isEERingGap[pho_n]/I");   

  //shower shape variables
  ana->Branch("pho_see",&pho_see,"pho_see[pho_n]/F");
  ana->Branch("pho_sieie",&pho_sieie,"pho_sieie[pho_n]/F");
  ana->Branch("pho_e1x5",&pho_e1x5,"pho_e1x5[pho_n]/F");
  ana->Branch("pho_e1x3",&pho_e1x3,"pho_e1x3[pho_n]/F");
  ana->Branch("pho_e2x2",&pho_e2x2,"pho_e2x2[pho_n]/F");
  ana->Branch("pho_e3x3",&pho_e3x3,"pho_e3x3[pho_n]/F");
  ana->Branch("pho_e5x5",&pho_e5x5,"pho_e5x5[pho_n]/F");
  ana->Branch("pho_emaxxtal",&pho_emaxxtal,"pho_emaxxtal[pho_n]/F");
  ana->Branch("pho_hoe",&pho_hoe,"pho_hoe[pho_n]/F");
  ana->Branch("pho_h1oe",&pho_h1oe,"pho_h1oe[pho_n]/F");
  ana->Branch("pho_h2oe",&pho_h2oe,"pho_h2oe[pho_n]/F");
  ana->Branch("pho_hoe_bc",&pho_hoe_bc,"pho_hoe_bc[pho_n]/F");
  ana->Branch("pho_h1oe_bc",&pho_h1oe_bc,"pho_h1oe_bc[pho_n]/F");
  ana->Branch("pho_h2oe_bc",&pho_h2oe_bc,"pho_h2oe_bc[pho_n]/F");
  ana->Branch("pho_r1x5", &pho_r1x5, "pho_r1x5[pho_n]/F");
  ana->Branch("pho_r2x5", &pho_r2x5, "pho_r2x5[pho_n]/F");
  ana->Branch("pho_r9", &pho_r9,"pho_r9[pho_n]/F");
  ana->Branch("pho_smaj", &pho_smaj, "pho_smaj[pho_n]/F");
  ana->Branch("pho_alpha", &pho_alpha, "pho_alpha[pho_n]/F");
  ana->Branch("pho_smin", &pho_smin, "pho_smin[pho_n]/F");

  ana->Branch("pho_r9_cleaned", &pho_r9_cleaned,"pho_r9_cleaned[pho_n]/F");
  ana->Branch("pho_sieie_cleaned",&pho_sieie_cleaned,"pho_sieie_cleaned[pho_n]/F");
  ana->Branch("pho_sieip_cleaned",&pho_sieip_cleaned,"pho_sieip_cleaned[pho_n]/F");
  ana->Branch("pho_e2x2_cleaned",&pho_e2x2_cleaned,"pho_e2x2_cleaned[pho_n]/F");
  ana->Branch("pho_e5x5_cleaned",&pho_e5x5_cleaned,"pho_e5x5_cleaned[pho_n]/F");

  ana->Branch("pho_r9_frac", &pho_r9_frac,"pho_r9_frac[pho_n]/F");
  ana->Branch("pho_sieie_frac",&pho_sieie_frac,"pho_sieie_frac[pho_n]/F");
  ana->Branch("pho_sieip_frac",&pho_sieip_frac,"pho_sieip_frac[pho_n]/F");
  ana->Branch("pho_sipip_frac",&pho_sipip_frac,"pho_sipip_frac[pho_n]/F");
  ana->Branch("pho_e2x2_frac",&pho_e2x2_frac,"pho_e2x2_frac[pho_n]/F");
  ana->Branch("pho_e5x5_frac",&pho_e5x5_frac,"pho_e5x5_frac[pho_n]/F");
  ana->Branch("pho_e3x3_frac",&pho_e3x3_frac,"pho_e3x3_frac[pho_n]/F");
  ana->Branch("pho_smaj_frac", &pho_smaj_frac, "pho_smaj_frac[pho_n]/F");
  ana->Branch("pho_alpha_frac", &pho_alpha_frac, "pho_alpha_frac[pho_n]/F");
  ana->Branch("pho_smin_frac", &pho_smin_frac, "pho_smin_frac[pho_n]/F");


  
  ana->Branch("pho_eseffsixix",&pho_eseffsixix,"pho_eseffsixix[pho_n]/F");
  ana->Branch("pho_eseffsiyiy",&pho_eseffsiyiy,"pho_eseffsiyiy[pho_n]/F");

  // NN variable
  ana->Branch("pho_r19", &pho_r19, "pho_r19[pho_n]/F");
  ana->Branch("pho_maxoraw", &pho_maxoraw, "pho_maxoraw[pho_n]/F");
  ana->Branch("pho_cep", &pho_cep, "pho_cep[pho_n]/F");
  ana->Branch("pho_lambdaratio", &pho_lambdaratio, "pho_lambdaratio[pho_n]/F");
  ana->Branch("pho_lambdadivcov", &pho_lambdadivcov, "pho_lambdadivcov[pho_n]/F");
  ana->Branch("pho_cep_global", &pho_cep_global, "pho_cep_global[pho_n]/F");
  ana->Branch("pho_lambdaratio_global", &pho_lambdaratio_global, "pho_lambdaratio_global[pho_n]/F");
  ana->Branch("pho_lambdadivcov_global", &pho_lambdadivcov_global, "pho_lambdadivcov_global[pho_n]/F");
  ana->Branch("pho_etawidth", &pho_etawidth, "pho_etawidth[pho_n]/F");
  ana->Branch("pho_brem", &pho_brem, "pho_brem[pho_n]/F");


  //isolation variables

  ana->Branch("pho_pfiso_myneutral01", &pho_pfiso_myneutral01, "pho_pfiso_myneutral01[pho_n]/F");
  ana->Branch("pho_pfiso_myphoton01", &pho_pfiso_myphoton01, "pho_pfiso_myphoton01[pho_n]/F");  
  ana->Branch("pho_pfiso_mycharged01", "std::vector<std::vector<float> >", &pho_pfiso_mycharged01);

  ana->Branch("pho_pfiso_myneutral02", &pho_pfiso_myneutral02, "pho_pfiso_myneutral02[pho_n]/F");
  ana->Branch("pho_pfiso_myphoton02", &pho_pfiso_myphoton02, "pho_pfiso_myphoton02[pho_n]/F");  
  ana->Branch("pho_pfiso_mycharged02", "std::vector<std::vector<float> >", &pho_pfiso_mycharged02);

  ana->Branch("pho_pfiso_myneutral03", &pho_pfiso_myneutral03, "pho_pfiso_myneutral03[pho_n]/F");
  ana->Branch("pho_pfiso_myphoton03", &pho_pfiso_myphoton03, "pho_pfiso_myphoton03[pho_n]/F");  
  ana->Branch("pho_pfiso_mycharged03", "std::vector<std::vector<float> >", &pho_pfiso_mycharged03);

  ana->Branch("pho_pfiso_myneutral04", &pho_pfiso_myneutral04, "pho_pfiso_myneutral04[pho_n]/F");
  ana->Branch("pho_pfiso_myphoton04", &pho_pfiso_myphoton04, "pho_pfiso_myphoton04[pho_n]/F");
  ana->Branch("pho_pfiso_mycharged04", "std::vector<std::vector<float> >", &pho_pfiso_mycharged04);

  ana->Branch("pho_pfiso_myneutral05", &pho_pfiso_myneutral05, "pho_pfiso_myneutral05[pho_n]/F");
  ana->Branch("pho_pfiso_myphoton05", &pho_pfiso_myphoton05, "pho_pfiso_myphoton05[pho_n]/F");
  ana->Branch("pho_pfiso_mycharged05", "std::vector<std::vector<float> >", &pho_pfiso_mycharged05);

  ana->Branch("pho_pfiso_myneutral06", &pho_pfiso_myneutral06, "pho_pfiso_myneutral06[pho_n]/F");
  ana->Branch("pho_pfiso_myphoton06", &pho_pfiso_myphoton06, "pho_pfiso_myphoton06[pho_n]/F");
  ana->Branch("pho_pfiso_mycharged06", "std::vector<std::vector<float> >", &pho_pfiso_mycharged06);

  ana->Branch("pho_schits", "std::vector<std::vector<UInt_t> >", &pho_schits);
  ana->Branch("pho_bchits", "std::vector<std::vector<UInt_t> >", &pho_bchits);

  //ana->Branch("pho_frixiso", "std::vector<std::vector<float> >", &pho_frixiso);  

  ana->Branch("pho_pfconvVtxZ", &pho_pfconvVtxZ, "pho_pfconvVtxZ[pho_n]/F");
  ana->Branch("pho_pfconvVtxZErr", &pho_pfconvVtxZErr, "pho_pfconvVtxZErr[pho_n]/F");
  ana->Branch("pho_hasConvPf", &pho_hasConvPf, "pho_hasConvPf[pho_n]/I");
  ana->Branch("pho_hasSLConvPf", &pho_hasSLConvPf, "pho_hasSLConvPf[pho_n]/I");

  ana->Branch("pho_must", &pho_must, "pho_must[pho_n]/F");
  ana->Branch("pho_mustnc", &pho_mustnc, "pho_mustnc[pho_n]/I");
  ana->Branch("pho_pfpresh1", &pho_pfpresh1, "pho_pfpresh1[pho_n]/F");
  ana->Branch("pho_pfpresh2", &pho_pfpresh2, "pho_pfpresh2[pho_n]/F");
  ana->Branch("pho_mustenergy", &pho_mustenergy, "pho_mustenergy[pho_n]/F");
  ana->Branch("pho_mustenergyout", &pho_mustenergyout, "pho_mustenergyout[pho_n]/F");
  
  ana->Branch("pho_mustEtout", &pho_mustEtout, "pho_mustEtout[pho_n]/F");
  ana->Branch("pho_pflowE", &pho_pflowE, "pho_pflowE[pho_n]/F");
  ana->Branch("pho_pfdeta", &pho_pfdeta, "pho_pfdeta[pho_n]/F");
  ana->Branch("pho_pfdphi", &pho_pfdphi, "pho_pfdphi[pho_n]/F");
  ana->Branch("pho_pfclusrms", &pho_pfclusrms, "pho_pfclusrms[pho_n]/F");
  ana->Branch("pho_pfclusrmsmust", &pho_pfclusrmsmust, "pho_pfclusrmsmust[pho_n]/F");
  ana->Branch("pho_pfClusECorr", &pho_pfClusECorr, "pho_pfClusECorr[pho_n]/F");  
  ana->Branch("pho_pfMatch", &pho_pfMatch, "pho_pfMatch[pho_n]/I");
  ana->Branch("pho_PfEleVeto", &pho_PfEleVeto, "pho_PfEleVeto[pho_n]/I");

  ana->Branch("pho_pfRawEnergy", &pho_pfRawEnergy, "pho_pfRawEnergy[pho_n]/F");
  ana->Branch("pho_pfe2x2", &pho_pfe2x2, "pho_pfe2x2[pho_n]/F");
  ana->Branch("pho_pfe3x3", &pho_pfe3x3, "pho_pfe3x3[pho_n]/F");
  ana->Branch("pho_pfe5x5", &pho_pfe5x5, "pho_pfe5x5[pho_n]/F");
  ana->Branch("pho_pfsieie", &pho_pfsieie, "pho_pfsieie[pho_n]/F");
  ana->Branch("pho_pfsieip", &pho_pfsieip, "pho_pfsieip[pho_n]/F");
  ana->Branch("pho_pfsipip", &pho_pfsipip, "pho_pfsipip[pho_n]/F");
  ana->Branch("pho_pfemaxxtal", &pho_pfemaxxtal, "pho_pfemaxxtal[pho_n]/F");
  ana->Branch("pho_pfe2nd", &pho_pfe2nd, "pho_pfe2nd[pho_n]/F");

  ana->Branch("pho_ecalsumetconedr04",&pho_ecalsumetconedr04,"pho_ecalsumetconedr04[pho_n]/F");
  ana->Branch("pho_hcalsumetconedr04",&pho_hcalsumetconedr04,"pho_hcalsumetconedr04[pho_n]/F");
  ana->Branch("pho_hcal1sumetconedr04",&pho_hcal1sumetconedr04,"pho_hcal1sumetconedr04[pho_n]/F");
  ana->Branch("pho_hcal2sumetconedr04",&pho_hcal2sumetconedr04,"pho_hcal2sumetconedr04[pho_n]/F");
  ana->Branch("pho_trksumptsolidconedr04",&pho_trksumptsolidconedr04,"pho_trksumptsolidconedr04[pho_n]/F");
  ana->Branch("pho_trksumpthollowconedr04",&pho_trksumpthollowconedr04,"pho_trksumpthollowconedr04[pho_n]/F");
  ana->Branch("pho_ntrksolidconedr04",&pho_ntrksolidconedr04,"pho_ntrksolidconedr04[pho_n]/F");
  ana->Branch("pho_ntrkhollowconedr04",&pho_ntrkhollowconedr04,"pho_ntrkhollowconedr04[pho_n]/F");
  ana->Branch("pho_ecalsumetconedr03",&pho_ecalsumetconedr03,"pho_ecalsumetconedr03[pho_n]/F");
  ana->Branch("pho_hcalsumetconedr03",&pho_hcalsumetconedr03,"pho_hcalsumetconedr03[pho_n]/F");
  ana->Branch("pho_hcal1sumetconedr03",&pho_hcal1sumetconedr03,"pho_hcal1sumetconedr03[pho_n]/F");
  ana->Branch("pho_hcal2sumetconedr03",&pho_hcal2sumetconedr03,"pho_hcal2sumetconedr03[pho_n]/F");
  ana->Branch("pho_trksumptsolidconedr03",&pho_trksumptsolidconedr03,"pho_trksumptsolidconedr03[pho_n]/F");
  ana->Branch("pho_trksumpthollowconedr03",&pho_trksumpthollowconedr03,"pho_trksumpthollowconedr03[pho_n]/F");
  ana->Branch("pho_ntrksolidconedr03",&pho_ntrksolidconedr03,"pho_ntrksolidconedr03[pho_n]/F");
  ana->Branch("pho_ntrkhollowconedr03",&pho_ntrkhollowconedr03,"pho_ntrkhollowconedr03[pho_n]/F");

  ana->Branch("pho_hcalbcsumetconedr03", &pho_hcalbcsumetconedr03, "pho_hcalbcsumetconedr03[pho_n]/F");
  ana->Branch("pho_hcalbc1sumetconedr03",&pho_hcalbc1sumetconedr03,"pho_hcalbc1sumetconedr03[pho_n]/F");
  ana->Branch("pho_hcalbc2sumetconedr03",&pho_hcalbc2sumetconedr03,"pho_hcalbc2sumetconedr03[pho_n]/F");
  ana->Branch("pho_hcalbcsumetconedr04", &pho_hcalbcsumetconedr04, "pho_hcalbcsumetconedr04[pho_n]/F");
  ana->Branch("pho_hcalbc1sumetconedr04",&pho_hcalbc1sumetconedr04,"pho_hcalbc1sumetconedr04[pho_n]/F");
  ana->Branch("pho_hcalbc2sumetconedr04",&pho_hcalbc2sumetconedr04,"pho_hcalbc2sumetconedr04[pho_n]/F");

  ana->Branch("pho_p4", "TClonesArray", &pho_p4, 32000, 0);
  ana->Branch("pho_calopos", "TClonesArray", &pho_calopos, 32000, 0);
  ana->Branch("pho_barrel", &pho_barrel, "pho_barrel[pho_n]/I");
  ana->Branch("pho_scind", &pho_scind, "pho_scind[pho_n]/I");
  
  ana->Branch("pho_haspixseed",&pho_haspixseed,"pho_haspixseed[pho_n]/I");
  
  ana->Branch("pho_hasconvtks",&pho_hasconvtks,"pho_hasconvtks[pho_n]/I");
  ana->Branch("pho_nconv",&pho_nconv,"pho_nconv[pho_n]/I");
  ana->Branch("pho_conv_ntracks",&pho_conv_ntracks,"pho_conv_ntracks[pho_n]/I");
  ana->Branch("pho_conv_pairinvmass",&pho_conv_pairinvmass,"pho_conv_pairinvmass[pho_n]/F");
  ana->Branch("pho_conv_paircotthetasep",&pho_conv_paircotthetasep,"pho_conv_paircotthetasep[pho_n]/F");
  ana->Branch("pho_conv_eoverp",&pho_conv_eoverp,"pho_conv_eoverp[pho_n]/F");
  ana->Branch("pho_conv_zofprimvtxfromtrks",&pho_conv_zofprimvtxfromtrks,"pho_conv_zofprimvtxfromtrks[pho_n]/F");
  ana->Branch("pho_conv_distofminapproach",&pho_conv_distofminapproach,"pho_conv_distofminapproach[pho_n]/F");
  ana->Branch("pho_conv_dphitrksatvtx",&pho_conv_dphitrksatvtx,"pho_conv_dphitrksatvtx[pho_n]/F");
  ana->Branch("pho_conv_dphitrksatecal",&pho_conv_dphitrksatecal,"pho_conv_dphitrksatecal[pho_n]/F");
  ana->Branch("pho_conv_detatrksatecal",&pho_conv_detatrksatecal,"pho_conv_detatrksatecal[pho_n]/F");
  ana->Branch("pho_conv_tk1_d0",&pho_conv_tk1_d0,"pho_conv_tk1_d0[pho_n]/F");
  ana->Branch("pho_conv_tk1_pout",&pho_conv_tk1_pout,"pho_conv_tk1_pout[pho_n]/F");
  ana->Branch("pho_conv_tk1_pin",&pho_conv_tk1_pin,"pho_conv_tk1_pin[pho_n]/F");
  ana->Branch("pho_conv_tk2_d0",&pho_conv_tk2_d0,"pho_conv_tk2_d0[pho_n]/F");
  ana->Branch("pho_conv_tk2_pout",&pho_conv_tk2_pout,"pho_conv_tk2_pout[pho_n]/F");
  ana->Branch("pho_conv_tk2_pin",&pho_conv_tk2_pin,"pho_conv_tk2_pin[pho_n]/F");
  
  //added by marco
  ana->Branch("pho_conv_tk1_dz",&pho_conv_tk1_dz,"pho_conv_tk1_dz[pho_n]/F");
  ana->Branch("pho_conv_tk2_dz",&pho_conv_tk2_dz,"pho_conv_tk2_dz[pho_n]/F");
  ana->Branch("pho_conv_tk1_dzerr",&pho_conv_tk1_dzerr,"pho_conv_tk1_dzerr[pho_n]/F");
  ana->Branch("pho_conv_tk2_dzerr",&pho_conv_tk2_dzerr,"pho_conv_tk2_dzerr[pho_n]/F");
  ana->Branch("pho_conv_tk1_nh",&pho_conv_tk1_nh,"pho_conv_tk1_nh[pho_n]/I");
  ana->Branch("pho_conv_tk2_nh",&pho_conv_tk2_nh,"pho_conv_tk2_nh[pho_n]/I");
  ana->Branch("pho_conv_chi2",&pho_conv_chi2,"pho_conv_chi2[pho_n]/F");
  ana->Branch("pho_conv_chi2_probability",&pho_conv_chi2_probability,"pho_conv_chi2_probability[pho_n]/F");
  ana->Branch("pho_conv_ch1ch2",&pho_conv_ch1ch2,"pho_conv_ch1ch2[pho_n]/I");
  ana->Branch("pho_conv_validvtx",&pho_conv_validvtx,"pho_conv_validvtx[pho_n]/I");
  ana->Branch("pho_conv_MVALikelihood",&pho_conv_MVALikelihood,"pho_conv_MVALikelihood[pho_n]/I");

  // added by pasquale
  ana->Branch("pho_sipip",&pho_sipip,"pho_sipip[pho_n]/F");
  ana->Branch("pho_sieip",&pho_sieip,"pho_sieip[pho_n]/F");
  ana->Branch("pho_zernike20",&pho_zernike20,"pho_zernike20[pho_n]/F");
  ana->Branch("pho_zernike42",&pho_zernike42,"pho_zernike42[pho_n]/F");
  ana->Branch("pho_e2nd",&pho_e2nd,"pho_e2nd[pho_n]/F");
  ana->Branch("pho_e2x5right",&pho_e2x5right,"pho_e2x5right[pho_n]/F");
  ana->Branch("pho_e2x5left",&pho_e2x5left,"pho_e2x5left[pho_n]/F");
  ana->Branch("pho_e2x5top",&pho_e2x5top,"pho_e2x5top[pho_n]/F");
  ana->Branch("pho_e2x5bottom",&pho_e2x5bottom,"pho_e2x5bottom[pho_n]/F");
  ana->Branch("pho_e2x5max",&pho_e2x5max,"pho_e2x5max[pho_n]/F");
  ana->Branch("pho_eright",&pho_eright,"pho_eright[pho_n]/F");
  ana->Branch("pho_eleft",&pho_eleft,"pho_eleft[pho_n]/F");
  ana->Branch("pho_etop",&pho_etop,"pho_etop[pho_n]/F");
  ana->Branch("pho_ebottom",&pho_ebottom,"pho_ebottom[pho_n]/F");

  ana->Branch("pho_biphi",&pho_biphi,"pho_biphi[pho_n]/I");
  ana->Branch("pho_bieta",&pho_bieta,"pho_bieta[pho_n]/I");
  ana->Branch("pho_betacry",&pho_betacry,"pho_betacry[pho_n]/F");
  ana->Branch("pho_phicry",&pho_bphicry,"pho_bphicry[pho_n]/F");
  //ana->Branch("pho_bthetatilt",&pho_bthetatilt,"pho_bthetatilt[pho_n]/F");
  //ana->Branch("pho_bphitilt",&pho_bphitilt,"pho_bphitilt[pho_n]/F");

  //ana->Branch("pho_e2overe9",&pho_e2overe9,"pho_e2overe9[pho_n]/F");
  ana->Branch("pho_seed_severity",&pho_seed_severity,"pho_seed_severity[pho_n]/F");
  ana->Branch("pho_seed_time",&pho_seed_time,"pho_seed_time[pho_n]/F");
  ana->Branch("pho_seed_outoftimechi2",&pho_seed_outoftimechi2,"pho_seed_outoftimechi2[pho_n]/F");
  ana->Branch("pho_seed_chi2",&pho_seed_chi2,"pho_seed_chi2[pho_n]/F");
  ana->Branch("pho_seed_recoflag",&pho_seed_recoflag,"pho_seed_recoflag[pho_n]/F");

  ana->Branch("pho_isconv", &pho_isconv, "pho_isconv[pho_n]/I");
  ana->Branch("pho_residCorrEnergy", &pho_residCorrEnergy, "pho_residCorrEnergy[pho_n]/F");
  ana->Branch("pho_residCorrResn", &pho_residCorrResn, "pho_residCorrResn[pho_n]/F");

  ana->Branch("pho_regr_energy", &pho_regr_energy, "pho_regr_energy[pho_n]/F");
  ana->Branch("pho_regr_energyerr", &pho_regr_energyerr, "pho_regr_energyerr[pho_n]/F");

  //ana->Branch("pho_id_4cat", &pho_id_4cat, "pho_id_4cat[pho_n][100]/I");
  //ana->Branch("pho_id_6cat", &pho_id_6cat, "pho_id_6cat[pho_n][100]/I");  
  //ana->Branch("pho_id_6catpf", &pho_id_6catpf, "pho_id_6catpf[pho_n][100]/I");
   
  pho_conv_vtx = new TClonesArray("TVector3", MAX_PHOTONS);
  ana->Branch("pho_conv_vtx", "TClonesArray", &pho_conv_vtx, 32000, 0);
  pho_conv_pair_momentum = new TClonesArray("TVector3", MAX_PHOTONS);
  ana->Branch("pho_conv_pair_momentum", "TClonesArray", &pho_conv_pair_momentum, 32000, 0);
  pho_conv_refitted_momentum = new TClonesArray("TVector3", MAX_PHOTONS);
  ana->Branch("pho_conv_refitted_momentum", "TClonesArray", &pho_conv_refitted_momentum, 32000, 0);
  pho_conv_vertexcorrected_p4 = new TClonesArray("TLorentzVector", MAX_PHOTONS);
  ana->Branch("pho_conv_vertexcorrected_p4", "TClonesArray", &pho_conv_vertexcorrected_p4, 32000, 0);
}

bool GlobePhotons::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  if (debug_level > 9) 
    std::cout << "GlobePhotons: Start analyze" << std::endl;

  //PhotonFix::initialiseGeometry(iSetup);
  checkSetup(iSetup);

  if (haveTowers) {
    hcalHelper->readEvent(const_cast<edm::Event &>(iEvent));
    hcalHelperPflow->readEvent(const_cast<edm::Event &>(iEvent));
  }

  // get collections
  edm::Handle<reco::PhotonCollection> phoH;
  iEvent.getByLabel(photonCollStd, phoH);

  edm::Handle<reco::PhotonCollection> phoHpf;
  iEvent.getByLabel(photonCollPf, phoHpf);

  // take the pi0 rejection info from RECO
  edm::Handle<reco::PhotonPi0DiscriminatorAssociationMap>  map;
  reco::PhotonPi0DiscriminatorAssociationMap::const_iterator mapIter;
  
  edm::Handle<reco::PhotonCollection> R_PhotonHandle;
  iEvent.getByLabel(photonCollStd, R_PhotonHandle);
  const reco::PhotonCollection R_photons = *(R_PhotonHandle.product());  
  
  //// Get the Out In CKF tracks from conversions 
  //bool validTrackInputs=true;
  edm::Handle<reco::TrackCollection> outInTrkHandle;

  if (!doFastSim) {
    iEvent.getByLabel("ckfOutInTracksFromConversions",  outInTrkHandle);
    
    if (!outInTrkHandle.isValid()) {
      std::cout << "Error! Can't get the conversionOITrack " << "\n";
      //validTrackInputs=false;
      if (debug_level > 9)
	std::cout  << "ConvertedPhotonProducer  outInTrack collection size " << (*outInTrkHandle).size() << "\n";
    }
  }

  //// Get the association map between CKF Out In tracks and the SC where they originated
  edm::Handle<reco::TrackCaloClusterPtrAssociation> outInTrkSCAssocHandle; 

  if (!doFastSim) {
    iEvent.getByLabel("ckfOutInTracksFromConversions" , "outInTrackSCAssociationCollection", outInTrkSCAssocHandle);
    if (!outInTrkSCAssocHandle.isValid()) {
      //  std::cout << "Error! Can't get the product " <<  outInTrackSCAssociationCollection_.c_str() <<"\n";
      //validTrackInputs=false;
    }
  }
  
  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByLabel(beamSpotColl, bsHandle);
  const reco::BeamSpot &thebs = *bsHandle.product();
  
  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel(convertedPhotonColl, hConversions);
  iEvent.getByLabel(electronColl, hElectrons);
  iEvent.getByLabel(vtxCollection, hVertex);
  iEvent.getByLabel(tkCollection, tkHandle);

  edm::Handle<double> rhoHandle;
  iEvent.getByLabel(rhoCollection, rhoHandle);
  rho = *(rhoHandle.product());
  
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  const CaloGeometry& geometry = *geoHandle;

  const CaloTopology *topology;
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  topology = pTopology.product();
  
  const CaloSubdetectorGeometry *geometryES = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
  CaloSubdetectorTopology *topology_p = 0;
  if (geometryES) 
    topology_p = new EcalPreshowerTopology(geoHandle);
     
  // FOR PF ISOLATION
  edm::Handle<reco::PFCandidateCollection> pfCollection;
  iEvent.getByLabel(pfColl, pfCollection);

  if (debug_level > 9) {
    std::cout << "GlobePhotons: Start analyze" << std::endl;
  }

  edm::Handle<reco::SuperClusterCollection> superClustersHybridH; 
  edm::Handle<reco::SuperClusterCollection> superClustersEndcapH; 
  edm::Handle<reco::BasicClusterShapeAssociationCollection> hybridClusterShapeBarrelH; 
  edm::Handle<reco::BasicClusterShapeAssociationCollection> basicClusterShapeEndcapH; 

  iEvent.getByLabel(barrelHybridClusterShapeColl, hybridClusterShapeBarrelH);
  iEvent.getByLabel(endcapBasicClusterShapeColl, basicClusterShapeEndcapH);

  iEvent.getByLabel(hybridSuperClusterColl, superClustersHybridH);
  iEvent.getByLabel(endcapSuperClusterColl, superClustersEndcapH);

  // pass the collection to the ID calculator
  cicPhotonId->configure(hVertex, tkHandle, hElectrons, pfCollection, rho); 

  if (debug_level > 9) {
    std::cout << "GlobePhotons: Photon collection size: "<< phoH->size() << std::endl;

    std::cout << "GlobePhotons: superClustersEndcap collection size: "<< superClustersEndcapH->size() << std::endl;
  }

  // now have collections
  pho_p4->Clear();
  pho_calopos->Clear();
  pho_conv_vtx->Clear();
  pho_conv_pair_momentum->Clear();
  pho_conv_refitted_momentum->Clear();
  pho_conv_vertexcorrected_p4->Clear();
  pho_pfiso_mycharged01->clear();
  pho_pfiso_mycharged02->clear();
  pho_pfiso_mycharged03->clear();
  pho_pfiso_mycharged04->clear();
  pho_pfiso_mycharged05->clear();
  pho_pfiso_mycharged06->clear();

  pho_schits->clear();
  pho_bchits->clear();

  pho_n = 0;

  if(debug_level>9)std::cout << "GlobePhotons: photons" << std::endl;

  for(unsigned int iPho = 0; iPho < phoH->size(); ++iPho) {
    if (pho_n >= MAX_PHOTONS) {
      std::cout << "GlobePhotons: WARNING TOO MANY PHOTONS: " << phoH->size() << " (allowed " << MAX_PHOTONS << ")" << std::endl;
      break;
    }

    reco::PhotonRef localPho(phoH, iPho);
    
    if(gCUT->cut(*localPho)) 
      continue;

    new ((*pho_p4)[pho_n]) TLorentzVector();
    new ((*pho_calopos)[pho_n]) TVector3();
    new ((*pho_conv_vtx)[pho_n]) TVector3();
    new ((*pho_conv_pair_momentum)[pho_n]) TVector3();
    new ((*pho_conv_refitted_momentum)[pho_n]) TVector3();
    new ((*pho_conv_vertexcorrected_p4)[pho_n]) TLorentzVector();

    if(debug_level>9)
      std::cout << "GlobePhotons: -21 "<< std::endl;

    ((TLorentzVector *)pho_p4->At(pho_n))->SetXYZT(localPho->px(), localPho->py(), localPho->pz(), localPho->energy());
    ((TVector3 *)pho_calopos->At(pho_n))->SetXYZ(localPho->caloPosition().x(), localPho->caloPosition().y(), localPho->caloPosition().z());

    reco::SuperClusterRef theClus = localPho->superCluster();

    std::vector<UInt_t> schits,bchits;
    pho_scind[pho_n] = -1;

    //PF info
    pho_isPFPhoton[pho_n] = 0;
    pho_isPFElectron[pho_n] = 0;

    pho_residCorrEnergy[pho_n] = 0;//ResidCorrector.fixedEnergy();
    pho_residCorrResn[pho_n] = 0;//ResidCorrector.sigmaEnergy();

    EcalClusterLazyTools lazyTool(iEvent, iSetup, ecalHitEBColl, ecalHitEEColl);   
    EcalClusterLazyTools lazyToolCleaned(iEvent, iSetup, ecalHitEBCleanedColl, ecalHitEECleanedColl);   
   
    pho_sc_time[pho_n] = lazyTool.SuperClusterTime(*theClus, iEvent);
	
    int index = 0;
    for(int isuperClusterType=0; isuperClusterType<3; ++isuperClusterType) {
      if (isuperClusterType == 0) {
        for(reco::SuperClusterCollection::size_type j = 0; j<superClustersHybridH->size(); ++j){

          reco::SuperClusterRef sc(superClustersHybridH, j);

          //apply cuts
          if(gCUT->cut(*sc))
	    continue;

          //passed cuts
          if (&(*localPho->superCluster()) == &(*sc)) {
            pho_scind[pho_n] = index;
            break;
          }
          index++;
        }
      }

      if (isuperClusterType == 2) {
        for(reco::SuperClusterCollection::size_type j = 0; j<superClustersEndcapH->size(); ++j){

          reco::SuperClusterRef sc(superClustersEndcapH, j);
          //apply cuts
          if(gCUT->cut(*sc))
	    continue;

          //passed cuts
          if (&(*(localPho->superCluster())) == &(*sc)) {
            pho_scind[pho_n] = index;
            break;
          }
          index++;
        }
      }
    }

    DetId id=localPho->superCluster()->seed()->hitsAndFractions()[0].first;

    bool isBarrel=(id.subdetId() == EcalBarrel);
    pho_barrel[pho_n]=(Int_t)isBarrel;

    // Rech-Hits related
     edm::Handle<EcalRecHitCollection> prechits;
    iEvent.getByLabel( (localPho->isEB() ? ecalHitEBColl : ecalHitEEColl) ,prechits);

    edm::Handle<EcalRecHitCollection> ESRecHits;
    iEvent.getByLabel(ecalHitESColl , ESRecHits);

    const reco::CaloClusterPtr  seed_clu = localPho->superCluster()->seed();
    EcalRecHitCollection::const_iterator seedcry_rh = prechits->find(id);

    //fiducial flags
    pho_isEB[pho_n] = localPho->isEB();
    pho_isEE[pho_n] = localPho->isEE();
    pho_isEBGap[pho_n] = localPho->isEBGap();
    pho_isEEGap[pho_n] = localPho->isEEGap();
    pho_isEBEEGap[pho_n] = localPho->isEBEEGap();
    pho_isEBEtaGap[pho_n] = localPho->isEBEtaGap();
    pho_isEBPhiGap[pho_n] = localPho->isEBPhiGap();
    pho_isEEDeeGap[pho_n] = localPho->isEEDeeGap();
    pho_isEERingGap[pho_n] = localPho->isEERingGap();
    
    

    //shower shape variables
    pho_see[pho_n] = localPho->sigmaEtaEta();
    std::vector<float> etaphimomentsnoZS = myClusterTools::noZSEcalClusterTools::localCovariances(*localPho->superCluster()->seed(), &(*prechits), &(*topology));
    pho_sieie[pho_n] = sqrt(etaphimomentsnoZS[0]);
    pho_sieip[pho_n] = etaphimomentsnoZS[1];
    pho_sipip[pho_n] = sqrt(etaphimomentsnoZS[2]);
    pho_e3x3[pho_n] = myClusterTools::noZSEcalClusterTools::e3x3(*localPho->superCluster()->seed(), &(*prechits), &(*topology));
    pho_e5x5[pho_n] = myClusterTools::noZSEcalClusterTools::e5x5(*localPho->superCluster()->seed(), &(*prechits), &(*topology));
    pho_r9[pho_n] =  pho_e3x3[pho_n]/localPho->superCluster()->rawEnergy();
    myClusterTools::Cluster2ndMoments momentsNoZS = myClusterTools::noZSEcalClusterTools::cluster2ndMoments(*localPho->superCluster()->seed(), *prechits);
    pho_smaj[pho_n] = momentsNoZS.sMaj;
    pho_smin[pho_n] = momentsNoZS.sMin;
    pho_alpha[pho_n] = momentsNoZS.alpha;
    
    std::vector<float> etaphimoments = myClusterTools::EcalClusterTools::localCovariances(*localPho->superCluster()->seed(), &(*prechits), &(*topology));
    pho_sieie_frac[pho_n] = sqrt(etaphimoments[0]);
    pho_sieip_frac[pho_n] = etaphimoments[1];
    pho_sipip_frac[pho_n] = sqrt(etaphimoments[2]);
    pho_e3x3_frac[pho_n] = myClusterTools::EcalClusterTools::e3x3(*localPho->superCluster()->seed(), &(*prechits), &(*topology));
    pho_e5x5_frac[pho_n] = myClusterTools::EcalClusterTools::e5x5(*localPho->superCluster()->seed(), &(*prechits), &(*topology));
    pho_r9_frac[pho_n] =  pho_e3x3_frac[pho_n]/localPho->superCluster()->rawEnergy();
    myClusterTools::Cluster2ndMoments moments = myClusterTools::EcalClusterTools::cluster2ndMoments(*localPho->superCluster()->seed(), *prechits);
    pho_smaj_frac[pho_n] = moments.sMaj;
    pho_smin_frac[pho_n] = moments.sMin;
    pho_alpha_frac[pho_n] = moments.alpha;
        

    pho_r9_cleaned[pho_n] = lazyToolCleaned.e3x3(*seed_clu)/localPho->superCluster()->rawEnergy();
    pho_e1x5[pho_n] = localPho->e1x5();
    pho_e2x2[pho_n] = lazyTool.e2x2(*seed_clu);
    pho_e1x3[pho_n] = lazyTool.e1x3(*seed_clu);
    pho_e2x2_cleaned[pho_n] = lazyToolCleaned.e2x2(*seed_clu);
    pho_e5x5_cleaned[pho_n] = lazyToolCleaned.e5x5(*seed_clu);
    pho_sieie_cleaned[pho_n] = sqrt(lazyToolCleaned.localCovariances(*seed_clu)[0]);
    pho_sieip_cleaned[pho_n] = lazyToolCleaned.localCovariances(*seed_clu)[1];
    pho_emaxxtal[pho_n] = localPho->maxEnergyXtal();
    pho_hoe[pho_n] = localPho->hadronicOverEm();
    pho_h1oe[pho_n] = localPho->hadronicDepth1OverEm();
    pho_h2oe[pho_n] = localPho->hadronicDepth2OverEm();
    
    pho_h1oe_bc[pho_n] = 0;//hcalHelper->hcalESumDepth1BehindClusters(caloTwId)/localPho->superCluster()->energy();
    pho_h2oe_bc[pho_n] = 0;//hcalHelper->hcalESumDepth2BehindClusters(caloTwId)/localPho->superCluster()->energy();
    pho_hoe_bc[pho_n]  = 0;//pho_h1oe_bc[pho_n] + pho_h2oe_bc[pho_n];
    if (haveTowers) {
      if (!(localPho->isStandardPhoton())) {
	std::vector<CaloTowerDetId> caloTwId = hcalHelperPflow->hcalTowersBehindClusters(*(localPho->superCluster()));
	pho_h1oe_bc[pho_n] = hcalHelperPflow->hcalESumDepth1BehindClusters(caloTwId)/localPho->superCluster()->energy();
	pho_h2oe_bc[pho_n] = hcalHelperPflow->hcalESumDepth2BehindClusters(caloTwId)/localPho->superCluster()->energy();
	pho_hoe_bc[pho_n]  = pho_h1oe_bc[pho_n] + pho_h2oe_bc[pho_n];
      } else {
	std::vector<CaloTowerDetId> caloTwId = hcalHelper->hcalTowersBehindClusters(*(localPho->superCluster()));
	pho_h1oe_bc[pho_n] = hcalHelper->hcalESumDepth1BehindClusters(caloTwId)/localPho->superCluster()->energy();
	pho_h2oe_bc[pho_n] = hcalHelper->hcalESumDepth2BehindClusters(caloTwId)/localPho->superCluster()->energy();
	pho_hoe_bc[pho_n]  = pho_h1oe_bc[pho_n] + pho_h2oe_bc[pho_n];
      }
    }

    pho_r1x5[pho_n] = localPho->r1x5();
    pho_r2x5[pho_n] = localPho->r2x5();

    pho_must[pho_n] = -9999.;
    pho_mustnc[pho_n] = -1;
    //reco::Mustache m;
    //m.MustacheID(*(localPho->superCluster()), pho_mustnc[pho_n], pho_must[pho_n]);

    edm::Handle<EcalRecHitCollection> EBReducedRecHits;
    edm::Handle<EcalRecHitCollection> EEReducedRecHits;
    iEvent.getByLabel(ecalHitEBColl, EBReducedRecHits);
    iEvent.getByLabel(ecalHitEEColl, EEReducedRecHits);

    //ggPFPhotons ggPFPhoton(*localPho, phoHpf,hElectrons,
    //			   pfCollection,
    //			   EBReducedRecHits,
    //			   EEReducedRecHits,
    //			   ESRecHits,
    //			   geomBar_,
    //			   geomEnd_,
    //			   bsHandle
    //  			   );
    
    pho_pfconvVtxZ[pho_n] = -9999.;
    pho_pfconvVtxZErr[pho_n] = -9999.;
    pho_hasConvPf[pho_n] = -9999.;
    pho_hasSLConvPf[pho_n] = -9999.;
    pho_pfpresh1[pho_n] = -9999.;
    pho_pfpresh2[pho_n] = -9999.;
    pho_mustenergy[pho_n] = -9999.;
    pho_mustenergyout[pho_n] = -9999.;
    pho_pflowE[pho_n] = -9999.;
    pho_pfdeta[pho_n] = -9999.;
    pho_pfdphi[pho_n] = -9999.;
    pho_pfclusrms[pho_n] = -9999.;
    pho_pfclusrmsmust[pho_n] = -9999.;

    pho_pfsieie[pho_n] = -9999.;
    pho_pfsieip[pho_n] = -9999.;
    pho_pfsipip[pho_n] = -9999.;
    pho_pfe2x2[pho_n] = -9999.;
    pho_pfe3x3[pho_n] = -9999.;
    pho_pfe5x5[pho_n] =-9999.;
    pho_pfemaxxtal[pho_n] = -9999.;
    pho_pfe2nd[pho_n] = -9999.;
    
    pho_pfRawEnergy[pho_n] = 0;

    // more cluster shapes from Lazy Tools
    std::vector<float> viCov;
    viCov = lazyTool.localCovariances(*seed_clu);
    std::vector<float> cov;
    cov = lazyTool.covariances(*localPho->superCluster());
    
    //pho_sipip[pho_n] = viCov[2];
    //pho_sieip[pho_n] = viCov[1];
    pho_zernike20[pho_n] = lazyTool.zernike20(*seed_clu);
    pho_zernike42[pho_n] = lazyTool.zernike42(*seed_clu);
    pho_e2nd[pho_n] = lazyTool.e2nd(*seed_clu);
    pho_e2x5right[pho_n] = lazyTool.e2x5Right(*seed_clu);
    pho_e2x5left[pho_n] = lazyTool.e2x5Left(*seed_clu);
    pho_e2x5top[pho_n] = lazyTool.e2x5Top(*seed_clu);
    pho_e2x5bottom[pho_n] = lazyTool.e2x5Bottom(*seed_clu);
    pho_e2x5max[pho_n] = lazyTool.e2x5Max(*seed_clu);
    pho_eright[pho_n] = lazyTool.eRight(*seed_clu);
    pho_eleft[pho_n] = lazyTool.eLeft(*seed_clu);
    pho_etop[pho_n] = lazyTool.eTop(*seed_clu);
    pho_ebottom[pho_n] = lazyTool.eBottom(*seed_clu);

    if (seed_clu->hitsAndFractions().at(0).first.subdetId() == EcalBarrel) {
      EcalClusterLocal _ecalLocal;
      _ecalLocal.localCoordsEB(*seed_clu, iSetup, pho_betacry[pho_n], pho_bphicry[pho_n], pho_bieta[pho_n], pho_biphi[pho_n], pho_bthetatilt[pho_n], pho_bphitilt[pho_n]);
    } else {
      pho_betacry[pho_n] = -999;
      pho_bphicry[pho_n] = -999;
      pho_bieta[pho_n] = -999;
      pho_biphi[pho_n] = -999;
      pho_bthetatilt[pho_n] = -999;
      pho_bphitilt[pho_n] = -999;
      
    }
    // NN variables
    pho_r19[pho_n]            = lazyTool.eMax(*seed_clu)/pho_e3x3[pho_n];
    pho_maxoraw[pho_n]        = lazyTool.eMax(*seed_clu)/localPho->superCluster()->rawEnergy();
    pho_cep[pho_n]            = viCov[1];
    float lambdaMinus         = (viCov[0] + viCov[2] - sqrt(pow(viCov[0] - viCov[2], 2) + 4*pow(viCov[1], 2)));
    float lambdaPlus          = (viCov[0] + viCov[2] + sqrt(pow(viCov[0] - viCov[2], 2) + 4*pow(viCov[1], 2)));
    pho_lambdaratio[pho_n]    = lambdaMinus/lambdaPlus;
    pho_lambdadivcov[pho_n]   = lambdaMinus/viCov[0];
    pho_etawidth[pho_n]       = localPho->superCluster()->etaWidth();
    pho_brem[pho_n]           = localPho->superCluster()->phiWidth()/localPho->superCluster()->etaWidth();
    
    // NN variables
    pho_cep_global[pho_n]          = cov[1];
    float lambdaMinus_glo          = (cov[0] + cov[2] - sqrt(pow(cov[0] - cov[2], 2) + 4*pow(cov[1], 2)));
    float lambdaPlus_glo           = (cov[0] + cov[2] + sqrt(pow(cov[0] - cov[2], 2) + 4*pow(cov[1], 2)));
    pho_lambdaratio_global[pho_n]  = lambdaMinus_glo/lambdaPlus_glo;
    pho_lambdadivcov_global[pho_n] = lambdaMinus_glo/cov[0];

    //spike-ID
    //pho_e2overe9[pho_n] = EcalSeverityLevelAlgo::E2overE9( id, *prechits, 5.0, 0.0);
    pho_seed_severity[pho_n] = 0;//sevLevel->severityLevel(id, *prechits);

    pho_seed_time[pho_n] = seedcry_rh != prechits->end() ? seedcry_rh->time() : 999.;
    pho_seed_outoftimechi2[pho_n] = seedcry_rh != prechits->end() ? seedcry_rh->outOfTimeChi2() : 999.;
    pho_seed_chi2[pho_n] = seedcry_rh != prechits->end() ? seedcry_rh->chi2() : 999.;
    pho_seed_recoflag[pho_n] = seedcry_rh != prechits->end() ? seedcry_rh->recoFlag() : 999.;

    pho_ecalsumetconedr04[pho_n] = localPho->ecalRecHitSumEtConeDR04();
    pho_hcalsumetconedr04[pho_n] = localPho->hcalTowerSumEtConeDR04();
    pho_hcal1sumetconedr04[pho_n] = localPho->hcalDepth1TowerSumEtConeDR04();
    pho_hcal2sumetconedr04[pho_n] = localPho->hcalDepth2TowerSumEtConeDR04();
    pho_trksumptsolidconedr04[pho_n] = localPho->trkSumPtSolidConeDR04();
    pho_trksumpthollowconedr04[pho_n] = localPho->trkSumPtHollowConeDR04();
    pho_ntrksolidconedr04[pho_n] = localPho->nTrkSolidConeDR04();
    pho_ntrkhollowconedr04[pho_n] = localPho->nTrkHollowConeDR04();
    pho_ecalsumetconedr03[pho_n] = localPho->ecalRecHitSumEtConeDR03();
    pho_hcalsumetconedr03[pho_n] = localPho->hcalTowerSumEtConeDR03();
    pho_hcal1sumetconedr03[pho_n] = localPho->hcalDepth1TowerSumEtConeDR03();
    pho_hcal2sumetconedr03[pho_n] = localPho->hcalDepth2TowerSumEtConeDR03();
    pho_trksumptsolidconedr03[pho_n] = localPho->trkSumPtSolidConeDR03();
    pho_trksumpthollowconedr03[pho_n] = localPho->trkSumPtHollowConeDR03();
    pho_ntrksolidconedr03[pho_n] = localPho->nTrkSolidConeDR03();
    pho_ntrkhollowconedr03[pho_n] = localPho->nTrkHollowConeDR03();

    bool passelectronveto = !ConversionTools::hasMatchedPromptElectron(localPho->superCluster(), hElectrons, hConversions, thebs.position());
    pho_isconv[pho_n] = int(passelectronveto);

    // ES variables
    pho_eseffsixix[pho_n] = 0.;
    pho_eseffsiyiy[pho_n] = 0.;
    if (ESRecHits.isValid() && (fabs(localPho->superCluster()->eta()) > 1.6 && fabs(localPho->superCluster()->eta()) < 3)) {
      std::vector<float> phoESHits0 = gES->getESHits(localPho->superCluster()->x(), localPho->superCluster()->y(), localPho->superCluster()->z(), rechits_map_, geometry, topology_p, 0);
      std::vector<float> phoESShape = gES->getESShape(phoESHits0);
      pho_eseffsixix[pho_n] = phoESShape[0];
      pho_eseffsiyiy[pho_n] = phoESShape[1];
    }

    //other variables
    pho_haspixseed[pho_n] = localPho->hasPixelSeed();
    pho_hasconvtks[pho_n] = localPho->hasConversionTracks();
    pho_nconv[pho_n] = localPho->conversions().size();
    pho_conv_ntracks[pho_n]=0;
    pho_conv_pairinvmass[pho_n]=-999.;
    pho_conv_paircotthetasep[pho_n]=-999.;
    pho_conv_eoverp[pho_n]=-999.;
    pho_conv_zofprimvtxfromtrks[pho_n]=-999.;
    pho_conv_distofminapproach[pho_n]=-999.;
    pho_conv_dphitrksatvtx[pho_n]=-999.;
    pho_conv_dphitrksatecal[pho_n]=-999.;
    pho_conv_detatrksatecal[pho_n]=-999.;
    pho_conv_tk1_d0[pho_n]=-999.;
    pho_conv_tk1_pout[pho_n]=-999.;
    pho_conv_tk1_pin[pho_n]=-999.;
    pho_conv_validvtx[pho_n]=0;

    ((TVector3 *)pho_conv_vtx->At(pho_n))->SetXYZ(-999, -999, -999);
    ((TVector3 *)pho_conv_pair_momentum->At(pho_n))->SetXYZ(-999, -999, -999);
    ((TVector3 *)pho_conv_refitted_momentum->At(pho_n))->SetXYZ(-999, -999, -999);

    pho_n++;
  }

  if(debug_level>9)
    std::cout << "End Photon" << std::endl;

  delete topology_p;

  return true;
}
