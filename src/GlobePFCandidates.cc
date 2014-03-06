#include "HiggsAnalysis/HiggsTo2photons/interface/GlobePFCandidates.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "RecoParticleFlow/PFClusterTools/interface/ClusterClusterMapping.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidateFwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidateFwd.h"

#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>

GlobePFCandidates::GlobePFCandidates(const edm::ParameterSet& iConfig) {
  
  pfColl = iConfig.getParameter<edm::InputTag>("PFCandidateColl");
  electronCollStd = iConfig.getParameter<edm::InputTag>("ElectronColl");
  photonCollStd =  iConfig.getParameter<edm::InputTag>("PhotonCollStd");
  PFIsoOuterConeSize = iConfig.getParameter<double>("PFIsoOuterCone");
  debug_level = iConfig.getParameter<int>("Debug_Level");
  
  // get cut thresholds
  gCUT = new GlobeCuts(iConfig);

  pfcand_hitind = new std::vector<std::vector<short>>;
}

void GlobePFCandidates::defineBranch(GlobeAnalyzer* ana) {
  
  pfcand_p4 = new TClonesArray("TLorentzVector", MAX_PFCANDS);
  pfcand_poscalo = new TClonesArray("TVector3", MAX_PFCANDS);
  pfcand_posvtx = new TClonesArray("TVector3", MAX_PFCANDS);
  
  ana->Branch("pfcand_p4", "TClonesArray", &pfcand_p4, 32000, 0);
  ana->Branch("pfcand_poscalo", "TClonesArray", &pfcand_poscalo, 32000, 0);
  ana->Branch("pfcand_posvtx", "TClonesArray", &pfcand_posvtx, 32000, 0);

  ana->Branch("pfcand_n", &pfcand_n, "pfcand_n/I");
  ana->Branch("pfcand_pdgid",&pfcand_pdgid,"pfcand_pdgid[pfcand_n]/I");
  ana->Branch("pfcand_ecalenergy", &pfcand_ecalEnergy, "pfcand_ecalenergy[pfcand_n]/F");
  ana->Branch("pfcand_hcalenergy", &pfcand_hcalEnergy, "pfcand_hcalenergy[pfcand_n]/F");
  ana->Branch("pfcand_rawecalenergy", &pfcand_rawEcalEnergy, "pfcand_rawecalenergy[pfcand_n]/F");
  ana->Branch("pfcand_rawhcalenergy", &pfcand_rawHcalEnergy, "pfcand_rawhcalenergy[pfcand_n]/F");
  ana->Branch("pfcand_ps1energy", &pfcand_ps1Energy, "pfcand_ps1energy[pfcand_n]/F");
  ana->Branch("pfcand_ps2energy", &pfcand_ps2Energy, "pfcand_ps2energy[pfcand_n]/F");
  ana->Branch("pfcand_momerr", &pfcand_momErr, "pfcand_momerr[pfcand_n]/F");
  //ana->Branch("pfcand_mva_e_pi", &pfcand_mva_e_pi, "pfcand_mva_e_pi[pfcand_n]/F");
  //ana->Branch("pfcand_mva_e_mu", &pfcand_mva_e_mu, "pfcand_mva_e_mu[pfcand_n]/F");
  //ana->Branch("pfcand_mva_pi_mu", &pfcand_mva_pi_mu, "pfcand_mva_pi_mu[pfcand_n]/F");
  //ana->Branch("pfcand_mva_nothing_gamma", &pfcand_mva_nothing_gamma, "pfcand_mva_nothing_gamma[pfcand_n]/F");
  //ana->Branch("pfcand_mva_nothing_nh", &pfcand_mva_nothing_nh, "pfcand_mva_nothing_nh[pfcand_n]/F");
  //ana->Branch("pfcand_mva_gamma_nh", &pfcand_mva_gamma_nh, "pfcand_mva_gamma_nh[pfcand_n]/F");
  //ana->Branch("pfcand_vz",&pfcand_vz,"pfcand_vz[pfcand_n]/F");
  //ana->Branch("pfcand_overlappho",&pfcand_overlappho,"pfcand_overlappho[pfcand_n]/i");
  ana->Branch("pfcand_ispu",&pfcand_ispu,"pfcand_ispu[pfcand_n]/i");
  ana->Branch("pfcand_hitind", "std::vector<std::vector<short> > ", &pfcand_hitind);
}

bool GlobePFCandidates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup, GlobeTracks* tks, GlobeMuons* mus, GlobePhotons* phos, GlobeEcalHits* hits) {

  pfcand_hitind->clear();
  pfcand_p4->Clear(); 
  pfcand_poscalo->Clear(); 
  pfcand_posvtx->Clear(); 
  pfcand_n = 0;
  //pfcandtimespho_n = 0;

  // All PF Candidate
  edm::Handle<reco::PFCandidateCollection> pfCandidatesH;
  iEvent.getByLabel(pfColl, pfCandidatesH);
  std::vector<reco::PFCandidate> candidates = (*pfCandidatesH.product());

  //PF Candidates from PileUp
  edm::Handle<std::vector< edm::FwdPtr<reco::PFCandidate> >> pfCandidatesPileUpH;
  iEvent.getByLabel("pfPileUp", pfCandidatesPileUpH);
  std::vector< edm::FwdPtr<reco::PFCandidate> > pucandidates = (*pfCandidatesPileUpH.product());

  edm::Handle < reco::GsfElectronCollection > theEGammaCollection;
  iEvent.getByLabel(electronCollStd, theEGammaCollection);

  edm::Handle<reco::PhotonCollection> phoH;
  iEvent.getByLabel(photonCollStd, phoH);

  for (std::vector<reco::PFCandidate>::iterator it = candidates.begin(); it != candidates.end(); ++it) {

    if (pfcand_n >= MAX_PFCANDS) {
      std::cout << "GlobePFCandidates: WARNING TOO MANY PFCandidates: " << pfcand_n << " (allowed " << MAX_PFCANDS << ")" << std::endl;
      break;
    }

    bool save = false;  
    bool isCandFromPU = false;
    unsigned npucandidates = pucandidates.size();
    for (unsigned npuit = 0; npuit < npucandidates; npuit++) {
      float dR = deltaR(pucandidates[npuit]->p4(), it->p4());
      if (dR < 1e-5) {
	isCandFromPU = true;
	break;
      }
    }
    
    for (unsigned int iele = 0; iele < theEGammaCollection->size(); iele++) {
      reco::GsfElectronRef electronRef(theEGammaCollection, iele);
      // do not consider the same gsf-pf electron in the iso cone. 
      if(it->particleId() == 2) {
	reco::GsfTrackRef egGsfTrackRef = electronRef->gsfTrack();
	reco::GsfTrackRef pfGsfTrackRef = (*it).gsfTrackRef();
	if (egGsfTrackRef == pfGsfTrackRef)
	  continue;
      }

      double dphi = deltaPhi( it->phi(), electronRef->phi() );
      double deta = it->eta() - electronRef->eta();
      double dR = sqrt(deta * deta + dphi * dphi);
      if (dR < 0.6) {
	if (it->particleId() == 1 ||//charged hadrons
	    it->particleId() == 2 ||// electrons do you want them?
	    it->particleId() == 3 ||// muons do you want them ?
	    it->particleId() == 4 ||// photons
	    it->particleId() == 5) {// netrual hadrons

	  save = true;
	} else {
	  save = false;
	}
      }
    }


    for (unsigned int ipho = 0; ipho < phoH->size(); ipho++) {
      reco::PhotonRef photonRef(phoH, ipho);

      // do not consider the same gsf-pf electron in the iso cone. 
      if(it->particleId() == 4) {
	reco::PhotonRef pfPhotonRef = (*it).photonRef();
	if (photonRef == pfPhotonRef)
	  continue;
      }

      double dphi = deltaPhi( it->phi(), photonRef->phi() );
      double deta = it->eta() - photonRef->eta();
      double dR = sqrt(deta * deta + dphi * dphi);
      if (dR < 0.6) {
	if (it->particleId() == 1 ||//charged hadrons
	    it->particleId() == 2 ||// electrons do you want them?
	    it->particleId() == 3 ||// muons do you want them ?
	    it->particleId() == 4 ||// photons
	    it->particleId() == 5) {// netrual hadrons

	  save = true;
	} else {
	  save = false;
	}
      }
    }

    if (save) {
      std::vector<short> ecalhits;
      reco::SuperClusterRef scRef = it->superClusterRef();
      std::vector< std::pair<DetId, float> >::const_iterator itSC;
      if (scRef.isNonnull()) {
	for (int i=0; i < hits->ecalhit_n; i++) {
	  for (itSC=scRef->hitsAndFractions().begin(); itSC != scRef->hitsAndFractions().end(); ++itSC) {
	    if (hits->ecalhit_detid[i] == itSC->first.rawId()) {
	      ecalhits.push_back(i);
	      break;
	    }
	  }
	}
      }
      pfcand_hitind->push_back(ecalhits);

      pfcand_pdgid[pfcand_n] = it->particleId();
      pfcand_ecalEnergy[pfcand_n] = it->ecalEnergy();
      pfcand_hcalEnergy[pfcand_n] = it->hcalEnergy();
      pfcand_rawEcalEnergy[pfcand_n] = it->rawEcalEnergy();
      pfcand_rawHcalEnergy[pfcand_n] = it->rawHcalEnergy();
      pfcand_ps1Energy[pfcand_n] = it->pS1Energy();
      pfcand_ps2Energy[pfcand_n] = it->pS2Energy();
      pfcand_momErr[pfcand_n] = it->deltaP();
      //pfcand_mva_e_pi[pfcand_n] = it->mva_e_pi();
      //pfcand_mva_e_mu[pfcand_n] = it->mva_e_mu();
      //pfcand_mva_pi_mu[pfcand_n] = it->mva_pi_mu();
      //pfcand_mva_nothing_gamma[pfcand_n] = it->mva_nothing_gamma();
      //pfcand_mva_nothing_nh[pfcand_n] = it->mva_nothing_nh();
      //pfcand_mva_gamma_nh[pfcand_n] = it->mva_gamma_nh();
      
      new ((*pfcand_p4)[pfcand_n]) TLorentzVector();
      ((TLorentzVector *)pfcand_p4->At(pfcand_n))->SetXYZT(it->px(), it->py(), it->pz(), it->energy());
      
      new ((*pfcand_poscalo)[pfcand_n]) TVector3();
      ((TVector3 *)pfcand_poscalo->At(pfcand_n))->SetXYZ(it->positionAtECALEntrance().x(), 
							 it->positionAtECALEntrance().y(), 
							 it->positionAtECALEntrance().z());
      
      new ((*pfcand_posvtx)[pfcand_n]) TVector3();
      ((TVector3 *)pfcand_posvtx->At(pfcand_n))->SetXYZ(it->vx(), 
							it->vy(), 
							it->vz());
      
      //pfcand_vz[pfcand_n] = 9999;
      //if(it->particleId() == 1 ) {
      //pfcand_vz[pfcand_n] = it->trackRef()->vz();
      //}

      pfcand_ispu[pfcand_n] = 0;
      if (isCandFromPU)
	pfcand_ispu[pfcand_n] = 1;    // 1 means candidates from PU, 0 means candidates from PV
      
      pfcand_n++;
    }      
  }

  /*
  // take collections
  edm::Handle<reco::PFCandidateCollection> pfH;
  iEvent.getByLabel(pfColl, pfH);

  */

  return true;
}

