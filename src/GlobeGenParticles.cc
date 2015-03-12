#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeGenParticles.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"

#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"
#include "SimDataFormats/JetMatching/interface/MatchedPartons.h"
#include "SimDataFormats/JetMatching/interface/JetMatchedPartons.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <iostream>

GlobeGenParticles::GlobeGenParticles(const edm::ParameterSet& iConfig) {
  
  genParticlesColl = iConfig.getParameter<edm::InputTag>("GenParticlesColl");
  jetFlavorColl = iConfig.getParameter<edm::InputTag>("JetFlavorColl");
  photonColl = iConfig.getParameter<edm::InputTag>("PhotonCollStd");
  debug_level = iConfig.getParameter<int>("Debug_Level");
  gCUT = new GlobeCuts(iConfig);
}

void GlobeGenParticles::defineBranch(GlobeAnalyzer* ana) {

  // think about changing branch names for duplicate collections
  gp_p4 = new TClonesArray("TLorentzVector", MAX_GENERATOR);
  gp_vtx = new TClonesArray("TVector3", MAX_GENERATOR);
  
  ana->Branch("gp_n", &gp_n, "gp_n/I");
  
  ana->Branch("gp_p4", "TClonesArray", &gp_p4, 32000, 0);
  ana->Branch("gp_vtx", "TClonesArray", &gp_vtx, 32000, 0);
  
  ana->Branch("gp_status", gp_status, "gp_status[gp_n]/S");
  ana->Branch("gp_pdgid", gp_pdgid, "gp_pdgid[gp_n]/S");
  ana->Branch("gp_mother", gp_mother, "gp_mother[gp_n]/S");

  ana->Branch("pho_flavor", flavor, "pho_flavor[pho_n]/F");
}

bool GlobeGenParticles::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // take collections
  edm::Handle<reco::GenParticleCollection> gpH;
  iEvent.getByLabel(genParticlesColl, gpH);

  gp_p4->Clear();
  gp_vtx->Clear();
  gp_n = 0;
  
  if(debug_level > 999)
    std::cout << "#"
	      << "\t pdgid"
	      << "\t status"
	      << "\t mother"
	      << "\t px \t py \t pz \t E"
	      << "\t vx \t vy \t vz"
	      << std::endl;
  
  
  for(size_t i = 0; i < gpH->size(); ++i) {

    const reco::GenParticleRef gp(gpH, i);

    if(!gCUT->cut(*gp))
      continue;

    gp_pdgid[gp_n] = gp->pdgId();
    gp_status[gp_n] = gp->status();
    if (gp->numberOfMothers() != 0)
      gp_mother[gp_n] = gp->motherRef().key();
    else
      gp_mother[gp_n] = -1;

    new ((*gp_p4)[gp_n]) TLorentzVector();
    ((TLorentzVector *)gp_p4->At(gp_n))->SetPtEtaPhiM(gp->pt(), gp->eta(), gp->phi(), gp->mass());

    new ((*gp_vtx)[gp_n]) TVector3();
    ((TVector3 *)gp_vtx->At(gp_n))->SetXYZ(gp->vx(), gp->vy(), gp->vz());
    
    if(debug_level > 999)
      std::cout << gp_n	<< "\t" << gp_pdgid[ gp_n]
                << "\t" << gp_status[gp_n] << "\t" << gp_mother[gp_n]
                << "\t" << gp->px() << "\t" << gp->py() 
                << "\t" << gp->pz() << "\t" << gp->energy()
                << "\t" << gp->vx() << "\t" << gp->vy() << "\t" << gp->vz()
                << std::endl;
    
    gp_n++;
  }

  edm::Handle<reco::PhotonCollection> fakePhotonH;
  iEvent.getByLabel(photonColl, fakePhotonH);
  
  edm::Handle<reco::JetFlavourMatchingCollection> theTagByValue;
  iEvent.getByLabel (jetFlavorColl, theTagByValue);
  
  for (unsigned int p=0; p<fakePhotonH->size(); p++) {
    reco::PhotonRef fakePhoton(fakePhotonH, p);
    flavor[p] = -1;

    std::vector< const reco::Candidate * > constituents;	
    for (reco::JetFlavourMatchingCollection::const_iterator j  = theTagByValue->begin(); j != theTagByValue->end(); j++) {
      edm::RefToBase<reco::Jet> aJet  = (*j).first; //jet  
      const reco::JetFlavour aFlav = (*j).second; //flavour
      constituents = aJet->getJetConstituentsQuick();
      for(unsigned int i=0;i<constituents.size();i++) {
	float dr = reco::deltaR(constituents.at(i)->eta(), constituents.at(i)->phi(), fakePhoton->eta(), fakePhoton->phi());
	if (dr < 0.05) {
	  //the fake photon is contained in this jet
	  if(aFlav.getFlavour()==21) 
	    flavor[p] = 1;
	  if(fabs(aFlav.getFlavour())<6) 
	    flavor[p] = 2;
	  break;
	}
      }
    }
  }

  return true;
}
