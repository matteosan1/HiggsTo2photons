#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeCuts.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"

GlobeCuts::GlobeCuts(const edm::ParameterSet& iConfig) {

  edm::ParameterSet pset = iConfig.getParameter<edm::ParameterSet>("Cuts");
  if (pset.exists("PdgId"))
    pdgidCut_  = pset.getParameter<std::vector<int> >("PdgId");
  if (pset.exists("Keep"))
    keepCut_   = pset.getParameter<bool>("Keep");
  if (pset.exists("EtCut"))
    etCut_     = pset.getParameter<double>("EtCut");
  if (pset.exists("EnergyCut"))
    energyCut_ = pset.getParameter<double>("EnergyCut");
  if (pset.exists("EtaCut"))
    etaCut_    = pset.getParameter<double>("EtaCut");
  if (pset.exists("TvpCut"))
    tvpCut_    = pset.getParameter<double>("TvpCut");
  if (pset.exists("LvpCut"))
    lvpCut_    = pset.getParameter<double>("LvpCut");
  if (pset.exists("HBHEEnergyCut"))
    hbHEEnergyCut_   = pset.getParameter<double>("HBHEEnergyCut");
  if (pset.exists("HFEnergyCut"))
    hfEnergyCut_     = pset.getParameter<double>("HFEnergyCut");
  if (pset.exists("KeepOutsideCone"))
    keepOutsideCone_ = pset.getParameter<bool>("KeepOutsideCone");
  if (pset.exists("HOEnergyCut"))
    hoEnergyCut_     = pset.getParameter<double>("HOEnergyCut");
  if (pset.exists("BarrelEnergyCut"))
    barrelECut_   = pset.getParameter<double>("BarrelEnergyCut");
  if (pset.exists("EndcapEnergyCut"))
    endcapECut_   = pset.getParameter<double>("EndcapEnergyCut");
  if (pset.exists("PreEnergyCut"))
    preECut_      = pset.getParameter<double>("PreEnergyCut");
  if (pset.exists("MaxDR"))
    maxDR_           = pset.getParameter<double>("MaxDR");
  if (pset.exists("InnerCone"))
    innerCone_        = pset.getParameter<double>("InnerCone");
  if (pset.exists("OuterCone"))
    outerCone_        = pset.getParameter<double>("OuterCone");
}

// The Functions return "true" if the object should be cut
// GenParticle
bool GlobeCuts::cut(const reco::GenParticle &gp) {

  if (gp.et() < etCut_)
    return 0;
  
  if (pdgidCut_.size() == 0)
    return 1;
  
  for(unsigned int i=0; i<pdgidCut_.size(); i++)
    if ((abs(gp.pdgId()) == pdgidCut_[i] and keepCut_) or 
	(abs(gp.pdgId()) != pdgidCut_[i] and !keepCut_))
      return 1;
  
  return 0;
}

//PFCands
bool GlobeCuts::cut(const reco::PFCandidate &pf) {
  return 0;
}

// Photons
bool GlobeCuts::cut(const reco::Photon &photon) { 
  return (photon.et() < etCut_); 
}


// Conversions
bool GlobeCuts::cut(const reco::Conversion &conv) { 
  return (sqrt(conv.refittedPairMomentum().perp2()) < etCut_); 
}

// Electrons
bool GlobeCuts::cut(const reco::GsfElectron& electron) { 
  return (electron.et() < etCut_); 
}
// Super Clusters
bool GlobeCuts::cut(const reco::SuperCluster &supercluster) { 
  return (supercluster.energy()/cosh(supercluster.eta()) < etCut_); 
}
// Basic Clusters
bool GlobeCuts::cut(const reco::BasicCluster &basiccluster) { 
  return (basiccluster.energy()/cosh(basiccluster.eta()) < etCut_); 
}
// CaloTowers
bool GlobeCuts::cut(const CaloTower &calotower) { 
  return (calotower.et() < etCut_); 
}

// Hcal RecHits Barrel
bool GlobeCuts::cut(const HBHERecHit &hcalrechit, double dR) { 
  
  if( dR < maxDR_ ) 
    return false; //keep if inside cone

  if (!keepOutsideCone_) 
    return true; //dont keep if user doesnt want

  return (hcalrechit.energy() < hbHEEnergyCut_ );
}

// Hcal RecHits Forward
bool GlobeCuts::cut(const HFRecHit &hcalrechit) { 
  return (hcalrechit.energy() < hfEnergyCut_); 
}

// Hcal RecHits Outer
bool GlobeCuts::cut(const HORecHit &hcalrechit) { 
  return (hcalrechit.energy() < hoEnergyCut_); 
}

// Tracks
bool GlobeCuts::cut(const reco::Track &track) { 
  return (track.pt() < etCut_); 
}

// TrackingParticle
bool GlobeCuts::cut(const TrackingParticle &tp) { 

  bool pdgidCut = true;
  for(unsigned int i=0; i<pdgidCut_.size(); i++) {
    if (abs(tp.pdgId()) == pdgidCut_[i]) {
      pdgidCut = false;
      break;
    }
  }

  return ( tp.pt() < etCut_ || 
           fabs(tp.eta()) > etaCut_ || 
           fabs(tp.vertex().z()) > lvpCut_ || 
           fabs(tp.vertex().rho()) > tvpCut_ ||
           pdgidCut); 
}


// Sim Hits
bool GlobeCuts::cut(const PSimHit &simhit) { 
  return (simhit.energyLoss() < energyCut_); 
}
// Sim Tracks
bool GlobeCuts::cut(const SimTrack &simtrack) { 
  return (simtrack.momentum().E() < energyCut_ || simtrack.noVertex());
}
// Jets
bool GlobeCuts::cut(const reco::CaloJet &jet) { 
  return (jet.energy() < energyCut_); 
}

bool GlobeCuts::cut(const reco::PFJet &jet) { 
  return (jet.energy() < energyCut_); 
}

bool GlobeCuts::cut(const reco::BasicJet &jet) { 
  return (jet.energy() < energyCut_); 
}

// Muons
bool GlobeCuts::cut(const reco::Muon &muon) { 
  return (muon.pt() < etCut_); 
}
// Gen Jets
bool GlobeCuts::cut(const reco::GenJet &genjet) { 
  return (genjet.et() < etCut_); 
}

//Cut to determine if RecHit falls within the lepton cone
//const EcalRecHit here is just used as a function identifier
bool GlobeCuts::cut(const EcalRecHit &ecalhit, int type, double dR) { 

  if (type == 2)  // preshower
    return (ecalhit.energy() < preECut_); 
   
  if( dR < maxDR_ ) 
    return false; //keep if inside cone
  
  if (!keepOutsideCone_) 
    return true; //dont keep if user doesnt want stuff outside cone

  if (type == 0)  // barrel
    return (fabs(ecalhit.energy()) < barrelECut_); 
  if (type == 1)  // endcap
    return (fabs(ecalhit.energy()) < endcapECut_); 
  return false;
}

bool GlobeCuts::isocut(const reco::Track &tk, const reco::Track &lep, const reco::Vertex &vtx) {
    
    using ROOT::Math::VectorUtil::DeltaR;
    
    //Inside the veto cone
    if(DeltaR(tk.innerMomentum(),lep.innerMomentum()) < innerCone_) return false;

    //Outside the isolation cone
    if(DeltaR(tk.innerMomentum(),lep.innerMomentum()) > outerCone_) return false;

    double d0 = sqrt((tk.vertex() - vtx.position()).perp2());
    double dZ = (tk.vertex() - vtx.position()).z();

    if(tk.numberOfValidHits() >= 10)
      return (tk.pt() > 1 && d0 < 1 && dZ < 5);
    else if(tk.numberOfValidHits() >= 8 && tk.numberOfValidHits() <= 9)
      return (tk.pt() > 1 && d0 < 0.2 && dZ < 2.0 && (d0 / tk.d0Error() < 10) && (dZ / tk.dzError() < 10));
    else if(tk.numberOfValidHits() >= 5 && tk.numberOfValidHits() <= 7)
      return (tk.pt() > 1 && d0 < 0.04 && dZ < 0.5 && (d0 / tk.d0Error() < 7) && (dZ / tk.dzError() < 7));
    else if(tk.numberOfValidHits() >= 3 && tk.numberOfValidHits() <= 4)
      return (tk.pt() > 1 && d0 < 0.02 && dZ < 0.2 && (d0 / tk.d0Error() < 3) && (dZ / tk.dzError() < 3));
    else return false;
}

bool GlobeCuts::isocut(const reco::Track &tk, const reco::Track &lep) {
    
    using ROOT::Math::VectorUtil::DeltaR;
    
    //Inside the veto cone
    if(DeltaR(tk.innerMomentum(),lep.innerMomentum()) < innerCone_) return false;

    //Outside the isolation cone
    if(DeltaR(tk.innerMomentum(),lep.innerMomentum()) > outerCone_) return false;

    double d0 = sqrt((tk.vertex() - lep.vertex()).perp2());
    double dZ = (tk.vertex() - lep.vertex()).z();
    
    if(tk.numberOfValidHits() >= 10)
      return (tk.pt() > 1 && d0 < 1 && dZ < 5);
    else if(tk.numberOfValidHits() >= 8 && tk.numberOfValidHits() <= 9)
        return (tk.pt() > 1 && d0 < 0.2 && dZ < 2.0 && d0 / tk.d0Error() < 10 && dZ / tk.dzError() < 10);
    else if(tk.numberOfValidHits() >= 5 && tk.numberOfValidHits() <= 7)
        return (tk.pt() > 1 && d0 < 0.04 && dZ < 0.5 && d0 / tk.d0Error() < 7 && dZ / tk.dzError() < 7);
    else if(tk.numberOfValidHits() >= 3 && tk.numberOfValidHits() <= 4)
        return (tk.pt() > 1 && d0 < 0.02 && dZ < 0.2 && d0 / tk.d0Error() < 3 && dZ / tk.dzError() < 3);
    else return false;
}

    




