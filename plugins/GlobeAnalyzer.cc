S#include "GlobeAnalyzer.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/Tools.h"

#include <stdio.h>
#include <time.h>

#define TIMINGDEBUG 0

static string memory_usage() {
  ostringstream mem;
  //PP("hi");
  ifstream proc("/proc/self/status");
  string s;
  while(getline(proc, s), !proc.fail()) {
    if(s.substr(0, 6) == "VmSize") {
      mem << s;
      return mem.str();
    }
  }
  return mem.str();
}

double diffclock1(clock_t clock1,clock_t clock2) {
  double diffticks=clock1-clock2;
  double diffms=(diffticks*10)/CLOCKS_PER_SEC;
  return diffms;
}

GlobeAnalyzer::GlobeAnalyzer(const edm::ParameterSet& iConfig) {

  fileName = iConfig.getParameter<std::string>("RootFileName");
  jobmaker = iConfig.getParameter<std::string>("JobMaker");
  branchesToSkim = iConfig.getParameter<std::vector<std::string> >("branchesToSkim");

  ExecCommand ex("");
  ex.setResults(iConfig.getParameter<std::string>("h2gAnalyzerVersion"));
  version = ex.getTag();

  globalCountersNames = iConfig.getParameter<std::vector<std::string> >("globalCounters");

  bool doGenerator=false, doGenParticles=false, doElectron=false, doMuon=false, doPhoton=false;
  std::vector<edm::ParameterSet> psets = iConfig.getParameter<std::vector<edm::ParameterSet>>("modules");
  for (unsigned int i=0; i<psets.size(); ++i) {
    std::string type = psets[i].getParameter<std::string>("type");

    if (type == "Electron") {
      container.push_back(new GlobeElectron(psets[i]));
      doElectron = true;
    } else if (type == "Common")
      container.push_back(new GlobeCommon(psets[i]));
    else if (type == "Generator") {
      container.push_back( new GlobeGenerator(iConfig));
      doGenerator = true;
    } else if (type == "GenParticle") {
      container.push_back(new GlobeGenParticles(iConfig));
      doGenParticles = true;
    } else if (type == "GenVertex")
      container.push_back(new GlobeGenVertices(iConfig));
    else if (type == "SimHit")
      container.push_back(new GlobeSimHits(iConfig));
    else if (type == "SimTrack")
      container.push_back(new GlobeSimTracks(iConfig));
    else if (type == "L1")
      container.push_back(new GlobeL1(iConfig)); 
    else if (type == "HKT")
      container.push_back(new GlobeHLT(iConfig)); 
    else if (type == "ECALCluster")
      container.push_back(new GlobeEcalClusters(iConfig));
    else if (type == "CaloTower")
      container.push_back(new GlobeCaloTowers(iConfig));
    else if (type == "Track")
      container.push_back(new GlobeTracks(iConfig));
    else if (type == "GsfTrack")
      container.push_back(new GlobeGsfTracks(iConfig));
    else if (type == "TrackingParticle")
      container.push_back(new GlobeTrackingParticles(iConfig));
    else if (type == "Vertex")
      container.push_back(new GlobeVertex(iConfig));
    else if (type == "Photon") {
      container.push_back(new GlobePhotons(iConfig));
      doPhoton = true;
    } else if (type == "Conversion")
      container.push_back(new GlobeConversions(iConfig));
    else if (type == "Muon") {
      container.push_back(new GlobeMuons(iConfig));
      doMuon = true;
    } else if (type == "MET")
      container.push_back(new GlobeMET(iConfig)); 
    else if (type == "GenJet")
      container.push_back(new GlobeGenJets(iConfig));
    else if (type == "Jet")
      container.push_back(new GlobeJets(iConfig));
    else if (type == "Jet")
      container.push_back(new GlobeSelector(iConfig));
    else if (type == "ReducedGen")
      container.push_back(new GlobeReducedGen(iConfig));
    else if (type == "PF")
      container.push_back(new GlobePFCandidates(iConfig));
    else if (type == "Rho")
      container.push_back(new GlobeRho(iConfig));
    else if (type == "PU")
      container.push_back(new GlobePileup(iConfig));
    else if (type == "PdfWeight")
      container.push_back(new GlobePdfWeights(iConfig));
    else if (type == "ECALHit")
      container.push_back(new GlobeEcalHits(iConfig)); 
    else if (type == "HCALCluster")
      container.push_back(new GlobeHcal(iConfig));
    else
      std::cout << "WARNING: unrecognized type " << type << " fix your configuation." <<std::endl;
  }

  debug_level = iConfig.getParameter<int>("Debug_Level");

  if(doGenerator and doGenParticles) {
    std::cout << "doGenerator and doGenParticles cannot be true at the same time." << std::endl;
    throw;
  }

  if (!doElectron and !doMuon and !doPhoton) {
    std::cout << "WARNING: EcalRecHits and HcalRecHits need Electrons, Muons or Photons." << std::endl;
    throw;
  }
  
  readConfiguration(iConfig);
  nProcessedEvents = 0;
}

GlobeAnalyzer::~GlobeAnalyzer() {}

void GlobeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  clock_t begin, t;
  if (TIMINGDEBUG)
    begin=clock();
  
  if(debug_level > 9) 
    std::cout << "GlobeAnalyzer: Begin" << std::endl;
  
  tot_events++;
  
  for (unsigned int i=0; i<container.size(); ++i) {
    if(debug_level > 2) 
      std::cout << "GlobeAnalyzer: " << container[i]->type() << ":" << container[i]->prefix() << std::endl;  
    container[i]->analyze(iEvent, iSetup);
    if (TIMINGDEBUG) {
      t=clock();
      std::cout << "Time elapsed " << i << ": " << double(diffclock1(t,begin)) << " ms"<< std::endl;
      std::cout << "Mem" << i << ": " << memory_usage() << std::endl;
    }
  }

  // fill the tree
  tree->Fill();
  nProcessedEvents++;
  
  if (nProcessedEvents % 500 == 0) {
    file->cd();
    tree->AutoSave();
    nProcessedEvents = 0;
  }

  
  //  tracks->GetAssociatedTrackingParticleIndex(iEvent,iSetup,trackingParticles);
  
  //  ecalrechits->analyze(iEvent, iSetup, std_electrons, muons, photons);
  //  hcalhits->analyze(iEvent, iSetup, std_electrons, muons, photons);
  
  //if (doGenParticles || doGenerator)
  //  selector_bits = selector->select(std_electrons, muons, photons, gen, leptons, reducedgen).to_ulong();
  //else
  //  selector_bits = selector->select(std_electrons, muons, photons).to_ulong();
  //if(debug_level > 2) 
  //  std::cout << "GlobeAnalyzer: selectorbits = " << selector_bits << std::endl;
  //if (selector_bits > 0) {
  //  sel_events++;
  
  if(debug_level > 9) 
    std::cout << "GlobeAnalyzer: End" << std::endl;
}

void GlobeAnalyzer::beginJob() { 

  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("event", "Event data");
  tree2 = new TTree("global_variables", "Global Parameters"); // need a different tree to fill once per job
  lumitree = new TTree("lumi", "Processed lumi sections");

  for (unsigned int i=0; i<container.size(); ++i) {
    container[i]->defineBranch(this);
    container[i]->defineLumiBranch(lumitree);
    container[i]->definePathBranch(tree2);
  }

  defineBranch();

  tot_events = 0;
  sel_events = 0;
}

int GlobeAnalyzer::findModule(std::string type, std::string prefix) {
  
  for (unsigned int i=0; i<container.size(); ++i)
    if (container[i]->type() == type) 
      if (prefix == "" or container[i]->prefix == prefix)
	return i;

  return -1;
}

void GlobeAnalyzer::endJob() { 

  fillTree();
  file->cd();
  
  int index = findModule("PU");
  if (index != -1) {
    TH1D* h = container[index]->getHisto(); 
    Int_t last_bin = h->GetNbinsX();
    h->SetBinContent(last_bin-1, h->GetBinContent(last_bin)+h->GetBinContent(last_bin-1));
    h->Write();  
    
    h = container[index]->getHistoTrue(); 
    last_bin = h->GetNbinsX();
    h->SetBinContent(last_bin-1, h->GetBinContent(last_bin)+h->GetBinContent(last_bin-1));
    h->Write();  
  }

  tree->Write(0, TObject::kWriteDelete);
  tree2->Write(0, TObject::kWriteDelete);
  lumitree->Write();

  file->Close();
}

void GlobeAnalyzer::endLuminosityBlock(const edm::LuminosityBlock & l, const edm::EventSetup & es) {

  int index = findModule("Common");
  if (index != -1) {
    container[index]->endLumiBlock(l, es);
    for(size_t ii=0; ii< globalCounters.size(); ++ii ) {
      edm::Handle<edm::MergeableCounter> ctrHandle;
      l.getByLabel(globalCountersNames[ii], ctrHandle);
      globalCounters[ii] += ctrHandle->value;
      globalCountersPerLumi[ii] = ctrHandle->value;
    }
    lumitree->Fill();
  }
}


void GlobeAnalyzer::readConfiguration(const edm::ParameterSet& iConfig) {
  parameters = new std::vector<std::string>;
  parameters->push_back(iConfig.dump());
}

void GlobeAnalyzer::defineBranch() {

  tree->Branch("selector_bits", &selector_bits, "selector_bits/I");

  tree2->Branch("version", &version, "version/I");
  tree2->Branch("type", &type, "type/I");
  tree2->Branch("tot_events", &tot_events, "tot_events/I");
  tree2->Branch("sel_events", &sel_events, "sel_events/I");
  tree2->Branch("parameters", "std::vector<std::string>", &parameters); 
  tree2->Branch("jobmaker", "std::string", &jobmaker); 
  
  globalCounters.clear();
  globalCounters.resize(globalCountersNames.size(),0);
  for(size_t ii=0; ii< globalCounters.size(); ++ii ) {
    tree2->Branch( globalCountersNames[ii].c_str(), &globalCounters[ii], (globalCountersNames[ii]+"/I").c_str() );
  }

  globalCountersPerLumi.clear();
  globalCountersPerLumi.resize(globalCountersNames.size(),0);
  for(size_t ii=0; ii< globalCountersPerLumi.size(); ++ii ) {
    lumitree->Branch( globalCountersNames[ii].c_str(), &globalCountersPerLumi[ii], (globalCountersNames[ii]+"/I").c_str() );
  }
  
}

void GlobeAnalyzer::fillTree() {
  
  type = 0 ;
  tree2->Fill();
}

DEFINE_FWK_MODULE(GlobeAnalyzer);
