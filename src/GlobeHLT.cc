#include "HiggsAnalysis/HiggsTo2photons/interface/Limits.h"
#include "HiggsAnalysis/HiggsTo2photons/plugins/GlobeAnalyzer.h"

#include "HiggsAnalysis/HiggsTo2photons/interface/GlobeHLT.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"

#include "TLorentzVector.h"

#include <bitset>

void set(unsigned int& x, int bit) {
  x |= (1 << bit);
}

void reset(unsigned int& x, int bit) {
  x &= ~ (1 << bit);
}

bool check(unsigned int x, int bit) {
  return (x && (1 << bit));
}

GlobeHLT::GlobeHLT(const edm::ParameterSet& iConfig) {
  
  GlobeBase::GlobeBase(iConfig);
  order = -1;

  edm::ParameterSet psetHLT = iConfig.getParameter<edm::ParameterSet>("HLTParameters");
  inputTag_ = psetHLT.getParameter<edm::InputTag>("TriggerResultsTag");
  hltTag_  = psetHLT.getParameter<edm::InputTag>("PrimaryTriggerResultsTag");

  hlt_filter_names = new std::vector<std::string>;
  hlt_path_names = new std::vector<std::vector<std::string>>;
  hlt_bit = new std::vector<unsigned short>;
  hlt_eta = new std::vector<std::vector<float>>;
  hlt_phi = new std::vector<std::vector<float>>;
  hlt_et =  new std::vector<std::vector<float>>;

  std::vector<std::string> temp_names = psetHLT.getParameter<std::vector<std::string>>("TriggerFilters");

  for (unsigned int i=0; i<temp_names.size(); i++)
    hlt_filter_names->push_back(temp_names[i]);
}

void GlobeHLT::definePathBranch(TTree* tree) {
  
  tree->Branch("hlt_path_names", "std::vector<std::vector<std::string>", &hlt_path_names);
  tree->Branch("hlt_filter_names", "std::vector<std::string>", &hlt_filter_names);
}

void GlobeHLT::defineBranch(GlobeAnalyzer* ana) {

  GlobeBase::defineBranch(ana);
  ana->Branch("hlt_path_n", &hlt_path_n, "hlt_path_n/I");
  ana->Branch("hlt_bit", "std::vector<unsigned short>", &hlt_bit);
  
  ana->Branch("hlt_eta", "std::vector<std::vector<float>>", &hlt_eta);
  ana->Branch("hlt_et",  "std::vector<std::vector<float>>", &hlt_et);
  ana->Branch("hlt_phi", "std::vector<std::vector<float>>", &hlt_phi);
}

bool GlobeHLT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  edm::Handle<trigger::TriggerEvent> triggerObj;
  iEvent.getByLabel(inputTag_, triggerObj);

  // HLT
  hlt_bit->clear();
  hlt_eta->clear();
  hlt_phi->clear();
  hlt_et->clear();

  bool changed = false;
  configProvider.init(iEvent.getRun(), iSetup, hltTag_.process(), changed);
  edm::Handle<edm::TriggerResults> h_triggerResults_HLT;
  iEvent.getByLabel(hltTag_, h_triggerResults_HLT);
  
  if (changed) {
    std::vector<std::string> temp_path_names;

    if(debug_level > 9) 
      std::cout << "Fill names HLT" << std::endl;

    for (size_t i = 0; i < configProvider.size(); ++i)
      temp_path_names.push_back(configProvider.triggerName(i));
    
    hlt_path_names->push_back(temp_path_names);
  }

  if (h_triggerResults_HLT.isValid()) {
    hlt_path_n = hlt_path_names->size()-1;
    
    // Trigger Results
    for (size_t i = 0; i < configProvider.size(); ++i) {
      if(h_triggerResults_HLT->accept(i))
	hlt_bit->push_back((unsigned short)(i));
    }
  }
  
  if(!triggerObj.isValid()) 
    throw(cms::Exception("Release Validation Error") << "RAW-type HLT results not found" );

  //added Trigger Objects
  edm::InputTag trigResultsTag("TriggerResults","","HLT");
  edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); 
  edm::Handle<trigger::TriggerEvent> trigEvent; 
  iEvent.getByLabel(trigEventTag,trigEvent);


  std::vector<trigger::TriggerObject> hltRefs;
  std::vector<std::string>::const_iterator itFilter;

  if ( trigEvent.isValid() ){
    const trigger::TriggerObjectCollection & triggerObjects = trigEvent -> getObjects();
    trigger::size_type n_filters   = trigEvent->sizeFilters();

    for (itFilter = hlt_filter_names->begin(); itFilter != hlt_filter_names->end(); ++itFilter) {
      trigger::size_type filter1_idx = trigEvent->filterIndex(edm::InputTag(*itFilter, "", hltTag_.process()));
      
      if (filter1_idx < n_filters) {
	const trigger::Keys & triggerKeys(trigEvent->filterKeys(filter1_idx));
	const int nkeys = triggerKeys.size();

	std::vector<float> temp_eta, temp_phi, temp_et;
	for (int ikey = 0; ikey < nkeys; ++ikey ) {
	  trigger::TriggerObject obj = triggerObjects[triggerKeys[ikey]];
	  temp_eta.push_back(obj.eta());
	  temp_et.push_back(obj.et());
	  temp_phi.push_back(obj.phi());
	  
	}

	hlt_eta->push_back(temp_eta);
	hlt_et->push_back(temp_et);
	hlt_phi->push_back(temp_phi);
      }
    }
  }

  return true;
}
