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

GlobeHLT::GlobeHLT(const edm::ParameterSet& iConfig, const char* n): nome(n) {
  
  edm::ParameterSet psetHLT = iConfig.getParameter<edm::ParameterSet>("HLTParameters");
  inputTag_ = psetHLT.getParameter<edm::InputTag>("TriggerResultsTag");
  hltTag_  = psetHLT.getParameter<edm::InputTag>("PrimaryTriggerResultsTag");

  hlt_filter_names = new std::vector<std::string>;
  hlt_path_names = new std::vector<std::vector<std::string>>;
  hlt_bit = new std::vector<unsigned short>;
  hlt_eta = new std::vector<std::vector<float>>;
  hlt_phi = new std::vector<std::vector<float>>;
  hlt_et =  new std::vector<std::vector<float>>;

  debug_level = iConfig.getParameter<int>("Debug_Level");
<<<<<<< HEAD
}

void GlobeHLT::defineBranch(GlobeAnalyzer* ana) {

  hlt_p4  = new TClonesArray("TLorentzVector", MAX_HLT);

  // Event Trigger
  ana->Branch("hlt_bit", "std::vector<unsigned short>", &hlt_bit);
  ana->Branch("hlt_path_names_HLT", "std::vector<std::string>", &hlt_path_names_HLT);

  // Trigger Candidates
  ana->Branch("hlt_n", &hlt_n, "hlt_n/I");
  ana->Branch("hlt_p4", "TClonesArray", &hlt_p4, 32000, 0);
  ana->Branch("hlt_candpath", "std::vector<std::vector<unsigned short> >", &hlt_candpath);
  ana->Branch("hlt_candpath2", "std::vector<std::vector<unsigned short> >", &hlt_candpath2);

  //filter passes
  ana->Branch("filter_names_HLT1", "std::vector<std::string>", &filter_names_HLT1);
  ana->Branch("filter_pass", "std::vector<unsigned int>", &filter_pass);

  //Mass filter decisions
  //RELATED to making sure OR is same as four separate paths
  /*
  ana->Branch("pass_Mass60_isoiso",&pass_Mass60_isoiso, "pass_Mass60_isoiso/I");
  ana->Branch("pass_Mass60_R9R9",&pass_Mass60_R9R9, "pass_Mass60_R9R9/I");
  ana->Branch("pass_Mass60_mix",&pass_Mass60_mix, "pass_Mass60_mix/I");
  ana->Branch("pass_Mass70_isoiso",&pass_Mass70_isoiso, "pass_Mass70_isoiso/I");
  ana->Branch("pass_Mass70_R9R9",&pass_Mass70_R9R9, "pass_Mass70_R9R9/I");
  ana->Branch("pass_Mass70_mix",&pass_Mass70_mix, "pass_Mass70_mix/I");
*/

  //Ele32_SC17 trigger objects
  ana->Branch("trg_SC_ele_n", &ElectronRefs0_n,"ElectronRefs0_n/I");
  ana->Branch("trg_SC_ele_eta", &ElectronRefs0_eta,"ElectronRefs0_eta[ElectronRefs0_n]/F");
  ana->Branch("trg_SC_ele_et", &ElectronRefs0_et,"ElectronRefs0_et[ElectronRefs0_n]/F");
  ana->Branch("trg_SC_ele_phi", &ElectronRefs0_phi,"ElectronRefs0_phi[ElectronRefs0_n]/F");
  ana->Branch("trg_ele_n", &ElectronRefs1_n,"ElectronRefs1_n/I");
  ana->Branch("trg_ele_eta", &ElectronRefs1_eta,"ElectronRefs1_eta[ElectronRefs1_n]/F");
  ana->Branch("trg_ele_et", &ElectronRefs1_et,"ElectronRefs1_et[ElectronRefs1_n]/F");
  ana->Branch("trg_ele_phi", &ElectronRefs1_phi,"ElectronRefs1_phi[ElectronRefs1_n]/F");


//Ele17_Ele8 trigger objects
  ana->Branch("trg8_ele_n", &ElectronRefs2_n,"ElectronRefs2_n/I");
  ana->Branch("trg8_ele_eta", &ElectronRefs2_eta,"ElectronRefs2_eta[ElectronRefs2_n]/F");
  ana->Branch("trg8_ele_et", &ElectronRefs2_et,"ElectronRefs2_et[ElectronRefs2_n]/F");
  ana->Branch("trg8_ele_phi", &ElectronRefs2_phi,"ElectronRefs2_phi[ElectronRefs2_n]/F");
  ana->Branch("trg17_ele_n", &ElectronRefs3_n,"ElectronRefs3_n/I");
  ana->Branch("trg17_ele_eta", &ElectronRefs3_eta,"ElectronRefs3_eta[ElectronRefs3_n]/F");
  ana->Branch("trg17_ele_et", &ElectronRefs3_et,"ElectronRefs3_et[ElectronRefs3_n]/F");
  ana->Branch("trg17_ele_phi", &ElectronRefs3_phi,"ElectronRefs3_phi[ElectronRefs3_n]/F");

//Ele17_Ele8_mass50 trigger objects
  ana->Branch("trg8_mass50_ele_n", &ElectronRefs4_n,"ElectronRefs4_n/I");
  ana->Branch("trg8_mass50_ele_eta", &ElectronRefs4_eta,"ElectronRefs4_eta[ElectronRefs4_n]/F");
  ana->Branch("trg8_mass50_ele_et", &ElectronRefs4_et,"ElectronRefs4_et[ElectronRefs4_n]/F");
  ana->Branch("trg8_mass50_ele_phi", &ElectronRefs4_phi,"ElectronRefs4_phi[ElectronRefs4_n]/F");
  ana->Branch("trg17_mass50_ele_n", &ElectronRefs5_n,"ElectronRefs3_n/I");
  ana->Branch("trg17_mass50_ele_eta", &ElectronRefs5_eta,"ElectronRefs3_eta[ElectronRefs3_n]/F");
  ana->Branch("trg17_mass50_ele_et", &ElectronRefs5_et,"ElectronRefs3_et[ElectronRefs3_n]/F");
  ana->Branch("trg17_mass50_ele_phi", &ElectronRefs5_phi,"ElectronRefs3_phi[ElectronRefs3_n]/F");


  //26_18 OR trigger objects
  ana->Branch("PhotonRefs0_n", &PhotonRefs0_n,"PhotonRefs0_n/I");
  ana->Branch("PhotonRefs0_eta", &PhotonRefs0_eta,"PhotonRefs0_eta[PhotonRefs0_n]/F");
  ana->Branch("PhotonRefs0_et", &PhotonRefs0_et,"PhotonRefs0_et[PhotonRefs0_n]/F");
  ana->Branch("PhotonRefs0_phi", &PhotonRefs0_phi,"PhotonRefs0_phi[PhotonRefs0_n]/F");
  ana->Branch("PhotonRefs1_n", &PhotonRefs1_n,"PhotonRefs1_n/I");
  ana->Branch("PhotonRefs1_eta", &PhotonRefs1_eta,"PhotonRefs1_eta[PhotonRefs1_n]/F");
  ana->Branch("PhotonRefs1_et", &PhotonRefs1_et,"PhotonRefs1_et[PhotonRefs1_n]/F");
  ana->Branch("PhotonRefs1_phi", &PhotonRefs1_phi,"PhotonRefs1_phi[PhotonRefs1_n]/F");
  ana->Branch("PhotonRefs3_n", &PhotonRefs3_n,"PhotonRefs3_n/I");
  ana->Branch("PhotonRefs3_eta", &PhotonRefs3_eta,"PhotonRefs3_eta[PhotonRefs3_n]/F");
  ana->Branch("PhotonRefs3_et", &PhotonRefs3_et,"PhotonRefs3_et[PhotonRefs3_n]/F");
  ana->Branch("PhotonRefs3_phi", &PhotonRefs3_phi,"PhotonRefs3_phi[PhotonRefs3_n]/F");
  ana->Branch("PhotonRefs4_n", &PhotonRefs4_n,"PhotonRefs4_n/I");
  ana->Branch("PhotonRefs4_eta", &PhotonRefs4_eta,"PhotonRefs4_eta[PhotonRefs4_n]/F");
  ana->Branch("PhotonRefs4_et", &PhotonRefs4_et,"PhotonRefs4_et[PhotonRefs4_n]/F");
  ana->Branch("PhotonRefs4_phi", &PhotonRefs4_phi,"PhotonRefs4_phi[PhotonRefs4_n]/F");

  //36_22 OR trigger objects
  ana->Branch("PhotonRefs5_n", &PhotonRefs5_n,"PhotonRefs5_n/I");
  ana->Branch("PhotonRefs5_eta", &PhotonRefs5_eta,"PhotonRefs5_eta[PhotonRefs5_n]/F");
  ana->Branch("PhotonRefs5_et", &PhotonRefs5_et,"PhotonRefs5_et[PhotonRefs5_n]/F");
  ana->Branch("PhotonRefs5_phi", &PhotonRefs5_phi,"PhotonRefs5_phi[PhotonRefs5_n]/F");
  ana->Branch("PhotonRefs6_n", &PhotonRefs6_n,"PhotonRefs6_n/I");
  ana->Branch("PhotonRefs6_eta", &PhotonRefs6_eta,"PhotonRefs6_eta[PhotonRefs6_n]/F");
  ana->Branch("PhotonRefs6_et", &PhotonRefs6_et,"PhotonRefs6_et[PhotonRefs6_n]/F");
  ana->Branch("PhotonRefs6_phi", &PhotonRefs6_phi,"PhotonRefs6_phi[PhotonRefs6_n]/F");
  ana->Branch("PhotonRefs8_n", &PhotonRefs8_n,"PhotonRefs8_n/I");
  ana->Branch("PhotonRefs8_eta", &PhotonRefs8_eta,"PhotonRefs8_eta[PhotonRefs8_n]/F");
  ana->Branch("PhotonRefs8_et", &PhotonRefs8_et,"PhotonRefs8_et[PhotonRefs8_n]/F");
  ana->Branch("PhotonRefs8_phi", &PhotonRefs8_phi,"PhotonRefs8_phi[PhotonRefs8_n]/F");
  ana->Branch("PhotonRefs9_n", &PhotonRefs9_n,"PhotonRefs9_n/I");
  ana->Branch("PhotonRefs9_eta", &PhotonRefs9_eta,"PhotonRefs9_eta[PhotonRefs9_n]/F");
  ana->Branch("PhotonRefs9_et", &PhotonRefs9_et,"PhotonRefs9_et[PhotonRefs9_n]/F");
  ana->Branch("PhotonRefs9_phi", &PhotonRefs9_phi,"PhotonRefs9_phi[PhotonRefs9_n]/F");

  //spin trigger
  ana->Branch("PhotonRefs10_n", &PhotonRefs10_n,"PhotonRefs10_n/I");
  ana->Branch("PhotonRefs10_eta", &PhotonRefs10_eta,"PhotonRefs10_eta[PhotonRefs10_n]/F");
  ana->Branch("PhotonRefs10_et", &PhotonRefs10_et,"PhotonRefs10_et[PhotonRefs10_n]/F");
  ana->Branch("PhotonRefs10_phi", &PhotonRefs10_phi,"PhotonRefs10_phi[PhotonRefs10_n]/F");

  ana->Branch("PhotonRefs11_n", &PhotonRefs11_n,"PhotonRefs11_n/I");
  ana->Branch("PhotonRefs11_eta", &PhotonRefs11_eta,"PhotonRefs11_eta[PhotonRefs11_n]/F");
  ana->Branch("PhotonRefs11_et", &PhotonRefs11_et,"PhotonRefs11_et[PhotonRefs11_n]/F");
  ana->Branch("PhotonRefs11_phi", &PhotonRefs11_phi,"PhotonRefs11_phi[PhotonRefs11_n]/F");
=======

  std::vector<std::string> temp_names = psetHLT.getParameter<std::vector<std::string>>("TriggerFilters");
>>>>>>> development

  for (unsigned int i=0; i<temp_names.size(); i++)
    hlt_filter_names->push_back(temp_names[i]);
  
 //temp_names.push_back("hltEG26CaloId10Iso50HcalIsoLastFilter");
 //temp_names.push_back("hltEG26R9Id85LastFilter");
 ////temp_names.push_back("HLTEG26R9Id85ORCaloId10Iso50LegCombLastFilter");
 //temp_names.push_back("hltEG18R9Id85LastFilterUnseeded");
 //temp_names.push_back("hltEG18CaloId10Iso50TrackIsoLastFilterUnseeded");
 //
 //temp_names.push_back("hltEG36CaloId10Iso50HcalIsoLastFilter");
 //temp_names.push_back("hltEG36R9Id85LastFilter");
 ////temp_names.push_back("HLTEG36R9Id85ORCaloId10Iso50LegCombLastFilter");
 //temp_names.push_back("hltEG22R9Id85LastFilterUnseeded");
 //temp_names.push_back("hltEG22CaloId10Iso50TrackIsoLastFilterUnseeded");
 //
 ////temp_names.push_back("hltEG18CaloId10Iso50TrackIsoDoubleLastFilterUnseeded");
 //temp_names.push_back("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter");
 //temp_names.push_back("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter");
 //
 ////spin trigger:
 //temp_names.push_back("hltEG10R9Id85LastFilterUnseeded");
 //temp_names.push_back("hltEG10CaloId10Iso50TrackIsoLastFilterUnseeded");
 //
 ///new electron triggers:
 //temp_names.push_back("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter");
 //temp_names.push_back("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter");
 //temp_names.push_back("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter");
 //temp_names.push_back("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter");
}

void GlobeHLT::definePathBranch(TTree* tree) {
  
  tree->Branch("hlt_path_names", "std::vector<std::vector<std::string>", &hlt_path_names);
  tree->Branch("hlt_filter_names", "std::vector<std::string>", &hlt_filter_names);
}

void GlobeHLT::defineBranch(GlobeAnalyzer* ana) {

  // Event Trigger
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
