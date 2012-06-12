#ifndef __MASSFACTORIZEDANALYSIS__
#define __MASSFACTORIZEDANALYSIS__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis.h"
#include "StatAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
#include <iostream>
#include <fstream>
#include "math.h"

// ------------------------------------------------------------------------------------
class MassFactorizedMvaAnalysis : public StatAnalysis 
{
 public:
    
    MassFactorizedMvaAnalysis();
    virtual ~MassFactorizedMvaAnalysis();
    
    virtual const std::string & name() const { return name_; };
    
    // LoopAll analysis interface implementation
    void Init(LoopAll&);
    void Term(LoopAll&);
    
    virtual void ResetAnalysis();
    //// virtual void Analysis(LoopAll&, Int_t); 

    virtual int GetBDTBoundaryCategory(float,bool,bool);

    void fillZeeControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, 
	                     const TLorentzVector & Higgs, float lead_r9, float sublead_r9,
			     float phoid_mvaout_lead, float phoid_mvaout_sublead, 
			     float diphobdt_output, float sigmaMrv, float sigmaMwv, float vtxProb,
			     int diphoton_id, int category, int selectioncategory, float evweight, LoopAll & l );

    bool doPhotonMvaIdSyst;
    bool doPhotonMvaIdSmear;
    bool doRegressionSmear, doRegressionSyst;

    std::string bdtTrainingPhilosophy;
    std::string photonLevelMvaUCSD  ;
    std::string eventLevelMvaUCSD   ;                    
    std::string photonLevelMvaMIT_EB;
    std::string photonLevelMvaMIT_EE;
    std::string eventLevelMvaMIT    ;
    std::string photonLevelNewIDMVA_EB;
    std::string photonLevelNewIDMVA_EE;

 protected:

    virtual bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, int & category, int & diphoton_id,
		      bool & isCorrectVertex, float &kinematic_bdtout,
		      bool isSyst=false, 
		      float syst_shift=0., bool skipSelection=false,
		      BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0); 

    EnergySmearer  *eRegressionSmearer ; 
    DiPhoEfficiencySmearer *photonMvaIdSmearer ;
    
    std::string name_;
    std::map<int,std::string> signalLabels;
    
    
    HggVertexAnalyzer vtxAna_;
    HggVertexFromConversions vtxConv_;
};

#endif


// Local Variables:
// mode: c++
// mode: sensitive
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4