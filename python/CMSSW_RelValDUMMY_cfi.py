
import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( (
        'file:/afs/cern.ch/user/s/sani/work/CMSSW_7_1_0_pre3/src/test.root',
        #'file:/tmp/sani/test.root',
#/store/relval/CMSSW_7_1_0_pre3/RelValZEE_13/GEN-SIM-RECO/PU50ns_POSTLS171_V2-v1/00000/E4821204-FCA0-E311-837E-02163E00EA2B.root',
    ) );


