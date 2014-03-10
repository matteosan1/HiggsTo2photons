
import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( (
        '/store/relval/CMSSW_6_2_5_cand1/RelValTTbar_13/GEN-SIM-RECO/PU_POSTLS162_V1_pu40bx50-v2/00000/8804A4DC-1F61-E311-900C-0025905A612E.root',
        ) );


