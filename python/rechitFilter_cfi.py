import FWCore.ParameterSet.Config as cms

rechitTimeFilter = cms.EDFilter("RecHitFilter",
                                EBCollection  = cms.InputTag("reducedEcalRecHitsEB"),
                                EECollection  = cms.InputTag("reducedEcalRecHitsEE"),
                                timeCut       =  cms.double(10.),
                                )

