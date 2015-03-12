import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import copy

#DATA TYPE
flagData = 'OFF'
flagMC = 'ON'
flagFastSim = 'OFF'
flagPG = 'OFF'

#SKIM TYPE
flagSkimDiphoton = 'OFF'
flagSkimPJet = 'OFF'
flagVLPreSelection = 'OFF'
flagMyPreSelection = 'OFF'
flagNoSkim = 'ON'
flagMuMuSkim = 'OFF'
flagMMgSkim = 'OFF'
flagSkimworz = 'OFF'
flagSkim1El = 'OFF'
flagAddPdfWeight = 'OFF'
flagSkimHmm = 'OFF'
flagSkimHee = 'OFF'
flagSkimMu = 'OFF'

#ADDITIONAL OPTIONS
flagAOD = 'ON'
jobMaker = 'job_maker --MC --AOD --noskim --skimBranches=smallest --scheduler=remoteGlidein --storage=eos --datasetpath=/GJet_Pt-15to3000_Tune4C_14TeV_pythia8/TP2023SHCALDR-SHCALJan23_PU140BX25_PH2_1K_FB_V6-v1/GEN-SIM-RECO --hlttag=HLT --disable=doHLT,doL1,doCaloTower,doMuon,doJet_algo1,doJet_algo2,doJet_algo3,doJet_algoPF1,doJet_algoPF2,doJet_algoPF3,doGenJet_algo1,doGenJet_algo2,doGenJet_algo3,doGsfTracks,doTracks --globaltag=POSTLS161_V15::All --events=-1 --outputdir=group/phys_higgs/cmshgg/processed/TP_14TeV_Phase2/'

if (not((flagNoSkim is 'ON') ^ (flagSkimDiphoton is 'ON') ^ (flagMMgSkim is 'ON') ^ (flagVLPreSelection is 'ON') ^ (flagSkim1El is 'ON') ^ (flagSkimworz is 'ON') ^ (flagMyPreSelection is 'ON') ^ (flagSkimPJet is 'ON')^ (flagSkimHmm is 'ON')^ (flagSkimHee is 'ON') ^ (flagSkimMu is 'ON'))):
  print "You must skim or not skim... these are your options"
  exit(-1)

process = cms.Process("Globe") 
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
process.load("Configuration.Geometry.GeometryExtended2023SHCalNoTaperReco_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.Services_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V6::All', '') 

process.load("HiggsAnalysis.HiggsTo2photons.h2ganalyzer_61X_cfi")

process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedElectrons = cms.PSet(
  initialSeed = cms.untracked.uint32(1),
  engineName = cms.untracked.string('TRandom3')
  ),
                                                   )

process.load("HiggsAnalysis.HiggsTo2photons.rechitFilter_cfi")

hltLabel = 'HLT'
#
process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  SkipEvent = cms.untracked.vstring('FatalRootError','InvalidReference')
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.dummySelector = cms.EDFilter("CandViewCountFilter",
                                     src = cms.InputTag("gsfElectrons"),
                                     minNumber = cms.uint32(0)
                                     )

process.eventFilter1 = cms.Sequence(process.dummySelector)
process.eventFilter2 = cms.Sequence(process.dummySelector)

process.h2ganalyzer.RootFileName = 'GJet_Pt-15to3000_Tune4C_14TeV_pythia8_TP2023SHCALDR-SHCALJan23_PU140BX25_PH2_1K_FB_V6-v1.root'
process.h2ganalyzer.Debug_Level = 0


##---------------------ELECTRON REGRESSION AND SMEARING ------------------------------
process.load("EgammaAnalysis.ElectronTools.calibratedElectrons_cfi")

## dataset to correct
if (flagMC == 'ON'):
  process.calibratedElectrons.isMC = cms.bool(True)
  process.calibratedElectrons.inputDataset = cms.string("Summer12_DR53X_HCP2012")
else:
  process.calibratedElectrons.isMC = cms.bool(False)
  process.calibratedElectrons.inputDataset = cms.string("Moriond2013")
  
process.calibratedElectrons.updateEnergyError = cms.bool(True)
process.calibratedElectrons.applyCorrections = cms.int32(1)
process.calibratedElectrons.smearingRatio = cms.double(0.607)
process.calibratedElectrons.verbose = cms.bool(False)
#process.calibratedElectrons.synchronization = cms.bool(True) 

process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('gsfElectrons')
process.eleRegressionEnergy.inputCollectionType = cms.uint32(0)
process.eleRegressionEnergy.useRecHitCollections = cms.bool(True)
process.eleRegressionEnergy.produceValueMaps = cms.bool(True)

process.calibratedElectronsSequence = cms.Sequence(process.eleRegressionEnergy*process.calibratedElectrons)

process.load("RecoMET.METFilters.ecalLaserCorrFilter_cfi")
process.load("RecoParticleFlow.PFProducer.particleFlowTmpPtrs_cfi")
process.particleFlowTmpPtrs.src = cms.InputTag("particleFlow")
process.load("CommonTools.ParticleFlow.pfPileUp_cfi")
process.load("CommonTools.ParticleFlow.pfNoPileUp_cff")
process.pfPileUp.PFCandidates = cms.InputTag("particleFlow")
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

# event counters
process.processedEvents = cms.EDProducer("EventCountProducer")
process.eventCounters = cms.Sequence(process.processedEvents)

process.eventCounters = cms.Sequence(process.processedEvents)
process.h2ganalyzer.globalCounters.extend(['processedEvents']) 
process.h2ganalyzerPath = cms.Sequence(process.h2ganalyzer)


# ----------------------------------------------------------------------
# plugins for JetFlavour
# ----------------------------------------------------------------------
process.myPartons = cms.EDProducer("PartonSelector",
    src = cms.InputTag("genParticles"),
    withLeptons = cms.bool(False)
)

process.flavourByRef = cms.EDProducer("JetPartonMatcher",
    jets = cms.InputTag("ak5GenJets"),
    coneSizeToAssociate = cms.double(0.3),
    partons = cms.InputTag("myPartons")
)

process.flavourByVal = cms.EDProducer("JetFlavourIdentifier",
    srcByReference = cms.InputTag("flavourByRef"),
    physicsDefinition = cms.bool(False)
)


#################################################
# Define path, first for AOD case then for RECO #
#################################################
process.p11 = cms.Path(process.eventCounters*process.eventFilter1*process.particleFlowTmpPtrs*process.pfNoPileUpSequence*process.calibratedElectronsSequence*process.myPartons*process.flavourByRef*process.flavourByVal)
process.p11 *= (process.rechitTimeFilter*process.h2ganalyzerPath)

process.p12 = copy.deepcopy(process.p11)
process.p12.replace(process.eventFilter1, process.eventFilter2)

#################################################
# End of Path definition                        #
#################################################

process.h2ganalyzer.JobMaker = jobMaker

if flagAOD is 'ON':
  process.h2ganalyzer.doAodSim = True
  process.h2ganalyzer.doHcal = False
  process.h2ganalyzer.doHFHcal = False
  process.h2ganalyzer.doPreshowerHits = False
else:
  process.h2ganalyzer.doAodSim = False
  process.h2ganalyzer.doHcal = True
  process.h2ganalyzer.doHFHcal = True
  process.h2ganalyzer.doPreshowerHits = True
  process.h2ganalyzer.EcalHitEBColl = cms.InputTag("ecalRecHit","EcalRecHitsEB")
  process.h2ganalyzer.EcalHitEEColl = cms.InputTag("ecalRecHit","EcalRecHitsEE")
  process.h2ganalyzer.EcalHitESColl = cms.InputTag("ecalPreshowerRecHit","EcalRecHitsES")
  process.h2ganalyzer.HcalHitsBEColl = cms.InputTag("hbhereco")
  process.h2ganalyzer.HcalHitsFColl = cms.InputTag("hfreco")
  process.h2ganalyzer.HcalHitsHoColl = cms.InputTag("horeco")
  process.h2ganalyzer.BarrelBasicClusterColl = cms.InputTag("")
  process.h2ganalyzer.BarrelBasicClusterShapeColl = cms.InputTag("multi5x5BasicClusters","multi5x5BarrelShapeAssoc")
  process.h2ganalyzer.JetTrackAssociationColl_algo3 = cms.InputTag("kt4JetTracksAssociatorAtVertex")

if (flagPG is 'OFF'):
  process.h2ganalyzer.doL1 = True
  process.h2ganalyzer.doHLT = True
else:
  process.h2ganalyzer.doL1 = False
  process.h2ganalyzer.doHLT = False
  process.h2ganalyzer.doJet_algoPF3 = False
  process.h2ganalyzer.doParticleGun = True


process.h2ganalyzer.HLTParameters.PrimaryTriggerResultsTag = cms.InputTag("TriggerResults","", hltLabel)
process.h2ganalyzer.HLTParameters.useSecondaryTrigger = cms.bool(False)
process.h2ganalyzer.HLTParameters.TriggerResultsTag = cms.InputTag("hltTriggerSummaryAOD","", hltLabel)

from HiggsAnalysis.HiggsTo2photons.h2ganalyzer_SkimBranch_smallest_cfi import *
process.h2ganalyzer.branchesToSkim = branchesToSkim


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        '/store/mc/TP2023SHCALDR//GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_14TeV-pythia6/GEN-SIM-RECO/SHCALJan23_PU140BX25_PH2_1K_FB_V6-v1/20000/F2584185-DCC1-E411-A34A-001E67398E62.root',
        '/store/mc/TP2023SHCALDR//GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_14TeV-pythia6/GEN-SIM-RECO/SHCALJan23_PU140BX25_PH2_1K_FB_V6-v1/20000/E6A2C1C5-65C1-E411-B0B5-002590A371AC.root',
        '/store/mc/TP2023SHCALDR//GJet_Pt-20to40_doubleEMEnriched_TuneZ2star_14TeV-pythia6/GEN-SIM-RECO/SHCALJan23_PU140BX25_PH2_1K_FB_V6-v1/20000/E6EBCFBA-F0C1-E411-9CF7-001E67398052.root',
 )
)  

process.h2ganalyzer.doHLT = False
process.h2ganalyzer.doL1 = False
process.h2ganalyzer.doCaloTower = False
process.h2ganalyzer.doMuon = False
process.h2ganalyzer.doJet_algo1 = False
process.h2ganalyzer.doJet_algo2 = False
process.h2ganalyzer.doJet_algo3 = False
process.h2ganalyzer.doJet_algoPF1 = False
process.h2ganalyzer.doJet_algoPF2 = False
process.h2ganalyzer.doJet_algoPF3 = False
process.h2ganalyzer.doGenJet_algo1 = False
process.h2ganalyzer.doGenJet_algo2 = False
process.h2ganalyzer.doGenJet_algo3 = False
process.h2ganalyzer.doGsfTracks = False
#process.h2ganalyzer.doTracks = False
process.h2ganalyzer.doVertices_nobs=False
process.h2ganalyzer.doMet=False
process.h2ganalyzer.dotcMet=False
process.h2ganalyzer.h2gAnalyzerVersion = 'V15_00_11'

process.h2ganalyzer.EcalHitEBColl = cms.InputTag("ecalRecHit:EcalRecHitsEB")
process.h2ganalyzer.EcalHitEEColl = cms.InputTag("ecalRecHit:EcalRecHitsEK")

process.rechitTimeFilter.EBCollection = cms.InputTag("ecalRecHit:EcalRecHitsEB")
process.rechitTimeFilter.EECollection = cms.InputTag("ecalRecHit:EcalRecHitsEK")

process.h2ganalyzer.PhotonCollStd = cms.InputTag("mustachePhotons")

process.h2ganalyzer.EndcapBasicClusterColl = cms.InputTag("particleFlowSuperClusterECAL:particleFlowBasicClusterECALEndcap")
process.h2ganalyzer.EndcapBasicClusterShapeColl = cms.InputTag("particleFlowSuperClusterECAL:PFClusterAssociationEBEE")

process.h2ganalyzer.BarrelHybridClusterColl = cms.InputTag("particleFlowSuperClusterECAL:particleFlowBasicClusterECALBarrel")
process.h2ganalyzer.BarrelHybridClusterShapeColl = cms.InputTag("particleFlowSuperClusterECAL:PFClusterAssociationEBEE")

process.h2ganalyzer.EndcapSuperClusterColl = cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALEndcapWithPreshower")
process.h2ganalyzer.HybridSuperClusterColl = cms.InputTag("particleFlowSuperClusterECAL:particleFlowSuperClusterECALBarrel")
   
