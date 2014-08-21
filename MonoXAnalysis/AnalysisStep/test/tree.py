import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("Tree")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START53_V23::All'

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())
process.source.fileNames = [
    'file:/hadoop/cms/store/user/avartak/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/MONOX_53X_V02_v1/4d6d2f0c3bfd3fade3e91c3ca6d2886c/monoxSkim_1_1_1Sy.root'
]

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.TFileService = cms.Service("TFileService", fileName = cms.string("tree.root"))

process.ak5PFL1Fastjet = cms.ESProducer('L1FastjetCorrectionESProducer',
    level       = cms.string('L1FastJet'),
    algorithm   = cms.string('AK5PF'),
    srcRho      = cms.InputTag('kt6PFJets','rho')
)
process.ak5PFL2Relative = cms.ESProducer('LXXXCorrectionESProducer',
    level     = cms.string('L2Relative'),
    algorithm = cms.string('AK5PF')
)
process.ak5PFL3Absolute = cms.ESProducer('LXXXCorrectionESProducer',
    level     = cms.string('L3Absolute'),
    algorithm = cms.string('AK5PF')
)
process.ak5PFL1L2 = cms.ESProducer('JetCorrectionESChain',
    correctors = cms.vstring('ak5PFL1Fastjet','ak5PFL2Relative')
)
process.ak5PFL1L2L3 = cms.ESProducer('JetCorrectionESChain',
    correctors = cms.vstring('ak5PFL1Fastjet','ak5PFL2Relative','ak5PFL3Absolute')
)

process.tree = cms.EDAnalyzer("MonoJetTreeMaker",
    pileup = cms.InputTag("addPileupInfo"),
    vertices = cms.InputTag("goodPrimaryVertices"),
    muons = cms.InputTag("cleanPatMuons"),
    electrons = cms.InputTag("cleanPatElectrons"),
    taus = cms.InputTag("vetoTaus"),
    jets = cms.InputTag("cleanPatJets"),
    gens = cms.InputTag("genParticles"),
    met = cms.InputTag("metNoMu"),
    jec = cms.string("ak5PFL1L2L3"),
    weight = cms.double(1.0),
    isControlSample = cms.bool(True)
)

process.treePath = cms.Path(
    process.tree 
)

process.schedule = cms.Schedule(process.treePath)
