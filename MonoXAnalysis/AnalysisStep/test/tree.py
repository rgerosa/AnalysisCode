import FWCore.ParameterSet.Config as cms

# Define the CMSSW process
process = cms.Process("TREE")

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True) 
)

# How many events to process
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1)
)

# Is this a simulation or real data
isMC = True

# Is this a Gamma+Jets Sample
isGammaSample = False

# Define the global tag depending on the sample type
if isMC:
    process.GlobalTag.globaltag = 'START53_V27::All'     # for Simulation
else:
    process.GlobalTag.globaltag = 'FT_53_V21_AN5::All'   # for Jan 22 2013 ReReco

# Define the input source
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('/store/mc/Summer12_DR53X/DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v2/0000/0006A6CC-D5D2-E111-9094-00266CFFBF84.root')
)

# Define the output -- Needed for PAT
process.out = cms.OutputModule("PoolOutputModule", 
    outputCommands =  cms.untracked.vstring(), 
    fileName = cms.untracked.string('monoxSkim.root') 
)

# Select good primary vertices
process.goodVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)

# Define rho to be used for photon isolation
from RecoJets.Configuration.RecoPFJets_cff import *
process.kt6PFJets25 = kt6PFJets.clone( doRhoFastjet = True )
process.kt6PFJets25.Rho_EtaMax = cms.double(2.5)

# Load the PAT sequence
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Define ak5PFJets without leptons and pileup
process.pfPileUp = cms.EDProducer("PFPileUp",
    checkClosestZVertex = cms.bool(True),
    Enable = cms.bool(True),
    PFCandidates = cms.InputTag("particleFlow"),
    verbose = cms.untracked.bool(False),
    Vertices = cms.InputTag("goodVertices")
)

process.particleFlowClean = cms.EDProducer("PFCleaner",
    src = cms.InputTag("particleFlow"),
    vertices = cms.InputTag("goodVertices"),
    pfpileup = cms.InputTag("pfPileUp"),
    rho = cms.InputTag("kt6PFJets25", "rho"),
    muselection = cms.string("pt > 10 && abs(eta) < 2.4 && (isGlobalMuon || isTrackerMuon)"),
    eleselection = cms.string("pt > 10 && abs(superCluster.eta) < 2.5 && gsfTrack.trackerExpectedHitsInner.numberOfHits <= 1 && " +
                "((isEB && sigmaIetaIeta < 0.01 && hadronicOverEm < 0.15 && abs(deltaPhiSuperClusterTrackAtVtx) < 0.8 && abs(deltaEtaSuperClusterTrackAtVtx) < 0.007) || " +
                " (isEE && sigmaIetaIeta < 0.03 && hadronicOverEm < 0.07 && abs(deltaPhiSuperClusterTrackAtVtx) < 0.7 && abs(deltaEtaSuperClusterTrackAtVtx) < 0.01))"),
    photonselection = cms.string("pt > 100 && abs(eta) < 2.5 && !hasPixelSeed && r9 > 0.9 && " +
                "((isEB && sigmaIetaIeta < 0.011) || (isEE && sigmaIetaIeta < 0.030))"),
    d0cut = cms.double(0.2),
    dzcut = cms.double(0.5),
    muisocut = cms.double(0.2),
    eleisocut = cms.double(0.2),
    vetophotons = cms.bool(isGammaSample)
)


from RecoJets.JetProducers.ak5PFJets_cfi import *
process.ak5PFJetsClean = ak5PFJets.clone(src = cms.InputTag("particleFlowClean", "pfcands"))

from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(
    process,
    cms.InputTag('ak5PFJetsClean'),
    doJTA            = True,
    doBTagging       = True,
    jetCorrLabel     = ('AK5PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute']),
    doType1MET       = True,
    genJetCollection = cms.InputTag('ak5GenJets'),
    doJetID          = True
)
process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
    cms.InputTag('pfMETcorrType0'),
    cms.InputTag('pfJetMETcorr', 'type1')        
)
process.pfJetMETcorr.src = "ak5PFJetsClean" 

# Use HPS taus
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

process.selectedPatTaus.cut = "pt > 20 && abs(eta) < 2.3 && tauID('decayModeFinding') > 0.5 && tauID('byLooseCombinedIsolationDeltaBetaCorr') > 0.5 && tauID('againstMuonTight2') > 0.5 && tauID('againstElectronLoose') > 0.5"

# Remove MC matching when running on data
if not isMC:
    runOnData(process)

# Use good primary vertices in the PAT sequence
massSearchReplaceAnyInputTag(process.patDefaultSequence, cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodVertices"), True)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tree.root"))

process.tree = cms.EDAnalyzer("MonoJetTreeMaker",
    pileup = cms.InputTag("addPileupInfo"),
    vertices = cms.InputTag("goodVertices"),
    muons = cms.InputTag("particleFlowClean", "muons"),
    electrons = cms.InputTag("particleFlowClean", "electrons"),
    taus = cms.InputTag("selectedPatTaus"),
    jets = cms.InputTag("patJets"),
    gens = cms.InputTag("genParticles"),
    pfmet = cms.InputTag("pfMet"),
    t1pfmet = cms.InputTag("pfType1CorrectedMet"),
    calomet = cms.InputTag("met"),
    pfpileup = cms.InputTag("pfPileUp"),
    photons = cms.InputTag("particleFlowClean", "photons"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    weight = cms.double(1000.0*34.1/2662137.0),
    isWorZMCSample = cms.bool(True),
    isPhotonSample = cms.bool(False)
)

# Select events passing the monojet triggers and having MET > 200 GeV
process.triggerfilter = cms.EDFilter("HLTCheckFilter",
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    triggerPaths   = cms.vstring('HLT_MET120_HBHENoiseCleaned', 'MonoCentralPFJet80_PFMETnoMu95_NHEF0p95', 'MonoCentralPFJet80_PFMETnoMu105_NHEF0p95')   
)

process.metfilter = cms.EDFilter("METFilter",
    muons = cms.InputTag("particleFlowClean", "muons"),
    type1pfmet = cms.InputTag("pfType1CorrectedMet"),
    metcut = cms.double(200.0)
)

process.treePath = cms.Path(process.goodVertices + process.tree)
#if not isGammaSample : 
#    process.treePath = cms.Path(process.goodVertices + process.triggerfilter + process.metfilter + process.tree)
#else :
#    process.treePath = cms.Path(process.goodVertices + process.tree)
