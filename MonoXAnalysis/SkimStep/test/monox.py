import FWCore.ParameterSet.Config as cms

# Define the CMSSW process
process = cms.Process("SKIM")

# Load the standard set of modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

# Message Logger Stuff
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

# Display summary at the end of the job
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Is this a simulation or real data
isMC = False

# Do you want to make a skim or get the trees directly
doSkim = True

# Define the global tag depending on the sample type
if isMC:
    process.GlobalTag.globaltag = 'START53_V21::All'     # for Simulation
else:
    process.GlobalTag.globaltag = 'FT_53_V21_AN5::All'   # for Jan 22 2013 ReReco

# Define the input source
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring('file:monoxCopy.root'),
)

# Define the output
process.out = cms.OutputModule("PoolOutputModule", 
    outputCommands =  cms.untracked.vstring(), 
    fileName = cms.untracked.string('monoxSkim.root') 
)

# Define a set of good vertices and require atleast one of them
process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)

# Add vertex selection to a 'prePatSequence' which is exactly what the name suggests
process.prePatSequence = cms.Sequence(process.goodPrimaryVertices)

# Load the PAT sequence ... This is where the action begins
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Electron IDs -- the old VBTF WPs
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi as vbtfid
process.eidVBTFRel95 = vbtfid.simpleCutBasedElectronID.clone(electronQuality = '95relIso')
process.eidVBTFRel80 = vbtfid.simpleCutBasedElectronID.clone(electronQuality = '80relIso')
process.eidVBTFCom95 = vbtfid.simpleCutBasedElectronID.clone(electronQuality = '95cIso')
process.eidVBTFCom80 = vbtfid.simpleCutBasedElectronID.clone(electronQuality = '80cIso')

process.eidSequence = cms.Sequence(
    process.eidVBTFRel95 +
    process.eidVBTFRel80 +
    process.eidVBTFCom95 +
    process.eidVBTFCom80
)

process.prePatSequence += process.eidSequence

process.patElectrons.electronIDSources = cms.PSet(
    eidVBTFRel95     = cms.InputTag("eidVBTFRel95"),
    eidVBTFRel80     = cms.InputTag("eidVBTFRel80"),
    eidVBTFCom95     = cms.InputTag("eidVBTFCom95"),
    eidVBTFCom80     = cms.InputTag("eidVBTFCom80"),
)

# Embed the tracker track in PAT muons
process.patMuons.embedTrack = True

# Jet energy corrections
if isMC:
    corrLabels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute')
else:
    corrLabels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual')

# Use ak5PFJets 
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(
    process,
    cms.InputTag('ak5PFJets'),
    doJTA        = True,
    doBTagging   = True,
    jetCorrLabel = ('AK5PF', corrLabels),
    doType1MET   = True,
    genJetCollection=cms.InputTag('ak5GenJets'),
    doJetID      = True
)

# Use HPS taus
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process)

# Compute MET without PF Muons
process.particleFlowNoMuons =  cms.EDFilter("CandViewSelector",
    src = cms.InputTag("particleFlow"),
    cut = cms.string("muonRef.isNull")
)

from RecoMET.METProducers.PFMET_cfi import *
process.metNoMu = pfMet.clone()
process.metNoMu.src = cms.InputTag("particleFlowNoMuons")
process.metNoMu.METType = cms.string('NoMuMET')
process.metNoMu.alias = cms.string('NoMuMET')
process.metNoMu.calculateSignificance = cms.bool(False)

process.metNoMuSequence = cms.Sequence(process.particleFlowNoMuons + process.metNoMu)
process.prePatSequence += process.metNoMuSequence 

# Physics object selections
process.vetoMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("cleanPatMuons"),
    cut = cms.string("pt > 10 & abs(eta) < 2.4 & (isGlobalMuon | isTrackerMuon) & isPFMuon"),
)

process.vetoElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("cleanPatElectrons"),
    cut = cms.string('electronID("eidVBTFCom95") == 7')
)

process.vetoTaus = cms.EDFilter("PATTauSelector",
    src = cms.InputTag("cleanPatTaus"),
    cut = cms.string('pt > 10 & abs(eta) < 2.3 & tauID("decayModeFinding") > 0.5 & tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5 & tauID("againstMuonTight2") > 0.5 & tauID("againstElectronLoose") > 0.5')
)

process.vetoJets = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("cleanPatJets"),    
    cut = cms.string("pt > 30 & abs(eta) < 4.5 & electronEnergyFraction < 0.5 & muonEnergyFraction < 0.5"),
)

process.goodMuons = cms.EDFilter("PATMuonSelector",
    src = cms.InputTag("cleanPatMuons"),
    cut = cms.string("pt > 20 & abs(eta) < 2.4 & isGlobalMuon & isTrackerMuon & globalTrack.hitPattern.numberOfValidMuonHits > 0 & innerTrack.hitPattern.numberOfValidPixelHits > 0 & numberOfMatchedStations > 1 & innerTrack.hitPattern.trackerLayersWithMeasurement > 5 & globalTrack.normalizedChi2 < 10 & dB < 0.02 & (pfIsolationR04.sumChargedHadronPt + pfIsolationR04.sumNeutralHadronEt + pfIsolationR04.sumPhotonEt)/pt < 0.12"),
)

process.goodJets = cms.EDFilter("PATJetSelector",
    src = cms.InputTag("cleanPatJets"),
    cut = cms.string("pt > 110 & abs(eta) < 2.4 & electronEnergyFraction < 0.5 & muonEnergyFraction < 0.5 & chargedHadronEnergyFraction > 0.2 & neutralHadronEnergyFraction < 0.7 & neutralEmEnergyFraction < 0.7"),
)

# Remove MC matching when running on data
if not isMC:
    runOnData(process)

# Decide what we want to keep in the skim
process.out.outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_genParticles_*_*',
    'keep GenEventInfoProduct_generator_*_*',
    'keep *_addPileupInfo_*_*',
    'keep *_TriggerResults_*_*',
    'keep *_goodPrimaryVertices_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_cleanPatMuons_*_*',
    'keep *_cleanPatElectrons_*_*',
    'keep *_cleanPatPhotons_*_*',
    'keep *_cleanPatTaus_*_*',
    'keep *_cleanPatJets_*_*',
    'keep *_vetoMuons_*_*',
    'keep *_vetoElectrons_*_*',
    'keep *_vetoTaus_*_*',
    'keep *_vetoJets_*_*',
    'keep *_goodMuons_*_*',
    'keep *_goodJets_*_*',
    'keep recoPFMETs_*_*_*',
    'keep *_metNoMu_*_*',
    'keep *_*_rho_*'
)

# Use good primary vertices in the PAT sequence
massSearchReplaceAnyInputTag(process.patDefaultSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodPrimaryVertices"),True)

# Sequences to be run after PAT 
process.postPatSequence = cms.Sequence(process.vetoMuons + process.vetoElectrons + process.vetoTaus + process.vetoJets + process.goodMuons + process.goodJets)
if doSkim :
    process.patPath = cms.Path(process.prePatSequence + process.patDefaultSequence + process.postPatSequence)
    process.outPath = cms.EndPath(process.out)
    process.out.SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('patPath'))
    process.schedule = cms.Schedule(process.patPath, process.outPath)
else :
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
        triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
        jec = cms.string("ak5PFL1L2L3"),
        weight = cms.double(1.0),
        isControlSample = cms.bool(False)
    )

    process.postPatSequence += process.tree
    process.patPath = cms.Path(process.prePatSequence + process.patDefaultSequence + process.postPatSequence)
    process.schedule = cms.Schedule(process.patPath)
