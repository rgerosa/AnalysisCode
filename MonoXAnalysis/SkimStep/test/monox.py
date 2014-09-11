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

# Do you want to make a skim or get the trees directly
doSkim = False

# Define the global tag depending on the sample type
if isMC:
    process.GlobalTag.globaltag = 'START53_V21::All'     # for Simulation
else:
    process.GlobalTag.globaltag = 'FT_53_V21_AN5::All'   # for Jan 22 2013 ReReco

# Define the input source
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('/store/mc/Summer12_DR53X/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v1/0000/1CDD7B73-E8E1-E111-9694-0030487D8661.root')
)

# Define the output
process.out = cms.OutputModule("PoolOutputModule", 
    outputCommands =  cms.untracked.vstring(), 
    fileName = cms.untracked.string('monoxSkim.root') 
)

# Define a set of good vertices and require atleast one of them
process.goodVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlinePrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)

# Load the PAT sequence ... This is where the action begins
process.load("PhysicsTools.PatAlgos.patSequences_cff")

# Electron IDs -- the old VBTF WPs
import ElectroWeakAnalysis.WENu.simpleCutBasedElectronIDSpring10_cfi as vbtfid
process.eidVBTFRel95 = vbtfid.simpleCutBasedElectronID.clone(electronQuality = '95relIso')
process.eidVBTFRel80 = vbtfid.simpleCutBasedElectronID.clone(electronQuality = '80relIso')
process.eidVBTFCom95 = vbtfid.simpleCutBasedElectronID.clone(electronQuality = '95cIso')
process.eidVBTFCom80 = vbtfid.simpleCutBasedElectronID.clone(electronQuality = '80cIso')

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
    'keep *_goodVertices_*_*',
    'keep *_offlineBeamSpot_*_*',
    'keep *_patMuons_*_*',
    'keep *_patElectrons_*_*',
    'keep *_patPhotons_*_*',
    'keep *_patTaus_*_*',
    'keep *_patJets_*_*',
    'keep recoPFMETs_*_*_*',
    'keep recoCaloMETs_met_*_*',
    'keep *_*_rho_*'
)

# Use good primary vertices in the PAT sequence
massSearchReplaceAnyInputTag(process.patDefaultSequence,cms.InputTag("offlinePrimaryVertices"), cms.InputTag("goodVertices"),True)

# Sequences to be run after PAT 
if doSkim :
    process.Flag_goodVertices = cms.Path(process.goodVertices)
    process.outPath = cms.EndPath(process.out)
    process.out.SelectEvents   = cms.untracked.PSet(SelectEvents = cms.vstring('Flag_goodVertices'))
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
    process.ak5PFResidual = cms.ESProducer('LXXXCorrectionESProducer',
        level     = cms.string('L2L3Residual'),
        algorithm = cms.string('AK5PF')
    )

    if isMC:
        corrections = cms.vstring('ak5PFL1Fastjet','ak5PFL2Relative','ak5PFL3Absolute')
    else:
        corrections = cms.vstring('ak5PFL1Fastjet','ak5PFL2Relative','ak5PFL3Absolute', 'ak5PFResidual')

    process.ak5PFJEC = cms.ESProducer('JetCorrectionESChain',
        correctors = corrections
    )

    process.tree = cms.EDAnalyzer("MonoJetTreeMaker",
        pileup = cms.InputTag("addPileupInfo"),
        vertices = cms.InputTag("goodVertices"),
        muons = cms.InputTag("patMuons"),
        electrons = cms.InputTag("patElectrons"),
        taus = cms.InputTag("patTaus"),
        jets = cms.InputTag("patJets"),
        gens = cms.InputTag("genParticles"),
        pfmet = cms.InputTag("pfMet"),
        type1pfmet = cms.InputTag("pfType1CorrectedMet"),
        calomet = cms.InputTag("met"),
        triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
        jec = cms.string("ak5PFJEC"),
        weight = cms.double(1.0),
        isControlSample = cms.bool(False)
    )

    process.treePath = cms.Path(process.goodVertices + process.tree)
