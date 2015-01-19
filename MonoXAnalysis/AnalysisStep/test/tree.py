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
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

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

# Filter on high MET events
filterHighMETEvents = True

# Define the global tag depending on the sample type
if isMC:
    process.GlobalTag.globaltag = 'PHYS14_25_V1::All'   # for Simulation
else:
    process.GlobalTag.globaltag = 'PHYS14_25_V1::All'   # for Data (for now set to the MC tag)

# Define the input source
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring([
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/08405777-2560-E411-BE59-0025902009B8.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/929341DB-2660-E411-8D11-002590200900.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/DC071E96-DF60-E411-A43C-001E67397D91.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/F8D96877-2560-E411-9418-0025902009B8.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/1E5B3288-DA60-E411-9891-002590A371C4.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/30B2D673-9061-E411-B97D-001E67396DF1.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/463CFD0B-9161-E411-B7E5-001E67398E12.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/6CC45722-9061-E411-A3B4-001E67396DF1.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/78ED8A43-0D60-E411-8B05-0025B3E05D8C.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/8498B043-0D60-E411-835E-002481E14F1E.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/9601DCE4-9161-E411-915D-001E67396D5B.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/B22166BD-0D60-E411-95E0-0025B3E0658E.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/DA1BF4BF-9061-E411-9BB8-001E67396DF1.root',
       '/store/mc/Phys14DR/DarkMatter_Monojet_M-10_AV_Tune4C_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/20000/EEFD4E0F-9161-E411-8556-001E67396DF1.root'         
    ])
)

# Select good primary vertices
process.goodVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)

process.selectedObjects = cms.EDProducer("PFCleaner",
    src = cms.InputTag("particleFlow"),
    beamspot = cms.InputTag("offlineBeamSpot"),
    vertices = cms.InputTag("goodVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
    pfcands = cms.InputTag("packedPFCandidates"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    photons = cms.InputTag("slimmedPhotons"),
    conversions = cms.InputTag("reducedEgamma", "reducedConversions"),
    d0cut = cms.double(0.2),
    dzcut = cms.double(0.5)
)

# Define all the METs corrected for lepton/photon momenta
process.mumet = cms.EDProducer("CandCorrectedMETProducer",
    met = cms.InputTag("slimmedMETs"),
    cands = cms.VInputTag(cms.InputTag("selectedObjects", "muons")) 
)

process.phmet = cms.EDProducer("CandCorrectedMETProducer",
    met = cms.InputTag("slimmedMETs"),
    cands = cms.VInputTag(cms.InputTag("selectedObjects", "photons")) 
)

from RecoMET.METProducers.PFMET_cfi import pfMet

process.pfMet = pfMet.clone(src = "packedPFCandidates")
process.pfMet.calculateSignificance = False

# Remove MC matching when running on data
if not isMC:
    runOnData(process)

process.TFileService = cms.Service("TFileService", fileName = cms.string("tree.root"))

process.tree = cms.EDAnalyzer("MonoJetTreeMaker",
    pileup = cms.InputTag("addPileupInfo"),
    vertices = cms.InputTag("goodVertices"),
    gens = cms.InputTag("prunedGenParticles"),
    photons = cms.InputTag("selectedObjects", "photons"),
    muons = cms.InputTag("selectedObjects", "muons"),
    electrons = cms.InputTag("selectedObjects", "electrons"),
    electronsnew = cms.InputTag("selectedObjects", "electronsnew"),
    tightmuons = cms.InputTag("selectedObjects", "tightmuons"),
    tightelectrons = cms.InputTag("selectedObjects", "tightelectrons"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJets"),
    fatjets = cms.InputTag("slimmedJetsAK8"),
    pfmet = cms.InputTag("pfMet"),
    t1pfmet = cms.InputTag("slimmedMETs"),
    mumet = cms.InputTag("mumet"),
    phmet = cms.InputTag("phmet"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    weight = cms.double(1000.0/191200.0),
    isWorZMCSample = cms.bool(False),
    isSignalSample = cms.bool(True)
)

# Select events passing the monojet triggers and having MET > 200 GeV
process.triggerfilter = cms.EDFilter("HLTCheckFilter",
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    triggerPaths   = cms.vstring('HLT_MonoCentralPFJet140_PFMETNoMu100_PFMHTNoMu140_NoiseCleaned', 'HLT_MonoCentralPFJet140_PFMETNoMu140_PFMHTNoMu140_NoiseCleaned', 'HLT_MonoCentralPFJet150_PFMETNoMu150_PFMHTNoMu150_NoiseCleaned')
)

process.metfilter = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("mumet"),
    cut = cms.string("et > 200"),
    filter = cms.bool(True)
)


if filterHighMETEvents: 
    #process.treePath = cms.Path(process.goodVertices + process.triggerfilter + process.metfilter + process.tree)
    process.treePath = cms.Path(process.goodVertices + process.metfilter + process.tree)
else :
    process.treePath = cms.Path(process.goodVertices + process.tree)

