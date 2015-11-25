import FWCore.ParameterSet.Config as cms

# Define the CMSSW process
process = cms.Process("TREE")

# Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False) 
)

# How many events to process
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1)
)

# Is this a simulation or real data
isMC = False

# Filter on high MET events
filterHighMETEvents = False

# Filter on triggered events
filterOnHLT = True

# Use private JECs since the GTs are not updated
usePrivateSQlite = True

# Apply L2L3 residual corrections
applyL2L3Residuals = True

# Process name used in MiniAOD -- needed to get the correct trigger results, and also for redoing the MET
miniAODProcess = "RECO"

# Define the input source
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring([
        '/store/data/Run2015D/MET/MINIAOD/PromptReco-v4/000/258/750/00000/5EE58B11-7572-E511-B952-02163E014378.root'
    ])
)

# Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", fileName = cms.string("tree.root"))

# Set the global tag depending on the sample type
from Configuration.AlCa.GlobalTag import GlobalTag
if isMC:
    process.GlobalTag.globaltag = '74X_mcRun2_asymptotic_v2'   # for Simulation
else:
    process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v4'            # for Data

# Setup the private SQLite -- Ripped from PhysicsTools/PatAlgos/test/corMETFromMiniAOD.py
if usePrivateSQlite:
    from CondCore.DBCommon.CondDBSetup_cfi import *
    import os
    era = "Summer15_25nsV6"
    if isMC : 
        era += "_MC"
    else :
        era += "_DATA"
    dBFile = os.path.expandvars(era+".db")
    process.jec = cms.ESSource("PoolDBESSource",
        CondDBSetup,
        connect = cms.string("sqlite_file:"+dBFile),
        toGet =  cms.VPSet(
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PF"),
                label= cms.untracked.string("AK4PF")
            ),
            cms.PSet(
                record = cms.string("JetCorrectionsRecord"),
                tag = cms.string("JetCorrectorParametersCollection_"+era+"_AK4PFchs"),
                label= cms.untracked.string("AK4PFchs")
            ),
        )
    )
    process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')


# JEC levels when redoing jets and MET
if isMC:
    JECLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
else :
    if not applyL2L3Residuals : 
        JECLevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
    else : 
        JECLevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']

# Re-run the HBHE Noise filter
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

# Select good primary vertices
process.goodVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True)
)

# Electron and Photon ValueMaps for identification
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer(process, dataFormat)

ele_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']
ph_id_modules  = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_50ns_V1_cff']
for idmod in ele_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
for idmod in ph_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)


# Re-running of jets and MET
process.load("JetMETCorrections.Configuration.JetCorrectors_cff")
    
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
process.patJetCorrFactorsReapplyJEC = patJetCorrFactorsUpdated.clone(
    src = cms.InputTag("slimmedJets"),
    levels = JECLevels,
    payload = 'AK4PFchs' 
)

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
process.slimmedJetsRecorrected = patJetsUpdated.clone(
    jetSource = cms.InputTag("slimmedJets"),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
)

process.pfCandsCHS = cms.EDFilter("CandPtrSelector", 
    src = cms.InputTag("packedPFCandidates"), 
    cut = cms.string("fromPV")
)

from RecoJets.Configuration.RecoPFJets_cff import ak4PFJetsCHS
process.ak4PFJetsCHS = ak4PFJetsCHS.clone(src = "pfCandsCHS")

from RecoMET.METProducers.PFMET_cfi import pfMet
process.pfMet = pfMet.clone(src = "packedPFCandidates")

process.load("JetMETCorrections.Type1MET.correctionTermsPfMetType1Type2_cff")
if isMC : 
    process.corrPfMetType1.jetCorrLabelRes = "ak4PFCHSL1FastL2L3Corrector"
else :
    if not applyL2L3Residuals :
        process.corrPfMetType1.jetCorrLabelRes = "ak4PFCHSL1FastL2L3Corrector"
    else :
        process.corrPfMetType1.jetCorrLabelRes = "ak4PFCHSL1FastL2L3ResidualCorrector"

from JetMETCorrections.Type1MET.correctedMet_cff import pfMetT1
process.pfMetT1 = pfMetT1.clone()

# Create a set of objects to read from
process.selectedObjects = cms.EDProducer("PFCleaner",
    vertices = cms.InputTag("goodVertices"),
    pfcands = cms.InputTag("packedPFCandidates"),
    muons = cms.InputTag("slimmedMuons"),
    electrons = cms.InputTag("slimmedElectrons"),
    photons = cms.InputTag("slimmedPhotons"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    jets = cms.InputTag("slimmedJetsRecorrected"),
    electronidveto = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    electronidmedium = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    photonidloose = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-loose"),
    photonsieie = cms.InputTag("photonIDValueMapProducer", "phoFull5x5SigmaIEtaIEta"),
    photonphiso = cms.InputTag("photonIDValueMapProducer", "phoPhotonIsolation"),
    photonchiso = cms.InputTag("photonIDValueMapProducer", "phoChargedIsolation")
)

# Define all the METs corrected for lepton/photon momenta
process.mumet = cms.EDProducer("MuonCorrectedRecoMETProducer",
    met = cms.InputTag("pfMet"),
    muons = cms.InputTag("selectedObjects", "muons")
)
process.t1mumet = cms.EDProducer("MuonCorrectedRecoMETProducer",
    met = cms.InputTag("pfMetT1"),
    muons = cms.InputTag("selectedObjects", "muons")
)
process.elmet = cms.EDProducer("CandCorrectedRecoMETProducer",
    met = cms.InputTag("pfMet"),
    cands = cms.VInputTag(cms.InputTag("selectedObjects", "electrons"))
)
process.t1elmet = cms.EDProducer("CandCorrectedRecoMETProducer",
    met = cms.InputTag("pfMetT1"),
    cands = cms.VInputTag(cms.InputTag("selectedObjects", "electrons")),
)
process.phmet = cms.EDProducer("CandCorrectedRecoMETProducer",
    met = cms.InputTag("pfMet"),
    cands = cms.VInputTag(cms.InputTag("selectedObjects", "photons"))

)
process.t1phmet = cms.EDProducer("CandCorrectedRecoMETProducer",
    met = cms.InputTag("pfMetT1"),
    cands = cms.VInputTag(cms.InputTag("selectedObjects", "photons")),
)

# Make the tree 
process.tree = cms.EDAnalyzer("MonoJetTreeMaker",
    pileup = cms.InputTag("addPileupInfo"),
    genevt = cms.InputTag("generator"),
    vertices = cms.InputTag("goodVertices"),
    gens = cms.InputTag("prunedGenParticles"),
    muons = cms.InputTag("selectedObjects", "muons"),
    electrons = cms.InputTag("selectedObjects", "electrons"),
    photons = cms.InputTag("selectedObjects", "photons"),
    tightmuons = cms.InputTag("selectedObjects", "tightmuons"),
    tightelectrons = cms.InputTag("selectedObjects", "tightelectrons"),
    tightphotons = cms.InputTag("selectedObjects", "tightphotons"),
    electronLooseId = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    photonLooseId = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-loose"),
    photonMediumId = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-medium"),
    photonTightId = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-tight"),
    photonHighPtId = cms.InputTag("selectedObjects", "photonHighPtId"),
    taus = cms.InputTag("slimmedTaus"),
    jets = cms.InputTag("slimmedJetsRecorrected"),
    pfmet = cms.InputTag("pfMet"),
    t1pfmet = cms.InputTag("pfMetT1"),
    mumet = cms.InputTag("mumet"),
    t1mumet = cms.InputTag("t1mumet"),
    elmet = cms.InputTag("elmet"),
    t1elmet = cms.InputTag("t1elmet"),
    phmet = cms.InputTag("phmet"),
    t1phmet = cms.InputTag("t1phmet"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
    filterResults = cms.InputTag("TriggerResults", "", miniAODProcess),
    hbheloose = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Loose"),
    hbhetight = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
    hbheiso   = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),
    xsec = cms.double(0.001),
    cleanMuonJet = cms.bool(True),
    cleanElectronJet = cms.bool(True),
    cleanPhotonJet = cms.bool(True),
    applyHLTFilter = cms.bool(filterOnHLT),
    uselheweights = cms.bool(False),
    isWorZMCSample = cms.bool(False)
)

# Tree for the generator weights
process.gentree = cms.EDAnalyzer("LHEWeightsTreeMaker",
    lheinfo = cms.InputTag("externalLHEProducer"),
    geninfo = cms.InputTag("generator"),
    uselheweights = cms.bool(False),
    addqcdpdfweights = cms.bool(False)
)

# MET filter
process.metfilter = cms.EDFilter("CandViewSelector",
    src = cms.InputTag("t1mumet"),
    cut = cms.string("et > 200"),
    filter = cms.bool(True)
)

# Set up the path
if filterHighMETEvents: 
    if (isMC):
        process.treePath = cms.Path(process.gentree + process.goodVertices + process.metfilter + process.tree)
    else :
        process.treePath = cms.Path(                  process.goodVertices + process.metfilter + process.tree)
else :
    if (isMC):
        process.treePath = cms.Path(process.gentree + process.goodVertices                     + process.tree)
    else :
        process.treePath = cms.Path(                  process.goodVertices                     + process.tree)
