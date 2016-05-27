import FWCore.ParameterSet.Config as cms

### CMSSW command line parameter parser                                                                                                                                       
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

## data or MC options                                                                                                                                                        
options.register (
        'isMC',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'flag to indicate data or MC');

options.register (
        'globalTag','80X_dataRun2_Prompt_v8',VarParsing.multiplicity.singleton,VarParsing.varType.string,
        'gloabl tag to be uses');


## parsing command line arguments                                                                                                                                            
options.parseArguments()

if options.isMC and 'dataRun2' in options.globalTag:
        options.globalTag = '80X_mcRun2_asymptotic_2016_miniAODv2';

print "##### Settings ######"
print "Running with isMC                = ",options.isMC
print "Running with globalTag           = ",options.globalTag

# Define the CMSSW process
process = cms.Process("TNP")

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
    wantSummary = cms.untracked.bool(True) 
)

# How many events to process
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(options.maxEvents)
)

# Define the input source
if options.inputFiles == []:

        process.source = cms.Source("PoolSource",
                 fileNames = cms.untracked.vstring())

        if options.isMC:
            process.source.fileNames.append(
                '/store/mc/RunIISpring16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext1-v1/50000/E2654E59-4D1C-E611-922C-0CC47A713A04.root'
                )
        else:
            process.source.fileNames.append(
                '/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/554/00000/AA246637-E61F-E611-A971-02163E01187E.root'
                )
            
else:
    process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(options.inputFiles))


# Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", fileName = cms.string("tnptree.root"))

# Set the global tag depending on the sample type
process.GlobalTag.globaltag = options.globalTag

# Select good primary vertices
process.goodVertices = cms.EDFilter("VertexSelector",
				    src = cms.InputTag("offlineSlimmedPrimaryVertices"),
				    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
				    filter = cms.bool(True)
				    )

# Probe muons --> definition ->even looser than the loose muon ID
process.probemuons = cms.EDFilter("PATMuonSelector",
				  src = cms.InputTag("slimmedMuons"),
				  cut = cms.string("pt > 10 && abs(eta) < 2.4 && (isTrackerMuon || isStandAloneMuon)"),
				  filter = cms.bool(True)  
				  )

# Probe electrons
process.probeelectrons = cms.EDFilter("PATElectronSelector",
				      src = cms.InputTag("slimmedElectrons"),
				      cut = cms.string("pt > 10 && abs(eta) < 2.5"),
				      filter = cms.bool(True)  
				      )

# Electron and Photon ValueMaps for identification
from AnalysisCode.MonoXAnalysis.ElectronTools_cff import ElectronTools
ElectronTools(process,False,options.isMC)
### change the input for the electron ID into the probe electorns
process.egmGsfElectronIDs.physicsObjectSrc = "probeelectrons"


##### set of single muons triggers --> no eta restricted path are neeeded 
tagmuontriggernames = cms.vstring([
		"HLT_IsoMu20_v*",
		"HLT_IsoMu22_v*",
		"HLT_IsoMu24_v*",
		"HLT_IsoTkMu20*",
		"HLT_IsoTkMu22*",
		"HLT_IsoTkMu24*"
])

tagelectrontriggernames = cms.vstring([
        "HLT_Ele23_WPLoose_Gsf_v*",
        "HLT_Ele27_WPLoose_Gsf_v*",
        "HLT_Ele27_WP85_Gsf_v*",
        "HLT_Ele25_WPTight_Gsf_v*",
        "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
        "HLT_Ele115_CaloIdVT_GsfTrkIdT_v*",
	"HLT_Ele23_CaloIdL_TrackIdL_IsoVL_v*"
    ])    

# Make the ValueMap for muon tight ID -- cannot pass it through a string selection due to the vertex argument
process.probeinfo = cms.EDProducer("LeptonTnPInfoProducer",
    muons     = cms.InputTag("probemuons"),     ## input muons
    electrons = cms.InputTag("probeelectrons"), ## input electrons
    geninfo   = cms.InputTag("generator"),
    vertices  = cms.InputTag("goodVertices"), 
    triggerobjects = cms.InputTag("selectedPatTrigger"),
    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),				   
    #### Muon information for identification --> pt cut and matching with trigger info				   
    loosemuisocut  = cms.double(0.25),
    tightmuisocut  = cms.double(0.15),
    tagmuonptcut   = cms.double(22),
    tagmuonetacut  = cms.double(2.4),
    tagmuontrigmatchdR = cms.double(0.3),
    requiremuonhlt = cms.bool(True),
    tagmuontriggers = tagmuontriggernames,
    #### Electron information for identification --> pt cut and matching with trigger info
    tagelectronptcut   = cms.double(35),
    tagelectronetacut  = cms.double(2.5),
    tagelectrontrigmatchdR = cms.double(0.3),
    requireelectronhlt = cms.bool(True),
    tagelectrontriggers = tagelectrontriggernames,
    electronvetoid   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    electronlooseid  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    electronmediumid = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    electrontightid  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
)

# Tag muons --> filter on the collection content --> at least one
process.tagmuons = cms.EDFilter("PATMuonSelector", 
    src = cms.InputTag("probeinfo", "tightmuons"),
    cut = cms.string(""),
    filter = cms.bool(True) 
)

# Tag electrons --> filter on the collection content --> at least one
process.tagelectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("probeinfo", "tightelectrons"),
    cut = cms.string(""),
    filter = cms.bool(True)
)

# Tag and Probe pairs --> invariant mass and charge selection for electorns and muons
process.muontnp = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagmuons@+ probemuons@-"),
    cut   = cms.string("60 < mass < 120 & charge=0"),
)

process.electrontnp = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagelectrons@+ probeelectrons@-"),
    cut   = cms.string("60 < mass < 120 & charge=0"),
)

# Make the TnP tree --> official tag and probe analyzer
process.muontnptree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("muontnp"), ## pairs for Z->ll candidate
    arbitration   = cms.string("BestMass"),
    massForArbitration = cms.double(91.186),
    variables = cms.PSet(
        pt   = cms.string("pt()"),
        eta  = cms.string("eta()"),
        abseta = cms.string("abs(eta())"),
        phi  = cms.string("phi()"),
        nvtx = cms.InputTag("probeinfo", "munvtxmap"), ## store nvtx
        wgt  = cms.InputTag("probeinfo", "muwgtmap")   ## store gen weight
    ),
    flags = cms.PSet(
        pfid      = cms.string("isPFMuon"), ## if is a PF muon
        hltmu20   = cms.InputTag("probeinfo", "hltmu20muonrefs"),   ## if it belongs to hltmu20
        hlttkmu20 = cms.InputTag("probeinfo", "hlttkmu20muonrefs"), ## if it belongs to isotk20
        looseid   = cms.InputTag("probeinfo", "loosemuonrefs"),     ## if pass the loose
        tightid   = cms.InputTag("probeinfo", "tightmuonrefs"),     ## if pass the tight
    ),
    isMC = cms.bool(options.isMC)
)

process.electrontnptree = cms.EDAnalyzer("TagProbeFitTreeProducer",
    tagProbePairs = cms.InputTag("electrontnp"),
    arbitration   = cms.string("BestMass"),
    massForArbitration = cms.double(91.186),
    variables = cms.PSet(
        pt  = cms.string("pt()"),
        eta = cms.string("superCluster().eta()"),
        abseta = cms.string("abs(superCluster().eta())"),
        phi = cms.string("phi()"),
        nvtx = cms.InputTag("probeinfo", "elnvtxmap"),
        wgt  = cms.InputTag("probeinfo", "elwgtmap")
    ),
    flags = cms.PSet(
        vetoid   = cms.InputTag("probeinfo", "vetoelectronrefs"),
        looseid  = cms.InputTag("probeinfo", "looseelectronrefs"),
        mediumid = cms.InputTag("probeinfo", "mediumelectronrefs"),
        tightid  = cms.InputTag("probeinfo", "tightelectronrefs"),
    ),
    isMC = cms.bool(options.isMC)
)

# MC Truth Matching
if options.isMC : 
    process.mcmutagmatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
        matchPDGId = cms.vint32(13),
        src     = cms.InputTag("tagmuons"),
        distMin = cms.double(0.1),
        matched = cms.InputTag("prunedGenParticles"),
        checkCharge = cms.bool(True)
    )

    process.mcmuprobematch = cms.EDProducer("MCTruthDeltaRMatcherNew",
        matchPDGId = cms.vint32(13),
        src        = cms.InputTag("probemuons"),
        distMin    = cms.double(0.1),
        matched    = cms.InputTag("prunedGenParticles"),
        checkCharge = cms.bool(True)
    )

    process.mceltagmatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
        matchPDGId  = cms.vint32(11),
        src         = cms.InputTag("tagelectrons"),
        distMin     = cms.double(0.1),
        matched     = cms.InputTag("prunedGenParticles"),
        checkCharge = cms.bool(True)
    )

    process.mcelprobematch = cms.EDProducer("MCTruthDeltaRMatcherNew",
        matchPDGId  = cms.vint32(11),
        src         = cms.InputTag("probeelectrons"),
        distMin     = cms.double(0.1),
        matched     = cms.InputTag("prunedGenParticles"),
        checkCharge = cms.bool(True)
    )

    #### muon gen info
    process.muontnptree.tagMatches   = cms.InputTag("mcmutagmatch")
    process.muontnptree.probeMatches = cms.InputTag("mcmuprobematch")
    process.muontnptree.motherPdgId  = cms.int32(23)
    process.muontnptree.makeMCUnbiasTree = cms.bool(True)
    process.muontnptree.checkMotherInUnbiasEff = cms.bool(True)
    process.muontnptree.allProbes = cms.InputTag("probemuons")
    
    process.electrontnptree.tagMatches = cms.InputTag("mceltagmatch")
    process.electrontnptree.probeMatches = cms.InputTag("mcelprobematch")
    process.electrontnptree.motherPdgId = cms.int32(23)
    process.electrontnptree.makeMCUnbiasTree = cms.bool(True)
    process.electrontnptree.checkMotherInUnbiasEff = cms.bool(True)
    process.electrontnptree.allProbes = cms.InputTag("probeelectrons")

process.mutreePath = cms.Path(process.muontnptree)
process.eltreePath = cms.Path(process.electrontnptree)
