# Define the CMSSW process
import FWCore.ParameterSet.Config as cms
process = cms.Process("KFAC")

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

# Define the input source
process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring('/store/mc/RunIIFall15MiniAODv2/ZJetsToNuNu_HT-100To200_13TeV-madgraph/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/060FC9A4-C8BD-E511-B138-000F530E46D0.root')
)

# Output file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string("gentree.root")
)

# Processing setup
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1)
)

# Make the tree 
process.gentree = cms.EDAnalyzer("GenTreeMaker",
                                 genevt = cms.InputTag("generator"),
                                 gens   = cms.InputTag("prunedGenParticles"),
                                 xsec   = cms.double(55.38*3.0),   
                                 jets   = cms.InputTag("slimmedGenJets"),
                                 met    = cms.InputTag("genMetTrue"),
                                 sample = cms.int32(23)
                                 )

process.gentreePath = cms.Path(process.gentree)
