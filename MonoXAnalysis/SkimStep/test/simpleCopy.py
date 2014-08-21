import FWCore.ParameterSet.Config as cms

# Define the CMSSW process
process = cms.Process("COPY")

# Message Logger Stuff
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 10

# Display summary at the end of the job
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# How many events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

# Define the input source
process.source = cms.Source("PoolSource", 
     fileNames = cms.untracked.vstring('/store/mc/Summer12_DR53X/WToENu_MSUB166_Pt-100to300_TuneZ2Star_8TeV-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/08F3738E-31ED-E111-BE72-001EC9D2887E.root'),
)

# Define the output
process.out = cms.OutputModule("PoolOutputModule", 
    fileName = cms.untracked.string('monoxCopy.root') 
)

# Schedule 
process.outpath = cms.EndPath(process.out)
process.schedule = cms.Schedule(process.outpath)
