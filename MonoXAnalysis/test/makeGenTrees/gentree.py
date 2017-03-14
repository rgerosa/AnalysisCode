### CMSSW command line parameter parser                                                                                                                                                               
from FWCore.ParameterSet.VarParsing import VarParsing
import os
options = VarParsing ('python')

## data or MC options                                                                                                                                                                                 
options.register (
        'outputName',"gentree.root",VarParsing.multiplicity.singleton,VarParsing.varType.string,
        'name for the outputfile');

options.register(
        'crossSection',-1.,VarParsing.multiplicity.singleton, VarParsing.varType.float,
        'external value for sample cross section, in case of data it is fixed to 0.001');

options.register(
        'sample',23,VarParsing.multiplicity.singleton, VarParsing.varType.int,
        'to identify which kind of sample');

options.register(
        'isMiniAOD',True,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
        'flag to tell if one is running on miniAOD or GEN files');

options.register(
        'isSignalSample',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
        'in case one want to run on signal files at gen level');

options.register(
        'minBosonPt',70,VarParsing.multiplicity.singleton, VarParsing.varType.float,
        'apply a minimum boson pt requirement, which is DM mediator in case of DM signals');
   
options.parseArguments()

# Define the CMSSW process
import FWCore.ParameterSet.Config as cms
process = cms.Process("KFAC")

# Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

# Define the input source
if options.inputFiles == []:
    process.source = cms.Source("PoolSource", 
                                fileNames = cms.untracked.vstring(
            '/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/620000/0429E5A7-80C3-E611-B89C-0242AC130003.root'))
else:
    process.source = cms.Source("PoolSource", 
                                fileNames = cms.untracked.vstring(options.inputFiles))

#process.source.eventsToProcess = cms.untracked.VEventRange('1:1450:1672340')
    
# Output file
process.TFileService = cms.Service("TFileService", 
    fileName = cms.string(options.outputName)
)

# Processing setup
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(options.maxEvents)
)

# Make the tree 
process.gentree = cms.EDAnalyzer("GenTreeMaker",
                                 lherun = cms.InputTag("externalLHEProducer"),
                                 lheevt = cms.InputTag("externalLHEProducer"),
                                 genevt = cms.InputTag("generator"),
                                 gens   = cms.InputTag("prunedGenParticles"),
                                 xsec   = cms.double(options.crossSection),   
                                 jets   = cms.InputTag("slimmedGenJets"),
                                 met    = cms.InputTag("slimmedMETs"),
                                 sample = cms.int32(options.sample),
                                 isSignalSample = cms.bool(options.isSignalSample),
                                 isMiniAOD = cms.bool(options.isMiniAOD),
                                 minBosonPt = cms.double(options.minBosonPt)                                 
                              )

if not options.isMiniAOD:
    process.gentree.gens   = cms.InputTag("genParticles");
    process.gentree.jets   = cms.InputTag("ak4GenJetsNoNu");
    process.gentree.met    = cms.InputTag("genMetTrue");
    process.gentree.lheevt = cms.InputTag("source");

process.gentreePath = cms.Path(process.gentree)

processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
