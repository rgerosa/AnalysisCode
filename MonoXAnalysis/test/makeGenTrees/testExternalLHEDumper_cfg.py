import FWCore.ParameterSet.Config as cms

### CMSSW command line parameter parser                                                                                                                                         
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.register (
        'outputFileName','file.lhe',VarParsing.multiplicity.singleton,VarParsing.varType.string,
        'output file name created by cmsRun');

options.parseArguments()

process = cms.Process("dumLHE")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/RunIIWinter15wmLHE/VectorMonoW_Mphi-500_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph/LHE/MCRUN2_71_V1-v4/10000/04E93F09-696F-E511-A316-B083FECF8ACE.root',
                                      '/store/mc/RunIIWinter15wmLHE/VectorMonoW_Mphi-500_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph/LHE/MCRUN2_71_V1-v4/10000/50F00C3E-7D6F-E511-A4F8-001E6757EAA4.root',
                                      '/store/mc/RunIIWinter15wmLHE/VectorMonoW_Mphi-500_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph/LHE/MCRUN2_71_V1-v4/10000/82239C23-9E6F-E511-BF2D-00304865C40C.root'),
    processingMode = cms.untracked.string('Runs'),
)

process.externalLHEAsciiDumper = cms.EDAnalyzer('LHEInfoReader',
                                                outputLHEFileName = cms.string(options.outputFileName))

process.p = cms.Path(process.externalLHEAsciiDumper)
