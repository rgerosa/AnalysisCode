import FWCore.ParameterSet.Config as cms
import os,sys
### CMSSW command line parameter parser                                                                                                                                                               
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

# data or MC options                                                                                                                                                                                  
options.register (
    'outputFileName',"",VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'outputFileName, should have extension .root');

options.register (
    'isAMCNLO',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'Tell whether the generation is at NLO');

options.register (
    'partonMultiplicity',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,
    'if > 0 load the MLM matching for LO or FXFX for NLO');

options.register (
    'hadronicDecaysVBoson',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'if true force hadronic decays of the vector bosons');

options.register (
    'partonInBorn',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,
    'number of partons at Born Level');

## parsing command line arguments                                                                                                                                                                     
options.parseArguments()

if options.partonInBorn > 0 and options.partonMultiplicity:
    sys.exit("partonMultiplicity and partonInBorn cannot be > 0 at the same time: oen is used for MLM or FxFx matching, the other for fixed order NLO shower --> exit");

if not options.outputFileName.find(".root"):
    sys.exit("outputFileName should have .root extension --> exit");

sampleName = options.outputFileName.replace(".root","");

###########
from Configuration.StandardSequences.Eras import eras
process = cms.Process('GEN',eras.Run2_25ns)

# import of standard configurations                                                                                                                                                                    
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

# Message Logger settings                                                                                                                                                                           
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

# Input source from LHE file                                                                                                                                                                      
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(False),
    wantSummary = cms.untracked.bool(True)
)

# Production Info                                                                                                                                                                                   
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('alpha'),
    annotation = cms.untracked.string('LHE Input'),
    name = cms.untracked.string(sampleName)
)


# Output definition                                                                                                                                                                                
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outputFileName),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),                                     
    outputCommands = cms.untracked.vstring("drop *",
                                           "keep GenEventInfoProduct_generator_*_*",
                                           "keep *GenParticle*_genParticles_*_*",
                                           "keep *_source_*_*",
                                           "keep *GenMET*_genMetTrue_*_*",
                                           "keep *GenJet*_ak4GenJetsNoNu_*_*",
                                           "keep *GenLumiInfoProduct*_*_*_*",
                                           "keep *GenFilterInfo*_*_*_*",
                                           )
                                  )


# Other statements                                                                                                                                                                                    
process.genstepfilter.triggerConditions = cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = "MCRUN2_71_V1"


if not options.isAMCNLO and options.partonMultiplicity <= 0:
    process.load("Configuration.Generator.Hadronizer_TuneCUETP8M1_13TeV_generic_LHE_pythia8_cff");
elif not options.isAMCNLO and options.partonMultiplicity > 0:
    process.load("Configuration.Generator.Hadronizer_TuneCUETP8M1_13TeV_MLM_5f_max4j_LHE_pythia8_cff");
    process.generator.PythiaParameters.JetMatchingParameters.remove('JetMatching:nJetMax = 4');
    process.generator.PythiaParameters.JetMatchingParameters += cms.vstring('JetMatching:nJetMax = '+str(options.partonMultiplicity));
elif options.isAMCNLO and options.partonMultiplicity <= 0:
    process.load('Configuration.Generator.Pythia8aMCatNLOSettings_cfi');
    from Configuration.Generator.Hadronizer_TuneCUETP8M1_13TeV_generic_LHE_pythia8_cff import generator, pythia8CommonSettingsBlock, pythia8CUEP8M1SettingsBlock
    process.generator = generator.clone( PythiaParameters = cms.PSet(
            pythia8CommonSettingsBlock,
            pythia8CUEP8M1SettingsBlock,
            process.pythia8aMCatNLOSettingsBlock,
            parameterSets = cms.vstring('pythia8CommonSettings',
                                        'pythia8CUEP8M1Settings',
                                        'Pythia8aMCatNLOSettings',
                                        'TimeShower:nPartonsInBorn = '+str(options.partonInBorn))
            )
                                         )
elif options.isAMCNLO and options.partonMultiplicity > 0:    
    process.load("Configuration.Generator.Hadronizer_TuneCUETP8M1_13TeV_aMCatNLO_FXFX_5f_max2j_max1p_LHE_pythia8_cff");
    process.generator.PythiaParameters.processParameters.remove('JetMatching:nJetMax = 2');
    process.generator.PythiaParameters.processParameters += cms.vstring('JetMatching:nJetMax = '+str(options.partonMultiplicity));

if options.hadronicDecaysVBoson:
    process.generator.PythiaParameters.processParameters += cms.vstring('24:onMode  = off', 
                                                                        '24:onIfAny = 1 2 3 4 5 -1 -2 -3 -4 -5',
                                                                        '23:onMode  = off',
                                                                        '23:onIfAny = 1 2 3 4 5');
    
    

# Path and EndPath definitions                                                                                                                                                                        
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.output_step = cms.EndPath(process.output)

# Schedule definition                                                                                                                                                                                 
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.output_step)

# filter all path with the production filter sequence                                                                                                                                                
for path in process.paths:
    getattr(process,path)._seq = process.generator * getattr(process,path)._seq


from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils                                                                                                        
process = addMonitoring(process)

processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
