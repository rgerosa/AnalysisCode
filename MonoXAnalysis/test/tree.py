import os
import FWCore.ParameterSet.Config as cms
  
### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

## data or MC
options.register (
	'isMC',False,
	VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate data or MC');

## filter or not high MET events
options.register (
	'filterHighMETEvents',False,
	VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate if apply or not MET filters');

## filter or not using HLT trigger path
options.register (
	'filterOnHLT',True,
	VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate if apply or not trigger requirements');

## private SQL file for JEC
options.register (
	'usePrivateSQlite',True,
    VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'if a private SQL file with JEC to be found in test directory');

## apply or not L2L3 Residual for data
options.register (
	'applyL2L3Residuals',True,
	VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'apply or not L2L3 Residual JEC on data');

## re-compute puppi MET with dedicated JEC
options.register (
	'doMETSystematics',True,
	VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'recompute Puppi MET propagating JEC from Jet + systematics');

## re-compute pileup-jet id for AK4 jets
options.register (
	'addPileupJetID',False,
	VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	're-compute pileup-jet id for AK4 jets');

options.register (
	'addQGLikelihood',True,
	VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'compute QGLikelihood for AK4 jets');

options.register (
	'addMVAMet',False,
	VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'compute MVAMet');
  	
## processName
options.register (
	'processName','TREE',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'process name to be considered');
    
## miniAOD process name    
options.register (
	'miniAODProcess','RECO',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'process name used for miniAOD production (target is miniAODv2)');

## outputFile Name
options.register (
	'outputFileName','tree.root',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'output file name created by cmsRun');

## GT to be used    
options.register (
	'globalTag','74X_dataRun2_Prompt_v4',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'gloabl tag to be uses');
  
## JEC    
options.register (
	'JECEra','Summer15_25nsV6',
    VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'JEC correction era');
                                    
## Debug options
options.register (
	'dropAnalyzerDumpEDM',False,
    VarParsing.multiplicity.singleton, VarParsing.varType.bool,
    'not run the analyzer and store an edm file');

options.register (
	'reportEvery',100,
    VarParsing.multiplicity.singleton, VarParsing.varType.int,
    'report message logger CMSSW');

options.register (
	'wantSummary',True,
     VarParsing.multiplicity.singleton, VarParsing.varType.bool, 
     'report message logger CMSSW');

## parsing command line arguments
options.parseArguments()

if options.isMC and 'dataRun2' in options.globalTag:
	options.globalTag = '74X_mcRun2_asymptotic_v2';

if options.isMC and options.applyL2L3Residuals:
	options.applyL2L3Residuals = False


print "##### Settings ######"
print "Running with isMC = ",options.isMC	
print "Running with filterHighMETEvents = ",options.filterHighMETEvents	
print "Running with filterOnHLT = ",options.filterOnHLT	
print "Running with usePrivateSQlite = ",options.usePrivateSQlite	
print "Running with applyL2L3Residuals = ",options.applyL2L3Residuals	
print "Running with doMETSystematics = ",options.doMETSystematics	
print "Running with processName = ",options.processName	
print "Running with miniAODProcess = ",options.miniAODProcess	
print "Running with outputFileName = ",options.outputFileName	
print "Running with globalTag = ",options.globalTag	
print "Running with JEC Era = ",options.JECEra	
print "Running with dropAnalyzerDumpEDM = ",options.dropAnalyzerDumpEDM	
print "Running with reportEvery = ",options.reportEvery	    
print "Running with wantSummary = ",options.wantSummary	
print "Running with addPileupJetID = ",options.addPileupJetID
print "Running with addQGLikelihood = ",options.addQGLikelihood
print "Running with addMVAMet = ",options.addMVAMet
print "#####################"

## Define the CMSSW process
process = cms.Process(options.processName)

## Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

## Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

## Define the input source
if options.inputFiles == []:

	process.source = cms.Source("PoolSource", 
   		 fileNames = cms.untracked.vstring()
   	)

	if not options.isMC :
		process.source.fileNames.append(
        	'/store/data/Run2015D/MET/MINIAOD/PromptReco-v4/000/258/750/00000/5EE58B11-7572-E511-B952-02163E014378.root'
    	)
	else:
		process.source.fileNames.append(     'root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/ZJetsToNuNu_HT-100To200_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/008902DD-9F6F-E511-BCE9-0025904C540C.root'
    	)    	
else:
   process.source = cms.Source("PoolSource",
   	  fileNames = cms.untracked.vstring(options.inputFiles)	
	)


## Set the process options -- Display summary at the end, enable unscheduled execution
process.options = cms.untracked.PSet( 
    allowUnscheduled = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(options.wantSummary) 
)

## How many events to process
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(options.maxEvents)
)

# Set the global tag depending on the sample type
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = options.globalTag  

## Setup the private SQLite -- Ripped from PhysicsTools/PatAlgos/test/corMETFromMiniAOD.py
from AnalysisCode.MonoXAnalysis.JECConfiguration_cff import JECConfiguration
## connect to a local SQLite file or take corrections from GT
JECConfiguration(process,options.usePrivateSQlite,options.JECEra,options.isMC,options.applyL2L3Residuals)

## Setup MET filters or not
process.load('AnalysisCode.MonoXAnalysis.METFilters_cff')

# run cut-based electron ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
from AnalysisCode.MonoXAnalysis.ElectronTools_cff import ElectronTools
ElectronTools(process)

# run cut-based photon ID 
from AnalysisCode.MonoXAnalysis.PhotonTools_cff import PhotonTools
PhotonTools(process)

# Apply JEC on jets and propagation on Missing ET
from AnalysisCode.MonoXAnalysis.JetMetCorrector_cff import JetMetCorrector
## apply JEC and propagation on MET for AK4PFchs
JetMetCorrector(process,"slimmedJets","slimmedMETs","AK4PFchs",options.isMC, options.applyL2L3Residuals, options.doMETSystematics)
## apply JEC and propagation on MET for AK4PFPuppi
JetMetCorrector(process,"slimmedJetsPuppi","slimmedMETsPuppi","AK4PFPuppi",options.isMC, options.applyL2L3Residuals, options.doMETSystematics)

# Create a set of objects to read from
process.selectedObjects = cms.EDProducer("PFCleaner",
     vertices  = cms.InputTag("goodVertices"),
     pfcands   = cms.InputTag("packedPFCandidates"),
     muons     = cms.InputTag("slimmedMuons"),
     electrons = cms.InputTag("slimmedElectrons"),
     photons   = cms.InputTag("slimmedPhotons"),
     rho       = cms.InputTag("fixedGridRhoFastjetAll"),
     jets      = cms.InputTag("slimmedJetsRecorrectedAK4PFchs"),
     jetsPuppi = cms.InputTag("slimmedJetsRecorrectedAK4PFPuppi"),
     electronidveto  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
     electronidtight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
     electronidheep  = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
     photonidloose = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
     photonsieie = cms.InputTag("photonIDValueMapProducer", "phoFull5x5SigmaIEtaIEta"),
     photonphiso = cms.InputTag("photonIDValueMapProducer", "phoPhotonIsolation"),
     photonchiso = cms.InputTag("photonIDValueMapProducer", "phoChargedIsolation")
)

## modify some existing jet collections adding pileup-jet id and QGLikelihood from GT
from AnalysisCode.MonoXAnalysis.JetTools_cff import addPileupJetID, addQGLikelihood

if options.addPileupJetID:  ## so far not working for slimmedJets, so re-clustering and update
	addPileupJetID(process, collection = "slimmedJetsRecorrectedAK4PFchs", postfix = "")
	addPileupJetID(process, collection = "slimmedJetsRecorrectedAK4PFPuppi", postfix = "Puppi")
	
if options.addQGLikelihood:

	inputJet = "slimmedJetsRecorrectedAK4PFchs";
	inputJetPuppi = "slimmedJetsRecorrectedAK4PFPuppi";
	if options.addPileupJetID:
		inputJet += "PUID";
		inputJetPuppi += "PUID";

	addQGLikelihood( process,collection = inputJet, postfix = "");
	addQGLikelihood( process,collection = inputJetPuppi, postfix = "Puppi");


## in case run the MVA met producer
from AnalysisCode.MonoXAnalysis.MVAMet_cff import runMVAMet

if options.addMVAMet:
	## to parse leptons we need a list of CandidateView not a value map with Refs
	leptons = ["PFCleaner:tightmuons","PFCleaner:tightelectrons"]
	runMVAMet(process,isMC = options.isMC,leptons = leptons )


# Define all the METs corrected for lepton/photon momenta
#process.mumet = cms.EDProducer("MuonCorrectedRecoMETProducer",
#    met = cms.InputTag("pfMet"),
#    muons = cms.InputTag("selectedObjects", "muons")
#)
#process.t1mumet = cms.EDProducer("MuonCorrectedRecoMETProducer",
#    met = cms.InputTag("pfMetT1"),
#    muons = cms.InputTag("selectedObjects", "muons")
#)
#process.elmet = cms.EDProducer("CandCorrectedRecoMETProducer",
#    met = cms.InputTag("pfMet"),
#    cands = cms.VInputTag(cms.InputTag("selectedObjects", "electrons"))
#)
#process.t1elmet = cms.EDProducer("CandCorrectedRecoMETProducer",
#    met = cms.InputTag("pfMetT1"),
#    cands = cms.VInputTag(cms.InputTag("selectedObjects", "electrons")),
#)
#process.phmet = cms.EDProducer("CandCorrectedRecoMETProducer",
#    met = cms.InputTag("pfMet"),
#    cands = cms.VInputTag(cms.InputTag("selectedObjects", "photons"))
#)
#process.t1phmet = cms.EDProducer("CandCorrectedRecoMETProducer",
#    met = cms.InputTag("pfMetT1"),
#    cands = cms.VInputTag(cms.InputTag("selectedObjects", "photons")),
#)

# Make the tree 
#process.tree = cms.EDAnalyzer("MonoJetTreeMaker",
#    pileup = cms.InputTag("addPileupInfo"),
#    genevt = cms.InputTag("generator"),
#    vertices = cms.InputTag("goodVertices"),
#    gens = cms.InputTag("prunedGenParticles"),
#    muons = cms.InputTag("selectedObjects", "muons"),
#    electrons = cms.InputTag("selectedObjects", "electrons"),
#    photons = cms.InputTag("selectedObjects", "photons"),
#    tightmuons = cms.InputTag("selectedObjects", "tightmuons"),
#    tightelectrons = cms.InputTag("selectedObjects", "tightelectrons"),
#    tightphotons = cms.InputTag("selectedObjects", "tightphotons"),
#    electronLooseId = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
#    photonLooseId = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-loose"),
#    photonMediumId = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-medium"),
#    photonTightId = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-tight"),
#    photonHighPtId = cms.InputTag("selectedObjects", "photonHighPtId"),
#    taus = cms.InputTag("slimmedTaus"),
#    jets = cms.InputTag("slimmedJetsRecorrected"),
#    pfmet = cms.InputTag("pfMet"),
#    t1pfmet = cms.InputTag("pfMetT1"),
#    mumet = cms.InputTag("mumet"),
#    t1mumet = cms.InputTag("t1mumet"),
#    elmet = cms.InputTag("elmet"),
#    t1elmet = cms.InputTag("t1elmet"),
#    phmet = cms.InputTag("phmet"),
#    t1phmet = cms.InputTag("t1phmet"),
#    triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
#    filterResults = cms.InputTag("TriggerResults", "", miniAODProcess),
#    hbheloose = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Loose"),
#    hbhetight = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
#    hbheiso   = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),
#    xsec = cms.double(0.001),
#    cleanMuonJet = cms.bool(True),
#    cleanElectronJet = cms.bool(True),
#    cleanPhotonJet = cms.bool(True),
#    applyHLTFilter = cms.bool(filterOnHLT),
#    uselheweights = cms.bool(False),
#    isWorZMCSample = cms.bool(False)
#)

# Tree for the generator weights
#process.gentree = cms.EDAnalyzer("LHEWeightsTreeMaker",
#    lheinfo = cms.InputTag("externalLHEProducer"),
#    geninfo = cms.InputTag("generator"),
#    uselheweights = cms.bool(False),
#    addqcdpdfweights = cms.bool(False)
#)

# MET filter
#process.metfilter = cms.EDFilter("CandViewSelector",
#    src = cms.InputTag("t1mumet"),
#    cut = cms.string("et > 200"),
#    filter = cms.bool(True)
#)

# Set up the path
#if filterHighMETEvents: 
#    if (isMC):
#        process.treePath = cms.Path(process.gentree + process.metFilters + process.metfilter + process.tree)
#    else :
#        process.treePath = cms.Path(                  process.metFilters + process.metfilter + process.tree)
#else :
#    if (isMC):
#       process.treePath = cms.Path(process.gentree + process.metFilters                     + process.tree)
#   else :
#        process.treePath = cms.Path(                  process.metFilters                     + process.tree)

## Create output file
if options.dropAnalyzerDumpEDM == False:	
   ## Setup the service to make a ROOT TTree
   process.TFileService = cms.Service("TFileService", 
		fileName = cms.string(options.outputFileName))
else:
   # Make edm File storing all the products from the processName (current process)
   process.out = cms.OutputModule("PoolOutputModule",
                                      fileName = cms.untracked.string(options.outputFileName),
                                      outputCommands = cms.untracked.vstring(
                                        'drop *',
                                      	'keep *_*pat*_*_*'+options.processName+'*',
                                      	'keep *_*Pat*_*_*'+options.processName+'*',
                                      	'keep *_*T1*_*_*'+options.processName+'*',
                                      	'keep *_*metSysProducer*_*_*'+options.processName+'*',
                                      	'drop *_*T0*_*_*'+options.processName+'*',
                                      	'drop *_*T2*_*_*'+options.processName+'*',
                                      	'keep *_*slimmed*_*_*'+options.processName+'*',
                                      	'keep *_*slimmedMETs*_*_*',
					'keep *_*selectedObjects*_*_*',
                                      	),
                                    )

   process.output = cms.EndPath(process.out)

processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()