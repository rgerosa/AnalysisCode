import os, re
import FWCore.ParameterSet.Config as cms
  
### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

## data or MC options
options.register (
	'isMC',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate data or MC');

options.register (
	'isFastSIM',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate full or fast SIM for MC');

## filter Events options
options.register (
	'filterHighMETEvents',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate if apply or not MET filters');

options.register (
	'filterOnHLT',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate if apply or not trigger requirements');

## JEC options
options.register (
	'usePrivateSQlite',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'if a private SQL file with JEC to be found in test directory');

options.register (
	'applyL2L3Residuals',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'apply or not L2L3 Residual JEC on data');

options.register (
	'addPileupJetID',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	're-compute pileup-jet id for AK4 jets');

options.register (
	'addQGLikelihood',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'compute QGLikelihood for AK4 jets');

## PUPPI Jets MET, MVA Met
options.register (
	'addPuppiJets',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'add PUPPI jets to the output');

options.register (
	'addPuppiMET',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'add PUPPI met to the output');

options.register (
	'addMVAMet',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'compute MVAMet');
  	
## MET options
options.register (
        'useMiniAODMet',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'use the default MET in minoAOD without re-applying corrections');
	
options.register (
	'addMETSystematics',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'recompute Puppi MET propagating JEC from Jet + systematics');

## do substructure for CHS or Puppi jets
options.register (
	'doSubstructureCHS',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'run substructure algo for AK8CHS jets (Pruning, softDrop)');

options.register (
	'doSubstructurePuppi',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'run substructure algo for AK8Puppi jets (Pruning, softDrop)');

## processName
options.register (
	'processName','TREE',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'process name to be considered');

options.register (
	'miniAODProcess','PAT',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'process name used for miniAOD production (target is miniAODv2)');

## outputFile Name
options.register (
	'outputFileName','tree.root',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'output file name created by cmsRun');

## GT to be used    
options.register (
	'globalTag','74X_dataRun2_v5',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'gloabl tag to be uses');

## JEC    
options.register (
	'JECEra','Summer15_25nsV6',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'JEC correction era');

## Dump Gen Level info
options.register (
	'useLHEWeights',True,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
	'Dump LHE weights');

options.register(
	'addQCDPDFWeights',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
	'Dump weights related to QCD scale and PDF variations');

options.register(
	'isWorZorSignalMCSample',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
	' dump gen info for W, Z , Photon and DM particles');

## input cross section in case you want to store a different value wrt to the LHE file
options.register(
	'crossSection',-1.,VarParsing.multiplicity.singleton, VarParsing.varType.float,
	'external value for sample cross section, in case of data it is fixed to 0.001');
## Debug options
options.register (
	'dropAnalyzerDumpEDM',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
	'not run the analyzer and store an edm file');

options.register ('nThreads',4,VarParsing.multiplicity.singleton, VarParsing.varType.int,
		  'default number of threads');

## parsing command line arguments
options.parseArguments()

### check consistentcy of basic options
if options.isMC and 'dataRun2' in options.globalTag:
		options.globalTag = '74X_mcRun2_asymptotic_v2';

if options.isMC and options.applyL2L3Residuals:
	options.applyL2L3Residuals = False
	
if not options.isMC:
	options.crossSection = -1.;

## set by default a right environment when 76X is used
CMSSW_VERSION = os.environ['CMSSW_VERSION'];   
if re.match("CMSSW_7_6_.*",CMSSW_VERSION):
	if "74X_mcRun2" in options.globalTag:
		options.globalTag = '76X_mcRun2_asymptotic_v12'
	elif "74X_dataRun2" in options.globalTag:
		options.globalTag = '76X_dataRun2_v15'

	if options.usePrivateSQlite == True:
		options.usePrivateSQlite = False

print "##### Settings ######"
print "Running with isMC                = ",options.isMC	
print "Running with filterHighMETEvents = ",options.filterHighMETEvents	
print "Running with filterOnHLT         = ",options.filterOnHLT	
print "Running with usePrivateSQlite    = ",options.usePrivateSQlite	
print "Running with applyL2L3Residuals  = ",options.applyL2L3Residuals	
print "Running with addPuppiJets        = ",options.addPuppiJets
print "Running with addPuppiMET         = ",options.addPuppiMET
print "Running with useMiniAODMet       = ",options.useMiniAODMet
print "Running with addMETSystematics    = ",options.addMETSystematics	
print "Running with processName         = ",options.processName	
print "Running with miniAODProcess      = ",options.miniAODProcess	
print "Running with outputFileName      = ",options.outputFileName	
print "Running with globalTag           = ",options.globalTag	
print "Running with JEC Era             = ",options.JECEra	
print "Running with dropAnalyzerDumpEDM = ",options.dropAnalyzerDumpEDM	
print "Running with addPileupJetID      = ",options.addPileupJetID
print "Running with addQGLikelihood     = ",options.addQGLikelihood
print "Running with addMVAMet           = ",options.addMVAMet
print "Running with doSubstructureCHS   = ",options.doSubstructureCHS
print "Running with doSubstructurePuppi = ",options.doSubstructurePuppi
print "Running with useLHEWeights       = ",options.useLHEWeights
print "Running with addQCDPDFWeights    = ",options.addQCDPDFWeights
print "Running with isWorZorSignalMCSample  = ",options.isWorZorSignalMCSample
print "Running with nThreads            = ",options.nThreads
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
process.MessageLogger.cerr.FwkReport.reportEvery = 100

## Define the input source
if options.inputFiles == []:

	process.source = cms.Source("PoolSource", 
   		 fileNames = cms.untracked.vstring())

	if not options.isMC :
		process.source.fileNames.append(
			'/store/data/Run2015D/MET/MINIAOD/PromptReco-v4/000/259/810/00000/DC35E3E2-297B-E511-B0E5-02163E011E2B.root')
#        	'/store/data/Run2015D/MET/MINIAOD/PromptReco-v4/000/258/750/00000/5EE58B11-7572-E511-B952-02163E014378.root')
	else:
		process.source.fileNames.append( 
#			'/store/mc/RunIISpring15MiniAODv2/VBF_HToInvisible_M110_13TeV_powheg_pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/084D5D5F-F36D-E511-8D98-02163E014DC0.root'
#			'root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/AxialMonoW_Mphi-1000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/D4A17AAE-9174-E511-B19C-002590FD5A78.root'
			'/store/mc/RunIISpring15MiniAODv2/DM_ScalarWH_Mphi-100_Mchi-50_gSM-1p0_gDM-1p0_13TeV-JHUGen/MINIAODSIM/Asympt25ns_74X_mcRun2_asymptotic_v2-v1/40000/9048CCCB-ADA4-E511-87C8-0025B3E015D2.root'
#			'/store/mc/RunIISpring15MiniAODv2/DMS_NNPDF30_Pseudoscalar_Mphi-50_Mchi-10_gSM-1p0_gDM-1p0_13TeV-powheg/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/50000/70D1C4EA-7174-E511-9B2D-441EA1615F8A.root'
#			'root://xrootd.unl.edu//store/user/emanuele/POWHEG_DMV_AV_NNPDF30_13TeV/POWHEG_DMV_AV_NNPDF30_13TeV_MINIAOD_V1/160222_152908/0000/miniaodsim_1.root'
			#'root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/ZJetsToNuNu_HT-100To200_13TeV-madgraph/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/40000/008902DD-9F6F-E511-BCE9-0025904C540C.root'
#			'root://xrootd.unl.edu//store/mc/RunIISpring15MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/10000/12608B5D-E66D-E511-B233-441EA173397A.root'			
#			'/store/mc/RunIISpring15MiniAODv2/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/74X_mcRun2_asymptotic_v2-v1/60000/0002EB7B-E76D-E511-8AD6-00269E95B17C.root'
			)    	
else:
   process.source = cms.Source("PoolSource",
   	  fileNames = cms.untracked.vstring(options.inputFiles))


## Set the process options -- Display summary at the end, enable unscheduled execution
if options.nThreads == 1 or options.nThreads == 0:
	process.options = cms.untracked.PSet( 
		allowUnscheduled = cms.untracked.bool(True),
		wantSummary = cms.untracked.bool(True))
else:
	process.options = cms.untracked.PSet( 
		allowUnscheduled = cms.untracked.bool(True),
		wantSummary = cms.untracked.bool(True),
		numberOfThreads = cms.untracked.uint32(options.nThreads),
		numberOfStreams = cms.untracked.uint32(options.nThreads))

## How many events to process
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(options.maxEvents))

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
JetMetCorrector(process,"slimmedJets","slimmedMETs","AK4PFchs",options.isMC, options.applyL2L3Residuals, options.useMiniAODMet, options.addMETSystematics)
## apply JEC and propagation on MET for AK4PFPuppi
if options.addPuppiJets or options.addPuppiMET:
	JetMetCorrector(process,"slimmedJetsPuppi","slimmedMETsPuppi","AK4PFPuppi",options.isMC, options.applyL2L3Residuals, options.useMiniAODMet, options.addMETSystematics)

# Create a set of objects to read from
process.selectedObjects = cms.EDProducer("PFCleaner",
     vertices  = cms.InputTag("goodVertices"),
     pfcands   = cms.InputTag("packedPFCandidates"),
     muons     = cms.InputTag("slimmedMuons"),
     electrons = cms.InputTag("slimmedElectrons"),
     photons   = cms.InputTag("slimmedPhotons"),
     rho       = cms.InputTag("fixedGridRhoFastjetAll"),
     jets      = cms.InputTag("slimmedJetsRecorrectedAK4PFchs"),
     electronidveto  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
     electronidloose = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
     electronidtight = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
     electronidheep  = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
     photonidloose  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
     photonidmedium = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
     photonidtight  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
     photonsieie = cms.InputTag("photonIDValueMapProducer", "phoFull5x5SigmaIEtaIEta"),
     photonphiso = cms.InputTag("photonIDValueMapProducer", "phoPhotonIsolation"),
     photonchiso = cms.InputTag("photonIDValueMapProducer", "phoChargedIsolation")
)

## modify some existing jet collections adding pileup-jet id and QGLikelihood from GT
from AnalysisCode.MonoXAnalysis.JetTools_cff import addPileupJetID, addQGLikelihood

if options.addPileupJetID:  ## so far not working for slimmedJets, so re-clustering and update
	addPileupJetID(process, collection = "slimmedJetsRecorrectedAK4PFchs", postfix = "")
	if options.addPuppiJets:
		addPileupJetID(process, collection = "slimmedJetsRecorrectedAK4PFPuppi", postfix = "Puppi")
	
if options.addQGLikelihood:
	inputJet = "slimmedJetsRecorrectedAK4PFchs";
	inputJetPuppi = "slimmedJetsRecorrectedAK4PFPuppi";
	if options.addPileupJetID:
		inputJet += "PUID";
		inputJetPuppi += "PUID";

	addQGLikelihood( process,collection = inputJet, postfix = "");
	if options.addPuppiJets:
		addQGLikelihood( process,collection = inputJetPuppi, postfix = "Puppi");


## in case run the MVA met producer
from AnalysisCode.MonoXAnalysis.MVAMet_cff import runMVAMet

if options.addMVAMet:
	## to parse leptons we need a list of CandidateView not a value map with Refs
	#leptons = ["PFCleaner:tightmuons","PFCleaner:tightelectrons"]
	leptons = [];
	runMVAMet(process,isMC = options.isMC,leptons = leptons )


# Define all the METs corrected for lepton/photon momenta
process.t1mumet = cms.EDProducer("MuonCorrectedMETProducer",
   met   = cms.InputTag("slimmedMETs","","TREE"),
   cands = cms.VInputTag(cms.InputTag("selectedObjects", "muons")),
   isPuppi = cms.bool(False),
   pfCandidates = cms.InputTag("packedPFCandidates"))

if options.useMiniAODMet:
	process.t1mumet.met = cms.InputTag("slimmedMETs","","PAT")

process.t1elmet = cms.EDProducer("ElectronCorrectedMETProducer",
   met = cms.InputTag("slimmedMETs","","TREE"),
   cands = cms.VInputTag(cms.InputTag("selectedObjects", "electrons")),
   isPuppi = cms.bool(False),
   pfCandidates = cms.InputTag("packedPFCandidates"))

if options.useMiniAODMet:
	process.t1elmet.met = cms.InputTag("slimmedMETs","","PAT")

process.t1phmet = cms.EDProducer("PhotonCorrectedMETProducer",
   met = cms.InputTag("slimmedMETs","","TREE"),
   cands = cms.VInputTag(cms.InputTag("selectedObjects", "photons")),
   isPuppi = cms.bool(False),
   pfCandidates = cms.InputTag("packedPFCandidates"))

if options.useMiniAODMet:
	process.t1phmet.met = cms.InputTag("slimmedMETs","","PAT")

if options.addPuppiMET:

	process.puppit1mumet = process.t1mumet.clone(
		met   = cms.InputTag("slimmedMETsPuppi","","TREE"),
		isPuppi = cms.bool(True))

	if options.useMiniAODMet:
		process.puppit1mumet.met = cms.InputTag("slimmedMETsPuppi","","PAT")

	process.puppit1elmet = process.t1elmet.clone(
		met = cms.InputTag("slimmedMETsPuppi","","TREE"),
		isPuppi = cms.bool(True))

	if options.useMiniAODMet:
		process.puppit1elmet.met = cms.InputTag("slimmedMETsPuppi","","PAT")

	process.puppit1phmet = process.t1phmet.clone(
		met = cms.InputTag("slimmedMETsPuppi","","TREE"),
		isPuppi = cms.bool(True))

	if options.useMiniAODMet:
		process.puppit1phmet.met = cms.InputTag("slimmedMETsPuppi","","PAT")
	

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
                                      	'keep *_*T1*_*_*'+options.processName+'*',
                                      	'keep *_*metSysProducer*_*_*'+options.processName+'*',
                                      	'drop *_*T0*_*_*'+options.processName+'*',
                                      	'drop *_*T2*_*_*'+options.processName+'*',
                                      	'keep *_*slimmed*_*_*'+options.processName+'*',
					'keep *_*slimmedJets*_*_*',
                                      	'keep *_*slimmedMETs*_*_*',
                                      	'keep *_patJetsAK8*_*_*',
                                      	'keep *_*Matched_*_*',
                                      	'keep *_*Packed_*_*',
					'keep *_*selectedObjects*_*_*',
					'keep *_*mvaMET*_*_*',
					'keep *_*t1mumet*_*_*',
					'keep *_*t1elmet*_*_*',
					'keep *_*t1phmet*_*_*',
					'keep *_*t1taumet*_*_*',
                                      	))

   process.output = cms.EndPath(process.out)

#### substructure sequence
from AnalysisCode.MonoXAnalysis.JetSubstructure_cff import JetSubstructure
if options.doSubstructureCHS:
	JetSubstructure(process,
			options.isMC,
			coneSize = 0.8, algo = "AK",
			pileupMethod = "chs", selection = "pt > 150 && abs(eta) < 2.5",
			addPruning = True, addSoftDrop = True, addTrimming = False, addFiltering = False,
			addNsubjettiness = True, addEnergyCorrelation = True, addQJets = False,
			addQGLikelihood = True);

if options.doSubstructurePuppi:

	JetSubstructure(process,
			options.isMC,
			coneSize = 0.8, algo = "AK",
			pileupMethod = "Puppi", selection = "pt > 150 && abs(eta) < 2.5",
			addPruning = True, addSoftDrop = True, addTrimming = False, addFiltering = False,
			addNsubjettiness = True, addEnergyCorrelation = True, addQJets = False,
			addQGLikelihood = True);

# MET filter --> making the or of different met
process.filterHighRecoil = cms.EDFilter("PATMETFilter",					
    metCollections = cms.VPSet(
		cms.PSet( srcMet = cms.InputTag("slimmedMETs","","TREE"),
			  metCut = cms.double(100)),
		cms.PSet(srcMet = cms.InputTag("t1mumet","","TREE"),
			 metCut = cms.double(100)),
		cms.PSet(srcMet = cms.InputTag("t1elmet","","TREE"),
			 metCut = cms.double(100)),
		cms.PSet(srcMet = cms.InputTag("t1phmet","","TREE"),
			 metCut = cms.double(100)),
		),
    filterEvents = cms.bool(options.filterHighMETEvents),
    graterThan = cms.bool(True),
    applyAndInsteadOfOr = cms.bool(False)
)

if options.useMiniAODMet:
	process.filterHighRecoil.metCollections[0].srcMet = cms.InputTag("slimmedMETs","","PAT")
     					
# Tree for the generator weights
process.gentree = cms.EDAnalyzer("LHEWeightsTreeMaker",
    lheinfo = cms.InputTag("externalLHEProducer"),
    lheRuninfo = cms.InputTag("externalLHEProducer"),
    geninfo = cms.InputTag("generator"),
    pileupinfo = cms.InputTag("slimmedAddPileupInfo"),				 
    uselheweights = cms.bool(options.useLHEWeights),
    addqcdpdfweights = cms.bool(options.addQCDPDFWeights))

# Make the tree 
process.tree = cms.EDAnalyzer("MonoJetTreeMaker",
   ## gen info			     
   isMC    = cms.bool(options.isMC),
   uselheweights  = cms.bool(options.useLHEWeights),
   isWorZorSignalMCSample = cms.bool(options.isWorZorSignalMCSample),
   lheinfo = cms.InputTag("externalLHEProducer"),			      
   lheRuninfo = cms.InputTag("externalLHEProducer"),			      
   pileup  = cms.InputTag("slimmedAddPileupInfo"),
   genevt  = cms.InputTag("generator"),
   gens    = cms.InputTag("prunedGenParticles"),
   xsec    = cms.double(options.crossSection),   
   ## trigger info
   triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
   prescales = cms.InputTag("patTrigger"),    
   filterResults  = cms.InputTag("TriggerResults", "", options.miniAODProcess),
   hbheloose = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Loose"),
   hbhetight = cms.InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResultRun2Tight"),
   hbheiso   = cms.InputTag("HBHENoiseFilterResultProducer","HBHEIsoNoiseFilterResult"),
   ## vertexes			    
   vertices = cms.InputTag("goodVertices"),
   ## muons    
   muons       = cms.InputTag("selectedObjects","muons"),
   tightmuons  = cms.InputTag("selectedObjects","tightmuons"),
   highptmuons = cms.InputTag("selectedObjects","highptmuons"),
   ## electrons
   electrons       = cms.InputTag("selectedObjects", "electrons"),
   tightelectrons  = cms.InputTag("selectedObjects", "tightelectrons"),
   heepelectrons   = cms.InputTag("selectedObjects", "heepelectrons"),
   electronLooseId = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
   ## photons
   photons        = cms.InputTag("selectedObjects", "photons"),
   tightphotons   = cms.InputTag("selectedObjects", "tightphotons"),
   photonHighPtId = cms.InputTag("selectedObjects", "photonHighPtId"),
   photonLooseId  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose"),
   photonMediumId = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium"),
   photonTightId  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"),
   ## taus
   taus = cms.InputTag("slimmedTaus"),
   ## jets AK4
   jets         = cms.InputTag("slimmedJetsRecorrectedAK4PFchs"),
   addPuppiJets = cms.bool(options.addPuppiJets),			      
   puppijets    = cms.InputTag("slimmedJetsRecorrectedAK4PFPuppi"),			      
   ## MET
   t1met = cms.InputTag("slimmedMETs","","TREE"),
   t1mumet = cms.InputTag("t1mumet"),
   t1elmet = cms.InputTag("t1elmet"),
   t1phmet = cms.InputTag("t1phmet"),
   ## Puppi MET
   addPuppiMET = cms.bool(options.addPuppiMET),
   puppit1met = cms.InputTag("slimmedMETsPuppi","","TREE"),
   puppit1mumet = cms.InputTag("puppit1mumet"),
   puppit1elmet = cms.InputTag("puppit1elmet"),
   puppit1phmet = cms.InputTag("puppit1phmet"),
   ## MET systematics
   addMETSystematics = cms.bool(options.addMETSystematics),    			      
   ## mvaMet
   addMVAMet = cms.bool(options.addMVAMet),			     
   mvaMET = cms.InputTag("mvaMET"),			      
   ## trigger filter
   applyHLTFilter = cms.bool(options.filterOnHLT),
   ## clean objects
   cleanMuonJet = cms.bool(True),
   cleanElectronJet = cms.bool(True),
   cleanPhotonJet = cms.bool(True),
   ## CHS jet substructure
   addSubstructureCHS = cms.bool(options.doSubstructureCHS),
   boostedJetsCHS     = cms.InputTag("packedPatJetsAK8PFJetsCHS"),
   boostedJetsOriginal     = cms.InputTag("slimmedJetsAK8"),
   addSubstructurePuppi = cms.bool(options.doSubstructurePuppi),
   boostedJetsPuppi     = cms.InputTag("packedPatJetsAK8PFJetsPuppi"),
   ## b-tag scale factors
   addBTagScaleFactor  = cms.bool(True),
   bTagScaleFactorFile = cms.FileInPath('AnalysisCode/MonoXAnalysis/data/BTagScaleFactors/pfCombinedInclusiveSecondaryVertexV2BJetTags_74X.csv') 			      
)

if options.useMiniAODMet:
	process.tree.t1met = cms.InputTag("slimmedMETs","","PAT")
	process.tree.puppit1met = cms.InputTag("slimmedMETsPuppi","","PAT")

jetColl      = "slimmedJetsRecorrectedAK4PFchs"
jetCollPuppi = "slimmedJetsRecorrectedAK4PFPuppi"

if options.addPileupJetID == True and options.addQGLikelihood == False:
	process.tree.jets = cms.InputTag("slimmedJetsRecorrectedAK4PFchsPUID");
	process.tree.puppijets = cms.InputTag("slimmedJetsRecorrectedAK4PFPuppiPUID");
	jetColl = "slimmedJetsRecorrectedAK4PFchsPUID"
	jetCollPuppi = "slimmedJetsRecorrectedAK4PFPuppiPUID"
	
if options.addQGLikelihood == True and options.addPileupJetID == False:
	process.tree.jets = cms.InputTag("slimmedJetsRecorrectedAK4PFchsQG");
	process.tree.puppijets = cms.InputTag("slimmedJetsRecorrectedAK4PFPuppiQG");
	jetColl = "slimmedJetsRecorrectedAK4PFchsQG"
	jetCollPuppi = "slimmedJetsRecorrectedAK4PFPuppiQG"
else:
	process.tree.jets = cms.InputTag("slimmedJetsRecorrectedAK4PFchsPUIDQG");
	process.tree.puppijets = cms.InputTag("slimmedJetsRecorrectedAK4PFPuppiPUIDQG");
	jetColl = "slimmedJetsRecorrectedAK4PFchsPUIDQG"
	jetCollPuppi = "slimmedJetsRecorrectedAK4PFPuppiPUIDQG"


# Histo for Btag efficiency
process.btageff = cms.EDAnalyzer("BTaggingEfficiencyTreeMaker",
				 directoryName = cms.string("btagEff"),
				 srcJets = cms.InputTag(jetColl),
				 dRClean = cms.double(0.4),
				 cleanMuonJet = cms.bool(True),
				 srcMuons = cms.InputTag("selectedObjects","muons"),
				 cleanElectronJet = cms.bool(True),
				 srcElectrons = cms.InputTag("selectedObjects","electrons"),
				 cleanPhotonJet = cms.bool(True),
				 srcPhotons = cms.InputTag("selectedObjects","photons"),
				 selection = cms.string('abs(eta)<2.4 && pt > 20'),
				 ptBins  = cms.vdouble(20,30,40,50,60,80,110,150,1000),
				 etaBins = cms.vdouble(0.,0.5,1.,1.5,2.0,2.4),
				 bDiscriminatorInfo = cms.VPSet(
		cms.PSet(discriminatorName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
			 wpLabel = cms.string("Loose"),
			 wpValue = cms.double(0.605)),

		cms.PSet(discriminatorName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
			 wpLabel = cms.string("Medium"),
			 wpValue = cms.double(0.89)),

		cms.PSet(discriminatorName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
			 wpLabel = cms.string("Tight"),
			 wpValue = cms.double(0.97))
		))

## some options for 76X				 
if re.match("CMSSW_7_6_.*",CMSSW_VERSION):
	process.btageff.bDiscriminatorInfo[0].wpValue = cms.double(0.460)
	process.btageff.bDiscriminatorInfo[1].wpValue = cms.double(0.800)
	process.btageff.bDiscriminatorInfo[2].wpValue = cms.double(0.935)
	
if options.addPuppiJets:
	setattr(process,"btageffPuppi",process.btageff.clone(
			srcJets = cms.InputTag(jetCollPuppi)))


## fast sim business
if options.isFastSIM:

	## fix LHE info
	process.tree.lheinfo = cms.InputTag("source")
	process.tree.lheRuninfo = cms.InputTag("source")

	process.gentree.lheinfo = cms.InputTag("source")
	process.gentree.lheRuninfo = cms.InputTag("source")
	
	## remove HCAL noise
	process.metFilters.remove(process.HBHENoiseFilterResultProducer)

	## remove HCAL noise info
	process.tree.hbheloose = cms.InputTag("");
	process.tree.hbhetight = cms.InputTag("");
	process.tree.hbheiso   = cms.InputTag("");

	process.tree.pileup = cms.InputTag("addPileupInfo");
	## tem fix
	process.slimmedMETs.runningOnMiniAOD = cms.bool(False);
	process.slimmedMETsPuppi.runningOnMiniAOD = cms.bool(False);
	
# Set up the path
if options.dropAnalyzerDumpEDM == False:
	if (options.isMC):
		process.treePath = cms.Path(process.gentree + 
					    process.metFilters + 
					    process.btageff+
					    process.filterHighRecoil + 
					    process.tree)
		if(options.addPuppiJets):
			process.treePath = cms.Path(process.gentree +
						    process.metFilters +
						    process.btageff+
						    process.btageffPuppi+
						    process.filterHighRecoil +
						    process.tree)
			
	else :
		process.treePath = cms.Path(process.metFilters + 
					    process.filterHighRecoil + 
					    process.tree)

processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
