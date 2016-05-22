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
	'metCut',125.,VarParsing.multiplicity.singleton,VarParsing.varType.float,
	'met/recoil cut to be applied if filterHighMETEvents is set to true');

options.register (
	'filterOnHLT',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate if apply or not trigger requirements');

## JEC options
options.register (
	'usePrivateSQliteJEC',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'if a private SQL file with JEC to be found in test directory');

options.register (
	'usePrivateSQliteJER',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'if a private SQL file with JER to be found in test directory');

options.register (
	'applyL2L3Residuals',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'apply or not L2L3 Residual JEC on data');

options.register (
	'addPileupJetID',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
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

#### Add scale and smear corrections for electrons and photons
options.register (
	'addEGMSmear',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'add e-gamma scale and resolution corrections for electrons and photons');

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

options.register (
	'useOfficialMETSystematics',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'run the official tool for met uncertainty --> does a lot of things but slow .. otherwise minimal home made validated code');

## do substructure for CHS or Puppi jets
options.register (
	'addSubstructureCHS',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'run substructure algo for AK8CHS jets (Pruning, softDrop)');

options.register (
	'addSubstructurePuppi',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'run substructure algo for AK8Puppi jets (Pruning, softDrop)');

## processName
options.register (
	'processName','TREE',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'process name to be considered');

options.register (
	'miniAODProcess','RECO',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'process name used for miniAOD production (target is miniAODv2)');

## outputFile Name
options.register (
	'outputFileName','tree.root',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'output file name created by cmsRun');

## GT to be used    
options.register (
	'globalTag','76X_dataRun2_16Dec2015_v0',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'gloabl tag to be uses');

## JEC    
options.register (
	'JECEra','Fall15_25nsV2',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'JEC correction era');

## Dump Gen Level info
options.register (
	'useLHEWeights',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
	'Dump LHE weights');

options.register(
	'addQCDPDFWeights',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
	'Dump weights related to QCD scale and PDF variations');

options.register(
	'addGenParticles',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
	'dump gen info for W,Z,top,photon');

options.register(
	'isSignalSample',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
	'dump DM particles and read LHE header');

## input cross section in case you want to store a different value wrt to the LHE file
options.register(
	'crossSection',-1.,VarParsing.multiplicity.singleton, VarParsing.varType.float,
	'external value for sample cross section, in case of data it is fixed to 0.001');
## Debug options
options.register (
	'dropAnalyzerDumpEDM',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
	'not run the analyzer and store an edm file');

options.register (
	'nThreads',4,VarParsing.multiplicity.singleton, VarParsing.varType.int,
	'default number of threads');

## to be used when running crab jobs with local files                                                                                                                           
options.register ('isCrab',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
		  'to be used to handle local files with crab');

## parsing command line arguments
options.parseArguments()

### check consistentcy of basic options
if options.isMC and 'dataRun2' in options.globalTag:
	options.globalTag = '76X_mcRun2_asymptotic_RunIIFall15DR76_v1';

if options.isMC and options.applyL2L3Residuals:
	options.applyL2L3Residuals = False
	
if not options.isMC:
	options.crossSection = -1.;

if options.isMC and options.miniAODProcess != 'PAT':
	options.miniAODProcess  = 'PAT'

print "##### Settings ######"
print "Running with isMC                = ",options.isMC	
print "Running with filterHighMETEvents = ",options.filterHighMETEvents	
if options.filterHighMETEvents:
	print "Running with metCut              = ",options.metCut
print "Running with filterOnHLT         = ",options.filterOnHLT	
print "Running with usePrivateSQliteJEC = ",options.usePrivateSQliteJEC	
print "Running with usePrivateSQliteJER = ",options.usePrivateSQliteJER	
print "Running with applyL2L3Residuals  = ",options.applyL2L3Residuals	
print "Running with addPuppiJets        = ",options.addPuppiJets
print "Running with addPuppiMET         = ",options.addPuppiMET
print "Running with addEGMSmear         = ",options.addEGMSmear
print "Running with useMiniAODMet       = ",options.useMiniAODMet
print "Running with addMETSystematics   = ",options.addMETSystematics	
print "Running with processName         = ",options.processName	
print "Running with miniAODProcess      = ",options.miniAODProcess	
print "Running with outputFileName      = ",options.outputFileName	
print "Running with globalTag           = ",options.globalTag	
print "Running with JEC Era             = ",options.JECEra	
print "Running with dropAnalyzerDumpEDM = ",options.dropAnalyzerDumpEDM	
print "Running with addPileupJetID      = ",options.addPileupJetID
print "Running with addQGLikelihood     = ",options.addQGLikelihood
print "Running with addMVAMet           = ",options.addMVAMet
print "Running with addSubstructureCHS  = ",options.addSubstructureCHS
print "Running with addSubstructurePuppi= ",options.addSubstructurePuppi
print "Running with useLHEWeights       = ",options.useLHEWeights
print "Running with addQCDPDFWeights    = ",options.addQCDPDFWeights
print "Running with isSignalSample      = ",options.isSignalSample
print "Running with addGenParticles     = ",options.addGenParticles
print "Running with nThreads            = ",options.nThreads
print "Running with isCrab              = ",options.isCrab
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
			'/store/data/Run2015D/SingleElectron/MINIAOD/16Dec2015-v1/20000/00050EF1-F9A6-E511-86B2-0025905A48D0.root')
	else:
		process.source.fileNames.append(
#			'file:pickevents.root'
#			'/store/mc/RunIIFall15MiniAODv2/ZJetsToNuNu_HT-100To200_13TeV-madgraph/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/060FC9A4-C8BD-E511-B138-000F530E46D0.root',
			'root://xrootd.unl.edu//store/mc/RunIIFall15MiniAODv2/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/70000/00761843-D4BD-E511-853E-000F53273498.root'		       
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
JECConfiguration(process,options.usePrivateSQliteJEC,options.JECEra,options.isMC,options.applyL2L3Residuals,options.isCrab)

from AnalysisCode.MonoXAnalysis.JERConfiguration_cff import JERConfiguration
## connect to a local SQLite file or take corrections from GT
JERConfiguration(process,options.usePrivateSQliteJER,options.JECEra,options.isMC,options.isCrab)

## Setup MET filters or not --> 76X everything inside miniAOD is already good
process.load('AnalysisCode.MonoXAnalysis.METFilters_cff')

# run cut-based electron ID https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
from AnalysisCode.MonoXAnalysis.ElectronTools_cff import ElectronTools
ElectronTools(process,options.addEGMSmear,options.isMC)

# run cut-based photon ID 
from AnalysisCode.MonoXAnalysis.PhotonTools_cff import PhotonTools
PhotonTools(process,options.addEGMSmear,options.isMC)

# Apply JEC on jets and update them
from AnalysisCode.MonoXAnalysis.JetTools_cff import JetCorrector
## apply JEC and propagation on MET for AK4PFchs
jetCollName      = "slimmedJets"
jetPuppiCollName = "slimmedJetsPuppi"

jetCollName = JetCorrector(process,jetCollName,"AK4PFchs",options.isMC, options.applyL2L3Residuals)
## apply JEC and propagation on MET for AK4PFPuppi
if options.addPuppiJets:
	jetPuppiCollName = JetCorrector(process,jetPuppiCollName,"AK4PFPuppi",options.isMC,options.applyL2L3Residuals)
	
# to apply analysis selections
process.load('AnalysisCode.MonoXAnalysis.selectionObjects_cfi')
process.selectedObjects.jets = cms.InputTag(jetCollName)
process.selectedObjects.useCalibratedElectrons = cms.bool(options.addEGMSmear)
process.selectedObjects.useCalibratedPhotons = cms.bool(options.addEGMSmear)

## modify some existing jet collections adding pileup-jet id and QGLikelihood from GT
from AnalysisCode.MonoXAnalysis.JetTools_cff import addPileupJetID, addQGLikelihood

if options.addPileupJetID:  ## so far not working for slimmedJets, so re-clustering and update	
	jetCollName = addPileupJetID(process,jetCollName,"",options.isMC)

if options.addQGLikelihood:
	jetCollName = addQGLikelihood(process,jetCollName,"");
	if options.addPuppiJets:
		jetPuppiCollName = addQGLikelihood(process,jetPuppiCollName,"Puppi");

### correct the MET
from AnalysisCode.MonoXAnalysis.metCorrector_cff import metCorrector
if not options.useMiniAODMet:
	metCollection = "slimmedMETs"
	metCorrector(process,jetCollName,metCollection,options.isMC,"AK4PFchs",options.applyL2L3Residuals,options.addMETSystematics,options.useOfficialMETSystematics);
	
	if options.addPuppiMET:
		metCollectionPuppi = "slimmedMETsPuppi"
		metCorrector(process,jetPuppiCollName,metCollectionPuppi,options.isMC,"AK4PFPuppi",options.applyL2L3Residuals,options.addMETSystematics,options.useOfficialMETSystematics);
		

## in case run the MVA met producer
from AnalysisCode.MonoXAnalysis.MVAMet_cff import runMVAMet

if options.addMVAMet:
	## to parse leptons we need a list of CandidateView not a value map with Refs
	#leptons = ["PFCleaner:tightmuons","PFCleaner:tightelectrons"]
	leptons = [];
	runMVAMet(process,isMC = options.isMC,leptons = leptons )


# Define all the METs corrected for lepton/photon momenta
from AnalysisCode.MonoXAnalysis.recoil_cff import recoilComputation
recoilComputation(process,options.processName,options.miniAODProcess,options.useMiniAODMet,False)
if options.addPuppiMET:
	recoilComputation(process,options.processName,options.miniAODProcess,options.useMiniAODMet,True)


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
boostedJetCollection = "";
boostedPuppiJetCollection = "";
if options.addSubstructureCHS:
	boostedJetCollection = JetSubstructure(process,
					       options.isMC,
					       coneSize = 0.8, algo = "AK",
					       pileupMethod = "chs", selection = "pt > 150 && abs(eta) < 2.5",
					       addPruning = True, addSoftDrop = True, addTrimming = False, addFiltering = False,
					       addNsubjettiness = True, addEnergyCorrelation = True, addQJets = False,
					       addQGLikelihood = True);

if options.addSubstructurePuppi:
	boostedPuppiJetCollection = JetSubstructure(process,
						    options.isMC,
						    coneSize = 0.8, algo = "AK",
						    pileupMethod = "Puppi", selection = "pt > 150 && abs(eta) < 2.5",
						    addPruning = True, addSoftDrop = True, addTrimming = False, addFiltering = False,
						    addNsubjettiness = True, addEnergyCorrelation = True, addQJets = False,
						    addQGLikelihood = True);

# MET filter --> making the or of different met
process.filterHighRecoil = cms.EDFilter("PATMETFilter",					
    metCollections = cms.VPSet(
		cms.PSet( srcMet = cms.InputTag("slimmedMETs","",options.processName),
			  metCut = cms.double(options.metCut)),
		cms.PSet(srcMet = cms.InputTag("t1mumet","",options.processName),
			 metCut = cms.double(options.metCut)),
		cms.PSet(srcMet = cms.InputTag("t1elmet","",options.processName),
			 metCut = cms.double(options.metCut)),
		cms.PSet(srcMet = cms.InputTag("t1phmet","",options.processName),
			 metCut = cms.double(options.metCut)),
		),
    filterEvents = cms.bool(options.filterHighMETEvents),
    graterThan = cms.bool(True),
    applyAndInsteadOfOr = cms.bool(False)
)

if options.useMiniAODMet:
	process.filterHighRecoil.metCollections[0].srcMet = cms.InputTag("slimmedMETs","",options.miniAODProcess)
     					
# Tree for the generator weights
process.gentree = cms.EDAnalyzer("LHEWeightsTreeMaker",
    lheinfo    = cms.InputTag("externalLHEProducer"),
    lheRuninfo = cms.InputTag("externalLHEProducer"),
    geninfo    = cms.InputTag("generator"),
    pileupinfo = cms.InputTag("slimmedAddPileupInfo"),				 
    genParticles  = cms.InputTag("prunedGenParticles"),				 
    uselheweights = cms.bool(options.useLHEWeights),
    addqcdpdfweights = cms.bool(options.addQCDPDFWeights),
    isSignalSample = cms.bool(options.isSignalSample))

# Make the tree 
process.tree = cms.EDAnalyzer("MonoJetTreeMaker",
			      ## gen info			     
			      isMC                   = cms.bool(options.isMC),
			      uselheweights          = cms.bool(options.useLHEWeights),
			      isSignalSample         = cms.bool(options.isSignalSample),			      
			      addGenParticles        = cms.bool(options.addGenParticles),			      
			      lheinfo                = cms.InputTag("externalLHEProducer"),			      
			      lheRuninfo             = cms.InputTag("externalLHEProducer"),			      
			      pileup                 = cms.InputTag("slimmedAddPileupInfo"),
			      genevt                 = cms.InputTag("generator"),
			      gens                   = cms.InputTag("prunedGenParticles"),
			      xsec                   = cms.double(options.crossSection),   
			      ## trigger info
			      triggerResults = cms.InputTag("TriggerResults", "", "HLT"),
			      prescales      = cms.InputTag("patTrigger"),    
			      filterResults  = cms.InputTag("TriggerResults", "", options.miniAODProcess),
			      ## vertexes			    
			      vertices       = cms.InputTag("goodVertices"),
			      ## muons    
			      muons          = cms.InputTag("selectedObjects","muons"),
			      tightmuons     = cms.InputTag("selectedObjects","tightmuons"),
			      highptmuons    = cms.InputTag("selectedObjects","highptmuons"),
			      ## electrons
			      electrons       = cms.InputTag("selectedObjects", "electrons"),
			      looseelectrons  = cms.InputTag("selectedObjects", "looseelectrons"),
			      tightelectrons  = cms.InputTag("selectedObjects", "tightelectrons"),
			      heepelectrons   = cms.InputTag("selectedObjects", "heepelectrons"),
			      ## photons --> can be matched by reference 
			      photons        = cms.InputTag("selectedObjects", "photons"),
			      mediumphotons   = cms.InputTag("selectedObjects", "mediumphotons"),
			      tightphotons   = cms.InputTag("selectedObjects", "tightphotons"),			     
			      photonHighPtId = cms.InputTag("selectedObjects", "photonHighPtId"),
			      photonLooseId = cms.InputTag("selectedObjects", "photonLooseId"),			     
			      ## taus
			      taus           = cms.InputTag("selectedObjects","taus"),
			      ## jets AK4
			      jets         = cms.InputTag(jetCollName),
			      addPuppiJets = cms.bool(options.addPuppiJets),			      
			      puppijets    = cms.InputTag(jetPuppiCollName),			      
			      ## MET
			      t1met     = cms.InputTag("slimmedMETs","",options.processName),
			      t1mumet   = cms.InputTag("t1mumet"),
			      t1elmet   = cms.InputTag("t1elmet"),
			      t1phmet   = cms.InputTag("t1phmet"),
			      t1taumet  = cms.InputTag("t1taumet"),
			      ## Puppi MET
			      addPuppiMET   = cms.bool(options.addPuppiMET),
			      puppit1met    = cms.InputTag("slimmedMETsPuppi","",options.processName),
			      puppit1mumet  = cms.InputTag("puppit1mumet"),
			      puppit1elmet  = cms.InputTag("puppit1elmet"),
			      puppit1phmet  = cms.InputTag("puppit1phmet"),
			      puppit1taumet = cms.InputTag("puppit1taumet"),
			      ## MET systematics
			      addMETSystematics = cms.bool(options.addMETSystematics),    			      
			      ## mvaMet
			      addMVAMet = cms.bool(options.addMVAMet),			     
			      mvaMET    = cms.InputTag("mvaMET"),			      
			      ## trigger filter
			      applyHLTFilter = cms.bool(options.filterOnHLT),
			      ## clean objects
			      cleanMuonJet     = cms.bool(True),
			      cleanElectronJet = cms.bool(True),
			      cleanPhotonJet   = cms.bool(True),
			      ## CHS jet substructure
			      addSubstructureCHS   = cms.bool(options.addSubstructureCHS),
			      boostedJetsCHS       = cms.InputTag(boostedJetCollection),
			      addSubstructurePuppi = cms.bool(options.addSubstructurePuppi),
			      boostedJetsPuppi     = cms.InputTag(boostedPuppiJetCollection),
			      ## b-tag scale factors
			      addBTagScaleFactor         = cms.bool(True),
			      bTagScaleFactorFileCSV     = cms.FileInPath('AnalysisCode/MonoXAnalysis/data/BTagScaleFactors/pfCombinedInclusiveSecondaryVertexV2BJetTags_76X.csv'), 	     
			      bTagScaleFactorFileMVA     = cms.FileInPath('AnalysisCode/MonoXAnalysis/data/BTagScaleFactors/pfCombinedMVAV2BJetTags_76X.csv'), 			      
			      bTagScaleFactorFileSubCSV  = cms.FileInPath('AnalysisCode/MonoXAnalysis/data/BTagScaleFactors/pfCombinedInclusiveSecondaryVertexV2BJetTags_76X_subjet.csv') 	  
			      )

if options.useMiniAODMet:
	process.tree.t1met = cms.InputTag("slimmedMETs","",options.miniAODProcess)
	process.tree.puppit1met = cms.InputTag("slimmedMETsPuppi","",options.miniAODProcess)


# Histo for Btag efficiency
process.btageff = cms.EDAnalyzer("BTaggingEfficiencyTreeMaker",
				 directoryName = cms.string("btagEff"),
				 srcJets       = cms.InputTag(jetCollName),
				 dRClean       = cms.double(0.4),
				 cleanMuonJet  = cms.bool(True),
				 srcMuons      = cms.InputTag("selectedObjects","muons"),
				 cleanElectronJet = cms.bool(True),
				 srcElectrons    = cms.InputTag("selectedObjects","electrons"),
				 cleanPhotonJet  = cms.bool(True),
				 srcPhotons      = cms.InputTag("selectedObjects","photons"),
				 selection       = cms.string('abs(eta)<2.4 && pt > 20'),
				 ptBins          = cms.vdouble(20,30,40,50,60,80,110,150,1000),
				 etaBins         = cms.vdouble(0.,0.5,1.,1.5,2.0,2.4),
				 ## CSV v2 wp in 76X
				 bDiscriminatorInfo = cms.VPSet(
		cms.PSet(discriminatorName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
			 wpLabel = cms.string("Loose"),
			 wpValue = cms.double(0.460)),

		cms.PSet(discriminatorName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
			 wpLabel = cms.string("Medium"),
			 wpValue = cms.double(0.800)),

		cms.PSet(discriminatorName = cms.string("pfCombinedMVAV2BJetTags"),
			 wpLabel = cms.string("Loose"),
			 wpValue = cms.double(-0.715)),

		cms.PSet(discriminatorName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
			 wpLabel = cms.string("Medium"),
			 wpValue = cms.double(0.185)),
		))

if options.addPuppiJets: ## make b-tagging efficiency maps also for puppi jets
	setattr(process,"btageffPuppi",process.btageff.clone(
			srcJets = cms.InputTag(jetPuppiCollName)))

## do the same of subjets of boosted jets
if options.addSubstructureCHS:
	setattr(process,"btageffBooostedJet",process.btageff.clone(
			srcJets = cms.InputTag(boostedJetCollection),
			useSubjets = cms.bool(True),
			etaBins         = cms.vdouble(0.,1.25,2.4),
			ptBins          = cms.vdouble(20,50,100,150,250,1000),
			bDiscriminatorInfo = cms.VPSet(
				cms.PSet(discriminatorName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
					 wpLabel = cms.string("Loose"),
					 wpValue = cms.double(0.460)),
				cms.PSet(discriminatorName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
					 wpLabel = cms.string("Medium"),
					 wpValue = cms.double(0.800)),
				)
			))

if options.addSubstructurePuppi:
	setattr(process,"btageffBooostedPuppiJet",process.btageff.clone(
			srcJets = cms.InputTag(boostedPuppiJetCollection),
			useSubjets = cms.bool(True),
			etaBins         = cms.vdouble(0.,1.25,2.4),
			ptBins          = cms.vdouble(20,50,100,150,250,1000),
			bDiscriminatorInfo = cms.VPSet(
				cms.PSet(discriminatorName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
					 wpLabel = cms.string("Loose"),
					 wpValue = cms.double(0.460)),
				cms.PSet(discriminatorName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
					 wpLabel = cms.string("Medium"),
					 wpValue = cms.double(0.800)),
				)
			))

## fast sim business
if options.isFastSIM:

	## fix LHE info
	process.tree.lheinfo = cms.InputTag("source")
	process.tree.lheRuninfo = cms.InputTag("source")

	process.gentree.lheinfo = cms.InputTag("source")
	process.gentree.lheRuninfo = cms.InputTag("source")
	
	process.tree.pileup = cms.InputTag("addPileupInfo");
	
# Set up the path
if options.dropAnalyzerDumpEDM == False:
	if (options.isMC):
		process.treePath = cms.Path(process.gentree + 					    
					    process.metFilters+
					    process.btageff+
					    process.filterHighRecoil + 
					    process.tree)
		if(options.addPuppiJets):
			process.treePath = cms.Path(process.gentree +
						    process.metFilters+
						    process.btageff+
						    process.btageffPuppi+
						    process.filterHighRecoil +
						    process.tree)

		if options.addSubstructureCHS:
			process.treePath.replace(process.btageff,process.btageff+process.btageffBooostedJet);
		if options.addSubstructurePuppi:
			process.treePath.replace(process.btageffPuppi,process.btageffPuppi+process.btageffBooostedPuppiJet);
			
	else :
		process.treePath = cms.Path(
			process.metFilters+
			process.filterHighRecoil + 
			process.tree)

processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
