import os, re
import sys
import FWCore.ParameterSet.Config as cms
  
### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

## data or MC options
options.register (
	'isMC',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate data or MC');

options.register (
	'isReMiniAOD',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate new re-miniAOD');

## MET filter options
options.register (
	'filterHighMETEvents',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate if apply or not MET filters');

options.register (
	'metCut',150.,VarParsing.multiplicity.singleton,VarParsing.varType.float,
	'met/recoil cut to be applied if filterHighMETEvents is set to true');

## HLT filter options
options.register (
	'filterOnHLT',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to indicate if apply or not trigger requirements');

options.register (
	'setHLTFilterFlag',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'flag to dump all the HLT flags to true');

## JET correction options
options.register (
	'usePrivateSQliteJEC',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'if a private SQL file with JEC to be found in test directory');

options.register (
	'usePrivateSQliteJER',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'if a private SQL file with JER to be found in test directory');

options.register (
	'applyL2L3Residuals',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'apply or not L2L3 Residual JEC on data');

options.register (
	'JECEra','Summer16_23Sep2016V4',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'JEC correction era');

## JET information
options.register (
	'addPileupJetID',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	're-compute pileup-jet id for AK4 jets');

options.register (
	'addQGLikelihood',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'compute QGLikelihood for AK4 jets');

## PUPPI Jets 
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
	'addEGMRegression',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'add new regression for photons/electrons');

## MET options
options.register (
        'useMiniAODMet',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'use the default MET in minoAOD without re-applying corrections');

options.register (
        'useMiniAODPuppiMet',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'use the default puppi MET in minoAOD without re-applying corrections');

options.register (
        'addMETSystematics',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'store met sys variation in the output tree');

options.register (
        'addPuppiMETSystematics',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'store met sys variation in the output tree');
	
options.register (
	'addBadMuonClean',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'Cleaning Bad muons from MET --> as done in 2016 re-miniAOD');
  	
options.register (
	'useOfficialMETSystematics',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'run the official tool for met uncertainty --> does a lot of things but slow .. otherwise minimal home made validated code');

options.register (
	'addMETBreakDown',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'produce the pf met breakdown in different components: pfMet, pfMetChargedHadrons, pfMetNeutralHadrons, pfMetPhotons ... etc');

## photon purity studies
options.register (
        'isPhotonPurity',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'photon purity studies --> add some more branches');

## Di-muon skim
options.register (
        'applyDiMuonFilter',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'select zmm events --> useful for response/resolution studies');

## Di-electron skim
options.register (
        'applyDiElectronFilter',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'select zee events --> useful for response/resolution studies');

## photon+jets skim
options.register (
        'applyPhotonJetsFilter',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'select gamma+jets events --> useful for response/resolution studies');

## QCD background studies
options.register (
        'isQCDTree',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'output tree layout for QCD background studies');

## Photon and electron id variables
options.register (
        'addPhotonIDVariables',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'photon id variables study --> dump id variables');

options.register (
        'addElectronIDVariables',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'electorn id variables study --> dump id variables');

## do substructure for CHS or Puppi jets
options.register (
	'addSubstructureCHS',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'run substructure algo for AK8CHS jets (Pruning, softDrop)');

options.register (
	'addSubstructurePuppi',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'run substructure algo for AK8Puppi jets (Pruning, softDrop)');

options.register (
	'useMiniAODSubstructure',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'use miniAOD AK8CHS and AK8Puppi jets instead of re-running');

## processName
options.register (
	'processName','TREE',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'process name to be considered');

options.register (
	'miniAODProcess','PAT',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'process name used for miniAOD production (target is miniAODv2)');

## specific to produce trees for trigger studies
options.register (
	'triggerName','HLT',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'process name used for miniAOD production (target is miniAODv2)');

options.register (
	'isTriggerTree',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'reduce the number of braches for trigger studies: photons and met triggers');

options.register (
	'addTriggerObjects',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
	'add L1 and HLT objects dump');

## outputFile Name
options.register (
	'outputFileName','tree.root',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'output file name created by cmsRun');

## GT to be used    
options.register (
	'globalTag','80X_dataRun2_2016SeptRepro_v7',VarParsing.multiplicity.singleton,VarParsing.varType.string,
	'gloabl tag to be uses');

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
	'nThreads',1,VarParsing.multiplicity.singleton, VarParsing.varType.int,
	'default number of threads');

## to be used when running crab jobs with local files                                                                                                                           
options.register ('isCrab',False,VarParsing.multiplicity.singleton, VarParsing.varType.bool,
		  'to be used to handle local files with crab');

## parsing command line arguments
options.parseArguments()

### check consistentcy of basic options
if options.isMC and 'dataRun2' in options.globalTag:
	options.globalTag = '80X_mcRun2_asymptotic_2016_TrancheIV_v8';
	options.JECEra  = 'Summer16_23Sep2016V4';
if options.isMC and options.applyL2L3Residuals:
	options.applyL2L3Residuals = False
	
if not options.isMC:
	options.crossSection = -1.;

if options.isMC and options.miniAODProcess != 'PAT':
	options.miniAODProcess  = 'PAT'

if options.isMC and options.isReMiniAOD :
	sys.exit("Cannot be isMC and isReMiniAOD at the same time");

print "##### Settings ######"
print "##### General #####"
print "Running with isMC                = ",options.isMC	
print "Running with isReMiniAOD         = ",options.isReMiniAOD	
print "Running with processName         = ",options.processName	
print "Running with miniAODProcess      = ",options.miniAODProcess	
print "Running with outputFileName      = ",options.outputFileName	
print "Running with globalTag           = ",options.globalTag	
print "Running with nThreads            = ",options.nThreads
print "Running with isCrab              = ",options.isCrab
print "Running with dropAnalyzerDumpEDM = ",options.dropAnalyzerDumpEDM	
print "##### Skim #####"
print "Running with filterHighMETEvents = ",options.filterHighMETEvents	
if options.filterHighMETEvents:
	print "Running with metCut              = ",options.metCut
print "Running with filterOnHLT         = ",options.filterOnHLT	
print "Running with setHLTFilterFlag    = ",options.setHLTFilterFlag
print "##### Trigger info #####"
print "Running with triggerName         = ",options.triggerName
print "Running with isTriggerTree       = ",options.isTriggerTree
print "Running with addTriggerObjects   = ",options.addTriggerObjects
print "##### QCD background trees #####"
print "Running with isQCDTree           = ",options.isQCDTree
print "##### Regular jets #####"
print "Running with JEC Era             = ",options.JECEra	
print "Running with usePrivateSQliteJEC = ",options.usePrivateSQliteJEC	
print "Running with usePrivateSQliteJER = ",options.usePrivateSQliteJER	
print "Running with applyL2L3Residuals  = ",options.applyL2L3Residuals	
print "Running with addPileupJetID      = ",options.addPileupJetID
print "Running with addQGLikelihood     = ",options.addQGLikelihood
print "Running with addPuppiJets        = ",options.addPuppiJets
print "##### Missing energy #####"
print "Running with addPuppiMET            = ",options.addPuppiMET
print "Running with useMiniAODMet          = ",options.useMiniAODMet
print "Running with useMiniAODPuppiMet     = ",options.useMiniAODPuppiMet
print "Running with addMETSystematics      = ",options.addMETSystematics
print "Running with addPuppiMETSystematics = ",options.addPuppiMETSystematics
print "Running with useOfficialMETSystematics = ",options.useOfficialMETSystematics
print "Running with addMETBreakDown        = ",options.addMETBreakDown	
print "##### Electrons/Photons #####"
print "Running with addEGMSmear            = ",options.addEGMSmear
print "Running with addEGMRegression       = ",options.addEGMRegression
print "Running with isPhotonPurity         = ",options.isPhotonPurity	
print "Running with addPhotonIDVariables   = ",options.addPhotonIDVariables
print "Running with applyPhotonJetsFilter  = ",options.applyPhotonJetsFilter
print "Running with addElectronIDVariables = ",options.addElectronIDVariables
print "Running with applyDiMuonFilter      = ",options.applyDiMuonFilter
print "Running with applyDiElectronFilter  = ",options.applyDiElectronFilter
print "##### Jet Substructure #####"
print "Running with addSubstructureCHS   = ",options.addSubstructureCHS
print "Running with addSubstructurePuppi = ",options.addSubstructurePuppi
print "Running with useMiniAODSubstructure = ",options.useMiniAODSubstructure
print "##### Generator info #####"
print "Running with useLHEWeights       = ",options.useLHEWeights
print "Running with addQCDPDFWeights    = ",options.addQCDPDFWeights
print "Running with isSignalSample      = ",options.isSignalSample
print "Running with addGenParticles     = ",options.addGenParticles
print "Running with crossSection        = ",options.crossSection
print "#####################"

if (options.addSubstructureCHS or options.addSubstructurePuppi) and options.useMiniAODSubstructure:
	sys.exit("cannot specifiy both useMiniAODSubstructure and one between addSubstructureCHS or addSubstructurePuppi --> exit");

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
process.MessageLogger.cerr.FwkReport.reportEvery = 200

## Define the input source
if options.inputFiles == []:

	process.source = cms.Source("PoolSource", 
   		 fileNames = cms.untracked.vstring())
	if not options.isMC :

		#process.source.fileNames.append('/store/data/Run2016C/MET/MINIAOD/PromptReco-v2/000/275/782/00000/327DB0D1-F93C-E611-8A11-02163E01450B.root')
		#process.source.fileNames.append('/store/data/Run2016C/MET/MINIAOD/PromptReco-v2/000/275/782/00000/768B07C9-F93C-E611-A52A-02163E0145D2.root')
		#process.source.fileNames.append('/store/data/Run2016B/SinglePhoton/MINIAOD/23Sep2016-v1/50000/004F3A63-2E84-E611-AEFB-00266CFFC43C.root')
		#process.source.fileNames.append('/store/data/Run2016G/MET/MINIAOD/23Sep2016-v1/90000/124A2693-B38A-E611-BC48-002590FC5ACC.root')
		process.source.fileNames.append('/store/data/Run2016C/JetHT/MINIAOD/03Feb2017-v1/110000/F6340B16-C1EB-E611-BE2F-008CFA197D60.root')

	else:
		#process.source.fileNames.append('file:pickevents.root')
		#process.source.fileNames.append('/store/mc/RunIISpring16MiniAODv2/DYJetsToNuNu_PtZ-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/027B63CF-D72B-E611-988C-002590A52B4A.root')
		#process.source.fileNames.append('/store/mc/RunIISpring16MiniAODv2/WJetsToLNu_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/F0025F27-AA2B-E611-9077-0CC47A4DED1A.root')
		process.source.fileNames.append('/store/mc/RunIISummer16MiniAODv2/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/4CAA7972-9BCB-E611-A395-0025905A6060.root')
		#process.source.fileNames.append('/store/mc/RunIISpring16MiniAODv2/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/7481FFE2-521A-E611-A18F-0025904C7B48.root')
		#process.source.fileNames.append('/store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/5E97F1F8-04D3-E611-9E11-549F3525DB98.root')
#		process.source.fileNames.append('/store/mc/RunIISummer16MiniAODv2/Scalar_MonoJ_NLO_Mphi-100_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/A2E89D89-52D6-E611-93DE-02163E011949.root')
		#process.source.fileNames.append('/store/mc/RunIISummer16MiniAODv2/SMM_MonoW_Mphi-1000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/8C1FE6FF-69D0-E611-91A5-FA163EF2F0A3.root')
		#process.source.fileNames.append('/store/mc/RunIISummer16MiniAODv2/Axial_MonoJ_NLO_Mphi-1000_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/80000/B6B3A7ED-05D6-E611-BA96-008CFA11113C.root')
#		process.source.fileNames.append('/store/mc/RunIISummer16MiniAODv2/Pseudo_MonoJ_NLO_Mphi-1000_Mchi-1_gSM-1p0_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/60000/D462BEC9-41DA-E611-83AD-02163E015C72.root')
#		process.source.fileNames.append('/store/mc/RunIISummer16MiniAODv2/Vector_MonoW_NLO_Mphi-1000_Mchi-300_gSM-0p25_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/747DC36C-A2D0-E611-B518-0CC47A4D75F0.root')
#		process.source.fileNames.append('/store/mc/RunIISummer16MiniAODv2/Scalar_MonoZ_NLO_Mphi-300_Mchi-50_gSM-1p0_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/9E496385-49D0-E611-948B-047D7B881DD4.root')
#		process.source.fileNames.append('/store/mc/RunIISummer16MiniAODv2/Axial_MonoW_NLO_Mphi-1750_Mchi-1_gSM-0p25_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/92A53279-1DD5-E611-A6EA-7845C4F92F87.root')
#		process.source.fileNames.append('/store/mc/RunIISummer16MiniAODv2/DMS_NNPDF30_Scalar_Mphi-350_Mchi-100_gSM-1p0_gDM-1p0_v2_13TeV-powheg/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/768DDA99-E8C5-E611-817B-A0000420FE80.root')
#		process.source.fileNames.append('/store/mc/RunIISummer16MiniAODv2/VectorMonoW_Mphi-1000_Mchi-1_gSM-0p25_gDM-1p0_v2_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/BEC73603-86CC-E611-ACA7-A0000420FE80.root')
		#process.source.fileNames.append('/store/mc/RunIISummer16MiniAODv2/SMM_MonoJ_Mphi-300_Mchi-100_gSM-1p0_gDM-1p0_13TeV-madgraph/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/1019B0E7-4ED0-E611-9C34-0025905A611C.root')

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


#process.source.eventsToProcess = cms.untracked.VEventRange("273725:305:443915886");
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#					ignoreTotal = cms.untracked.int32(1),
#					moduleMemorySummary = cms.untracked.bool(True)
#					)

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
ElectronTools(process,options.addEGMSmear,options.isMC,options.addEGMRegression, addElectronCorrection = False)

# run cut-based photon ID 
from AnalysisCode.MonoXAnalysis.PhotonTools_cff import PhotonTools
PhotonTools(process,options.addEGMSmear,options.isMC, options.addEGMRegression, addPhotonCorrection = False)

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
### gain correction always applied for the time-being
process.selectedObjects.useCalibratedElectrons = cms.bool(True)
process.selectedObjects.useCalibratedPhotons = cms.bool(True)
process.selectedObjects.addPhotonPurity = cms.bool(options.isPhotonPurity)
if hasattr(process,"correctedElectrons"):
	process.selectedObjects.calibratedElectrons = cms.InputTag("correctedElectrons")
if hasattr(process,"correctedPhotons"):
	process.selectedObjects.calibratedPhotons = cms.InputTag("correctedPhotons")
	
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
if not options.isReMiniAOD: #### don't 
	if not options.useMiniAODMet:
		metCollection = "slimmedMETs"
		metCorrector(process,jetCollName,metCollection,options.isMC,"AK4PFchs",options.applyL2L3Residuals,options.useOfficialMETSystematics,options.addMETSystematics,options.addBadMuonClean);
	if options.addPuppiMET and not options.useMiniAODPuppiMet:
		metCollectionPuppi = "slimmedMETsPuppi"
		metCorrector(process,jetPuppiCollName,metCollectionPuppi,options.isMC,"AK4PFPuppi",options.applyL2L3Residuals,options.useOfficialMETSystematics,options.addPuppiMETSystematics,options.addBadMuonClean);
#### e-gamma fix in met for re-miniAOD
elif options.isReMiniAOD and not options.useMiniAODMet:
	from PhysicsTools.PatUtils.tools.corMETFromMuonAndEG import corMETFromMuonAndEG
	from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
	runMetCorAndUncFromMiniAOD(process,
				   isData=True,
				   )
	corMETFromMuonAndEG(process,pfCandCollection="",electronCollection="slimmedElectronsBeforeGSFix",photonCollection="slimmedPhotonsBeforeGSFix",
			    corElectronCollection="slimmedElectrons", corPhotonCollection="slimmedPhotons",
			    allMETEGCorrected=True,muCorrection=False,eGCorrection=True,runOnMiniAOD=True,postfix="MuEGClean")

	process.slimmedMETsMuEGClean = process.slimmedMETs.clone()
	process.slimmedMETsMuEGClean.src = cms.InputTag("patPFMetT1MuEGClean")
	process.slimmedMETsMuEGClean.rawVariation =  cms.InputTag("patPFMetRawMuEGClean")
	process.slimmedMETsMuEGClean.t1Uncertainties = cms.InputTag("patPFMetT1%sMuEGClean")
	del process.slimmedMETsMuEGClean.caloMET


# Define all the METs corrected for lepton/photon momenta
from AnalysisCode.MonoXAnalysis.recoil_cff import recoilComputation
recoilComputation(process,options.processName,options.miniAODProcess,options.useMiniAODMet,False,options.isReMiniAOD,options.addBadMuonClean)
if options.addPuppiMET:
	recoilComputation(process,options.processName,options.miniAODProcess,options.useMiniAODPuppiMet,True,options.isReMiniAOD,options.addBadMuonClean)

### met breakdown
if options.addMETBreakDown:
	process.METBreakDown = cms.EDProducer("PATMETBreakDownProducer",
					      pfcands = cms.InputTag("packedPFCandidates"))
	
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
                                      	'keep *_*selectedObjects*_*_*',
                                      	'keep *_*t1*_*_*',
                                      	))

   process.output = cms.EndPath(process.out)

#### substructure sequence
boostedJetCollection = "";
boostedPuppiJetCollection = "";
if not options.useMiniAODSubstructure:
	from AnalysisCode.MonoXAnalysis.JetSubstructure_cff import JetSubstructure
	if options.addSubstructureCHS:
		boostedJetCollection = JetSubstructure(process,
						       options.isMC,
						       coneSize = 0.8, 
						       algo = "AK",
						       pileupMethod = "chs", 
						       selection = "pt > 190 && abs(eta) < 2.5",
						       addPruning   = True, 
						       addSoftDrop  = True, 
						       addTrimming  = False, 
						       addFiltering = False,
						       addNsubjettiness = True, 
						       addEnergyCorrelation = False, 
						       addQJets        = False,
						       addQGLikelihood = False);
		
	if options.addSubstructurePuppi:
			boostedPuppiJetCollection = JetSubstructure(process,
								    options.isMC,
								    coneSize = 0.8, 
								    algo = "AK",
								    pileupMethod = "Puppi", 
								    selection = "pt > 190 && abs(eta) < 2.5",
								    addPruning  = True, 
								    addSoftDrop = True, 
								    addTrimming = False,
								    addFiltering = False,
								    addNsubjettiness = True, 
								    addEnergyCorrelation = False, 
								    addQJets = False,
								    addQGLikelihood = False);
			
else:	
	boostedJetCollection = JetCorrector(process,"slimmedJetsAK8","AK8PFchs",options.isMC, options.applyL2L3Residuals);
	
### apply event selections
from AnalysisCode.MonoXAnalysis.applyEventFilters_cff import applyEventFilters
looseMuonPt = 10. ; tightMuonPt = 20.; 
looseElectronPt = 10.; tightElectronPt = 35.; useMVAElectronID = False;
photonPt = 50.; useMVAPhotonID = False;

applyEventFilters(process,
		  options.processName,
		  options.miniAODProcess,
		  options.filterHighMETEvents, ### if apply or not filter on MET
		  options.metCut, ## value for reoil selection
		  options.isPhotonPurity, ### in case one wants to make a specific filter
		  options.applyDiMuonFilter, ### in case one wants to apply di-muon Zmm filter
		  looseMuonPt,
		  tightMuonPt,	
		  options.applyDiElectronFilter, ### in case one wants to apply di-electron Zee filter
		  looseElectronPt,
		  tightElectronPt,
		  useMVAElectronID,
		  options.applyPhotonJetsFilter, ### in case one wants to apply single-photon filter
		  photonPt,
		  useMVAPhotonID,
		  options.isReMiniAOD,
		  options.addBadMuonClean)
	
if options.useMiniAODMet and not options.isReMiniAOD:
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

### fix the trigger label for data and MC
triggerLabel = "RECO";
if options.isMC:
	triggerLabel = "PAT";

# Make the tree 
process.tree = cms.EDAnalyzer("MonoJetTreeMaker",
			      ## gen info			     
			      isMC                   = cms.bool(options.isMC),
			      isReMiniAOD            = cms.bool(options.isReMiniAOD),
			      addBadMuonClean        = cms.bool(options.addBadMuonClean),
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
			      isTriggerTree  = cms.bool(options.isTriggerTree),
			      addTriggerObjects = cms.bool(options.addTriggerObjects),
			      triggerResults = cms.InputTag("TriggerResults", "",options.triggerName),
			      prescales      = cms.InputTag("patTrigger"),    
			      filterResults  = cms.InputTag("TriggerResults", "",triggerLabel),
			      triggerObjects = cms.InputTag("selectedPatTrigger"),
			      triggerL1EG    = cms.InputTag("caloStage2Digis"   , "EGamma"),
			      triggerL1Jet   = cms.InputTag("caloStage2Digis"   , "Jet"   ),
			      triggerL1Mu    = cms.InputTag("gmtStage2Digis"    , "Muon"  ),
			      triggerL1Sums  = cms.InputTag("caloStage2Digis"   , "EtSum" ),
			      badChargedCandidate = cms.InputTag("BadChargedCandidateFilter"),
			      badPFMuon           = cms.InputTag("BadPFMuonFilter"),
			      ## flag for QCD background studies
			      isQCDTree      = cms.bool(options.isQCDTree),
			      ## vertexes			    
			      vertices       = cms.InputTag("goodVertices"),
			      ## muons    
			      applyDiMuonFilter = cms.bool(options.applyDiMuonFilter),
			      muons          = cms.InputTag("selectedObjects","muons"),
			      tightmuons     = cms.InputTag("selectedObjects","tightmuons"),
			      highptmuons    = cms.InputTag("selectedObjects","highptmuons"),			      
			      ## electrons
			      applyDiElectronFilter = cms.bool(options.applyDiElectronFilter),
			      electrons       = cms.InputTag("selectedObjects", "electrons"),
			      looseelectrons  = cms.InputTag("selectedObjects", "looseelectrons"),
			      tightelectrons  = cms.InputTag("selectedObjects", "tightelectrons"),
			      triggerelectrons = cms.InputTag("selectedObjects","triggerelectrons"),
			      heepelectrons   = cms.InputTag("selectedObjects", "heepelectrons"),
			      mvalooseelectrons   = cms.InputTag("selectedObjects", "mvalooseelectrons"),
			      mvatightelectrons   = cms.InputTag("selectedObjects", "mvatightelectrons"),
			      ## photons --> can be matched by reference 
			      applyPhotonJetsFilter = cms.bool(options.applyPhotonJetsFilter),
			      rho             = cms.InputTag("fixedGridRhoFastjetAll"),
			      photons         = cms.InputTag("selectedObjects", "photons"),
			      mediumphotons   = cms.InputTag("selectedObjects", "mediumphotons"),
			      tightphotons    = cms.InputTag("selectedObjects", "tightphotons"),			     
			      mvaloosephotons = cms.InputTag("selectedObjects", "mvaloosephotons"),			     
			      mvatightphotons = cms.InputTag("selectedObjects", "mvatightphotons"),			     
			      photonHighPtId  = cms.InputTag("selectedObjects", "photonHighPtId"),
			      ## photon purity
			      isPhotonPurity  = cms.bool(options.isPhotonPurity),
			      photonsPurity   = cms.InputTag("selectedObjects", "photonsPurity"),
			      rndgammaiso04   = cms.InputTag("selectedObjects", "rndgammaiso04"),
			      rndgammaiso08   = cms.InputTag("selectedObjects", "rndgammaiso08"),
			      rndchhadiso04   = cms.InputTag("selectedObjects", "rndchhadiso04"),
			      rndchhadiso08   = cms.InputTag("selectedObjects", "rndchhadiso08"),
			      photonsieie = cms.InputTag("selectedObjects", "sigmaietaieta"),
			      photonPHiso = cms.InputTag("selectedObjects", "gammaiso"),
			      photonCHiso = cms.InputTag("selectedObjects", "chhadiso"),
			      photonNHiso = cms.InputTag("selectedObjects", "nhhadiso"),
			      ## taus
			      tausVLNew       = cms.InputTag("selectedObjects","tausVLNew"),
			      tausVLOld       = cms.InputTag("selectedObjects","tausVLOld"),
			      tausRawNew      = cms.InputTag("selectedObjects","tausRawNew"),
			      tausRawOld      = cms.InputTag("selectedObjects","tausRawOld"),
			      tausTightNew    = cms.InputTag("selectedObjects","tausTightNew"),
			      tausTightOld    = cms.InputTag("selectedObjects","tausTightOld"),
			      ## jets AK4
			      jets            = cms.InputTag(jetCollName),
			      addPuppiJets    = cms.bool(options.addPuppiJets),			      
			      puppijets       = cms.InputTag(jetPuppiCollName),			      
			      ## MET
			      t1met     = cms.InputTag("slimmedMETs","",options.processName),
			      t1mumet   = cms.InputTag("t1mumet"),
			      t1elmet   = cms.InputTag("t1elmet"),
			      t1phmet   = cms.InputTag("t1phmet"),
			      t1taumet  = cms.InputTag("t1taumet"),
			      ## MET
			      addMETBreakDown    = cms.bool(options.addMETBreakDown),
			      pfMetHadronHF      = cms.InputTag("METBreakDown","pfMetHadronHF"),
			      pfMetEgammaHF      = cms.InputTag("METBreakDown","pfMetEgammaHF"),
			      pfMetChargedHadron = cms.InputTag("METBreakDown","pfMetChargedHadron"),
			      pfMetNeutralHadron = cms.InputTag("METBreakDown","pfMetNeutralHadron"),
			      pfMetElectrons     = cms.InputTag("METBreakDown","pfMetElectrons"),
 			      pfMetPhotons       = cms.InputTag("METBreakDown","pfMetPhotons"),
 			      pfMetMuons         = cms.InputTag("METBreakDown","pfMetMuons"),
			      pfMetUnclustered   = cms.InputTag("METBreakDown","pfMetUnclustered"),
			      ## Puppi MET
			      addPuppiMET   = cms.bool(options.addPuppiMET),
			      puppit1met    = cms.InputTag("slimmedMETsPuppi","",options.processName),
			      puppit1mumet  = cms.InputTag("puppit1mumet"),
			      puppit1elmet  = cms.InputTag("puppit1elmet"),
			      puppit1phmet  = cms.InputTag("puppit1phmet"),
			      puppit1taumet = cms.InputTag("puppit1taumet"),
			      ## MET systematics
			      addMETSystematics = cms.bool(options.addMETSystematics),    			      
			      ## trigger filter
			      applyHLTFilter = cms.bool(options.filterOnHLT),
			      setHLTFilterFlag = cms.bool(options.setHLTFilterFlag),
			      ## clean objects
			      cleanMuonJet     = cms.bool(True),
			      cleanElectronJet = cms.bool(True),
			      cleanPhotonJet   = cms.bool(True),
			      dRCleaningAK4    = cms.double(0.4),
                              dRCleaningAK8    = cms.double(0.8),
			      jetidwp          = cms.string("loose"),
			      applypileupjetid = cms.bool(False),
			      pileupjetidwp    = cms.string("medium"),
			      btaggingCSVWP    = cms.double(0.8484),
			      btaggingMVAWP    = cms.double(0.4432),
			      minJetPtCountAK4     = cms.double(30),
			      minJetPtBveto        = cms.double(20),
			      minJetPtAK4Store     = cms.double(25),
			      ## CHS jet substructure
			      useMiniAODSubstructure = cms.bool(options.useMiniAODSubstructure),
			      addSubstructureCHS   = cms.bool(options.addSubstructureCHS),
			      boostedJetsCHS       = cms.InputTag(boostedJetCollection),
			      addSubstructurePuppi = cms.bool(options.addSubstructurePuppi),
			      boostedJetsPuppi     = cms.InputTag(boostedPuppiJetCollection),
			      ## b-tag scale factors
			      addBTagScaleFactor         = cms.bool(True),
			      bTagScaleFactorFileCSV     = cms.FileInPath('AnalysisCode/MonoXAnalysis/data/BTagScaleFactors/CSVv2_Moriond17_B_H.csv'), 	     
			      bTagScaleFactorFileMVA     = cms.FileInPath('AnalysisCode/MonoXAnalysis/data/BTagScaleFactors/cMVAv2_Moriond17_B_H.csv'), 	 
			      bTagScaleFactorFileSubCSV  = cms.FileInPath('AnalysisCode/MonoXAnalysis/data/BTagScaleFactors/subjet_CSVv2_Moriond17_B_H.csv'),
			      ## photon id
			      addPhotonIDVariables = cms.bool(options.addPhotonIDVariables),
			      photonIDCollection   = cms.InputTag("slimmedPhotons"),
			      addElectronIDVariables = cms.bool(options.addElectronIDVariables),
			      electronIDCollection   = cms.InputTag("slimmedElectrons")
			      )

if options.useMiniAODMet:
	process.tree.t1met = cms.InputTag("slimmedMETs","",options.miniAODProcess)
	process.tree.puppit1met = cms.InputTag("slimmedMETsPuppi","",options.miniAODProcess)
if options.useMiniAODPuppiMet:
	process.tree.puppit1met = cms.InputTag("slimmedMETsPuppi","",options.miniAODProcess)

if options.isReMiniAOD:
	process.tree.fakeMuonCandidates = cms.InputTag("packedPFCandidatesDiscarded");
	## re-run on the fly EGMUclean met
	if not options.useMiniAODMet:
		process.tree.t1met = cms.InputTag("slimmedMETsMuEGClean","",options.processName);
	else:
		process.tree.t1met = cms.InputTag("slimmedMETsMuEGClean","",options.miniAODProcess);

	process.tree.t1metEGClean = cms.InputTag("slimmedMETsEGClean","",options.miniAODProcess);
	## whether bad muons done again or taken from miniAOD
	if not options.addBadMuonClean:
		process.tree.t1metMuClean = cms.InputTag("slimmedMETs","",options.miniAODProcess);
	else:
		process.tree.t1metMuClean = cms.InputTag("slimmedMETsMuClean");
	process.tree.t1metOriginal = cms.InputTag("slimmedMETsUncorrected","",options.miniAODProcess);

	process.tree.t1mumet = cms.InputTag("t1mumet");
	process.tree.t1mumetEGClean = cms.InputTag("t1mumetEGClean");
	process.tree.t1mumetMuClean = cms.InputTag("t1mumetMuClean");

	process.tree.t1elmet = cms.InputTag("t1elmet");
	process.tree.t1elmetEGClean = cms.InputTag("t1elmetEGClean");
	process.tree.t1elmetMuClean = cms.InputTag("t1elmetMuClean");

	process.tree.t1phmet = cms.InputTag("t1phmet");
	process.tree.t1phmetEGClean = cms.InputTag("t1phmetEGClean");
	process.tree.t1phmetMuClean = cms.InputTag("t1phmetMuClean");

	process.tree.t1taumet = cms.InputTag("t1taumet");
	process.tree.t1taumetEGClean = cms.InputTag("t1taumetEGClean");
	process.tree.t1taumetMuClean = cms.InputTag("t1taumetMuClean");	

if options.isMC and options.addBadMuonClean:
	process.tree.t1met = cms.InputTag("slimmedMETsMuClean");
	process.tree.t1metOriginal = cms.InputTag("slimmedMETs");

if options.addBadMuonClean:
	process.patCaloMet.metSource = cms.InputTag("metrawCaloMuClean")

###### sys scaling objects 
process.tree.jetsJESUp = cms.InputTag("");
process.tree.jetsJESDw = cms.InputTag("");
process.tree.jetsJER   = cms.InputTag("");
process.tree.puppijetsJESUp =  cms.InputTag("");
process.tree.puppijetsJESDw =  cms.InputTag("");
process.tree.puppijetsJER   =  cms.InputTag("");

if options.addMETSystematics:
	if options.useOfficialMETSystematics :
		process.tree.jetsJESUp = cms.InputTag("shiftedPatJetEnUp")
		process.tree.jetsJESDw = cms.InputTag("shiftedPatJetEnDown")
		process.tree.jetsJER   = cms.InputTag("patSmearedJets")
	else:
		process.tree.jetsJESUp = cms.InputTag("metSysProducer",jetCollName+"EnUp")
		process.tree.jetsJESDw = cms.InputTag("metSysProducer",jetCollName+"EnDown")
		process.tree.jetsJER   = cms.InputTag("metSysProducer",jetCollName+"Smear")

if options.addPuppiMETSystematics:
	
	if options.addPuppiJets and options.addPuppiMET and options.useOfficialMETSystematics:
		process.tree.puppijetsJESUp = cms.InputTag("shiftedPatJetEnUpPuppi")
		process.tree.puppijetsJESDw = cms.InputTag("shiftedPatJetEnDownPuppi")
		process.tree.puppijetsJER   = cms.InputTag("patSmearedJetsPuppi")
		
	if options.addPuppiJets and options.addPuppiMET and not options.useOfficialMETSystematics:
		process.tree.puppijetsJESUp = cms.InputTag("metSysProducerPuppi",jetPuppiCollName+"EnUp")
		process.tree.puppijetsJESDw = cms.InputTag("metSysProducerPuppi",jetPuppiCollName+"EnDown")
		process.tree.puppijetsJER   = cms.InputTag("metSysProducerPuppi",jetPuppiCollName+"Smear")



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
				 ptBins          = cms.vdouble(20,25,30,40,50,60,80,120,160,200,250,350,450,600,800),
				 etaBins         = cms.vdouble(0.,0.5,1.,1.5,2.0,2.4),
				 ## CSV v2 wp in 76X
				 bDiscriminatorInfo = cms.VPSet(
		cms.PSet(discriminatorName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
			 wpLabel = cms.string("Loose"),
			 wpValue = cms.double(0.5426)),

		cms.PSet(discriminatorName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
			 wpLabel = cms.string("Medium"),
			 wpValue = cms.double(0.8484)),

		cms.PSet(discriminatorName = cms.string("pfCombinedMVAV2BJetTags"),
			 wpLabel = cms.string("Loose"),
			 wpValue = cms.double(-0.5884)),

		cms.PSet(discriminatorName = cms.string("pfCombinedMVAV2BJetTags"),
			 wpLabel = cms.string("Medium"),
			 wpValue = cms.double(0.4432)),
		))



process.taueff = cms.EDAnalyzer("TauTaggingEfficiencyTreeMaker",
				 directoryName = cms.string("tautagEff"),
				 srcTaus       = cms.InputTag("slimmedTaus"),
				 dRClean       = cms.double(0.4),
				 cleanMuonJet  = cms.bool(True),
				 srcMuons      = cms.InputTag("selectedObjects","muons"),
				 cleanElectronJet = cms.bool(True),
				 srcElectrons    = cms.InputTag("selectedObjects","electrons"),
				 cleanPhotonJet  = cms.bool(True),
				 srcPhotons      = cms.InputTag("selectedObjects","photons"),
				 selection       = cms.string('abs(eta) < 2.3 && pt > 18'),
				 ptBins          = cms.vdouble(18,21,25,30,40,50,60,80,120,160,200,250,350,450,600,800),
				 etaBins         = cms.vdouble(0.,0.5,1.,1.5,2.0,2.3),
				 ## CSV v2 wp in 76X
				 tauDiscriminatorInfo = cms.VPSet(
		cms.PSet(discriminatorName = cms.string("VLooseTauOldDM"),
			 wpLabel = cms.string("byVLooseIsolationMVArun2v1DBoldDMwLT"),
			 decayModeFinding = cms.string("decayModeFinding")),

		cms.PSet(discriminatorName = cms.string("VLooseTauNewDM"),
			 wpLabel = cms.string("byVLooseIsolationMVArun2v1DBnewDMwLT"),
			 decayModeFinding = cms.string("decayModeFindingNewDMs"))
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

# Set up the path
if options.dropAnalyzerDumpEDM == False:
	if options.isMC:
		process.treePath = cms.Path(process.gentree + 					    
					    process.metFilters+
					    process.btageff+
					    process.taueff+
					    getattr(process,"eventFilters")+
					    process.tree)
		if(options.addPuppiJets):
			process.treePath = cms.Path(process.gentree +
						    process.metFilters+
						    process.btageff+
						    process.taueff+
						    process.btageffPuppi+
						    getattr(process,"eventFilters")+
						    process.tree)		

		if options.addSubstructureCHS:
			process.treePath.replace(process.btageff,process.btageff+process.btageffBooostedJet);
		if options.addSubstructurePuppi:
			process.treePath.replace(process.btageffPuppi,process.btageffPuppi+process.btageffBooostedPuppiJet);
		if options.addMETBreakDown:
			process.treePath.replace(getattr(process,"process.filterHighRecoil"),getattr(process,"process.filterHighRecoil")+ process.METBreakDown)
			
	else :
		process.treePath = cms.Path(
			process.metFilters+
			getattr(process,"eventFilters")+
			process.tree)

		if options.addMETBreakDown:
			process.treePath.replace(getattr(process,"process.filterHighRecoil"), getattr(process,"process.filterHighRecoil")+process.METBreakDown)


processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
