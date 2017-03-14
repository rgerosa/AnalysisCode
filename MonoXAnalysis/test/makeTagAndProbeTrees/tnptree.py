import FWCore.ParameterSet.Config as cms

### CMSSW command line parameter parser                                                                                                                                       
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

## data or MC options                                                                                                                                                        
options.register (
        'isMC',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'flag to indicate data or MC');

options.register (
        'globalTag','80X_dataRun2_Prompt_v14',VarParsing.multiplicity.singleton,VarParsing.varType.string,
        'gloabl tag to be uses');

options.register (
        'usePrivateSQliteJEC',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'if a private SQL file with JEC to be found in test directory');

options.register (
        'JECEra','Summer16_23Sep2016V3',VarParsing.multiplicity.singleton,VarParsing.varType.string,
        'JEC correction era');

options.register (
        'isCrab',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
        'if a private SQL file with JEC to be used via crab');


## parsing command line arguments                                                                                                                                            
options.parseArguments()

if options.isMC and 'dataRun2' in options.globalTag:
        options.globalTag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6';

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
                '/store/mc/RunIISpring16MiniAODv2/DYToEE_NNPDF30_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/88BF817F-FA1A-E611-B23E-14187733AD81.root'
                )
        else:
            process.source.fileNames.append(
#		    'root://xrootd-cms.infn.it:1194//store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/274/094/00000/5C319205-6425-E611-BBF4-02163E011F60.root',
		    '/store/data/Run2016B/DoubleMuon/MINIAOD/PromptReco-v2/000/273/554/00000/AA246637-E61F-E611-A971-02163E01187E.root'
#		    '/store/data/Run2016B/SingleElectron/MINIAOD/PromptReco-v2/000/273/728/00000/221C84FC-F620-E611-8A0A-02163E013752.root'
                )
            
else:
    process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(options.inputFiles))

# Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", fileName = cms.string("tnptree.root"))

# Set the global tag depending on the sample type
process.GlobalTag.globaltag = options.globalTag

## Setup the private SQLite -- Ripped from PhysicsTools/PatAlgos/test/corMETFromMiniAOD.py                                                                                                            
from AnalysisCode.MonoXAnalysis.JECConfiguration_cff import JECConfiguration
## connect to a local SQLite file or take corrections from GT                                                                                                                                         
JECConfiguration(process,options.usePrivateSQliteJEC,options.JECEra,options.isMC,True,options.isCrab)


# Select good primary vertices
process.goodVertices = cms.EDFilter("VertexSelector",
				    src = cms.InputTag("offlineSlimmedPrimaryVertices"),
				    cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
				    filter = cms.bool(True)
				    )

# Probe muons --> definition ->even looser than the loose muon ID
process.probemuons = cms.EDFilter("PATMuonSelector",
				  src = cms.InputTag("slimmedMuons"),
				  cut = cms.string("pt > 10 && abs(eta) < 2.4 && (isStandAloneMuon || isTrackerMuon)"),
				  filter = cms.bool(True)  
				  )
	

# Probe electron is just a reconstructed electron
process.probeelectrons = cms.EDFilter("PATElectronSelector",
				      src = cms.InputTag("slimmedElectrons"),
				      cut = cms.string("pt > 10 && abs(eta) < 2.5"),
				      filter = cms.bool(True)  
				      )
	

# Electron ValueMaps for identification
from AnalysisCode.MonoXAnalysis.ElectronTools_cff import ElectronTools
ElectronTools(process,False,options.isMC,False)
process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag("probeelectrons")
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag("probeelectrons"); 
process.electronRegressionValueMapProducer.srcMiniAOD = cms.InputTag("probeelectrons");

# Photon ValueMaps for identification
from AnalysisCode.MonoXAnalysis.PhotonTools_cff import PhotonTools
PhotonTools(process,False,options.isMC,False)

##### set of single muons triggers --> no eta restricted path are neeeded 
tagmuontriggernames = cms.vstring([
		"HLT_IsoMu20_v*",
		"HLT_IsoMu22_v*",
		"HLT_IsoMu24_v*",
		"HLT_IsoTkMu20*",
		"HLT_IsoTkMu22*",
		"HLT_IsoTkMu24*",
		"HLT_Mu50_v*",
		"HLT_TkMu50_v*",
		"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
		"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
		"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
		"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"
		])

tagelectrontriggernames = cms.vstring([
		"HLT_Ele24_eta2p1_WPLoose_Gsf_v*",
		"HLT_Ele25_eta2p1_WPTight_Gsf_v*",
		"HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
		"HLT_Ele27_eta2p1_WPTight_Gsf_v*",
		"HLT_Ele27_WPTight_Gsf_v*",
		"HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
		"HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"
		])    

# Make the ValueMap for muon tight ID -- cannot pass it through a string selection due to the vertex argument
process.probeinfo = cms.EDProducer("LeptonTnPInfoProducer",
				   ## probes
				   muons     = cms.InputTag("probemuons"),     
				   electrons = cms.InputTag("probeelectrons"), 
				   photons   = cms.InputTag("slimmedPhotons"), 
				   ## in case of reco electron efficiency used to match probe-electorns with reco gsf ones
				   electronsFullCollection = cms.InputTag("slimmedElectrons"),
				   ## additional event info
				   geninfo   = cms.InputTag("generator"),
				   vertices  = cms.InputTag("goodVertices"), 
				   ## trigger collections
				   triggerobjects = cms.InputTag("selectedPatTrigger"),
				   triggerResults = cms.InputTag("TriggerResults", "", "HLT"),				   
				   #### Muon information for identification --> pt cut and matching with trigger info				   
				   tagloosemuons = cms.PSet(
		isocut  = cms.double(0.25),
		ptcut   = cms.double(6), ### small pt threhsold to include low pt muons for double muon trigger measurements
		etacut  = cms.double(2.4)
		),
				   tagtightmuons = cms.PSet(
		isocut  = cms.double(0.15),
		ptcut   = cms.double(6), ### small pt threhsold to include low pt muons for double muon trigger measurements
		etacut  = cms.double(2.4)
		),
				   tagmuontriggermatch = cms.PSet(
		tagmuontrigmatchdR = cms.double(0.3),
		requiremuonhlt  = cms.bool(True),
		tagmuontriggers = tagmuontriggernames
		),
				   #### Electron information for identification --> pt cut and matching with trigger info
				   tagelectrons = cms.PSet(
		tagelectronptcut   = cms.double(35),
		tagelectronetacut  = cms.double(2.5),
		applyPVSelection   = cms.bool(True),
		d0Barrel = cms.double(0.05),
                d0Endcap = cms.double(0.10),
                dzBarrel = cms.double(0.10),
                dzEndcap = cms.double(0.20),   		
		),
				   tagelectrontriggermatch = cms.PSet(
		tagelectrontrigmatchdR = cms.double(0.3),
		requireelectronhlt  = cms.bool(True),
		tagelectrontriggers = tagelectrontriggernames,
		),
				   ### electorn ID
				   electronvetoid   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"),
				   electronlooseid  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose"),
				   electronmediumid = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium"),
				   electrontightid  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight"),
				   electronhltsafeid  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1"),
				   electronmvalooseid = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90"),
				   electronmvatightid = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"),
				   ### photon id
				   photonlooseid  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose"),
				   photonmediumid = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium"),
				   photontightid  = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight"),
    				   photonmvalooseid = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90"),
    				   photonmvatightid = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp80"),				   
				   )

if options.isMC:
	process.probeinfo.tagmuontriggermatch.requiremuonhlt = cms.bool(False)
	process.probeinfo.tagelectrontriggermatch.requireelectronhlt = cms.bool(False)

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
				 checkCharge = cms.bool(True)
				 )

process.electrontnp = cms.EDProducer("CandViewShallowCloneCombiner",
				     decay = cms.string("tagelectrons@+ probeelectrons@-"),
				     cut   = cms.string("60 < mass < 120 & charge=0"),
				     checkCharge = cms.bool(True)
				     )


process.photontnp = cms.EDProducer("CandViewShallowCloneCombiner",
				   decay = cms.string("tagelectrons@+ slimmedPhotons@-"),
				   cut   = cms.string("60 < mass < 120"),
				   checkCharge = cms.bool(False)
				   )

# Make the TnP tree --> official tag and probe analyzer
process.muontnptree = cms.EDAnalyzer("TagProbeFitTreeProducer",
				     tagProbePairs = cms.InputTag("muontnp"), ## pairs for Z->ll candidate
				     arbitration   = cms.string("BestMass"),
				     massForArbitration = cms.double(91.186),				     
				     addRunLumiInfo        = cms.bool (True),
				     variables = cms.PSet(## probe info
		pt   = cms.string("pt()"),
		eta  = cms.string("eta()"),
		abseta = cms.string("abs(eta())"),
		phi  = cms.string("phi()"),
		probe_mass = cms.string("mass()"),
		charge = cms.string("charge()"),
		nvtx      = cms.InputTag("probeinfo", "munvtxmap"), ## store nvtx
		wgt       = cms.InputTag("probeinfo", "muwgtmap"),   ## store gen weight		
		nstation  = cms.string("numberOfMatchedStations()"),#		
		chi2      = cms.InputTag("probeinfo","muchi2map"),
		nvalidhit = cms.InputTag("probeinfo","munvalidhitmap"),
		npixelhit = cms.InputTag("probeinfo","munpixelhitmap"),
		ntrackerlayerhit = cms.InputTag("probeinfo","muntrackerlayermap"),
		dxy     = cms.InputTag("probeinfo","mudxymap"),
		dz      = cms.InputTag("probeinfo","mudzmap"),
		nhiso   = cms.string("pfIsolationR04().sumNeutralHadronEt"),
		emiso   = cms.string("pfIsolationR04().sumPhotonEt"),
		puiso   = cms.string("pfIsolationR04().sumPUPt"),
		chiso   = cms.string("pfIsolationR04().sumChargedHadronPt"),
		),
				     flags = cms.PSet( ## flags 
		pfid      = cms.string("isPFMuon"),
		globalid  = cms.string("isGlobalMuon"),
		standaloneid = cms.string("isStandAloneMuon"),
		trackerid    = cms.string("isTrackerMuon"),
		hltmu20   = cms.InputTag("probeinfo", "hltmu20muonrefs"),   ## if it belongs to hltmu20
		hlttkmu20 = cms.InputTag("probeinfo", "hlttkmu20muonrefs"), ## if it belongs to isotk20
		hltmu22   = cms.InputTag("probeinfo", "hltmu22muonrefs"),   ## if it belongs to hltmu22
		hlttkmu22 = cms.InputTag("probeinfo", "hlttkmu22muonrefs"), ## if it belongs to isotk22
		hltmu24   = cms.InputTag("probeinfo", "hltmu24muonrefs"),   ## if it belongs to hltmu22
		hlttkmu24 = cms.InputTag("probeinfo", "hlttkmu24muonrefs"), ## if it belongs to isotk22
		hltmu50   = cms.InputTag("probeinfo", "hltmu50muonrefs"),   ## if it belongs to hltmu22
		hlttkmu50 = cms.InputTag("probeinfo", "hlttkmu50muonrefs"), ## if it belongs to isotk22
		hltmu     = cms.InputTag("probeinfo", "hltmumuonrefs"),     ## if it belongs to hltmu22		
		hlttkmu   = cms.InputTag("probeinfo", "hlttkmumuonrefs"),   ## if it belongs to isotk22
		hltmu17mu8_leg17 = cms.InputTag("probeinfo", "hltmu17mu8Leg17"),
		hltmu17mu8_leg8 = cms.InputTag("probeinfo", "hltmu17mu8Leg8"),
		hltmu17tkmu8_leg17 = cms.InputTag("probeinfo", "hltmu17tkmu8Leg17"),
		hltmu17tkmu8_leg8 = cms.InputTag("probeinfo", "hltmu17tkmu8Leg8"),
		hltmu17mu8dz_leg17 = cms.InputTag("probeinfo", "hltmu17mu8dzLeg17"),
		hltmu17mu8dz_leg8 = cms.InputTag("probeinfo", "hltmu17mu8dzLeg8"),
		hltmu17tkmu8dz_leg17 = cms.InputTag("probeinfo", "hltmu17tkmu8dzLeg17"),
		hltmu17tkmu8dz_leg8 = cms.InputTag("probeinfo", "hltmu17tkmu8dzLeg8"),
		looseid   = cms.InputTag("probeinfo", "loosemuonrefs"),     ## if pass the loose
		tightid   = cms.InputTag("probeinfo", "tightmuonrefs")     ## if pass the tight
		),
				     tagVariables   =  cms.PSet(### tag info
		pt   = cms.string("pt()"),
		eta  = cms.string("eta()"),
		abseta = cms.string("abs(eta())"),
		phi  = cms.string("phi()"),
		mass = cms.string("mass()"),
		charge = cms.string("charge()"),
		nstation  = cms.string("numberOfMatchedStations()"),#                                                                                                                  
		nhiso   = cms.string("pfIsolationR04().sumNeutralHadronEt"),
		emiso   = cms.string("pfIsolationR04().sumPhotonEt"),
		puiso   = cms.string("pfIsolationR04().sumPUPt"),
		chiso   = cms.string("pfIsolationR04().sumChargedHadronPt"),
		),
				     tagFlags = cms.PSet( ## flags 
		pfid      = cms.string("isPFMuon"),
		globalid  = cms.string("isGlobalMuon"),
		standaloneid = cms.string("isStandAloneMuon"),
		trackerid    = cms.string("isTrackerMuon"),
		),
				     
				     pairVariables = cms.PSet(
		zmass = cms.string("mass()"),
		zpt  = cms.string("pt()"),
		zeta = cms.string("eta()"),
		zphi = cms.string("phi()")
		),
				     pairFlags = cms.PSet(),				     
				     isMC = cms.bool(options.isMC),
				     allProbes = cms.InputTag("probemuons")		
				     )

				     
process.electrontnptree = cms.EDAnalyzer("TagProbeFitTreeProducer",
					 tagProbePairs = cms.InputTag("electrontnp"),
					 arbitration   = cms.string("BestMass"),
					 massForArbitration = cms.double(91.186),
					 addRunLumiInfo        = cms.bool (True),
					 variables = cms.PSet(
		pt  = cms.string("pt()"),
		eta = cms.string("superCluster().eta()"),
		abseta    = cms.string("abs(superCluster().eta())"),
		phi       = cms.string("phi()"),
		probe_mass      = cms.string("mass()"),
		charge    = cms.string("charge()"),
		nvtx      = cms.InputTag("probeinfo", "elnvtxmap"),
		wgt       = cms.InputTag("probeinfo", "elwgtmap"),
		energy    = cms.string("energy()"),
		rawenergy = cms.string("superCluster().rawEnergy()"),
		hovere    = cms.string("hadronicOverEm()"),
		sigietaieta = cms.string("full5x5_sigmaIetaIeta()"),
		chiso     = cms.string("pfIsolationVariables().sumChargedHadronPt"),
		nhiso     = cms.string("pfIsolationVariables().sumNeutralHadronEt"),
		emiso     = cms.string("pfIsolationVariables().sumPhotonEt"),
		trackpt   = cms.string("gsfTrack().pt()"),
		eop       = cms.string("abs(1-eSuperClusterOverP())/ecalEnergy()"),
		dphiIn    = cms.string("deltaPhiSuperClusterTrackAtVtx()"),
		detaIn    = cms.string("deltaEtaSuperClusterTrackAtVtx()"),
		missHit   = cms.string("gsfTrack().hitPattern().numberOfHits(\'MISSING_INNER_HITS\')"),
		conversion = cms.string("passConversionVeto()"),
		dxy  = cms.InputTag("probeinfo","eldxymap"),
		dz   = cms.InputTag("probeinfo","eldzmap"),
		),
					 flags = cms.PSet(
		hltele24eta2p1wpl = cms.InputTag("probeinfo", "hltele24eta2p1wplooseelectronrefs"),       
		hltele25eta2p1wpt = cms.InputTag("probeinfo", "hltele25eta2p1wptightelectronrefs"),       
		hltele27eta2p1wpl = cms.InputTag("probeinfo", "hltele27eta2p1wplooseelectronrefs"),       
		hltele27eta2p1wpt = cms.InputTag("probeinfo", "hltele27eta2p1wptightelectronrefs"),       
		hltele27wpt       = cms.InputTag("probeinfo", "hltele27wptightelectronrefs"),       
		hltele105 = cms.InputTag("probeinfo", "hltele105electronrefs"),       
		hltele115 = cms.InputTag("probeinfo", "hltele115electronrefs"),       
		hltele    = cms.InputTag("probeinfo", "hltelelectronrefs"),       
		vetoid    = cms.InputTag("probeinfo", "vetoelectronrefs"),
		looseid   = cms.InputTag("probeinfo", "looseelectronrefs"),
		mediumid  = cms.InputTag("probeinfo", "mediumelectronrefs"),
		tightid   = cms.InputTag("probeinfo", "tightelectronrefs"),
		hltsafeid = cms.InputTag("probeinfo", "hltsafeelectronrefs"),   
		mvalooseid  = cms.InputTag("probeinfo", "mvalooseelectronrefs"),
		mvatightid  = cms.InputTag("probeinfo", "mvatightelectronrefs")
		),
					 tagVariables   =  cms.PSet(### tag info
		pt   = cms.string("pt()"),
		eta  = cms.string("eta()"),
		abseta = cms.string("abs(eta())"),
		phi  = cms.string("phi()"),
		mass = cms.string("mass()"),
		charge = cms.string("charge()"),
		energy    = cms.string("energy()"),
		rawenergy = cms.string("superCluster().rawEnergy()"),
		hovere    = cms.string("hadronicOverEm()"),
		sigietaieta = cms.string("full5x5_sigmaIetaIeta()"),
		chiso     = cms.string("pfIsolationVariables().sumChargedHadronPt"),
		nhiso     = cms.string("pfIsolationVariables().sumNeutralHadronEt"),
		emiso     = cms.string("pfIsolationVariables().sumPhotonEt"),
		trackpt   = cms.string("gsfTrack().pt()"),
		eop       = cms.string("abs(1-eSuperClusterOverP())/ecalEnergy()"),
		dphiIn    = cms.string("deltaPhiSuperClusterTrackAtVtx()"),
		detaIn    = cms.string("deltaEtaSuperClusterTrackAtVtx()"),
		missHit   = cms.string("gsfTrack().hitPattern().numberOfHits(\'MISSING_INNER_HITS\')"),
		conversion = cms.string("passConversionVeto()"),
		),
					 tagFlags = cms.PSet(),
					 pairVariables = cms.PSet(
		zmass = cms.string("mass()"),
		zpt  = cms.string("pt()"),
		zeta = cms.string("eta()"),
		zphi = cms.string("phi()")
		),
					 pairFlags  = cms.PSet(),		
					 isMC = cms.bool(options.isMC),
					 allProbes = cms.InputTag("probeelectrons")
					 )


process.photontnptree = cms.EDAnalyzer("TagProbeFitTreeProducer",
				       tagProbePairs = cms.InputTag("photontnp"),
				       arbitration   = cms.string("BestMass"),
				       massForArbitration = cms.double(91.186),
				       addRunLumiInfo        = cms.bool (True),
				       variables = cms.PSet(
		pt  = cms.string("pt()"),
		eta = cms.string("superCluster().eta()"),
		abseta = cms.string("abs(superCluster().eta())"),
		phi       = cms.string("phi()"),
		probe_mass  = cms.string("mass()"),
		charge      = cms.string("charge()"),
		nvtx      = cms.InputTag("probeinfo", "phnvtxmap"),
		wgt       = cms.InputTag("probeinfo", "phwgtmap"),
		energy    = cms.string("energy()"),
		rawenergy = cms.string("superCluster().rawEnergy()"),
		hovere    = cms.string("hadTowOverEm()"),
		eleveto   = cms.string("passElectronVeto()"),
		sigietaieta = cms.InputTag("photonIDValueMapProducer", "phoFull5x5SigmaIEtaIEta"),
		chiso     = cms.InputTag("photonIDValueMapProducer", "phoPhotonIsolation"),
		nhiso     = cms.InputTag("photonIDValueMapProducer", "phoChargedIsolation"),
		emiso     = cms.InputTag("photonIDValueMapProducer", "phoNeutralHadronIsolation"),
		),
				       flags = cms.PSet(
		hltele24eta2p1wpl = cms.InputTag("probeinfo", "hltele24eta2p1wplooseelectronrefs"),       
		hltele25eta2p1wpt = cms.InputTag("probeinfo", "hltele25eta2p1wptightelectronrefs"),       
		hltele27eta2p1wpl = cms.InputTag("probeinfo", "hltele27eta2p1wplooseelectronrefs"),       
		hltele27eta2p1wpt = cms.InputTag("probeinfo", "hltele27eta2p1wptightelectronrefs"),       
		hltele27wpt       = cms.InputTag("probeinfo", "hltele27wptightelectronrefs"),       
		hltele105 = cms.InputTag("probeinfo", "hltele105electronrefs"),       
		hltele115 = cms.InputTag("probeinfo", "hltele115electronrefs"),       
		hltele    = cms.InputTag("probeinfo", "hltelelectronrefs"),       
		looseid   = cms.InputTag("probeinfo", "loosephotonrefs"),
		mediumid  = cms.InputTag("probeinfo", "mediumphotonrefs"),
		tightid   = cms.InputTag("probeinfo", "tightphotonrefs"),
		mvalooseid  = cms.InputTag("probeinfo", "mvaloosephotonrefs"),
		mvatightid  = cms.InputTag("probeinfo", "mvatightphotonrefs"),
		recoelectronmatch = cms.InputTag("probeinfo","recoelectronmatch")
		),
				       tagVariables   =  cms.PSet(### tag info
		pt   = cms.string("pt()"),
		eta  = cms.string("eta()"),
		abseta = cms.string("abs(eta())"),
		phi  = cms.string("phi()"),
		mass = cms.string("mass()"),
		charge = cms.string("charge()"),
		energy    = cms.string("energy()"),
		rawenergy = cms.string("superCluster().rawEnergy()"),
		hovere    = cms.string("hadronicOverEm()"),
		sigietaieta = cms.string("full5x5_sigmaIetaIeta()"),
		chiso     = cms.string("pfIsolationVariables().sumChargedHadronPt"),
		nhiso     = cms.string("pfIsolationVariables().sumNeutralHadronEt"),
		emiso     = cms.string("pfIsolationVariables().sumPhotonEt"),
		),
				       tagFlags = cms.PSet(),
				       pairVariables = cms.PSet(
		zmass = cms.string("mass()"),
		zpt  = cms.string("pt()"),
		zeta = cms.string("eta()"),
		zphi = cms.string("phi()")
		),				     
				       pairFlags = cms.PSet(),
				       isMC = cms.bool(options.isMC),
				       allProbes = cms.InputTag("slimmedPhotons")
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

    process.mcphprobematch = cms.EDProducer("MCTruthDeltaRMatcherNew",
        matchPDGId  = cms.vint32(11),
        src         = cms.InputTag("slimmedPhotons"),
        distMin     = cms.double(0.1),
        matched     = cms.InputTag("prunedGenParticles"),
        checkCharge = cms.bool(False)
    )

    #### muon gen info
    process.muontnptree.tagMatches   = cms.InputTag("mcmutagmatch")
    process.muontnptree.probeMatches = cms.InputTag("mcmuprobematch")
    process.muontnptree.motherPdgId  = cms.int32(23)
    process.muontnptree.makeMCUnbiasTree = cms.bool(True)
    process.muontnptree.checkMotherInUnbiasEff = cms.bool(True)
    
    process.electrontnptree.tagMatches = cms.InputTag("mceltagmatch")
    process.electrontnptree.probeMatches = cms.InputTag("mcelprobematch")
    process.electrontnptree.motherPdgId = cms.int32(23)
    process.electrontnptree.makeMCUnbiasTree = cms.bool(True)
    process.electrontnptree.checkMotherInUnbiasEff = cms.bool(True)

    process.photontnptree.tagMatches = cms.InputTag("mceltagmatch")
    process.photontnptree.probeMatches = cms.InputTag("mcphprobematch")
    process.photontnptree.motherPdgId = cms.int32(23)
    process.photontnptree.makeMCUnbiasTree = cms.bool(True)
    process.photontnptree.checkMotherInUnbiasEff = cms.bool(True)


process.muonPath     = cms.Path(
	process.muontnp*
	process.muontnptree)

process.electronPath     = cms.Path(
	process.electrontnp*
	process.electrontnptree)

process.photonPath = cms.Path(
        process.photontnp*
        process.photontnptree)

processDumpFile = open('processDump.py', 'w')
print >> processDumpFile, process.dumpPython()
