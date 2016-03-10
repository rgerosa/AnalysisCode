import os, copy
import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
from RecoMET.METProducers.PFMET_cfi import pfMet

## generic function that corrects jet and MET given a JEC
def metCorrector(process,jetCollection,metCollection,isMC,payloadName,applyL2L3Residuals,addMETSystematics):

	## propagation on missing energy + full systematics
	if "Puppi" in metCollection or "PUPPI" in metCollection:
		postfix = "Puppi"
	else:
		postfix = ""
			
	## extract genMet		
	if not hasattr(process,"genMetExtractor") and isMC:
		setattr(process,"genMetExtractor",
			cms.EDProducer("GenMETExtractor",
				       metSource = cms.InputTag(metCollection,"","@skipCurrentProcess")))			   
		
	## redo raw PF met for both collections
	if not hasattr(process,"pfMet"+postfix): 
		setattr(process,"pfMet"+postfix, pfMet.clone( 
				src = cms.InputTag("packedPFCandidates"),
				alias = cms.string('pfMet'+postfix)))

	## re-cast PFMets into PAT objects
	addMETCollection(process, labelName='patPFMet'+postfix, metSource='pfMet'+postfix)

	if isMC:
		getattr(process,"patPFMet"+postfix).addGenMET = cms.bool(True)
		getattr(process,"patPFMet"+postfix).genMETSource = cms.InputTag("genMetExtractor")
	else:
		getattr(process,"patPFMet"+postfix).addGenMET = cms.bool(False)

		 
	## derive type-I corrector object
	setattr(process,"patPFMetT1Corr"+postfix, 
		cms.EDProducer("PATPFJetMETcorrInputProducer",
			       isMC = cms.bool(isMC),
			       offsetCorrLabel = cms.InputTag("L1FastJet"),
			       jetCorrLabel = cms.InputTag("L3Absolute"), ## info embedded in the jet object
			       jetCorrLabelRes = cms.InputTag("L2L3Residual"), ## info embedded in the jet object 
			       skipEM = cms.bool(True),
			       skipEMfractionThreshold = cms.double(0.9),
			       skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
			       skipMuons = cms.bool(True),
			       src = cms.InputTag(jetCollection),
			       type1JetPtThreshold = cms.double(15.0),
			       type2ExtraCorrFactor = cms.double(1.0),
			       type2ResidualCorrEtaMax = cms.double(9.9),
			       type2ResidualCorrLabel = cms.InputTag(""),
			       type2ResidualCorrOffset = cms.double(0.0)
			       ))
	

	## if not residuals
	if applyL2L3Residuals == False :
		getattr(process,"patPFMetT1Corr"+postfix).isMC = cms.bool(True) ## avoid to apply residual JEC

	## no L1 correction for Puppi
	if postfix == "Puppi":
		getattr(process,"patPFMetT1Corr"+postfix).offsetCorrLabel = cms.InputTag("") ## info embedded in the jet object 
		  
	## apply type-I corrections
	setattr(process,"patPFMetT1"+postfix,cms.EDProducer("CorrectedPATMETProducer",
							    src = cms.InputTag("patPFMet"+postfix),
							    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1Corr"+postfix,"type1"))))
			

	## re-compute all MET systematics
	if addMETSystematics:
			
		setattr(process,"metSysProducer"+postfix,
			cms.EDProducer("METSystematicsProducer",
				       inputMET    = cms.InputTag("patPFMetT1"+postfix),
				       rho         = cms.InputTag("fixedGridRhoFastjetAll"),
				       pfCandidate = cms.InputTag("packedPFCandidates"),
				       ## skip candidates
				       storeSmearedShiftedCollections = cms.bool(False),
				       skipMuon    = cms.bool(False),
				       skipElectron    = cms.bool(False),
				       skipTau     = cms.bool(False),
				       skipPhoton  = cms.bool(False),
				       skipJet     = cms.bool(False),
				       ## muons
				       muon = cms.PSet(
					src = cms.InputTag("slimmedMuons"),
					useExternalUncertainty = cms.bool(True),
					binning = cms.VPSet(
						cms.PSet(binSelection = cms.string("pt < 100"),
							 uncertainty = cms.double(0.002)),
						cms.PSet(binSelection = cms.string("pt >= 100"),
							 uncertainty = cms.double(0.05))
						)
					),
				       ## electrons
				       electron = cms.PSet(
					src = cms.InputTag("slimmedElectrons"),
					useExternalUncertainty = cms.bool(True),
					binning = cms.VPSet(
						cms.PSet(binSelection = cms.string("isEB"),
							 uncertainty = cms.double(0.006)),
						cms.PSet(binSelection = cms.string("!isEB"),
							 uncertainty = cms.double(0.015))
						)
					),
				       ## electrons
				       photon = cms.PSet(
					src = cms.InputTag("slimmedPhotons"),
					useExternalUncertainty = cms.bool(True),
					binning = cms.VPSet(
						cms.PSet(binSelection = cms.string('isEB'),
							 uncertainty = cms.double(0.01)),
						cms.PSet(binSelection = cms.string('!isEB'),
							 uncertainty = cms.double(0.025))
						)
					),
				       ## taus
				       tau = cms.PSet(
					src = cms.InputTag("slimmedTaus"),
					useExternalUncertainty = cms.bool(True),
					binning = cms.VPSet(
						cms.PSet(binSelection = cms.string("abs(eta) < 2.5 && pt > 18. && tauID(\'decayModeFindingNewDMs\')> 0.5"),
							 uncertainty = cms.double(0.03)),
						)
					),
				       jet = cms.PSet(
					## input collection
					src = cms.InputTag(jetCollection),
					## information for jet energy correction
					payloadName = cms.string(payloadName),
					useExternalJECUncertainty = cms.bool(False),
					JECUncFile = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JEC/Fall15_25nsV2_DATA_Uncertainty_"+payloadName+".txt"),
					#https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
					JERFile    = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JER/Fall15_25nsV2_DATA_PtResolution_"+payloadName+".txt"),
					JERSFFile  = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JER/Fall15_25nsV2_DATA_SF_"+payloadName+".txt"),
					useExternalJERSF = cms.bool(False),
					useExternalJER   = cms.bool(False)),
				       ## unclustered component
				       unclustered = cms.PSet(
					useExternalUncertainty = cms.bool(True),
					binning = cms.VPSet(
						cms.PSet( binSelection = cms.string("abs(eta) < 9.9"),
							  uncertainty = cms.double(0.1)
							  )
						))
				       )
			)
		
		if isMC:
			getattr(process,"metSysProducer"+postfix).jet.JECUncFile = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JEC/Fall15_25nsV2_MC_Uncertainty_"+payloadName+".txt");
			getattr(process,"metSysProducer"+postfix).jet.JERFile = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JER/Fall15_25nsV2_MC_PtResolution_"+payloadName+".txt");
					 
		## final slimmed MET			
		setattr(process,metCollection, cms.EDProducer("PATMETSlimmer",
							      caloMET = cms.InputTag("patPFMet"+postfix),
							      rawVariation = cms.InputTag("patPFMet"+postfix),
							      runningOnMiniAOD = cms.bool(True),
							      src = cms.InputTag("patPFMetT1"+postfix),
							      t01Variation = cms.InputTag(metCollection,"","@skipCurrentProcess"),
							      t1Uncertainties = cms.InputTag("metSysProducer"+postfix,"patPFMetT1"+postfix+"%s"),
							      tXYUncForRaw = cms.InputTag(metCollection,"","@skipCurrentProcess"),
							      tXYUncForT1 = cms.InputTag(metCollection,"","@skipCurrentProcess"),
							      t1SmearedVarsAndUncs = cms.InputTag("metSysProducer"+postfix,"patPFMetT1"+postfix+"Smear%s")
							      ))
	else:
		setattr(process,metCollection, cms.EDProducer("PATMETSlimmer",
							      caloMET = cms.InputTag("patPFMet"+postfix),
							      rawVariation = cms.InputTag("patPFMet"+postfix),
							      runningOnMiniAOD = cms.bool(True),
							      src = cms.InputTag("patPFMetT1"+postfix),
							      t01Variation = cms.InputTag(metCollection,"","@skipCurrentProcess"),
							      tXYUncForRaw = cms.InputTag(metCollection,"","@skipCurrentProcess"),
							      tXYUncForT1 = cms.InputTag(metCollection,"","@skipCurrentProcess"),
							      t1Uncertainties = cms.InputTag(metCollection,"","@skipCurrentProcess"),
							      t1SmearedVarsAndUncs = cms.InputTag(metCollection,"","@skipCurrentProcess")
							      ))
				 
