import os, copy
import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetCorrFactorsUpdated
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
from RecoMET.METProducers.PFMET_cfi import pfMet

## generic function that corrects jet and MET given a JEC
def JetMetCorrector(process, jetCollection, metCollection, payloadName, isMC, applyL2L3Residuals, useMiniAODMet, doMETSystematics):

	## apply corrections on jets
	if not hasattr(process,"patJetCorrFactorsReapplyJEC"+payloadName):
		setattr(process,"patJetCorrFactorsReapplyJEC"+payloadName, patJetCorrFactorsUpdated.clone(
				src     = cms.InputTag(jetCollection),
				levels  = process.JECLevels.labels,
				payload = payloadName
				));

		if "Puppi" in metCollection or "PUPPI" in metCollection:
			puppiJEC = copy.deepcopy(process.JECLevels.labels)
			puppiJEC.remove('L1FastJet')
			getattr(process,"patJetCorrFactorsReapplyJEC"+payloadName).levels = puppiJEC
			getattr(process,"patJetCorrFactorsReapplyJEC"+payloadName).useRho = False
  
	if not hasattr(process,"slimmedJetsRecorrected"+payloadName):
		setattr(process,"slimmedJetsRecorrected"+payloadName,
			patJetsUpdated.clone(jetSource = cms.InputTag(jetCollection),
					     jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"+payloadName))))

	## in case don't want to rely on miniAOD met, propagate corrections and redo-systeamtics
	if not useMiniAODMet: 
		
  	        ## propagation on missing energy + full systematics
		if "Puppi" in metCollection or "PUPPI" in metCollection:
			postfix = "Puppi"
		else:
			postfix = ""
			
	        ## extract genMet		
		if not hasattr(process,"genMetExtractor") and isMC:
			setattr(process,"genMetExtractor",cms.EDProducer("GenMETExtractor",
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
			setattr(process,"patPFMetT1Corr"+postfix, cms.EDProducer("PATPFJetMETcorrInputProducer",
										 isMC = cms.bool(isMC),
										 offsetCorrLabel = cms.InputTag("L1FastJet"),
										 jetCorrLabel = cms.InputTag("L3Absolute"), ## info embedded in the jet object
										 jetCorrLabelRes = cms.InputTag("L2L3Residual"), ## info embedded in the jet object 
										 skipEM = cms.bool(True),
										 skipEMfractionThreshold = cms.double(0.9),
										 skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
										 skipMuons = cms.bool(True),
										 src = cms.InputTag("slimmedJetsRecorrected"+payloadName),
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
			if doMETSystematics:

			      ### global systematics
				 if "Puppi" in payloadName:
					 payloadNameUnc = "AK4PFchs";
				 else:
					 payloadNameUnc = payloadName;
			

				 setattr(process,"metSysProducer"+postfix,
					 cms.EDProducer("METSystematicsProducer",
							rho = cms.InputTag("fixedGridRhoFastjetAll"),
							inputMET = cms.InputTag("patPFMetT1"+postfix),
							pfCandidate = cms.InputTag("packedPFCandidates"),
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
							tau = cms.PSet(
							 src = cms.InputTag("slimmedTaus"),
							 useExternalUncertainty = cms.bool(True),
							 binning = cms.VPSet(
								 cms.PSet(binSelection = cms.string("abs(eta) < 2.5 && pt > 18. && tauID(\'decayModeFindingNewDMs\')> 0.5"),
									  uncertainty = cms.double(0.03)),
								 )
							 ),
							jet = cms.PSet(
							 payloadName = cms.string(payloadNameUnc),
							 src = cms.InputTag(jetCollection),
							 useExternalJECUncertainty = cms.bool(False),
							 JECUncFile = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JEC/Summer15_25nsV6_DATA_Uncertainty_AK4PFchs.txt"),
							 #https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
							 JERFile    = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JEC/Summer15_25nsV6_MC_PtResolution_AK4PFchs.txt"),
							 JERFormula = cms.string("TMath::Sqrt([0]*TMath::Abs([0])/(x*x)+[1]*[1]*TMath::Power(x,[3])+[2]*[2])"),
							 binningJERSF = cms.VPSet(
								 #https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
								 cms.PSet(
									 binSelection = cms.string("abs(eta) < 0.8"),
									 scaleFactor     = cms.double(1.061),
									 scaleFactorUnc  = cms.double(0.023),
									 ),
								 cms.PSet(
									 binSelection = cms.string("abs(eta) >= 0.8 && abs(eta)< 1.3"),
									 scaleFactor     = cms.double(1.088),
									 scaleFactorUnc  = cms.double(0.029),
									 ),
								 cms.PSet(
									 binSelection = cms.string("abs(eta) >= 1.3 && abs(eta) < 1.9"),
									 scaleFactor     = cms.double(1.106),
									 scaleFactorUnc  = cms.double(0.030),
									 ),
								 cms.PSet(
									 binSelection = cms.string("abs(eta) >= 1.9 && abs(eta) < 2.5"),
									 scaleFactor     = cms.double(1.126),
									 scaleFactorUnc  = cms.double(0.094),
									 ),
								 cms.PSet(
									 binSelection = cms.string("abs(eta) >= 2.5 && abs(eta) < 3.0"),
									 scaleFactor     = cms.double(1.343),
									 scaleFactorUnc  = cms.double(0.123),
									 ),
								 cms.PSet(
									 binSelection = cms.string("abs(eta) >= 3.0 && abs(eta) < 3.2"),
									 scaleFactor     = cms.double(1.303),
									 scaleFactorUnc  = cms.double(0.111),
									 ),
								 cms.PSet(
									 binSelection = cms.string("abs(eta) >= 3.2 && abs(eta) < 5"),
									 scaleFactor     = cms.double(1.320),
									 scaleFactorUnc  = cms.double(0.286),
									 )
								 
								 )
							 ),
							unclustered = cms.PSet(
							 useExternalUncertainty = cms.bool(True),
							 binning = cms.VPSet(
								 cms.PSet( binSelection = cms.string("abs(eta) < 9.9"),
									   uncertainty = cms.double(0.1)
									   )
								 )
							 )
					       )
					 )
				 
				 if isMC:
					 getattr(process,"metSysProducer"+postfix).jet.JECUncFile = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JEC/Summer15_25nsV6_MC_Uncertainty_AK4PFchs.txt");

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
				 
			         ## final slimmed MET
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
				 
