import os, copy
import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.metTools import addMETCollection
from RecoMET.METProducers.PFMET_cfi import pfMet
from PhysicsTools.PatUtils.patPFMETCorrections_cff import patPFMetT1T2Corr
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
from JetMETCorrections.Type1MET.pfMETmultShiftCorrections_cfi import pfMEtMultShiftCorr

## generic function that corrects jet and MET given a JEC
def metCorrector(process,jetCollection,metCollection,isMC,payloadName,applyL2L3Residuals,addMETSystematics,useOfficialMETSystematics):
	
	## propagation on missing energy + full systematics
	if "Puppi" in metCollection or "PUPPI" in metCollection:
		postfix = "Puppi"
	else:
		postfix = ""

	if postfix == "Puppi" and not hasattr(process,"puppi"):
		from PhysicsTools.PatAlgos.slimming.puppiForMET_cff import makePuppiesFromMiniAOD
		makePuppiesFromMiniAOD( process, False);
		process.puppi.useExistingWeights = cms.bool(False)
		process.puppiNoLep.useExistingWeights = cms.bool(False)
		
		
	######################
	if useOfficialMETSystematics and addMETSystematics:
		## use the official jet-MET tool
		## re-run for standard met
		if postfix == "Puppi" :
			if isMC:
				runMetCorAndUncFromMiniAOD(process,isData=False,pfCandColl=cms.InputTag("puppiForMET"),metType=postfix,postfix=postfix,jetFlavor="AK4PFPuppi")		 
			else:
				runMetCorAndUncFromMiniAOD(process,isData=True,pfCandColl=cms.InputTag("puppiForMET"),metType=postfix,postfix=postfix,jetFlavor="AK4PFPuppi")		       
	
				
			if applyL2L3Residuals == False and not isMC:
				process.patPFMetT1T2CorrPuppi.jetCorrLabelRes = cms.InputTag("L3Absolute")
				process.patPFMetT1T2SmearCorrPuppi.jetCorrLabelRes = cms.InputTag("L3Absolute")
				process.patPFMetT2CorrPuppi.jetCorrLabelRes = cms.InputTag("L3Absolute")
				process.patPFMetT2SmearCorrPuppi.jetCorrLabelRes = cms.InputTag("L3Absolute")
				process.shiftedPatJetEnDownPuppi.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
				process.shiftedPatJetEnUpPuppi.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
		else:
			if isMC:
				runMetCorAndUncFromMiniAOD(process,isData=False)
			else:
				runMetCorAndUncFromMiniAOD(process,isData=True)
				
			if applyL2L3Residuals == False and not isMC:
				process.patPFMetT1T2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
				process.patPFMetT1T2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
				process.patPFMetT2Corr.jetCorrLabelRes = cms.InputTag("L3Absolute")
				process.patPFMetT2SmearCorr.jetCorrLabelRes = cms.InputTag("L3Absolute")
				process.shiftedPatJetEnDown.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
				process.shiftedPatJetEnUp.jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector")
				

        else:		
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
		if postfix == "Puppi":
			getattr(process,"pfMet"+postfix).src = cms.InputTag('puppiForMET')

  	         ## re-cast PFMets into PAT objects
		addMETCollection(process, labelName='patPFMet'+postfix, metSource='pfMet'+postfix)
		
		if isMC:
			getattr(process,"patPFMet"+postfix).addGenMET = cms.bool(True)
			getattr(process,"patPFMet"+postfix).genMETSource = cms.InputTag("genMetExtractor")
		else:
			getattr(process,"patPFMet"+postfix).addGenMET = cms.bool(False)

		 
 	        ## derive type-I corrector object
		setattr(process,"patPFMetT1Corr"+postfix, 
			patPFMetT1T2Corr.clone(
				src = cms.InputTag(jetCollection),
				type1JetPtThreshold = cms.double(15.0),
				))
	


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
					       storeSmearedShiftedCollections = cms.bool(True),
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
						selection = cms.string('pt > 15 && abs(eta) < 9.9 && (chargedEmEnergyFraction+neutralEmEnergyFraction) < 0.9'),
						## information for jet energy correction
						payloadName = cms.string(payloadName),
						useExternalJECUncertainty = cms.bool(False),
						JECUncFile = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JEC/Fall15_25nsV2_DATA_Uncertainty_"+payloadName+".txt"),
						#https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
						JERFile      = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JER/Fall15_25nsV2_DATA_PtResolution_"+payloadName+".txt"),
						JERSFFile    = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JER/Fall15_25nsV2_DATA_SF_"+payloadName+".txt"),
						jetCorrLabel    = cms.InputTag("L3Absolute"),
						jetCorrLabelRes = cms.InputTag("L2L3Residual"), 
						useExternalJERSF = cms.bool(False),
						useExternalJER   = cms.bool(False)),
					       ## unclustered component
					       unclustered = cms.PSet(
						useExternalUncertainty = cms.bool(True),
						binning = cms.VPSet(
							cms.PSet( binSelection   = cms.string("charge!=0"),
								  binUncertainty = cms.string("sqrt(pow(0.00009*pt,2)+pow(0.0085/sqrt(sin(2*atan(exp(eta)))),2))")
								  ),
							cms.PSet( binSelection   = cms.string("pdgId==130"),
								  binUncertainty = cms.string("((abs(eta)<1.3)?(max(0.25,sqrt(pow(0.8/sqrt(energy), 2)+0.05*0.05))):(max(0.30,sqrt(pow(1.0/sqrt(energy),2)+0.04*0.04))))")
								  ),
							cms.PSet( binSelection   = cms.string("pdgId==22"),
								  binUncertainty = cms.string("sqrt(pow(0.00009*energy,2)+pow(0.0085/sqrt(sin(2*atan(exp(eta)))),2))")
                                                                  ),
							
							cms.PSet(
								binSelection = cms.string('pdgId==1 || pdgId==2'),
								binUncertainty = cms.string('sqrt(pow(1/sqrt(energy),2)+0.05*0.05)+0*eta'),
								)),
						
						))
				)
			
		
			if isMC:
				getattr(process,"metSysProducer"+postfix).jet.JECUncFile = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JEC/Fall15_25nsV2_MC_Uncertainty_"+payloadName+".txt");
				getattr(process,"metSysProducer"+postfix).jet.JERFile = cms.FileInPath("AnalysisCode/MonoXAnalysis/data/JER/Fall15_25nsV2_MC_PtResolution_"+payloadName+".txt");
			
			### add xy corrections
			if not hasattr(process,"patPFMetTxyCorr"+postfix):
				setattr(process,"patPFMetTxyCorr"+postfix,pfMEtMultShiftCorr.clone(
						srcPFlow = cms.InputTag('packedPFCandidates'),
						vertexCollection = cms.InputTag('offlineSlimmedPrimaryVertices')
						))
				
				setattr(process,'patPFMetT1Txy'+postfix,
					cms.EDProducer("CorrectedPATMETProducer",
						       src = cms.InputTag("patPFMet"+postfix),
						       srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1Corr"+postfix,"type1"), cms.InputTag("patPFMetTxyCorr"+postfix))))
						       
					 
			## final slimmed MET			
		       	setattr(process,metCollection, cms.EDProducer("PATMETSlimmer",
								      caloMET = cms.InputTag("patPFMet"+postfix),
								      rawVariation = cms.InputTag("patPFMet"+postfix),
								      runningOnMiniAOD = cms.bool(True),
								      src = cms.InputTag("patPFMetT1"+postfix),
								      t01Variation = cms.InputTag(metCollection,"","@skipCurrentProcess"),
								      t1Uncertainties = cms.InputTag("metSysProducer"+postfix,"patPFMetT1"+postfix+"%s"),
								      tXYUncForRaw = cms.InputTag(metCollection,"","@skipCurrentProcess"),
								      tXYUncForT1 = cms.InputTag('patPFMetT1Txy'+postfix),
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
				 
