import os
import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
import random

def ElectronTools(process,addEGMSmear,isMC):


	# Electron ValueMaps for identification --> takes into account both id+iso
	dataFormat = DataFormat.MiniAOD
	switchOnVIDElectronIdProducer(process, dataFormat);
	ele_id_modules = [];
	ele_id_modules.append('RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff');
	#ele_id_modules.append('RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff');
	ele_id_modules.append('RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff');
	ele_id_modules.append('RecoEgamma.ElectronIdentification.Identification.cutBasedElectronHLTPreselecition_Summer16_V1_cff');

	for idmod in ele_id_modules:
	   	setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

		
	# add scale + smearing energy corrections
	if addEGMSmear:

		process.selectedElectrons = cms.EDFilter("PATElectronSelector", 
							 src = cms.InputTag("slimmedElectrons"), 
							 cut = cms.string("pt > 5 && abs(eta) < 3.0")) 
		
		from EgammaAnalysis.ElectronTools.calibratedElectronsRun2_cfi import calibratedPatElectrons,files
		setattr(process,"calibratedElectrons",calibratedPatElectrons.clone(
				isMC = cms.bool(isMC),				
				electrons = cms.InputTag('selectedElectrons'),
				correctionFile = cms.string(files["80Xapproval"])
				))

		
	
		process.RandomNumberGeneratorService.calibratedElectrons = cms.PSet(
			initialSeed = cms.untracked.uint32(int(random.uniform(0,1000000))),
			engineName = cms.untracked.string('TRandom3')
			)
		
		
	#### apply gain corrections
	setattr(process,"correctedElectrons",cms.EDProducer("PATElectronCorrector",
							    src = cms.InputTag("slimmedElectrons"),
							    isMC = cms.bool(isMC),
							    correction = cms.VPSet(
				cms.PSet(
					eMin = cms.double(200),
					eMax = cms.double(300),
					value = cms.double(1.0199)),
				cms.PSet(
					eMin = cms.double(300),
					eMax = cms.double(400),
					value = cms.double(1.0520)),
				cms.PSet(
					eMin = cms.double(400),
					eMax = cms.double(500),
					value = cms.double(1.0150)),
				cms.PSet(
					eMin = cms.double(500),
					eMax = cms.double(10000),
					value = cms.double(1.0150))),
							    recHitEB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
							    recHitEE = cms.InputTag("reducedEgamma","reducedEERecHits")
							    ));
	if addEGMSmear:
		getattr(process,"correctedElectrons").src = cms.InputTag("calibratedElectrons");
						    
			
			
