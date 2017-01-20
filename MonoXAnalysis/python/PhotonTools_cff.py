import os
import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
import random

def PhotonTools(process,addEGMSmear,isMC):

	# Photon ValueMaps for identification
	dataFormat = DataFormat.MiniAOD;
	switchOnVIDPhotonIdProducer(process, dataFormat);
	ph_id_modules = []
	ph_id_modules.append('RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff')
	ph_id_modules.append('RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring16_nonTrig_V1_cff')
		
	for idmod in ph_id_modules:
		setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)


	if addEGMSmear:
                from EgammaAnalysis.ElectronTools.calibratedPhotonsRun2_cfi import calibratedPatPhotons,files
                setattr(process,"calibratedPhotons",calibratedPatPhotons.clone(
                                isMC = cms.bool(isMC),
				correctionFile = cms.string(files["80Xapproval"]),
				
                                ))

		process.RandomNumberGeneratorService.calibratedPhotons = cms.PSet(
			initialSeed = cms.untracked.uint32(int(random.uniform(0,100000))),
			engineName = cms.untracked.string('TRandom3')
			)

	#### apply gain corrections                                                                                                                                                                   
	setattr(process,"correctedPhotons",cms.EDProducer("PATPhotonCorrector",
							  src = cms.InputTag("slimmedPhotons"),
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
							  ))
	if addEGMSmear:
		getattr(process,"correctedPhotons").src = cms.InputTag("calibratedPhotons");
