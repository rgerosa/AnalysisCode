import os
import FWCore.ParameterSet.Config as cms
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
import random

def PhotonTools(process,addEGMSmear,isMC):

	# Photon ValueMaps for identification
	dataFormat = DataFormat.MiniAOD;
	switchOnVIDPhotonIdProducer(process, dataFormat);
	ph_id_modules = []
	ph_id_modules.append('RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff')
		
	for idmod in ph_id_modules:
		setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)


	if addEGMSmear:
                from EgammaAnalysis.ElectronTools.calibratedPhotonsRun2_cfi import calibratedPatPhotons,files
                setattr(process,"calibratedPhotons",calibratedPatPhotons.clone(
                                isMC = cms.bool(isMC),
				correctionFile = cms.string(files["76XReReco"]),
				
                                ))

		process.RandomNumberGeneratorService.calibratedPhotons = cms.PSet(
			initialSeed = cms.untracked.uint32(int(random.uniform(0,100000))),
			engineName = cms.untracked.string('TRandom3')
			)
