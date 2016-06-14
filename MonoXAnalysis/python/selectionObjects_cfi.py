import os
import FWCore.ParameterSet.Config as cms

selectedObjects = cms.EDProducer("PFCleaner",
                                 ### vertex collection
                                 vertices  = cms.InputTag("goodVertices"),
                                 ### pf candidates
                                 pfcands   = cms.InputTag("packedPFCandidates"),
                                 #### rho
                                 rho       = cms.InputTag("fixedGridRhoFastjetAll"),
                                 ### jets 
                                 jets      = cms.InputTag("slimmedAK4Jets"),
                                 ### muon information
                                 muons     = cms.InputTag("slimmedMuons"),
                                 muonSelection = cms.VPSet(
        cms.PSet(
            idType = cms.string("loose"),
            muonCollectionName = cms.string("muons"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.4),
            deltaBeta = cms.double(0.5),
            isolation = cms.double(0.25)),
        cms.PSet(
            idType = cms.string("tight"),
            muonCollectionName = cms.string("tightmuons"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.4),
            deltaBeta = cms.double(0.5),
            isolation = cms.double(0.15)),
        cms.PSet(
            idType = cms.string("highPt"),
            muonCollectionName = cms.string("highptmuons"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.4),
            deltaBeta = cms.double(0.5),
            isolation = cms.double(0.15))),
                                 ### electrons
                                 electrons = cms.InputTag("slimmedElectrons"),
                                 calibratedElectrons = cms.InputTag("calibratedElectrons"),
                                 useCalibratedElectrons = cms.bool(False),
                                 electronSelection = cms.VPSet(
        cms.PSet(
            electronCollectionName = cms.string("electrons"),
            idType = cms.string("veto"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.5),
            eleValueMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto")),
        cms.PSet(
            electronCollectionName = cms.string("looseelectrons"),
            idType = cms.string("loose"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.5),
            eleValueMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose")),
        cms.PSet(
            electronCollectionName = cms.string("tightelectrons"),
            idType = cms.string("tight"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.5),
            eleValueMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight")),
        cms.PSet(
            electronCollectionName = cms.string("heepelectrons"),
            idType = cms.string("heep"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.5),
            eleValueMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"))),
                                 ### taus
                                 taus      = cms.InputTag("slimmedTaus"),
                                 tauSelection = cms.VPSet(
        cms.PSet(
            tauCollectionName = cms.string("taus"),
            dRCleaning = cms.double(0.4),
            tauIDName  = cms.string("byLooseCombinedIsolationDeltaBetaCorr3Hits"),
            ptMin  = cms.double(18),
            absEta = cms.double(2.3),
            decayModeFinding = cms.double(0.5),
            isolation = cms.double(0.5),
            graterThan = cms.bool(True)),
        cms.PSet(
            tauCollectionName = cms.string("tausRaw"),
            tauIDName  = cms.string("byCombinedIsolationDeltaBetaCorrRaw3Hits"),
            dRCleaning = cms.double(0.4),
            ptMin  = cms.double(18),
            absEta = cms.double(2.3),
            decayModeFinding = cms.double(0.5),
            isolation = cms.double(5),
            graterThan = cms.bool(False))
        ),
                                 #### photons
                                 photons   = cms.InputTag("slimmedPhotons"),
                                 calibratedPhotons = cms.InputTag("calibratedPhotons"),
                                 useCalibratedPhotons = cms.bool(False),
                                 addPhotonPurity = cms.bool(False),
                                 photonSelection = cms.VPSet(
        cms.PSet(
            photonCollectionName = cms.string("photons"),
            ptMin  = cms.double(15),
            absEta = cms.double(2.5),
            photonValueMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-loose")),
        cms.PSet(
            photonCollectionName = cms.string("mediumphotons"),
            ptMin  = cms.double(15),
            absEta = cms.double(2.5),
            photonValueMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-medium")),        
        cms.PSet(
            photonCollectionName = cms.string("tightphotons"),
            ptMin  = cms.double(15),
            absEta = cms.double(2.5),
            photonValueMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-25ns-V1-standalone-tight"))),        
                                 ### cluster shape and isolation value maps
                                 photonsieie = cms.InputTag("photonIDValueMapProducer", "phoFull5x5SigmaIEtaIEta"),
                                 photonphiso = cms.InputTag("photonIDValueMapProducer", "phoPhotonIsolation"),
                                 photonchiso = cms.InputTag("photonIDValueMapProducer", "phoChargedIsolation"),
                                 ### my loose photon id
                                 loosePhotonID = cms.PSet(
        ptMin  = cms.double(15),
        absEta = cms.double(2.5),
        R9min = cms.double(0.8),
        chIso = cms.double(20),
        chIsoFrac = cms.double(0.3)),
                                 ### high pt photon id: AN (2016/079)
                                 highPtPhotonID = cms.PSet(
        absEta = cms.double(1.4442),
        ptMin  = cms.double(15),
        chIso  = cms.double(5),
        sigmaIetaIeta = cms.double(0.0105),
        HOverE = cms.double(0.05),
        isolation = cms.double(2.75),
        alpha  = cms.double(2.5),
        k      = cms.double(0.0045))                                                                  
)
