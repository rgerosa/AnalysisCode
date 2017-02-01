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
        #### loose muons
        cms.PSet(
            idType = cms.string("loose"),
            muonCollectionName = cms.string("muons"),
            ptMin     = cms.double(10),
            absEta    = cms.double(2.4),
            deltaBeta = cms.double(0.5),
            isolation = cms.double(0.25)),
        #### tight muons
        cms.PSet(
            idType = cms.string("tight"),
            muonCollectionName = cms.string("tightmuons"),
            ptMin     = cms.double(10),
            absEta    = cms.double(2.4),
            deltaBeta = cms.double(0.5),
            isolation = cms.double(0.15)),
        #### high pt muons
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
                                 useCalibratedElectrons = cms.bool(True),
                                 electronSelection = cms.VPSet(
        #### veto electrons
        cms.PSet(
            electronCollectionName = cms.string("electrons"),
            idType = cms.string("veto"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.5),
            applyPVSelection = cms.bool(True),
            PVSelection = cms.PSet(   
                d0Barrel = cms.double(0.05),
                d0Endcap = cms.double(0.10),
                dzBarrel = cms.double(0.10),
                dzEndcap = cms.double(0.20),
                ),            
            eleValueMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto")),
            #### loose electrons
        cms.PSet(
            electronCollectionName = cms.string("looseelectrons"),
            idType = cms.string("loose"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.5),
            applyPVSelection = cms.bool(True),
            PVSelection = cms.PSet(   
                d0Barrel = cms.double(0.05),
                d0Endcap = cms.double(0.10),
                dzBarrel = cms.double(0.10),
                dzEndcap = cms.double(0.20),
                ),            
            eleValueMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose")),
        #### tight electrons
        cms.PSet(
            electronCollectionName = cms.string("tightelectrons"),
            idType = cms.string("tight"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.5),
            applyPVSelection = cms.bool(True),
            PVSelection = cms.PSet(   
                d0Barrel = cms.double(0.05),
                d0Endcap = cms.double(0.10),
                dzBarrel = cms.double(0.10),
                dzEndcap = cms.double(0.20),
                ),            
            eleValueMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight")),
        ##### trigger hlt safe
        cms.PSet(
            electronCollectionName = cms.string("triggerelectrons"),
            idType = cms.string("trigger"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.5),
            applyPVSelection = cms.bool(False),
            PVSelection = cms.PSet(   
                d0Barrel = cms.double(0.05),
                d0Endcap = cms.double(0.10),
                dzBarrel = cms.double(0.10),
                dzEndcap = cms.double(0.20),
                ),            
            eleValueMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronHLTPreselection-Summer16-V1")),
        ####### heep id
        cms.PSet(
            electronCollectionName = cms.string("heepelectrons"),
            idType = cms.string("heep"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.5),
            applyPVSelection = cms.bool(False),
            PVSelection = cms.PSet(   
                d0Barrel = cms.double(0.05),
                d0Endcap = cms.double(0.10),
                dzBarrel = cms.double(0.10),
                dzEndcap = cms.double(0.20),
                ),            
            eleValueMap = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60")),
        ##### mvaloose
        cms.PSet(
            electronCollectionName = cms.string("mvalooseelectrons"),
            idType = cms.string("mvaloose"),
            ptMin  = cms.double(10),            
            absEta = cms.double(2.5),
            applyPVSelection = cms.bool(False),
            PVSelection = cms.PSet(
                d0Barrel = cms.double(0.05),
                d0Endcap = cms.double(0.10),
                dzBarrel = cms.double(0.10),
                dzEndcap = cms.double(0.20),
                ),
            eleValueMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90")),
        ##### mvatight
        cms.PSet(
            electronCollectionName = cms.string("mvatightelectrons"),
            idType = cms.string("mvatight"),
            ptMin  = cms.double(10),
            absEta = cms.double(2.5),
            applyPVSelection = cms.bool(False),
            PVSelection = cms.PSet(
                d0Barrel = cms.double(0.05),
                d0Endcap = cms.double(0.10),
                dzBarrel = cms.double(0.10),
                dzEndcap = cms.double(0.20),
                ),
            eleValueMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80"))),
                                 ### taus
                                 taus      = cms.InputTag("slimmedTaus"),
                                 tauSelection = cms.VPSet(
        cms.PSet(
            #### very loose MVA with new decay mode
            tauCollectionName = cms.string("tausVLNew"),
            dRCleaning = cms.double(0.4),
            tauIDName  = cms.string("byVLooseIsolationMVArun2v1DBnewDMwLT"),
            ptMin  = cms.double(18),
            absEta = cms.double(2.3),
            decayModeFinding = cms.double(0.5),
            isolation        = cms.double(0.5),
            useNewDecayMode  = cms.bool(True),
            graterThan       = cms.bool(True)),
            #### very loose MVA with old decay mode
        cms.PSet(
            tauCollectionName = cms.string("tausVLOld"),
            dRCleaning = cms.double(0.4),
            tauIDName  = cms.string("byVLooseIsolationMVArun2v1DBoldDMwLT"),
            ptMin  = cms.double(18),
            absEta = cms.double(2.3),
            decayModeFinding  = cms.double(0.5),
            isolation         = cms.double(0.5),
            useNewDecayMode   = cms.bool(False),
            graterThan        = cms.bool(True)),
        cms.PSet(
            #### raw isolation with new decay mode
            tauCollectionName = cms.string("tausRawNew"),
            tauIDName  = cms.string("byCombinedIsolationDeltaBetaCorrRaw3Hits"),
            dRCleaning = cms.double(0.4),
            ptMin  = cms.double(18),
            absEta = cms.double(2.3),
            decayModeFinding = cms.double(0.5),
            isolation        = cms.double(5),
            useNewDecayMode  = cms.bool(True),
            graterThan       = cms.bool(False)),            
        #### raw isolation with old decay mode
        cms.PSet(
            tauCollectionName = cms.string("tausRawOld"),
            tauIDName  = cms.string("byCombinedIsolationDeltaBetaCorrRaw3Hits"),
            dRCleaning = cms.double(0.4),
            ptMin  = cms.double(18),
            absEta = cms.double(2.3),
            decayModeFinding = cms.double(0.5),
            isolation        = cms.double(5),
            useNewDecayMode  = cms.bool(False),
            graterThan       = cms.bool(False)),
        #### tight tau-id with new decay mode
        cms.PSet(
            tauCollectionName = cms.string("tausTightNew"),
            tauIDName  = cms.string("byTightIsolationMVArun2v1DBnewDMwLT"),
            dRCleaning = cms.double(0.4),
            ptMin  = cms.double(18),
            absEta = cms.double(2.3),
            decayModeFinding = cms.double(0.5),
            isolation        = cms.double(0.5),
            useNewDecayMode  = cms.bool(True),
            graterThan       = cms.bool(True)),
        #### tight tau-id with old decay mode
        cms.PSet(
            tauCollectionName = cms.string("tausTightOld"),
            tauIDName  = cms.string("byTightIsolationMVArun2v1DBoldDMwLT"),
            dRCleaning = cms.double(0.4),
            ptMin  = cms.double(18),
            absEta = cms.double(2.3),
            decayModeFinding = cms.double(0.5),
            isolation        = cms.double(0.5),
            useNewDecayMode  = cms.bool(False),
            graterThan       = cms.bool(True)),
        ),
                                 #### photons
                                 photons   = cms.InputTag("slimmedPhotons"),
                                 calibratedPhotons = cms.InputTag("calibratedPhotons"),
                                 useCalibratedPhotons = cms.bool(True),
                                 addPhotonPurity = cms.bool(False),
                                 photonSelection = cms.VPSet(
        #### loose photons
        cms.PSet(
            photonCollectionName = cms.string("photons"),
            ptMin  = cms.double(15),
            absEta = cms.double(2.5),
            photonValueMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-loose")),
        ##### medium photons
        cms.PSet(
            photonCollectionName = cms.string("mediumphotons"),
            ptMin  = cms.double(15),
            absEta = cms.double(2.5),
            photonValueMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium")),        
        ##### tight phtotons
        cms.PSet(
            photonCollectionName = cms.string("tightphotons"),
            ptMin  = cms.double(15),
            absEta = cms.double(2.5),
            photonValueMap = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-tight")),
        cms.PSet(
            photonCollectionName = cms.string("mvaloosephotons"),
            ptMin  = cms.double(15),
            absEta = cms.double(2.5),
            photonValueMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp90")),            
        cms.PSet(
            photonCollectionName = cms.string("mvatightphotons"),
            ptMin  = cms.double(15),
            absEta = cms.double(2.5),
            photonValueMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring16-nonTrig-V1-wp80")),            
        ),
                                 ### cluster shape and isolation value maps
                                 photonsieie = cms.InputTag("photonIDValueMapProducer", "phoFull5x5SigmaIEtaIEta"),
                                 photonphiso = cms.InputTag("photonIDValueMapProducer", "phoPhotonIsolation"),
                                 photonchiso = cms.InputTag("photonIDValueMapProducer", "phoChargedIsolation"),
                                 photonnhiso = cms.InputTag("photonIDValueMapProducer", "phoNeutralHadronIsolation"),
                                 ### my loose photon id
                                 photonPurityID = cms.PSet(
                                     ptMin  = cms.double(150),
                                     etaMax = cms.double(2.4),
                                     HOverE = cms.double(0.07),
                                     chIso  = cms.double(3.0), ### using effective area correction
                                     nhIso  = cms.double(20.),
                                     sigmaIetaIeta  = cms.double(0.05))
                                 )
