import os
import FWCore.ParameterSet.Config as cms
from RecoMET.METPUSubtraction.mvaPFMET_cff import calibratedAK4PFJetsForPFMVAMEt, pfMVAMEt

def runMVAMet(process,isMC,leptons):

    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets

    if not hasattr(process,"ak4PFJets"):
        setattr(process,"ak4PFJets", ak4PFJets.clone(
                src = cms.InputTag('packedPFCandidates'))
                )

    if not hasattr(process,"ak4PFL1FastL2L3Corrector"):
        process.load('JetMETCorrections.Configuration.JetCorrectors_cff')


    process.calibratedAK4PFJetsForPFMVAMEt = calibratedAK4PFJetsForPFMVAMEt.clone()
    if not isMC:
        process.calibratedAK4PFJetsForPFMVAMEt.correctors = cms.vstring("ak4PFL1FastL2L3Residual")
    
    from RecoJets.JetProducers.PileupJetID_cfi import pileupJetIdEvaluator
    from RecoJets.JetProducers.PileupJetIDParams_cfi import JetIdParams

    process.puJetIdForPFMVAMEt = pileupJetIdEvaluator.clone(
        produceJetIds = cms.bool(True),
        jets = cms.InputTag("calibratedAK4PFJetsForPFMVAMEt"),
        vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
        jec     = cms.string("AK4PF")
        )

    process.mvaMET = pfMVAMEt.clone(
        srcPFCandidates = cms.InputTag('packedPFCandidates'),
        srcVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
        srcLeptons = cms.VInputTag(leptons),
        inputFileNames = cms.PSet(
            U     = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru_7_4_X_miniAOD_50NS_July2015.root'),
            DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_4_X_miniAOD_50NS_July2015.root'),
            CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_4_X_miniAOD_50NS_July2015.root'),
            CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_4_X_miniAOD_50NS_July2015.root')
            ),
        )
        
                
        
        
        
    
