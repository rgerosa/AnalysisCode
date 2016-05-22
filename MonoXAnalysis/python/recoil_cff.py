import os
import FWCore.ParameterSet.Config as cms

def recoilComputation(process,processName,miniAODProcess,useMiniAODMet,isPuppi):

    if not isPuppi:
        metName = "slimmedMETs"
    else:
        metName = "slimmedMETsPuppi"

    process.t1mumet = cms.EDProducer("MuonCorrectedMETProducer",
                                     met     = cms.InputTag(metName,"",processName),
                                     cands   = cms.VInputTag(cms.InputTag("selectedObjects", "muons")),
                                     isPuppi = cms.bool(isPuppi),
                                     pfCandidates = cms.InputTag("packedPFCandidates"))
    if useMiniAODMet:
        process.t1mumet.met = cms.InputTag(metName,"",miniAODProcess)
        
    process.t1elmet = cms.EDProducer("ElectronCorrectedMETProducer",
                                     met     = cms.InputTag(metName,"",processName),
                                     cands   = cms.VInputTag(cms.InputTag("selectedObjects", "electrons")),
                                     isPuppi = cms.bool(isPuppi),
                                     pfCandidates = cms.InputTag("packedPFCandidates"))
    
    if useMiniAODMet:
        process.t1elmet.met = cms.InputTag(metName,"",miniAODProcess)

    process.t1phmet = cms.EDProducer("PhotonCorrectedMETProducer",
                                     met     = cms.InputTag(metName,"",processName),
                                     cands   = cms.VInputTag(cms.InputTag("selectedObjects", "photons")),
                                     isPuppi = cms.bool(isPuppi),
                                     pfCandidates = cms.InputTag("packedPFCandidates"))
    
    if useMiniAODMet:
        process.t1phmet.met = cms.InputTag(metName,"",miniAODProcess)

    process.t1taumet = cms.EDProducer("TauCorrectedMETProducer",
                                      met     = cms.InputTag(metName,"",processName),
                                      cands   = cms.VInputTag(cms.InputTag("selectedObjects", "taus")),
                                      isPuppi = cms.bool(isPuppi),
                                      pfCandidates = cms.InputTag("packedPFCandidates"))
    
    if useMiniAODMet:
        process.t1phmet.met = cms.InputTag(metName,"",miniAODProcess)

