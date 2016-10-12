import os
import FWCore.ParameterSet.Config as cms

def recoilComputation(process,processName,miniAODProcess,useMiniAODMet,isPuppi):

    if not isPuppi:
        metName = "slimmedMETs"
    else:
        metName = "slimmedMETsPuppi"

    if(isPuppi):
        postfix = "puppi"
    else:
        postfix = ""
    
    if not hasattr(process,postfix+"t1mumet"):
        setattr(process, postfix+"t1mumet",cms.EDProducer("MuonCorrectedMETProducer",
                                                          met     = cms.InputTag(metName,"",processName),
                                                          cands   = cms.VInputTag(cms.InputTag("selectedObjects", "muons")),
                                                          isPuppi = cms.bool(isPuppi),
                                                          pfCandidates = cms.InputTag("packedPFCandidates")))
        if useMiniAODMet:
            getattr(process,postfix+"t1mumet").met = cms.InputTag(metName,"",miniAODProcess)
        
    if not hasattr(process,postfix+"t1elmet"):
        setattr(process, postfix+"t1elmet",cms.EDProducer("ElectronCorrectedMETProducer",
                                                          met     = cms.InputTag(metName,"",processName),
                                                          cands   = cms.VInputTag(cms.InputTag("selectedObjects", "electrons")),
                                                          isPuppi = cms.bool(isPuppi),
                                                          pfCandidates = cms.InputTag("packedPFCandidates")))
                
        if useMiniAODMet:
           getattr(process,postfix+"t1elmet").met = cms.InputTag(metName,"",miniAODProcess)

    if not hasattr(process,postfix+"t1phmet"):
        setattr(process, postfix+"t1phmet",cms.EDProducer("PhotonCorrectedMETProducer",
                                                          met     = cms.InputTag(metName,"",processName),
                                                          cands   = cms.VInputTag(cms.InputTag("selectedObjects", "photons")),
                                                          isPuppi = cms.bool(isPuppi),
                                                          pfCandidates = cms.InputTag("packedPFCandidates")))
    
        if useMiniAODMet:
            getattr(process,postfix+"t1phmet").met  = cms.InputTag(metName,"",miniAODProcess)

    if not hasattr(process,postfix+"t1taumet"):
        setattr(process, postfix+"t1taumet",cms.EDProducer("TauCorrectedMETProducer",
                                                           met     = cms.InputTag(metName,"",processName),
                                                           cands   = cms.VInputTag(cms.InputTag("selectedObjects", "tausNew")),
                                                           isPuppi = cms.bool(isPuppi),
                                                           pfCandidates = cms.InputTag("packedPFCandidates")))
        
        if useMiniAODMet:
           getattr(process,postfix+"t1taumet").met = cms.InputTag(metName,"",miniAODProcess)

