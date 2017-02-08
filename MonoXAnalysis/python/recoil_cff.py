import os
import FWCore.ParameterSet.Config as cms

def recoilComputation(process,processName,miniAODProcess,useMiniAODMet,isPuppi,isReMiniAOD):
    
    if not isPuppi:
        metName = "slimmedMETs";
        postfix = "";
    else:
        metName = "slimmedMETsPuppi";
        postfix = "puppi";

    ### standard applications for MET
    if not isReMiniAOD:

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
                                                               cands   = cms.VInputTag(cms.InputTag("selectedObjects", "tausVLNew")),
                                                               isPuppi = cms.bool(isPuppi),
                                                               pfCandidates = cms.InputTag("packedPFCandidates")))
        
        if useMiniAODMet:
            getattr(process,postfix+"t1taumet").met = cms.InputTag(metName,"",miniAODProcess)

    else: ### only for new re-miniAOD in data

        if not isPuppi:
            if not hasattr(process,postfix+"t1mumet"):
                setattr(process, postfix+"t1mumet",cms.EDProducer("MuonCorrectedMETProducer",
                                                                  met     = cms.InputTag(metName+"MuEGClean","",miniAODProcess),
                                                                  cands   = cms.VInputTag(cms.InputTag("selectedObjects", "muons")),
                                                                  isPuppi = cms.bool(isPuppi),
                                                                  pfCandidates = cms.InputTag("packedPFCandidates")))
            if not hasattr(process,postfix+"t1mumetEGClean"):
                setattr(process, postfix+"t1mumetEGClean",cms.EDProducer("MuonCorrectedMETProducer",
                                                                         met     = cms.InputTag(metName+"EGClean","",miniAODProcess),
                                                                         cands   = cms.VInputTag(cms.InputTag("selectedObjects", "muons")),
                                                                         isPuppi = cms.bool(isPuppi),
                                                                         pfCandidates = cms.InputTag("packedPFCandidates")))                
            if not hasattr(process,postfix+"t1mumetMuClean"):
                setattr(process, postfix+"t1mumetMuClean",cms.EDProducer("MuonCorrectedMETProducer",
                                                                         met     = cms.InputTag(metName,"",miniAODProcess),
                                                                         cands   = cms.VInputTag(cms.InputTag("selectedObjects", "muons")),
                                                                         isPuppi = cms.bool(isPuppi),
                                                                         pfCandidates = cms.InputTag("packedPFCandidates")))

            if not hasattr(process,postfix+"t1elmet"):
                setattr(process, postfix+"t1elmet",cms.EDProducer("MuonCorrectedMETProducer",
                                                                  met     = cms.InputTag(metName+"MuEGClean","",miniAODProcess),
                                                                  cands   = cms.VInputTag(cms.InputTag("selectedObjects", "electrons")),
                                                                  isPuppi = cms.bool(isPuppi),
                                                                  pfCandidates = cms.InputTag("packedPFCandidates")))

            if not hasattr(process,postfix+"t1elmetEGClean"):
                setattr(process, postfix+"t1elmetEGClean",cms.EDProducer("MuonCorrectedMETProducer",
                                                                         met     = cms.InputTag(metName+"EGClean","",miniAODProcess),
                                                                         cands   = cms.VInputTag(cms.InputTag("selectedObjects", "electrons")),
                                                                         isPuppi = cms.bool(isPuppi),
                                                                         pfCandidates = cms.InputTag("packedPFCandidates")))
                
            if not hasattr(process,postfix+"t1elmetMuClean"):
                setattr(process, postfix+"t1elmetMuClean",cms.EDProducer("MuonCorrectedMETProducer",
                                                                         met     = cms.InputTag(metName,"",miniAODProcess),
                                                                         cands   = cms.VInputTag(cms.InputTag("selectedObjects", "electrons")),
                                                                         isPuppi = cms.bool(isPuppi),
                                                                         pfCandidates = cms.InputTag("packedPFCandidates")))

            
            if not hasattr(process,postfix+"t1phmet"):
                setattr(process, postfix+"t1phmet",cms.EDProducer("MuonCorrectedMETProducer",
                                                                  met     = cms.InputTag(metName+"MuEGClean","",miniAODProcess),
                                                                  cands   = cms.VInputTag(cms.InputTag("selectedObjects", "photons")),
                                                                  isPuppi = cms.bool(isPuppi),
                                                                  pfCandidates = cms.InputTag("packedPFCandidates")))


            if not hasattr(process,postfix+"t1phmetEGClean"):
                setattr(process, postfix+"t1phmetEGClean",cms.EDProducer("MuonCorrectedMETProducer",
                                                                         met     = cms.InputTag(metName+"EGClean","",miniAODProcess),
                                                                         cands   = cms.VInputTag(cms.InputTag("selectedObjects", "photons")),
                                                                         isPuppi = cms.bool(isPuppi),
                                                                         pfCandidates = cms.InputTag("packedPFCandidates")))

            if not hasattr(process,postfix+"t1phmetMuClean"):
                setattr(process, postfix+"t1phmetMuEGClean",cms.EDProducer("MuonCorrectedMETProducer",
                                                                           met     = cms.InputTag(metName,"",miniAODProcess),
                                                                           cands   = cms.VInputTag(cms.InputTag("selectedObjects", "photons")),
                                                                           isPuppi = cms.bool(isPuppi),
                                                                           pfCandidates = cms.InputTag("packedPFCandidates")))

                
            if not hasattr(process,postfix+"t1taumet"):
                setattr(process, postfix+"t1taumet",cms.EDProducer("MuonCorrectedMETProducer",
                                                                   met     = cms.InputTag(metName+"MuEGClean","",miniAODProcess),
                                                                   cands   = cms.VInputTag(cms.InputTag("selectedObjects", "tausVLNew")),
                                                                   isPuppi = cms.bool(isPuppi),
                                                                   pfCandidates = cms.InputTag("packedPFCandidates")))


            if not hasattr(process,postfix+"t1taumetEGClean"):
                setattr(process, postfix+"t1taumetEGClean",cms.EDProducer("MuonCorrectedMETProducer",
                                                                          met     = cms.InputTag(metName+"EGClean","",miniAODProcess),
                                                                          cands   = cms.VInputTag(cms.InputTag("selectedObjects", "tausVLNew")),
                                                                          isPuppi = cms.bool(isPuppi),
                                                                          pfCandidates = cms.InputTag("packedPFCandidates")))
                
            if not hasattr(process,postfix+"t1taumetMuClean"):
                setattr(process, postfix+"t1taumetMuClean",cms.EDProducer("MuonCorrectedMETProducer",
                                                                          met     = cms.InputTag(metName,"",miniAODProcess),
                                                                          cands   = cms.VInputTag(cms.InputTag("selectedObjects", "tausVLNew")),
                                                                          isPuppi = cms.bool(isPuppi),
                                                                          pfCandidates = cms.InputTag("packedPFCandidates")))
        else:

            if not hasattr(process,postfix+"t1mumet"):
                setattr(process, postfix+"t1mumet",cms.EDProducer("MuonCorrectedMETProducer",
                                                                  met     = cms.InputTag(metName,"",miniAODProcess),
                                                                  cands   = cms.VInputTag(cms.InputTag("selectedObjects", "muons")),
                                                                  isPuppi = cms.bool(isPuppi),
                                                                  pfCandidates = cms.InputTag("packedPFCandidates")))
                
            if not hasattr(process,postfix+"t1elmet"):
                setattr(process, postfix+"t1elmet",cms.EDProducer("ElectronCorrectedMETProducer",
                                                                  met     = cms.InputTag(metName,"",miniAODProcess),
                                                                  cands   = cms.VInputTag(cms.InputTag("selectedObjects", "electrons")),
                                                                  isPuppi = cms.bool(isPuppi),
                                                                  pfCandidates = cms.InputTag("packedPFCandidates")))
                

            if not hasattr(process,postfix+"t1phmet"):
                setattr(process, postfix+"t1phmet",cms.EDProducer("PhotonCorrectedMETProducer",
                                                                  met     = cms.InputTag(metName,"",miniAODProcess),
                                                                  cands   = cms.VInputTag(cms.InputTag("selectedObjects", "photons")),
                                                                  isPuppi = cms.bool(isPuppi),
                                                                  pfCandidates = cms.InputTag("packedPFCandidates")))
    
            if not hasattr(process,postfix+"t1taumet"):
                setattr(process, postfix+"t1taumet",cms.EDProducer("TauCorrectedMETProducer",
                                                                   met     = cms.InputTag(metName,"",miniAODProcess),
                                                                   cands   = cms.VInputTag(cms.InputTag("selectedObjects", "tausVLNew")),
                                                                   isPuppi = cms.bool(isPuppi),
                                                                   pfCandidates = cms.InputTag("packedPFCandidates")))
        
