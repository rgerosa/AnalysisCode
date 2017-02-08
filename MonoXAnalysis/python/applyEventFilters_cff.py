import os, copy
import FWCore.ParameterSet.Config as cms

########################################    
def applyEventFilters(process,
                      processName,
                      miniAODProcess,
                      filterHighMETEvents = True,
                      metCut = 0,
                      isPhotonPurity = False,
                      applyDiMuonFilter = False,
                      looseMuonPt = 10,
                      tightMuonPt = 20,
                      applyDiElectronFilter = False,
                      looseElectronPt = 10,
                      tightElectronPt = 35,
                      useMVAElectronID = False,
                      applyPhotonJetsFilter = False,
                      photonPt = 50,
                      useMVAPhotonID = False,
                      isReMiniAOD = False):


    ########################################    
    ########################################    
    ########################################    

    if isPhotonPurity and applyDiMuonFilter:
        print "Please check --> you can either selection isPhotonPurity or applyDiMuonFilter";
    if isPhotonPurity and applyDiElectronFilter:
        print "Please check --> you can either selection isPhotonPurity or applyDiElectronFilter";
    if isPhotonPurity and applyPhotonJetsFilter:
        print "Please check --> you can either selection isPhotonPurity or applyPhotonJetsFilter";
    if applyDiMuonFilter and applyDiElectronFilter:    
        print "Please check --> you can either selection applyDiMuonFilter or applyDiElectronFilter";
    if applyDiMuonFilter and applyPhotonJetsFilter:    
        print "Please check --> you can either selection applyDiMuonFilter or applyPhotonJetsFilter";
    if applyDiMuonFilter and applyDiElectronFilter:    
        print "Please check --> you can either selection applyDiMuonFilter or applyDiElectronFilter";
    if applyDiElectronFilter and applyPhotonJetsFilter:    
        print "Please check --> you can either selection applyDiElectronFilter or applyPhotonJetsFilter";

        
    ########################################    
    ########################################    
    if not isReMiniAOD:
        setattr(process,"filterHighRecoil", cms.EDFilter("PATMETFilter",
                                                         metCollections = cms.VPSet(
                    cms.PSet( srcMet = cms.InputTag("slimmedMETs","",processName),
                              metCut = cms.double(metCut)),
                    cms.PSet(srcMet = cms.InputTag("t1mumet","",processName),
                             metCut = cms.double(metCut)),
                    cms.PSet(srcMet = cms.InputTag("t1elmet","",processName),
                             metCut = cms.double(metCut)),
                    cms.PSet(srcMet = cms.InputTag("t1phmet","",processName),
                         metCut = cms.double(metCut)),
                    ),
                                                         filterEvents = cms.bool(filterHighMETEvents),
                                                         graterThan = cms.bool(True),
                                                         applyAndInsteadOfOr = cms.bool(False)
                                                         ))
    else:
        setattr(process,"filterHighRecoil", cms.EDFilter("PATMETFilter",
                                                         metCollections = cms.VPSet(
                    cms.PSet( srcMet = cms.InputTag("slimmedMETs","",miniAODProcess),
                              metCut = cms.double(metCut)),
                    cms.PSet( srcMet = cms.InputTag("slimmedMETsEGClean","",miniAODProcess),
                              metCut = cms.double(metCut)),
                    cms.PSet( srcMet = cms.InputTag("slimmedMETsEGMuClean","",miniAODProcess),
                              metCut = cms.double(metCut)),
                    cms.PSet(srcMet = cms.InputTag("t1mumet","",processName),
                             metCut = cms.double(metCut)),
                    cms.PSet(srcMet = cms.InputTag("t1mumetEGClean","",processName),
                             metCut = cms.double(metCut)),
                    cms.PSet(srcMet = cms.InputTag("t1mumetMuClean","",processName),
                             metCut = cms.double(metCut)),
                    cms.PSet(srcMet = cms.InputTag("t1elmet","",processName),
                             metCut = cms.double(metCut)),
                    cms.PSet(srcMet = cms.InputTag("t1elmetEGClean","",processName),
                             metCut = cms.double(metCut)),
                    cms.PSet(srcMet = cms.InputTag("t1elmetMuClean","",processName),
                             metCut = cms.double(metCut)),
                    cms.PSet(srcMet = cms.InputTag("t1phmet","",processName),
                         metCut = cms.double(metCut)),
                    cms.PSet(srcMet = cms.InputTag("t1phmetEGClean","",processName),
                         metCut = cms.double(metCut)),
                    cms.PSet(srcMet = cms.InputTag("t1phmetMuClean","",processName),
                         metCut = cms.double(metCut)),
                    ),
                                                         filterEvents = cms.bool(filterHighMETEvents),
                                                         graterThan = cms.bool(True),
                                                         applyAndInsteadOfOr = cms.bool(False)
                                                         ))

    

    recoilSequence = cms.Sequence(getattr(process,"filterHighRecoil"))
    puritySequence = cms.Sequence();
    diMuonSequence = cms.Sequence();
    diElectronSequence = cms.Sequence();
    photonSequence = cms.Sequence();

    ########################################    
    ### recoil cut can bias the measurement                                                                                                                                                       
    ########################################    
    if isPhotonPurity:
        if filterHighMETEvents and metCut != 0:
            for j in range(getattr(process,"filterHighRecoil").metCollections):
                getattr(process,"filterHighRecoil").metCollections[j].metCut = cms.double(0);
            
        setattr(process,"filterPhotonCandidates",cms.EDFilter("PhotonRefCountFilter",
                                                              src = cms.InputTag("selectedObjects","photonsPurity"),
                                                              minNumber = cms.int32(1),
                                                              maxNumber = cms.int32(999)
                                                              ))
        
        puritySequence += getattr(process,"filterPhotonCandidates")
        
    ########################################    
    ### ask two leptons with same flavor and opposite sign, at least one tight and Zmass                                                                                                    
    ########################################    
    if applyDiMuonFilter:

        setattr(process,"muonLooseCandidates",cms.EDFilter("MuonRefCountFilter",
                                                           src = cms.InputTag("selectedObjects","muons"),
                                                           minNumber = cms.int32(2),
                                                           maxNumber = cms.int32(2),
                                                           selection = cms.string("pt > "+str(looseMuonPt)+" && abs(eta) < 2.4"),
                                                           produceOutputCollection = cms.bool(True),
                                                           ))

        setattr(process,"muonTightCandidates",cms.EDFilter("MuonRefCountFilter",
                                                           src = cms.InputTag("selectedObjects","tightmuons"),
                                                           selection  = cms.string("pt > "+str(tightMuonPt)+" && (abs(eta) < 2.4)"),
                                                           minNumber = cms.int32(1),
                                                           maxNumber = cms.int32(2),
                                                           produceOutputCollection = cms.bool(True)
                                                           ))
        
        setattr(process,"diMuonCandidates",cms.EDProducer("CandViewShallowCloneCombiner",
                                                          decay = cms.string("muonTightCandidates@+ muonLooseCandidates@-"),
                                                          cut  = cms.string("60 < mass < 120 && charge=0"),
                                                          checkCharge = cms.bool(True),                                                          
                                                          ))
        
        setattr(process,"filterDiMuonCandidates",cms.EDFilter("CandViewCountFilter",
                                                              src = cms.InputTag("diMuonCandidates"),
                                                              minNumber = cms.uint32(1),
                                                              maxNumber = cms.uint32(999)
                                                              ))
        
        diMuonSequence += getattr(process,"muonLooseCandidates");
        diMuonSequence += getattr(process,"muonTightCandidates");
        diMuonSequence += getattr(process,"diMuonCandidates");
        diMuonSequence += getattr(process,"filterDiMuonCandidates");


    ########################################    
    ######## ask for two electrons, same flavor and opposite sign, at least one tight and Zmass
    ########################################    
    if applyDiElectronFilter:

        setattr(process,"electronLooseCandidates",cms.EDFilter("ElectronRefCountFilter",
                                                               src = cms.InputTag("selectedObjects","electrons"),
                                                               minNumber = cms.int32(2),
                                                               maxNumber = cms.int32(2),
                                                               selection = cms.string("pt > "+str(looseElectronPt)+" && abs(eta) < 2.5"),
                                                               produceOutputCollection = cms.bool(True)                                                               
                                                               ))
        
        setattr(process,"electronTightCandidates",cms.EDFilter("ElectronRefCountFilter",
                                                               src = cms.InputTag("selectedObjects","tightelectrons"),
                                                               minNumber = cms.int32(1),
                                                               maxNumber = cms.int32(2),
                                                               selection = cms.string("pt > "+str(tightElectronPt)+" && abs(eta) < 2.5"),
                                                               produceOutputCollection = cms.bool(True)
                                                               ))

        if useMVAElectronID:
            getattr(process,"electronTightCandidates").src = cms.InputTag("selectedObjects","mvatightelectrons"),
            

        setattr(process,"diElectronCandidates",cms.EDProducer("CandViewShallowCloneCombiner",
                                                              decay = cms.string("electronTightCandidates@+ electronLooseCandidates@-"),
                                                              cut  = cms.string("60 < mass < 120 && charge=0"),
                                                              checkCharge = cms.bool(True)
                                                              ))
        
        setattr(process,"filterDiElectronCandidates",cms.EDFilter("CandViewCountFilter",
                                                                  src = cms.InputTag("diElectronCandidates"),
                                                                  minNumber = cms.uint32(1),
                                                                  maxNumber = cms.uint32(999)
                                                                  ))

        diElectronSequence += getattr(process,"electronLooseCandidates");
        diElectronSequence += getattr(process,"electronTightCandidates");
        diElectronSequence += getattr(process,"diElectronCandidates");
        diElectronSequence += getattr(process,"filterDiElectronCandidates");
        
    ########################################    
    ########################################    

    if applyPhotonJetsFilter:
    
        setattr(process,"filterPhotonCandidates",cms.EDFilter("PhotonRefCountFilter",
                                                              src = cms.InputTag("selectedObjects","mediumphotons"),
                                                              minNumber = cms.int32(1),
                                                              maxNumber = cms.int32(1),
                                                              selection = cms.string("pt > "+str(photonPt)+" && abs(eta) < 2.5")
                                                              ))
        
        if useMVAElectronID:
            getattr(process,"filterPhotonCandidates").src = cms.InputTag("selectedObjects","mvatightphotons"),
            
        photonSequence +=  getattr(process,"filterPhotonCandidates");


    ########################################    
    ###### set final sequence
    ########################################    

    if not isPhotonPurity and not applyPhotonJetsFilter and not isPhotonPurity and not applyDiElectronFilter and not applyDiMuonFilter:
        setattr(process,"eventFilters", cms.Sequence(recoilSequence));
    elif isPhotonPurity:
        setattr(process,"eventFilters", cms.Sequence(recoilSequence+puritySequence));
    elif applyPhotonJetsFilter:
        setattr(process,"eventFilters", cms.Sequence(recoilSequence+photonSequence));
    elif applyDiMuonFilter:
        setattr(process,"eventFilters", cms.Sequence(recoilSequence+diMuonSequence));
    elif applyDiElectronFilter:
        setattr(process,"eventFilters", cms.Sequence(recoilSequence+diElectronSequence));
