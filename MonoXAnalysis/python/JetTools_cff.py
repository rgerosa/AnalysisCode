import os, copy, re
import FWCore.ParameterSet.Config as cms
from RecoJets.JetProducers.pileupjetidproducer_cfi import pileupJetId
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import patJetCorrFactors
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
from CondCore.DBCommon.CondDBSetup_cfi import *

## generic function that corrects jet and MET given a JEC                                                                                                                    
def JetCorrector(process,jetCollection,payloadName,isMC,applyL2L3Residuals):
    
    ## apply corrections on jets                                                                                                                                             
    if not hasattr(process,"patJetCorrFactorsReapplyJEC"+payloadName):
        setattr(process,"patJetCorrFactorsReapplyJEC"+payloadName, updatedPatJetCorrFactors.clone(
                src     = cms.InputTag(jetCollection),
                levels  = process.JECLevels.labels,
                payload = payloadName
                ));

    if "Puppi" in jetCollection or "PUPPI" in jetCollection: ## fix corrections for puppi jets removing L1                                                           
        puppiJEC = copy.deepcopy(process.JECLevels.labels)
        puppiJEC.remove('L1FastJet')
        getattr(process,"patJetCorrFactorsReapplyJEC"+payloadName).levels = puppiJEC
        getattr(process,"patJetCorrFactorsReapplyJEC"+payloadName).useRho = False

    if not hasattr(process,"slimmedJetsRecorrected"+payloadName):
        setattr(process,"slimmedJetsRecorrected"+payloadName,
                updatedPatJets.clone(jetSource = cms.InputTag(jetCollection),
                                     jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"+payloadName))))


    return "slimmedJetsRecorrected"+payloadName;

def addPileupJetID(process,jetCollection,postfix,isMC):

    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import patJets 
        
    ## add pileup jet id modules for the jetCollection
    setattr(process,'puid'+postfix,
            pileupJetId.clone(
            jets     = cms.InputTag(jetCollection),
            rho      = cms.InputTag("fixedGridRhoFastjetAll"),
            vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
            applyJec = cms.bool(True),
            inputIsCorrected = cms.bool(True))) ## if already calibrated jets are parsed, there is no need to re-calibrate

    if "Puppi" in postfix or "puppi" in postfix or "PUPPI" in postfix:
        getattr(process,'puid'+postfix).jec = cms.string("AK4PFPuppi")
        

    ## modify jets                                                                                                                                                           
    if not hasattr(process,jetCollection+"PUID"):
        setattr(process,jetCollection+"PUID",
                updatedPatJets.clone(jetSource = cms.InputTag(jetCollection),
                                     addJetCorrFactors = cms.bool(False),
                                     jetCorrFactorsSource = cms.VInputTag()))


        ## add info inside PAT jets
        getattr(process,jetCollection+"PUID").userData.userFloats.src = cms.VInputTag("puid"+postfix+":fullDiscriminant");

    return jetCollection+"PUID";

### QGLikelihood adder
def addQGLikelihood(process,jetCollection,postfix):

    CMSSW_VERSION = os.environ['CMSSW_VERSION'];

    ## connect to the DB
    if not hasattr(process,"QGPoolDBESSource"):

        process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
                                                CondDBSetup,
                                                toGet = cms.VPSet(),
                                                connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'))
        qgDatabaseVersion = 'v1'

        for type in ['AK4PFchs']:
            process.QGPoolDBESSource.toGet.extend(cms.VPSet(cms.PSet(
                        record = cms.string('QGLikelihoodRcd'),
                        tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_'+type),
                        label  = cms.untracked.string('QGL_'+type)
                        )))
        

        if re.match("CMSSW_7_6_.*",CMSSW_VERSION) and not hasattr(process,"es_prefer_QGL"):
            process.es_prefer_QGL = cms.ESPrefer("PoolDBESSource",'QGPoolDBESSource')   

    ## run evaluator
    from RecoJets.JetProducers.QGTagger_cfi import QGTagger

    if not hasattr(process,'QGTagger'+postfix):
        setattr(process,'QGTagger'+postfix, QGTagger.clone(
                srcJets = cms.InputTag(jetCollection),
                jetsLabel = cms.string('QGL_AK4PFchs'),
                srcVertexCollection   = cms.InputTag('offlineSlimmedPrimaryVertices')))

        ## modify jets
        setattr(process,jetCollection+"QG",
                updatedPatJets.clone(jetSource = cms.InputTag(jetCollection),
                                     addJetCorrFactors = cms.bool(False),
                                     jetCorrFactorsSource = cms.VInputTag()))

        
        getattr(process,jetCollection+"QG").userData.userFloats = cms.PSet(
            src = cms.VInputTag('QGTagger'+postfix+':qgLikelihood'))
        

    return jetCollection+"QG";

    
