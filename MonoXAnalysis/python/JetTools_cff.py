import os, copy, re
import FWCore.ParameterSet.Config as cms
from RecoJets.JetProducers.pileupjetidproducer_cfi import pileupJetIdCalculator,pileupJetIdEvaluator
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import patJetCorrFactors
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import patJetsUpdated
from CondCore.DBCommon.CondDBSetup_cfi import *

def addPileupJetID(process,collection, postfix, isMC = True):

    from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
    from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import patJets 

    ## load corrections
    process.load('JetMETCorrections.Configuration.JetCorrectors_cff')

    ## in case of puppi
    if "Puppi" in postfix or "puppi" in postfix or "PUPPI" in postfix:
        ## filter out puppi particles
        if not hasattr(process,'puppi'):
            setattr(process, 'puppi', cms.EDFilter('CandPtrSelector', 
                                                   src = cms.InputTag("packedPFCandidates"), 
                                                   cut = cms.string('puppiWeight > 0')) )

        ## cluster particles
        if not hasattr(process,'ak4PFJetsPuppi'):

            setattr(process, 'ak4PFJetsPuppi', ak4PFJets.clone(
                    src = cms.InputTag('puppi'),
                    doAreaFastjet = cms.bool(True),
                    rParam = cms.double(0.4),
                    jetAlgorithm = cms.string("AntiKt")))

            jetCollection = 'ak4PFJetsPuppi'
            
        ## generate objects for JEC
        if not hasattr(process,'ak4PFPuppiL1FastjetCorrector'):
            process.ak4PFPuppiL1FastjetCorrector = process.ak4PFCHSL1FastjetCorrector.clone(algorithm   = cms.string('AK4PFPuppi'))
        if not hasattr(process,'ak4PFPuppiL2RelativeCorrector'):
           process.ak4PFPuppiL2RelativeCorrector = process.ak4PFCHSL2RelativeCorrector.clone(algorithm   = cms.string('AK4PFPuppi'))
        if not hasattr(process,'ak4PFPuppiL3AbsoluteCorrector'):
            process.ak4PFPuppiL3AbsoluteCorrector = process.ak4PFCHSL3AbsoluteCorrector.clone(algorithm   = cms.string('AK4PFPuppi'))
        if not hasattr(process,'ak4PFPuppiL1FastL2L3Corrector'):
            process.ak4PFPuppiL1FastL2L3Corrector = process.ak4PFCHSL1FastL2L3Corrector.clone(
                correctors = cms.VInputTag("ak4PFPuppiL1FastjetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector"))
        if not hasattr(process,'ak4PFPuppiResidualCorrector'):
            process.ak4PFPuppiResidualCorrector = process.ak4PFResidualCorrector.clone( algorithm = 'AK4PFPuppi' )
        if not hasattr(process,'ak4PFPuppiL1FastL2L3ResidualCorrector'):
            process.ak4PFPuppiL1FastL2L3ResidualCorrector = process.ak4PFCHSL1FastL2L3ResidualCorrector.clone( 
                correctors = cms.VInputTag("ak4PFPuppiL1FastjetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector", "ak4PFPuppiResidualCorrector"))

        ## move from reco to patjets
        if not hasattr(process, 'patJetsAK4PFPuppi'):

            ## corrector
            puppiJEC = copy.deepcopy(process.JECLevels.labels)
            puppiJEC.remove('L1FastJet')
            setattr(process,"patJetCorrFactorsAK4Puppi", patJetCorrFactors.clone(
                    src     = cms.InputTag(jetCollection),
                    levels  = puppiJEC,
                    primaryVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                    payload = "AK4PFPuppi"));

            getattr(process,"patJetCorrFactorsAK4Puppi").useRho = False

            ## producer
            setattr(process,'patJetsAK4PFPuppi', patJets.clone(
                    jetSource = cms.InputTag(jetCollection),
                    jetCorrFactorsSource = cms.VInputTag(cms.InputTag('patJetCorrFactorsAK4Puppi')),
                    addBTagInfo          = cms.bool(False),
                    addDiscriminators    = cms.bool(False),
                    discriminatorSources  = cms.VInputTag('None'),
                    addAssociatedTracks    = cms.bool(False),
                    addJetCharge  = cms.bool(False),
                    addGenPartonMatch   = cms.bool(False),
                    embedGenPartonMatch = cms.bool(False),
                    addGenJetMatch      = cms.bool(False),
                    embedGenJetMatch    = cms.bool(False),
                    getJetMCFlavour    = cms.bool(False),
                    addJetFlavourInfo  = cms.bool(False)))
            
            jetCollectionPAT = 'patJetsAK4PFPuppi'

    else:
        ## filter out chs particles
        if not hasattr(process,'chs'):
            setattr(process, 'chs', cms.EDFilter('CandPtrSelector', src = cms.InputTag("packedPFCandidates"), cut = cms.string('fromPV')) )
            
        ## cluster particles
        if not hasattr(process, 'ak4PFJetsCHS'):
            setattr(process, 'ak4PFJetsCHS', ak4PFJets.clone(
                    src = cms.InputTag('chs'),
                    doAreaFastjet = cms.bool(True),
                    rParam = cms.double(0.4),
                    jetAlgorithm = cms.string("AntiKt"),
                    ))

            jetCollection = 'ak4PFJetsCHS'

            ## corrector
            setattr(process,"patJetCorrFactorsAK4CHS", patJetCorrFactors.clone(
                    src     = cms.InputTag(jetCollection),
                    levels  = process.JECLevels.labels,
                    primaryVertices = cms.InputTag('offlineSlimmedPrimaryVertices'),
                    payload = "AK4PFchs"));
            
            ## producer
            setattr(process,'patJetsAK4PFCHS', patJets.clone(
                    jetSource = cms.InputTag(jetCollection),
                    jetCorrFactorsSource = cms.VInputTag(cms.InputTag('patJetCorrFactorsAK4CHS')),
                    addBTagInfo          = cms.bool(False),
                    addDiscriminators    = cms.bool(False),
                    discriminatorSources  = cms.VInputTag('None'),
                    addAssociatedTracks    = cms.bool(False),
                    addJetCharge  = cms.bool(False),
                    addGenPartonMatch   = cms.bool(False),
                    embedGenPartonMatch = cms.bool(False),
                    addGenJetMatch      = cms.bool(False),
                    embedGenJetMatch    = cms.bool(False),
                    getJetMCFlavour    = cms.bool(False),
                    addJetFlavourInfo  = cms.bool(False)))
                    
            jetCollectionPAT = 'patJetsAK4PFCHS'
        
    ## add pileup jet id modules
    setattr(process,'pileupJetIdCalculator'+postfix,
            pileupJetIdCalculator.clone(
            jets = cms.InputTag(jetCollection),
            rho  = cms.InputTag("fixedGridRhoFastjetAll"),
            vertexes = cms.InputTag("offlineSlimmedPrimaryVertices"),
            applyJec = cms.bool(True),
            inputIsCorrected = cms.bool(False)))

    setattr(process,'pileupJetIdEvaluator'+postfix,
            pileupJetIdEvaluator.clone(
            jetids = cms.InputTag('pileupJetIdCalculator'+postfix),
            jets = cms.InputTag(jetCollection),
            rho = cms.InputTag("fixedGridRhoFastjetAll"),
            vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")))

    if "Puppi" in postfix or "puppi" in postfix or "PUPPI" in postfix:
        getattr(process,'pileupJetIdCalculator'+postfix).jec = cms.string("AK4PFPuppi")
        getattr(process,'pileupJetIdEvaluator'+postfix).jec = cms.string("AK4PFPuppi")


    ## add info inside PAT jets
    getattr(process,jetCollectionPAT).userData.userFloats.src = cms.VInputTag("pileupJetIdEvaluator"+postfix+":fullDiscriminant");

    ## matched orignal slimmed with PAT only for fullDiscriminant since RecoJetDeltaRValueMapProducer produces a value map of floats
    setattr(process,"puid"+postfix,cms.EDProducer("PATJetDeltaRValueMapProducer",
                                                  src = cms.InputTag(collection),
                                                  matched = cms.InputTag(jetCollectionPAT),
                                                  distMax = cms.double(0.4),
                                                  values = cms.vstring("userFloat('pileupJetIdEvaluator"+postfix+":fullDiscriminant')"),
                                                  valueLabels = cms.vstring("fullDiscriminant")
                                                  ))
    
    ## update original jets
    setattr(process,collection+"PUID",
            patJetsUpdated.clone(jetSource = cms.InputTag(collection),
                                 addJetCorrFactors = cms.bool(False),
                                 jetCorrFactorsSource = cms.VInputTag())
            )

        
    getattr(process,collection+"PUID").userData.userFloats = cms.PSet(
        src = cms.VInputTag("puid"+postfix+":fullDiscriminant"))
                            

### QGLikelihood adder
def addQGLikelihood(process,collection,postfix):

    CMSSW_VERSION = os.environ['CMSSW_VERSION'];

    ## connect to the DB
    if not hasattr(process,"QGPoolDBESSource"):

        process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
                                                CondDBSetup,
                                                toGet = cms.VPSet(),
                                                connect = cms.string('frontier://FrontierProd/CMS_COND_PAT_000'))
        qgDatabaseVersion = 'v1'

        for type in ['AK4PFchs','AK4PFchs_antib']:
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
                srcJets = cms.InputTag(collection),
                jetsLabel = cms.string('QGL_AK4PFchs'),
                srcVertexCollection   = cms.InputTag('offlineSlimmedPrimaryVertices')
                )
                )

        ## modify jets
        setattr(process,collection+"QG",
                patJetsUpdated.clone(jetSource = cms.InputTag(collection),
                                     addJetCorrFactors = cms.bool(False),
                                     jetCorrFactorsSource = cms.VInputTag()))


        getattr(process,collection+"QG").userData.userFloats = cms.PSet(
            src = cms.VInputTag('QGTagger'+postfix+':qgLikelihood'))
        

        

    
