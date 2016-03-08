import os, sys, copy, re
import FWCore.ParameterSet.Config as cms
from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets
from RecoJets.Configuration.RecoPFJets_cff import*
from RecoJets.JetProducers.ak4PFJetsPuppi_cfi import ak4PFJetsPuppi
from PhysicsTools.PatAlgos.tools.jetTools import addJetCollection
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
from RecoJets.JetProducers.QGTagger_cfi import QGTagger
from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
from RecoJets.JetProducers.ECF_cfi import ECF 
from RecoJets.JetProducers.qjetsadder_cfi import QJetsAdder
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import patJets
from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching

def runGroomedMethod(process, isMC,
                     jetCollection,
                     coneSize, algo,
                     payloadName, 
                     payloadNameSubjet,
                     JECLevel,
                     pfCand, btagDiscriminators,
                     addSubJets = True,
                     addQGLikelihood = True,
                     isPruning   = True,
                     isSoftDrop  = True,
                     isTrimming  = True, 
                     isFiltering = True):

    ### name for cone size and algo
    coneSizeStr = str(coneSize).replace("0","").replace(".","")
    jetAlgo = algo + coneSizeStr

    if algo == "AK" :
        ALGO = "AntiKt"
    elif algo == "CA":
        ALGO = "CambridgeAachen"

    btagSubjets = ['pfCombinedInclusiveSecondaryVertexV2BJetTags']

    ### form the right postfix
    if not isPruning and not isSoftDrop and not isTrimming and not isFiltering:
        sys.exit("runGroomedMethod: all bools are false --> please check")

    if isPruning and not isSoftDrop and not isTrimming and not isFiltering:
        postfix = "Pruned"
    elif not isPruning and isSoftDrop and not isTrimming and not isFiltering:
        postfix = "SoftDrop"
    elif not isPruning and not isSoftDrop and isTrimming and not isFiltering:
        postfix = "Trimmed"
    elif not isPruning and not isSoftDrop and not isTrimming and isFiltering:
        postfix = "Filtered"
    else:
        sys.exit("wrong setup, only one groomer is allowed each time --> please check");

    ## build default reco groomed jets                                                                                                                                   
    if not hasattr(process,jetCollection+postfix):
        if isPruning:
            setattr( process, jetCollection+postfix,ak8PFJetsCHSPruned.clone())
        elif isSoftDrop:
            setattr( process, jetCollection+postfix,ak8PFJetsCHSSoftDrop.clone())
        elif isTrimming:
            setattr( process, jetCollection+postfix,ak8PFJetsCHSTrimmed.clone())
        elif isFiltering:
            setattr( process, jetCollection+postfix,ak8PFJetsCHSFiltered.clone())


        getattr(process, jetCollection+postfix).src    = cms.InputTag(jetCollection+"Reduced", 'constituents')
        getattr(process, jetCollection+postfix).rParam = coneSize
        getattr(process, jetCollection+postfix).jetAlgorithm = ALGO
        ## write compounds means that subjets are stored as reco::Jet while the groomed jet is stored as base::Jet
        getattr(process, jetCollection+postfix).writeCompound = cms.bool(True)
        getattr(process, jetCollection+postfix).doAreaFastjet = cms.bool(True)

    ## gen groomed jets                                                                                                                                                       
    if isMC:
        if not hasattr(process,"genJets"+jetAlgo+postfix):
            ## store AK8 gen::Jet
            setattr( process, "genJets"+jetAlgo+postfix,
                     getattr(process, jetCollection+postfix).clone(
                    src = cms.InputTag('genParticlesForJetsNoNu'),
                    writeCompound = cms.bool(False),
                    useExplicitGhosts = cms.bool(True),
                    jetPtMin = cms.double(50.),
                    jetType = cms.string('GenJet'),
                    jetCollInstanceName = cms.string("")
                    ))


    ## pruned pat jets with embeded gen-jet info                                                                                                                             
    if not hasattr(process,"patJets"+jetCollection+postfix):
        addJetCollection(
            process,
            labelName = jetCollection+postfix,
            jetSource = cms.InputTag(jetCollection+postfix),
            algo      = algo,
            pfCandidates = cms.InputTag(pfCand),
            rParam    = coneSize,
            jetCorrections = (payloadName, JECLevel, 'None'),
            svSource       = cms.InputTag('slimmedSecondaryVertices'),
            genJetCollection = cms.InputTag("genJets"+jetAlgo+postfix),
            pvSource         = cms.InputTag('offlineSlimmedPrimaryVertices'),
            btagDiscriminators = btagDiscriminators, ## no b-tag info for pruned jets                                                                                     
            getJetMCFlavour  = isMC, ## no flavor info                                                                                                                       
            genParticles     = cms.InputTag("prunedGenParticles")
            )

        getattr(process,"patJetCorrFactors"+jetCollection+postfix).useRho = cms.bool(False)


    ## matched fat jet with groomed one adding info as user float                                                                                                            
    if not hasattr(process,jetCollection+postfix+'Matched'):                                                                                                         
        setattr(process,jetCollection+postfix+'Matched',                                                                                                  
                cms.EDProducer("RecoPATJetDeltaRValueMapProducer",                                                                                       
                               ## basic reco::jet ungroomed
                               src = cms.InputTag(jetCollection+"Reduced"),                                                                                       
                               ## mathched groomed pat jet
                               matched = cms.InputTag('patJets'+jetCollection+postfix),                                                                                     
                               distMax = cms.double(coneSize),                                                                                                    
                               values  = cms.vstring("mass","pt","eta","phi",
                                                     "bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",
                                                     "bDiscriminator('pfBoostedDoubleSecondaryVertexAK8BJetTags')",
                                                     "correctedP4(0).mass()","correctedP4(0).pt()"),                               
                               valueLabels = cms.vstring("mass","pt","eta","phi",
                                                         "pfCombinedInclusiveSecondaryVertexV2BJetTags",
                                                         "pfBoostedDoubleSecondaryVertexAK8BJetTags",
                                                         "rawmass","rawpt")))            
            
        if isMC:
            getattr(process,jetCollection+postfix+'Matched').valueLabels += ["hadronFlavour","partonFlavour","genMass","genPt","genEta","genPhi"]
            getattr(process,jetCollection+postfix+'Matched').values      += ["hadronFlavour","partonFlavour","genJet().mass","genJet().pt","genJet().eta","genJet().phi"]
 
        getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:mass']                                                         
        getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:pt']                                                          
        getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:pfCombinedInclusiveSecondaryVertexV2BJetTags']                  
        getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:pfBoostedDoubleSecondaryVertexAK8BJetTags']                  
        getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:eta']                                                          
        getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:phi']                                                          
        getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:rawmass']                                                          
        getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:rawpt']                                                          


        if isMC:
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:hadronFlavour']                                       
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:partonFlavour']                     
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:genMass']                                              
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:genPt']                                                     
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:genEta']                                                     
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:genPhi']                                                     

        
    ## add QGL --> some tricks are needed
    if addQGLikelihood:
        ## redo groomed jet with a special postfix (ForQGL)
        if not hasattr(process,jetCollection+postfix+"ForQGL"):
            setattr( process, jetCollection+postfix+"ForQGL",
                     getattr(process, jetCollection+postfix).clone(
                    writeCompound = cms.bool(False),
                    jetCollInstanceName = cms.string("")))

        ## run QGL evaluator    
        if not hasattr(process,jetCollection+postfix+"QGL"):
                setattr(process,jetCollection+postfix+"QGL", QGTagger.clone(
                        srcJets = cms.InputTag(jetCollection+postfix+"ForQGL"),
                        jetsLabel = cms.string('QGL_AK4PFchs'),
                        srcVertexCollection   = cms.InputTag('offlineSlimmedPrimaryVertices')))

        ## pattify jets on the fly
        if not hasattr(process,'patJets'+jetCollection+postfix+"QGL"):
            setattr(process,'patJets'+jetCollection+postfix+"QGL", patJets.clone(
                    jetSource = cms.InputTag(jetCollection+postfix+"ForQGL"),
                    addJetCorrFactors     = cms.bool(False), 
                    addBTagInfo           = cms.bool(False),
                    addDiscriminators     = cms.bool(False),
                    discriminatorSources  = cms.VInputTag('None'),
                    addAssociatedTracks   = cms.bool(False),
                    addJetCharge          = cms.bool(False),
                    addGenPartonMatch     = cms.bool(False),
                    embedGenPartonMatch   = cms.bool(False),
                    addGenJetMatch        = cms.bool(False),
                    embedGenJetMatch      = cms.bool(False),
                    getJetMCFlavour       = cms.bool(False),
                    addJetFlavourInfo     = cms.bool(False)))

            getattr(process,'patJets'+jetCollection+postfix+"QGL").userData.userFloats.src += [jetCollection+postfix+"QGL:qgLikelihood"]

        ## match QGL value with the original jet
        if not hasattr(process,jetCollection+postfix+"QGLMatched"):
            setattr(process,jetCollection+postfix+'QGLMatched',
                    cms.EDProducer("RecoPATJetDeltaRValueMapProducer",
                                   src = cms.InputTag(jetCollection+"Reduced"),
                                   matched = cms.InputTag('patJets'+jetCollection+postfix+"QGL"),
                                   distMax = cms.double(coneSize),
                                   values = cms.vstring("userFloat('"+jetCollection+postfix+"QGL:qgLikelihood')"),
                                   valueLabels = cms.vstring("qgLikelihood")))

            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'QGLMatched:qgLikelihood']


    ## add subjet information
    if addSubJets:

        ## gen groomed sub-jets --> star with default clustering
        if isMC and not hasattr(process,"genJets"+jetAlgo+postfix+"SubJets"):           
                setattr( process, "genJets"+jetAlgo+postfix+"SubJets",
                         getattr(process,"genJets"+jetAlgo+postfix).clone(
                        writeCompound = cms.bool(True),
                        jetCollInstanceName=cms.string('SubJets')
                        ))

        if not hasattr(process,"patJets"+jetCollection+postfix+"SubJets"):

            addJetCollection(
                process,
                labelName = jetCollection+postfix+'SubJets',
                jetSource = cms.InputTag(jetCollection+postfix, 'SubJets'),
                algo   = algo,  # needed for subjet b tagging
                rParam = coneSize,  # needed for subjet b tagging
                pfCandidates = cms.InputTag(pfCand),
                jetCorrections = (payloadNameSubjet, JECLevel, 'None'),
                pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'), 
                svSource = cms.InputTag('slimmedSecondaryVertices'), 
                getJetMCFlavour = isMC,
                genParticles = cms.InputTag("prunedGenParticles"),
                btagDiscriminators = btagSubjets,
                genJetCollection = cms.InputTag("genJets"+jetAlgo+postfix+"SubJets","SubJets"),
                explicitJTA  = True,  # needed for subjet b tagging
                svClustering = True, # needed for subjet b tagging
                fatJets = cms.InputTag(jetCollection+"Reduced"),             # needed for subjet flavor clustering
                groomedFatJets=cms.InputTag(jetCollection+postfix), # needed for subjet flavor clustering
                ) 

            getattr(process,"patJetCorrFactors"+jetCollection+postfix+"SubJets").useRho = cms.bool(False)

        ## adding sub-jet QGL
        if addQGLikelihood:
            if not hasattr(process,jetCollection+postfix+"SubJetsQGL"):
                setattr(process,jetCollection+postfix+"SubJetsQGL", QGTagger.clone(
                        srcJets = cms.InputTag(jetCollection+postfix,"SubJets"),
                        jetsLabel = cms.string('QGL_AK4PFchs'),
                        srcVertexCollection   = cms.InputTag('offlineSlimmedPrimaryVertices')
                        ))
                
                
                getattr(process,'patJets'+jetCollection+postfix+"SubJets").userData.userFloats.src += [jetCollection+postfix+'SubJetsQGL:qgLikelihood']
        

        ## Establish references between PATified fat jets and subjets using the BoostedJetMerger
        if not hasattr(process,'patJets'+jetCollection+postfix+'Packed'):
            setattr( process, 'patJets'+jetCollection+postfix+'Packed', 
                     cms.EDProducer("BoostedJetMerger",
                                    jetSrc = cms.InputTag('selectedPatJets'+jetCollection+postfix),
                                    subjetSrc = cms.InputTag("selectedPatJets"+jetCollection+postfix+"SubJets")
                                    ))


def JetSubstructure(process,
                    isMC, 
                    coneSize = 0.8, algo = "AK", 
                    pileupMethod = "chs", 
                    selection = "pt > 150 && abs(eta) < 2.5", 
                    addPruning   = True,   addPrunedSubjets   = True,
                    addSoftDrop  = True,   addSoftDropSubjets = True,
                    addTrimming  = False,  addTrimmedSubjets  = False,
                    addFiltering = False,  addFilteredSubjets = False, 
                    addNsubjettiness = True, addEnergyCorrelation = True, addQJets = False,
                    addQGLikelihood  = True):
    
    print "############################"
    print "add substructure information"
    print "isMC             = ",isMC;
    print "coneSize         = ",coneSize;
    print "algo             = ",algo;
    print "pileupMethod     = ",pileupMethod;
    print "selection        = ",selection;
    print "addPruning       = ",addPruning;
    print "addPrunedSubjets = ",addPrunedSubjets;
    print "addSoftDrop      = ",addSoftDrop;
    print "addSoftDropSubjets = ",addSoftDropSubjets;
    print "addTrimming        = ",addTrimming;
    print "addTrimmedSubjets  = ",addTrimmedSubjets;
    print "addFiltering       = ",addFiltering;
    print "addFilteredSubjets = ",addFilteredSubjets;
    print "addNsubjettiness   = ",addNsubjettiness;
    print "addEnergyCorrelation = ",addEnergyCorrelation;
    print "addQJets             = ",addQJets;
    print "addQGLikelihood      = ",addQGLikelihood;
    print "############################"

    ## build jet algo name for jet clustering
    coneSizeStr = str(coneSize).replace("0","").replace(".","")
    jetAlgo = algo + coneSizeStr
    if algo == "AK" :
        ALGO = "AntiKt"
    elif algo == "CA":
        ALGO = "CambridgeAachen"
    
    if isMC : 
        ## run gen jet clustering from genParticles
        if not hasattr(process,"genParticlesForJetsNoNu"):
            setattr( process, 'genParticlesForJetsNoNu', 
                     cms.EDFilter("CandPtrSelector", 
                                  src = cms.InputTag("packedGenParticles"), 
                                  cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16 && abs(pdgId) != 1000022 && abs(pdgId) != 1000012 &&"+
                                                   "abs(pdgId) != 1000014 && abs(pdgId) != 1000016 && abs(pdgId) != 2000012 && abs(pdgId) != 2000014 &&"+
                                                   "abs(pdgId) != 2000016 && abs(pdgId) != 1000039 && abs(pdgId) != 5100039 && abs(pdgId) != 4000012 &&"+
                                                   "abs(pdgId) != 4000014 && abs(pdgId) != 4000016 && abs(pdgId) != 9900012 && abs(pdgId) != 9900014 &&"+
                                                   "abs(pdgId) != 9900016 && abs(pdgId) != 39 && abs(pdgId) != 9100012"))) ## skip all the possible invisible particles

        if not hasattr(process,"genJets"+jetAlgo):
            setattr(process,"genJets"+jetAlgo, ak4GenJets.clone( src = 'genParticlesForJetsNoNu', 
                                                                 rParam       = coneSize, 
                                                                 jetAlgorithm = ALGO,
                                                                 jetPtMin = cms.double(50.))) ## fix a lower reasonable threshold

        ## filter only hadronically decaying W/Z and Higgs at generator level
        if not hasattr(process,"genBosons"):
            setattr(process,"genBosons",cms.EDFilter("CandPtrSelector",
                                                     src = cms.InputTag("prunedGenParticles"),
                                                     cut = cms.string("(abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25) &&"+ 
                                                                      "numberOfDaughters > 1 && abs(daughter(0).pdgId) <= 5 && abs(daughter(1).pdgId) <= 5")))

    ## b-tag discriminators to be considered
    bTagDiscriminators = [
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'pfBoostedDoubleSecondaryVertexAK8BJetTags' ## new tag for Boosted double b-tagging
        ]


    ## jet energy corrections already loaded in JECLevels
    if pileupMethod != "chs" and pileupMethod != "Puppi" and pileupMethod != "":
        sys.exit("Invalid pileupMethod setup --> only recognized options are 'chs' and 'Puppi'");

    ## JEC
    JECLevel = copy.deepcopy(process.JECLevels.labels)
    if pileupMethod == "Puppi" and 'L1FastJet' in JECLevel:
        JECLevel.remove('L1FastJet')
        
    payloadName       = ""
    payloadNameSubJet = ""
    jetCollection = ""
    pfCand        = ""

    ## in case of CHS select particles and re-cluster jets
    if pileupMethod == "chs":
        if not hasattr(process,"chs"):
            setattr( process, 'chs', 
                     cms.EDFilter('CandPtrSelector', src = cms.InputTag("packedPFCandidates"), cut = cms.string('fromPV')) )

        if not hasattr(process,jetAlgo+'PFJets'+pileupMethod.upper()):
            setattr(process,jetAlgo+'PFJets'+pileupMethod.upper(),
                    ak4PFJetsCHS.clone( src = cms.InputTag('chs'), 
                                        doAreaFastjet = True, 
                                        rParam = coneSize, 
                                        jetAlgorithm = ALGO ) ) 

        payloadName       = jetAlgo+'PF'+pileupMethod
        payloadNameSubjet = 'AK4PF'+pileupMethod
        jetCollection     = jetAlgo+'PFJets'+pileupMethod.upper()
        pfCand            = "chs"

    ## for puppi jets
    elif pileupMethod == "Puppi":
        if not hasattr(process,"puppi"):
            setattr( process, 'puppi', 
                     cms.EDFilter('CandPtrSelector', 
                                  src = cms.InputTag("packedPFCandidates"), 
                                  cut = cms.string('puppiWeight > 0')))
        
            
        if not hasattr(process,jetAlgo+'PFJets'+pileupMethod):
            setattr(process,jetAlgo+'PFJets'+pileupMethod,
                    ak4PFJetsPuppi.clone( src = cms.InputTag('puppi'),
                                          doAreaFastjet = True, 
                                          rParam = coneSize, 
                                          jetAlgorithm = ALGO ) )  

        payloadName       = jetAlgo+'PF'+pileupMethod
        payloadNameSubjet = 'AK4PF'+pileupMethod
        jetCollection     = jetAlgo+'PFJets'+pileupMethod;
        pfCand            = "puppi"

    ## standard PF jets
    elif pileupMethod == "":

        if not hasattr(process,jetAlgo+'PFJets'):
            setattr(process,jetAlgo+'PFJets',
                    ak4PFJetsPuppi.clone( src = cms.InputTag('packedPFCandidates'),
                                          doAreaFastjet = True, 
                                          rParam = coneSize, 
                                          jetAlgorithm = ALGO ) )  
        payloadName       = jetAlgo+'PF'
        payloadNameSubjet = 'AK4PF'
        jetCollection     = jetAlgo+'PFJets'
        pfCand            = "packedPFCandidates"

    ## apply selection and produce a restricted set of consituents only for jets passig the selection
    if not hasattr(process,jetCollection+"Reduced"):
        setattr(process,jetCollection+"Reduced",
                cms.EDFilter("MiniAODJetConstituentSelector",
                             src = cms.InputTag(jetCollection), 
                             cut = cms.string(selection)))

    ## build pat-jets from this skimmed collection: example
    if not hasattr(process,"patJets"+jetCollection):

        addJetCollection(
            process,
            labelName = jetCollection,
            jetSource = cms.InputTag(jetCollection+"Reduced"),
            algo      = algo,
            rParam    = coneSize,
            pfCandidates   = cms.InputTag(pfCand),
            jetCorrections = (payloadName, JECLevel, 'None'), 
            svSource       = cms.InputTag('slimmedSecondaryVertices'),  
            genJetCollection = cms.InputTag("genJets"+jetAlgo),
            pvSource           = cms.InputTag('offlineSlimmedPrimaryVertices'), 
            btagDiscriminators = bTagDiscriminators,
            getJetMCFlavour    = isMC,
            genParticles       = cms.InputTag("prunedGenParticles"),
            ) 

        if "Puppi" in pfCand or "puppi" in pfCand:
            getattr(process,"patJetCorrFactors"+jetCollection).useRho = cms.bool(False)


    ## match reco-jets with hadronically decaying genBosons (W,Z,H)
    if isMC:
        if not hasattr(process,jetCollection+'GenBosonMatched'):                                                                                                
            setattr(process,jetCollection+'GenBosonMatched',                                                                                                  
                    cms.EDProducer("RecoJetCandDeltaRValueMapProducer",                                                                                       
                                   ## basic reco::jet ungroomed
                                   src = cms.InputTag(jetCollection+"Reduced"),                                                                                       
                                   ## mathched groomed pat jet
                                   matched = cms.InputTag('genBosons'),                                                                                     
                                   distMax = cms.double(coneSize),                                                                                                    
                                   values  = cms.vstring("pt","eta","phi","mass"),
                                   valueLabels = cms.vstring("pt","eta","phi","mass")))            

            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+'GenBosonMatched:pt']
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+'GenBosonMatched:eta']
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+'GenBosonMatched:phi']
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+'GenBosonMatched:mass']


    ## add QGLikelihood on fat jet
    if addQGLikelihood:
     if not hasattr(process,jetCollection+"QGL"):
            setattr(process,jetCollection+"QGL", QGTagger.clone(
                    srcJets = cms.InputTag(jetCollection+"Reduced"),
                    jetsLabel = cms.string('QGL_AK4PFchs'),
                    srcVertexCollection   = cms.InputTag('offlineSlimmedPrimaryVertices')))


            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+'QGL:qgLikelihood']


     ## addNsubjettiness
     if addNsubjettiness:
         if not hasattr(process,'Njettiness'+jetCollection):
             setattr(process,'Njettiness'+jetCollection, 
                     Njettiness.clone(
                     src = cms.InputTag(jetCollection+"Reduced"),
                     cone = cms.double(coneSize)));
         
             getattr(process,'patJets'+jetCollection).userData.userFloats.src += ['Njettiness'+jetCollection+':tau1',
                                                                                  'Njettiness'+jetCollection+':tau2',
                                                                                  'Njettiness'+jetCollection+':tau3',
                                                                                  'Njettiness'+jetCollection+':tau4']

         ## on gen jets
         if isMC:
             if not hasattr(process,"NjettinessGenJets"+jetAlgo):
                 setattr(process,"NjettinessGenJets"+jetAlgo,
                         Njettiness.clone(
                         src = cms.InputTag("genJets"+jetAlgo),
                         cone = cms.double(coneSize)));

                 ## pattify gen jets --> temp
                 if not hasattr(process,'patGenJets'+jetAlgo):
                     setattr(process,'patGenJets'+jetAlgo, patJets.clone(
                             jetSource = cms.InputTag("genJets"+jetAlgo),
                             addJetCorrFactors     = cms.bool(False),
                             addBTagInfo           = cms.bool(False),
                             addDiscriminators     = cms.bool(False),
                             discriminatorSources  = cms.VInputTag('None'),
                             addAssociatedTracks   = cms.bool(False),
                             addJetCharge          = cms.bool(False),
                             addGenPartonMatch     = cms.bool(False),
                             embedGenPartonMatch   = cms.bool(False),
                             addGenJetMatch        = cms.bool(False),
                             embedGenJetMatch      = cms.bool(False),
                             getJetMCFlavour       = cms.bool(False),
                             addJetFlavourInfo     = cms.bool(False)))

                     getattr(process,'patGenJets'+jetAlgo).userData.userFloats.src += ["NjettinessGenJets"+jetAlgo+':tau1']
                     getattr(process,'patGenJets'+jetAlgo).userData.userFloats.src += ["NjettinessGenJets"+jetAlgo+':tau2']
                     getattr(process,'patGenJets'+jetAlgo).userData.userFloats.src += ["NjettinessGenJets"+jetAlgo+':tau3']
                     getattr(process,'patGenJets'+jetAlgo).userData.userFloats.src += ["NjettinessGenJets"+jetAlgo+':tau4']
                     

             if not hasattr(process,jetCollection+'GenNjettinessMatched'):
                 setattr(process,jetCollection+'GenNjettinessMatched',
                         cms.EDProducer("RecoPATJetDeltaRValueMapProducer",
                                        ## basic reco::jet ungroomed                                                                                                         
                                        src = cms.InputTag(jetCollection+"Reduced"),
                                        ## mathched groomed pat jet                                                                                                       
                                        matched = cms.InputTag('patGenJets'+jetAlgo),
                                        distMax = cms.double(coneSize),
                                        values  = cms.vstring("userFloat('NjettinessGenJets"+jetAlgo+":tau1')",
                                                              "userFloat('NjettinessGenJets"+jetAlgo+":tau2')",
                                                              "userFloat('NjettinessGenJets"+jetAlgo+":tau3')",
                                                              "userFloat('NjettinessGenJets"+jetAlgo+":tau4')"),
                                        valueLabels = cms.vstring("tau1","tau2","tau3","tau4")))

                  
                 getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+'GenNjettinessMatched:tau1',
                                                                                      jetCollection+'GenNjettinessMatched:tau2',
                                                                                      jetCollection+'GenNjettinessMatched:tau3',
                                                                                      jetCollection+'GenNjettinessMatched:tau4']
                                                                                      
     

     ## add ECF
     if addEnergyCorrelation:
         if not hasattr(process,'ecf'+jetCollection):
             setattr(process,'ecf'+jetCollection,ECF.clone(
                 src = cms.InputTag(jetCollection+"Reduced"),
                 Njets = cms.vuint32(1,2,3)))

             getattr(process,'patJets'+jetCollection).userData.userFloats.src += ['ecf'+jetCollection+':ecf1'];
             getattr(process,'patJets'+jetCollection).userData.userFloats.src += ['ecf'+jetCollection+':ecf2'];
             getattr(process,'patJets'+jetCollection).userData.userFloats.src += ['ecf'+jetCollection+':ecf3'];

     ## add QJets
     if addQJets:
         if not hasattr(process,"qjets"+jetCollection):
             setattr(process,"qjets"+jetCollection,QJetsAdder.clone(
                 src = cms.InputTag(jetCollection+"Reduced"),
                 jeRad = cms.double(coneSize),
                 jetAlgo = cms.string(algo)))

             getattr(process,'patJets'+jetCollection).userData.userFloats.src += ["qjets"+jetCollection+":QjetsVolatility"]

            
    ### start with substructure: Pruning (run it on both reco and gen jets):
    Tags   = [];
    Labels = [];
    JECLevelTemp = copy.deepcopy(JECLevel)
    if 'L1FastJet' in JECLevelTemp:
        JECLevelTemp.remove('L1FastJet') ## in any case groomers are removing already pileup

    if addPruning:        
        runGroomedMethod(process, 
                         isMC = isMC,
                         jetCollection = jetCollection,
                         coneSize = coneSize, 
                         algo = algo,
                         payloadName = payloadName, 
                         payloadNameSubjet = payloadNameSubjet,
                         JECLevel =  JECLevelTemp,
                         pfCand = pfCand, 
                         btagDiscriminators = bTagDiscriminators,
                         addSubJets      = addPrunedSubjets,
                         addQGLikelihood = addQGLikelihood,
                         isPruning   = True,
                         isSoftDrop  = False,
                         isTrimming  = False,
                         isFiltering = False);

        if addPrunedSubjets:
            Tags   += ['patJets'+jetCollection+'PrunedPacked']
            Labels += ['Pruned']

    if addSoftDrop:        
        runGroomedMethod(process, 
                         isMC = isMC,
                         jetCollection = jetCollection,
                         coneSize = coneSize, 
                         algo = algo,
                         payloadName = payloadName, 
                         payloadNameSubjet = payloadNameSubjet,
                         JECLevel =  JECLevelTemp,
                         pfCand = pfCand, 
                         btagDiscriminators = bTagDiscriminators,
                         addSubJets      = addSoftDropSubjets,
                         addQGLikelihood = addQGLikelihood,
                         isPruning   = False,
                         isSoftDrop  = True,
                         isTrimming  = False,
                         isFiltering = False);

        if addSoftDropSubjets :
            Tags   += ['patJets'+jetCollection+'SoftDropPacked']
            Labels += ['SoftDrop']

    if addTrimming:        
        runGroomedMethod(process, 
                         isMC = isMC,
                         jetCollection = jetCollection,
                         coneSize = coneSize, 
                         algo = algo,
                         payloadName = payloadName, 
                         payloadNameSubjet = payloadNameSubjet,
                         JECLevel =  JECLevelTemp,
                         pfCand = pfCand, 
                         btagDiscriminators = bTagDiscriminators,
                         addSubJets      = addTrimmedSubjets,
                         addQGLikelihood = addQGLikelihood,
                         isPruning   = False,
                         isSoftDrop  = False,
                         isTrimming  = True,
                         isFiltering = False);

        if addTrimmedSubjets:
            Tags   += ['patJets'+jetCollection+'TrimmedPacked']
            Labels += ['Trimmed']


    if addFiltering:        
        runGroomedMethod(process, 
                         isMC = isMC,
                         jetCollection = jetCollection,
                         coneSize = coneSize, 
                         algo = algo,
                         payloadName = payloadName, 
                         payloadNameSubjet =  payloadNameSubjet,
                         JECLevel =  JECLevelTemp,
                         pfCand = pfCand, 
                         btagDiscriminators = bTagDiscriminators,
                         addSubJets      = addFilteredSubjets,
                         addQGLikelihood = addQGLikelihood,
                         isPruning   = False,
                         isSoftDrop  = False,
                         isTrimming  = False,
                         isFiltering = True);

        if addFilteredSubjets:
            Tags   += ['patJets'+jetCollection+'FilteredPacked']
            Labels += ['Filtered']

    ## finally fix the pat fat jet collection
    if not hasattr(process,"packedPatJets"+jetCollection):
        setattr(process,"packedPatJets"+jetCollection,
                cms.EDProducer("JetSubstructurePacker",
                               jetSrc = cms.InputTag("selectedPatJets"+jetCollection),
                               distMax = cms.double(coneSize),
                               fixDaughters = cms.bool(False),
                               algoTags = cms.VInputTag(),
                               algoLabels = cms.vstring()))

        getattr(process,"packedPatJets"+jetCollection).algoTags   = Tags;
        getattr(process,"packedPatJets"+jetCollection).algoLabels = Labels;


    if not isMC:
        removeMCMatching(process, names=['Jets'], outputModules=[])

    return "packedPatJets"+jetCollection
