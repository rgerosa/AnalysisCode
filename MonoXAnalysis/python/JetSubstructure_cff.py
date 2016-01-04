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


    coneSizeStr = str(coneSize).replace("0","").replace(".","")
    jetAlgo = algo + coneSizeStr

    if algo == "AK" :
        ALGO = "AntiKt"
    elif algo == "CA":
        ALGO = "CambridgeAachen"


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


        getattr(process, jetCollection+postfix).src = cms.InputTag(jetCollection+"Reduced", 'constituents')
        getattr(process, jetCollection+postfix).rParam = coneSize
        getattr(process, jetCollection+postfix).jetAlgorithm = ALGO
        getattr(process, jetCollection+postfix).writeCompound = cms.bool(False)
        getattr(process, jetCollection+postfix).doAreaFastjet = cms.bool(True)
        getattr(process, jetCollection+postfix).jetCollInstanceName = cms.string("")

    ## gen groomed jets                                                                                                                                                       
    if isMC:
        if not hasattr(process,"genJets"+jetAlgo+postfix):
            setattr( process, "genJets"+jetAlgo+postfix,
                     getattr(process, jetCollection+postfix).clone(
                    src = cms.InputTag('genParticlesForJetsNoNu'),
                    writeCompound = cms.bool(False),
                    useExplicitGhosts = cms.bool(True),
                    jetPtMin = cms.double(100),
                    jetType = cms.string('GenJet')
                    ))


    ## pruned pat jets with embeded gen-jet info                                                                                                                             
    if not hasattr(process,"patJets"+jetCollection+postfix):
        addJetCollection(
            process,
            labelName = jetCollection+postfix,
            jetSource = cms.InputTag(jetCollection+postfix),
            algo = algo,
            rParam = coneSize,
            jetCorrections = (payloadName, JECLevel, 'None'),
            pfCandidates = cms.InputTag( pfCand ),
            svSource = cms.InputTag('slimmedSecondaryVertices'),
            genJetCollection = cms.InputTag("genJets"+jetAlgo+postfix),
            pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
            btagDiscriminators = ['None'], ## no b-tag info for pruned jets                                                                                                 
            getJetMCFlavour = isMC, ## no flavor info                                                                                                                       
            genParticles = cms.InputTag("prunedGenParticles")
            )

        if "Puppi" in pfCand or "puppi" in pfCand:
            getattr(process,"patJetCorrFactor"+jetCollection+postfix).useRho = cms.bool(False)
    
    if addQGLikelihood:
        if not hasattr(process,jetCollection+postfix+"QGL"):
                setattr(process,jetCollection+postfix+"QGL", QGTagger.clone(
                        srcJets = cms.InputTag(jetCollection+postfix),
                        jetsLabel = cms.string('QGL_AK4PFchs'),
                        srcVertexCollection   = cms.InputTag('offlineSlimmedPrimaryVertices')))


        getattr(process,'patJets'+jetCollection+postfix).userData.userFloats.src += [jetCollection+postfix+'QGL:qgLikelihood']



    ## matched fat jet with groomed one adding info as user float                                                                                                            
    if not hasattr(process,jetCollection+postfix+'Matched'):                                                                                                         
            setattr(process,jetCollection+postfix+'Matched',                                                                                                  
                    cms.EDProducer("PATRecoJetDeltaRValueMapProducer",                                                                                       
                                   src = cms.InputTag(jetCollection+"Reduced"),                                                                                       
                                   matched = cms.InputTag('patJets'+jetCollection+postfix),                                                                                     
                                   distMax = cms.double(coneSize),                                                                                                    
                                   values = cms.vstring("mass","pt","userFloat('"+jetCollection+"PrunedQGL:qgLikelihood')","bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')"),                                                                                           
                                   valueLabels = cms.vstring("mass","pt","qgLikelihood","bTagCSVIVFV2")))            
            
            if isMC:
                getattr(process,jetCollection+postfix+'Matched').valueLabels += ["hadronFlavour","partonFlavor","genMass","genPt"]
                getattr(process,jetCollection+postfix+'Matched').values += ["hadronFlavour","partonFlavour","genJet().mass","genJet().pt"]
 
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:mass']                                                         
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:pt']                                                          
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:qgLikelihood']                                                  
            getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:bTagCSVIVFV2']                                                  

            if isMC:
                getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:hadronFlavour']                                       
                getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:partonFlavor']                     
                getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:genMass']                                              
                getattr(process,'patJets'+jetCollection).userData.userFloats.src += [jetCollection+postfix+'Matched:genPt']                                                     

    if addSubJets:

        ## groomed sub-jets --> star with default clustering
        if not hasattr(process,jetCollection+postfix+'SubJets'):        
            setattr( process, jetCollection+postfix+'SubJets', 
                     getattr(process, jetCollection+postfix).clone(                    
                    src = cms.InputTag(jetCollection+"Reduced", 'constituents'),
                    writeCompound = cms.bool(True),
                    jetCollInstanceName=cms.string('SubJets')
                    ))

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
                jetSource = cms.InputTag(jetCollection+postfix+'SubJets', 'SubJets'),
                algo   = algo,  # needed for subjet b tagging
                rParam = coneSize,  # needed for subjet b tagging
                jetCorrections = (payloadNameSubjet, JECLevel, 'None'),
                pfCandidates = cms.InputTag( pfCand ),  
                pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'), 
                svSource = cms.InputTag('slimmedSecondaryVertices'), 
                getJetMCFlavour = isMC,
                genParticles = cms.InputTag("prunedGenParticles"),
                btagDiscriminators = btagDiscriminators,
                genJetCollection = cms.InputTag("genJets"+jetAlgo+postfix+"SubJets","SubJets"),
                explicitJTA  = True,  # needed for subjet b tagging
                svClustering = True, # needed for subjet b tagging
                fatJets = cms.InputTag(jetCollection+"Reduced"),             # needed for subjet flavor clustering
                groomedFatJets=cms.InputTag(jetCollection+postfix+"SubJets"), # needed for subjet flavor clustering
                ) 

            if "Puppi" in pfCand or "puppi" in pfCand:
                getattr(process,"patJetCorrFactor"+jetCollection+postfix+"SubJets").useRho = cms.bool(False)

        ## adding sub-jet QGL
        if addQGLikelihood:
            if not hasattr(process,jetCollection+postfix+"SubJetsQGL"):
                setattr(process,jetCollection+postfix+"SubJetsQGL", QGTagger.clone(
                        srcJets = cms.InputTag(jetCollection+postfix+"SubJets","SubJets"),
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
                    pileupMethod = "chs", selection = "pt > 150 && abs(eta) < 2.5", 
                    addPruning = True,   addPrunedSubjets = True,
                    addSoftDrop = True,  addSoftDropSubjets = True,
                    addTrimming = False, addTrimmedSubjets = False,
                    addFiltering = False,addFilteredSubjets = False, 
                    addNsubjettiness = True, addEnergyCorrelation = True, addQJets = False,
                    addQGLikelihood = True):
    
    print "############################"
    print "add substructure information"
    print "isMC = ",isMC;
    print "coneSize = ",coneSize;
    print "algo = ",algo;
    print "pileupMethod = ",pileupMethod;
    print "selection = ",selection;
    print "addPruning = ",addPruning;
    print "addPrunedSubjets = ",addPrunedSubjets;
    print "addSoftDrop = ",addSoftDrop;
    print "addSoftDropSubjets = ",addSoftDropSubjets;
    print "addTrimming = ",addTrimming;
    print "addTrimmedSubjets = ",addTrimmedSubjets;
    print "addFiltering = ",addFiltering;
    print "addFilteredSubjets = ",addFilteredSubjets;
    print "addNsubjettiness = ",addNsubjettiness;
    print "addEnergyCorrelation = ",addEnergyCorrelation;
    print "addQJets = ",addQJets;
    print "addQGLikelihood = ",addQGLikelihood;
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
                                  cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16")))

        if not hasattr(process,"genJets"+jetAlgo):
            setattr(process,"genJets"+jetAlgo, ak4GenJets.clone( src = 'genParticlesForJetsNoNu', 
                                                                 rParam       = coneSize, 
                                                                 jetAlgorithm = ALGO,
                                                                 jetPtMin = cms.double(100)))

        ## filter only W/Z and Higgs at generator level
        if not hasattr(process,"genBosons"):
            setattr(process,"genBosons",cms.EDFilter("CandPtrSelector",
                                                     src = cms.InputTag("prunedGenParticles"),
                                                     cut = cms.string("abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25 ")))
    
    ## b-tag discriminators to be considered
    bTagDiscriminators = [
        'pfTrackCountingHighEffBJetTags',
        'pfTrackCountingHighPurBJetTags',
        'pfSimpleSecondaryVertexHighEffBJetTags',
        'pfSimpleSecondaryVertexHighPurBJetTags',
        'pfCombinedSecondaryVertexV2BJetTags',
        'pfCombinedInclusiveSecondaryVertexV2BJetTags',
        'pfBoostedDoubleSecondaryVertexAK8BJetTags' ## new tag for Boosted double b-tagging
        ]

    CMSSW_VERSION = os.environ['CMSSW_VERSION'];
    if re.match("CMSSW_7_4_.*",CMSSW_VERSION) or re.match("CMSSW_7_5_.*",CMSSW_VERSION):
        bTagDiscriminators.remove('pfBoostedDoubleSecondaryVertexAK8BJetTags')

    ## jet energy corrections already loaded in JECLevels
    if pileupMethod != "chs" and pileupMethod != "Puppi" and pileupMethod != "":
        sys.exit("Invalid pileupMethod setup --> only recognized options are 'chs' and 'Puppi'");

    ## JEC
    JECLevel = copy.deepcopy(process.JECLevels.labels)
    if pileupMethod == "Puppi":
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

        payloadName = jetAlgo+'PF'+pileupMethod
        payloadNameSubjet = 'AK4PF'+pileupMethod
        jetCollection = jetAlgo+'PFJets'+pileupMethod.upper()
        pfCand = "chs"

    elif pileupMethod == "Puppi":
        if not hasattr(process,"puppi"):
            setattr( process, 'puppi', 
                     cms.EDFilter('CandPtrSelector', src = cms.InputTag("packedPFCandidates"), cut = cms.string('fromPV')), cut = cms.string('puppiWeight > 0'))
        
            
        if not hasattr(process,jetAlgo+'PFJets'+pileupMethod):
            setattr(process,jetAlgo+'PFJets'+pileupMethod,
                    ak4PFJetsPuppi.clone( src = cms.InputTag('puppi'),
                                          doAreaFastjet = True, 
                                          rParam = coneSize, 
                                          jetAlgorithm = ALGO ) )  

        payloadName = jetAlgo+'PF'+pileupMethod
        payloadNameSubjet = 'AK4PF'+pileupMethod
        jetCollection = jetAlgo+'PFJets'+pileupMethod;
        pfCand = "puppi"

    elif pileupMethod == "":

        if not hasattr(process,jetAlgo+'PFJets'):
            setattr(process,jetAlgo+'PFJets',
                    ak4PFJetsPuppi.clone( src = cms.InputTag('packedPFCandidates'),
                                          doAreaFastjet = True, 
                                          rParam = coneSize, 
                                          jetAlgorithm = ALGO ) )  
        payloadName = jetAlgo+'PF'
        payloadNameSubjet = 'AK4PF'
        jetCollection = jetAlgo+'PFJets'
        pfCand = "packedPFCandidates"
    
    
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
            algo   = algo,
            rParam = coneSize,
            jetCorrections = (payloadName, JECLevel, 'None'), 
            pfCandidates = cms.InputTag( pfCand ),  
            svSource = cms.InputTag('slimmedSecondaryVertices'),  
            genJetCollection = cms.InputTag("genJets"+jetAlgo),
            pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'), 
            btagDiscriminators = bTagDiscriminators,
            getJetMCFlavour = isMC,
            genParticles = cms.InputTag("prunedGenParticles"),
            ) 

        if "Puppi" in pfCand or "puppi" in pfCand:
            getattr(process,"patJetCorrFactor"+jetCollection).useRho = cms.bool(False)
        
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
                                                                              'Njettiness'+jetCollection+':tau4',
                                                                              ]

     ## add ECF
     if addEnergyCorrelation:
         if not hasattr(process,'ecf'+jetCollection):
             setattr(process,'ecf'+jetCollection,ECF.clone(
                 src = cms.InputTag(jetCollection+"Reduced"),
                 Njets = cms.vuint32(1,2,3),
                 ))

         getattr(process,'patJets'+jetCollection).userData.userFloats.src += ['ecf'+jetCollection+':ecf1'];
         getattr(process,'patJets'+jetCollection).userData.userFloats.src += ['ecf'+jetCollection+':ecf2'];
         getattr(process,'patJets'+jetCollection).userData.userFloats.src += ['ecf'+jetCollection+':ecf3'];

     ## add QJets
     if addQJets:
         if not hasattr(process,"qjets"+jetCollection):
             setattr(process,"qjets"+jetCollection,QJetsAdder.clone(
                 src = cms.InputTag(jetCollection+"Reduced"),
                 jeRad = cms.double(coneSize),
                 jetAlgo = cms.string(algo)
                 ))
         getattr(process,'patJets'+jetCollection).userData.userFloats.src += ["qjets"+jetCollection+":QjetsVolatility"]

            
    ### start with substructure: Pruning (run it on both reco and gen jets):
    Tags = [];
    Labels = [];
    if addPruning:        
        runGroomedMethod(process, 
                         isMC = isMC,
                         jetCollection = jetCollection,
                         coneSize = coneSize, 
                         algo = algo,
                         payloadName = payloadName, 
                         payloadNameSubjet = payloadNameSubjet,
                         JECLevel =  JECLevel,
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
                         JECLevel =  JECLevel,
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
                         JECLevel =  JECLevel,
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
                         JECLevel =  JECLevel,
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

        getattr(process,"packedPatJets"+jetCollection).algoTags = Tags;
        getattr(process,"packedPatJets"+jetCollection).algoLabels = Labels;


    if not isMC:

        from PhysicsTools.PatAlgos.tools.coreTools import removeMCMatching
        removeMCMatching(process, names=['Jets'], outputModules=['out'])
