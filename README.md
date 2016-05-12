============
AnalysisCode
============

Repository for analysis code

Recipe for 76X (temp fixes for pileup jet id):

       cmsrel CMSSW_7_6_4
       cd CMSSW_7_6_4/src
       cmsenv
       git cms-init
       #echo /CommonTools/PileupAlgos/ > .git/info/sparse-checkout
       #echo /CommonTools/Utils/ >> .git/info/sparse-checkout
       #echo /JetMETCorrections/Configuration/ >> .git/info/sparse-checkout
       #echo /JetMETCorrections/Modules/ >> .git/info/sparse-checkout
       #echo /JetMETCorrections/Type1MET/ >> .git/info/sparse-checkout
       #echo /PhysicsTools/PatAlgos/ >> .git/info/sparse-checkout
       #echo /PhysicsTools/PatUtils/ >> .git/info/sparse-checkout
       #echo /RecoMET/METAlgorithms/ >> .git/info/sparse-checkout
       #echo /RecoMET/METProducers/ >> .git/info/sparse-checkout
       git cms-merge-topic cms-met:metTool76X
       git cms-addpkg RecoJets/JetProducers
       git-cms-merge-topic jbrands:pileupJetId76X
       cd RecoJets/JetProducers/data/
       wget https://github.com/jbrands/RecoJets-JetProducers/raw/3dad903ed25d025f68be94d6f781ca957d6f86ac/pileupJetId_76x_Eta0to2p5_BDT.weights.xml.gz
       wget https://github.com/jbrands/RecoJets-JetProducers/raw/3dad903ed25d025f68be94d6f781ca957d6f86ac/pileupJetId_76x_Eta2p5to2p75_BDT.weights.xml.gz
       wget https://github.com/jbrands/RecoJets-JetProducers/raw/3dad903ed25d025f68be94d6f781ca957d6f86ac/pileupJetId_76x_Eta2p75to3_BDT.weights.xml.gz
       wget https://github.com/jbrands/RecoJets-JetProducers/raw/3dad903ed25d025f68be94d6f781ca957d6f86ac/pileupJetId_76x_Eta3to5_BDT.weights.xml.gz	   
       cd ../../..	
       git clone git@github.com:avartak/AnalysisCode.git -b Raffaele_76X
       git clone https://github.com/rfriese/RecoMET-METPUSubtraction data -b 74X-13TeV-Summer15-July2015
       rm -rf data/.git  
       rm JetMETCorrections/Modules/src/SealModule.cc
	   
How to Run the ntuple production (for analysis):

       cd CMSSW_7_6_4/src/AnalysisCode/MonoXAnalysis/test
       cmsRun tree.py <list of options>

Example:

	cmsRun tree.py isMC=True 	  


Options:   

	   isMC                : set to True when running on simulated events
	   isFastSIM           : set to True when running on fast sim miniAOD (still temp and relying on Emanuele miniAOD)
	   filterHighMETEvents : apply metCut GeV seletion on recoil/met (mumet or elmet or phmet or met)
	   metCut              : met threshold
	   filterOnHLT         : require the OR of all the trigger path definied in the MonoJetTreeMaker (beginRun)
	   usePrivateSQliteJEC : use a local db file for jet energy correction
	   usePrivateSQliteJER : use a local db file for jet energy smearings (by default true since no records in existing GT)
	   applyL2L3Residuals  : apply or not L2L3 residuals on data
	   addPileupJetID      : re-run the pileup-jet id on AK4PFchs jets (new tune for 76X)
	   addQGLikelihood     : add QGL value for AK4PFchs and AK4PFPuppi in case
	   addPuppiJets        : dump Puppi jet info in the MonoJetTreeMaker and add all the related modules in the process
	   addPuppiMET         : dump Puppi met info in the MonoJetTreeMaker and add all the related modules in the process
	   addMVAMet           : add the MVA met (still not working in 76X)
	   useMiniAODMet       : set to True when using miniAOD computed met and related systematics
	   addMETSystematics   : dump MET systematic variations in the tree and, when useMiniAODMet = False, re-evaluate thorugh METSystematicProducer (for PUPPI too)
	   useOfficialMETSystematics : use the official MET tool for systematics and not the private code
	   addSubstructureCHS  : compute all substructure info for PFCHS jets (different parameters can be changed directly in tree.py)
	   addSubstructurePuppi: compute all substructure info for PFPuppi jets 
	   processName         : name of the cms.Process
	   miniAODProcess      : name of the cms.Process used in miniAOD production --> useful for MET filters..etc
	   outputFileName      : name of the output root file from TFileService
	   globalTag           : GT to be used--> default for MC = 76X_mcRun2_asymptotic_RunIIFall15DR76_v1, default for DATA = 76X_dataRun2_16Dec2015_v0
	   JECEra              : JEC era string
	   useLHEWeights       : dump LHE weights and info (only for MC samples)
	   addQCDPDFWeights    : dump scalar and PDF variations (only for MC samples)
	   addGenParticles     : dump gen level info like W/Z boson and related decaying particles, gen photons, top and anti-top
	   isSignalSample      : dump gen level info for DM particles and mediator, as well as reading masses from the LHE header
	   crossSection        : fix the cross section branch in fb-1 from an external value
	   dropAnalyzerDumpEDM : to debug, avoid analyzers and dump edm collections produced in the event
	   nThreads            : number of threads
	   isCrab              : to be used to handle correctly local files with Crab

	   
Some useful information:

     Global Tag:
     	    https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideFrontierConditions

     Recipe for pileup-reweight:
     	    https://twiki.cern.ch/twiki/bin/view/CMS/PileupJSONFileforData#2015_Pileup_JSON_Files

	    Example:
	    pileupCalc.py -i Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt --inputLumiJSON pileup_JSON_11-19-2015.txt  --calcMode true --minBiasXsec 69000 --maxPileupBin 50 --numPileupBins 50 Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.root

     MiniAOD content:
     	    https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2015

     Recipe for muon ID:
     	    https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2

     Recipe for electron ID (cut based):
     	    https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2

     Recipe for photon ID (cut based):
     	    https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2

     Recipe for jet ID:
     	    https://twiki.cern.ch/twiki/bin/view/CMS/JetID

     Recipe for pileup jet ID:
     	    https://twiki.cern.ch/twiki/bin/view/CMS/PileupJetID

     Recipe for Puppi and Puppi MET:
     	    https://twiki.cern.ch/twiki/bin/viewauth/CMS/PUPPI

     Recipe for Quark-Gluon likelihood ID:
     	    https://twiki.cern.ch/twiki/bin/view/CMS/QuarkGluonLikelihood

     Recipe for B-tagging and scale factor:
     	    https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation

     Recipe for V-tagging:
     	    https://twiki.cern.ch/twiki/bin/view/CMS/JetWtagging

     Recipe for DoubleB-tagging:
     	    https://twiki.cern.ch/twiki/bin/view/CMS/Hbbtagging

     Recipe for JEC variations:
     	    https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC

     Recipe for JER smearings:
     	    https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution

     Recipe for MET filters:
     	    https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2

     Recipe for MVA MET:
     	    https://twiki.cern.ch/twiki/bin/viewauth/CMS/MVAMet
	 
     Recipe for Recoil Correction:
     	    https://twiki.cern.ch/twiki/bin/view/CMS/MissingETRun2RecoilCorrection

==========================
Combine Package and ROOT 6
==========================

Please have a look to the following twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideHiggsAnalysisCombinedLimit#ROOT6_SLC6_release_CMSSW_7_4_X

=================================
Run Crab Jobs for a given dataset
=================================

Python multicrab script can be used to submit on multiple dataset in one shot

       cd AnalysisCode/MonoXAnalysis/test/crab
       python <multicrab.py> <multicrab cfg> 

Examples on how to setup the python for data and MC are crab/multicrab_MC_76X.py crab/multicrab_DATA_76X.py:

	 pyCfgParams : list of all the parameters to be given to the cmsRun application when the job is executed on the grid
	 config.JobType.numCores : to allow multi-thread process .. the value should be the same of nThread
	 All the other parameters are the standard crab3 ones: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile

The config file which defines the dataset and additional parmaters for the cmsRun application can be found in SampleList_DATA_76X and SampleList_MC_76X. 

    samples[<name of crab di>] = [<dataset name>, [<additional parameters>]]
    
Example: 

    samples['ZJetsToNuNu_HT-100To200']  = ['/ZJetsToNuNu_HT-100To200_13TeV-madgraph/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM',
    					   ['useLHEWeights=True','addQCDPDFWeights=True','isSignalSample=False','addGenParticles=True','crossSection=280.5']]

=================================
Run Filter step for a dataset
=================================

==========================================
Run Filter step for all the crab outputs
==========================================

=================================
Copy Files to a local laptop
=================================

=================================
Merge directories in a clever way
=================================

=======================================
Run Basic makeTemplate for the analysis
======================================
