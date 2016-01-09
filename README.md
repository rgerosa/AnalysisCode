AnalysisCode
============

Repository for analysis code

Recipe for 74X:

	   cmsrel CMSSW_7_4_16
	   cd CMSSW_7_4_16/src
	   cmsenv
	   git cms-merge-topic matteosan1:egm_tnp_v7
	   git clone git@github.com:avartak/AnalysisCode.git -b Raffaele_74X
	   git clone https://github.com/rfriese/RecoMET-METPUSubtraction data -b 74X-13TeV-Summer15-July2015
	   rm -rf data/.git


Recipe for 76X:

	   cmsrel CMSSW_7_6_3
	   cd CMSSW_7_6_3/src
	   cmsenv
	   git clone git@github.com:avartak/AnalysisCode.git -b Raffaele_74X
	   git clone https://github.com/rfriese/RecoMET-METPUSubtraction data -b 74X-13TeV-Summer15-July2015
	   rm -rf data/.git


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