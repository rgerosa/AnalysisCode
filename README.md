AnalysisCode
============

Repository for analysis code

Recipe for 74X

	   cmsrel CMSSW_7_4_16
	   cd CMSSW_7_4_16/src
	   cmsenv
	   git cms-merge-topic matteosan1:egm_tnp_v7
	   git clone git@github.com:avartak/AnalysisCode.git -b Raffaele_74X
	   git clone https://github.com/rfriese/RecoMET-METPUSubtraction data -b 74X-13TeV-Summer15-July2015
	   rm -rf data/.git


Recipe for 76X

	   cmsrel CMSSW_7_6_3
	   cd CMSSW_7_6_3/src
	   cmsenv
	   git clone git@github.com:avartak/AnalysisCode.git -b Raffaele_74X
	   git clone https://github.com/rfriese/RecoMET-METPUSubtraction data -b 74X-13TeV-Summer15-July2015
	   rm -rf data/.git
