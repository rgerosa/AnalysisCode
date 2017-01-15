#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#include "makehist.h"
#include "makeCorrHistograms.C"
#include "makeDataHistograms.C"
#include "makeSignalHistograms.C"
#include "makeTemplatesUtils.h"

using namespace std;

// Run the final analysis:
// 1) Store all corrections templates from input files (complient to combine)
// 2) Make data and expected yields templates for all the other processes
static bool makeResonantSelection = false; // split top in resonant and non resonant
static int  typeOfDMSignal        = 0;     // 0 means both mono-j and mono-V, 1 is mono-j, 2 is mono-V
static bool runHiggsInvisible     = false; // run Higgs invisible analysis
static bool addTop                = false;
static bool addQCD                = false;
static bool addWgamma             = true; 
static bool addZgamma             = true;
static bool skipTFsystematics     = false;
static bool skipDataAnalysis      = false;
static SamplesNLO nloSamples (false,false,false,false);
static bool useTheoriestKFactors  = false;
static bool useNewTheoryUncertainty = false;

void makeTemplates(bool doCorrectionHistograms   = false,  // calculate transfer factors and sys
		   bool skipCorrectionHistograms = false,  // skip to open and dump transfer factors
		   Category category             = Category::monojet,  // 0 = inclusive mono-j, 1 = exclsuive mono-j, 2 V-tag HP ..
		   double lumi                   = 36.46, // 
		   string outDir                 = "", // output dir for template file
		   string templateSuffix         = "",  // suffix for the output file
		   vector<string> observables    = {"met"}, // 1D histo
		   vector<string> observables_2D = {},  // 2D histo
		   bool doShapeSystematics    = false, // run all the met, b-tag shape variations
		   bool runOnlySignal         = false, // produce a file with only signal templates
		   bool runOnlyBackground     = false, // produce a file with only background templates
		   bool applyPostFitWeights   = false,
		   bool addHistoForCutAndCount= false) {

  system(("mkdir -p "+outDir).c_str());

  // to initialize the binning map
  initializeBinning();

  // find all possible mass pont to use in the analysis for each Model: Vector, Axial, Scalar and Pseudoscalar .. if onlyMonoJetSignal is true just use all the available mono-j signal
  vector<signalSample> signalMassPoint;
  if(not skipDataAnalysis){
    if(not runHiggsInvisible and not runOnlyBackground){
      findAllPossibleMassPoints(signalMassPoint,"Vector",typeOfDMSignal);  
      findAllPossibleMassPoints(signalMassPoint,"Axial",typeOfDMSignal);
      findAllPossibleMassPoints(signalMassPoint,"Scalar",typeOfDMSignal);
      findAllPossibleMassPoints(signalMassPoint,"Pseudoscalar",typeOfDMSignal);
    }
  }

  //////////////////////////
  nloSamples.WJetsDIR = "WJets";
  nloSamples.ZJetsDIR = "ZJets";
  nloSamples.DYJetsDIR = "DYJets";
  nloSamples.PhotonJetsDIR = "PhotonJets";
  
  if(nloSamples.useWJetsNLO)
    nloSamples.WJetsDIR = "WJetsNLO";
  if(nloSamples.useZJetsNLO)
    nloSamples.ZJetsDIR = "ZJetsNLO";
  if(nloSamples.useDYJetsNLO)
    nloSamples.DYJetsDIR = "DYJetsNLO";
  if(nloSamples.usePhotonJetsNLO)
    nloSamples.PhotonJetsDIR = "PhotonJetsNLO";
  ////////////////////////////

  if(doCorrectionHistograms){    

    if (not skipTFsystematics and category != Category::monojet and useNewTheoryUncertainty){
      cerr<<"Protection --> new theory uncertainty can be used only for monojet category"<<endl;
      return;
    }

    cout<<"make correction histogram for Zmm to Znn"<<endl;      
    makezmmcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		   baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/",
		   category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors); 

    cout<<"make correction histogram for Zee to Znn"<<endl;
    makezeecorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		   baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zeefilter/",
		   category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors); 

    cout<<"make correction histogram for Wmn to WJets"<<endl;
    makewmncorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		   baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/",
		   category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors); 
    
    cout<<"make correction histogram for Wen to WJets"<<endl;
    makewencorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		   baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wenfilter/",
		   category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors); 

    if(category != Category::VBF or (category == Category::VBF and addZgamma)){
      
      cout<<"make correction histogram for Gam+jets to Znn"<<endl;
      makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		     baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/", 		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false);
          
      cout<<"systematics on Z/gamma ratio --> NLO QCD "<<endl;
      makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		     baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"qcd",1);
      
      cout<<"systematics on Z/gamma ratio --> NLO EWK "<<endl;
      makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		     baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"ewk",2);

      if(not skipTFsystematics){
	
	cout<<"systematics on Z/gamma ratio --> RE 1 "<<endl;
	makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
		       "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"re1",3);
	
	cout<<"systematics on Z/gamma ratio --> FA 1 "<<endl;
	makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
		       "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"fa1",4);
      
	
	cout<<"systematics on Z/gamma ratio --> RE 2 "<<endl;
	makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
		       "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"re2",5);
	
	cout<<"systematics on Z/gamma ratio --> FA 2 "<<endl;
	makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
		       "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"fa2",6);
	
	cout<<"systematics on Z/gamma ratio --> PDF "<<endl;
	makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
		       "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"pdf",7);
	
	cout<<"systematics on Z/gamma ratio --> FP "<<endl;
	makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
		       "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"fpc",8);
      }
    }
    
    cout<<"make Z/W ratio"<<endl;
    makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		   baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		   category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty); 
    
    cout<<"systematics on Z/W ratio --> NLO QCD"<<endl;
    makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		   baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		   category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcd",1);
      
    cout<<"systematics on Z/W ratio --> NLO EWK"<<endl;
    makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		   baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		   category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"ewk",2);

    if(not skipTFsystematics){      

      if(not useNewTheoryUncertainty){
	cout<<"systematics on Z/W ratio --> RE 1"<<endl;
	makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",		   
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"re1",3);
      }
      else{// only for mono-jet
	cout<<"systematics on Z/W ratio --> QCD scale up"<<endl;
	makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",		   
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcd_scaleup",3);
      }

      if(not useNewTheoryUncertainty){
	cout<<"systematics on Z/W ratio --> FA 1"<<endl;
	makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"fa1",4);
      }
      else{
	cout<<"systematics on Z/W ratio --> QCD scale dw"<<endl;
	makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcd_scaledw",4);
      }

      if(not useNewTheoryUncertainty){
	cout<<"systematics on Z/W ratio --> RE 2"<<endl;
	makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"re2",5);
      }
      else{ //only for mono-jet
	cout<<"systematics on Z/W ratio --> NLO-EWK up"<<endl;
	makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nloewk_up",5);
      }
      
      if(not useNewTheoryUncertainty){
	cout<<"systematics on Z/W ratio --> FA 2"<<endl;
	makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"fa2",6);
      }
      else{
	cout<<"systematics on Z/W ratio --> NLO-EWK dw"<<endl;
	makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nloewk_dw",6);
      }

      if(not useNewTheoryUncertainty){
	cout<<"systematics on Z/W ratio --> PDF"<<endl;
	makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"pdf",7);
      }
      else{
	cout<<"systematics on Z/W ratio --> SUD-EWK up"<<endl;
	makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudewk_up",7);	
      }
      
      if(useNewTheoryUncertainty){

	cout<<"systematics on Z/W ratio --> SUD-EWK dw"<<endl;
        makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
                       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
                       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudewk_dw",8);

	cout<<"systematics on Z/W ratio --> MIX QCD-EWK up"<<endl;
        makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
                       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
                       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"ewkqcd_up",9);
   

	cout<<"systematics on Z/W ratio --> MIX QCD-EWK dw"<<endl;
        makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
                       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
                       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"ewkqcd_dw",10);	
      }
    }

    if(addWgamma and category != Category::VBF){
      cout<<"make W/gamma ratio "<<endl;
      makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		      baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false);
	
      cout<<"systematics W/gamma ratio --> NLO QCD "<<endl;
      makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		      baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"qcd",1);
      
      cout<<"systematics W/gamma ratio --> NLO QCD+EWK "<<endl;
      makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		      baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"ewk",2);
      
      if(not skipTFsystematics){
	cout<<"systematics W/gamma ratio --> RE 1 "<<endl;
	makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			"$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"re1",3);
	
	cout<<"systematics W/gamma ratio --> FA 1 "<<endl;
	makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			"$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"fa1",4);
	
	cout<<"systematics W/gamma ratio --> RE 2 "<<endl;
	makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			"$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"re2",5);
	
	cout<<"systematics W/gamma ratio --> FA 2 "<<endl;
	makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			"$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"fa2",6);
	  
	cout<<"systematics W/gamma ratio --> PDF "<<endl;
	makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			"$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"pdf",7);
	
	cout<<"systematics on W/gamma ratio --> FP "<<endl;
	makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
			"$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"fpc",8);
      }
    }  

    if(addTop){
      cout<<"make TOP+MU ratio"<<endl;
      maketopmucorhist(baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       category,observables,observables_2D,lumi,
		       baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       outDir,"",runHiggsInvisible);
      
      cout<<"systematics on TOP+MU ratio --> bUp"<<endl;
      maketopmucorhist(baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       category,observables,observables_2D,lumi,
		       baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       outDir,"btagUp",runHiggsInvisible,"bUp");
      
      
      cout<<"systematics on TOP+MU ratio --> bDw"<<endl;
      maketopmucorhist(baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       category,observables,observables_2D,lumi,
		       baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       outDir,"btagDown",runHiggsInvisible,"bDown");
      
      cout<<"make TOP+EL ratio"<<endl;
      maketopelcorhist(baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       category,observables,observables_2D,lumi,
		       baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       outDir,"",runHiggsInvisible);
      
      
      cout<<"systematics on TOP+EL ratio --> bUp"<<endl;
      maketopelcorhist(baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       category,observables,observables_2D,lumi,
		       baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       outDir,"btagUp",runHiggsInvisible,"bUp");
      
      cout<<"systematics on TOP+EL ratio --> bDw"<<endl;
      maketopelcorhist(baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       category,observables,observables_2D,lumi,
		       baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       outDir,"btagDown",runHiggsInvisible,"bDown");
    }

    // need to add EWK V-jet TFs
    if(category == Category::VBF){
      
      cout<<"make correction histogram for Zmm EWK to Znn EWK"<<endl;      
      makezmmcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/ZJetsToLLEWK/zmmfilter/",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 
      
      cout<<"make correction histogram for Zee EWK to Znn EWK"<<endl;
      makezeecorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/ZJetsToLLEWK/zeefilter/",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 
      
      cout<<"make correction histogram for Wmn EWK to WJets EWK"<<endl;
      makewmncorhist(baseInputTreePath+"/WJetsEWK/sigfilter/",
		     baseInputTreePath+"/WJetsEWK/wmnfilter/",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 
      
      cout<<"make correction histogram for Wen EWK to WJets EWK"<<endl;
      makewencorhist(baseInputTreePath+"/WJetsEWK/sigfilter/",
		     baseInputTreePath+"/WJetsEWK/wenfilter/",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 
      
      
      cout<<"make Z-EWK/W-EWK ratio"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/WJetsEWK/sigfilter/",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,useNewTheoryUncertainty); 
    }
  }

  TFile outfile((outDir+"/templates_"+templateSuffix+".root").c_str(), "RECREATE");  

  if(not skipCorrectionHistograms){
    fillAndSaveCorrQCDHistograms(observables,outfile,outDir,category,addZgamma,addWgamma,addTop,"",addHistoForCutAndCount,useNewTheoryUncertainty);
    if(not observables_2D.empty())
      fillAndSaveCorrQCDHistograms(observables_2D,outfile,outDir,category,addZgamma,addWgamma,addTop,"",addHistoForCutAndCount,useNewTheoryUncertainty);
    if(category == Category::VBF){
      fillAndSaveCorrEWKHistograms(observables,outfile,outDir,category,false,false,"",addHistoForCutAndCount);
      if(not observables_2D.empty())
	fillAndSaveCorrEWKHistograms(observables_2D,outfile,outDir,category,false,false,"",addHistoForCutAndCount);
    }
  }

  // signal region templates
  if(not skipDataAnalysis){
    cout<<"start signal region shapes for signal"<<endl;
    if(not runHiggsInvisible and not runOnlyBackground){
      signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Vector",lumi,doShapeSystematics);
      signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Axial",lumi,doShapeSystematics);
      signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Scalar",lumi,doShapeSystematics);
      signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Pseudoscalar",lumi,doShapeSystematics);
    }
    else if(runHiggsInvisible and not runOnlyBackground){
      signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"110",{5.507E+04,4.434E+03,2.194E+03,1.309E+03},1);
      signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"125",{4.414E+04,3.782E+03,1.373E+03,7.612E+02,1.227E+02},0);
      signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"150",{3.210E+04,3.239E+03,8.154E+02,5.279E+02},1);
      signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"200",{1.812E+04,2.282E+03,3.023E+02,2.054E+02},1);
      signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"300",{9.823E+03,1.256E+03,6.724E+01,4.132E+01},1);
      signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"400",{9.516E+04,7.580E+02,2.163E+01,1.273E+01},1);
      signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"500",{4.538E+03,4.872E+02,8.621E+00,5.256E+00},1);
    }
  }

  if(not skipDataAnalysis and not runOnlySignal){
    cout<<"start signal region data"<<endl;
    sigdatamchist(&outfile,category,observables,observables_2D,lumi,nloSamples,doShapeSystematics,false,false,runHiggsInvisible,applyPostFitWeights,useTheoriestKFactors);
    // gamma + jets
    cout<<"start gamma+jets region data"<<endl;
    gamdatamchist(&outfile,category,observables,observables_2D,nloSamples,lumi,runHiggsInvisible,true,applyPostFitWeights,useTheoriestKFactors);
    // lepton control regions
    cout<<"start zmumu region data"<<endl;
    lepdatamchist(&outfile,Sample::zmm,category,observables,observables_2D,lumi,nloSamples,doShapeSystematics,runHiggsInvisible,false,false,applyPostFitWeights,useTheoriestKFactors); 
    cout<<"start wmunu region data"<<endl;
    lepdatamchist(&outfile,Sample::wmn,category,observables,observables_2D,lumi,nloSamples,doShapeSystematics,runHiggsInvisible,false,false,applyPostFitWeights,useTheoriestKFactors); 
    cout<<"start zee region data"<<endl;
    lepdatamchist(&outfile,Sample::zee,category,observables,observables_2D,lumi,nloSamples,doShapeSystematics,runHiggsInvisible,false,true,applyPostFitWeights,useTheoriestKFactors); 
    cout<<"start wenu region data"<<endl;
    lepdatamchist(&outfile,Sample::wen,category,observables,observables_2D,lumi,nloSamples,doShapeSystematics,runHiggsInvisible,false,true,applyPostFitWeights,useTheoriestKFactors);     
    // top control regions
    if(addTop){
      cout<<"start top+mu region data"<<endl;
      topdatamchist(&outfile,Sample::topmu,category,observables,observables_2D,lumi,nloSamples,makeResonantSelection,doShapeSystematics,runHiggsInvisible,false,applyPostFitWeights);
      cout<<"start Top+el region data"<<endl;
      topdatamchist(&outfile,Sample::topel,category,observables,observables_2D,lumi,nloSamples,makeResonantSelection,doShapeSystematics,runHiggsInvisible,false,applyPostFitWeights);
    }

    if(addQCD){
      cout<<"start QCD region data"<<endl;
      qcddatamchist(&outfile,category,observables,observables_2D,lumi,nloSamples,false,runHiggsInvisible,applyPostFitWeights);
    }

    //add qcd data templates
    TFile* qcdfile_data = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/QCD/templates_2016_12p9.root");
    if(qcdfile_data and (category == Category::monojet or category == Category::monoV)){
      cout<<"Take templates QCD from data"<<endl;
      vector<double> met_bins = selectBinning("met",category);
      TH1F*  qcd_nominal    = new TH1F("qbkghistDD_met","",int(met_bins.size()-1),&met_bins[0]);
      TH1F*  qcd_nominal_up = new TH1F("qbkghistDD_shapeUp_met","",int(met_bins.size()-1),&met_bins[0]);
      TH1F*  qcd_nominal_dw = new TH1F("qbkghistDD_shapeDw_met","",int(met_bins.size()-1),&met_bins[0]);
      
      TH1F* temp = NULL;
      if(category ==  Category::monojet)
	temp = (TH1F*) qcdfile_data->FindObjectAny("hQCD_MonoJ_nominal");
      else if(category ==  Category::monoV)
	temp = (TH1F*) qcdfile_data->FindObjectAny("hQCD_MonoV_nominal");
      
      for(int iBinX = 0; iBinX < qcd_nominal->GetNbinsX(); iBinX++)   
	qcd_nominal->SetBinContent(iBinX+1,temp->GetBinContent(temp->FindBin(qcd_nominal->GetBinCenter(iBinX+1))));
      
      if(category == Category::monojet)
	temp = (TH1F*) qcdfile_data->FindObjectAny("hQCD_MonoJ_AllUp");
      else if(category ==  Category::monoV)
	temp = (TH1F*) qcdfile_data->FindObjectAny("hQCD_MonoV_AllUp");
      
      for(int iBinX = 0; iBinX < qcd_nominal->GetNbinsX(); iBinX++)   
	qcd_nominal_up->SetBinContent(iBinX+1,temp->GetBinContent(temp->FindBin(qcd_nominal->GetBinCenter(iBinX+1))));
      
      if(category == Category::monojet)
	temp = (TH1F*) qcdfile_data->FindObjectAny("hQCD_MonoJ_AllDown");
      else if(category ==  Category::monoV)
	temp = (TH1F*) qcdfile_data->FindObjectAny("hQCD_MonoV_AllDown");
      
      for(int iBinX = 0; iBinX < qcd_nominal->GetNbinsX(); iBinX++)   
	qcd_nominal_dw->SetBinContent(iBinX+1,temp->GetBinContent(temp->FindBin(qcd_nominal->GetBinCenter(iBinX+1))));
      
      outfile.cd();
      outfile.cd("SR");
      qcd_nominal->Write();
      qcd_nominal_up->Write();
      qcd_nominal_dw->Write();      
      outfile.cd();
    }
  }  
  outfile.Close();
}

