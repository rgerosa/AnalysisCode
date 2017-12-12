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
static bool addTauCR              = false;
static bool addWgamma             = true; 
static bool addZgamma             = true;
static bool addZWratio            = true;
static bool skipTFsystematics     = false;
static bool skipDataAnalysis      = false;
static SamplesNLO nloSamples (false,false,false,false);
static bool useTheoriestKFactors  = false;
static bool useNewTheoryUncertainty = true;

void makeTemplates(const bool & doCorrectionHistograms   = false,  // calculate transfer factors and sys
		   const bool & skipCorrectionHistograms = false,  // skip to open and dump transfer factors
		   Category category             = Category::monojet,  // 0 = inclusive mono-j, 1 = exclsuive mono-j, 2 V-tag HP ..
		   const double & lumi           = 35.9, // 
		   const string & outDir         = "", // output dir for template file
		   const string & templateSuffix = "",  // suffix for the output file
		   vector<string> observables    = {"met"}, // 1D histo
		   vector<string> observables_2D = {},  // 2D histo
		   const bool & doShapeSystematics = false, // run all the met, b-tag shape variations
		   const bool & runOnlySignal      = false, // produce a file with only signal templates
		   const bool & runOnlyBackground  = false, // produce a file with only background templates
		   const bool & applyPostFitWeights = false,		   
		   const bool & addHistoForCutAndCount = false,
		   const bool & doVBFOptimization = false, /// to change the VBF selection values from main function call
		   const float & mjjVal = 0,
		   const float & detajjVal = 0,
		   const float & dphijjVal = 3.14) {


  ///// ------
  if(category == Category::VBF or category == Category::VBFrelaxed or category == Category::twojet){
    addZgamma = false;
    addWgamma = false;
  }


  ///// ------
  if(category != Category::VBFrelaxed and category != Category::VBF and doVBFOptimization){
    cerr<<"In order to make the VBFOptimization the category must be VBF or VBFrelaxed --> continue"<<endl;
    return;
  }

  ///// ------
  if(doVBFOptimization){
    if(category == Category::VBFrelaxed){
      detajjrelaxed = detajjVal;
      mjjrelaxed    = mjjVal;
      dphijjrelaxed = dphijjVal;
    }
    else if(category == Category::VBF){
      detajj = detajjVal;
      mjj    = mjjVal;
      dphijj = dphijjrelaxed;
    }
  }

  ///// ------
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

  if (not skipTFsystematics and (category != Category::monojet and category != Category::monoV) and useNewTheoryUncertainty){
    cerr<<"Protection --> new theory uncertainty can be used only for monojet category --> switched to false"<<endl;
    useNewTheoryUncertainty = false;
  }


  ////////// Transfer factors
  if(doCorrectionHistograms){    

    // NLO QCD + NLO EWK
    cout<<"make correction histogram for Zmm to Znn"<<endl;      
    makezmmcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		   baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/",
		   category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors); 

    if(not skipTFsystematics and (category == Category::VBF or category == Category::VBFrelaxed)){      
      if(not useNewTheoryUncertainty){

	cout<<"systematics on Z/Zmm ratio --> FA"<<endl;
	makezmmcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/",		   
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,"fa",1);

	cout<<"systematics on Z/Zmm ratio --> RE"<<endl;
	makezmmcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,"re",2);

	cout<<"systematics on Z/Zmm ratio --> PDF"<<endl;
	makezmmcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,"pdf",3);
      }
    }
    
    // NLO QCD + NLO EWK
    cout<<"make correction histogram for Zee to Znn"<<endl;
    makezeecorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		   baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zeefilter/",
		   category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors); 

    if(not skipTFsystematics and (category == Category::VBF or category == Category::VBFrelaxed)){      
      
      if(not useNewTheoryUncertainty){
	cout<<"systematics on Z/Zee ratio --> FA"<<endl;
	makezeecorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zeefilter/",		   
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,"fa",1);

	cout<<"systematics on Z/Zee ratio --> RE"<<endl;
	makezeecorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zeefilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,"re",2);

	cout<<"systematics on Z/Zee ratio --> PDF"<<endl;
	makezeecorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zeefilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,"pdf",3);
      }
    }    

    // NLO QCD + NLO EWK
    cout<<"make correction histogram for Wmn to WJets"<<endl;
    makewmncorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		   baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/",
		   category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors); 

    if(not skipTFsystematics and (category == Category::VBF or category == Category::VBFrelaxed)){      
      
      if(not useNewTheoryUncertainty){
	cout<<"systematics on W/Wmn ratio --> FA"<<endl;
	makewmncorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/",		   
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,"fa",1);

	cout<<"systematics on W/Wmn ratio --> RE"<<endl;
	makewmncorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,"re",2);

	cout<<"systematics on W/Wmn ratio --> PDF"<<endl;
	makewmncorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,"pdf",3);
      }
    }    
      
    // NLO QCD + NLO EWK
    cout<<"make correction histogram for Wen to WJets"<<endl;
    makewencorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		   baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wenfilter/",
		   category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors); 

    
    if(not skipTFsystematics and (category == Category::VBF or category == Category::VBFrelaxed)){      

      if(not useNewTheoryUncertainty){
	cout<<"systematics on W/Wen ratio --> FA"<<endl;
	makewencorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wenfilter/",		   
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,"fa",1);

	cout<<"systematics on W/Wen ratio --> RE"<<endl;
	makewencorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wenfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,"re",2);

	cout<<"systematics on W/Wen ratio --> PDF"<<endl;
	makewencorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR +"/sigfilter/",
		       baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wenfilter/",
		       category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,"pdf",3);
      }
    }    

    //////////////////////
    /// Add Zgamma TFs ///
    //////////////////////

    if(addZgamma and category != Category::VBF and category != Category::VBFrelaxed and category != Category::twojet){      

      ///////// no re-weight at all
      cout<<"make correction histogram for Gam+jets to Znn"<<endl;
      makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		     baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/", 		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty);
          
      ///////// NLO QCD
      cout<<"systematics on Z/gamma ratio --> NLO QCD "<<endl;
      makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		     baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcd",1);
      
      ////////// NLO QCD + EWK
      cout<<"systematics on Z/gamma ratio --> NLO EWK "<<endl;
      makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		     baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"ewk",2);

      ////////// SYS uncertainties
      if(not skipTFsystematics){	

	if(not useNewTheoryUncertainty){
	  cout<<"systematics on Z/gamma ratio --> RE 1 "<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"re1",3);
	}
	else{
	  cout<<"systematics on Z/gamma ratio --> QCD-Scale up "<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdscale_up",3);
	}
	
	if(not useNewTheoryUncertainty){
	  cout<<"systematics on Z/gamma ratio --> FA 1 "<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"fa1",4);
	}
	else{
	  cout<<"systematics on Z/gamma ratio --> QCD-Scale dw "<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdscale_dw",4);
	}
	
	if(not useNewTheoryUncertainty){
	  cout<<"systematics on Z/gamma ratio --> RE 2 "<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"re2",5);
	}
	else{
	  cout<<"systematics on Z/gamma ratio --> QCD-Shape up "<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdshape_up",5);
	}
	
	if(not useNewTheoryUncertainty){
	  cout<<"systematics on Z/gamma ratio --> FA 2 "<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"fa2",6);
	}
	else{
	  cout<<"systematics on Z/gamma ratio --> QCD shape dw "<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdshape_dw",6);
	}

	if(not useNewTheoryUncertainty){
	  cout<<"systematics on Z/gamma ratio --> PDF"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"pdf",7);
	}
	else{
	  cout<<"systematics on Z/gamma ratio --> QCD-Process up"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdproc_up",7);	
	}
	
	if(not useNewTheoryUncertainty){
	  cout<<"systematics on Z/gamma ratio --> FP "<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"fpc",8);
	}
	else{
	  cout<<"systematics on Z/gamma ratio --> QCD-Process dw"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdproc_dw",8);
	}
	if(useNewTheoryUncertainty){

	  cout<<"systematics on Z/gamma ratio --> NNLO EWK up"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnloewk_up",9);
	  
	  cout<<"systematics on Z/gamma ratio --> NNLO EWK dw"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnloewk_dw",10);	
	  
	  cout<<"systematics on Z/gamma ratio --> NNLO Miss EWK up 1"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnlomiss_up_1",11);	
	  
	  cout<<"systematics on Z/gamma ratio --> NNLO Miss EWK dw 1"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnlomiss_dw_1",12);	
	  
	  cout<<"systematics on Z/gamma ratio --> NNLO Miss EWK up 2"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnlomiss_up_2",13);	
	  
	  cout<<"systematics on Z/gamma ratio --> NNLO Miss EWK dw 2"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnlomiss_dw_2",14);	
	  
	  cout<<"systematics on Z/gamma ratio --> NNLO Sud EWK up 1"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudakov_up_1",15);	
	  
	  cout<<"systematics on Z/gamma ratio --> NNLO Sud EWK dw 1"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudakov_dw_1",16);	
	  
	  cout<<"systematics on Z/gamma ratio --> NNLO Sud EWK up 2"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudakov_up_2",17);	
	  
	  cout<<"systematics on Z/gamma ratio --> NNLO Sud EWK dw 2"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudakov_dw_2",18);	
	  
	  cout<<"systematics on Z/gamma ratio --> NNLO MIX QCD-EWK up"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"mix_up",19);	
	  
	  cout<<"systematics on Z/gamma ratio --> NNLO MIX QCD-EWK dw"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"mix_dw",20);	

	  cout<<"systematics on Z/gamma ratio --> PDF"<<endl;
	  makegamcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			 "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"pdf",21);		  
	  
	}
      }
    }
    
    /////////////////
    /// Z/W ratio ///
    /////////////////
    if(addZWratio){

      string ext ;

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
	
	/////
	if(not useNewTheoryUncertainty){
	  if(decorrelateZWUncertainties){
	    cout<<"systematics on Z/W ratio --> RE Z+jets"<<endl;
	    ext = "rez";
	  }
	  else{
	    cout<<"systematics on Z/W ratio --> RE 1"<<endl;
	    ext = "re1";
	  }
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",		   
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,ext,3);
	}
	else{// only for mono-jet
	  cout<<"systematics on Z/W ratio --> QCD scale up"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",		   
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcd_scaleup",3);
	}
	
	/////
	if(not useNewTheoryUncertainty){
	  if(decorrelateZWUncertainties){
	    cout<<"systematics on Z/W ratio --> FA Z+jets"<<endl;
            ext = "faz";
	  }
          else{
	    cout<<"systematics on Z/W ratio --> FA 1"<<endl;
            ext = "fa1";
	  }
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,ext,4);
	}
	else{
	  cout<<"systematics on Z/W ratio --> QCD scale dw"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcd_scaledw",4);
	}

	//////
	if(not useNewTheoryUncertainty){
	  if(decorrelateZWUncertainties){
	    cout<<"systematics on Z/W ratio --> PDF Z+jets"<<endl;
            ext = "pdfz";
	  }
          else{
	    cout<<"systematics on Z/W ratio --> RE 2"<<endl;
            ext = "re2";
	  }

	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,ext,5);
	}
	else{ //only for mono-jet
	  cout<<"systematics on Z/W ratio --> QCD-Shape up"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdshape_up",5);
	}

	////////
	if(not useNewTheoryUncertainty){
	  if(decorrelateZWUncertainties){
	    cout<<"systematics on Z/W ratio --> RE W+jets"<<endl;
            ext = "rew";
	  }
          else{
	    cout<<"systematics on Z/W ratio --> FA 2"<<endl;
            ext = "fa2";
	  }

	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,ext,6);
	}
	else{
	  cout<<"systematics on Z/W ratio --> QCD-Shape dw"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdshape_dw",6);
	}
	
	///////////
	if(not useNewTheoryUncertainty){
	  if(decorrelateZWUncertainties){
            cout<<"systematics on Z/W ratio --> FA W+jets"<<endl;
            ext = "faw";
          }
          else{
            cout<<"systematics on Z/W ratio --> PDF"<<endl;
            ext = "pdf";
          }

	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,ext,7);
	}
	else{
	  cout<<"systematics on Z/W ratio --> QCD-Proc up"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdproc_up",7);	
	}

	///////
	if(not useNewTheoryUncertainty and decorrelateZWUncertainties){
	  cout<<"systematics on Z/W ratio --> PDF W+jets"<<endl;
	  ext = "pdfw";

          makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
                         baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
                         category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,ext,8);
	}
	else if(useNewTheoryUncertainty){
	  cout<<"systematics on Z/W ratio --> QCD-Proc dw"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdproc_dw",8);
	}

	if(useNewTheoryUncertainty){
	  cout<<"systematics on Z/W ratio --> NNLO EWK up"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnloewk_up",9);
	  
	  cout<<"systematics on Z/W ratio --> NNLO EWK dw"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnloewk_dw",10);	
	  
	  cout<<"systematics on Z/W ratio --> NNLO Miss EWK up 1"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnlomiss_up_1",11);	
	  
	  cout<<"systematics on Z/W ratio --> NNLO Miss EWK dw 1"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnlomiss_dw_1",12);	
	  
	  cout<<"systematics on Z/W ratio --> NNLO Miss EWK up 2"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnlomiss_up_2",13);	
	  
	  cout<<"systematics on Z/W ratio --> NNLO Miss EWK dw 2"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnlomiss_dw_2",14);	
	  
	  cout<<"systematics on Z/W ratio --> NNLO Sud EWK up 1"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudakov_up_1",15);	
	  
	  cout<<"systematics on Z/W ratio --> NNLO Sud EWK dw 1"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudakov_dw_1",16);	
	  
	  cout<<"systematics on Z/W ratio --> NNLO Sud EWK up 2"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudakov_up_2",17);	
	  
	  cout<<"systematics on Z/W ratio --> NNLO Sud EWK dw 2"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudakov_dw_2",18);	
	  
	  cout<<"systematics on Z/W ratio --> NNLO MIX QCD-EWK up"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"mix_up",19);	
	  
	  cout<<"systematics on Z/W ratio --> NNLO MIX QCD-EWK dw"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"mix_dw",20);	
	  
	  cout<<"systematics on Z/W ratio --> PDF"<<endl;
	  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR+"/sigfilter/",
			 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"pdf",21);	
	  
	}
      }
    }
    
    // last block
    if(addWgamma and category != Category::VBF and category != Category::VBFrelaxed and category != Category::twojet){
      cout<<"make W/gamma ratio "<<endl;
      makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		      baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty);
	
      cout<<"systematics W/gamma ratio --> NLO QCD "<<endl;
      makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		      baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcd",1);
      
      cout<<"systematics W/gamma ratio --> NLO QCD+EWK "<<endl;
      makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		      baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"ewk",2);
      
      if(not skipTFsystematics){
	
	if(not useNewTheoryUncertainty){
	  cout<<"systematics W/gamma ratio --> RE 1 "<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"re1",3);
	}
	else{	  
	  cout<<"systematics W/gamma ratio --> QCD scale up "<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdscale_up",3);
	}

	if(not useNewTheoryUncertainty){
	  cout<<"systematics W/gamma ratio --> FA 1 "<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"fa1",4);
	}
	else{
	  cout<<"systematics W/gamma ratio --> QCD scale dw "<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdscale_dw",4);
	}

	if(not useNewTheoryUncertainty){
	  cout<<"systematics W/gamma ratio --> RE 2 "<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"re2",5);
	}
	else{
	  cout<<"systematics W/gamma ratio --> QCD Shape up "<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdshape_up",5);
	}


	if(not useNewTheoryUncertainty){	  
	  cout<<"systematics W/gamma ratio --> FA 2 "<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"fa2",6);
	}
	else{
	  cout<<"systematics W/gamma ratio --> QCD Shape dw "<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdshape_dw",6);
	}

	if(not useNewTheoryUncertainty){
	  cout<<"systematics W/gamma ratio --> PDF "<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"pdf",7);
	}
	else{
	  cout<<"systematics W/gamma ratio --> QCD proc up "<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdproc_up",7);
	}

	if(not useNewTheoryUncertainty){
	  cout<<"systematics on W/gamma ratio --> FP "<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",		   
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"fpc",8);
	}
	else{
	  cout<<"systematics W/gamma ratio --> QCD proc dw "<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"qcdproc_dw",8);
	}

	if(useNewTheoryUncertainty){

	  cout<<"systematics on W/gamma ratio --> NNLO EWK up"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnloewk_up",9);
	  
	  cout<<"systematics on W/gamma ratio --> NNLO EWK dw"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnloewk_dw",10);	
	  
	  cout<<"systematics on W/gamma ratio --> NNLO Miss EWK up 1"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnlomiss_up_1",11);	
	  
	  cout<<"systematics on W/gamma ratio --> NNLO Miss EWK dw 1"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnlomiss_dw_1",12);	
	  
	  cout<<"systematics on W/gamma ratio --> NNLO Miss EWK up 2"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnlomiss_up_2",13);	
	  
	  cout<<"systematics on W/gamma ratio --> NNLO Miss EWK dw 2"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"nnlomiss_dw_2",14);	
	  
	  cout<<"systematics on W/gamma ratio --> NNLO Sud EWK up 1"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudakov_up_1",15);	
	  
	  cout<<"systematics on W/gamma ratio --> NNLO Sud EWK dw 1"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudakov_dw_1",16);	
	  
	  cout<<"systematics on W/gamma ratio --> NNLO Sud EWK up 2"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudakov_up_2",17);	
	  
	  cout<<"systematics on W/gamma ratio --> NNLO Sud EWK dw 2"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"sudakov_dw_2",18);	
	  
	  cout<<"systematics on W/gamma ratio --> NNLO MIX QCD-EWK up"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"mix_up",19);	
	  
	  cout<<"systematics on W/gamma ratio --> NNLO MIX QCD-EWK dw"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"mix_dw",20);	
	  
	  cout<<"systematics on W/gamma ratio --> PDF"<<endl;
	  makewgamcorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
			  baseInputTreePath+"/"+nloSamples.PhotonJetsDIR+"/gamfilter/",
			  "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
			  category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,useTheoriestKFactors,useNewTheoryUncertainty,"pdf",21);		  	  
	 
	}
      }    
    }  

    // need to add EWK V-jet TFs
    if(category == Category::VBF or category == Category::VBFrelaxed){

      cout<<"make correction histogram for Zmm EWK to Znn EWK"<<endl;      
      makezmmcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/ZJetsToLLEWK/zmmfilter/",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 

      
      if(not skipTFsystematics and (category == Category::VBF or category == Category::VBFrelaxed)){      
	if(not useNewTheoryUncertainty){
	  
	  cout<<"systematics on Z/Zmm EWK ratio --> FA"<<endl;
	  makezmmcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
			 baseInputTreePath+"/ZJetsToLLEWK/zmmfilter/",		   
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,"fa",1);
	  
	  cout<<"systematics on Z/Zmm EWK ratio --> RE"<<endl;
	  makezmmcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
			 baseInputTreePath+"/ZJetsToLLEWK/zmmfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,"re",2);
	  
	  cout<<"systematics on Z/Zmm EWK ratio --> PDF"<<endl;
	  makezmmcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
			 baseInputTreePath+"/ZJetsToLLEWK/zmmfilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,"pdf",3);
	}
      }
    
      cout<<"make correction histogram for Zee EWK to Znn EWK"<<endl;
      makezeecorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/ZJetsToLLEWK/zeefilter/",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 

      if(not skipTFsystematics and (category == Category::VBF or category == Category::VBFrelaxed)){      
	if(not useNewTheoryUncertainty){
	  
	  cout<<"systematics on Z/Zee EWK ratio --> FA"<<endl;
	  makezeecorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
			 baseInputTreePath+"/ZJetsToLLEWK/zeefilter/",		   
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,"fa",1);
	  
	  cout<<"systematics on Z/Zee EWK ratio --> RE"<<endl;
	  makezeecorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
			 baseInputTreePath+"/ZJetsToLLEWK/zeefilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,"re",2);
	  
	  cout<<"systematics on Z/Zee EWK ratio --> PDF"<<endl;
	  makezeecorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
			 baseInputTreePath+"/ZJetsToLLEWK/zeefilter/",
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,"pdf",3);
	}
      }
    
      
      cout<<"make correction histogram for Wmn EWK to WJets EWK"<<endl;
      makewmncorhist(baseInputTreePath+"/WJetsEWK/sigfilter/",
		     baseInputTreePath+"/WJetsEWK/wmnfilter/",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 

      if(not skipTFsystematics and (category == Category::VBF or category == Category::VBFrelaxed)){      
	if(not useNewTheoryUncertainty){

	  cout<<"systematics on W/Wmn ratio --> FA"<<endl;
	  makewmncorhist(baseInputTreePath+"/WJetsEWK/sigfilter/",
			 baseInputTreePath+"/WJetsEWK/wmnfilter/",		   
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,"fa",1);
	  
	  cout<<"systematics on W/Wmn ratio --> RE"<<endl;
	  makewmncorhist(baseInputTreePath+"/WJetsEWK/sigfilter/",
			 baseInputTreePath+"/WJetsEWK/wmnfilter/",		   
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,"re",2);
	  
	  cout<<"systematics on W/Wmn ratio --> PDF"<<endl;
	  makewmncorhist(baseInputTreePath+"/WJetsEWK/sigfilter/",
			 baseInputTreePath+"/WJetsEWK/wmnfilter/",		   
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,"pdf",3);
	  
	}
      }
      
      cout<<"make correction histogram for Wen EWK to WJets EWK"<<endl;
      makewencorhist(baseInputTreePath+"/WJetsEWK/sigfilter/",
		     baseInputTreePath+"/WJetsEWK/wenfilter/",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 
      
      
      if(not skipTFsystematics and (category == Category::VBF or category == Category::VBFrelaxed)){      
	if(not useNewTheoryUncertainty){

	  cout<<"systematics on W/Wen ratio --> FA"<<endl;
	  makewencorhist(baseInputTreePath+"/WJetsEWK/sigfilter/",
			 baseInputTreePath+"/WJetsEWK/wenfilter/",		   
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,"fa",1);
	  
	  cout<<"systematics on W/Wen ratio --> RE"<<endl;
	  makewencorhist(baseInputTreePath+"/WJetsEWK/sigfilter/",
			 baseInputTreePath+"/WJetsEWK/wenfilter/",		   
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,"re",2);
	  
	  cout<<"systematics on W/Wen ratio --> PDF"<<endl;
	  makewencorhist(baseInputTreePath+"/WJetsEWK/sigfilter/",
			 baseInputTreePath+"/WJetsEWK/wenfilter/",		   
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,"pdf",3);
	  
	}
      }

      cout<<"make Z-EWK/W-EWK ratio"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/WJetsEWK/sigfilter/",
		     category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,useNewTheoryUncertainty); 


      if(not skipTFsystematics and (category == Category::VBF or category == Category::VBFrelaxed)){ //////////////

	string ext;
	if(not useNewTheoryUncertainty){
	  if(decorrelateZWUncertainties){
	    cout<<"systematics on Z/W ratio EW --> Z+jets RE"<<endl;                                                                                                                                 
	    ext = "rez";
	  }
	  else{
	    cout<<"systematics on Z/W ratio EW --> RE 1"<<endl;                                                                                                                                 
	    ext = "re1";
	  }
	    
	  makezwjcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",                                                                                                                  
			 baseInputTreePath+"/WJetsEWK/sigfilter/",  
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,useNewTheoryUncertainty,ext,3);                           
	  
	  if(decorrelateZWUncertainties){
            cout<<"systematics on Z/W ratio EW --> Z+jets FA"<<endl;
            ext= "faz";
          }
          else{
            cout<<"systematics on Z/W ratio EW --> FA 1"<<endl;
            ext= "fa1";
          }
	  
	  makezwjcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",                                                                                                                  
			 baseInputTreePath+"/WJetsEWK/sigfilter/",  
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,useNewTheoryUncertainty,ext,4);                           

	  if(decorrelateZWUncertainties){
            cout<<"systematics on Z/W ratio EW --> Z+jets PDF"<<endl;
            ext= "pdfz";
          }
          else{
            cout<<"systematics on Z/W ratio EW --> RE 2"<<endl;
            ext= "re2";
          }

	  makezwjcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",                                                                                                                  
			 baseInputTreePath+"/WJetsEWK/sigfilter/",  
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,useNewTheoryUncertainty,ext,5);                           

	  if(decorrelateZWUncertainties){
            cout<<"systematics on Z/W ratio EW --> W+jets RE"<<endl;
            ext= "rew";
          }
          else{
            cout<<"systematics on Z/W ratio EW --> FA 2"<<endl;
            ext= "fa2";
          }
	  	  
	  makezwjcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",                                                                                                                  
			 baseInputTreePath+"/WJetsEWK/sigfilter/",  
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,useNewTheoryUncertainty,ext,6);                           

	  if(decorrelateZWUncertainties){
            cout<<"systematics on Z/W ratio EW --> W+jets FA"<<endl;
            ext= "faw";
          }
          else{
            cout<<"systematics on Z/W ratio EW --> PDF"<<endl;
            ext= "pdf";
          }

	  makezwjcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",                                                                                                                  
			 baseInputTreePath+"/WJetsEWK/sigfilter/",  
			 category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,useNewTheoryUncertainty,ext,7);                           
	  
	  if(decorrelateZWUncertainties){
            cout<<"systematics on Z/W ratio EW --> W+jets PDF"<<endl;
            ext= "pdfw";
	    makezwjcorhist(baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",                                                                                                                  
			   baseInputTreePath+"/WJetsEWK/sigfilter/",  
			   category,nloSamples,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,useTheoriestKFactors,useNewTheoryUncertainty,ext,8);                           
	    
          }	  
	}                                                                                                                                                                                             
      }
    }
  }

  TFile outfile ((outDir+"/templates_"+templateSuffix+".root").c_str(), "RECREATE");  

  if(not skipCorrectionHistograms){
    fillAndSaveCorrQCDHistograms(observables,outfile,outDir,category,addZWratio,addZgamma,addWgamma,"",addHistoForCutAndCount,useNewTheoryUncertainty);
    if(not observables_2D.empty())
      fillAndSaveCorrQCDHistograms(observables_2D,outfile,outDir,category,addZWratio,addZgamma,addWgamma,"",addHistoForCutAndCount,useNewTheoryUncertainty);
    if(category == Category::VBF or category == Category::VBFrelaxed){
      fillAndSaveCorrQCDHistograms(observables,outfile,outDir,category,addZWratio,addZgamma,addWgamma,"ewk",addHistoForCutAndCount,useNewTheoryUncertainty); 
      //fillAndSaveCorrEWKHistograms(observables,outfile,outDir,category,addZWratio,"",addHistoForCutAndCount);
      if(not observables_2D.empty())
	fillAndSaveCorrQCDHistograms(observables_2D,outfile,outDir,category,addZWratio,addZgamma,addWgamma,"ewk",addHistoForCutAndCount,useNewTheoryUncertainty);
      //fillAndSaveCorrEWKHistograms(observables_2D,outfile,outDir,category,addZWratio,"",addHistoForCutAndCount);
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
    if(category != Category::VBF and category != Category::VBFrelaxed){
      cout<<"start gamma+jets region data"<<endl;
      gamdatamchist(&outfile,category,observables,observables_2D,nloSamples,lumi,runHiggsInvisible,false,applyPostFitWeights,useTheoriestKFactors);
    }
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

    if(addTauCR){
      cout<<"start Tau region data"<<endl;
      taudatamchist(&outfile,category,observables,observables_2D,lumi,nloSamples,false,runHiggsInvisible,applyPostFitWeights);
    }

    //add qcd data templates
    if(category == Category::monojet or category == Category::monoV){

      TFile* qcdfile_data = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/QCD/MonojetTemplates/templates_2016_12p9.root");

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
      

      /// Scaling for lumi
      qcd_nominal->Scale(lumi/12.9);
      qcd_nominal_up->Scale(lumi/12.9);
      qcd_nominal_dw->Scale(lumi/12.9);
      
      outfile.cd();
      outfile.cd("SR");
      qcd_nominal->Write();
      qcd_nominal_up->Write();
      qcd_nominal_dw->Write();      
      outfile.cd();
    }
    else if(category == Category::VBFrelaxed){

      TFile* qcdfile_data = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/QCD/VBFTemplates/templates_QCD_DD_v2.root");
      cout<<"Take templates QCD from data"<<endl;      
      TH1F* qcd_nominal = (TH1F*) qcdfile_data->Get("template_QCD_SR_fromDD");
      // set the bin errors to reflect the bin-by-bin variations in the TFs
      for(int iBin = 0; iBin < qcd_nominal->GetNbinsX(); iBin++){
	TH1F* qcd_temp_up = (TH1F*) qcdfile_data->Get(Form("BinByBin/template_QCD_SR_fromDD_bin_%d_statUp",iBin));
	TH1F* qcd_temp_dw = (TH1F*) qcdfile_data->Get(Form("BinByBin/template_QCD_SR_fromDD_bin_%d_statDown",iBin));
	qcd_nominal->SetBinError(iBin+1,(qcd_temp_up->GetBinContent(iBin+1)-qcd_temp_dw->GetBinContent(iBin+1))/2);
      }

      outfile.cd();
      outfile.cd("SR");
      qcd_nominal->Write("qbkghistDD_mjj");
      outfile.cd();      
    }
    else if(category == Category::VBF){
      TFile* qcdfile_data = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/QCD/VBFTemplates/templates_QCD_DD_v2.root");
      cout<<"Take templates QCD from data"<<endl;      

      vector<double> met_bins = selectBinning("met_onebin",category);
      TH1F* qbkghistDD_met_onebin = new TH1F("qbkghistDD_met_onebin","",met_bins.size()-1,&met_bins[0]);
      
      TH1F* qcd_nominal = (TH1F*) qcdfile_data->Get("template_QCD_SR_fromDD"); 
      qbkghistDD_met_onebin->SetBinContent(1,qcd_nominal->Integral(qcd_nominal->FindBin(1300),qcd_nominal->GetNbinsX())); // special trick for the cut and count case

      // set the bin errors to reflect the bin-by-bin variations in the TFs
      float error = 0;
      for(int iBin = qcd_nominal->FindBin(1300); iBin < qcd_nominal->GetNbinsX(); iBin++){
	TH1F* qcd_temp_up = (TH1F*) qcdfile_data->Get(Form("BinByBin/template_QCD_SR_fromDD_bin_%d_statUp",iBin-1));
	TH1F* qcd_temp_dw = (TH1F*) qcdfile_data->Get(Form("BinByBin/template_QCD_SR_fromDD_bin_%d_statDown",iBin-1));
	error += pow((qcd_temp_up->GetBinContent(iBin)-qcd_temp_dw->GetBinContent(iBin))/2,2);
      }      
      qbkghistDD_met_onebin->SetBinError(1,sqrt(error));
      outfile.cd();
      outfile.cd("SR");
      qbkghistDD_met_onebin->Write("qbkghistDD_met_onebin");
      outfile.cd();      
    }      
  }  
  outfile.Close();
}

