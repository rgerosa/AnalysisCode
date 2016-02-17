#include "makehist.h"
#include "makeCorrHistograms.C"
#include "makeDataHistograms.C"

using namespace std;

// Run the final analysis:
// 1) Store all corrections templates from input files (complient to combine)
// 2) Make data and expected yields templates for all the other processes

void makeTemplates(bool doCorrectionHistograms = false, 
		   bool skipCorrectionHistograms = false, 
		   int category  = 0, 
		   double lumi   = 2.24, 
		   string outDir = "", 
		   string templateSuffix = "", 
		   vector<string> observables    = {"met"}, 
		   vector<string> observables_2D = {}, 
		   bool applyQGLReweight      = false,
		   bool doShapeSystematics    = false,
		   bool makeResonantSelection = false,
		   string ext ="") {

  system(("mkdir -p "+outDir).c_str());

  vector<signalSample> signalMassPoint;
  signalMassPoint.push_back(signalSample("Vector","500","1"));
  signalMassPoint.push_back(signalSample("Vector","500","10"));
  signalMassPoint.push_back(signalSample("Vector","1000","50"));
  signalMassPoint.push_back(signalSample("Vector","2000","10"));

  signalMassPoint.push_back(signalSample("Axial","500","1"));
  signalMassPoint.push_back(signalSample("Axial","500","10"));
  signalMassPoint.push_back(signalSample("Axial","1000","10"));
  signalMassPoint.push_back(signalSample("Axial","2000","1"));

  signalMassPoint.push_back(signalSample("Scalar","100","50"));
  signalMassPoint.push_back(signalSample("Scalar","1000","10"));
  signalMassPoint.push_back(signalSample("Scalar","2000","10"));

  signalMassPoint.push_back(signalSample("Pseudoscalar","100","1"));
  signalMassPoint.push_back(signalSample("Pseudoscalar","200","1"));
  signalMassPoint.push_back(signalSample("Pseudoscalar","500","1"));
  signalMassPoint.push_back(signalSample("Pseudoscalar","1000","1"));

  if(doCorrectionHistograms){

    cout<<"make correction histogram for Zmm to Znn"<<endl;
    // make central values
    makezmmcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
    		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/zmmfilter/zmm_tree_DYJetsToLL_M-50.root",
    		   category,observables,lumi,applyQGLReweight,outDir,"",ext); 

    cout<<"make correction histogram for Zee to Znn"<<endl;
    makezeecorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/zeefilter/zee_tree_DYJetsToLL_M-50.root",
    		   category,observables,lumi,applyQGLReweight,outDir,"",ext); 

    cout<<"make correction histogram for Wmn to WJets"<<endl;
    makewmncorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/wmnfilter/wmn_tree_WJetsToLNu.root",
    		   category,observables,lumi,applyQGLReweight,outDir,"",ext); 

    cout<<"make correction histogram for Wen to WJets"<<endl;
    makewencorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/wenfilter/wen_tree_WJetsToLNu.root",
    		   category,observables,lumi,applyQGLReweight,outDir,"",ext); 
   
    cout<<"make correction histogram for Gam+jets to Znn"<<endl;
    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root", 		   
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,observables,lumi,applyQGLReweight,outDir,"",ext);


    cout<<"systematics on Z/gamma ratio --> NLO QCD "<<endl;
    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",		   
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","qcd"+ext,1);

    cout<<"systematics on Z/gamma ratio --> NLO EWK "<<endl;
    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",		   
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","ewk"+ext,2);
    
    cout<<"systematics on Z/gamma ratio --> RE 1 "<<endl;
    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",		   
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","re1"+ext,3);

    cout<<"systematics on Z/gamma ratio --> FA 1 "<<endl;
    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",		   
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","fa1"+ext,4);


    cout<<"systematics on Z/gamma ratio --> RE 2 "<<endl;
    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",		   
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","re2"+ext,5);

    cout<<"systematics on Z/gamma ratio --> FA 2 "<<endl;
    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",		   
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","fa2"+ext,6);

    cout<<"systematics on Z/gamma ratio --> PDF "<<endl;
    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",		   
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","pdf"+ext,7);

    cout<<"systematics on Z/gamma ratio --> FP "<<endl;
    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",		   
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","fpc"+ext,8);


    cout<<"make Z/W ratio"<<endl;
    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
    		   category,observables,lumi,applyQGLReweight,outDir,"",ext); 

    cout<<"systematics on Z/W ratio --> NLO QCD"<<endl;
    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","qcd"+ext,1);

    cout<<"systematics on Z/W ratio --> NLO EWK"<<endl;
    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","ewk"+ext,2);

    cout<<"systematics on Z/W ratio --> RE 1"<<endl;
    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",		   
		   category,observables,lumi,applyQGLReweight,outDir,"","re1"+ext,3);


    cout<<"systematics on Z/W ratio --> FA 1"<<endl;
    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","fa1"+ext,4);


    cout<<"systematics on Z/W ratio --> RE 2"<<endl;
    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","re2"+ext,5);


    cout<<"systematics on Z/W ratio --> FA 2"<<endl;
    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","fa2"+ext,6);

    cout<<"systematics on Z/W ratio --> PDF"<<endl;
    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   category,observables,lumi,applyQGLReweight,outDir,"","pdf"+ext,7);

    cout<<"make TOP+MU ratio"<<endl;
    maketopmucorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top_amc.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top_amc.root",
		     category,observables,lumi,
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root",
		     applyQGLReweight,outDir,"",ext);

    cout<<"systematics on TOP+MU ratio --> bUp"<<endl;
    maketopmucorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top_amc.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top_amc.root",
		     category,observables,lumi,
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root",
		     applyQGLReweight,outDir,"btagUp",ext+"bUp");


    cout<<"systematics on TOP+MU ratio --> bDw"<<endl;
    maketopmucorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top_amc.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top_amc.root",
		     category,observables,lumi,
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root",
		     applyQGLReweight,outDir,"btagDown",ext+"bDown");

    cout<<"make TOP+EL ratio"<<endl;
    maketopelcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top_amc.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top_amc.root",
		     category,observables,lumi,
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root",
		     applyQGLReweight,outDir,"",ext);


    cout<<"systematics on TOP+MU ratio --> bUp"<<endl;
    maketopelcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top_amc.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top_amc.root",
		     category,observables,lumi,
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root",
		     applyQGLReweight,outDir,"btagUp",ext+"bUp");

    cout<<"systematics on TOP+MU ratio --> bDw"<<endl;
    maketopelcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top_amc.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top_amc.root",
		     category,observables,lumi,
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root",
		     applyQGLReweight,outDir,"btagDown",ext+"bDown");

    /*
    if(category == 2 or category == 3){

      cout<<"make sideband ratio for Zvv "<<std::endl;      
      makesidebandcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
			  "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
			  category,category+2,observables,lumi,outDir,ext+"Z");

      cout<<"make sideband ratio for W+jets "<<std::endl;      
      makesidebandcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
			  "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
			  category,category+2,observables,lumi,outDir,ext+"W");
    }
    */

  }

  TFile outfile((outDir+"/templates_"+templateSuffix+".root").c_str(), "RECREATE");

  if(not skipCorrectionHistograms){

    // take correction files --> central value
    cout<<"Re-open file for correction histo"<<endl;
    TFile* zmmcorfile = TFile::Open((outDir+"/zmmcor"+ext+".root").c_str());
    TFile* zeecorfile = TFile::Open((outDir+"/zeecor"+ext+".root").c_str());
    TFile* wmncorfile = TFile::Open((outDir+"/wmncor"+ext+".root").c_str());
    TFile* wencorfile = TFile::Open((outDir+"/wencor"+ext+".root").c_str());
    TFile* zwjcorfile = TFile::Open((outDir+"/zwjcor"+ext+".root").c_str());
    TFile* gamcorfile = TFile::Open((outDir+"/gamcor"+ext+".root").c_str());
    TFile* topmucorfile = TFile::Open((outDir+"/topmucor"+ext+".root").c_str());
    TFile* topelcorfile = TFile::Open((outDir+"/topelcor"+ext+".root").c_str());

    /*
    TFile* sidebandfileZ = NULL;
    TFile* sidebandfileW = NULL;
    if(category == 2 or category == 3){
      sidebandfileZ = TFile::Open((outDir+"/sidebandcor"+ext+"Z.root").c_str());
      sidebandfileW = TFile::Open((outDir+"/sidebandcor"+ext+"W.root").c_str());
    }
    */

    // QCD, EWK, factm re and footprint on Z/gamma
    TFile* gamcorqcdfile = TFile::Open((outDir+"/gamcorqcd"+ext+".root").c_str());
    TFile* gamcorewkfile = TFile::Open((outDir+"/gamcorewk"+ext+".root").c_str());
    TFile* gamcorre1file = TFile::Open((outDir+"/gamcorre1"+ext+".root").c_str());
    TFile* gamcorfa1file = TFile::Open((outDir+"/gamcorfa1"+ext+".root").c_str());
    TFile* gamcorre2file = TFile::Open((outDir+"/gamcorre2"+ext+".root").c_str());
    TFile* gamcorfa2file = TFile::Open((outDir+"/gamcorfa2"+ext+".root").c_str());
    TFile* gamcorpdffile = TFile::Open((outDir+"/gamcorpdf"+ext+".root").c_str());
    TFile* gamcorfpcfile = TFile::Open((outDir+"/gamcorfpc"+ext+".root").c_str());

    // QCD, EWK, factm re and footprint on Z/W
    TFile* zwjcorqcdfile = TFile::Open((outDir+"/zwjcorqcd"+ext+".root").c_str());
    TFile* zwjcorewkfile = TFile::Open((outDir+"/zwjcorewk"+ext+".root").c_str());
    TFile* zwjcorre1file = TFile::Open((outDir+"/zwjcorre1"+ext+".root").c_str());
    TFile* zwjcorfa1file = TFile::Open((outDir+"/zwjcorfa1"+ext+".root").c_str());
    TFile* zwjcorre2file = TFile::Open((outDir+"/zwjcorre2"+ext+".root").c_str());
    TFile* zwjcorfa2file = TFile::Open((outDir+"/zwjcorfa2"+ext+".root").c_str());
    TFile* zwjcorpdffile = TFile::Open((outDir+"/zwjcorpdf"+ext+".root").c_str());
    
    // Top btag up and down
    TFile* topmucorbupfile   = TFile::Open((outDir+"/topmucor"+ext+"bUp.root").c_str());
    TFile* topmucorbdownfile = TFile::Open((outDir+"/topmucor"+ext+"bDown.root").c_str());
    TFile* topelcorbupfile   = TFile::Open((outDir+"/topelcor"+ext+"bUp.root").c_str());
    TFile* topelcorbdownfile = TFile::Open((outDir+"/topelcor"+ext+"bDown.root").c_str());
    
    // get histograms  
    vector<TH1*> zmmcorhist;
    vector<TH1*> zeecorhist;
    vector<TH1*> wmncorhist;
    vector<TH1*> wencorhist;
    vector<TH1*> zwjcorhist;
    vector<TH1*> gamcorhist;
    vector<TH1*> topmucorhist;
    vector<TH1*> topelcorhist;
    
    vector<TH1*> gamcorewkhist;
    vector<TH1*> gamcorqcdhist;
    vector<TH1*> gamcorre1hist;
    vector<TH1*> gamcorfa1hist;
    vector<TH1*> gamcorre2hist;
    vector<TH1*> gamcorfa2hist;
    vector<TH1*> gamcorpdfhist;
    vector<TH1*> gamcorfpchist;
    
    vector<TH1*> zwjcorewkhist;
    vector<TH1*> zwjcorqcdhist;
    vector<TH1*> zwjcorre1hist;
    vector<TH1*> zwjcorre2hist;
    vector<TH1*> zwjcorfa1hist;
    vector<TH1*> zwjcorfa2hist;
    vector<TH1*> zwjcorpdfhist;
    
    vector<TH1F*> sidebandZhist;
    vector<TH1F*> sidebandWhist;
    
    vector<TH1*> topmucorbuphist;
    vector<TH1*> topmucorbdownhist;
    vector<TH1*> topelcorbuphist;
    vector<TH1*> topelcorbdownhist;

    // output file
    for(auto obs : observables){

      cout<<"Get histograms for observable "<<obs<<endl;
      
      zmmcorhist.push_back( (TH1*)zmmcorfile->Get(("zmmcor"+ext+"hist_"+obs).c_str()));    
      zeecorhist.push_back( (TH1*)zeecorfile->Get(("zeecor"+ext+"hist_"+obs).c_str()));    
      wmncorhist.push_back( (TH1*)wmncorfile->Get(("wmncor"+ext+"hist_"+obs).c_str()));    
      wencorhist.push_back( (TH1*)wencorfile->Get(("wencor"+ext+"hist_"+obs).c_str()));    
      zwjcorhist.push_back( (TH1*)zwjcorfile->Get(("zwjcor"+ext+"hist_"+obs).c_str()));    
      gamcorhist.push_back( (TH1*)gamcorfile->Get(("gamcor"+ext+"hist_"+obs).c_str()));    
      topmucorhist.push_back( (TH1*)topmucorfile->Get(("topmucor"+ext+"hist_"+obs).c_str()));    
      topelcorhist.push_back( (TH1*)topelcorfile->Get(("topelcor"+ext+"hist_"+obs).c_str()));    

      /*
      if(category == 2 or category == 3){
	sidebandZhist.push_back((TH1F*) sidebandfileZ->Get(("sidebandcor"+ext+"Zhist_"+obs).c_str()));
	sidebandWhist.push_back((TH1F*) sidebandfileW->Get(("sidebandcor"+ext+"Whist_"+obs).c_str()));
      }
      */

      // get histograms Z/gamma
      cout<<"Make Z/gamma sys histograms"<<endl;
      gamcorewkhist.push_back( (TH1*)gamcorewkfile->Get(("gamcor"+ext+"ewkhist_"+obs).c_str()));    
      gamcorqcdhist.push_back( (TH1*)gamcorqcdfile->Get(("gamcor"+ext+"qcdhist_"+obs).c_str()));    
      gamcorre1hist.push_back( (TH1*)gamcorre1file->Get(("gamcor"+ext+"re1hist_"+obs).c_str()));    
      gamcorfa1hist.push_back( (TH1*)gamcorfa1file->Get(("gamcor"+ext+"fa1hist_"+obs).c_str()));    
      gamcorre2hist.push_back( (TH1*)gamcorre2file->Get(("gamcor"+ext+"re2hist_"+obs).c_str()));    
      gamcorfa2hist.push_back( (TH1*)gamcorfa2file->Get(("gamcor"+ext+"fa2hist_"+obs).c_str()));    
      gamcorpdfhist.push_back( (TH1*)gamcorpdffile->Get(("gamcor"+ext+"pdfhist_"+obs).c_str()));    
      gamcorfpchist.push_back( (TH1*)gamcorfpcfile->Get(("gamcor"+ext+"fpchist_"+obs).c_str()));    
      
      // uncertainty histogram for combine
      TH1* gamuncewkhist = (TH1*)gamcorewkhist.back()->Clone(("gamuncewk"+ext+"hist_"+obs).c_str());    
      gamuncewkhist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncewkhist->GetNbinsX(); i++) 
      gamuncewkhist->SetBinContent(i, fabs(gamuncewkhist->GetBinContent(i)-1.0));
      gamuncewkhist->SetName(("ZG_EWK_"+obs).c_str());
      
      TH1* gamuncre1hist = (TH1*)gamcorre1hist.back()->Clone(("gamuncre1"+ext+"hist_"+obs).c_str());    
      gamuncre1hist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncre1hist->GetNbinsX(); i++) 
	gamuncre1hist->SetBinContent(i, fabs(gamuncre1hist->GetBinContent(i)-1.0));
      gamuncre1hist->SetName(("ZG_RenScale1_"+obs).c_str());
      
      TH1* gamuncfa1hist = (TH1*)gamcorfa1hist.back()->Clone(("gamuncfa1"+ext+"hist_"+obs).c_str());    
      gamuncfa1hist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncfa1hist->GetNbinsX(); i++) 
      gamuncfa1hist->SetBinContent(i, fabs(gamuncfa1hist->GetBinContent(i)-1.0));
      gamuncfa1hist->SetName(("ZG_FactScale1_"+obs).c_str());
    
      TH1* gamuncre2hist = (TH1*)gamcorre2hist.back()->Clone(("gamuncre2"+ext+"hist_"+obs).c_str());    
      gamuncre2hist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncre2hist->GetNbinsX(); i++) 
	gamuncre2hist->SetBinContent(i, fabs(gamuncre2hist->GetBinContent(i)-1.0));
      gamuncre2hist->SetName(("ZG_RenScale2_"+obs).c_str());
      
      TH1* gamuncfa2hist = (TH1*)gamcorfa2hist.back()->Clone(("gamuncfa2"+ext+"hist_"+obs).c_str());    
      gamuncfa2hist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncfa2hist->GetNbinsX(); i++) 
	gamuncfa2hist->SetBinContent(i, fabs(gamuncfa2hist->GetBinContent(i)-1.0));
      gamuncfa2hist->SetName(("ZG_FactScale2_"+obs).c_str());
      
      TH1* gamuncpdfhist = (TH1*)gamcorpdfhist.back()->Clone(("gamuncpdf"+ext+"hist_"+obs).c_str());    
      gamuncpdfhist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncpdfhist->GetNbinsX(); i++) 
	gamuncpdfhist->SetBinContent(i, fabs(gamuncpdfhist->GetBinContent(i)-1.0));
      gamuncpdfhist->SetName(("ZG_PDF_"+obs).c_str());
    
      TH1* gamuncfpchist = (TH1*)gamcorfpchist.back()->Clone(("gamuncfpc"+ext+"hist_"+obs).c_str());    
      gamuncfpchist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncfpchist->GetNbinsX(); i++) 
	gamuncfpchist->SetBinContent(i, fabs(gamuncfpchist->GetBinContent(i)-1.0));
      gamuncfpchist->SetName(("ZG_Footprint_"+obs).c_str());
      
      // Same thing for Z/W ratio
      cout<<"Make Z/W sys histograms"<<endl;
      zwjcorewkhist.push_back( (TH1*)zwjcorewkfile->Get(("zwjcorewk"+ext+"hist_"+obs).c_str()));    
      zwjcorqcdhist.push_back( (TH1*)zwjcorqcdfile->Get(("zwjcorqcd"+ext+"hist_"+obs).c_str()));    
      zwjcorre1hist.push_back( (TH1*)zwjcorre1file->Get(("zwjcorre1"+ext+"hist_"+obs).c_str()));    
      zwjcorfa1hist.push_back( (TH1*)zwjcorfa1file->Get(("zwjcorfa1"+ext+"hist_"+obs).c_str()));    
      zwjcorre2hist.push_back( (TH1*)zwjcorre2file->Get(("zwjcorre2"+ext+"hist_"+obs).c_str()));    
      zwjcorfa2hist.push_back( (TH1*)zwjcorfa2file->Get(("zwjcorfa2"+ext+"hist_"+obs).c_str()));    
      zwjcorpdfhist.push_back( (TH1*)zwjcorpdffile->Get(("zwjcorpdf"+ext+"hist_"+obs).c_str()));    
      
      TH1* zwjuncewkhist = (TH1*)zwjcorewkhist.back()->Clone(("zwjuncewk"+ext+"hist_"+obs).c_str());    
      zwjuncewkhist->Divide(zwjcorqcdhist.back());
      for (int i = 1; i <= zwjuncewkhist->GetNbinsX(); i++) 
	zwjuncewkhist->SetBinContent(i, fabs(zwjuncewkhist->GetBinContent(i)-1.0));
      zwjuncewkhist->SetName(("ZW_EWK_"+obs).c_str());
      
      TH1* zwjuncre1hist = (TH1*)zwjcorre1hist.back()->Clone(("zwjuncre1"+ext+"hist_"+obs).c_str());
      zwjuncre1hist->Divide(zwjcorqcdhist.back());
      for (int i = 1; i <= zwjuncre1hist->GetNbinsX(); i++) 
	zwjuncre1hist->SetBinContent(i, fabs(zwjuncre1hist->GetBinContent(i)-1.0));
      zwjuncre1hist->SetName(("ZW_RenScale1_"+obs).c_str());

      TH1* zwjuncfa1hist = (TH1*)zwjcorfa1hist.back()->Clone(("zwjuncfa1"+ext+"hist_"+obs).c_str());
      zwjuncfa1hist->Divide(zwjcorqcdhist.back());
      for (int i = 1; i <= zwjuncfa1hist->GetNbinsX(); i++) 
	zwjuncfa1hist->SetBinContent(i, fabs(zwjuncfa1hist->GetBinContent(i)-1.0));
      zwjuncfa1hist->SetName(("ZW_FactScale1_"+obs).c_str());
    
      TH1* zwjuncre2hist = (TH1*)zwjcorre2hist.back()->Clone(("zwjuncre2"+ext+"hist_"+obs).c_str());
      zwjuncre2hist->Divide(zwjcorqcdhist.back());
      for (int i = 1; i <= zwjuncre2hist->GetNbinsX(); i++) 
	zwjuncre2hist->SetBinContent(i, fabs(zwjuncre2hist->GetBinContent(i)-1.0));
      zwjuncre2hist->SetName(("ZW_RenScale2_"+obs).c_str());
      
      TH1* zwjuncfa2hist = (TH1*)zwjcorfa2hist.back()->Clone(("zwjuncfa2"+ext+"hist_"+obs).c_str());
      zwjuncfa2hist->Divide(zwjcorqcdhist.back());
      for (int i = 1; i <= zwjuncfa2hist->GetNbinsX(); i++) 
	zwjuncfa2hist->SetBinContent(i, fabs(zwjuncfa2hist->GetBinContent(i)-1.0));
      zwjuncfa2hist->SetName(("ZW_FactScale2_"+obs).c_str());
    
      TH1* zwjuncpdfhist = (TH1*)zwjcorpdfhist.back()->Clone(("zwjuncpdf"+ext+"hist_"+obs).c_str());
      zwjuncpdfhist->Divide(zwjcorqcdhist.back());
      for (int i = 1; i <= zwjuncpdfhist->GetNbinsX(); i++) 
	zwjuncpdfhist->SetBinContent(i, fabs(zwjuncpdfhist->GetBinContent(i)-1.0));
      zwjuncpdfhist->SetName(("ZW_PDF_"+obs).c_str());

      // make b-tagging top
      cout<<"Make top sys histograms"<<endl;
      topmucorbuphist.push_back( (TH1*) topmucorbupfile->Get(("topmucor"+ext+"bUphist_"+obs).c_str()));    
      topmucorbdownhist.push_back( (TH1*) topmucorbdownfile->Get(("topmucor"+ext+"bDownhist_"+obs).c_str()));    
      topelcorbuphist.push_back( (TH1*) topelcorbupfile->Get(("topelcor"+ext+"bUphist_"+obs).c_str()));    
      topelcorbdownhist.push_back( (TH1*) topelcorbdownfile->Get(("topelcor"+ext+"bDownhist_"+obs).c_str()));    
      
      // make symmetrization
      TH1* topmucorbuphist_tmp = (TH1*) topmucorbuphist.back()->Clone(("topmucorbup_tmp"+ext+"hist_"+obs).c_str());    
      topmucorbuphist_tmp->Divide(topmucorhist.back());
      TH1* topmucorbdownhist_tmp = (TH1*) topmucorbdownhist.back()->Clone(("topmucorbdown_tmp"+ext+"hist_"+obs).c_str());    
      topmucorbdownhist_tmp->Divide(topmucorhist.back());
      
      TH1* topmucoruncbhist = (TH1*) topmucorhist.back()->Clone(("topmucoruncbhist"+ext+"hist_"+obs).c_str());    
      for (int i = 1; i <= topmucoruncbhist->GetNbinsX(); i++) {      
	topmucoruncbhist->SetBinContent(i, fabs(fabs(topmucorbuphist_tmp->GetBinContent(i)+topmucorbdownhist_tmp->GetBinContent(i))/2-1.0));
      }
      topmucoruncbhist->SetName(("TOP_MU_B_"+obs).c_str());
    
      // make symmetrization
      TH1* topelcorbuphist_tmp = (TH1*) topelcorbuphist.back()->Clone(("topelcorbup_tmp"+ext+"hist_"+obs).c_str());    
      topelcorbuphist_tmp->Divide(topelcorhist.back());
      TH1* topelcorbdownhist_tmp = (TH1*) topelcorbdownhist.back()->Clone(("topelcorbdown_tmp"+ext+"hist_"+obs).c_str());    
      topelcorbdownhist_tmp->Divide(topelcorhist.back());
      
      TH1* topelcoruncbhist = (TH1*) topelcorhist.back()->Clone(("topelcoruncbhist"+ext+"hist_"+obs).c_str());    
      for (int i = 1; i <= topelcoruncbhist->GetNbinsX(); i++) {      
	topelcoruncbhist->SetBinContent(i, fabs(fabs(topelcorbuphist_tmp->GetBinContent(i)+topelcorbdownhist_tmp->GetBinContent(i))/2-1.0));
      }
      topelcoruncbhist->SetName(("TOP_EL_B_"+obs).c_str());
      
      outfile.cd();
      
      cout<<"Save transfer factor"<<endl;
      zmmcorhist.back()->Write();
      zeecorhist.back()->Write();
      wmncorhist.back()->Write();
      wencorhist.back()->Write();
      zwjcorhist.back()->Write();
      gamcorhist.back()->Write();
      topmucorhist.back()->Write();
      topelcorhist.back()->Write();

      /*      
      if(category == 2 or category == 3){
	sidebandZhist.back()->Write();
	sidebandWhist.back()->Write();
      }
      */

      gamcorqcdhist.back()->Write();
      gamcorewkhist.back()->Write();
      gamcorre1hist.back()->Write();
      gamcorfa1hist.back()->Write();
      gamcorre2hist.back()->Write();
      gamcorfa2hist.back()->Write();
      gamcorpdfhist.back()->Write();
      gamcorfpchist.back()->Write();
      
      gamuncewkhist->Write();
      gamuncre1hist->Write();
      gamuncfa1hist->Write();
      gamuncre2hist->Write();
      gamuncfa2hist->Write();
      gamuncpdfhist->Write();
      gamuncfpchist->Write();
      
      zwjcorqcdhist.back()->Write();
      zwjcorewkhist.back()->Write();
      zwjcorre1hist.back()->Write();
      zwjcorfa1hist.back()->Write();
      zwjcorre2hist.back()->Write();
      zwjcorfa2hist.back()->Write();
      zwjcorpdfhist.back()->Write();
      
      zwjuncewkhist->Write();
      zwjuncre1hist->Write();
      zwjuncfa1hist->Write();
      zwjuncre2hist->Write();
      zwjuncfa2hist->Write();
      zwjuncpdfhist->Write();
      
      topmucoruncbhist->Write();
      topelcoruncbhist->Write();      
    }
  }

  // signal region templates
  cout<<"start signal region shapes for signal"<<endl;
  signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Vector",lumi,doShapeSystematics);
  signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Axial",lumi,doShapeSystematics);
  signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Scalar",lumi,doShapeSystematics);
  signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Pseudoscalar",lumi,doShapeSystematics);
  cout<<"start signal region data"<<endl;
  sigdatamchist(&outfile,category,observables,observables_2D,lumi,applyQGLReweight,doShapeSystematics,true,false);
  // gamma + jets
  cout<<"start gamma+jets region data"<<endl;
  gamdatamchist(&outfile,category,observables,observables_2D,lumi,applyQGLReweight);
  // lepton control regions
  cout<<"start zmumu region data"<<endl;
  lepdatamchist(&outfile,1,category,observables,observables_2D,lumi,applyQGLReweight,doShapeSystematics); 
  cout<<"start wmunu region data"<<endl;
  lepdatamchist(&outfile,2,category,observables,observables_2D,lumi,applyQGLReweight,doShapeSystematics); 
  cout<<"start zee region data"<<endl;
  lepdatamchist(&outfile,3,category,observables,observables_2D,lumi,applyQGLReweight,doShapeSystematics); 
  cout<<"start wenu region data"<<endl;
  lepdatamchist(&outfile,4,category,observables,observables_2D,lumi,applyQGLReweight,doShapeSystematics);     
  // top control regions
  cout<<"start top+mu region data"<<endl;
  topdatamchist(&outfile,7,category,observables,observables_2D,lumi,applyQGLReweight,makeResonantSelection,doShapeSystematics);
  cout<<"start Top+el region data"<<endl;
  topdatamchist(&outfile,8,category,observables,observables_2D,lumi,applyQGLReweight,makeResonantSelection,doShapeSystematics);

  //add qcd data templates
  TFile* qcdfile_data = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/QCD/templates.root");
  if(qcdfile_data){
    cout<<"Take templates QCD from data"<<endl;
    vector<float> met_bins = selectBinning("met",category);
    TH1F*  qcd_nominal    = new TH1F("qbkghistDD_met","",int(met_bins.size()-1),&met_bins[0]);
    TH1F*  qcd_nominal_up = new TH1F("qbkghistDD_shapeUp_met","",int(met_bins.size()-1),&met_bins[0]);
    TH1F*  qcd_nominal_dw = new TH1F("qbkghistDD_shapeDw_met","",int(met_bins.size()-1),&met_bins[0]);

    
    TH1F* temp = NULL;
    if(category <= 1)
      temp = (TH1F*) qcdfile_data->Get("hQCD_MonoJ_nominal");
    else
      temp = (TH1F*) qcdfile_data->Get("hQCD_MonoV_nominal");

    for(int iBinX = 0; iBinX < qcd_nominal->GetNbinsX(); iBinX++)   
      qcd_nominal->SetBinContent(iBinX+1,temp->GetBinContent(temp->FindBin(qcd_nominal->GetBinCenter(iBinX+1))));
    
    if(category <= 1)
      temp = (TH1F*) qcdfile_data->Get("hQCD_MonoJ_AllUp");
    else
      temp = (TH1F*) qcdfile_data->Get("hQCD_MonoV_AllUp");

    for(int iBinX = 0; iBinX < qcd_nominal->GetNbinsX(); iBinX++)   
      qcd_nominal_up->SetBinContent(iBinX+1,temp->GetBinContent(temp->FindBin(qcd_nominal->GetBinCenter(iBinX+1))));

    if(category <= 1)
      temp = (TH1F*) qcdfile_data->Get("hQCD_MonoJ_AllDown");
    else
      temp = (TH1F*) qcdfile_data->Get("hQCD_MonoV_AllDown");

    for(int iBinX = 0; iBinX < qcd_nominal->GetNbinsX(); iBinX++)   
      qcd_nominal_dw->SetBinContent(iBinX+1,temp->GetBinContent(temp->FindBin(qcd_nominal->GetBinCenter(iBinX+1))));

    outfile.cd();
    qcd_nominal->Write();
    qcd_nominal_up->Write();
    qcd_nominal_dw->Write();

  }

  outfile.Close();
}
