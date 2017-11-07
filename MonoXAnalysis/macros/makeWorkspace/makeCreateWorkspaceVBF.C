#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <utility>
#include <sstream>

#include "TSystem.h"
#include "workspaceUtils.h"

using namespace std;

/// some basic options to be set manually
static float scaleQCD      = 2;      // scale QCD prediction in the signal region in case a data driven template is not found
static bool  useQCDDataDriven = true; // wether use MC or DD
static bool  connectWZ     = true;   // apply a Z/W ratio
static bool  correlateEWK  = true;   // to correlate EWK uncertainties across bins on the Z/gamma Z/W ratio
static bool  connectEWKQCD = true;   // Connect EWK-Z and Z-QCD TFs in case of VBF analysis
static bool  runOnlySignal     = false;    // run only on signal templates --> workspace with only signals
static bool  runOnlyBackground = true;    // run only on signal templates --> workspace with only background
static bool  freezeZQCDOverZEWKNuisances = true; // set to constant the Z-QCD / Z-EWK constraint in VBF
static bool  addBinByBinMCUncertainty = true;
static float scaleWZUncertainty     = 1.0; // scale up/down the size of theory uncertanty on the Z/W-QCD ratio;                                                                                      
static float scaleWZEWKUncertainty  = 1.0; // scale up/down the size of theory uncertanty on the Z/W-EWK ratio;                                                                                       
static float scaleSignal = 1;
static bool  applyUncertaintyOnNumerator = true; // write the perturbation nuisance on the numerator

// function to create workspace, to be run from a release which has the combine package
void makeCreateWorkspaceVBF(string   inputName,                        // input template file
			    Category category,
			    string   outputName    = "workspace.root", // output workspace name
			    string   observable    = "met",            // observable 1D or 2D
			    bool     splitEWKQCD   = true, // split in two independent set of TFs
			    bool     RunOnlySignal     = false,  // run only on signal templates --> workspace with only signals 
			    bool     RunOnlyBackground = false,  // run only on signal templates --> workspace with only background 
			    string   mediatorMass  = "125"   // Med mass
			    ){

  // parsing
  runOnlySignal     = RunOnlySignal;
  runOnlyBackground = RunOnlyBackground;

  // basic loads
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  RooMsgService::instance().setSilentMode(kTRUE); 
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  // to load all the variables information
  initializeBinning();

  // for templates and sys naming
  string suffix = "VBF";
  
  // create the output workspac
  cout<<"Create output file ..."<<endl;
  TFile *outfile = new TFile(outputName.c_str(),"RECREATE");
  outfile->cd();

  // Select observable and binning
  cout<<"Load binning and observable ..."<<endl;
  double xMin = 0., xMax = 0.;
  vector<double> bins = selectBinning(observable,category);  
  RooBinning *binning = NULL;

  if(not bins.empty()){ // non empty
    xMin = bins.at(0);
    xMax = bins.back();
    binning = new RooBinning(bins.size()-1,&bins[0],(observable+"_"+suffix+"_binning").c_str()); 
  }
  else{
    
    bins.clear();
    bin2D bin = selectBinning2D(observable,category);    
    if(not bin.binX.empty() and not bin.binY.empty()){ // in case of 2D analysis --> unrolled histo
      xMin = 0.;
      xMax = double((bin.binX.size()-1)*(bin.binY.size()-1));
      for(size_t iBin = xMin; iBin <= xMax ; iBin++)
	bins.push_back(double(iBin));
      binning = new RooBinning(bins.size()-1,&bins[0],(observable+"_"+suffix+"_binning").c_str()); 
    }      
    else
      cout<<"Binning not implemented for the observable "<<observable<<" --> please define it "<<endl;
  }
  
  
  // Build the observable
  RooRealVar* met = new RooRealVar((observable+"_"+suffix).c_str(),"",xMin,xMax);
    
  met->setBinning(*binning);
  RooArgList vars(*met);

  // Templates
  cout<<"Open inputFile ..."<<endl;
  TFile* templatesfile = TFile::Open(inputName.c_str());

  ///////////////////////////////////////
  // -------- SIGNAL REGION  -------- ///
  ///////////////////////////////////////

  cout<<"Make SR templates ..."<<endl;
  // create a workspace for the signal region
  RooWorkspace wspace_SR(("SR_"+suffix).c_str(),(suffix+"_SR").c_str());

  // Signal shape
  if(runOnlySignal or not runOnlyBackground){

    // different Higgs invisible signals
    TH1F* ggH   = (TH1F*)templatesfile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+observable).c_str());
    TH1F* qqH   = (TH1F*)templatesfile->FindObjectAny(("vbfHhist_"+mediatorMass+"_"+observable).c_str());
    TH1F* wH    = (TH1F*)templatesfile->FindObjectAny(("wHhist_"+mediatorMass+"_"+observable).c_str());
    TH1F* zH    = (TH1F*)templatesfile->FindObjectAny(("zHhist_"+mediatorMass+"_"+observable).c_str());
    TH1F* ggZH  = (TH1F*)templatesfile->FindObjectAny(("ggZHhist_"+mediatorMass+"_"+observable).c_str());

    ggH->Scale(scaleSignal);
    qqH->Scale(scaleSignal);
    wH->Scale(scaleSignal);
    zH->Scale(scaleSignal);
    ggZH->Scale(scaleSignal);

    addTemplate("ggH_SR_"+suffix,vars,wspace_SR,ggH);
    addTemplate("qqH_SR_"+suffix,vars,wspace_SR,qqH);
    addTemplate("WH_SR_"+suffix, vars,wspace_SR,wH);
    addTemplate("ZH_SR_"+suffix, vars,wspace_SR,zH);
    addTemplate("ggZH_SR_"+suffix, vars,wspace_SR,ggZH);
    
    // generate MC statist bin-by-bin variations
    if(addBinByBinMCUncertainty){
      generateStatTemplate("ggH_SR_"+suffix,vars,wspace_SR,ggH,1);
      generateStatTemplate("qqH_SR_"+suffix,vars,wspace_SR,qqH,1);
      generateStatTemplate("WH_SR_"+suffix,vars,wspace_SR,wH,1);
      generateStatTemplate("ZH_SR_"+suffix,vars,wspace_SR,zH,1);
      generateStatTemplate("ggZH_SR_"+suffix,vars,wspace_SR,ggZH,1);      
    }    
  }

  
  //// Background part in the signal region
  if(not runOnlySignal or runOnlyBackground){

    // Add Data
    addTemplate("data_obs_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("datahist_"+observable).c_str()));
    

    /////////////////////
    // Zvv QCD background --> to be extracted from CRs
    TH1F* znn_SR_hist = (TH1F*) templatesfile->FindObjectAny(("zinvhist_"+observable).c_str());
    TH1F* znn_ewk_SR_hist = (TH1F*) templatesfile->FindObjectAny(("ewkbkgzhist_"+observable).c_str());
    TH1F* znn_SR_total_hist = NULL;

    RooArgList znn_SR_bins; 
    RooArgList znn_ewk_SR_bins;

    if(not splitEWKQCD){// in case one wants to keep them independently
      znn_SR_total_hist = (TH1F*) znn_SR_hist->Clone(Form("%s_total",znn_SR_hist->GetName()));
      znn_SR_total_hist->Add(znn_ewk_SR_hist);
      makeBinList("Znunu_SR_"+suffix,*met,wspace_SR,znn_SR_total_hist,znn_SR_bins,false);    
    }
    else{	

      // create a RooParametric hist with one RooRealVar per bin 
      makeBinList("Znunu_SR_"+suffix,*met,wspace_SR,znn_SR_hist,znn_SR_bins,false);    
      if(not connectEWKQCD) // without explicit link between the two processes
	makeBinList("Znunu_EWK_SR_"+suffix,*met,wspace_SR,znn_ewk_SR_hist,znn_ewk_SR_bins,false);
      else{
	// make the Z-QCD / Z-EWK ratio and propagate stat uncertainty --> theory ones goes only on Z/W-QCD and Z/W-EWK ratios
	TH1F* znn_qcd_SR = (TH1F*) znn_SR_hist->Clone("z_qcd_over_z_ewk");
	znn_qcd_SR->Divide(znn_ewk_SR_hist);
	// create Z-QCD/Z-EWK link                                                                                                                                                               
	vector<pair<RooRealVar*,TH1*> > znn_ewk_SR_syst;
	makeConnectedBinList("Znunu_EWK_SR_"+suffix,*met,wspace_SR,
			     znn_qcd_SR,
			     znn_ewk_SR_syst, //list of systematic variations for the TFs                                                                                                             
			     znn_SR_bins,     //bins for Znunu                                                                                                                                        
			     &znn_ewk_SR_bins, 
			     observable,
			     applyUncertaintyOnNumerator,
			     not freezeZQCDOverZEWKNuisances // decide if the stat uncertainty has to be frozen or not
			     );
	
      }      
    }

    /////////////////////
    // Look at the W+jets

    TH1F* wln_SR_hist = (TH1F*) templatesfile->FindObjectAny(("wjethist_"+observable).c_str());
    TH1F* wln_ewk_SR_hist = (TH1F*) templatesfile->FindObjectAny(("ewkbkgwhist_"+observable).c_str());
    TH1F* wln_SR_total_hist = NULL;

    RooArgList wln_SR_bins;
    RooArgList wln_ewk_SR_bins;

    RooRealVar* wln_SR_re1 = new RooRealVar("ZW_QCD_SR_RenScale1",""  ,0.,-5.,5.);
    RooRealVar* wln_SR_fa1 = new RooRealVar("ZW_QCD_SR_FactScale1","" ,0.,-5.,5.);
    RooRealVar* wln_SR_re2 = new RooRealVar("ZW_QCD_SR_RenScale2",""  ,0.,-5.,5.);
    RooRealVar* wln_SR_fa2 = new RooRealVar("ZW_QCD_SR_FactScale2","" ,0.,-5.,5.);
    RooRealVar* wln_SR_pdf = new RooRealVar("ZW_QCD_SR_PDF",""        ,0.,-5.,5.);
    RooRealVar* wln_SR_ewk = new RooRealVar(("ZW_QCD_SR_"+suffix+"_EWK").c_str(),"",0.,-5.,5.);
    
    RooRealVar* wln_ewk_SR_re1 = new RooRealVar("ZW_EWK_SR_RenScale1",""  ,0.,-5.,5.);
    RooRealVar* wln_ewk_SR_fa1 = new RooRealVar("ZW_EWK_SR_FactScale1","" ,0.,-5.,5.);
    RooRealVar* wln_ewk_SR_re2 = new RooRealVar("ZW_EWK_SR_RenScale2",""  ,0.,-5.,5.);
    RooRealVar* wln_ewk_SR_fa2 = new RooRealVar("ZW_EWK_SR_FactScale2","" ,0.,-5.,5.);
    RooRealVar* wln_ewk_SR_pdf = new RooRealVar("ZW_EWK_SR_PDF",""        ,0.,-5.,5.);
    RooRealVar* wln_ewk_SR_ewk = new RooRealVar(("ZW_EWK_SR_"+suffix+"_EWK").c_str(),"",0.,-5.,5.);
    

    if(not splitEWKQCD){ // summing them up

      wln_SR_total_hist = (TH1F*) wln_SR_hist->Clone(Form("%s_total",wln_SR_hist->GetName()));
      wln_SR_total_hist->Add(wln_ewk_SR_hist);

      if (!connectWZ)
	makeBinList("WJets_SR_"+suffix,*met,wspace_SR,wln_SR_total_hist,wln_SR_bins,true);      
      else{
      
	vector<pair<RooRealVar*,TH1*> > wln_SR_syst;

	// V-EWK part
	TH1* uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale1_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_RenScale1_"+observable);
	uncertaintyPartial->Multiply(wln_ewk_SR_hist);
	uncertaintyPartial->Divide(wln_SR_total_hist);
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_re1,uncertaintyPartial));

	uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale1_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_FactScale1_"+observable);
	uncertaintyPartial->Multiply(wln_ewk_SR_hist);
	uncertaintyPartial->Divide(wln_SR_total_hist);
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_fa1,uncertaintyPartial));

	uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale2_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_RenScale2_"+observable);
	uncertaintyPartial->Multiply(wln_ewk_SR_hist);
	uncertaintyPartial->Divide(wln_SR_total_hist);
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_re2,uncertaintyPartial));

	uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale2_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_FactScale2_"+observable);
	uncertaintyPartial->Multiply(wln_ewk_SR_hist);
	uncertaintyPartial->Divide(wln_SR_total_hist);
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_fa2,uncertaintyPartial));

	uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_PDF_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_PDF_"+observable);
	uncertaintyPartial->Multiply(wln_ewk_SR_hist);
	uncertaintyPartial->Divide(wln_SR_total_hist);
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_pdf,uncertaintyPartial));
	
	uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_EWK_"+observable);
	uncertaintyPartial->Multiply(wln_ewk_SR_hist);
	uncertaintyPartial->Divide(wln_SR_total_hist);

	if(not correlateEWK)
	  wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(NULL,uncertaintyPartial));
	else
	  wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_ewk,uncertaintyPartial));


	// V-QCD part
	uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale1_"+observable).c_str()),scaleWZUncertainty,"");
	uncertaintyPartial->Multiply(wln_SR_hist);
	uncertaintyPartial->Divide(wln_SR_total_hist);
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_re1,uncertaintyPartial));

	uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale1_"+observable).c_str()),scaleWZUncertainty,"");
	uncertaintyPartial->Multiply(wln_SR_hist);
	uncertaintyPartial->Divide(wln_SR_total_hist);
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_fa1,uncertaintyPartial));

	uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale2_"+observable).c_str()),scaleWZUncertainty,"");
	uncertaintyPartial->Multiply(wln_SR_hist);
	uncertaintyPartial->Divide(wln_SR_total_hist);
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_re2,uncertaintyPartial));

	uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale2_"+observable).c_str()),scaleWZUncertainty,"");
	uncertaintyPartial->Multiply(wln_SR_hist);
	uncertaintyPartial->Divide(wln_SR_total_hist);
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_fa2,uncertaintyPartial));

	uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_PDF_"+observable).c_str()),scaleWZUncertainty,"");
	uncertaintyPartial->Multiply(wln_SR_hist);
	uncertaintyPartial->Divide(wln_SR_total_hist);
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_pdf,uncertaintyPartial));

	uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str()),scaleWZUncertainty,"");
	uncertaintyPartial->Multiply(wln_SR_hist);
	uncertaintyPartial->Divide(wln_SR_total_hist);

	if(not correlateEWK)
	  wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(NULL,uncertaintyPartial));
	else
	  wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_ewk,uncertaintyPartial));

	TH1F* zw_ratio_total = (TH1F*) znn_SR_total_hist->Clone("zw_ratio_total");
	zw_ratio_total->Divide(wln_SR_total_hist);
	
	// create Z/W link QCD
	makeConnectedBinList("WJets_SR_"+suffix,*met,wspace_SR,
			     zw_ratio_total, // ratio
			     wln_SR_syst,  //list of systematic variations for the TFs
			     znn_SR_bins,  //bins for Znunu
			     &wln_SR_bins, // W+jets -> empty list
			     observable,
			     applyUncertaintyOnNumerator
			     );

      }
    }
    else{
      
      if (!connectWZ){ // independent links also here
	makeBinList("WJets_SR_"+suffix,*met,wspace_SR,wln_SR_hist,wln_SR_bins,true);
	makeBinList("WJets_EWK_SR_"+suffix,*met,wspace_SR,wln_ewk_SR_hist,wln_ewk_SR_bins,true);
      }
      else{

	// set of correlated systematic uncertainties for the Z/W ratio
	vector<pair<RooRealVar*,TH1*> > wln_ewk_SR_syst;
	wln_ewk_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_re1,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale1_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_RenScale1_"+observable)));
	wln_ewk_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_fa1,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale1_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_FactScale1_"+observable)));
	
	wln_ewk_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_re2,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale2_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_RenScale2_"+observable)));
	
	wln_ewk_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_fa2,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale2_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_FactScale2_"+observable)));
	
	wln_ewk_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_pdf,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_PDF_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_PDF_"+observable)));
	
	if(not correlateEWK)
	  wln_ewk_SR_syst.push_back(pair<RooRealVar*,TH1*>(NULL,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str()),scaleWZUncertainty,"")));
	else
	  wln_ewk_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_ewk,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str()),scaleWZUncertainty,"")));
	
	
	////
	vector<pair<RooRealVar*,TH1*> > wln_SR_syst;    
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_re1,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale1_"+observable).c_str()),scaleWZUncertainty,"")));
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_fa1,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale1_"+observable).c_str()),scaleWZUncertainty,"")));
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_re2,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale2_"+observable).c_str()),scaleWZUncertainty,"")));
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_fa2,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale2_"+observable).c_str()),scaleWZUncertainty,"")));
	wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_pdf,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_PDF_"+observable).c_str()),scaleWZUncertainty,"")));
	if(not correlateEWK)
	  wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(NULL,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str()),scaleWZUncertainty,"")));
	else
	  wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_ewk,cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str()),scaleWZUncertainty,"")));
	
       
	// create Z/W link QCD
	makeConnectedBinList("WJets_SR_"+suffix,*met,wspace_SR,
			     (TH1F*)templatesfile->FindObjectAny(("zwjcorewkhist_"+observable).c_str()), //Z/W ratio --> central value + stat unc.
			     wln_SR_syst, //list of systematic variations for the TFs
			     znn_SR_bins, //bins for Znunu
			     &wln_SR_bins, // W+jets -> empty list
			     observable,
			     applyUncertaintyOnNumerator
			     );
	
	makeConnectedBinList("WJets_EWK_SR_"+suffix,*met,wspace_SR,
			     (TH1F*)templatesfile->FindObjectAny(("zwjewkcorhist_"+observable).c_str()), //Z/W ratio --> central value + stat unc.
			     wln_ewk_SR_syst, //list of systematic variations for the TFs
			     znn_ewk_SR_bins, //bins for Znunu
			     &wln_ewk_SR_bins, // W+jets -> empty list
			     observable,
			     applyUncertaintyOnNumerator
			     );
	
      }
    }

    // Other MC backgrounds
    addTemplate("Top_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("tbkghist_"+observable).c_str()));
    addTemplate("ZJets_SR_"+suffix     ,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("zjethist_"+observable).c_str()));
    addTemplate("Dibosons_SR_"+suffix  ,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("dbkghist_"+observable).c_str()));
    addTemplate("GJets_SR_"+suffix     ,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("gbkghist_"+observable).c_str()));
    addTemplate("VGamma_SR_"+suffix     ,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("vgbkghist_"+observable).c_str()));
    
    if(addBinByBinMCUncertainty){
      generateStatTemplate("Top_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("tbkghist_"+observable).c_str()),1);      
      generateStatTemplate("ZJets_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("zjethist_"+observable).c_str()),1);
      generateStatTemplate("Dibosons_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("dbkghist_"+observable).c_str()),1);
      generateStatTemplate("GJets_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("gbkghist_"+observable).c_str()),1);
      generateStatTemplate("VGamma_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("vgbkghist_"+observable).c_str()),1);
    }
    
    // look for DD qcd background,otherwise MC scaled by a factor 2
    TH1F* qcdhist = (TH1F*)templatesfile->FindObjectAny(("qbkghistDD_"+observable).c_str());
    if(qcdhist and useQCDDataDriven){
      addTemplate("QCD_SR_"+suffix,vars,wspace_SR,qcdhist);
      generateStatTemplate("QCD_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("qbkghistDD_"+observable).c_str()),1); // to take into account TF uncertainties	
    }
    else{
      qcdhist = (TH1F*)templatesfile->FindObjectAny(("qbkghist_"+observable).c_str());
      qcdhist->Scale(scaleQCD);  
      addTemplate("QCD_SR_"+suffix,vars,wspace_SR,qcdhist);
      generateStatTemplate("QCD_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("qbkghist_"+observable).c_str()),1);
    }
  
    /// start with control regions
    RooWorkspace* wspace_ZM = NULL;
    RooWorkspace* wspace_ZE = NULL;
    RooWorkspace* wspace_WM = NULL;
    RooWorkspace* wspace_WE = NULL;

        
    ////////////////////////////////////
    // -------- CR Di-Muon  -------- //
    /////////////////////////////////// 
    
    cout<<"Make CR Di-Muon  templates ..."<<endl;    
    wspace_ZM = new RooWorkspace(("ZM_"+suffix).c_str(),("ZM_"+suffix).c_str());
    
    addTemplate("data_obs_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("datahistzmm_"+observable).c_str()));
    
    // Other MC backgrounds in dimuon control region
    addTemplate("WJets_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_"+observable).c_str()));
    addTemplate("Top_ZM_"+suffix ,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_"+observable).c_str()));
    addTemplate("QCD_ZM_"+suffix ,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("qbkghistzmm_"+observable).c_str()));
    addTemplate("Dibosons_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_"+observable).c_str()));
    addTemplate("VGamma_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("vgbkghistzmm_"+observable).c_str()));
    addTemplate("WJets_EWK_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("ewkwbkghistzmm_"+observable).c_str()));
    
    if(addBinByBinMCUncertainty){
      generateStatTemplate("WJets_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_"+observable).c_str()),1);
      generateStatTemplate("Top_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_"+observable).c_str()),1);
      generateStatTemplate("QCD_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("qbkghistzmm_"+observable).c_str()),1);
      generateStatTemplate("Dibosons_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_"+observable).c_str()),1);
      generateStatTemplate("VGamma_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("vgbkghistzmm_"+observable).c_str()),1);
      generateStatTemplate("WJets_EWK_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("ewkwbkghistzmm_"+observable).c_str()),1);
    }
    
    
    if(splitEWKQCD){

      vector<pair<RooRealVar*,TH1*> > znn_ZM_syst;
      vector<pair<RooRealVar*,TH1*> > znn_ewk_ZM_syst;    

      float rescale = 1;
      if(category == Category::VBFrelaxed)
	rescale = 1.5;

      makeConnectedBinList("Znunu_ZM_"+suffix,
			   *met,
			   *wspace_ZM,
			   (TH1F*)templatesfile->FindObjectAny(("zmmcorhist_"+observable).c_str()),
			   znn_ZM_syst,znn_SR_bins,
			   NULL,observable,applyUncertaintyOnNumerator,true,rescale);

      makeConnectedBinList("Znunu_EWK_ZM_"+suffix,
			   *met,*wspace_ZM,
			   (TH1F*)templatesfile->FindObjectAny(("zewkmmcorhist_"+observable).c_str()),
			   znn_ewk_ZM_syst,znn_ewk_SR_bins,
			   NULL,observable,applyUncertaintyOnNumerator,true,rescale);
    }
    
    else{

      vector<pair<RooRealVar*,TH1*> > znn_ZM_syst;
      TH1F* numerator = (TH1F*) znn_SR_total_hist->Clone(("zmmcorhist_"+observable).c_str());
      TH1F* denominator = (TH1F*) ((TH1F*) templatesfile->Get(("ZM/vllbkghistzmm_"+observable).c_str()))->Clone(Form("denzmm_%s",observable.c_str()));
      denominator->Add((TH1F*) templatesfile->Get(("ZM/ewkzbkghistzmm_"+observable).c_str()));
      numerator->Divide(denominator);

      makeConnectedBinList("Znunu_ZM_"+suffix,
			   *met,
			   *wspace_ZM,numerator,
			   znn_ZM_syst,znn_SR_bins,
			   NULL,observable,applyUncertaintyOnNumerator);
      
    }

    ////////////////////////////////////
    // -------- CR Di-Electron  -------- //
    /////////////////////////////////// 
    
    cout<<"Make CR Di-Electron  templates ..."<<endl;    
    wspace_ZE = new RooWorkspace(("ZE_"+suffix).c_str(),("ZE_"+suffix).c_str());
    
    addTemplate("data_obs_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("datahistzee_"+observable).c_str()));
    
    // Other MC backgrounds in dimuon control region
    addTemplate("WJets_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_"+observable).c_str()));
    addTemplate("Top_ZE_"+suffix ,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("tbkghistzee_"+observable).c_str()));
    addTemplate("QCD_ZE_"+suffix ,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("qbkghistzee_"+observable).c_str()));
    addTemplate("Dibosons_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("dbkghistzee_"+observable).c_str()));
    addTemplate("VGamma_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("vgbkghistzee_"+observable).c_str()));
    addTemplate("WJets_EWK_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("ewkwbkghistzee_"+observable).c_str()));
    
    if(addBinByBinMCUncertainty){
      generateStatTemplate("WJets_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_"+observable).c_str()),1);
      generateStatTemplate("Top_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("tbkghistzee_"+observable).c_str()),1);
      generateStatTemplate("QCD_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("qbkghistzee_"+observable).c_str()),1);
      generateStatTemplate("Dibosons_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("dbkghistzee_"+observable).c_str()),1);
      generateStatTemplate("VGamma_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("vgbkghistzee_"+observable).c_str()),1);
      generateStatTemplate("WJets_EWK_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("ewkwbkghistzee_"+observable).c_str()),1);
    }
    
    
    if(splitEWKQCD){

      vector<pair<RooRealVar*,TH1*> > znn_ZE_syst;
      vector<pair<RooRealVar*,TH1*> > znn_ewk_ZE_syst;    

      float rescale = 1;
      if(category == Category::VBFrelaxed)
	rescale = 1.5;

      makeConnectedBinList("Znunu_ZE_"+suffix,
			   *met,
			   *wspace_ZE,
			   (TH1F*)templatesfile->FindObjectAny(("zeecorhist_"+observable).c_str()),
			   znn_ZE_syst,znn_SR_bins,NULL,observable,applyUncertaintyOnNumerator,true,rescale);

      makeConnectedBinList("Znunu_EWK_ZE_"+suffix,
			   *met,
			   *wspace_ZE,
			   (TH1F*)templatesfile->FindObjectAny(("zewkeecorhist_"+observable).c_str()),
			   znn_ewk_ZE_syst,znn_ewk_SR_bins,NULL,observable,applyUncertaintyOnNumerator,true,rescale);

    }    
    else{
      
      vector<pair<RooRealVar*,TH1*> > znn_ZE_syst;
      TH1F* numerator = (TH1F*) znn_SR_total_hist->Clone(("zeecorhist_"+observable).c_str());
      TH1F* denominator = (TH1F*) ((TH1F*) templatesfile->Get(("ZE/vllbkghistzee_"+observable).c_str()))->Clone(Form("denzee_%s",observable.c_str()));
      denominator->Add((TH1F*) templatesfile->Get(("ZE/ewkzbkghistzee_"+observable).c_str()));
      numerator->Divide(denominator);
      makeConnectedBinList("Znunu_ZE_"+suffix,*met,*wspace_ZE,numerator,znn_ZE_syst,znn_SR_bins,NULL,observable,applyUncertaintyOnNumerator);
      
    }

    ///////////////////////////////////////
    // -------- CR Single-Muon  -------- //
    //////////////////////////////////////
    
    cout<<"Make CR Single-Mu  templates ..."<<endl;
    
    wspace_WM = new RooWorkspace(("WM_"+suffix).c_str(),("WM_"+suffix).c_str());
    addTemplate("data_obs_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("datahistwmn_"+observable).c_str()));
      
    // Other MC backgrounds in single muon control region
    addTemplate("ZJets_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_"+observable).c_str()));
    addTemplate("Top_WM_"+suffix  ,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_"+observable).c_str()));
    addTemplate("QCD_WM_"+suffix  ,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("qbkghistwmn_"+observable).c_str()));
    addTemplate("Dibosons_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_"+observable).c_str()));
    addTemplate("VGamma_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("vgbkghistwmn_"+observable).c_str()));
    addTemplate("ZJets_EWK_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("ewkzbkghistwmn_"+observable).c_str()));
    
    if(addBinByBinMCUncertainty){
      generateStatTemplate("ZJets_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_"+observable).c_str()),1);
      generateStatTemplate("Top_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_"+observable).c_str()),1);
      generateStatTemplate("QCD_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("qbkghistwmn_"+observable).c_str()),1);
      generateStatTemplate("Dibosons_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_"+observable).c_str()),1);
      generateStatTemplate("VGamma_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("vgbkghistwmn_"+observable).c_str()),1);
      generateStatTemplate("ZJets_EWK_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("ewkzbkghistwmn_"+observable).c_str()),1);
    }
    
    
    if(splitEWKQCD){

      vector<pair<RooRealVar*,TH1*> > wln_WM_syst;
      vector<pair<RooRealVar*,TH1*> > wln_ewk_WM_syst;
	
      makeConnectedBinList("WJets_WM_"+suffix,
			   *met,
			   *wspace_WM,
			   (TH1F*)templatesfile->FindObjectAny(("wmncorhist_"+observable).c_str()),
			   wln_WM_syst,wln_SR_bins,NULL,observable,applyUncertaintyOnNumerator);
      
      makeConnectedBinList("WJets_EWK_WM_"+suffix,
			   *met,
			   *wspace_WM,
			   (TH1F*)templatesfile->FindObjectAny(("wewkmncorhist_"+observable).c_str()),
			   wln_ewk_WM_syst,wln_ewk_SR_bins,NULL,observable,applyUncertaintyOnNumerator);
    }
    else{

      vector<pair<RooRealVar*,TH1*> > wln_WM_syst;
      TH1F* numerator = (TH1F*) wln_SR_total_hist->Clone(("wmncorhist_"+observable).c_str());
      TH1F* denominator = (TH1F*) ((TH1F*) templatesfile->Get(("WM/vlbkghistwmn_"+observable).c_str()))->Clone(Form("denwmn_%s",observable.c_str()));
      denominator->Add((TH1F*) templatesfile->Get(("WM/ewkwbkghistwmn_"+observable).c_str()));
      numerator->Divide(denominator);

      // V-EWK part                                                                                                                                                                                   
      TH1* uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale1_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_RenScale1_"+observable);
      uncertaintyPartial->Multiply(numerator);
      wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_re1,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale1_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_FactScale1_"+observable);
      uncertaintyPartial->Multiply(numerator);
      wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_fa1,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale2_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_RenScale2_"+observable);
      uncertaintyPartial->Multiply(numerator);
      wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_re2,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale2_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_FactScale2_"+observable);
      uncertaintyPartial->Multiply(numerator);
      wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_fa2,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_PDF_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_PDF_"+observable);
      uncertaintyPartial->Multiply(numerator);
      wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_pdf,uncertaintyPartial));
	
      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_EWK_"+observable);
      uncertaintyPartial->Multiply(numerator);

      if(not correlateEWK)
	wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(NULL,uncertaintyPartial));
      else
	wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_ewk,uncertaintyPartial));
      

      // V-QCD part
      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale1_"+observable).c_str()),scaleWZUncertainty,"");
      uncertaintyPartial->Multiply(numerator);
      wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_re1,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale1_"+observable).c_str()),scaleWZUncertainty,"");
      uncertaintyPartial->Multiply(numerator);
      wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_fa1,uncertaintyPartial));
      
      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale2_"+observable).c_str()),scaleWZUncertainty,"");
      uncertaintyPartial->Multiply(numerator);
      wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_re2,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale2_"+observable).c_str()),scaleWZUncertainty,"");
      uncertaintyPartial->Multiply(numerator);
      wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_fa2,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_PDF_"+observable).c_str()),scaleWZUncertainty,"");
      uncertaintyPartial->Multiply(numerator);
      wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_pdf,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str()),scaleWZUncertainty,"");
      uncertaintyPartial->Multiply(numerator);
      
      if(not correlateEWK)
	wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(NULL,uncertaintyPartial));
      else
	wln_WM_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_ewk,uncertaintyPartial));
      
      makeConnectedBinList("WJets_WM_"+suffix,*met,*wspace_WM,numerator,wln_WM_syst,wln_SR_bins,NULL,observable,applyUncertaintyOnNumerator);
      
    }

    ///////////////////////////////////////////
    // -------- CR Single-Electron  -------- //
    //////////////////////////////////////////
    cout<<"Make CR Single-El  templates ..."<<endl;
    
    wspace_WE = new RooWorkspace(("WE_"+suffix).c_str(),("WE_"+suffix).c_str());
      
    addTemplate("data_obs_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("datahistwen_"+observable).c_str()));
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_"+observable).c_str()));
    addTemplate("Top_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("tbkghistwen_"+observable).c_str()));
    addTemplate("QCD_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("qbkghistwen_"+observable).c_str()));
    addTemplate("Dibosons_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("dbkghistwen_"+observable).c_str()));
    addTemplate("VGamma_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("vgbkghistwen_"+observable).c_str()));
    addTemplate("ZJets_EWK_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("ewkzbkghistwen_"+observable).c_str()));
    
    
    if(addBinByBinMCUncertainty){
      generateStatTemplate("ZJets_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_"+observable).c_str()),1);
      generateStatTemplate("Top_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("tbkghistwen_"+observable).c_str()),1);
      generateStatTemplate("QCD_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("qbkghistwen_"+observable).c_str()),1);
      generateStatTemplate("Dibosons_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("dbkghistwen_"+observable).c_str()),1);
      generateStatTemplate("VGamma_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("vgbkghistwen_"+observable).c_str()),1);
      generateStatTemplate("ZJets_EWK_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("ewkzbkghistwen_"+observable).c_str()),1);
    }

    if(splitEWKQCD){

      vector<pair<RooRealVar*,TH1*> > wln_WE_syst;
      vector<pair<RooRealVar*,TH1*> > wln_ewk_WE_syst;
	
      makeConnectedBinList("WJets_WE_"+suffix,
			   *met,
			   *wspace_WE,
			   (TH1F*)templatesfile->FindObjectAny(("wencorhist_"+observable).c_str()),
			   wln_WE_syst,wln_SR_bins,NULL,observable,applyUncertaintyOnNumerator);

      makeConnectedBinList("WJets_EWK_WE_"+suffix,
			   *met,
			   *wspace_WE,
			   (TH1F*)templatesfile->FindObjectAny(("wewkencorhist_"+observable).c_str()),
			   wln_ewk_WE_syst,wln_ewk_SR_bins,NULL,observable,applyUncertaintyOnNumerator);

    }
    else{

      vector<pair<RooRealVar*,TH1*> > wln_WE_syst;
      TH1F* numerator = (TH1F*) wln_SR_total_hist->Clone(("wmncorhist_"+observable).c_str());
      TH1F* denominator = (TH1F*) ((TH1F*) templatesfile->Get(("WE/vlbkghistwmn_"+observable).c_str()))->Clone(Form("denwmn_%s",observable.c_str()));
      denominator->Add((TH1F*) templatesfile->Get(("WE/ewkwbkghistwmn_"+observable).c_str()));
      numerator->Divide(denominator);

      // V-EWK part                                                                                                                                                                                   
      TH1* uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale1_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_RenScale1_"+observable);
      uncertaintyPartial->Multiply(numerator);
      wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_re1,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale1_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_FactScale1_"+observable);
      uncertaintyPartial->Multiply(numerator);
      wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_fa1,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale2_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_RenScale2_"+observable);
      uncertaintyPartial->Multiply(numerator);
      wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_re2,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale2_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_FactScale2_"+observable);
      uncertaintyPartial->Multiply(numerator);
      wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_fa2,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_PDF_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_PDF_"+observable);
      uncertaintyPartial->Multiply(numerator);
      wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_pdf,uncertaintyPartial));
	
      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str()),scaleWZEWKUncertainty,"ZW_ewk_EWK_"+observable);
      uncertaintyPartial->Multiply(numerator);

      if(not correlateEWK)
	wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(NULL,uncertaintyPartial));
      else
	wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(wln_ewk_SR_ewk,uncertaintyPartial));
      

      // V-QCD part
      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale1_"+observable).c_str()),scaleWZUncertainty,"");
      uncertaintyPartial->Multiply(numerator);
      wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_re1,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale1_"+observable).c_str()),scaleWZUncertainty,"");
      uncertaintyPartial->Multiply(numerator);
      wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_fa1,uncertaintyPartial));
      
      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_RenScale2_"+observable).c_str()),scaleWZUncertainty,"");
      uncertaintyPartial->Multiply(numerator);
      wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_re2,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_FactScale2_"+observable).c_str()),scaleWZUncertainty,"");
      uncertaintyPartial->Multiply(numerator);
      wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_fa2,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_PDF_"+observable).c_str()),scaleWZUncertainty,"");
      uncertaintyPartial->Multiply(numerator);
      wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_pdf,uncertaintyPartial));

      uncertaintyPartial = cloneAndRescale((TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str()),scaleWZUncertainty,"");
      uncertaintyPartial->Multiply(numerator);
      
      if(not correlateEWK)
	wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(NULL,uncertaintyPartial));
      else
	wln_WE_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_ewk,uncertaintyPartial));
      
      makeConnectedBinList("WJets_WE_"+suffix,*met,*wspace_WE,numerator,wln_WE_syst,wln_SR_bins,NULL,observable,applyUncertaintyOnNumerator);
      
    }

    outfile->cd();
    wspace_ZM->Write();
    wspace_ZE->Write();
    wspace_WM->Write();
    wspace_WE->Write();
  }
      
  outfile->cd();
  wspace_SR.Write();    
  outfile->Close();
  return;  
  
}
