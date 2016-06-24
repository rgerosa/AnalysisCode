#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <utility>
#include <sstream>

#include "TSystem.h"
#include "workspaceUtils.h"

using namespace std;

// function to create workspace, to be run from a release which has the combine package
void createWorkspace(string inputName,  // input template file
		     Category category,  // analysis category
		     string outputName    = "workspace.root", // output workspace name
		     string observable    = "met",    // observable 1D or 2D
		     bool   isHiggsInvisible = false, // Higgs invisible or DM analsysis
		     float  scaleQCD      = 1,    // scale for QCD MC in SR
		     bool   connectWZ     = true, // connect W and Z in SR 
		     bool   connectTop    = false, // connect top CR --> fit also top eneriched samples
		     bool   addShapeSystematics = false, // add shapeN2
		     bool   mergeLeptons  = false, // merge mm and ee final sates
		     bool   isCombination = false,    // convention for HIG-16 invisible combo
		     string interaction   = "Vector", // DM interaction 
		     string mediatorMass  = "1000",   // Med mass
		     string DMMass        = "50",     // DM mass
		     bool   isCutAndCount = false,    // to produce a workspace for a cut and count analysis
		     float  normalizeSignal = -99,    // to scale signal templates to a fixed rate
		     std::pair<float,float> xAxisSelection = {-10000,10000}, // define bins for cut and count
		     std::pair<float,float> yAxisSelection = {-10000,10000}  // define bins for cut and count --> only when 2D templates are considered
){
  
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  RooMsgService::instance().setSilentMode(kTRUE); 
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  initializeBinning();

  // for templates and sys naming
  string suffix;
  if(category == Category::monojet) suffix = "MJ";
  else if(category == Category::monoV) suffix = "MV";
  else if(category == Category::VBF) suffix = "VBF";
    
  // create the output workspace
  cout<<"Create output file ..."<<endl;
  TFile *outfile = new TFile(outputName.c_str(),"RECREATE");
  outfile->cd();
  cout<<"Load binning and observable ..."<<endl;

  // Select observable and binning
  double xMin = 0., xMax = 0.;
  double yMin = 0., yMax = 0.;
  vector<double> bins = selectBinning(observable,category);  
  RooBinning *binning = NULL;

  if(not isCutAndCount){ //shape analysis propagate the binning to the RooRealVar

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
  }
  else{  // for cut and count

    // in case of a 1D analysis
    if(not bins.empty()){
      if(xAxisSelection.first > bins.at(0))
	xMin = xAxisSelection.first;
      else
	xMin = bins.at(0);
      if(xAxisSelection.second < bins.back())
	xMax = xAxisSelection.second;
      else
	xMax = bins.back() ;

      bins.clear();
      bins.push_back(xMin);
      bins.push_back(xMax);
      binning = new RooBinning(bins.size()-1,&bins[0],(observable+"_"+suffix+"_binning").c_str());
    }
    else{ 

      bins.clear();
      bin2D bin = selectBinning2D(observable,category);
      if(not bin.binX.empty() and not bin.binY.empty()){ // in case of 2D analysis --> unrolled histo                                               

	// range for the xAxis
	if(xAxisSelection.first > bin.binX.at(0))
	  xMin = xAxisSelection.first;
	else
	  xMin = bin.binX.at(0);
	
	if(xAxisSelection.second < bin.binX.back())
	  xMax = xAxisSelection.second;
	else
	  xMax = bin.binX.back() ;
	
	// range for the yAxis
	if(yAxisSelection.first > bin.binY.at(0))
	  yMin = yAxisSelection.first;
	else
	  yMin = bin.binY.at(0);
	if(yAxisSelection.second < bin.binY.back())
	  yMax = yAxisSelection.second;
	else
	  yMax = bin.binY.back() ;
	
	// find the right location in un-rolled histograms
	bins.clear();
	int xMinPos = std::distance(bin.binX.begin(),std::find(bin.binX.begin(),bin.binX.end(),xMin));
	int xMaxPos = std::distance(bin.binX.begin(),std::find(bin.binX.begin(),bin.binX.end(),xMax));
	if(xMinPos == xMaxPos){
	  cerr<<" 2D cut and count --> extremes are not matching any histo X-axis bin --> please check"<<endl;
	  return;
	}
	int yMinPos = std::distance(bin.binX.begin(),std::find(bin.binY.begin(),bin.binY.end(),yMin));
	int yMaxPos = std::distance(bin.binX.begin(),std::find(bin.binY.begin(),bin.binY.end(),yMax));
	if(yMinPos == yMaxPos){
	  cerr<<" 2D cut and count --> extremes are not matching any histo Y-axis bin --> please check"<<endl;
	  return;
	}

	if(xMinPos %2 == 0 and xMinPos != 0) xMinPos --;
	if(xMaxPos %2 == 0 and xMaxPos != 0) xMaxPos --;
	  
	// unrolling along the X-axis perfomed when filled hte template file
	bins.push_back(xMinPos*(bin.binY.size()-1)+yMinPos);
	bins.push_back(xMaxPos*(bin.binY.size()-1)+yMaxPos);
	binning = new RooBinning(bins.size()-1,&bins[0],(observable+"_"+suffix+"_binning").c_str());
	
      }
    }
  }
  // Build the RooRealVar
  RooRealVar met((observable+"_"+suffix).c_str(),"",xMin,xMax);
  met.setBinning(*binning);
  RooArgList vars(met);

  // Templates
  cout<<"Open inputFile ..."<<endl;
  TFile* templatesfile = TFile::Open(inputName.c_str());

  ///////////////////////////////////////
  // -------- SIGNAL REGION  -------- //
  ///////////////////////////////////////
  cout<<"Make SR templates ..."<<endl;
  // create a workspace for the signal region
  RooWorkspace wspace_SR(("SR_"+suffix).c_str(),(suffix+"_SR").c_str());
  // Add Data
  addTemplate("data_obs_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("datahist_"+observable).c_str()),isCutAndCount);
  // Signal shape
  if(!isHiggsInvisible){

    TH1F* monoJ = (TH1F*)templatesfile->FindObjectAny(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
    TH1F* monoW = (TH1F*)templatesfile->FindObjectAny(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
    TH1F* monoZ = (TH1F*)templatesfile->FindObjectAny(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());

    if(isCutAndCount and normalizeSignal > 0){
      if(monoJ and monoJ->Integral(monoJ->FindBin(met.getMin()),monoJ->FindBin(met.getMax())) != 0)
	monoJ->Scale(normalizeSignal/monoJ->Integral(monoJ->FindBin(met.getMin()),monoJ->FindBin(met.getMax())));
      if(monoW and monoW->Integral(monoW->FindBin(met.getMin()),monoW->FindBin(met.getMax())) != 0)
	monoW->Scale(normalizeSignal/monoW->Integral(monoW->FindBin(met.getMin()),monoW->FindBin(met.getMax())));
      if(monoZ and monoZ->Integral(monoZ->FindBin(met.getMin()),monoZ->FindBin(met.getMax())) != 0)
	monoZ->Scale(normalizeSignal/monoZ->Integral(monoZ->FindBin(met.getMin()),monoZ->FindBin(met.getMax())));      
    }

    addTemplate("MonoJ_SR_"+suffix,vars,wspace_SR,monoJ,isCutAndCount);
    addTemplate("MonoW_SR_"+suffix,vars,wspace_SR,monoW,isCutAndCount);
    addTemplate("MonoZ_SR_"+suffix,vars,wspace_SR,monoZ,isCutAndCount);
    
    if(addShapeSystematics){

      addShapeVariations("monoJhist","MonoJ_SR",suffix,observable,vars,wspace_SR,templatesfile,interaction+"_"+mediatorMass+"_"+DMMass,isCombination,isCutAndCount,normalizeSignal);
      addShapeVariations("monoWhist","MonoW_SR",suffix,observable,vars,wspace_SR,templatesfile,interaction+"_"+mediatorMass+"_"+DMMass,isCombination,isCutAndCount,normalizeSignal);
      addShapeVariations("monoZhist","MonoZ_SR",suffix,observable,vars,wspace_SR,templatesfile,interaction+"_"+mediatorMass+"_"+DMMass,isCombination,isCutAndCount,normalizeSignal);
      
      // statistics
      generateStatTemplate("MonoJ_SR_"+suffix,vars,wspace_SR,monoJ,1,isCutAndCount);
      generateStatTemplate("MonoW_SR_"+suffix,vars,wspace_SR,monoW,1,isCutAndCount);
      generateStatTemplate("MonoZ_SR_"+suffix,vars,wspace_SR,monoZ,1,isCutAndCount);
      
    }
  }
  else{

    TH1F* ggH = (TH1F*)templatesfile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+observable).c_str());
    TH1F* qqH = (TH1F*)templatesfile->FindObjectAny(("vbfHhist_"+mediatorMass+"_"+observable).c_str());
    TH1F* wH  = (TH1F*)templatesfile->FindObjectAny(("wHhist_"+mediatorMass+"_"+observable).c_str());
    TH1F* zH  = (TH1F*)templatesfile->FindObjectAny(("zHhist_"+mediatorMass+"_"+observable).c_str());
    TH1F* ggZH  = (TH1F*)templatesfile->FindObjectAny(("ggZHhist_"+mediatorMass+"_"+observable).c_str());

    if(isCutAndCount and normalizeSignal > 0){

      if(ggH and ggH->Integral(ggH->FindBin(met.getMin()),ggH->FindBin(met.getMax())) != 0)
	ggH->Scale(normalizeSignal/ggH->Integral(ggH->FindBin(met.getMin()),ggH->FindBin(met.getMax())));      

      if(qqH and qqH->Integral(qqH->FindBin(met.getMin()),qqH->FindBin(met.getMax())) != 0)
	qqH->Scale(normalizeSignal/qqH->Integral(qqH->FindBin(met.getMin()),qqH->FindBin(met.getMax())));

      if(wH and wH->Integral(wH->FindBin(met.getMin()),wH->FindBin(met.getMax())) != 0)
	wH->Scale(normalizeSignal/wH->Integral(wH->FindBin(met.getMin()),wH->FindBin(met.getMax())));

      if(zH and zH->Integral(zH->FindBin(met.getMin()),zH->FindBin(met.getMax())) != 0)
	zH->Scale(normalizeSignal/zH->Integral(zH->FindBin(met.getMin()),zH->FindBin(met.getMax())));

      if(ggZH and ggZH->Integral(ggZH->FindBin(met.getMin()),ggZH->FindBin(met.getMax())) != 0)
	ggZH->Scale(normalizeSignal/ggZH->Integral(ggZH->FindBin(met.getMin()),ggZH->FindBin(met.getMax())));
    }
    
    addTemplate("ggH_SR_"+suffix,vars,wspace_SR,ggH,isCutAndCount);
    addTemplate("qqH_SR_"+suffix,vars,wspace_SR,qqH,isCutAndCount);
    addTemplate("WH_SR_"+suffix, vars,wspace_SR,wH,isCutAndCount);
    addTemplate("ZH_SR_"+suffix, vars,wspace_SR,zH,isCutAndCount);
    addTemplate("ggZH_SR_"+suffix, vars,wspace_SR,ggZH,isCutAndCount);

    if(addShapeSystematics){

      addShapeVariations("ggHhist","ggH_SR",suffix,observable,vars,wspace_SR,templatesfile,mediatorMass,isCombination,isCutAndCount,normalizeSignal);
      addShapeVariations("vbfHhist","qqH_SR",suffix,observable,vars,wspace_SR,templatesfile,mediatorMass,isCombination,isCutAndCount,normalizeSignal);
      addShapeVariations("wHhist","WH_SR",suffix,observable,vars,wspace_SR,templatesfile,mediatorMass,isCombination,isCutAndCount,normalizeSignal);
      addShapeVariations("zHhist","ZH_SR",suffix,observable,vars,wspace_SR,templatesfile,mediatorMass,isCombination,isCutAndCount,normalizeSignal);
      addShapeVariations("ggZHhist","ggZH_SR",suffix,observable,vars,wspace_SR,templatesfile,mediatorMass,isCombination,isCutAndCount,normalizeSignal);

      // ggH higgs pT uncertainties 
      TH1F* histoRenUp = (TH1F*)templatesfile->FindObjectAny(("ggHhist_renUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoRenDw = (TH1F*)templatesfile->FindObjectAny(("ggHhist_renDw_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoFacUp = (TH1F*)templatesfile->FindObjectAny(("ggHhist_facUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoFacDw = (TH1F*)templatesfile->FindObjectAny(("ggHhist_facDw_"+mediatorMass+"_"+observable).c_str());

      
      if(isCutAndCount and normalizeSignal > 0){
	if(histoRenUp)
	  histoRenUp->Scale(normalizeSignal/histoRenUp->Integral(histoRenUp->FindBin(met.getMin()),histoRenUp->FindBin(met.getMax())));
	if(histoRenDw)
	  histoRenDw->Scale(normalizeSignal/histoRenDw->Integral(histoRenDw->FindBin(met.getMin()),histoRenDw->FindBin(met.getMax())));
	if(histoFacUp)
	  histoFacUp->Scale(normalizeSignal/histoFacUp->Integral(histoFacUp->FindBin(met.getMin()),histoFacUp->FindBin(met.getMax())));
	if(histoFacDw)
	  histoFacDw->Scale(normalizeSignal/histoFacDw->Integral(histoFacDw->FindBin(met.getMin()),histoFacDw->FindBin(met.getMax())));
      }

      if(ggH){
	vector<TH1F*> histoVec;
	if(histoRenUp != 0)
	  histoVec.push_back(histoRenUp);
	if(histoRenDw != 0)
	  histoVec.push_back(histoRenDw);
	if(histoFacUp != 0)
	  histoVec.push_back(histoFacUp);
	if(histoFacDw != 0)
	  histoVec.push_back(histoFacDw);		

	addTemplate("ggH_SR_"+suffix+"_hptUp",vars,wspace_SR,generateEnvelopeMax(histoVec,"ggH_SR_"+suffix),isCutAndCount);
	addTemplate("ggH_SR_"+suffix+"_hptDown",vars,wspace_SR,generateEnvelopeMin(histoVec,"ggH_SR_"+suffix),isCutAndCount);
      }
      
      // statistics
      generateStatTemplate("ggH_SR_"+suffix,vars,wspace_SR,ggH,0.5,isCutAndCount);
      generateStatTemplate("qqH_SR_"+suffix,vars,wspace_SR,qqH,0.5,isCutAndCount);
      generateStatTemplate("WH_SR_"+suffix,vars,wspace_SR,wH,0.5,isCutAndCount);
      generateStatTemplate("ZH_SR_"+suffix,vars,wspace_SR,zH,0.5,isCutAndCount);
      generateStatTemplate("ggZH_SR_"+suffix,vars,wspace_SR,ggZH,0.5,isCutAndCount);      
    }    
  }

  // Zvv background --> to be extracted from CRs
  TH1F* znn_SR_hist = (TH1F*) templatesfile->FindObjectAny(("zinvhist_"+observable).c_str());
  RooArgList znn_SR_bins; 
  // create a RooParametric hist with one RooRealVar per bin 
  makeBinList("Znunu_SR_"+suffix,met,wspace_SR,znn_SR_hist,znn_SR_bins,false,isCutAndCount);
  // Top background --> to be extracted from CRs
  RooArgList top_SR_bins;
  TH1F* top_SR_hist = NULL;

  // for data driven top estimation
  if(connectTop){
    top_SR_hist = (TH1F*) templatesfile->FindObjectAny(("tbkghist_"+observable).c_str());
    RooArgList top_SR_bins; 
    makeBinList("Top_SR_"+suffix,met,wspace_SR,top_SR_hist,top_SR_bins,false,isCutAndCount);
  }
  else{ // rely on MC + systematics
    addTemplate("Top_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("tbkghist_"+observable).c_str()),isCutAndCount);
    if(addShapeSystematics){
      addShapeVariations("tbkghist","Top_SR",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination,isCutAndCount);
      generateStatTemplate("Top_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("tbkghist_"+observable).c_str()),1,isCutAndCount);
    }
  }

  // WJets background --> to be extracted from CRs,with connection to Z->nunu
  TH1F* wln_SR_hist = (TH1F*) templatesfile->FindObjectAny(("wjethist_"+observable).c_str());
  RooArgList wln_SR_bins;
  if (!connectWZ) 
    makeBinList("WJets_SR_"+suffix,met,wspace_SR,wln_SR_hist,wln_SR_bins,true,isCutAndCount);
  else{

    RooRealVar* wln_SR_re1 = new RooRealVar("WJets_SR_RenScale1",""  ,0.,-5.,5.);
    RooRealVar* wln_SR_fa1 = new RooRealVar("WJets_SR_FactScale1","" ,0.,-5.,5.);
    RooRealVar* wln_SR_re2 = new RooRealVar("WJets_SR_RenScale2",""  ,0.,-5.,5.);
    RooRealVar* wln_SR_fa2 = new RooRealVar("WJets_SR_FactScale2","" ,0.,-5.,5.);
    RooRealVar* wln_SR_pdf = new RooRealVar("WJets_SR_PDF",""        ,0.,-5.,5.);
    
    if(not isCutAndCount){
      
      // set of correlated systematic uncertainties for the Z/W ratio
      vector<pair<RooRealVar*,TH1*> > wln_SR_syst;

      // NULL means bin-by-bin
      wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(NULL      ,(TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str())));
      wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_re1,(TH1F*)templatesfile->FindObjectAny(("ZW_RenScale1_"+observable).c_str())));
      wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_fa1,(TH1F*)templatesfile->FindObjectAny(("ZW_FactScale1_"+observable).c_str())));
      wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_re2,(TH1F*)templatesfile->FindObjectAny(("ZW_RenScale2_"+observable).c_str())));
      wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_fa2,(TH1F*)templatesfile->FindObjectAny(("ZW_FactScale2_"+observable).c_str())));
      wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_pdf,(TH1F*)templatesfile->FindObjectAny(("ZW_PDF_"+observable).c_str())));

      // create Z/W link
      makeConnectedBinList("WJets_SR_"+suffix,met,wspace_SR,
			   (TH1F*)templatesfile->FindObjectAny(("zwjcorewkhist_"+observable).c_str()), //Z/W ratio --> central value + stat unc.
			   wln_SR_syst, //list of systematic variations for the TFs
			   znn_SR_bins, //bins for Znunu
			   &wln_SR_bins, // W+jets -> empty list
			   observable);
      
    }
    
    else{
      
      vector<pair<RooRealVar*,systematicCutAndCount> > wln_SR_syst;
      
      RooRealVar* wln_SR_ewk = new RooRealVar(("WJets_SR_"+suffix+"_ZW_EWK").c_str(),""  ,0.,-5.,5.);
      systematicCutAndCount wln_SR_ewk_sys; // single bin everything here
      wln_SR_ewk_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_ewk_"+observable).c_str());
      wln_SR_ewk_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_ewk_"+observable).c_str());
      wln_SR_ewk_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_qcd_"+observable).c_str());
      wln_SR_ewk_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_qcd_"+observable).c_str());
      wln_SR_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(wln_SR_ewk,wln_SR_ewk_sys));
      
      systematicCutAndCount wln_SR_re1_sys;
      wln_SR_re1_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_re1_"+observable).c_str());
      wln_SR_re1_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_re1_"+observable).c_str());
      wln_SR_re1_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_qcd_"+observable).c_str());
      wln_SR_re1_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_qcd_"+observable).c_str());
      wln_SR_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(wln_SR_re1,wln_SR_re1_sys));
      
      systematicCutAndCount wln_SR_re2_sys;
      wln_SR_re2_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_re2_"+observable).c_str());
      wln_SR_re2_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_re2_"+observable).c_str());
      wln_SR_re2_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_qcd_"+observable).c_str());
      wln_SR_re2_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_qcd_"+observable).c_str());
      wln_SR_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(wln_SR_re2,wln_SR_re2_sys));
      
      systematicCutAndCount wln_SR_fa1_sys;
      wln_SR_fa1_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_fa1_"+observable).c_str());
      wln_SR_fa1_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_fa1_"+observable).c_str());
      wln_SR_fa1_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_qcd_"+observable).c_str());
      wln_SR_fa1_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_qcd_"+observable).c_str());
      wln_SR_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(wln_SR_fa1,wln_SR_fa1_sys));
      
      systematicCutAndCount wln_SR_fa2_sys;
      wln_SR_fa2_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_fa2_"+observable).c_str());
      wln_SR_fa2_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_fa2_"+observable).c_str());
      wln_SR_fa2_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_qcd_"+observable).c_str());
      wln_SR_fa2_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_qcd_"+observable).c_str());
      wln_SR_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(wln_SR_fa2,wln_SR_fa2_sys));
    
      systematicCutAndCount wln_SR_pdf_sys;
      wln_SR_pdf_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_pdf_"+observable).c_str());
      wln_SR_pdf_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_pdf_"+observable).c_str());
      wln_SR_pdf_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_qcd_"+observable).c_str());
      wln_SR_pdf_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_qcd_"+observable).c_str());
      wln_SR_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(wln_SR_pdf,wln_SR_pdf_sys));
      
      makeConnectedBinListCutAndCount("WJets_SR_"+suffix,met,wspace_SR,
				      (TH1F*)templatesfile->FindObjectAny(("nhist_zwj_ewk_"+observable).c_str()),
				      (TH1F*)templatesfile->FindObjectAny(("dhist_zwj_ewk_"+observable).c_str()),
				      wln_SR_syst,znn_SR_bins,&wln_SR_bins,observable);
    }
  }

  // Other MC backgrounds
  addTemplate("ZJets_SR_"+suffix     ,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("zjethist_"+observable).c_str()),isCutAndCount);
  addTemplate("Dibosons_SR_"+suffix  ,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("dbkghist_"+observable).c_str()),isCutAndCount);
  addTemplate("GJets_SR_"+suffix     ,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("gbkghist_"+observable).c_str()),isCutAndCount);
  addTemplate("EWKZ_SR_"+suffix      ,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("ewkbkgzhist_"+observable).c_str()),isCutAndCount);
  addTemplate("EWKW_SR_"+suffix      ,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("ewkbkgwhist_"+observable).c_str()),isCutAndCount);


  if(addShapeSystematics){
    addShapeVariations("zjethist","ZJets_SR",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination,isCutAndCount);
    generateStatTemplate("ZJets_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("zjethist_"+observable).c_str()),1,isCutAndCount);

    addShapeVariations("dbkghist","Dibosons_SR",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination,isCutAndCount);
    generateStatTemplate("Dibosons_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("dbkghist_"+observable).c_str()),1,isCutAndCount);

    addShapeVariations("gbkghist","GJets_SR",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination,isCutAndCount);
    generateStatTemplate("GJets_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("gbkghist_"+observable).c_str()),1,isCutAndCount);

    addShapeVariations("ewkbkgzhist","EWKZ_SR",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination,isCutAndCount);
    generateStatTemplate("EWKZ_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("ewkbkgzhist_"+observable).c_str()),1,isCutAndCount);

    addShapeVariations("ewkbkgwhist","EWKW_SR",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination,isCutAndCount);
    generateStatTemplate("EWKW_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("ewkbkgwhist_"+observable).c_str()),1,isCutAndCount);
  }

  // look for DD qcd background,otherwise MC scaled by a factor 2
  TH1F* qcdhist = (TH1F*)templatesfile->FindObjectAny(("qbkghistDD_"+observable).c_str());
  if(qcdhist){
    addTemplate("QCD_SR_"+suffix,vars,wspace_SR,qcdhist,isCutAndCount);
    addTemplate("QCD_SR_"+suffix+"_CMS_QCD_SRUp",vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("qbkghistDD_shapeUp_"+observable).c_str()),isCutAndCount);
    addTemplate("QCD_SR_"+suffix+"_CMS_QCD_SRDown",vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("qbkghistDD_shapeDw_"+observable).c_str()),isCutAndCount);
  }
  else{
    qcdhist = (TH1F*)templatesfile->FindObjectAny(("qbkghist_"+observable).c_str());
    qcdhist->Scale(scaleQCD);  
    addTemplate("QCD_SR_"+suffix,vars,wspace_SR,qcdhist,isCutAndCount);
  }

  RooWorkspace* wspace_ZM = NULL;
  RooWorkspace* wspace_ZE = NULL;
  RooWorkspace* wspace_WM = NULL;
  RooWorkspace* wspace_WE = NULL;
  RooWorkspace* wspace_ZL = NULL;
  RooWorkspace* wspace_WL = NULL;

  if(not mergeLeptons){

    ////////////////////////////////////
    // -------- CR Di-Muon  -------- //
    /////////////////////////////////// 

    cout<<"Make CR Di-Muon  templates ..."<<endl;    
    wspace_ZM = new RooWorkspace(("ZM_"+suffix).c_str(),("ZM_"+suffix).c_str());
  
    addTemplate("data_obs_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("datahistzmm_"+observable).c_str()),isCutAndCount);
  
    // Other MC backgrounds in dimuon control region
    addTemplate("WJets_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_"+observable).c_str()),isCutAndCount);
    addTemplate("Top_ZM_"+suffix ,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_"+observable).c_str()),isCutAndCount);
    addTemplate("QCD_ZM_"+suffix ,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("qbkghistzmm_"+observable).c_str()),isCutAndCount);
    addTemplate("Dibosons_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_"+observable).c_str()),isCutAndCount);
    addTemplate("EWKW_ZM_"+suffix  ,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("ewkwbkghistzmm_"+observable).c_str()),isCutAndCount);
    addTemplate("EWKZ_ZM_"+suffix  ,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("ewkzbkghistzmm_"+observable).c_str()),isCutAndCount);

    if(addShapeSystematics){      
      addShapeVariations("vlbkghistzmm","WJets_ZM",suffix,observable,vars,*wspace_ZM,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("tbkghistzmm","Top_ZM",suffix,observable,vars,*wspace_ZM,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("dbkghistzmm","Dibosons_ZM",suffix,observable,vars,*wspace_ZM,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("ewkzbkghistzmm","EWKZ_ZM",suffix,observable,vars,*wspace_ZM,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("ewkwbkghistzmm","EWKW_ZM",suffix,observable,vars,*wspace_ZM,templatesfile,"",isCombination,isCutAndCount);
    }
    
    if(not isCutAndCount){
      // Z->mumu connected with Z->nunu SR
      vector<pair<RooRealVar*,TH1*> >   znn_ZM_syst;
      makeConnectedBinList("Znunu_ZM_"+suffix,met,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("zmmcorhist_"+observable).c_str()),znn_ZM_syst,znn_SR_bins,NULL,observable);
    }
    else{

      vector<pair<RooRealVar*,systematicCutAndCount> > znn_ZM_syst;
      makeConnectedBinListCutAndCount("Znunu_ZM_"+suffix,met,*wspace_ZM,
				      (TH1F*)templatesfile->FindObjectAny(("nhist_zmm_"+observable).c_str()),
				      (TH1F*)templatesfile->FindObjectAny(("dhist_zmm_"+observable).c_str()),
				      znn_ZM_syst,znn_SR_bins,NULL,observable);
      
    }
    
    ////////////////////////////////////////
    // -------- CR Di-Electron  -------- //
    ///////////////////////////////////////

    cout<<"Make CR Di-Electron  templates ..."<<endl;      
    wspace_ZE = new RooWorkspace(("ZE_"+suffix).c_str(),("ZE_"+suffix).c_str());

    addTemplate("data_obs_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("datahistzee_"+observable).c_str()),isCutAndCount);
    
    // Other MC backgrounds in dielectron control region
    addTemplate("WJets_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_"+observable).c_str()),isCutAndCount);
    addTemplate("Top_ZE_"+suffix  ,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("tbkghistzee_"+observable).c_str()),isCutAndCount);
    addTemplate("QCD_ZE_"+suffix  ,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("qbkghistzee_"+observable).c_str()),isCutAndCount);
    addTemplate("Dibosons_ZE_"+suffix  ,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("dbkghistzee_"+observable).c_str()),isCutAndCount);
    addTemplate("EWKW_ZE_"+suffix  ,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("ewkwbkghistzee_"+observable).c_str()),isCutAndCount);
    addTemplate("EWKZ_ZE_"+suffix  ,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("ewkzbkghistzee_"+observable).c_str()),isCutAndCount);

    if( addShapeSystematics){      
      addShapeVariations("vlbkghistzee","WJets_ZE",suffix,observable,vars,*wspace_ZE,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("tbkghistzee","Top_ZE",suffix,observable,vars,*wspace_ZE,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("dbkghistzee","Dibosons_ZE",suffix,observable,vars,*wspace_ZE,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("ewkzbkghistzee","EWKZ_ZE",suffix,observable,vars,*wspace_ZE,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("ewkwbkghistzee","EWKW_ZE",suffix,observable,vars,*wspace_ZE,templatesfile,"",isCombination,isCutAndCount);
    }

    if(not isCutAndCount){
      // Z->ee connected with Z->nunu SR
      vector<pair<RooRealVar*,TH1*> > znn_ZE_syst;
      makeConnectedBinList("Znunu_ZE_"+suffix,met,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("zeecorhist_"+observable).c_str()),znn_ZE_syst,znn_SR_bins,NULL,observable);

    }
    else{
      vector<pair<RooRealVar*,systematicCutAndCount> > znn_ZE_syst;
      makeConnectedBinListCutAndCount("Znunu_ZE_"+suffix,met,*wspace_ZE,
				      (TH1F*)templatesfile->FindObjectAny(("nhist_zee_"+observable).c_str()),
				      (TH1F*)templatesfile->FindObjectAny(("dhist_zee_"+observable).c_str()),
				      znn_ZE_syst,znn_SR_bins,NULL,observable);


    }
  }
  else{
    
    cout<<"Make CR Di-Lepton  templates ..."<<endl;    
    wspace_ZL = new RooWorkspace(("ZL_"+suffix).c_str(),("ZL_"+suffix).c_str());
  
    addTemplate("data_obs_ZL_"+suffix,vars,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("datahistzll_"+observable).c_str()),isCutAndCount);
  
    // Other MC backgrounds in dimuon control region
    addTemplate("WJets_ZL_"+suffix,vars,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("vlbkghistzll_"+observable).c_str()),isCutAndCount);
    addTemplate("Top_ZL_"+suffix ,vars,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("tbkghistzll_"+observable).c_str()),isCutAndCount);
    addTemplate("QCD_ZL_"+suffix ,vars,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("qbkghistzll_"+observable).c_str()),isCutAndCount);
    addTemplate("Dibosons_ZL_"+suffix,vars,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("dbkghistzll_"+observable).c_str()),isCutAndCount);
    addTemplate("EWKW_ZL_"+suffix  ,vars,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("ewkwbkghistzll_"+observable).c_str()),isCutAndCount);
    addTemplate("EWKZ_ZL_"+suffix  ,vars,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("ewkzbkghistzll_"+observable).c_str()),isCutAndCount);

    if(addShapeSystematics){
      addShapeVariations("vlbkghistzll","WJets_ZL",suffix,observable,vars,*wspace_ZL,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("tbkghistzll","Top_ZL",suffix,observable,vars,*wspace_ZL,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("dbkghistzll","Dibosons_ZL",suffix,observable,vars,*wspace_ZL,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("ewkzbkghistzll","EWKZ_ZL",suffix,observable,vars,*wspace_ZL,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("ewkwbkghistzll","EWKW_ZL",suffix,observable,vars,*wspace_ZL,templatesfile,"",isCombination,isCutAndCount);
    }            

    if(not isCutAndCount){
      // Z->mumu connected with Z->nunu SR
      vector<pair<RooRealVar*,TH1*> >   znn_ZL_syst;
      makeConnectedBinList("Znunu_ZL_"+suffix,met,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("zllcorhist_"+observable).c_str()),znn_ZL_syst,znn_SR_bins,NULL,observable);
    }
    else{
      vector<pair<RooRealVar*,systematicCutAndCount> > znn_ZL_syst;
      makeConnectedBinListCutAndCount("Znunu_ZL_"+suffix,met,*wspace_ZL,
				      (TH1F*)templatesfile->FindObjectAny(("nhist_zll_"+observable).c_str()),
				      (TH1F*)templatesfile->FindObjectAny(("dhist_zll_"+observable).c_str()),
				      znn_ZL_syst,znn_SR_bins,NULL,observable);
    }
  }


  if(not mergeLeptons){
    ///////////////////////////////////////
    // -------- CR Single-Muon  -------- //
    //////////////////////////////////////
    cout<<"Make CR Single-Mu  templates ..."<<endl;
    
    wspace_WM = new RooWorkspace(("WM_"+suffix).c_str(),("WM_"+suffix).c_str());
    addTemplate("data_obs_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("datahistwmn_"+observable).c_str()),isCutAndCount);
    
    // Other MC backgrounds in single muon control region
    addTemplate("ZJets_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_"+observable).c_str()),isCutAndCount);
    addTemplate("Top_WM_"+suffix  ,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_"+observable).c_str()),isCutAndCount);
    addTemplate("QCD_WM_"+suffix  ,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("qbkghistwmn_"+observable).c_str()),isCutAndCount);
    addTemplate("Dibosons_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_"+observable).c_str()),isCutAndCount);
    addTemplate("EWKZ_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("ewkzbkghistwmn_"+observable).c_str()),isCutAndCount);
    addTemplate("EWKW_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("ewkwbkghistwmn_"+observable).c_str()),isCutAndCount);
    
    if(addShapeSystematics){      
      addShapeVariations("vllbkghistwmn","ZJets_WM",suffix,observable,vars,*wspace_WM,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("tbkghistwmn","Top_WM",suffix,observable,vars,*wspace_WM,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("dbkghistwmn","Dibosons_WM",suffix,observable,vars,*wspace_WM,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("ewkzbkghistwmn","EWKZ_WM",suffix,observable,vars,*wspace_WM,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("ewkwbkghistwmn","EWKW_WM",suffix,observable,vars,*wspace_WM,templatesfile,"",isCombination,isCutAndCount);
    }

    if(not isCutAndCount){
      // connected W->munu with W+jets SR
      vector<pair<RooRealVar*,TH1*> > wln_WM_syst;
      makeConnectedBinList("WJets_WM_"+suffix,met,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("wmncorhist_"+observable).c_str()),wln_WM_syst,wln_SR_bins,NULL,observable);
    }
    else{
      vector<pair<RooRealVar*,systematicCutAndCount> > wln_WM_syst;
      makeConnectedBinListCutAndCount("WJets_WM_"+suffix,met,*wspace_WM,
				      (TH1F*)templatesfile->FindObjectAny(("nhist_wmn_"+observable).c_str()),
				      (TH1F*)templatesfile->FindObjectAny(("dhist_wmn_"+observable).c_str()),
				      wln_WM_syst,wln_SR_bins,NULL,observable);

    }

    //////////////////////////////////..../////
    // -------- CR Single-Electron  -------- //
    //////////////////////////////////////////
    cout<<"Make CR Single-El  templates ..."<<endl;
    
    wspace_WE = new RooWorkspace(("WE_"+suffix).c_str(),("WE_"+suffix).c_str());
    
    addTemplate("data_obs_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("datahistwen_"+observable).c_str()),isCutAndCount);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_"+observable).c_str()),isCutAndCount);
    addTemplate("Top_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("tbkghistwen_"+observable).c_str()),isCutAndCount);
    addTemplate("QCD_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("qbkghistwen_"+observable).c_str()),isCutAndCount);
    addTemplate("Dibosons_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("dbkghistwen_"+observable).c_str()),isCutAndCount);
    addTemplate("EWKZ_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("ewkzbkghistwen_"+observable).c_str()),isCutAndCount);
    addTemplate("EWKW_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("ewkwbkghistwen_"+observable).c_str()),isCutAndCount);

    if(addShapeSystematics){
      addShapeVariations("vllbkghistwen","ZJets_WE",suffix,observable,vars,*wspace_WE,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("tbkghistwen","Top_WE",suffix,observable,vars,*wspace_WE,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("dbkghistwen","Dibosons_WE",suffix,observable,vars,*wspace_WE,templatesfile,"",isCombination,isCutAndCount);      
      addShapeVariations("ewkzbkghistwen","EWKZ_WE",suffix,observable,vars,*wspace_WE,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("ewkwbkghistwen","EWKW_WE",suffix,observable,vars,*wspace_WE,templatesfile,"",isCombination,isCutAndCount);
    }

    if(not isCutAndCount){
      // connected W->enu with W+jets SR 
      vector<pair<RooRealVar*,TH1*> > wln_WE_syst;
      makeConnectedBinList("WJets_WE_"+suffix,met,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("wencorhist_"+observable).c_str()),wln_WE_syst,wln_SR_bins,NULL,observable);

    }
    else{

      vector<pair<RooRealVar*,systematicCutAndCount> > wln_WE_syst;
      makeConnectedBinListCutAndCount("WJets_WE_"+suffix,met,*wspace_WE,
                                      (TH1F*)templatesfile->FindObjectAny(("nhist_wen_"+observable).c_str()),
                                      (TH1F*)templatesfile->FindObjectAny(("dhist_wen_"+observable).c_str()),
                                      wln_WE_syst,wln_SR_bins,NULL,observable);

    }

  }
  else{

    cout<<"Make CR Single-Lepton templates ..."<<endl;
    
    wspace_WL = new RooWorkspace(("WL_"+suffix).c_str(),("WL_"+suffix).c_str());
    
    addTemplate("data_obs_WL_"+suffix,vars,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("datahistwln_"+observable).c_str()),isCutAndCount);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_WL_"+suffix,vars,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("vllbkghistwln_"+observable).c_str()),isCutAndCount);
    addTemplate("Top_WL_"+suffix,vars,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("tbkghistwln_"+observable).c_str()),isCutAndCount);
    addTemplate("QCD_WL_"+suffix,vars,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("qbkghistwln_"+observable).c_str()),isCutAndCount);
    addTemplate("Dibosons_WL_"+suffix,vars,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("dbkghistwln_"+observable).c_str()),isCutAndCount);
    addTemplate("EWKZ_WL_"+suffix,vars,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("ewkzbkghistwln_"+observable).c_str()),isCutAndCount);
    addTemplate("EWKW_WL_"+suffix,vars,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("ewkwbkghistwln_"+observable).c_str()),isCutAndCount);
    
    if(addShapeSystematics){
      addShapeVariations("vllbkghistwln","ZJets_WL",suffix,observable,vars,*wspace_WL,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("tbkghistwln","Top_WL",suffix,observable,vars,*wspace_WL,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("dbkghistwln","Dibosons_WL",suffix,observable,vars,*wspace_WL,templatesfile,"",isCombination,isCutAndCount);      
      addShapeVariations("ewkzbkghistwln","EWKZ_WL",suffix,observable,vars,*wspace_WL,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("ewkwbkghistwln","EWKW_WL",suffix,observable,vars,*wspace_WL,templatesfile,"",isCombination,isCutAndCount);
    }

    if(not isCutAndCount){
      // connected W->enu with W+jets SR 
      vector<pair<RooRealVar*,TH1*> > wln_WL_syst;
      makeConnectedBinList("WJets_WL_"+suffix,met,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("wlncorhist_"+observable).c_str()),wln_WL_syst,wln_SR_bins,NULL,observable);
    }
    else{
      vector<pair<RooRealVar*,systematicCutAndCount> > wln_WL_syst;
      makeConnectedBinListCutAndCount("WJets_WL_"+suffix,met,*wspace_WL,
                                      (TH1F*)templatesfile->FindObjectAny(("nhist_wln_"+observable).c_str()),
                                      (TH1F*)templatesfile->FindObjectAny(("dhist_wln_"+observable).c_str()),
                                      wln_WL_syst,wln_SR_bins,NULL,observable);
      
    }
    
  }


  ///////////////////////////////////////
  // -------- CR Gamma+jets  -------- //
  //////////////////////////////////////
  cout<<"Make CR Gamma+jets  templates ..."<<endl;
  RooWorkspace wspace_GJ(("GJ_"+suffix).c_str(),("GJ_"+suffix).c_str());
  
  addTemplate("data_obs_GJ_"+suffix,vars,wspace_GJ,(TH1F*)templatesfile->FindObjectAny(("datahistgam_"+observable).c_str()),isCutAndCount);
  // Gamma+jets --> connected with Z->nunu

  RooRealVar* znn_GJ_re1 = 0;
  RooRealVar* znn_GJ_fa1 = 0;
  RooRealVar* znn_GJ_re2 = 0;
  RooRealVar* znn_GJ_fa2 = 0;
  RooRealVar* znn_GJ_pdf = 0;
  RooRealVar* znn_GJ_fpc = 0;
  if(not isCombination){
    znn_GJ_re1 = new RooRealVar("Znunu_GJ_RenScale1"  ,"",0.,-5.,5.);
    znn_GJ_fa1 = new RooRealVar("Znunu_GJ_FactScale1" ,"",0.,-5.,5.);
    znn_GJ_re2 = new RooRealVar("Znunu_GJ_RenScale2"  ,"",0.,-5.,5.);
    znn_GJ_fa2 = new RooRealVar("Znunu_GJ_FactScale2" ,"",0.,-5.,5.);
    znn_GJ_pdf = new RooRealVar("Znunu_GJ_PDF"        ,"",0.,-5.,5.);
    znn_GJ_fpc = new RooRealVar("Znunu_GJ_Footprint"  ,"",0.,-5.,5.);
  }
  else{
    znn_GJ_re1 = new RooRealVar("mr"  ,"",0.,-5.,5.);
    znn_GJ_re2 = new RooRealVar("mr2"  ,"",0.,-5.,5.);
    znn_GJ_fa1 = new RooRealVar("mf"  ,"",0.,-5.,5.);
    znn_GJ_fa2 = new RooRealVar("mf2"  ,"",0.,-5.,5.);
    znn_GJ_pdf = new RooRealVar("pdf"  ,"",0.,-5.,5.);
    znn_GJ_fpc = new RooRealVar("fp"  ,"",0.,-5.,5.);
  }

  if(not isCutAndCount){

  vector<pair<RooRealVar*,TH1*> > znn_GJ_syst;

  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(NULL      ,(TH1F*)templatesfile->FindObjectAny(("ZG_EWK_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(znn_GJ_re1,(TH1F*)templatesfile->FindObjectAny(("ZG_RenScale1_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(znn_GJ_fa1,(TH1F*)templatesfile->FindObjectAny(("ZG_FactScale1_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(znn_GJ_re2,(TH1F*)templatesfile->FindObjectAny(("ZG_RenScale2_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(znn_GJ_fa2,(TH1F*)templatesfile->FindObjectAny(("ZG_FactScale2_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(znn_GJ_pdf,(TH1F*)templatesfile->FindObjectAny(("ZG_PDF_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(znn_GJ_fpc,(TH1F*)templatesfile->FindObjectAny(("ZG_Footprint_"+observable).c_str())));

  makeConnectedBinList("Znunu_GJ_"+suffix,met,wspace_GJ,(TH1F*)templatesfile->FindObjectAny(("gamcorewkhist_"+observable).c_str()),znn_GJ_syst,znn_SR_bins,NULL,observable);  

  }
  else{

      vector<pair<RooRealVar*,systematicCutAndCount> > znn_GJ_syst;
      
      RooRealVar* znn_GJ_ewk = new RooRealVar(("Znunu_GJ_"+suffix+"_ZG_EWK").c_str(),""  ,0.,-5.,5.);
      systematicCutAndCount znn_GJ_ewk_sys; // single bin everything here
      znn_GJ_ewk_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_ewk_"+observable).c_str());
      znn_GJ_ewk_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_ewk_"+observable).c_str());
      znn_GJ_ewk_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_qcd_"+observable).c_str());
      znn_GJ_ewk_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_qcd_"+observable).c_str());
      znn_GJ_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(znn_GJ_ewk,znn_GJ_ewk_sys));
      
      systematicCutAndCount znn_GJ_re1_sys;
      znn_GJ_re1_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_re1_"+observable).c_str());
      znn_GJ_re1_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_re1_"+observable).c_str());
      znn_GJ_re1_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_qcd_"+observable).c_str());
      znn_GJ_re1_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_qcd_"+observable).c_str());
      znn_GJ_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(znn_GJ_re1,znn_GJ_re1_sys));
      
      systematicCutAndCount znn_GJ_re2_sys;
      znn_GJ_re2_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_re2_"+observable).c_str());
      znn_GJ_re2_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_re2_"+observable).c_str());
      znn_GJ_re2_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_qcd_"+observable).c_str());
      znn_GJ_re2_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_qcd_"+observable).c_str());
      znn_GJ_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(znn_GJ_re2,znn_GJ_re2_sys));
      
      systematicCutAndCount znn_GJ_fa1_sys;
      znn_GJ_fa1_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_fa1_"+observable).c_str());
      znn_GJ_fa1_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_fa1_"+observable).c_str());
      znn_GJ_fa1_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_qcd_"+observable).c_str());
      znn_GJ_fa1_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_qcd_"+observable).c_str());
      znn_GJ_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(znn_GJ_fa1,znn_GJ_fa1_sys));
      
      systematicCutAndCount znn_GJ_fa2_sys;
      znn_GJ_fa2_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_fa2_"+observable).c_str());
      znn_GJ_fa2_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_fa2_"+observable).c_str());
      znn_GJ_fa2_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_qcd_"+observable).c_str());
      znn_GJ_fa2_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_qcd_"+observable).c_str());
      znn_GJ_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(znn_GJ_fa2,znn_GJ_fa2_sys));
    
      systematicCutAndCount znn_GJ_pdf_sys;
      znn_GJ_pdf_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_pdf_"+observable).c_str());
      znn_GJ_pdf_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_pdf_"+observable).c_str());
      znn_GJ_pdf_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_qcd_"+observable).c_str());
      znn_GJ_pdf_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_qcd_"+observable).c_str());
      znn_GJ_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(znn_GJ_pdf,znn_GJ_pdf_sys));

      systematicCutAndCount znn_GJ_fpc_sys;
      znn_GJ_fpc_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_fpc_"+observable).c_str());
      znn_GJ_fpc_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_fpc_"+observable).c_str());
      znn_GJ_fpc_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_gam_qcd_"+observable).c_str());
      znn_GJ_fpc_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_gam_qcd_"+observable).c_str());
      znn_GJ_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(znn_GJ_fpc,znn_GJ_fpc_sys));
      
      makeConnectedBinListCutAndCount("Znunu_GJ_"+suffix,met,wspace_GJ,
				      (TH1F*)templatesfile->FindObjectAny(("nhist_gam_ewk_"+observable).c_str()),
				      (TH1F*)templatesfile->FindObjectAny(("dhist_gam_ewk_"+observable).c_str()),
				      znn_GJ_syst,znn_SR_bins,NULL,observable);

  }
  
  // Other MC backgrounds photon+jets control region
  addTemplate("QCD_GJ_"+suffix,vars,wspace_GJ,(TH1F*)templatesfile->FindObjectAny(("qbkghistgam_"+observable).c_str()),isCutAndCount);

  RooWorkspace* wspace_TM = NULL;
  RooWorkspace* wspace_TE = NULL;

  if(connectTop){
    
    /////////////////////////////////////
    // -------- CR Top-Muon  -------- //
    ////////////////////////////////////
    cout<<"Make CR Top-mu  templates ..."<<endl;
    wspace_TM = new RooWorkspace(("TM_"+suffix).c_str(),("TM_"+suffix).c_str());

    addTemplate("data_obs_TM_"+suffix,vars,*wspace_TM,(TH1F*)templatesfile->FindObjectAny(("datahisttopmu_"+observable).c_str()),isCutAndCount);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_TM_"+suffix,vars,*wspace_TM,(TH1F*)templatesfile->FindObjectAny(("vllbkghisttopmu_"+observable).c_str()),isCutAndCount);
    addTemplate("WJets_TM_"+suffix,vars,*wspace_TM,(TH1F*)templatesfile->FindObjectAny(("vlbkghisttopmu_"+observable).c_str()),isCutAndCount);
    addTemplate("QCD_TM_"+suffix,vars,*wspace_TM,(TH1F*)templatesfile->FindObjectAny(("qbkghisttopmu_"+observable).c_str()),isCutAndCount);
    addTemplate("Dibosons_TM_"+suffix,vars,*wspace_TM,(TH1F*)templatesfile->FindObjectAny(("dbkghisttopmu_"+observable).c_str()),isCutAndCount);

    if(addShapeSystematics){      
      addShapeVariations("vllbkghisttopmu","ZJets_TM",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("vlbkghisttopmu","WJets_TM",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("dbkghisttopmu","Dibosons_TM",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination,isCutAndCount);
    }
    
    if(not isCutAndCount){

      // connect tt->mu+b with tt SR
      vector<pair<RooRealVar*,TH1*> > top_TM_syst;
      RooRealVar* top_btag   = new RooRealVar("Top_btag","",0.,-5.,5.);
      top_TM_syst.push_back(pair<RooRealVar*,TH1*>(top_btag,(TH1F*)templatesfile->FindObjectAny(("TOP_MU_B_"+observable).c_str())));    
      makeConnectedBinList("Top_TM_"+suffix,met,*wspace_TM,(TH1F*)templatesfile->FindObjectAny(("topmucorhist_"+observable).c_str()),top_TM_syst,top_SR_bins,NULL,observable);

    }
    else{
      vector<pair<RooRealVar*,systematicCutAndCount> > top_TM_syst;
      RooRealVar* top_btag   = new RooRealVar("Top_btag","",0.,-5.,5.);
      systematicCutAndCount top_TM_btag_sys;
      top_TM_btag_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_topmu_btagup_"+observable).c_str());
      top_TM_btag_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_topmu_btagup_"+observable).c_str());
      top_TM_btag_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_topmu_"+observable).c_str());
      top_TM_btag_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_topmu_"+observable).c_str());
      top_TM_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(top_btag,top_TM_btag_sys));

      makeConnectedBinListCutAndCount("TOP_TM_"+suffix,met,*wspace_TM,
                                      (TH1F*)templatesfile->FindObjectAny(("nhist_topmu_"+observable).c_str()),
                                      (TH1F*)templatesfile->FindObjectAny(("dhist_topmu_"+observable).c_str()),
                                      top_TM_syst,top_SR_bins,NULL,observable);
      
    }


    /////////////////////////////////////
    // -------- CR Top-Electron-------- //
    ////////////////////////////////////
    cout<<"Make CR Top-el  templates ..."<<endl;
    wspace_TE = new RooWorkspace(("TE_"+suffix).c_str(),("TE_"+suffix).c_str());

    addTemplate("data_obs_TE_"+suffix,vars,*wspace_TE,(TH1F*)templatesfile->FindObjectAny(("datahisttopel_"+observable).c_str()),isCutAndCount);

    if(not isCutAndCount){
      
      // connect tt->mu+b with tt SR
      vector<pair<RooRealVar*,TH1*> > top_TE_syst;
      RooRealVar* top_btag   = new RooRealVar("Top_btag","",0.,-5.,5.);
      top_TE_syst.push_back(pair<RooRealVar*,TH1*>(top_btag,(TH1F*)templatesfile->FindObjectAny(("TOP_MU_B_"+observable).c_str())));    
      makeConnectedBinList("Top_TE_"+suffix,met,*wspace_TE,(TH1F*)templatesfile->FindObjectAny(("topelcorhist_"+observable).c_str()),top_TE_syst,top_SR_bins,NULL,observable);

    }
    else{
      vector<pair<RooRealVar*,systematicCutAndCount> > top_TE_syst;
      RooRealVar* top_btag   = new RooRealVar("Top_btag","",0.,-5.,5.);
      systematicCutAndCount top_TE_btag_sys;
      top_TE_btag_sys.num_1 = (TH1F*)templatesfile->FindObjectAny(("nhist_topel_btagup_"+observable).c_str());
      top_TE_btag_sys.den_1 = (TH1F*)templatesfile->FindObjectAny(("dhist_topel_btagup_"+observable).c_str());
      top_TE_btag_sys.num_2 = (TH1F*)templatesfile->FindObjectAny(("nhist_topel_"+observable).c_str());
      top_TE_btag_sys.den_2 = (TH1F*)templatesfile->FindObjectAny(("dhist_topel_"+observable).c_str());
      top_TE_syst.push_back(pair<RooRealVar*,systematicCutAndCount>(top_btag,top_TE_btag_sys));

      makeConnectedBinListCutAndCount("TOP_TE_"+suffix,met,*wspace_TE,
                                      (TH1F*)templatesfile->FindObjectAny(("nhist_topel_"+observable).c_str()),
                                      (TH1F*)templatesfile->FindObjectAny(("dhist_topel_"+observable).c_str()),
                                      top_TE_syst,top_SR_bins,NULL,observable);
      
    }
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_TE_"+suffix ,vars,*wspace_TE,(TH1F*)templatesfile->FindObjectAny(("vllbkghisttopel_"+observable).c_str()),isCutAndCount);
    addTemplate("WJets_TE_"+suffix ,vars,*wspace_TE,(TH1F*)templatesfile->FindObjectAny(("vlbkghisttopel_"+observable).c_str()),isCutAndCount);
    addTemplate("QCD_TE_"+suffix   ,vars,*wspace_TE,(TH1F*)templatesfile->FindObjectAny(("qbkghisttopel_"+observable).c_str()),isCutAndCount);
    addTemplate("Dibosons_TE_"+suffix,vars,*wspace_TE,(TH1F*)templatesfile->FindObjectAny(("dbkghisttopel_"+observable).c_str()),isCutAndCount);

    if(addShapeSystematics){
      addShapeVariations("vllbkghisttopel","ZJets_TE",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("vlbkghisttopel","WJets_TE",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination,isCutAndCount);
      addShapeVariations("dbkghisttopel","Dibosons_TE",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination,isCutAndCount);                 
    }
    
  }
  
  // ---------------------------- Write out the workspace -----------------------------------------------------------------//
  outfile->cd();
  wspace_SR.Write();  
  if(not mergeLeptons){
    wspace_ZM->Write();
    wspace_ZE->Write();
    wspace_WM->Write();
    wspace_WE->Write();
  }
  else{
    wspace_ZL->Write();
    wspace_WL->Write();
  }

  wspace_GJ.Write();
  if(connectTop){
    wspace_TM->Write();
    wspace_TE->Write();
  }
  
  outfile->Close();
  return;  
}

//  LocalWords:  isCutAndCount
