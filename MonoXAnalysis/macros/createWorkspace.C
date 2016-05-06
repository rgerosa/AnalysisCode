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
		     int category,  // analysis category
		     string outputName    = "workspace.root", // output workspace name
		     string observable    = "met",    // observable 1D or 2D
		     bool   isHiggsInvisible = false, // Higgs invisible or DM analsysis
		     float  scaleQCD      = 2,    // scale for QCD MC in SR
		     bool   connectWZ     = true, // connect W and Z in SR 
		     bool   connectTop    = false, // connect top CR --> fit also top eneriched samples
		     bool   addShapeSystematics = false, // add shapeN2
		     bool   mergeLeptons  = false, // merge mm and ee final sates
		     bool   isCombination = false, // convention for HIG-16 invisible combo
		     string interaction   = "Vector", // DM interaction 
		     string mediatorMass  = "1000",   // Med mass
		     string DMMass        = "50"      // DM mass
){
  
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  RooMsgService::instance().setSilentMode(kTRUE); 
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  // for templates and sys naming
  string suffix;
  if(category <=1) suffix = "MJ";
  else suffix = "MV";
    
  // create the output workspace
  cout<<"Create output file ..."<<endl;
  TFile *outfile = new TFile(outputName.c_str(),"RECREATE");
  outfile->cd();
  cout<<"Load binning and observable ..."<<endl;

  // Select observable and binning
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
  addTemplate("data_obs_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("datahist_"+observable).c_str()));

  // Signal shape
  if(!isHiggsInvisible){

    addTemplate("MonoJ_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
    addTemplate("MonoW_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
    addTemplate("MonoZ_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));

    if(addShapeSystematics){

      addShapeVariations("monoJhist","MonoJ_SR",suffix,observable,vars,wspace_SR,templatesfile,interaction+"_"+mediatorMass+"_"+DMMass,isCombination);
      addShapeVariations("monoWhist","MonoW_SR",suffix,observable,vars,wspace_SR,templatesfile,interaction+"_"+mediatorMass+"_"+DMMass,isCombination);
      addShapeVariations("monoZhist","MonoZ_SR",suffix,observable,vars,wspace_SR,templatesfile,interaction+"_"+mediatorMass+"_"+DMMass,isCombination);

      // statistics
      generateStatTemplate("MonoJ_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
      generateStatTemplate("MonoW_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
      generateStatTemplate("MonoZ_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
    }
  }
  else{

    addTemplate("ggH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+observable).c_str()));
    addTemplate("qqH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("vbfHhist_"+mediatorMass+"_"+observable).c_str()));
    addTemplate("WH_SR_"+suffix, vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("wHhist_"+mediatorMass+"_"+observable).c_str()));
    addTemplate("ZH_SR_"+suffix, vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("zHhist_"+mediatorMass+"_"+observable).c_str()));
    addTemplate("ggZH_SR_"+suffix, vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("ggzHhist_"+mediatorMass+"_"+observable).c_str()));

    if(addShapeSystematics){

      addShapeVariations("ggHhist","ggH_SR",suffix,observable,vars,wspace_SR,templatesfile,mediatorMass,isCombination);
      addShapeVariations("vbfHhist","qqH_SR",suffix,observable,vars,wspace_SR,templatesfile,mediatorMass,isCombination);
      addShapeVariations("wHhist","WH_SR",suffix,observable,vars,wspace_SR,templatesfile,mediatorMass,isCombination);
      addShapeVariations("zHhist","ZH_SR",suffix,observable,vars,wspace_SR,templatesfile,mediatorMass,isCombination);
      addShapeVariations("ggZHhist","ggZH_SR",suffix,observable,vars,wspace_SR,templatesfile,mediatorMass,isCombination);

      // ggH higgs pT uncertainties 
      TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoRenUp = (TH1F*)templatesfile->FindObjectAny(("ggHhist_renUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoRenDw = (TH1F*)templatesfile->FindObjectAny(("ggHhist_renDw_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoFacUp = (TH1F*)templatesfile->FindObjectAny(("ggHhist_facUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoFacDw = (TH1F*)templatesfile->FindObjectAny(("ggHhist_facDw_"+mediatorMass+"_"+observable).c_str());

      if(nominalHisto){
	vector<TH1F*> histoVec;
	histoVec.push_back(histoRenUp);
	histoVec.push_back(histoRenDw);
	histoVec.push_back(histoFacUp);
	histoVec.push_back(histoFacDw);		
	addTemplate("ggH_SR_"+suffix+"_hptUp",vars,wspace_SR,generateEnvelopeMax(histoVec));
	addTemplate("ggH_SR_"+suffix+"_hptDown",vars,wspace_SR,generateEnvelopeMin(histoVec));
      }
      
      // statistics
      generateStatTemplate("ggH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
      generateStatTemplate("qqH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("vbfHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
      generateStatTemplate("WH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("wHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
      generateStatTemplate("ZH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("zHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
      generateStatTemplate("ggZH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("zHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
    }
  }
  
  // Zvv background --> to be extracted from CRs
  TH1F* znn_SR_hist = (TH1F*) templatesfile->FindObjectAny(("zinvhist_"+observable).c_str());
  RooArgList znn_SR_bins; 
  // create a RooParametric hist with one RooRealVar per bin 
  makeBinList("Znunu_SR_"+suffix,met,wspace_SR,znn_SR_hist,znn_SR_bins);

  // Top background --> to be extracted from CRs
  RooArgList top_SR_bins;
  TH1F* top_SR_hist = NULL;

  // for data driven top estimation
  if(connectTop){
    top_SR_hist = (TH1F*) templatesfile->FindObjectAny(("tbkghist_"+observable).c_str());
    RooArgList top_SR_bins; 
    makeBinList("Top_SR_"+suffix,met,wspace_SR,top_SR_hist,top_SR_bins);
  }
  else{ // rely on MC + systematics

    addTemplate("Top_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("tbkghist_"+observable).c_str()));

    if(addShapeSystematics){
      addShapeVariations("tbkghist","Top_SR",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      generateStatTemplate("Top_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("tbkghist_"+observable).c_str()));
    }
  }

  // WJets background --> to be extracted from CRs,with connection to Z->nunu
  TH1F* wln_SR_hist = (TH1F*) templatesfile->FindObjectAny(("wjethist_"+observable).c_str());
  RooArgList wln_SR_bins;
  // set of correlated systematic uncertainties for the Z/W ratio
  vector<pair<RooRealVar*,TH1*> > wln_SR_syst;

  RooRealVar* wln_SR_re1 = new RooRealVar("WJets_SR_RenScale1",""  ,0.,-5.,5.);
  RooRealVar* wln_SR_fa1 = new RooRealVar("WJets_SR_FactScale1","" ,0.,-5.,5.);
  RooRealVar* wln_SR_re2 = new RooRealVar("WJets_SR_RenScale2",""  ,0.,-5.,5.);
  RooRealVar* wln_SR_fa2 = new RooRealVar("WJets_SR_FactScale2","" ,0.,-5.,5.);
  RooRealVar* wln_SR_pdf = new RooRealVar("WJets_SR_PDF",""        ,0.,-5.,5.);

  // NULL means bin-by-bin
  wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(NULL      ,(TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_re1,(TH1F*)templatesfile->FindObjectAny(("ZW_RenScale1_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_fa1,(TH1F*)templatesfile->FindObjectAny(("ZW_FactScale1_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_re2,(TH1F*)templatesfile->FindObjectAny(("ZW_RenScale2_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_fa2,(TH1F*)templatesfile->FindObjectAny(("ZW_FactScale2_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*,TH1*>(wln_SR_pdf,(TH1F*)templatesfile->FindObjectAny(("ZW_PDF_"+observable).c_str())));

  if (!connectWZ) 
    makeBinList("WJets_SR_"+suffix,met,wspace_SR,wln_SR_hist,wln_SR_bins);
  else   
    makeConnectedBinList("WJets_SR_"+suffix,met,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("zwjcorewkhist_"+observable).c_str()),wln_SR_syst,znn_SR_bins,&wln_SR_bins,observable);

  // Other MC backgrounds
  addTemplate("ZJets_SR_"+suffix     ,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("zjethist_"+observable).c_str()));
  addTemplate("Dibosons_SR_"+suffix  ,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("dbkghist_"+observable).c_str()));
  addTemplate("GJets_SR_"+suffix     ,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("gbkghist_"+observable).c_str()));

  if(addShapeSystematics){
    addShapeVariations("zjethist","ZJets_SR",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
    generateStatTemplate("ZJets_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("zjethist_"+observable).c_str()));

    addShapeVariations("dbkghist","Dibosons_SR",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
    generateStatTemplate("Dibosons_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("dbkghist_"+observable).c_str()));

    addShapeVariations("gbkghist","GJets_SR",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
    generateStatTemplate("GJets_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("gbkghist_"+observable).c_str()));
  }

  // look for DD qcd background,otherwise MC scaled by a factor 2
  TH1F* qcdhist = (TH1F*)templatesfile->FindObjectAny(("qbkghistDD_"+observable).c_str());
  if(qcdhist){
    addTemplate("QCD_SR_"+suffix,vars,wspace_SR,qcdhist);
    addTemplate("QCD_SR_"+suffix+"_CMS_QCD_SRUp",vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("qbkghistDD_shapeUp_"+observable).c_str()));
    addTemplate("QCD_SR_"+suffix+"_CMS_QCD_SRDown",vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("qbkghistDD_shapeDw_"+observable).c_str()));
  }
  else{
    qcdhist = (TH1F*)templatesfile->FindObjectAny(("qbkghist_"+observable).c_str());
    qcdhist->Scale(scaleQCD);  
    addTemplate("QCD_SR_"+suffix,vars,wspace_SR,qcdhist);
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
  
    addTemplate("data_obs_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("datahistzmm_"+observable).c_str()));
    // Z->mumu connected with Z->nunu SR
    vector<pair<RooRealVar*,TH1*> >   znn_ZM_syst;
  
    // Other MC backgrounds in dimuon control region
    addTemplate("WJets_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_"+observable).c_str()));
    addTemplate("Top_ZM_"+suffix ,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_"+observable).c_str()));
    addTemplate("QCD_ZM_"+suffix ,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("qbkghistzmm_"+observable).c_str()));
    addTemplate("Dibosons_ZM_"+suffix,vars,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_"+observable).c_str()));

    if(addShapeSystematics){      
      addShapeVariations("vlbkghist","WJets_ZM",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("tbkghist","Top_ZM",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("dbkghist","Dibosons_ZM",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
    }

    makeConnectedBinList("Znunu_ZM_"+suffix,met,*wspace_ZM,(TH1F*)templatesfile->FindObjectAny(("zmmcorhist_"+observable).c_str()),znn_ZM_syst,znn_SR_bins,NULL,observable);
    
    ////////////////////////////////////////
    // -------- CR Di-Electron  -------- //
    ///////////////////////////////////////

    cout<<"Make CR Di-Electron  templates ..."<<endl;      
    wspace_ZE = new RooWorkspace(("ZE_"+suffix).c_str(),("ZE_"+suffix).c_str());

    addTemplate("data_obs_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("datahistzee_"+observable).c_str()));
    // Z->ee connected with Z->nunu SR
    vector<pair<RooRealVar*,TH1*> > znn_ZE_syst;
    makeConnectedBinList("Znunu_ZE_"+suffix,met,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("zeecorhist_"+observable).c_str()),znn_ZE_syst,znn_SR_bins,NULL,observable);
    
    // Other MC backgrounds in dielectron control region
    addTemplate("WJets_ZE_"+suffix,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_"+observable).c_str()));
    addTemplate("Top_ZE_"+suffix  ,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("tbkghistzee_"+observable).c_str()));
    addTemplate("QCD_ZE_"+suffix  ,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("qbkghistzee_"+observable).c_str()));
    addTemplate("Dibosons_ZE_"+suffix  ,vars,*wspace_ZE,(TH1F*)templatesfile->FindObjectAny(("dbkghistzee_"+observable).c_str()));

    if( addShapeSystematics){      
      addShapeVariations("vlbkghist","WJets_ZE",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("tbkghist","Top_ZE",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("dbkghist","Dibosons_ZE",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
    }
  }
  else{
    
    cout<<"Make CR Di-Lepton  templates ..."<<endl;    
    wspace_ZL = new RooWorkspace(("ZL_"+suffix).c_str(),("ZL_"+suffix).c_str());
  
    addTemplate("data_obs_ZL_"+suffix,vars,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("datahistzll_"+observable).c_str()));
    // Z->mumu connected with Z->nunu SR
    vector<pair<RooRealVar*,TH1*> >   znn_syst;
    makeConnectedBinList("Znunu_ZL_"+suffix,met,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("zllcorhist_"+observable).c_str()),znn_syst,znn_SR_bins,NULL,observable);
  
    // Other MC backgrounds in dimuon control region
    addTemplate("WJets_ZL_"+suffix,vars,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("vlbkghistzll_"+observable).c_str()));
    addTemplate("Top_ZL_"+suffix ,vars,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("tbkghistzll_"+observable).c_str()));
    addTemplate("QCD_ZL_"+suffix ,vars,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("qbkghistzll_"+observable).c_str()));
    addTemplate("Dibosons_ZL_"+suffix,vars,*wspace_ZL,(TH1F*)templatesfile->FindObjectAny(("dbkghistzll_"+observable).c_str()));

    if(addShapeSystematics){
      addShapeVariations("vlbkghist","WJets_ZL",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("tbkghist","Top_ZL",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("dbkghist","Dibosons_ZL",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
    }            
  }

  ///////////////////////////////////////
  // -------- CR Gamma+jets  -------- //
  //////////////////////////////////////
  cout<<"Make CR Gamma+jets  templates ..."<<endl;
  RooWorkspace wspace_GJ(("GJ_"+suffix).c_str(),("GJ_"+suffix).c_str());

  addTemplate("data_obs_GJ_"+suffix,vars,wspace_GJ,(TH1F*)templatesfile->FindObjectAny(("datahistgam_"+observable).c_str()));
  // Gamma+jets --> connected with Z->nunu
  vector<pair<RooRealVar*,TH1*> > znn_GJ_syst;
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

  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(NULL      ,(TH1F*)templatesfile->FindObjectAny(("ZG_EWK_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(znn_GJ_re1,(TH1F*)templatesfile->FindObjectAny(("ZG_RenScale1_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(znn_GJ_fa1,(TH1F*)templatesfile->FindObjectAny(("ZG_FactScale1_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(znn_GJ_re2,(TH1F*)templatesfile->FindObjectAny(("ZG_RenScale2_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(znn_GJ_fa2,(TH1F*)templatesfile->FindObjectAny(("ZG_FactScale2_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(znn_GJ_pdf,(TH1F*)templatesfile->FindObjectAny(("ZG_PDF_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*,TH1*>(znn_GJ_fpc,(TH1F*)templatesfile->FindObjectAny(("ZG_Footprint_"+observable).c_str())));

  makeConnectedBinList("Znunu_GJ_"+suffix,met,wspace_GJ,(TH1F*)templatesfile->FindObjectAny(("gamcorewkhist_"+observable).c_str()),znn_GJ_syst,znn_SR_bins,NULL,observable);  
  // Other MC backgrounds photon+jets control region
  addTemplate("QCD_GJ_"+suffix,vars,wspace_GJ,(TH1F*)templatesfile->FindObjectAny(("qbkghistgam_"+observable).c_str()));

  
  if(not mergeLeptons){
    ///////////////////////////////////////
    // -------- CR Single-Muon  -------- //
    //////////////////////////////////////
    cout<<"Make CR Single-Mu  templates ..."<<endl;
    
    wspace_WM = new RooWorkspace(("WM_"+suffix).c_str(),("WM_"+suffix).c_str());

    addTemplate("data_obs_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("datahistwmn_"+observable).c_str()));
    // connected W->munu with W+jets SR
    vector<pair<RooRealVar*,TH1*> > wln_WM_syst;
    makeConnectedBinList("WJets_WM_"+suffix,met,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("wmncorhist_"+observable).c_str()),wln_WM_syst,wln_SR_bins,NULL,observable);
    
  // Other MC backgrounds in single muon control region
    addTemplate("ZJets_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_"+observable).c_str()));
    addTemplate("Top_WM_"+suffix  ,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_"+observable).c_str()));
    addTemplate("QCD_WM_"+suffix  ,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("qbkghistwmn_"+observable).c_str()));
    addTemplate("Dibosons_WM_"+suffix,vars,*wspace_WM,(TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_"+observable).c_str()));
    
    if(addShapeSystematics){      
      addShapeVariations("vllbkghist","ZJets_WM",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("tbkghist","Top_WM",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("dbkghist","Dibosons_WM",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
    }
  
    //////////////////////////////////..../////
    // -------- CR Single-Electron  -------- //
    //////////////////////////////////////////
    cout<<"Make CR Single-El  templates ..."<<endl;
    
    wspace_WE = new RooWorkspace(("WE_"+suffix).c_str(),("WE_"+suffix).c_str());
    
    addTemplate("data_obs_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("datahistwen_"+observable).c_str()));
    // connected W->enu with W+jets SR 
    vector<pair<RooRealVar*,TH1*> > wln_WE_syst;
    makeConnectedBinList("WJets_WE_"+suffix,met,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("wencorhist_"+observable).c_str()),wln_WE_syst,wln_SR_bins,NULL,observable);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_"+observable).c_str()));
    addTemplate("Top_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("tbkghistwen_"+observable).c_str()));
    addTemplate("QCD_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("qbkghistwen_"+observable).c_str()));
    addTemplate("Dibosons_WE_"+suffix,vars,*wspace_WE,(TH1F*)templatesfile->FindObjectAny(("dbkghistwen_"+observable).c_str()));
    
    if(addShapeSystematics){
      addShapeVariations("vllbkghist","ZJets_WE",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("tbkghist","Top_WE",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("dbkghist","Dibosons_WE",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);      
    }
  }
  else{

    cout<<"Make CR Single-Lepton templates ..."<<endl;
    
    wspace_WL = new RooWorkspace(("WL_"+suffix).c_str(),("WL_"+suffix).c_str());
    
    addTemplate("data_obs_WL_"+suffix,vars,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("datahistwln_"+observable).c_str()));
    // connected W->enu with W+jets SR 
    vector<pair<RooRealVar*,TH1*> > wln_WL_syst;
    makeConnectedBinList("WJets_WL_"+suffix,met,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("wlncorhist_"+observable).c_str()),wln_WL_syst,wln_SR_bins,NULL,observable);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_WL_"+suffix,vars,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("vllbkghistwln_"+observable).c_str()));
    addTemplate("Top_WL_"+suffix,vars,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("tbkghistwln_"+observable).c_str()));
    addTemplate("QCD_WL_"+suffix,vars,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("qbkghistwln_"+observable).c_str()));
    addTemplate("Dibosons_WL_"+suffix,vars,*wspace_WL,(TH1F*)templatesfile->FindObjectAny(("dbkghistwln_"+observable).c_str()));
    
    if(addShapeSystematics){
      addShapeVariations("vllbkghist","ZJets_WL",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("tbkghist","Top_WL",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("dbkghist","Dibosons_WL",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);      
    }
  }

  RooWorkspace* wspace_TM = NULL;
  RooWorkspace* wspace_TE = NULL;

  if(connectTop){
    
    /////////////////////////////////////
    // -------- CR Top-Muon  -------- //
    ////////////////////////////////////
    cout<<"Make CR Top-mu  templates ..."<<endl;
    wspace_TM = new RooWorkspace(("TM_"+suffix).c_str(),("TM_"+suffix).c_str());

    addTemplate("data_obs_TM_"+suffix,vars,*wspace_TM,(TH1F*)templatesfile->FindObjectAny(("datahisttopmu_"+observable).c_str()));
    // connect tt->mu+b with tt SR
    vector<pair<RooRealVar*,TH1*> > top_TM_syst;
    RooRealVar* top_btag   = new RooRealVar("Top_btag","",0.,-5.,5.);
    top_TM_syst.push_back(pair<RooRealVar*,TH1*>(top_btag,(TH1F*)templatesfile->FindObjectAny(("TOP_MU_B_"+observable).c_str())));
    
    makeConnectedBinList("Top_TM_"+suffix,met,*wspace_TM,(TH1F*)templatesfile->FindObjectAny(("topmucorhist_"+observable).c_str()),top_TM_syst,top_SR_bins,NULL,observable);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_TM_"+suffix,vars,*wspace_TM,(TH1F*)templatesfile->FindObjectAny(("vllbkghisttopmu_"+observable).c_str()));
    addTemplate("WJets_TM_"+suffix,vars,*wspace_TM,(TH1F*)templatesfile->FindObjectAny(("vlbkghisttopmu_"+observable).c_str()));
    addTemplate("QCD_TM_"+suffix,vars,*wspace_TM,(TH1F*)templatesfile->FindObjectAny(("qbkghisttopmu_"+observable).c_str()));
    addTemplate("Dibosons_TM_"+suffix,vars,*wspace_TM,(TH1F*)templatesfile->FindObjectAny(("dbkghisttopmu_"+observable).c_str()));

    if(addShapeSystematics){      
      addShapeVariations("vllbkghist","ZJets_TM",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("vlbkghist","WJets_TM",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("dbkghist","Dibosons_TM",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
    }

    /////////////////////////////////////
    // -------- CR Top-Electron-------- //
    ////////////////////////////////////
    cout<<"Make CR Top-el  templates ..."<<endl;
    wspace_TE = new RooWorkspace(("TE_"+suffix).c_str(),("TE_"+suffix).c_str());

    addTemplate("data_obs_TE_"+suffix,vars,*wspace_TE,(TH1F*)templatesfile->FindObjectAny(("datahisttopel_"+observable).c_str()));
    // connect tt->mu+b with tt SR
    vector<pair<RooRealVar*,TH1*> > top_TE_syst;
    top_TE_syst.push_back(pair<RooRealVar*,TH1*>(top_btag,(TH1F*)templatesfile->FindObjectAny(("TOP_EL_B_"+observable).c_str())));
    
    makeConnectedBinList("Top_TE_"+suffix,met,*wspace_TE,(TH1F*)templatesfile->FindObjectAny(("topelcorhist_"+observable).c_str()),top_TE_syst,top_SR_bins,NULL,observable);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_TE_"+suffix ,vars,*wspace_TE,(TH1F*)templatesfile->FindObjectAny(("vllbkghisttopel_"+observable).c_str()));
    addTemplate("WJets_TE_"+suffix ,vars,*wspace_TE,(TH1F*)templatesfile->FindObjectAny(("vlbkghisttopel_"+observable).c_str()));
    addTemplate("QCD_TE_"+suffix   ,vars,*wspace_TE,(TH1F*)templatesfile->FindObjectAny(("qbkghisttopel_"+observable).c_str()));
    addTemplate("Dibosons_TE_"+suffix,vars,*wspace_TE,(TH1F*)templatesfile->FindObjectAny(("dbkghisttopel_"+observable).c_str()));

    if(addShapeSystematics){
      addShapeVariations("vllbkghist","ZJets_TE",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("vlbkghist","WJets_TE",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);
      addShapeVariations("dbkghist","Dibosons_TE",suffix,observable,vars,wspace_SR,templatesfile,"",isCombination);                 
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
