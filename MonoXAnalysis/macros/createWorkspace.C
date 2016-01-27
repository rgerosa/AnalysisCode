#include <string>
#include <vector>
#include <utility>

using namespace std;

double defaultVal = 1e-10;

// function to create a RooDataHist from TH1F and import it in a workspace
void addTemplate(string procname, RooArgList& varlist, RooWorkspace& ws, TH1F* hist) {
  RooDataHist rhist(procname.c_str(), "", varlist, hist);
  ws.import(rhist);
}

// Make list of bins of a TH1F as RooArgList of RooRealVar and building the RooParametricHist (to be Run in the release with combine)
void makeBinList(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* hist, RooArgList& binlist, bool setConst = false) {

  // loop over histo bins
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    stringstream binss;
    binss << procname << "_bin" << i;
    RooRealVar* binvar;
    // make a RooRealVar for each bin [0,2*binContent]
    if (!setConst) 
      binvar = new RooRealVar(binss.str().c_str(), "", hist->GetBinContent(i), 0., hist->GetBinContent(i)*2.0);
    else           
      binvar = new RooRealVar(binss.str().c_str(), "", hist->GetBinContent(i));

    if(hist->GetBinContent(i) == 0)
      binvar->setVal(defaultVal);

    binlist.add(*binvar);
  }
  
  stringstream normss;
  normss << procname << "_norm";
  
  // build a RooParamtric hist using the bin list and the histogram
  RooParametricHist phist(procname.c_str(), "", var, binlist, *hist);
  RooAddition norm(normss.str().c_str(), "", binlist);

  ws.import(phist,RooFit::RecycleConflictNodes());
  ws.import(norm, RooFit::RecycleConflictNodes());
}

// make connections betweem signal region and control region
void makeConnectedBinList(string procname, RooRealVar& var, 
			  RooWorkspace& ws, TH1F* rhist, 
			  vector<pair<RooRealVar*, TH1*> > syst, const RooArgList& srbinlist, RooArgList* crbinlist=NULL) {
  
  // bin list for the CR
  if (crbinlist == NULL) 
    crbinlist = new RooArgList();

  // Loop on ratio hist
  for (int i = 1; i <= rhist->GetNbinsX(); i++) {
    stringstream rbinss;
    rbinss << "r_" << procname << "_bin" << i;
    // Fixed value for each bin of the ratio
    RooRealVar* rbinvar = new RooRealVar(rbinss.str().c_str(), "", rhist->GetBinContent(i));

    if(rhist->GetBinContent(i) == 0)
      rbinvar->setVal(defaultVal);
    
    // uncertainty histograms for systematics
    stringstream rerrbinss;
    rerrbinss << procname << "_bin" << i << "_Runc";
    // Nuisance for the Final fit for each bin (bin-by-bin unc)
    RooRealVar* rerrbinvar = new RooRealVar(rerrbinss.str().c_str(), "", 0., -5., 5.);
    
    stringstream binss;
    binss << procname << "_bin" << i;

    // list of bins
    RooArgList fobinlist;
    // signal region bins
    fobinlist.add(srbinlist[i-1]);
    // connection bins from rHisto (ratio histo)
    fobinlist.add(*rbinvar);
    // uncertainty
    fobinlist.add(*rerrbinvar);
    
    // bin [i] (CR) = bin [i] (SR) /( Rbin [i] *(1+RhistError[i]/Rbin[i]*Rbin_Err[i] )) --> statstical uncertainty
    double err_stat = 10e-6;
    if(rhist->GetBinContent(i) != 0)
      err_stat = rhist->GetBinError(i)/rhist->GetBinContent(i);

    stringstream formss;
    formss << "@0/";
    formss << "(";
    formss << "@1";
    formss << "*(1+" << err_stat << "*@2)";
    // systemaitc uncertainty
    for (int j = 0; j < syst.size(); j++) {
      stringstream systbinss;
      if (syst[j].first == NULL) { // add bin by bin
	systbinss << procname << "_bin" << i << "_" << syst[j].second->GetName();
	RooRealVar* systbinvar = new RooRealVar(systbinss.str().c_str(), "", 0., -5., 5.);
	// Add all the systeamtics as new Multiplicative Nuisance for each bin
	fobinlist.add(*systbinvar);
      }
      else 
	fobinlist.add(*syst[j].first);
      
      formss << "*(1+" << syst[j].second->GetBinContent(i) << "*@" << j+3 << ")";
    }
    formss << ")";
    
    // create a single RooFormulaVar
    RooFormulaVar* binvar = new RooFormulaVar(binss.str().c_str(), "", formss.str().c_str(), RooArgList(fobinlist));
    crbinlist->add(*binvar);
  }
  
  stringstream normss;
  normss << procname << "_norm";

  // Make the parametric Histograms for the control regions \mu*SR[i]/(R(...))
  RooParametricHist phist(procname.c_str(), "", var, *crbinlist, *rhist);
  RooAddition norm(normss.str().c_str(),"", *crbinlist);
  
  ws.import(phist,RooFit::RecycleConflictNodes());
  ws.import(norm, RooFit::RecycleConflictNodes());
}


// function to create workspace, to be run from a release which has the combine package
void createWorkspace(string inputName, int category, 
		     string outputName = "workspace.root", string observable = "met", float scaleQCD = 2, bool connectWZ = true){

  gSystem->Load("libHiggsAnalysisCombinedLimit.so");  

  // create the output workspace
  cout<<"Create output file ..."<<endl;
  TFile *outfile = new TFile(outputName.c_str(),"RECREATE");
  RooWorkspace wspace("w","w");

  cout<<"Load binning and observable ..."<<endl;

  // to be fix one combine will switch to ROOT 6
  double xMin = 0.;
  double xMax = 0.;
  if(observable == "met" && category <=1){
    xMin = 200.;
    xMax = 1090.;
  }
  else if(observable == "met" && category >1){
    xMin = 200.;
    xMax = 1000.;
  }
  else
    cout<<"Binning not implemented for the observable "<<observable<<" --> please define it "<<endl;
  
  RooRealVar met(observable.c_str(),"",xMin,xMax);
  RooArgList vars(met);

  // Templates
  cout<<"Open inputFile ..."<<endl;
  TFile* templatesfile = new TFile(inputName.c_str());

  ///////////////////////////////////////
  // -------- SIGNAL REGION  -------- //
  ///////////////////////////////////////
  cout<<"Make SR templates ..."<<endl;

  // Data
  addTemplate("data_obs_SR", 
	      vars, wspace, 
	      (TH1F*)templatesfile->Get(("datahist_"+observable).c_str()));

  // Signal shape
  addTemplate("MonoJ_SR", 
	      vars, wspace, 
	      (TH1F*)templatesfile->Get(("monoJhist_"+observable).c_str()));

  addTemplate("MonoW_SR", 
	      vars, wspace, 
	      (TH1F*)templatesfile->Get(("monoWhist_"+observable).c_str()));

  addTemplate("MonoZ_SR", 
	      vars, wspace, 
	      (TH1F*)templatesfile->Get(("monoZhist_"+observable).c_str()));

  // Zvv background --> to be extracted from CRs
  TH1F* znn_SR_hist = (TH1F*) templatesfile->Get(("zinvhist_"+observable).c_str());
  RooArgList znn_SR_bins; 
  // create a RooParametric hist with one RooRealVar per bin 
  makeBinList("Znunu_SR", met, wspace, znn_SR_hist, znn_SR_bins);

  // Top background --> to be extracted from CRs
  TH1F* top_SR_hist = (TH1F*) templatesfile->Get(("tbkghist_"+observable).c_str());
  RooArgList top_SR_bins; 
  makeBinList("Top_SR", met, wspace, top_SR_hist, top_SR_bins);

  // WJets background --> to be extracted from CRs, with connection to Z->nunu
  TH1F* wln_SR_hist = (TH1F*) templatesfile->Get(("wjethist_"+observable).c_str());
  RooArgList wln_SR_bins;
  // set of correlated systematic uncertainties for the Z/W ratio
  vector<pair<RooRealVar*, TH1*> > wln_SR_syst;
  RooRealVar* wln_SR_re1 = new RooRealVar("WJets_SR_RenScale1" , "", 0., -5., 5.);
  RooRealVar* wln_SR_fa1 = new RooRealVar("WJets_SR_FactScale1", "", 0., -5., 5.);
  RooRealVar* wln_SR_re2 = new RooRealVar("WJets_SR_RenScale2" , "", 0., -5., 5.);
  RooRealVar* wln_SR_fa2 = new RooRealVar("WJets_SR_FactScale2", "", 0., -5., 5.);
  RooRealVar* wln_SR_pdf = new RooRealVar("WJets_SR_PDF"       , "", 0., -5., 5.);
  // NULL means bin-by-bin
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(NULL      , (TH1F*)templatesfile->Get("ZW_EWK")));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_re1, (TH1F*)templatesfile->Get("ZW_RenScale1")));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_fa1, (TH1F*)templatesfile->Get("ZW_FactScale1")));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_re2, (TH1F*)templatesfile->Get("ZW_RenScale2")));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_fa2, (TH1F*)templatesfile->Get("ZW_FactScale2")));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_pdf, (TH1F*)templatesfile->Get("ZW_PDF")));
  if (!connectWZ) 
    makeBinList("WJets_SR", met, wspace, wln_SR_hist, wln_SR_bins);
  else   
    makeConnectedBinList("WJets_SR", met, wspace, (TH1F*)templatesfile->Get(("zwjcorewkhist_"+observable).c_str()), wln_SR_syst, znn_SR_bins, &wln_SR_bins);
  
  // Other MC backgrounds
  addTemplate("ZJets_SR"     , vars, wspace, (TH1F*)templatesfile->Get(("zjethist_"+observable).c_str()));
  addTemplate("Dibosons_SR"  , vars, wspace, (TH1F*)templatesfile->Get(("dbkghist_"+observable).c_str()));
  // QCD background
  TH1F* qcdhist = (TH1F*)templatesfile->Get(("qbkghist_"+observable).c_str());
  qcdhist->Scale(2.0);  
  addTemplate("QCD_SR"       , vars, wspace, qcdhist);

  ////////////////////////////////////
  // -------- CR Di-Muon  -------- //
  ///////////////////////////////////
  cout<<"Make CR Di-Muon  templates ..."<<endl;

  addTemplate("data_obs_ZM", vars, wspace, (TH1F*)templatesfile->Get(("datahistzmm_"+observable).c_str()));
  // Z->mumu connected with Z->nunu SR
  vector<pair<RooRealVar*, TH1*> >   znn_ZM_syst;
  makeConnectedBinList("Znunu_ZM", met, wspace, (TH1F*)templatesfile->Get(("zmmcorhist_"+observable).c_str()), znn_ZM_syst, znn_SR_bins);
  
  // Other MC backgrounds in dimuon control region
  addTemplate("WJets_ZM"     , vars, wspace, (TH1F*)templatesfile->Get(("vlbkghistzmm_"+observable).c_str()));
  addTemplate("Top_ZM"       , vars, wspace, (TH1F*)templatesfile->Get(("tbkghistzmm_"+observable).c_str()));
  addTemplate("QCD_ZM"       , vars, wspace, (TH1F*)templatesfile->Get(("qbkghistzmm_"+observable).c_str()));
  addTemplate("Dibosons_ZM"  , vars, wspace, (TH1F*)templatesfile->Get(("dbkghistzmm_"+observable).c_str()));

  ////////////////////////////////////////
  // -------- CR Di-Electron  -------- //
  ///////////////////////////////////////
  cout<<"Make CR Di-Electron  templates ..."<<endl;

  addTemplate("data_obs_ZE"  , vars, wspace, (TH1F*)templatesfile->Get(("datahistzee_"+observable).c_str()));
  // Z->ee connected with Z->nunu SR
  vector<pair<RooRealVar*, TH1*> > znn_ZE_syst;
  makeConnectedBinList("Znunu_ZE", met, wspace, (TH1F*)templatesfile->Get(("zeecorhist_"+observable).c_str()), znn_ZE_syst, znn_SR_bins);
  
  // Other MC backgrounds in dielectron control region
  addTemplate("WJets_ZE"     , vars, wspace, (TH1F*)templatesfile->Get(("vlbkghistzee_"+observable).c_str()));
  addTemplate("Top_ZE"       , vars, wspace, (TH1F*)templatesfile->Get(("tbkghistzee_"+observable).c_str()));
  addTemplate("QCD_ZE"       , vars, wspace, (TH1F*)templatesfile->Get(("qbkghistzee_"+observable).c_str()));
  addTemplate("Dibosons_ZE"  , vars, wspace, (TH1F*)templatesfile->Get(("dbkghistzee_"+observable).c_str()));

  ///////////////////////////////////////
  // -------- CR Gamma+jets  -------- //
  //////////////////////////////////////
  cout<<"Make CR Gamma+jets  templates ..."<<endl;

  addTemplate("data_obs_GJ"  , vars, wspace, (TH1F*)templatesfile->Get(("datahistgam_"+observable).c_str()));
  // Gamma+jets --> connected with Z->nunu
  vector<pair<RooRealVar*, TH1*> > znn_GJ_syst;
  RooRealVar* znn_GJ_re1 = new RooRealVar("Znunu_GJ_RenScale1" , "", 0., -5., 5.);
  RooRealVar* znn_GJ_fa1 = new RooRealVar("Znunu_GJ_FactScale1", "", 0., -5., 5.);
  RooRealVar* znn_GJ_re2 = new RooRealVar("Znunu_GJ_RenScale2" , "", 0., -5., 5.);
  RooRealVar* znn_GJ_fa2 = new RooRealVar("Znunu_GJ_FactScale2", "", 0., -5., 5.);
  RooRealVar* znn_GJ_pdf = new RooRealVar("Znunu_GJ_PDF"       , "", 0., -5., 5.);
  RooRealVar* znn_GJ_fpc = new RooRealVar("Znunu_GJ_Footprint" , "", 0., -5., 5.);
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(NULL      , (TH1F*)templatesfile->Get("ZG_EWK")));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_re1, (TH1F*)templatesfile->Get("ZG_RenScale1")));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_fa1, (TH1F*)templatesfile->Get("ZG_FactScale1")));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_re2, (TH1F*)templatesfile->Get("ZG_RenScale2")));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_fa2, (TH1F*)templatesfile->Get("ZG_FactScale2")));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_pdf, (TH1F*)templatesfile->Get("ZG_PDF")));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_fpc, (TH1F*)templatesfile->Get("ZG_Footprint")));

  makeConnectedBinList("Znunu_GJ", met, wspace, (TH1F*)templatesfile->Get(("gamcorewkhist_"+observable).c_str()), znn_GJ_syst, znn_SR_bins);
  
  // Other MC backgrounds photon+jets control region
  addTemplate("QCD_GJ"     , vars, wspace, (TH1F*)templatesfile->Get(("qbkghistgam_"+observable).c_str()));
  

  ///////////////////////////////////////
  // -------- CR Single-Muon  -------- //
  //////////////////////////////////////
  cout<<"Make CR Single-Mu  templates ..."<<endl;

  addTemplate("data_obs_WM"  , vars, wspace, (TH1F*)templatesfile->Get(("datahistwmn_"+observable).c_str()));
  // connected W->munu with W+jets SR
  vector<pair<RooRealVar*, TH1*> > wln_WM_syst;
  makeConnectedBinList("WJets_WM", met, wspace, (TH1F*)templatesfile->Get(("wmncorhist_"+observable).c_str()), wln_WM_syst, wln_SR_bins);
		       
  // Other MC backgrounds in single muon control region
  addTemplate("ZJets_WM"     , vars, wspace, (TH1F*)templatesfile->Get(("vllbkghistwmn_"+observable).c_str()));
  addTemplate("Top_WM"       , vars, wspace, (TH1F*)templatesfile->Get(("tbkghistwmn_"+observable).c_str()));
  addTemplate("QCD_WM"       , vars, wspace, (TH1F*)templatesfile->Get(("qbkghistwmn_"+observable).c_str()));
  addTemplate("Dibosons_WM"  , vars, wspace, (TH1F*)templatesfile->Get(("dbkghistwmn_"+observable).c_str()));
  
  //////////////////////////////////..../////
  // -------- CR Single-Electron  -------- //
  //////////////////////////////////////////
  cout<<"Make CR Single-El  templates ..."<<endl;

  addTemplate("data_obs_WE"  , vars, wspace, (TH1F*)templatesfile->Get(("datahistwen_"+observable).c_str()));
  // connected W->enu with W+jets SR 
  vector<pair<RooRealVar*, TH1*> > wln_WE_syst;
  makeConnectedBinList("WJets_WE", met, wspace, (TH1F*)templatesfile->Get(("wencorhist_"+observable).c_str()), wln_WE_syst, wln_SR_bins);
  
  // Other MC backgrounds in single electron control region
  addTemplate("ZJets_WE"     , vars, wspace, (TH1F*)templatesfile->Get(("vllbkghistwen_"+observable).c_str()));
  addTemplate("Top_WE"       , vars, wspace, (TH1F*)templatesfile->Get(("tbkghistwen_"+observable).c_str()));
  addTemplate("QCD_WE"       , vars, wspace, (TH1F*)templatesfile->Get(("qbkghistwen_"+observable).c_str()));
  addTemplate("Dibosons_WE"  , vars, wspace, (TH1F*)templatesfile->Get(("dbkghistwen_"+observable).c_str()));

  /////////////////////////////////////
  // -------- CR Top-Muon  -------- //
  ////////////////////////////////////
  cout<<"Make CR Top-mu  templates ..."<<endl;

  addTemplate("data_obs_TM"  , vars, wspace, (TH1F*)templatesfile->Get(("datahisttopmu_"+observable).c_str()));
  // connect tt->mu+b with tt SR
  vector<pair<RooRealVar*, TH1*> > top_TM_syst;
  RooRealVar* top_TM_btag   = new RooRealVar("top_TM_btag" , "", 0., -5., 5.);
  top_TM_syst.push_back(pair<RooRealVar*, TH1*>(top_TM_btag, (TH1F*)templatesfile->Get("TOP_MU_B")));

  makeConnectedBinList("Top_TM", met, wspace, (TH1F*)templatesfile->Get(("topmucorhist_"+observable).c_str()), top_TM_syst, top_SR_bins);
  
  // Other MC backgrounds in single electron control region
  //  addTemplate("ZJets_TM"     , vars, wspace, (TH1F*)templatesfile->Get(("vllbkghisttopmu_"+observable).c_str()));
  //  addTemplate("WJets_TM"     , vars, wspace, (TH1F*)templatesfile->Get(("vlbkghisttopmu_"+observable).c_str()));
  //  addTemplate("QCD_TM"       , vars, wspace, (TH1F*)templatesfile->Get(("qbkghisttopmu_"+observable).c_str()));
  //  addTemplate("Dibosons_TM"  , vars, wspace, (TH1F*)templatesfile->Get(("dbkghisttopmu_"+observable).c_str()));
  /*
  /////////////////////////////////////
  // -------- CR Top-Electron-------- //
  ////////////////////////////////////
  cout<<"Make CR Top-el  templates ..."<<endl;

  addTemplate("data_obs_TE"  , vars, wspace, (TH1F*)templatesfile->Get(("datahisttopel_"+observable).c_str()));
  // connect tt->mu+b with tt SR
  vector<pair<RooRealVar*, TH1*> > top_TE_syst;
  RooRealVar* top_TE_btag   = new RooRealVar("top_TE_btag" , "", 0., -5., 5.);
  top_TE_syst.push_back(pair<RooRealVar*, TH1*>(top_TM_btag, (TH1F*)templatesfile->Get("TOP_EL_B")));

  makeConnectedBinList("Top_TE", met, wspace, (TH1F*)templatesfile->Get(("topelcorhist_"+observable).c_str()), top_TE_syst, top_SR_bins);
  
  // Other MC backgrounds in single electron control region
  addTemplate("ZJets_TE"     , vars, wspace, (TH1F*)templatesfile->Get(("vllbkghisttopel_"+observable).c_str()));
  addTemplate("WJets_TE"     , vars, wspace, (TH1F*)templatesfile->Get(("vlbkghisttopel_"+observable).c_str()));
  addTemplate("QCD_TE"       , vars, wspace, (TH1F*)templatesfile->Get(("qbkghisttopel_"+observable).c_str()));
  addTemplate("Dibosons_TE"  , vars, wspace, (TH1F*)templatesfile->Get(("dbkghisttopel_"+observable).c_str()));
  */    
  // ---------------------------- Write out the workspace -----------------------------------------------------------------//
  outfile->cd();
  wspace.Write();
  outfile->Close();

  return;
}
