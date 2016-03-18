#include <string>
#include <vector>
#include <utility>

using namespace std;

void addDummyBinContent(TH1* histo){

  for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
    histo->SetBinContent(iBin+1,10.e-6);
    histo->SetBinError(iBin+1,10.e-6);
  }

}

void fixShapeUncertainty(TH1* nominalHisto, TH1* sysHisto, float xPoint, float xValue){
  if(sysHisto == 0 || sysHisto == NULL) return;
  for(int iBin = 0; iBin < sysHisto->GetNbinsX(); iBin++){
    if(iBin >= nominalHisto->FindBin(xPoint))
      sysHisto->SetBinContent(iBin+1,nominalHisto->GetBinContent(iBin+1)*xValue);
  }
}


// function to create a RooDataHist from TH1F and import it in a workspace
void addTemplate(string procname, RooArgList& varlist, RooWorkspace& ws, TH1F* hist) {  
  if(hist == 0 || hist == NULL) return;
  if(hist->Integral() == 0)
    addDummyBinContent(hist); // avoind empty histograms in the workspace

  RooDataHist rhist((procname).c_str(), "", varlist, hist);
  ws.import(rhist);
}

void generateStatTemplate(string procname, RooArgList& varlist, RooWorkspace& ws, TH1* histo, float scaleUncertainty){

  vector<TH1F*> histStatUp;
  vector<TH1F*> histStatDw;

  if(histo == 0 || histo == NULL) return;
  if(histo->Integral() == 0){
    addDummyBinContent(hist);
  }

  for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
    histStatUp.push_back((TH1F*) histo->Clone(Form("%s_%s_bin_%i_statUp_tmp",procname.c_str(),procname.c_str(),iBin+1)));
    histStatDw.push_back((TH1F*) histo->Clone(Form("%s_%s_bin_%i_statDown_tmp",procname.c_str(),procname.c_str(),iBin+1)));
  }

  for( size_t iHisto =0; iHisto < histStatUp.size(); iHisto++){
    histStatUp.at(iHisto)->SetBinContent(iHisto+1,histo->GetBinContent(iHisto+1)+histo->GetBinError(iHisto+1)*scaleUncertainty);
    histStatDw.at(iHisto)->SetBinContent(iHisto+1,histo->GetBinContent(iHisto+1)-histo->GetBinError(iHisto+1)*scaleUncertainty);
  }

  for(size_t iHisto =0; iHisto < histStatUp.size(); iHisto++) {
    stringstream name;
    // add two times proc name in order to obtain correct lines in the datacard without name overlapping
    name << procname << "_"<< procname << "_CMS_bin" << iHisto+1 << "_statUp";
    addTemplate(name.str(),varlist,ws,histStatUp[iHisto]); 
  }
  for(size_t iHisto =0; iHisto < histStatDw.size(); iHisto++){
    stringstream name;
    // add two times proc name in order to obtain correct lines in the datacard without name overlapping
    name << procname << "_" << procname << "_CMS_bin" << iHisto+1 << "_statDown";
    addTemplate(name.str(),varlist,ws,histStatDw[iHisto]); 
  }
}


// Make list of bins of a TH1F as RooArgList of RooRealVar and building the RooParametricHist (to be Run in the release with combine)
void makeBinList(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* hist, RooArgList& binlist, bool setConst = false) {

  if(hist == 0 || hist == NULL) return;
  if(hist->Integral() == 0){
    addDummyBinContent(hist);
  }

  // loop over histo bins
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    stringstream binss;
    binss << procname << "_bin" << i;
    RooRealVar* binvar;

    // make a RooRealVar for each bin [0,2*binContent]
    if (!setConst)
      binvar = new RooRealVar(binss.str().c_str(), "", hist->GetBinContent(i), 0., hist->GetBinContent(i)*10.0);
    else         
      binvar = new RooRealVar(binss.str().c_str(), "", hist->GetBinContent(i));
    
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
			  vector<pair<RooRealVar*, TH1*> > syst, const RooArgList& srbinlist, 
			  RooArgList* crbinlist=NULL,
			  string observable = "met") {
  

  if(rhist == 0 || rhist == NULL) return;
  if(rhist->Integral() == 0){
    addDummyBinContent(hist);
  }

  // bin list for the CR
  if (crbinlist == NULL) 
    crbinlist = new RooArgList();

  // Loop on ratio hist
  float extreme_tmp = 5;
  for (int i = 1; i <= rhist->GetNbinsX(); i++) {
    stringstream rbinss;
    rbinss << "r_" << procname << "_bin" << i;
    // Fixed value for each bin of the ratio
    RooRealVar* rbinvar = new RooRealVar(rbinss.str().c_str(), "", rhist->GetBinContent(i));

    // uncertainty histograms for systematics
    stringstream rerrbinss;
    rerrbinss << procname << "_bin" << i << "_Runc";
    // Nuisance for the Final fit for each bin (bin-by-bin unc) --> avoid negative values
    float extreme = std::min(5,0.9*rhist->GetBinContent(i)/rhist->GetBinError(i));
    RooRealVar* rerrbinvar = new RooRealVar(rerrbinss.str().c_str(), "", 0., -extreme, extreme);
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
    stringstream formss;
    formss << "@0/";
    formss << "(";
    formss << "@1";
    formss << "*(abs(1+" << rhist->GetBinError(i)/rhist->GetBinContent(i) << "*@2))";
    // systemaitc uncertainty
    for (int j = 0; j < syst.size(); j++) {
      stringstream systbinss;
      if (syst[j].first == NULL) { // add bin by bin
	systbinss << procname << "_bin" << i << "_" << syst[j].second->GetName();
	TString nameSys (systbinss.str());
	nameSys.ReplaceAll(("_"+observable).c_str(),"");
	float extreme = std::min(5,0.9*rhist->GetBinContent(i)/syst[j].second->GetBinContent(i)); 
	RooRealVar* systbinvar = new RooRealVar(nameSys.Data(), "", 0., -extreme, extreme);
	// Add all the systeamtics as new Multiplicative Nuisance for each bin
	fobinlist.add(*systbinvar);
      }
      else{ 
	float extreme = std::min(5,0.9*rhist->GetBinContent(i)/syst[j].second->GetBinContent(i));
	if( extreme < extreme_tmp){
	  syst[j].first->setMin(-extreme);
	  syst[j].first->setMax(extreme);
	  extreme_tmp = extreme;
	}	
	fobinlist.add(*syst[j].first);
      }
      formss << "*(abs(1+" << syst[j].second->GetBinContent(i) << "*@" << j+3 << "))";
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
void createWorkspace(string inputName, 
		     int category, 
		     string outputName    = "workspace.root",
		     string observable    = "met", 
		     bool   isHiggsInvisible = false,
		     float  scaleQCD      = 2, 
		     bool   connectWZ     = true, 
		     bool   connectTop    = false,
		     bool   addShapeSystematics = false,
		     string interaction   = "Vector",
		     string mediatorMass  = "1000", 
OB		     string DMMass        = "50", 
){
  
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");  
  RooMsgService::instance().setSilentMode(kTRUE); 
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  // for templates and sys naming
  string suffix;
  if(category <=1)
    suffix = "MJ";
  else
    suffix = "MV";
    
  // create the output workspace
  cout<<"Create output file ..."<<endl;
  TFile *outfile = new TFile(outputName.c_str(),"RECREATE");
  outfile->cd();
  cout<<"Load binning and observable ..."<<endl;

  // to be fix one combine will switch to ROOT 6
  double xMin = 0.;
  double xMax = 0.;

  if(TString(observable).Contains("met") && category <=1){
    xMin = 200.;
    xMax = 1250.;
  }
  else if(TString(observable).Contains("met") && category >1){
    xMin = 250.;
    xMax = 1000.;
  }
  else
    cout<<"Binning not implemented for the observable "<<observable<<" --> please define it "<<endl;
  
  RooRealVar met((observable+"_"+suffix).c_str(),"",xMin,xMax);
  RooArgList vars(met);
  
  // Templates
  cout<<"Open inputFile ..."<<endl;
  TFile* templatesfile = new TFile(inputName.c_str());

  ///////////////////////////////////////
  // -------- SIGNAL REGION  -------- //
  ///////////////////////////////////////
  cout<<"Make SR templates ..."<<endl;

  // create a workspace for the signal region
  RooWorkspace wspace_SR(("SR_"+suffix).c_str(),(suffix+"_SR").c_str());

  // Add Data
  addTemplate("data_obs_SR_"+suffix, vars, wspace_SR, (TH1F*)templatesfile->Get(("datahist_"+observable).c_str()));

  // Signal shape
  if(!isHiggsInvisible){

    addTemplate("MonoJ_SR_"+suffix, vars, wspace_SR, (TH1F*)templatesfile->Get(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
    addTemplate("MonoW_SR_"+suffix, vars, wspace_SR, (TH1F*)templatesfile->Get(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
    addTemplate("MonoZ_SR_"+suffix, vars, wspace_SR, (TH1F*)templatesfile->Get(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));

    if(addShapeSystematics){

      // monoJ
      TH1F* nominalHisto = (TH1F*)templatesfile->Get(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histobUp  = (TH1F*)templatesfile->Get(("monoJhist_bUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histobDw  = (TH1F*)templatesfile->Get(("monoJhist_bDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histoJesUp  = (TH1F*)templatesfile->Get(("monoJhist_metJetUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histoJesDw  = (TH1F*)templatesfile->Get(("monoJhist_metJetDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histoJerUp  = (TH1F*)templatesfile->Get(("monoJhist_metResUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histoJerDw  = (TH1F*)templatesfile->Get(("monoJhist_metResDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histoUncUp  = (TH1F*)templatesfile->Get(("monoJhist_metUncUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histoUncDw  = (TH1F*)templatesfile->Get(("monoJhist_metUncDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      
      if(nominalHisto){
	
      if(nominalHisto->GetNbinsX() > 10){
	fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
      }
      
      addTemplate("MonoJ_SR_"+suffix+"_CMS_btagUp",   vars, wspace_SR, histobUp);
      addTemplate("MonoJ_SR_"+suffix+"_CMS_btagDown", vars, wspace_SR, histobDw);
      addTemplate("MonoJ_SR_"+suffix+"_CMS_scale_jUp",    vars, wspace_SR, histoJesUp);
      addTemplate("MonoJ_SR_"+suffix+"_CMS_scale_jDown",  vars, wspace_SR, histoJesDw);
      addTemplate("MonoJ_SR_"+suffix+"_CMS_res_jUp",    vars, wspace_SR, histoJerUp);
      addTemplate("MonoJ_SR_"+suffix+"_CMS_res_jDown",  vars, wspace_SR, histoJerDw);
      addTemplate("MonoJ_SR_"+suffix+"_CMS_uncUp",    vars, wspace_SR, histoUncUp);
      addTemplate("MonoJ_SR_"+suffix+"_CMS_uncDown",  vars, wspace_SR, histoUncDw);
    }
      
      nominalHisto = (TH1F*)templatesfile->Get(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histobUp  = (TH1F*)templatesfile->Get(("monoWhist_bUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histobDw  = (TH1F*)templatesfile->Get(("monoWhist_bDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJesUp  = (TH1F*)templatesfile->Get(("monoWhist_metJetUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJesDw  = (TH1F*)templatesfile->Get(("monoWhist_metJetDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJerUp  = (TH1F*)templatesfile->Get(("monoWhist_metResUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJerDw  = (TH1F*)templatesfile->Get(("monoWhist_metResDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoUncUp  = (TH1F*)templatesfile->Get(("monoWhist_metUncUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoUncDw  = (TH1F*)templatesfile->Get(("monoWhist_metUncDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      
      if(nominalHisto){
	
	if(nominalHisto->GetNbinsX() > 10){
	  fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	  fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	  fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	  fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	  fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	  fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	  fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	  fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
	}
	
	addTemplate("MonoW_SR_"+suffix+"_CMS_btagUp",   vars, wspace_SR, histobUp);
	addTemplate("MonoW_SR_"+suffix+"_CMS_btagDown", vars, wspace_SR, histobDw);
	addTemplate("MonoW_SR_"+suffix+"_CMS_scale_jUp",    vars, wspace_SR, histoJesUp);
	addTemplate("MonoW_SR_"+suffix+"_CMS_scale_jDown",  vars, wspace_SR, histoJesDw);
	addTemplate("MonoW_SR_"+suffix+"_CMS_res_jUp",    vars, wspace_SR, histoJerUp);
	addTemplate("MonoW_SR_"+suffix+"_CMS_res_jDown",  vars, wspace_SR, histoJerDw);
	addTemplate("MonoW_SR_"+suffix+"_CMS_uncUp",    vars, wspace_SR, histoUncUp);
	addTemplate("MonoW_SR_"+suffix+"_CMS_uncDown",  vars, wspace_SR, histoUncDw);
      }

      // monoW
      nominalHisto = (TH1F*)templatesfile->Get(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histobUp  = (TH1F*)templatesfile->Get(("monoZhist_bUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histobDw  = (TH1F*)templatesfile->Get(("monoZhist_bDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJesUp  = (TH1F*)templatesfile->Get(("monoZhist_metJetUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJesDw  = (TH1F*)templatesfile->Get(("monoZhist_metJetDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJerUp  = (TH1F*)templatesfile->Get(("monoZhist_metResUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJerDw  = (TH1F*)templatesfile->Get(("monoZhist_metResDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoUncUp  = (TH1F*)templatesfile->Get(("monoZhist_metUncUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoUncDw  = (TH1F*)templatesfile->Get(("monoZhist_metUncDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      
      if(nominalHisto){
	
	if(nominalHisto->GetNbinsX() > 10){
	  fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	  fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	  fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	  fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	  fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	  fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	  fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	  fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
	}
	
	addTemplate("MonoZ_SR_"+suffix+"_CMS_btagUp",   vars, wspace_SR, histobUp);
	addTemplate("MonoZ_SR_"+suffix+"_CMS_btagDown", vars, wspace_SR, histobDw);
	addTemplate("MonoZ_SR_"+suffix+"_CMS_scale_jUp",    vars, wspace_SR, histoJesUp);
	addTemplate("MonoZ_SR_"+suffix+"_CMS_scale_jDown",  vars, wspace_SR, histoJesDw);
	addTemplate("MonoZ_SR_"+suffix+"_CMS_res_jUp",    vars, wspace_SR, histoJerUp);
	addTemplate("MonoZ_SR_"+suffix+"_CMS_res_jDown",  vars, wspace_SR, histoJerDw);
	addTemplate("MonoZ_SR_"+suffix+"_CMS_uncUp",    vars, wspace_SR, histoUncUp);
	addTemplate("MonoZ_SR_"+suffix+"_CMS_uncDown",  vars, wspace_SR, histoUncDw);
      }
      // statistics
      generateStatTemplate("MonoJ_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->Get(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
      generateStatTemplate("MonoW_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->Get(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
      generateStatTemplate("MonoZ_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->Get(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
    }
  }
  else{

    addTemplate("ggH_SR_"+suffix, vars, wspace_SR, (TH1F*)templatesfile->Get(("ggHhist_"+mediatorMass+"_"+observable).c_str()));
    addTemplate("vbfH_SR_"+suffix,vars, wspace_SR, (TH1F*)templatesfile->Get(("vbfHhist_"+mediatorMass+"_"+observable).c_str()));
    addTemplate("wH_SR_"+suffix,  vars, wspace_SR, (TH1F*)templatesfile->Get(("wHhist_"+mediatorMass+"_"+observable).c_str()));
    addTemplate("zH_SR_"+suffix,  vars, wspace_SR, (TH1F*)templatesfile->Get(("zHhist_"+mediatorMass+"_"+observable).c_str()));

    if(addShapeSystematics){

      // ggH
      TH1F* nominalHisto = (TH1F*)templatesfile->Get(("ggHhist_"+mediatorMass+"_"+observable).c_str());
      TH1F* histobUp  = (TH1F*)templatesfile->Get(("ggHhist_bUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histobDw  = (TH1F*)templatesfile->Get(("ggHhist_bDw_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoJesUp  = (TH1F*)templatesfile->Get(("ggHhist_metJetUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoJesDw  = (TH1F*)templatesfile->Get(("ggHhist_metJetDw_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoJerUp  = (TH1F*)templatesfile->Get(("ggHhist_metResUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoJerDw  = (TH1F*)templatesfile->Get(("ggHhist_metResDw_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoUncUp  = (TH1F*)templatesfile->Get(("ggHhist_metUncUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoUncDw  = (TH1F*)templatesfile->Get(("ggHhist_metUncDw_"+mediatorMass+"_"+observable).c_str());

      TH1F* histoRenUp = (TH1F*)templatesfile->Get(("ggHhist_renUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoRenDw = (TH1F*)templatesfile->Get(("ggHhist_renDw_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoFacUp = (TH1F*)templatesfile->Get(("ggHhist_facUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoFacDw = (TH1F*)templatesfile->Get(("ggHhist_facDw_"+mediatorMass+"_"+observable).c_str());

      if(nominalHisto){
	
	if(nominalHisto->GetNbinsX() > 10){
	  fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	  fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	  fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	  fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	  fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	  fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	  fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	  fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
	}
	
	addTemplate("ggH_SR_"+suffix+"_QCDScale_ggH_ren_acceptUp", vars, wspace_SR, histoRenUp);
	addTemplate("ggH_SR_"+suffix+"_QCDScale_ggH_ren_acceptDown", vars, wspace_SR, histoRenDw);
	addTemplate("ggH_SR_"+suffix+"_QCDScale_ggH_fac_acceptUp", vars, wspace_SR, histoFacUp);
	addTemplate("ggH_SR_"+suffix+"_QCDScale_ggH_fac_acceptDown", vars, wspace_SR, histoFacDw);
	
	addTemplate("ggH_SR_"+suffix+"_CMS_btagUp",   vars, wspace_SR, histobUp);
	addTemplate("ggH_SR_"+suffix+"_CMS_btagDown", vars, wspace_SR, histobDw);
	addTemplate("ggH_SR_"+suffix+"_CMS_scale_jUp",    vars, wspace_SR, histoJesUp);
	addTemplate("ggH_SR_"+suffix+"_CMS_scale_jDown",  vars, wspace_SR, histoJesDw);
	addTemplate("ggH_SR_"+suffix+"_CMS_res_jUp",    vars, wspace_SR, histoJerUp);
	addTemplate("ggH_SR_"+suffix+"_CMS_res_jDown",  vars, wspace_SR, histoJerDw);
	addTemplate("ggH_SR_"+suffix+"_CMS_uncUp",    vars, wspace_SR, histoUncUp);
	addTemplate("ggH_SR_"+suffix+"_CMS_uncDown",  vars, wspace_SR, histoUncDw);
      }
      
      nominalHisto = (TH1F*)templatesfile->Get(("vbfHhist_"+mediatorMass+"_"+observable).c_str());
      histobUp  = (TH1F*)templatesfile->Get(("vbfHhist_bUp_"+mediatorMass+"_"+observable).c_str());
      histobDw  = (TH1F*)templatesfile->Get(("vbfHhist_bDw_"+mediatorMass+"_"+observable).c_str());
      histoJesUp  = (TH1F*)templatesfile->Get(("vbfHhist_metJetUp_"+mediatorMass+"_"+observable).c_str());
      histoJesDw  = (TH1F*)templatesfile->Get(("vbfHhist_metJetDw_"+mediatorMass+"_"+observable).c_str());
      histoJerUp  = (TH1F*)templatesfile->Get(("vbfHhist_metResUp_"+mediatorMass+"_"+observable).c_str());
      histoJerDw  = (TH1F*)templatesfile->Get(("vbfHhist_metResDw_"+mediatorMass+"_"+observable).c_str());
      histoUncUp  = (TH1F*)templatesfile->Get(("vbfHhist_metUncUp_"+mediatorMass+"_"+observable).c_str());
      histoUncDw  = (TH1F*)templatesfile->Get(("vbfHhist_metUncDw_"+mediatorMass+"_"+observable).c_str());
      
      if(nominalHisto){
	
	if(nominalHisto->GetNbinsX() > 10){
	  fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	  fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	  fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	  fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	  fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	  fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	  fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	  fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
	}
	
	addTemplate("vbfH_SR_"+suffix+"_CMS_btagUp",   vars, wspace_SR, histobUp);
	addTemplate("vbfH_SR_"+suffix+"_CMS_btagDown", vars, wspace_SR, histobDw);
	addTemplate("vbfH_SR_"+suffix+"_CMS_scale_jUp",    vars, wspace_SR, histoJesUp);
	addTemplate("vbfH_SR_"+suffix+"_CMS_scale_jDown",  vars, wspace_SR, histoJesDw);
	addTemplate("vbfH_SR_"+suffix+"_CMS_res_jUp",    vars, wspace_SR, histoJerUp);
	addTemplate("vbfH_SR_"+suffix+"_CMS_res_jDown",  vars, wspace_SR, histoJerDw);
	addTemplate("vbfH_SR_"+suffix+"_CMS_uncUp",    vars, wspace_SR, histoUncUp);
	addTemplate("vbfH_SR_"+suffix+"_CMS_uncDown",  vars, wspace_SR, histoUncDw);
      }

      // vbfH
      nominalHisto = (TH1F*)templatesfile->Get(("wHhist_"+mediatorMass+"_"+observable).c_str());
      histobUp  = (TH1F*)templatesfile->Get(("wHhist_bUp_"+mediatorMass+"_"+observable).c_str());
      histobDw  = (TH1F*)templatesfile->Get(("wHhist_bDw_"+mediatorMass+"_"+observable).c_str());
      histoJesUp  = (TH1F*)templatesfile->Get(("wHhist_metJetUp_"+mediatorMass+"_"+observable).c_str());
      histoJesDw  = (TH1F*)templatesfile->Get(("wHhist_metJetDw_"+mediatorMass+"_"+observable).c_str());
      histoJerUp  = (TH1F*)templatesfile->Get(("wHhist_metResUp_"+mediatorMass+"_"+observable).c_str());
      histoJerDw  = (TH1F*)templatesfile->Get(("wHhist_metResDw_"+mediatorMass+"_"+observable).c_str());
      histoUncUp  = (TH1F*)templatesfile->Get(("wHhist_metUncUp_"+mediatorMass+"_"+observable).c_str());
      histoUncDw  = (TH1F*)templatesfile->Get(("wHhist_metUncDw_"+mediatorMass+"_"+observable).c_str());
      
      if(nominalHisto){
	
	if(nominalHisto->GetNbinsX() > 10){
	  fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	  fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	  fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	  fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	  fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	  fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	  fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	  fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
	}
	
	addTemplate("wH_SR_"+suffix+"_CMS_btagUp",   vars, wspace_SR, histobUp);
	addTemplate("wH_SR_"+suffix+"_CMS_btagDown", vars, wspace_SR, histobDw);
	addTemplate("wH_SR_"+suffix+"_CMS_scale_jUp",    vars, wspace_SR, histoJesUp);
	addTemplate("wH_SR_"+suffix+"_CMS_scale_jDown",  vars, wspace_SR, histoJesDw);
	addTemplate("wH_SR_"+suffix+"_CMS_res_jUp",    vars, wspace_SR, histoJerUp);
	addTemplate("wH_SR_"+suffix+"_CMS_res_jDown",  vars, wspace_SR, histoJerDw);
	addTemplate("wH_SR_"+suffix+"_CMS_uncUp",    vars, wspace_SR, histoUncUp);
	addTemplate("wH_SR_"+suffix+"_CMS_uncDown",  vars, wspace_SR, histoUncDw);

      }

      // vbfH
      nominalHisto = (TH1F*)templatesfile->Get(("zHhist_"+mediatorMass+"_"+observable).c_str());
      histobUp  = (TH1F*)templatesfile->Get(("zHhist_bUp_"+mediatorMass+"_"+observable).c_str());
      histobDw  = (TH1F*)templatesfile->Get(("zHhist_bDw_"+mediatorMass+"_"+observable).c_str());
      histoJesUp  = (TH1F*)templatesfile->Get(("zHhist_metJetUp_"+mediatorMass+"_"+observable).c_str());
      histoJesDw  = (TH1F*)templatesfile->Get(("zHhist_metJetDw_"+mediatorMass+"_"+observable).c_str());
      histoJerUp  = (TH1F*)templatesfile->Get(("zHhist_metResUp_"+mediatorMass+"_"+observable).c_str());
      histoJerDw  = (TH1F*)templatesfile->Get(("zHhist_metResDw_"+mediatorMass+"_"+observable).c_str());
      histoUncUp  = (TH1F*)templatesfile->Get(("zHhist_metUncUp_"+mediatorMass+"_"+observable).c_str());
      histoUncDw  = (TH1F*)templatesfile->Get(("zHhist_metUncDw_"+mediatorMass+"_"+observable).c_str());
      
      if(nominalHisto){
	
	if(nominalHisto->GetNbinsX() > 10){
	  fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	  fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	  fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	  fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	  fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	  fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	  fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	  fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
	}
	
	addTemplate("zH_SR_"+suffix+"_CMS_btagUp",   vars, wspace_SR, histobUp);
	addTemplate("zH_SR_"+suffix+"_CMS_btagDown", vars, wspace_SR, histobDw);
	addTemplate("zH_SR_"+suffix+"_CMS_scale_jUp",    vars, wspace_SR, histoJesUp);
	addTemplate("zH_SR_"+suffix+"_CMS_scale_jDown",  vars, wspace_SR, histoJesDw);
	addTemplate("zH_SR_"+suffix+"_CMS_res_jUp",    vars, wspace_SR, histoJerUp);
	addTemplate("zH_SR_"+suffix+"_CMS_res_jDown",  vars, wspace_SR, histoJerDw);
	addTemplate("zH_SR_"+suffix+"_CMS_uncUp",    vars, wspace_SR, histoUncUp);
	addTemplate("zH_SR_"+suffix+"_CMS_uncDown",  vars, wspace_SR, histoUncDw);
      }
    }
    // statistics
    generateStatTemplate("ggH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->Get(("ggHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
    generateStatTemplate("vbfH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->Get(("vbfHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
    generateStatTemplate("wH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->Get(("wHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
    generateStatTemplate("zH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->Get(("zHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
  }
  
  // Zvv background --> to be extracted from CRs
  TH1F* znn_SR_hist = (TH1F*) templatesfile->Get(("zinvhist_"+observable).c_str());
  RooArgList znn_SR_bins; 
  // create a RooParametric hist with one RooRealVar per bin 
  makeBinList("Znunu_SR_"+suffix, met, wspace_SR, znn_SR_hist, znn_SR_bins);

  // Top background --> to be extracted from CRs
  RooArgList top_SR_bins;
  TH1F* top_SR_hist = NULL;

  // for data driven top estimation
  if(connectTop){
    top_SR_hist = (TH1F*) templatesfile->Get(("tbkghist_"+observable).c_str());
    RooArgList top_SR_bins; 
    makeBinList("Top_SR_"+suffix, met, wspace_SR, top_SR_hist, top_SR_bins);
  }
  else{ // rely on MC + systematics

    addTemplate("Top_SR_"+suffix,vars, wspace_SR, (TH1F*)templatesfile->Get(("tbkghist_"+observable).c_str()));

    if(addShapeSystematics){

      TH1F* nominalHisto = (TH1F*)templatesfile->Get(("tbkghist_"+observable).c_str());
      TH1F* histobUp = (TH1F*)templatesfile->Get(("tbkghist_bUp_"+observable).c_str());
      TH1F* histobDw = (TH1F*)templatesfile->Get(("tbkghist_bDw_"+observable).c_str());
      TH1F* histoJesUp = (TH1F*)templatesfile->Get(("tbkghist_metJetUp_"+observable).c_str());
      TH1F* histoJesDw = (TH1F*)templatesfile->Get(("tbkghist_metJetDw_"+observable).c_str());
      TH1F* histoJerUp = (TH1F*)templatesfile->Get(("tbkghist_metResUp_"+observable).c_str());
      TH1F* histoJerDw = (TH1F*)templatesfile->Get(("tbkghist_metResDw_"+observable).c_str());
      TH1F* histoUncUp = (TH1F*)templatesfile->Get(("tbkghist_metUncUp_"+observable).c_str());
      TH1F* histoUncDw = (TH1F*)templatesfile->Get(("tbkghist_metUncDw_"+observable).c_str());

      if(nominalHisto->GetNbinsX() > 10){
	fixShapeUncertainty(nominalHisto,histobUp,500.,1.06);
	fixShapeUncertainty(nominalHisto,histobDw,500.,0.94);
	fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.10);
	fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.90);
	fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.03);
	fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.97);
	fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
      }

      addTemplate("Top_SR_"+suffix+"_CMS_btagUp",vars, wspace_SR, histobUp);
      addTemplate("Top_SR_"+suffix+"_CMS_btagDown",vars, wspace_SR, histobDw);
      addTemplate("Top_SR_"+suffix+"_CMS_scale_jUp",vars, wspace_SR, histoJesUp);
      addTemplate("Top_SR_"+suffix+"_CMS_scale_jDown",vars, wspace_SR, histoJesDw);
      addTemplate("Top_SR_"+suffix+"_CMS_res_jUp",vars, wspace_SR, histoJerUp);
      addTemplate("Top_SR_"+suffix+"_CMS_res_jDown",vars, wspace_SR, histoJerDw);
      addTemplate("Top_SR_"+suffix+"_CMS_uncUp",vars, wspace_SR, histoUncUp);
      addTemplate("Top_SR_"+suffix+"_CMS_uncDown",vars, wspace_SR, histoUncDw);
      
      generateStatTemplate("Top_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->Get(("tbkghist_"+observable).c_str()));
    }
  }

  // WJets background --> to be extracted from CRs, with connection to Z->nunu
  TH1F* wln_SR_hist = (TH1F*) templatesfile->Get(("wjethist_"+observable).c_str());
  RooArgList wln_SR_bins;
  // set of correlated systematic uncertainties for the Z/W ratio
  vector<pair<RooRealVar*, TH1*> > wln_SR_syst;

  RooRealVar* wln_SR_re1 = new RooRealVar("WJets_SR_RenScale1", ""  , 0., -5., 5.);
  RooRealVar* wln_SR_fa1 = new RooRealVar("WJets_SR_FactScale1", "" , 0., -5., 5.);
  RooRealVar* wln_SR_re2 = new RooRealVar("WJets_SR_RenScale2", ""  , 0., -5., 5.);
  RooRealVar* wln_SR_fa2 = new RooRealVar("WJets_SR_FactScale2", "" , 0., -5., 5.);
  RooRealVar* wln_SR_pdf = new RooRealVar("WJets_SR_PDF", ""        , 0., -5., 5.);

  // NULL means bin-by-bin
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(NULL      , (TH1F*)templatesfile->Get(("ZW_EWK_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_re1, (TH1F*)templatesfile->Get(("ZW_RenScale1_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_fa1, (TH1F*)templatesfile->Get(("ZW_FactScale1_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_re2, (TH1F*)templatesfile->Get(("ZW_RenScale2_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_fa2, (TH1F*)templatesfile->Get(("ZW_FactScale2_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_pdf, (TH1F*)templatesfile->Get(("ZW_PDF_"+observable).c_str())));

  if (!connectWZ) 
    makeBinList("WJets_SR_"+suffix, met, wspace_SR, wln_SR_hist, wln_SR_bins);
  else   
    makeConnectedBinList("WJets_SR_"+suffix, met, wspace_SR, (TH1F*)templatesfile->Get(("zwjcorewkhist_"+observable).c_str()), wln_SR_syst, znn_SR_bins, &wln_SR_bins, observable);

  // Other MC backgrounds
  addTemplate("ZJets_SR_"+suffix     , vars, wspace_SR, (TH1F*)templatesfile->Get(("zjethist_"+observable).c_str()));
  addTemplate("Dibosons_SR_"+suffix  , vars, wspace_SR, (TH1F*)templatesfile->Get(("dbkghist_"+observable).c_str()));
  addTemplate("GJets_SR_"+suffix     , vars, wspace_SR, (TH1F*)templatesfile->Get(("gbkghist_"+observable).c_str()));

  if(addShapeSystematics){

    TH1F* nominalHisto = (TH1F*)templatesfile->Get(("zjethist_"+observable).c_str());
    TH1F* histobUp = (TH1F*)templatesfile->Get(("zjethist_bUp_"+observable).c_str());
    TH1F* histobDw = (TH1F*)templatesfile->Get(("zjethist_bDw_"+observable).c_str());
    TH1F* histoJesUp = (TH1F*)templatesfile->Get(("zjethist_metJetUp_"+observable).c_str());
    TH1F* histoJesDw = (TH1F*)templatesfile->Get(("zjethist_metJetDw_"+observable).c_str());
    TH1F* histoJerUp = (TH1F*)templatesfile->Get(("zjethist_metResUp_"+observable).c_str());
    TH1F* histoJerDw = (TH1F*)templatesfile->Get(("zjethist_metResDw_"+observable).c_str());
    TH1F* histoUncUp = (TH1F*)templatesfile->Get(("zjethist_metUncUp_"+observable).c_str());
    TH1F* histoUncDw = (TH1F*)templatesfile->Get(("zjethist_metUncDw_"+observable).c_str());
    
    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }
    
    addTemplate("ZJets_SR_"+suffix+"_CMS_btagUp", vars, wspace_SR, histobUp);
    addTemplate("ZJets_SR_"+suffix+"_CMS_btagDown", vars, wspace_SR, histobDw);
    addTemplate("ZJets_SR_"+suffix+"_CMS_scale_jUp", vars, wspace_SR, histoJesUp);
    addTemplate("ZJets_SR_"+suffix+"_CMS_scale_jDown", vars, wspace_SR, histoJesDw);
    addTemplate("ZJets_SR_"+suffix+"_CMS_res_jUp", vars, wspace_SR, histoJerUp);
    addTemplate("ZJets_SR_"+suffix+"_CMS_res_jDown", vars, wspace_SR, histoJerDw);
    addTemplate("ZJets_SR_"+suffix+"_CMS_uncUp", vars, wspace_SR, histoUncUp);
    addTemplate("ZJets_SR_"+suffix+"_CMS_uncDown", vars, wspace_SR, histoUncDw);

    generateStatTemplate("ZJets_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->Get(("zjethist_"+observable).c_str()));

    nominalHisto = (TH1F*)templatesfile->Get(("dbkghist_"+observable).c_str());
    histobUp = (TH1F*)templatesfile->Get(("dbkghist_bUp_"+observable).c_str());
    histobDw = (TH1F*)templatesfile->Get(("dbkghist_bDw_"+observable).c_str());
    histoJesUp = (TH1F*)templatesfile->Get(("dbkghist_metJetUp_"+observable).c_str());
    histoJesDw = (TH1F*)templatesfile->Get(("dbkghist_metJetDw_"+observable).c_str());
    histoJerUp = (TH1F*)templatesfile->Get(("dbkghist_metResUp_"+observable).c_str());
    histoJerDw = (TH1F*)templatesfile->Get(("dbkghist_metResDw_"+observable).c_str());
    histoUncUp = (TH1F*)templatesfile->Get(("dbkghist_metUncUp_"+observable).c_str());
    histoUncDw = (TH1F*)templatesfile->Get(("dbkghist_metUncDw_"+observable).c_str());
    
    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }
    
    addTemplate("Dibosons_SR_"+suffix+"_CMS_btagUp", vars, wspace_SR, histobUp);
    addTemplate("Dibosons_SR_"+suffix+"_CMS_btagDown", vars, wspace_SR, histobDw);
    addTemplate("Dibosons_SR_"+suffix+"_CMS_scale_jUp", vars, wspace_SR, histoJesUp);
    addTemplate("Dibosons_SR_"+suffix+"_CMS_scale_jDown", vars, wspace_SR, histoJesDw);
    addTemplate("Dibosons_SR_"+suffix+"_CMS_res_jUp", vars, wspace_SR, histoJerUp);
    addTemplate("Dibosons_SR_"+suffix+"_CMS_res_jDown", vars, wspace_SR, histoJerDw);
    addTemplate("Dibosons_SR_"+suffix+"_CMS_uncUp", vars, wspace_SR, histoUncUp);
    addTemplate("Dibosons_SR_"+suffix+"_CMS_uncDown", vars, wspace_SR, histoUncDw);

    generateStatTemplate("Dibosons_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->Get(("dbkghist_"+observable).c_str()));


    nominalHisto = (TH1F*)templatesfile->Get(("gbkghist_"+observable).c_str());
    histobUp = (TH1F*)templatesfile->Get(("gbkghist_bUp_"+observable).c_str());
    histobDw = (TH1F*)templatesfile->Get(("gbkghist_bDw_"+observable).c_str());
    histoJesUp = (TH1F*)templatesfile->Get(("gbkghist_metJetUp_"+observable).c_str());
    histoJesDw = (TH1F*)templatesfile->Get(("gbkghist_metJetDw_"+observable).c_str());
    histoJerUp = (TH1F*)templatesfile->Get(("gbkghist_metResUp_"+observable).c_str());
    histoJerDw = (TH1F*)templatesfile->Get(("gbkghist_metResDw_"+observable).c_str());
    histoUncUp = (TH1F*)templatesfile->Get(("gbkghist_metUncUp_"+observable).c_str());
    histoUncDw = (TH1F*)templatesfile->Get(("gbkghist_metUncDw_"+observable).c_str());
    
    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }
    
    addTemplate("GJets_SR_"+suffix+"_CMS_btagUp", vars, wspace_SR, histobUp);
    addTemplate("GJets_SR_"+suffix+"_CMS_btagDown", vars, wspace_SR, histobDw);
    addTemplate("GJets_SR_"+suffix+"_CMS_scale_jUp", vars, wspace_SR, histoJesUp);
    addTemplate("GJets_SR_"+suffix+"_CMS_scale_jDown", vars, wspace_SR, histoJesDw);
    addTemplate("GJets_SR_"+suffix+"_CMS_res_jUp", vars, wspace_SR, histoJerUp);
    addTemplate("GJets_SR_"+suffix+"_CMS_res_jDown", vars, wspace_SR, histoJerDw);
    addTemplate("GJets_SR_"+suffix+"_CMS_uncUp", vars, wspace_SR, histoUncUp);
    addTemplate("GJets_SR_"+suffix+"_CMS_uncDown", vars, wspace_SR, histoUncDw);

    generateStatTemplate("GJets_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->Get(("gbkghist_"+observable).c_str()));

  }

  // look for DD qcd background, otherwise MC scaled by a factor 2
  TH1F* qcdhist = (TH1F*)templatesfile->Get(("qbkghistDD_"+observable).c_str());
  if(qcdhist){
    addTemplate("QCD_SR_"+suffix, vars, wspace_SR, qcdhist);
    addTemplate("QCD_SR_"+suffix+"_CMS_qcdUp", vars, wspace_SR, (TH1F*)templatesfile->Get(("qbkghistDD_shapeUp_"+observable).c_str()));
    addTemplate("QCD_SR_"+suffix+"_CMS_qcdDown", vars, wspace_SR, (TH1F*)templatesfile->Get(("qbkghistDD_shapeDw_"+observable).c_str()));
  }
  else{
    qcdhist = (TH1F*)templatesfile->Get(("qbkghist_"+observable).c_str());
    qcdhist->Scale(scaleQCD);  
    addTemplate("QCD_SR_"+suffix, vars, wspace_SR, qcdhist);
  }


  ////////////////////////////////////
  // -------- CR Di-Muon  -------- //
  ///////////////////////////////////
  cout<<"Make CR Di-Muon  templates ..."<<endl;
 
  RooWorkspace wspace_ZM(("ZM_"+suffix).c_str(),("ZM_"+suffix).c_str());
  
  addTemplate("data_obs_ZM_"+suffix, vars, wspace_ZM, (TH1F*)templatesfile->Get(("datahistzmm_"+observable).c_str()));
  // Z->mumu connected with Z->nunu SR
  vector<pair<RooRealVar*, TH1*> >   znn_syst;
  makeConnectedBinList("Znunu_ZM_"+suffix, met, wspace_ZM, (TH1F*)templatesfile->Get(("zmmcorhist_"+observable).c_str()), znn_syst, znn_SR_bins,NULL, observable);
  
  // Other MC backgrounds in dimuon control region
  addTemplate("WJets_ZM_"+suffix, vars, wspace_ZM, (TH1F*)templatesfile->Get(("vlbkghistzmm_"+observable).c_str()));
  addTemplate("Top_ZM_"+suffix , vars, wspace_ZM, (TH1F*)templatesfile->Get(("tbkghistzmm_"+observable).c_str()));
  addTemplate("QCD_ZM_"+suffix , vars, wspace_ZM, (TH1F*)templatesfile->Get(("qbkghistzmm_"+observable).c_str()));
  addTemplate("Dibosons_ZM_"+suffix, vars, wspace_ZM, (TH1F*)templatesfile->Get(("dbkghistzmm_"+observable).c_str()));

  if(addShapeSystematics){

    TH1F* nominalHisto = (TH1F*)templatesfile->Get(("vlbkghistzmm_"+observable).c_str());
    TH1F* histobUp = (TH1F*)templatesfile->Get(("vlbkghistzmm_bUp_"+observable).c_str());
    TH1F* histobDw = (TH1F*)templatesfile->Get(("vlbkghistzmm_bDw_"+observable).c_str());
    TH1F* histoJesUp = (TH1F*)templatesfile->Get(("vlbkghistzmm_metJetUp_"+observable).c_str());
    TH1F* histoJesDw = (TH1F*)templatesfile->Get(("vlbkghistzmm_metJetDw_"+observable).c_str());
    TH1F* histoJerUp = (TH1F*)templatesfile->Get(("vlbkghistzmm_metResUp_"+observable).c_str());
    TH1F* histoJerDw = (TH1F*)templatesfile->Get(("vlbkghistzmm_metResDw_"+observable).c_str());
    TH1F* histoUncUp = (TH1F*)templatesfile->Get(("vlbkghistzmm_metUncUp_"+observable).c_str());
    TH1F* histoUncDw = (TH1F*)templatesfile->Get(("vlbkghistzmm_metUncDw_"+observable).c_str());

    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }

    addTemplate("WJets_ZM_"+suffix+"_CMS_btagUp", vars, wspace_ZM, histobUp);
    addTemplate("WJets_ZM_"+suffix+"_CMS_btagDown", vars, wspace_ZM, histobDw);
    addTemplate("WJets_ZM_"+suffix+"_CMS_scale_jUp", vars, wspace_ZM, histoJesUp);
    addTemplate("WJets_ZM_"+suffix+"_CMS_scale_jDown", vars, wspace_ZM, histoJesDw);
    addTemplate("WJets_ZM_"+suffix+"_CMS_res_jUp", vars, wspace_ZM, histoJerUp);
    addTemplate("WJets_ZM_"+suffix+"_CMS_res_jDown", vars, wspace_ZM, histoJerDw);
    addTemplate("WJets_ZM_"+suffix+"_CMS_uncUp", vars, wspace_ZM, histoUncUp);
    addTemplate("WJets_ZM_"+suffix+"_CMS_uncDown", vars, wspace_ZM, histoUncDw);

    ///
    nominalHisto = (TH1F*)templatesfile->Get(("tbkghistzmm_"+observable).c_str());
    histobUp = (TH1F*)templatesfile->Get(("tbkghistzmm_bUp_"+observable).c_str());
    histobDw = (TH1F*)templatesfile->Get(("tbkghistzmm_bDw_"+observable).c_str());
    histoJesUp = (TH1F*)templatesfile->Get(("tbkghistzmm_metJetUp_"+observable).c_str());
    histoJesDw = (TH1F*)templatesfile->Get(("tbkghistzmm_metJetDw_"+observable).c_str());
    histoJerUp = (TH1F*)templatesfile->Get(("tbkghistzmm_metResUp_"+observable).c_str());
    histoJerDw = (TH1F*)templatesfile->Get(("tbkghistzmm_metResDw_"+observable).c_str());
    histoUncUp = (TH1F*)templatesfile->Get(("tbkghistzmm_metUncUp_"+observable).c_str());
    histoUncDw = (TH1F*)templatesfile->Get(("tbkghistzmm_metUncDw_"+observable).c_str());

    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.10);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.90);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.03);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.97);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }

    addTemplate("Top_ZM_"+suffix+"_CMS_btagUp", vars, wspace_ZM, histobUp);
    addTemplate("Top_ZM_"+suffix+"_CMS_btagDown", vars, wspace_ZM, histobDw);
    addTemplate("Top_ZM_"+suffix+"_CMS_scale_jUp", vars, wspace_ZM, histoJesUp);
    addTemplate("Top_ZM_"+suffix+"_CMS_scale_jDown", vars, wspace_ZM, histoJesDw);
    addTemplate("Top_ZM_"+suffix+"_CMS_res_jUp", vars, wspace_ZM, histoJerUp);
    addTemplate("Top_ZM_"+suffix+"_CMS_res_jDown", vars, wspace_ZM, histoJerDw);
    addTemplate("Top_ZM_"+suffix+"_CMS_uncUp", vars, wspace_ZM, histoUncUp);
    addTemplate("Top_ZM_"+suffix+"_CMS_uncDown", vars, wspace_ZM, histoUncDw);

    ///
    nominalHisto = (TH1F*)templatesfile->Get(("dbkghistzmm_"+observable).c_str());
    histobUp = (TH1F*)templatesfile->Get(("dbkghistzmm_bUp_"+observable).c_str());
    histobDw = (TH1F*)templatesfile->Get(("dbkghistzmm_bDw_"+observable).c_str());
    histoJesUp = (TH1F*)templatesfile->Get(("dbkghistzmm_metJetUp_"+observable).c_str());
    histoJesDw = (TH1F*)templatesfile->Get(("dbkghistzmm_metJetDw_"+observable).c_str());
    histoJerUp = (TH1F*)templatesfile->Get(("dbkghistzmm_metResUp_"+observable).c_str());
    histoJerDw = (TH1F*)templatesfile->Get(("dbkghistzmm_metResDw_"+observable).c_str());
    histoUncUp = (TH1F*)templatesfile->Get(("dbkghistzmm_metUncUp_"+observable).c_str());
    histoUncDw = (TH1F*)templatesfile->Get(("dbkghistzmm_metUncDw_"+observable).c_str());

    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }


    addTemplate("Dibosons_ZM_"+suffix+"_CMS_btagUp", vars, wspace_ZM, histobUp);
    addTemplate("Dibosons_ZM_"+suffix+"_CMS_btagDown", vars, wspace_ZM, histobDw);
    addTemplate("Dibosons_ZM_"+suffix+"_CMS_scale_jUp", vars, wspace_ZM, histoJesUp);
    addTemplate("Dibosons_ZM_"+suffix+"_CMS_scale_jDown", vars, wspace_ZM, histoJesDw);
    addTemplate("Dibosons_ZM_"+suffix+"_CMS_res_jUp", vars, wspace_ZM, histoJerUp);
    addTemplate("Dibosons_ZM_"+suffix+"_CMS_res_jDown", vars, wspace_ZM, histoJerDw);
    addTemplate("Dibosons_ZM_"+suffix+"_CMS_uncUp", vars, wspace_ZM, histoUncUp);
    addTemplate("Dibosons_ZM_"+suffix+"_CMS_uncDown", vars, wspace_ZM, histoUncDw);

  }

  ////////////////////////////////////////
  // -------- CR Di-Electron  -------- //
  ///////////////////////////////////////
  cout<<"Make CR Di-Electron  templates ..."<<endl;

  RooWorkspace wspace_ZE(("ZE_"+suffix).c_str(),("ZE_"+suffix).c_str());

  addTemplate("data_obs_ZE_"+suffix, vars, wspace_ZE, (TH1F*)templatesfile->Get(("datahistzee_"+observable).c_str()));
  // Z->ee connected with Z->nunu SR
  vector<pair<RooRealVar*, TH1*> > znn_ZE_syst;
  makeConnectedBinList("Znunu_ZE_"+suffix, met, wspace_ZE, (TH1F*)templatesfile->Get(("zeecorhist_"+observable).c_str()), znn_ZE_syst, znn_SR_bins,NULL,observable);
  
  // Other MC backgrounds in dielectron control region
  addTemplate("WJets_ZE_"+suffix, vars, wspace_ZE, (TH1F*)templatesfile->Get(("vlbkghistzee_"+observable).c_str()));
  addTemplate("Top_ZE_"+suffix  , vars, wspace_ZE, (TH1F*)templatesfile->Get(("tbkghistzee_"+observable).c_str()));
  addTemplate("QCD_ZE_"+suffix  , vars, wspace_ZE, (TH1F*)templatesfile->Get(("qbkghistzee_"+observable).c_str()));
  addTemplate("Dibosons_ZE_"+suffix  , vars, wspace_ZE, (TH1F*)templatesfile->Get(("dbkghistzee_"+observable).c_str()));

  if( addShapeSystematics){

    TH1F* nominalHisto = (TH1F*)templatesfile->Get(("vlbkghistzee_"+observable).c_str());
    TH1F* histobUp = (TH1F*)templatesfile->Get(("vlbkghistzee_bUp_"+observable).c_str());
    TH1F* histobDw = (TH1F*)templatesfile->Get(("vlbkghistzee_bDw_"+observable).c_str());
    TH1F* histoJesUp = (TH1F*)templatesfile->Get(("vlbkghistzee_metJetUp_"+observable).c_str());
    TH1F* histoJesDw = (TH1F*)templatesfile->Get(("vlbkghistzee_metJetDw_"+observable).c_str());
    TH1F* histoJerUp = (TH1F*)templatesfile->Get(("vlbkghistzee_metResUp_"+observable).c_str());
    TH1F* histoJerDw = (TH1F*)templatesfile->Get(("vlbkghistzee_metResDw_"+observable).c_str());
    TH1F* histoUncUp = (TH1F*)templatesfile->Get(("vlbkghistzee_metUncUp_"+observable).c_str());
    TH1F* histoUncDw = (TH1F*)templatesfile->Get(("vlbkghistzee_metUncDw_"+observable).c_str());

    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }

    addTemplate("WJets_ZE_"+suffix+"_CMS_btagUp", vars, wspace_ZE, histobUp);
    addTemplate("WJets_ZE_"+suffix+"_CMS_btagDown", vars, wspace_ZE, histobDw);
    addTemplate("WJets_ZE_"+suffix+"_CMS_scale_jUp", vars, wspace_ZE, histoJesUp);
    addTemplate("WJets_ZE_"+suffix+"_CMS_scale_jDown", vars, wspace_ZE, histoJesDw);
    addTemplate("WJets_ZE_"+suffix+"_CMS_res_jUp", vars, wspace_ZE, histoJerUp);
    addTemplate("WJets_ZE_"+suffix+"_CMS_res_jDown", vars, wspace_ZE, histoJerDw);
    addTemplate("WJets_ZE_"+suffix+"_CMS_uncUp", vars, wspace_ZE, histoUncUp);
    addTemplate("WJets_ZE_"+suffix+"_CMS_uncDown", vars, wspace_ZE, histoUncDw);

    ///
    nominalHisto = (TH1F*)templatesfile->Get(("tbkghistzee_"+observable).c_str());
    histobUp = (TH1F*)templatesfile->Get(("tbkghistzee_bUp_"+observable).c_str());
    histobDw = (TH1F*)templatesfile->Get(("tbkghistzee_bDw_"+observable).c_str());
    histoJesUp = (TH1F*)templatesfile->Get(("tbkghistzee_metJetUp_"+observable).c_str());
    histoJesDw = (TH1F*)templatesfile->Get(("tbkghistzee_metJetDw_"+observable).c_str());
    histoJerUp = (TH1F*)templatesfile->Get(("tbkghistzee_metResUp_"+observable).c_str());
    histoJerDw = (TH1F*)templatesfile->Get(("tbkghistzee_metResDw_"+observable).c_str());
    histoUncUp = (TH1F*)templatesfile->Get(("tbkghistzee_metUncUp_"+observable).c_str());
    histoUncDw = (TH1F*)templatesfile->Get(("tbkghistzee_metUncDw_"+observable).c_str());

    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.10);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.90);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.03);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.97);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }

    addTemplate("Top_ZE_"+suffix+"_CMS_btagUp", vars, wspace_ZE, histobUp);
    addTemplate("Top_ZE_"+suffix+"_CMS_btagDown", vars, wspace_ZE, histobDw);
    addTemplate("Top_ZE_"+suffix+"_CMS_scale_jUp", vars, wspace_ZE, histoJesUp);
    addTemplate("Top_ZE_"+suffix+"_CMS_scale_jDown", vars, wspace_ZE, histoJesDw);
    addTemplate("Top_ZE_"+suffix+"_CMS_res_jUp", vars, wspace_ZE, histoJerUp);
    addTemplate("Top_ZE_"+suffix+"_CMS_res_jDown", vars, wspace_ZE, histoJerDw);
    addTemplate("Top_ZE_"+suffix+"_CMS_uncUp", vars, wspace_ZE, histoUncUp);
    addTemplate("Top_ZE_"+suffix+"_CMS_uncDown", vars, wspace_ZE, histoUncDw);

    ///
    nominalHisto = (TH1F*)templatesfile->Get(("dbkghistzee_"+observable).c_str());
    histobUp = (TH1F*)templatesfile->Get(("dbkghistzee_bUp_"+observable).c_str());
    histobDw = (TH1F*)templatesfile->Get(("dbkghistzee_bDw_"+observable).c_str());
    histoJesUp = (TH1F*)templatesfile->Get(("dbkghistzee_metJetUp_"+observable).c_str());
    histoJesDw = (TH1F*)templatesfile->Get(("dbkghistzee_metJetDw_"+observable).c_str());
    histoJerUp = (TH1F*)templatesfile->Get(("dbkghistzee_metResUp_"+observable).c_str());
    histoJerDw = (TH1F*)templatesfile->Get(("dbkghistzee_metResDw_"+observable).c_str());
    histoUncUp = (TH1F*)templatesfile->Get(("dbkghistzee_metUncUp_"+observable).c_str());
    histoUncDw = (TH1F*)templatesfile->Get(("dbkghistzee_metUncDw_"+observable).c_str());

    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }


    addTemplate("Dibosons_ZE_"+suffix+"_CMS_btagUp", vars, wspace_ZE, histobUp);
    addTemplate("Dibosons_ZE_"+suffix+"_CMS_btagDown", vars, wspace_ZE, histobDw);
    addTemplate("Dibosons_ZE_"+suffix+"_CMS_scale_jUp", vars, wspace_ZE, histoJesUp);
    addTemplate("Dibosons_ZE_"+suffix+"_CMS_scale_jDown", vars, wspace_ZE, histoJesDw);
    addTemplate("Dibosons_ZE_"+suffix+"_CMS_res_jUp", vars, wspace_ZE, histoJerUp);
    addTemplate("Dibosons_ZE_"+suffix+"_CMS_res_jDown", vars, wspace_ZE, histoJerDw);
    addTemplate("Dibosons_ZE_"+suffix+"_CMS_uncUp", vars, wspace_ZE, histoUncUp);
    addTemplate("Dibosons_ZE_"+suffix+"_CMS_uncDown", vars, wspace_ZE, histoUncDw);


  }

  ///////////////////////////////////////
  // -------- CR Gamma+jets  -------- //
  //////////////////////////////////////
  cout<<"Make CR Gamma+jets  templates ..."<<endl;
  RooWorkspace wspace_GJ(("GJ_"+suffix).c_str(),("GJ_"+suffix).c_str());

  addTemplate("data_obs_GJ_"+suffix, vars, wspace_GJ, (TH1F*)templatesfile->Get(("datahistgam_"+observable).c_str()));
  // Gamma+jets --> connected with Z->nunu
  vector<pair<RooRealVar*, TH1*> > znn_GJ_syst;
  RooRealVar* znn_GJ_re1 = new RooRealVar("Znunu_GJ_RenScale1"  , "", 0., -5., 5.);
  RooRealVar* znn_GJ_fa1 = new RooRealVar("Znunu_GJ_FactScale1" , "", 0., -5., 5.);
  RooRealVar* znn_GJ_re2 = new RooRealVar("Znunu_GJ_RenScale2"  , "", 0., -5., 5.);
  RooRealVar* znn_GJ_fa2 = new RooRealVar("Znunu_GJ_FactScale2" , "", 0., -5., 5.);
  RooRealVar* znn_GJ_pdf = new RooRealVar("Znunu_GJ_PDF"        , "", 0., -5., 5.);
  RooRealVar* znn_GJ_fpc = new RooRealVar("Znunu_GJ_Footprint"  , "", 0., -5., 5.);
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(NULL      , (TH1F*)templatesfile->Get(("ZG_EWK_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_re1, (TH1F*)templatesfile->Get(("ZG_RenScale1_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_fa1, (TH1F*)templatesfile->Get(("ZG_FactScale1_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_re2, (TH1F*)templatesfile->Get(("ZG_RenScale2_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_fa2, (TH1F*)templatesfile->Get(("ZG_FactScale2_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_pdf, (TH1F*)templatesfile->Get(("ZG_PDF_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_fpc, (TH1F*)templatesfile->Get(("ZG_Footprint_"+observable).c_str())));

  makeConnectedBinList("Znunu_GJ_"+suffix, met, wspace_GJ, (TH1F*)templatesfile->Get(("gamcorewkhist_"+observable).c_str()), znn_GJ_syst, znn_SR_bins,NULL,observable);  
  // Other MC backgrounds photon+jets control region
  addTemplate("QCD_GJ_"+suffix, vars, wspace_GJ, (TH1F*)templatesfile->Get(("qbkghistgam_"+observable).c_str()));
  
  ///////////////////////////////////////
  // -------- CR Single-Muon  -------- //
  //////////////////////////////////////
  cout<<"Make CR Single-Mu  templates ..."<<endl;

  RooWorkspace wspace_WM(("WM_"+suffix).c_str(),("WM_"+suffix).c_str());

  addTemplate("data_obs_WM_"+suffix, vars, wspace_WM, (TH1F*)templatesfile->Get(("datahistwmn_"+observable).c_str()));
  // connected W->munu with W+jets SR
  vector<pair<RooRealVar*, TH1*> > wln_WM_syst;
  makeConnectedBinList("WJets_WM_"+suffix, met, wspace_WM, (TH1F*)templatesfile->Get(("wmncorhist_"+observable).c_str()), wln_WM_syst, wln_SR_bins,NULL,observable);
		       
  // Other MC backgrounds in single muon control region
  addTemplate("ZJets_WM_"+suffix, vars, wspace_WM, (TH1F*)templatesfile->Get(("vllbkghistwmn_"+observable).c_str()));
  addTemplate("Top_WM_"+suffix  , vars, wspace_WM, (TH1F*)templatesfile->Get(("tbkghistwmn_"+observable).c_str()));
  addTemplate("QCD_WM_"+suffix  , vars, wspace_WM, (TH1F*)templatesfile->Get(("qbkghistwmn_"+observable).c_str()));
  addTemplate("Dibosons_WM_"+suffix, vars, wspace_WM, (TH1F*)templatesfile->Get(("dbkghistwmn_"+observable).c_str()));

  if(addShapeSystematics){

    TH1F* nominalHisto = (TH1F*)templatesfile->Get(("vllbkghistwmn_"+observable).c_str());
    TH1F* histobUp = (TH1F*)templatesfile->Get(("vllbkghistwmn_bUp_"+observable).c_str());
    TH1F* histobDw = (TH1F*)templatesfile->Get(("vllbkghistwmn_bDw_"+observable).c_str());
    TH1F* histoJesUp = (TH1F*)templatesfile->Get(("vllbkghistwmn_metJetUp_"+observable).c_str());
    TH1F* histoJesDw = (TH1F*)templatesfile->Get(("vllbkghistwmn_metJetDw_"+observable).c_str());
    TH1F* histoJerUp = (TH1F*)templatesfile->Get(("vllbkghistwmn_metResUp_"+observable).c_str());
    TH1F* histoJerDw = (TH1F*)templatesfile->Get(("vllbkghistwmn_metResDw_"+observable).c_str());
    TH1F* histoUncUp = (TH1F*)templatesfile->Get(("vllbkghistwmn_metUncUp_"+observable).c_str());
    TH1F* histoUncDw = (TH1F*)templatesfile->Get(("vllbkghistwmn_metUncDw_"+observable).c_str());

    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }
 
    addTemplate("ZJets_WM_"+suffix+"_CMS_btagUp", vars, wspace_WM, histobUp);
    addTemplate("ZJets_WM_"+suffix+"_CMS_btagDown", vars, wspace_WM, histobDw);
    addTemplate("ZJets_WM_"+suffix+"_CMS_scale_jUp", vars, wspace_WM, histoJesUp);
    addTemplate("ZJets_WM_"+suffix+"_CMS_scale_jDown", vars, wspace_WM, histoJesDw);
    addTemplate("ZJets_WM_"+suffix+"_CMS_res_jUp", vars, wspace_WM, histoJerUp);
    addTemplate("ZJets_WM_"+suffix+"_CMS_res_jDown", vars, wspace_WM, histoJerDw);
    addTemplate("ZJets_WM_"+suffix+"_CMS_uncUp", vars, wspace_WM, histoUncUp);
    addTemplate("ZJets_WM_"+suffix+"_CMS_uncDown", vars, wspace_WM, histoUncDw);

    ///
    nominalHisto = (TH1F*)templatesfile->Get(("tbkghistwmn_"+observable).c_str());
    histobUp = (TH1F*)templatesfile->Get(("tbkghistwmn_bUp_"+observable).c_str());
    histobDw = (TH1F*)templatesfile->Get(("tbkghistwmn_bDw_"+observable).c_str());
    histoJesUp = (TH1F*)templatesfile->Get(("tbkghistwmn_metJetUp_"+observable).c_str());
    histoJesDw = (TH1F*)templatesfile->Get(("tbkghistwmn_metJetDw_"+observable).c_str());
    histoJerUp = (TH1F*)templatesfile->Get(("tbkghistwmn_metResUp_"+observable).c_str());
    histoJerDw = (TH1F*)templatesfile->Get(("tbkghistwmn_metResDw_"+observable).c_str());
    histoUncUp = (TH1F*)templatesfile->Get(("tbkghistwmn_metUncUp_"+observable).c_str());
    histoUncDw = (TH1F*)templatesfile->Get(("tbkghistwmn_metUncDw_"+observable).c_str());

    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.10);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.90);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.03);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.97);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }

    addTemplate("Top_WM_"+suffix+"_CMS_btagUp", vars, wspace_WM, histobUp);
    addTemplate("Top_WM_"+suffix+"_CMS_btagDown", vars, wspace_WM, histobDw);
    addTemplate("Top_WM_"+suffix+"_CMS_scale_jUp", vars, wspace_WM, histoJesUp);
    addTemplate("Top_WM_"+suffix+"_CMS_scale_jDown", vars, wspace_WM, histoJesDw);
    addTemplate("Top_WM_"+suffix+"_CMS_res_jUp", vars, wspace_WM, histoJerUp);
    addTemplate("Top_WM_"+suffix+"_CMS_res_jDown", vars, wspace_WM, histoJerDw);
    addTemplate("Top_WM_"+suffix+"_CMS_uncUp", vars, wspace_WM, histoUncUp);
    addTemplate("Top_WM_"+suffix+"_CMS_uncDown", vars, wspace_WM, histoUncDw);

    ///
    nominalHisto = (TH1F*)templatesfile->Get(("dbkghistwmn_"+observable).c_str());
    histobUp = (TH1F*)templatesfile->Get(("dbkghistwmn_bUp_"+observable).c_str());
    histobDw = (TH1F*)templatesfile->Get(("dbkghistwmn_bDw_"+observable).c_str());
    histoJesUp = (TH1F*)templatesfile->Get(("dbkghistwmn_metJetUp_"+observable).c_str());
    histoJesDw = (TH1F*)templatesfile->Get(("dbkghistwmn_metJetDw_"+observable).c_str());
    histoJerUp = (TH1F*)templatesfile->Get(("dbkghistwmn_metResUp_"+observable).c_str());
    histoJerDw = (TH1F*)templatesfile->Get(("dbkghistwmn_metResDw_"+observable).c_str());
    histoUncUp = (TH1F*)templatesfile->Get(("dbkghistwmn_metUncUp_"+observable).c_str());
    histoUncDw = (TH1F*)templatesfile->Get(("dbkghistwmn_metUncDw_"+observable).c_str());

    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }


    addTemplate("Dibosons_WM_"+suffix+"_CMS_btagUp", vars, wspace_WM, histobUp);
    addTemplate("Dibosons_WM_"+suffix+"_CMS_btagDown", vars, wspace_WM, histobDw);
    addTemplate("Dibosons_WM_"+suffix+"_CMS_scale_jUp", vars, wspace_WM, histoJesUp);
    addTemplate("Dibosons_WM_"+suffix+"_CMS_scale_jDown", vars, wspace_WM, histoJesDw);
    addTemplate("Dibosons_WM_"+suffix+"_CMS_res_jUp", vars, wspace_WM, histoJerUp);
    addTemplate("Dibosons_WM_"+suffix+"_CMS_res_jDown", vars, wspace_WM, histoJerDw);
    addTemplate("Dibosons_WM_"+suffix+"_CMS_uncUp", vars, wspace_WM, histoUncUp);
    addTemplate("Dibosons_WM_"+suffix+"_CMS_uncDown", vars, wspace_WM, histoUncDw);

  }
  
  //////////////////////////////////..../////
  // -------- CR Single-Electron  -------- //
  //////////////////////////////////////////
  cout<<"Make CR Single-El  templates ..."<<endl;

  RooWorkspace wspace_WE(("WE_"+suffix).c_str(),("WE_"+suffix).c_str());

  addTemplate("data_obs_WE_"+suffix, vars, wspace_WE, (TH1F*)templatesfile->Get(("datahistwen_"+observable).c_str()));
  // connected W->enu with W+jets SR 
  vector<pair<RooRealVar*, TH1*> > wln_WE_syst;
  makeConnectedBinList("WJets_WE_"+suffix, met, wspace_WE, (TH1F*)templatesfile->Get(("wencorhist_"+observable).c_str()), wln_WE_syst, wln_SR_bins,NULL,observable);

  // Other MC backgrounds in single electron control region
  addTemplate("ZJets_WE_"+suffix, vars, wspace_WE, (TH1F*)templatesfile->Get(("vllbkghistwen_"+observable).c_str()));
  addTemplate("Top_WE_"+suffix, vars, wspace_WE, (TH1F*)templatesfile->Get(("tbkghistwen_"+observable).c_str()));
  addTemplate("QCD_WE_"+suffix, vars, wspace_WE, (TH1F*)templatesfile->Get(("qbkghistwen_"+observable).c_str()));
  addTemplate("Dibosons_WE_"+suffix, vars, wspace_WE, (TH1F*)templatesfile->Get(("dbkghistwen_"+observable).c_str()));

  if(addShapeSystematics){

    TH1F* nominalHisto = (TH1F*)templatesfile->Get(("vllbkghistwen_"+observable).c_str());
    TH1F* histobUp = (TH1F*)templatesfile->Get(("vllbkghistwen_bUp_"+observable).c_str());
    TH1F* histobDw = (TH1F*)templatesfile->Get(("vllbkghistwen_bDw_"+observable).c_str());
    TH1F* histoJesUp = (TH1F*)templatesfile->Get(("vllbkghistwen_metJetUp_"+observable).c_str());
    TH1F* histoJesDw = (TH1F*)templatesfile->Get(("vllbkghistwen_metJetDw_"+observable).c_str());
    TH1F* histoJerUp = (TH1F*)templatesfile->Get(("vllbkghistwen_metResUp_"+observable).c_str());
    TH1F* histoJerDw = (TH1F*)templatesfile->Get(("vllbkghistwen_metResDw_"+observable).c_str());
    TH1F* histoUncUp = (TH1F*)templatesfile->Get(("vllbkghistwen_metUncUp_"+observable).c_str());
    TH1F* histoUncDw = (TH1F*)templatesfile->Get(("vllbkghistwen_metUncDw_"+observable).c_str());

    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }

    addTemplate("ZJets_WE_"+suffix+"_CMS_btagUp", vars, wspace_WE, histobUp);
    addTemplate("ZJets_WE_"+suffix+"_CMS_btagDown", vars, wspace_WE, histobDw);
    addTemplate("ZJets_WE_"+suffix+"_CMS_scale_jUp", vars, wspace_WE, histoJesUp);
    addTemplate("ZJets_WE_"+suffix+"_CMS_scale_jDown", vars, wspace_WE, histoJesDw);
    addTemplate("ZJets_WE_"+suffix+"_CMS_res_jUp", vars, wspace_WE, histoJerUp);
    addTemplate("ZJets_WE_"+suffix+"_CMS_res_jDown", vars, wspace_WE, histoJerDw);
    addTemplate("ZJets_WE_"+suffix+"_CMS_uncUp", vars, wspace_WE, histoUncUp);
    addTemplate("ZJets_WE_"+suffix+"_CMS_uncDown", vars, wspace_WE, histoUncDw);

    ///
    nominalHisto = (TH1F*)templatesfile->Get(("tbkghistwen_"+observable).c_str());
    histobUp = (TH1F*)templatesfile->Get(("tbkghistwen_bUp_"+observable).c_str());
    histobDw = (TH1F*)templatesfile->Get(("tbkghistwen_bDw_"+observable).c_str());
    histoJesUp = (TH1F*)templatesfile->Get(("tbkghistwen_metJetUp_"+observable).c_str());
    histoJesDw = (TH1F*)templatesfile->Get(("tbkghistwen_metJetDw_"+observable).c_str());
    histoJerUp = (TH1F*)templatesfile->Get(("tbkghistwen_metResUp_"+observable).c_str());
    histoJerDw = (TH1F*)templatesfile->Get(("tbkghistwen_metResDw_"+observable).c_str());
    histoUncUp = (TH1F*)templatesfile->Get(("tbkghistwen_metUncUp_"+observable).c_str());
    histoUncDw = (TH1F*)templatesfile->Get(("tbkghistwen_metUncDw_"+observable).c_str());

    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.10);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.90);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.03);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.97);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }

    addTemplate("Top_WE_"+suffix+"_CMS_btagUp", vars, wspace_WE, histobUp);
    addTemplate("Top_WE_"+suffix+"_CMS_btagDown", vars, wspace_WE, histobDw);
    addTemplate("Top_WE_"+suffix+"_CMS_scale_jUp", vars, wspace_WE, histoJesUp);
    addTemplate("Top_WE_"+suffix+"_CMS_scale_jDown", vars, wspace_WE, histoJesDw);
    addTemplate("Top_WE_"+suffix+"_CMS_res_jUp", vars, wspace_WE, histoJerUp);
    addTemplate("Top_WE_"+suffix+"_CMS_res_jDown", vars, wspace_WE, histoJerDw);
    addTemplate("Top_WE_"+suffix+"_CMS_uncUp", vars, wspace_WE, histoUncUp);
    addTemplate("Top_WE_"+suffix+"_CMS_uncDown", vars, wspace_WE, histoUncDw);

    ///
    nominalHisto = (TH1F*)templatesfile->Get(("dbkghistwen_"+observable).c_str());
    histobUp = (TH1F*)templatesfile->Get(("dbkghistwen_bUp_"+observable).c_str());
    histobDw = (TH1F*)templatesfile->Get(("dbkghistwen_bDw_"+observable).c_str());
    histoJesUp = (TH1F*)templatesfile->Get(("dbkghistwen_metJetUp_"+observable).c_str());
    histoJesDw = (TH1F*)templatesfile->Get(("dbkghistwen_metJetDw_"+observable).c_str());
    histoJerUp = (TH1F*)templatesfile->Get(("dbkghistwen_metResUp_"+observable).c_str());
    histoJerDw = (TH1F*)templatesfile->Get(("dbkghistwen_metResDw_"+observable).c_str());
    histoUncUp = (TH1F*)templatesfile->Get(("dbkghistwen_metUncUp_"+observable).c_str());
    histoUncDw = (TH1F*)templatesfile->Get(("dbkghistwen_metUncDw_"+observable).c_str());

    if(nominalHisto->GetNbinsX() > 10){
      fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
    }

    addTemplate("Dibosons_WE_"+suffix+"_CMS_btagUp", vars, wspace_WE, histobUp);
    addTemplate("Dibosons_WE_"+suffix+"_CMS_btagDown", vars, wspace_WE, histobDw);
    addTemplate("Dibosons_WE_"+suffix+"_CMS_scale_jUp", vars, wspace_WE, histoJesUp);
    addTemplate("Dibosons_WE_"+suffix+"_CMS_scale_jDown", vars, wspace_WE, histoJesDw);
    addTemplate("Dibosons_WE_"+suffix+"_CMS_res_jUp", vars, wspace_WE, histoJerUp);
    addTemplate("Dibosons_WE_"+suffix+"_CMS_res_jDown", vars, wspace_WE, histoJerDw);
    addTemplate("Dibosons_WE_"+suffix+"_CMS_uncUp", vars, wspace_WE, histoUncUp);
    addTemplate("Dibosons_WE_"+suffix+"_CMS_uncDown", vars, wspace_WE, histoUncDw);

  }

  if(connectTop){
    
    /////////////////////////////////////
    // -------- CR Top-Muon  -------- //
    ////////////////////////////////////
    cout<<"Make CR Top-mu  templates ..."<<endl;
    RooWorkspace wspace_TM(("TM_"+suffix).c_str(),("TM_"+suffix).c_str());

    addTemplate("data_obs_TM_"+suffix, vars, wspace_TM, (TH1F*)templatesfile->Get(("datahisttopmu_"+observable).c_str()));
    // connect tt->mu+b with tt SR
    vector<pair<RooRealVar*, TH1*> > top_TM_syst;
    RooRealVar* top_btag   = new RooRealVar("Top_btag", "", 0., -5., 5.);
    top_TM_syst.push_back(pair<RooRealVar*, TH1*>(top_btag, (TH1F*)templatesfile->Get(("TOP_MU_B_"+observable).c_str())));
    
    makeConnectedBinList("Top_TM_"+suffix, met, wspace_TM, (TH1F*)templatesfile->Get(("topmucorhist_"+observable).c_str()), top_TM_syst, top_SR_bins,NULL,observable);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_TM_"+suffix, vars, wspace_TM, (TH1F*)templatesfile->Get(("vllbkghisttopmu_"+observable).c_str()));
    addTemplate("WJets_TM_"+suffix, vars, wspace_TM, (TH1F*)templatesfile->Get(("vlbkghisttopmu_"+observable).c_str()));
    addTemplate("QCD_TM_"+suffix, vars, wspace_TM, (TH1F*)templatesfile->Get(("qbkghisttopmu_"+observable).c_str()));
    addTemplate("Dibosons_TM_"+suffix, vars, wspace_TM, (TH1F*)templatesfile->Get(("dbkghisttopmu_"+observable).c_str()));

    if(addShapeSystematics){
      
      TH1F* nominalHisto = (TH1F*)templatesfile->Get(("vllbkghisttopmu_"+observable).c_str());
      TH1F* histobUp = (TH1F*)templatesfile->Get(("vllbkghisttopmu_bUp_"+observable).c_str());
      TH1F* histobDw = (TH1F*)templatesfile->Get(("vllbkghisttopmu_bDw_"+observable).c_str());
      TH1F* histoJesUp = (TH1F*)templatesfile->Get(("vllbkghisttopmu_metJetUp_"+observable).c_str());
      TH1F* histoJesDw = (TH1F*)templatesfile->Get(("vllbkghisttopmu_metJetDw_"+observable).c_str());
      TH1F* histoJerUp = (TH1F*)templatesfile->Get(("vllbkghisttopmu_metResUp_"+observable).c_str());
      TH1F* histoJerDw = (TH1F*)templatesfile->Get(("vllbkghisttopmu_metResDw_"+observable).c_str());
      TH1F* histoUncUp = (TH1F*)templatesfile->Get(("vllbkghisttopmu_metUncUp_"+observable).c_str());
      TH1F* histoUncDw = (TH1F*)templatesfile->Get(("vllbkghisttopmu_metUncDw_"+observable).c_str());
      
      if(nominalHisto->GetNbinsX() > 10){
	fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
      }
      
      addTemplate("ZJets_TM_"+suffix+"_CMS_btagUp", vars, wspace_TM, histobUp);
      addTemplate("ZJets_TM_"+suffix+"_CMS_btagDown", vars, wspace_TM, histobDw);
      addTemplate("ZJets_TM_"+suffix+"_CMS_scale_jUp", vars, wspace_TM, histoJesUp);
      addTemplate("ZJets_TM_"+suffix+"_CMS_scale_jDown", vars, wspace_TM, histoJesDw);
      addTemplate("ZJets_TM_"+suffix+"_CMS_res_jUp", vars, wspace_TM, histoJerUp);
      addTemplate("ZJets_TM_"+suffix+"_CMS_res_jDown", vars, wspace_TM, histoJerDw);
      addTemplate("ZJets_TM_"+suffix+"_CMS_uncUp", vars, wspace_TM, histoUncUp);
      addTemplate("ZJets_TM_"+suffix+"_CMS_uncDown", vars, wspace_TM, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->Get(("vlbkghisttopmu_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->Get(("vlbkghisttopmu_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->Get(("vlbkghisttopmu_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->Get(("vlbkghisttopmu_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->Get(("vlbkghisttopmu_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->Get(("vlbkghisttopmu_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->Get(("vlbkghisttopmu_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->Get(("vlbkghisttopmu_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->Get(("vlbkghisttopmu_metUncDw_"+observable).c_str());

      if(nominalHisto->GetNbinsX() > 10){
	fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
      }
      
      addTemplate("WJets_TM_"+suffix+"_CMS_btagUp", vars, wspace_TM, histobUp);
      addTemplate("WJets_TM_"+suffix+"_CMS_btagDown", vars, wspace_TM, histobDw);
      addTemplate("WJets_TM_"+suffix+"_CMS_scale_jUp", vars, wspace_TM, histoJesUp);
      addTemplate("WJets_TM_"+suffix+"_CMS_scale_jDown", vars, wspace_TM, histoJesDw);
      addTemplate("WJets_TM_"+suffix+"_CMS_res_jUp", vars, wspace_TM, histoJerUp);
      addTemplate("WJets_TM_"+suffix+"_CMS_res_jDown", vars, wspace_TM, histoJerDw);
      addTemplate("WJets_TM_"+suffix+"_CMS_uncUp", vars, wspace_TM, histoUncUp);
      addTemplate("WJets_TM_"+suffix+"_CMS_uncDown", vars, wspace_TM, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->Get(("dbkghisttopmu_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->Get(("dbkghisttopmu_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->Get(("dbkghisttopmu_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->Get(("dbkghisttopmu_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->Get(("dbkghisttopmu_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->Get(("dbkghisttopmu_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->Get(("dbkghisttopmu_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->Get(("dbkghisttopmu_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->Get(("dbkghisttopmu_metUncDw_"+observable).c_str());
      
      if(nominalHisto->GetNbinsX() > 10){
	fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
      }
      
      
      addTemplate("Dibosons_TM_"+suffix+"_CMS_btagUp", vars, wspace_TM, histobUp);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_btagDown", vars, wspace_TM, histobDw);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_scale_jUp", vars, wspace_TM, histoJesUp);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_scale_jDown", vars, wspace_TM, histoJesDw);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_res_jUp", vars, wspace_TM, histoJerUp);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_res_jDown", vars, wspace_TM, histoJerDw);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_uncUp", vars, wspace_TM, histoUncUp);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_uncDown", vars, wspace_TM, histoUncDw);
           
    }

    /////////////////////////////////////
    // -------- CR Top-Electron-------- //
    ////////////////////////////////////
    cout<<"Make CR Top-el  templates ..."<<endl;
    RooWorkspace wspace_TE(("TE_"+suffix).c_str(),("TE_"+suffix).c_str());

    addTemplate("data_obs_TE_"+suffix, vars, wspace_TE, (TH1F*)templatesfile->Get(("datahisttopel_"+observable).c_str()));
    // connect tt->mu+b with tt SR
    vector<pair<RooRealVar*, TH1*> > top_TE_syst;
    top_TE_syst.push_back(pair<RooRealVar*, TH1*>(top_btag, (TH1F*)templatesfile->Get(("TOP_EL_B_"+observable).c_str())));
    
    makeConnectedBinList("Top_TE_"+suffix, met, wspace_TE, (TH1F*)templatesfile->Get(("topelcorhist_"+observable).c_str()), top_TE_syst, top_SR_bins,NULL,observable);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_TE_"+suffix , vars, wspace_TE, (TH1F*)templatesfile->Get(("vllbkghisttopel_"+observable).c_str()));
    addTemplate("WJets_TE_"+suffix , vars, wspace_TE, (TH1F*)templatesfile->Get(("vlbkghisttopel_"+observable).c_str()));
    addTemplate("QCD_TE_"+suffix   , vars, wspace_TE, (TH1F*)templatesfile->Get(("qbkghisttopel_"+observable).c_str()));
    addTemplate("Dibosons_TE_"+suffix, vars, wspace_TE, (TH1F*)templatesfile->Get(("dbkghisttopel_"+observable).c_str()));

    if(addShapeSystematics){
      
      TH1F* nominalHisto = (TH1F*)templatesfile->Get(("vllbkghisttopel_"+observable).c_str());
      TH1F* histobUp = (TH1F*)templatesfile->Get(("vllbkghisttopel_bUp_"+observable).c_str());
      TH1F* histobDw = (TH1F*)templatesfile->Get(("vllbkghisttopel_bDw_"+observable).c_str());
      TH1F* histoJesUp = (TH1F*)templatesfile->Get(("vllbkghisttopel_metJetUp_"+observable).c_str());
      TH1F* histoJesDw = (TH1F*)templatesfile->Get(("vllbkghisttopel_metJetDw_"+observable).c_str());
      TH1F* histoJerUp = (TH1F*)templatesfile->Get(("vllbkghisttopel_metResUp_"+observable).c_str());
      TH1F* histoJerDw = (TH1F*)templatesfile->Get(("vllbkghisttopel_metResDw_"+observable).c_str());
      TH1F* histoUncUp = (TH1F*)templatesfile->Get(("vllbkghisttopel_metUncUp_"+observable).c_str());
      TH1F* histoUncDw = (TH1F*)templatesfile->Get(("vllbkghisttopel_metUncDw_"+observable).c_str());
      
      if(nominalHisto->GetNbinsX() > 10){
	fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
      }
      
      addTemplate("ZJets_TE_"+suffix+"_CMS_btagUp", vars, wspace_TE, histobUp);
      addTemplate("ZJets_TE_"+suffix+"_CMS_btagDown", vars, wspace_TE, histobDw);
      addTemplate("ZJets_TE_"+suffix+"_CMS_scale_jUp", vars, wspace_TE, histoJesUp);
      addTemplate("ZJets_TE_"+suffix+"_CMS_scale_jDown", vars, wspace_TE, histoJesDw);
      addTemplate("ZJets_TE_"+suffix+"_CMS_res_jUp", vars, wspace_TE, histoJerUp);
      addTemplate("ZJets_TE_"+suffix+"_CMS_res_jDown", vars, wspace_TE, histoJerDw);
      addTemplate("ZJets_TE_"+suffix+"_CMS_uncUp", vars, wspace_TE, histoUncUp);
      addTemplate("ZJets_TE_"+suffix+"_CMS_uncDown", vars, wspace_TE, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->Get(("vlbkghisttopel_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->Get(("vlbkghisttopel_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->Get(("vlbkghisttopel_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->Get(("vlbkghisttopel_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->Get(("vlbkghisttopel_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->Get(("vlbkghisttopel_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->Get(("vlbkghisttopel_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->Get(("vlbkghisttopel_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->Get(("vlbkghisttopel_metUncDw_"+observable).c_str());

      if(nominalHisto->GetNbinsX() > 10){
	fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
      }
      
      addTemplate("WJets_TE_"+suffix+"_CMS_btagUp", vars, wspace_TE, histobUp);
      addTemplate("WJets_TE_"+suffix+"_CMS_btagDown", vars, wspace_TE, histobDw);
      addTemplate("WJets_TE_"+suffix+"_CMS_scale_jUp", vars, wspace_TE, histoJesUp);
      addTemplate("WJets_TE_"+suffix+"_CMS_scale_jDown", vars, wspace_TE, histoJesDw);
      addTemplate("WJets_TE_"+suffix+"_CMS_res_jUp", vars, wspace_TE, histoJerUp);
      addTemplate("WJets_TE_"+suffix+"_CMS_res_jDown", vars, wspace_TE, histoJerDw);
      addTemplate("WJets_TE_"+suffix+"_CMS_uncUp", vars, wspace_TE, histoUncUp);
      addTemplate("WJets_TE_"+suffix+"_CMS_uncDown", vars, wspace_TE, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->Get(("dbkghisttopel_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->Get(("dbkghisttopel_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->Get(("dbkghisttopel_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->Get(("dbkghisttopel_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->Get(("dbkghisttopel_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->Get(("dbkghisttopel_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->Get(("dbkghisttopel_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->Get(("dbkghisttopel_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->Get(("dbkghisttopel_metUncDw_"+observable).c_str());
      
      if(nominalHisto->GetNbinsX() > 10){
	fixShapeUncertainty(nominalHisto,histobUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histobDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoJesUp,500.,1.06);
	fixShapeUncertainty(nominalHisto,histoJesDw,500.,0.94);
	fixShapeUncertainty(nominalHisto,histoJerUp,500.,1.02);
	fixShapeUncertainty(nominalHisto,histoJerDw,500.,0.98);
	fixShapeUncertainty(nominalHisto,histoUncUp,500.,1.01);
	fixShapeUncertainty(nominalHisto,histoUncDw,500.,0.99);
      }
      
      
      addTemplate("Dibosons_TE_"+suffix+"_CMS_btagUp", vars, wspace_TE, histobUp);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_btagDown", vars, wspace_TE, histobDw);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_scale_jUp", vars, wspace_TE, histoJesUp);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_scale_jDown", vars, wspace_TE, histoJesDw);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_res_jUp", vars, wspace_TE, histoJerUp);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_res_jDown", vars, wspace_TE, histoJerDw);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_uncUp", vars, wspace_TE, histoUncUp);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_uncDown", vars, wspace_TE, histoUncDw);
      
      
    }


  }

  // ---------------------------- Write out the workspace -----------------------------------------------------------------//
  outfile->cd();
  wspace_SR.Write();
  wspace_ZM.Write();
  wspace_ZE.Write();
  wspace_WM.Write();
  wspace_WE.Write();
  wspace_GJ.Write();
  if(connectTop){
    wspace_TE.Write();
    wspace_TM.Write();
  }
  outfile->Close();
  return;
}
