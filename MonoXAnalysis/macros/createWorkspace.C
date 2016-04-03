#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <utility>
#include <sstream>

#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooAddition.h"

using namespace std;
#include "histoUtils.h"
#include "../../../HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#include "TSystem.h"

void addDummyBinContent(TH1* histo){

  cerr<<"addDummyBinContent: called for histo "<<histo->GetName()<<endl;
  for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
    histo->SetBinContent(iBin+1,10.e-6);
    histo->SetBinError(iBin+1,10.e-6);
  }

}


// function to create a RooDataHist from TH1F and import it in a workspace
void addTemplate(string procname, RooArgList& varlist, RooWorkspace& ws, TH1F* hist) {  
  
  if(hist == 0 || hist == NULL){
    cerr<<"addTemplate: found an null pointer --> check "<<endl;
    return;
  }
  if(hist->Integral() == 0)
    addDummyBinContent(hist); // avoind empty histograms in the workspace

  RooDataHist rhist((procname).c_str(), "", varlist, hist);
  ws.import(rhist);
}

void generateStatTemplate(string procname, RooArgList& varlist, RooWorkspace& ws, TH1* histo, float scaleUncertainty = 1.){

  vector<TH1F*> histStatUp;
  vector<TH1F*> histStatDw;

  if(histo == 0 || histo == NULL){
    cerr<<"generateStatTemplate: found an null pointer --> check "<<endl;
    return;
  }
  if(histo->Integral() == 0){
    addDummyBinContent(histo);
  }

  for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
    histStatUp.push_back((TH1F*) histo->Clone(Form("%s_%s_bin_%i_statUp_tmp",procname.c_str(),procname.c_str(),iBin+1)));
    histStatDw.push_back((TH1F*) histo->Clone(Form("%s_%s_bin_%i_statDown_tmp",procname.c_str(),procname.c_str(),iBin+1)));
  }

  for( size_t iHisto =0; iHisto < histStatUp.size(); iHisto++){
    histStatUp.at(iHisto)->SetBinContent(iHisto+1,histo->GetBinContent(iHisto+1)+histo->GetBinError(iHisto+1)*scaleUncertainty);
    if(histo->GetBinContent(iHisto+1)-histo->GetBinError(iHisto+1)*scaleUncertainty >= 0.)
      histStatDw.at(iHisto)->SetBinContent(iHisto+1,histo->GetBinContent(iHisto+1)-histo->GetBinError(iHisto+1)*scaleUncertainty);
    else
      histStatDw.at(iHisto)->SetBinContent(iHisto+1,0.);
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

  if(hist == 0 || hist == NULL){
    cerr<<"makeBinList: found an null pointer --> check "<<endl;
    return;
  }
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
  

  if(rhist == 0 || rhist == NULL){
    cerr<<"makeConnectedBinList: found an null pointer --> check "<<endl;
    return;
  }
  if(rhist->Integral() == 0){
    addDummyBinContent(rhist);
  }

  // bin list for the CR
  if (crbinlist == NULL){
    cerr<<"makeConnectedBinList: crbinlist empty --> create "<<endl;
    crbinlist = new RooArgList();
  }

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
    float extreme = min(5.,0.9*rhist->GetBinContent(i)/rhist->GetBinError(i));
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
    for (size_t j = 0; j < syst.size(); j++) {
      stringstream systbinss;
      if (syst[j].first == NULL) { // add bin by bin
	systbinss << procname << "_bin" << i << "_" << syst[j].second->GetName();
	TString nameSys (systbinss.str());
	nameSys.ReplaceAll(("_"+observable).c_str(),"");
	float extreme = min(5.,0.9*rhist->GetBinContent(i)/syst[j].second->GetBinContent(i)); 
	RooRealVar* systbinvar = new RooRealVar(nameSys.Data(), "", 0., -extreme, extreme);
	// Add all the systeamtics as new Multiplicative Nuisance for each bin
	fobinlist.add(*systbinvar);
      }
      else{ 
	float extreme = min(5.,0.9*rhist->GetBinContent(i)/syst[j].second->GetBinContent(i));
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
		     bool   mergeLeptons  = false,
		     string interaction   = "Vector",
		     string mediatorMass  = "1000", 
		     string DMMass        = "50"
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
  addTemplate("data_obs_SR_"+suffix, vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("datahist_"+observable).c_str()));

  // Signal shape
  if(!isHiggsInvisible){

    addTemplate("MonoJ_SR_"+suffix, vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
    addTemplate("MonoW_SR_"+suffix, vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
    addTemplate("MonoZ_SR_"+suffix, vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));

    if(addShapeSystematics){

      // monoJ
      TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histobUp  = (TH1F*)templatesfile->FindObjectAny(("monoJhist_bUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histobDw  = (TH1F*)templatesfile->FindObjectAny(("monoJhist_bDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histoJesUp  = (TH1F*)templatesfile->FindObjectAny(("monoJhist_metJetUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histoJesDw  = (TH1F*)templatesfile->FindObjectAny(("monoJhist_metJetDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histoJerUp  = (TH1F*)templatesfile->FindObjectAny(("monoJhist_metResUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histoJerDw  = (TH1F*)templatesfile->FindObjectAny(("monoJhist_metResDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histoUncUp  = (TH1F*)templatesfile->FindObjectAny(("monoJhist_metUncUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      TH1F* histoUncDw  = (TH1F*)templatesfile->FindObjectAny(("monoJhist_metUncDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      
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
      
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histobUp  = (TH1F*)templatesfile->FindObjectAny(("monoWhist_bUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histobDw  = (TH1F*)templatesfile->FindObjectAny(("monoWhist_bDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJesUp  = (TH1F*)templatesfile->FindObjectAny(("monoWhist_metJetUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJesDw  = (TH1F*)templatesfile->FindObjectAny(("monoWhist_metJetDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJerUp  = (TH1F*)templatesfile->FindObjectAny(("monoWhist_metResUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJerDw  = (TH1F*)templatesfile->FindObjectAny(("monoWhist_metResDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoUncUp  = (TH1F*)templatesfile->FindObjectAny(("monoWhist_metUncUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoUncDw  = (TH1F*)templatesfile->FindObjectAny(("monoWhist_metUncDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      
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
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histobUp  = (TH1F*)templatesfile->FindObjectAny(("monoZhist_bUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histobDw  = (TH1F*)templatesfile->FindObjectAny(("monoZhist_bDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJesUp  = (TH1F*)templatesfile->FindObjectAny(("monoZhist_metJetUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJesDw  = (TH1F*)templatesfile->FindObjectAny(("monoZhist_metJetDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJerUp  = (TH1F*)templatesfile->FindObjectAny(("monoZhist_metResUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoJerDw  = (TH1F*)templatesfile->FindObjectAny(("monoZhist_metResDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoUncUp  = (TH1F*)templatesfile->FindObjectAny(("monoZhist_metUncUp_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      histoUncDw  = (TH1F*)templatesfile->FindObjectAny(("monoZhist_metUncDw_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      
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
      generateStatTemplate("MonoJ_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
      generateStatTemplate("MonoW_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
      generateStatTemplate("MonoZ_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()));
    }
  }
  else{

    addTemplate("ggH_SR_"+suffix, vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+observable).c_str()));
    addTemplate("vbfH_SR_"+suffix,vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("vbfHhist_"+mediatorMass+"_"+observable).c_str()));
    addTemplate("wH_SR_"+suffix,  vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("wHhist_"+mediatorMass+"_"+observable).c_str()));
    addTemplate("zH_SR_"+suffix,  vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("zHhist_"+mediatorMass+"_"+observable).c_str()));

    if(addShapeSystematics){

      // ggH
      TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+observable).c_str());
      TH1F* histobUp  = (TH1F*)templatesfile->FindObjectAny(("ggHhist_bUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histobDw  = (TH1F*)templatesfile->FindObjectAny(("ggHhist_bDw_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoJesUp  = (TH1F*)templatesfile->FindObjectAny(("ggHhist_metJetUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoJesDw  = (TH1F*)templatesfile->FindObjectAny(("ggHhist_metJetDw_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoJerUp  = (TH1F*)templatesfile->FindObjectAny(("ggHhist_metResUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoJerDw  = (TH1F*)templatesfile->FindObjectAny(("ggHhist_metResDw_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoUncUp  = (TH1F*)templatesfile->FindObjectAny(("ggHhist_metUncUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoUncDw  = (TH1F*)templatesfile->FindObjectAny(("ggHhist_metUncDw_"+mediatorMass+"_"+observable).c_str());

      TH1F* histoRenUp = (TH1F*)templatesfile->FindObjectAny(("ggHhist_renUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoRenDw = (TH1F*)templatesfile->FindObjectAny(("ggHhist_renDw_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoFacUp = (TH1F*)templatesfile->FindObjectAny(("ggHhist_facUp_"+mediatorMass+"_"+observable).c_str());
      TH1F* histoFacDw = (TH1F*)templatesfile->FindObjectAny(("ggHhist_facDw_"+mediatorMass+"_"+observable).c_str());

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
      
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("vbfHhist_"+mediatorMass+"_"+observable).c_str());
      histobUp  = (TH1F*)templatesfile->FindObjectAny(("vbfHhist_bUp_"+mediatorMass+"_"+observable).c_str());
      histobDw  = (TH1F*)templatesfile->FindObjectAny(("vbfHhist_bDw_"+mediatorMass+"_"+observable).c_str());
      histoJesUp  = (TH1F*)templatesfile->FindObjectAny(("vbfHhist_metJetUp_"+mediatorMass+"_"+observable).c_str());
      histoJesDw  = (TH1F*)templatesfile->FindObjectAny(("vbfHhist_metJetDw_"+mediatorMass+"_"+observable).c_str());
      histoJerUp  = (TH1F*)templatesfile->FindObjectAny(("vbfHhist_metResUp_"+mediatorMass+"_"+observable).c_str());
      histoJerDw  = (TH1F*)templatesfile->FindObjectAny(("vbfHhist_metResDw_"+mediatorMass+"_"+observable).c_str());
      histoUncUp  = (TH1F*)templatesfile->FindObjectAny(("vbfHhist_metUncUp_"+mediatorMass+"_"+observable).c_str());
      histoUncDw  = (TH1F*)templatesfile->FindObjectAny(("vbfHhist_metUncDw_"+mediatorMass+"_"+observable).c_str());
      
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
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("wHhist_"+mediatorMass+"_"+observable).c_str());
      histobUp  = (TH1F*)templatesfile->FindObjectAny(("wHhist_bUp_"+mediatorMass+"_"+observable).c_str());
      histobDw  = (TH1F*)templatesfile->FindObjectAny(("wHhist_bDw_"+mediatorMass+"_"+observable).c_str());
      histoJesUp  = (TH1F*)templatesfile->FindObjectAny(("wHhist_metJetUp_"+mediatorMass+"_"+observable).c_str());
      histoJesDw  = (TH1F*)templatesfile->FindObjectAny(("wHhist_metJetDw_"+mediatorMass+"_"+observable).c_str());
      histoJerUp  = (TH1F*)templatesfile->FindObjectAny(("wHhist_metResUp_"+mediatorMass+"_"+observable).c_str());
      histoJerDw  = (TH1F*)templatesfile->FindObjectAny(("wHhist_metResDw_"+mediatorMass+"_"+observable).c_str());
      histoUncUp  = (TH1F*)templatesfile->FindObjectAny(("wHhist_metUncUp_"+mediatorMass+"_"+observable).c_str());
      histoUncDw  = (TH1F*)templatesfile->FindObjectAny(("wHhist_metUncDw_"+mediatorMass+"_"+observable).c_str());
      
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
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("zHhist_"+mediatorMass+"_"+observable).c_str());
      histobUp  = (TH1F*)templatesfile->FindObjectAny(("zHhist_bUp_"+mediatorMass+"_"+observable).c_str());
      histobDw  = (TH1F*)templatesfile->FindObjectAny(("zHhist_bDw_"+mediatorMass+"_"+observable).c_str());
      histoJesUp  = (TH1F*)templatesfile->FindObjectAny(("zHhist_metJetUp_"+mediatorMass+"_"+observable).c_str());
      histoJesDw  = (TH1F*)templatesfile->FindObjectAny(("zHhist_metJetDw_"+mediatorMass+"_"+observable).c_str());
      histoJerUp  = (TH1F*)templatesfile->FindObjectAny(("zHhist_metResUp_"+mediatorMass+"_"+observable).c_str());
      histoJerDw  = (TH1F*)templatesfile->FindObjectAny(("zHhist_metResDw_"+mediatorMass+"_"+observable).c_str());
      histoUncUp  = (TH1F*)templatesfile->FindObjectAny(("zHhist_metUncUp_"+mediatorMass+"_"+observable).c_str());
      histoUncDw  = (TH1F*)templatesfile->FindObjectAny(("zHhist_metUncDw_"+mediatorMass+"_"+observable).c_str());
      
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
    generateStatTemplate("ggH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
    generateStatTemplate("vbfH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("vbfHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
    generateStatTemplate("wH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("wHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
    generateStatTemplate("zH_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("zHhist_"+mediatorMass+"_"+observable).c_str()),0.5);
  }
  
  // Zvv background --> to be extracted from CRs
  TH1F* znn_SR_hist = (TH1F*) templatesfile->FindObjectAny(("zinvhist_"+observable).c_str());
  RooArgList znn_SR_bins; 
  // create a RooParametric hist with one RooRealVar per bin 
  makeBinList("Znunu_SR_"+suffix, met, wspace_SR, znn_SR_hist, znn_SR_bins);

  // Top background --> to be extracted from CRs
  RooArgList top_SR_bins;
  TH1F* top_SR_hist = NULL;

  // for data driven top estimation
  if(connectTop){
    top_SR_hist = (TH1F*) templatesfile->FindObjectAny(("tbkghist_"+observable).c_str());
    RooArgList top_SR_bins; 
    makeBinList("Top_SR_"+suffix, met, wspace_SR, top_SR_hist, top_SR_bins);
  }
  else{ // rely on MC + systematics

    addTemplate("Top_SR_"+suffix,vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("tbkghist_"+observable).c_str()));

    if(addShapeSystematics){

      TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("tbkghist_"+observable).c_str());
      TH1F* histobUp = (TH1F*)templatesfile->FindObjectAny(("tbkghist_bUp_"+observable).c_str());
      TH1F* histobDw = (TH1F*)templatesfile->FindObjectAny(("tbkghist_bDw_"+observable).c_str());
      TH1F* histoJesUp = (TH1F*)templatesfile->FindObjectAny(("tbkghist_metJetUp_"+observable).c_str());
      TH1F* histoJesDw = (TH1F*)templatesfile->FindObjectAny(("tbkghist_metJetDw_"+observable).c_str());
      TH1F* histoJerUp = (TH1F*)templatesfile->FindObjectAny(("tbkghist_metResUp_"+observable).c_str());
      TH1F* histoJerDw = (TH1F*)templatesfile->FindObjectAny(("tbkghist_metResDw_"+observable).c_str());
      TH1F* histoUncUp = (TH1F*)templatesfile->FindObjectAny(("tbkghist_metUncUp_"+observable).c_str());
      TH1F* histoUncDw = (TH1F*)templatesfile->FindObjectAny(("tbkghist_metUncDw_"+observable).c_str());

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
      
      generateStatTemplate("Top_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("tbkghist_"+observable).c_str()));
    }
  }

  // WJets background --> to be extracted from CRs, with connection to Z->nunu
  TH1F* wln_SR_hist = (TH1F*) templatesfile->FindObjectAny(("wjethist_"+observable).c_str());
  RooArgList wln_SR_bins;
  // set of correlated systematic uncertainties for the Z/W ratio
  vector<pair<RooRealVar*, TH1*> > wln_SR_syst;

  RooRealVar* wln_SR_re1 = new RooRealVar("WJets_SR_RenScale1", ""  , 0., -5., 5.);
  RooRealVar* wln_SR_fa1 = new RooRealVar("WJets_SR_FactScale1", "" , 0., -5., 5.);
  RooRealVar* wln_SR_re2 = new RooRealVar("WJets_SR_RenScale2", ""  , 0., -5., 5.);
  RooRealVar* wln_SR_fa2 = new RooRealVar("WJets_SR_FactScale2", "" , 0., -5., 5.);
  RooRealVar* wln_SR_pdf = new RooRealVar("WJets_SR_PDF", ""        , 0., -5., 5.);

  // NULL means bin-by-bin
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(NULL      , (TH1F*)templatesfile->FindObjectAny(("ZW_EWK_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_re1, (TH1F*)templatesfile->FindObjectAny(("ZW_RenScale1_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_fa1, (TH1F*)templatesfile->FindObjectAny(("ZW_FactScale1_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_re2, (TH1F*)templatesfile->FindObjectAny(("ZW_RenScale2_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_fa2, (TH1F*)templatesfile->FindObjectAny(("ZW_FactScale2_"+observable).c_str())));
  wln_SR_syst.push_back(pair<RooRealVar*, TH1*>(wln_SR_pdf, (TH1F*)templatesfile->FindObjectAny(("ZW_PDF_"+observable).c_str())));

  if (!connectWZ) 
    makeBinList("WJets_SR_"+suffix, met, wspace_SR, wln_SR_hist, wln_SR_bins);
  else   
    makeConnectedBinList("WJets_SR_"+suffix, met, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("zwjcorewkhist_"+observable).c_str()), wln_SR_syst, znn_SR_bins, &wln_SR_bins, observable);

  // Other MC backgrounds
  addTemplate("ZJets_SR_"+suffix     , vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("zjethist_"+observable).c_str()));
  addTemplate("Dibosons_SR_"+suffix  , vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("dbkghist_"+observable).c_str()));
  addTemplate("GJets_SR_"+suffix     , vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("gbkghist_"+observable).c_str()));

  if(addShapeSystematics){

    TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("zjethist_"+observable).c_str());
    TH1F* histobUp = (TH1F*)templatesfile->FindObjectAny(("zjethist_bUp_"+observable).c_str());
    TH1F* histobDw = (TH1F*)templatesfile->FindObjectAny(("zjethist_bDw_"+observable).c_str());
    TH1F* histoJesUp = (TH1F*)templatesfile->FindObjectAny(("zjethist_metJetUp_"+observable).c_str());
    TH1F* histoJesDw = (TH1F*)templatesfile->FindObjectAny(("zjethist_metJetDw_"+observable).c_str());
    TH1F* histoJerUp = (TH1F*)templatesfile->FindObjectAny(("zjethist_metResUp_"+observable).c_str());
    TH1F* histoJerDw = (TH1F*)templatesfile->FindObjectAny(("zjethist_metResDw_"+observable).c_str());
    TH1F* histoUncUp = (TH1F*)templatesfile->FindObjectAny(("zjethist_metUncUp_"+observable).c_str());
    TH1F* histoUncDw = (TH1F*)templatesfile->FindObjectAny(("zjethist_metUncDw_"+observable).c_str());
    
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

    generateStatTemplate("ZJets_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("zjethist_"+observable).c_str()));

    nominalHisto = (TH1F*)templatesfile->FindObjectAny(("dbkghist_"+observable).c_str());
    histobUp = (TH1F*)templatesfile->FindObjectAny(("dbkghist_bUp_"+observable).c_str());
    histobDw = (TH1F*)templatesfile->FindObjectAny(("dbkghist_bDw_"+observable).c_str());
    histoJesUp = (TH1F*)templatesfile->FindObjectAny(("dbkghist_metJetUp_"+observable).c_str());
    histoJesDw = (TH1F*)templatesfile->FindObjectAny(("dbkghist_metJetDw_"+observable).c_str());
    histoJerUp = (TH1F*)templatesfile->FindObjectAny(("dbkghist_metResUp_"+observable).c_str());
    histoJerDw = (TH1F*)templatesfile->FindObjectAny(("dbkghist_metResDw_"+observable).c_str());
    histoUncUp = (TH1F*)templatesfile->FindObjectAny(("dbkghist_metUncUp_"+observable).c_str());
    histoUncDw = (TH1F*)templatesfile->FindObjectAny(("dbkghist_metUncDw_"+observable).c_str());
    
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

    generateStatTemplate("Dibosons_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("dbkghist_"+observable).c_str()));


    nominalHisto = (TH1F*)templatesfile->FindObjectAny(("gbkghist_"+observable).c_str());
    histobUp = (TH1F*)templatesfile->FindObjectAny(("gbkghist_bUp_"+observable).c_str());
    histobDw = (TH1F*)templatesfile->FindObjectAny(("gbkghist_bDw_"+observable).c_str());
    histoJesUp = (TH1F*)templatesfile->FindObjectAny(("gbkghist_metJetUp_"+observable).c_str());
    histoJesDw = (TH1F*)templatesfile->FindObjectAny(("gbkghist_metJetDw_"+observable).c_str());
    histoJerUp = (TH1F*)templatesfile->FindObjectAny(("gbkghist_metResUp_"+observable).c_str());
    histoJerDw = (TH1F*)templatesfile->FindObjectAny(("gbkghist_metResDw_"+observable).c_str());
    histoUncUp = (TH1F*)templatesfile->FindObjectAny(("gbkghist_metUncUp_"+observable).c_str());
    histoUncDw = (TH1F*)templatesfile->FindObjectAny(("gbkghist_metUncDw_"+observable).c_str());
    
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

    generateStatTemplate("GJets_SR_"+suffix,vars,wspace_SR,(TH1F*)templatesfile->FindObjectAny(("gbkghist_"+observable).c_str()));

  }

  // look for DD qcd background, otherwise MC scaled by a factor 2
  TH1F* qcdhist = (TH1F*)templatesfile->FindObjectAny(("qbkghistDD_"+observable).c_str());
  if(qcdhist){
    addTemplate("QCD_SR_"+suffix, vars, wspace_SR, qcdhist);
    addTemplate("QCD_SR_"+suffix+"_CMS_qcdUp", vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("qbkghistDD_shapeUp_"+observable).c_str()));
    addTemplate("QCD_SR_"+suffix+"_CMS_qcdDown", vars, wspace_SR, (TH1F*)templatesfile->FindObjectAny(("qbkghistDD_shapeDw_"+observable).c_str()));
  }
  else{
    qcdhist = (TH1F*)templatesfile->FindObjectAny(("qbkghist_"+observable).c_str());
    qcdhist->Scale(scaleQCD);  
    addTemplate("QCD_SR_"+suffix, vars, wspace_SR, qcdhist);
  }


  ////////////////////////////////////
  // -------- CR Di-Muon  -------- //
  /////////////////////////////////// 
  RooWorkspace* wspace_ZM = NULL;
  RooWorkspace* wspace_ZE = NULL;
  RooWorkspace* wspace_WM = NULL;
  RooWorkspace* wspace_WE = NULL;
  RooWorkspace* wspace_ZL = NULL;
  RooWorkspace* wspace_WL = NULL;
  
  if(not mergeLeptons){

    cout<<"Make CR Di-Muon  templates ..."<<endl;    
    wspace_ZM = new RooWorkspace(("ZM_"+suffix).c_str(),("ZM_"+suffix).c_str());
  
    addTemplate("data_obs_ZM_"+suffix, vars, *wspace_ZM, (TH1F*)templatesfile->FindObjectAny(("datahistzmm_"+observable).c_str()));
    // Z->mumu connected with Z->nunu SR
    vector<pair<RooRealVar*, TH1*> >   znn_ZM_syst;
  
    // Other MC backgrounds in dimuon control region
    addTemplate("WJets_ZM_"+suffix, vars, *wspace_ZM, (TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_"+observable).c_str()));
    addTemplate("Top_ZM_"+suffix , vars, *wspace_ZM, (TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_"+observable).c_str()));
    addTemplate("QCD_ZM_"+suffix , vars, *wspace_ZM, (TH1F*)templatesfile->FindObjectAny(("qbkghistzmm_"+observable).c_str()));
    addTemplate("Dibosons_ZM_"+suffix, vars, *wspace_ZM, (TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_"+observable).c_str()));

    if(addShapeSystematics){
      
      TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_"+observable).c_str());
      TH1F* histobUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_bUp_"+observable).c_str());
      TH1F* histobDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_bDw_"+observable).c_str());
      TH1F* histoJesUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_metJetUp_"+observable).c_str());
      TH1F* histoJesDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_metJetDw_"+observable).c_str());
      TH1F* histoJerUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_metResUp_"+observable).c_str());
      TH1F* histoJerDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_metResDw_"+observable).c_str());
      TH1F* histoUncUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_metUncUp_"+observable).c_str());
      TH1F* histoUncDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzmm_metUncDw_"+observable).c_str());

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
      
      addTemplate("WJets_ZM_"+suffix+"_CMS_btagUp", vars, *wspace_ZM, histobUp);
      addTemplate("WJets_ZM_"+suffix+"_CMS_btagDown", vars, *wspace_ZM, histobDw);
      addTemplate("WJets_ZM_"+suffix+"_CMS_scale_jUp", vars, *wspace_ZM, histoJesUp);
      addTemplate("WJets_ZM_"+suffix+"_CMS_scale_jDown", vars, *wspace_ZM, histoJesDw);
      addTemplate("WJets_ZM_"+suffix+"_CMS_res_jUp", vars, *wspace_ZM, histoJerUp);
      addTemplate("WJets_ZM_"+suffix+"_CMS_res_jDown", vars, *wspace_ZM, histoJerDw);
      addTemplate("WJets_ZM_"+suffix+"_CMS_uncUp", vars, *wspace_ZM, histoUncUp);
      addTemplate("WJets_ZM_"+suffix+"_CMS_uncDown", vars, *wspace_ZM, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistzmm_metUncDw_"+observable).c_str());
      
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
      
      addTemplate("Top_ZM_"+suffix+"_CMS_btagUp", vars, *wspace_ZM, histobUp);
      addTemplate("Top_ZM_"+suffix+"_CMS_btagDown", vars, *wspace_ZM, histobDw);
      addTemplate("Top_ZM_"+suffix+"_CMS_scale_jUp", vars, *wspace_ZM, histoJesUp);
      addTemplate("Top_ZM_"+suffix+"_CMS_scale_jDown", vars, *wspace_ZM, histoJesDw);
      addTemplate("Top_ZM_"+suffix+"_CMS_res_jUp", vars, *wspace_ZM, histoJerUp);
      addTemplate("Top_ZM_"+suffix+"_CMS_res_jDown", vars, *wspace_ZM, histoJerDw);
      addTemplate("Top_ZM_"+suffix+"_CMS_uncUp", vars, *wspace_ZM, histoUncUp);
      addTemplate("Top_ZM_"+suffix+"_CMS_uncDown", vars, *wspace_ZM, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistzmm_metUncDw_"+observable).c_str());

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
      
      
      addTemplate("Dibosons_ZM_"+suffix+"_CMS_btagUp", vars, *wspace_ZM, histobUp);
      addTemplate("Dibosons_ZM_"+suffix+"_CMS_btagDown", vars, *wspace_ZM, histobDw);
      addTemplate("Dibosons_ZM_"+suffix+"_CMS_scale_jUp", vars, *wspace_ZM, histoJesUp);
      addTemplate("Dibosons_ZM_"+suffix+"_CMS_scale_jDown", vars, *wspace_ZM, histoJesDw);
      addTemplate("Dibosons_ZM_"+suffix+"_CMS_res_jUp", vars, *wspace_ZM, histoJerUp);
      addTemplate("Dibosons_ZM_"+suffix+"_CMS_res_jDown", vars, *wspace_ZM, histoJerDw);
      addTemplate("Dibosons_ZM_"+suffix+"_CMS_uncUp", vars, *wspace_ZM, histoUncUp);
      addTemplate("Dibosons_ZM_"+suffix+"_CMS_uncDown", vars, *wspace_ZM, histoUncDw);

      // temp Znunu shape systematics
      /*
      TFile* file_tmp = new TFile("monoV_hinv/postfit_weights_ZM.root");
      TH1F*   data     = (TH1F*) file_tmp->Get("data");
      TH1F*   prefit   = (TH1F*) file_tmp->Get("pre_fit_post_fit");
      data->Divide(prefit);
      TH1F*   zjets_pre_fit = (TH1F*) file_tmp->Get("zjets_pre_fit");
      TH1F*   wjets_pre_fit = (TH1F*) file_tmp->Get("wjets_pre_fit");
      TH1F*   top_pre_fit = (TH1F*) file_tmp->Get("top_pre_fit");
      TH1F*   dibos_pre_fit = (TH1F*) file_tmp->Get("diboson_pre_fit");
      
      TH1F*   zjets_pre_fit_shapeUp = (TH1F*) zjets_pre_fit->Clone("ZM_DMC");
      TH1F*   wjets_pre_fit_shapeUp = (TH1F*) wjets_pre_fit->Clone("wjets_pre_fit_shapeUp");
      TH1F*   wjets_pre_fit_shapeDw = (TH1F*) wjets_pre_fit->Clone("wjets_pre_fit_shapeDw");
      TH1F*   top_pre_fit_shapeUp = (TH1F*) top_pre_fit->Clone("top_pre_fit_shapeUp");
      TH1F*   top_pre_fit_shapeDw = (TH1F*) top_pre_fit->Clone("top_pre_fit_shapeDw");
      TH1F*   dibos_pre_fit_shapeUp = (TH1F*) dibos_pre_fit->Clone("diboson_pre_fit_shapeUp");
      TH1F*   dibos_pre_fit_shapeDw = (TH1F*) dibos_pre_fit->Clone("diboson_pre_fit_shapeDw");
      
      for(int iBin = 1 ; iBin <= zjets_pre_fit_shapeUp->GetNbinsX(); iBin++){
	if(data->GetBinContent(iBin) != 0)
	  zjets_pre_fit_shapeUp->SetBinContent(iBin,fabs(1-data->GetBinContent(iBin)));
	else 
	  zjets_pre_fit_shapeUp->SetBinContent(iBin,0.);
	wjets_pre_fit_shapeUp->SetBinContent(iBin,wjets_pre_fit->GetBinContent(iBin)*data->GetBinContent(iBin));
	top_pre_fit_shapeUp->SetBinContent(iBin,top_pre_fit->GetBinContent(iBin)*data->GetBinContent(iBin));
	dibos_pre_fit_shapeUp->SetBinContent(iBin,dibos_pre_fit->GetBinContent(iBin)*data->GetBinContent(iBin));
      }
      
      for(int iBin = 1 ; iBin <= wjets_pre_fit_shapeDw->GetNbinsX(); iBin++){
	if(data->GetBinContent(iBin) != 0){
	  top_pre_fit_shapeDw->SetBinContent(iBin,top_pre_fit->GetBinContent(iBin)*1./data->GetBinContent(iBin));
	  dibos_pre_fit_shapeDw->SetBinContent(iBin,dibos_pre_fit->GetBinContent(iBin)*1./data->GetBinContent(iBin));
	  wjets_pre_fit_shapeDw->SetBinContent(iBin,wjets_pre_fit->GetBinContent(iBin)*1./data->GetBinContent(iBin));
	}
	else{
	  wjets_pre_fit_shapeDw->SetBinContent(iBin,0);
	  top_pre_fit_shapeDw->SetBinContent(iBin,0);
	  dibos_pre_fit_shapeDw->SetBinContent(iBin,0);
	}
      }

      addTemplate("WJets_ZM_"+suffix+"_CMS_datamcUp", vars, *wspace_ZM, wjets_pre_fit_shapeUp);
      addTemplate("WJets_ZM_"+suffix+"_CMS_datamcDown", vars, *wspace_ZM, wjets_pre_fit_shapeDw);
      addTemplate("Top_ZM_"+suffix+"_CMS_datamcUp", vars, *wspace_ZM, top_pre_fit_shapeUp);
      addTemplate("Top_ZM_"+suffix+"_CMS_datamcDown", vars, *wspace_ZM, top_pre_fit_shapeDw);
      addTemplate("Dibosons_ZM_"+suffix+"_CMS_datamcUp", vars, *wspace_ZM, dibos_pre_fit_shapeUp);
      addTemplate("Dibosons_ZM_"+suffix+"_CMS_datamcDown", vars, *wspace_ZM, dibos_pre_fit_shapeDw);

      RooRealVar* znn_ZM_re1 = new RooRealVar("Znunu_ZM_Shape"  , "", 0., -5., 5.);
      znn_ZM_syst.push_back(pair<RooRealVar*, TH1*>(znn_ZM_re1 , zjets_pre_fit_shapeUp));      
      */
    }

    makeConnectedBinList("Znunu_ZM_"+suffix, met, *wspace_ZM, (TH1F*)templatesfile->FindObjectAny(("zmmcorhist_"+observable).c_str()), znn_ZM_syst, znn_SR_bins,NULL, observable);
 
    
    ////////////////////////////////////////
    // -------- CR Di-Electron  -------- //
    ///////////////////////////////////////
    cout<<"Make CR Di-Electron  templates ..."<<endl;
  
    
    wspace_ZE = new RooWorkspace(("ZE_"+suffix).c_str(),("ZE_"+suffix).c_str());

    addTemplate("data_obs_ZE_"+suffix, vars, *wspace_ZE, (TH1F*)templatesfile->FindObjectAny(("datahistzee_"+observable).c_str()));
    // Z->ee connected with Z->nunu SR
    vector<pair<RooRealVar*, TH1*> > znn_ZE_syst;
    makeConnectedBinList("Znunu_ZE_"+suffix, met, *wspace_ZE, (TH1F*)templatesfile->FindObjectAny(("zeecorhist_"+observable).c_str()), znn_ZE_syst, znn_SR_bins,NULL,observable);
    
    // Other MC backgrounds in dielectron control region
    addTemplate("WJets_ZE_"+suffix, vars, *wspace_ZE, (TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_"+observable).c_str()));
    addTemplate("Top_ZE_"+suffix  , vars, *wspace_ZE, (TH1F*)templatesfile->FindObjectAny(("tbkghistzee_"+observable).c_str()));
    addTemplate("QCD_ZE_"+suffix  , vars, *wspace_ZE, (TH1F*)templatesfile->FindObjectAny(("qbkghistzee_"+observable).c_str()));
    addTemplate("Dibosons_ZE_"+suffix  , vars, *wspace_ZE, (TH1F*)templatesfile->FindObjectAny(("dbkghistzee_"+observable).c_str()));

    if( addShapeSystematics){
      
      TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_"+observable).c_str());
      TH1F* histobUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_bUp_"+observable).c_str());
      TH1F* histobDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_bDw_"+observable).c_str());
      TH1F* histoJesUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_metJetUp_"+observable).c_str());
      TH1F* histoJesDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_metJetDw_"+observable).c_str());
      TH1F* histoJerUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_metResUp_"+observable).c_str());
      TH1F* histoJerDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_metResDw_"+observable).c_str());
      TH1F* histoUncUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_metUncUp_"+observable).c_str());
      TH1F* histoUncDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzee_metUncDw_"+observable).c_str());

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

      addTemplate("WJets_ZE_"+suffix+"_CMS_btagUp", vars, *wspace_ZE, histobUp);
      addTemplate("WJets_ZE_"+suffix+"_CMS_btagDown", vars, *wspace_ZE, histobDw);
      addTemplate("WJets_ZE_"+suffix+"_CMS_scale_jUp", vars, *wspace_ZE, histoJesUp);
      addTemplate("WJets_ZE_"+suffix+"_CMS_scale_jDown", vars, *wspace_ZE, histoJesDw);
      addTemplate("WJets_ZE_"+suffix+"_CMS_res_jUp", vars, *wspace_ZE, histoJerUp);
      addTemplate("WJets_ZE_"+suffix+"_CMS_res_jDown", vars, *wspace_ZE, histoJerDw);
      addTemplate("WJets_ZE_"+suffix+"_CMS_uncUp", vars, *wspace_ZE, histoUncUp);
      addTemplate("WJets_ZE_"+suffix+"_CMS_uncDown", vars, *wspace_ZE, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("tbkghistzee_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistzee_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistzee_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistzee_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistzee_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistzee_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistzee_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistzee_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistzee_metUncDw_"+observable).c_str());

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
      
      addTemplate("Top_ZE_"+suffix+"_CMS_btagUp", vars, *wspace_ZE, histobUp);
      addTemplate("Top_ZE_"+suffix+"_CMS_btagDown", vars, *wspace_ZE, histobDw);
      addTemplate("Top_ZE_"+suffix+"_CMS_scale_jUp", vars, *wspace_ZE, histoJesUp);
      addTemplate("Top_ZE_"+suffix+"_CMS_scale_jDown", vars, *wspace_ZE, histoJesDw);
      addTemplate("Top_ZE_"+suffix+"_CMS_res_jUp", vars, *wspace_ZE, histoJerUp);
      addTemplate("Top_ZE_"+suffix+"_CMS_res_jDown", vars, *wspace_ZE, histoJerDw);
      addTemplate("Top_ZE_"+suffix+"_CMS_uncUp", vars, *wspace_ZE, histoUncUp);
      addTemplate("Top_ZE_"+suffix+"_CMS_uncDown", vars, *wspace_ZE, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("dbkghistzee_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistzee_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistzee_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistzee_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistzee_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistzee_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistzee_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistzee_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistzee_metUncDw_"+observable).c_str());
      
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
      
      
      addTemplate("Dibosons_ZE_"+suffix+"_CMS_btagUp", vars, *wspace_ZE, histobUp);
      addTemplate("Dibosons_ZE_"+suffix+"_CMS_btagDown", vars, *wspace_ZE, histobDw);
      addTemplate("Dibosons_ZE_"+suffix+"_CMS_scale_jUp", vars, *wspace_ZE, histoJesUp);
      addTemplate("Dibosons_ZE_"+suffix+"_CMS_scale_jDown", vars, *wspace_ZE, histoJesDw);
      addTemplate("Dibosons_ZE_"+suffix+"_CMS_res_jUp", vars, *wspace_ZE, histoJerUp);
      addTemplate("Dibosons_ZE_"+suffix+"_CMS_res_jDown", vars, *wspace_ZE, histoJerDw);
      addTemplate("Dibosons_ZE_"+suffix+"_CMS_uncUp", vars, *wspace_ZE, histoUncUp);
      addTemplate("Dibosons_ZE_"+suffix+"_CMS_uncDown", vars, *wspace_ZE, histoUncDw);          
    }
  }
  else{
    
    cout<<"Make CR Di-Lepton  templates ..."<<endl;    
    wspace_ZL = new RooWorkspace(("ZL_"+suffix).c_str(),("ZL_"+suffix).c_str());
  
    addTemplate("data_obs_ZL_"+suffix, vars, *wspace_ZL, (TH1F*)templatesfile->FindObjectAny(("datahistzll_"+observable).c_str()));
    // Z->mumu connected with Z->nunu SR
    vector<pair<RooRealVar*, TH1*> >   znn_syst;
    makeConnectedBinList("Znunu_ZL_"+suffix, met, *wspace_ZL, (TH1F*)templatesfile->FindObjectAny(("zllcorhist_"+observable).c_str()), znn_syst, znn_SR_bins,NULL, observable);
  
    // Other MC backgrounds in dimuon control region
    addTemplate("WJets_ZL_"+suffix, vars, *wspace_ZL, (TH1F*)templatesfile->FindObjectAny(("vlbkghistzll_"+observable).c_str()));
    addTemplate("Top_ZL_"+suffix , vars, *wspace_ZL, (TH1F*)templatesfile->FindObjectAny(("tbkghistzll_"+observable).c_str()));
    addTemplate("QCD_ZL_"+suffix , vars, *wspace_ZL, (TH1F*)templatesfile->FindObjectAny(("qbkghistzll_"+observable).c_str()));
    addTemplate("Dibosons_ZL_"+suffix, vars, *wspace_ZL, (TH1F*)templatesfile->FindObjectAny(("dbkghistzll_"+observable).c_str()));

    if(addShapeSystematics){
      
      TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzll_"+observable).c_str());
      TH1F* histobUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzll_bUp_"+observable).c_str());
      TH1F* histobDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzll_bDw_"+observable).c_str());
      TH1F* histoJesUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzll_metJetUp_"+observable).c_str());
      TH1F* histoJesDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzll_metJetDw_"+observable).c_str());
      TH1F* histoJerUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzll_metResUp_"+observable).c_str());
      TH1F* histoJerDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzll_metResDw_"+observable).c_str());
      TH1F* histoUncUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzll_metUncUp_"+observable).c_str());
      TH1F* histoUncDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghistzll_metUncDw_"+observable).c_str());

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
      
      addTemplate("WJets_ZL_"+suffix+"_CMS_btagUp", vars, *wspace_ZL, histobUp);
      addTemplate("WJets_ZL_"+suffix+"_CMS_btagDown", vars, *wspace_ZL, histobDw);
      addTemplate("WJets_ZL_"+suffix+"_CMS_scale_jUp", vars, *wspace_ZL, histoJesUp);
      addTemplate("WJets_ZL_"+suffix+"_CMS_scale_jDown", vars, *wspace_ZL, histoJesDw);
      addTemplate("WJets_ZL_"+suffix+"_CMS_res_jUp", vars, *wspace_ZL, histoJerUp);
      addTemplate("WJets_ZL_"+suffix+"_CMS_res_jDown", vars, *wspace_ZL, histoJerDw);
      addTemplate("WJets_ZL_"+suffix+"_CMS_uncUp", vars, *wspace_ZL, histoUncUp);
      addTemplate("WJets_ZL_"+suffix+"_CMS_uncDown", vars, *wspace_ZL, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("tbkghistzll_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistzll_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistzll_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistzll_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistzll_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistzll_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistzll_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistzll_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistzll_metUncDw_"+observable).c_str());
      
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
      
      addTemplate("Top_ZL_"+suffix+"_CMS_btagUp", vars, *wspace_ZL, histobUp);
      addTemplate("Top_ZL_"+suffix+"_CMS_btagDown", vars, *wspace_ZL, histobDw);
      addTemplate("Top_ZL_"+suffix+"_CMS_scale_jUp", vars, *wspace_ZL, histoJesUp);
      addTemplate("Top_ZL_"+suffix+"_CMS_scale_jDown", vars, *wspace_ZL, histoJesDw);
      addTemplate("Top_ZL_"+suffix+"_CMS_res_jUp", vars, *wspace_ZL, histoJerUp);
      addTemplate("Top_ZL_"+suffix+"_CMS_res_jDown", vars, *wspace_ZL, histoJerDw);
      addTemplate("Top_ZL_"+suffix+"_CMS_uncUp", vars, *wspace_ZL, histoUncUp);
      addTemplate("Top_ZL_"+suffix+"_CMS_uncDown", vars, *wspace_ZL, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("dbkghistzll_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistzll_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistzll_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistzll_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistzll_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistzll_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistzll_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistzll_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistzll_metUncDw_"+observable).c_str());

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
      
      
      addTemplate("Dibosons_ZL_"+suffix+"_CMS_btagUp", vars, *wspace_ZL, histobUp);
      addTemplate("Dibosons_ZL_"+suffix+"_CMS_btagDown", vars, *wspace_ZL, histobDw);
      addTemplate("Dibosons_ZL_"+suffix+"_CMS_scale_jUp", vars, *wspace_ZL, histoJesUp);
      addTemplate("Dibosons_ZL_"+suffix+"_CMS_scale_jDown", vars, *wspace_ZL, histoJesDw);
      addTemplate("Dibosons_ZL_"+suffix+"_CMS_res_jUp", vars, *wspace_ZL, histoJerUp);
      addTemplate("Dibosons_ZL_"+suffix+"_CMS_res_jDown", vars, *wspace_ZL, histoJerDw);
      addTemplate("Dibosons_ZL_"+suffix+"_CMS_uncUp", vars, *wspace_ZL, histoUncUp);
      addTemplate("Dibosons_ZL_"+suffix+"_CMS_uncDown", vars, *wspace_ZL, histoUncDw);
      
    }            
  }
  ///////////////////////////////////////
  // -------- CR Gamma+jets  -------- //
  //////////////////////////////////////
  cout<<"Make CR Gamma+jets  templates ..."<<endl;
  RooWorkspace wspace_GJ(("GJ_"+suffix).c_str(),("GJ_"+suffix).c_str());

  addTemplate("data_obs_GJ_"+suffix, vars, wspace_GJ, (TH1F*)templatesfile->FindObjectAny(("datahistgam_"+observable).c_str()));
  // Gamma+jets --> connected with Z->nunu
  vector<pair<RooRealVar*, TH1*> > znn_GJ_syst;
  RooRealVar* znn_GJ_re1 = new RooRealVar("Znunu_GJ_RenScale1"  , "", 0., -5., 5.);
  RooRealVar* znn_GJ_fa1 = new RooRealVar("Znunu_GJ_FactScale1" , "", 0., -5., 5.);
  RooRealVar* znn_GJ_re2 = new RooRealVar("Znunu_GJ_RenScale2"  , "", 0., -5., 5.);
  RooRealVar* znn_GJ_fa2 = new RooRealVar("Znunu_GJ_FactScale2" , "", 0., -5., 5.);
  RooRealVar* znn_GJ_pdf = new RooRealVar("Znunu_GJ_PDF"        , "", 0., -5., 5.);
  RooRealVar* znn_GJ_fpc = new RooRealVar("Znunu_GJ_Footprint"  , "", 0., -5., 5.);
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(NULL      , (TH1F*)templatesfile->FindObjectAny(("ZG_EWK_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_re1, (TH1F*)templatesfile->FindObjectAny(("ZG_RenScale1_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_fa1, (TH1F*)templatesfile->FindObjectAny(("ZG_FactScale1_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_re2, (TH1F*)templatesfile->FindObjectAny(("ZG_RenScale2_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_fa2, (TH1F*)templatesfile->FindObjectAny(("ZG_FactScale2_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_pdf, (TH1F*)templatesfile->FindObjectAny(("ZG_PDF_"+observable).c_str())));
  znn_GJ_syst.push_back(pair<RooRealVar*, TH1*>(znn_GJ_fpc, (TH1F*)templatesfile->FindObjectAny(("ZG_Footprint_"+observable).c_str())));

  makeConnectedBinList("Znunu_GJ_"+suffix, met, wspace_GJ, (TH1F*)templatesfile->FindObjectAny(("gamcorewkhist_"+observable).c_str()), znn_GJ_syst, znn_SR_bins,NULL,observable);  
  // Other MC backgrounds photon+jets control region
  addTemplate("QCD_GJ_"+suffix, vars, wspace_GJ, (TH1F*)templatesfile->FindObjectAny(("qbkghistgam_"+observable).c_str()));

  
  if(not mergeLeptons){
    ///////////////////////////////////////
    // -------- CR Single-Muon  -------- //
    //////////////////////////////////////
    cout<<"Make CR Single-Mu  templates ..."<<endl;
    
    wspace_WM = new RooWorkspace(("WM_"+suffix).c_str(),("WM_"+suffix).c_str());

    addTemplate("data_obs_WM_"+suffix, vars, *wspace_WM, (TH1F*)templatesfile->FindObjectAny(("datahistwmn_"+observable).c_str()));
    // connected W->munu with W+jets SR
    vector<pair<RooRealVar*, TH1*> > wln_WM_syst;
    makeConnectedBinList("WJets_WM_"+suffix, met, *wspace_WM, (TH1F*)templatesfile->FindObjectAny(("wmncorhist_"+observable).c_str()), wln_WM_syst, wln_SR_bins,NULL,observable);
    
  // Other MC backgrounds in single muon control region
    addTemplate("ZJets_WM_"+suffix, vars, *wspace_WM, (TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_"+observable).c_str()));
    addTemplate("Top_WM_"+suffix  , vars, *wspace_WM, (TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_"+observable).c_str()));
    addTemplate("QCD_WM_"+suffix  , vars, *wspace_WM, (TH1F*)templatesfile->FindObjectAny(("qbkghistwmn_"+observable).c_str()));
    addTemplate("Dibosons_WM_"+suffix, vars, *wspace_WM, (TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_"+observable).c_str()));
    
    if(addShapeSystematics){
      
      TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_"+observable).c_str());
      TH1F* histobUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_bUp_"+observable).c_str());
      TH1F* histobDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_bDw_"+observable).c_str());
      TH1F* histoJesUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_metJetUp_"+observable).c_str());
      TH1F* histoJesDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_metJetDw_"+observable).c_str());
      TH1F* histoJerUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_metResUp_"+observable).c_str());
      TH1F* histoJerDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_metResDw_"+observable).c_str());
      TH1F* histoUncUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_metUncUp_"+observable).c_str());
      TH1F* histoUncDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwmn_metUncDw_"+observable).c_str());
      
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
      
      addTemplate("ZJets_WM_"+suffix+"_CMS_btagUp", vars, *wspace_WM, histobUp);
      addTemplate("ZJets_WM_"+suffix+"_CMS_btagDown", vars, *wspace_WM, histobDw);
      addTemplate("ZJets_WM_"+suffix+"_CMS_scale_jUp", vars, *wspace_WM, histoJesUp);
      addTemplate("ZJets_WM_"+suffix+"_CMS_scale_jDown", vars, *wspace_WM, histoJesDw);
      addTemplate("ZJets_WM_"+suffix+"_CMS_res_jUp", vars, *wspace_WM, histoJerUp);
      addTemplate("ZJets_WM_"+suffix+"_CMS_res_jDown", vars, *wspace_WM, histoJerDw);
      addTemplate("ZJets_WM_"+suffix+"_CMS_uncUp", vars, *wspace_WM, histoUncUp);
      addTemplate("ZJets_WM_"+suffix+"_CMS_uncDown", vars, *wspace_WM, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistwmn_metUncDw_"+observable).c_str());
      
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
      
      addTemplate("Top_WM_"+suffix+"_CMS_btagUp", vars, *wspace_WM, histobUp);
      addTemplate("Top_WM_"+suffix+"_CMS_btagDown", vars, *wspace_WM, histobDw);
      addTemplate("Top_WM_"+suffix+"_CMS_scale_jUp", vars, *wspace_WM, histoJesUp);
      addTemplate("Top_WM_"+suffix+"_CMS_scale_jDown", vars, *wspace_WM, histoJesDw);
      addTemplate("Top_WM_"+suffix+"_CMS_res_jUp", vars, *wspace_WM, histoJerUp);
      addTemplate("Top_WM_"+suffix+"_CMS_res_jDown", vars, *wspace_WM, histoJerDw);
      addTemplate("Top_WM_"+suffix+"_CMS_uncUp", vars, *wspace_WM, histoUncUp);
      addTemplate("Top_WM_"+suffix+"_CMS_uncDown", vars, *wspace_WM, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistwmn_metUncDw_"+observable).c_str());

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
      
      
      addTemplate("Dibosons_WM_"+suffix+"_CMS_btagUp", vars, *wspace_WM, histobUp);
      addTemplate("Dibosons_WM_"+suffix+"_CMS_btagDown", vars, *wspace_WM, histobDw);
      addTemplate("Dibosons_WM_"+suffix+"_CMS_scale_jUp", vars, *wspace_WM, histoJesUp);
      addTemplate("Dibosons_WM_"+suffix+"_CMS_scale_jDown", vars, *wspace_WM, histoJesDw);
      addTemplate("Dibosons_WM_"+suffix+"_CMS_res_jUp", vars, *wspace_WM, histoJerUp);
      addTemplate("Dibosons_WM_"+suffix+"_CMS_res_jDown", vars, *wspace_WM, histoJerDw);
      addTemplate("Dibosons_WM_"+suffix+"_CMS_uncUp", vars, *wspace_WM, histoUncUp);
      addTemplate("Dibosons_WM_"+suffix+"_CMS_uncDown", vars, *wspace_WM, histoUncDw);      
    }
  
    //////////////////////////////////..../////
    // -------- CR Single-Electron  -------- //
    //////////////////////////////////////////
    cout<<"Make CR Single-El  templates ..."<<endl;
    
    wspace_WE = new RooWorkspace(("WE_"+suffix).c_str(),("WE_"+suffix).c_str());
    
    addTemplate("data_obs_WE_"+suffix, vars, *wspace_WE, (TH1F*)templatesfile->FindObjectAny(("datahistwen_"+observable).c_str()));
    // connected W->enu with W+jets SR 
    vector<pair<RooRealVar*, TH1*> > wln_WE_syst;
    makeConnectedBinList("WJets_WE_"+suffix, met, *wspace_WE, (TH1F*)templatesfile->FindObjectAny(("wencorhist_"+observable).c_str()), wln_WE_syst, wln_SR_bins,NULL,observable);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_WE_"+suffix, vars, *wspace_WE, (TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_"+observable).c_str()));
    addTemplate("Top_WE_"+suffix, vars, *wspace_WE, (TH1F*)templatesfile->FindObjectAny(("tbkghistwen_"+observable).c_str()));
    addTemplate("QCD_WE_"+suffix, vars, *wspace_WE, (TH1F*)templatesfile->FindObjectAny(("qbkghistwen_"+observable).c_str()));
    addTemplate("Dibosons_WE_"+suffix, vars, *wspace_WE, (TH1F*)templatesfile->FindObjectAny(("dbkghistwen_"+observable).c_str()));
    
    if(addShapeSystematics){
      
      TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_"+observable).c_str());
      TH1F* histobUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_bUp_"+observable).c_str());
      TH1F* histobDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_bDw_"+observable).c_str());
      TH1F* histoJesUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_metJetUp_"+observable).c_str());
      TH1F* histoJesDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_metJetDw_"+observable).c_str());
      TH1F* histoJerUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_metResUp_"+observable).c_str());
      TH1F* histoJerDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_metResDw_"+observable).c_str());
      TH1F* histoUncUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_metUncUp_"+observable).c_str());
      TH1F* histoUncDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwen_metUncDw_"+observable).c_str());
      
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

      addTemplate("ZJets_WE_"+suffix+"_CMS_btagUp", vars, *wspace_WE, histobUp);
      addTemplate("ZJets_WE_"+suffix+"_CMS_btagDown", vars, *wspace_WE, histobDw);
      addTemplate("ZJets_WE_"+suffix+"_CMS_scale_jUp", vars, *wspace_WE, histoJesUp);
      addTemplate("ZJets_WE_"+suffix+"_CMS_scale_jDown", vars, *wspace_WE, histoJesDw);
      addTemplate("ZJets_WE_"+suffix+"_CMS_res_jUp", vars, *wspace_WE, histoJerUp);
      addTemplate("ZJets_WE_"+suffix+"_CMS_res_jDown", vars, *wspace_WE, histoJerDw);
      addTemplate("ZJets_WE_"+suffix+"_CMS_uncUp", vars, *wspace_WE, histoUncUp);
      addTemplate("ZJets_WE_"+suffix+"_CMS_uncDown", vars, *wspace_WE, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("tbkghistwen_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistwen_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistwen_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistwen_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistwen_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistwen_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistwen_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistwen_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistwen_metUncDw_"+observable).c_str());
      
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
      
      addTemplate("Top_WE_"+suffix+"_CMS_btagUp", vars, *wspace_WE, histobUp);
      addTemplate("Top_WE_"+suffix+"_CMS_btagDown", vars, *wspace_WE, histobDw);
      addTemplate("Top_WE_"+suffix+"_CMS_scale_jUp", vars, *wspace_WE, histoJesUp);
      addTemplate("Top_WE_"+suffix+"_CMS_scale_jDown", vars, *wspace_WE, histoJesDw);
      addTemplate("Top_WE_"+suffix+"_CMS_res_jUp", vars, *wspace_WE, histoJerUp);
      addTemplate("Top_WE_"+suffix+"_CMS_res_jDown", vars, *wspace_WE, histoJerDw);
      addTemplate("Top_WE_"+suffix+"_CMS_uncUp", vars, *wspace_WE, histoUncUp);
      addTemplate("Top_WE_"+suffix+"_CMS_uncDown", vars, *wspace_WE, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("dbkghistwen_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistwen_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistwen_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistwen_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistwen_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistwen_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistwen_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistwen_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistwen_metUncDw_"+observable).c_str());
      
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
      
      addTemplate("Dibosons_WE_"+suffix+"_CMS_btagUp", vars, *wspace_WE, histobUp);
      addTemplate("Dibosons_WE_"+suffix+"_CMS_btagDown", vars, *wspace_WE, histobDw);
      addTemplate("Dibosons_WE_"+suffix+"_CMS_scale_jUp", vars, *wspace_WE, histoJesUp);
      addTemplate("Dibosons_WE_"+suffix+"_CMS_scale_jDown", vars, *wspace_WE, histoJesDw);
      addTemplate("Dibosons_WE_"+suffix+"_CMS_res_jUp", vars, *wspace_WE, histoJerUp);
      addTemplate("Dibosons_WE_"+suffix+"_CMS_res_jDown", vars, *wspace_WE, histoJerDw);
      addTemplate("Dibosons_WE_"+suffix+"_CMS_uncUp", vars, *wspace_WE, histoUncUp);
      addTemplate("Dibosons_WE_"+suffix+"_CMS_uncDown", vars, *wspace_WE, histoUncDw);

    }
  }
  else{

    cout<<"Make CR Single-Lepton templates ..."<<endl;
    
    wspace_WL = new RooWorkspace(("WL_"+suffix).c_str(),("WL_"+suffix).c_str());
    
    addTemplate("data_obs_WL_"+suffix, vars, *wspace_WL, (TH1F*)templatesfile->FindObjectAny(("datahistwln_"+observable).c_str()));
    // connected W->enu with W+jets SR 
    vector<pair<RooRealVar*, TH1*> > wln_WL_syst;
    makeConnectedBinList("WJets_WL_"+suffix, met, *wspace_WL, (TH1F*)templatesfile->FindObjectAny(("wlncorhist_"+observable).c_str()), wln_WL_syst, wln_SR_bins,NULL,observable);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_WL_"+suffix, vars, *wspace_WL, (TH1F*)templatesfile->FindObjectAny(("vllbkghistwln_"+observable).c_str()));
    addTemplate("Top_WL_"+suffix, vars, *wspace_WL, (TH1F*)templatesfile->FindObjectAny(("tbkghistwln_"+observable).c_str()));
    addTemplate("QCD_WL_"+suffix, vars, *wspace_WL, (TH1F*)templatesfile->FindObjectAny(("qbkghistwln_"+observable).c_str()));
    addTemplate("Dibosons_WL_"+suffix, vars, *wspace_WL, (TH1F*)templatesfile->FindObjectAny(("dbkghistwln_"+observable).c_str()));
    
    if(addShapeSystematics){
      
      TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwln_"+observable).c_str());
      TH1F* histobUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwln_bUp_"+observable).c_str());
      TH1F* histobDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwln_bDw_"+observable).c_str());
      TH1F* histoJesUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwln_metJetUp_"+observable).c_str());
      TH1F* histoJesDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwln_metJetDw_"+observable).c_str());
      TH1F* histoJerUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwln_metResUp_"+observable).c_str());
      TH1F* histoJerDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwln_metResDw_"+observable).c_str());
      TH1F* histoUncUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwln_metUncUp_"+observable).c_str());
      TH1F* histoUncDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghistwln_metUncDw_"+observable).c_str());
      
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

      addTemplate("ZJets_WL_"+suffix+"_CMS_btagUp", vars, *wspace_WL, histobUp);
      addTemplate("ZJets_WL_"+suffix+"_CMS_btagDown", vars, *wspace_WL, histobDw);
      addTemplate("ZJets_WL_"+suffix+"_CMS_scale_jUp", vars, *wspace_WL, histoJesUp);
      addTemplate("ZJets_WL_"+suffix+"_CMS_scale_jDown", vars, *wspace_WL, histoJesDw);
      addTemplate("ZJets_WL_"+suffix+"_CMS_res_jUp", vars, *wspace_WL, histoJerUp);
      addTemplate("ZJets_WL_"+suffix+"_CMS_res_jDown", vars, *wspace_WL, histoJerDw);
      addTemplate("ZJets_WL_"+suffix+"_CMS_uncUp", vars, *wspace_WL, histoUncUp);
      addTemplate("ZJets_WL_"+suffix+"_CMS_uncDown", vars, *wspace_WL, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("tbkghistwln_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistwln_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistwln_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistwln_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistwln_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistwln_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistwln_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("tbkghistwln_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("tbkghistwln_metUncDw_"+observable).c_str());
      
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
      
      addTemplate("Top_WL_"+suffix+"_CMS_btagUp", vars, *wspace_WL, histobUp);
      addTemplate("Top_WL_"+suffix+"_CMS_btagDown", vars, *wspace_WL, histobDw);
      addTemplate("Top_WL_"+suffix+"_CMS_scale_jUp", vars, *wspace_WL, histoJesUp);
      addTemplate("Top_WL_"+suffix+"_CMS_scale_jDown", vars, *wspace_WL, histoJesDw);
      addTemplate("Top_WL_"+suffix+"_CMS_res_jUp", vars, *wspace_WL, histoJerUp);
      addTemplate("Top_WL_"+suffix+"_CMS_res_jDown", vars, *wspace_WL, histoJerDw);
      addTemplate("Top_WL_"+suffix+"_CMS_uncUp", vars, *wspace_WL, histoUncUp);
      addTemplate("Top_WL_"+suffix+"_CMS_uncDown", vars, *wspace_WL, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("dbkghistwln_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistwln_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistwln_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistwln_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistwln_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistwln_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistwln_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("dbkghistwln_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("dbkghistwln_metUncDw_"+observable).c_str());
      
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
      
      addTemplate("Dibosons_WL_"+suffix+"_CMS_btagUp", vars, *wspace_WL, histobUp);
      addTemplate("Dibosons_WL_"+suffix+"_CMS_btagDown", vars, *wspace_WL, histobDw);
      addTemplate("Dibosons_WL_"+suffix+"_CMS_scale_jUp", vars, *wspace_WL, histoJesUp);
      addTemplate("Dibosons_WL_"+suffix+"_CMS_scale_jDown", vars, *wspace_WL, histoJesDw);
      addTemplate("Dibosons_WL_"+suffix+"_CMS_res_jUp", vars, *wspace_WL, histoJerUp);
      addTemplate("Dibosons_WL_"+suffix+"_CMS_res_jDown", vars, *wspace_WL, histoJerDw);
      addTemplate("Dibosons_WL_"+suffix+"_CMS_uncUp", vars, *wspace_WL, histoUncUp);
      addTemplate("Dibosons_WL_"+suffix+"_CMS_uncDown", vars, *wspace_WL, histoUncDw);

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

    addTemplate("data_obs_TM_"+suffix, vars, *wspace_TM, (TH1F*)templatesfile->FindObjectAny(("datahisttopmu_"+observable).c_str()));
    // connect tt->mu+b with tt SR
    vector<pair<RooRealVar*, TH1*> > top_TM_syst;
    RooRealVar* top_btag   = new RooRealVar("Top_btag", "", 0., -5., 5.);
    top_TM_syst.push_back(pair<RooRealVar*, TH1*>(top_btag, (TH1F*)templatesfile->FindObjectAny(("TOP_MU_B_"+observable).c_str())));
    
    makeConnectedBinList("Top_TM_"+suffix, met, *wspace_TM, (TH1F*)templatesfile->FindObjectAny(("topmucorhist_"+observable).c_str()), top_TM_syst, top_SR_bins,NULL,observable);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_TM_"+suffix, vars, *wspace_TM, (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopmu_"+observable).c_str()));
    addTemplate("WJets_TM_"+suffix, vars, *wspace_TM, (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopmu_"+observable).c_str()));
    addTemplate("QCD_TM_"+suffix, vars, *wspace_TM, (TH1F*)templatesfile->FindObjectAny(("qbkghisttopmu_"+observable).c_str()));
    addTemplate("Dibosons_TM_"+suffix, vars, *wspace_TM, (TH1F*)templatesfile->FindObjectAny(("dbkghisttopmu_"+observable).c_str()));

    if(addShapeSystematics){
      
      TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopmu_"+observable).c_str());
      TH1F* histobUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopmu_bUp_"+observable).c_str());
      TH1F* histobDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopmu_bDw_"+observable).c_str());
      TH1F* histoJesUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopmu_metJetUp_"+observable).c_str());
      TH1F* histoJesDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopmu_metJetDw_"+observable).c_str());
      TH1F* histoJerUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopmu_metResUp_"+observable).c_str());
      TH1F* histoJerDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopmu_metResDw_"+observable).c_str());
      TH1F* histoUncUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopmu_metUncUp_"+observable).c_str());
      TH1F* histoUncDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopmu_metUncDw_"+observable).c_str());
      
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
      
      addTemplate("ZJets_TM_"+suffix+"_CMS_btagUp", vars, *wspace_TM, histobUp);
      addTemplate("ZJets_TM_"+suffix+"_CMS_btagDown", vars, *wspace_TM, histobDw);
      addTemplate("ZJets_TM_"+suffix+"_CMS_scale_jUp", vars, *wspace_TM, histoJesUp);
      addTemplate("ZJets_TM_"+suffix+"_CMS_scale_jDown", vars, *wspace_TM, histoJesDw);
      addTemplate("ZJets_TM_"+suffix+"_CMS_res_jUp", vars, *wspace_TM, histoJerUp);
      addTemplate("ZJets_TM_"+suffix+"_CMS_res_jDown", vars, *wspace_TM, histoJerDw);
      addTemplate("ZJets_TM_"+suffix+"_CMS_uncUp", vars, *wspace_TM, histoUncUp);
      addTemplate("ZJets_TM_"+suffix+"_CMS_uncDown", vars, *wspace_TM, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopmu_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopmu_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopmu_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopmu_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopmu_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopmu_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopmu_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopmu_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopmu_metUncDw_"+observable).c_str());

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
      
      addTemplate("WJets_TM_"+suffix+"_CMS_btagUp", vars, *wspace_TM, histobUp);
      addTemplate("WJets_TM_"+suffix+"_CMS_btagDown", vars, *wspace_TM, histobDw);
      addTemplate("WJets_TM_"+suffix+"_CMS_scale_jUp", vars, *wspace_TM, histoJesUp);
      addTemplate("WJets_TM_"+suffix+"_CMS_scale_jDown", vars, *wspace_TM, histoJesDw);
      addTemplate("WJets_TM_"+suffix+"_CMS_res_jUp", vars, *wspace_TM, histoJerUp);
      addTemplate("WJets_TM_"+suffix+"_CMS_res_jDown", vars, *wspace_TM, histoJerDw);
      addTemplate("WJets_TM_"+suffix+"_CMS_uncUp", vars, *wspace_TM, histoUncUp);
      addTemplate("WJets_TM_"+suffix+"_CMS_uncDown", vars, *wspace_TM, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopmu_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopmu_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopmu_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopmu_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopmu_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopmu_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopmu_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopmu_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopmu_metUncDw_"+observable).c_str());
      
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
      
      
      addTemplate("Dibosons_TM_"+suffix+"_CMS_btagUp", vars, *wspace_TM, histobUp);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_btagDown", vars, *wspace_TM, histobDw);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_scale_jUp", vars, *wspace_TM, histoJesUp);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_scale_jDown", vars, *wspace_TM, histoJesDw);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_res_jUp", vars, *wspace_TM, histoJerUp);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_res_jDown", vars, *wspace_TM, histoJerDw);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_uncUp", vars, *wspace_TM, histoUncUp);
      addTemplate("Dibosons_TM_"+suffix+"_CMS_uncDown", vars, *wspace_TM, histoUncDw);
           
    }

    /////////////////////////////////////
    // -------- CR Top-Electron-------- //
    ////////////////////////////////////
    cout<<"Make CR Top-el  templates ..."<<endl;
    wspace_TE = new RooWorkspace(("TE_"+suffix).c_str(),("TE_"+suffix).c_str());

    addTemplate("data_obs_TE_"+suffix, vars, *wspace_TE, (TH1F*)templatesfile->FindObjectAny(("datahisttopel_"+observable).c_str()));
    // connect tt->mu+b with tt SR
    vector<pair<RooRealVar*, TH1*> > top_TE_syst;
    top_TE_syst.push_back(pair<RooRealVar*, TH1*>(top_btag, (TH1F*)templatesfile->FindObjectAny(("TOP_EL_B_"+observable).c_str())));
    
    makeConnectedBinList("Top_TE_"+suffix, met, *wspace_TE, (TH1F*)templatesfile->FindObjectAny(("topelcorhist_"+observable).c_str()), top_TE_syst, top_SR_bins,NULL,observable);
    
    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_TE_"+suffix , vars, *wspace_TE, (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopel_"+observable).c_str()));
    addTemplate("WJets_TE_"+suffix , vars, *wspace_TE, (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopel_"+observable).c_str()));
    addTemplate("QCD_TE_"+suffix   , vars, *wspace_TE, (TH1F*)templatesfile->FindObjectAny(("qbkghisttopel_"+observable).c_str()));
    addTemplate("Dibosons_TE_"+suffix, vars, *wspace_TE, (TH1F*)templatesfile->FindObjectAny(("dbkghisttopel_"+observable).c_str()));

    if(addShapeSystematics){
      
      TH1F* nominalHisto = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopel_"+observable).c_str());
      TH1F* histobUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopel_bUp_"+observable).c_str());
      TH1F* histobDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopel_bDw_"+observable).c_str());
      TH1F* histoJesUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopel_metJetUp_"+observable).c_str());
      TH1F* histoJesDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopel_metJetDw_"+observable).c_str());
      TH1F* histoJerUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopel_metResUp_"+observable).c_str());
      TH1F* histoJerDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopel_metResDw_"+observable).c_str());
      TH1F* histoUncUp = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopel_metUncUp_"+observable).c_str());
      TH1F* histoUncDw = (TH1F*)templatesfile->FindObjectAny(("vllbkghisttopel_metUncDw_"+observable).c_str());
      
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
      
      addTemplate("ZJets_TE_"+suffix+"_CMS_btagUp", vars, *wspace_TE, histobUp);
      addTemplate("ZJets_TE_"+suffix+"_CMS_btagDown", vars, *wspace_TE, histobDw);
      addTemplate("ZJets_TE_"+suffix+"_CMS_scale_jUp", vars, *wspace_TE, histoJesUp);
      addTemplate("ZJets_TE_"+suffix+"_CMS_scale_jDown", vars, *wspace_TE, histoJesDw);
      addTemplate("ZJets_TE_"+suffix+"_CMS_res_jUp", vars, *wspace_TE, histoJerUp);
      addTemplate("ZJets_TE_"+suffix+"_CMS_res_jDown", vars, *wspace_TE, histoJerDw);
      addTemplate("ZJets_TE_"+suffix+"_CMS_uncUp", vars, *wspace_TE, histoUncUp);
      addTemplate("ZJets_TE_"+suffix+"_CMS_uncDown", vars, *wspace_TE, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopel_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopel_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopel_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopel_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopel_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopel_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopel_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopel_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("vlbkghisttopel_metUncDw_"+observable).c_str());

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
      
      addTemplate("WJets_TE_"+suffix+"_CMS_btagUp", vars, *wspace_TE, histobUp);
      addTemplate("WJets_TE_"+suffix+"_CMS_btagDown", vars, *wspace_TE, histobDw);
      addTemplate("WJets_TE_"+suffix+"_CMS_scale_jUp", vars, *wspace_TE, histoJesUp);
      addTemplate("WJets_TE_"+suffix+"_CMS_scale_jDown", vars, *wspace_TE, histoJesDw);
      addTemplate("WJets_TE_"+suffix+"_CMS_res_jUp", vars, *wspace_TE, histoJerUp);
      addTemplate("WJets_TE_"+suffix+"_CMS_res_jDown", vars, *wspace_TE, histoJerDw);
      addTemplate("WJets_TE_"+suffix+"_CMS_uncUp", vars, *wspace_TE, histoUncUp);
      addTemplate("WJets_TE_"+suffix+"_CMS_uncDown", vars, *wspace_TE, histoUncDw);
      
      ///
      nominalHisto = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopel_"+observable).c_str());
      histobUp = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopel_bUp_"+observable).c_str());
      histobDw = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopel_bDw_"+observable).c_str());
      histoJesUp = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopel_metJetUp_"+observable).c_str());
      histoJesDw = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopel_metJetDw_"+observable).c_str());
      histoJerUp = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopel_metResUp_"+observable).c_str());
      histoJerDw = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopel_metResDw_"+observable).c_str());
      histoUncUp = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopel_metUncUp_"+observable).c_str());
      histoUncDw = (TH1F*)templatesfile->FindObjectAny(("dbkghisttopel_metUncDw_"+observable).c_str());
      
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
      
      
      addTemplate("Dibosons_TE_"+suffix+"_CMS_btagUp", vars, *wspace_TE, histobUp);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_btagDown", vars, *wspace_TE, histobDw);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_scale_jUp", vars, *wspace_TE, histoJesUp);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_scale_jDown", vars, *wspace_TE, histoJesDw);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_res_jUp", vars, *wspace_TE, histoJerUp);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_res_jDown", vars, *wspace_TE, histoJerDw);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_uncUp", vars, *wspace_TE, histoUncUp);
      addTemplate("Dibosons_TE_"+suffix+"_CMS_uncDown", vars, *wspace_TE, histoUncDw);
      
      
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
