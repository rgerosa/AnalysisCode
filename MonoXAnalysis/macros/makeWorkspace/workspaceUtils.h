#ifndef WORKSPACEUTILS_H
#define WORKSPACEUTILS_H

#include <vector>
#include <fstream>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "TString.h"

#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooAddition.h"
#include "../../../../HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#include "../makeTemplates/histoUtils.h"
#include "../makeTemplates/histoUtils2D.h"

const int binMaxForShapeSys = 100;

using namespace std;


class systematicCutAndCount {

 public:

 systematicCutAndCount(const string & sysName, TH1* num_1, TH1* den_1, TH1* num_2, TH1* den_2):
  sysName(sysName),
  num_1(num_1),
    den_1(den_1),
    num_2(num_2),
    den_2(den_2){
    };

  systematicCutAndCount(){};

  string sysName;
  TH1* num_1;
  TH1* den_1;
  TH1* num_2;
  TH1* den_2;  

};

void checkNegativeBin(TH1* histo){
  for(int iBin = 0; iBin < histo->GetNbinsX()+1; iBin++){
    if(histo->GetBinContent(iBin+1) < 0)
      histo->SetBinContent(iBin+1,0);
  }
}

// function to create a RooDataHist from TH1F and import it in a workspace                                                                                                   
void addTemplate(string procname, 
		 const RooArgList& varlist, 
		 RooWorkspace& ws, 
		 TH1F* hist, 
		 const bool & isCoutAndCount = false) {

  if(hist == 0 || hist == NULL){
    cerr<<"addTemplate: found an null pointer --> check "<<endl;
    return;
  }

  if(hist->Integral() == 0)
    addDummyBinContent(hist); // avoind empty histograms in the workspace                                                                                                     

  checkNegativeBin(hist);

  if(not isCoutAndCount){
    RooDataHist rhist((procname).c_str(), "", varlist, hist);
    ws.import(rhist);
  }
  else{
    RooRealVar *var = dynamic_cast<RooRealVar*>(varlist.at(0));
    vector<double> bin;
    bin.push_back(var->getMin());
    bin.push_back(var->getMax());
    TH1F* hist_temp = new TH1F(Form("%s_cutAndCount",hist->GetName()),"",bin.size()-1,&bin[0]);
    hist_temp->SetBinContent(var->getBins(),hist->Integral(hist->FindBin(var->getMin()),hist->FindBin(var->getMax())));
    RooDataHist rhist((procname).c_str(), "", varlist,hist_temp);
    ws.import(rhist);
  }
}

// function to generate stat variations bin-by-bin and put them into a Worksapce
void generateStatTemplate(string procname, 
			  const RooArgList& varlist, 
			  RooWorkspace& ws, 
			  TH1* histo, 
			  float scaleUncertainty = 1., 
			  const bool & isCutAndCount = false){
  vector<TH1F*> histStatUp;
  vector<TH1F*> histStatDw;

  if(histo == 0 || histo == NULL){
    cerr<<"generateStatTemplate: found an null pointer --> check "<<endl;
    return;
  }
  if(histo->Integral() == 0)
    addDummyBinContent(histo);

  if(not isCutAndCount){    
    for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
      histStatUp.push_back((TH1F*) histo->Clone(Form("%s_%s_bin_%i_statUp_tmp",procname.c_str(),procname.c_str(),iBin+1)));
      histStatDw.push_back((TH1F*) histo->Clone(Form("%s_%s_bin_%i_statDown_tmp",procname.c_str(),procname.c_str(),iBin+1)));
    }
    
    for( size_t iHisto =0; iHisto < histStatUp.size(); iHisto++){
      histStatUp.at(iHisto)->SetBinContent(iHisto+1,histo->GetBinContent(iHisto+1)+histo->GetBinError(iHisto+1)*scaleUncertainty);
      if(histo->GetBinContent(iHisto+1)-histo->GetBinError(iHisto+1)*scaleUncertainty >= 0.)
	histStatDw.at(iHisto)->SetBinContent(iHisto+1,histo->GetBinContent(iHisto+1)-histo->GetBinError(iHisto+1)*scaleUncertainty);
      else
	histStatDw.at(iHisto)->SetBinContent(iHisto+1,histo->GetBinContent(iHisto+1)/10.);
    }
  }
  else{ // for cut and count stat uncertainty
    histStatUp.push_back((TH1F*) histo->Clone(Form("%s_%s_statUp_tmp",procname.c_str(),procname.c_str())));
    histStatDw.push_back((TH1F*) histo->Clone(Form("%s_%s_statDown_tmp",procname.c_str(),procname.c_str())));
    for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
      histStatUp.at(0)->SetBinContent(iBin+1,histo->GetBinContent(iBin+1)+histo->GetBinError(iBin+1)*scaleUncertainty);
      histStatDw.at(0)->SetBinContent(iBin+1,histo->GetBinContent(iBin+1)-histo->GetBinError(iBin+1)*scaleUncertainty);
    }
  }

  for(size_t iHisto =0; iHisto < histStatUp.size(); iHisto++) {
    stringstream name;
    // add two times proc name in order to obtain correct lines in the datacard without name overlapping                                                                      
    if(not isCutAndCount)
      name << procname << "_"<< procname << "_CMS_bin" << iHisto+1 << "_statUp";
    else
      name << procname << "_"<< procname << "_CMS_bin_statUp";
    
    addTemplate(name.str(),varlist,ws,histStatUp[iHisto],isCutAndCount);
  }
  for(size_t iHisto =0; iHisto < histStatDw.size(); iHisto++){
    stringstream name;
    // add two times proc name in order to obtain correct lines in the datacard without name overlapping                                                                      
    if(not isCutAndCount)
      name << procname << "_" << procname << "_CMS_bin" << iHisto+1 << "_statDown";
    else
      name << procname << "_" << procname << "_CMS_bin_statDown";
    
    addTemplate(name.str(),varlist,ws,histStatDw[iHisto],isCutAndCount);
  }
}

// cloneAndRescale
TH1* cloneAndRescale(TH1* input, const float & rescale, const string & postfix){
  if(postfix == ""){
    input->Scale(rescale);
    return input;
  }
  else{
    TH1* output = (TH1*) input->Clone(postfix.c_str());
    output->Scale(rescale);
    return output;
  }
}


// to add shapeVariation to a given process
void addShapeVariations(const string & inputName, 
			const string & workspaceName,
			const string & suffix, 
			const string & observable,
			const RooArgList& varlist, 
			RooWorkspace& workspace, 
			TFile* templateFile, 
			const string & postfix = "",
			bool  isCombination    = false,
			const bool & isCoutAndCount = false,
			const float & normalizeSignal = -99){

  //Access to the different systematic variations
  TH1F* nominalHisto = NULL;
  TH1F* histoJesUp   = NULL;
  TH1F* histoJesDw   = NULL;
  TH1F* histoJerUp   = NULL;
  TH1F* histoJerDw   = NULL;
  TH1F* histoUncUp   = NULL;
  TH1F* histoUncDw   = NULL;

  RooRealVar* var = dynamic_cast<RooRealVar*>(varlist.at(0));

  if(postfix != ""){
    nominalHisto = (TH1F*)templateFile->FindObjectAny((inputName+"_"+postfix+"_"+observable).c_str());
    histoJesUp  = (TH1F*)templateFile->FindObjectAny((inputName+"_metJetUp_"+postfix+"_"+observable).c_str());
    histoJesDw  = (TH1F*)templateFile->FindObjectAny((inputName+"_metJetDw_"+postfix+"_"+observable).c_str());
    histoJerUp  = (TH1F*)templateFile->FindObjectAny((inputName+"_metResUp_"+postfix+"_"+observable).c_str());
    histoJerDw  = (TH1F*)templateFile->FindObjectAny((inputName+"_metResDw_"+postfix+"_"+observable).c_str());
    histoUncUp  = (TH1F*)templateFile->FindObjectAny((inputName+"_metUncUp_"+postfix+"_"+observable).c_str());
    histoUncDw  = (TH1F*)templateFile->FindObjectAny((inputName+"_metUncDw_"+postfix+"_"+observable).c_str());
  }
  else{
    nominalHisto = (TH1F*)templateFile->FindObjectAny((inputName+"_"+observable).c_str());
    histoJesUp  = (TH1F*)templateFile->FindObjectAny((inputName+"_metJetUp_"+observable).c_str());
    histoJesDw  = (TH1F*)templateFile->FindObjectAny((inputName+"_metJetDw_"+observable).c_str());
    histoJerUp  = (TH1F*)templateFile->FindObjectAny((inputName+"_metResUp_"+observable).c_str());
    histoJerDw  = (TH1F*)templateFile->FindObjectAny((inputName+"_metResDw_"+observable).c_str());
    histoUncUp  = (TH1F*)templateFile->FindObjectAny((inputName+"_metUncUp_"+observable).c_str());
    histoUncDw  = (TH1F*)templateFile->FindObjectAny((inputName+"_metUncDw_"+observable).c_str());
  }

  if(nominalHisto){

    if(normalizeSignal > 0 and isCoutAndCount){
      if(histoJesUp and histoJesUp->Integral(histoJesUp->FindBin(var->getMin()),histoJesUp->FindBin(var->getMax())) != 0)
	histoJesUp->Scale(normalizeSignal/histoJesUp->Integral(histoJesUp->FindBin(var->getMin()),histoJesUp->FindBin(var->getMax())));
      if(histoJesDw and histoJesDw->Integral(histoJesDw->FindBin(var->getMin()),histoJesDw->FindBin(var->getMax())) != 0)
	histoJesDw->Scale(normalizeSignal/histoJesDw->Integral(histoJesDw->FindBin(var->getMin()),histoJesDw->FindBin(var->getMax())));
      if(histoJerUp and histoJerUp->Integral(histoJerUp->FindBin(var->getMin()),histoJerUp->FindBin(var->getMax())) != 0)
	histoJerUp->Scale(normalizeSignal/histoJerUp->Integral(histoJerUp->FindBin(var->getMin()),histoJerUp->FindBin(var->getMax())));
      if(histoJerDw and histoJerDw->Integral(histoJerDw->FindBin(var->getMin()),histoJerDw->FindBin(var->getMax())) != 0)
	histoJerDw->Scale(normalizeSignal/histoJerDw->Integral(histoJerDw->FindBin(var->getMin()),histoJerDw->FindBin(var->getMax())));
      if(histoUncUp and histoUncUp->Integral(histoUncUp->FindBin(var->getMin()),histoUncUp->FindBin(var->getMax())) != 0)
	histoUncUp->Scale(normalizeSignal/histoUncUp->Integral(histoUncUp->FindBin(var->getMin()),histoUncUp->FindBin(var->getMax())));
      if(histoUncDw and histoUncDw->Integral(histoUncDw->FindBin(var->getMin()),histoUncDw->FindBin(var->getMax())) != 0)
	histoUncDw->Scale(normalizeSignal/histoUncDw->Integral(histoUncDw->FindBin(var->getMin()),histoUncDw->FindBin(var->getMax())));      
    }

    if(nominalHisto->GetNbinsX() > binMaxForShapeSys){
      fixShapeUncertainty(nominalHisto,histoJesUp,int(nominalHisto->GetNbinsX()/2),1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,int(nominalHisto->GetNbinsX()/2),0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,int(nominalHisto->GetNbinsX()/2),1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,int(nominalHisto->GetNbinsX()/2),0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,int(nominalHisto->GetNbinsX()/2),1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,int(nominalHisto->GetNbinsX()/2),0.99);
    }

    if(not isCombination){
      addTemplate(workspaceName+"_"+suffix+"_CMS_scale_j_13TeVUp",varlist,workspace,histoJesUp,isCoutAndCount);
      addTemplate(workspaceName+"_"+suffix+"_CMS_scale_j_13TeVDown",varlist,workspace,histoJesDw,isCoutAndCount);
      addTemplate(workspaceName+"_"+suffix+"_CMS_res_j_13TeVUp",varlist,workspace,histoJerUp,isCoutAndCount);
      addTemplate(workspaceName+"_"+suffix+"_CMS_res_j_13TeVDown",varlist,workspace,histoJerDw,isCoutAndCount);
      addTemplate(workspaceName+"_"+suffix+"_CMS_scale_met_13TeVUp",varlist,workspace,histoUncUp,isCoutAndCount);
      addTemplate(workspaceName+"_"+suffix+"_CMS_scale_met_13TeVDown",varlist,workspace,histoUncDw,isCoutAndCount);
    }
    else{
      addTemplate(workspaceName+"_"+suffix+"_CMS_mono_scale_j_13TeVUp",varlist,workspace,histoJesUp,isCoutAndCount);
      addTemplate(workspaceName+"_"+suffix+"_CMS_mono_scale_j_13TeVDown",varlist,workspace,histoJesDw,isCoutAndCount);
      addTemplate(workspaceName+"_"+suffix+"_CMS_mono_res_j_13TeVUp",varlist,workspace,histoJerUp,isCoutAndCount);
      addTemplate(workspaceName+"_"+suffix+"_CMS_mono_res_j_13TeVDown",varlist,workspace,histoJerDw,isCoutAndCount);
      addTemplate(workspaceName+"_"+suffix+"_CMS_scale_metUp",varlist,workspace,histoUncUp,isCoutAndCount);
      addTemplate(workspaceName+"_"+suffix+"_CMS_scale_metDown",varlist,workspace,histoUncDw,isCoutAndCount);
    }
  }   
} 


// Create a basic RooParamtricHist
// Make list of bins of a TH1F as RooArgList of RooRealVar and building the RooParametricHist (to be Run in the release with combine)                                       
void makeBinList(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* hist, RooArgList& binlist, bool setConst = false, const bool & isCoutAndCount = false) {

  if(hist == 0 || hist == NULL){
    cerr<<"makeBinList: found an null pointer --> check "<<endl;
    return;
  }
  if(hist->Integral() == 0){
    addDummyBinContent(hist);
  }
  
  if(not isCoutAndCount){

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
  else{

    // name is simply the process, any reference to the bin number 
    stringstream binss;
    binss << procname <<"_bin";
    RooRealVar* binvar;

    // make a RooRealVar for each bin [0,2*binContent]                                                                                                                        
    if (!setConst)
      binvar = new RooRealVar(binss.str().c_str(), "", hist->Integral(hist->FindBin(var.getMin()),hist->FindBin(var.getMax())), 
			      0.,hist->Integral(hist->FindBin(var.getMin()),hist->FindBin(var.getMax()))*10.0);
    else
      binvar = new RooRealVar(binss.str().c_str(), "", hist->Integral(hist->FindBin(var.getMin()),hist->FindBin(var.getMax())));
    
    binlist.add(*binvar);

    vector<double> bin;
    bin.push_back(var.getMin());
    bin.push_back(var.getMax());
    TH1F* hist_temp = new TH1F(Form("%s_cutAndCount",hist->GetName()),"",bin.size()-1,&bin[0]);
    hist_temp->SetBinContent(var.getBins(),hist->Integral(hist->FindBin(var.getMin()),hist->FindBin(var.getMax())));

    stringstream normss;
    normss << procname << "_norm";
    // build a RooParamtric hist using the bin list and the histogram                                                                                                           
    RooParametricHist phist(procname.c_str(), "", var, binlist, *hist_temp);
    RooAddition norm(normss.str().c_str(), "", binlist);
    ws.import(phist,RooFit::RecycleConflictNodes());
    ws.import(norm, RooFit::RecycleConflictNodes());
  }
}

// make connections betweem signal region and control region                                                                                                                
void makeConnectedBinList(string procname,  // name to be used to fill the workspace 
			  RooRealVar& var,  // observable RooRealVar
                          RooWorkspace& ws, // RooWorkspace
			  TH1F* rhist,      // TF histogram
                          vector<pair<RooRealVar*, TH1*> > syst, // set of systematics to be considered
			  const RooArgList& srbinlist,  // Binning into the signal region
                          RooArgList* crbinlist = NULL, // constrol region bins
                          string observable     = "met", // observable name
			  bool   applyUncertaintyOnNumerator = false,
			  bool   addStatUncertainty = true,
			  float  reductionFactorLastBin = 1
			  ) {


  if(rhist == 0 || rhist == NULL){
    cerr<<"makeConnectedBinList: found an null pointer --> check "<<endl;
    return;
  }

  if(rhist->Integral() == 0)
    addDummyBinContent(rhist);  

  // bin list for the CR                                                                                                                                               
  if (crbinlist == NULL){
    cerr<<"makeConnectedBinList: crbinlist empty --> create "<<endl;
    crbinlist = new RooArgList();
  }


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
    if(addStatUncertainty == false)
      rerrbinvar->setConstant(kTRUE);
    
    stringstream binss;
    binss << procname << "_bin" << i ;
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

    if(applyUncertaintyOnNumerator)
      formss << ")";

    if(i != rhist->GetNbinsX()){
      if(rhist->GetBinContent(i) == 0)
	formss << "*(TMath::Max(0,1+1*@2))";
      else if(rhist->GetBinError(i)/rhist->GetBinContent(i) >= 0)
	formss << "*(TMath::Max(0,1+" << rhist->GetBinError(i)/rhist->GetBinContent(i) << "*@2))";
      else
	formss << "*(TMath::Max(0,1-" << fabs(rhist->GetBinError(i)/rhist->GetBinContent(i)) << "*@2))";
    }
    else{
      if(rhist->GetBinContent(i) == 0)
        formss << "*(TMath::Max(0,1+1*@2))";	
      else if(rhist->GetBinError(i)/rhist->GetBinContent(i) >= 0)
        formss << "*(TMath::Max(0,1+" << (rhist->GetBinError(i)/reductionFactorLastBin)/rhist->GetBinContent(i) << "*@2))";
      else
        formss << "*(TMath::Max(0,1-" << (fabs(rhist->GetBinError(i)/reductionFactorLastBin)/rhist->GetBinContent(i)) << "*@2))"; 
    }

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
	  syst[j].first->setMin(-fabs(extreme));
	  syst[j].first->setMax(fabs(extreme));
	  extreme_tmp = extreme;
	}
	fobinlist.add(*syst[j].first);
      }
      if(syst[j].second->GetBinContent(i) >= 0)
	formss << "*(TMath::Max(0,1+" << syst[j].second->GetBinContent(i) << "*@" << j+3 << "))";
      else
	formss << "*(TMath::Max(0,1-" << fabs(syst[j].second->GetBinContent(i)) << "*@" << j+3 << "))";
      
    }

    if(not applyUncertaintyOnNumerator)
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


// make connections betweem signal region and control region                                                                                                                
void makeConnectedBinListCutAndCount(string procname,  // name to be used to fill the workspace 
				     RooRealVar& var,  // observable RooRealVar
				     RooWorkspace& ws, // RooWorkspace
				     TH1F* nhist,      // TF numerator histogram
				     TH1F* dhist,      // TF denominator histogram
				     vector<pair<RooRealVar*, systematicCutAndCount> > syst, // set of systematics to be considered
				     const RooArgList& srbinlist,  // Binning into the signal region
				     RooArgList* crbinlist = NULL, // constrol region bins
				     string observable = "met",     // observable name
				     bool   addStatUncertainty = true,
				     bool   applyUncertaintyOnNumerator = false
				     ) {


  if(nhist == 0 || nhist == NULL){
    cerr<<"makeConnectedBinListCutAndCount: found an null pointer  for numerator --> check "<<endl;
    return;
  }

  if(dhist == 0 || dhist == NULL){
    cerr<<"makeConnectedBinListCutAndCount: found an null pointer  for denominator --> check "<<endl;
    return;
  }

  if(nhist->Integral() == 0)
    addDummyBinContent(nhist);  

  if(dhist->Integral() == 0)
    addDummyBinContent(dhist);  

  // bin list for the CR                                                                                                                                               
  if (crbinlist == NULL){
    cerr<<"makeConnectedBinListCutAndCount: crbinlist empty --> create "<<endl;
    crbinlist = new RooArgList();
  }

  // determi the transfer factor value for the cut and count  --> ratio of integrals
  stringstream rbinss;
  rbinss << "r_" << procname;

  double integral_num_err = 0.;
  double integral_num = nhist->IntegralAndError(nhist->FindBin(var.getMin()),nhist->FindBin(var.getMax()),integral_num_err);
  double integral_den_err = 0.;
  double integral_den = dhist->IntegralAndError(dhist->FindBin(var.getMin()),dhist->FindBin(var.getMax()),integral_den_err);
  // z = Int(n)/Int(d) --> error propagation --> sigmaZ = sqrt((1/Int(d))^2*error(Int(n))^2+(Int(n)/Int(d)^2)^2*error(Int(d)^2)
  double ratio_err = sqrt((1/(pow(integral_den,2)))*pow(integral_num_err,2)+(pow(integral_num,2)/pow(integral_den,4))*pow(integral_den_err,2));
  RooRealVar* rbinvar = new RooRealVar(rbinss.str().c_str(), "",integral_num/integral_den);
  // stat uncertainty 
  stringstream rerrbinss;
  rerrbinss << procname << "_Runc";
  double extreme = min(5.,0.9*rbinvar->getVal()/ratio_err);
  RooRealVar* rerrbinvar = new RooRealVar(rerrbinss.str().c_str(), "", 0., -extreme, extreme);
  if(not addStatUncertainty)
    rerrbinvar->setConstant(kTRUE);
  
  stringstream binss;
  binss << procname+"_bin";
  // list of bins                                                                                                                                                          
  RooArgList fobinlist;
  // signal region bins    
  if(srbinlist.getSize() > 1) {
    cerr<<"makeConnectedBinListCutAndCount : more than one bin defined in the signal region, please check ";
    return;
  }

  fobinlist.add(srbinlist[0]);
  // connection bins from rHisto (ratio histo)                                                                                                                             
  fobinlist.add(*rbinvar);
  // uncertainty                                                                                                                                                          
  fobinlist.add(*rerrbinvar);
    
  // CR = SR /( R *(1+RhistError/Rbin*Rbin_Err)) --> statstical uncertainty                                                         
  stringstream formss;
  formss << "@0/";
  formss << "(";
  formss << "@1";
  if(applyUncertaintyOnNumerator)
    formss << ")";
    
  if(ratio_err/rbinvar->getVal() >= 0)
    formss << "*(TMath::Max(0,1+" << ratio_err/rbinvar->getVal() << "*@2))";
  else
    formss << "*(TMath::Max(0,1-" <<fabs(ratio_err/rbinvar->getVal())<< "*@2))";

  float extreme_tmp = 5;

  // systemaitc uncertainty                                                                                                                                               
  for (size_t j = 0; j < syst.size(); j++) {
    stringstream systbinss;
    if (syst[j].first == NULL) { 
      cerr<<"makeConnectedBinListCutAndCount: problem --> no way to add bin by bin uncertainty in a cut and count analysis --> return and check"<<endl;
      return;
    }
    else{

      // ratio for the first histo
      double integral_num_sys_1 = syst[j].second.num_1->Integral(syst[j].second.num_1->FindBin(var.getMin()),syst[j].second.num_1->FindBin(var.getMax()));
      double integral_den_sys_1 = syst[j].second.den_1->Integral(syst[j].second.den_1->FindBin(var.getMin()),syst[j].second.den_1->FindBin(var.getMax()));
      double ratio_1 = integral_num_sys_1/integral_den_sys_1;

      // ratio for the second histo
      double integral_num_sys_2 = syst[j].second.num_2->Integral(syst[j].second.num_2->FindBin(var.getMin()),syst[j].second.num_2->FindBin(var.getMax()));
      double integral_den_sys_2 = syst[j].second.den_2->Integral(syst[j].second.den_2->FindBin(var.getMin()),syst[j].second.den_2->FindBin(var.getMax()));
      double ratio_2 = integral_num_sys_2/integral_den_sys_2;

      double sys_value = fabs(ratio_1/ratio_2-1.);      
      double extreme = min(5.,0.9*rbinvar->getVal()/sys_value);

      if(extreme < extreme_tmp){
	syst[j].first->setMin(-extreme);
	syst[j].first->setMax(extreme);      
	extreme_tmp = extreme;
      }

      fobinlist.add(*syst[j].first);
      if(sys_value >= 0)
	formss << "*(TMath::Max(0,1+" << sys_value << "*@" << j+3 << "))";      
      else
	formss << "*(TMath::Max(0,1-" << fabs(sys_value) << "*@" << j+3 << "))";      

    }
  }

  if(not applyUncertaintyOnNumerator)
    formss << ")";
    
  // create a single RooFormulaVar                                                                                                                                         
  RooFormulaVar* binvar = new RooFormulaVar(binss.str().c_str(), "", formss.str().c_str(), RooArgList(fobinlist));
  crbinlist->add(*binvar);

  stringstream normss;
  normss << procname << "_norm";
  
  vector<double> bin;
  bin.push_back(var.getMin());
  bin.push_back(var.getMax());
  TH1F* hist_temp = new TH1F(Form("%s_cutAndCount",nhist->GetName()),"",bin.size()-1,&bin[0]);
  RooParametricHist phist(procname.c_str(), "", var, *crbinlist, *hist_temp);
  RooAddition norm(normss.str().c_str(),"", *crbinlist);

  ws.import(phist,RooFit::RecycleConflictNodes());
  ws.import(norm, RooFit::RecycleConflictNodes());
  
}

#endif
