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
#include "../../../HiggsAnalysis/CombinedLimit/interface/RooParametricHist.h"
#include "histoUtils.h"
#include "histoUtils2D.h"

const int binMaxForShapeSys = 10;

using namespace std;

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

// function to generate stat variations bin-by-bin and put them into a Worksapce
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


// to add shapeVariation to a given process
void addShapeVariations(const string & inputName, 
			const string & workspaceName,
			const string & suffix, 
			const string & observable,
			RooArgList& varlist, RooWorkspace& workspace, TFile* templateFile, 
			const string & postfix = "",
			bool  isCombination    = false){

  //Access to the different systematic variations
  TH1F* nominalHisto = (TH1F*)templateFile->FindObjectAny((inputName+"_"+postfix+"_"+observable).c_str());
  TH1F* histobUp  = (TH1F*)templateFile->FindObjectAny((inputName+"_bUp_"+postfix+"_"+observable).c_str());
  TH1F* histobDw  = (TH1F*)templateFile->FindObjectAny((inputName+"_bDw_"+postfix+"_"+observable).c_str());
  TH1F* histoJesUp  = (TH1F*)templateFile->FindObjectAny((inputName+"_metJetUp_"+postfix+"_"+observable).c_str());
  TH1F* histoJesDw  = (TH1F*)templateFile->FindObjectAny((inputName+"_metJetDw_"+postfix+"_"+observable).c_str());
  TH1F* histoJerUp  = (TH1F*)templateFile->FindObjectAny((inputName+"_metResUp_"+postfix+"_"+observable).c_str());
  TH1F* histoJerDw  = (TH1F*)templateFile->FindObjectAny((inputName+"_metResDw_"+postfix+"_"+observable).c_str());
  TH1F* histoUncUp  = (TH1F*)templateFile->FindObjectAny((inputName+"_metUncUp_"+postfix+"_"+observable).c_str());
  TH1F* histoUncDw  = (TH1F*)templateFile->FindObjectAny((inputName+"_metUncDw_"+postfix+"_"+observable).c_str());

  if(nominalHisto){

    if(nominalHisto->GetNbinsX() > binMaxForShapeSys){
      fixShapeUncertainty(nominalHisto,histobUp,int(nominalHisto->GetNbinsX()/2),1.02);
      fixShapeUncertainty(nominalHisto,histobDw,int(nominalHisto->GetNbinsX()/2),0.98);
      fixShapeUncertainty(nominalHisto,histoJesUp,int(nominalHisto->GetNbinsX()/2),1.06);
      fixShapeUncertainty(nominalHisto,histoJesDw,int(nominalHisto->GetNbinsX()/2),0.94);
      fixShapeUncertainty(nominalHisto,histoJerUp,int(nominalHisto->GetNbinsX()/2),1.02);
      fixShapeUncertainty(nominalHisto,histoJerDw,int(nominalHisto->GetNbinsX()/2),0.98);
      fixShapeUncertainty(nominalHisto,histoUncUp,int(nominalHisto->GetNbinsX()/2),1.01);
      fixShapeUncertainty(nominalHisto,histoUncDw,int(nominalHisto->GetNbinsX()/2),0.99);
    }

    addTemplate(workspaceName+"_"+suffix+"_CMS_btag_13TeVUp",varlist,workspace,histobUp);
    addTemplate(workspaceName+"_"+suffix+"_CMS_btag_13TeVDown",varlist,workspace,histobDw);
    if(not isCombination){
      addTemplate(workspaceName+"_"+suffix+"_CMS_scale_j_13TeVUp",varlist,workspace,histoJesUp);
      addTemplate(workspaceName+"_"+suffix+"_CMS_scale_j_13TeVDown",varlist,workspace,histoJesDw);
      addTemplate(workspaceName+"_"+suffix+"_CMS_res_j_13TeVUp",varlist,workspace,histoJerUp);
      addTemplate(workspaceName+"_"+suffix+"_CMS_res_j_13TeVDown",varlist,workspace,histoJerDw);
      addTemplate(workspaceName+"_"+suffix+"_CMS_scale_met_13TeVUp",varlist,workspace,histoUncUp);
      addTemplate(workspaceName+"_"+suffix+"_CMS_scale_met_13TeVDown",varlist,workspace,histoUncDw);
    }
    else{
      addTemplate(workspaceName+"_"+suffix+"_CMS_mono_scale_j_13TeVUp",varlist,workspace,histoJesUp);
      addTemplate(workspaceName+"_"+suffix+"_CMS_mono_scale_j_13TeVDown",varlist,workspace,histoJesDw);
      addTemplate(workspaceName+"_"+suffix+"_CMS_mono_res_j_13TeVUp",varlist,workspace,histoJerUp);
      addTemplate(workspaceName+"_"+suffix+"_CMS_mono_res_j_13TeVDown",varlist,workspace,histoJerDw);
      addTemplate(workspaceName+"_"+suffix+"_CMS_scale_metUp",varlist,workspace,histoUncUp);
      addTemplate(workspaceName+"_"+suffix+"_CMS_scale_metDown",varlist,workspace,histoUncDw);
    }
  }   
} 


// Create a basic RooParamtricHist
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

  if(rhist->Integral() == 0)
    addDummyBinContent(rhist);  

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



#endif
