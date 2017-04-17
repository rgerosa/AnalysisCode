#include "../CMS_lumi.h"

enum class Nuisance {theory, ratios, experimental, stat, otherbkg};

void plotNuisance(TCanvas* canvas, map<TString,RooRealVar*> & listParamPreFit, map<TString,RooRealVar*> & listParamPostFit,string outputDIR, string postfix){
  
  TH1F* histopull = new TH1F(("histopull"+postfix).c_str(),"",listParamPostFit.size(),0,listParamPostFit.size()+1);
  TH1F* histopull_band_1s = new TH1F(("histopull_band_1s"+postfix).c_str(),"",listParamPostFit.size(),0,listParamPostFit.size()+1);
  TH1F* histopull_band_2s = new TH1F(("histopull_band_2s"+postfix).c_str(),"",listParamPostFit.size(),0,listParamPostFit.size()+1);
  histopull->Sumw2();
  histopull_band_1s->Sumw2();
  histopull_band_2s->Sumw2();

  int iBin= 0;
  for(auto& postfit : listParamPostFit){
    histopull->SetBinContent(iBin+1,postfit.second->getVal()-listParamPreFit[postfit.first]->getVal());
    histopull->SetBinError(iBin+1,postfit.second->getError());
    histopull_band_1s->SetBinContent(iBin+1,listParamPreFit[postfit.first]->getVal());
    histopull_band_1s->SetBinError(iBin+1,listParamPreFit[postfit.first]->getError());
    histopull_band_2s->SetBinContent(iBin+1,listParamPreFit[postfit.first]->getVal());
    histopull_band_2s->SetBinError(iBin+1,2*listParamPreFit[postfit.first]->getError());
    histopull_band_2s->GetXaxis()->SetBinLabel(iBin+1,listParamPreFit[postfit.first]->GetName());      
    iBin++;
  }
  
  histopull_band_2s->GetXaxis()->SetTitle("");
  histopull_band_2s->GetYaxis()->SetTitle("(#theta_{post}-#theta_{pre})/#sigma_{post}");
  histopull_band_2s->SetFillColor(kOrange);
  histopull_band_2s->SetMarkerColor(kOrange);
  histopull_band_2s->SetLineColor(kOrange);
  histopull_band_2s->GetYaxis()->SetRangeUser(-3,3);
  histopull_band_2s->Draw("E2");
  histopull_band_2s->GetXaxis()->LabelsOption("v");
  histopull_band_2s->GetXaxis()->SetNdivisions(-414);

  histopull_band_1s->SetFillColor(kGreen);
  histopull_band_1s->SetMarkerColor(kGreen);
  histopull_band_1s->SetLineColor(kGreen);
  histopull_band_1s->Draw("E2same");
  
  histopull->SetMarkerColor(kBlack);
  histopull->SetMarkerStyle(20);
  histopull->SetLineColor(kBlack);
  histopull->SetMarkerSize(1.2);
  histopull->SetLineWidth(2);
  histopull->SetFillColor(kBlack);
  histopull->SetFillStyle(3001);
  histopull->Draw("E2same");
  histopull->Draw("E1Psame");

  
  CMS_lumi(canvas,"35.9");

  canvas->RedrawAxis("sameaxis");
  
  canvas->SaveAs((outputDIR+"/pull_nuisance_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/pull_nuisance_"+postfix+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/pull_nuisance_"+postfix+".root").c_str(),"root");


}

void makePostFitNuisancePull(string inputFileName, string outputDIR, Nuisance nuisancetype, bool plotSBFit = false){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

 
  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  RooFitResult* res = NULL;
  if(plotSBFit)
    res = (RooFitResult*) inputFile->Get("fit_s");
  else
    res = (RooFitResult*) inputFile->Get("fit_b");


  RooArgList parlist_final = res->floatParsFinal();
  RooArgList parlist_init  = RooArgList(*((RooArgSet*) inputFile->Get("nuisances_prefit")));
  
  map<TString,RooRealVar*> listParamPreFit;
  map<TString,RooRealVar*> listParamPostFit;

  map<TString,RooRealVar*> listParamPreFit_ZG;
  map<TString,RooRealVar*> listParamPreFit_ZW;
  map<TString,RooRealVar*> listParamPreFit_WW;
  map<TString,RooRealVar*> listParamPreFit_ZZ;

  map<TString,RooRealVar*> listParamPostFit_ZG;
  map<TString,RooRealVar*> listParamPostFit_ZW;
  map<TString,RooRealVar*> listParamPostFit_WW;
  map<TString,RooRealVar*> listParamPostFit_ZZ;
  
  for(int isize = 0; isize < parlist_init.getSize(); isize++){
    TString name (parlist_init.at(isize)->GetName());
    if(nuisancetype == Nuisance::theory){
      if((name.Contains("QCDScale") and (name.Contains("ZW") or name.Contains("ZG"))) or name.Contains("QCDShape") or name.Contains("QCDProcess") or name.Contains("PDF") or name.Contains("NNLOMiss") or name.Contains("Sudakov") or name.Contains("QCDEWKMIX") or name.Contains("NNLOEWK")){
	listParamPreFit[name] = (RooRealVar*) parlist_init.at(isize);
      }
    }
    else if(nuisancetype == Nuisance::experimental){
      if(name.Contains("CMS_eff") or name.Contains("CMS_trig") or name.Contains("CMS_met") or name.Contains("CMS_btag") or name.Contains("CMS_reco") or name.Contains("veto") or 
	 name.Contains("CMS_Vtag"))
	listParamPreFit[name] = (RooRealVar*) parlist_init.at(isize);
    }
    else if(nuisancetype == Nuisance::stat){
      if(name.Contains("bin") and name.Contains("Runc"))
	listParamPreFit[name] = (RooRealVar*) parlist_init.at(isize);
    }
    else if(nuisancetype == Nuisance::ratios){

      if((name.Contains("ZM") and name.Contains("Runc")) or (name.Contains("ZE") and name.Contains("Runc")) or name.Contains("CMS_eff_m") or name.Contains("CMS_reco_m") or name.Contains("CMS_eff_e") or name.Contains("CMS_reco_e") or name.Contains("CMS_met_trig"))
	listParamPreFit_ZZ[name] = (RooRealVar*) parlist_init.at(isize);

      if((name.Contains("WM") and name.Contains("Runc")) or (name.Contains("WE") and name.Contains("Runc")) or name.Contains("CMS_eff_m") or name.Contains("CMS_reco_m") or name.Contains("CMS_eff_e") or name.Contains("CMS_reco_e") or name.Contains("CMS_met_trig") or name.Contains("CMS_eff_trig_e") or name.Contains("WtoWPDF") or name.Contains("muon_veto") or name.Contains("ele_veto") or name.Contains("tau_veto"))
	listParamPreFit_WW[name] = (RooRealVar*) parlist_init.at(isize);
      
      if((name.Contains("WJets_SR") and name.Contains("Runc"))  or name.Contains("ZW_SR") or (name.Contains("Znunu_SR") and not name.Contains("Znunu_SR_MJ") and not name.Contains("Znunu_SR_MV")) or name.Contains("WJets_SR"))
	listParamPreFit_ZW[name] = (RooRealVar*) parlist_init.at(isize);

      if((name.Contains("Znunu_GJ") and name.Contains("Runc"))  or name.Contains("ZG_GJ") or (name.Contains("Znunu_SR") and not name.Contains("Znunu_SR_MJ") and not name.Contains("Znunu_SR_MV")) or name.Contains("Gamma_GJ"))
	listParamPreFit_ZG[name] = (RooRealVar*) parlist_init.at(isize);
    }
    else if(nuisancetype == Nuisance::otherbkg){
      if((name.Contains("Top") or name.Contains("Diboson") or name.Contains("QCD_") or name.Contains("Purity") or name.Contains("GJets")) and not name.Contains("bin") and not name.Contains("stat"))
	listParamPreFit[name] = (RooRealVar*) parlist_init.at(isize);
    }
  }

  ///////
  for(int isize = 0; isize < parlist_final.getSize(); isize++){
    TString name (parlist_final.at(isize)->GetName());
    if(nuisancetype == Nuisance::theory){
      if((name.Contains("QCDScale") and (name.Contains("ZW") or name.Contains("ZG")))  or name.Contains("QCDShape") or name.Contains("QCDProcess") or name.Contains("PDF") or name.Contains("NNLOMiss") or name.Contains("Sudakov") or name.Contains("QCDEWKMIX") or name.Contains("NNLOEWK"))
	listParamPostFit[name] = (RooRealVar*) parlist_final.at(isize);
    }
    else if(nuisancetype == Nuisance::experimental){
      if(name.Contains("CMS_eff") or name.Contains("CMS_trig") or name.Contains("CMS_met") or name.Contains("CMS_btag") or name.Contains("CMS_reco") or name.Contains("veto") or 
	 name.Contains("CMS_Vtag"))
	listParamPostFit[name] = (RooRealVar*) parlist_final.at(isize);
    }
    else if(nuisancetype == Nuisance::stat){
      if(name.Contains("bin") and name.Contains("Runc"))
	listParamPostFit[name] = (RooRealVar*) parlist_final.at(isize);
    }
    else if(nuisancetype == Nuisance::ratios){

      if((name.Contains("ZM") and name.Contains("Runc")) or (name.Contains("ZE") and name.Contains("Runc")) or name.Contains("CMS_eff_m") or name.Contains("CMS_reco_m") or name.Contains("CMS_eff_e") or name.Contains("CMS_reco_e") or name.Contains("CMS_met_trig"))
	listParamPostFit_ZZ[name] = (RooRealVar*) parlist_final.at(isize);

      if((name.Contains("WM") and name.Contains("Runc")) or (name.Contains("WE") and name.Contains("Runc")) or name.Contains("CMS_eff_m") or name.Contains("CMS_reco_m") or name.Contains("CMS_eff_e") or name.Contains("CMS_reco_e") or name.Contains("CMS_met_trig") or name.Contains("CMS_eff_trig_e") or name.Contains("WtoWPDF") or name.Contains("muon_veto") or name.Contains("ele_veto") or name.Contains("tau_veto"))
	listParamPostFit_WW[name] = (RooRealVar*) parlist_final.at(isize);


      if((name.Contains("WJets_SR") and name.Contains("Runc"))  or name.Contains("ZW_SR") or (name.Contains("Znunu_SR") and not name.Contains("Znunu_SR_MJ") and not name.Contains("Znunu_SR_MV")) or name.Contains("WJets_SR"))
	listParamPostFit_ZW[name] = (RooRealVar*) parlist_final.at(isize);

      if((name.Contains("Znunu_GJ") and name.Contains("Runc"))  or name.Contains("ZG_GJ") or (name.Contains("Znunu_SR") and not name.Contains("Znunu_SR_MJ") and not name.Contains("Znunu_SR_MV")) or name.Contains("Gamma_GJ"))
	listParamPostFit_ZG[name] = (RooRealVar*) parlist_final.at(isize);
    }
    else if(nuisancetype == Nuisance::otherbkg){
      if((name.Contains("Top") or name.Contains("Diboson") or name.Contains("QCD_") or name.Contains("Purity") or name.Contains("GJets")) and not name.Contains("bin") and not name.Contains("stat"))
	listParamPostFit[name] = (RooRealVar*) parlist_final.at(isize);
    }
  }

  /////////////// ------ 
  if(nuisancetype != Nuisance::ratios){    
    
    TCanvas* canvas = NULL;
    if(nuisancetype != Nuisance::stat) {
      canvas = new TCanvas ("canvas","canvas",800,600);
      canvas->SetBottomMargin(0.3);
    }
    else{
      canvas = new TCanvas ("canvas","canvas",900,600);
      canvas->SetBottomMargin(0.35);
    }


    string postfix;
    if(nuisancetype == Nuisance::theory) postfix = "theory";
    else if(nuisancetype == Nuisance::experimental) postfix = "experimental";
    else if(nuisancetype == Nuisance::stat) postfix = "stat";    
    else if(nuisancetype == Nuisance::otherbkg) postfix = "otherbkgs";    
    plotNuisance(canvas,listParamPreFit,listParamPostFit,outputDIR,postfix);
    
  }
  else{
    
    TCanvas* canvas = new TCanvas ("canvas","canvas",900,600);
    canvas->SetBottomMargin(0.35);
    
    plotNuisance(canvas,listParamPreFit_ZZ,listParamPostFit_ZZ,outputDIR,"ZZ_ratio");
    plotNuisance(canvas,listParamPreFit_WW,listParamPostFit_WW,outputDIR,"WW_ratio");
    plotNuisance(canvas,listParamPreFit_ZW,listParamPostFit_ZW,outputDIR,"ZW_ratio");
    plotNuisance(canvas,listParamPreFit_ZG,listParamPostFit_ZG,outputDIR,"ZG_ratio");    
  }
}
			
