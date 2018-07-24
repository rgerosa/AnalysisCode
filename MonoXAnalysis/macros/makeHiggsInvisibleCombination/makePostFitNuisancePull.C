#include "../CMS_lumi.h"

enum class Category {monojet, monoV, monoZ, VBF, combined};

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
    TString name (listParamPreFit[postfit.first]->GetName());
    name.ReplaceAll("StatBounding","");
    name.ReplaceAll("_13TeV2016","");
    name.ReplaceAll("2016","");
    name.ReplaceAll("MVA","");
    name.ReplaceAll("_Runc","");
    name.ReplaceAll("stat_error","");
    name.ReplaceAll("__","_");
    name.ReplaceAll("__","_");
    histopull_band_2s->GetXaxis()->SetBinLabel(iBin+1,name);      
    iBin++;
  }

  histopull_band_2s->GetXaxis()->SetTitle("");
  histopull_band_2s->GetYaxis()->SetTitle("(#theta_{post}-#theta_{pre})/#sigma_{pre}");
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

//////
void makePostFitNuisancePull(string inputFileName, string outputDIR, Category category, bool plotSBFit = false){

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

  //// Look at experimental nuisances
  map<TString,RooRealVar*> listParamPreFit_experimental;
  map<TString,RooRealVar*> listParamPostFit_experimental;

  //// Look at theory nuisances
  map<TString,RooRealVar*> listParamPreFit_theory;
  map<TString,RooRealVar*> listParamPostFit_theory;

  //// Look at the stat bin-by-bin uncertainties
  map<TString,RooRealVar*> listParamPreFit_ZoverG;
  map<TString,RooRealVar*> listParamPreFit_ZoverW;
  map<TString,RooRealVar*> listParamPreFit_WoverW;
  map<TString,RooRealVar*> listParamPreFit_ZoverZ;
  map<TString,RooRealVar*> listParamPreFit_WZ3l;
  map<TString,RooRealVar*> listParamPreFit_ZZ4l;
  map<TString,RooRealVar*> listParamPreFit_ZZSR;
  map<TString,RooRealVar*> listParamPostFit_ZoverG;
  map<TString,RooRealVar*> listParamPostFit_ZoverW;
  map<TString,RooRealVar*> listParamPostFit_WoverW;
  map<TString,RooRealVar*> listParamPostFit_ZoverZ;
  map<TString,RooRealVar*> listParamPostFit_WZ3l;
  map<TString,RooRealVar*> listParamPostFit_ZZ4l;
  map<TString,RooRealVar*> listParamPostFit_ZZSR;

  for(int isize = 0; isize < parlist_init.getSize(); isize++){
    TString name (parlist_init.at(isize)->GetName());
    bool accepted = false;
  
    //// experimental uncertainties
    if(name.Contains("CMS_eff") or 
       name.Contains("CMS_fake") or
       name.Contains("CMS_reco") or
       name.Contains("CMS_trigger") or
       name.Contains("CMS_scale") or
       name.Contains("CMS_BDT_scale") or
       name.Contains("CMS_met") or
       name.Contains("FakeE") or
       name.Contains("FakeM") or
       name.Contains("lep2016") or
       name.Contains("tau2016") or
       name.Contains("CMS_veto") or
       name.Contains("lumi") or
       name.Contains("QCD_") or
       name.Contains("CMS_pu")){

      if(not (name.Contains("Bin") or name.Contains("bin") or name.Contains("ZnunuWJets"))){ 
	accepted = true;
	listParamPreFit_experimental[name] = (RooRealVar*) parlist_init.at(isize);
      }
    }

    ///// theory uncertainties
    if(name.Contains("QCDscale") or
       name.Contains("Top_Reweight") or
       name.Contains("pdf_") or
       name.Contains("_pdf_") or
       name.Contains("Higgs") or
       name.Contains("UEPS") or
       name.Contains("GJets") or
       name.Contains("ZH") or
       name.Contains("ggH") or
       name.Contains("EMSyst2016") or
       name.Contains("ZZWZ_EWKCorr") or
       name.Contains("gZZCorr")){
      
      if(not plotSBFit and  // not interested in signal nuisances
	 (name.Contains("ggH") or
	  name.Contains("qqH") or
	  name.Contains("ZH") or
	  name.Contains("VH") or
	  name.Contains("Higgs"))) continue;

      if(plotSBFit and name.Contains("stat")) continue;
      if(plotSBFit and name.Contains("Bin")) continue;

      accepted = true;
      listParamPreFit_theory[name] = (RooRealVar*) parlist_init.at(isize);
    }

    ///////////
    if((category == Category::monojet or category == Category::combined) and 
       (name.Contains("ZnunuZmm_MJ") or name.Contains("ZnunuZee_MJ"))){
      accepted = true;
      listParamPreFit_ZoverZ[name] = (RooRealVar*) parlist_init.at(isize);
    }
    else if((category == Category::monoV or category == Category::combined) and 
	    (name.Contains("ZnunuZmm_MV") or name.Contains("ZnunuZee_MV"))){
      listParamPreFit_ZoverZ[name] = (RooRealVar*) parlist_init.at(isize);
      accepted = true;
    }
    else if((category == Category::VBF or category == Category::combined) and 
	    (name.Contains("ZnunuZmm_QCD_VBF") or name.Contains("ZnunuZee_QCD_VBF") or name.Contains("dimuon") or name.Contains("dielectron") )){
      accepted = true;
      listParamPreFit_ZoverZ[name] = (RooRealVar*) parlist_init.at(isize);
    }
    else if((category == Category::VBF or category == Category::combined) and 
	    (name.Contains("ZnunuZmm_EWK_VBF") or name.Contains("ZnunuZee_EWK_VBF"))){
      accepted = true;
      listParamPreFit_ZoverZ[name] = (RooRealVar*) parlist_init.at(isize);
    }
    if((category == Category::monojet or category == Category::combined) and 
       (name.Contains("WJetsWmn_MJ") or name.Contains("WJetsWen_MJ") or name.Contains("WJetsWln"))){
      accepted = true;
      listParamPreFit_WoverW[name] = (RooRealVar*) parlist_init.at(isize);
    }
    else if((category == Category::monoV or category == Category::combined) and 
	    (name.Contains("WJetsWmn_MV") or name.Contains("WJetsWen_MV") or name.Contains("WJetsWln"))){
      accepted = true;
      listParamPreFit_WoverW[name] = (RooRealVar*) parlist_init.at(isize);
    }
    else if((category == Category::VBF or category == Category::combined) and 
	    (name.Contains("WJetsWmn_QCD_VBF") or name.Contains("WJetsWen_QCD_VBF") or name.Contains("singlemuon") or name.Contains("singleelectron") )){
      listParamPreFit_WoverW[name] = (RooRealVar*) parlist_init.at(isize);
      accepted = true;
    }
    else if((category == Category::VBF or category == Category::combined) and 
	    (name.Contains("WJetsWmn_EWK_VBF") or name.Contains("WJetsWen_EWK_VBF"))){
      accepted = true;
      listParamPreFit_WoverW[name] = (RooRealVar*) parlist_init.at(isize);
    }
    
    if((category == Category::monojet or category == Category::combined) and 
       (name.Contains("ZnunuWJets_MJ") or name.Contains("Znunu_MJ")  or name.Contains("WJets_MJ"))){
      accepted = true;
      listParamPreFit_ZoverW[name] = (RooRealVar*) parlist_init.at(isize);
    }
    else if((category == Category::monoV or category == Category::combined) and 
	    (name.Contains("ZnunuWJets_MV") or name.Contains("Znunu_MJ")  or name.Contains("WJets_MJ"))){
      accepted = true;
      listParamPreFit_ZoverW[name] = (RooRealVar*) parlist_init.at(isize);
    }
    else if((category == Category::VBF or category == Category::combined) and 
	    (name.Contains("ZnunuWJets_QCD_VBF") or name.Contains("ZnunuWJets") or name.Contains("wzCR") or name.Contains("qcd_ewk"))){
      accepted = true;
      listParamPreFit_ZoverW[name] = (RooRealVar*) parlist_init.at(isize);
    }
    else if((category == Category::VBF or category == Category::combined) and 
	    (name.Contains("ZnunuWJets_EWK_VBF") or name.Contains("ewk_ewk"))){
      accepted = true;
      listParamPreFit_ZoverW[name] = (RooRealVar*) parlist_init.at(isize);
    }

    if((category == Category::monojet or category == Category::combined) and 
       (name.Contains("ZnunuGamma_MJ") or name.Contains("Znunu_MJ") or name.Contains("Gamma_MJ"))){
      accepted = true;
      listParamPreFit_ZoverG[name] = (RooRealVar*) parlist_init.at(isize);
    }
    else if((category == Category::monoV or category == Category::combined) and 
	    (name.Contains("ZnunuGamma_MV") or name.Contains("Znunu_MJ") or name.Contains("Gamma_MJ"))){
      accepted = true;
      listParamPreFit_ZoverG[name] = (RooRealVar*) parlist_init.at(isize);
    }    

    if((category == Category::monoZ or category == Category::combined) and name.Contains("CMS_wz3l3l")){
      accepted = true;
      listParamPreFit_WZ3l[name] = (RooRealVar*) parlist_init.at(isize);
    }
    if((category == Category::monoZ or category == Category::combined) and name.Contains("CMS_zllhinvllll1j")){
      accepted = true;
      listParamPreFit_ZZ4l[name] = (RooRealVar*) parlist_init.at(isize);
    }
    if((category == Category::monoZ or category == Category::combined) and name.Contains("CMS_zllhinvll1j")){
      accepted = true;
      listParamPreFit_ZZSR[name] = (RooRealVar*) parlist_init.at(isize);
    }

  }

  ////
  for(int isize = 0; isize < parlist_final.getSize(); isize++){
    TString name (parlist_final.at(isize)->GetName());
    bool accepted = false;

    //// experimental uncertainties
    if(name.Contains("CMS_eff") or 
       name.Contains("CMS_fake") or
       name.Contains("CMS_reco") or
       name.Contains("CMS_trigger") or
       name.Contains("CMS_scale") or
       name.Contains("CMS_BDT_scale") or
       name.Contains("CMS_met") or
       name.Contains("FakeE") or
       name.Contains("FakeM") or
       name.Contains("lep2016") or
       name.Contains("tau2016") or
       name.Contains("CMS_veto") or
       name.Contains("lumi") or
       name.Contains("QCD_") or
       name.Contains("CMS_pu")){

      if(not (name.Contains("Bin") or name.Contains("bin") or name.Contains("ZnunuWJets"))){
      accepted = true;
      listParamPostFit_experimental[name] = (RooRealVar*) parlist_final.at(isize);
      }
    }

    ///// theory uncertainties
    if(name.Contains("QCDscale") or
       name.Contains("Top_Reweight") or
       name.Contains("pdf_") or
       name.Contains("_pdf_") or
       name.Contains("Higgs") or
       name.Contains("UEPS") or
       name.Contains("GJets") or
       name.Contains("ZH") or
       name.Contains("ggH") or
       name.Contains("EMSyst2016") or
       name.Contains("ZZWZ_EWKCorr") or
       name.Contains("gZZCorr")){
      
      if(not plotSBFit and  // not interested in signal nuisances
	 (name.Contains("ggH") or
	  name.Contains("qqH") or
	  name.Contains("ZH") or
	  name.Contains("VH") or
	  name.Contains("Higgs"))) continue;

      if(plotSBFit and name.Contains("stat")) continue;
      if(plotSBFit and name.Contains("Bin")) continue;

      accepted = true;
      listParamPostFit_theory[name] = (RooRealVar*) parlist_final.at(isize);
    }

    ///////////
    if((category == Category::monojet or category == Category::combined) and 
       (name.Contains("ZnunuZmm_MJ") or name.Contains("ZnunuZee_MJ"))){
      accepted = true;
      listParamPostFit_ZoverZ[name] = (RooRealVar*) parlist_final.at(isize);
    }
    else if((category == Category::monoV or category == Category::combined) and 
	    (name.Contains("ZnunuZmm_MV") or name.Contains("ZnunuZee_MV"))){
      listParamPostFit_ZoverZ[name] = (RooRealVar*) parlist_final.at(isize);
      accepted = true;
    }
    else if((category == Category::VBF or category == Category::combined) and 
	    (name.Contains("ZnunuZmm_QCD_VBF") or name.Contains("ZnunuZee_QCD_VBF") or name.Contains("dimuon") or name.Contains("dielectron") )){
      accepted = true;
      listParamPostFit_ZoverZ[name] = (RooRealVar*) parlist_final.at(isize);
    }
    else if((category == Category::VBF or category == Category::combined) and 
	    (name.Contains("ZnunuZmm_EWK_VBF") or name.Contains("ZnunuZee_EWK_VBF"))){
      accepted = true;
      listParamPostFit_ZoverZ[name] = (RooRealVar*) parlist_final.at(isize);
    }
    if((category == Category::monojet or category == Category::combined) and 
       (name.Contains("WJetsWmn_MJ") or name.Contains("WJetsWen_MJ") or name.Contains("WJetsWln"))){
      accepted = true;
      listParamPostFit_WoverW[name] = (RooRealVar*) parlist_final.at(isize);
    }
    else if((category == Category::monoV or category == Category::combined) and 
	    (name.Contains("WJetsWmn_MV") or name.Contains("WJetsWen_MV") or name.Contains("WJetsWln"))){
      accepted = true;
      listParamPostFit_WoverW[name] = (RooRealVar*) parlist_final.at(isize);
    }
    else if((category == Category::VBF or category == Category::combined) and 
	    (name.Contains("WJetsWmn_QCD_VBF") or name.Contains("WJetsWen_QCD_VBF") or name.Contains("singlemuon") or name.Contains("singleelectron") )){
      listParamPostFit_WoverW[name] = (RooRealVar*) parlist_final.at(isize);
      accepted = true;
    }
    else if((category == Category::VBF or category == Category::combined) and 
	    (name.Contains("WJetsWmn_EWK_VBF") or name.Contains("WJetsWen_EWK_VBF"))){
      accepted = true;
      listParamPostFit_WoverW[name] = (RooRealVar*) parlist_final.at(isize);
    }
    
    if((category == Category::monojet or category == Category::combined) and 
       (name.Contains("ZnunuWJets_MJ") or name.Contains("Znunu_MJ")  or name.Contains("WJets_MJ"))){
      accepted = true;
      listParamPostFit_ZoverW[name] = (RooRealVar*) parlist_final.at(isize);
    }
    else if((category == Category::monoV or category == Category::combined) and 
	    (name.Contains("ZnunuWJets_MV") or name.Contains("Znunu_MJ")  or name.Contains("WJets_MJ"))){
      accepted = true;
      listParamPostFit_ZoverW[name] = (RooRealVar*) parlist_final.at(isize);
    }
    else if((category == Category::VBF or category == Category::combined) and 
	    (name.Contains("ZnunuWJets_QCD_VBF") or name.Contains("ZnunuWJets") or name.Contains("wzCR") or name.Contains("qcd_ewk"))){
      accepted = true;
      listParamPostFit_ZoverW[name] = (RooRealVar*) parlist_final.at(isize);
    }
    else if((category == Category::VBF or category == Category::combined) and 
	    (name.Contains("ZnunuWJets_EWK_VBF") or name.Contains("ewk_ewk"))){
      accepted = true;
      listParamPostFit_ZoverW[name] = (RooRealVar*) parlist_final.at(isize);
    }

    if((category == Category::monojet or category == Category::combined) and 
       (name.Contains("ZnunuGamma_MJ") or name.Contains("Znunu_MJ") or name.Contains("Gamma_MJ"))){
      accepted = true;
      listParamPostFit_ZoverG[name] = (RooRealVar*) parlist_final.at(isize);
    }
    else if((category == Category::monoV or category == Category::combined) and 
	    (name.Contains("ZnunuGamma_MV") or name.Contains("Znunu_MJ") or name.Contains("Gamma_MJ"))){
      accepted = true;
      listParamPostFit_ZoverG[name] = (RooRealVar*) parlist_final.at(isize);
    }    

    if((category == Category::monoZ or category == Category::combined) and name.Contains("CMS_wz3l3l")){
      accepted = true;
      listParamPostFit_WZ3l[name] = (RooRealVar*) parlist_final.at(isize);
    }
    if((category == Category::monoZ or category == Category::combined) and name.Contains("CMS_zllhinvllll1j")){
      accepted = true;
      listParamPostFit_ZZ4l[name] = (RooRealVar*) parlist_final.at(isize);
    }
    if((category == Category::monoZ or category == Category::combined) and name.Contains("CMS_zllhinvll1j")){
      accepted = true;
      listParamPostFit_ZZSR[name] = (RooRealVar*) parlist_final.at(isize);
    }
  }
  
  /////////////// ------ 
  TCanvas* canvas = new TCanvas ("canvas","canvas",900,600);
  canvas->SetBottomMargin(0.35);
  plotNuisance(canvas,listParamPreFit_experimental,listParamPostFit_experimental,outputDIR,"experimental");
  plotNuisance(canvas,listParamPreFit_theory,listParamPostFit_theory,outputDIR,"theory");

  delete canvas;
  canvas = new TCanvas ("canvas","canvas",1500,700);
  canvas->SetBottomMargin(0.35);

  if(category == Category::monojet or category == Category::monoV or category == Category::combined or category == Category::VBF)
    plotNuisance(canvas,listParamPreFit_ZoverZ,listParamPostFit_ZoverZ,outputDIR,"ZZ");
  if(category == Category::monojet or category == Category::monoV or category == Category::combined or category == Category::VBF)
    plotNuisance(canvas,listParamPreFit_WoverW,listParamPostFit_WoverW,outputDIR,"WW");
  if(category == Category::monojet or category == Category::monoV or category == Category::combined or category == Category::VBF)
    plotNuisance(canvas,listParamPreFit_ZoverW,listParamPostFit_ZoverW,outputDIR,"ZW");

  delete canvas;
  canvas = new TCanvas ("canvas","canvas",1100,600);
  canvas->SetBottomMargin(0.35);

  if(category == Category::monojet or category == Category::monoV or category == Category::combined)
    plotNuisance(canvas,listParamPreFit_ZoverG,listParamPostFit_ZoverG,outputDIR,"ZG");
  if(category == Category::monoZ or category == Category::combined)
    plotNuisance(canvas,listParamPreFit_WZ3l,listParamPostFit_WZ3l,outputDIR,"WZ3l");
  if(category == Category::monoZ or category == Category::combined)
    plotNuisance(canvas,listParamPreFit_ZZ4l,listParamPostFit_ZZ4l,outputDIR,"ZZ4l");
  if(category == Category::monoZ or category == Category::combined)
    plotNuisance(canvas,listParamPreFit_ZZSR,listParamPostFit_ZZSR,outputDIR,"ZZSR");
}

