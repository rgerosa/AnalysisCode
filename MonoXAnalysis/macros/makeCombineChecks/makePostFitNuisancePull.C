#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

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
  
  map<TString,RooRealVar*> listParamPreFit_theory;
  map<TString,RooRealVar*> listParamPostFit_theory;
  map<TString,RooRealVar*> listParamPreFit_experimental;
  map<TString,RooRealVar*> listParamPostFit_experimental;

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
    if(name.Contains("QCDscale") or
       name.Contains("Top_Reweight") or
       name.Contains("pdf_") or
       name.Contains("QCDprocess") or
       name.Contains("QCDscale") or
       name.Contains("QCDshape") or
       name.Contains("sudakov") or
       name.Contains("nnlomiss") or
       name.Contains("qcdewkmix") or
       name.Contains("nnloewk") or
       name.Contains("_pdf_"))	  
      listParamPreFit_theory[name] = (RooRealVar*) parlist_init.at(isize);    

    if(name.Contains("CMS_eff") or
       name.Contains("CMS_fake") or
       name.Contains("CMS_reco") or
       name.Contains("CMS_trigger") or
       name.Contains("CMS_scale") or
       name.Contains("CMS_BDT_scale") or
       name.Contains("CMS_met") or
       name.Contains("CMS_veto") or
       name.Contains("lumi") or
       name.Contains("CMS_pu"))
      listParamPreFit_experimental[name] = (RooRealVar*) parlist_init.at(isize);

    if((category == Category::monojet and (name.Contains("ZnunuZmm_MJ") or name.Contains("ZnunuZee_MJ"))) or
       (category == Category::monoV and (name.Contains("ZnunuZmm_MV") or name.Contains("ZnunuZee_MV"))) or
       ((category == Category::VBF or category == Category::VBFrelaxed) and (name.Contains("ZnunuZmm_QCD_VBF") or name.Contains("ZnunuZee_QCD_VBF"))) or
       ((category == Category::VBF or category == Category::VBFrelaxed) and (name.Contains("ZnunuZmm_EWK_VBF") or name.Contains("ZnunuZee_EWK_VBF"))))      
      listParamPreFit_ZZ[name] = (RooRealVar*) parlist_init.at(isize);

    if((category == Category::monojet and (name.Contains("WJetsWmn_MJ") or name.Contains("WJetsWen_MJ"))) or
       (category == Category::monoV and (name.Contains("WJetsWmn_MV") or name.Contains("WJetsWen_MV"))) or
       ((category == Category::VBF or category == Category::VBFrelaxed) and (name.Contains("WJetsWmn_QCD_VBF") or name.Contains("WJetsWen_QCD_VBF"))) or
       ((category == Category::VBF or category == Category::VBFrelaxed) and (name.Contains("WJetsWmn_EWK_VBF") or name.Contains("WJetsWen_EWK_VBF"))))      
      listParamPreFit_WW[name] = (RooRealVar*) parlist_init.at(isize);
    

    if((category == Category::monojet and name.Contains("ZnunuGamma_MJ")) or
       (category == Category::monoV and name.Contains("ZnunuGamma_MV")))       
      listParamPreFit_ZG[name] = (RooRealVar*) parlist_init.at(isize);

    if((category == Category::monojet and name.Contains("ZnunuWJets_MJ")) or
       (category == Category::monoV and name.Contains("ZnunuWJets_MV")) or
       ((category == Category::VBF or category == Category::VBFrelaxed) and name.Contains("ZnunuWJets_QCD_VBF")) or
       ((category == Category::VBF or category == Category::VBFrelaxed) and name.Contains("ZnunuWJets_EWK_VBF")))
      listParamPreFit_ZW[name] = (RooRealVar*) parlist_init.at(isize);
    
  }
  
  ///////
  for(int isize = 0; isize < parlist_final.getSize(); isize++){
    TString name (parlist_final.at(isize)->GetName());

    if(name.Contains("QCDscale") or
       name.Contains("Top_Reweight") or
       name.Contains("pdf_") or
       name.Contains("QCDprocess") or
       name.Contains("QCDscale") or
       name.Contains("QCDshape") or
       name.Contains("sudakov") or
       name.Contains("nnlomiss") or
       name.Contains("qcdewkmix") or
       name.Contains("nnloewk") or
       name.Contains("_pdf_"))	  
      listParamPostFit_theory[name] = (RooRealVar*) parlist_final.at(isize);    

    if(name.Contains("CMS_eff") or
       name.Contains("CMS_fake") or
       name.Contains("CMS_reco") or
       name.Contains("CMS_trigger") or
       name.Contains("CMS_scale") or
       name.Contains("CMS_BDT_scale") or
       name.Contains("CMS_met") or
       name.Contains("CMS_veto") or
       name.Contains("lumi") or
       name.Contains("CMS_pu"))
      listParamPostFit_experimental[name] = (RooRealVar*) parlist_final.at(isize);
    
    if((category == Category::monojet and (name.Contains("ZnunuZmm_MJ") or name.Contains("ZnunuZee_MJ"))) or
       (category == Category::monoV and (name.Contains("ZnunuZmm_MV") or name.Contains("ZnunuZee_MV"))) or
       ((category == Category::VBF or category == Category::VBFrelaxed) and (name.Contains("ZnunuZmm_QCD_VBF") or name.Contains("ZnunuZee_QCD_VBF"))) or
       ((category == Category::VBF or category == Category::VBFrelaxed) and (name.Contains("ZnunuZmm_EWK_VBF") or name.Contains("ZnunuZee_EWK_VBF"))))      
      listParamPostFit_ZZ[name] = (RooRealVar*) parlist_final.at(isize);
    
    if((category == Category::monojet and (name.Contains("WJetsWmn_MJ") or name.Contains("WJetsWen_MJ"))) or
       (category == Category::monoV and (name.Contains("WJetsWmn_MV") or name.Contains("WJetsWen_MV"))) or
       ((category == Category::VBF or category == Category::VBFrelaxed) and (name.Contains("WJetsWmn_QCD_VBF") or name.Contains("WJetsWen_QCD_VBF"))) or
       ((category == Category::VBF or category == Category::VBFrelaxed) and (name.Contains("WJetsWmn_EWK_VBF") or name.Contains("WJetsWen_EWK_VBF"))))      
      listParamPostFit_WW[name] = (RooRealVar*) parlist_final.at(isize);
    
    
    if((category == Category::monojet and name.Contains("ZnunuGamma_MJ")) or
       (category == Category::monoV and name.Contains("ZnunuGamma_MV")))       
      listParamPostFit_ZG[name] = (RooRealVar*) parlist_final.at(isize);
    
    if((category == Category::monojet and name.Contains("ZnunuWJets_MJ")) or
       (category == Category::monoV and name.Contains("ZnunuWJets_MV")) or
       ((category == Category::VBF or category == Category::VBFrelaxed) and name.Contains("ZnunuWJets_QCD_VBF")) or
       ((category == Category::VBF or category == Category::VBFrelaxed) and name.Contains("ZnunuWJets_EWK_VBF")))
      listParamPostFit_ZW[name] = (RooRealVar*) parlist_final.at(isize);
  }
  
  /////////////// ------                                                                                                                                                                                
  TCanvas* canvas = new TCanvas ("canvas","canvas",900,600);
  canvas->SetBottomMargin(0.35);
  plotNuisance(canvas,listParamPreFit_experimental,listParamPostFit_experimental,outputDIR,"experimental");
  plotNuisance(canvas,listParamPreFit_theory,listParamPostFit_theory,outputDIR,"theory");

  delete canvas;
  canvas = new TCanvas ("canvas","canvas",1100,600);
  canvas->SetBottomMargin(0.35);

  if(category == Category::monojet or category == Category::monoV or category == Category::VBF or category == Category::VBFrelaxed)
    plotNuisance(canvas,listParamPreFit_ZZ,listParamPostFit_ZZ,outputDIR,"ZZ");
  if(category == Category::monojet or category == Category::monoV or category == Category::VBFrelaxed or category == Category::VBF)
    plotNuisance(canvas,listParamPreFit_WW,listParamPostFit_WW,outputDIR,"WW");
  if(category == Category::monojet or category == Category::monoV or category == Category::VBFrelaxed or category == Category::VBF)
    plotNuisance(canvas,listParamPreFit_ZW,listParamPostFit_ZW,outputDIR,"ZW");
  if(category == Category::monojet or category == Category::monoV)
    plotNuisance(canvas,listParamPreFit_ZG,listParamPostFit_ZG,outputDIR,"ZG");
       
}
			
