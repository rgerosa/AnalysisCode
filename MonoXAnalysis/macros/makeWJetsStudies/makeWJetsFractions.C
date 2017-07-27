#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static float luminosity = 35.9;
static int reductionFactor = 1;

void plotFraction(TH1* histo_total, TH1* histo_tau, TH1* histo_mu, TH1* histo_el, const string & outputDIR, const string & observable){

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);

  //                                                                                                                                                                                                  
  TH1* ratio_tau =  (TH1*) histo_tau->Clone("ratio_tau");
  TH1* ratio_mu  =  (TH1*) histo_mu->Clone("ratio_mu");
  TH1* ratio_el  =  (TH1*) histo_el->Clone("ratio_el");

  ratio_tau->Divide(histo_total);
  ratio_mu->Divide(histo_total);
  ratio_el->Divide(histo_total);

  histo_total->GetXaxis()->SetLabelSize(0);
  histo_total->GetXaxis()->SetTitleSize(0);
  
  if(observable == "mjj")
    ratio_tau->GetXaxis()->SetTitle("M_{jj} [GeV]");
  else if(observable == "met")
    ratio_tau->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");

  histo_total->GetYaxis()->SetRangeUser(min(histo_total->GetMinimum(),min(histo_tau->GetMinimum(),min(histo_mu->GetMinimum(),histo_el->GetMinimum())))*0.1,
					max(histo_total->GetMaximum(),max(histo_tau->GetMaximum(),max(histo_mu->GetMaximum(),histo_el->GetMaximum())))*100);

  histo_total->GetYaxis()->SetTitle("Events");
  histo_total->GetYaxis()->SetTitleOffset(1.1);
  histo_total->SetLineColor(kBlack);
  histo_total->SetLineWidth(2);
  histo_mu->SetLineColor(kBlue);
  histo_mu->SetLineWidth(2);
  histo_el->SetLineColor(kRed);
  histo_el->SetLineWidth(2);
  histo_tau->SetLineColor(kGreen+1);
  histo_tau->SetLineWidth(2);

  histo_total->Draw("hist");
  histo_mu->Draw("hist same");
  histo_el->Draw("hist same");
  histo_tau->Draw("hist same");

  ///                                                                                                                                                                                                  
  CMS_lumi(canvas,Form("%.1f",luminosity));
  
  TLegend* leg = new TLegend(0.55,0.6,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histo_total,"W+jets inclusive","L");
  leg->AddEntry(histo_mu,"W #rightarrow #mu#nu","L");
  leg->AddEntry(histo_el,"W #rightarrow e#nu","L");
  leg->AddEntry(histo_tau,"W #rightarrow #tau#nu","L");
  leg->Draw("same");

  canvas->SetLogy();

  pad2->Draw();
  pad2->cd();
  ratio_tau->GetYaxis()->SetTitle("Ratio");
  ratio_tau->GetYaxis()->SetTitleOffset(1.20);
  ratio_tau->GetYaxis()->SetTitleSize(0.04);
  ratio_tau->GetYaxis()->SetLabelSize(0.03);
  ratio_tau->GetYaxis()->SetNdivisions(505);
  ratio_tau->GetXaxis()->SetTitleOffset(1.10);
  ratio_tau->GetXaxis()->SetNdivisions(505);
  ratio_tau->SetLineColor(kGreen+1);
  ratio_tau->SetLineWidth(2);
  ratio_mu->SetLineColor(kBlue);
  ratio_mu->SetLineWidth(2);
  ratio_el->SetLineColor(kRed);
  ratio_el->SetLineWidth(2);
  ratio_tau->Draw("hist");
  ratio_mu->Draw("hist same");
  ratio_el->Draw("hist same");
  ratio_tau->GetYaxis()->SetRangeUser(min(ratio_tau->GetMinimum(),min(ratio_mu->GetMinimum(),ratio_el->GetMinimum()))*0.8,
  				      max(ratio_tau->GetMaximum(),max(ratio_mu->GetMaximum(),ratio_el->GetMaximum()))*1.2);

  canvas->SaveAs((outputDIR+"/"+"wjet_fraction_"+observable+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+"wjet_fraction_"+observable+".pdf").c_str(),"pdf");

  if(canvas) delete canvas;


}

void plotAcceptance(TH1* histo1, TH1* histo2, TH1* histo3, const string & outputDIR, const string & plot, const string & observable){

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);

  //                                                                                                                                                                                                  
  TH1* ratio_1 =  (TH1*) histo2->Clone("ratio_1");
  TH1* ratio_2 =  (TH1*) histo3->Clone("ratio_2");
  ratio_1->Divide(histo1);
  ratio_2->Divide(histo1);
  histo1->GetXaxis()->SetLabelSize(0);
  histo1->GetXaxis()->SetTitleSize(0);
  
  if(observable == "mjj")
    ratio_1->GetXaxis()->SetTitle("M_{jj} [GeV]");
  else if(observable == "met")
    ratio_1->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");

  histo1->GetYaxis()->SetRangeUser(min(histo1->GetMinimum(),min(histo2->GetMinimum(),histo3->GetMinimum()))*0.1,
				   max(histo1->GetMaximum(),max(histo2->GetMaximum(),histo3->GetMaximum()))*100);
  
  histo1->GetYaxis()->SetTitle("Events");
  histo1->GetYaxis()->SetTitleOffset(1.1);
  histo1->SetLineColor(kBlack);
  histo1->SetLineWidth(2);
  histo2->SetLineColor(kRed);
  histo2->SetLineWidth(2);
  histo2->SetLineStyle(7);
  histo3->SetLineColor(kBlue);
  histo3->SetLineWidth(2);
  histo3->SetLineStyle(7);

  histo1->Draw("hist");
  histo2->Draw("hist same");
  histo3->Draw("hist same");

  ///                                                                                                                                                                                                 
  CMS_lumi(canvas,Form("%.1f",luminosity));
  
  TLegend* leg = new TLegend(0.55,0.6,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  if(plot == "muonAcc"){
    leg->AddEntry(histo1,"W #rightarrow #mu#nu","L");
    leg->AddEntry(histo2,"W #rightarrow #mu#nu IN","L");
    leg->AddEntry(histo3,"W #rightarrow #mu#nu OUT","L");    
  }
  else if(plot == "eleAcc"){
    leg->AddEntry(histo1,"W #rightarrow e#nu","L");
    leg->AddEntry(histo2,"W #rightarrow e#nu IN","L");
    leg->AddEntry(histo3,"W #rightarrow e#nu OUT","L");    
  }
  else if(plot == "tauAcc"){
    leg->AddEntry(histo1,"W #rightarrow #tau#nu","L");
    leg->AddEntry(histo2,"W #rightarrow #tau#nu IN","L");
    leg->AddEntry(histo3,"W #rightarrow #tau#nu OUT","L");    
  }
  leg->Draw("same");

  canvas->SetLogy();

  pad2->Draw();
  pad2->cd();
  ratio_1->GetYaxis()->SetTitle("Ratio");
  ratio_1->GetYaxis()->SetTitleOffset(1.20);
  ratio_1->GetYaxis()->SetTitleSize(0.04);
  ratio_1->GetYaxis()->SetLabelSize(0.03);
  ratio_1->GetYaxis()->SetNdivisions(505);
  ratio_1->GetXaxis()->SetTitleOffset(1.10);
  ratio_1->GetXaxis()->SetNdivisions(505);
  ratio_1->SetLineColor(kRed);
  ratio_1->SetLineWidth(2);
  ratio_1->SetLineStyle(7);
  ratio_2->SetLineColor(kBlue);
  ratio_2->SetLineWidth(2);
  ratio_2->SetLineStyle(7);
  ratio_1->Draw("hist");
  ratio_2->Draw("hist same");
  ratio_1->GetYaxis()->SetRangeUser(min(ratio_1->GetMinimum(),ratio_2->GetMinimum())*0.8,
				    max(ratio_1->GetMaximum(),ratio_2->GetMaximum())*1.2);
  canvas->SaveAs((outputDIR+"/"+plot+"_"+observable+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+plot+"_"+observable+".pdf").c_str(),"pdf");

  if(canvas) delete canvas;

}



void plotUncertainty(TH1* histo_nominal, TH1* histo_up, TH1* histo_dw, const string & outputDIR, const string & plot, const string & observable){

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);

  //                                                                                                                                                                                                  
  TH1* ratio_up =  (TH1*) histo_up->Clone("ratio_up");
  TH1* ratio_dw =  (TH1*) histo_dw->Clone("ratio_dw");
  ratio_up->Divide(histo_nominal);
  ratio_dw->Divide(histo_nominal);
  histo_nominal->GetXaxis()->SetLabelSize(0);
  histo_nominal->GetXaxis()->SetTitleSize(0);
  
  if(observable == "mjj")
    ratio_up->GetXaxis()->SetTitle("M_{jj} [GeV]");
  else if(observable == "met")
    ratio_up->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");

  histo_nominal->GetYaxis()->SetRangeUser(min(histo_nominal->GetMinimum(),min(histo_up->GetMinimum(),histo_dw->GetMinimum()))*0.1,
					  max(histo_nominal->GetMaximum(),max(histo_up->GetMaximum(),histo_dw->GetMaximum()))*100);
  
  histo_nominal->GetYaxis()->SetTitle("Events");
  histo_nominal->GetYaxis()->SetTitleOffset(1.1);
  histo_nominal->SetLineColor(kBlack);
  histo_nominal->SetLineWidth(2);
  histo_up->SetLineColor(kBlack);
  histo_up->SetLineWidth(2);
  histo_up->SetLineStyle(7);
  histo_dw->SetLineColor(kBlack);
  histo_dw->SetLineWidth(2);
  histo_dw->SetLineStyle(7);
  
  histo_nominal->Draw("hist");
  histo_up->Draw("hist same");
  histo_dw->Draw("hist same");

  ///                                                                                                                                                                                                 
  CMS_lumi(canvas,Form("%.1f",luminosity));
  
  TLegend* leg = new TLegend(0.55,0.6,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  if(TString(plot).Contains("muon")){
    if(not TString(plot).Contains("total"))
      leg->AddEntry(histo_nominal,"W #rightarrow #mu#nu nominal","L");
    else
      leg->AddEntry(histo_nominal,"W+jets nominal","L");
    leg->AddEntry(histo_up,"W #rightarrow #mu#nu uncertainty","L");
  }
  if(TString(plot).Contains("ele")){
    if(not TString(plot).Contains("total"))
      leg->AddEntry(histo_nominal,"W #rightarrow e#nu nominal","L");
    else
      leg->AddEntry(histo_nominal,"W+jets nominal","L");
    leg->AddEntry(histo_up,"W #rightarrow e#nu uncertainty","L");
  }
  if(TString(plot).Contains("tau")){
    if(not TString(plot).Contains("total"))
      leg->AddEntry(histo_nominal,"W #rightarrow #tau#nu nominal","L");
    else
      leg->AddEntry(histo_nominal,"W+jets nominal","L");
    leg->AddEntry(histo_up,"W #rightarrow #tau#nu uncertainty","L");
  }
  
  leg->Draw("same");  
  canvas->SetLogy();

  pad2->Draw();
  pad2->cd();
  ratio_up->GetYaxis()->SetTitle("Ratio");
  ratio_up->GetYaxis()->SetTitleOffset(1.20);
  ratio_up->GetYaxis()->SetTitleSize(0.04);
  ratio_up->GetYaxis()->SetLabelSize(0.03);
  ratio_up->GetYaxis()->SetNdivisions(505);
  ratio_up->GetXaxis()->SetTitleOffset(1.10);
  ratio_up->GetXaxis()->SetNdivisions(505);
  ratio_up->SetLineColor(kBlack);
  ratio_up->SetLineWidth(2);
  ratio_up->SetLineStyle(7);
  ratio_dw->SetLineColor(kBlack);
  ratio_dw->SetLineWidth(2);
  ratio_dw->SetLineStyle(7);
  ratio_up->Draw("hist");
  ratio_dw->Draw("hist same");
  if(not TString(plot).Contains("total"))
    ratio_up->GetYaxis()->SetRangeUser(min(ratio_up->GetMinimum(),ratio_dw->GetMinimum())*0.9,
				       max(ratio_up->GetMaximum(),ratio_dw->GetMaximum())*1.1);
  else
    ratio_up->GetYaxis()->SetRangeUser(min(ratio_up->GetMinimum(),ratio_dw->GetMinimum())*0.95,
				       max(ratio_up->GetMaximum(),ratio_dw->GetMaximum())*1.05);

  canvas->SaveAs((outputDIR+"/"+plot+"_"+observable+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+plot+"_"+observable+".pdf").c_str(),"pdf");
  
  if(canvas) delete canvas;

}

static float sfel_unc  = 0.015;
static float sftau_unc = 0.06;
static float sftau_val = 0.99;

void makeWJetsFractions(string inputDIR, string outputDIR, string kfactorFile, Category category){

  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  initializeBinning();

  // k-factors
  cout<<"Load k-factors"<<endl;
  TFile* kffile = TFile::Open(kfactorFile.c_str());
  TH1* hist_nloqcdewk = NULL;
  TH1* hist_nloqcd    = NULL;
  TH1* hist_loqcd     = NULL;

  hist_nloqcdewk = (TH1*) kffile->Get("EWKcorr/W");
  hist_nloqcd    = (TH1*) kffile->Get("WJets_012j_NLO/nominal");
  hist_loqcd     = (TH1*) kffile->Get("WJets_LO/inv_pt");

  hist_nloqcdewk->Divide(hist_nloqcd);
  hist_nloqcd->Divide(hist_loqcd);

  TFile* kfactwjet_vbf = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_VBF_wjets_v2.root");
  ///////////////////                                                                                                                                                                                
  TH1* wjet_nlo_vbf = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_vbf");
  if(category == Category::VBFrelaxed)
    wjet_nlo_vbf = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_vbf_relaxed");
  if(category == Category::VBF)
    wjet_nlo_vbf->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_vbf"));
  else if(category == Category::VBFrelaxed)
    wjet_nlo_vbf->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_vbf_relaxed"));

  TH1* wjet_nlo_mj  = (TH1*) kfactwjet_vbf->Get("bosonPt_NLO_monojet");
  wjet_nlo_mj->Divide((TH1*) kfactwjet_vbf->Get("bosonPt_LO_monojet"));

  wjet_nlo_vbf->Divide(wjet_nlo_mj);

  // VBF k-factor
  vector<TH1*> khists; khists.push_back(hist_nloqcd); khists.push_back(hist_nloqcdewk); khists.push_back(wjet_nlo_vbf);

  // histograms
  cout<<"Book histograms"<<endl;
  vector<double> mjj_bin;
  vector<double> met_bin;
  if(category == Category::VBFrelaxed){
    mjj_bin = selectBinning("mjj",category);
    met_bin = selectBinning("met",category);
  }
  else if(category == Category::VBF){
    mjj_bin = {1300,2000,3500};
    met_bin = {250,350,500,1000};
  }

  TH1F* histo_mjj_total = new TH1F("histo_mjj_total","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_total_muon_up = new TH1F("histo_mjj_total_muon_up","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_total_muon_dw = new TH1F("histo_mjj_total_muon_dw","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_total_ele_up = new TH1F("histo_mjj_total_ele_up","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_total_ele_dw = new TH1F("histo_mjj_total_ele_dw","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_total_tau_up = new TH1F("histo_mjj_total_tau_up","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_total_tau_dw = new TH1F("histo_mjj_total_tau_dw","",mjj_bin.size()-1,&mjj_bin[0]);

  TH1F* histo_mjj_muon  = new TH1F("histo_mjj_muon","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_muon_up  = new TH1F("histo_mjj_muon_up","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_muon_dw  = new TH1F("histo_mjj_muon_dw","",mjj_bin.size()-1,&mjj_bin[0]);

  TH1F* histo_mjj_ele   = new TH1F("histo_mjj_ele","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_ele_up   = new TH1F("histo_mjj_ele_up","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_ele_dw   = new TH1F("histo_mjj_ele_dw","",mjj_bin.size()-1,&mjj_bin[0]);

  TH1F* histo_mjj_tau   = new TH1F("histo_mjj_tau","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_tau_up   = new TH1F("histo_mjj_tau_up","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_tau_dw   = new TH1F("histo_mjj_tau_dw","",mjj_bin.size()-1,&mjj_bin[0]);

  TH1F* histo_mjj_muon_inaccept  = new TH1F("histo_mjj_muon_inaccept","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_muon_outaccept = new TH1F("histo_mjj_muon_outaccept","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_ele_inaccept  = new TH1F("histo_mjj_ele_inaccept","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_ele_outaccept = new TH1F("histo_mjj_ele_outaccept","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_tau_inaccept  = new TH1F("histo_mjj_tau_inaccept","",mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* histo_mjj_tau_outaccept = new TH1F("histo_mjj_tau_outaccept","",mjj_bin.size()-1,&mjj_bin[0]);

  histo_mjj_total->Sumw2();
  histo_mjj_total_muon_up->Sumw2();
  histo_mjj_total_muon_dw->Sumw2();
  histo_mjj_total_ele_up->Sumw2();
  histo_mjj_total_ele_dw->Sumw2();
  histo_mjj_total_tau_up->Sumw2();
  histo_mjj_total_tau_dw->Sumw2();

  histo_mjj_muon->Sumw2();
  histo_mjj_muon_up->Sumw2();
  histo_mjj_muon_dw->Sumw2();

  histo_mjj_ele->Sumw2();
  histo_mjj_ele_up->Sumw2();
  histo_mjj_ele_dw->Sumw2();

  histo_mjj_tau->Sumw2();
  histo_mjj_tau_up->Sumw2();
  histo_mjj_tau_dw->Sumw2();

  histo_mjj_muon_inaccept->Sumw2();
  histo_mjj_ele_inaccept->Sumw2();
  histo_mjj_tau_inaccept->Sumw2();

  histo_mjj_muon_outaccept->Sumw2();
  histo_mjj_ele_outaccept->Sumw2();
  histo_mjj_tau_outaccept->Sumw2();

  //////////-
  TH1F* histo_met_total = new TH1F("histo_met_total","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_total_muon_up = new TH1F("histo_met_total_muon_up","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_total_muon_dw = new TH1F("histo_met_total_muon_dw","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_total_ele_up = new TH1F("histo_met_total_ele_up","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_total_ele_dw = new TH1F("histo_met_total_ele_dw","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_total_tau_up = new TH1F("histo_met_total_tau_up","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_total_tau_dw = new TH1F("histo_met_total_tau_dw","",met_bin.size()-1,&met_bin[0]);

  TH1F* histo_met_muon  = new TH1F("histo_met_muon","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_muon_up  = new TH1F("histo_met_muon_up","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_muon_dw  = new TH1F("histo_met_muon_dw","",met_bin.size()-1,&met_bin[0]);

  TH1F* histo_met_ele   = new TH1F("histo_met_ele","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_ele_up   = new TH1F("histo_met_ele_up","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_ele_dw   = new TH1F("histo_met_ele_dw","",met_bin.size()-1,&met_bin[0]);

  TH1F* histo_met_tau   = new TH1F("histo_met_tau","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_tau_up   = new TH1F("histo_met_tau_up","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_tau_dw   = new TH1F("histo_met_tau_dw","",met_bin.size()-1,&met_bin[0]);

  TH1F* histo_met_muon_inaccept  = new TH1F("histo_met_muon_inaccept","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_muon_outaccept = new TH1F("histo_met_muon_outaccept","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_ele_inaccept  = new TH1F("histo_met_ele_inaccept","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_ele_outaccept = new TH1F("histo_met_ele_outaccept","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_tau_inaccept  = new TH1F("histo_met_tau_inaccept","",met_bin.size()-1,&met_bin[0]);
  TH1F* histo_met_tau_outaccept = new TH1F("histo_met_tau_outaccept","",met_bin.size()-1,&met_bin[0]);

  histo_met_total->Sumw2();
  histo_met_total_muon_up->Sumw2();
  histo_met_total_muon_dw->Sumw2();
  histo_met_total_ele_up->Sumw2();
  histo_met_total_ele_dw->Sumw2();
  histo_met_total_tau_up->Sumw2();
  histo_met_total_tau_dw->Sumw2();

  histo_met_muon->Sumw2();
  histo_met_muon_up->Sumw2();
  histo_met_muon_dw->Sumw2();

  histo_met_ele->Sumw2();
  histo_met_ele_up->Sumw2();
  histo_met_ele_dw->Sumw2();

  histo_met_tau->Sumw2();
  histo_met_tau_up->Sumw2();
  histo_met_tau_dw->Sumw2();

  histo_met_muon_inaccept->Sumw2();
  histo_met_ele_inaccept->Sumw2();
  histo_met_tau_inaccept->Sumw2();

  histo_met_muon_outaccept->Sumw2();
  histo_met_ele_outaccept->Sumw2();
  histo_met_tau_outaccept->Sumw2();


  
  // tree reader
  TChain* tree_wjet = new TChain("tree/tree");
  tree_wjet->Add((inputDIR+"/*root").c_str());

  TTreeReader reader (tree_wjet);
  TTreeReaderValue<int>    putrue  (reader,"putrue");
  TTreeReaderValue<float>  xsec    (reader,"xsec");
  TTreeReaderValue<float>  wgt     (reader,"wgt");
  TTreeReaderValue<double> wgtsum  (reader,"wgtsum");
  TTreeReaderValue<float>  wgtpu   (reader,"wgtpileup");
  TTreeReaderValue<float>  wgtbtag (reader,"wgtbtag");
  TTreeReaderValue<UChar_t> hltm90      (reader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100     (reader,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110     (reader,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120     (reader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm90    (reader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hltmwm120   (reader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170   (reader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300   (reader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");
  TTreeReaderValue<unsigned int> nbjets     (reader,"nbjetslowpt");
  TTreeReaderValue<unsigned int> ntaus      (reader,"ntausold");
  TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
  TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
  TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
  TTreeReaderValue<vector<float> > chfrac  (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac  (reader,"combinejetNHfrac");
  TTreeReaderValue<unsigned int> nincjets  (reader,"njetsinc");
  TTreeReaderValue<float> met         (reader,"t1pfmet");
  TTreeReaderValue<float> metcalo     (reader,"calomet");
  TTreeReaderValue<float> jmmdphi (reader,"incjetmumetdphimin4");
  TTreeReaderValue<float> l1eta   (reader,"l1eta");
  TTreeReaderValue<float> l1pt    (reader,"l1pt");
  TTreeReaderValue<int>   l1id    (reader,"l1id");
  TTreeReaderValue<float> l2eta   (reader,"l2eta");
  TTreeReaderValue<float> l2pt    (reader,"l2pt");
  TTreeReaderValue<int>   l2id    (reader,"l2id");
  TTreeReaderValue<float> wzpt    (reader,"wzpt");
  


  // Met trigger efficiency                                                                                                                                                                         
  vector<TFile*> triggerfile_MET_binned_Wmn;
  vector<TF1*> triggermet_func_binned_Wmn;
  if(category == Category::VBF or category == Category::VBFrelaxed){ // monojet                                                                                    
    triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_0.0_800.0.root"));
    triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_800.0_1200.0.root"));
    triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_1200.0_1700.0.root"));
    triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_1700.0_3000.0.root"));
    if(triggerfile_MET_binned_Wmn.size() != 0){
      for(auto ifile : triggerfile_MET_binned_Wmn)
	triggermet_func_binned_Wmn.push_back((TF1*) ifile->Get("efficiency_func"));
    }
  }

  /// muon and electron loose scale factor
  TFile* muon_looseSF_file = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/muon_scalefactors.root");
  TFile* ele_vetoSF_file   = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/scalefactors_80x_egpog_37ifb.root");

  TH2F* msfloose_id = (TH2F*) muon_looseSF_file->Get("scalefactors_MuonLooseId_Muon");
  TH2F* msfloose_iso = (TH2F*) muon_looseSF_file->Get("scalefactors_Iso_MuonLooseId");
  TH2F* esfveto     = (TH2F*) ele_vetoSF_file->Get("scalefactors_Veto_Electron");

  /// muon, electron and tau efficiencies
  TFile* muon_efficiency_file = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonEfficiency_MC_2016/muon_efficiency_MC.root");
  TFile* ele_efficiency_file  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonEfficiency_MC_2016/egamma_efficiency_MC.root");
  TFile* tau_efficiency_file  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonEfficiency_MC_2016/efficiency_tau_MC.root");

  TH2F* effmuloose_id = (TH2F*) muon_efficiency_file->Get("efficiency_MC_loose_id");
  TH2F* effmuloose_iso = (TH2F*) muon_efficiency_file->Get("efficiency_MC_loose_iso");
  TH2F* effeleveto_id  = (TH2F*) ele_efficiency_file->Get("efficiency_MC_vetoid");
  TH2F* efftauloose_id =  (TH2F*) tau_efficiency_file->Get("efficiency_tau_MC");

  // loop on events                                                                                                                                                                                   
  cout<<"Event loop"<<endl;
  long int nEvents = tree_wjet->GetEntries();
  long int iEvent = 0;

  while(reader.Next()){

    cout.flush();
    if(iEvent %10000 == 0) cout<<"\r"<<"iEvent "<<100*double(iEvent)/(double(nEvents)/reductionFactor)<<" % ";
    if(iEvent > double(nEvents)/reductionFactor) break;
    iEvent++;

    // check trigger depending on the sample                                                                                                                                                          
    Double_t hlt = *hltm90+*hltm100+*hltm110+*hltm120+*hltmwm90+*hltmwm120+*hltmwm170+*hltmwm300;
    if(hlt == 0) continue;
    ///----
    if(*fhbhe == 0 or *fhbiso == 0 or *feeb == 0 or *fetp == 0 or *fvtx == 0 or *fcsc == 0 or *fbadmu == 0 or *fbadch == 0) continue;
    ///----
    if(*nbjets > 0) continue;
    if(*ntaus > 0) continue;
    ///----
    if(fabs(*met-*metcalo)/(*met) > 0.5) continue;
    if(*met < 250) continue;
    ///----
    if(*nincjets < 2) continue;
    if(jetpt->size() < 2) continue;
    if(jetpt->at(0) < 80) continue;
    if(jetpt->at(1) < 40) continue;
    if(fabs(jeteta->at(0)) > 4.7) continue;
    if(fabs(jeteta->at(1)) > 4.7) continue;
    if(jeteta->at(0)*jeteta->at(1) > 0) continue;
    ///----
    if(fabs(jeteta->at(0)) < 2.4 and chfrac->at(0) < 0.1) continue;
    if(fabs(jeteta->at(0)) < 2.4 and nhfrac->at(0) > 0.8) continue;
    if(*jmmdphi < 0.5) continue;

    TLorentzVector jet1,jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));

    if(category == Category::VBF){
      if(fabs(jeteta->at(0)-jeteta->at(1)) < 4) continue;
      if((jet1+jet2).M() < 1300) continue;
      if(fabs(jet1.DeltaPhi(jet2)) > 1.5) continue;
    }
    else if(category == Category::VBFrelaxed){
      if(fabs(jeteta->at(0)-jeteta->at(1)) < 1) continue;
      if((jet1+jet2).M() < 200) continue;
      if(fabs(jet1.DeltaPhi(jet2)) > 1.3) continue;
    }

    // met trigger scale factor                                                                                                                                                                    
    Double_t trig_wgt = 1.;
    if(category == Category::VBF or category == Category::twojet or category == Category::VBFrelaxed){
      double pfmet = *met;
      if((jet1+jet2).M() < 800)
	trig_wgt *= triggermet_func_binned_Wmn.at(0)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(0)->GetXaxis()->GetXmax()));
      else if((jet1+jet2).M() >= 800 and (jet1+jet2).M() < 1200)
	trig_wgt *= triggermet_func_binned_Wmn.at(1)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(1)->GetXaxis()->GetXmax()));
      else if((jet1+jet2).M() >= 1200 and (jet1+jet2).M() < 1700)
	trig_wgt *= triggermet_func_binned_Wmn.at(2)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(2)->GetXaxis()->GetXmax()));
      else if((jet1+jet2).M() >= 1700)
	trig_wgt *= triggermet_func_binned_Wmn.at(3)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(3)->GetXaxis()->GetXmax()));
    }

    //Gen level info --> NLO re-weight                                                                                                                                                                
    Double_t kwgt = 1.0;
    double genpt = *wzpt;
    for (size_t i = 0; i < khists.size(); i++) {
      if (khists[i]) {
        if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
        if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
        kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
      }
    }

    // pileup weight
    double puwgt = 1;
    if(fabs(*wgtpu) > 100)  puwgt = 1;
    else if(fabs(*wgtpu) < 0.01) puwgt = 1;
    else puwgt = *wgtpu;

    // lepton type event
    int leptonPDG = -1;
    bool inAcceptance = false;
    float ptLep   = 0;
    float etaLep  = 0;

    float sfwgt   = 1.; // scale factor for lepton veto
    float sfwgt_up  = 1.; // scale factor for lepton veto
    float sfwgt_dw  = 1.; // scale factor for lepton veto
    
      
    if(fabs(*l1id) == 13 or fabs(*l2id) == 13){ // muon event either from direct W-decay or leptonic tau-decay

      leptonPDG = 13;      
      
      if(fabs(*l1id) == 13){
	ptLep = *l1pt;
	etaLep = *l1eta;
      }
      else if(fabs(*l2id) == 13){
	ptLep = *l2pt;
	etaLep = *l2eta;
      }
      
      if(ptLep > 20 and fabs(etaLep) < 2.4){ // in acceptance --> read efficiency and SF
	inAcceptance = true;
	float ptVal = ptLep;
	if(ptVal > effmuloose_id->GetYaxis()->GetXmax()) ptVal = effmuloose_id->GetYaxis()->GetXmax()-1;
	if(ptVal < effmuloose_id->GetYaxis()->GetXmin()) ptVal = effmuloose_id->GetYaxis()->GetXmin()+1;
	float efficiency = effmuloose_id->GetBinContent(effmuloose_id->GetXaxis()->FindBin(fabs(etaLep)),effmuloose_id->GetYaxis()->FindBin(ptVal))*
	  effmuloose_iso->GetBinContent(effmuloose_iso->GetXaxis()->FindBin(fabs(etaLep)),effmuloose_iso->GetYaxis()->FindBin(ptVal));
	float sf = msfloose_id->GetBinContent(msfloose_id->GetXaxis()->FindBin(fabs(etaLep)),msfloose_id->GetYaxis()->FindBin(ptVal))*
	  msfloose_iso->GetBinContent(msfloose_iso->GetXaxis()->FindBin(fabs(etaLep)),msfloose_iso->GetYaxis()->FindBin(ptVal));	

	float sfmu_up = (msfloose_id->GetBinContent(msfloose_id->GetXaxis()->FindBin(fabs(etaLep)),msfloose_id->GetYaxis()->FindBin(ptVal)) + msfloose_id->GetBinError(msfloose_id->GetXaxis()->FindBin(fabs(etaLep)),msfloose_id->GetYaxis()->FindBin(ptVal)))*(msfloose_iso->GetBinContent(msfloose_iso->GetXaxis()->FindBin(fabs(etaLep)),msfloose_iso->GetYaxis()->FindBin(ptVal))+msfloose_iso->GetBinError(msfloose_iso->GetXaxis()->FindBin(fabs(etaLep)),msfloose_iso->GetYaxis()->FindBin(ptVal)));

	float sfmu_dw = (msfloose_id->GetBinContent(msfloose_id->GetXaxis()->FindBin(fabs(etaLep)),msfloose_id->GetYaxis()->FindBin(ptVal)) - msfloose_id->GetBinError(msfloose_id->GetXaxis()->FindBin(fabs(etaLep)),msfloose_id->GetYaxis()->FindBin(ptVal)))*(msfloose_iso->GetBinContent(msfloose_iso->GetXaxis()->FindBin(fabs(etaLep)),msfloose_iso->GetYaxis()->FindBin(ptVal))-msfloose_iso->GetBinError(msfloose_iso->GetXaxis()->FindBin(fabs(etaLep)),msfloose_iso->GetYaxis()->FindBin(ptVal)));

	if(efficiency == 0 or sf == 0) cout<<"muon issue: eta "<<etaLep<<" pt "<<ptVal<<endl;

	sfwgt *= (1-efficiency*sf)/(1-efficiency);	
	sfwgt_up *= (1-efficiency*sfmu_up)/(1-efficiency);
	sfwgt_dw *= (1-efficiency*sfmu_dw)/(1-efficiency);

      }
      else{
	sfwgt *= 1;
	sfwgt_up *= 1;
	sfwgt_dw *= 1;
      }
    }

    else if(fabs(*l1id) == 11 or fabs(*l2id) == 11){ // electron event either from direct W-decay or leptonic tau-decay 

      leptonPDG = 11;      
      if(fabs(*l1id) == 11){
	ptLep = *l1pt;
	etaLep = *l1eta;
      }
      else if(fabs(*l2id) == 11){
	ptLep = *l2pt;
	etaLep = *l2eta;
      }

      if(ptLep > 10 and fabs(etaLep) < 2.5){  // in acceptance --> read efficiency and SF   
	inAcceptance = true;
	float ptVal = ptLep; 
	if(ptVal > effeleveto_id->GetYaxis()->GetXmax()) ptVal = effeleveto_id->GetYaxis()->GetXmax()-1;
        if(ptVal < effeleveto_id->GetYaxis()->GetXmin()) ptVal = effeleveto_id->GetYaxis()->GetXmin()+1;
	float efficiency = effeleveto_id->GetBinContent(effeleveto_id->GetXaxis()->FindBin(etaLep),effeleveto_id->GetYaxis()->FindBin(ptVal));
	float sf = esfveto->GetBinContent(esfveto->GetXaxis()->FindBin(etaLep),esfveto->GetYaxis()->FindBin(ptVal));
	sfwgt *= (1-efficiency*sf)/(1-efficiency);
	sfwgt_up *= (1-efficiency*(sf+sfel_unc/2))/(1-efficiency);
	sfwgt_dw *= (1-efficiency*(sf-sfel_unc/2))/(1-efficiency);

	if(efficiency == 0 or sf == 0) cout<<"electron issue: eta "<<etaLep<<" pt "<<ptVal<<endl;

      }
      else{
	sfwgt *= 1;
	sfwgt_up *= 1;
	sfwgt_dw *= 1;
      }
    }

    else if(fabs(*l1id) == 15 or fabs(*l2id) == 15){ // hadronic taus

      leptonPDG = 15;      
      if(fabs(*l1id) == 15){
	ptLep = *l1pt;
	etaLep = *l1eta;
      }
      else if(fabs(*l2id) == 15){
	ptLep = *l2pt;
	etaLep = *l2eta;
      }

      if(ptLep > 18 and fabs(etaLep) < 2.3){  // in acceptance --> read efficiency and SF   
	inAcceptance = true;
	float ptVal = ptLep; 
	if(ptVal > efftauloose_id->GetXaxis()->GetXmax()) ptVal = efftauloose_id->GetXaxis()->GetXmax()-1;
        if(ptVal < efftauloose_id->GetXaxis()->GetXmin()) ptVal = efftauloose_id->GetXaxis()->GetXmin()+1;
	float efficiency = efftauloose_id->GetBinContent(efftauloose_id->GetXaxis()->FindBin(ptVal),efftauloose_id->GetYaxis()->FindBin(fabs(etaLep)));
	float sf = sftau_val;
	sfwgt *= (1-efficiency*sf)/(1-efficiency);
	sfwgt_up *= (1-efficiency*(sf+sftau_unc))/(1-efficiency);
	sfwgt_dw *= (1-efficiency*(sf-sftau_unc))/(1-efficiency);

	if(efficiency == 0 or sf == 0) cout<<"tau issue: eta "<<etaLep<<" pt "<<ptVal<<endl;

      }
      else{
	sfwgt *= 1;
	sfwgt_up *= 1;
	sfwgt_dw *= 1;
      }
    }
    
    // Fill histograms
    histo_mjj_total->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
    histo_met_total->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));

    if(leptonPDG == 11){
      histo_mjj_total_ele_up->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_up)/(*wgtsum));
      histo_mjj_total_ele_dw->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_dw)/(*wgtsum));
      histo_met_total_ele_up->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_up)/(*wgtsum));
      histo_met_total_ele_dw->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_dw)/(*wgtsum));
    }
    else{
      histo_mjj_total_ele_up->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_mjj_total_ele_dw->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_met_total_ele_up->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_met_total_ele_dw->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
    }

    if(leptonPDG == 13){
      histo_mjj_total_muon_up->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_up)/(*wgtsum));
      histo_mjj_total_muon_dw->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_dw)/(*wgtsum));
      histo_met_total_muon_up->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_up)/(*wgtsum));
      histo_met_total_muon_dw->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_dw)/(*wgtsum));
    }
    else{
      histo_mjj_total_muon_up->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_mjj_total_muon_dw->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_met_total_muon_up->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_met_total_muon_dw->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
    }

    if(leptonPDG == 15){
      histo_mjj_total_tau_up->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_up)/(*wgtsum));
      histo_mjj_total_tau_dw->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_dw)/(*wgtsum));
      histo_met_total_tau_up->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_up)/(*wgtsum));
      histo_met_total_tau_dw->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_dw)/(*wgtsum));
    }
    else{
      histo_mjj_total_tau_up->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_mjj_total_tau_dw->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_met_total_tau_up->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_met_total_tau_dw->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
    }

    if(leptonPDG == 11){// electron case
      histo_mjj_ele->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_mjj_ele_up->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_up)/(*wgtsum));
      histo_mjj_ele_dw->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_dw)/(*wgtsum));
      histo_met_ele->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_met_ele_up->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_up)/(*wgtsum));
      histo_met_ele_dw->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_dw)/(*wgtsum));
    }      
    else if(leptonPDG == 13){ //muon case
      histo_mjj_muon->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_mjj_muon_up->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_up)/(*wgtsum));
      histo_mjj_muon_dw->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_dw)/(*wgtsum));
      histo_met_muon->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_met_muon_up->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_up)/(*wgtsum));
      histo_met_muon_dw->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_dw)/(*wgtsum));
    }
      
    else if(leptonPDG == 15){ // hadronic tau
      histo_mjj_tau->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_mjj_tau_up->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_up)/(*wgtsum));
      histo_mjj_tau_dw->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_dw)/(*wgtsum));
      histo_met_tau->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum));
      histo_met_tau_up->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_up)/(*wgtsum));
      histo_met_tau_dw->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt_dw)/(*wgtsum));
    }
    
    if(leptonPDG == 11 and inAcceptance){      
      histo_mjj_ele_inaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum)); // IN of acceptance electrons	
      histo_met_ele_inaccept->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum)); // IN of acceptance electrons	
    }
    else if(leptonPDG == 11 and not inAcceptance){
      histo_mjj_ele_outaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum)); // OUT of acceptance electrons
      histo_met_ele_outaccept->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum)); // OUT of acceptance electrons
    }
    else if(leptonPDG == 13 and inAcceptance){      
      histo_mjj_muon_inaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum)); // IN of acceptance muons
      histo_met_muon_inaccept->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum)); // IN of acceptance muons
    }
    else if(leptonPDG == 13 and not inAcceptance){
      histo_mjj_muon_outaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum)); // OUT of acceptance muons
      histo_met_muon_outaccept->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum)); // OUT of acceptance muons
    }
    else if(leptonPDG == 15 and inAcceptance){      
      histo_mjj_tau_inaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum)); // IN of acceptance taus
      histo_met_tau_inaccept->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum)); // IN of acceptance taus
    }
    else if(leptonPDG == 15 and not inAcceptance){
      histo_mjj_tau_outaccept->Fill((jet1+jet2).M(),(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum)); // OUT of acceptance taus
      histo_met_tau_outaccept->Fill(*met,(*xsec)*luminosity*(*wgt)*(trig_wgt)*(kwgt)*(puwgt)*(*wgtbtag)*(sfwgt)/(*wgtsum)); // OUT of acceptance taus
    }
  }
  cout<<endl;

  //// ----- ////
  cout<<"######################"<<endl;
  cout<<"W+jets events in SR "<<histo_mjj_total->Integral()<<endl;
  cout<<"Wmn events in SR "<<histo_mjj_muon->Integral()<<" fraction "<<histo_mjj_muon->Integral()/histo_mjj_total->Integral()<<endl;
  cout<<"Wen events in SR "<<histo_mjj_ele->Integral()<<" fraction "<<histo_mjj_ele->Integral()/histo_mjj_total->Integral()<<endl;
  cout<<"Wtau events in SR "<<histo_mjj_tau->Integral()<<" fraction "<<histo_mjj_tau->Integral()/histo_mjj_total->Integral()<<endl;
  cout<<"######################"<<endl;

  cout<<"######################"<<endl;
  cout<<"Wmn events in SR: in acceptance "<<histo_mjj_muon_inaccept->Integral()<<"  out acceptance "<<histo_mjj_muon_outaccept->Integral()<<" fraction "<<histo_mjj_muon_inaccept->Integral()/(histo_mjj_muon_inaccept->Integral()+histo_mjj_muon_outaccept->Integral())<<endl;
  cout<<"Wen events in SR: in acceptance "<<histo_mjj_ele_inaccept->Integral()<<"  out acceptance "<<histo_mjj_ele_outaccept->Integral()<<" fraction "<<histo_mjj_ele_inaccept->Integral()/(histo_mjj_ele_inaccept->Integral()+histo_mjj_ele_outaccept->Integral())<<endl;
  cout<<"Wtau events in SR: in acceptance "<<histo_mjj_tau_inaccept->Integral()<<"  out acceptance "<<histo_mjj_tau_outaccept->Integral()<<" fraction "<<histo_mjj_tau_inaccept->Integral()/(histo_mjj_tau_inaccept->Integral()+histo_mjj_tau_outaccept->Integral())<<endl;
  cout<<"######################"<<endl;

  ////// ------------ //////
  TH1F* histo_mjj_muon_v2 = (TH1F*) histo_mjj_muon->Clone("histo_mjj_muon_v2");
  TH1F* histo_mjj_ele_v2 = (TH1F*) histo_mjj_ele->Clone("histo_mjj_ele_v2");
  TH1F* histo_mjj_tau_v2 = (TH1F*) histo_mjj_tau->Clone("histo_mjj_tau_v2");

  TH1F* histo_met_muon_v2 = (TH1F*) histo_met_muon->Clone("histo_met_muon_v2");
  TH1F* histo_met_ele_v2 = (TH1F*) histo_met_ele->Clone("histo_met_ele_v2");
  TH1F* histo_met_tau_v2 = (TH1F*) histo_met_tau->Clone("histo_met_tau_v2");

  plotFraction(histo_mjj_total,histo_mjj_tau,histo_mjj_muon,histo_mjj_ele,outputDIR,"mjj");
  plotAcceptance(histo_mjj_muon,histo_mjj_muon_inaccept,histo_mjj_muon_outaccept,outputDIR,"muonAcc","mjj");
  plotAcceptance(histo_mjj_ele,histo_mjj_ele_inaccept,histo_mjj_ele_outaccept,outputDIR,"eleAcc","mjj");
  plotAcceptance(histo_mjj_tau,histo_mjj_tau_inaccept,histo_mjj_tau_outaccept,outputDIR,"tauAcc","mjj");  

  plotFraction(histo_met_total,histo_met_tau,histo_met_muon,histo_met_ele,outputDIR,"met");
  plotAcceptance(histo_met_muon,histo_met_muon_inaccept,histo_met_muon_outaccept,outputDIR,"muonAcc","met");
  plotAcceptance(histo_met_ele,histo_met_ele_inaccept,histo_met_ele_outaccept,outputDIR,"eleAcc","met");
  plotAcceptance(histo_met_tau,histo_met_tau_inaccept,histo_met_tau_outaccept,outputDIR,"tauAcc","met");  
   
  ///// ---------- ///
  plotUncertainty(histo_mjj_muon_v2,histo_mjj_muon_up,histo_mjj_muon_dw,outputDIR,"muon_uncertainty","mjj");
  plotUncertainty(histo_mjj_ele_v2,histo_mjj_ele_up,histo_mjj_ele_dw,outputDIR,"electron_uncertainty","mjj");
  plotUncertainty(histo_mjj_tau_v2,histo_mjj_tau_up,histo_mjj_tau_dw,outputDIR,"tau_uncertainty","mjj");

  ///// ---------- ///
  plotUncertainty(histo_met_muon_v2,histo_met_muon_up,histo_met_muon_dw,outputDIR,"muon_uncertainty","mjt");
  plotUncertainty(histo_met_ele_v2,histo_met_ele_up,histo_met_ele_dw,outputDIR,"electron_uncertainty","met");
  plotUncertainty(histo_met_tau_v2,histo_met_tau_up,histo_met_tau_dw,outputDIR,"tau_uncertainty","met");

  ///// ---------- final plot for the uncertianty ///  
  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);

  //                                                                                                                                                                                                  
  TH1* ratio_tau_up =  (TH1*) histo_mjj_total_tau_up->Clone("ratio_tau_up");
  TH1* ratio_tau_dw =  (TH1*) histo_mjj_total_tau_dw->Clone("ratio_tau_dw");
  TH1* ratio_muon_up =  (TH1*) histo_mjj_total_muon_up->Clone("ratio_muon_up");
  TH1* ratio_muon_dw =  (TH1*) histo_mjj_total_muon_dw->Clone("ratio_muon_dw");
  TH1* ratio_ele_up =  (TH1*) histo_mjj_total_ele_up->Clone("ratio_ele_up");
  TH1* ratio_ele_dw =  (TH1*) histo_mjj_total_ele_dw->Clone("ratio_ele_dw");

  ratio_tau_up->Divide(histo_mjj_total);
  ratio_tau_dw->Divide(histo_mjj_total);
  ratio_muon_up->Divide(histo_mjj_total);
  ratio_muon_dw->Divide(histo_mjj_total);
  ratio_ele_up->Divide(histo_mjj_total);
  ratio_ele_dw->Divide(histo_mjj_total);

  histo_mjj_total->GetXaxis()->SetLabelSize(0);
  histo_mjj_total->GetXaxis()->SetTitleSize(0);
  
  histo_mjj_total->GetYaxis()->SetTitle("Events");
  histo_mjj_total->GetYaxis()->SetTitleOffset(1.1);
  histo_mjj_total->SetLineColor(kBlack);
  histo_mjj_total->SetLineWidth(2);
  
  histo_mjj_total_tau_up->SetLineColor(kGreen+1);
  histo_mjj_total_tau_up->SetLineWidth(2);
  histo_mjj_total_tau_up->SetLineStyle(7);
  histo_mjj_total_tau_dw->SetLineColor(kGreen+1);
  histo_mjj_total_tau_dw->SetLineWidth(2);
  histo_mjj_total_tau_dw->SetLineStyle(7);
  
  histo_mjj_total_muon_up->SetLineColor(kBlue);
  histo_mjj_total_muon_up->SetLineWidth(2);
  histo_mjj_total_muon_up->SetLineStyle(7);
  histo_mjj_total_muon_dw->SetLineColor(kBlue);
  histo_mjj_total_muon_dw->SetLineWidth(2);
  histo_mjj_total_muon_dw->SetLineStyle(7);

  histo_mjj_total_ele_up->SetLineColor(kRed);
  histo_mjj_total_ele_up->SetLineWidth(2);
  histo_mjj_total_ele_up->SetLineStyle(7);
  histo_mjj_total_ele_dw->SetLineColor(kRed);
  histo_mjj_total_ele_dw->SetLineWidth(2);
  histo_mjj_total_ele_dw->SetLineStyle(7);
  
  histo_mjj_total->Draw("hist");
  histo_mjj_total_ele_up->Draw("hist same");
  histo_mjj_total_ele_dw->Draw("hist same");
  histo_mjj_total_muon_up->Draw("hist same");
  histo_mjj_total_muon_dw->Draw("hist same");
  histo_mjj_total_tau_up->Draw("hist same");
  histo_mjj_total_tau_dw->Draw("hist same");
  histo_mjj_total->Draw("hist same");

  ///                                                                                                                                                                                                 
  CMS_lumi(canvas,Form("%.1f",luminosity));
  
  TLegend* leg = new TLegend(0.55,0.6,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histo_mjj_total,"W+jets SR","L");
  leg->AddEntry(histo_mjj_total_tau_up,"W+jets #tau-uncertainty","L");
  leg->AddEntry(histo_mjj_total_muon_up,"W+jets #mu-uncertainty","L");
  leg->AddEntry(histo_mjj_total_ele_up,"W+jets ele-uncertainty","L");
  leg->Draw("same");  
  canvas->SetLogy();

  pad2->Draw();
  pad2->cd();
  ratio_tau_up->GetXaxis()->SetTitle("M_{jj} [GeV]");
  ratio_tau_up->GetYaxis()->SetTitle("Ratio");
  ratio_tau_up->GetYaxis()->SetTitleOffset(1.20);
  ratio_tau_up->GetYaxis()->SetTitleSize(0.04);
  ratio_tau_up->GetYaxis()->SetLabelSize(0.03);
  ratio_tau_up->GetYaxis()->SetNdivisions(505);
  ratio_tau_up->GetXaxis()->SetTitleOffset(1.10);
  ratio_tau_up->GetXaxis()->SetNdivisions(505);

  ratio_tau_up->SetLineColor(kGreen+1);
  ratio_tau_up->SetLineWidth(2);
  ratio_tau_dw->SetLineColor(kGreen+1);
  ratio_tau_dw->SetLineWidth(2);

  ratio_muon_up->SetLineColor(kBlue);
  ratio_muon_up->SetLineWidth(2);
  ratio_muon_dw->SetLineColor(kBlue);
  ratio_muon_dw->SetLineWidth(2);

  ratio_ele_up->SetLineColor(kRed);
  ratio_ele_up->SetLineWidth(2);
  ratio_ele_dw->SetLineColor(kRed);
  ratio_ele_dw->SetLineWidth(2);

  ratio_tau_up->Draw("hist");
  ratio_tau_dw->Draw("hist same");
  ratio_muon_up->Draw("hist same");
  ratio_muon_dw->Draw("hist same");
  ratio_ele_up->Draw("hist same");
  ratio_ele_dw->Draw("hist same");

  ratio_tau_up->GetYaxis()->SetRangeUser(0.95,1.05);
  
  canvas->SaveAs((outputDIR+"/wjets_leptonveto_unc.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/wjets_leptonveto_unc.pdf").c_str(),"pdf");
  
  if(canvas) delete canvas;
}
