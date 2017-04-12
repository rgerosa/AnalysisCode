#include "../CMS_lumi.h"

void drawPlot(TCanvas* canvas, TH1* histo1, TH1* histo2, string outputDIR, string yaxisTitle, string postfix,
	      string legend1 = "", string legend2 = ""){

  histo1->GetXaxis()->SetTitle("boson p_{T} [GeV]");
  histo1->GetYaxis()->SetTitle(yaxisTitle.c_str());
  histo1->GetXaxis()->SetTitleSize(0.045);
  histo1->GetYaxis()->SetTitleSize(0.045);
  histo1->GetXaxis()->SetTitleOffset(1.1);
  histo1->GetYaxis()->SetTitleOffset(1.1);
  
  histo1->SetLineColor(kRed);
  histo2->SetLineColor(kBlue);
  histo1->SetLineWidth(2);
  histo2->SetLineWidth(2);

  histo1->GetXaxis()->SetRangeUser(150,1250);

  histo1->GetYaxis()->SetRangeUser(min(histo1->GetMinimum(),histo2->GetMinimum())*0.9,max(histo1->GetMaximum(),histo2->GetMaximum())*1.1);
  histo1->Draw("hist");
  histo2->Draw("hist same");

  TLegend leg (0.2,0.8,0.4,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  if(legend1 == "")
    leg.AddEntry(histo1,"CMS","L");
  else
    leg.AddEntry(histo1,legend1.c_str(),"L");
  if(legend2 == "")
     leg.AddEntry(histo2,"Theorist","L");
  else
    leg.AddEntry(histo2,legend2.c_str(),"L");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");

}

void makeKFactorComparisonWithTheory(string outputDIR){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  TFile* zvv_theory = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors_theorist_v3/vvj.root");
  TFile* wln_theory = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors_theorist_v3/evj.root");
  TFile* zll_theory = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors_theorist_v3/eej.root");
  
  TFile* kfactor_cms = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root");

  TH1* zvv_kfactor_qcd = (TH1*) zvv_theory->Get("vvj_pTV_K_NLO");
  TH1* wln_kfactor_qcd = (TH1*) wln_theory->Get("evj_pTV_K_NLO");
  TH1* zll_kfactor_qcd = (TH1*) zll_theory->Get("eej_pTV_K_NLO");

  TH1* zvv_kfactor_ewk = (TH1*) ((TH1*) zvv_theory->Get("vvj_pTV_kappa_NLO_EW"))->Clone("zvv_kfactor_ewk");
  zvv_kfactor_ewk->Add((TH1*) zvv_theory->Get("vvj_pTV_kappa_NNLO_Sud"));
  for(int iBin =1; iBin < zvv_kfactor_ewk->GetNbinsX()+1; iBin++)
    zvv_kfactor_ewk->AddBinContent(iBin);

  TH1* wln_kfactor_ewk = (TH1*) ((TH1*) wln_theory->Get("evj_pTV_kappa_NLO_EW"))->Clone("wln_kfactor_ewk");
  wln_kfactor_ewk->Add((TH1*) wln_theory->Get("evj_pTV_kappa_NNLO_Sud"));
  for(int iBin =1; iBin < wln_kfactor_ewk->GetNbinsX()+1; iBin++)
    wln_kfactor_ewk->AddBinContent(iBin);

  TH1* zll_kfactor_ewk = (TH1*) ((TH1*) zll_theory->Get("eej_pTV_kappa_NLO_EW"))->Clone("zll_kfactor_ewk");
  zll_kfactor_ewk->Add((TH1*) zll_theory->Get("eej_pTV_kappa_NNLO_Sud"));
  for(int iBin =1; iBin < zll_kfactor_ewk->GetNbinsX()+1; iBin++)
    zll_kfactor_ewk->AddBinContent(iBin);

  TH1* zvv_kfactor_ewk_cms = (TH1*) kfactor_cms->Get("EWKcorr/Z");
  zvv_kfactor_ewk_cms->Divide((TH1*) kfactor_cms->Get("ZJets_012j_NLO/nominal"));

  TH1* wln_kfactor_ewk_cms = (TH1*) kfactor_cms->Get("EWKcorr/W");
  wln_kfactor_ewk_cms->Divide((TH1*) kfactor_cms->Get("WJets_012j_NLO/nominal"));

  TH1* zvv_kfactor_qcd_cms = (TH1*) kfactor_cms->Get("ZJets_012j_NLO/nominal");
  zvv_kfactor_qcd_cms->Divide((TH1*) kfactor_cms->Get("ZJets_LO/inv_pt"));

  TH1* wln_kfactor_qcd_cms = (TH1*) kfactor_cms->Get("WJets_012j_NLO/nominal");
  wln_kfactor_qcd_cms->Divide((TH1*) kfactor_cms->Get("WJets_LO/inv_pt"));

  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();

  // ---
  drawPlot(canvas,zvv_kfactor_qcd_cms,zvv_kfactor_qcd,outputDIR,"NLO QCD K-factor","kfactor_qcd_zvv");
  drawPlot(canvas,wln_kfactor_qcd_cms,wln_kfactor_qcd,outputDIR,"NLO QCD k-factor","kfactor_qcd_wln");
  drawPlot(canvas,zvv_kfactor_ewk_cms,zvv_kfactor_ewk,outputDIR,"NLO EWK K-factor","kfactor_ewk_zvv");
  drawPlot(canvas,wln_kfactor_ewk_cms,wln_kfactor_ewk,outputDIR,"NLO EWK k-factor","kfactor_ewk_wln");

  // zll vs zvv
  drawPlot(canvas,zll_kfactor_ewk,zvv_kfactor_ewk,outputDIR,"NLO EWK k-factor","kfactor_ewk_zll_zvv","Zvv","Zll");
  drawPlot(canvas,zll_kfactor_qcd,zvv_kfactor_qcd,outputDIR,"NLO QCD k-factor","kfactor_qcd_zll_zvv","Zvv","Zll");
  
}
