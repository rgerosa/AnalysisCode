#include "../CMS_lumi.h"

TH1* flipHisto(TH1* histo){

  histo_flip = (TH1*) histo->Clone(Form("%s_flip",histo->GetName()));
  for(int iBin=0; iBin < histo->GetNbinsX()+1; iBin++){
    histo_flip->SetBinContent(iBin,-histo->GetBinContent(iBin));
  }

  return histo_flip;

}

void drawPlot(TCanvas* canvas, TH1* histo1, TH1* histo2, string outputDIR, string yaxisTitle, string postfix){

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

  histo1->GetYaxis()->SetRangeUser(min(histo1->GetMinimum(),histo2->GetMinimum())*0.9,max(histo1->GetMaximum(),histo2->GetMaximum())*1.1);
  histo1->Draw("hist");
  histo2->Draw("hist same");

  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");

}

void makeUncertaintyComparisonWithTheory(string outputDIR){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  TFile* zwjewkfile_cms = TFile::Open("/afs/cern.ch/user/r/rgerosa/work/MONOJET_ANALYSIS/CMSSW_8_0_20_patch1/src/AnalysisCode/MonoXAnalysis/macros/makeTemplates/templates_monojet_moriond_AN/zwjcorewk.root");
  TFile* zwjqcdfile_cms = TFile::Open("/afs/cern.ch/user/r/rgerosa/work/MONOJET_ANALYSIS/CMSSW_8_0_20_patch1/src/AnalysisCode/MonoXAnalysis/macros/makeTemplates/templates_monojet_moriond_AN/zwjcorqcd.root");

  TH1* zvv_ewk_cms = (TH1*) zwjewkfile_cms->Get("nhist_zwj_ewk_met");
  TH1* zvv_qcd_cms = (TH1*) zwjqcdfile_cms->Get("nhist_zwj_qcd_met");
  TH1* wln_ewk_cms = (TH1*) zwjewkfile_cms->Get("dhist_zwj_ewk_met");
  TH1* wln_qcd_cms = (TH1*) zwjqcdfile_cms->Get("dhist_zwj_qcd_met");

  TH1* zvv_ewk_unc_cms = (TH1*) zvv_ewk_cms->Clone("zvv_ewk_unc_cms");
  for(int iBin = 0; iBin < zvv_ewk_unc_cms->GetNbinsX()+1; iBin++)
    zvv_ewk_unc_cms->SetBinContent(iBin,fabs(zvv_ewk_cms->GetBinContent(iBin)-zvv_qcd_cms->GetBinContent(iBin))/zvv_ewk_cms->GetBinContent(iBin));  
  TH1* zvv_ewk_unc_cms_flip = flipHisto(zvv_ewk_unc_cms);

  TH1* wln_ewk_unc_cms = (TH1*) wln_ewk_cms->Clone("wln_ewk_unc_cms");
  for(int iBin = 0; iBin < wln_ewk_unc_cms->GetNbinsX()+1; iBin++)
    wln_ewk_unc_cms->SetBinContent(iBin,fabs(wln_ewk_cms->GetBinContent(iBin)-wln_qcd_cms->GetBinContent(iBin))/wln_ewk_cms->GetBinContent(iBin));  
  TH1* wln_ewk_unc_cms_flip = flipHisto(wln_ewk_unc_cms);


  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();

  drawPlot(canvas,zvv_ewk_unc_cms,zvv_ewk_unc_cms_flip,outputDIR,"NLO EWK uncertainty","zvv_nloewk_unc_cms");
  drawPlot(canvas,wln_ewk_unc_cms,wln_ewk_unc_cms_flip,outputDIR,"NLO EWK uncertainty","wln_nloewk_unc_cms");

  TFile* zwjewkfile_nloewk_up = TFile::Open("../../macros/makeTemplates/templates_theorist_kfactors_test/zwjcornloewk_up.root");
  TFile* zwjewkfile_nloewk_dw = TFile::Open("../../macros/makeTemplates/templates_theorist_kfactors_test/zwjcornloewk_dw.root");

  TH1* zvv_nloewk_th_up = (TH1*) zwjewkfile_nloewk_up->Get("nhist_zwj_nloewk_up_met");
  TH1* zvv_nloewk_th_dw = (TH1*) zwjewkfile_nloewk_dw->Get("nhist_zwj_nloewk_dw_met");
  TH1* wln_nloewk_th_up = (TH1*) zwjewkfile_nloewk_up->Get("dhist_zwj_nloewk_up_met");
  TH1* wln_nloewk_th_dw = (TH1*) zwjewkfile_nloewk_dw->Get("dhist_zwj_nloewk_dw_met");
  
  TH1* zvv_nloewk_unc_th = (TH1*) zvv_ewk_cms->Clone("zvv_nloewk_unc_th");
  for(int iBin = 0; iBin < zvv_nloewk_unc_th->GetNbinsX()+1; iBin++)
    zvv_nloewk_unc_th->SetBinContent(iBin,fabs(zvv_nloewk_th_up->GetBinContent(iBin)-zvv_nloewk_th_dw->GetBinContent(iBin))/(2*zvv_ewk_cms->GetBinContent(iBin)));
  TH1* zvv_nloewk_unc_th_flip = flipHisto(zvv_nloewk_unc_th);

  TH1* wln_nloewk_unc_th = (TH1*) wln_ewk_cms->Clone("wln_nloewk_unc_th");
  for(int iBin = 0; iBin < wln_nloewk_unc_th->GetNbinsX()+1; iBin++)
    wln_nloewk_unc_th->SetBinContent(iBin,fabs(wln_nloewk_th_up->GetBinContent(iBin)-wln_nloewk_th_dw->GetBinContent(iBin))/(2*wln_ewk_cms->GetBinContent(iBin)));
  TH1* wln_nloewk_unc_th_flip = flipHisto(wln_nloewk_unc_th);

  drawPlot(canvas,zvv_nloewk_unc_th,zvv_nloewk_unc_th_flip,outputDIR,"NLO EWK uncertainty","zvv_nloewk_unc_th");
  drawPlot(canvas,wln_nloewk_unc_th,wln_nloewk_unc_th_flip,outputDIR,"NLO EWK uncertainty","wln_nloewk_unc_th");

  TFile* zwjewkfile_sudewk_up = TFile::Open("../../macros/makeTemplates/templates_theorist_kfactors_test/zwjcorsudewk_up.root");
  TFile* zwjewkfile_sudewk_dw = TFile::Open("../../macros/makeTemplates/templates_theorist_kfactors_test/zwjcorsudewk_dw.root");

  TH1* zvv_sudewk_th_up = (TH1*) zwjewkfile_sudewk_up->Get("nhist_zwj_sudewk_up_met");
  TH1* zvv_sudewk_th_dw = (TH1*) zwjewkfile_sudewk_dw->Get("nhist_zwj_sudewk_dw_met");
  TH1* wln_sudewk_th_up = (TH1*) zwjewkfile_sudewk_up->Get("dhist_zwj_sudewk_up_met");
  TH1* wln_sudewk_th_dw = (TH1*) zwjewkfile_sudewk_dw->Get("dhist_zwj_sudewk_dw_met");
  
  TH1* zvv_sudewk_unc_th = (TH1*) zvv_ewk_cms->Clone("zvv_sudewk_unc_th");
  for(int iBin = 0; iBin < zvv_sudewk_unc_th->GetNbinsX()+1; iBin++)
    zvv_sudewk_unc_th->SetBinContent(iBin,fabs(zvv_sudewk_th_up->GetBinContent(iBin)-zvv_sudewk_th_dw->GetBinContent(iBin))/(2*zvv_ewk_cms->GetBinContent(iBin)));
  TH1* zvv_sudewk_unc_th_flip = flipHisto(zvv_sudewk_unc_th);

  TH1* wln_sudewk_unc_th = (TH1*) wln_ewk_cms->Clone("wln_sudewk_unc_th");
  for(int iBin = 0; iBin < wln_sudewk_unc_th->GetNbinsX()+1; iBin++)
    wln_sudewk_unc_th->SetBinContent(iBin,fabs(wln_sudewk_th_up->GetBinContent(iBin)-wln_sudewk_th_dw->GetBinContent(iBin))/(2*wln_ewk_cms->GetBinContent(iBin)));
  TH1* wln_sudewk_unc_th_flip = flipHisto(wln_sudewk_unc_th);

  drawPlot(canvas,zvv_sudewk_unc_th,zvv_sudewk_unc_th_flip,outputDIR,"NNLL EWK uncertainty","zvv_sudewk_unc_th");
  drawPlot(canvas,wln_sudewk_unc_th,wln_sudewk_unc_th_flip,outputDIR,"NNLL EWK uncertainty","wln_sudewk_unc_th");

  TFile* zwjewkfile_qcdscale_up = TFile::Open("../../macros/makeTemplates/templates_theorist_kfactors_test/zwjcorqcd_scaleup.root");
  TFile* zwjewkfile_qcdscale_dw = TFile::Open("../../macros/makeTemplates/templates_theorist_kfactors_test/zwjcorqcd_scaledw.root");

  TH1* zvv_qcd_scale_th_up = (TH1*) zwjewkfile_qcdscale_up->Get("nhist_zwj_qcd_scaleup_met");
  TH1* zvv_qcd_scale_th_dw = (TH1*) zwjewkfile_qcdscale_dw->Get("nhist_zwj_qcd_scaledw_met");
  TH1* wln_qcd_scale_th_up = (TH1*) zwjewkfile_qcdscale_up->Get("dhist_zwj_qcd_scaleup_met");
  TH1* wln_qcd_scale_th_dw = (TH1*) zwjewkfile_qcdscale_dw->Get("dhist_zwj_qcd_scaledw_met");
  
  TH1* zvv_qcd_scale_unc_th = (TH1*) zvv_ewk_cms->Clone("zvv_qcd_scale_unc_th");
  for(int iBin = 0; iBin < zvv_qcd_scale_unc_th->GetNbinsX()+1; iBin++)
    zvv_qcd_scale_unc_th->SetBinContent(iBin,fabs(zvv_qcd_scale_th_up->GetBinContent(iBin)-zvv_qcd_scale_th_dw->GetBinContent(iBin))/(2*zvv_ewk_cms->GetBinContent(iBin)));
  TH1* zvv_qcd_scale_unc_th_flip = flipHisto(zvv_qcd_scale_unc_th);

  TH1* wln_qcd_scale_unc_th = (TH1*) wln_ewk_cms->Clone("wln_qcd_scale_unc_th");
  for(int iBin = 0; iBin < wln_qcd_scale_unc_th->GetNbinsX()+1; iBin++)
    wln_qcd_scale_unc_th->SetBinContent(iBin,fabs(wln_qcd_scale_th_up->GetBinContent(iBin)-wln_qcd_scale_th_dw->GetBinContent(iBin))/(2*wln_ewk_cms->GetBinContent(iBin)));
  TH1* wln_qcd_scale_unc_th_flip = flipHisto(wln_qcd_scale_unc_th);

  drawPlot(canvas,zvv_qcd_scale_unc_th,zvv_qcd_scale_unc_th_flip,outputDIR,"QCD-scale uncertainty","zvv_qcd_scale_unc_th");
  drawPlot(canvas,wln_qcd_scale_unc_th,wln_qcd_scale_unc_th_flip,outputDIR,"QCD-scale uncertainty","wln_qcd_scale_unc_th");

}
