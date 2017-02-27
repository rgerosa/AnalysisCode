#include "../CMS_lumi.h"

TH1* flipHisto(TH1* histo){

  TH1* histo_flip = (TH1*) histo->Clone(Form("%s_flip",histo->GetName()));
  for(int iBin=0; iBin < histo->GetNbinsX()+1; iBin++){
    histo_flip->SetBinContent(iBin,-histo->GetBinContent(iBin));
  }

  return histo_flip;

}

void drawPlot(TCanvas* canvas, pair<TH1*,TH1*> zjet, pair<TH1*,TH1*> wjet, pair<TH1*,TH1*> gamma, string outputDIR, string yaxisTitle, string postfix){
  

  zjet.first->GetXaxis()->SetTitle("Recoil [GeV]");
  zjet.first->GetYaxis()->SetTitle(yaxisTitle.c_str());
  zjet.first->GetXaxis()->SetTitleSize(0.042);
  zjet.first->GetYaxis()->SetTitleSize(0.042);
  zjet.first->GetXaxis()->SetLabelSize(0.032);
  zjet.first->GetYaxis()->SetLabelSize(0.032);
  zjet.first->GetXaxis()->SetTitleOffset(1.1);
  zjet.first->GetYaxis()->SetTitleOffset(1.2);

  zjet.first->SetLineColor(kRed);
  zjet.second->SetLineColor(kRed);
  zjet.first->SetLineWidth(2);
  zjet.second->SetLineWidth(2);
  if(TString(postfix).Contains("qcdproc"))
    zjet.second->SetLineStyle(2);

  wjet.first->SetLineColor(kBlue);
  wjet.second->SetLineColor(kBlue);
  wjet.first->SetLineWidth(2);
  wjet.second->SetLineWidth(2);
  if(TString(postfix).Contains("qcdproc"))
    wjet.second->SetLineStyle(2);

  gamma.first->SetLineColor(kBlack);
  gamma.second->SetLineColor(kBlack);
  gamma.first->SetLineWidth(2);
  gamma.second->SetLineWidth(2);
  if(TString(postfix).Contains("qcdproc"))
    gamma.second->SetLineStyle(2);

  if(TString(postfix).Contains("qcdproc"))
    zjet.first->GetYaxis()->SetRangeUser(-0.1,0.1);
  else
    zjet.first->GetYaxis()->SetRangeUser(min(zjet.second->GetMinimum(),min(wjet.second->GetMinimum(),gamma.second->GetMinimum()))*0.9,
					 max(zjet.first->GetMaximum(),max(wjet.first->GetMaximum(),gamma.first->GetMaximum()))*1.1);
  zjet.first->Draw("hist");
  zjet.second->Draw("hist same");
  wjet.first->Draw("hist same");
  wjet.second->Draw("hist same");
  gamma.first->Draw("hist same");
  gamma.second->Draw("hist same");

  TLegend* leg = NULL;
  if(TString(postfix).Contains("qcdproc"))
    leg = new TLegend(0.75,0.75,0.90,0.9);
  else if(TString(postfix).Contains("ewk") or TString(postfix).Contains("sudakov") or TString(postfix).Contains("miss") or TString(postfix).Contains("mix"))
    leg = new TLegend(0.25,0.70,0.40,0.85);
  else
    leg = new TLegend(0.75,0.50,0.90,0.65);

  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(zjet.first,"Z+jets","L");
  leg->AddEntry(wjet.first,"W+jets","L");
  leg->AddEntry(gamma.first,"#gamma+jets","L");
  leg->Draw("same");

  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");

}

void makeUncertaintySingleProcess(string inputDIR, string outputDIR){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  TCanvas* canvas = new TCanvas("canvas","",600,600);

  // central values
  TFile* zw_central = TFile::Open((inputDIR+"/zwjcorewk.root").c_str(),"READ");
  TFile* zg_central = TFile::Open((inputDIR+"/gamcorewk.root").c_str(),"READ");
  
  TH1* zvv_central = (TH1*) zw_central->Get("nhist_zwj_ewk_met");
  TH1* wln_central = (TH1*) zw_central->Get("dhist_zwj_ewk_met");
  TH1* gam_central = (TH1*) zg_central->Get("dhist_gam_ewk_met");

  // qcd scale  
  TFile* zwewkfile_qcdscale_up = TFile::Open((inputDIR+"/zwjcorqcd_scaleup.root").c_str(),"READ");
  TFile* zwewkfile_qcdscale_dw = TFile::Open((inputDIR+"/zwjcorqcd_scaledw.root").c_str(),"READ");
  TFile* zgewkfile_qcdscale_up = TFile::Open((inputDIR+"/gamcorqcdscale_up.root").c_str(),"READ");
  TFile* zgewkfile_qcdscale_dw = TFile::Open((inputDIR+"/gamcorqcdscale_dw.root").c_str(),"READ");

  TH1* zvv_qcdscale_up = (TH1*) zwewkfile_qcdscale_up->Get("nhist_zwj_qcd_scaleup_met");
  TH1* zvv_qcdscale_dw = (TH1*) zwewkfile_qcdscale_dw->Get("nhist_zwj_qcd_scaledw_met");
  TH1* wln_qcdscale_up = (TH1*) zwewkfile_qcdscale_up->Get("dhist_zwj_qcd_scaleup_met");
  TH1* wln_qcdscale_dw = (TH1*) zwewkfile_qcdscale_dw->Get("dhist_zwj_qcd_scaledw_met");
  TH1* gam_qcdscale_up = (TH1*) zgewkfile_qcdscale_up->Get("dhist_gam_qcdscale_up_met");
  TH1* gam_qcdscale_dw = (TH1*) zgewkfile_qcdscale_dw->Get("dhist_gam_qcdscale_dw_met");

  TH1* zvv_qcdscale_unc = (TH1*) zvv_qcdscale_up->Clone("zvv_qcdscale_unc");
  for(int iBin = 0; iBin < zvv_qcdscale_unc->GetNbinsX()+1; iBin++)
    zvv_qcdscale_unc->SetBinContent(iBin,(zvv_qcdscale_up->GetBinContent(iBin)-zvv_qcdscale_dw->GetBinContent(iBin))/(2*zvv_central->GetBinContent(iBin)));
  TH1* zvv_qcdscale_unc_flip = flipHisto(zvv_qcdscale_unc);

  TH1* wln_qcdscale_unc = (TH1*) wln_qcdscale_up->Clone("wln_qcdscale_unc");
  for(int iBin = 0; iBin < wln_qcdscale_unc->GetNbinsX()+1; iBin++)
    wln_qcdscale_unc->SetBinContent(iBin,(wln_qcdscale_up->GetBinContent(iBin)-wln_qcdscale_dw->GetBinContent(iBin))/(2*wln_central->GetBinContent(iBin)));
  TH1* wln_qcdscale_unc_flip = flipHisto(wln_qcdscale_unc);

  TH1* gam_qcdscale_unc = (TH1*) gam_qcdscale_up->Clone("gam_qcdscale_unc");
  for(int iBin = 0; iBin < gam_qcdscale_unc->GetNbinsX()+1; iBin++)
    gam_qcdscale_unc->SetBinContent(iBin,(gam_qcdscale_up->GetBinContent(iBin)-gam_qcdscale_dw->GetBinContent(iBin))/(2*gam_central->GetBinContent(iBin)));
  TH1* gam_qcdscale_unc_flip = flipHisto(gam_qcdscale_unc);

  pair<TH1*,TH1*> zjet_qcdscale (zvv_qcdscale_unc,zvv_qcdscale_unc_flip);
  pair<TH1*,TH1*> wjet_qcdscale (wln_qcdscale_unc,wln_qcdscale_unc_flip);
  pair<TH1*,TH1*> gamma_qcdscale (gam_qcdscale_unc,gam_qcdscale_unc_flip);
  
  drawPlot(canvas,zjet_qcdscale,wjet_qcdscale,gamma_qcdscale,outputDIR,"QCD Scale uncertainty","qcdscale_unc");

  // qcd shape
  TFile* zwewkfile_qcdshape_up = TFile::Open((inputDIR+"/zwjcorqcdshape_up.root").c_str(),"READ");
  TFile* zwewkfile_qcdshape_dw = TFile::Open((inputDIR+"/zwjcorqcdshape_dw.root").c_str(),"READ");
  TFile* zgewkfile_qcdshape_up = TFile::Open((inputDIR+"/gamcorqcdshape_up.root").c_str(),"READ");
  TFile* zgewkfile_qcdshape_dw = TFile::Open((inputDIR+"/gamcorqcdshape_dw.root").c_str(),"READ");

  TH1* zvv_qcdshape_up = (TH1*) zwewkfile_qcdshape_up->Get("nhist_zwj_qcdshape_up_met");
  TH1* zvv_qcdshape_dw = (TH1*) zwewkfile_qcdshape_dw->Get("nhist_zwj_qcdshape_dw_met");
  TH1* wln_qcdshape_up = (TH1*) zwewkfile_qcdshape_up->Get("dhist_zwj_qcdshape_up_met");
  TH1* wln_qcdshape_dw = (TH1*) zwewkfile_qcdshape_dw->Get("dhist_zwj_qcdshape_dw_met");
  TH1* gam_qcdshape_up = (TH1*) zgewkfile_qcdshape_up->Get("dhist_gam_qcdshape_up_met");
  TH1* gam_qcdshape_dw = (TH1*) zgewkfile_qcdshape_dw->Get("dhist_gam_qcdshape_dw_met");

  TH1* zvv_qcdshape_unc = (TH1*) zvv_qcdshape_up->Clone("zvv_qcdshape_unc");
  for(int iBin = 0; iBin < zvv_qcdshape_unc->GetNbinsX()+1; iBin++)
    zvv_qcdshape_unc->SetBinContent(iBin,(zvv_qcdshape_up->GetBinContent(iBin)-zvv_qcdshape_dw->GetBinContent(iBin))/(2*zvv_central->GetBinContent(iBin)));
  TH1* zvv_qcdshape_unc_flip = flipHisto(zvv_qcdshape_unc);

  TH1* wln_qcdshape_unc = (TH1*) wln_qcdshape_up->Clone("wln_qcdshape_unc");
  for(int iBin = 0; iBin < wln_qcdshape_unc->GetNbinsX()+1; iBin++)
    wln_qcdshape_unc->SetBinContent(iBin,(wln_qcdshape_up->GetBinContent(iBin)-wln_qcdshape_dw->GetBinContent(iBin))/(2*wln_central->GetBinContent(iBin)));
  TH1* wln_qcdshape_unc_flip = flipHisto(wln_qcdshape_unc);

  TH1* gam_qcdshape_unc = (TH1*) gam_qcdshape_up->Clone("gam_qcdshape_unc");
  for(int iBin = 0; iBin < gam_qcdshape_unc->GetNbinsX()+1; iBin++)
    gam_qcdshape_unc->SetBinContent(iBin,(gam_qcdshape_up->GetBinContent(iBin)-gam_qcdshape_dw->GetBinContent(iBin))/(2*gam_central->GetBinContent(iBin)));
  TH1* gam_qcdshape_unc_flip = flipHisto(gam_qcdshape_unc);

  pair<TH1*,TH1*> zjet_qcdshape (zvv_qcdshape_unc,zvv_qcdshape_unc_flip);
  pair<TH1*,TH1*> wjet_qcdshape (wln_qcdshape_unc,wln_qcdshape_unc_flip);
  pair<TH1*,TH1*> gamma_qcdshape (gam_qcdshape_unc,gam_qcdshape_unc_flip);
  
  drawPlot(canvas,zjet_qcdshape,wjet_qcdshape,gamma_qcdshape,outputDIR,"QCD Shape uncertainty","qcdshape_unc");

  // qcd proc
  TFile* zwewkfile_qcdproc_up = TFile::Open((inputDIR+"/zwjcorqcdproc_up.root").c_str(),"READ");
  TFile* zwewkfile_qcdproc_dw = TFile::Open((inputDIR+"/zwjcorqcdproc_dw.root").c_str(),"READ");
  TFile* zgewkfile_qcdproc_up = TFile::Open((inputDIR+"/gamcorqcdproc_up.root").c_str(),"READ");
  TFile* zgewkfile_qcdproc_dw = TFile::Open((inputDIR+"/gamcorqcdproc_dw.root").c_str(),"READ");

  TH1* zvv_qcdproc_up = (TH1*) zwewkfile_qcdproc_up->Get("nhist_zwj_qcdproc_up_met");
  TH1* zvv_qcdproc_dw = (TH1*) zwewkfile_qcdproc_dw->Get("nhist_zwj_qcdproc_dw_met");
  TH1* wln_qcdproc_up = (TH1*) zwewkfile_qcdproc_up->Get("dhist_zwj_qcdproc_up_met");
  TH1* wln_qcdproc_dw = (TH1*) zwewkfile_qcdproc_dw->Get("dhist_zwj_qcdproc_dw_met");
  TH1* gam_qcdproc_up = (TH1*) zgewkfile_qcdproc_up->Get("dhist_gam_qcdproc_up_met");
  TH1* gam_qcdproc_dw = (TH1*) zgewkfile_qcdproc_dw->Get("dhist_gam_qcdproc_dw_met");

  TH1* zvv_qcdproc_unc = (TH1*) zvv_qcdproc_up->Clone("zvv_qcdproc_unc");
  for(int iBin = 0; iBin < zvv_qcdproc_unc->GetNbinsX()+1; iBin++)
    zvv_qcdproc_unc->SetBinContent(iBin,(zvv_qcdproc_up->GetBinContent(iBin)-zvv_qcdproc_dw->GetBinContent(iBin))/(2*zvv_central->GetBinContent(iBin)));
  TH1* zvv_qcdproc_unc_flip = flipHisto(zvv_qcdproc_unc);

  TH1* wln_qcdproc_unc = (TH1*) wln_qcdproc_up->Clone("wln_qcdproc_unc");
  for(int iBin = 0; iBin < wln_qcdproc_unc->GetNbinsX()+1; iBin++)
    wln_qcdproc_unc->SetBinContent(iBin,(wln_qcdproc_up->GetBinContent(iBin)-wln_qcdproc_dw->GetBinContent(iBin))/(2*wln_central->GetBinContent(iBin)));
  TH1* wln_qcdproc_unc_flip = flipHisto(wln_qcdproc_unc);

  TH1* gam_qcdproc_unc = (TH1*) gam_qcdproc_up->Clone("gam_qcdproc_unc");
  for(int iBin = 0; iBin < gam_qcdproc_unc->GetNbinsX()+1; iBin++)
    gam_qcdproc_unc->SetBinContent(iBin,(gam_qcdproc_up->GetBinContent(iBin)-gam_qcdproc_dw->GetBinContent(iBin))/(2*gam_central->GetBinContent(iBin)));
  TH1* gam_qcdproc_unc_flip = flipHisto(gam_qcdproc_unc);

  pair<TH1*,TH1*> zjet_qcdproc (zvv_qcdproc_unc,zvv_qcdproc_unc_flip);
  pair<TH1*,TH1*> wjet_qcdproc (wln_qcdproc_unc,wln_qcdproc_unc_flip);
  pair<TH1*,TH1*> gamma_qcdproc (gam_qcdproc_unc,gam_qcdproc_unc_flip);
  
  drawPlot(canvas,zjet_qcdproc,wjet_qcdproc,gamma_qcdproc,outputDIR,"QCD Proc uncertainty","qcdproc_unc");

  // NNLO EWK
  TFile* zwewkfile_nnloewk_up = TFile::Open((inputDIR+"/zwjcornnloewk_up.root").c_str(),"READ");
  TFile* zwewkfile_nnloewk_dw = TFile::Open((inputDIR+"/zwjcornnloewk_dw.root").c_str(),"READ");
  TFile* zgewkfile_nnloewk_up = TFile::Open((inputDIR+"/gamcornnloewk_up.root").c_str(),"READ");
  TFile* zgewkfile_nnloewk_dw = TFile::Open((inputDIR+"/gamcornnloewk_dw.root").c_str(),"READ");

  TH1* zvv_nnloewk_up = (TH1*) zwewkfile_nnloewk_up->Get("nhist_zwj_nnloewk_up_met");
  TH1* zvv_nnloewk_dw = (TH1*) zwewkfile_nnloewk_dw->Get("nhist_zwj_nnloewk_dw_met");
  TH1* wln_nnloewk_up = (TH1*) zwewkfile_nnloewk_up->Get("dhist_zwj_nnloewk_up_met");
  TH1* wln_nnloewk_dw = (TH1*) zwewkfile_nnloewk_dw->Get("dhist_zwj_nnloewk_dw_met");
  TH1* gam_nnloewk_up = (TH1*) zgewkfile_nnloewk_up->Get("dhist_gam_nnloewk_up_met");
  TH1* gam_nnloewk_dw = (TH1*) zgewkfile_nnloewk_dw->Get("dhist_gam_nnloewk_dw_met");

  TH1* zvv_nnloewk_unc = (TH1*) zvv_nnloewk_up->Clone("zvv_nnloewk_unc");
  for(int iBin = 0; iBin < zvv_nnloewk_unc->GetNbinsX()+1; iBin++)
    zvv_nnloewk_unc->SetBinContent(iBin,(zvv_nnloewk_up->GetBinContent(iBin)-zvv_nnloewk_dw->GetBinContent(iBin))/(2*zvv_central->GetBinContent(iBin)));
  TH1* zvv_nnloewk_unc_flip = flipHisto(zvv_nnloewk_unc);

  TH1* wln_nnloewk_unc = (TH1*) wln_nnloewk_up->Clone("wln_nnloewk_unc");
  for(int iBin = 0; iBin < wln_nnloewk_unc->GetNbinsX()+1; iBin++)
    wln_nnloewk_unc->SetBinContent(iBin,(wln_nnloewk_up->GetBinContent(iBin)-wln_nnloewk_dw->GetBinContent(iBin))/(2*wln_central->GetBinContent(iBin)));
  TH1* wln_nnloewk_unc_flip = flipHisto(wln_nnloewk_unc);

  TH1* gam_nnloewk_unc = (TH1*) gam_nnloewk_up->Clone("gam_nnloewk_unc");
  for(int iBin = 0; iBin < gam_nnloewk_unc->GetNbinsX()+1; iBin++)
    gam_nnloewk_unc->SetBinContent(iBin,(gam_nnloewk_up->GetBinContent(iBin)-gam_nnloewk_dw->GetBinContent(iBin))/(2*gam_central->GetBinContent(iBin)));
  TH1* gam_nnloewk_unc_flip = flipHisto(gam_nnloewk_unc);
 
  pair<TH1*,TH1*> zjet_nnloewk (zvv_nnloewk_unc,zvv_nnloewk_unc_flip);
  pair<TH1*,TH1*> wjet_nnloewk (wln_nnloewk_unc,wln_nnloewk_unc_flip);
  pair<TH1*,TH1*> gamma_nnloewk (gam_nnloewk_unc,gam_nnloewk_unc_flip);
  
  drawPlot(canvas,zjet_nnloewk,wjet_nnloewk,gamma_nnloewk,outputDIR,"NNLO EWK uncertainty","nnloewk_unc");

  // QCD-EWK
  TFile* zwewkfile_mix_up = TFile::Open((inputDIR+"/zwjcormix_up.root").c_str(),"READ");
  TFile* zwewkfile_mix_dw = TFile::Open((inputDIR+"/zwjcormix_dw.root").c_str(),"READ");
  TFile* zgewkfile_mix_up = TFile::Open((inputDIR+"/gamcormix_up.root").c_str(),"READ");
  TFile* zgewkfile_mix_dw = TFile::Open((inputDIR+"/gamcormix_dw.root").c_str(),"READ");

  TH1* zvv_mix_up = (TH1*) zwewkfile_mix_up->Get("nhist_zwj_mix_up_met");
  TH1* zvv_mix_dw = (TH1*) zwewkfile_mix_dw->Get("nhist_zwj_mix_dw_met");
  TH1* wln_mix_up = (TH1*) zwewkfile_mix_up->Get("dhist_zwj_mix_up_met");
  TH1* wln_mix_dw = (TH1*) zwewkfile_mix_dw->Get("dhist_zwj_mix_dw_met");
  TH1* gam_mix_up = (TH1*) zgewkfile_mix_up->Get("dhist_gam_mix_up_met");
  TH1* gam_mix_dw = (TH1*) zgewkfile_mix_dw->Get("dhist_gam_mix_dw_met");

  TH1* zvv_mix_unc = (TH1*) zvv_mix_up->Clone("zvv_mix_unc");
  for(int iBin = 0; iBin < zvv_mix_unc->GetNbinsX()+1; iBin++)
    zvv_mix_unc->SetBinContent(iBin,(zvv_mix_up->GetBinContent(iBin)-zvv_mix_dw->GetBinContent(iBin))/(2*zvv_central->GetBinContent(iBin)));
  TH1* zvv_mix_unc_flip = flipHisto(zvv_mix_unc);

  TH1* wln_mix_unc = (TH1*) wln_mix_up->Clone("wln_mix_unc");
  for(int iBin = 0; iBin < wln_mix_unc->GetNbinsX()+1; iBin++)
    wln_mix_unc->SetBinContent(iBin,(wln_mix_up->GetBinContent(iBin)-wln_mix_dw->GetBinContent(iBin))/(2*wln_central->GetBinContent(iBin)));
  TH1* wln_mix_unc_flip = flipHisto(wln_mix_unc);

  TH1* gam_mix_unc = (TH1*) gam_mix_up->Clone("gam_mix_unc");
  for(int iBin = 0; iBin < gam_mix_unc->GetNbinsX()+1; iBin++)
    gam_mix_unc->SetBinContent(iBin,(gam_mix_up->GetBinContent(iBin)-gam_mix_dw->GetBinContent(iBin))/(2*gam_central->GetBinContent(iBin)));
  TH1* gam_mix_unc_flip = flipHisto(gam_mix_unc);
 
  pair<TH1*,TH1*> zjet_mix (zvv_mix_unc,zvv_mix_unc_flip);
  pair<TH1*,TH1*> wjet_mix (wln_mix_unc,wln_mix_unc_flip);
  pair<TH1*,TH1*> gamma_mix (gam_mix_unc,gam_mix_unc_flip);
  
  drawPlot(canvas,zjet_mix,wjet_mix,gamma_mix,outputDIR,"NNLO EWK uncertainty","mix_unc");

  // SUDAKOV 

  TFile* zvvfile_sudakov_up   = TFile::Open((inputDIR+"/zwjcorsudakov_up_1.root").c_str(),"READ");
  TFile* zvvfile_sudakov_dw   = TFile::Open((inputDIR+"/zwjcorsudakov_dw_1.root").c_str(),"READ");
  TFile* wlnfile_sudakov_up   = TFile::Open((inputDIR+"/zwjcorsudakov_up_2.root").c_str(),"READ");
  TFile* wlnfile_sudakov_dw   = TFile::Open((inputDIR+"/zwjcorsudakov_dw_2.root").c_str(),"READ");
  TFile* gamfile_sudakov_up = TFile::Open((inputDIR+"/gamcorsudakov_up_2.root").c_str(),"READ");
  TFile* gamfile_sudakov_dw = TFile::Open((inputDIR+"/gamcorsudakov_dw_2.root").c_str(),"READ");

  TH1* zvv_sudakov_up = (TH1*) zvvfile_sudakov_up->Get("nhist_zwj_sudakov_up_1_met");
  TH1* zvv_sudakov_dw = (TH1*) zvvfile_sudakov_dw->Get("nhist_zwj_sudakov_dw_1_met");
  TH1* wln_sudakov_up = (TH1*) wlnfile_sudakov_up->Get("dhist_zwj_sudakov_up_2_met");
  TH1* wln_sudakov_dw = (TH1*) wlnfile_sudakov_dw->Get("dhist_zwj_sudakov_dw_2_met");
  TH1* gam_sudakov_up = (TH1*) gamfile_sudakov_up->Get("dhist_gam_sudakov_up_2_met");
  TH1* gam_sudakov_dw = (TH1*) gamfile_sudakov_dw->Get("dhist_gam_sudakov_dw_2_met");

  TH1* zvv_sudakov_unc = (TH1*) zvv_sudakov_up->Clone("zvv_sudakov_unc");
  for(int iBin = 0; iBin < zvv_sudakov_unc->GetNbinsX()+1; iBin++)
    zvv_sudakov_unc->SetBinContent(iBin,(zvv_sudakov_up->GetBinContent(iBin)-zvv_sudakov_dw->GetBinContent(iBin))/(2*zvv_central->GetBinContent(iBin)));
  TH1* zvv_sudakov_unc_flip = flipHisto(zvv_sudakov_unc);

  TH1* wln_sudakov_unc = (TH1*) wln_sudakov_up->Clone("wln_sudakov_unc");
  for(int iBin = 0; iBin < wln_sudakov_unc->GetNbinsX()+1; iBin++)
    wln_sudakov_unc->SetBinContent(iBin,(wln_sudakov_up->GetBinContent(iBin)-wln_sudakov_dw->GetBinContent(iBin))/(2*wln_central->GetBinContent(iBin)));
  TH1* wln_sudakov_unc_flip = flipHisto(wln_sudakov_unc);

  TH1* gam_sudakov_unc = (TH1*) gam_sudakov_up->Clone("gam_sudakov_unc");
  for(int iBin = 0; iBin < gam_sudakov_unc->GetNbinsX()+1; iBin++)
    gam_sudakov_unc->SetBinContent(iBin,(gam_sudakov_up->GetBinContent(iBin)-gam_sudakov_dw->GetBinContent(iBin))/(2*gam_central->GetBinContent(iBin)));
  TH1* gam_sudakov_unc_flip = flipHisto(gam_sudakov_unc);
 
  pair<TH1*,TH1*> zjet_sudakov (zvv_sudakov_unc,zvv_sudakov_unc_flip);
  pair<TH1*,TH1*> wjet_sudakov (wln_sudakov_unc,wln_sudakov_unc_flip);
  pair<TH1*,TH1*> gamma_sudakov (gam_sudakov_unc,gam_sudakov_unc_flip);
  
  drawPlot(canvas,zjet_sudakov,wjet_sudakov,gamma_sudakov,outputDIR,"NNLO Sudakov uncertainty","sudakov_unc");

  // SUDAKOV 

  TFile* zvvfile_nnlomiss_up   = TFile::Open((inputDIR+"/zwjcornnlomiss_up_1.root").c_str(),"READ");
  TFile* zvvfile_nnlomiss_dw   = TFile::Open((inputDIR+"/zwjcornnlomiss_dw_1.root").c_str(),"READ");
  TFile* wlnfile_nnlomiss_up   = TFile::Open((inputDIR+"/zwjcornnlomiss_up_2.root").c_str(),"READ");
  TFile* wlnfile_nnlomiss_dw   = TFile::Open((inputDIR+"/zwjcornnlomiss_dw_2.root").c_str(),"READ");
  TFile* gamfile_nnlomiss_up = TFile::Open((inputDIR+"/gamcornnlomiss_up_2.root").c_str(),"READ");
  TFile* gamfile_nnlomiss_dw = TFile::Open((inputDIR+"/gamcornnlomiss_dw_2.root").c_str(),"READ");

  TH1* zvv_nnlomiss_up = (TH1*) zvvfile_nnlomiss_up->Get("nhist_zwj_nnlomiss_up_1_met");
  TH1* zvv_nnlomiss_dw = (TH1*) zvvfile_nnlomiss_dw->Get("nhist_zwj_nnlomiss_dw_1_met");
  TH1* wln_nnlomiss_up = (TH1*) wlnfile_nnlomiss_up->Get("dhist_zwj_nnlomiss_up_2_met");
  TH1* wln_nnlomiss_dw = (TH1*) wlnfile_nnlomiss_dw->Get("dhist_zwj_nnlomiss_dw_2_met");
  TH1* gam_nnlomiss_up = (TH1*) gamfile_nnlomiss_up->Get("dhist_gam_nnlomiss_up_2_met");
  TH1* gam_nnlomiss_dw = (TH1*) gamfile_nnlomiss_dw->Get("dhist_gam_nnlomiss_dw_2_met");

  TH1* zvv_nnlomiss_unc = (TH1*) zvv_nnlomiss_up->Clone("zvv_nnlomiss_unc");
  for(int iBin = 0; iBin < zvv_nnlomiss_unc->GetNbinsX()+1; iBin++)
    zvv_nnlomiss_unc->SetBinContent(iBin,(zvv_nnlomiss_up->GetBinContent(iBin)-zvv_nnlomiss_dw->GetBinContent(iBin))/(2*zvv_central->GetBinContent(iBin)));
  TH1* zvv_nnlomiss_unc_flip = flipHisto(zvv_nnlomiss_unc);

  TH1* wln_nnlomiss_unc = (TH1*) wln_nnlomiss_up->Clone("wln_nnlomiss_unc");
  for(int iBin = 0; iBin < wln_nnlomiss_unc->GetNbinsX()+1; iBin++)
    wln_nnlomiss_unc->SetBinContent(iBin,(wln_nnlomiss_up->GetBinContent(iBin)-wln_nnlomiss_dw->GetBinContent(iBin))/(2*wln_central->GetBinContent(iBin)));
  TH1* wln_nnlomiss_unc_flip = flipHisto(wln_nnlomiss_unc);

  TH1* gam_nnlomiss_unc = (TH1*) gam_nnlomiss_up->Clone("gam_nnlomiss_unc");
  for(int iBin = 0; iBin < gam_nnlomiss_unc->GetNbinsX()+1; iBin++)
    gam_nnlomiss_unc->SetBinContent(iBin,(gam_nnlomiss_up->GetBinContent(iBin)-gam_nnlomiss_dw->GetBinContent(iBin))/(2*gam_central->GetBinContent(iBin)));
  TH1* gam_nnlomiss_unc_flip = flipHisto(gam_nnlomiss_unc);
 
  pair<TH1*,TH1*> zjet_nnlomiss (zvv_nnlomiss_unc,zvv_nnlomiss_unc_flip);
  pair<TH1*,TH1*> wjet_nnlomiss (wln_nnlomiss_unc,wln_nnlomiss_unc_flip);
  pair<TH1*,TH1*> gamma_nnlomiss (gam_nnlomiss_unc,gam_nnlomiss_unc_flip);
  
  drawPlot(canvas,zjet_nnlomiss,wjet_nnlomiss,gamma_nnlomiss,outputDIR,"NNLO Miss uncertainty","nnlomiss_unc");

}
