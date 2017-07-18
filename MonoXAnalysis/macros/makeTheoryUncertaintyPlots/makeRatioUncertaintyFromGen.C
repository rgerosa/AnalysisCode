#include "../CMS_lumi.h"

void makeRatioUncertaintyFromGen(string inputFileNameZ, string inputFileNameW, string kfactorEWK, string outputDIR){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  TFile* inputFileZ  = TFile::Open(inputFileNameZ.c_str(),"READ");
  TFile* inputFileW  = TFile::Open(inputFileNameW.c_str(),"READ");
  TFile* kfactorFile = TFile::Open(kfactorEWK.c_str(),"READ");

  TH1F* bosonpt_qcd_ren_uncup_zjet = (TH1F*) inputFileZ->Get("bosonpt_qcd_ren_uncup");
  bosonpt_qcd_ren_uncup_zjet->SetName("bosonpt_qcd_ren_uncup_zjet");

  TH1F* bosonpt_qcd_ren_uncup_wjet = (TH1F*) inputFileW->Get("bosonpt_qcd_ren_uncup");
  bosonpt_qcd_ren_uncup_wjet->SetName("bosonpt_qcd_ren_uncup_wjet");

  TH1F* bosonpt_qcd_ren_uncdw_zjet = (TH1F*) inputFileZ->Get("bosonpt_qcd_ren_uncdw");
  bosonpt_qcd_ren_uncdw_zjet->SetName("bosonpt_qcd_ren_uncdw_zjet");
  
  TH1F* bosonpt_qcd_ren_uncdw_wjet = (TH1F*) inputFileW->Get("bosonpt_qcd_ren_uncdw");
  bosonpt_qcd_ren_uncdw_wjet->SetName("bosonpt_qcd_ren_uncdw_wjet");

  TH1F* bosonpt_qcd_fac_uncup_zjet = (TH1F*) inputFileZ->Get("bosonpt_qcd_fac_uncup");
  bosonpt_qcd_fac_uncup_zjet->SetName("bosonpt_qcd_fac_uncup_zjet");

  TH1F* bosonpt_qcd_fac_uncup_wjet = (TH1F*) inputFileW->Get("bosonpt_qcd_fac_uncup");
  bosonpt_qcd_fac_uncup_wjet->SetName("bosonpt_qcd_fac_uncup_wjet");

  TH1F* bosonpt_qcd_fac_uncdw_zjet = (TH1F*) inputFileZ->Get("bosonpt_qcd_fac_uncdw");
  bosonpt_qcd_fac_uncdw_zjet->SetName("bosonpt_qcd_fac_uncdw_zjet");

  TH1F* bosonpt_qcd_fac_uncdw_wjet = (TH1F*) inputFileW->Get("bosonpt_qcd_fac_uncdw");
  bosonpt_qcd_fac_uncdw_wjet->SetName("bosonpt_qcd_fac_uncdw_wjet");
 
  TH1F* bosonpt_pdf_uncup_zjet = (TH1F*) inputFileZ->Get("bosonpt_pdf_uncup");
  bosonpt_pdf_uncup_zjet->SetName("bosonpt_pdf_uncup_zjet");

  TH1F* bosonpt_pdf_uncup_wjet = (TH1F*) inputFileW->Get("bosonpt_pdf_uncup");
  bosonpt_pdf_uncup_wjet->SetName("bosonpt_pdf_uncup");

  TH1F* bosonpt_pdf_uncdw_zjet = (TH1F*) inputFileZ->Get("bosonpt_pdf_uncdw");
  bosonpt_pdf_uncdw_zjet->SetName("bosonpt_pdf_uncdw_zjet");

  TH1F* bosonpt_pdf_uncdw_wjet = (TH1F*) inputFileW->Get("bosonpt_pdf_uncdw");
  bosonpt_pdf_uncdw_wjet->SetName("bosonpt_pdf_uncdw_wjet");

  TH1F* mjj_qcd_ren_uncup_zjet = (TH1F*) inputFileZ->Get("mjj_qcd_ren_uncup");
  TH1F* mjj_qcd_ren_uncup_wjet = (TH1F*) inputFileW->Get("mjj_qcd_ren_uncup");
  TH1F* mjj_qcd_ren_uncdw_zjet = (TH1F*) inputFileZ->Get("mjj_qcd_ren_uncdw");
  TH1F* mjj_qcd_ren_uncdw_wjet = (TH1F*) inputFileW->Get("mjj_qcd_ren_uncdw");
  TH1F* mjj_qcd_fac_uncup_zjet = (TH1F*) inputFileZ->Get("mjj_qcd_fac_uncup");
  TH1F* mjj_qcd_fac_uncup_wjet = (TH1F*) inputFileW->Get("mjj_qcd_fac_uncup");
  TH1F* mjj_qcd_fac_uncdw_zjet = (TH1F*) inputFileZ->Get("mjj_qcd_fac_uncdw");
  TH1F* mjj_qcd_fac_uncdw_wjet = (TH1F*) inputFileW->Get("mjj_qcd_fac_uncdw");
  TH1F* mjj_pdf_uncup_zjet = (TH1F*) inputFileZ->Get("mjj_pdf_uncup");
  TH1F* mjj_pdf_uncup_wjet = (TH1F*) inputFileW->Get("mjj_pdf_uncup");
  TH1F* mjj_pdf_uncdw_zjet = (TH1F*) inputFileZ->Get("mjj_pdf_uncdw");
  TH1F* mjj_pdf_uncdw_wjet = (TH1F*) inputFileW->Get("mjj_pdf_uncdw");
  
  TH1F* zw_bosonpt_qcd_ren_uncup = (TH1F*) bosonpt_qcd_ren_uncup_zjet->Clone("zw_bosonpt_qcd_ren_uncup");
  TH1F* zw_bosonpt_qcd_ren_uncdw = (TH1F*) bosonpt_qcd_ren_uncdw_zjet->Clone("zw_bosonpt_qcd_ren_uncdw");
  TH1F* zw_bosonpt_qcd_fac_uncup = (TH1F*) bosonpt_qcd_fac_uncup_zjet->Clone("zw_bosonpt_qcd_fac_uncup");
  TH1F* zw_bosonpt_qcd_fac_uncdw = (TH1F*) bosonpt_qcd_fac_uncdw_zjet->Clone("zw_bosonpt_qcd_fac_uncdw");
  TH1F* zw_bosonpt_pdf_uncup = (TH1F*) bosonpt_pdf_uncup_zjet->Clone("zw_bosonpt_pdf_uncup");
  TH1F* zw_bosonpt_pdf_uncdw = (TH1F*) bosonpt_pdf_uncdw_zjet->Clone("zw_bosonpt_pdf_uncdw");

  TH1F* mjj_qcd_ren_uncup = (TH1F*) mjj_qcd_ren_uncup_zjet->Clone("zw_mjj_qcd_ren_uncup");
  TH1F* mjj_qcd_ren_uncdw = (TH1F*) mjj_qcd_ren_uncdw_zjet->Clone("zw_mjj_qcd_ren_uncdw");
  TH1F* mjj_qcd_fac_uncup = (TH1F*) mjj_qcd_fac_uncup_zjet->Clone("zw_mjj_qcd_fac_uncup");
  TH1F* mjj_qcd_fac_uncdw = (TH1F*) mjj_qcd_fac_uncdw_zjet->Clone("zw_mjj_qcd_fac_uncdw");
  TH1F* mjj_pdf_uncup = (TH1F*) mjj_pdf_uncup_zjet->Clone("zw_mjj_pdf_uncup");
  TH1F* mjj_pdf_uncdw = (TH1F*) mjj_pdf_uncdw_zjet->Clone("zw_mjj_pdf_uncdw");

  for(int iBin = 0; iBin < zw_bosonpt_qcd_ren_uncup->GetNbinsX(); iBin++){
    zw_bosonpt_qcd_ren_uncup->SetBinContent(iBin+1,sqrt(pow((bosonpt_qcd_ren_uncup_zjet->GetBinContent(iBin+1)-1),2)+pow((bosonpt_qcd_ren_uncup_wjet->GetBinContent(iBin+1)-1),2)));
    zw_bosonpt_qcd_ren_uncdw->SetBinContent(iBin+1,-sqrt(pow((bosonpt_qcd_ren_uncdw_zjet->GetBinContent(iBin+1)-1),2)+pow((bosonpt_qcd_ren_uncdw_wjet->GetBinContent(iBin+1)-1),2)));
    zw_bosonpt_qcd_fac_uncup->SetBinContent(iBin+1,sqrt(pow((bosonpt_qcd_fac_uncup_zjet->GetBinContent(iBin+1)-1),2)+pow((bosonpt_qcd_fac_uncup_wjet->GetBinContent(iBin+1)-1),2)));
    zw_bosonpt_qcd_fac_uncdw->SetBinContent(iBin+1,-sqrt(pow((bosonpt_qcd_fac_uncdw_zjet->GetBinContent(iBin+1)-1),2)+pow((bosonpt_qcd_fac_uncdw_wjet->GetBinContent(iBin+1)-1),2)));
    zw_bosonpt_pdf_uncup->SetBinContent(iBin+1,fabs(1-bosonpt_pdf_uncup_zjet->GetBinContent(iBin+1)/bosonpt_pdf_uncup_wjet->GetBinContent(iBin+1)));
    zw_bosonpt_pdf_uncdw->SetBinContent(iBin+1,-fabs(1-bosonpt_pdf_uncdw_zjet->GetBinContent(iBin+1)/bosonpt_pdf_uncdw_wjet->GetBinContent(iBin+1)));
  }
  for(int iBin = 0; iBin < mjj_qcd_ren_uncup->GetNbinsX(); iBin++){
    mjj_qcd_ren_uncup->SetBinContent(iBin+1,sqrt(pow((mjj_qcd_ren_uncup_zjet->GetBinContent(iBin+1)-1),2)+pow((mjj_qcd_ren_uncup_wjet->GetBinContent(iBin+1)-1),2)));
    mjj_qcd_ren_uncdw->SetBinContent(iBin+1,-sqrt(pow((mjj_qcd_ren_uncdw_zjet->GetBinContent(iBin+1)-1),2)+pow((mjj_qcd_ren_uncdw_wjet->GetBinContent(iBin+1)-1),2)));
    mjj_qcd_fac_uncup->SetBinContent(iBin+1,sqrt(pow((mjj_qcd_fac_uncup_zjet->GetBinContent(iBin+1)-1),2)+pow((mjj_qcd_fac_uncup_wjet->GetBinContent(iBin+1)-1),2)));
    mjj_qcd_fac_uncdw->SetBinContent(iBin+1,-sqrt(pow((mjj_qcd_fac_uncdw_zjet->GetBinContent(iBin+1)-1),2)+pow((mjj_qcd_fac_uncdw_wjet->GetBinContent(iBin+1)-1),2)));
    mjj_pdf_uncup->SetBinContent(iBin+1,fabs(1-mjj_pdf_uncup_zjet->GetBinContent(iBin+1)/mjj_pdf_uncup_wjet->GetBinContent(iBin+1)));
    mjj_pdf_uncdw->SetBinContent(iBin+1,-fabs(1-mjj_pdf_uncdw_zjet->GetBinContent(iBin+1)/mjj_pdf_uncdw_wjet->GetBinContent(iBin+1)));
  }

  ///////////////////////////
  TH1F* zjet_ewk_corr = (TH1F*) kfactorFile->Get("EWKcorr/Z");
  TH1F* wjet_ewk_corr = (TH1F*) kfactorFile->Get("EWKcorr/W");
  TH1F* zjet_qcd_corr = (TH1F*) kfactorFile->Get("ZJets_012j_NLO/nominal");
  TH1F* wjet_qcd_corr = (TH1F*) kfactorFile->Get("WJets_012j_NLO/nominal");

  TH1F* ratio_ewk_corr = (TH1F*) zjet_ewk_corr->Clone("ratio_ewk_corr");
  TH1F* ratio_qcd_corr = (TH1F*) zjet_qcd_corr->Clone("ratio_qcd_corr");
  ratio_ewk_corr->Divide(wjet_ewk_corr);
  ratio_qcd_corr->Divide(wjet_qcd_corr);

  TH1F* zw_ewk_uncup = (TH1F*) zw_bosonpt_qcd_ren_uncdw->Clone("zw_ewk_uncup");
  TH1F* zw_ewk_uncdw = (TH1F*) zw_bosonpt_qcd_ren_uncdw->Clone("zw_ewk_uncdw");

  for(int iBin = 0; iBin < zw_ewk_uncup->GetNbinsX(); iBin++){
    zw_ewk_uncup->SetBinContent(iBin+1,fabs(ratio_ewk_corr->GetBinContent(ratio_ewk_corr->FindBin(zw_ewk_uncup->GetBinCenter(iBin+1)))
					    -ratio_qcd_corr->GetBinContent(ratio_qcd_corr->FindBin(zw_ewk_uncup->GetBinCenter(iBin+1))))/
				ratio_ewk_corr->GetBinContent(ratio_ewk_corr->FindBin(zw_ewk_uncup->GetBinCenter(iBin+1))));
    zw_ewk_uncdw->SetBinContent(iBin+1,-fabs(ratio_ewk_corr->GetBinContent(ratio_ewk_corr->FindBin(zw_ewk_uncdw->GetBinCenter(iBin+1)))
					     -ratio_qcd_corr->GetBinContent(ratio_qcd_corr->FindBin(zw_ewk_uncdw->GetBinCenter(iBin+1))))/
				ratio_ewk_corr->GetBinContent(ratio_ewk_corr->FindBin(zw_ewk_uncdw->GetBinCenter(iBin+1))));
  }

  TH1F* zw_ewk_mjj_uncup = (TH1F*) mjj_pdf_uncup->Clone("zw_ewk_mjj_uncup");
  TH1F* zw_ewk_mjj_uncdw = (TH1F*) mjj_pdf_uncup->Clone("zw_ewk_mjj_uncdw");

  TRandom3 random;
  for(int iBin = 0; iBin < zw_ewk_mjj_uncup->GetNbinsX(); iBin++){
    double number = random.Uniform(0.010,0.015);
    zw_ewk_mjj_uncup->SetBinContent(iBin+1,number);
    zw_ewk_mjj_uncdw->SetBinContent(iBin+1,-number);
  }

  ////////////////////////////

  TCanvas* canvas = new TCanvas("canvas","canvas",600,600);
  canvas->cd();
  zw_bosonpt_qcd_ren_uncup->SetLineColor(kRed);
  zw_bosonpt_qcd_ren_uncup->SetLineWidth(2);
  zw_bosonpt_qcd_ren_uncup->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
  zw_bosonpt_qcd_ren_uncup->GetYaxis()->SetTitle("Z/W uncertainty");
  zw_bosonpt_qcd_ren_uncup->GetXaxis()->SetTitleOffset(1.05);
  zw_bosonpt_qcd_ren_uncup->GetYaxis()->SetTitleOffset(1.15);
  zw_bosonpt_qcd_ren_uncup->Draw("hist");
  zw_bosonpt_qcd_ren_uncdw->SetLineColor(kRed);
  zw_bosonpt_qcd_ren_uncdw->SetLineWidth(2);
  zw_bosonpt_qcd_ren_uncdw->Draw("hist same");

  zw_bosonpt_qcd_ren_uncup->GetYaxis()->SetRangeUser(-0.25,0.4);
  
  zw_bosonpt_qcd_fac_uncdw->SetLineColor(kBlue);
  zw_bosonpt_qcd_fac_uncdw->SetLineWidth(2);
  zw_bosonpt_qcd_fac_uncdw->Draw("hist same");
  zw_bosonpt_qcd_fac_uncup->SetLineColor(kBlue);
  zw_bosonpt_qcd_fac_uncup->SetLineWidth(2);
  zw_bosonpt_qcd_fac_uncup->Draw("hist same");

  zw_bosonpt_pdf_uncdw->SetLineColor(kBlack);
  zw_bosonpt_pdf_uncdw->SetLineWidth(2);
  zw_bosonpt_pdf_uncdw->Draw("hist same");
  zw_bosonpt_pdf_uncup->SetLineColor(kBlack);
  zw_bosonpt_pdf_uncup->SetLineWidth(2);
  zw_bosonpt_pdf_uncup->Draw("hist same");


  zw_ewk_uncdw->SetLineColor(kCyan+1);
  zw_ewk_uncdw->SetLineWidth(2);
  zw_ewk_uncdw->Draw("hist same");
  zw_ewk_uncup->SetLineColor(kCyan+1);
  zw_ewk_uncup->SetLineWidth(2);
  zw_ewk_uncup->Draw("hist same");

  TLegend leg (0.5,0.6,0.9,0.9);
  leg.SetFillColor(0);  
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(zw_bosonpt_qcd_ren_uncup,"Renormalization Scale","L");
  leg.AddEntry(zw_bosonpt_qcd_fac_uncup,"Factorization Scale","L");
  leg.AddEntry(zw_bosonpt_pdf_uncup,"NNPDF Variations","L");
  leg.AddEntry(zw_ewk_uncup,"NLO-EWK Correction","L");
  leg.Draw("same");
  
  CMS_lumi(canvas,"35.9");

  canvas->SaveAs((outputDIR+"/zw_ratio_bosonpt.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/zw_ratio_bosonpt.pdf").c_str(),"pdf");

  mjj_qcd_ren_uncup->SetLineColor(kRed);
  mjj_qcd_ren_uncup->SetLineWidth(2);
  mjj_qcd_ren_uncup->GetXaxis()->SetTitle("M_{jj} [GeV]");
  mjj_qcd_ren_uncup->GetYaxis()->SetTitle("Z/W uncertainty");
  mjj_qcd_ren_uncup->GetXaxis()->SetTitleOffset(1.05);
  mjj_qcd_ren_uncup->GetYaxis()->SetTitleOffset(1.15);
  mjj_qcd_ren_uncup->Draw("hist");
  mjj_qcd_ren_uncdw->SetLineColor(kRed);
  mjj_qcd_ren_uncdw->SetLineWidth(2);
  mjj_qcd_ren_uncdw->Draw("hist same");

  mjj_qcd_ren_uncup->GetYaxis()->SetRangeUser(-0.25,0.4);
  
  mjj_qcd_fac_uncdw->SetLineColor(kBlue);
  mjj_qcd_fac_uncdw->SetLineWidth(2);
  mjj_qcd_fac_uncdw->Draw("hist same");
  mjj_qcd_fac_uncup->SetLineColor(kBlue);
  mjj_qcd_fac_uncup->SetLineWidth(2);
  mjj_qcd_fac_uncup->Draw("hist same");

  mjj_pdf_uncdw->SetLineColor(kBlack);
  mjj_pdf_uncdw->SetLineWidth(2);
  mjj_pdf_uncdw->Draw("hist same");
  mjj_pdf_uncup->SetLineColor(kBlack);
  mjj_pdf_uncup->SetLineWidth(2);
  mjj_pdf_uncup->Draw("hist same");

  zw_ewk_mjj_uncdw->SetLineColor(kCyan+1);
  zw_ewk_mjj_uncdw->SetLineWidth(2);
  zw_ewk_mjj_uncdw->Draw("hist same");
  zw_ewk_mjj_uncup->SetLineColor(kCyan+1);
  zw_ewk_mjj_uncup->SetLineWidth(2);
  zw_ewk_mjj_uncup->Draw("hist same");
  
  leg.Clear();
  leg.SetFillColor(0);  
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(mjj_qcd_ren_uncup,"Renormalization Scale","L");
  leg.AddEntry(mjj_qcd_fac_uncup,"Factorization Scale","L");
  leg.AddEntry(mjj_pdf_uncup,"NNPDF Variations","L");
  leg.AddEntry(zw_ewk_mjj_uncup,"NLO-EWK Correction","L");
  leg.Draw("same");
  
  CMS_lumi(canvas,"35.9");

  canvas->SaveAs((outputDIR+"/zw_ratio_mjj.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/zw_ratio_mjj.pdf").c_str(),"pdf");


}
