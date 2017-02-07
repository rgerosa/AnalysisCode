#include "../CMS_lumi.h"

TH1* flipHisto(TH1* histo){
  TH1* histo2 = (TH1*) histo->Clone(Form("%s_flip",histo->GetName()));
  for(int iBin =1; iBin <= histo2->GetNbinsX(); iBin++)
    histo2->SetBinContent(iBin,-histo->GetBinContent(iBin));
  return histo2;
}


void makeTheoryUncertaintyBand(string inputFileName, string outputDIR , bool useNewTheoryUncertainties, string observable){

  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  TFile* inputFile = TFile::Open(inputFileName.c_str());
  
  TH1* zw_fac1_unc = NULL;
  TH1* zw_fac2_unc = NULL;
  TH1* zw_ren1_unc = NULL;
  TH1* zw_ren2_unc = NULL;
  TH1* zw_pdf_unc = NULL;
  TH1* zw_ewk_unc = NULL;
  TH1* zw_qcdscale_unc = NULL;
  TH1* zw_nloewk_unc = NULL;
  TH1* zw_sudewk_unc = NULL;
  TH1* zw_ewkqcd_unc = NULL;

  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();

  if(not useNewTheoryUncertainties){
    zw_fac1_unc = (TH1*) inputFile->FindObjectAny(("ZW_FactScale1_"+observable).c_str());
    zw_fac2_unc = (TH1*) inputFile->FindObjectAny(("ZW_FactScale2_"+observable).c_str());
    zw_ren1_unc = (TH1*) inputFile->FindObjectAny(("ZW_RenScale1_"+observable).c_str());
    zw_ren2_unc = (TH1*) inputFile->FindObjectAny(("ZW_RenScale2_"+observable).c_str());
    zw_pdf_unc  = (TH1*) inputFile->FindObjectAny(("ZW_PDF_"+observable).c_str());
    zw_ewk_unc  = (TH1*) inputFile->FindObjectAny(("ZW_EWK_"+observable).c_str());

    TH1* zw_fac_unc = (TH1*) zw_fac1_unc->Clone("zw_fac_unc");
    for(int iBin = 1; iBin <= zw_fac_unc->GetNbinsX(); iBin++)
      zw_fac_unc->SetBinContent(iBin,sqrt(zw_fac1_unc->GetBinContent(iBin)*zw_fac1_unc->GetBinContent(iBin)+zw_fac2_unc->GetBinContent(iBin)*zw_fac2_unc->GetBinContent(iBin)));

    TH1* zw_ren_unc = (TH1*) zw_ren1_unc->Clone("zw_ren_unc");
    for(int iBin = 1; iBin <= zw_ren_unc->GetNbinsX(); iBin++)
      zw_ren_unc->SetBinContent(iBin,sqrt(zw_ren1_unc->GetBinContent(iBin)*zw_ren1_unc->GetBinContent(iBin)+zw_ren2_unc->GetBinContent(iBin)*zw_ren2_unc->GetBinContent(iBin)));

    TH1* zw_fac_unc_flip = flipHisto(zw_fac_unc);
    TH1* zw_ren_unc_flip = flipHisto(zw_ren_unc);
    TH1* zw_pdf_unc_flip = flipHisto(zw_pdf_unc);
    TH1* zw_ewk_unc_flip = flipHisto(zw_ewk_unc);

    zw_ren_unc->SetLineColor(kBlack);
    zw_ren_unc->SetLineWidth(2);
    zw_ren_unc_flip->SetLineColor(kBlack);
    zw_ren_unc_flip->SetLineWidth(2);
    zw_ren_unc->GetXaxis()->SetTitle("m_{jj} [GeV]");
    zw_ren_unc->GetXaxis()->SetTitleOffset(1.15);
    zw_ren_unc->GetXaxis()->SetTitleSize(0.05);
    zw_ren_unc->GetYaxis()->SetTitle("Z/W Variation/Nominal");
    zw_ren_unc->GetYaxis()->SetTitleSize(0.042);
    zw_ren_unc->GetYaxis()->SetTitleOffset(1.25);

    zw_ren_unc->GetYaxis()->SetRangeUser(-0.15,0.25);
    zw_ren_unc->Draw("hist");
    zw_ren_unc_flip->Draw("hist same");
    
    zw_fac_unc_flip->SetLineColor(kRed);
    zw_fac_unc_flip->SetLineWidth(2);
    zw_fac_unc->SetLineColor(kRed);
    zw_fac_unc->SetLineWidth(2);
    
    zw_fac_unc->Draw("hist same");
    zw_fac_unc_flip->Draw("hist same");

    zw_pdf_unc_flip->SetLineColor(kGreen+1);
    zw_pdf_unc_flip->SetLineWidth(2);
    zw_pdf_unc->SetLineColor(kGreen+1);
    zw_pdf_unc->SetLineWidth(2);

    zw_pdf_unc->Draw("hist same");
    zw_pdf_unc_flip->Draw("hist same");

    zw_ewk_unc_flip->SetLineColor(kBlue);
    zw_ewk_unc_flip->SetLineWidth(2);
    zw_ewk_unc->SetLineColor(kBlue);
    zw_ewk_unc->SetLineWidth(2);

    zw_ewk_unc->Draw("hist same");
    zw_ewk_unc_flip->Draw("hist same");
    
    TLegend leg (0.7,0.7,0.9,0.9);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(zw_ren_unc,"Ren-Scale #mu_{r}","L");
    leg.AddEntry(zw_fac_unc,"Fact-Scale #mu_{f}","L");
    leg.AddEntry(zw_pdf_unc,"PDF","L");
    leg.AddEntry(zw_ewk_unc,"NLO-EWK","L");
    leg.Draw("same");

    canvas->SaveAs((outputDIR+"/thoeryUnc_ZW_ratio.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/thoeryUnc_ZW_ratio.pdf").c_str(),"pdf");

  }
  else{
    zw_qcdscale_unc = (TH1*) inputFile->FindObjectAny(("ZW_QCDScale_"+observable).c_str());
    zw_nloewk_unc = (TH1*) inputFile->FindObjectAny(("ZW_NLOEWK_"+observable).c_str());
    zw_sudewk_unc = (TH1*) inputFile->FindObjectAny(("ZW_EWKSudakov_"+observable).c_str());
    zw_ewkqcd_unc = (TH1*) inputFile->FindObjectAny(("ZW_QCDEWK_"+observable).c_str());

    TH1* zw_qcdscale_unc_flip = flipHisto(zw_qcdscale_unc);
    TH1* zw_nloewk_unc_flip = flipHisto(zw_nloewk_unc);
    TH1* zw_sudewk_unc_flip = flipHisto(zw_sudewk_unc);
    TH1* zw_ewkqcd_unc_flip = flipHisto(zw_ewkqcd_unc);

    zw_qcdscale_unc->SetLineColor(kGreen+1);
    zw_qcdscale_unc->SetLineWidth(2);
    zw_qcdscale_unc_flip->SetLineColor(kGreen+1);
    zw_qcdscale_unc_flip->SetLineWidth(2);
    zw_qcdscale_unc->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
    zw_qcdscale_unc->GetXaxis()->SetTitleOffset(1.1);
    zw_qcdscale_unc->GetXaxis()->SetTitleSize(0.05);
    zw_qcdscale_unc->GetYaxis()->SetTitle("Z/W Variation/Nominal");
    zw_qcdscale_unc->GetYaxis()->SetTitleSize(0.045);
    zw_qcdscale_unc->GetYaxis()->SetTitleOffset(1.15);

    zw_qcdscale_unc->GetYaxis()->SetRangeUser(-0.06,0.06);
    zw_qcdscale_unc->Draw("hist");
    zw_qcdscale_unc_flip->Draw("hist same");
    
    zw_nloewk_unc_flip->SetLineColor(kRed);
    zw_nloewk_unc_flip->SetLineWidth(2);
    zw_nloewk_unc->SetLineColor(kRed);
    zw_nloewk_unc->SetLineWidth(2);
    
    zw_nloewk_unc->Draw("hist same");
    zw_nloewk_unc_flip->Draw("hist same");

    zw_sudewk_unc_flip->SetLineColor(kBlack);
    zw_sudewk_unc_flip->SetLineWidth(2);
    zw_sudewk_unc->SetLineColor(kBlack);
    zw_sudewk_unc->SetLineWidth(2);

    zw_sudewk_unc->Draw("hist same");
    zw_sudewk_unc_flip->Draw("hist same");

    zw_ewkqcd_unc_flip->SetLineColor(kBlue);
    zw_ewkqcd_unc_flip->SetLineWidth(2);
    zw_ewkqcd_unc->SetLineColor(kBlue);
    zw_ewkqcd_unc->SetLineWidth(2);

    zw_ewkqcd_unc->Draw("hist same");
    zw_ewkqcd_unc_flip->Draw("hist same");
    
    TLegend leg (0.2,0.8,0.7,0.9);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.SetNColumns(2);
    leg.AddEntry(zw_qcdscale_unc,"QCD #mu_{r},#mu_{f}","L");
    leg.AddEntry(zw_nloewk_unc,"NLO-EWK","L");
    leg.AddEntry(zw_sudewk_unc,"NNLO Sudakov","L");
    leg.AddEntry(zw_ewkqcd_unc,"QCD-EWK Mix","L");
    leg.Draw("same");

    canvas->SaveAs((outputDIR+"/thoeryUnc_ZW_ratio.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/thoeryUnc_ZW_ratio.pdf").c_str(),"pdf");

    TFile* theory_unc = new TFile("theory_unc_ZW.root","RECREATE");
    zw_qcdscale_unc->Write();
    zw_nloewk_unc->Write();
    zw_sudewk_unc->Write();
    zw_ewkqcd_unc->Write();
      
  }


  TH1* zg_fac1_unc = (TH1*) inputFile->FindObjectAny(("ZG_FactScale1_"+observable).c_str());
  TH1* zg_fac2_unc = (TH1*) inputFile->FindObjectAny(("ZG_FactScale2_"+observable).c_str());
  TH1* zg_ren1_unc = (TH1*) inputFile->FindObjectAny(("ZG_RenScale1_"+observable).c_str());
  TH1* zg_ren2_unc = (TH1*) inputFile->FindObjectAny(("ZG_RenScale2_"+observable).c_str());
  TH1* zg_pdf_unc = (TH1*) inputFile->FindObjectAny(("ZG_PDF_"+observable).c_str());
  TH1* zg_pfp_unc = (TH1*) inputFile->FindObjectAny(("ZG_Footprint_"+observable).c_str());
  TH1* zg_ewk_unc = (TH1*) inputFile->FindObjectAny(("ZG_EWK_"+observable).c_str());

  
  TH1* zg_fac_unc = (TH1*) zg_fac1_unc->Clone("zg_fac_unc");
  for(int iBin = 1; iBin <= zg_fac_unc->GetNbinsX(); iBin++)
    zg_fac_unc->SetBinContent(iBin,sqrt(zg_fac1_unc->GetBinContent(iBin)*zg_fac1_unc->GetBinContent(iBin)+zg_fac2_unc->GetBinContent(iBin)*zg_fac2_unc->GetBinContent(iBin)));
  
  TH1* zg_ren_unc = (TH1*) zg_ren1_unc->Clone("zg_ren_unc");
  for(int iBin = 1; iBin <= zg_ren_unc->GetNbinsX(); iBin++)
    zg_ren_unc->SetBinContent(iBin,sqrt(zg_ren1_unc->GetBinContent(iBin)*zg_ren1_unc->GetBinContent(iBin)+zg_ren2_unc->GetBinContent(iBin)*zg_ren2_unc->GetBinContent(iBin)));
  
  TH1* zg_fac_unc_flip = flipHisto(zg_fac_unc);
  TH1* zg_ren_unc_flip = flipHisto(zg_ren_unc);
  TH1* zg_pdf_unc_flip = flipHisto(zg_pdf_unc);
  TH1* zg_ewk_unc_flip = flipHisto(zg_ewk_unc);
  
  zg_ren_unc->SetLineColor(kBlack);
  zg_ren_unc->SetLineWidth(2);
  zg_ren_unc_flip->SetLineColor(kBlack);
  zg_ren_unc_flip->SetLineWidth(2);
  zg_ren_unc->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
  zg_ren_unc->GetXaxis()->SetTitleOffset(1.1);
  zg_ren_unc->GetXaxis()->SetTitleSize(0.05);
  zg_ren_unc->GetYaxis()->SetTitle("Z/#gamma Variation/Nominal");
  zg_ren_unc->GetYaxis()->SetTitleSize(0.045);
  zg_ren_unc->GetYaxis()->SetTitleOffset(1.15);
  
  zg_ren_unc->GetYaxis()->SetRangeUser(-0.25,0.25);
  zg_ren_unc->Draw("hist");
  zg_ren_unc_flip->Draw("hist same");
  
  zg_fac_unc_flip->SetLineColor(kRed);
  zg_fac_unc_flip->SetLineWidth(2);
  zg_fac_unc->SetLineColor(kRed);
  zg_fac_unc->SetLineWidth(2);
  
  zg_fac_unc->Draw("hist same");
  zg_fac_unc_flip->Draw("hist same");
  
  zg_pdf_unc_flip->SetLineColor(kGreen+1);
  zg_pdf_unc_flip->SetLineWidth(2);
  zg_pdf_unc->SetLineColor(kGreen+1);
  zg_pdf_unc->SetLineWidth(2);
  
  zg_pdf_unc->Draw("hist same");
  zg_pdf_unc_flip->Draw("hist same");
  
  zg_ewk_unc_flip->SetLineColor(kBlue);
  zg_ewk_unc_flip->SetLineWidth(2);
  zg_ewk_unc->SetLineColor(kBlue);
  zg_ewk_unc->SetLineWidth(2);
  
  zg_ewk_unc->Draw("hist same");
  zg_ewk_unc_flip->Draw("hist same");
  
  TLegend leg (0.2,0.8,0.6,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.SetNColumns(2);  
  leg.AddEntry(zg_ren_unc,"QCD #mu_{r}","L");
  leg.AddEntry(zg_fac_unc,"QCD #mu_{f}","L");
  leg.AddEntry(zg_pdf_unc,"PDF","L");
  leg.AddEntry(zg_ewk_unc,"NLO-EWK","L");
  leg.Draw("same");
  
  canvas->SaveAs((outputDIR+"/thoeryUnc_ZG_ratio.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/thoeryUnc_ZG_ratio.pdf").c_str(),"pdf");
    

  

}
