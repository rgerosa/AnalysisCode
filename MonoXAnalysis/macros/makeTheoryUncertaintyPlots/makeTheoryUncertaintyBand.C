#include "../CMS_lumi.h"

TH1* flipHisto(TH1* histo){
  TH1* histo2 = (TH1*) histo->Clone(Form("%s_flip",histo->GetName()));
  for(int iBin =1; iBin <= histo2->GetNbinsX(); iBin++)
    histo2->SetBinContent(iBin,-histo->GetBinContent(iBin));
  return histo2;
}


void makeTheoryUncertaintyBand(string inputFileName, string outputDIR , bool useNewTheoryUncertainties, string observable, string observableLatex, bool addWgamma = false, bool twoSided = true){

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
  TH1* zw_qcdshape_unc = NULL;
  TH1* zw_qcdproc_unc = NULL;
  TH1* zw_nnloewk_unc = NULL;
  TH1* zw_nnlomiss1_unc = NULL;
  TH1* zw_nnlomiss2_unc = NULL;
  TH1* zw_sudakov1_unc = NULL;
  TH1* zw_sudakov2_unc = NULL;
  TH1* zw_mix_unc = NULL;

  TCanvas* canvas = new TCanvas("canvas","",600,600);
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
    zw_ren_unc->GetXaxis()->SetTitle(observableLatex.c_str());
    zw_ren_unc->GetXaxis()->SetTitleOffset(1.15);
    zw_ren_unc->GetXaxis()->SetTitleSize(0.05);
    zw_ren_unc->GetYaxis()->SetTitle("Z/W Variation/Nominal");
    zw_ren_unc->GetYaxis()->SetTitleSize(0.042);
    zw_ren_unc->GetYaxis()->SetTitleOffset(1.25);

    zw_ren_unc->GetYaxis()->SetRangeUser(-0.15,0.25);
    zw_ren_unc->Draw("hist");
    if(twoSided)
      zw_ren_unc_flip->Draw("hist same");
    
    zw_fac_unc_flip->SetLineColor(kRed);
    zw_fac_unc_flip->SetLineWidth(2);
    zw_fac_unc->SetLineColor(kRed);
    zw_fac_unc->SetLineWidth(2);
    
    zw_fac_unc->Draw("hist same");
    if(twoSided)
      zw_fac_unc_flip->Draw("hist same");

    zw_pdf_unc_flip->SetLineColor(kGreen+1);
    zw_pdf_unc_flip->SetLineWidth(2);
    zw_pdf_unc->SetLineColor(kGreen+1);
    zw_pdf_unc->SetLineWidth(2);

    zw_pdf_unc->Draw("hist same");
    if(twoSided)
      zw_pdf_unc_flip->Draw("hist same");

    zw_ewk_unc_flip->SetLineColor(kBlue);
    zw_ewk_unc_flip->SetLineWidth(2);
    zw_ewk_unc->SetLineColor(kBlue);
    zw_ewk_unc->SetLineWidth(2);

    zw_ewk_unc->Draw("hist same");
    if(twoSided)
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

    canvas->SaveAs((outputDIR+"/theoryUnc_ZW_ratio.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/theoryUnc_ZW_ratio.pdf").c_str(),"pdf");

  }
  else{

    zw_qcdscale_unc  = (TH1*) inputFile->FindObjectAny(("ZW_QCDScale_"+observable).c_str());
    zw_qcdshape_unc  = (TH1*) inputFile->FindObjectAny(("ZW_QCDShape_"+observable).c_str());
    zw_qcdproc_unc   = (TH1*) inputFile->FindObjectAny(("ZW_QCDProcess_"+observable).c_str());
    zw_nnloewk_unc   = (TH1*) inputFile->FindObjectAny(("ZW_NNLOEWK_"+observable).c_str());
    zw_nnlomiss1_unc = (TH1*) inputFile->FindObjectAny(("ZW_NNLOMiss1_"+observable).c_str());
    zw_nnlomiss2_unc = (TH1*) inputFile->FindObjectAny(("ZW_NNLOMiss2_"+observable).c_str());
    zw_sudakov1_unc  = (TH1*) inputFile->FindObjectAny(("ZW_Sudakov1_"+observable).c_str());
    zw_sudakov2_unc  = (TH1*) inputFile->FindObjectAny(("ZW_Sudakov2_"+observable).c_str());
    zw_mix_unc       = (TH1*) inputFile->FindObjectAny(("ZW_MIX_"+observable).c_str());
    zw_pdf_unc  = (TH1*) inputFile->FindObjectAny(("ZW_PDF_"+observable).c_str());

    TH1* zw_qcdscale_unc_flip = flipHisto(zw_qcdscale_unc);
    TH1* zw_qcdshape_unc_flip = flipHisto(zw_qcdshape_unc);
    TH1* zw_qcdproc_unc_flip = flipHisto(zw_qcdproc_unc);
    TH1* zw_nnloewk_unc_flip = flipHisto(zw_nnloewk_unc);
    TH1* zw_nnlomiss1_unc_flip = flipHisto(zw_nnlomiss1_unc);
    TH1* zw_nnlomiss2_unc_flip = flipHisto(zw_nnlomiss2_unc);
    TH1* zw_sudakov1_unc_flip = flipHisto(zw_sudakov1_unc);
    TH1* zw_sudakov2_unc_flip = flipHisto(zw_sudakov2_unc);
    TH1* zw_mix_unc_flip = flipHisto(zw_mix_unc);
    TH1* zw_pdf_unc_flip = flipHisto(zw_pdf_unc);

    zw_qcdscale_unc->GetXaxis()->SetTitle(observableLatex.c_str());
    zw_qcdscale_unc->GetXaxis()->SetTitleOffset(1.15);
    zw_qcdscale_unc->GetXaxis()->SetTitleSize(0.05);
    zw_qcdscale_unc->GetYaxis()->SetTitle("Z/W Variation/Nominal");
    zw_qcdscale_unc->GetYaxis()->SetTitleSize(0.045);
    zw_qcdscale_unc->GetYaxis()->SetTitleOffset(1.15);
    
    zw_qcdscale_unc->GetYaxis()->SetRangeUser(-0.03,0.05);
    zw_qcdscale_unc->SetLineColor(kBlack);
    zw_qcdscale_unc->SetLineWidth(2);
    zw_qcdscale_unc_flip->SetLineColor(kBlack);
    zw_qcdscale_unc_flip->SetLineWidth(2);

    zw_qcdscale_unc->Draw("hist");
    if(twoSided)
      zw_qcdscale_unc_flip->Draw("hist same");
    CMS_lumi(canvas,"");

    zw_qcdshape_unc->SetLineColor(kRed);
    zw_qcdshape_unc->SetLineWidth(2);
    zw_qcdshape_unc_flip->SetLineColor(kRed);
    zw_qcdshape_unc_flip->SetLineWidth(2);
    zw_qcdshape_unc->Draw("hist same");
    if(twoSided)
      zw_qcdshape_unc_flip->Draw("hist same");

    zw_qcdproc_unc->SetLineColor(kBlue);
    zw_qcdproc_unc->SetLineWidth(2);
    zw_qcdproc_unc_flip->SetLineColor(kBlue);
    zw_qcdproc_unc_flip->SetLineWidth(2);
    zw_qcdproc_unc->Draw("hist same");
    if(twoSided)
      zw_qcdproc_unc_flip->Draw("hist same");

    zw_pdf_unc_flip->SetLineColor(kCyan);
    zw_pdf_unc_flip->SetLineWidth(2);
    zw_pdf_unc->SetLineColor(kCyan);
    zw_pdf_unc->SetLineWidth(2);
    zw_pdf_unc->Draw("hist same");
    if(twoSided)
      zw_pdf_unc_flip->Draw("hist same");

    zw_mix_unc_flip->SetLineColor(kOrange+1);
    zw_mix_unc_flip->SetLineWidth(2);
    zw_mix_unc->SetLineColor(kOrange+1);
    zw_mix_unc->SetLineWidth(2);
    zw_mix_unc->Draw("hist same");
    if(twoSided)
      zw_mix_unc_flip->Draw("hist same");

    TLegend leg (0.4,0.6,0.7,0.9);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(zw_qcdscale_unc,"QCD #mu_{r},#mu_{f}","L");
    leg.AddEntry(zw_qcdshape_unc,"QCD Shape","L");
    leg.AddEntry(zw_qcdproc_unc,"QCD Process","L");
    leg.AddEntry(zw_mix_unc,"QCD-EWK Mix","L");
    leg.AddEntry(zw_pdf_unc,"PDF","L");
    leg.Draw("same");

    canvas->SaveAs((outputDIR+"/theoryUnc_ZW_ratio_qcdPart.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/theoryUnc_ZW_ratio_qcdPart.pdf").c_str(),"pdf");


    zw_nnloewk_unc->GetXaxis()->SetTitle(observableLatex.c_str());
    zw_nnloewk_unc->GetXaxis()->SetTitleOffset(1.15);
    zw_nnloewk_unc->GetXaxis()->SetTitleSize(0.05);
    zw_nnloewk_unc->GetYaxis()->SetTitle("Z/W Variation/Nominal");
    zw_nnloewk_unc->GetYaxis()->SetTitleSize(0.045);
    zw_nnloewk_unc->GetYaxis()->SetTitleOffset(1.15);
    zw_nnloewk_unc->GetYaxis()->SetRangeUser(-0.03,0.05);
    zw_nnloewk_unc->Draw("hist");
    zw_nnloewk_unc_flip->SetLineColor(kBlack);
    zw_nnloewk_unc_flip->SetLineWidth(2);
    zw_nnloewk_unc->SetLineColor(kBlack);
    zw_nnloewk_unc->SetLineWidth(2);    
    if(twoSided)
      zw_nnloewk_unc_flip->Draw("hist same");
    CMS_lumi(canvas,"");
    

    zw_sudakov1_unc_flip->SetLineColor(kRed);
    zw_sudakov1_unc_flip->SetLineWidth(2);
    zw_sudakov1_unc->SetLineColor(kRed);
    zw_sudakov1_unc->SetLineWidth(2);
    zw_sudakov1_unc->Draw("hist same");
    if(twoSided)
      zw_sudakov1_unc_flip->Draw("hist same");

    zw_sudakov2_unc_flip->SetLineColor(kCyan);
    zw_sudakov2_unc_flip->SetLineWidth(2);
    zw_sudakov2_unc->SetLineColor(kCyan);
    zw_sudakov2_unc->SetLineWidth(2);
    zw_sudakov2_unc->Draw("hist same");
    if(twoSided)
      zw_sudakov2_unc_flip->Draw("hist same");

    zw_nnlomiss1_unc_flip->SetLineColor(kBlue);
    zw_nnlomiss1_unc_flip->SetLineWidth(2);
    zw_nnlomiss1_unc->SetLineColor(kBlue);
    zw_nnlomiss1_unc->SetLineWidth(2);
    zw_nnlomiss1_unc->Draw("hist same");
    if(twoSided)
      zw_nnlomiss1_unc_flip->Draw("hist same");

    zw_nnlomiss2_unc_flip->SetLineColor(kOrange+1);
    zw_nnlomiss2_unc_flip->SetLineWidth(2);
    zw_nnlomiss2_unc->SetLineColor(kOrange+1);
    zw_nnlomiss2_unc->SetLineWidth(2);
    zw_nnlomiss2_unc->Draw("hist same");
    if(twoSided)
      zw_nnlomiss2_unc_flip->Draw("hist same");

    leg.Clear();
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(zw_nnloewk_unc,"N^{3}LO","L");
    leg.AddEntry(zw_sudakov1_unc,"Z+jets Sudakov","L");
    leg.AddEntry(zw_sudakov2_unc,"W+jets Sudakov","L");
    leg.AddEntry(zw_nnlomiss1_unc,"Z+jets NNLO Miss","L");
    leg.AddEntry(zw_nnlomiss2_unc,"W+jets NNLO Miss","L");
    leg.Draw("same");

    canvas->SaveAs((outputDIR+"/theoryUnc_ZW_ratio_ewkPart.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/theoryUnc_ZW_ratio_ewkPart.pdf").c_str(),"pdf");


    TFile* theory_unc = new TFile((outputDIR+"/theory_unc_ZW.root").c_str(),"RECREATE");
    zw_qcdscale_unc->Write();
    zw_qcdshape_unc->Write();
    zw_qcdproc_unc->Write();
    zw_nnloewk_unc->Write();
    zw_sudakov1_unc->Write();
    zw_sudakov2_unc->Write();
    zw_nnlomiss1_unc->Write();
    zw_nnlomiss2_unc->Write();
    zw_mix_unc->Write();
    zw_pdf_unc->Write();
    theory_unc->Close();
      
  }

  /////// Z-gamma part

  TH1* zg_fac1_unc = NULL;
  TH1* zg_fac2_unc = NULL;
  TH1* zg_ren1_unc = NULL;
  TH1* zg_ren2_unc = NULL;
  TH1* zg_pdf_unc = NULL;
  TH1* zg_ewk_unc = NULL;
  TH1* zg_qcdscale_unc = NULL;
  TH1* zg_qcdshape_unc = NULL;
  TH1* zg_qcdproc_unc = NULL;
  TH1* zg_nnloewk_unc = NULL;
  TH1* zg_nnlomiss1_unc = NULL;
  TH1* zg_nnlomiss2_unc = NULL;
  TH1* zg_sudakov1_unc = NULL;
  TH1* zg_sudakov2_unc = NULL;
  TH1* zg_mix_unc = NULL;

  if(not useNewTheoryUncertainties){
    zg_fac1_unc = (TH1*) inputFile->FindObjectAny(("ZG_FactScale1_"+observable).c_str());
    zg_fac2_unc = (TH1*) inputFile->FindObjectAny(("ZG_FactScale2_"+observable).c_str());
    zg_ren1_unc = (TH1*) inputFile->FindObjectAny(("ZG_RenScale1_"+observable).c_str());
    zg_ren2_unc = (TH1*) inputFile->FindObjectAny(("ZG_RenScale2_"+observable).c_str());
    zg_pdf_unc  = (TH1*) inputFile->FindObjectAny(("ZG_PDF_"+observable).c_str());
    zg_ewk_unc  = (TH1*) inputFile->FindObjectAny(("ZG_EWK_"+observable).c_str());

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
    zg_ren_unc->GetXaxis()->SetTitle(observableLatex.c_str());
    zg_ren_unc->GetXaxis()->SetTitleOffset(1.15);
    zg_ren_unc->GetXaxis()->SetTitleSize(0.05);
    zg_ren_unc->GetYaxis()->SetTitle("Z/#gamma Variation/Nominal");
    zg_ren_unc->GetYaxis()->SetTitleSize(0.042);
    zg_ren_unc->GetYaxis()->SetTitleOffset(1.25);

    zg_ren_unc->GetYaxis()->SetRangeUser(-0.20,0.35);
    zg_ren_unc->Draw("hist");
    if(twoSided)
      zg_ren_unc_flip->Draw("hist same");
    
    zg_fac_unc_flip->SetLineColor(kRed);
    zg_fac_unc_flip->SetLineWidth(2);
    zg_fac_unc->SetLineColor(kRed);
    zg_fac_unc->SetLineWidth(2);
    
    zg_fac_unc->Draw("hist same");
    if(twoSided)
      zg_fac_unc_flip->Draw("hist same");
      
    zg_pdf_unc_flip->SetLineColor(kGreen+1);
    zg_pdf_unc_flip->SetLineWidth(2);
    zg_pdf_unc->SetLineColor(kGreen+1);
    zg_pdf_unc->SetLineWidth(2);

    zg_pdf_unc->Draw("hist same");
    if(twoSided)
      zg_pdf_unc_flip->Draw("hist same");

    zg_ewk_unc_flip->SetLineColor(kBlue);
    zg_ewk_unc_flip->SetLineWidth(2);
    zg_ewk_unc->SetLineColor(kBlue);
    zg_ewk_unc->SetLineWidth(2);

    zg_ewk_unc->Draw("hist same");
    if(twoSided)
      zg_ewk_unc_flip->Draw("hist same");
    
    TLegend leg (0.7,0.7,0.9,0.9);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(zg_ren_unc,"Ren-Scale #mu_{r}","L");
    leg.AddEntry(zg_fac_unc,"Fact-Scale #mu_{f}","L");
    leg.AddEntry(zg_pdf_unc,"PDF","L");
    leg.AddEntry(zg_ewk_unc,"NLO-EWK","L");
    leg.Draw("same");

    canvas->SaveAs((outputDIR+"/theoryUnc_ZG_ratio.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/theoryUnc_ZG_ratio.pdf").c_str(),"pdf");

  }
  else{

    zg_qcdscale_unc  = (TH1*) inputFile->FindObjectAny(("ZG_QCDScale_"+observable).c_str());
    zg_qcdshape_unc  = (TH1*) inputFile->FindObjectAny(("ZG_QCDShape_"+observable).c_str());
    zg_qcdproc_unc   = (TH1*) inputFile->FindObjectAny(("ZG_QCDProcess_"+observable).c_str());
    zg_nnloewk_unc   = (TH1*) inputFile->FindObjectAny(("ZG_NNLOEWK_"+observable).c_str());
    zg_nnlomiss1_unc = (TH1*) inputFile->FindObjectAny(("ZG_NNLOMiss1_"+observable).c_str());
    zg_nnlomiss2_unc = (TH1*) inputFile->FindObjectAny(("ZG_NNLOMiss2_"+observable).c_str());
    zg_sudakov1_unc  = (TH1*) inputFile->FindObjectAny(("ZG_Sudakov1_"+observable).c_str());
    zg_sudakov2_unc  = (TH1*) inputFile->FindObjectAny(("ZG_Sudakov2_"+observable).c_str());
    zg_mix_unc       = (TH1*) inputFile->FindObjectAny(("ZG_MIX_"+observable).c_str());
    zg_pdf_unc       = (TH1*) inputFile->FindObjectAny(("ZG_PDF_"+observable).c_str());

    TH1* zg_qcdscale_unc_flip = flipHisto(zg_qcdscale_unc);
    TH1* zg_qcdshape_unc_flip = flipHisto(zg_qcdshape_unc);
    TH1* zg_qcdproc_unc_flip = flipHisto(zg_qcdproc_unc);
    TH1* zg_nnloewk_unc_flip = flipHisto(zg_nnloewk_unc);
    TH1* zg_nnlomiss1_unc_flip = flipHisto(zg_nnlomiss1_unc);
    TH1* zg_nnlomiss2_unc_flip = flipHisto(zg_nnlomiss2_unc);
    TH1* zg_sudakov1_unc_flip = flipHisto(zg_sudakov1_unc);
    TH1* zg_sudakov2_unc_flip = flipHisto(zg_sudakov2_unc);
    TH1* zg_mix_unc_flip = flipHisto(zg_mix_unc);
    TH1* zg_pdf_unc_flip = flipHisto(zg_pdf_unc);

    zg_qcdscale_unc->GetXaxis()->SetTitle(observableLatex.c_str());
    zg_qcdscale_unc->GetXaxis()->SetTitleOffset(1.15);
    zg_qcdscale_unc->GetXaxis()->SetTitleSize(0.05);
    zg_qcdscale_unc->GetYaxis()->SetTitle("Z/#gamma Variation/Nominal");
    zg_qcdscale_unc->GetYaxis()->SetTitleSize(0.045);
    zg_qcdscale_unc->GetYaxis()->SetTitleOffset(1.15);
    
    zg_qcdscale_unc->GetYaxis()->SetRangeUser(-0.04,0.10);
    zg_qcdscale_unc->SetLineColor(kBlack);
    zg_qcdscale_unc->SetLineWidth(2);
    zg_qcdscale_unc_flip->SetLineColor(kBlack);
    zg_qcdscale_unc_flip->SetLineWidth(2);

    zg_qcdscale_unc->Draw("hist");
    if(twoSided){
      zg_qcdscale_unc_flip->SetLineStyle(2);
      zg_qcdscale_unc_flip->Draw("hist same");
    }
    CMS_lumi(canvas,"");

    zg_qcdshape_unc->SetLineColor(kRed);
    zg_qcdshape_unc->SetLineWidth(2);
    zg_qcdshape_unc_flip->SetLineColor(kRed);
    zg_qcdshape_unc_flip->SetLineWidth(2);
    zg_qcdshape_unc->Draw("hist same");
    if(twoSided)
      zg_qcdshape_unc_flip->Draw("hist same");

    zg_qcdproc_unc->SetLineColor(kBlue);
    zg_qcdproc_unc->SetLineWidth(2);
    zg_qcdproc_unc_flip->SetLineColor(kBlue);
    zg_qcdproc_unc_flip->SetLineWidth(2);
    zg_qcdproc_unc->Draw("hist same");
    if(twoSided){
      zg_qcdproc_unc_flip->SetLineStyle(2);
      zg_qcdproc_unc_flip->Draw("hist same");
    }

    zg_pdf_unc_flip->SetLineColor(kCyan);
    zg_pdf_unc_flip->SetLineWidth(2);
    zg_pdf_unc->SetLineColor(kCyan);
    zg_pdf_unc->SetLineWidth(2);
    zg_pdf_unc->Draw("hist same");
    if(twoSided)
      zg_pdf_unc_flip->Draw("hist same");

    zg_mix_unc_flip->SetLineColor(kOrange+1);
    zg_mix_unc_flip->SetLineWidth(2);
    zg_mix_unc->SetLineColor(kOrange+1);
    zg_mix_unc->SetLineWidth(2);
    zg_mix_unc->Draw("hist same");
    if(twoSided)
      zg_mix_unc_flip->Draw("hist same");

    TLegend leg (0.4,0.6,0.7,0.9);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(zg_qcdscale_unc,"QCD #mu_{r},#mu_{f}","L");
    leg.AddEntry(zg_qcdshape_unc,"QCD Shape","L");
    leg.AddEntry(zg_qcdproc_unc,"QCD Process","L");
    leg.AddEntry(zg_mix_unc,"QCD-EWK Mix","L");
    leg.AddEntry(zg_pdf_unc,"PDF","L");
    leg.Draw("same");

    canvas->SaveAs((outputDIR+"/theoryUnc_ZG_ratio_qcdPart.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/theoryUnc_ZG_ratio_qcdPart.pdf").c_str(),"pdf");


    zg_nnloewk_unc->GetXaxis()->SetTitle(observableLatex.c_str());
    zg_nnloewk_unc->GetXaxis()->SetTitleOffset(1.15);
    zg_nnloewk_unc->GetXaxis()->SetTitleSize(0.05);
    zg_nnloewk_unc->GetYaxis()->SetTitle("Z/#gamma Variation/Nominal");
    zg_nnloewk_unc->GetYaxis()->SetTitleSize(0.045);
    zg_nnloewk_unc->GetYaxis()->SetTitleOffset(1.15);
    zg_nnloewk_unc->GetYaxis()->SetRangeUser(-0.04,0.07);
    zg_nnloewk_unc->Draw("hist");
    zg_nnloewk_unc_flip->SetLineColor(kBlack);
    zg_nnloewk_unc_flip->SetLineWidth(2);
    zg_nnloewk_unc->SetLineColor(kBlack);
    zg_nnloewk_unc->SetLineWidth(2);    
    if(twoSided)
      zg_nnloewk_unc_flip->Draw("hist same");
    CMS_lumi(canvas,"");
    

    zg_sudakov1_unc_flip->SetLineColor(kRed);
    zg_sudakov1_unc_flip->SetLineWidth(2);
    zg_sudakov1_unc->SetLineColor(kRed);
    zg_sudakov1_unc->SetLineWidth(2);
    zg_sudakov1_unc->Draw("hist same");
    if(twoSided)
      zg_sudakov1_unc_flip->Draw("hist same");

    zg_sudakov2_unc_flip->SetLineColor(kCyan);
    zg_sudakov2_unc_flip->SetLineWidth(2);
    zg_sudakov2_unc->SetLineColor(kCyan);
    zg_sudakov2_unc->SetLineWidth(2);
    zg_sudakov2_unc->Draw("hist same");
    if(twoSided)
      zg_sudakov2_unc_flip->Draw("hist same");

    zg_nnlomiss1_unc_flip->SetLineColor(kBlue);
    zg_nnlomiss1_unc_flip->SetLineWidth(2);
    zg_nnlomiss1_unc->SetLineColor(kBlue);
    zg_nnlomiss1_unc->SetLineWidth(2);
    zg_nnlomiss1_unc->Draw("hist same");
    if(twoSided)
      zg_nnlomiss1_unc_flip->Draw("hist same");

    zg_nnlomiss2_unc_flip->SetLineColor(kOrange+1);
    zg_nnlomiss2_unc_flip->SetLineWidth(2);
    zg_nnlomiss2_unc->SetLineColor(kOrange+1);
    zg_nnlomiss2_unc->SetLineWidth(2);
    zg_nnlomiss2_unc->Draw("hist same");
    if(twoSided)
      zg_nnlomiss2_unc_flip->Draw("hist same");

    leg.Clear();
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry(zg_nnloewk_unc,"N^{3}LO","L");
    leg.AddEntry(zg_sudakov1_unc,"Z+jets Sudakov","L");
    leg.AddEntry(zg_sudakov2_unc,"#gamma+jets Sudakov","L");
    leg.AddEntry(zg_nnlomiss1_unc,"Z+jets NNLO Miss","L");
    leg.AddEntry(zg_nnlomiss2_unc,"#gamma+jets NNLO Miss","L");
    leg.Draw("same");

    canvas->SaveAs((outputDIR+"/theoryUnc_ZG_ratio_ewkPart.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/theoryUnc_ZG_ratio_ewkPart.pdf").c_str(),"pdf");


    TFile* theory_unc = new TFile((outputDIR+"/theory_unc_ZG.root").c_str(),"RECREATE");
    zg_qcdscale_unc->Write();
    zg_qcdshape_unc->Write();
    zg_qcdproc_unc->Write();
    zg_nnloewk_unc->Write();
    zg_sudakov1_unc->Write();
    zg_sudakov2_unc->Write();
    zg_nnlomiss1_unc->Write();
    zg_nnlomiss2_unc->Write();
    zg_mix_unc->Write();
    zg_pdf_unc->Write();
    theory_unc->Close();
      
  }

  if(addWgamma){
    /////// W-gamma part
    TH1* wg_fac1_unc = NULL;
    TH1* wg_fac2_unc = NULL;
    TH1* wg_ren1_unc = NULL;
    TH1* wg_ren2_unc = NULL;
    TH1* wg_pdf_unc = NULL;
    TH1* wg_ewk_unc = NULL;
    TH1* wg_qcdscale_unc = NULL;
    TH1* wg_qcdshape_unc = NULL;
    TH1* wg_qcdproc_unc = NULL;
    TH1* wg_nnloewk_unc = NULL;
    TH1* wg_nnlomiss1_unc = NULL;
    TH1* wg_nnlomiss2_unc = NULL;
    TH1* wg_sudakov1_unc = NULL;
    TH1* wg_sudakov2_unc = NULL;
    TH1* wg_mix_unc = NULL;
    
    if(not useNewTheoryUncertainties){
      wg_fac1_unc = (TH1*) inputFile->FindObjectAny(("WG_FactScale1_"+observable).c_str());
      wg_fac2_unc = (TH1*) inputFile->FindObjectAny(("WG_FactScale2_"+observable).c_str());
      wg_ren1_unc = (TH1*) inputFile->FindObjectAny(("WG_RenScale1_"+observable).c_str());
      wg_ren2_unc = (TH1*) inputFile->FindObjectAny(("WG_RenScale2_"+observable).c_str());
      wg_pdf_unc  = (TH1*) inputFile->FindObjectAny(("WG_PDF_"+observable).c_str());
      wg_ewk_unc  = (TH1*) inputFile->FindObjectAny(("WG_EWK_"+observable).c_str());
      
      TH1* wg_fac_unc = (TH1*) wg_fac1_unc->Clone("wg_fac_unc");
      for(int iBin = 1; iBin <= wg_fac_unc->GetNbinsX(); iBin++)
	wg_fac_unc->SetBinContent(iBin,sqrt(wg_fac1_unc->GetBinContent(iBin)*wg_fac1_unc->GetBinContent(iBin)+wg_fac2_unc->GetBinContent(iBin)*wg_fac2_unc->GetBinContent(iBin)));
      
      TH1* wg_ren_unc = (TH1*) wg_ren1_unc->Clone("wg_ren_unc");
      for(int iBin = 1; iBin <= wg_ren_unc->GetNbinsX(); iBin++)
	wg_ren_unc->SetBinContent(iBin,sqrt(wg_ren1_unc->GetBinContent(iBin)*wg_ren1_unc->GetBinContent(iBin)+wg_ren2_unc->GetBinContent(iBin)*wg_ren2_unc->GetBinContent(iBin)));
      
      TH1* wg_fac_unc_flip = flipHisto(wg_fac_unc);
      TH1* wg_ren_unc_flip = flipHisto(wg_ren_unc);
      TH1* wg_pdf_unc_flip = flipHisto(wg_pdf_unc);
      TH1* wg_ewk_unc_flip = flipHisto(wg_ewk_unc);
      
      wg_ren_unc->SetLineColor(kBlack);
      wg_ren_unc->SetLineWidth(2);
      wg_ren_unc_flip->SetLineColor(kBlack);
      wg_ren_unc_flip->SetLineWidth(2);
      wg_ren_unc->GetXaxis()->SetTitle(observableLatex.c_str());
      wg_ren_unc->GetXaxis()->SetTitleOffset(1.15);
      wg_ren_unc->GetXaxis()->SetTitleSize(0.05);
      wg_ren_unc->GetYaxis()->SetTitle("W/#gamma Variation/Nominal");
      wg_ren_unc->GetYaxis()->SetTitleSize(0.042);
      wg_ren_unc->GetYaxis()->SetTitleOffset(1.25);
      
      wg_ren_unc->GetYaxis()->SetRangeUser(-0.20,0.35);
      wg_ren_unc->Draw("hist");
      if(twoSided)
	wg_ren_unc_flip->Draw("hist same");
    
      wg_fac_unc_flip->SetLineColor(kRed);
      wg_fac_unc_flip->SetLineWidth(2);
      wg_fac_unc->SetLineColor(kRed);
      wg_fac_unc->SetLineWidth(2);
      
      wg_fac_unc->Draw("hist same");
      if(twoSided)
	wg_fac_unc_flip->Draw("hist same");
      
      wg_pdf_unc_flip->SetLineColor(kGreen+1);
      wg_pdf_unc_flip->SetLineWidth(2);
      wg_pdf_unc->SetLineColor(kGreen+1);
      wg_pdf_unc->SetLineWidth(2);
      
      wg_pdf_unc->Draw("hist same");
      if(twoSided)
	wg_pdf_unc_flip->Draw("hist same");
      
      wg_ewk_unc_flip->SetLineColor(kBlue);
      wg_ewk_unc_flip->SetLineWidth(2);
      wg_ewk_unc->SetLineColor(kBlue);
      wg_ewk_unc->SetLineWidth(2);
      
      wg_ewk_unc->Draw("hist same");
      if(twoSided)
	wg_ewk_unc_flip->Draw("hist same");
      
      TLegend leg (0.7,0.7,0.9,0.9);
      leg.SetFillColor(0);
      leg.SetFillStyle(0);
      leg.SetBorderSize(0);
      leg.AddEntry(wg_ren_unc,"Ren-Scale #mu_{r}","L");
      leg.AddEntry(wg_fac_unc,"Fact-Scale #mu_{f}","L");
      leg.AddEntry(wg_pdf_unc,"PDF","L");
      leg.AddEntry(wg_ewk_unc,"NLO-EWK","L");
      leg.Draw("same");
      
      canvas->SaveAs((outputDIR+"/theoryUnc_WG_ratio.png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/theoryUnc_WG_ratio.pdf").c_str(),"pdf");    
    }
    else{
      
      wg_qcdscale_unc  = (TH1*) inputFile->FindObjectAny(("WG_QCDScale_"+observable).c_str());
      wg_qcdshape_unc  = (TH1*) inputFile->FindObjectAny(("WG_QCDShape_"+observable).c_str());
      wg_qcdproc_unc   = (TH1*) inputFile->FindObjectAny(("WG_QCDProcess_"+observable).c_str());
      wg_nnloewk_unc   = (TH1*) inputFile->FindObjectAny(("WG_NNLOEWK_"+observable).c_str());
      wg_nnlomiss1_unc = (TH1*) inputFile->FindObjectAny(("WG_NNLOMiss1_"+observable).c_str());
      wg_nnlomiss2_unc = (TH1*) inputFile->FindObjectAny(("WG_NNLOMiss2_"+observable).c_str());
      wg_sudakov1_unc  = (TH1*) inputFile->FindObjectAny(("WG_Sudakov1_"+observable).c_str());
      wg_sudakov2_unc  = (TH1*) inputFile->FindObjectAny(("WG_Sudakov2_"+observable).c_str());
      wg_mix_unc       = (TH1*) inputFile->FindObjectAny(("WG_MIX_"+observable).c_str());
      wg_pdf_unc       = (TH1*) inputFile->FindObjectAny(("WG_PDF_"+observable).c_str());
      
      TH1* wg_qcdscale_unc_flip = flipHisto(wg_qcdscale_unc);
      TH1* wg_qcdshape_unc_flip = flipHisto(wg_qcdshape_unc);
      TH1* wg_qcdproc_unc_flip = flipHisto(wg_qcdproc_unc);
      TH1* wg_nnloewk_unc_flip = flipHisto(wg_nnloewk_unc);
      TH1* wg_nnlomiss1_unc_flip = flipHisto(wg_nnlomiss1_unc);
      TH1* wg_nnlomiss2_unc_flip = flipHisto(wg_nnlomiss2_unc);
      TH1* wg_sudakov1_unc_flip = flipHisto(wg_sudakov1_unc);
      TH1* wg_sudakov2_unc_flip = flipHisto(wg_sudakov2_unc);
      TH1* wg_mix_unc_flip = flipHisto(wg_mix_unc);
      TH1* wg_pdf_unc_flip = flipHisto(wg_pdf_unc);

      wg_qcdscale_unc->GetXaxis()->SetTitle(observableLatex.c_str());
      wg_qcdscale_unc->GetXaxis()->SetTitleOffset(1.15);
      wg_qcdscale_unc->GetXaxis()->SetTitleSize(0.05);
      wg_qcdscale_unc->GetYaxis()->SetTitle("W/#gamma Variation/Nominal");
      wg_qcdscale_unc->GetYaxis()->SetTitleSize(0.045);
      wg_qcdscale_unc->GetYaxis()->SetTitleOffset(1.15);
      
      wg_qcdscale_unc->GetYaxis()->SetRangeUser(-0.04,0.07);
      wg_qcdscale_unc->SetLineColor(kBlack);
      wg_qcdscale_unc->SetLineWidth(2);
      wg_qcdscale_unc_flip->SetLineColor(kBlack);
      wg_qcdscale_unc_flip->SetLineWidth(2);
      
      wg_qcdscale_unc->Draw("hist");
      if(twoSided)
	wg_qcdscale_unc_flip->Draw("hist same");
      CMS_lumi(canvas,"");
      
      wg_qcdshape_unc->SetLineColor(kRed);
      wg_qcdshape_unc->SetLineWidth(2);
      wg_qcdshape_unc_flip->SetLineColor(kRed);
      wg_qcdshape_unc_flip->SetLineWidth(2);
      wg_qcdshape_unc->Draw("hist same");
      if(twoSided)
	wg_qcdshape_unc_flip->Draw("hist same");

      wg_qcdproc_unc->SetLineColor(kBlue);
      wg_qcdproc_unc->SetLineWidth(2);
      wg_qcdproc_unc_flip->SetLineColor(kBlue);
      wg_qcdproc_unc_flip->SetLineWidth(2);
      wg_qcdproc_unc->Draw("hist same");
      if(twoSided)
	wg_qcdproc_unc_flip->Draw("hist same");
      
      wg_pdf_unc_flip->SetLineColor(kCyan);
      wg_pdf_unc_flip->SetLineWidth(2);
      wg_pdf_unc->SetLineColor(kCyan);
      wg_pdf_unc->SetLineWidth(2);
      wg_pdf_unc->Draw("hist same");
      if(twoSided)
	wg_pdf_unc_flip->Draw("hist same");
      
      wg_mix_unc_flip->SetLineColor(kOrange+1);
      wg_mix_unc_flip->SetLineWidth(2);
      wg_mix_unc->SetLineColor(kOrange+1);
      wg_mix_unc->SetLineWidth(2);
      wg_mix_unc->Draw("hist same");
      if(twoSided)
	wg_mix_unc_flip->Draw("hist same");

      TLegend leg (0.4,0.6,0.7,0.9);
      leg.SetFillColor(0);
      leg.SetFillStyle(0);
      leg.SetBorderSize(0);
      leg.AddEntry(wg_qcdscale_unc,"QCD #mu_{r},#mu_{f}","L");
      leg.AddEntry(wg_qcdshape_unc,"QCD Shape","L");
      leg.AddEntry(wg_qcdproc_unc,"QCD Process","L");
      leg.AddEntry(wg_mix_unc,"QCD-EWK Mix","L");
      leg.AddEntry(wg_pdf_unc,"PDF","L");
      leg.Draw("same");
      
      canvas->SaveAs((outputDIR+"/theoryUnc_WG_ratio_qcdPart.png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/theoryUnc_WG_ratio_qcdPart.pdf").c_str(),"pdf");

      
      wg_nnloewk_unc->GetXaxis()->SetTitle(observableLatex.c_str());
      wg_nnloewk_unc->GetXaxis()->SetTitleOffset(1.15);
      wg_nnloewk_unc->GetXaxis()->SetTitleSize(0.05);
      wg_nnloewk_unc->GetYaxis()->SetTitle("W/#gamma Variation/Nominal");
      wg_nnloewk_unc->GetYaxis()->SetTitleSize(0.045);
      wg_nnloewk_unc->GetYaxis()->SetTitleOffset(1.15);
      wg_nnloewk_unc->GetYaxis()->SetRangeUser(-0.04,0.07);
      wg_nnloewk_unc->Draw("hist");
      wg_nnloewk_unc_flip->SetLineColor(kBlack);
      wg_nnloewk_unc_flip->SetLineWidth(2);
      wg_nnloewk_unc->SetLineColor(kBlack);
      wg_nnloewk_unc->SetLineWidth(2);    
      if(twoSided)
	wg_nnloewk_unc_flip->Draw("hist same");
      CMS_lumi(canvas,"");
      

      wg_sudakov1_unc_flip->SetLineColor(kRed);
      wg_sudakov1_unc_flip->SetLineWidth(2);
      wg_sudakov1_unc->SetLineColor(kRed);
      wg_sudakov1_unc->SetLineWidth(2);
      wg_sudakov1_unc->Draw("hist same");
      if(twoSided)
	wg_sudakov1_unc_flip->Draw("hist same");
      
      wg_sudakov2_unc_flip->SetLineColor(kCyan);
      wg_sudakov2_unc_flip->SetLineWidth(2);
      wg_sudakov2_unc->SetLineColor(kCyan);
      wg_sudakov2_unc->SetLineWidth(2);
      wg_sudakov2_unc->Draw("hist same");
      if(twoSided)
	wg_sudakov2_unc_flip->Draw("hist same");
      
      wg_nnlomiss1_unc_flip->SetLineColor(kBlue);
      wg_nnlomiss1_unc_flip->SetLineWidth(2);
      wg_nnlomiss1_unc->SetLineColor(kBlue);
      wg_nnlomiss1_unc->SetLineWidth(2);
      wg_nnlomiss1_unc->Draw("hist same");
      if(twoSided)
	wg_nnlomiss1_unc_flip->Draw("hist same");
      
      wg_nnlomiss2_unc_flip->SetLineColor(kOrange+1);
      wg_nnlomiss2_unc_flip->SetLineWidth(2);
      wg_nnlomiss2_unc->SetLineColor(kOrange+1);
      wg_nnlomiss2_unc->SetLineWidth(2);
      wg_nnlomiss2_unc->Draw("hist same");
      if(twoSided)
      wg_nnlomiss2_unc_flip->Draw("hist same");
      
      leg.Clear();
      leg.SetFillColor(0);
      leg.SetFillStyle(0);
      leg.SetBorderSize(0);
      leg.AddEntry(wg_nnloewk_unc,"N^{3}LO","L");
      leg.AddEntry(wg_sudakov1_unc,"Sudakov 1","L");
      leg.AddEntry(wg_sudakov2_unc,"Sudakov 2","L");
      leg.AddEntry(wg_nnlomiss1_unc,"NNLO miss 1","L");
      leg.AddEntry(wg_nnlomiss2_unc,"NNLO miss 2","L");
      leg.Draw("same");
      
      canvas->SaveAs((outputDIR+"/theoryUnc_WG_ratio_ewkPart.png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/theoryUnc_WG_ratio_ewkPart.pdf").c_str(),"pdf");
      
      
      TFile* theory_unc = new TFile((outputDIR+"/theory_unc_WG.root").c_str(),"RECREATE");
      wg_qcdscale_unc->Write();
      wg_qcdshape_unc->Write();
      wg_qcdproc_unc->Write();
      wg_nnloewk_unc->Write();
      wg_sudakov1_unc->Write();
      wg_sudakov2_unc->Write();
      wg_nnlomiss1_unc->Write();
      wg_nnlomiss2_unc->Write();
      wg_mix_unc->Write();
      wg_pdf_unc->Write();
      theory_unc->Close();
      
    }
  }
}
