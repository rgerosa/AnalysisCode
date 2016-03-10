#include "histoUtils.h"
#include "CMS_lumi.h"

void makeDataValidationPlots(string inputFileName, int category, string observable, string observableLatex){

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);

  // open the input file
  TFile* inputFile = TFile::Open(inputFileName.c_str());
  
  TH1* data_zmm = (TH1*) inputFile->Get(("datahistzmm_"+observable).c_str());
  TH1* data_zee = (TH1*) inputFile->Get(("datahistzee_"+observable).c_str());
  TH1* data_wen = (TH1*) inputFile->Get(("datahistwen_"+observable).c_str());
  TH1* data_wmn = (TH1*) inputFile->Get(("datahistwmn_"+observable).c_str());
  TH1* data_gam = (TH1*) inputFile->Get(("datahistgam_"+observable).c_str());

  // MC
  TH1* vllbkg_zmm = (TH1*) inputFile->Get(("vllbkghistzmm_"+observable).c_str());
  TH1* vlbkg_zmm = (TH1*) inputFile->Get(("vlbkghistzmm_"+observable).c_str());
  TH1* dbbkg_zmm = (TH1*) inputFile->Get(("dbkghistzmm_"+observable).c_str());
  TH1* ttbkg_zmm = (TH1*) inputFile->Get(("tbkghistzmm_"+observable).c_str());
  vllbkg_zmm->Add(vlbkg_zmm);
  vllbkg_zmm->Add(dbbkg_zmm);
  vllbkg_zmm->Add(ttbkg_zmm);
  
  TH1* vllbkg_zee = (TH1*) inputFile->Get(("vllbkghistzee_"+observable).c_str());
  TH1* vlbkg_zee = (TH1*) inputFile->Get(("vlbkghistzee_"+observable).c_str());
  TH1* dbbkg_zee = (TH1*) inputFile->Get(("dbkghistzee_"+observable).c_str());
  TH1* ttbkg_zee = (TH1*) inputFile->Get(("tbkghistzee_"+observable).c_str());
  vllbkg_zee->Add(vlbkg_zee);
  vllbkg_zee->Add(dbbkg_zee);
  vllbkg_zee->Add(ttbkg_zee);

  TH1* vlbkg_wen   = (TH1*) inputFile->Get(("vlbkghistwen_"+observable).c_str());
  TH1* vllbkg_wen  = (TH1*) inputFile->Get(("vllbkghistwen_"+observable).c_str());
  TH1* dbbkg_wen   = (TH1*) inputFile->Get(("dbkghistwen_"+observable).c_str());
  TH1* ttbkg_wen   = (TH1*) inputFile->Get(("tbkghistwen_"+observable).c_str());
  TH1* qbkg_wen    = (TH1*) inputFile->Get(("qbkghistwen_"+observable).c_str());
  vlbkg_wen->Add(vllbkg_wen);
  vlbkg_wen->Add(dbbkg_wen);
  vlbkg_wen->Add(ttbkg_wen);
  vlbkg_wen->Add(qbkg_wen);

  TH1* vlbkg_wmn   = (TH1*) inputFile->Get(("vlbkghistwmn_"+observable).c_str());
  TH1* vllbkg_wmn  = (TH1*) inputFile->Get(("vllbkghistwmn_"+observable).c_str());
  TH1* dbbkg_wmn   = (TH1*) inputFile->Get(("dbkghistwmn_"+observable).c_str());
  TH1* ttbkg_wmn   = (TH1*) inputFile->Get(("tbkghistwmn_"+observable).c_str());
  TH1* qbkg_wmn    = (TH1*) inputFile->Get(("qbkghistwmn_"+observable).c_str());
  vlbkg_wmn->Add(vllbkg_wmn);
  vlbkg_wmn->Add(dbbkg_wmn);
  vlbkg_wmn->Add(ttbkg_wmn);
  vlbkg_wmn->Add(qbkg_wmn);

  TH1* gbkg_gam   = (TH1*) inputFile->Get(("gbkghistgam_"+observable).c_str());
  TH1* qbkg_gam   = (TH1*) inputFile->Get(("qbkghistgam_"+observable).c_str());
  gbkg_gam->Add(qbkg_gam);

  //SYS Unc
  TH1*  ZG_ewk = (TH1*)inputFile->Get(("ZG_EWK_"+observable).c_str());
  TH1*  ZG_re1 = (TH1*)inputFile->Get(("ZG_RenScale1_"+observable).c_str());
  TH1*  ZG_re2 = (TH1*)inputFile->Get(("ZG_RenScale2_"+observable).c_str());
  TH1*  ZG_fa1 = (TH1*)inputFile->Get(("ZG_FactScale1_"+observable).c_str());
  TH1*  ZG_fa2 = (TH1*)inputFile->Get(("ZG_FactScale2_"+observable).c_str());
  TH1*  ZG_pdf = (TH1*)inputFile->Get(("ZG_PDF_"+observable).c_str());
  TH1*  ZG_fp  = (TH1*)inputFile->Get(("ZG_Footprint_"+observable).c_str());

  TH1*  ZW_ewk = (TH1*)inputFile->Get(("ZW_EWK_"+observable).c_str());
  TH1*  ZW_re1 = (TH1*)inputFile->Get(("ZW_RenScale1_"+observable).c_str());
  TH1*  ZW_re2 = (TH1*)inputFile->Get(("ZW_RenScale2_"+observable).c_str());
  TH1*  ZW_fa1 = (TH1*)inputFile->Get(("ZW_FactScale1_"+observable).c_str());
  TH1*  ZW_fa2 = (TH1*)inputFile->Get(("ZW_FactScale2_"+observable).c_str());
  TH1*  ZW_pdf = (TH1*)inputFile->Get(("ZW_PDF_"+observable).c_str());

  //Ratios Data
  TH1* ZGData_mm = (TH1*) data_zmm->Clone("ZGData_mm");
  ZGData_mm->Divide(data_gam);
  TH1* ZGData_ee = (TH1*) data_zee->Clone("ZGData_ee");
  ZGData_ee->Divide(data_gam);
  TH1* ZGData_ll = (TH1*) data_zmm->Clone("ZGData_ll");
  ZGData_ll->Add(data_zee);
  ZGData_ll->Divide(data_gam);
  
  TH1* ZWData_ee = (TH1*) data_zee->Clone("ZWData_e");
  ZWData_ee->Divide(data_wen);
  TH1* ZWData_mm = (TH1*) data_zmm->Clone("ZWData_m");
  ZWData_mm->Divide(data_wmn);
  TH1* ZWData_ll = (TH1*) data_zmm->Clone("ZWData_l");
  ZWData_ll->Add(data_zee);
  TH1* temp = (TH1*) data_wmn->Clone("temp"); 
  temp->Add(data_wen);
  ZWData_ll->Divide(temp);

  //Ratios MC
  TH1* ZGMC_mm = (TH1*) vllbkg_zmm->Clone("ZGMC_mm");
  ZGMC_mm->Divide(gbkg_gam);
  TH1* ZGMC_ee = (TH1*) vllbkg_zee->Clone("ZGMC_ee");
  ZGMC_ee->Divide(gbkg_gam);
  TH1* ZGMC_ll = (TH1*) vllbkg_zmm->Clone("ZGMC_ll");
  ZGMC_ll->Add(vllbkg_zee);
  ZGMC_ll->Divide(gbkg_gam);
  
  TH1* ZWMC_ee = (TH1*) vllbkg_zee->Clone("ZWMC_e");
  ZWMC_ee->Divide(vlbkg_wen);
  TH1* ZWMC_mm = (TH1*) vllbkg_zmm->Clone("ZWMC_m");
  ZWMC_mm->Divide(vlbkg_wmn);
  TH1* ZWMC_ll = (TH1*) vllbkg_zmm->Clone("ZWMC_l");
  ZWMC_ll->Add(vllbkg_zee);
  temp = (TH1*) vlbkg_wmn->Clone("temp"); 
  temp->Add(vlbkg_wen);
  ZWMC_ll->Divide(temp);

  //Add systematic uncertainties
  for(int iBin = 0; iBin < ZGMC_mm->GetNbinsX(); iBin++){
    double err = 0.;
    err += ZGMC_mm->GetBinError(iBin+1)*ZGMC_mm->GetBinError(iBin+1);
    err += pow(ZG_ewk->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZG_re1->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZG_re2->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZG_fa1->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZG_fa2->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZG_pdf->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZG_fp->GetBinContent(iBin+1)*ZGMC_mm->GetBinContent(iBin+1), 2);
    ZGMC_mm->SetBinError(iBin+1,sqrt(err));
  }
  
  for(int iBin = 0; iBin < ZGMC_ee->GetNbinsX(); iBin++){
    double err = 0.;
    err += ZGMC_ee->GetBinError(iBin+1)*ZGMC_ee->GetBinError(iBin+1);
    err += pow(ZG_ewk->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZG_re1->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZG_re2->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZG_fa1->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZG_fa2->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZG_pdf->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZG_fp->GetBinContent(iBin+1)*ZGMC_ee->GetBinContent(iBin+1), 2);
    ZGMC_ee->SetBinError(iBin+1,sqrt(err));
  }

  for(int iBin = 0; iBin < ZGMC_ll->GetNbinsX(); iBin++){
    double err = 0.;
    err += ZGMC_ll->GetBinError(iBin+1)*ZGMC_ll->GetBinError(iBin+1);
    err += pow(ZG_ewk->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZG_re1->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZG_re2->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZG_fa1->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZG_fa2->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZG_pdf->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZG_fp->GetBinContent(iBin+1)*ZGMC_ll->GetBinContent(iBin+1), 2);
    ZGMC_ll->SetBinError(iBin+1,sqrt(err));
  }

  for(int iBin = 0; iBin < ZWMC_mm->GetNbinsX(); iBin++){
    double err = 0.;
    err += ZWMC_mm->GetBinError(iBin+1)*ZWMC_mm->GetBinError(iBin+1);
    err += pow(ZW_ewk->GetBinContent(iBin+1)*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZW_re1->GetBinContent(iBin+1)*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZW_re2->GetBinContent(iBin+1)*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZW_fa1->GetBinContent(iBin+1)*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZW_fa2->GetBinContent(iBin+1)*ZWMC_mm->GetBinContent(iBin+1), 2);
    err += pow(ZW_pdf->GetBinContent(iBin+1)*ZWMC_mm->GetBinContent(iBin+1), 2);
    ZWMC_mm->SetBinError(iBin+1,sqrt(err));
  }
  
  for(int iBin = 0; iBin < ZWMC_ee->GetNbinsX(); iBin++){
    double err = 0.;
    err += ZWMC_ee->GetBinError(iBin+1)*ZWMC_ee->GetBinError(iBin+1);
    err += pow(ZW_ewk->GetBinContent(iBin+1)*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZW_re1->GetBinContent(iBin+1)*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZW_re2->GetBinContent(iBin+1)*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZW_fa1->GetBinContent(iBin+1)*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZW_fa2->GetBinContent(iBin+1)*ZWMC_ee->GetBinContent(iBin+1), 2);
    err += pow(ZW_pdf->GetBinContent(iBin+1)*ZWMC_ee->GetBinContent(iBin+1), 2);
    ZWMC_ee->SetBinError(iBin+1,sqrt(err));
  }

  for(int iBin = 0; iBin < ZWMC_ll->GetNbinsX(); iBin++){
    double err = 0.;
    err += ZWMC_ll->GetBinError(iBin+1)*ZWMC_ll->GetBinError(iBin+1);
    err += pow(ZW_ewk->GetBinContent(iBin+1)*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZW_re1->GetBinContent(iBin+1)*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZW_re2->GetBinContent(iBin+1)*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZW_fa1->GetBinContent(iBin+1)*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZW_fa2->GetBinContent(iBin+1)*ZWMC_ll->GetBinContent(iBin+1), 2);
    err += pow(ZW_pdf->GetBinContent(iBin+1)*ZWMC_ll->GetBinContent(iBin+1), 2);
    ZWMC_ll->SetBinError(iBin+1,sqrt(err));
  }

  // final plot
  TCanvas* canvas_ZG = new TCanvas("canvas_ZG","",600,700);
  canvas_ZG->SetTickx();
  canvas_ZG->SetTicky();
  canvas_ZG->cd();
  canvas_ZG->SetLeftMargin(0.11);
  TPad *pad1_ZG = new TPad("pad1_ZG","pad1_ZG",0,0.3,1,1);
  pad1_ZG->SetTickx();
  pad1_ZG->SetTicky();
  TPad *pad2_ZG = new TPad("pad2_ZG","pad2_ZG",0,0.,1,0.27);
  pad2_ZG->SetTickx();
  pad2_ZG->SetTicky();

  // Draw Pad1
  pad1_ZG->SetRightMargin(0.075);
  pad1_ZG->SetTopMargin(0.06);
  pad1_ZG->SetBottomMargin(0.0);
  pad1_ZG->Draw();
  pad1_ZG->cd();

  vector<float> bins = selectBinning(observable,category);

  TH1* frame_ZG  = pad1_ZG->DrawFrame(bins.front(),0.,bins.back(),0.2, "");
  frame_ZG->GetXaxis()->SetTitle(observableLatex.c_str());
  frame_ZG->GetYaxis()->SetTitle("Ratio Z/#gamma");
  frame_ZG->GetYaxis()->CenterTitle();
  frame_ZG->GetXaxis()->SetLabelSize(0.);
  frame_ZG->GetXaxis()->SetLabelOffset(1.10);
  frame_ZG->GetXaxis()->SetTitleSize(0.);
  frame_ZG->GetYaxis()->SetTitleSize(0.050);

  frame_ZG->Draw();
  CMS_lumi(pad1_ZG,"2.30",true);
 
  canvas_ZG->cd();
  pad2_ZG->SetTopMargin(0.04);
  pad2_ZG->SetBottomMargin(0.35);
  pad2_ZG->SetRightMargin(0.075);
  pad2_ZG->Draw();
  pad2_ZG->cd();

  TH1* frame2_ZG = NULL;
  if(category <= 1)
    frame2_ZG = pad2_ZG->DrawFrame(bins.front(), 0.5, bins.back(), 1.5, "");
  else
    frame2_ZG = pad2_ZG->DrawFrame(bins.front(), 0.0, bins.back(), 2.0, "");

  frame2_ZG->GetXaxis()->SetLabelSize(0.10);
  frame2_ZG->GetXaxis()->SetLabelOffset(0.03);
  frame2_ZG->GetXaxis()->SetTitleSize(0.13);
  frame2_ZG->GetXaxis()->SetTitleOffset(1.05);
  frame2_ZG->GetYaxis()->SetLabelSize(0.08);
  frame2_ZG->GetYaxis()->SetTitleSize(0.10);
  frame2_ZG->GetXaxis()->SetTitle(observableLatex.c_str());
  frame2_ZG->GetYaxis()->SetNdivisions(504, false);
  frame2_ZG->GetYaxis()->SetTitle("Data/Pred.");
  frame2_ZG->GetYaxis()->SetTitleOffset(0.5);
  frame2_ZG->Draw();
 
  // histo style
  ZGData_mm->SetLineColor(kBlack);
  ZGData_mm->SetLineWidth(2);
  ZGData_mm->SetMarkerColor(kBlack);
  ZGData_mm->SetMarkerStyle(20);
  ZGData_mm->SetMarkerSize(1.);

  ZGMC_mm->SetLineColor(kRed);
  ZGMC_mm->SetLineWidth(2);
  ZGMC_mm->SetMarkerSize(0);

  ZGData_ee->SetLineColor(kBlack);
  ZGData_ee->SetLineWidth(2);
  ZGData_ee->SetMarkerColor(kBlack);
  ZGData_ee->SetMarkerStyle(20);
  ZGData_ee->SetMarkerSize(1.);

  ZGMC_ee->SetLineColor(kRed);
  ZGMC_ee->SetLineWidth(2);
  ZGMC_ee->SetMarkerSize(0);

  ZGData_ll->SetLineColor(kBlack);
  ZGData_ll->SetLineWidth(2);
  ZGData_ll->SetMarkerColor(kBlack);
  ZGData_ll->SetMarkerStyle(20);
  ZGData_ll->SetMarkerSize(1.);

  ZGMC_ll->SetLineColor(kRed);
  ZGMC_ll->SetLineWidth(2);
  ZGMC_ll->SetMarkerSize(0);

  // Draw things
  pad1_ZG->cd();
  CMS_lumi(pad1_ZG, "2.30",true);
  TH1* ZGMC_mm_band = (TH1*) ZGMC_mm->Clone("ZGMC_mm_band");
  ZGMC_mm_band->SetFillColor(kGray);
  ZGMC_mm_band->Draw("E2same");
  ZGMC_mm->Draw("HIST same");
  ZGData_mm->Draw("PESAME");

  TLegend* leg = new TLegend(0.18, 0.66, 0.45, 0.92);
  leg->AddEntry(ZGMC_mm_band,"Z(#mu#mu)/#gamma MC","FL");
  leg->AddEntry(ZGData_mm,"Z(#mu#mu)/#gamma Data","PL");    
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("same");
  pad1_ZG->RedrawAxis("sameaxis");

  pad2_ZG->cd();
  TH1* ratioZG_mm       = (TH1*) ZGData_mm->Clone("ratioZG_mm");
  TH1* ratioZGd_mm      = (TH1*) ZGMC_mm->Clone("ratioZGd_mm");
  TH1* ratioZGD_mm      = (TH1*) ZGMC_mm->Clone("ratioZGD_mm");
  TH1* ratioZGD_mm_band = (TH1*) ZGMC_mm_band->Clone("ratioZGD_mm_band");

  for (int i = 1; i <= ratioZGd_mm->GetNbinsX(); i++) ratioZGd_mm->SetBinError(i, 0);
  
  ratioZG_mm->Divide(ratioZGd_mm);
  ratioZGD_mm->Divide(ratioZGd_mm);
  ratioZGD_mm_band->Divide(ratioZGd_mm);

  ratioZGD_mm_band->Draw("E2 same");
  ratioZGD_mm->Draw("HIST same");
  ratioZG_mm->Draw("PE same");
  pad2_ZG->RedrawAxis("sameaxis");

  canvas_ZG->SaveAs("ZG_mumu.png","png");
  canvas_ZG->SaveAs("ZG_mumu.pdf","pdf");

  ////////////
  pad1_ZG->cd(); 
  frame_ZG->Draw();
  TH1* ZGMC_ee_band = (TH1*) ZGMC_ee->Clone("ZGMC_ee_band");
  CMS_lumi(pad1_ZG, "2.30",true);
  ZGMC_ee_band->SetFillColor(kGray);
  ZGMC_ee_band->Draw("E2same");
  ZGMC_ee->Draw("HIST same");
  ZGData_ee->Draw("PESAME");
  
  leg->Clear();
  leg->AddEntry(ZGMC_ee_band,"Z(ee)/#gamma MC","FL");
  leg->AddEntry(ZGData_ee,"Z(ee)/#gamma Data","PL");    
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("same");
  pad1_ZG->RedrawAxis("sameaxis");

  pad2_ZG->cd();
  frame2_ZG->Draw();
  TH1* ratioZG_ee       = (TH1*) ZGData_ee->Clone("ratioZG_ee");
  TH1* ratioZGd_ee      = (TH1*) ZGMC_ee->Clone("ratioZGd_ee");
  TH1* ratioZGD_ee      = (TH1*) ZGMC_ee->Clone("ratioZGD_ee");
  TH1* ratioZGD_ee_band = (TH1*) ZGMC_ee_band->Clone("ratioZGD_ee_band");

  for (int i = 1; i <= ratioZGd_ee->GetNbinsX(); i++) ratioZGd_ee->SetBinError(i, 0);
  
  ratioZG_ee->Divide(ratioZGd_ee);
  ratioZGD_ee->Divide(ratioZGd_ee);
  ratioZGD_ee_band->Divide(ratioZGd_ee);

  ratioZGD_ee_band->Draw("E2 same");
  ratioZGD_ee->Draw("HIST same");
  ratioZG_ee->Draw("PE same");
  pad2_ZG->RedrawAxis("sameaxis");

  canvas_ZG->SaveAs("ZG_ee.png","png");
  canvas_ZG->SaveAs("ZG_ee.pdf","pdf");


  ////////////
  pad1_ZG->cd();
  frame_ZG->Draw();
  TH1* ZGMC_ll_band = (TH1*) ZGMC_ll->Clone("ZGMC_ll_band");
  ZGMC_ll_band->SetFillColor(kGray);
  CMS_lumi(pad1_ZG, "2.30",true);
  ZGMC_ll_band->Draw("E2same");
  ZGMC_ll->Draw("HIST same");
  ZGData_ll->Draw("PESAME");
  
  leg->Clear();
  leg->AddEntry(ZGMC_ll_band,"Z(ll)/#gamma MC","FL");
  leg->AddEntry(ZGData_ll,"Z(ll)/#gamma Data","PL");    
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("same");
  pad1_ZG->RedrawAxis("sameaxis");

  pad2_ZG->cd();
  frame2_ZG->Draw();
  TH1* ratioZG_ll       = (TH1*) ZGData_ll->Clone("ratioZG_ll");
  TH1* ratioZGd_ll      = (TH1*) ZGMC_ll->Clone("ratioZGd_ll");
  TH1* ratioZGD_ll      = (TH1*) ZGMC_ll->Clone("ratioZGD_ll");
  TH1* ratioZGD_ll_band = (TH1*) ZGMC_ll_band->Clone("ratioZGD_ll_band");

  for (int i = 1; i <= ratioZGd_ll->GetNbinsX(); i++) ratioZGd_ll->SetBinError(i, 0);
  
  ratioZG_ll->Divide(ratioZGd_ll);
  ratioZGD_ll->Divide(ratioZGd_ll);
  ratioZGD_ll_band->Divide(ratioZGd_ll);

  ratioZGD_ll_band->Draw("E2 same");
  ratioZGD_ll->Draw("HIST same");
  ratioZG_ll->Draw("PE same");
  pad2_ZG->RedrawAxis("sameaxis");

  canvas_ZG->SaveAs("ZG_ll.png","png");
  canvas_ZG->SaveAs("ZG_ll.pdf","pdf");


  // final plot
  TCanvas* canvas_ZW = new TCanvas("canvas_ZW","",600,700);
  canvas_ZW->SetTickx();
  canvas_ZW->SetTicky();
  canvas_ZW->cd();
  canvas_ZW->SetLeftMargin(0.11);
  TPad *pad1_ZW = new TPad("pad1_ZW","pad1_ZW",0,0.3,1,1);
  pad1_ZW->SetTickx();
  pad1_ZW->SetTicky();
  TPad *pad2_ZW = new TPad("pad2_ZW","pad2_ZW",0,0.,1,0.28);
  pad2_ZW->SetTickx();
  pad2_ZW->SetTicky();

  // Draw Pad1
  pad1_ZW->SetRightMargin(0.075);
  pad1_ZW->SetTopMargin(0.06);
  pad1_ZW->SetBottomMargin(0.0);
  pad1_ZW->Draw();
  pad1_ZW->cd();

  TH1* frame_ZW  = pad1_ZW->DrawFrame(bins.front(),0.,bins.back(),0.30, "");
  frame_ZW->GetXaxis()->SetTitle(observableLatex.c_str());
  frame_ZW->GetYaxis()->SetTitle("Ratio Z/W");
  frame_ZW->GetYaxis()->CenterTitle();
  frame_ZW->GetXaxis()->SetLabelSize(0.);
  frame_ZW->GetXaxis()->SetLabelOffset(1.10);
  frame_ZW->GetXaxis()->SetTitleSize(0.);
  frame_ZW->GetYaxis()->SetTitleSize(0.050);

  frame_ZW->Draw();
  CMS_lumi(pad1_ZW, "2.30",true);
 
  canvas_ZW->cd();
  pad2_ZW->SetTopMargin(0.04);
  pad2_ZW->SetBottomMargin(0.35);
  pad2_ZW->SetRightMargin(0.075);
  pad2_ZW->Draw();
  pad2_ZW->cd();

  TH1* frame2_ZW = NULL;
  if(category <= 1)
    frame2_ZW = pad2_ZW->DrawFrame(bins.front(), 0.5, bins.back(), 1.5, "");
  else
    frame2_ZW = pad2_ZW->DrawFrame(bins.front(), 0.0, bins.back(), 2.0, "");

  frame2_ZW->GetXaxis()->SetLabelSize(0.10);
  frame2_ZW->GetXaxis()->SetLabelOffset(0.03);
  frame2_ZW->GetXaxis()->SetTitleSize(0.13);
  frame2_ZW->GetXaxis()->SetTitleOffset(1.05);
  frame2_ZW->GetYaxis()->SetLabelSize(0.08);
  frame2_ZW->GetYaxis()->SetTitleSize(0.10);
  frame2_ZW->GetXaxis()->SetTitle(observableLatex.c_str());
  frame2_ZW->GetYaxis()->SetNdivisions(504, false);
  frame2_ZW->GetYaxis()->SetTitle("Data/Pred.");
  frame2_ZW->GetYaxis()->SetTitleOffset(0.5);
  frame2_ZW->Draw();
 
  // histo style
  ZWData_mm->SetLineColor(kBlack);
  ZWData_mm->SetLineWidth(2);
  ZWData_mm->SetMarkerColor(kBlack);
  ZWData_mm->SetMarkerStyle(20);
  ZWData_mm->SetMarkerSize(1.);

  ZWMC_mm->SetLineColor(kRed);
  ZWMC_mm->SetLineWidth(2);
  ZWMC_mm->SetMarkerSize(0);

  ZWData_ee->SetLineColor(kBlack);
  ZWData_ee->SetLineWidth(2);
  ZWData_ee->SetMarkerColor(kBlack);
  ZWData_ee->SetMarkerStyle(20);
  ZWData_ee->SetMarkerSize(1.);

  ZWMC_ee->SetLineColor(kRed);
  ZWMC_ee->SetLineWidth(2);
  ZWMC_ee->SetMarkerSize(0);

  ZWData_ll->SetLineColor(kBlack);
  ZWData_ll->SetLineWidth(2);
  ZWData_ll->SetMarkerColor(kBlack);
  ZWData_ll->SetMarkerStyle(20);
  ZWData_ll->SetMarkerSize(1.);

  ZWMC_ll->SetLineColor(kRed);
  ZWMC_ll->SetLineWidth(2);
  ZWMC_ll->SetMarkerSize(0);

  // Draw things
  pad1_ZW->cd();
  TH1* ZWMC_mm_band = (TH1*) ZWMC_mm->Clone("ZWMC_mm_band");
  CMS_lumi(pad1_ZW, "2.30",true);
  ZWMC_mm_band->SetFillColor(kGray);
  ZWMC_mm_band->Draw("E2same");
  ZWMC_mm->Draw("HIST same");
  ZWData_mm->Draw("PESAME");

  leg->Clear();
  leg->AddEntry(ZWMC_mm_band,"Z(#mu#mu)/W(#mu#nu) MC","FL");
  leg->AddEntry(ZWData_mm,"Z(#mu#mu)/W(#mu#nu) Data","PL");    
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("same");
  pad1_ZW->RedrawAxis("sameaxis");

  pad2_ZW->cd();
  TH1* ratioZW_mm       = (TH1*) ZWData_mm->Clone("ratioZW_mm");
  TH1* ratioZWd_mm      = (TH1*) ZWMC_mm->Clone("ratioZWd_mm");
  TH1* ratioZWD_mm      = (TH1*) ZWMC_mm->Clone("ratioZWD_mm");
  TH1* ratioZWD_mm_band = (TH1*) ZWMC_mm_band->Clone("ratioZWD_mm_band");

  for (int i = 1; i <= ratioZWd_mm->GetNbinsX(); i++) ratioZWd_mm->SetBinError(i, 0);
  
  ratioZW_mm->Divide(ratioZWd_mm);
  ratioZWD_mm->Divide(ratioZWd_mm);
  ratioZWD_mm_band->Divide(ratioZWd_mm);

  ratioZWD_mm_band->Draw("E2 same");
  ratioZWD_mm->Draw("HIST same");
  ratioZW_mm->Draw("PE same");
  pad2_ZW->RedrawAxis("sameaxis");

  canvas_ZW->SaveAs("ZW_mumu.png","png");
  canvas_ZW->SaveAs("ZW_mumu.pdf","pdf");

  ////////////
  pad1_ZW->cd(); 
  frame_ZW->Draw();
  TH1* ZWMC_ee_band = (TH1*) ZWMC_ee->Clone("ZWMC_ee_band");
  CMS_lumi(pad1_ZW, "2.30",true);
  ZWMC_ee_band->SetFillColor(kGray);
  ZWMC_ee_band->Draw("E2same");
  ZWMC_ee->Draw("HIST same");
  ZWData_ee->Draw("PESAME");
  
  leg->Clear();
  leg->AddEntry(ZWMC_ee_band,"Z(ee)/W(e#nu) MC","FL");
  leg->AddEntry(ZWData_ee,"Z(ee)/W(e#nu) Data","PL");    
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("same");
  pad1_ZW->RedrawAxis("sameaxis");

  pad2_ZW->cd();
  frame2_ZW->Draw();
  TH1* ratioZW_ee       = (TH1*) ZWData_ee->Clone("ratioZW_ee");
  TH1* ratioZWd_ee      = (TH1*) ZWMC_ee->Clone("ratioZWd_ee");
  TH1* ratioZWD_ee      = (TH1*) ZWMC_ee->Clone("ratioZWD_ee");
  TH1* ratioZWD_ee_band = (TH1*) ZWMC_ee_band->Clone("ratioZWD_ee_band");

  for (int i = 1; i <= ratioZWd_ee->GetNbinsX(); i++) ratioZWd_ee->SetBinError(i, 0);
  
  ratioZW_ee->Divide(ratioZWd_ee);
  ratioZWD_ee->Divide(ratioZWd_ee);
  ratioZWD_ee_band->Divide(ratioZWd_ee);

  ratioZWD_ee_band->Draw("E2 same");
  ratioZWD_ee->Draw("HIST same");
  ratioZW_ee->Draw("PE same");
  pad2_ZW->RedrawAxis("sameaxis");

  canvas_ZW->SaveAs("ZW_ee.png","png");
  canvas_ZW->SaveAs("ZW_ee.pdf","pdf");


  ////////////
  pad1_ZW->cd();
  frame_ZW->Draw();
  TH1* ZWMC_ll_band = (TH1*) ZWMC_ll->Clone("ZWMC_ll_band");
  CMS_lumi(pad1_ZW, "2.30",true);
  ZWMC_ll_band->SetFillColor(kGray);
  ZWMC_ll_band->Draw("E2same");
  ZWMC_ll->Draw("HIST same");
  ZWData_ll->Draw("PESAME");
  
  leg->Clear();
  leg->AddEntry(ZWMC_ll_band,"Z(ll)/W(l#nu) MC","FL");
  leg->AddEntry(ZWData_ll,"Z(ll)/W(l#nu) Data","PL");    
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("same");
  pad1_ZW->RedrawAxis("sameaxis");

  pad2_ZW->cd();
  frame2_ZW->Draw();
  TH1* ratioZW_ll       = (TH1*) ZWData_ll->Clone("ratioZW_ll");
  TH1* ratioZWd_ll      = (TH1*) ZWMC_ll->Clone("ratioZWd_ll");
  TH1* ratioZWD_ll      = (TH1*) ZWMC_ll->Clone("ratioZWD_ll");
  TH1* ratioZWD_ll_band = (TH1*) ZWMC_ll_band->Clone("ratioZWD_ll_band");

  for (int i = 1; i <= ratioZWd_ll->GetNbinsX(); i++) ratioZWd_ll->SetBinError(i, 0);
  
  ratioZW_ll->Divide(ratioZWd_ll);
  ratioZWD_ll->Divide(ratioZWd_ll);
  ratioZWD_ll_band->Divide(ratioZWd_ll);

  ratioZWD_ll_band->Draw("E2 same");
  ratioZWD_ll->Draw("HIST same");
  ratioZW_ll->Draw("PE same");
  pad2_ZW->RedrawAxis("sameaxis");

  canvas_ZW->SaveAs("ZW_ll.png","png");
  canvas_ZW->SaveAs("ZW_ll.pdf","pdf");


}
