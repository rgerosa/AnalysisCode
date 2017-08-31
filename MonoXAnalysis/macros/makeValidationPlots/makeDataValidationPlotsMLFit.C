#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void makePlot(TH1* histoData, TH1* histoMC, const string & observable, const Category & category, const string & observableLatex, const string & postfix){
  
  TCanvas* canvas = new TCanvas(("canvas_"+postfix).c_str(), "", 600, 700);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad* pad2 = new TPad(("pad2_"+postfix).c_str(),"",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);
  
  canvas->cd();
  // take the observable binning
  vector<double> bins = selectBinning(observable,category);

  // make the frame
  TH1* frame  = canvas->DrawFrame(bins.front(),0.,bins.back(),0.3, "");
  frame->GetXaxis()->SetTitle(observableLatex.c_str());

  // y-axis title
  if(TString(postfix).Contains("WG") and TString(postfix).Contains("mm"))
    frame->GetYaxis()->SetTitle("W #rightarrow #mu#mu / #gamma");
  if(TString(postfix).Contains("WG") and TString(postfix).Contains("ee"))
    frame->GetYaxis()->SetTitle("W #rightarrow ee / #gamma");
  if(TString(postfix).Contains("WG") and TString(postfix).Contains("ll"))
    frame->GetYaxis()->SetTitle("W #rightarrow ll / #gamma");

  if(TString(postfix).Contains("WG") and TString(postfix).Contains("m"))
    frame->GetYaxis()->SetTitle("W #rightarrow #mu#nu / #gamma");
  if(TString(postfix).Contains("WG") and TString(postfix).Contains("e"))
    frame->GetYaxis()->SetTitle("W #rightarrow e#nu / #gamma");
  if(TString(postfix).Contains("WG") and TString(postfix).Contains("l"))
    frame->GetYaxis()->SetTitle("W #rightarrow l#nu / #gamma");

  if(TString(postfix).Contains("ZW") and TString(postfix).Contains("mm"))
    frame->GetYaxis()->SetTitle("Z #rightarrow #mu#mu / W #rightarrow #mu#nu");
  if(TString(postfix).Contains("ZW") and TString(postfix).Contains("ee"))
    frame->GetYaxis()->SetTitle("Z #rightarrow ee / W #rightarrow e#nu");
  if(TString(postfix).Contains("ZW") and TString(postfix).Contains("ll"))
    frame->GetYaxis()->SetTitle("Z #rightarrow ll / W #rightarrow l#nu");

  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitleOffset(1.25);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleSize(0.045);

  if(category == Category::monojet)
    frame->GetXaxis()->SetNdivisions(510);
  else
    frame->GetXaxis()->SetNdivisions(506);

  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
  frame->Draw();

  CMS_lumi(canvas,"35.9");

  // Set the axis ranges -->
  if(category == Category::monojet){
    if(TString(postfix).Contains("ZG") and not TString(postfix).Contains("ll"))
      frame->GetYaxis()->SetRangeUser(0.0,0.12);
    if(TString(postfix).Contains("ZG") and TString(postfix).Contains("ll"))
      frame->GetYaxis()->SetRangeUser(0.02,0.18);
    if(TString(postfix).Contains("ZW") and not TString(postfix).Contains("ll"))
      frame->GetYaxis()->SetRangeUser(0.0,0.22);
    if(TString(postfix).Contains("ZW") and TString(postfix).Contains("ll"))
      frame->GetYaxis()->SetRangeUser(0.0,0.20);
    if(TString(postfix).Contains("WG") and TString(postfix).Contains("_m"))
      frame->GetYaxis()->SetRangeUser(0.3,1.4);
    if(TString(postfix).Contains("WG") and TString(postfix).Contains("_e"))
      frame->GetYaxis()->SetRangeUser(0.2,0.7);
    if(TString(postfix).Contains("WG") and TString(postfix).Contains("_l"))
      frame->GetYaxis()->SetRangeUser(0.4,1.7);
  }
  else if(category == Category::monoV){
    if(TString(postfix).Contains("ZG") and not TString(postfix).Contains("ll"))
      frame->GetYaxis()->SetRangeUser(0.0,0.14);
    if(TString(postfix).Contains("ZG") and TString(postfix).Contains("ll"))
      frame->GetYaxis()->SetRangeUser(0.02,0.22);
    if(TString(postfix).Contains("ZW") and not TString(postfix).Contains("ll"))
      frame->GetYaxis()->SetRangeUser(0.0,0.22);
    if(TString(postfix).Contains("ZW") and TString(postfix).Contains("ll"))
      frame->GetYaxis()->SetRangeUser(0.0,0.20);
    if(TString(postfix).Contains("WG") and TString(postfix).Contains("_m"))
      frame->GetYaxis()->SetRangeUser(0.2,1.4);
    if(TString(postfix).Contains("WG") and TString(postfix).Contains("_e"))
      frame->GetYaxis()->SetRangeUser(0.1,1.1);
    if(TString(postfix).Contains("WG") and TString(postfix).Contains("_l"))
      frame->GetYaxis()->SetRangeUser(0.6,2.2);
  }

  // histo style
  histoData->SetLineColor(kBlack);
  histoData->SetLineWidth(2);
  histoData->SetMarkerColor(kBlack);
  histoData->SetMarkerStyle(20);
  histoData->SetMarkerSize(1.);

  histoMC->SetLineColor(kRed);
  histoMC->SetLineWidth(2);
  histoMC->SetMarkerSize(0);

  TH1* histoMCband = (TH1*) histoMC->Clone("histoMCband");
  histoMCband->SetFillColor(kGray);
  histoMCband->SetLineColor(kGray);
  histoMCband->Draw("E2same");
  histoMC->Draw("HIST same");
  histoData->Draw("PESAME");

  TLegend* leg = new TLegend(0.55, 0.7, 0.9, 0.9);
  leg->AddEntry(histoData,"Data","PEL");
  leg->AddEntry(histoMC,"MC Prediction","L");
  leg->AddEntry(histoMCband,"Stat + Sys Unc.","F");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("same");
  canvas->RedrawAxis("sameaxis");

  canvas->cd();
  pad2->Draw();
  pad2->cd();
  
  TH1* frame2 = NULL;
  if(category == Category::monojet)
    frame2 = pad2->DrawFrame(bins.front(), 0.5, bins.back(), 1.5, "");
  else if(category == Category::monoV)
    frame2 = pad2->DrawFrame(bins.front(), 0.5, bins.back(), 1.5, "");
  else if(category == Category::VBF or category == Category::VBFrelaxed)
    frame2 = pad2->DrawFrame(bins.front(), 0.5, bins.back(), 1.5, "");

  frame2->GetYaxis()->SetNdivisions(5);
  frame2->GetXaxis()->SetTitle(observableLatex.c_str());
  frame2->GetYaxis()->SetTitle("Data / Pred.");
  frame2->GetYaxis()->CenterTitle();
  frame2->GetYaxis()->SetTitleOffset(1.5);
  frame2->GetYaxis()->SetLabelSize(0.04);
  frame2->GetYaxis()->SetTitleSize(0.04);
  frame2->GetXaxis()->SetLabelSize(0.04);
  frame2->GetXaxis()->SetTitleSize(0.05);
  frame2->GetXaxis()->SetTitleOffset(1.1);
  if(category == Category::monojet)
    frame2->GetXaxis()->SetNdivisions(510);
  else
    frame2->GetXaxis()->SetNdivisions(506);

 
  TH1* ratio         = (TH1*) histoData->Clone(("ratio_"+postfix).c_str());
  TH1* ratiod        = (TH1*) histoMC->Clone(("ratiod_"+postfix).c_str());
  TH1* ratioD        = (TH1*) histoMC->Clone(("ratioD_"+postfix).c_str());
  TH1* ratioD_band   = (TH1*) histoMCband->Clone(("ratioD_band_"+postfix).c_str());

  for (int i = 1; i <= ratiod->GetNbinsX(); i++) ratiod->SetBinError(i, 0);
  
  ratio->Divide(ratiod);
  ratioD->Divide(ratiod);
  ratioD_band->Divide(ratiod);

  ratioD_band->Draw("E2 same");
  ratioD->Draw("HIST same");
  ratio->Draw("PE same");
  pad2->RedrawAxis("sameaxis");

  canvas->SaveAs((postfix+".png").c_str(),"png");
  canvas->SaveAs((postfix+".pdf").c_str(),"pdf");
}


// make the data-validation plot from the mlfit.root given by combine --> using the prefit shapes
void makeDataValidationPlotsMLFit(string inputFileName, Category category, string observable, string observableLatex, int rebinFactor = 1, bool isCombination = false){

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  setTDRStyle();
  initializeBinning();

  // open the input file with all the templates
  TFile* inputFile = TFile::Open(inputFileName.c_str());
  
  string dir_zmm = "ch2";
  string dir_wmn = "ch3";
  string dir_gam = "ch4";
  string dir_zee = "ch5";
  string dir_wen = "ch6";

  bool useGammaJets = true;

  if(isCombination and category == Category::monojet){
    dir_zmm = "ch1_ch2";
    dir_wmn = "ch1_ch3";
    dir_gam = "ch1_ch4";
    dir_zee = "ch1_ch5";
    dir_wen = "ch1_ch6";
  }
  else if(isCombination and category == Category::monoV){
    dir_zmm = "ch2_ch2";
    dir_wmn = "ch2_ch3";
    dir_gam = "ch2_ch4";
    dir_zee = "ch2_ch5";
    dir_wen = "ch2_ch6";
  }
  
  if(category == Category::VBF or category == Category::VBFrelaxed){
    useGammaJets = false;
    dir_zmm = "ch2";
    dir_wmn = "ch3";
    dir_zee = "ch4";
    dir_wen = "ch5";
  }

  TGraphAsymmErrors* data_zmm = (TGraphAsymmErrors*) inputFile->Get(("shapes_prefit/"+dir_zmm+"/data").c_str());
  TGraphAsymmErrors* data_zee = (TGraphAsymmErrors*) inputFile->Get(("shapes_prefit/"+dir_zee+"/data").c_str());
  TGraphAsymmErrors* data_wmn = (TGraphAsymmErrors*) inputFile->Get(("shapes_prefit/"+dir_wmn+"/data").c_str());
  TGraphAsymmErrors* data_wen = (TGraphAsymmErrors*) inputFile->Get(("shapes_prefit/"+dir_wen+"/data").c_str());
  TGraphAsymmErrors* data_gam = NULL;
  if(useGammaJets)
    data_gam = (TGraphAsymmErrors*) inputFile->Get(("shapes_prefit/"+dir_gam+"/data").c_str());

  TH1* total_background_zmm = (TH1*) inputFile->Get(("shapes_prefit/"+dir_zmm+"/total_background").c_str());
  TH1* total_background_zee = (TH1*) inputFile->Get(("shapes_prefit/"+dir_zee+"/total_background").c_str());
  TH1* total_background_wmn = (TH1*) inputFile->Get(("shapes_prefit/"+dir_wmn+"/total_background").c_str());
  TH1* total_background_wen = (TH1*) inputFile->Get(("shapes_prefit/"+dir_wen+"/total_background").c_str());
  TH1* total_background_gam = NULL;
  if(useGammaJets)
    total_background_gam = (TH1*) inputFile->Get(("shapes_prefit/"+dir_gam+"/total_background").c_str());

  // convert TGraphAsymmErrors in histograms
  TH1* data_zmm_hist = (TH1*) total_background_zmm->Clone("data_zmm_hist");
  data_zmm_hist->Reset();
  TH1* data_zee_hist = (TH1*) total_background_zee->Clone("data_zee_hist");
  data_zee_hist->Reset();
  TH1* data_wmn_hist = (TH1*) total_background_wmn->Clone("data_wmn_hist");
  data_wmn_hist->Reset();
  TH1* data_wen_hist = (TH1*) total_background_wen->Clone("data_wen_hist");
  data_wen_hist->Reset();
  TH1* data_gam_hist = NULL;
  if(useGammaJets){
    data_gam_hist = (TH1*) total_background_gam->Clone("data_gam_hist");
    data_gam_hist->Reset();
  }

  for(int iPoint = 0; iPoint < data_zmm->GetN(); iPoint++){
    double x,y;
    data_zmm->GetPoint(iPoint,x,y);
    data_zmm_hist->SetBinContent(iPoint+1,y);
    data_zmm_hist->SetBinError(iPoint+1,(data_zmm->GetErrorYlow(iPoint)+data_zmm->GetErrorYhigh(iPoint))/2);

    data_zee->GetPoint(iPoint,x,y);
    data_zee_hist->SetBinContent(iPoint+1,y);
    data_zee_hist->SetBinError(iPoint+1,(data_zee->GetErrorYlow(iPoint)+data_zee->GetErrorYhigh(iPoint))/2);

    data_wmn->GetPoint(iPoint,x,y);
    data_wmn_hist->SetBinContent(iPoint+1,y);
    data_wmn_hist->SetBinError(iPoint+1,(data_wmn->GetErrorYlow(iPoint)+data_wmn->GetErrorYhigh(iPoint))/2);

    data_wen->GetPoint(iPoint,x,y);
    data_wen_hist->SetBinContent(iPoint+1,y);
    data_wen_hist->SetBinError(iPoint+1,(data_wen->GetErrorYlow(iPoint)+data_wen->GetErrorYhigh(iPoint))/2);

    if(data_gam_hist){
      data_gam->GetPoint(iPoint,x,y);
      data_gam_hist->SetBinContent(iPoint+1,y);
      data_gam_hist->SetBinError(iPoint+1,(data_gam->GetErrorYlow(iPoint)+data_gam->GetErrorYhigh(iPoint))/2);
    }
  }

  // Rebin if needed
  data_zmm_hist->Rebin(rebinFactor);
  data_zee_hist->Rebin(rebinFactor);
  data_wmn_hist->Rebin(rebinFactor);
  data_wen_hist->Rebin(rebinFactor);
  if(data_gam_hist)
    data_gam_hist->Rebin(rebinFactor);

  total_background_zmm->Rebin(rebinFactor);
  total_background_zee->Rebin(rebinFactor);
  total_background_wmn->Rebin(rebinFactor);
  total_background_wen->Rebin(rebinFactor);
  if(total_background_gam)
    total_background_gam->Rebin(rebinFactor);

  // Make the ratios Z/W
  TH1* ZWData_mm = (TH1*) data_zmm_hist->Clone("ZWData_mm");
  ZWData_mm->Divide(data_wmn_hist);
  TH1* ZWData_ee = (TH1*) data_zee_hist->Clone("ZWData_ee");
  ZWData_ee->Divide(data_wen_hist);
  TH1* ZWData_ll = (TH1*) data_zmm_hist->Clone("ZWData_ll");
  ZWData_ll->Add(data_zee_hist);
  TH1* ZWData_ll_den = (TH1*) data_wmn_hist->Clone("ZWData_ll_den");
  ZWData_ll_den->Add(data_wen_hist);
  ZWData_ll->Divide(ZWData_ll_den);

  TH1* ZWMC_mm = (TH1*) total_background_zmm->Clone("ZWMC_mm");
  ZWMC_mm->Divide(total_background_wmn);
  TH1* ZWMC_ee = (TH1*) total_background_zee->Clone("ZWMC_ee");
  ZWMC_ee->Divide(total_background_wen);
  TH1* ZWMC_ll = (TH1*) total_background_zmm->Clone("ZWMC_ll");
  ZWMC_ll->Add(total_background_zee);
  TH1* ZWMC_ll_den = (TH1*) total_background_wmn->Clone("ZWMC_ll_den");
  ZWMC_ll_den->Add(total_background_wen);
  ZWMC_ll->Divide(ZWMC_ll_den);

  // make the plots 
  makePlot(ZWData_mm,ZWMC_mm,observable,category,observableLatex,"ZW_mm");
  makePlot(ZWData_ee,ZWMC_ee,observable,category,observableLatex,"ZW_ee");
  makePlot(ZWData_ll,ZWMC_ll,observable,category,observableLatex,"ZW_ll");

  if(useGammaJets){

    TH1* ZGData_mm = (TH1*) data_zmm_hist->Clone("ZGData_mm");
    ZGData_mm->Divide(data_gam_hist);
    TH1* ZGData_ee = (TH1*) data_zee_hist->Clone("ZGData_ee");
    ZGData_ee->Divide(data_gam_hist);
    TH1* ZGData_ll = (TH1*) data_zmm_hist->Clone("ZGData_ll");
    ZGData_ll->Add(data_zee_hist);
    ZGData_ll->Divide(data_gam_hist);

    TH1* ZGMC_mm = (TH1*) total_background_zmm->Clone("ZGMC_mm");
    ZGMC_mm->Divide(total_background_gam);
    TH1* ZGMC_ee = (TH1*) total_background_zee->Clone("ZGMC_ee");
    ZGMC_ee->Divide(total_background_gam);
    TH1* ZGMC_ll = (TH1*) total_background_zmm->Clone("ZGMC_ll");
    ZGMC_ll->Add(total_background_zee);
    ZGMC_ll->Divide(total_background_gam);

    makePlot(ZGData_mm,ZGMC_mm,observable,category,observableLatex,"ZG_mm");
    makePlot(ZGData_ee,ZGMC_ee,observable,category,observableLatex,"ZG_ee");
    makePlot(ZGData_ll,ZGMC_ll,observable,category,observableLatex,"ZG_ll");
    

    TH1* WGData_m = (TH1*) data_wmn_hist->Clone("WGData_m");
    WGData_m->Divide(data_gam_hist);
    TH1* WGData_e = (TH1*) data_wen_hist->Clone("WGData_e");
    WGData_e->Divide(data_gam_hist);
    TH1* WGData_l = (TH1*) data_wmn_hist->Clone("WGData_l");
    WGData_l->Add(data_wen_hist);
    WGData_l->Divide(data_gam_hist);

    TH1* WGMC_m = (TH1*) total_background_wmn->Clone("WGMC_m");
    WGMC_m->Divide(total_background_gam);
    TH1* WGMC_e = (TH1*) total_background_wen->Clone("WGMC_e");
    WGMC_e->Divide(total_background_gam);
    TH1* WGMC_l = (TH1*) total_background_wmn->Clone("WGMC_l");
    WGMC_l->Add(total_background_wen);
    WGMC_l->Divide(total_background_gam);

    makePlot(WGData_m,WGMC_m,observable,category,observableLatex,"WG_m");
    makePlot(WGData_e,WGMC_e,observable,category,observableLatex,"WG_e");
    makePlot(WGData_l,WGMC_l,observable,category,observableLatex,"WG_l");
    
  }
}
