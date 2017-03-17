#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void makePlot(TH1* histoData, TH1* histoMC,const string & observable, const Category & category, const string & observableLatex, string postfix){

  TCanvas* canvas = new TCanvas(("canvas"+postfix).c_str(), "canvas", 600, 700);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad* pad2 = new TPad(("pad2"+postfix).c_str(),"pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);

  // Draw Pad1                                                                                                                                                                                        
  canvas->cd();
  vector<double> bins = selectBinning(observable,category);
  TH1* frame  = canvas->DrawFrame(bins.front(),0.,bins.back(),0.3, "");

  frame->GetXaxis()->SetTitle(observableLatex.c_str());
  if(TString(postfix).Contains("ZG"))
    frame->GetYaxis()->SetTitle("Ratio Z/#gamma");
  else if(TString(postfix).Contains("ZW"))
    frame->GetYaxis()->SetTitle("Ratio Z/W");
  else if(TString(postfix).Contains("WG"))
    frame->GetYaxis()->SetTitle("Ratio W/#gamma");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitleOffset(1.25);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleSize(0.045);
  if(category == Category::monojet)
    frame->GetXaxis()->SetNdivisions(510);
  else
    frame->GetXaxis()->SetNdivisions(504);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
  frame->Draw();
  CMS_lumi(canvas,"35.9");

  float maxdata  = -1;
  float mindata  = 1000000;
  for(int iBin = 0; iBin < histoData->GetNbinsX(); iBin++){
    if(histoData->GetBinContent(iBin+1)+histoData->GetBinError(iBin+1) >= maxdata)
      maxdata = histoData->GetBinContent(iBin+1)+histoData->GetBinError(iBin+1);
    if(histoData->GetBinContent(iBin+1)-histoData->GetBinError(iBin+1) <= mindata)
      mindata = histoData->GetBinContent(iBin+1)-histoData->GetBinError(iBin+1);
  }
  float maxmc  = -1;
  float minmc  = 1000000;
  for(int iBin = 0; iBin < histoMC->GetNbinsX(); iBin++){
    if(histoMC->GetBinContent(iBin+1)+histoMC->GetBinError(iBin+1) >= maxmc)
      maxmc = histoMC->GetBinContent(iBin+1)+histoMC->GetBinError(iBin+1);
    if(histoMC->GetBinContent(iBin+1)-histoMC->GetBinError(iBin+1) <= maxmc)
      minmc = histoMC->GetBinContent(iBin+1)-histoMC->GetBinError(iBin+1);
  }

  if(TString(postfix).Contains("ZW"))
    frame->GetYaxis()->SetRangeUser(0.0,0.25);
  if(TString(postfix).Contains("ZG"))
    frame->GetYaxis()->SetRangeUser(0.0,0.20);
  else
    if(TString(postfix).Contains("WG"))
      frame->GetYaxis()->SetRangeUser(0.0,2.2);

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
  histoMCband->Draw("E2same");
  histoMC->Draw("HIST same");
  histoData->Draw("PESAME");

  TLegend* leg = new TLegend(0.55, 0.7, 0.9, 0.9);
  if(TString(postfix).Contains("ZG_mm")){
    leg->AddEntry(histoData,"Z(#mu#mu)/#gamma Data","PL");
    leg->AddEntry(histoMCband,"Z(#mu#mu)/#gamma MC","FL");
  }
  else if(TString(postfix).Contains("ZG_ee")){
    leg->AddEntry(histoData,"Z(ee)/#gamma Data","PL");
    leg->AddEntry(histoMCband,"Z(ee)/#gamma MC","FL");
  }
  else if(TString(postfix).Contains("ZW_mm")){
    leg->AddEntry(histoData,"Z(#mu#mu)/W(#mu#nu) Data","PL");
    leg->AddEntry(histoMCband,"Z(#mu#mu)/W(#mu#nu) MC","FL");
  }
  else if(TString(postfix).Contains("ZW_ee")){
    leg->AddEntry(histoData,"Z(ee)/W(e#nu) Data","PL");
    leg->AddEntry(histoMCband,"Z(ee)/W(e#nu) MC","FL");
  }
  else if(TString(postfix).Contains("WG_m")){
    leg->AddEntry(histoData,"W(#mu#nu)/#gamma Data","PL");
    leg->AddEntry(histoMCband,"W(#mu#nu)/#gamma MC","FL");
  }
  else if(TString(postfix).Contains("WG_e")){
    leg->AddEntry(histoData,"W(e#nu)/#gamma Data","PL");
    leg->AddEntry(histoMCband,"W(e#nu)/#gamma MC","FL");
  }

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
    frame2 = pad2->DrawFrame(bins.front(), 0.25, bins.back(), 1.75, "");
  else if(category == Category::monoV)
    frame2 = pad2->DrawFrame(bins.front(), 0.25, bins.back(), 1.75, "");
  else if(category == Category::VBF)
    frame2 = pad2->DrawFrame(bins.front(), 0.25, bins.back(), 1.75, "");

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


void makePostFitValidationPlot(string fitFilename, string observable, Category category, bool isCombinedFit = false,bool plotSBFit = false){
  
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  initializeBinning();
  
  TCanvas* canvas = NULL;
  TPad *pad2 = NULL;

  canvas = new TCanvas("canvas", "canvas", 600, 700);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);


  TFile* pfile = new TFile(fitFilename.c_str());

  TGraphAsymmErrors* dthist_gam = NULL;
  TGraphAsymmErrors* dthist_wmn = NULL;
  TGraphAsymmErrors* dthist_wen = NULL;
  TGraphAsymmErrors* dthist_zmm = NULL;
  TGraphAsymmErrors* dthist_zee = NULL;
  TH1* totalbkg_gam = NULL;
  TH1* totalbkg_wmn = NULL;
  TH1* totalbkg_wen = NULL;
  TH1* totalbkg_zmm = NULL;
  TH1* totalbkg_zee = NULL;

  string fit_dir = "shapes_fit_b";
  if(plotSBFit)
    fit_dir = "shapes_fit_s";

  string dir_gam;
  string dir_wen;
  string dir_wmn;
  string dir_zmm;
  string dir_zee;
  if(isCombinedFit){
    if(category == Category::monojet){
      dir_gam = "ch1_ch4";
      dir_wen = "ch1_ch6";
      dir_wmn = "ch1_ch3";
      dir_zmm = "ch1_ch2";
      dir_zee = "ch1_ch5";
    }
    else if(category == Category::monoV){
      dir_gam = "ch2_ch4";
      dir_wen = "ch2_ch6";
      dir_wmn = "ch2_ch3";
      dir_zmm = "ch2_ch2";
      dir_zee = "ch2_ch5";
    }
    else if(category == Category::VBF){
      dir_wen = "ch3_ch5";
      dir_wmn = "ch3_ch3";
      dir_zmm = "ch3_ch2";
      dir_zee = "ch3_ch4";
    }
  }
  else{    
    dir_gam = "ch4";
    dir_wen = "ch6";
    dir_wmn = "ch3";
    dir_zmm = "ch2";
    dir_zee = "ch5";
  }

  string postfix = "_MJ";
  if(category == Category::monoV)
    postfix = "_MV";
  else if(category == Category::VBF)
    postfix = "_VBF";

  dthist_gam = (TGraphAsymmErrors*)pfile->Get((fit_dir+"/"+dir_gam+"/data").c_str());
  dthist_wen = (TGraphAsymmErrors*)pfile->Get((fit_dir+"/"+dir_wen+"/data").c_str());
  dthist_wmn = (TGraphAsymmErrors*)pfile->Get((fit_dir+"/"+dir_wmn+"/data").c_str());
  dthist_zmm = (TGraphAsymmErrors*)pfile->Get((fit_dir+"/"+dir_zmm+"/data").c_str());
  dthist_zee = (TGraphAsymmErrors*)pfile->Get((fit_dir+"/"+dir_zee+"/data").c_str());

  totalbkg_gam = (TH1*)pfile->Get((fit_dir+"/"+dir_gam+"/total_background").c_str());
  totalbkg_wen = (TH1*)pfile->Get((fit_dir+"/"+dir_wen+"/total_background").c_str());
  totalbkg_wmn = (TH1*)pfile->Get((fit_dir+"/"+dir_wmn+"/total_background").c_str());
  totalbkg_zmm = (TH1*)pfile->Get((fit_dir+"/"+dir_zmm+"/total_background").c_str());
  totalbkg_zee = (TH1*)pfile->Get((fit_dir+"/"+dir_zee+"/total_background").c_str());

  // Zmm-gamma
  TH1* zmmgam = (TH1*) totalbkg_zmm->Clone("zmmgam");
  zmmgam->Divide(totalbkg_gam);
  // Zee-gamma
  TH1* zeegam = (TH1*) totalbkg_zee->Clone("zeegam");
  zeegam->Divide(totalbkg_gam);
  // Zmm-Wmn
  TH1* zmmwmn = (TH1*) totalbkg_zmm->Clone("zmmwmn");
  zmmwmn->Divide(totalbkg_wmn);
  // Zww-Wen
  TH1* zeewen = (TH1*) totalbkg_zee->Clone("zeewen");
  zeewen->Divide(totalbkg_wen);
  // Wmn-gam
  TH1* wmngam = (TH1*) totalbkg_wmn->Clone("wmngam");
  wmngam->Divide(totalbkg_gam);
  // Wen-gam
  TH1* wengam = (TH1*) totalbkg_wen->Clone("wengam");
  wengam->Divide(totalbkg_gam);

  TH1* zll = (TH1*) totalbkg_zmm->Clone("zll");
  zll->Add(totalbkg_zee);
  TH1* wln = (TH1*) totalbkg_wmn->Clone("wln");
  wln->Add(totalbkg_wen);

  TH1* zllgam = (TH1*) zll->Clone("zllgam");
  zllgam->Divide(totalbkg_gam);
  TH1* zllwln = (TH1*) zll->Clone("zllwln");
  zllwln->Divide(wln);
  TH1* wlngam = (TH1*) wln->Clone("wlngam");
  wlngam->Divide(totalbkg_gam);

  TH1* zmm_data = (TH1*) totalbkg_zmm->Clone("zmm_data");
  TH1* zee_data = (TH1*) totalbkg_zee->Clone("zee_data");
  TH1* wmn_data = (TH1*) totalbkg_wmn->Clone("wmn_data");
  TH1* wen_data = (TH1*) totalbkg_wen->Clone("wen_data");
  TH1* gam_data = (TH1*) totalbkg_gam->Clone("gam_data");

  for(int iBin = 0; iBin < zmm_data->GetNbinsX()+1; iBin++){
    double x,y;
    dthist_zmm->GetPoint(iBin,x,y);
    zmm_data->SetBinContent(iBin+1,y);
    zmm_data->SetBinError(iBin+1,(fabs(dthist_zmm->GetErrorYlow(iBin))+fabs(dthist_zmm->GetErrorYhigh(iBin)))/2);
    dthist_zee->GetPoint(iBin,x,y);
    zee_data->SetBinContent(iBin+1,y);
    zee_data->SetBinError(iBin+1,(fabs(dthist_zee->GetErrorYlow(iBin))+fabs(dthist_zee->GetErrorYhigh(iBin)))/2);
    dthist_wen->GetPoint(iBin,x,y);
    wen_data->SetBinContent(iBin+1,y);
    wen_data->SetBinError(iBin+1,(fabs(dthist_wen->GetErrorYlow(iBin))+fabs(dthist_wen->GetErrorYhigh(iBin)))/2);
    dthist_wmn->GetPoint(iBin,x,y);
    wmn_data->SetBinContent(iBin+1,y);
    wmn_data->SetBinError(iBin+1,(fabs(dthist_wmn->GetErrorYlow(iBin))+fabs(dthist_wmn->GetErrorYhigh(iBin)))/2);
    dthist_gam->GetPoint(iBin,x,y);
    gam_data->SetBinContent(iBin+1,y);
    gam_data->SetBinError(iBin+1,(fabs(dthist_gam->GetErrorYlow(iBin))+fabs(dthist_gam->GetErrorYhigh(iBin)))/2);
  }

  TH1* zll_data = (TH1*) zmm_data->Clone("zll_data");
  zll_data->Add(zee_data);
  TH1* wln_data = (TH1*) wmn_data->Clone("wln_data");
  wln_data->Add(wen_data);

  // Zmm-gamma
  TH1* zmmgam_data = (TH1*) zmm_data->Clone("zmmgam_data");
  zmmgam_data->Divide(gam_data);
  // zee-gamma
  TH1* zeegam_data = (TH1*) zee_data->Clone("zeegam_data");
  zeegam_data->Divide(gam_data);
  // Zmm-wmnma
  TH1* zmmwmn_data = (TH1*) zmm_data->Clone("zmmwmn_data");
  zmmwmn_data->Divide(wmn_data);
  // zee-wenma
  TH1* zeewen_data = (TH1*) zee_data->Clone("zeewen_data");
  zeewen_data->Divide(wen_data);
  // Wmn-gam
  TH1* wmngam_data = (TH1*) wmn_data->Clone("wmngam_data");
  wmngam_data->Divide(gam_data);
  // Wen-gam
  TH1* wengam_data = (TH1*) wen_data->Clone("wengam_data");
  wengam_data->Divide(gam_data);

  // Wln-gam
  TH1* wlngam_data = (TH1*) wln_data->Clone("wlngam_data");
  wlngam_data->Divide(gam_data);
  // Zll-gam
  TH1* zllgam_data = (TH1*) zll_data->Clone("zllgam_data");
  zllgam_data->Divide(gam_data);
  // Zll-Wln
  TH1* zllwln_data = (TH1*) zll_data->Clone("zllwln_data");
  zllwln_data->Divide(wln_data);
  

  makePlot(zmmgam_data,zmmgam,observable,category,"Recoil [GeV]","ZG_mm");
  makePlot(zeegam_data,zeegam,observable,category,"Recoil [GeV]","ZG_ee");
  makePlot(zmmwmn_data,zmmwmn,observable,category,"Recoil [GeV]","ZW_mm");
  makePlot(zeewen_data,zeewen,observable,category,"Recoil [GeV]","ZW_ee");
  makePlot(wmngam_data,wmngam,observable,category,"Recoil [GeV]","WG_m");
  makePlot(wengam_data,wengam,observable,category,"Recoil [GeV]","WG_e");

  makePlot(zllgam_data,zllgam,observable,category,"Recoil [GeV]","ZG_ll");
  makePlot(wlngam_data,wlngam,observable,category,"Recoil [GeV]","WG_l");
  makePlot(zllwln_data,zllwln,observable,category,"Recoil [GeV]","ZW_ll");

}

  
