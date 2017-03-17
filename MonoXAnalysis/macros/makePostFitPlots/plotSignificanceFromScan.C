#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static float lumiScale_DM = 2.78;

//// significance wrt the CR-only fit
void plotSignificanceFromScan(string   fitFilename, 
			      Category category, 
			      bool     isCombinedFit = false,
			      string   postfix = ""
			      ){
  
  
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  if(postfix != "Vector" and postfix != "Scalar"){
    cerr<<"Problem --> only available options are Vector and Scalar "<<endl;
    return;
  }

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();

  TFile* pfile = new TFile(fitFilename.c_str());
  string fit_dir = "shapes_fit_s";

  //////////////////
  TFile*monoj_av = NULL, *monow_av = NULL, *monoz_av = NULL;  
  if(category == Category::monoV and postfix == "Vector"){
    monoj_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_800_0.25_catmonov_13TeV_v1.root","READ");
    monow_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_800_0.25_catmonov_13TeV_v1.root","READ");
    monoz_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_800_0.25_catmonov_13TeV_v1.root","READ");
  }
  else if(category == Category::monojet and postfix == "Vector"){
    monoj_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_800_0.25_catmonojet_13TeV_v1.root","READ");
    monow_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_800_0.25_catmonojet_13TeV_v1.root","READ");
    monoz_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_800_0.25_catmonojet_13TeV_v1.root","READ");
  }
  else if(category == Category::monojet and postfix == "Scalar"){
    monoj_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_805_1.0_catmonojet_13TeV_v1.root","READ");
    monow_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_805_1.0_catmonojet_13TeV_v1.root","READ");
    monoz_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_805_1.0_catmonojet_13TeV_v1.root","READ");
  }
  else if(category == Category::monoV and postfix == "Scalar"){
    monoj_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_805_1.0_catmonov_13TeV_v1.root","READ");
    monow_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_805_1.0_catmonov_13TeV_v1.root","READ");
    monoz_av = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_805_1.0_catmonov_13TeV_v1.root","READ");
  }

  // in case of b-only fit just dispaly three possible signal on the stack
  TH1* mjhist_av_1 = NULL;
  TH1* mwhist_av_1 = NULL;
  TH1* mzhist_av_1 = NULL;
  TH1* mjhist_av_2 = NULL;
  TH1* mwhist_av_2 = NULL;
  TH1* mzhist_av_2 = NULL;
  TH1* mjhist_av_3 = NULL;
  TH1* mwhist_av_3 = NULL;
  TH1* mzhist_av_3 = NULL;

  // signals for axial vector model 
  if(postfix == "Vector"){
    mjhist_av_1 = (TH1*) monoj_av->FindObjectAny("signal_signal_80015250001");
    mwhist_av_1 = (TH1*) monow_av->FindObjectAny("signal_signal_80015250001");
    mzhist_av_1 = (TH1*) monoz_av->FindObjectAny("signal_signal_80015250001");
    mjhist_av_2 = (TH1*) monoj_av->FindObjectAny("signal_signal_80020000001");
    mwhist_av_2 = (TH1*) monow_av->FindObjectAny("signal_signal_80020000001");
    mzhist_av_2 = (TH1*) monoz_av->FindObjectAny("signal_signal_80020000001");
    mjhist_av_3 = (TH1*) monoj_av->FindObjectAny("signal_signal_80025000001");
    mwhist_av_3 = (TH1*) monow_av->FindObjectAny("signal_signal_80025000001");
    mzhist_av_3 = (TH1*) monoz_av->FindObjectAny("signal_signal_80025000001");
  }
  else if(postfix == "Scalar"){
    mjhist_av_1 = (TH1*) monoj_av->FindObjectAny("signal_signal_80501000001");
    mwhist_av_1 = (TH1*) monow_av->FindObjectAny("signal_signal_80501000001");
    mzhist_av_1 = (TH1*) monoz_av->FindObjectAny("signal_signal_80501000001");
    mjhist_av_2 = (TH1*) monoj_av->FindObjectAny("signal_signal_80504000001");
    mwhist_av_2 = (TH1*) monow_av->FindObjectAny("signal_signal_80504000001");
    mzhist_av_2 = (TH1*) monoz_av->FindObjectAny("signal_signal_80504000001");
    mjhist_av_3 = (TH1*) monoj_av->FindObjectAny("signal_signal_80506000001");
    mwhist_av_3 = (TH1*) monow_av->FindObjectAny("signal_signal_80506000001");
    mzhist_av_3 = (TH1*) monoz_av->FindObjectAny("signal_signal_80506000001");
  }

  mjhist_av_1->Scale(1.,"width");
  if(mwhist_av_1 != NULL) mwhist_av_1->Scale(1.,"width");
  if(mzhist_av_1 != NULL) mzhist_av_1->Scale(1.,"width");
  mjhist_av_2->Scale(1.,"width");
  if(mwhist_av_2 != NULL) mwhist_av_2->Scale(1.,"width");
  if(mzhist_av_2 != NULL) mzhist_av_2->Scale(1.,"width");
  mjhist_av_3->Scale(1.,"width");
  if(mwhist_av_3 != NULL) mwhist_av_3->Scale(1.,"width");
  if(mzhist_av_3 != NULL) mzhist_av_3->Scale(1.,"width");
  
  // summing all signals together
  mjhist_av_1->Add(mwhist_av_1);
  mjhist_av_1->Add(mzhist_av_1);
  mjhist_av_1->Scale(lumiScale_DM);
  mjhist_av_2->Add(mwhist_av_2);
  mjhist_av_2->Add(mzhist_av_2);
  mjhist_av_2->Scale(lumiScale_DM);
  mjhist_av_3->Add(mwhist_av_3);
  mjhist_av_3->Add(mzhist_av_3);
  mjhist_av_3->Scale(lumiScale_DM);

  /////////////////
  string dir;
  if(isCombinedFit){
    if(category == Category::monojet)
      dir = "ch1_ch1";
    else if(category == Category::monoV)
      dir = "ch2_ch1";
    else if(category == Category::VBF)
      dir = "ch3_ch1";
  }
  else if( category != Category::VBF)
    dir = "ch1";
  else
    dir = "ch1";

  /////////////////
  if(category == Category::monojet)
    postfix += "_MJ";
  else if(category == Category::monoV)
    postfix += "_MV";
  else if(category == Category::VBF)
    postfix += "_VBF";

  // background
  TH1* tohist = (TH1*) ((TH1*)pfile->Get((fit_dir+"/"+dir+"/total_background").c_str()))->Clone("tohist");    

  TH1* frame = (TH1*) tohist->Clone("frame");
  frame->Reset();
  frame->SetLineColor(kBlack);
  frame->SetLineWidth(1);
  frame->GetYaxis()->SetTitle("(S+B)/B");
  frame->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.050);
  if(category == Category::monojet)
    frame->GetXaxis()->SetNdivisions(510);
  else
    frame->GetXaxis()->SetNdivisions(504);
  frame->Draw();
  CMS_lumi(canvas,"35.9");

  TLatex* categoryLabel = new TLatex();
  categoryLabel->SetNDC();
  categoryLabel->SetTextSize(0.5*canvas->GetTopMargin());
  categoryLabel->SetTextFont(42);
  categoryLabel->SetTextAlign(11);
  if(category == Category::monojet)
    categoryLabel ->DrawLatex(0.175,0.80,"monojet");
  else if(category == Category::monoV)
    categoryLabel ->DrawLatex(0.175,0.80,"mono-V");
  else if(category == Category::VBF)
    categoryLabel ->DrawLatex(0.175,0.80,"VBF");
  categoryLabel->Draw("same");

  // make signal uncertainty zero
  for(int iBin = 0; iBin < mjhist_av_1->GetNbinsX()+1; iBin++){
    mjhist_av_1->SetBinError(iBin+1,0);
    mjhist_av_2->SetBinError(iBin+1,0);
    mjhist_av_3->SetBinError(iBin+1,0);
  }


  TLegend* leg = new TLegend(0.35, 0.45, 0.62, 0.72);
  // leg->SetFillColor(0);
  //  leg->SetFillStyle(3001);
  //leg->SetBorderSize(0);


  mjhist_av_1->SetLineColor(kBlack);
  mjhist_av_1->SetLineWidth(2);
  mjhist_av_2->SetLineColor(kRed);
  mjhist_av_2->SetLineWidth(2);
  mjhist_av_3->SetLineColor(kBlue);
  mjhist_av_3->SetLineWidth(2);
  TH1* mjhist_av_1_norm = (TH1*) mjhist_av_1->Clone("mjhist_av_1_norm");
  mjhist_av_1_norm->Scale(1./mjhist_av_1_norm->Integral());
  TH1* mjhist_av_2_norm = (TH1*) mjhist_av_2->Clone("mjhist_av_2_norm");
  mjhist_av_2_norm->Scale(1./mjhist_av_2_norm->Integral());
  TH1* mjhist_av_3_norm = (TH1*) mjhist_av_3->Clone("mjhist_av_3_norm");
  mjhist_av_3_norm->Scale(1./mjhist_av_3_norm->Integral());
  TH1* bkgtot_norm = (TH1*) tohist->Clone("bkgtot_norm");
  bkgtot_norm->Scale(1./bkgtot_norm->Integral());

  mjhist_av_1->Add(tohist);
  mjhist_av_2->Add(tohist);
  mjhist_av_3->Add(tohist);
  mjhist_av_1->Divide(tohist);
  mjhist_av_2->Divide(tohist);
  mjhist_av_3->Divide(tohist);

  mjhist_av_1_norm->Divide(bkgtot_norm);
  mjhist_av_2_norm->Divide(bkgtot_norm);
  mjhist_av_3_norm->Divide(bkgtot_norm);


  frame->GetYaxis()->SetRangeUser(1.,max(mjhist_av_1->GetMaximum(),max(mjhist_av_2->GetMaximum(),mjhist_av_3->GetMaximum()))*1.05);

  TH1* postfit_band = (TH1*) mjhist_av_3->Clone("postfit_band");
  for(int iBin = 0; iBin < postfit_band->GetNbinsX(); iBin++)
    postfit_band->SetBinContent(iBin+1,1);
  postfit_band->SetLineColor(0);
    
  postfit_band->SetMarkerColor(0);
  postfit_band->SetMarkerSize(0);
  postfit_band->SetFillColor(kGray);
  postfit_band->SetFillStyle(1001);
  postfit_band->Draw("E2 SAME");

  TH1* unhist = (TH1*) mjhist_av_3->Clone("unhist");
  unhist->Reset();
  for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinContent(i, 1);
  for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinError(i, 0);
  unhist->SetMarkerSize(0);
  unhist->SetLineColor(kBlack);
  unhist->SetLineStyle(2);
  unhist->SetLineWidth(2);
  unhist->SetFillColor(0);  
  //unhist->Draw("SAME");

  mjhist_av_1->Draw("hist same");
  mjhist_av_2->Draw("hist same");
  mjhist_av_3->Draw("hist same");

  if(TString(postfix).Contains("Vector")){
    leg->AddEntry(mjhist_av_1,"Vector m_{med} = 1.5 TeV","L");
    leg->AddEntry(mjhist_av_2,"Vector m_{med} = 2.0 TeV","L");
    leg->AddEntry(mjhist_av_3,"Vector m_{med} = 2.5 TeV","L");
  }
  else if(TString(postfix).Contains("Scalar")){
    leg->AddEntry(mjhist_av_1,"Scalar m_{med} = 125 GeV","L");
    leg->AddEntry(mjhist_av_2,"Scalar m_{med} = 400 GeV","L");
    leg->AddEntry(mjhist_av_3,"Scalar m_{med} = 600 GeV","L");
  }

  leg->Draw("same");
  canvas->RedrawAxis("sameaxis");
  canvas->SaveAs(("significance_"+postfix+".pdf").c_str());
  canvas->SaveAs(("significance_"+postfix+".png").c_str());

  frame->Draw();
  CMS_lumi(canvas,"35.9");
  categoryLabel->Draw("same");
  frame->GetYaxis()->SetTitle("S/B a.u.");
  frame->GetYaxis()->SetRangeUser(1.,max(mjhist_av_1_norm->GetMaximum(),max(mjhist_av_2_norm->GetMaximum(),mjhist_av_3_norm->GetMaximum()))*1.05);
  mjhist_av_1_norm->Draw("hist same");
  mjhist_av_2_norm->Draw("hist same");
  mjhist_av_3_norm->Draw("hist same");
  leg->Draw("same");
  canvas->RedrawAxis("sameaxis");
  canvas->SaveAs(("significance_norm_"+postfix+".pdf").c_str());
  canvas->SaveAs(("significance_norm_"+postfix+".png").c_str());

}



