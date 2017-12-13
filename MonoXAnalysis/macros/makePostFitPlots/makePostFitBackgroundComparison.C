#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void plotComparison(TH1* histo_1, TH1* histo_2, const string & observable, const Category & category, const string & postfix){

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


  canvas->cd();
  TH1* frame = (TH1*) histo_1->Clone("frame");
  frame->Reset();
  frame->SetLineColor(kBlack);
  frame->SetLineWidth(1);

  if(category == Category::monojet)
    frame->GetYaxis()->SetRangeUser(0.002,histo_1->GetMaximum()*500);
  else if(category == Category::monoV)
    frame->GetYaxis()->SetRangeUser(0.01,histo_1->GetMaximum()*500);
  else if(category == Category::VBF)
    frame->GetYaxis()->SetRangeUser(0.015,histo_1->GetMaximum()*500);
  else if(category == Category::VBFrelaxed)
    frame->GetYaxis()->SetRangeUser(0.0007,histo_1->GetMaximum()*500);

  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("Events / GeV");
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
  else if(category == Category::VBF or category == Category::VBFrelaxed)
    categoryLabel ->DrawLatex(0.175,0.80,"VBF");
  categoryLabel->Draw("same");

  histo_1->Draw("hist same");
  histo_2->Draw("hist same");

  TLegend* leg =  new TLegend(0.6,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histo_1,(postfix+" CR").c_str(),"L");
  leg->AddEntry(histo_2,(postfix+" CR+SR").c_str(),"L");
  leg->Draw("SAME");    
  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();
  canvas->cd();

  pad2->Draw();
  pad2->cd();

  ////
  TH1* frame2 =  (TH1*) histo_1->Clone("frame");
  frame2->Reset();
  frame2->SetLineColor(kBlack);
  frame2->SetLineWidth(1);

  if(category == Category::monojet)
    frame2->GetYaxis()->SetRangeUser(0.7,1.3);
  else if(category == Category::monoV)
    frame2->GetYaxis()->SetRangeUser(0.7,1.3);
  else
    frame2->GetYaxis()->SetRangeUser(0.7,1.3);

  if(category == Category::monojet)
    frame2->GetXaxis()->SetNdivisions(510);
  else if(category == Category::VBFrelaxed)
    frame2->GetXaxis()->SetNdivisions(505);
  else
    frame2->GetXaxis()->SetNdivisions(210);
  frame2->GetYaxis()->SetNdivisions(5);

  frame2->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  if((category == Category::VBF or category == Category::VBFrelaxed) and TString(observable).Contains("mjj"))
    frame2->GetXaxis()->SetTitle("M_{jj} [GeV]");
  frame2->GetYaxis()->SetTitle("CR-only/(CR+SR)");
  frame2->GetYaxis()->CenterTitle();
  frame2->GetYaxis()->SetTitleOffset(1.5);
  frame2->GetYaxis()->SetLabelSize(0.04);
  frame2->GetYaxis()->SetTitleSize(0.04);
  frame2->GetXaxis()->SetLabelSize(0.04);
  frame2->GetXaxis()->SetTitleSize(0.05);
  frame2->GetXaxis()->SetTickLength(0.025);
  frame2->Draw();

  TH1F* ratio = (TH1F*) histo_1->Clone("ratio");
  ratio->Divide(histo_2);
  ratio->SetLineColor(kBlack);
  ratio->SetMarkerColor(kBlack);
  ratio->SetMarkerSize(1);
  ratio->SetMarkerStyle(20);
  ratio->Draw("PEsame");

  canvas->SaveAs(("comparison_"+postfix+".pdf").c_str());
  canvas->SaveAs(("comparison_"+postfix+".png").c_str());
}

void makePostFitBackgroundComparison(string   fileName_crOnly,
				     string   fileName_bOnly,
				     string   observable, 
				     Category category, 
				     bool     plotSBFit = false,
				     bool     isCombinedFit = false
				     ){
  
  
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  TColor *color; // for color definition with alpha                                                                                                                             

  TFile* file_crOnly = new TFile(fileName_crOnly.c_str());
  TFile* file_bOnly = new TFile(fileName_bOnly.c_str());

  string fit_dir = "shapes_fit_b";
  if(plotSBFit)
    fit_dir = "shapes_fit_s";

  string dir;
  if(isCombinedFit){
    if(category == Category::monojet)
      dir = "ch1_ch1";
    else if(category == Category::monoV)
      dir = "ch2_ch1";
    else if(category == Category::VBF or category == Category::VBFrelaxed)
      dir = "ch3_ch1";
  }
  else if( category != Category::VBF and category != Category::VBFrelaxed)
    dir = "ch1";
  else
    dir = "ch1";

  string postfix = "_MJ";
  if(category == Category::monoV)
    postfix = "_MV";
  else if(category == Category::VBF or category == Category::VBFrelaxed)
    postfix = "_VBF";

  TH1* zvvhist_1 = NULL;
  TH1* zvvhist_2 = NULL;
  TH1* zvvewkhist_1 = NULL;
  TH1* zvvewkhist_2 = NULL;
  TH1* wjethist_1 = NULL;
  TH1* wjethist_2 = NULL;
  TH1* wjetewkhist_1 = NULL;
  TH1* wjetewkhist_2 = NULL;

  zvvhist_1 = (TH1*) file_crOnly->Get((fit_dir+"/"+dir+"/qcd_znunu").c_str());
  zvvewkhist_1 = (TH1*) file_crOnly->Get((fit_dir+"/"+dir+"/ewk_znunu").c_str());
  wjethist_1 = (TH1*) file_crOnly->Get((fit_dir+"/"+dir+"/qcd_wjets").c_str());
  wjetewkhist_1 = (TH1*) file_crOnly->Get((fit_dir+"/"+dir+"/ewk_wjets").c_str());

  zvvhist_2 = (TH1*) file_bOnly->Get((fit_dir+"/"+dir+"/qcd_znunu").c_str());
  zvvewkhist_2 = (TH1*) file_bOnly->Get((fit_dir+"/"+dir+"/ewk_znunu").c_str());
  wjethist_2 = (TH1*) file_bOnly->Get((fit_dir+"/"+dir+"/qcd_wjets").c_str());
  wjetewkhist_2 = (TH1*) file_bOnly->Get((fit_dir+"/"+dir+"/ewk_wjets").c_str());


  zvvhist_1->SetLineColor(kRed);
  zvvewkhist_1->SetLineColor(kRed);
  wjethist_1->SetLineColor(kRed);
  wjetewkhist_1->SetLineColor(kRed);

  zvvhist_2->SetLineColor(kBlue);
  zvvewkhist_2->SetLineColor(kBlue);
  wjethist_2->SetLineColor(kBlue);
  wjetewkhist_2->SetLineColor(kBlue);

  zvvhist_1->SetLineWidth(2);
  zvvewkhist_1->SetLineWidth(2);
  wjethist_1->SetLineWidth(2);
  wjetewkhist_1->SetLineWidth(2);

  zvvhist_2->SetLineWidth(2);
  zvvewkhist_2->SetLineWidth(2);
  wjethist_2->SetLineWidth(2);
  wjetewkhist_2->SetLineWidth(2);

  plotComparison(zvvhist_1,zvvhist_2,observable,category,"Zvv-QCD");
  plotComparison(zvvewkhist_1,zvvewkhist_2,observable,category,"Zvv-EW");
  plotComparison(wjethist_1,wjethist_2,observable,category,"Wjets-QCD");
  plotComparison(wjetewkhist_1,wjetewkhist_2,observable,category,"Wjets-EW");
  
}

