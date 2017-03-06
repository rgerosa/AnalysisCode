#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void drawPlot(TH1* template_mj_1, TH1* template_mj_2, TH1* template_mj_3, string outputDIR, float lumi, string model, Category category, string postfix);

void fixZeroBins (TH1* histo){
  for(int iBin = 0; iBin < histo->GetNbinsX()+1; iBin++){
    if(histo->GetBinContent(iBin+1) == 0) histo->SetBinContent(iBin+1,0.001);
  }
}

void plotSignalPrefit (string model, Category category, string outputDIR, float lumi){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputDIR).c_str());

  TFile* monoJfile = NULL;
  TFile* monoWfile = NULL;
  TFile* monoZfile = NULL;

  TH1* template_mj_1 = NULL;
  TH1* template_mj_2 = NULL;
  TH1* template_mj_3 = NULL;
  TH1* template_mv_1 = NULL;
  TH1* template_mv_2 = NULL;
  TH1* template_mv_3 = NULL;

  if(model == "Vector" and category ==  Category::monojet){
    monoJfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_800_0.25_catmonojet_13TeV_v1.root","READ");
    monoWfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_800_0.25_catmonojet_13TeV_v1.root","READ");
    monoZfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_800_0.25_catmonojet_13TeV_v1.root","READ");

    template_mj_1 = (TH1*) monoJfile->FindObjectAny("signal_signal_80015250001"); 
    template_mj_2 = (TH1*) monoJfile->FindObjectAny("signal_signal_80018000001"); 
    template_mj_3 = (TH1*) monoJfile->FindObjectAny("signal_signal_80020000001"); 

    fixZeroBins(template_mj_1);
    fixZeroBins(template_mj_2);
    fixZeroBins(template_mj_3);

    template_mv_1 = (TH1*) monoWfile->FindObjectAny("signal_signal_80015250001"); 
    template_mv_2 = (TH1*) monoWfile->FindObjectAny("signal_signal_80018000001"); 
    template_mv_3 = (TH1*) monoWfile->FindObjectAny("signal_signal_80020000001"); 

    template_mv_1->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80015250001"));
    template_mv_2->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80018000001"));
    template_mv_3->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80020000001"));

    fixZeroBins(template_mv_1);
    fixZeroBins(template_mv_2);
    fixZeroBins(template_mv_3);

  }
  else if(model == "Vector" and category ==  Category::monoV){
    monoJfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_800_0.25_catmonov_13TeV_v1.root","READ");
    monoWfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_800_0.25_catmonov_13TeV_v1.root","READ");
    monoZfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_800_0.25_catmonov_13TeV_v1.root","READ");

    template_mj_1 = (TH1*) monoJfile->FindObjectAny("signal_signal_80015250001"); 
    template_mj_2 = (TH1*) monoJfile->FindObjectAny("signal_signal_80018000001"); 
    template_mj_3 = (TH1*) monoJfile->FindObjectAny("signal_signal_80020000001"); 

    fixZeroBins(template_mj_1);
    fixZeroBins(template_mj_2);
    fixZeroBins(template_mj_3);

    template_mv_1 = (TH1*) monoWfile->FindObjectAny("signal_signal_80015250001"); 
    template_mv_2 = (TH1*) monoWfile->FindObjectAny("signal_signal_80018000001"); 
    template_mv_3 = (TH1*) monoWfile->FindObjectAny("signal_signal_80020000001"); 

    template_mv_1->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80015250001"));
    template_mv_2->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80018000001"));
    template_mv_3->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80020000001"));

    fixZeroBins(template_mv_1);
    fixZeroBins(template_mv_2);
    fixZeroBins(template_mv_3);

  }


  if(model == "Axial" and category ==  Category::monojet){
    monoJfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_801_0.25_catmonojet_13TeV_v1.root","READ");
    monoWfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_801_0.25_catmonojet_13TeV_v1.root","READ");
    monoZfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_801_0.25_catmonojet_13TeV_v1.root","READ");

    template_mj_1 = (TH1*) monoJfile->FindObjectAny("signal_signal_80115250001"); 
    template_mj_2 = (TH1*) monoJfile->FindObjectAny("signal_signal_80118000001"); 
    template_mj_3 = (TH1*) monoJfile->FindObjectAny("signal_signal_80120000001"); 

    fixZeroBins(template_mj_1);
    fixZeroBins(template_mj_2);
    fixZeroBins(template_mj_3);

    template_mv_1 = (TH1*) monoWfile->FindObjectAny("signal_signal_80115250001"); 
    template_mv_2 = (TH1*) monoWfile->FindObjectAny("signal_signal_80118000001"); 
    template_mv_3 = (TH1*) monoWfile->FindObjectAny("signal_signal_80120000001"); 

    template_mv_1->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80115250001"));
    template_mv_2->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80118000001"));
    template_mv_3->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80120000001"));

    fixZeroBins(template_mv_1);
    fixZeroBins(template_mv_2);
    fixZeroBins(template_mv_3);

  }
  else if(model == "Axial" and category ==  Category::monoV){
    monoJfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_801_0.25_catmonov_13TeV_v1.root","READ");
    monoWfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_801_0.25_catmonov_13TeV_v1.root","READ");
    monoZfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_801_0.25_catmonov_13TeV_v1.root","READ");

    template_mj_1 = (TH1*) monoJfile->FindObjectAny("signal_signal_80115250001"); 
    template_mj_2 = (TH1*) monoJfile->FindObjectAny("signal_signal_80118000001"); 
    template_mj_3 = (TH1*) monoJfile->FindObjectAny("signal_signal_80120000001"); 

    fixZeroBins(template_mj_1);
    fixZeroBins(template_mj_2);
    fixZeroBins(template_mj_3);

    template_mv_1 = (TH1*) monoWfile->FindObjectAny("signal_signal_80115250001"); 
    template_mv_2 = (TH1*) monoWfile->FindObjectAny("signal_signal_80118000001"); 
    template_mv_3 = (TH1*) monoWfile->FindObjectAny("signal_signal_80120000001"); 

    template_mv_1->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80115250001"));
    template_mv_2->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80118000001"));
    template_mv_3->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80120000001"));

    fixZeroBins(template_mv_1);
    fixZeroBins(template_mv_2);
    fixZeroBins(template_mv_3);

  }

  ////
  if(model == "Scalar" and category ==  Category::monojet){
    monoJfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_805_1.0_catmonojet_13TeV_v1.root","READ");
    monoWfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_805_1.0_catmonojet_13TeV_v1.root","READ");
    monoZfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_805_1.0_catmonojet_13TeV_v1.root","READ");

    template_mj_1 = (TH1*) monoJfile->FindObjectAny("signal_signal_80500100001"); 
    template_mj_2 = (TH1*) monoJfile->FindObjectAny("signal_signal_80500500001"); 
    template_mj_3 = (TH1*) monoJfile->FindObjectAny("signal_signal_80501000001"); 

    fixZeroBins(template_mj_1);
    fixZeroBins(template_mj_2);
    fixZeroBins(template_mj_3);

    template_mv_1 = (TH1*) monoWfile->FindObjectAny("signal_signal_80500100001"); 
    template_mv_2 = (TH1*) monoWfile->FindObjectAny("signal_signal_80500500001"); 
    template_mv_3 = (TH1*) monoWfile->FindObjectAny("signal_signal_80501000001"); 

    template_mv_1->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80500100001"));
    template_mv_2->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80500500001"));
    template_mv_3->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80501000001"));

    fixZeroBins(template_mv_1);
    fixZeroBins(template_mv_2);
    fixZeroBins(template_mv_3);

  }
  else if(model == "Scalar" and category ==  Category::monoV){
    monoJfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_805_1.0_catmonov_13TeV_v1.root","READ");
    monoWfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoW_805_1.0_catmonov_13TeV_v1.root","READ");
    monoZfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoZ_805_1.0_catmonov_13TeV_v1.root","READ");

    template_mj_1 = (TH1*) monoJfile->FindObjectAny("signal_signal_80500100001"); 
    template_mj_2 = (TH1*) monoJfile->FindObjectAny("signal_signal_80500500001"); 
    template_mj_3 = (TH1*) monoJfile->FindObjectAny("signal_signal_80501000001"); 

    fixZeroBins(template_mj_1);
    fixZeroBins(template_mj_2);
    fixZeroBins(template_mj_3);

    template_mv_1 = (TH1*) monoWfile->FindObjectAny("signal_signal_80500100001"); 
    template_mv_2 = (TH1*) monoWfile->FindObjectAny("signal_signal_80500500001"); 
    template_mv_3 = (TH1*) monoWfile->FindObjectAny("signal_signal_80501000001"); 

    template_mv_1->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80500100001"));
    template_mv_2->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80500500001"));
    template_mv_3->Add((TH1*) monoZfile->FindObjectAny("signal_signal_80501000001"));

    fixZeroBins(template_mv_1);
    fixZeroBins(template_mv_2);
    fixZeroBins(template_mv_3);

  }


  ////
  if(model == "PseudoScalar" and category ==  Category::monojet){

    monoJfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_806_1.0_catmonojet_13TeV_v1.root","READ");

    template_mj_1 = (TH1*) monoJfile->FindObjectAny("signal_signal_80600100001"); 
    template_mj_2 = (TH1*) monoJfile->FindObjectAny("signal_signal_80601000001"); 
    template_mj_3 = (TH1*) monoJfile->FindObjectAny("signal_signal_80602000001"); 

    fixZeroBins(template_mj_1);
    fixZeroBins(template_mj_2);
    fixZeroBins(template_mj_3);

  }
  else if(model == "PseudoScalar" and category ==  Category::monoV){
    monoJfile = new TFile("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Signal_v3/MonoJ_806_1.0_catmonov_13TeV_v1.root","READ");
    template_mj_1 = (TH1*) monoJfile->FindObjectAny("signal_signal_80600100001"); 
    template_mj_2 = (TH1*) monoJfile->FindObjectAny("signal_signal_80601000001"); 
    template_mj_3 = (TH1*) monoJfile->FindObjectAny("signal_signal_80602000001"); 

    fixZeroBins(template_mj_1);
    fixZeroBins(template_mj_2);
    fixZeroBins(template_mj_3);

  }

  
  drawPlot(template_mj_1,template_mj_2,template_mj_3,outputDIR,lumi,model,category,"monoJ");
  if(model != "PseudoScalar")
    drawPlot(template_mv_1,template_mv_2,template_mv_3,outputDIR,lumi,model,category,"monoV");

}


/////////////////
void drawPlot(TH1* template_mj_1, TH1* template_mj_2, TH1* template_mj_3, string outputDIR, float lumi, string model, Category category, string postfix){

  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->SetBottomMargin(0.3);

  TH1F* frame = (TH1F*) template_mj_1->Clone("frame");
  frame->Reset();
  TH1* frame2 =  (TH1*) frame->Clone("frame2");
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("Events");
  frame->SetFillColor(0);
  frame->SetFillStyle(0);
  frame->GetYaxis()->SetRangeUser(template_mj_3->GetMinimum()*0.5,template_mj_1->GetMaximum()*100);
  frame->Draw();
  CMS_lumi(canvas,Form("%1.f",lumi));
  
  TLegend leg (0.65,0.65,0.9,0.9);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);

  if(template_mj_1 != NULL){
    template_mj_1->Scale(lumi/12.9);
    template_mj_1->SetLineColor(kBlack);
    template_mj_1->SetLineWidth(2);
    template_mj_1->Draw("hist same");
    if(model == "Vector" or model == "Axial")
      leg.AddEntry(template_mj_1,"V/AV m_{med} = 1.5 TeV","L");
    else if (model == "Scalar" or model == "PseudoScalar")
      leg.AddEntry(template_mj_1,"S/PS m_{med} = 10 GeV","L");
  }
  if(template_mj_2 != NULL){
    template_mj_2->Scale(lumi/12.9);
    template_mj_2->SetLineColor(kRed);
    template_mj_2->SetLineWidth(2);
    template_mj_2->Draw("hist same");
    if(model == "Vector" or model == "Axial")
      leg.AddEntry(template_mj_2,"V/AV m_{med} = 1.8 TeV","L");
    else if (model == "Scalar" or model == "PseudoScalar")
      leg.AddEntry(template_mj_2,"S/PS m_{med} = 50 GeV","L");
  }
  if(template_mj_3 != NULL){
    template_mj_3->Scale(lumi/12.9);
    template_mj_3->SetLineColor(kBlue);
    template_mj_3->SetLineWidth(2);
    template_mj_3->Draw("hist same");
    if(model == "Vector" or model == "Axial")
      leg.AddEntry(template_mj_3,"V/AV m_{med} = 2.0 TeV","L");
    else if (model == "Scalar" or model == "PseudoScalar")
      leg.AddEntry(template_mj_3,"S/PS m_{med} = 100 GeV","L");
  }
  canvas->SetLogy();
  leg.Draw("same");

  canvas->cd();
  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,1);
  pad2->SetTopMargin(0.72);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);
  pad2->Draw();
  pad2->cd();
  
  frame2->GetYaxis()->SetTitle("Ratio");
  frame2->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  frame2->GetXaxis()->SetTitleOffset(1.1);
  frame2->GetYaxis()->SetTitleOffset(1.2);
  frame2->GetYaxis()->SetLabelSize(0.04);
  frame2->GetYaxis()->SetTitleSize(0.04);
  frame2->GetYaxis()->SetNdivisions(5);
  frame2->Draw();

  TH1* ratio_1 = (TH1*) template_mj_2->Clone("ratio_1");
  TH1* ratio_2 = (TH1*) template_mj_3->Clone("ratio_2");
  ratio_1->Divide(template_mj_1);
  ratio_2->Divide(template_mj_1);
  ratio_1->Draw("hist same");
  ratio_2->Draw("hist same");
  frame2->GetYaxis()->SetRangeUser(0.5,1.2);
  
  if(category == Category::monoV) postfix += "_catMV";
  else if(category == Category::monojet) postfix += "_catMJ";

  canvas->SaveAs((outputDIR+"/spectrum_"+model+"_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/spectrum_"+model+"_"+postfix+".pdf").c_str(),"pdf");

  cout<<"integral "<<template_mj_1->Integral()<<" "<<template_mj_2->Integral()<<" "<<template_mj_3->Integral()<<endl;

}
