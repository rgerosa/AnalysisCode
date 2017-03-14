#include "../CMS_lumi.h"

void makeTriggerMETPlotComparison(string inputFileDataWmn, string inputFileDataZmm, string inputFileMC,string outputDIR, bool isWen){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  TFile* inputFile1 = TFile::Open(inputFileDataWmn.c_str());
  TFile* inputFile2 = TFile::Open(inputFileDataZmm.c_str());
  TFile* inputFile3 = NULL;
  if(inputFileMC != "")
    inputFile3 = TFile::Open(inputFileMC.c_str());
    
  TGraphAsymmErrors* graph_file1 = (TGraphAsymmErrors*) ((TEfficiency*) inputFile1->Get("efficiency"))->CreateGraph();
  TGraphAsymmErrors* graph_file2 = (TGraphAsymmErrors*) ((TEfficiency*) inputFile2->Get("efficiency"))->CreateGraph();

  TGraphAsymmErrors* graph_file3 = NULL;
  TGraphAsymmErrors* graph_file4 = NULL;
  if(inputFile3)
    graph_file3 = (TGraphAsymmErrors*) ((TEfficiency*) inputFile3->Get("efficiency_wmn"))->CreateGraph();
  if(inputFile3 and not isWen)
    graph_file4 = (TGraphAsymmErrors*) ((TEfficiency*) inputFile3->Get("efficiency_zmm"))->CreateGraph();
  else if(inputFile3 and isWen)
    graph_file4 = (TGraphAsymmErrors*) ((TEfficiency*) inputFile3->Get("efficiency_wen"))->CreateGraph();

  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();
  graph_file1->GetXaxis()->SetTitle("Recoil [GeV]");
  graph_file1->GetYaxis()->SetTitle("Trigger Efficiency");
  graph_file1->SetMarkerColor(kBlack);
  graph_file1->SetLineColor(kBlack);
  graph_file1->Draw("AEPL");
 
  graph_file2->SetMarkerColor(kRed);
  graph_file2->SetLineColor(kRed);
  graph_file2->Draw("EPLsame");

  if(graph_file3){
    graph_file3->SetMarkerColor(kBlue);
    graph_file3->SetLineColor(kBlue);
    graph_file3->Draw("EPLsame");
  }

  if(graph_file4){
    graph_file4->SetMarkerColor(kCyan);
    graph_file4->SetLineColor(kCyan);
    graph_file4->Draw("EPLsame");
  }

  CMS_lumi(canvas,"35.9");
  TLegend leg (0.5,0.35,0.8,0.55);
  leg.AddEntry(graph_file1,"W #rightarrow #mu#nu Data","EPL");
  if(isWen)
    leg.AddEntry(graph_file2,"W #rightarrow e#nu Data","EPL");
  else
    leg.AddEntry(graph_file2,"Z #rightarrow #mu#mu Data","EPL");

  if(graph_file3) 
    leg.AddEntry(graph_file3,"W #rightarrow #mu#nu MC","EPL");
  if(graph_file4 and not isWen) 
    leg.AddEntry(graph_file4,"Z #rightarrow #mu#mu MC","EPL");
  else if(graph_file4 and isWen) 
    leg.AddEntry(graph_file4,"W #rightarrow e#nu MC","EPL");


  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/comparisonTurnOn.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/comparisonTurnOn.pdf").c_str(),"pdf");

  graph_file1->GetXaxis()->SetRangeUser(150,700);
  graph_file1->GetYaxis()->SetRangeUser(0.85,1.05);

  canvas->SaveAs((outputDIR+"/comparisonTurnOn_zoom.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/comparisonTurnOn_zoom.pdf").c_str(),"pdf");

  TGraphAsymmErrors* ratio_band_up_data = new TGraphAsymmErrors();
  TGraphAsymmErrors* ratio_band_dw_data = new TGraphAsymmErrors();
  for(int iPoint = 0; iPoint < graph_file3->GetN(); iPoint++){
    double x,y;
    graph_file3->GetPoint(iPoint,x,y);
    ratio_band_up_data->SetPoint(iPoint,x,graph_file1->Eval(x)/graph_file3->Eval(x));
    ratio_band_dw_data->SetPoint(iPoint,x,graph_file1->Eval(x)/graph_file3->Eval(x));
  }
  ratio_band_up_data->GetXaxis()->SetTitle("Recoil [GeV]");
  ratio_band_up_data->GetYaxis()->SetTitle("Scale Factor");
  ratio_band_up_data->GetXaxis()->SetRangeUser(200,600);
  ratio_band_up_data->GetYaxis()->SetRangeUser(0.9,1.05);
  ratio_band_up_data->SetMarkerColor(kBlack);
  ratio_band_up_data->SetLineColor(kBlack);
  ratio_band_up_data->SetLineWidth(2);
  ratio_band_dw_data->SetLineWidth(2);
  ratio_band_up_data->SetMarkerSize(1);
  ratio_band_up_data->SetMarkerStyle(20);
  ratio_band_up_data->Draw("AEPL");
  //ratio_band_dw_data->Draw("EPLsame");
  
  TGraphAsymmErrors* ratio_band_up_mc = new TGraphAsymmErrors();
  TGraphAsymmErrors* ratio_band_dw_mc = new TGraphAsymmErrors();
  if(graph_file3 and graph_file4){
    for(int iPoint = 0; iPoint < graph_file4->GetN(); iPoint++){
      double x,y;
      graph_file4->GetPoint(iPoint,x,y);
      ratio_band_up_mc->SetPoint(iPoint,x,graph_file2->Eval(x)/graph_file4->Eval(x));
      ratio_band_dw_mc->SetPoint(iPoint,x,graph_file2->Eval(x)/graph_file4->Eval(x));
    }
    ratio_band_up_mc->SetMarkerSize(1);
    ratio_band_up_mc->SetMarkerStyle(20);
    ratio_band_up_mc->Draw("EPLsame");
    //ratio_band_dw_mc->Draw("EPLsame");
  }

  TLegend leg2 (0.2,0.25,0.55,0.5);
  leg2.AddEntry(ratio_band_up_data,"W #rightarrow #mu#nu (Data) / W #rightarrow #mu#nu (MC) ","L");
  if(isWen)
    leg2.AddEntry(ratio_band_up_mc,"W #rightarrow e#nu (Data) / W #rightarrow e#nu (MC) ","L");
  else
    leg2.AddEntry(ratio_band_up_mc,"Z #rightarrow #mu#mu (Data) / Z #rightarrow #mu#mu (MC) ","L");

  leg2.Draw("same");

  ratio_band_up_mc->GetYaxis()->SetRangeUser(0.9,1.05);
  ratio_band_up_mc->SetMarkerColor(kRed);
  ratio_band_up_mc->SetLineColor(kRed);
  ratio_band_up_mc->SetLineWidth(2);
  ratio_band_dw_mc->SetMarkerColor(kRed);
  ratio_band_dw_mc->SetLineColor(kRed);
  ratio_band_dw_mc->SetLineWidth(2);
  CMS_lumi(canvas,"35.9");

  canvas->SaveAs((outputDIR+"/ratio_turnon.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/ratio_turnon.pdf").c_str(),"pdf");

  TFile* outputFile = new TFile((outputDIR+"/metTriggerScaleFactor.root").c_str(),"RECREATE");
  outputFile->cd();
  ratio_band_up_data->GetXaxis()->SetRangeUser(0.,1500);
  ratio_band_up_data->GetYaxis()->SetRangeUser(0.,1.05);
  ratio_band_up_data->Write("scalefactor_wmn");
  ratio_band_up_mc->GetXaxis()->SetRangeUser(0.,1500);
  ratio_band_up_mc->GetYaxis()->SetRangeUser(0.,1.05);
  if(not isWen)
    ratio_band_up_mc->Write("scalefactor_zmm");
  else
    ratio_band_up_mc->Write("scalefactor_wen");

  TGraphAsymmErrors* eff_zvv = (TGraphAsymmErrors*) ((TEfficiency*) inputFile3->Get("efficiency_zvv"))->CreateGraph();
  TGraphAsymmErrors* eff_wjet = (TGraphAsymmErrors*) ((TEfficiency*) inputFile3->Get("efficiency_wjet"))->CreateGraph();
  eff_zvv->GetYaxis()->SetRangeUser(0.,1.05);
  eff_zvv->Write("efficiency_zvv");
  eff_wjet->GetYaxis()->SetRangeUser(0.,1.05);
  eff_wjet->Write("efficiency_wjet");
  graph_file1->GetXaxis()->SetRangeUser(0.,1500);
  graph_file1->GetYaxis()->SetRangeUser(0.,1.05);
  graph_file1->Write("efficiency_wmn_data");
  graph_file2->GetXaxis()->SetRangeUser(0.,1500);
  graph_file2->GetYaxis()->SetRangeUser(0.,1.05);
  if(not isWen)
    graph_file2->Write("efficiency_zmm_data");
  else
    graph_file2->Write("efficiency_wen_data");
  outputFile->Close();
}
