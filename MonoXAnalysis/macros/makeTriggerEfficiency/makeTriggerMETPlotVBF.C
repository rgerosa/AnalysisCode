#include "../CMS_lumi.h"

void makeTriggerPlotVBF (string observable){

  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  TFile* file_mjj_bin1 = NULL;
  TFile* file_mjj_bin2 = NULL;
  TFile* file_mjj_bin3 = NULL;
  TFile* file_mjj_bin4 = NULL;
    
  if(observable == "mjj"){
    file_mjj_bin1 = TFile::Open(("triggerEfficiencyMET_muon/metTriggerEfficiency_"+observable+"_vbf_0.0_800.0.root").c_str());
    file_mjj_bin2 = TFile::Open(("triggerEfficiencyMET_muon/metTriggerEfficiency_"+observable+"_vbf_800.0_1200.0.root").c_str());
    file_mjj_bin3 = TFile::Open(("triggerEfficiencyMET_muon/metTriggerEfficiency_"+observable+"_vbf_1200.0_1700.0.root").c_str());
    file_mjj_bin4 = TFile::Open(("triggerEfficiencyMET_muon/metTriggerEfficiency_"+observable+"_vbf_1700.0_3000.0.root").c_str());
  }
  else if(observable == "detajj"){
    file_mjj_bin1 = TFile::Open(("triggerEfficiencyMET_muon/metTriggerEfficiency_"+observable+"_vbf_0.0_1.5.root").c_str());
    file_mjj_bin2 = TFile::Open(("triggerEfficiencyMET_muon/metTriggerEfficiency_"+observable+"_vbf_1.5_3.0.root").c_str());
    file_mjj_bin3 = TFile::Open(("triggerEfficiencyMET_muon/metTriggerEfficiency_"+observable+"_vbf_3.0_5.0.root").c_str());
    file_mjj_bin4 = TFile::Open(("triggerEfficiencyMET_muon/metTriggerEfficiency_"+observable+"_vbf_5.0_9.0.root").c_str());
  }

  TEfficiency* mjj_bin1 = (TEfficiency*) file_mjj_bin1->Get("efficiency");
  TEfficiency* mjj_bin2 = (TEfficiency*) file_mjj_bin2->Get("efficiency");
  TEfficiency* mjj_bin3 = (TEfficiency*) file_mjj_bin3->Get("efficiency");
  TEfficiency* mjj_bin4 = (TEfficiency*) file_mjj_bin4->Get("efficiency");

  mjj_bin1->SetMarkerColor(kBlack);
  mjj_bin2->SetMarkerColor(kBlue);
  mjj_bin3->SetMarkerColor(kRed);
  mjj_bin4->SetMarkerColor(kGreen+1);

  mjj_bin1->SetLineColor(kBlack);
  mjj_bin2->SetLineColor(kBlue);
  mjj_bin3->SetLineColor(kRed);
  mjj_bin4->SetLineColor(kGreen+1);

  TF1* fit_mjj_bin1 = (TF1*) file_mjj_bin1->Get("efficiency_func");
  TF1* fit_mjj_bin2 = (TF1*) file_mjj_bin2->Get("efficiency_func");
  TF1* fit_mjj_bin3 = (TF1*) file_mjj_bin3->Get("efficiency_func");
  TF1* fit_mjj_bin4 = (TF1*) file_mjj_bin4->Get("efficiency_func");

  fit_mjj_bin1->SetLineColor(kBlack);
  fit_mjj_bin2->SetLineColor(kBlue);
  fit_mjj_bin3->SetLineColor(kRed);
  fit_mjj_bin4->SetLineColor(kGreen+1);
  
  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();

  TH1* frame = canvas->DrawFrame(fit_mjj_bin1->GetXaxis()->GetXmin(),0.,fit_mjj_bin1->GetXaxis()->GetXmax(), 1.1, "");
  frame->GetXaxis()->SetTitle("Recoil [GeV]");
  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetYaxis()->SetLabelSize(0.8*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetTitleSize(0.8*frame->GetYaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleOffset(1.0);
  frame->Draw();

  mjj_bin1->Draw("EPsame");  
  fit_mjj_bin1->Draw("Lsame");
  mjj_bin2->Draw("EPsame");
  fit_mjj_bin2->Draw("Lsame");
  mjj_bin3->Draw("EPsame");
  fit_mjj_bin3->Draw("Lsame");
  mjj_bin4->Draw("EPsame");
  fit_mjj_bin4->Draw("Lsame");

  CMS_lumi(canvas,"35.9",true);

  TLegend* leg = new TLegend(0.5,0.4,0.8,0.6);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  if(observable == "mjj"){
    leg->AddEntry(mjj_bin1,"0.0 < m_{jj} < 0.8 GeV","EP");
    leg->AddEntry(mjj_bin2,"0.8 < m_{jj} < 1.2 GeV","EP");
    leg->AddEntry(mjj_bin3,"1.2 < m_{jj} < 1.7 GeV","EP");
    leg->AddEntry(mjj_bin4,"1.7 < m_{jj} < 3.0 GeV","EP");
  }
  if(observable == "detajj"){
    leg->AddEntry(mjj_bin1,"0.0 < #Delta#eta_{jj} < 1.5","EP");
    leg->AddEntry(mjj_bin2,"1.5 < #Delta#eta_{jj} < 3.0","EP");
    leg->AddEntry(mjj_bin3,"3.0 < #Delta#eta_{jj} < 5.0","EP");
    leg->AddEntry(mjj_bin4,"5.0 < #Delta#eta_{jj} < 9.0","EP");
  }
  leg->Draw("same");

  canvas->SaveAs(("triggerPlotVBF"+observable+".pdf").c_str(),"pdf");
  canvas->SaveAs(("triggerPlotVBF"+observable+".png").c_str(),"png");


}
