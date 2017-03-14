#include "../CMS_lumi.h"

void makeTriggerMETDistributionComparison(string inputFileName, string outputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TFile* inputFile = TFile::Open(inputFileName.c_str());
  
  TH1* caloMet_SR_pas = (TH1*) inputFile->FindObjectAny("metCalo_SR_pas");
  TH1* caloMet_WM_pas = (TH1*) inputFile->FindObjectAny("metCalo_WM_pas");
  TH1* caloMet_ZM_pas = (TH1*) inputFile->FindObjectAny("metCalo_ZM_pas");

  caloMet_SR_pas->Scale(1./caloMet_SR_pas->Integral());
  caloMet_WM_pas->Scale(1./caloMet_WM_pas->Integral());
  caloMet_ZM_pas->Scale(1./caloMet_ZM_pas->Integral());

  TH1* pfMet_SR_pas = (TH1*) inputFile->FindObjectAny("metPF_SR_pas");
  TH1* pfMet_WM_pas = (TH1*) inputFile->FindObjectAny("metPF_WM_pas");
  TH1* pfMet_ZM_pas = (TH1*) inputFile->FindObjectAny("metPF_ZM_pas");

  pfMet_SR_pas->Scale(1./pfMet_SR_pas->Integral());
  pfMet_WM_pas->Scale(1./pfMet_WM_pas->Integral());
  pfMet_ZM_pas->Scale(1./pfMet_ZM_pas->Integral());

  TH1* caloMet_SR_fail = (TH1*) inputFile->FindObjectAny("metCalo_SR_fail");
  TH1* caloMet_WM_fail = (TH1*) inputFile->FindObjectAny("metCalo_WM_fail");
  TH1* caloMet_ZM_fail = (TH1*) inputFile->FindObjectAny("metCalo_ZM_fail");

  caloMet_SR_fail->Scale(1./caloMet_SR_fail->Integral());
  caloMet_WM_fail->Scale(1./caloMet_WM_fail->Integral());
  caloMet_ZM_fail->Scale(1./caloMet_ZM_fail->Integral());

  TH1* pfMet_SR_fail = (TH1*) inputFile->FindObjectAny("metPF_SR_fail");
  TH1* pfMet_WM_fail = (TH1*) inputFile->FindObjectAny("metPF_WM_fail");
  TH1* pfMet_ZM_fail = (TH1*) inputFile->FindObjectAny("metPF_ZM_fail");

  pfMet_SR_fail->Scale(1./pfMet_SR_fail->Integral());
  pfMet_WM_fail->Scale(1./pfMet_WM_fail->Integral());
  pfMet_ZM_fail->Scale(1./pfMet_ZM_fail->Integral());

  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();

  TH1F* frame = (TH1F*) caloMet_SR_pas->Clone("frame");
  frame->Reset();
  frame->GetXaxis()->SetTitle("Calo-MET [GeV]");
  frame->GetYaxis()->SetTitle("a.u.");
  frame->Draw();				
  CMS_lumi(canvas,"35.9");

  caloMet_SR_pas->SetLineColor(kBlack);
  caloMet_SR_pas->SetLineWidth(2);
  caloMet_SR_pas->Draw("hist same");

  caloMet_SR_fail->SetLineColor(kRed);
  caloMet_SR_fail->SetLineWidth(2);
  caloMet_SR_fail->Draw("hist same");

  float minim = 0;
  if(min(caloMet_SR_pas->GetMinimum(),caloMet_SR_fail->GetMinimum()) == 0) minim = 0.000001;
  else minim = min(caloMet_SR_pas->GetMinimum(),caloMet_SR_fail->GetMinimum())*0.1;

  frame->GetYaxis()->SetRangeUser(minim, max(caloMet_SR_pas->GetMaximum(),caloMet_SR_fail->GetMaximum())*100);

  TLegend leg (0.7,0.7,0.9,0.9);
  leg.SetBorderSize(0);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.AddEntry(caloMet_SR_pas,"Zvv SR pass","L");
  leg.AddEntry(caloMet_SR_fail,"Zvv SR fail","L");
  leg.Draw("same");

  canvas->SetLogy();

  canvas->SaveAs((outputDIR+"/calomet_zvv_SR.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/calomet_zvv_SR.pdf").c_str(),"pdf");

  frame->Draw();
  CMS_lumi(canvas,"35.9");

  caloMet_WM_pas->SetLineColor(kBlack);
  caloMet_WM_pas->SetLineWidth(2);
  caloMet_WM_pas->Draw("hist same");

  caloMet_WM_fail->SetLineColor(kRed);
  caloMet_WM_fail->SetLineWidth(2);
  caloMet_WM_fail->Draw("hist same");

  minim = 0;
  if(min(caloMet_WM_pas->GetMinimum(),caloMet_WM_fail->GetMinimum()) == 0) minim = 0.000001;
  else minim = min(caloMet_WM_pas->GetMinimum(),caloMet_WM_fail->GetMinimum())*0.01;

  frame->GetYaxis()->SetRangeUser(minim, max(caloMet_WM_pas->GetMaximum(),caloMet_WM_fail->GetMaximum())*100);

  leg.Clear();
  leg.AddEntry(caloMet_WM_pas,"Wmn SR pass","L");
  leg.AddEntry(caloMet_WM_fail,"Wmn SR fail","L");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/calomet_wjet_WM.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/calomet_wjet_WM.pdf").c_str(),"pdf");

  frame->Draw();
  CMS_lumi(canvas,"35.9");

  caloMet_ZM_pas->SetLineColor(kBlack);
  caloMet_ZM_pas->SetLineWidth(2);
  caloMet_ZM_pas->Draw("hist same");

  caloMet_ZM_fail->SetLineColor(kRed);
  caloMet_ZM_fail->SetLineWidth(2);
  caloMet_ZM_fail->Draw("hist same");

  minim = 0;
  if(min(caloMet_ZM_pas->GetMinimum(),caloMet_ZM_fail->GetMinimum()) == 0) minim = 0.000001;
  else minim = min(caloMet_ZM_pas->GetMinimum(),caloMet_ZM_fail->GetMinimum())*0.01;

  frame->GetYaxis()->SetRangeUser(minim, max(caloMet_ZM_pas->GetMaximum(),caloMet_ZM_fail->GetMaximum())*100);

  leg.Clear();
  leg.AddEntry(caloMet_ZM_pas,"Zmm SR pass","L");
  leg.AddEntry(caloMet_ZM_fail,"Zmm SR fail","L");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/calomet_zjet_ZM.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/calomet_zjet_ZM.pdf").c_str(),"pdf");


  //////////////////////////////////////
  frame->GetXaxis()->SetTitle("PF-MET [GeV]");
  frame->Draw();				
  CMS_lumi(canvas,"35.9");

  pfMet_SR_pas->SetLineColor(kBlack);
  pfMet_SR_pas->SetLineWidth(2);
  pfMet_SR_pas->Draw("hist same");

  pfMet_SR_fail->SetLineColor(kRed);
  pfMet_SR_fail->SetLineWidth(2);
  pfMet_SR_fail->Draw("hist same");

  minim = 0;
  if(min(pfMet_SR_pas->GetMinimum(),pfMet_SR_fail->GetMinimum()) == 0) minim = 0.000001;
  else minim = min(pfMet_SR_pas->GetMinimum(),pfMet_SR_fail->GetMinimum())*0.01;

  frame->GetYaxis()->SetRangeUser(minim, max(pfMet_SR_pas->GetMaximum(),pfMet_SR_fail->GetMaximum())*100);

  leg.Clear();
  leg.AddEntry(pfMet_SR_pas,"Zvv SR pass","L");
  leg.AddEntry(pfMet_SR_fail,"Zvv SR fail","L");
  leg.Draw("same");

  canvas->SetLogy();

  canvas->SaveAs((outputDIR+"/pfmet_zvv_SR.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/pfmet_zvv_SR.pdf").c_str(),"pdf");

  frame->Draw();
  CMS_lumi(canvas,"35.9");

  pfMet_WM_pas->SetLineColor(kBlack);
  pfMet_WM_pas->SetLineWidth(2);
  pfMet_WM_pas->Draw("hist same");

  pfMet_WM_fail->SetLineColor(kRed);
  pfMet_WM_fail->SetLineWidth(2);
  pfMet_WM_fail->Draw("hist same");

  minim = 0;
  if(min(pfMet_WM_pas->GetMinimum(),pfMet_WM_fail->GetMinimum()) == 0) minim = 0.000001;
  else minim = min(pfMet_WM_pas->GetMinimum(),pfMet_WM_fail->GetMinimum())*0.01;

  frame->GetYaxis()->SetRangeUser(minim, max(pfMet_WM_pas->GetMaximum(),pfMet_WM_fail->GetMaximum())*100);

  leg.Clear();
  leg.AddEntry(pfMet_WM_pas,"Wmn SR pass","L");
  leg.AddEntry(pfMet_WM_fail,"Wmn SR fail","L");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/pfmet_wjet_WM.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/pfmet_wjet_WM.pdf").c_str(),"pdf");

  frame->Draw();
  CMS_lumi(canvas,"35.9");

  pfMet_ZM_pas->SetLineColor(kBlack);
  pfMet_ZM_pas->SetLineWidth(2);
  pfMet_ZM_pas->Draw("hist same");

  pfMet_ZM_fail->SetLineColor(kRed);
  pfMet_ZM_fail->SetLineWidth(2);
  pfMet_ZM_fail->Draw("hist same");

  minim = 0;
  if(min(pfMet_ZM_pas->GetMinimum(),pfMet_ZM_fail->GetMinimum()) == 0) minim = 0.000001;
  else minim = min(pfMet_ZM_pas->GetMinimum(),pfMet_ZM_fail->GetMinimum())*0.01;

  frame->GetYaxis()->SetRangeUser(minim, max(pfMet_ZM_pas->GetMaximum(),pfMet_ZM_fail->GetMaximum())*100);

  leg.Clear();
  leg.AddEntry(pfMet_ZM_pas,"Zmm SR pass","L");
  leg.AddEntry(pfMet_ZM_fail,"Zmm SR fail","L");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/pfmet_zjet_ZM.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/pfmet_zjet_ZM.pdf").c_str(),"pdf");

}
