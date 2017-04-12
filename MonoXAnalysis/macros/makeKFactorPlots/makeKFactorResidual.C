#include "../CMS_lumi.h"

void makeKFactorResidual(string inputFile,  string outputDIR){

  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  TFile* input = TFile::Open(inputFile.c_str());
  
  TH1* kfactorZ = (TH1*) ((TH1*) input->Get("ZJets_012j_NLO/nominal"))->Clone("kfactorZ");
  kfactorZ->Divide((TH1*) input->Get("ZJets_LO/inv_pt"));
  TH1* kfactorW = (TH1*) ((TH1*) input->Get("WJets_012j_NLO/nominal"))->Clone("kfactorW");
  kfactorW->Divide((TH1*) input->Get("WJets_LO/inv_pt"));  
  TH1* kfactorG = (TH1*) ((TH1*) input->Get("GJets_1j_NLO/nominal_G"))->Clone("kfactorG");
  kfactorG->Divide((TH1*) input->Get("GJets_LO/inv_pt_G"));

  TH1* kfactor_mean = (TH1*) kfactorZ->Clone("kfactor");  
  kfactorZ->Rebin(2);
  kfactorW->Rebin(2);
  kfactorG->Rebin(2);
  kfactor_mean->Rebin(2);

  kfactorZ->Add(kfactor_mean,-1);
  kfactorW->Add(kfactor_mean,-1);
  kfactorG->Add(kfactor_mean,-1);

  kfactorZ->Divide(kfactor_mean);
  kfactorW->Divide(kfactor_mean);
  kfactorG->Divide(kfactor_mean);

  TCanvas* canvas = new TCanvas("canvas","canvas",600,600);
  kfactorZ->GetXaxis()->SetTitle("boson p_{T} [GeV]");
  kfactorZ->GetYaxis()->SetTitle("#Delta K_{NLO-QCD}/K_{mean}");
  kfactorZ->GetYaxis()->SetTitleOffset(1.3);
  kfactorZ->SetLineColor(kBlack);
  kfactorZ->SetLineWidth(2);
  kfactorW->SetLineColor(kBlue);
  kfactorW->SetLineWidth(2);
  kfactorG->SetLineColor(kRed);
  kfactorG->SetLineWidth(2);

  kfactorZ->GetXaxis()->SetRangeUser(250,1250);
  kfactorZ->GetYaxis()->SetRangeUser(min(kfactorZ->GetMinimum(),min(kfactorW->GetMinimum(),kfactorG->GetMinimum()))*1.5,
				     max(kfactorZ->GetMaximum(),max(kfactorW->GetMaximum(),kfactorG->GetMaximum()))*2);
  
  kfactorZ->Draw("hist");
  kfactorW->Draw("hist same");
  kfactorG->Draw("hist same");
  CMS_lumi(canvas,"");

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetBorderSize(0);
  leg->AddEntry(kfactorZ,"#DeltaK Z-jets","L");
  leg->AddEntry(kfactorW,"#DeltaK W-jets","L");
  leg->AddEntry(kfactorG,"#DeltaK #gamma-jets","L");
  leg->Draw("same");

  canvas->SaveAs((outputDIR+"/kfactor_difference.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/kfactor_difference.png").c_str(),"png");
}
