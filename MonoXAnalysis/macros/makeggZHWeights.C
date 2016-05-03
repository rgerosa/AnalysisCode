#include "CMS_lumi.h"

void makeggZHWeights(){

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.2f");

  TFile* qqZH_MM = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-19-2-2016/ggZH/tree_qqZH_MM.root");
  TFile* qqZH_EE = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-19-2-2016/ggZH/tree_qqZH_EE.root");
  TFile* ggZH = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-19-2-2016/ggZH/tree_ggZH.root");

  TTree* qqZHTree_MM = (TTree*) qqZH_MM->Get("tree/tree");
  TTree* qqZHTree_EE = (TTree*) qqZH_EE->Get("tree/tree");
  TTree* ggZHTree = (TTree*) ggZH->Get("tree/tree");

  TH1F* qqZH_hpt_MM = new TH1F("qqZH_hpt_MM","",100,0,500);
  TH1F* qqZH_hpt_EE = new TH1F("qqZH_hpt_EE","",100,0,500);
  TH1F* qqZH_hpt    = NULL;
  TH1F* ggZH_hpt = new TH1F("ggZH_hpt","",100,0,500);
  qqZH_hpt_MM->Sumw2();
  qqZH_hpt_EE->Sumw2();
  ggZH_hpt->Sumw2();
  
  TH1F* qqZH_zpt_MM = new TH1F("qqZH_zpt_MM","",100,0,500);
  TH1F* qqZH_zpt_EE = new TH1F("qqZH_zpt_EE","",100,0,500);
  TH1F* qqZH_zpt = NULL;
  TH1F* ggZH_zpt = new TH1F("ggZH_zpt","",100,0,500);
  qqZH_zpt_MM->Sumw2();
  qqZH_zpt_EE->Sumw2();
  ggZH_zpt->Sumw2();
  
  qqZHTree_MM->Draw("wzpt >> qqZH_zpt_MM","","goff");
  qqZHTree_EE->Draw("wzpt >> qqZH_zpt_EE","","goff");
  ggZHTree->Draw("wzpt >> ggZH_zpt","","goff");

  qqZHTree_MM->Draw("dmpt >> qqZH_hpt_MM","","goff");
  qqZHTree_EE->Draw("dmpt >> qqZH_hpt_EE","","goff");
  ggZHTree->Draw("dmpt >> ggZH_hpt","","goff");

  qqZH_hpt = (TH1F*) qqZH_hpt_MM->Clone("qqZH_hpt");
  qqZH_hpt->Add(qqZH_hpt_EE);

  qqZH_zpt = (TH1F*) qqZH_zpt_MM->Clone("qqZH_zpt");
  qqZH_zpt->Add(qqZH_zpt_EE);

  vector<float> higgsPT  = {50,100,150,200,250,350.,550};
  vector<float> bosonPT  = {50,100,150,200,275,550};

  TH2F* qqZH_MM_2D = new TH2F("qqZH_MM_2D","",higgsPT.size()-1,&higgsPT[0],bosonPT.size()-1,&bosonPT[0]);
  TH2F* qqZH_EE_2D = new TH2F("qqZH_EE_2D","",higgsPT.size()-1,&higgsPT[0],bosonPT.size()-1,&bosonPT[0]);
  TH2F* qqZH_2D = NULL;
  TH2F* ggZH_2D = new TH2F("ggZH_2D","",higgsPT.size()-1,&higgsPT[0],bosonPT.size()-1,&bosonPT[0]);
  qqZH_MM_2D->Sumw2();
  qqZH_EE_2D->Sumw2();
  ggZH_2D->Sumw2();

  qqZHTree_MM->Draw("dmpt:wzpt >> qqZH_MM_2D","","goff");
  qqZHTree_EE->Draw("dmpt:wzpt >> qqZH_EE_2D","","goff");
  ggZHTree->Draw("dmpt:wzpt >> ggZH_2D","","goff");

  qqZH_2D = (TH2F*) qqZH_MM_2D->Clone("qqZH_2D");
  qqZH_2D->Add(qqZH_EE_2D);
  
  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->SetLeftMargin(0.15);
  
  canvas->SetTickx();
  canvas->SetTicky();

  qqZH_hpt->Scale(1./qqZH_hpt->Integral());
  ggZH_hpt->Scale(1./ggZH_hpt->Integral());

  qqZH_hpt->GetXaxis()->SetTitle("Higgs p_{T}");
  qqZH_hpt->GetYaxis()->SetTitle("a.u.");
  qqZH_hpt->SetLineColor(kRed);
  qqZH_hpt->SetLineWidth(2);
  ggZH_hpt->SetLineColor(kBlue);
  ggZH_hpt->SetLineWidth(2);

  qqZH_hpt->Draw("hist");
  ggZH_hpt->Draw("hist same");

  CMS_lumi(canvas,"2.30");

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.8);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(qqZH_hpt,"qqZH","l");
  leg->AddEntry(ggZH_hpt,"ggZH","l");
  leg->Draw("same");

  TLine* line = new TLine(200.,0.,200.,qqZH_hpt->GetMaximum());
  line->SetLineColor(kBlack);
  line->SetLineWidth(2);
  line->Draw("same");
  canvas->SetLogy();
  canvas->SaveAs("higgs_pt.pdf","pdf");
  canvas->SaveAs("higgs_pt.png","png");

  qqZH_zpt->Scale(1./qqZH_zpt->Integral());
  ggZH_zpt->Scale(1./ggZH_zpt->Integral());

  qqZH_zpt->GetXaxis()->SetTitle("Z-boson p_{T}");
  qqZH_zpt->GetYaxis()->SetTitle("a.u.");
  qqZH_zpt->SetLineColor(kRed);
  qqZH_zpt->SetLineWidth(2);
  ggZH_zpt->SetLineColor(kBlue);
  ggZH_zpt->SetLineWidth(2);

  qqZH_zpt->Draw("hist");
  ggZH_zpt->Draw("hist same");
  leg->Draw("same");

  TLine* line2 = new TLine(250.,0.,250.,qqZH_hpt->GetMaximum());
  line2->SetLineColor(kBlack);
  line2->SetLineWidth(2);
  line2->Draw("same");

  TLine* line3 = new TLine(100.,0.,100.,qqZH_hpt->GetMaximum());
  line3->SetLineColor(kBlack);
  line3->SetLineWidth(2);
  line3->Draw("same");

  CMS_lumi(canvas,"2.30");

  canvas->SetLogy();
  canvas->SaveAs("z_pt.pdf","pdf");
  canvas->SaveAs("z_pt.png","png");

  canvas->SetRightMargin(0.18);

  qqZH_2D->Scale(1./qqZH_2D->Integral());
  ggZH_2D->Scale(1./ggZH_2D->Integral());

  qqZH_2D->GetXaxis()->SetTitle("Higgs p_{T}");
  qqZH_2D->GetYaxis()->SetTitle("Z-boson p_{T}");
  qqZH_2D->GetZaxis()->SetTitle("a.u.");
  qqZH_2D->Draw("colz");

  CMS_lumi(canvas,"2.30");
  canvas->SetLogz();
  canvas->SetLogy(0);
  canvas->SaveAs("qqZH_2D.pdf","pdf");
  canvas->SaveAs("qqZH_2D.png","png");
  
  ggZH_2D->GetXaxis()->SetTitle("Higgs p_{T}");
  ggZH_2D->GetYaxis()->SetTitle("Z-boson p_{T}");
  ggZH_2D->GetZaxis()->SetTitle("a.u.");
  ggZH_2D->Draw("colz");

  CMS_lumi(canvas,"2.30");

  canvas->SaveAs("ggZH_2D.pdf","pdf");
  canvas->SaveAs("ggZH_2D.png","png");

  TH2F* weight = (TH2F*) ggZH_2D->Clone("weight");
  weight->Reset();
  for(int iBin = 1; iBin <= ggZH_2D->GetNbinsX(); iBin++){
    for(int jBin = 1; jBin <= ggZH_2D->GetNbinsY(); jBin++){
      if(qqZH_2D->GetBinContent(iBin,jBin) != 0)
	weight->SetBinContent(iBin,jBin,ggZH_2D->GetBinContent(iBin,jBin)/qqZH_2D->GetBinContent(iBin,jBin));
      else
	weight->SetBinContent(iBin,jBin,0);
    } 
  }
  weight->GetZaxis()->SetTitle("weight");
  weight->Draw("colz text2"); 
  CMS_lumi(canvas,"2.30");

  canvas->SetLogz(0);
  canvas->SaveAs("weight_2D.pdf","pdf");
  canvas->SaveAs("weight_2D.png","png");

  TGraph2D* graph = new TGraph2D(weight);
  graph->SetNpx(100);
  graph->SetNpy(100);
  for(int iBinX = 0; iBinX <= weight->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY <= weight->GetNbinsY(); iBinY++){
      graph->Interpolate(weight->GetXaxis()->GetBinCenter(iBinX),weight->GetYaxis()->GetBinCenter(iBinY));
    }
  }

  weight = (TH2F*) graph->GetHistogram();  
  weight->GetXaxis()->SetTitle("Higgs p_{T}");
  weight->GetYaxis()->SetTitle("Z-boson p_{T}");
  weight->GetZaxis()->SetTitle("weight");
  weight->Draw("colz"); 

  CMS_lumi(canvas,"2.30");
  canvas->SaveAs("weight_2D_smooth.pdf","pdf");
  canvas->SaveAs("weight_2D_smooth.png","png");

  TFile* output = new TFile("ggZH_weight.root","RECREATE");
  output->cd();
  weight->Write();
  output->Close();

}
