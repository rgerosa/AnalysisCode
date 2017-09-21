#include "../CMS_lumi.h"

void plotHisto(TCanvas* canvas, TH1* nominal, TH1* up, TH1* dw, string outputDIR, string postfix){

  if(TString(postfix).Contains("jes")){
    for(int iBin = 1; iBin < nominal->GetNbinsX()+1; iBin ++){
      if(up->GetBinContent(iBin)/nominal->GetBinContent(iBin) > 1.25) up->SetBinContent(iBin,nominal->GetBinContent(iBin)*1.25); 
      if(dw->GetBinContent(iBin)/nominal->GetBinContent(iBin) > 1.25) dw->SetBinContent(iBin,nominal->GetBinContent(iBin)*1.25); 
      if(up->GetBinContent(iBin)/nominal->GetBinContent(iBin) < 0.8) up->SetBinContent(iBin,nominal->GetBinContent(iBin)*0.8); 
      if(dw->GetBinContent(iBin)/nominal->GetBinContent(iBin) < 0.8) dw->SetBinContent(iBin,nominal->GetBinContent(iBin)*0.8); 
    }
  }
  else if(TString(postfix).Contains("jer")){
    for(int iBin = 1; iBin < nominal->GetNbinsX()+1; iBin ++){
      if(up->GetBinContent(iBin)/nominal->GetBinContent(iBin) < 1) up->SetBinContent(iBin,nominal->GetBinContent(iBin)*1.02); 
      if(dw->GetBinContent(iBin)/nominal->GetBinContent(iBin) > 1) dw->SetBinContent(iBin,nominal->GetBinContent(iBin)*1.025); 
    }
  }
  else if(TString(postfix).Contains("unc")){
    for(int iBin = 1; iBin < nominal->GetNbinsX()+1; iBin ++){
      if(up->GetBinContent(iBin)/nominal->GetBinContent(iBin) < 1) up->SetBinContent(iBin,nominal->GetBinContent(iBin)*1.007); 
      if(dw->GetBinContent(iBin)/nominal->GetBinContent(iBin) > 1) dw->SetBinContent(iBin,nominal->GetBinContent(iBin)*1.005); 
    }
  }

  canvas->cd();
  canvas->SetBottomMargin(0.3);
  nominal->SetLineColor(kBlack);
  nominal->SetMarkerColor(kBlack);
  nominal->SetMarkerStyle(20);
  nominal->SetMarkerSize(1);
  nominal->GetYaxis()->SetRangeUser(dw->GetMinimum()*0.8,up->GetMaximum()*1.2);

  nominal->GetXaxis()->SetTitle("M_{jj} [GeV]");
  nominal->GetXaxis()->SetTitleSize(0);
  nominal->GetXaxis()->SetLabelSize(0);
  nominal->GetYaxis()->SetTitle("Events / GeV");
  nominal->Draw("EP");

  up->SetLineColor(kRed);
  up->SetMarkerColor(kRed);
  up->SetMarkerStyle(20);
  up->SetMarkerSize(1);
  up->Draw("hist same");

  dw->SetLineColor(kRed);
  dw->SetMarkerColor(kRed);
  dw->SetMarkerStyle(20);
  dw->SetMarkerSize(1);
  dw->Draw("hist same");
  
  CMS_lumi(canvas,"35.9");
  
  TLegend leg (0.7,0.7,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(nominal,"nominal","EP");
  leg.AddEntry(up,"up/dw","L");
  leg.Draw("same");

  //// ----                                                                                                                                                                                             
  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);
  pad2->Draw();
  pad2->cd();

  TH1* ratio_up = (TH1*) up->Clone("ratio_up");
  TH1* ratio_dw = (TH1*) dw->Clone("ratio_dw");
  ratio_up->Divide(nominal);
  ratio_dw->Divide(nominal);

  ratio_up->Draw("hist");
  ratio_dw->Draw("hist same");

  ratio_up->GetYaxis()->SetRangeUser(0.7,1.3);
  ratio_up->GetYaxis()->SetNdivisions(504);
  ratio_up->GetYaxis()->SetTitle("Ratio");
  ratio_up->GetXaxis()->SetTitle("M_{jj} [GeV]");

  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  
  if(pad2) delete pad2;

}

void makeJetEnergyUncertainties(string inputFileName, string outputDIR){
  
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  
  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  
  TH1F* ggh_nominal = (TH1F*) inputFile->Get("ggH/ggHhist_125_mjj");
  TH1F* ggh_jesUp = (TH1F*) inputFile->Get("ggH/sysShape/ggHhist_metJetUp_125_mjj");
  TH1F* ggh_jesDw = (TH1F*) inputFile->Get("ggH/sysShape/ggHhist_metJetDw_125_mjj");
  TH1F* ggh_jerUp = (TH1F*) inputFile->Get("ggH/sysShape/ggHhist_metResUp_125_mjj");
  TH1F* ggh_jerDw = (TH1F*) inputFile->Get("ggH/sysShape/ggHhist_metResDw_125_mjj");
  TH1F* ggh_uncUp = (TH1F*) inputFile->Get("ggH/sysShape/ggHhist_metUncUp_125_mjj");
  TH1F* ggh_uncDw = (TH1F*) inputFile->Get("ggH/sysShape/ggHhist_metUncDw_125_mjj");

  TH1F* vbf_nominal = (TH1F*) inputFile->Get("vbfH/vbfHhist_125_mjj");
  TH1F* vbf_jesUp = (TH1F*) inputFile->Get("vbfH/sysShape/vbfHhist_metJetUp_125_mjj");
  TH1F* vbf_jesDw = (TH1F*) inputFile->Get("vbfH/sysShape/vbfHhist_metJetDw_125_mjj");
  TH1F* vbf_jerUp = (TH1F*) inputFile->Get("vbfH/sysShape/vbfHhist_metResUp_125_mjj");
  TH1F* vbf_jerDw = (TH1F*) inputFile->Get("vbfH/sysShape/vbfHhist_metResDw_125_mjj");
  TH1F* vbf_uncUp = (TH1F*) inputFile->Get("vbfH/sysShape/vbfHhist_metUncUp_125_mjj");
  TH1F* vbf_uncDw = (TH1F*) inputFile->Get("vbfH/sysShape/vbfHhist_metUncDw_125_mjj");

  vbf_nominal->Scale(1,"width");
  vbf_jesUp->Scale(1,"width");
  vbf_jesDw->Scale(1,"width");
  vbf_jerUp->Scale(1,"width");
  vbf_jerDw->Scale(1,"width");
  vbf_uncUp->Scale(1,"width");
  vbf_uncDw->Scale(1,"width");

  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();

  plotHisto(canvas,ggh_nominal,ggh_jesUp,ggh_jesDw,outputDIR,"ggh_jes");
  plotHisto(canvas,ggh_nominal,ggh_jerUp,ggh_jerDw,outputDIR,"ggh_jer");
  plotHisto(canvas,ggh_nominal,ggh_uncUp,ggh_uncDw,outputDIR,"ggh_unc");

  plotHisto(canvas,vbf_nominal,vbf_jesUp,vbf_jesDw,outputDIR,"vbf_jes");
  plotHisto(canvas,vbf_nominal,vbf_jerUp,vbf_jerDw,outputDIR,"vbf_jer");
  plotHisto(canvas,vbf_nominal,vbf_uncUp,vbf_uncDw,outputDIR,"vbf_unc");

}
