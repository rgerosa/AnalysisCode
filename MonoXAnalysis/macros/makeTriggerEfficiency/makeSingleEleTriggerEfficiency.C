#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

// calculation from the jetHT dataset
void makeSingleElectronTriggerEfficiency(string inputDIR, string outputDIR, float lumi = 0.86, bool doFit = false) {

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  TF1 *fitfunc = new TF1("fitfunc", ErfCB, 20, 1200, 5);
  fitfunc->SetParameters(30., 5., 5., 4., 1.);

  vector<float> binsPt  = {100,125,150,200,250,300,350,400,450,500,550,600,700,800,1200};
  vector<float> binsEta = {0,1.5,2.5};  
  fitfunc->SetRange(binsPt.front(),binsPt.back());
  
  TChain* tree = new TChain("tree/tree");
  // should use the wmnu events triggered by single muon
  tree->Add((inputDIR+"/*root").c_str());

  TH2F* hnum = new TH2F("hnum", "", binsEta.size()-1, &binsEta[0],binsPt.size()-1, &binsPt[0]);
  TH2F* hden = new TH2F("hden", "", binsEta.size()-1, &binsEta[0],binsPt.size()-1, &binsPt[0]);
  hnum->Sumw2();
  hden->Sumw2();

  TH2F* hnum_recover = new TH2F("hnum_recover", "", binsEta.size()-1, &binsEta[0],binsPt.size()-1, &binsPt[0]);
  TH2F* hden_recover = new TH2F("hden_recover", "", binsEta.size()-1, &binsEta[0],binsPt.size()-1, &binsPt[0]);
  hnum_recover->Sumw2();
  hden_recover->Sumw2();

  // define numerator as event with a medium photon + trigger requirement
  tree->Draw("el1pt:abs(el1eta) >> hnum",Form("el1id == 1 && nelectrons == 1 && abs(el1eta) < 2.5 && (hltPFHT400 || hltPFHT650 || hltPFHT800) && (hltsingleel || hltelnoiso) && t1pfmet > 50"));
  tree->Draw("el1pt:abs(el1eta) >> hden",Form("el1id == 1 && nelectrons == 1 && abs(el1eta) < 2.5 && (hltPFHT400 || hltPFHT650 || hltPFHT800) && t1pfmet > 50"));
  tree->Draw("el1pt:abs(el1eta) >> hnum_recover",Form("el1id == 1 && nelectrons == 1 && abs(el1eta) < 2.5 && (hltPFHT400 || hltPFHT650) && (hltsingleel || hltelnoiso || hltPFHT800) && t1pfmet > 50"));
  tree->Draw("el1pt:abs(el1eta) >> hden_recover",Form("el1id == 1 && nelectrons == 1 && abs(el1eta) < 2.5 && (hltPFHT400 || hltPFHT650) && t1pfmet > 50"));
  
  TEfficiency* efficiency = new TEfficiency(*hnum,*hden);
  TEfficiency* efficiency_recover = new TEfficiency(*hnum_recover,*hden_recover);
  TH2* histoEff = efficiency->CreateHistogram();
  TH2* histoEff_recover = efficiency_recover->CreateHistogram();
  // in order to plot the 2D histo
  fitfunc->SetLineColor(kBlue);
  fitfunc->SetLineWidth(2);
  
  gPad->SetBottomMargin(0.10);
  gPad->SetRightMargin(0.17);
  TH1* frame = canvas->DrawFrame(binsEta.front(),binsPt.front(), binsEta.back(), binsPt.back(), "");  
  frame->GetYaxis()->SetTitle("Electron p_{T} [GeV]");
  frame->GetXaxis()->SetTitle("Electron #eta");
  frame->GetZaxis()->SetTitle("Trigger Efficiency");
  frame->GetYaxis()->SetLabelSize(0.8*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame->GetZaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetTitleSize(0.8*frame->GetYaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame->GetZaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleOffset(1.0);  
  canvas->SetTopMargin(0.06);
  canvas->Draw();
  canvas->cd();
  canvas->SetLogy();
  gStyle->SetPaintTextFormat(".2f");
  efficiency->Draw("colztext same");
  canvas->RedrawAxis();
  CMS_lumi(canvas,string(Form("%.1f",lumi)),true);
  canvas->SaveAs((outputDIR+"/electronTriggerEff.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/electronTriggerEff.pdf").c_str(),"pdf");
  canvas->SetLogy(0);

  canvas->Draw();
  canvas->cd();
  canvas->SetLogy();
  gStyle->SetPaintTextFormat(".2f");
  efficiency_recover->Draw("colztext same");
  canvas->RedrawAxis();
  CMS_lumi(canvas,string(Form("%.1f",lumi)),true);
  canvas->SaveAs((outputDIR+"/electronTriggerEff_recover.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/electronTriggerEff_recover.pdf").c_str(),"pdf");
  canvas->SetLogy(0);

  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 600, 625);
  canvas2->cd();
  gPad->SetRightMargin(0.06);
  for(int iBinX = 0; iBinX < histoEff->GetNbinsX(); iBinX++){
    TGraphAsymmErrors* projection_pt = new TGraphAsymmErrors();
    TGraphAsymmErrors* projection_pt_recover = new TGraphAsymmErrors();
    for(int iBinY = 0; iBinY < histoEff->GetNbinsY(); iBinY++){
      int globalBin = histoEff->GetBin(iBinX+1,iBinY+1);
      projection_pt->SetPoint(iBinY+1,histoEff->GetYaxis()->GetBinCenter(iBinY+1),efficiency->GetEfficiency(globalBin));
      projection_pt->SetPointError(iBinY+1,histoEff->GetYaxis()->GetBinWidth(iBinY+1)/2,histoEff->GetYaxis()->GetBinWidth(iBinY+1)/2,efficiency->GetEfficiencyErrorLow(globalBin),efficiency->GetEfficiencyErrorUp(globalBin));
    }

    for(int iBinY = 0; iBinY < histoEff_recover->GetNbinsY(); iBinY++){
      int globalBin = histoEff_recover->GetBin(iBinX+1,iBinY+1);
      projection_pt_recover->SetPoint(iBinY+1,histoEff_recover->GetYaxis()->GetBinCenter(iBinY+1),efficiency_recover->GetEfficiency(globalBin));
      projection_pt_recover->SetPointError(iBinY+1,histoEff_recover->GetYaxis()->GetBinWidth(iBinY+1)/2,histoEff_recover->GetYaxis()->GetBinWidth(iBinY+1)/2,efficiency_recover->GetEfficiencyErrorLow(globalBin),efficiency_recover->GetEfficiencyErrorUp(globalBin));
    }

    projection_pt->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
    projection_pt->GetXaxis()->SetRangeUser(binsPt.front(),binsPt.back());
    projection_pt->GetYaxis()->SetTitle("Trigger EfficiencyHisto");
    projection_pt->SetMarkerSize(1);
    projection_pt->SetMarkerStyle(20);
    projection_pt->SetMarkerColor(kRed);
    projection_pt->SetLineColor(kRed);
    projection_pt->Draw("AE1P");
    projection_pt_recover->SetMarkerSize(1);
    projection_pt_recover->SetMarkerStyle(20);
    projection_pt_recover->SetMarkerColor(kBlue);
    projection_pt_recover->SetLineColor(kBlue);
    projection_pt_recover->Draw("E1Psame");

    TLegend* leg = new TLegend(0.6,0.25,0.9,0.45);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(projection_pt,"HLT_Ele27 || HLT_Ele_105","PE");
    leg->AddEntry(projection_pt_recover,"HLT_Ele27 || HLT_Ele_105 || HLT_PFHT800","PE");
      
    
    if(doFit)
      projection_pt->Fit(fitfunc);
    CMS_lumi(canvas2,string(Form("%.1f",lumi)),true);
    canvas2->SaveAs((outputDIR+"/"+string(histoEff->GetName())+string(Form("_projection_pt_eta_%.1f_%.1f.png",histoEff->GetXaxis()->GetBinLowEdge(iBinX+1),
									   histoEff->GetXaxis()->GetBinLowEdge(iBinX+2)))).c_str(),"png");
    canvas2->SaveAs((outputDIR+"/"+string(histoEff->GetName())+string(Form("_projection_pt_eta_%.1f_%.1f.pdf",histoEff->GetXaxis()->GetBinLowEdge(iBinX+1),
										 histoEff->GetXaxis()->GetBinLowEdge(iBinX+2)))).c_str(),"pdf");
  }
  
  TFile* outputFile = new TFile((outputDIR+"/triggerEfficiencyHisto_DATA_SingleElectron.root").c_str(),"RECREATE");
  outputFile->cd();
  efficiency->Write("efficiency");
}

