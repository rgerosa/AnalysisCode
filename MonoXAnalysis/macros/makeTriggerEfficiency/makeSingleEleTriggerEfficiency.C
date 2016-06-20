#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

// calculation from the jetHT dataset
void makeSingleElectronTriggerEfficiency(string inputDIR, string outputDIR, float lumi = 0.86) {

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  TF1 *fitfunc = new TF1("fitfunc", ErfCB, 20, 1200, 5);
  fitfunc->SetParameters(30., 5., 5., 4., 1.);

  vector<float> binsPt  = {20,30,40,50,75,100,125,150,200,250,300,400,500,600,1000};
  vector<float> binsEta = {0,0.5,1.0,1.444,1.56,2.1,2.5};  
  fitfunc->SetRange(binsPt.front(),binsPt.back());
  
  TChain* tree = new TChain("tree/tree");
  // should use the wmnu events triggered by single muon
  tree->Add((inputDIR+"/*root").c_str());

  TH2F* hnum = new TH2F("hnum", "", binsEta.size()-1, &binsEta[0],binsPt.size()-1, &binsPt[0]);
  TH2F* hden = new TH2F("hden", "", binsEta.size()-1, &binsEta[0],binsPt.size()-1, &binsPt[0]);
  hnum->Sumw2();
  hden->Sumw2();

  // define numerator as event with a medium photon + trigger requirement
  //  tree->Draw("el1pt:abs(el1eta) >> hnum",Form("el1id == 1 && nelectrons == 1 && abs(el1eta) < 2.5 && (hltPFHT800 || hltPFHT650) && (hltsingleel) && t1pfmet > 50"));
  tree->Draw("el1pt:abs(el1eta) >> hnum",Form("el1id == 1 && nelectrons == 1 && abs(el1eta) < 2.5 && (hltPFHT800 || hltPFHT650) && (hltsingleel || hltelnoiso) && t1pfmet > 50"));
  tree->Draw("el1pt:abs(el1eta) >> hden",Form("el1id == 1 && nelectrons == 1 && abs(el1eta) < 2.5 && (hltPFHT650 || hltPFHT800) && t1pfmet > 50"));
  
  TH2F* efficiency = (TH2F*) hnum->Clone("efficiency");
  efficiency->Reset();
  efficiency->Divide(hnum,hden,1,1,"B");
  fitfunc->SetLineColor(kBlue);
  fitfunc->SetLineWidth(2);
  
  gPad->SetBottomMargin(0.10);
  gPad->SetRightMargin(0.17);
  efficiency->GetYaxis()->SetTitle("Electron p_{T} [GeV]");
  efficiency->GetXaxis()->SetTitle("Electron #eta");
  efficiency->GetZaxis()->SetTitle("Trigger Efficiency");
  efficiency->GetYaxis()->SetLabelSize(0.8*efficiency->GetYaxis()->GetLabelSize());
  efficiency->GetXaxis()->SetLabelSize(0.8*efficiency->GetXaxis()->GetLabelSize());
  efficiency->GetZaxis()->SetLabelSize(0.8*efficiency->GetXaxis()->GetLabelSize());
  efficiency->GetYaxis()->SetTitleSize(0.8*efficiency->GetYaxis()->GetTitleSize());
  efficiency->GetXaxis()->SetTitleSize(0.8*efficiency->GetXaxis()->GetTitleSize());
  efficiency->GetZaxis()->SetTitleSize(0.8*efficiency->GetXaxis()->GetTitleSize());
  efficiency->GetXaxis()->SetTitleOffset(1.0);  
  canvas->SetTopMargin(0.06);
  canvas->Draw();
  canvas->cd();
  canvas->SetLogy();
  gStyle->SetPaintTextFormat(".2f");
  efficiency->Draw("colztext");
  canvas->RedrawAxis();
  CMS_lumi(canvas,string(Form("%.2f",lumi)),true);
  canvas->SaveAs((outputDIR+"/electronTriggerEff.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/electronTriggerEff.pdf").c_str(),"pdf");
  canvas->SetLogy(0);

  gPad->SetRightMargin(0.06);
  TH1F projection_pt ("projection_pt","",efficiency->GetYaxis()->GetXbins()->GetSize()-1,efficiency->GetYaxis()->GetXbins()->GetArray());
  for(int iBinX = 0; iBinX < efficiency->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < efficiency->GetNbinsY(); iBinY++){
      projection_pt.SetBinContent(iBinY+1,efficiency->GetBinContent(iBinX+1,iBinY+1));
      projection_pt.SetBinError(iBinY+1,efficiency->GetBinError(iBinX+1,iBinY+1));
    }
    projection_pt.GetXaxis()->SetTitle("Electron p_{T} [GeV]");
    projection_pt.GetYaxis()->SetTitle("Trigger Efficiency");
    projection_pt.SetMarkerSize(1);
    projection_pt.SetMarkerStyle(20);
    projection_pt.SetMarkerColor(kBlack);
    projection_pt.SetLineColor(kBlack);
    projection_pt.Draw("E1P");
    projection_pt.Fit(fitfunc);
    CMS_lumi(canvas,string(Form("%.2f",lumi)),true);
    canvas->SaveAs((outputDIR+"/"+string(efficiency->GetName())+string(Form("_projection_pt_eta_%.1f_%.1f.png",efficiency->GetXaxis()->GetBinLowEdge(iBinX+1),
                                                                       efficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))).c_str(),"png");
    canvas->SaveAs((outputDIR+"/"+string(efficiency->GetName())+string(Form("_projection_pt_eta_%.1f_%.1f.pdf",efficiency->GetXaxis()->GetBinLowEdge(iBinX+1),
                                                                       efficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))).c_str(),"pdf");
  }
  
  TFile* outputFile = new TFile((outputDIR+"/triggerEfficiency_DATA_SingleElectron.root").c_str(),"RECREATE");
  outputFile->cd();
  efficiency->Write("efficiency");
}

