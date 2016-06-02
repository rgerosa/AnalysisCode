#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

void makeSinglePhotonTriggerEfficiency(string inputDIR, string ouputDIR, float lumi = 0.59) {

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+ouputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  TF1 *fitfunc = new TF1("fitfunc", ErfCB, 150, 300, 5);
  fitfunc->SetParameters(160., 5., 5., 4., 1.);
  vector<float> bins = {160,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,190,200};
  
  TChain* tree = new TChain("tree/tree");
  // should use the wmnu events triggered by single muon
  tree->Add((inputDIR+"/*root").c_str());
  
  TH1F* hnum = new TH1F("hnum", "", bins.size()-1, &bins[0]);
  TH1F* hden = new TH1F("hden", "", bins.size()-1, &bins[0]);
  hnum->Sumw2();
  hden->Sumw2();

  // define numerator as event with a medium photon + trigger requirement
  tree->Draw("phpt>>hnum","(hltphoton165 || hltphoton175) && phidm == 1 && abs(pheta) < 1.4442");    
  // define denominator as an event with a tight muon passing single muon trigger
  tree->Draw("phpt>>hden","phidm == 1 && abs(pheta) < 1.4442");    
  
  TEfficiency* eff = new TEfficiency(*hnum,*hden);
  eff->Fit(fitfunc);
  eff->SetMarkerColor(kBlack);
  eff->SetLineColor(kBlack);
  eff->SetMarkerStyle(20);
  eff->SetMarkerSize(1);
  fitfunc->SetLineColor(kBlue);
  fitfunc->SetLineWidth(2);
  
  TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 1.1, "");
  frame->GetXaxis()->SetTitle("Photon p_{T} [GeV]");
  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetYaxis()->SetLabelSize(0.8*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetTitleSize(0.8*frame->GetYaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleOffset(1.0);
  
  canvas->SetRightMargin(0.075);
  canvas->SetTopMargin(0.06);
  canvas->Draw();
  canvas->cd();
  frame->Draw();
  eff->Draw("E1PSAME");
  fitfunc->Draw("SAME");
  canvas->RedrawAxis();
  CMS_lumi(canvas,string(Form("%.2f",lumi)),true);

  canvas->SaveAs((ouputDIR+"/photonTriggerEff.png").c_str(),"png");
  canvas->SaveAs((ouputDIR+"/photonTriggerEff.pdf").c_str(),"pdf");
  TFile* outputFile = new TFile((ouputDIR+"/photonTriggerEfficiency.root").c_str(),"RECREATE");
  outputFile->cd();
  eff->Write("efficiency");
  fitfunc->Write("efficiency_func");
}

