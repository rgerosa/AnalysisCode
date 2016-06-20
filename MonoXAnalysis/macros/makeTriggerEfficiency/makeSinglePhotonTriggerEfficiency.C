#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

void makeSinglePhotonTriggerEfficiency(string inputDIR, string ouputDIR, float lumi = 0.86, bool useJetHT = false, int runCut = 999999999) {

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+ouputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  TF1 *fitfunc = new TF1("fitfunc", ErfCB, 150, 300, 5);
  fitfunc->SetParameters(165., 5., 5., 4., 1.);
  vector<float> bins = {160,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,180,190,200,210,220,235,250};
  if(useJetHT)
    bins = {150,160,165,170,175,180,190,200,250,300,400,500,600,700,800,1000,1200,1400};  
  fitfunc->SetRange(bins.front(),bins.back());
  
  TChain* tree = new TChain("tree/tree");
  // should use the wmnu events triggered by single muon
  tree->Add((inputDIR+"/*root").c_str());

  TH1F* hnum   = new TH1F("hnum", "", bins.size()-1, &bins[0]);
  TH1F* hnum_2 = new TH1F("hnum_2", "", bins.size()-1, &bins[0]);
  TH1F* hden   = new TH1F("hden", "", bins.size()-1, &bins[0]);
  TH1F* hden_2 = new TH1F("hden_2", "", bins.size()-1, &bins[0]);
  hnum->Sumw2();
  hnum_2->Sumw2();
  hden->Sumw2();
  hden_2->Sumw2();

  if(not useJetHT){
    // define numerator as event with a medium photon + trigger requirement
    tree->Draw("phpt>>hnum",Form("(hltphoton50 || hltphoton75 || hltphoton90 || hltphoton120) && (hltphoton175 || hltphoton165) && phidm == 1 && abs(pheta) < 1.4442 && run <= %d",runCut));    
    // define denominator as an event with a tight muon passing single muon trigger
    tree->Draw("phpt>>hden",Form("phidm == 1 && abs(pheta) < 1.4442 && (hltphoton50 || hltphoton75 || hltphoton90 || hltphoton120) && run <= %d",runCut));    
  }
  else{
    tree->Draw("phpt>>hnum",Form("phidm == 1 && abs(pheta) < 1.4442 && (hltPFHT800) && (hltphoton175 || hltphoton165) && run <= %d",runCut));
    tree->Draw("phpt>>hden",Form("phidm == 1 && abs(pheta) < 1.4442 && (hltPFHT800) && run <= %d",runCut));

    tree->Draw("phpt>>hnum_2",Form("phidm == 1 && abs(pheta) < 1.4442 && (hltPFHT400 || hltPFHT125 || hltPFHT200 || hltPFHT250 || hltPFHT300 || hltPFHT350 || hltPFHT475 || hltPFHT600 || hltPFHT650) && (hltphoton175 || hltphoton165 || hltPFHT800) && run <= %d",runCut));
    tree->Draw("phpt>>hden_2",Form("phidm == 1 && abs(pheta) < 1.4442 && (hltPFHT400 || hltPFHT125 || hltPFHT200 || hltPFHT250 || hltPFHT300 || hltPFHT350 || hltPFHT475 || hltPFHT600 || hltPFHT650) && run <= %d",runCut));
  }

  TEfficiency* eff = new TEfficiency(*hnum,*hden);  
  eff->Fit(fitfunc);
  eff->SetMarkerColor(kBlack);
  eff->SetLineColor(kBlack);
  eff->SetMarkerStyle(20);
  eff->SetMarkerSize(1);
  fitfunc->SetLineColor(kBlue);
  fitfunc->SetLineWidth(2);
  
  TH1* frame = canvas->DrawFrame(bins.front(), 0.02, bins.back(), 1.1, "");
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

  
  TEfficiency* eff_2 = NULL;
  TF1* fitfunc_2 = NULL;

  if(useJetHT){

    eff_2 = new TEfficiency(*hnum_2,*hden_2);  
    fitfunc_2 = (TF1*) fitfunc->Clone("fitfunc_2");
    eff_2->Fit(fitfunc_2);
    eff->SetMarkerColor(kBlue);
    eff_2->SetMarkerColor(kRed);
    eff->SetLineColor(kBlue);
    eff_2->SetLineColor(kRed);
    eff_2->SetMarkerStyle(20);
    eff_2->SetMarkerSize(1);
    fitfunc_2->SetLineColor(kRed);
    fitfunc_2->SetLineWidth(2);
    
    eff_2->Draw("E1PSAME");
    fitfunc_2->Draw("SAME");
    canvas->RedrawAxis();
    
    TLegend* leg = new TLegend(0.65,0.3,0.9,0.55);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(eff,"HLT Photon 165/175","PLE");
    leg->AddEntry(eff_2,"HLT Photon 165/175 or HLT Photon+PFHT ","PLE");
    leg->Draw("same");
  }
 
  canvas->SaveAs((ouputDIR+"/photonTriggerEff.png").c_str(),"png");
  canvas->SaveAs((ouputDIR+"/photonTriggerEff.pdf").c_str(),"pdf");
  TFile* outputFile = new TFile((ouputDIR+"/photonTriggerEfficiency.root").c_str(),"RECREATE");
  outputFile->cd();
  eff->Write("efficiency_photon");
  if(eff_2 and eff_2 != NULL)
    eff_2->Write("efficiency_photon_pfht");
  fitfunc->Write("efficiency_func");
  if(fitfunc_2 and fitfunc_2 != NULL)
    fitfunc_2->Write("efficiency_func_pfht");

}

