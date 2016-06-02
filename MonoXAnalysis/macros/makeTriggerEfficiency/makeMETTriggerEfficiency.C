#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

void makeMETTriggerEfficiency(string inputDIR, string ouputDIR, float lumi = 0.59) {

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+ouputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  TF1 *fitfunc = new TF1("fitfunc", ErfCB, 50, 1000, 5);
  fitfunc->SetParameters(117., 25., 30., 4., 1.);
  vector<float> bins = {50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100., 105., 110., 115., 120., 130., 140., 150., 160., 180., 200., 225, 250., 275., 300., 350., 400., 500., 600., 1000.};
  
  TChain* tree = new TChain("tree/tree");
  // should use the wmnu events triggered by single muon
  tree->Add((inputDIR+"/*root").c_str());
  
  TH1F* hnum = new TH1F("hnum", "", bins.size()-1, &bins[0]);
  TH1F* hden = new TH1F("hden", "", bins.size()-1, &bins[0]);
  hnum->Sumw2();
  hden->Sumw2();

  // define numerator as event with tight muon + trigger requirement
  tree->Draw("t1mumet>>hnum","(hltmet90==1 || hltmetwithmu170==1) && hltsinglemu==1 && mu1pt>20 && mu1id==1 && centraljetpt[0]>100. && centraljetCHfrac[0]>0.1");    
  // define denominator as an event with a tight muon passing single muon trigger
  tree->Draw("t1mumet>>hden","hltsinglemu==1 && mu1pt>20 && mu1id==1 && centraljetpt[0]>100. && centraljetCHfrac[0]>0.1");    

  TEfficiency* eff = new TEfficiency(*hnum,*hden);
  eff->Fit(fitfunc);
  eff->SetMarkerColor(kBlack);
  eff->SetLineColor(kBlack);
  eff->SetMarkerStyle(20);
  eff->SetMarkerSize(1);
  fitfunc->SetLineColor(kBlue);
  fitfunc->SetLineWidth(2);
  
  TH1* frame = canvas->DrawFrame(bins.front(), 0., bins.back(), 1.1, "");
  frame->GetXaxis()->SetTitle("E_{T#mu}^{miss} [GeV]    ");
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

  canvas->SaveAs((ouputDIR+"/metTriggerEff.png").c_str(),"png");
  canvas->SaveAs((ouputDIR+"/metTriggerEff.pdf").c_str(),"pdf");
  TFile* outputFile = new TFile((ouputDIR+"/metTriggerEfficiency.root").c_str(),"RECREATE");
  outputFile->cd();
  eff->Write("efficiency");
  fitfunc->Write("efficiency_func");
}

