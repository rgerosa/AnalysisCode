#include "../CMS_lumi.h"

void makeWJetsSplitSR(string fileName, string observable, string observableLatex, string outputDIR, bool isEWK){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TFile* input = TFile::Open(fileName.c_str());
  TH1* wjet_total = NULL;
  TH1* wjet_muon  = NULL;
  TH1* wjet_ele   = NULL;
  TH1* wjet_tau   = NULL;

  if(not isEWK){
    wjet_total = (TH1*) input->FindObjectAny(("wjethist_"+observable).c_str());
    wjet_muon  = (TH1*) input->FindObjectAny(("wjethist_mu_"+observable).c_str());
    wjet_ele   = (TH1*) input->FindObjectAny(("wjethist_el_"+observable).c_str());
    wjet_tau   = (TH1*) input->FindObjectAny(("wjethist_ta_"+observable).c_str());
  }
  else{
    wjet_total = (TH1*) input->FindObjectAny(("ewkbkgwhist_"+observable).c_str());
    wjet_muon  = (TH1*) input->FindObjectAny(("ewkwjethist_mu_"+observable).c_str());
    wjet_ele   = (TH1*) input->FindObjectAny(("ewkwjethist_el_"+observable).c_str());
    wjet_tau   = (TH1*) input->FindObjectAny(("ewkwjethist_ta_"+observable).c_str());
  }

  wjet_total->Scale(1.,"width");
  wjet_muon->Scale(1.,"width");
  wjet_ele->Scale(1.,"width");
  wjet_tau->Scale(1.,"width");
  /*
  // as a check
  for(int iBin = 0; iBin < wjet_total->GetNbinsX()+1; iBin++){
    cout<<"iBin "<<iBin+1<<" muon "<<wjet_muon->GetBinContent(iBin+1)<<" ele "<<wjet_ele->GetBinContent(iBin+1)<<" tau "<<wjet_tau->GetBinContent(iBin+1)<<" total "<<
      wjet_muon->GetBinContent(iBin+1)+wjet_ele->GetBinContent(iBin+1)+wjet_tau->GetBinContent(iBin+1)<<" inclusive "<<wjet_total->GetBinContent(iBin+1)<<endl;
  }
  */

  // make the plot in absolute yields
  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  wjet_total->SetLineColor(kBlack);
  wjet_total->SetLineWidth(2);
  wjet_total->GetYaxis()->SetTitle("Events / GeV");
  wjet_total->GetYaxis()->SetTitleOffset(1.25);
  wjet_total->GetYaxis()->SetLabelSize(0.035);
  wjet_total->GetYaxis()->SetTitleSize(0.045);
  wjet_total->GetXaxis()->SetLabelSize(0);
  wjet_total->GetXaxis()->SetTitleSize(0);

  wjet_total->Draw("hist");
  
  wjet_muon->SetLineColor(kBlue);
  wjet_muon->SetLineWidth(2);
  wjet_muon->Draw("hist same");

  wjet_ele->SetLineColor(kRed);
  wjet_ele->SetLineWidth(2);
  wjet_ele->Draw("hist same");

  wjet_tau->SetLineColor(kCyan);
  wjet_tau->SetLineWidth(2);
  wjet_tau->Draw("hist same");

  wjet_total->GetYaxis()->SetRangeUser(min(wjet_total->GetMinimum(),min(wjet_muon->GetMinimum(),min(wjet_ele->GetMinimum(),wjet_tau->GetMinimum())))*0.1,
				       max(wjet_total->GetMaximum(),max(wjet_muon->GetMaximum(),max(wjet_ele->GetMaximum(),wjet_tau->GetMaximum())))*100);
  CMS_lumi(canvas,"35.9");
  canvas->SetLogy();

  TLegend leg (0.6,0.7,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  if(not isEWK){
    leg.AddEntry(wjet_total,"W+jets inclusive","L");
    leg.AddEntry(wjet_muon,"W #rightarrow #mu#nu ","L");
    leg.AddEntry(wjet_ele,"W #rightarrow e#nu ","L");
    leg.AddEntry(wjet_tau,"W #rightarrow #tau#nu ","L");
  }
  else{
    leg.AddEntry(wjet_total,"EWK-W inclusive","L");
    leg.AddEntry(wjet_muon,"EWK-W #rightarrow #mu#nu ","L");
    leg.AddEntry(wjet_ele,"EWK-W #rightarrow e#nu ","L");
    leg.AddEntry(wjet_tau,"EWK-W #rightarrow #tau#nu ","L");
  }

  leg.Draw("same");

  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);
  canvas->cd();
  pad2->Draw();
  pad2->cd();
  
  TH1* ratio_tau = (TH1*) wjet_tau->Clone("ratio_tau");
  TH1* ratio_muon = (TH1*) wjet_muon->Clone("ratio_muon");
  TH1* ratio_ele = (TH1*) wjet_ele->Clone("ratio_ele");

  ratio_tau->Divide(wjet_total);
  ratio_muon->Divide(wjet_total);
  ratio_ele->Divide(wjet_total);

  ratio_muon->GetYaxis()->SetNdivisions(5);
  ratio_muon->GetYaxis()->CenterTitle();
  ratio_muon->GetYaxis()->SetTitleOffset(1.5);
  ratio_muon->GetYaxis()->SetLabelSize(0.04);
  ratio_muon->GetYaxis()->SetTitleSize(0.04);
  ratio_muon->GetXaxis()->SetLabelSize(0.04);
  ratio_muon->GetXaxis()->SetTitleSize(0.05);
  ratio_muon->GetXaxis()->SetTitleOffset(1.1);
  ratio_muon->GetXaxis()->SetTitle(observableLatex.c_str());
  ratio_muon->GetYaxis()->SetTitle("Split / Incl.");
  ratio_muon->Draw("hist");
  ratio_ele->Draw("hist same");
  ratio_tau->Draw("hist same");
  ratio_muon->GetYaxis()->SetRangeUser(0.1,0.8);

  if(not isEWK){
    canvas->SaveAs((outputDIR+"/wjets_comparison_normalized.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/wjets_comparison_normalized.pdf").c_str(),"pdf");
  }
  else{
    canvas->SaveAs((outputDIR+"/wjets_comparison_normalized_ewk.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/wjets_comparison_normalized_ewk.pdf").c_str(),"pdf");
  }
}
