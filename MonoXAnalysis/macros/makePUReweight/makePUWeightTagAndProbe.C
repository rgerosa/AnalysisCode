#include "../CMS_lumi.h"

void makePUWeightTagAndProbe(string tagAndProbeDIR_data, string tagAndProbeDIR_MC, bool isSingleMuon, float lumi){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  TGaxis::SetMaxDigits(3);

  TChain* chain_data = NULL;
  TChain* chain_mc = NULL;
  if(isSingleMuon){
    chain_data = new TChain("muontnptree/fitter_tree");
    chain_mc = new TChain("muontnptree/fitter_tree");
  }
  else{
    chain_data = new TChain("electrontnptree/fitter_tree");
    chain_mc = new TChain("electrontnptree/fitter_tree");
  }

  chain_data->Add((tagAndProbeDIR_data+"/*root").c_str());
  chain_mc->Add((tagAndProbeDIR_MC+"/*root").c_str());
  
  TH1F* histo_data = new TH1F("histo_data","",40,0,40);
  TH1F* histo_mc   = new TH1F("histo_mc","",40,0,40);
  histo_data->Sumw2();
  histo_mc->Sumw2();
  
  if(isSingleMuon){
    chain_data->Draw("nvtx >> histo_data","(looseid && mass > 60 && mass < 120 && eta < 2.4)","goff");
    chain_mc->Draw("nvtx >> histo_mc","(wgt)*(looseid && mass > 60 && mass < 120 && eta < 2.4)","goff");
  }
  else{
    chain_data->Draw("nvtx >> histo_data","(vetoid && mass > 60 && mass < 120 && eta < 2.5)","goff");
    chain_mc->Draw("nvtx >> histo_mc","(wgt)*(vetoid && mass > 60 && mass < 120 && eta < 2.5)","goff");
  }

  
  TCanvas* canvas = new TCanvas("canvas","canvas",600,625);
  canvas->cd();
  
  // weight to the same area
  histo_mc->Scale(1.*histo_data->Integral()/histo_mc->Integral());
  histo_mc->SetLineColor(kBlack);
  histo_mc->GetXaxis()->SetTitle("N_{PV}");
  histo_mc->GetYaxis()->SetTitle("Entries");
  histo_mc->SetFillColor(kGray);
  histo_mc->GetYaxis()->SetRangeUser(0.,max(histo_mc->GetMaximum(),histo_data->GetMaximum())*1.2);
  histo_mc->Draw("hist");
  histo_data->SetMarkerColor(kBlack);
  histo_data->SetLineColor(kBlack);
  histo_data->SetMarkerStyle(20);
  histo_data->Draw("E1Psame");

  CMS_lumi(canvas,Form("%.2f",lumi),true);
  TLegend* legend = new TLegend(0.65,0.65,0.9,0.9);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(histo_data,"Data","EP");
  legend->AddEntry(histo_mc,"Z #rightarrow #mu #mu","F");
  legend->Draw("same");

  canvas->SaveAs("npv_before.pdf","pdf");
  canvas->SaveAs("npv_before.png","png");

  TH1F* weight = (TH1F*) histo_data->Clone("weight");
  weight->Divide(histo_mc);
  TH1F* histo_mc_afterWeight = (TH1F*) histo_mc->Clone("histo_mc_afterWeight");
  for(int iBin = 0; iBin < histo_mc_afterWeight->GetNbinsX(); iBin++){
    histo_mc_afterWeight->SetBinContent(iBin+1,histo_mc_afterWeight->GetBinContent(iBin+1)*weight->GetBinContent(iBin+1));
  }

  histo_mc_afterWeight->Draw("hist");
  histo_data->Draw("E1Psame");
  CMS_lumi(canvas,Form("%.2f",lumi),true);
  legend->Clear();
  legend->AddEntry(histo_data,"Data","EP");
  legend->AddEntry(histo_mc_afterWeight,"Z #rightarrow #mu #mu","F");
  legend->Draw("same");

  canvas->SaveAs("npv_after.pdf","pdf");
  canvas->SaveAs("npv_after.png","png");

  TFile* output = new TFile(Form("purwt_%.2f.root",lumi),"RECREATE");
  output->cd();
  weight->Write("puhist");
  output->Close();

}
