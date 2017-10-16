#include "../CMS_lumi.h"

void drawPlot(const vector<TH1F*> & histos,
	      string postfix,
	      string outputDIR,
	      string yAxisLabel,
	      const vector<string> & label,
	      TH1F* total_unc = NULL){


  TCanvas* canvas = new TCanvas(("canvas_"+postfix).c_str(),"",600,625);
  canvas->cd();

  int icolor = 1;
  double min = 1000;
  double max = -1000;
  for(auto hist: histos){    
    if(icolor == 3) icolor++;
    if(icolor == 5) icolor++;
    if(icolor == 6) icolor++;
    hist->SetLineWidth(2);
    hist->SetLineColor(icolor);
    hist->SetMarkerColor(icolor);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1);
    icolor++;
    if(hist->GetMinimum() < min) min = hist->GetMinimum();
    if(hist->GetMaximum() > max) max = hist->GetMaximum();
    
  }

  histos.at(0)->GetYaxis()->SetRangeUser(0.5,2.5);
  histos.at(0)->GetYaxis()->SetTitle(yAxisLabel.c_str());
  histos.at(0)->GetXaxis()->SetTitle("Boson p_{T} [GeV]"); 
  histos.at(0)->Draw("hist");

  if(total_unc != NULL){
    total_unc->SetFillColor(kGray);
    total_unc->SetLineColor(kGray);
    total_unc->Draw("E2 same");
    histos.at(0)->Draw("hist same");
  }
  
  int ihist = 0;
  for(auto hist: histos){
    if(ihist == 0){
      ihist++;
      continue;
    }
    hist->Draw("hist same");
    hist->Draw("P same");
  }
    
  CMS_lumi(canvas,"");

  TLegend leg (0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  for(size_t lab = 0; lab < label.size(); lab++)
    leg.AddEntry(histos.at(lab),label.at(lab).c_str(),"PL");
  leg.Draw("same");

  canvas->RedrawAxis("sameaxis");
  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");
    
}

void makeKFactorVBFPlots(string inputFileZjet, string inputFileWjet, string outputDIR, bool useSmoothed = false){
  
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  
  TFile* inputFile_zjets = TFile::Open(inputFileZjet.c_str());
  TFile* inputFile_wjets = TFile::Open(inputFileWjet.c_str());

  vector<TH1F*> zjets_kfactor;
  if(not useSmoothed){
    zjets_kfactor.push_back((TH1F*) inputFile_zjets->Get("kfactors_shape/kfactor_vbf_mjj_200_500"));
    zjets_kfactor.push_back((TH1F*) inputFile_zjets->Get("kfactors_shape/kfactor_vbf_mjj_500_1000"));
    zjets_kfactor.push_back((TH1F*) inputFile_zjets->Get("kfactors_shape/kfactor_vbf_mjj_1000_1500"));
    zjets_kfactor.push_back((TH1F*) inputFile_zjets->Get("kfactors_shape/kfactor_vbf_mjj_1500_5000"));
  }
  else{
    zjets_kfactor.push_back((TH1F*) inputFile_zjets->Get("kfactors_shape/kfactor_vbf_mjj_200_500_smoothed"));
    zjets_kfactor.push_back((TH1F*) inputFile_zjets->Get("kfactors_shape/kfactor_vbf_mjj_500_1000_smoothed"));
    zjets_kfactor.push_back((TH1F*) inputFile_zjets->Get("kfactors_shape/kfactor_vbf_mjj_1000_1500_smoothed"));
    zjets_kfactor.push_back((TH1F*) inputFile_zjets->Get("kfactors_shape/kfactor_vbf_mjj_1500_5000_smoothed"));
  }
  vector<TH1F*> wjets_kfactor;
  if(not useSmoothed){
    wjets_kfactor.push_back((TH1F*) inputFile_wjets->Get("kfactors_shape/kfactor_vbf_mjj_200_500"));
    wjets_kfactor.push_back((TH1F*) inputFile_wjets->Get("kfactors_shape/kfactor_vbf_mjj_500_1000"));
    wjets_kfactor.push_back((TH1F*) inputFile_wjets->Get("kfactors_shape/kfactor_vbf_mjj_1000_1500"));
    wjets_kfactor.push_back((TH1F*) inputFile_wjets->Get("kfactors_shape/kfactor_vbf_mjj_1500_5000"));
  }
  else{
    wjets_kfactor.push_back((TH1F*) inputFile_wjets->Get("kfactors_shape/kfactor_vbf_mjj_200_500_smoothed"));
    wjets_kfactor.push_back((TH1F*) inputFile_wjets->Get("kfactors_shape/kfactor_vbf_mjj_500_1000_smoothed"));
    wjets_kfactor.push_back((TH1F*) inputFile_wjets->Get("kfactors_shape/kfactor_vbf_mjj_1000_1500_smoothed"));
    wjets_kfactor.push_back((TH1F*) inputFile_wjets->Get("kfactors_shape/kfactor_vbf_mjj_1500_5000_smoothed"));
  }

  vector<TH1F*> zw_ratios;
  for(size_t ihist = 0; ihist < zjets_kfactor.size(); ihist++){
    zw_ratios.push_back((TH1F*) zjets_kfactor.at(ihist)->Clone(Form("zwratio_%s",zjets_kfactor.at(ihist)->GetName())));
    zw_ratios.back()->Divide(wjets_kfactor.at(ihist));    
  }

  vector<string> label;
  label.push_back("200 < M_{jj} < 500 GeV");
  label.push_back("500 < M_{jj} < 1000 GeV");
  label.push_back("1000 < M_{jj} < 1500 GeV");
  label.push_back("1500 < M_{jj} < 5000 GeV");

  TH1F* uncertaintiy_zw = new TH1F("uncertaintiy_zw","",zw_ratios.at(0)->GetNbinsX(),zw_ratios.at(0)->GetXaxis()->GetXbins()->GetArray());
  for(int iBin = 0; iBin <= uncertaintiy_zw->GetNbinsX()+1; iBin++){
    if(uncertaintiy_zw->GetXaxis()->GetBinCenter(iBin) > 100 and uncertaintiy_zw->GetXaxis()->GetBinCenter(iBin) <= 400) uncertaintiy_zw->SetBinContent(iBin,0.12);
    else if(uncertaintiy_zw->GetXaxis()->GetBinCenter(iBin) > 400 and uncertaintiy_zw->GetXaxis()->GetBinCenter(iBin) <= 600) uncertaintiy_zw->SetBinContent(iBin,0.14);
    else if(uncertaintiy_zw->GetXaxis()->GetBinCenter(iBin) > 600 and uncertaintiy_zw->GetXaxis()->GetBinCenter(iBin) <= 800) uncertaintiy_zw->SetBinContent(iBin,0.15);
    else if(uncertaintiy_zw->GetXaxis()->GetBinCenter(iBin) > 800 and uncertaintiy_zw->GetXaxis()->GetBinCenter(iBin) <= 1000) uncertaintiy_zw->SetBinContent(iBin,0.16);
  }

  // to draw the band
  for(int iBin = 0; iBin < zw_ratios.at(0)->GetNbinsX(); iBin++){
    float value_err = sqrt(uncertaintiy_zw->GetBinError(iBin+1)*uncertaintiy_zw->GetBinError(iBin+1)+pow(uncertaintiy_zw->GetBinContent(iBin+1)*zw_ratios.at(0)->GetBinContent(iBin+1),2));
    uncertaintiy_zw->SetBinError(iBin+1,uncertaintiy_zw->GetBinContent(iBin+1)*zw_ratios.at(0)->GetBinContent(iBin+1));    
    uncertaintiy_zw->SetBinContent(iBin+1,zw_ratios.at(0)->GetBinContent(iBin+1));
  }

  if(not useSmoothed)
    drawPlot(zw_ratios,"zw_ratio_comparison",outputDIR,"Z/W k-factor ratio",label,uncertaintiy_zw);
  else
    drawPlot(zw_ratios,"zw_ratio_comparison_smooth",outputDIR,"Z/W k-factor ratio",label,uncertaintiy_zw);

}
