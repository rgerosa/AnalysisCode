#include "../CMS_lumi.h"

void drawPlot(TH1F* histo1,
	      TH1F* histo2,
	      TH1F* histo3,
	      string postfix,
	      string outputDIR,
	      string yAxisLabel,
	      TH1F* total_unc = NULL){


  TCanvas* canvas = new TCanvas(("canvas_"+postfix).c_str(),"",600,625);
  canvas->cd();
  histo1->SetLineWidth(2);
  histo2->SetLineWidth(2);
  histo3->SetLineWidth(2);

  histo1->SetLineColor(kBlack);
  histo2->SetLineColor(kRed);
  histo3->SetLineColor(kBlue);

  histo1->SetMarkerColor(kBlack);
  histo2->SetMarkerColor(kRed);
  histo3->SetMarkerColor(kBlue);

  histo1->SetMarkerStyle(20);
  histo2->SetMarkerStyle(20);
  histo3->SetMarkerStyle(20);

  histo1->SetMarkerSize(1);
  histo2->SetMarkerSize(1);
  histo3->SetMarkerSize(1);

  histo1->GetYaxis()->SetRangeUser(min(histo1->GetMinimum(),min(histo2->GetMinimum(),histo3->GetMinimum()))*0.8,
				   max(histo1->GetMaximum(),max(histo2->GetMaximum(),histo3->GetMaximum()))*1.1);
  
  histo1->GetYaxis()->SetTitle(yAxisLabel.c_str());
  histo1->GetXaxis()->SetTitle("Boson p_{T} [GeV]");

 
  histo1->Draw("hist");

  if(total_unc != NULL){
    total_unc->SetFillColor(kGray);
    total_unc->SetLineColor(kGray);
    total_unc->Draw("E2 same");
    histo1->Draw("hist same");
  }

  histo1->Draw("P same");
  histo2->Draw("hist same");
  histo2->Draw("P same");
  histo3->Draw("hist same");
  histo3->Draw("P same");

 

  CMS_lumi(canvas,"");

  TLegend leg (0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  
  leg.AddEntry(histo1,"Two jet","PL");
  leg.AddEntry(histo2,"VBF","PL");
  leg.AddEntry(histo3,"VBF + #Delta#phi(j_{1},j_{2})","PL");
  leg.Draw("same");

  canvas->RedrawAxis("sameaxis");
  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");
    
}

void makeKFactorVBFPlots(string outputDIR){

  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  gROOT->SetBatch(kTRUE);
  
  TFile* inputFile_zjets_twojet = TFile::Open("kfactor_vbf_zjets_nodphijj/kfactor_znn.root");
  TFile* inputFile_zjets_vbf    = TFile::Open("kfactor_vbf_zjets_nodphijj/kfactor_znn.root");
  TFile* inputFile_zjets_vbf_dphijj = TFile::Open("kfactor_vbf_zjets_dphijj0p6/kfactor_znn.root");

  TFile* inputFile_wjets_twojet = TFile::Open("kfactor_vbf_wjets_nodphijj/kfactor_wjet.root");
  TFile* inputFile_wjets_vbf    = TFile::Open("kfactor_vbf_wjets_nodphijj/kfactor_wjet.root");
  TFile* inputFile_wjets_vbf_dphijj = TFile::Open("kfactor_vbf_wjets_dphijj0p6/kfactor_wjet.root");

  TH1F* kfactor_zjets_twojet     = (TH1F*) ((TH1F*) inputFile_zjets_twojet->Get("bosonPt_NLO_twojet"))->Clone("kfactor_zjets_twojet");
  TH1F* kfactor_zjets_vbf        = (TH1F*) ((TH1F*) inputFile_zjets_vbf->Get("bosonPt_NLO_vbf"))->Clone("kfactor_zjets_vbf");
  TH1F* kfactor_zjets_vbf_dphijj = (TH1F*) ((TH1F*) inputFile_zjets_vbf_dphijj->Get("bosonPt_NLO_vbf"))->Clone("kfactor_zjets_vbf_dphijj");

  kfactor_zjets_twojet->Divide((TH1F*) inputFile_zjets_twojet->Get("bosonPt_LO_twojet"));
  kfactor_zjets_vbf->Divide((TH1F*) inputFile_zjets_vbf->Get("bosonPt_LO_vbf"));
  kfactor_zjets_vbf_dphijj->Divide((TH1F*) inputFile_zjets_vbf_dphijj->Get("bosonPt_LO_vbf"));


  TH1F* kfactor_wjets_twojet     = (TH1F*) ((TH1F*) inputFile_wjets_twojet->Get("bosonPt_NLO_twojet"))->Clone("kfactor_wjets_twojet");
  TH1F* kfactor_wjets_vbf        = (TH1F*) ((TH1F*) inputFile_wjets_vbf->Get("bosonPt_NLO_vbf"))->Clone("kfactor_wjets_vbf");
  TH1F* kfactor_wjets_vbf_dphijj = (TH1F*) ((TH1F*) inputFile_wjets_vbf_dphijj->Get("bosonPt_NLO_vbf"))->Clone("kfactor_wjets_vbf_dphijj");

  kfactor_wjets_twojet->Divide((TH1F*) inputFile_wjets_twojet->Get("bosonPt_LO_twojet"));
  kfactor_wjets_vbf->Divide((TH1F*) inputFile_wjets_vbf->Get("bosonPt_LO_vbf"));
  kfactor_wjets_vbf_dphijj->Divide((TH1F*) inputFile_wjets_vbf_dphijj->Get("bosonPt_LO_vbf"));

  // draw k-factors
  drawPlot(kfactor_zjets_twojet,kfactor_zjets_vbf,kfactor_zjets_vbf_dphijj,"kfact_zjets",outputDIR,"k-factor");
  drawPlot(kfactor_wjets_twojet,kfactor_wjets_vbf,kfactor_wjets_vbf_dphijj,"kfact_wjets",outputDIR,"k-factor");

  // make Z/W ratio

  TH1F* ZW_ratio_twojet = (TH1F*) kfactor_zjets_twojet->Clone("ZW_ratio_twojet");
  ZW_ratio_twojet->Divide(kfactor_wjets_twojet);
  TH1F* ZW_ratio_vbf = (TH1F*) kfactor_zjets_vbf->Clone("ZW_ratio_vbf");
  ZW_ratio_vbf->Divide(kfactor_wjets_vbf);
  TH1F* ZW_ratio_vbf_dphijj = (TH1F*) kfactor_zjets_vbf_dphijj->Clone("ZW_ratio_vbf_dphijj");
  ZW_ratio_vbf_dphijj->Divide(kfactor_wjets_vbf_dphijj);

  //make unc band
  TFile* uncertainty_zw = TFile::Open("/afs/cern.ch/user/r/rgerosa/work/MONOJET_ANALYSIS/CMSSW_8_0_20_patch1/src/AnalysisCode/MonoXAnalysis/macros/makeTemplates/templates_monojet_moriond/templates_shape_met.root");
  TH1F*  ewk_band = (TH1F*) uncertainty_zw->FindObjectAny("ZW_EWK_met");
  TH1F*  factscale_band = (TH1F*) uncertainty_zw->FindObjectAny("ZW_FactScale2_met");
  TH1F*  renscale_band = (TH1F*) uncertainty_zw->FindObjectAny("ZW_RenScale2_met");
  TH1F*  pdf_band = (TH1F*) uncertainty_zw->FindObjectAny("ZW_PDF_met");

  TH1F* total_unc = (TH1F*) ZW_ratio_vbf_dphijj->Clone("total_unc");
  for(int iBin = 0; iBin < total_unc->GetNbinsX()+1; iBin++){
    if(ZW_ratio_vbf_dphijj->GetBinCenter(iBin+1) >= 200){
      total_unc->SetBinError(iBin+1,sqrt(pow(total_unc->GetBinContent(iBin+1)*ewk_band->GetBinContent(ewk_band->FindBin(ZW_ratio_vbf_dphijj->GetBinCenter(iBin+1))),2)+
					 pow(total_unc->GetBinContent(iBin+1)*factscale_band->GetBinContent(factscale_band->FindBin(ZW_ratio_vbf_dphijj->GetBinCenter(iBin+1))),2)+
					 pow(total_unc->GetBinContent(iBin+1)*renscale_band->GetBinContent(renscale_band->FindBin(ZW_ratio_vbf_dphijj->GetBinCenter(iBin+1))),2)+
					 pow(total_unc->GetBinContent(iBin+1)*pdf_band->GetBinContent(pdf_band->FindBin(ZW_ratio_vbf_dphijj->GetBinCenter(iBin+1))),2)));
    }
    else
      total_unc->SetBinError(iBin+1,sqrt(pow(total_unc->GetBinContent(iBin+1)*ewk_band->GetBinContent(ewk_band->FindBin(ZW_ratio_vbf_dphijj->GetBinCenter(iBin+2))),2)+
					 pow(total_unc->GetBinContent(iBin+1)*factscale_band->GetBinContent(factscale_band->FindBin(ZW_ratio_vbf_dphijj->GetBinCenter(iBin+2))),2)+
					 pow(total_unc->GetBinContent(iBin+1)*renscale_band->GetBinContent(renscale_band->FindBin(ZW_ratio_vbf_dphijj->GetBinCenter(iBin+2))),2)+
					 pow(total_unc->GetBinContent(iBin+1)*pdf_band->GetBinContent(pdf_band->FindBin(ZW_ratio_vbf_dphijj->GetBinCenter(iBin+2))),2)));
  }
  
  drawPlot(ZW_ratio_twojet,ZW_ratio_vbf,ZW_ratio_vbf_dphijj,"ZW_ratio",outputDIR,"Z/W",total_unc);

}
