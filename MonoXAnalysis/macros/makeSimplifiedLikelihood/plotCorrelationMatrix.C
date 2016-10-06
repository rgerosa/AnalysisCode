#include "../makeTemplates/histoUtils.h"
#include "../CMS_lumi.h"

void plotCorrelationMatrix(string inputFile, Category category, bool isZeynep, string outputDIR){

  system(("mkdir -p "+outputDIR).c_str());

  gROOT->SetBatch(kTRUE);
  setTDRStyle();  

  TFile *file = new TFile(inputFile.c_str(),"READ");
  string dir = "ch1_ch1";
  if(category == Category::monoV and not isZeynep)
    dir = "ch2_ch1";
  else if(category == Category::monojet and isZeynep)
    dir = "monojet_signal";
  else if(category == Category::monoV and isZeynep)
    dir = "monov_signal";

  
  TCanvas* canvas = NULL;
  if(category == Category::monoV)
    canvas = new TCanvas("canvas","",900,900);
  else
    canvas = new TCanvas("canvas","",1150,900);

  canvas->cd();  
  canvas->SetGrid();
  canvas->SetRightMargin(0.15);
  canvas->SetBottomMargin(0.24);
  canvas->SetLeftMargin(0.2);

  TH1D* bkg   = (TH1D*)file->Get(("shapes_fit_b/"+dir+"/total_background").c_str());
  TH2D* covar = (TH2D*)file->Get(("shapes_fit_b/"+dir+"/total_covar").c_str());

  TH2D* corr = (TH2D*) covar->Clone("test");
  int nbins = bkg->GetNbinsX();
  // loop on the bin --> square matrix --> multiply by the bin width to get real yields
  for (int b=1;b<=nbins;b++){
    double bw = bkg->GetBinWidth(b);
    for (int j=1;j<=nbins;j++){
      double bj = bkg->GetBinWidth(j);
      covar->SetBinContent(b,j,covar->GetBinContent(b,j)*bw*bj);
    }
  }

  for (int b=1;b<=nbins;b++){
    for (int j=1;j<=nbins;j++){
      // calculate the standard deviation (diagonal term)
      double sigb = TMath::Sqrt(covar->GetBinContent(b,b));
      double sigj = TMath::Sqrt(covar->GetBinContent(j,j));
      corr->SetBinContent(b,j,covar->GetBinContent(b,j)/(sigb*sigj));
      if(b == 1)
	corr->GetYaxis()->SetBinLabel(j,Form("%d < E_{T}^{miss} < %d (GeV)",int(bkg->GetXaxis()->GetBinLowEdge(j)),int(bkg->GetXaxis()->GetBinLowEdge(j+1))));
      if(j == 1) 
	corr->GetXaxis()->SetBinLabel(b,Form("%d < E_{T}^{miss} < %d (GeV)",int(bkg->GetXaxis()->GetBinLowEdge(b)),int(bkg->GetXaxis()->GetBinLowEdge(b+1))));
    }
  }

  corr->GetXaxis()->SetTitle("");
  corr->GetYaxis()->SetTitle("");
  corr->GetZaxis()->SetTitle("Correlation");
  corr->GetXaxis()->LabelsOption("v");
  corr->GetYaxis()->LabelsOption("v");
  corr->GetXaxis()->SetLabelSize(0.027);
  corr->GetYaxis()->SetLabelSize(0.027);
  corr->GetZaxis()->SetTitleOffset(1.25);

  gStyle->SetPaintTextFormat("1.2f");
  corr->Draw("COLZTEXT");
  CMS_lumi(canvas,"12.9",true,true,true,0.05,-0.09);

  if(category == Category::monojet){
    canvas->SaveAs((outputDIR+"/correlation_monojet.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_monojet.png").c_str());
    canvas->SetLogz();
    canvas->SaveAs((outputDIR+"/correlation_monojet_log.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_monojet_log.png").c_str());
    corr->Draw("TEXT");
    canvas->SetLogz(0);
    canvas->SetRightMargin(0.07);
    CMS_lumi(canvas,"12.9",true,true,true,0.05,-0.01);
    canvas->SaveAs((outputDIR+"/correlation_monojet_text.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_monojet_text.png").c_str());
  }
  else{
    canvas->SaveAs((outputDIR+"/correlation_monov.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_monov.png").c_str());
    canvas->SetLogz();
    canvas->SaveAs((outputDIR+"/correlation_monov_log.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_monov_log.png").c_str());
    corr->Draw("TEXT");
    canvas->SetLogz(0);
    canvas->SetRightMargin(0.07);
    CMS_lumi(canvas,"12.9",true,true,true,0.05,-0.01);
    canvas->SaveAs((outputDIR+"/correlation_monov_text.pdf").c_str());
    canvas->SaveAs((outputDIR+"/correlation_monov_text.png").c_str());
  }

}
