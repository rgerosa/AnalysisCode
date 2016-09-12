#include "../makeTemplates/histoUtils.h"
#include "../CMS_lumi.h"

void plotCorrelationMatrix(string inputFile, Category category, bool isCombo){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();  

  TFile *file = new TFile(inputFile.c_str(),"READ");
  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();  
  string dir = "ch1_ch1";
  if(category == Category::monoV)
    dir = "ch2_ch1";

  if(isCombo == false)
    dir = "ch1";


  canvas->SetRightMargin(0.15);

  TH1D* bkg   = (TH1D*)file->Get(("shapes_fit_b/"+dir+"/total_background").c_str());
  TH2D* covar = (TH2D*)file->Get(("shapes_fit_b/"+dir+"/total_covar").c_str());

  TH2D* corr = (TH2D*) covar->Clone("test");
  corr->GetZaxis()->SetTitle("Correlation");
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
    }
  }

  gStyle->SetPaintTextFormat("1.2f");
  corr->Draw("COLZTEXT");
  CMS_lumi(canvas,"7.6",true);
  canvas->SaveAs("correlation.pdf");
  canvas->SaveAs("correlation.png");
  canvas->SetLogz();
  canvas->SaveAs("correlation_log.pdf");
  canvas->SaveAs("correlation_log.png");

}
