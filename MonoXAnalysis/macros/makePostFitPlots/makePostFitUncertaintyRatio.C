#include "../CMS_lumi.h"

void makePostFitUncertaintyRatio(string fileName1, string fileName2, string channel){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  

  TFile* inputFile_1 = TFile::Open(fileName1.c_str());
  TFile* inputFile_2 = TFile::Open(fileName2.c_str());

  TH1F* background_1 = (TH1F*) inputFile_1->Get(("shapes_fit_b/"+channel+"/total_background").c_str());
  TH1F* background_2 = (TH1F*) inputFile_2->Get(("shapes_fit_b/"+channel+"/total_background").c_str());

  string postfix;
  if(channel == "ch1") postfix = "_sig";
  else if(channel == "ch2") postfix = "_zmm";
  else if(channel == "ch3") postfix = "_wmn";
  else if(channel == "ch4") postfix = "_zee";
  else if(channel == "ch5") postfix = "_wen";


  TH1F* uncertaintyRatio = (TH1F*) background_1->Clone("uncertaintyRatio");
  uncertaintyRatio->Reset();
  
  for(int iBin = 0; iBin < uncertaintyRatio->GetNbinsX()+1; iBin++){
    uncertaintyRatio->SetBinContent(iBin+1,background_1->GetBinError(iBin+1)/background_2->GetBinError(iBin+1));
  }
  
  TCanvas* canvas = new TCanvas("canvas","canvas",600,600);
  canvas->cd();

  uncertaintyRatio->GetXaxis()->SetTitle("M_{jj} [GeV]");
  uncertaintyRatio->GetYaxis()->SetTitle("Uncertainty ratio");
  uncertaintyRatio->GetYaxis()->SetRangeUser(uncertaintyRatio->GetMinimum()*0.85,uncertaintyRatio->GetMaximum()*1.15);
  uncertaintyRatio->SetLineColor(kRed);
  uncertaintyRatio->SetLineWidth(2);

  uncertaintyRatio->Draw("hist");

  canvas->SaveAs(("postfit_uncertainty"+postfix+".png").c_str(),"png");
  canvas->SaveAs(("postfit_uncertainty"+postfix+".pdf").c_str(),"pdf");
  
}
