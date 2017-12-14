#include "../CMS_lumi.h"

void makeGoodnessOfFit(string observedFileName,
		       string expectedFileName, 
		       string postfix){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();


  TFile* expectedFile = TFile::Open((expectedFileName).c_str(),"READ");
  TFile* observedFile = TFile::Open((observedFileName).c_str(),"READ");

  TTree* limit_exp = (TTree*) expectedFile->Get("limit");
  TTree* limit_obs = (TTree*) observedFile->Get("limit");

  TH1F* expected_statistics = new TH1F("expected_statistics","",250,0,100);
  limit_exp->Draw("limit >> expected_statistics","","goff");

  TH1F* observed_statistics = new TH1F("observed_statistics","",250,0,100);
  limit_obs->Draw("limit >> observed_statistics","","goff");

  float uppervalue = 100;
  for(int iBin = 0; iBin < expected_statistics->GetNbinsX(); iBin++){
    if(expected_statistics->GetBinContent(iBin+1) != 0) uppervalue = expected_statistics->GetBinCenter(iBin+2);
  }
    

  TCanvas* canvas = new TCanvas ("canvas","",600,600);
  canvas->cd();

  expected_statistics->SetLineColor(kBlack);
  expected_statistics->SetLineWidth(2);
  
  TLine* line = new TLine(observed_statistics->GetMean(),0,observed_statistics->GetMean(),expected_statistics->GetMaximum());
  line->SetLineColor(kRed);
  line->SetLineWidth(2);

  expected_statistics->GetXaxis()->SetTitle("Test statistics");
  expected_statistics->GetYaxis()->SetTitle("Entries");
  expected_statistics->GetXaxis()->SetRangeUser(0.,uppervalue);
  expected_statistics->GetYaxis()->SetRangeUser(0.,expected_statistics->GetMaximum()*1.2);
  expected_statistics->Draw("hist");
  line->Draw("same");

  CMS_lumi(canvas,"35.9");

  cout<<"integral above observed "<<expected_statistics->Integral(expected_statistics->FindBin(observed_statistics->GetMean()),expected_statistics->GetNbinsX())/expected_statistics->Integral()<<endl;

  canvas->SaveAs((postfix+".png").c_str(),"png");
  canvas->SaveAs((postfix+".pdf").c_str(),"pdf");
  
}
