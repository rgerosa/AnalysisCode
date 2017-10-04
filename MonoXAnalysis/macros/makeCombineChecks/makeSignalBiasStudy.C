#include "../CMS_lumi.h"

void makeSignalBiasStudy(string inputDIR, string nameToGrep, string outputDIR, string postfix, float signalStrenghtInjected = 0.5){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TCanvas* canvas = new TCanvas ("canvas","",600,600);
  canvas->cd();

  vector<float> pullVal;
  double min = 10000;
  double max = -10000;

  TChain* chain = new TChain("tree_fit_sb","tree_fit_sb");
  chain->Add((inputDIR+"/*"+nameToGrep+"*root").c_str());

  TTreeReader reader (chain);
  TTreeReaderValue<int> fit_status (reader,"fit_status");
  TTreeReaderValue<double> r (reader,"r");
  TTreeReaderValue<double> rLoErr (reader,"rLoErr");
  TTreeReaderValue<double> rHiErr (reader,"rHiErr");
  
  
  while(reader.Next()){

    if(*fit_status != 0) continue;
    double error = (*rHiErr+*rLoErr)/2;
    if((*r-signalStrenghtInjected)/error < -6 or (*r-signalStrenghtInjected)/error > 6) continue;
    pullVal.push_back((*r-signalStrenghtInjected)/error);    

    if(pullVal.back() > max) max = pullVal.back();
    else if(pullVal.back() < min) min = pullVal.back();

  }

  TH1F* biasDistribution = new TH1F("biasDistribution","",int((max-min)/0.2),-4,4);
  biasDistribution->Sumw2();

  for(auto val : pullVal) biasDistribution->Fill(val);

  canvas->cd();  
  biasDistribution->SetLineColor(kBlack);
  biasDistribution->SetLineWidth(2);
  biasDistribution->SetMarkerColor(kBlack);
  biasDistribution->SetMarkerStyle(20);
  biasDistribution->SetMarkerSize(1);
  biasDistribution->GetXaxis()->SetTitle("(#mu_{fit}-#mu_{inj})/#sigma_{#mu}");
  biasDistribution->GetYaxis()->SetTitle("N_{Toys}");
  biasDistribution->Draw("hist");
  TF1* gaus = new TF1("gaus","gaus(0)",min,max);
  gaus->SetLineColor(kRed);
  gaus->SetLineWidth(2);
  biasDistribution->Fit(gaus,"RSM");
  gaus->Draw("L same");
  biasDistribution->GetYaxis()->SetRangeUser(0.,biasDistribution->GetMaximum()*1.7);

  CMS_lumi(canvas,"35.9");

  TPaveText *pt = new TPaveText(0.55,0.62,.9,.9,"NDC");
  pt->SetFillColor(0);
  pt->SetFillStyle(0);
  pt->SetBorderSize(1);
  pt->SetTextAlign(11);
  pt->SetTextFont(42);
  pt->AddText(Form("Integral = %d",int(biasDistribution->Integral())));
  pt->AddText(Form("Mean  = %.2f #pm %.2f",biasDistribution->GetMean(),biasDistribution->GetMeanError()));
  pt->AddText(Form("Fit Mean = %.2f #pm %.2f",gaus->GetParameter(1),gaus->GetParError(1)));
  pt->AddText(Form("RMS = %.2f #pm %.2f",biasDistribution->GetRMS(),biasDistribution->GetRMSError()));
  pt->AddText(Form("Fit RMS = %.2f #pm %.2f",gaus->GetParameter(2),gaus->GetParError(2)));
  pt->Draw();

  canvas->SaveAs((outputDIR+"/biasSignalStrenght_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/biasSignalStrenght_"+postfix+".pdf").c_str(),"pdf");
}
