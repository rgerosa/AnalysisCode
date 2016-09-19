#include "../CMS_lumi.h"

void makeLikelihoodScan(string inputDirectory, string massID, string outputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputDIR).c_str());
  TChain* chain = new TChain("limit","limit");
  chain->Add((inputDirectory+"/*"+massID+"*root").c_str());

  TTreeReader* reader = new TTreeReader (chain);
  TTreeReaderValue<vector<float> >* muVal = new TTreeReaderValue<vector<float> >(*reader,"muVal");
  TTreeReaderValue<vector<float> >* nllObserved = new TTreeReaderValue<vector<float> >(*reader,"nllObserved");
  TTreeReaderValue<vector<float> >* nllExpected = new TTreeReaderValue<vector<float> >(*reader,"nllExpected");
   
  TGraph* observedLikelihood = new TGraph();
  TGraph* expectedLikelihood = new TGraph();
  while(reader->Next()){
    long int N = 0;
    for(size_t iMu = 0; iMu < (*muVal)->size(); iMu++){      
      if(isnan((*nllObserved)->at(iMu)) or isnan((*nllExpected)->at(iMu))) continue;
      observedLikelihood->SetPoint(N,(*muVal)->at(iMu),(*nllObserved)->at(iMu));
      expectedLikelihood->SetPoint(N,(*muVal)->at(iMu),(*nllExpected)->at(iMu));
      N++;
    }    
  }
  
  //Min of observed and expected likelihood
  float minNllObs = TMath::MinElement(observedLikelihood->GetN(),observedLikelihood->GetY());
  float minNllExp = TMath::MinElement(expectedLikelihood->GetN(),expectedLikelihood->GetY());
  float minMu = TMath::MinElement(observedLikelihood->GetN(),observedLikelihood->GetX());
  float maxMu = TMath::MaxElement(observedLikelihood->GetN(),observedLikelihood->GetX());

  float muRange = min(fabs(minMu),maxMu);

  // plotting
  TCanvas* canvas = new TCanvas("canvas","canvas",600,650);
  canvas->cd();
  TH1F* frame = canvas->DrawFrame(-muRange,0,muRange,max(TMath::MaxElement(observedLikelihood->GetN(),observedLikelihood->GetY()),TMath::MaxElement(expectedLikelihood->GetN(),expectedLikelihood->GetY())));
  frame->GetXaxis()->SetTitle("signal strength #mu");
  frame->GetYaxis()->SetTitle("-2*Log(L)");
  frame->Draw();
  CMS_lumi(canvas,"12.9");
  observedLikelihood->SetLineColor(kBlack);
  observedLikelihood->SetLineWidth(2);
  observedLikelihood->Draw("Csame");
  expectedLikelihood->SetLineColor(kRed);
  expectedLikelihood->SetLineWidth(2);
  expectedLikelihood->Draw("Csame");

  TLegend leg(0.7,0.75,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(observedLikelihood,"Observed","L");
  leg.AddEntry(expectedLikelihood,"Expected","L");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/likelihoodScan_massid_"+massID+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/likelihoodScan_massid_"+massID+".pdf").c_str(),"pdf");
}
