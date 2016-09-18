#include "../CMS_lumi.h"

void makeLikelihoodScan(string inputDirectory, string massID, string outputDIR, string inputDirectoryNoCorr = ""){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputDIR).c_str());
  TChain* chain = new TChain("limit","limit");
  chain->Add((inputDirectory+"/*"+massID+"*root").c_str());

  TChain* chainNoCorr = NULL;
  if(inputDirectoryNoCorr != ""){
    chainNoCorr = new TChain("limit","limit");
    chainNoCorr->Add((inputDirectoryNoCorr+"/*"+massID+"*root").c_str());
  }

  TTreeReader* reader = new TTreeReader (chain);
  TTreeReaderValue<vector<float> >* muVal = new TTreeReaderValue<vector<float> >(*reader,"muVal");
  TTreeReaderValue<vector<float> >* nllObserved = new TTreeReaderValue<vector<float> >(*reader,"nllObserved");
  TTreeReaderValue<vector<float> >* nllExpected = new TTreeReaderValue<vector<float> >(*reader,"nllExpected");
   
  TTreeReader* readerNoCorr = NULL;
  TTreeReaderValue<vector<float> >* muValNoCorr = NULL;
  TTreeReaderValue<vector<float> >* nllObservedNoCorr = NULL;
  TTreeReaderValue<vector<float> >* nllExpectedNoCorr = NULL;

  if(chainNoCorr != NULL){
    readerNoCorr = new TTreeReader (chainNoCorr);
    muValNoCorr = new TTreeReaderValue<vector<float> >(*readerNoCorr,"muVal");
    nllObservedNoCorr = new TTreeReaderValue<vector<float> >(*readerNoCorr,"nllObserved");
    nllExpectedNoCorr = new TTreeReaderValue<vector<float> >(*readerNoCorr,"nllExpected");
  }

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

  TGraph* observedLikelihoodNoCorr = NULL;
  TGraph* expectedLikelihoodNoCorr = NULL;
  if(chainNoCorr != NULL){
    observedLikelihoodNoCorr = new TGraph();
    expectedLikelihoodNoCorr = new TGraph();
    while(readerNoCorr->Next()){
      long int N = 0;
      for(size_t iMu = 0; iMu < (*muVal)->size(); iMu++){
	if(isnan((*nllObservedNoCorr)->at(iMu)) or isnan((*nllExpectedNoCorr)->at(iMu))) continue;
	observedLikelihoodNoCorr->SetPoint(N,(*muValNoCorr)->at(iMu),(*nllObservedNoCorr)->at(iMu));
	expectedLikelihoodNoCorr->SetPoint(N,(*muValNoCorr)->at(iMu),(*nllExpectedNoCorr)->at(iMu));
	N++;
      }
    }    
  }
  
  //Min of observed and expected likelihood
  float minNllObs = TMath::MinElement(observedLikelihood->GetN(),observedLikelihood->GetY());
  float minNllExp = TMath::MinElement(expectedLikelihood->GetN(),expectedLikelihood->GetY());
  float minMu = TMath::MinElement(observedLikelihood->GetN(),observedLikelihood->GetX());
  float maxMu = TMath::MaxElement(observedLikelihood->GetN(),observedLikelihood->GetX());

  float minNllObsNoCorr = 0;
  float minNllExpNoCorr = 0;
  if(observedLikelihoodNoCorr!= NULL and expectedLikelihoodNoCorr != NULL){
    minNllObsNoCorr = TMath::MinElement(observedLikelihoodNoCorr->GetN(),observedLikelihoodNoCorr->GetY());
    minNllExpNoCorr = TMath::MinElement(expectedLikelihoodNoCorr->GetN(),expectedLikelihoodNoCorr->GetY());    
  }


  // plotting
  TCanvas* canvas = new TCanvas("canvas","canvas",600,650);
  canvas->cd();
  TH1F* frame = canvas->DrawFrame(minMu,0,maxMu,max(TMath::MaxElement(observedLikelihood->GetN(),observedLikelihood->GetY()),TMath::MaxElement(expectedLikelihood->GetN(),expectedLikelihood->GetY())));
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

  if(observedLikelihoodNoCorr != NULL){
    observedLikelihoodNoCorr->SetLineColor(kBlue);
    observedLikelihoodNoCorr->SetLineWidth(2);
    observedLikelihoodNoCorr->Draw("Csame");
  }

  if(expectedLikelihoodNoCorr != NULL){
    expectedLikelihoodNoCorr->SetLineColor(kGreen+1);
    expectedLikelihoodNoCorr->SetLineWidth(2);
    expectedLikelihoodNoCorr->Draw("Csame");
  }

  TLegend leg(0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(observedLikelihood,"Observed","L");
  leg.AddEntry(expectedLikelihood,"Expected","L");
  if(observedLikelihoodNoCorr != NULL)
    leg.AddEntry(observedLikelihoodNoCorr,"Observed no correlation","L");
  if(observedLikelihoodNoCorr != NULL)
    leg.AddEntry(expectedLikelihoodNoCorr,"Expected no correlation","L");
  leg.Draw("same");
  
  canvas->SaveAs((outputDIR+"/likelihoodScan_massid_"+massID+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/likelihoodScan_massid_"+massID+".pdf").c_str(),"pdf");
}
