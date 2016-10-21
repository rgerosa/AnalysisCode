#include "../CMS_lumi.h"

void makeLikelihoodScanComparison(string outputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  // combine file for full model
  TFile* combineObs = TFile::Open("ResultsSimplifiedLikelihood/forLikelihoodScan/higgsCombinefullLikelihood_Obs.MultiDimFit.mH120.000000.root");
  TTree*  limitCombineObs = (TTree*) combineObs->Get("limit");

  TFile* combineExp = TFile::Open("ResultsSimplifiedLikelihood/forLikelihoodScan/higgsCombinefullLikelihood_Exp.MultiDimFit.mH120.000000.root");
  TTree*  limitCombineExp = (TTree*) combineExp->Get("limit");
  
  // file for simplified likelihood
  TFile* simpLikelihoodFullCorr = TFile::Open("ResultsSimplifiedLikelihood/forLikelihoodScan/simplifiedLikelihood_combined_80018000001_mask_fullcorr.root");
  TTree* limitSimpLikelihoodFullCorr = (TTree*) simpLikelihoodFullCorr->Get("limit");
  
  TFile* simpLikelihood = TFile::Open("ResultsSimplifiedLikelihood/forLikelihoodScan/simplifiedLikelihood_combined_80018000001_mask_corr.root");
  TTree* limitSimpLikelihood = (TTree*) simpLikelihood->Get("limit");

  TFile* simpLikelihoodNoCorr = TFile::Open("ResultsSimplifiedLikelihood/forLikelihoodScan/simplifiedLikelihood_combined_80018000001_mask_nocorr.root");
  TTree* limitSimpLikelihoodNoCorr = (TTree*) simpLikelihoodNoCorr->Get("limit");

  // file for simplified likelihood with SR
  TFile* simpLikelihoodSR = TFile::Open("ResultsSimplifiedLikelihood/forLikelihoodScan/simplifiedLikelihood_combined_80018000001_SR.root");
  TTree* limitSimpLikelihoodSR = (TTree*) simpLikelihoodSR->Get("limit");

  TFile* simpLikelihoodSRNoCorr = TFile::Open("ResultsSimplifiedLikelihood/forLikelihoodScan/simplifiedLikelihood_combined_80018000001_SR_nocorr.root");
  TTree* limitSimpLikelihoodSRNoCorr = (TTree*) simpLikelihoodSRNoCorr->Get("limit");

  system(("mkdir -p "+outputDIR).c_str());

  // Combine
  TTreeReader* readerCombineObs = new TTreeReader (limitCombineObs);
  TTreeReader* readerCombineExp = new TTreeReader (limitCombineExp);
  TTreeReaderValue<float >* muValObs = new TTreeReaderValue<float >(*readerCombineObs,"r");
  TTreeReaderValue<float >* muValExp = new TTreeReaderValue<float >(*readerCombineExp,"r");
  TTreeReaderValue<float >* nllObserved = new TTreeReaderValue<float >(*readerCombineObs,"deltaNLL");
  TTreeReaderValue<float >* nllExpected = new TTreeReaderValue<float >(*readerCombineExp,"deltaNLL");
   
  TGraph* observedLikelihoodComb = new TGraph();
  TGraph* expectedLikelihoodComb = new TGraph();

  long int N = 0;
  while(readerCombineObs->Next()){
    if(**nllObserved == 0) continue;
    if(isnan(**nllObserved)) continue;
    observedLikelihoodComb->SetPoint(N,**muValObs,**nllObserved);
    N++;
  }
  N = 0;
  while(readerCombineExp->Next()){
    if(**nllExpected == 0) continue;
    if(isnan(**nllExpected)) continue;
    expectedLikelihoodComb->SetPoint(N,**muValExp,**nllExpected);
    N++;
  }

  //SL Likelihood full correlation
  TTreeReader* readerSimpLikelihoodFullCorr = new TTreeReader  (limitSimpLikelihoodFullCorr);
  TTreeReaderValue<vector<float> >* muVal = new TTreeReaderValue<vector<float> >(*readerSimpLikelihoodFullCorr,"muVal");
  TTreeReaderValue<vector<float> >* nllObs = new TTreeReaderValue<vector<float> >(*readerSimpLikelihoodFullCorr,"nllObserved");
  TTreeReaderValue<vector<float> >* nllExp = new TTreeReaderValue<vector<float> >(*readerSimpLikelihoodFullCorr,"nllExpected");
  
  TGraph* observedLikelihoodSLFullCorr = new TGraph();
  TGraph* expectedLikelihoodSLFullCorr = new TGraph();

  while(readerSimpLikelihoodFullCorr->Next()){
    N = 0;
    for(size_t imu = 0; imu < (*muVal)->size(); imu++){
      if(isnan((*nllObs)->at(imu)) or isnan((*nllExp)->at(imu)))  continue;
      observedLikelihoodSLFullCorr->SetPoint(N,(*muVal)->at(imu),(*nllObs)->at(imu));
      expectedLikelihoodSLFullCorr->SetPoint(N,(*muVal)->at(imu),(*nllExp)->at(imu));
      N++;
    }
  }

  
  // SL Likelihood
  TTreeReader* readerSimpLikelihood = new TTreeReader (limitSimpLikelihood);
  muVal = new TTreeReaderValue<vector<float> >(*readerSimpLikelihood,"muVal");
  nllObs = new TTreeReaderValue<vector<float> >(*readerSimpLikelihood,"nllObserved");
  nllExp = new TTreeReaderValue<vector<float> >(*readerSimpLikelihood,"nllExpected");
   
  TGraph* observedLikelihoodSL = new TGraph();
  TGraph* expectedLikelihoodSL = new TGraph();

  while(readerSimpLikelihood->Next()){
    N = 0;
    for(size_t imu = 0; imu < (*muVal)->size(); imu++){
      if(isnan((*nllObs)->at(imu)) or isnan((*nllExp)->at(imu)))  continue;
      observedLikelihoodSL->SetPoint(N,(*muVal)->at(imu),(*nllObs)->at(imu));
      expectedLikelihoodSL->SetPoint(N,(*muVal)->at(imu),(*nllExp)->at(imu));
      N++;
    }
  }

  // SL Likelihood no corr
  TTreeReader* readerSimpLikelihoodNoCorr = new TTreeReader (limitSimpLikelihoodNoCorr);
  muVal = new TTreeReaderValue<vector<float> >(*readerSimpLikelihoodNoCorr,"muVal");
  nllObs = new TTreeReaderValue<vector<float> >(*readerSimpLikelihoodNoCorr,"nllObserved");
  nllExp = new TTreeReaderValue<vector<float> >(*readerSimpLikelihoodNoCorr,"nllExpected");
  
  TGraph* observedLikelihoodSLNoCorr = new TGraph();
  TGraph* expectedLikelihoodSLNoCorr = new TGraph();

  while(readerSimpLikelihoodNoCorr->Next()){
    N = 0;
    for(size_t imu = 0; imu < (*muVal)->size(); imu++){
      if(isnan((*nllObs)->at(imu)) or isnan((*nllExp)->at(imu)))  continue;
      observedLikelihoodSLNoCorr->SetPoint(N,(*muVal)->at(imu),(*nllObs)->at(imu));
      expectedLikelihoodSLNoCorr->SetPoint(N,(*muVal)->at(imu),(*nllExp)->at(imu));
      N++;
    }
  }

  // SL Likelihood SR
  TTreeReader* readerSimpLikelihoodSR = new TTreeReader (limitSimpLikelihoodSR);
  muVal = new TTreeReaderValue<vector<float> >(*readerSimpLikelihoodSR,"muVal");
  nllObs = new TTreeReaderValue<vector<float> >(*readerSimpLikelihoodSR,"nllObserved");
  nllExp = new TTreeReaderValue<vector<float> >(*readerSimpLikelihoodSR,"nllExpected");
  
  TGraph* observedLikelihoodSLSR = new TGraph();
  TGraph* expectedLikelihoodSLSR = new TGraph();

  while(readerSimpLikelihoodSR->Next()){
    N = 0;
    for(size_t imu = 0; imu < (*muVal)->size(); imu++){
      if(isnan((*nllObs)->at(imu)) or isnan((*nllExp)->at(imu)))  continue;
      observedLikelihoodSLSR->SetPoint(N,(*muVal)->at(imu),(*nllObs)->at(imu));
      expectedLikelihoodSLSR->SetPoint(N,(*muVal)->at(imu),(*nllExp)->at(imu));
      N++;
    }
  }

  // SL Likelihood no corr
  TTreeReader* readerSimpLikelihoodSRNoCorr = new TTreeReader (limitSimpLikelihoodSRNoCorr);
  muVal = new TTreeReaderValue<vector<float> >(*readerSimpLikelihoodSRNoCorr,"muVal");
  nllObs = new TTreeReaderValue<vector<float> >(*readerSimpLikelihoodSRNoCorr,"nllObserved");
  nllExp = new TTreeReaderValue<vector<float> >(*readerSimpLikelihoodSRNoCorr,"nllExpected");
  
  TGraph* observedLikelihoodSLSRNoCorr = new TGraph();
  TGraph* expectedLikelihoodSLSRNoCorr = new TGraph();

  while(readerSimpLikelihoodSRNoCorr->Next()){
    N = 0;
    for(size_t imu = 0; imu < (*muVal)->size(); imu++){
      if(isnan((*nllObs)->at(imu)) or isnan((*nllExp)->at(imu)))  continue;
      observedLikelihoodSLSRNoCorr->SetPoint(N,(*muVal)->at(imu),(*nllObs)->at(imu));
      expectedLikelihoodSLSRNoCorr->SetPoint(N,(*muVal)->at(imu),(*nllExp)->at(imu));
      N++;
    }
  }

  
  
  //Min of observed and expected likelihood
  float minMuObsComb = TMath::MinElement(observedLikelihoodComb->GetN(),observedLikelihoodComb->GetX());
  float maxMuObsComb = TMath::MaxElement(observedLikelihoodComb->GetN(),observedLikelihoodComb->GetX());
  float minMuExpComb = TMath::MinElement(expectedLikelihoodComb->GetN(),expectedLikelihoodComb->GetX());
  float maxMuExpComb = TMath::MaxElement(expectedLikelihoodComb->GetN(),expectedLikelihoodComb->GetX());
  float rangeComb = min(min(fabs(minMuObsComb),maxMuObsComb),min(fabs(minMuExpComb),maxMuExpComb));

  float minMuObsSLFullCorr = TMath::MinElement(observedLikelihoodSLFullCorr->GetN(),observedLikelihoodSLFullCorr->GetX());
  float maxMuObsSLFullCorr = TMath::MaxElement(observedLikelihoodSLFullCorr->GetN(),observedLikelihoodSLFullCorr->GetX());
  float minMuExpSLFullCorr = TMath::MinElement(expectedLikelihoodSLFullCorr->GetN(),expectedLikelihoodSLFullCorr->GetX());
  float maxMuExpSLFullCorr = TMath::MaxElement(expectedLikelihoodSLFullCorr->GetN(),expectedLikelihoodSLFullCorr->GetX()); 
  float rangeSLFullCorr = min(min(fabs(minMuObsSLFullCorr),maxMuObsSLFullCorr),min(fabs(minMuExpSLFullCorr),maxMuExpSLFullCorr));

  float minMuObsSL = TMath::MinElement(observedLikelihoodSL->GetN(),observedLikelihoodSL->GetX());
  float maxMuObsSL = TMath::MaxElement(observedLikelihoodSL->GetN(),observedLikelihoodSL->GetX());
  float minMuExpSL = TMath::MinElement(expectedLikelihoodSL->GetN(),expectedLikelihoodSL->GetX());
  float maxMuExpSL = TMath::MaxElement(expectedLikelihoodSL->GetN(),expectedLikelihoodSL->GetX()); 
  float rangeSL = min(min(fabs(minMuObsSL),maxMuObsSL),min(fabs(minMuExpSL),maxMuExpSL));

  float minMuObsSLNoCorr = TMath::MinElement(observedLikelihoodSLNoCorr->GetN(),observedLikelihoodSLNoCorr->GetX());
  float maxMuObsSLNoCorr = TMath::MaxElement(observedLikelihoodSLNoCorr->GetN(),observedLikelihoodSLNoCorr->GetX());
  float minMuExpSLNoCorr = TMath::MinElement(expectedLikelihoodSLNoCorr->GetN(),expectedLikelihoodSLNoCorr->GetX());
  float maxMuExpSLNoCorr = TMath::MaxElement(expectedLikelihoodSLNoCorr->GetN(),expectedLikelihoodSLNoCorr->GetX());
  float rangeSLNoCorr = min(min(fabs(minMuObsSLNoCorr),maxMuObsSLNoCorr),min(fabs(minMuExpSLNoCorr),maxMuExpSLNoCorr));

  float minMuObsSLSR = TMath::MinElement(observedLikelihoodSLSR->GetN(),observedLikelihoodSLSR->GetX());
  float maxMuObsSLSR = TMath::MaxElement(observedLikelihoodSLSR->GetN(),observedLikelihoodSLSR->GetX());
  float minMuExpSLSR = TMath::MinElement(expectedLikelihoodSLSR->GetN(),expectedLikelihoodSLSR->GetX());
  float maxMuExpSLSR = TMath::MaxElement(expectedLikelihoodSLSR->GetN(),expectedLikelihoodSLSR->GetX());
  float rangeSLSR = min(min(fabs(minMuObsSLSR),maxMuObsSLSR),min(fabs(minMuExpSLSR),maxMuExpSLSR));

  float minMuObsSLSRNoCorr = TMath::MinElement(observedLikelihoodSLSRNoCorr->GetN(),observedLikelihoodSLSRNoCorr->GetX());
  float maxMuObsSLSRNoCorr = TMath::MaxElement(observedLikelihoodSLSRNoCorr->GetN(),observedLikelihoodSLSRNoCorr->GetX());
  float minMuExpSLSRNoCorr = TMath::MinElement(expectedLikelihoodSLSRNoCorr->GetN(),expectedLikelihoodSLSRNoCorr->GetX());
  float maxMuExpSLSRNoCorr = TMath::MaxElement(expectedLikelihoodSLSRNoCorr->GetN(),expectedLikelihoodSLSRNoCorr->GetX());
  float rangeSLSRNoCorr = min(min(fabs(minMuObsSLSRNoCorr),maxMuObsSLSRNoCorr),min(fabs(minMuExpSLSRNoCorr),maxMuExpSLSRNoCorr));

  float muRange = min(rangeComb,min(rangeSLFullCorr,min(rangeSL,min(rangeSLNoCorr,min(rangeSLSR,rangeSLSRNoCorr)))));

  // plotting
  TCanvas* canvas = new TCanvas("canvas","canvas",600,650);
  canvas->cd();
  TH1F* frame = canvas->DrawFrame(-muRange,0,muRange,TMath::MaxElement(observedLikelihoodComb->GetN(),observedLikelihoodComb->GetY()));
  frame->GetXaxis()->SetTitle("signal strength #mu");
  frame->GetYaxis()->SetTitle("-2*Log(L)");
  frame->Draw();
  CMS_lumi(canvas,"12.9");
  observedLikelihoodComb->SetLineColor(kBlack);
  observedLikelihoodComb->SetLineWidth(2);
  observedLikelihoodComb->Draw("Csame");
  observedLikelihoodSLFullCorr->SetLineColor(kYellow);
  observedLikelihoodSLFullCorr->SetLineWidth(2);
  observedLikelihoodSLFullCorr->Draw("Csame");
  observedLikelihoodSL->SetLineColor(kBlue);
  observedLikelihoodSL->SetLineWidth(2);
  observedLikelihoodSL->Draw("Csame");
  observedLikelihoodSLNoCorr->SetLineColor(kRed);
  observedLikelihoodSLNoCorr->SetLineWidth(2);
  observedLikelihoodSLNoCorr->Draw("Csame");
  observedLikelihoodSLSR->SetLineColor(kGreen+1);
  observedLikelihoodSLSR->SetLineWidth(2);
  observedLikelihoodSLSR->Draw("Csame");
  observedLikelihoodSLSRNoCorr->SetLineColor(kOrange+1);
  observedLikelihoodSLSRNoCorr->SetLineWidth(2);
  observedLikelihoodSLSRNoCorr->Draw("Csame");

  TLegend leg(0.45,0.5,0.9,0.9);
  leg.SetBorderSize(0);
  leg.AddEntry(observedLikelihoodComb,"Observed full Likelihood","L");
  leg.AddEntry(observedLikelihoodSLFullCorr,"Observed Simp. Likelihood full correlation","L");
  leg.AddEntry(observedLikelihoodSL,"Observed Simp. Likelihood","L");
  leg.AddEntry(observedLikelihoodSLNoCorr,"Observed Simp. Likelihood w/o corr.","L");
  leg.AddEntry(observedLikelihoodSLSR,"Observed Fit with SR","L");
  leg.AddEntry(observedLikelihoodSLSRNoCorr,"Observed Fit with SR w/o corr.","L");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/likelihoodComparison_observed.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/likelihoodComparison_observed.pdf").c_str(),"pdf");

  frame->Draw();
  CMS_lumi(canvas,"12.9");
  expectedLikelihoodComb->SetLineColor(kBlack);
  expectedLikelihoodComb->SetLineWidth(2);
  expectedLikelihoodComb->Draw("Csame");
  expectedLikelihoodSLFullCorr->SetLineColor(kYellow);
  expectedLikelihoodSLFullCorr->SetLineWidth(2);
  expectedLikelihoodSLFullCorr->Draw("Csame");
  expectedLikelihoodSL->SetLineColor(kBlue);
  expectedLikelihoodSL->SetLineWidth(2);
  expectedLikelihoodSL->Draw("Csame");
  expectedLikelihoodSLNoCorr->SetLineColor(kRed+1);
  expectedLikelihoodSLNoCorr->SetLineWidth(2);
  expectedLikelihoodSLNoCorr->Draw("Csame");
  expectedLikelihoodSLSR->SetLineColor(kGreen+1);
  expectedLikelihoodSLSR->SetLineWidth(2);
  expectedLikelihoodSLSR->Draw("Csame");
  expectedLikelihoodSLSRNoCorr->SetLineColor(kOrange+1);
  expectedLikelihoodSLSRNoCorr->SetLineWidth(2);
  expectedLikelihoodSLSRNoCorr->Draw("Csame");
  
  leg.Clear();
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(expectedLikelihoodComb,"Expected full Likelihood","L");
  leg.AddEntry(expectedLikelihoodSLFullCorr,"Expected Simp. Likelihood full correlation","L");
  leg.AddEntry(expectedLikelihoodSL,"Expected Simp. Likelihood","L");
  leg.AddEntry(expectedLikelihoodSLNoCorr,"Expected Simp. Likelihood w/o corr.","L");
  leg.AddEntry(expectedLikelihoodSLSR,"Expected Fit with SR","L");
  leg.AddEntry(expectedLikelihoodSLSRNoCorr,"Expected Fit with SR w/o corr.","L");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/likelihoodComparison_expected.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/likelihoodComparison_expected.pdf").c_str(),"pdf");

  TFile* outputFile = new TFile((outputDIR+"/likelihoodComparison.root").c_str(),"RECREATE");
  outputFile->cd();
  observedLikelihoodComb->Write("observed_full_likelihood");
  expectedLikelihoodComb->Write("expected_full_likelihood");
  observedLikelihoodSLFullCorr->Write("observed_simp_likelihood_full_covariance");
  expectedLikelihoodSLFullCorr->Write("expected_simp_likelihood_full_covariance");
  observedLikelihoodSL->Write("observed_simp_likelihood");
  expectedLikelihoodSL->Write("expected_simp_likelihood");
  observedLikelihoodSLNoCorr->Write("observed_simp_likelihood_no_corr");
  expectedLikelihoodSLNoCorr->Write("expected_simp_likelihood_no_corr");
  observedLikelihoodSLSR->Write("observed_simp_likelihood_SR");
  expectedLikelihoodSLSR->Write("expected_simp_likelihood_SR");
  observedLikelihoodSLSRNoCorr->Write("observed_simp_likelihood_SR_no_corr");
  expectedLikelihoodSLSRNoCorr->Write("expected_simp_likelihood_SR_no_corr");
  outputFile->Close();
}

//  LocalWords:  maxMuExpComb
