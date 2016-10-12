#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

vector<float> bins_singlePhoton = {160,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,180,190,200,210,220,235,250};
vector<float> bins_jetHT = {175,180,190,200,250,300,350,400,450,500,550,600,650,700,800,1000,1200,1400};

//vector<string> RunEra = {"Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G"};                                                                                                        
vector<string> RunEra = {"Run2016B","Run2016C","Run2016D"};

void makeSinglePhotonTriggerEfficiency(string inputDIR, string ouputDIR, float lumi = 0.86, bool useJetHT = false, bool drawFitFunctions =  false, int runCut = 999999999) {

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+ouputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  TF1 *fitfunc = new TF1("fitfunc", ErfCB, 150, 300, 5);
  fitfunc->SetParameters(165., 5., 5., 4., 1.);

  TF1 *fitfunc_recover = new TF1("fitfunc_recover", ErfCB, 150, 300, 5);
  fitfunc_recover->SetParameters(165., 5., 5., 4., 1.);

  vector<float> bins;
  if(useJetHT)
    bins = bins_jetHT;
  else
    bins = bins_singlePhoton;  
  fitfunc->SetRange(bins.front(),bins.back());
  
  TChain* tree = new TChain("tree/tree");
  // should use the wmnu events triggered by single muon
  for(auto era : RunEra)
    tree->Add((inputDIR+"/*"+era+"*/*root").c_str());
    
  TH1F* hnum   = new TH1F("hnum", "",   bins.size()-1, &bins[0]);
  TH1F* hden   = new TH1F("hden", "",   bins.size()-1, &bins[0]);
  TH1F* hnum_recover = new TH1F("hnum_recover", "", bins.size()-1, &bins[0]);
  TH1F* hden_recover = new TH1F("hden_recover", "", bins.size()-1, &bins[0]);
  hnum->Sumw2();
  hden->Sumw2();
  hnum_recover->Sumw2();
  hden_recover->Sumw2();

  TTreeReader reader(tree);
  TTreeReaderValue<UChar_t> hltp50    (reader,"hltphoton50");
  TTreeReaderValue<UChar_t> hltp75    (reader,"hltphoton75");
  TTreeReaderValue<UChar_t> hltp90    (reader,"hltphoton90");
  TTreeReaderValue<UChar_t> hltp120   (reader,"hltphoton120");
  TTreeReaderValue<UChar_t> hltp165   (reader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175   (reader,"hltphoton175");

  TTreeReaderValue<UChar_t> hltht400    (reader,"hltPFHT400");
  TTreeReaderValue<UChar_t> hltht475    (reader,"hltPFHT475");
  TTreeReaderValue<UChar_t> hltht600    (reader,"hltPFHT600");
  TTreeReaderValue<UChar_t> hltht650    (reader,"hltPFHT650");
  TTreeReaderValue<UChar_t> hltht800   (reader,"hltPFHT800");

  TTreeReaderValue<double>  phpt     (reader,"phpt");
  TTreeReaderValue<double>  pheta    (reader,"pheta");
  TTreeReaderValue<int>     phidm    (reader,"phidm");

  TTreeReaderValue<int>     nphotons (reader,"nphotons");
  TTreeReaderValue<int>     nmuons   (reader,"nmuons");
  TTreeReaderValue<int>     nelectrons (reader,"nelectrons");
  TTreeReaderValue<int>     nbjets   (reader,"nbjetslowpt");
  TTreeReaderValue<int>     ntaus    (reader,"ntausraw");

  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");

  cout<<"Total number of events: "<<tree->GetEntries()<<endl;
  long int nEvents = 0;
  while(reader.Next()){

    cout.flush();
    if(nEvents%100000) cout<<" Analyzing events "<<double(nEvents)/tree->GetEntries()*100<<" % ";
    nEvents++;
    if(*nmuons != 0)   continue;
    if(*nelectrons != 0) continue;
    if(*nbjets != 0)   continue;
    if(*ntaus != 0)    continue;
    if(*nphotons == 1) continue;
    if(fabs(*pheta) > 1.4442) continue;
    if(*phidm != 1) continue;
    if(not *fcsc or not *fcsct or not *feeb or not *fetp or not *fvtx or not *fbadmu or not *fbadch) continue;

    if(not useJetHT){
      if(*hltp50 or *hltp75 or *hltp90 or *hltp120){
	hden->Fill(*phpt);
	if(*hltp165 or *hltp175)
	  hnum->Fill(*phpt);
      }
    }
    else{
      if(*hltht400 or *hltht475 or *hltht600 or *hltht650){
	hden->Fill(*phpt);
	hden_recover->Fill(*phpt);
	if(*hltp165 or *hltp175)
	  hnum->Fill(*phpt);
	if(*hltp165 or *hltp175 or *hltht800)
	  hnum_recover->Fill(*phpt);	  
      }      	
    }
  }
  cout<<endl;

  TEfficiency* eff = new TEfficiency(*hnum,*hden);  
  eff->SetMarkerColor(kBlack);
  eff->SetLineColor(kBlack);
  eff->SetMarkerStyle(20);
  eff->SetMarkerSize(1);
  fitfunc->SetLineColor(kBlack);
  fitfunc->SetLineWidth(2);

  
  TH1* frame = canvas->DrawFrame(bins.front(), 0.1, bins.back(), 1.1, "");
  frame->GetXaxis()->SetTitle("Photon p_{T} [GeV]");
  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetYaxis()->SetLabelSize(0.8*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetTitleSize(0.8*frame->GetYaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleOffset(1.0);
  
  canvas->SetRightMargin(0.075);
  canvas->SetTopMargin(0.06);
  canvas->Draw();
  canvas->cd();
  frame->Draw();
  eff->Draw("E1PSAME");
  if(drawFitFunctions){
    eff->Fit(fitfunc);
    fitfunc->Draw("SAME");
  }
  canvas->RedrawAxis();
  CMS_lumi(canvas,string(Form("%.2f",lumi)),true);

  TEfficiency* eff_recover = NULL;

  if(useJetHT){

    eff_recover = new TEfficiency(*hnum_recover,*hden_recover);  
    eff_recover->SetMarkerColor(kRed);
    eff_recover->SetLineColor(kRed);
    eff_recover->SetMarkerStyle(20);
    eff_recover->SetMarkerSize(1);
    fitfunc_recover->SetLineColor(kRed);
    fitfunc_recover->SetLineWidth(2);
    
    eff_recover->Draw("E1PSAME");
    if(drawFitFunctions){
      eff_recover->Fit(fitfunc_recover);
      fitfunc_recover->Draw("SAME");
    }
    canvas->RedrawAxis();
    
    TLegend* leg = new TLegend(0.35,0.25,0.9,0.45);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->AddEntry(eff,"Photon 165/175","PLE");
    leg->AddEntry(eff_recover,"Photon || PFHT800 ","PLE");
    leg->Draw("same");
  }
 
  canvas->SaveAs((ouputDIR+"/photonTriggerEff.png").c_str(),"png");
  canvas->SaveAs((ouputDIR+"/photonTriggerEff.pdf").c_str(),"pdf");

  TFile* outputFile = new TFile((ouputDIR+"/photonTriggerEfficiency.root").c_str(),"RECREATE");
  outputFile->cd();
  eff->Write("efficiency_photon");
  if(eff_recover and eff_recover != NULL)
    eff_recover->Write("efficiency_photon_recover");
  fitfunc->Write("efficiency_func");
  if(fitfunc_recover and fitfunc_recover != NULL)
    fitfunc_recover->Write("efficiency_func_recover");
  
}

