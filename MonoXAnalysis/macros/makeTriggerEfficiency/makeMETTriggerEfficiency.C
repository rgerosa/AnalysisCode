#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

vector <float> bins_monojet_muon = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 180., 200., 225, 250., 275., 300., 350., 400., 450., 500., 550., 650., 800., 1000.};
vector <float> bins_monojet_elec = {0.,50.,70.,90.,100,110.,120.,140.,160., 180., 200., 250., 350.,450,550,650,1000.};

vector <float> bins_vbf_muon = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 180., 200., 225, 250., 275., 300., 350., 400., 450., 500., 550., 650., 800., 1000.};
vector <float> bins_vbf_elec = {0.,50.,70.,90.,100,110.,120.,140.,160., 180., 200., 250., 350.,450,550,650,1000.};

void makeMETTriggerEfficiency(string inputDIR, string ouputDIR, float lumi = 0.81, bool singleMuon = true) {

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+ouputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();

  // fitting function for the turn-on 
  TF1 *fitfunc_monojet = new TF1("fitfunc_monojet", ErfCB, 50, 1000, 5);
  fitfunc_monojet->SetParameters(117., 25., 30., 4., 1.);
  TF1 *fitfunc_vbf = new TF1("fitfunc_vbf", ErfCB, 50, 1000, 5);
  fitfunc_vbf->SetParameters(117., 25., 30., 4., 1.);

  // input tree
  TChain* tree = new TChain("tree/tree");
  tree->Add((inputDIR+"/*root").c_str());

  vector<float> bins_monojet;
  vector<float> bins_vbf;
  if(singleMuon){
    bins_monojet = bins_monojet_muon;
    bins_vbf = bins_vbf_muon;
  }
  else{
    bins_monojet = bins_monojet_muon;
    bins_vbf = bins_vbf_muon;
  }

  TH1F* hnum_monojet = new TH1F("hnum_monojet", "", bins_monojet.size()-1, &bins_monojet[0]);
  TH1F* hden_monojet = new TH1F("hden_monojet", "", bins_monojet.size()-1, &bins_monojet[0]);
  hnum_monojet->Sumw2();
  hden_monojet->Sumw2();
  TH1F* hnum_vbf = new TH1F("hnum_vbf", "", bins_vbf.size()-1, &bins_vbf[0]);
  TH1F* hden_vbf = new TH1F("hden_vbf", "", bins_vbf.size()-1, &bins_vbf[0]);
  hnum_vbf->Sumw2();
  hden_vbf->Sumw2();
  
  // define numerator as event with tight muon + trigger requirement
  // define denominator as an event with a tight muon passing single muon trigger
  TTreeReader reader(tree);
  TTreeReaderValue<UChar_t> hltm90     (reader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100    (reader,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110    (reader,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120    (reader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120  (reader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170  (reader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300  (reader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90   (reader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hltsinglemu (reader,"hltsinglemu");
  TTreeReaderValue<UChar_t> hlte      (reader,"hltsingleel");
  TTreeReaderValue<double>  mu1pt     (reader,"mu1pt");
  TTreeReaderValue<double>  mu1eta    (reader,"mu1eta");
  TTreeReaderValue<int>     mu1id     (reader,"mu1id");
  TTreeReaderValue<double>  el1pt     (reader,"el1pt");
  TTreeReaderValue<double>  el1eta    (reader,"el1eta");
  TTreeReaderValue<int>     el1id     (reader,"el1id");

  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");                                                                                                                                       
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");                                                                                                                                       

  TTreeReaderValue<unsigned int> ntausraw    (reader,"ntausraw");
  TTreeReaderValue<unsigned int> nmuons      (reader,"nmuons");
  TTreeReaderValue<unsigned int> nelectrons  (reader,"nelectrons");
  TTreeReaderValue<unsigned int> nphotons    (reader,"nphotons");
  TTreeReaderValue<unsigned int> nincjets    (reader,"njetsinc");
  TTreeReaderValue<unsigned int> nbjets      (reader,"nbjetslowpt");
  
  TTreeReaderValue<vector<double> > jetpt   (reader,"combinejetpt");
  TTreeReaderValue<vector<double> > jetphi  (reader,"combinejetphi");
  TTreeReaderValue<vector<double> > jeteta  (reader,"combinejeteta");
  TTreeReaderValue<vector<double> > jetm    (reader,"combinejetm");  
  TTreeReaderValue<vector<double> > jetchfrac  (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<double> > jetnhfrac  (reader,"combinejetNHfrac");
  TTreeReaderValue<double> met         (reader,"t1pfmet");
  TTreeReaderValue<double> metphi      (reader,"t1pfmetphi");
  TTreeReaderValue<double> mmet        (reader,"t1mumet");
  TTreeReaderValue<double> mmetphi     (reader,"t1mumetphi");
  TTreeReaderValue<double> emet        (reader,"t1elmet");
  TTreeReaderValue<double> emetphi     (reader,"t1elmetphi");
  
  TTreeReaderValue<double> metpf       (reader,"pfmet");
  TTreeReaderValue<double> metcalo     (reader,"calomet");
  TTreeReaderValue<double> jmmdphi (reader,"incjetmumetdphimin4");
  TTreeReaderValue<double> jemdphi (reader,"incjetelmetdphimin4");

  //////////////////
  while(reader.Next()){
    // define the denominator
    if(*nbjets   != 0) continue;
    if(*ntausraw != 0) continue;
    if(*nphotons  != 0) continue;
    if(not *fcsc or not *fcsct or not *feeb or not *fetp or not *fvtx or not *fbadmu or not *fbadch) continue;

    if(singleMuon){
      if(not *hltsinglemu) continue;    
      if(*mu1pt < 20) continue;
      if(fabs(*mu1eta) > 2.4) continue;
      if(*mu1id != 1) continue;
      if(*nelectrons > 0 ) continue;
      if(fabs(*mmet-*metcalo)/(*metpf) > 0.5) continue;
      if(*jmmdphi < 0.5) continue;      
      // denominator monojet
      if(jetpt->at(0) > 100 and fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) > 0.1 and jetnhfrac->at(0) < 0.8){
	hden_monojet->Fill(*mmet);
	// numerator monojet
	if((*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300))
	  hnum_monojet->Fill(*mmet);	
      }

      // denominator VBF
      if(*nincjets > 1 and jetpt->at(0) > 70 and jetpt->at(1) > 50 and fabs(jeteta->at(0)-jeteta->at(1)) > 2 and 
	 *jmmdphi > 1 and jeteta->at(0)*jetpt->at(1) < 0){
	TLorentzVector jet1, jet2;
	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	if((jet1+jet2).M() < 450) continue;	
	hden_vbf->Fill(*mmet);
	if((*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300))
	  hnum_vbf->Fill(*mmet);	

      }
    }
    else{// electron case
      if(not *hlte) continue;    
      if(*el1pt < 40) continue;
      if(fabs(*el1eta) > 2.5) continue;
      if(*el1id != 1) continue;
      if(*nmuons > 0 ) continue;
      if(fabs(*emet-*metcalo)/(*metpf) > 0.5) continue;
      if(*jemdphi < 0.5) continue;
      // denominator monojet
      if(jetpt->at(0) > 100 and fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) > 0.1 and jetnhfrac->at(0) < 0.8){
	hden_monojet->Fill(*met);
	// numerator monojet
	if((*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300))
	  hnum_monojet->Fill(*met);	
      }

      // denominator VBF
      if(*nincjets > 1 and jetpt->at(0) > 70 and jetpt->at(1) > 50 and fabs(jeteta->at(0)-jeteta->at(1)) > 2 and 
	 *jemdphi > 1 and jeteta->at(0)*jetpt->at(1) < 0){
	TLorentzVector jet1, jet2;
	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	if((jet1+jet2).M() < 450) continue;	
	hden_vbf->Fill(*met);
	if((*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300))
	  hnum_vbf->Fill(*met);	
      }      
    }
  }
  

  TEfficiency* eff_monojet = new TEfficiency(*hnum_monojet,*hden_monojet);
  eff_monojet->Fit(fitfunc_monojet);
  eff_monojet->SetMarkerColor(kBlack);
  eff_monojet->SetLineColor(kBlack);
  eff_monojet->SetMarkerStyle(20);
  eff_monojet->SetMarkerSize(1);
  fitfunc_monojet->SetLineColor(kBlack);
  fitfunc_monojet->SetLineWidth(2);
  
  TEfficiency* eff_vbf = new TEfficiency(*hnum_vbf,*hden_vbf);
  eff_vbf->Fit(fitfunc_vbf);
  eff_vbf->SetMarkerColor(kRed);
  eff_vbf->SetLineColor(kRed);
  eff_vbf->SetMarkerStyle(20);
  eff_vbf->SetMarkerSize(1);
  fitfunc_vbf->SetLineColor(kRed);
  fitfunc_vbf->SetLineWidth(2);
  
  TH1* frame = canvas->DrawFrame(bins_monojet.front(),0.,bins_monojet.back(), 1.1, "");
  if(singleMuon)
    frame->GetXaxis()->SetTitle("E_{T#mu}^{miss} [GeV]");
  else
    frame->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  
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
  eff_monojet->Draw("E1PSAME");
  eff_vbf->Draw("E1PSAME");
  fitfunc_monojet->Draw("SAME");
  fitfunc_vbf->Draw("SAME");
  canvas->RedrawAxis();
  CMS_lumi(canvas,string(Form("%.2f",lumi)),true);

  if(singleMuon){
    canvas->SaveAs((ouputDIR+"/metTriggerEff.png").c_str(),"png");
    canvas->SaveAs((ouputDIR+"/metTriggerEff.pdf").c_str(),"pdf");
    TFile* outputFile = new TFile((ouputDIR+"/metTriggerEfficiency.root").c_str(),"RECREATE");
    outputFile->cd();
    eff_monojet->Write("efficiency_monojet");
    fitfunc_monojet->Write("efficiency_monojet_func");
    eff_vbf->Write("efficiency_vbf");
    fitfunc_vbf->Write("efficiency_vbf_func");
  }
  else{
    canvas->SaveAs((ouputDIR+"/metTriggerEff_el.png").c_str(),"png");
    canvas->SaveAs((ouputDIR+"/metTriggerEff_el.pdf").c_str(),"pdf");
    TFile* outputFile = new TFile((ouputDIR+"/metTriggerEfficiency_el.root").c_str(),"RECREATE");
    outputFile->cd();
    eff_monojet->Write("efficiency_monojet");
    fitfunc_monojet->Write("efficiency_monojet_func");
    eff_vbf->Write("efficiency_vbf");
    fitfunc_vbf->Write("efficiency_vbf_func");

  }
}

