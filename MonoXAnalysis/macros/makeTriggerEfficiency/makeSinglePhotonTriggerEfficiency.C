#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

// boson pt and recoil when using single photon
vector<float> bins_singlePhoton            = {160,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,180,190,200,210,220,235,250,300};
vector<float> bins_singlePhoton_vbf        = {160,165,170,175,180,190,200,225,250,300};
vector<float> bins_singlePhoton_recoil     = {100,110,120,130,140,150,160,170,180,190,200,210,220,230,250,300,350,450,550,650,850,1000,1200,1450};
vector<float> bins_singlePhoton_vbf_recoil = {100,125,150,175,200,225,250,350,450,650,1000,1450};

// boson pt and recoil when using jetHT data
vector<float> bins_jetHT            = {160,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,180,190,200,250,300,350,400,450,500,550,600,650,700,800,1000,1200,1400};
vector<float> bins_jetHT_vbf        = {160,165,170,175,180,190,200,225,250,300};
vector<float> bins_jetHT_recoil     = {100,110,120,130,140,150,160,170,180,190,200,210,220,230,250,300,350,450,550,650,850,1000,1200,1450};
vector<float> bins_jetHT_vbf_recoil = {100,125,150,175,200,225,250,350,450,650,1000,1450};

vector<string> RunEra = {"Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G","Run2016H"};

static float leadingJetVBF     = 80;
static float trailingJetVBF    = 40;
static float detajj            = 3.5;
static float mjj               = 1000;
static float jetmetdphi        = 0.5;
static float dphijj            = 1.5;
static float recoilSelection   = 150;
static float photonPtSelection = 120;
static bool  drawFitFunctions_ = false;
static float lumi_             = 12.9;

void makePlot(TCanvas* canvas, TEfficiency* efficiency, TF1* funz, TEfficiency* efficiency_recover, TF1* funz_recover, string labelX, string ouputDIR, string postfix){

  TGraphAsymmErrors* eff = efficiency->CreateGraph();
  TGraphAsymmErrors* eff_recover = NULL;
  if(efficiency_recover != NULL)
    eff_recover = efficiency_recover->CreateGraph();
  
  eff->GetXaxis()->SetTitle(labelX.c_str());
  eff->GetYaxis()->SetTitle("Trigger Efficiency");
  eff->GetYaxis()->CenterTitle();
  eff->GetXaxis()->SetTitleOffset(1.0);
  
  canvas->SetRightMargin(0.075);
  canvas->SetTopMargin(0.06);
  canvas->cd();
  eff->Draw("AP");
  if(eff_recover != NULL)
    eff_recover->Draw("PSAME");
  
  if(drawFitFunctions_){
    eff->Fit(funz);
    funz->Draw("SAME");
    if(eff_recover != NULL and funz_recover != NULL){
      eff_recover->Fit(funz_recover);
      funz_recover->Draw("SAME");
    }
  }

  canvas->RedrawAxis();
  CMS_lumi(canvas,string(Form("%.1f",lumi_)),true);


  TLegend* leg = new TLegend(0.25,0.15,0.55,0.35);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  if(eff != NULL)
    leg->AddEntry(eff,"Photon 165/175","PLE");
  if(eff_recover != NULL)
    leg->AddEntry(eff_recover,"Photon or PFHT800 ","PLE");
  leg->Draw("same");
  
  canvas->SaveAs((ouputDIR+"/photonTriggerEfficiency_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((ouputDIR+"/photonTriggerEfficiency_"+postfix+".pdf").c_str(),"pdf");
}


/////
void makeSinglePhotonTriggerEfficiency(string inputDIR, string ouputDIR, float lumi = 12.9, bool useJetHT = false, bool drawFitFunctions =  false, int runCut = 999999999) {
  
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  drawFitFunctions_ = drawFitFunctions;
  lumi_ = lumi;

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+ouputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  // input tree                                                                                                                                                                                       
  TChain* tree = new TChain("tree/tree");
  // use only a subset of directories                                                                                                                                                         
  if(not useJetHT)
    system(("ls "+inputDIR+"  | grep SinglePhoton > list_dir.txt").c_str());
  else
    system(("ls "+inputDIR+"  | grep JetHT > list_dir.txt").c_str());
  ifstream dirlist ("list_dir.txt");
  string dirname;
  if(dirlist.is_open()){
    while(not dirlist.eof()){
      getline(dirlist,dirname);
      bool found = false;
      for(auto era : RunEra){
        if(dirname.find(era) != string::npos)
          found = true;
      }
      if(found == false) continue;
      system(("find "+inputDIR+"/"+dirname+" -name  \"*.root\" > list.txt").c_str());
      ifstream file("list.txt");
      if(file.is_open()){
        string line;
        while(!file.eof()){
          getline(file,line);
          if(TString(line).Contains("failed")) continue;
          if(line == "" or not TString(line).Contains("root")) continue;
	  cout<<line<<endl;
          tree->Add(line.c_str());
        }
      }
      system("rm list.txt");
    }
  }
  system("rm list_dir.txt");

  // select the binning depending on the sample
  vector<float> bins_photonpt;
  vector<float> bins_vbf_photonpt;
  vector<float> bins_recoil;
  vector<float> bins_vbf_recoil;
  if(useJetHT){
    bins_photonpt     = bins_jetHT;
    bins_vbf_photonpt = bins_jetHT_vbf;
    bins_recoil       = bins_jetHT_recoil;
    bins_vbf_recoil   = bins_jetHT_vbf_recoil;
  }
  else{
    bins_photonpt     = bins_singlePhoton;
    bins_vbf_photonpt = bins_singlePhoton_vbf;
    bins_recoil       = bins_singlePhoton_recoil;
    bins_vbf_recoil   = bins_singlePhoton_vbf_recoil;
  }
 
  
  /// makeing turn-on functions
  TF1 *fitfunc_photonpt = new TF1("fitfunc_photonpt", ErfCB, bins_photonpt.front(), bins_photonpt.back(), 5);
  TF1 *fitfunc_recover_photonpt = new TF1("fitfunc_recover_photonpt", ErfCB, bins_photonpt.front(), bins_photonpt.back(), 5);
  TF1 *fitfunc_vbf_photonpt     = new TF1("fitfunc_vbf_photonpt", ErfCB, bins_vbf_photonpt.front(), bins_vbf_photonpt.back(), 5);
  TF1 *fitfunc_recover_vbf_photonpt = new TF1("fitfunc_recover_vbf_photonpt", ErfCB, bins_vbf_photonpt.front(), bins_vbf_photonpt.back(), 5);
  fitfunc_photonpt->SetParameters(165., 5., 5., 4., 1.);
  fitfunc_recover_photonpt->SetParameters(165., 5., 5., 4., 1.);
  fitfunc_vbf_photonpt->SetParameters(165., 5., 5., 4., 1.);
  fitfunc_recover_vbf_photonpt->SetParameters(165., 5., 5., 4., 1.);

  TF1 *fitfunc_recoil  = new TF1("fitfunc_recoil", ErfCB, bins_recoil.front(), bins_recoil.back(), 5);
  TF1 *fitfunc_recover_recoil = new TF1("fitfunc_recover_recoil", ErfCB, bins_recoil.front(), bins_recoil.back(), 5);
  TF1 *fitfunc_vbf_recoil     = new TF1("fitfunc_vbf_recoil", ErfCB, bins_vbf_recoil.front(), bins_vbf_recoil.back(), 5);
  TF1 *fitfunc_recover_vbf_recoil = new TF1("fitfunc_recover_vbf_recoil", ErfCB, bins_vbf_recoil.front(), bins_vbf_recoil.back(), 5);
  fitfunc_recoil->SetParameters(165., 5., 5., 4., 1.);
  fitfunc_recover_recoil->SetParameters(165., 5., 5., 4., 1.);
  fitfunc_vbf_recoil->SetParameters(165., 5., 5., 4., 1.);
  fitfunc_recover_vbf_recoil->SetParameters(165., 5., 5., 4., 1.);
  
   
  // efficiency vs photon pt
  TH1F* hnum_photonpt   = new TH1F("hnum_photonpt", "",   bins_photonpt.size()-1, &bins_photonpt[0]);
  TH1F* hden_photonpt   = new TH1F("hden_photonpt", "",   bins_photonpt.size()-1, &bins_photonpt[0]);
  TH1F* hnum_recover_photonpt = new TH1F("hnum_recover_photonpt", "", bins_photonpt.size()-1, &bins_photonpt[0]);
  TH1F* hden_recover_photonpt = new TH1F("hden_recover_photonpt", "", bins_photonpt.size()-1, &bins_photonpt[0]);
  hnum_photonpt->Sumw2();
  hden_photonpt->Sumw2();
  hnum_recover_photonpt->Sumw2();
  hden_recover_photonpt->Sumw2();

  TH1F* hnum_vbf_photonpt   = new TH1F("hnum_vbf_photonpt", "",   bins_vbf_photonpt.size()-1, &bins_vbf_photonpt[0]);
  TH1F* hden_vbf_photonpt   = new TH1F("hden_vbf_photonpt", "",   bins_vbf_photonpt.size()-1, &bins_vbf_photonpt[0]);
  TH1F* hnum_recover_vbf_photonpt = new TH1F("hnum_recover_vbf_photonpt", "", bins_vbf_photonpt.size()-1, &bins_vbf_photonpt[0]);
  TH1F* hden_recover_vbf_photonpt = new TH1F("hden_recover_vbf_photonpt", "", bins_vbf_photonpt.size()-1, &bins_vbf_photonpt[0]);
  hnum_vbf_photonpt->Sumw2();
  hden_vbf_photonpt->Sumw2();
  hnum_recover_vbf_photonpt->Sumw2();
  hden_recover_vbf_photonpt->Sumw2();

  // efficiency vs recoil
  TH1F* hnum_recoil   = new TH1F("hnum_recoil", "",   bins_recoil.size()-1, &bins_recoil[0]);
  TH1F* hden_recoil   = new TH1F("hden_recoil", "",   bins_recoil.size()-1, &bins_recoil[0]);
  TH1F* hnum_recover_recoil = new TH1F("hnum_recover_recoil", "", bins_recoil.size()-1, &bins_recoil[0]);
  TH1F* hden_recover_recoil = new TH1F("hden_recover_recoil", "", bins_recoil.size()-1, &bins_recoil[0]);
  hnum_recoil->Sumw2();
  hden_recoil->Sumw2();
  hnum_recover_recoil->Sumw2();
  hden_recover_recoil->Sumw2();

  TH1F* hnum_vbf_recoil   = new TH1F("hnum_vbf_recoil", "",   bins_vbf_recoil.size()-1, &bins_vbf_recoil[0]);
  TH1F* hden_vbf_recoil   = new TH1F("hden_vbf_recoil", "",   bins_vbf_recoil.size()-1, &bins_vbf_recoil[0]);
  TH1F* hnum_recover_vbf_recoil = new TH1F("hnum_recover_vbf_recoil", "", bins_vbf_recoil.size()-1, &bins_vbf_recoil[0]);
  TH1F* hden_recover_vbf_recoil = new TH1F("hden_recover_vbf_recoil", "", bins_vbf_recoil.size()-1, &bins_vbf_recoil[0]);
  hnum_vbf_recoil->Sumw2();
  hden_vbf_recoil->Sumw2();
  hnum_recover_vbf_recoil->Sumw2();
  hden_recover_vbf_recoil->Sumw2();
  
  TTreeReader reader(tree);
  // basic triggers
  TTreeReaderValue<UChar_t> hltp90    (reader,"hltphoton90");
  TTreeReaderValue<UChar_t> hltp90PFHT(reader,"hltphoton90PFHT");
  TTreeReaderValue<UChar_t> hltp120   (reader,"hltphoton120");
  TTreeReaderValue<UChar_t> hltp120vbf(reader,"hltphoton120vbf");
  TTreeReaderValue<UChar_t> hltp165   (reader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175   (reader,"hltphoton175");
  TTreeReaderValue<UChar_t> hltht400  (reader,"hltPFHT400");
  TTreeReaderValue<UChar_t> hltht475  (reader,"hltPFHT475");
  TTreeReaderValue<UChar_t> hltht600  (reader,"hltPFHT600");
  TTreeReaderValue<UChar_t> hltht650  (reader,"hltPFHT650");
  TTreeReaderValue<UChar_t> hltht800  (reader,"hltPFHT800");
  TTreeReaderValue<UChar_t> hltecalht800  (reader,"hltEcalHT800");

  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");

  TTreeReaderValue<unsigned int> ntausraw    (reader,"ntausrawold");
  TTreeReaderValue<unsigned int> nmuons      (reader,"nmuons");
  TTreeReaderValue<unsigned int> nelectrons  (reader,"nelectrons");
  TTreeReaderValue<unsigned int> nphotons    (reader,"nphotons");
  TTreeReaderValue<unsigned int> nincjets    (reader,"njetsinc");
  TTreeReaderValue<unsigned int> nbjets      (reader,"nbjetslowpt");
  TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
  TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
  TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
  TTreeReaderValue<vector<float> > jetchfrac  (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > jetnhfrac  (reader,"combinejetNHfrac");

  TTreeReaderValue<float> phpt   (reader,"phpt");
  TTreeReaderValue<float> pheta  (reader,"pheta");
  TTreeReaderValue<int>   phidm  (reader,"phidm");

  TTreeReaderValue<float> met       (reader,"t1pfmet");
  TTreeReaderValue<float> metphi    (reader,"t1pfmetphi");
  TTreeReaderValue<float> pmet      (reader,"t1phmet");
  TTreeReaderValue<float> pmetphi   (reader,"t1phmetphi");
  TTreeReaderValue<float> metpf     (reader,"pfmet");
  TTreeReaderValue<float> metcalo   (reader,"calomet");
  TTreeReaderValue<float> jpmdphi   (reader,"incjetphmetdphimin4");

  
  cout<<"Total number of events: "<<tree->GetEntries()<<endl;

  long int nEvents = 0;
  while(reader.Next()){
    cout.flush();
    if(nEvents%100000 == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/tree->GetEntries()*100<<" % ";
    nEvents++;

    if(*nmuons     != 0) continue;
    if(*nelectrons != 0) continue;
    if(*nbjets   != 0)   continue;
    if(*ntausraw != 0)   continue;
    if(*nphotons != 1)   continue;
    if(fabs(*pheta) > 1.4442) continue;
    if(*phidm != 1)     continue;
    if(not *fcsc)   continue;
    if(not *fcsct)  continue;
    if(not *feeb)   continue;
    if(not *fetp)   continue;
    if(not *fvtx) continue;
    if(not *fbadmu) continue;
    if(not *fbadch) continue;
    if(not *fhbhe)  continue;
    if(not *fhbiso) continue;
    if(jetpt->size() <=0) continue;
    /// single photon
    if(not useJetHT){      
      if(*hltp90 or *hltp120 or *hltp90PFHT){

	if(jetpt->at(0) > 100 and fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) > 0.1 and jetnhfrac->at(0) < 0.8 and *jpmdphi > 0.5){
	  if(*pmet > recoilSelection)
	    hden_photonpt->Fill(*phpt);
	  if(*phpt > photonPtSelection)
	    hden_recoil->Fill(*pmet);
	  if(*hltp165 or *hltp175){
	    if(*pmet > recoilSelection)
	      hnum_photonpt->Fill(*phpt);
	    if(*phpt > photonPtSelection)
	      hnum_recoil->Fill(*pmet);
	  } 
	}
      
	if(*nincjets >= 2){
	  if(jetpt->at(0) > leadingJetVBF and jetpt->at(1) > trailingJetVBF and fabs(jeteta->at(0)-jeteta->at(1)) > detajj and *jpmdphi > jetmetdphi){
	    
	    if(fabs(jetpt->at(0)) < 2.5 and jetchfrac->at(0) < 0.1) continue;
	    if(fabs(jetpt->at(0)) < 2.5 and jetnhfrac->at(0) > 0.8) continue;
	    if(fabs(jetpt->at(0)) < 3.2 and fabs(jetpt->at(0)) > 3.0 and jetnhfrac->at(0) > 0.96) continue;

	    float deltaPhi = fabs(jetphi->at(0)-jetphi->at(1));
	    if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi()-deltaPhi;
	    if(deltaPhi > dphijj) continue;

 	    TLorentzVector jet1, jet2;
	    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	    if((jet1+jet2).M() < mjj) continue;

	    if(*pmet > recoilSelection)
	      hden_vbf_photonpt->Fill(*phpt);
	    if(*phpt > photonPtSelection)
	      hden_vbf_recoil->Fill(*pmet);
	    
	    if(*hltp165 or *hltp175){
	      if(*pmet > recoilSelection)
		hnum_vbf_photonpt->Fill(*phpt);
	      if(*phpt > photonPtSelection)
		hnum_vbf_recoil->Fill(*pmet);
	    }	
	  }
	}
      }
    }
    // jet ht
    else{
      if(*hltht400 or *hltht475 or *hltht600 or *hltht650){

	if(jetpt->at(0) > 100 and fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) > 0.1 and jetnhfrac->at(0) < 0.8 and *jpmdphi > 0.5){	  
	  if(*pmet > recoilSelection)
	    hden_photonpt->Fill(*phpt);
	  if(*phpt > photonPtSelection)
	    hden_recoil->Fill(*pmet);
	  if(*pmet > recoilSelection)
	    hden_recover_photonpt->Fill(*phpt);
	  if(*phpt > photonPtSelection)
	    hden_recover_recoil->Fill(*pmet);

	  if(*hltp165 or *hltp175){
	    if(*pmet > recoilSelection)
	      hnum_photonpt->Fill(*phpt);
	    if(*phpt > photonPtSelection)
	      hnum_recoil->Fill(*pmet);
	  }
	  if(*hltp165 or *hltp175 or *hltht800){
	    if(*pmet > recoilSelection)
	      hnum_recover_photonpt->Fill(*phpt);	  
	    if(*phpt > photonPtSelection)
	      hnum_recover_recoil->Fill(*pmet);	  
	  }      	
	}

      	if(*nincjets >= 2){
	  if(jetpt->at(0) > leadingJetVBF and jetpt->at(1) > trailingJetVBF and fabs(jeteta->at(0)-jeteta->at(1)) > detajj and *jpmdphi > jetmetdphi){
	    
	    if(fabs(jetpt->at(0)) < 2.5 and jetchfrac->at(0) < 0.1) continue;
	    if(fabs(jetpt->at(0)) < 2.5 and jetnhfrac->at(0) > 0.8) continue;
	    if(fabs(jetpt->at(0)) < 3.2 and fabs(jetpt->at(0)) > 3.0 and jetnhfrac->at(0) > 0.96) continue;

	    float deltaPhi = fabs(jetphi->at(0)-jetphi->at(1));
	    if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi()-deltaPhi;
	    if(deltaPhi > dphijj) continue;

	    TLorentzVector jet1, jet2;
	    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	    if((jet1+jet2).M() < mjj) continue;
	    
	    if(*pmet > recoilSelection)
	      hden_vbf_photonpt->Fill(*phpt);
	    if(*phpt > photonPtSelection)
	      hden_vbf_recoil->Fill(*pmet);
	    if(*pmet > recoilSelection)
	      hden_recover_vbf_photonpt->Fill(*phpt);
	    if(*phpt > photonPtSelection)
	      hden_recover_vbf_recoil->Fill(*pmet);
	    
	    if(*hltp165 or *hltp175){
	      if(*pmet > recoilSelection)
		hnum_vbf_photonpt->Fill(*phpt);
	      if(*phpt > photonPtSelection)
		hnum_vbf_recoil->Fill(*pmet);   
	    }
	    
	    if(*hltp165 or *hltp175 or *hltht800){
	      if(*pmet > recoilSelection)
		hnum_recover_vbf_photonpt->Fill(*phpt);
	      if(*phpt > photonPtSelection)
		hnum_recover_vbf_recoil->Fill(*pmet);
	    }
	  }
	}
      }
    }
  }
  cout<<endl;

  //////////////////////
  TEfficiency* eff_photonpt = new TEfficiency(*hnum_photonpt,*hden_photonpt);  
  eff_photonpt->SetMarkerColor(kBlack);
  eff_photonpt->SetLineColor(kBlack);
  eff_photonpt->SetMarkerStyle(20);
  eff_photonpt->SetMarkerSize(1);
  fitfunc_photonpt->SetLineColor(kBlack);
  fitfunc_photonpt->SetLineWidth(2);

  TEfficiency* eff_vbf_photonpt = new TEfficiency(*hnum_vbf_photonpt,*hden_vbf_photonpt);  
  eff_vbf_photonpt->SetMarkerColor(kBlack);
  eff_vbf_photonpt->SetLineColor(kBlack);
  eff_vbf_photonpt->SetMarkerStyle(20);
  eff_vbf_photonpt->SetMarkerSize(1);
  fitfunc_vbf_photonpt->SetLineColor(kBlack);
  fitfunc_vbf_photonpt->SetLineWidth(2);

  TEfficiency* eff_recoil = new TEfficiency(*hnum_recoil,*hden_recoil);  
  eff_recoil->SetMarkerColor(kBlack);
  eff_recoil->SetLineColor(kBlack);
  eff_recoil->SetMarkerStyle(20);
  eff_recoil->SetMarkerSize(1);
  fitfunc_recoil->SetLineColor(kBlack);
  fitfunc_recoil->SetLineWidth(2);

  TEfficiency* eff_vbf_recoil = new TEfficiency(*hnum_vbf_recoil,*hden_vbf_recoil);  
  eff_vbf_recoil->SetMarkerColor(kBlack);
  eff_vbf_recoil->SetLineColor(kBlack);
  eff_vbf_recoil->SetMarkerStyle(20);
  eff_vbf_recoil->SetMarkerSize(1);
  fitfunc_vbf_recoil->SetLineColor(kBlack);
  fitfunc_vbf_recoil->SetLineWidth(2);

  TEfficiency* eff_recover_photonpt = NULL;
  TEfficiency* eff_recover_vbf_photonpt = NULL;
  TEfficiency* eff_recover_recoil = NULL;
  TEfficiency* eff_recover_vbf_recoil = NULL;
  
  if(not useJetHT){
    makePlot(canvas,eff_photonpt,fitfunc_photonpt,NULL,NULL,"Photon p_{T} [GeV]",ouputDIR,"photonpt");
    makePlot(canvas,eff_vbf_photonpt,fitfunc_vbf_photonpt,NULL,NULL,"Photon p_{T} [GeV] [GeV]",ouputDIR,"vbf_photonpt");
    makePlot(canvas,eff_recoil,fitfunc_recoil,NULL,NULL,"Recoil [GeV]",ouputDIR,"recoil");
    makePlot(canvas,eff_vbf_recoil,fitfunc_vbf_recoil,NULL,NULL,"Recoil [GeV]",ouputDIR,"vbf_recoil");
  }
  else{
        
    eff_recover_photonpt = new TEfficiency(*hnum_recover_photonpt,*hden_recover_photonpt);  
    eff_recover_photonpt->SetMarkerColor(kRed);
    eff_recover_photonpt->SetLineColor(kRed);
    eff_recover_photonpt->SetMarkerStyle(20);
    eff_recover_photonpt->SetMarkerSize(1);
    fitfunc_recover_photonpt->SetLineColor(kRed);
    fitfunc_recover_photonpt->SetLineWidth(2);

    eff_recover_vbf_photonpt = new TEfficiency(*hnum_recover_vbf_photonpt,*hden_recover_vbf_photonpt);  
    eff_recover_vbf_photonpt->SetMarkerColor(kRed);
    eff_recover_vbf_photonpt->SetLineColor(kRed);
    eff_recover_vbf_photonpt->SetMarkerStyle(20);
    eff_recover_vbf_photonpt->SetMarkerSize(1);
    fitfunc_recover_vbf_photonpt->SetLineColor(kRed);
    fitfunc_recover_vbf_photonpt->SetLineWidth(2);

    eff_recover_recoil = new TEfficiency(*hnum_recover_recoil,*hden_recover_recoil);  
    eff_recover_recoil->SetMarkerColor(kRed);
    eff_recover_recoil->SetLineColor(kRed);
    eff_recover_recoil->SetMarkerStyle(20);
    eff_recover_recoil->SetMarkerSize(1);
    fitfunc_recover_recoil->SetLineColor(kRed);
    fitfunc_recover_recoil->SetLineWidth(2);

    eff_recover_vbf_recoil = new TEfficiency(*hnum_recover_vbf_recoil,*hden_recover_vbf_recoil);  
    eff_recover_vbf_recoil->SetMarkerColor(kRed);
    eff_recover_vbf_recoil->SetLineColor(kRed);
    eff_recover_vbf_recoil->SetMarkerStyle(20);
    eff_recover_vbf_recoil->SetMarkerSize(1);
    fitfunc_recover_vbf_recoil->SetLineColor(kRed);
    fitfunc_recover_vbf_recoil->SetLineWidth(2);


    makePlot(canvas,eff_photonpt,fitfunc_photonpt,eff_recover_photonpt,fitfunc_recover_photonpt,"Photon p_{T} [GeV]",ouputDIR,"photonpt");
    makePlot(canvas,eff_vbf_photonpt,fitfunc_vbf_photonpt,eff_recover_vbf_photonpt,fitfunc_recover_vbf_photonpt,"Photon p_{T} [GeV] [GeV]",ouputDIR,"vbf_photonpt");
    makePlot(canvas,eff_recoil,fitfunc_recoil,eff_recover_recoil,fitfunc_recover_recoil,"Recoil [GeV]",ouputDIR,"recoil");
    makePlot(canvas,eff_vbf_recoil,fitfunc_vbf_recoil,eff_recover_vbf_recoil,fitfunc_recover_vbf_recoil,"Recoil [GeV]",ouputDIR,"vbf_recoil");
  }

  TFile* outputFile = new TFile((ouputDIR+"/photonTriggerEfficiency.root").c_str(),"RECREATE");
  outputFile->cd();
  if(eff_photonpt)
    eff_photonpt->Write("eff_photonpt");
  if(eff_vbf_photonpt)
    eff_vbf_photonpt->Write("eff_vbf_photonpt");
  if(eff_recoil)
    eff_recoil->Write("eff_recoil");
  if(eff_vbf_recoil)
    eff_vbf_recoil->Write("eff_vbf_recoil");
  if(eff_recover_photonpt)
    eff_recover_photonpt->Write("eff_recover_photonpt");
  if(eff_recover_vbf_photonpt)
    eff_recover_vbf_photonpt->Write("eff_recover_vbf_photonpt");
  if(eff_recover_recoil)
    eff_recover_recoil->Write("eff_recover_recoil");
  if(eff_recover_vbf_recoil)
    eff_recover_vbf_recoil->Write("eff_recover_vbf_recoil");


}

