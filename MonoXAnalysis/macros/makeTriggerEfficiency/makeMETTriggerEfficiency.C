#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

// recoil binning
vector <float> bins_monojet_recoil = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 225., 250., 275., 300., 325., 350., 400., 450., 500., 550., 650., 800., 1000., 1250};
vector <float> bins_vbf_recoil     = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 180., 200., 250., 300., 350., 400., 450., 500., 550., 650., 800., 1000., 1500};
vector <float> bins_monoV_recoil   = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 180., 200., 250., 300., 350., 400., 450., 500., 550., 650., 800., 1000., 1500};
// mjj
vector <float> bins_vbf_recoilvsmjj = {0.,50.,75.,100.,125.,150.,175.,200., 225, 250., 300., 400., 550., 800., 1500};
vector <float> bins_vbf_mjj         = {0.,800.,1200.,1700.,3000.};
// detajj
vector <float> bins_vbf_recoilvsdetajj = {0.,50.,75.,100.,125.,150.,175.,200., 225, 250., 300., 400., 550., 800., 1500};
vector <float> bins_vbf_detajj  = {0.,1.5,3.0,5.0,9};

// eras
vector<string> RunEra = {"Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G","Run2016H"};

static float leadingVBF  = 80;
static float trailingVBF = 40;
static float detajj      = 3.5; 
static float mjj         = 1000; 
static float jetmetdphi  = 0.5; 
static float dphijj      = 1.5;
static float recoil      = 200;
static bool  drawUncertaintyBand = false;

void plotTurnOn(TCanvas* canvas, TEfficiency* eff, TF1* fitfunc, const string & axisLabel, const TString & postfix, const string & ouputDIR, const float  & luminosity, 
		const bool & singleMuon, const TString & banner = "");

void GetConfidenceIntervals(TF1* funz, TH1F* obj, Double_t cl, const TFitResultPtr & fitResult);

void makeMETTriggerEfficiency(string inputDIR, string ouputDIR, float luminosity = 0.81, bool singleMuon = true, bool doubleMuon = false) {

  if(singleMuon and doubleMuon){
    cerr<<"Problem with the options --> decide to compute turn on in Wmn or Zmm events "<<endl;
    return;
  }
    

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+ouputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();

  // input tree
  TChain* tree = new TChain("tree/tree");
  // use only a subset of directories
  if(singleMuon or doubleMuon)
    system(("ls "+inputDIR+"  | grep SingleMu > list_dir.txt").c_str());
  else
    system(("ls "+inputDIR+"  | grep SingleEle > list_dir.txt").c_str());
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
	  cout<<"adding following file: "<<line<<endl;
	  tree->Add(line.c_str());
	}
      }
      system("rm list.txt");
    }
  }
  system("rm list_dir.txt");

  // fitting function for the turn-on  as a function of recoil
  TF1 *fitfunc_monojet_recoil = new TF1("fitfunc_monojet_recoil",ErfCB,bins_monojet_recoil.front(), bins_monojet_recoil.back(),5);
  fitfunc_monojet_recoil->SetParameters(120., 25., 30., 4., 1.);
  TF1 *fitfunc_monoV_recoil   = new TF1("fitfunc_monoV_recoil", ErfCB, bins_monoV_recoil.front(), bins_monoV_recoil.back(), 5);
  fitfunc_monoV_recoil->SetParameters(120., 25., 30., 4., 1.);
  TF1* fitfunc_vbf_recoil = new TF1("fitfunc_vbf_recoil", ErfCB, bins_vbf_recoil.front(), bins_vbf_recoil.back(), 5);
  fitfunc_vbf_recoil->SetParameters(120., 25., 30., 4., 1.);
  //////
  vector<TF1*> fitfunc_vbf_mjj;
  for(size_t ibin = 0; ibin < bins_vbf_mjj.size()-1; ibin++){
    fitfunc_vbf_mjj.push_back(new TF1(Form("fitfunc_vbf_mjj_%1f_%1f",bins_vbf_mjj.at(ibin),bins_vbf_mjj.at(ibin+1)), ErfCB,bins_vbf_recoilvsmjj.front(),bins_vbf_recoilvsmjj.back(), 5));
    fitfunc_vbf_mjj.back()->SetParameters(120., 25., 30., 4., 1.);
  }
  ///
  vector<TF1*> fitfunc_vbf_detajj;
  for(size_t ibin = 0; ibin < bins_vbf_detajj.size()-1; ibin++){
    fitfunc_vbf_detajj.push_back(new TF1(Form("fitfunc_vbf_detajj_%1f_%1f",bins_vbf_detajj.at(ibin),bins_vbf_detajj.at(ibin+1)), ErfCB, bins_vbf_recoilvsdetajj.front(),bins_vbf_recoilvsdetajj.back(), 5));
    fitfunc_vbf_detajj.back()->SetParameters(120., 25., 30., 4., 1.);
  }
  
  TH1F* hnum_monojet_recoil = new TH1F("hnum_monojet_recoil", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  TH1F* hden_monojet_recoil = new TH1F("hden_monojet_recoil", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  hnum_monojet_recoil->Sumw2();
  hden_monojet_recoil->Sumw2();
  TH1F* hnum_monoV_recoil = new TH1F("hnum_monoV_recoil", "", bins_monoV_recoil.size()-1, &bins_monoV_recoil[0]);
  TH1F* hden_monoV_recoil = new TH1F("hden_monoV_recoil", "", bins_monoV_recoil.size()-1, &bins_monoV_recoil[0]);
  hnum_monoV_recoil->Sumw2();
  hden_monoV_recoil->Sumw2();
  TH1F* hnum_vbf_recoil = new TH1F("hnum_vbf_recoil", "", bins_vbf_recoil.size()-1, &bins_vbf_recoil[0]);
  TH1F* hden_vbf_recoil = new TH1F("hden_vbf_recoil", "", bins_vbf_recoil.size()-1, &bins_vbf_recoil[0]);
  hnum_vbf_recoil->Sumw2();
  hden_vbf_recoil->Sumw2();

  vector<TH1F*> hnum_vbf_mjj;
  vector<TH1F*> hden_vbf_mjj;
  for(size_t ibin = 0; ibin < bins_vbf_mjj.size()-1; ibin++){
    hnum_vbf_mjj.push_back(new TH1F(Form("hnum_vbf_mjj_%1f_%1f",bins_vbf_mjj.at(ibin),bins_vbf_mjj.at(ibin+1)),"",bins_vbf_recoilvsmjj.size()-1, &bins_vbf_recoilvsmjj[0]));
    hden_vbf_mjj.push_back(new TH1F(Form("hden_vbf_mjj_%1f_%1f",bins_vbf_mjj.at(ibin),bins_vbf_mjj.at(ibin+1)),"",bins_vbf_recoilvsmjj.size()-1, &bins_vbf_recoilvsmjj[0]));
    hnum_vbf_mjj.back()->Sumw2();
    hden_vbf_mjj.back()->Sumw2();
  }

  vector<TH1F*> hnum_vbf_detajj;
  vector<TH1F*> hden_vbf_detajj;
  for(size_t ibin = 0; ibin < bins_vbf_detajj.size()-1; ibin++){
    hnum_vbf_detajj.push_back(new TH1F(Form("hnum_vbf_detajj_%1f_%1f",bins_vbf_detajj.at(ibin),bins_vbf_detajj.at(ibin+1)),"",bins_vbf_recoilvsdetajj.size()-1, &bins_vbf_recoilvsdetajj[0]));
    hden_vbf_detajj.push_back(new TH1F(Form("hden_vbf_detajj_%1f_%1f",bins_vbf_detajj.at(ibin),bins_vbf_detajj.at(ibin+1)),"",bins_vbf_recoilvsdetajj.size()-1, &bins_vbf_recoilvsdetajj[0]));
    hnum_vbf_detajj.back()->Sumw2();
    hden_vbf_detajj.back()->Sumw2();
  }
  
  // define numerator as event with tight muon + trigger requirement
  // define denominator as an event with a tight muon passing single muon trigger
  TTreeReader reader(tree);
  TTreeReaderValue<unsigned int> run    (reader,"run");
  TTreeReaderValue<unsigned int> lumi   (reader,"lumi");
  TTreeReaderValue<unsigned int> event  (reader,"event");
  TTreeReaderValue<UChar_t> hltm90     (reader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100    (reader,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110    (reader,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120    (reader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120  (reader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm100  (reader,"hltmetwithmu100");
  TTreeReaderValue<UChar_t> hltmwm110  (reader,"hltmetwithmu110");
  TTreeReaderValue<UChar_t> hltmwm170  (reader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300  (reader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90   (reader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hltjm      (reader,"hltjetmet");
  TTreeReaderValue<UChar_t> hltsinglemu (reader,"hltsinglemu");
  TTreeReaderValue<UChar_t> hlte      (reader,"hltsingleel");
  TTreeReaderValue<float>   mu1pt     (reader,"mu1pt");
  TTreeReaderValue<float>   mu1eta    (reader,"mu1eta");
  TTreeReaderValue<float>   mu1phi    (reader,"mu1phi");
  TTreeReaderValue<int>     mu1id     (reader,"mu1id");
  TTreeReaderValue<int>     mu1pid     (reader,"mu1pid");
  TTreeReaderValue<float>   mu2pt     (reader,"mu2pt");
  TTreeReaderValue<float>   mu2eta    (reader,"mu2eta");
  TTreeReaderValue<float>   mu2phi    (reader,"mu2phi");
  TTreeReaderValue<int>     mu2id     (reader,"mu2id");
  TTreeReaderValue<int>     mu2pid     (reader,"mu2pid");
  TTreeReaderValue<float>   el1pt     (reader,"el1pt");
  TTreeReaderValue<float>   el1eta    (reader,"el1eta");
  TTreeReaderValue<float>   el1phi    (reader,"el1phi");
  TTreeReaderValue<int>     el1id     (reader,"el1id");
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");                                                                                                                                       
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");                                                                                                                                       
  TTreeReaderValue<unsigned int> ntausraw    (reader,"ntausold");
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
  TTreeReaderValue<float> met         (reader,"t1pfmet");
  TTreeReaderValue<float> metphi      (reader,"t1pfmetphi");
  TTreeReaderValue<float> mmet        (reader,"t1mumet");
  TTreeReaderValue<float> mmetphi     (reader,"t1mumetphi");
  TTreeReaderValue<float> emet        (reader,"t1elmet");
  TTreeReaderValue<float> emetphi     (reader,"t1elmetphi");  
  TTreeReaderValue<float> metpf       (reader,"pfmet");
  TTreeReaderValue<float> metcalo     (reader,"calomet");
  TTreeReaderValue<float> jmmdphi (reader,"incjetmumetdphimin4");
  TTreeReaderValue<float> jemdphi (reader,"incjetelmetdphimin4");
 
  //////////////////
  long int nTotal = tree->GetEntries();
  cout<<"Total number of events: "<<nTotal<<endl;
  long int nEvents = 0;

  TH1D* efficiencyMonojetSelections = new TH1D("efficiencyMonojetSelections","efficiencyMonojetSelections",16,0,17);
  efficiencyMonojetSelections->Sumw2();
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(1,"Total");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(2,"B-veto");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(3,"Tau-veto");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(4,"Photon-veto");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(5,"MET filters");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(6,"HLT single muon");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(7,"Lepton pT and eta");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(8,"Lepton ID");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(9,"Lepton Veto");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(10,"Calo-PF MET");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(11,"Jet met dphi");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(12,"Jet cuts");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(13,"Trigger");

  TH1D* efficiencyMonoVSelections = (TH1D*) efficiencyMonojetSelections->Clone("efficiencyMonoVSelections");
  efficiencyMonoVSelections->GetXaxis()->SetBinLabel(13,"V-tagging");
  efficiencyMonoVSelections->GetXaxis()->SetBinLabel(14,"Trigger");

  TH1D* efficiencyVBFSelections = (TH1D*) efficiencyMonojetSelections->Clone("efficiencyVBFSelections");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(13,"Deta cut");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(14,"Mjj cut");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(15,"Trigger");

  long int nPart = 100000;
  while(reader.Next()){
    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    // define the denominator
    efficiencyMonojetSelections->SetBinContent(1,efficiencyMonojetSelections->GetBinContent(1)+1);
    efficiencyMonoVSelections->SetBinContent(1,efficiencyMonoVSelections->GetBinContent(1)+1);
    efficiencyVBFSelections->SetBinContent(1,efficiencyVBFSelections->GetBinContent(1)+1);
    if(*nbjets   != 0) continue;
    efficiencyMonojetSelections->SetBinContent(2,efficiencyMonojetSelections->GetBinContent(2)+1);
    efficiencyMonoVSelections->SetBinContent(2,efficiencyMonoVSelections->GetBinContent(2)+1);
    efficiencyVBFSelections->SetBinContent(2,efficiencyVBFSelections->GetBinContent(2)+1);
    ///
    if(*ntausraw != 0) continue;
    ///
    efficiencyMonojetSelections->SetBinContent(3,efficiencyMonojetSelections->GetBinContent(3)+1);
    efficiencyMonoVSelections->SetBinContent(3,efficiencyMonoVSelections->GetBinContent(3)+1);
    efficiencyVBFSelections->SetBinContent(3,efficiencyVBFSelections->GetBinContent(3)+1);
    if(*nphotons  != 0) continue;
    efficiencyMonojetSelections->SetBinContent(4,efficiencyMonojetSelections->GetBinContent(4)+1);
    efficiencyMonoVSelections->SetBinContent(4,efficiencyMonoVSelections->GetBinContent(4)+1);
    efficiencyVBFSelections->SetBinContent(4,efficiencyVBFSelections->GetBinContent(4)+1);
    //
    if(not *fcsc)  continue;
    if(not *fcsct) continue;
    if(not *feeb)  continue;
    if(not *fetp)  continue;
    if(not *fvtx)  continue;
    if(not *fbadmu) continue;
    if(not *fbadch) continue; 
    if(not *fhbhe)  continue;
    if(not *fhbiso) continue;

    efficiencyMonojetSelections->SetBinContent(5,efficiencyMonojetSelections->GetBinContent(5)+1);
    efficiencyMonoVSelections->SetBinContent(5,efficiencyMonoVSelections->GetBinContent(5)+1);
    efficiencyVBFSelections->SetBinContent(5,efficiencyVBFSelections->GetBinContent(5)+1);
    
    if(singleMuon or doubleMuon){
      if(not *hltsinglemu) continue;    
      efficiencyMonojetSelections->SetBinContent(6,efficiencyMonojetSelections->GetBinContent(6)+1);
      efficiencyMonoVSelections->SetBinContent(6,efficiencyMonoVSelections->GetBinContent(6)+1);
      efficiencyVBFSelections->SetBinContent(6,efficiencyVBFSelections->GetBinContent(6)+1);
      if(*mu1pt < 20) continue;
      if(fabs(*mu1eta) > 2.4) continue;
      efficiencyMonojetSelections->SetBinContent(7,efficiencyMonojetSelections->GetBinContent(7)+1);
      efficiencyMonoVSelections->SetBinContent(7,efficiencyMonoVSelections->GetBinContent(7)+1);
      efficiencyVBFSelections->SetBinContent(7,efficiencyVBFSelections->GetBinContent(7)+1);
      if(not doubleMuon  and *mu1id != 1) continue;      
      if(not doubleMuon  and *nmuons != 1) continue;      
      else if(doubleMuon and *nmuons != 2) continue;      
      if(doubleMuon){
	if(*mu1pid == *mu2pid) continue; //opposite charge
	TLorentzVector mu1, mu2;
	mu1.SetPtEtaPhiM(*mu1pt,*mu1eta,*mu1phi,0.);
	mu2.SetPtEtaPhiM(*mu2pt,*mu2eta,*mu2phi,0.);
	if((mu1+mu2).M() < 60 or (mu1+mu2).M() > 120) continue;
      }
      efficiencyMonojetSelections->SetBinContent(8,efficiencyMonojetSelections->GetBinContent(8)+1);
      efficiencyMonoVSelections->SetBinContent(8,efficiencyMonoVSelections->GetBinContent(8)+1);
      efficiencyVBFSelections->SetBinContent(8,efficiencyVBFSelections->GetBinContent(8)+1);
      if(*nelectrons > 0 ) continue;
      efficiencyMonojetSelections->SetBinContent(9,efficiencyMonojetSelections->GetBinContent(9)+1);
      efficiencyMonoVSelections->SetBinContent(9,efficiencyMonoVSelections->GetBinContent(9)+1);
      efficiencyVBFSelections->SetBinContent(9,efficiencyVBFSelections->GetBinContent(9)+1);
      if(fabs(*metpf-*metcalo)/(*metpf) > 0.5) continue;
      efficiencyMonojetSelections->SetBinContent(10,efficiencyMonojetSelections->GetBinContent(10)+1);
      efficiencyMonoVSelections->SetBinContent(10,efficiencyMonoVSelections->GetBinContent(10)+1);
      efficiencyVBFSelections->SetBinContent(10,efficiencyVBFSelections->GetBinContent(10)+1);
      if(*jmmdphi < 0.5) continue;      
      efficiencyMonojetSelections->SetBinContent(11,efficiencyMonojetSelections->GetBinContent(11)+1);
      efficiencyMonoVSelections->SetBinContent(11,efficiencyMonoVSelections->GetBinContent(11)+1);
      efficiencyVBFSelections->SetBinContent(11,efficiencyVBFSelections->GetBinContent(11)+1);
      // transverse mass cut
      if(not doubleMuon){
	float dphi = fabs(*mu1phi-*metphi);
	if(dphi > TMath::Pi())
	  dphi = 2*TMath::Pi()-dphi;
	float mtw = sqrt(2*(*mu1pt)*(*met)*(1-cos(dphi)));
	if(mtw > 160) continue;
      }
      
      // MONOJET
      if(jetpt->at(0) > 100 and fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) > 0.1 and jetnhfrac->at(0) < 0.8 and jetchfrac->at(0) < 0.997){
	
	// monojet vs recoil and met
	hden_monojet_recoil->Fill(*mmet);

	efficiencyMonojetSelections->SetBinContent(12,efficiencyMonojetSelections->GetBinContent(12)+1);
	// numerator monojet
	if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300){
	  hnum_monojet_recoil->Fill(*mmet);	
	  efficiencyMonojetSelections->SetBinContent(13,efficiencyMonojetSelections->GetBinContent(13)+1);	 
	}
	
	if(*mmet > 800 and not (*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300 or *hltmwm100 or *hltmwm110)){
	  cout<<"Monojet analysis: run "<<*run<<" lumi "<<*lumi<<" event "<<*event<<" t1mumet "<<*mmet<<" muon pt "<<*mu1pt<<" muon eta "<<*mu1eta<<" jet pt "<<jetpt->at(0)<<" eta "<<jeteta->at(0)<<" njets "<<*nincjets<<" pfmet "<<*met<<" calomet "<<*metcalo<<endl;
	}
      }

      // denominator VBF
      if(*nincjets > 1 and jetpt->at(0) > leadingVBF and jetpt->at(1) > trailingVBF and *jmmdphi > jetmetdphi){

	// charge fraction cut on the central jet
	if(fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) < 0.1) continue;
	if(fabs(jeteta->at(0)) < 2.5 and jetnhfrac->at(0) > 0.8) continue;
	if(fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) > 0.997) continue;
	if(fabs(jeteta->at(0)) < 3.2 and fabs(jeteta->at(0)) > 3.0 and jetnhfrac->at(0) > 0.96) continue; 

	float deltaPhi = fabs(jetphi->at(0)-jetphi->at(1));
	if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi()-deltaPhi;
	if(deltaPhi > dphijj) continue;

	efficiencyVBFSelections->SetBinContent(12,efficiencyVBFSelections->GetBinContent(12)+1);

	TLorentzVector jet1, jet2;
	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	if(fabs(jeteta->at(0)-jeteta->at(1)) > detajj){
	  efficiencyVBFSelections->SetBinContent(13,efficiencyVBFSelections->GetBinContent(13)+1);
	  if((jet1+jet2).M() > mjj){
	    efficiencyVBFSelections->SetBinContent(14,efficiencyVBFSelections->GetBinContent(14)+1);
	    hden_vbf_recoil->Fill(*mmet);
	    if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300){
	      hnum_vbf_recoil->Fill(*mmet);
	      efficiencyVBFSelections->SetBinContent(15,efficiencyVBFSelections->GetBinContent(15)+1);
	    }
	  }
	}
	  
	if(fabs(jeteta->at(0)-jeteta->at(1)) > detajj){
	  for(size_t ibin = 0; ibin < bins_vbf_mjj.size()-1; ibin++){
	    if((jet1+jet2).M() > bins_vbf_mjj.at(ibin) and (jet1+jet2).M() <= bins_vbf_mjj.at(ibin+1)){
	      hden_vbf_mjj.at(ibin)->Fill(*mmet);
	      if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300){
		hnum_vbf_mjj.at(ibin)->Fill(*mmet);
	      }
	    }
	  }
	  
	  if((jet1+jet2).M() > mjj){
	    for(size_t ibin = 0; ibin < bins_vbf_detajj.size()-1; ibin++){
	      if(fabs(jeteta->at(0)-jeteta->at(1)) > bins_vbf_detajj.at(ibin) and fabs(jeteta->at(0)-jeteta->at(1)) <= bins_vbf_mjj.at(ibin+1)){
		hden_vbf_detajj.at(ibin)->Fill(*mmet);
		if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300){
		  hnum_vbf_detajj.at(ibin)->Fill(*mmet);
		}
	      }
	    }
	  }
	  
	  if(*mmet > 800 and not (*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300 or *hltmwm100 or *hltmwm110)){
	    cout<<"VBF analysis: run "<<*run<<" lumi "<<*lumi<<" event "<<*event<<" t1pfmet "<<*mmet<<" muon pt "<<*mu1pt<<" muon eta "<<*mu1eta<<" jet pt "<<jetpt->at(0)<<" eta "<<jeteta->at(0)<<" njets "<<*nincjets<<endl;
	    
	  }
	}
      }	
    }
    else{// electron case
      if(not *hlte) continue;    
      efficiencyMonojetSelections->SetBinContent(6,efficiencyMonojetSelections->GetBinContent(6)+1);
      efficiencyMonoVSelections->SetBinContent(6,efficiencyMonoVSelections->GetBinContent(6)+1);
      efficiencyVBFSelections->SetBinContent(6,efficiencyVBFSelections->GetBinContent(6)+1);
      if(*el1pt < 40) continue;
      if(fabs(*el1eta) > 2.5) continue;
      efficiencyMonojetSelections->SetBinContent(7,efficiencyMonojetSelections->GetBinContent(7)+1);
      efficiencyMonoVSelections->SetBinContent(7,efficiencyMonoVSelections->GetBinContent(7)+1);
      efficiencyVBFSelections->SetBinContent(7,efficiencyVBFSelections->GetBinContent(7)+1);
  
      if(*el1id != 1) continue;
      if(*nelectrons > 1) continue;
      efficiencyMonojetSelections->SetBinContent(8,efficiencyMonojetSelections->GetBinContent(8)+1);
      efficiencyMonoVSelections->SetBinContent(8,efficiencyMonoVSelections->GetBinContent(8)+1);
      efficiencyVBFSelections->SetBinContent(8,efficiencyVBFSelections->GetBinContent(8)+1);
  
      if(*nmuons > 0 ) continue;
      efficiencyMonojetSelections->SetBinContent(9,efficiencyMonojetSelections->GetBinContent(9)+1);
      efficiencyMonoVSelections->SetBinContent(9,efficiencyMonoVSelections->GetBinContent(9)+1);
      efficiencyVBFSelections->SetBinContent(9,efficiencyVBFSelections->GetBinContent(9)+1);
  
      if(fabs(*metpf-*metcalo)/(*emet) > 0.5) continue;
      efficiencyMonojetSelections->SetBinContent(10,efficiencyMonojetSelections->GetBinContent(10)+1);
      efficiencyMonoVSelections->SetBinContent(10,efficiencyMonoVSelections->GetBinContent(10)+1);
      efficiencyVBFSelections->SetBinContent(10,efficiencyVBFSelections->GetBinContent(10)+1);
  
      if(*jemdphi < 0.5) continue;
      efficiencyMonojetSelections->SetBinContent(11,efficiencyMonojetSelections->GetBinContent(11)+1);
      efficiencyMonoVSelections->SetBinContent(11,efficiencyMonoVSelections->GetBinContent(11)+1);
      efficiencyVBFSelections->SetBinContent(11,efficiencyVBFSelections->GetBinContent(11)+1);

      // transverse mass cut
      float dphi = fabs(*el1phi-*metphi);
      if(dphi > TMath::Pi())
	dphi = 2*TMath::Pi()-dphi;
      float mtw = sqrt(2*(*el1pt)*(*met)*(1-cos(dphi)));
      if(mtw > 160) continue;
      if(*met < 50) continue;
      
      // denominator monojet
      if(jetpt->at(0) > 100 and fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) > 0.1 and jetnhfrac->at(0) < 0.8 and jetchfrac->at(0) < 0.997){
	hden_monojet_recoil->Fill(*met);
	efficiencyMonojetSelections->SetBinContent(12,efficiencyMonojetSelections->GetBinContent(12)+1);
	// numerator monojet
	if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300 or *hltmwm100 or *hltmwm110){
	  efficiencyMonojetSelections->SetBinContent(13,efficiencyMonojetSelections->GetBinContent(13)+1);
	  hnum_monojet_recoil->Fill(*met);	
	}
	
	if(*met > 800 and not (*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300 or *hltmwm100 or *hltmwm110)){
	  cout<<"Monojet analysis: run "<<*run<<" lumi "<<*lumi<<" event "<<*event<<" t1pfmet "<<*emet<<" ele pt "<<*el1pt<<" ele eta "<<*el1eta<<" jet pt "<<jetpt->at(0)<<" eta "<<jeteta->at(0)<<" njets "<<*nincjets<<" met "<<*met<<" calomet "<<*metcalo<<endl;
	}
      }

      // denominator VBF                                                                                                                                                                               
      if(*nincjets > 1 and jetpt->at(0) > leadingVBF and jetpt->at(1) > trailingVBF and *jmmdphi > jetmetdphi){

        // charge fraction cut on the central jet                                                                                                                                                     
        if(fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) < 0.1) continue;
        if(fabs(jeteta->at(0)) < 2.5 and jetnhfrac->at(0) > 0.8) continue;
        if(fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) > 0.997) continue;
        if(fabs(jeteta->at(0)) < 3.2 and fabs(jeteta->at(0)) > 3.0 and jetnhfrac->at(0) > 0.96) continue;
	float deltaPhi = fabs(jetphi->at(0)-jetphi->at(1));
	if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi()-deltaPhi;
	if(deltaPhi > dphijj) continue;

        efficiencyVBFSelections->SetBinContent(12,efficiencyVBFSelections->GetBinContent(12)+1);

        TLorentzVector jet1, jet2;
        jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
        jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
        if(fabs(jeteta->at(0)-jeteta->at(1)) > detajj){
          efficiencyVBFSelections->SetBinContent(13,efficiencyVBFSelections->GetBinContent(13)+1);
          if((jet1+jet2).M() > mjj){
            efficiencyVBFSelections->SetBinContent(14,efficiencyVBFSelections->GetBinContent(14)+1);
            hden_vbf_recoil->Fill(*met);
            if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300){
              hnum_vbf_recoil->Fill(*met);
	      efficiencyVBFSelections->SetBinContent(15,efficiencyVBFSelections->GetBinContent(15)+1);
            }
            if(*mmet > 800 and not (*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300 or *hltmwm100 or *hltmwm110)){
              cout<<"VBF analysis: run "<<*run<<" lumi "<<*lumi<<" event "<<*event<<" t1pfmet "<<*mmet<<" muon pt "<<*mu1pt<<" muon eta "<<*mu1eta<<" jet pt "<<jetpt->at(0)<<" eta "<<jeteta->at(0)<<" njets "<<*nincjets<<endl;
	      
            }
          }
        }

	if(fabs(jeteta->at(0)-jeteta->at(1)) > detajj){
	  for(size_t ibin = 0; ibin < bins_vbf_mjj.size()-1; ibin++){
	    if((jet1+jet2).M() > bins_vbf_mjj.at(ibin) and (jet1+jet2).M() <= bins_vbf_mjj.at(ibin+1)){
	      hden_vbf_mjj.at(ibin)->Fill(*met);
	      if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300){
		hnum_vbf_mjj.at(ibin)->Fill(*met);
	      }
	    }
	  }
	  
	  if((jet1+jet2).M() > mjj){
	    for(size_t ibin = 0; ibin < bins_vbf_detajj.size()-1; ibin++){
	      if(fabs(jeteta->at(0)-jeteta->at(1)) > bins_vbf_detajj.at(ibin) and fabs(jeteta->at(0)-jeteta->at(1)) <= bins_vbf_mjj.at(ibin+1)){
		hden_vbf_detajj.at(ibin)->Fill(*met);
		if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300){
		  hnum_vbf_detajj.at(ibin)->Fill(*met);
		}
	      }
	    }
	  }	  
        }
      }
    }
  }
    
  cout<<endl;
  
  //////////////
  for(int iBin = 0; iBin < efficiencyMonojetSelections->GetNbinsX(); iBin++){
    efficiencyMonojetSelections->SetBinContent(iBin+1,efficiencyMonojetSelections->GetBinContent(iBin+1)/double(nEvents));
    cout<<"Monojet selection efficiency: "<<efficiencyMonojetSelections->GetXaxis()->GetBinLabel(iBin+1)<<" eff "<<efficiencyMonojetSelections->GetBinContent(iBin+1)<<endl;
  }

  for(int iBin = 0; iBin < efficiencyVBFSelections->GetNbinsX(); iBin++){
    efficiencyVBFSelections->SetBinContent(iBin+1,efficiencyVBFSelections->GetBinContent(iBin+1)/double(nEvents));
    cout<<"VBF selection efficiency: "<<efficiencyVBFSelections->GetXaxis()->GetBinLabel(iBin+1)<<" eff "<<efficiencyVBFSelections->GetBinContent(iBin+1)<<endl;
  }


  TEfficiency* eff_monojet_recoil = new TEfficiency(*hnum_monojet_recoil,*hden_monojet_recoil);
  eff_monojet_recoil->SetMarkerColor(kBlack);
  eff_monojet_recoil->SetLineColor(kBlack);
  eff_monojet_recoil->SetMarkerStyle(20);
  eff_monojet_recoil->SetMarkerSize(1);
  fitfunc_monojet_recoil->SetLineColor(kBlue);
  fitfunc_monojet_recoil->SetLineWidth(2);
  
  TEfficiency* eff_vbf_recoil = new TEfficiency(*hnum_vbf_recoil,*hden_vbf_recoil);
  eff_vbf_recoil->SetMarkerColor(kBlack);
  eff_vbf_recoil->SetLineColor(kBlack);
  eff_vbf_recoil->SetMarkerStyle(20);
  eff_vbf_recoil->SetMarkerSize(1);
  fitfunc_vbf_recoil->SetLineColor(kBlue);
  fitfunc_vbf_recoil->SetLineWidth(2);

  
  vector<TEfficiency*> eff_vbf_mjj;
  for(size_t ibin = 0; ibin <  hnum_vbf_mjj.size(); ibin++){
    eff_vbf_mjj.push_back(new TEfficiency(*hnum_vbf_mjj.at(ibin),*hden_vbf_mjj.at(ibin)));
    eff_vbf_mjj.back()->SetMarkerColor(kBlack);
    eff_vbf_mjj.back()->SetLineColor(kBlack);
    eff_vbf_mjj.back()->SetMarkerStyle(20);
    eff_vbf_mjj.back()->SetMarkerSize(1);
    fitfunc_vbf_mjj.at(ibin)->SetLineColor(kBlue);
    fitfunc_vbf_mjj.at(ibin)->SetLineWidth(2);
  }

  vector<TEfficiency*> eff_vbf_detajj;
  for(size_t ibin = 0; ibin <  hnum_vbf_detajj.size(); ibin++){
    eff_vbf_detajj.push_back(new TEfficiency(*hnum_vbf_detajj.at(ibin),*hden_vbf_detajj.at(ibin)));
    eff_vbf_detajj.back()->SetMarkerColor(kBlack);
    eff_vbf_detajj.back()->SetLineColor(kBlack);
    eff_vbf_detajj.back()->SetMarkerStyle(20);
    eff_vbf_detajj.back()->SetMarkerSize(1);
    fitfunc_vbf_detajj.at(ibin)->SetLineColor(kBlue);
    fitfunc_vbf_detajj.at(ibin)->SetLineWidth(2);
  }

  bool forplot = false;
  if(singleMuon) forplot = true;
  if(doubleMuon) forplot = true;

  plotTurnOn(canvas,eff_monojet_recoil,fitfunc_monojet_recoil,"Recoil [GeV]","recoil_monojet",ouputDIR,luminosity,forplot);
  plotTurnOn(canvas,eff_vbf_recoil,fitfunc_vbf_recoil,"Recoil [GeV]","recoil_vbf",ouputDIR,luminosity,forplot);

  for(size_t ibin = 0; ibin < eff_vbf_mjj.size(); ibin++){
    plotTurnOn(canvas,eff_vbf_mjj.at(ibin),fitfunc_vbf_mjj.at(ibin),"Recoil [GeV]",
	       Form("mjj_vbf_%.1f_%.1f",bins_vbf_mjj.at(ibin),bins_vbf_mjj.at(ibin+1)),
	       ouputDIR,luminosity,forplot,Form("%.1f < m_{jj} < %.1f",bins_vbf_mjj.at(ibin),bins_vbf_mjj.at(ibin+1)));
  }
  for(size_t ibin = 0; ibin < eff_vbf_detajj.size(); ibin++){
    plotTurnOn(canvas,eff_vbf_detajj.at(ibin),fitfunc_vbf_detajj.at(ibin),"Recoil [GeV]",
	       Form("detajj_vbf_%.1f_%.1f",bins_vbf_detajj.at(ibin),bins_vbf_detajj.at(ibin+1)),
	       ouputDIR,luminosity,forplot,Form("%.1f < #Delta#eta_{jj} < %.1f",bins_vbf_detajj.at(ibin),bins_vbf_detajj.at(ibin+1)));
  }
  
}

void plotTurnOn(TCanvas* canvas, TEfficiency* eff, TF1* fitfunc, const string & axisLabel, const TString & postfix, const string & ouputDIR, const float  & luminosity, const bool & singleMuon, const TString & banner){


  TH1* frame = canvas->DrawFrame(fitfunc->GetXaxis()->GetXmin(),0.,fitfunc->GetXaxis()->GetXmax(), 1.1, "");
  frame->GetXaxis()->SetTitle(axisLabel.c_str());  
  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetYaxis()->SetLabelSize(0.8*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetTitleSize(0.8*frame->GetYaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleOffset(1.0);

  // make the fit
  TGraphAsymmErrors* graph = eff->CreateGraph();
  TFitResultPtr fitResult = graph->Fit(fitfunc,"RS");
  int npoints       = 350;                                                                                                                                                                   
  TH1F* error_band  = new TH1F(Form("%s_error_band",fitfunc->GetName()),"",npoints,fitfunc->GetXaxis()->GetXmin(),fitfunc->GetXaxis()->GetXmax());
  if(drawUncertaintyBand)
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(error_band,0.68);
  
  canvas->SetRightMargin(0.075);
  canvas->SetTopMargin(0.06);
  canvas->Draw();
  canvas->cd();
  frame->Draw();
  if(drawUncertaintyBand){    
    error_band->SetFillColor(kBlue);
    error_band->SetFillStyle(3001);
    error_band->Draw("e4 same");  
  }
  eff->Draw("E1PSAME");
  fitfunc->Draw("SAME");
  
  canvas->RedrawAxis();
  CMS_lumi(canvas,string(Form("%.2f",luminosity)),true);

  TLegend leg (0.6,0.3,0.9,0.4);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  leg.AddEntry((TObject*)0,banner,"");
  leg.Draw("same");

  if(singleMuon){
    canvas->SaveAs((ouputDIR+"/metTriggerEfficiency_"+string(postfix)+".png").c_str(),"png");
    canvas->SaveAs((ouputDIR+"/metTriggerEfficiency_"+string(postfix)+".pdf").c_str(),"pdf");
    TFile* outputFile = new TFile((ouputDIR+"/metTriggerEfficiency_"+string(postfix)+".root").c_str(),"RECREATE");
    outputFile->cd();
    // efficiency
    eff->Write("efficiency");
    fitfunc->Write("efficiency_func");
  }
  else{
    canvas->SaveAs((ouputDIR+"/metTriggerEfficiency_ele_"+string(postfix)+".png").c_str(),"png");
    canvas->SaveAs((ouputDIR+"/metTriggerEfficiency_ele_"+string(postfix)+".pdf").c_str(),"pdf");
    TFile* outputFile = new TFile((ouputDIR+"/metTriggerEfficiency_ele_"+string(postfix)+".root").c_str(),"RECREATE");
    outputFile->cd();
    // efficiency
    eff->Write("efficiency");
    fitfunc->Write("efficiency_func");
  }

}

