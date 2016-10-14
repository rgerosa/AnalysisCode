#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

vector <float> bins_monojet_muon = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 180., 200., 225, 250., 275., 300., 350., 400., 450., 500., 550., 650., 800., 1000., 1250};
vector <float> bins_vbf_muon     = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 180., 200., 225, 250., 275., 300., 350., 400., 450., 500., 550., 650., 800., 1000., 1500};
vector <float> bins_vbf_elec     = {0.,50.,60.,70.,80.,85.,90.,95.,100,110.,120.,140.,160., 180., 200., 225.,250., 300.,350.,450,550,650,800,1000.,1250,1500};
vector <float> bins_monojet_elec = {0.,50.,60.,70.,80.,85.,90.,95.,100,110.,120.,140.,160., 180., 200., 225.,250., 300.,350.,450,550,650,800,1000.,1250,1500};

//vector<string> RunEra = {"Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G"};
//vector<string> RunEra = {"Run2016B","Run2016C","Run2016D"};
vector<string> RunEra = {"Run2016B","Run2016C","Run2016D"};

static float leadingVBF  = 60;
static float trailingVBF = 50;
static float leadingVBFtight  = 80;
static float trailingVBFtight = 60;
static float detajj = 2.5; 
static float mjj    = 450; 
static float jetmetdphi = 1.0; 

void makeMETTriggerEfficiency(string inputDIR, string ouputDIR, float luminosity = 0.81, bool singleMuon = true) {

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+ouputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();

  // input tree
  TChain* tree = new TChain("tree/tree");
  // use only a subset of directories
  if(singleMuon)
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
	  tree->Add(line.c_str());
	}
      }
      system("rm list.txt");
    }
  }
  system("rm list_dir.txt");

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

  // fitting function for the turn-on 
  TF1 *fitfunc_monojet = new TF1("fitfunc_monojet", ErfCB, bins_monojet.front(), bins_monojet.back(), 5);
  fitfunc_monojet->SetParameters(117., 25., 30., 4., 1.);
  TF1 *fitfunc_vbf_loose = new TF1("fitfunc_vbf_loose", ErfCB, bins_vbf.front(), bins_vbf.back(), 5);
  fitfunc_vbf_loose->SetParameters(117., 25., 30., 4., 1.);
  TF1 *fitfunc_vbf_tight = new TF1("fitfunc_vbf_tight", ErfCB, bins_vbf.front(), bins_vbf.back(), 5);
  fitfunc_vbf_tight->SetParameters(117., 25., 30., 4., 1.);
  TF1 *fitfunc_vbf_central = new TF1("fitfunc_vbf_central", ErfCB, bins_vbf.front(), bins_vbf.back(), 5);
  fitfunc_vbf_central->SetParameters(117., 25., 30., 4., 1.);

  
  TH1F* hnum_monojet = new TH1F("hnum_monojet", "", bins_monojet.size()-1, &bins_monojet[0]);
  TH1F* hden_monojet = new TH1F("hden_monojet", "", bins_monojet.size()-1, &bins_monojet[0]);
  hnum_monojet->Sumw2();
  hden_monojet->Sumw2();
  TH1F* hnum_vbf_loose = new TH1F("hnum_vbf_loose", "", bins_vbf.size()-1, &bins_vbf[0]);
  TH1F* hden_vbf_loose = new TH1F("hden_vbf_loose", "", bins_vbf.size()-1, &bins_vbf[0]);
  hnum_vbf_loose->Sumw2();
  hden_vbf_loose->Sumw2();
  TH1F* hnum_vbf_tight = new TH1F("hnum_vbf_tight", "", bins_vbf.size()-1, &bins_vbf[0]);
  TH1F* hden_vbf_tight = new TH1F("hden_vbf_tight", "", bins_vbf.size()-1, &bins_vbf[0]);
  hnum_vbf_tight->Sumw2();
  hden_vbf_tight->Sumw2();
  TH1F* hnum_vbf_central = new TH1F("hnum_vbf_central", "", bins_vbf.size()-1, &bins_vbf[0]);
  TH1F* hden_vbf_central = new TH1F("hden_vbf_central", "", bins_vbf.size()-1, &bins_vbf[0]);
  hnum_vbf_central->Sumw2();
  hden_vbf_central->Sumw2();
  
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
  TTreeReaderValue<UChar_t> hltjm90   (reader,"hltjetmet90");
  TTreeReaderValue<UChar_t> hltjm100  (reader,"hltjetmet100");
  TTreeReaderValue<UChar_t> hltjm110  (reader,"hltjetmet110");
  TTreeReaderValue<UChar_t> hltsinglemu (reader,"hltsinglemu");
  TTreeReaderValue<UChar_t> hlte      (reader,"hltsingleel");
  TTreeReaderValue<double>  mu1pt     (reader,"mu1pt");
  TTreeReaderValue<double>  mu1eta    (reader,"mu1eta");
  TTreeReaderValue<double>  mu1phi    (reader,"mu1phi");
  TTreeReaderValue<int>     mu1id     (reader,"mu1id");
  TTreeReaderValue<double>  el1pt     (reader,"el1pt");
  TTreeReaderValue<double>  el1eta    (reader,"el1eta");
  TTreeReaderValue<double>  el1phi    (reader,"el1phi");
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
  long int nTotal = tree->GetEntries();
  cout<<"Total number of events: "<<nTotal<<endl;
  long int nEvents = 0;

  TH1D* efficiencyMonojetSelections = new TH1D("efficiencyMonojetSelections","efficiencyMonojetSelections",15,0,16);
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

  TH1D* efficiencyVBFSelections_loose = (TH1D*) efficiencyMonojetSelections->Clone("efficiencyVBFSelections_loose");
  efficiencyVBFSelections_loose->GetXaxis()->SetBinLabel(13,"Deta cut");
  efficiencyVBFSelections_loose->GetXaxis()->SetBinLabel(14,"Mjj cut");
  efficiencyVBFSelections_loose->GetXaxis()->SetBinLabel(15,"Trigger");

  TH1D* efficiencyVBFSelections_tight = (TH1D*) efficiencyMonojetSelections->Clone("efficiencyVBFSelections_tight");
  efficiencyVBFSelections_tight->GetXaxis()->SetBinLabel(13,"Deta cut");
  efficiencyVBFSelections_tight->GetXaxis()->SetBinLabel(14,"Mjj cut");
  efficiencyVBFSelections_tight->GetXaxis()->SetBinLabel(15,"Trigger");

  TH1D* efficiencyVBFSelections_central = (TH1D*) efficiencyMonojetSelections->Clone("efficiencyVBFSelections_central");
  efficiencyVBFSelections_central->GetXaxis()->SetBinLabel(13,"Deta cut");
  efficiencyVBFSelections_central->GetXaxis()->SetBinLabel(14,"Mjj cut");
  efficiencyVBFSelections_central->GetXaxis()->SetBinLabel(15,"Trigger");

  long int nPart = 100000;
  while(reader.Next()){
    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

   // define the denominator
    efficiencyMonojetSelections->SetBinContent(1,efficiencyMonojetSelections->GetBinContent(1)+1);
    efficiencyVBFSelections_loose->SetBinContent(1,efficiencyVBFSelections_loose->GetBinContent(1)+1);
    efficiencyVBFSelections_tight->SetBinContent(1,efficiencyVBFSelections_tight->GetBinContent(1)+1);
    efficiencyVBFSelections_central->SetBinContent(1,efficiencyVBFSelections_central->GetBinContent(1)+1);
    if(*nbjets   != 0) continue;
    ///
    efficiencyMonojetSelections->SetBinContent(2,efficiencyMonojetSelections->GetBinContent(2)+1);
    efficiencyVBFSelections_loose->SetBinContent(2,efficiencyVBFSelections_loose->GetBinContent(2)+1);
    efficiencyVBFSelections_tight->SetBinContent(2,efficiencyVBFSelections_tight->GetBinContent(2)+1);
    efficiencyVBFSelections_central->SetBinContent(2,efficiencyVBFSelections_central->GetBinContent(2)+1);
    if(*ntausraw != 0) continue;
    ///
    efficiencyMonojetSelections->SetBinContent(3,efficiencyMonojetSelections->GetBinContent(3)+1);
    efficiencyVBFSelections_loose->SetBinContent(3,efficiencyVBFSelections_loose->GetBinContent(3)+1);
    efficiencyVBFSelections_tight->SetBinContent(3,efficiencyVBFSelections_tight->GetBinContent(3)+1);
    efficiencyVBFSelections_central->SetBinContent(3,efficiencyVBFSelections_central->GetBinContent(3)+1);
    if(*nphotons  != 0) continue;
    efficiencyMonojetSelections->SetBinContent(4,efficiencyMonojetSelections->GetBinContent(4)+1);
    efficiencyVBFSelections_loose->SetBinContent(4,efficiencyVBFSelections_loose->GetBinContent(4)+1);
    efficiencyVBFSelections_tight->SetBinContent(4,efficiencyVBFSelections_tight->GetBinContent(4)+1);
    efficiencyVBFSelections_central->SetBinContent(4,efficiencyVBFSelections_central->GetBinContent(4)+1);
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
    efficiencyVBFSelections_loose->SetBinContent(5,efficiencyVBFSelections_loose->GetBinContent(5)+1);
    efficiencyVBFSelections_tight->SetBinContent(5,efficiencyVBFSelections_tight->GetBinContent(5)+1);
    efficiencyVBFSelections_central->SetBinContent(5,efficiencyVBFSelections_central->GetBinContent(5)+1);
    
    if(singleMuon){
      if(not *hltsinglemu) continue;    
      efficiencyMonojetSelections->SetBinContent(6,efficiencyMonojetSelections->GetBinContent(6)+1);
      efficiencyVBFSelections_loose->SetBinContent(6,efficiencyVBFSelections_loose->GetBinContent(6)+1);
      efficiencyVBFSelections_tight->SetBinContent(6,efficiencyVBFSelections_tight->GetBinContent(6)+1);
      efficiencyVBFSelections_central->SetBinContent(6,efficiencyVBFSelections_central->GetBinContent(6)+1);
      if(*mu1pt < 20) continue;
      if(fabs(*mu1eta) > 2.4) continue;
      efficiencyMonojetSelections->SetBinContent(7,efficiencyMonojetSelections->GetBinContent(7)+1);
      efficiencyVBFSelections_loose->SetBinContent(7,efficiencyVBFSelections_loose->GetBinContent(7)+1);
      efficiencyVBFSelections_tight->SetBinContent(7,efficiencyVBFSelections_tight->GetBinContent(7)+1);
      efficiencyVBFSelections_central->SetBinContent(7,efficiencyVBFSelections_central->GetBinContent(7)+1);
      if(*mu1id != 1) continue;
      if(*nmuons > 1) continue;      
      efficiencyMonojetSelections->SetBinContent(8,efficiencyMonojetSelections->GetBinContent(8)+1);
      efficiencyVBFSelections_loose->SetBinContent(8,efficiencyVBFSelections_loose->GetBinContent(8)+1);
      efficiencyVBFSelections_tight->SetBinContent(8,efficiencyVBFSelections_tight->GetBinContent(8)+1);
      efficiencyVBFSelections_central->SetBinContent(8,efficiencyVBFSelections_central->GetBinContent(8)+1);
      if(*nelectrons > 0 ) continue;
      efficiencyMonojetSelections->SetBinContent(9,efficiencyMonojetSelections->GetBinContent(9)+1);
      efficiencyVBFSelections_loose->SetBinContent(9,efficiencyVBFSelections_loose->GetBinContent(9)+1);
      efficiencyVBFSelections_tight->SetBinContent(9,efficiencyVBFSelections_tight->GetBinContent(9)+1);
      efficiencyVBFSelections_central->SetBinContent(9,efficiencyVBFSelections_central->GetBinContent(9)+1);
      if(fabs(*metpf-*metcalo)/(*metpf) > 0.5) continue;
      efficiencyMonojetSelections->SetBinContent(10,efficiencyMonojetSelections->GetBinContent(10)+1);
      efficiencyVBFSelections_loose->SetBinContent(10,efficiencyVBFSelections_loose->GetBinContent(10)+1);
      efficiencyVBFSelections_tight->SetBinContent(10,efficiencyVBFSelections_tight->GetBinContent(10)+1);
      efficiencyVBFSelections_central->SetBinContent(10,efficiencyVBFSelections_central->GetBinContent(10)+1);
      if(*jmmdphi < 0.5) continue;      
      efficiencyMonojetSelections->SetBinContent(11,efficiencyMonojetSelections->GetBinContent(11)+1);
      efficiencyVBFSelections_loose->SetBinContent(11,efficiencyVBFSelections_loose->GetBinContent(11)+1);
      efficiencyVBFSelections_tight->SetBinContent(11,efficiencyVBFSelections_tight->GetBinContent(11)+1);
      efficiencyVBFSelections_central->SetBinContent(11,efficiencyVBFSelections_central->GetBinContent(11)+1);
      // transverse mass cut
      float dphi = fabs(*mu1phi-*metphi);
      if(dphi > TMath::Pi())
	dphi = 2*TMath::Pi()-dphi;
      float mtw = sqrt(2*(*mu1pt)*(*met)*(1-cos(dphi)));
      if(mtw > 160) continue;
      // denominator monojet
      if(jetpt->at(0) > 100 and fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) > 0.1 and jetnhfrac->at(0) < 0.8){
	hden_monojet->Fill(*mmet);
	efficiencyMonojetSelections->SetBinContent(12,efficiencyMonojetSelections->GetBinContent(12)+1);
	// numerator monojet
	if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300){
	  hnum_monojet->Fill(*mmet);	
	  efficiencyMonojetSelections->SetBinContent(13,efficiencyMonojetSelections->GetBinContent(13)+1);	 
	}
	
	if(*mmet > 800 and not (*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300 or *hltmwm100 or *hltmwm110)){
	  cout<<"Monojet analysis: run "<<*run<<" lumi "<<*lumi<<" event "<<*event<<" t1mumet "<<*mmet<<" muon pt "<<*mu1pt<<" muon eta "<<*mu1eta<<" jet pt "<<jetpt->at(0)<<" eta "<<jeteta->at(0)<<" njets "<<*nincjets<<" pfmet "<<*met<<" calomet "<<*metcalo<<endl;
	}
      }
      
      // denominator VBF
      if(*nincjets > 1 and jetpt->at(0) > leadingVBF and jetpt->at(1) > trailingVBF and fabs(jeteta->at(0)-jeteta->at(1)) > detajj and *jmmdphi > jetmetdphi){

	// charge fraction cut on the central jet
	if(fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) < 0.1) continue;
	if(fabs(jeteta->at(0)) < 2.5 and jetnhfrac->at(0) > 0.8) continue;

	efficiencyVBFSelections_loose->SetBinContent(12,efficiencyVBFSelections_loose->GetBinContent(12)+1);
	if(jetpt->at(0) > leadingVBFtight and jetpt->at(1) > trailingVBFtight)
	  efficiencyVBFSelections_tight->SetBinContent(12,efficiencyVBFSelections_tight->GetBinContent(12)+1);
	if(fabs(jeteta->at(0)) < 2.5 and jetpt->at(0) > leadingVBFtight )
	  efficiencyVBFSelections_central->SetBinContent(12,efficiencyVBFSelections_central->GetBinContent(12)+1);

	TLorentzVector jet1, jet2;
	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	if((jet1+jet2).M() < mjj) continue;	
	efficiencyVBFSelections_loose->SetBinContent(13,efficiencyVBFSelections_loose->GetBinContent(13)+1);
	if(jetpt->at(0) > leadingVBFtight and jetpt->at(1) > trailingVBFtight)
	  efficiencyVBFSelections_tight->SetBinContent(13,efficiencyVBFSelections_tight->GetBinContent(13)+1);
	if(fabs(jeteta->at(0)) < 2.5 and jetpt->at(0) > leadingVBFtight )
	  efficiencyVBFSelections_central->SetBinContent(13,efficiencyVBFSelections_central->GetBinContent(13)+1);
	
	hden_vbf_loose->Fill(*mmet);
	if(jetpt->at(0) > leadingVBFtight and jetpt->at(1) > trailingVBFtight)
	  hden_vbf_tight->Fill(*mmet);
	if(fabs(jeteta->at(0)) < 2.5 and jetpt->at(0) > leadingVBFtight)
	  hden_vbf_central->Fill(*mmet);

	if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300){
	  hnum_vbf_loose->Fill(*mmet);	
	  efficiencyVBFSelections_loose->SetBinContent(14,efficiencyVBFSelections_loose->GetBinContent(14)+1);
	  if(jetpt->at(0) > leadingVBFtight and jetpt->at(1) > trailingVBFtight){
	    hnum_vbf_tight->Fill(*mmet);
	    efficiencyVBFSelections_tight->SetBinContent(14,efficiencyVBFSelections_tight->GetBinContent(14)+1);
	  }
	  if(fabs(jeteta->at(0)) < 2.5 and jetpt->at(0) > leadingVBFtight){
	    hnum_vbf_central->Fill(*mmet);
	    efficiencyVBFSelections_central->SetBinContent(14,efficiencyVBFSelections_central->GetBinContent(14)+1);
	  }
	}
	if(*mmet > 800 and not (*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300 or *hltmwm100 or *hltmwm110)){
	  cout<<"VBF analysis: run "<<*run<<" lumi "<<*lumi<<" event "<<*event<<" t1pfmet "<<*mmet<<" muon pt "<<*mu1pt<<" muon eta "<<*mu1eta<<" jet pt "<<jetpt->at(0)<<" eta "<<jeteta->at(0)<<" njets "<<*nincjets<<endl;
	}
      }
    }
    else{// electron case
      if(not *hlte) continue;    
      efficiencyMonojetSelections->SetBinContent(6,efficiencyMonojetSelections->GetBinContent(6)+1);
      efficiencyVBFSelections_loose->SetBinContent(6,efficiencyVBFSelections_loose->GetBinContent(6)+1);
      efficiencyVBFSelections_tight->SetBinContent(6,efficiencyVBFSelections_tight->GetBinContent(6)+1);
      efficiencyVBFSelections_central->SetBinContent(6,efficiencyVBFSelections_central->GetBinContent(6)+1);
      if(*el1pt < 40) continue;
      if(fabs(*el1eta) > 2.5) continue;
      efficiencyMonojetSelections->SetBinContent(7,efficiencyMonojetSelections->GetBinContent(7)+1);
      efficiencyVBFSelections_loose->SetBinContent(7,efficiencyVBFSelections_loose->GetBinContent(7)+1);
      efficiencyVBFSelections_tight->SetBinContent(7,efficiencyVBFSelections_tight->GetBinContent(7)+1);
      efficiencyVBFSelections_central->SetBinContent(7,efficiencyVBFSelections_central->GetBinContent(7)+1);
      if(*el1id != 1) continue;
      if(*nelectrons > 1) continue;
      efficiencyMonojetSelections->SetBinContent(8,efficiencyMonojetSelections->GetBinContent(8)+1);
      efficiencyVBFSelections_loose->SetBinContent(8,efficiencyVBFSelections_loose->GetBinContent(8)+1);
      efficiencyVBFSelections_tight->SetBinContent(8,efficiencyVBFSelections_tight->GetBinContent(8)+1);
      efficiencyVBFSelections_central->SetBinContent(8,efficiencyVBFSelections_central->GetBinContent(8)+1);
      if(*nmuons > 0 ) continue;
      efficiencyMonojetSelections->SetBinContent(9,efficiencyMonojetSelections->GetBinContent(9)+1);
      efficiencyVBFSelections_loose->SetBinContent(9,efficiencyVBFSelections_loose->GetBinContent(9)+1);
      efficiencyVBFSelections_tight->SetBinContent(9,efficiencyVBFSelections_tight->GetBinContent(9)+1);
      efficiencyVBFSelections_central->SetBinContent(9,efficiencyVBFSelections_central->GetBinContent(9)+1);
      if(fabs(*metpf-*metcalo)/(*emet) > 0.5) continue;
      efficiencyMonojetSelections->SetBinContent(10,efficiencyMonojetSelections->GetBinContent(10)+1);
      efficiencyVBFSelections_loose->SetBinContent(10,efficiencyVBFSelections_loose->GetBinContent(10)+1);
      efficiencyVBFSelections_tight->SetBinContent(10,efficiencyVBFSelections_tight->GetBinContent(10)+1);
      efficiencyVBFSelections_central->SetBinContent(10,efficiencyVBFSelections_central->GetBinContent(10)+1);
      if(*jemdphi < 0.5) continue;
      efficiencyMonojetSelections->SetBinContent(11,efficiencyMonojetSelections->GetBinContent(11)+1);
      efficiencyVBFSelections_loose->SetBinContent(11,efficiencyVBFSelections_loose->GetBinContent(11)+1);
      efficiencyVBFSelections_tight->SetBinContent(11,efficiencyVBFSelections_tight->GetBinContent(11)+1);
      efficiencyVBFSelections_central->SetBinContent(11,efficiencyVBFSelections_central->GetBinContent(11)+1);
      // transverse mass cut
      float dphi = fabs(*el1phi-*metphi);
      if(dphi > TMath::Pi())
	dphi = 2*TMath::Pi()-dphi;
      float mtw = sqrt(2*(*el1pt)*(*met)*(1-cos(dphi)));
      if(mtw > 160) continue;
      
      // denominator monojet
      if(jetpt->at(0) > 100 and fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) > 0.1 and jetnhfrac->at(0) < 0.8){
	hden_monojet->Fill(*met);
	efficiencyMonojetSelections->SetBinContent(12,efficiencyMonojetSelections->GetBinContent(12)+1);
	// numerator monojet
	if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300 or *hltmwm100 or *hltmwm110){
	  efficiencyMonojetSelections->SetBinContent(13,efficiencyMonojetSelections->GetBinContent(13)+1);
	  hnum_monojet->Fill(*met);	
	}
	
	if(*met > 800 and not (*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300 or *hltmwm100 or *hltmwm110)){
	  cout<<"Monojet analysis: run "<<*run<<" lumi "<<*lumi<<" event "<<*event<<" t1pfmet "<<*emet<<" ele pt "<<*el1pt<<" ele eta "<<*el1eta<<" jet pt "<<jetpt->at(0)<<" eta "<<jeteta->at(0)<<" njets "<<*nincjets<<" met "<<*met<<" calomet "<<*metcalo<<endl;
	}
      }

      // denominator VBF
      if(*nincjets > 1 and jetpt->at(0) > leadingVBF and jetpt->at(1) > trailingVBF and fabs(jeteta->at(0)-jeteta->at(1)) > detajj and *jemdphi > jetmetdphi){

	// charge fraction cut on the central jet
	if(fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) < 0.1) continue;
	if(fabs(jeteta->at(0)) < 2.5 and jetnhfrac->at(0) > 0.8) continue;

	efficiencyVBFSelections_loose->SetBinContent(12,efficiencyVBFSelections_loose->GetBinContent(12)+1);
	if(jetpt->at(0) > leadingVBFtight and jetpt->at(1) > trailingVBFtight)
	  efficiencyVBFSelections_tight->SetBinContent(12,efficiencyVBFSelections_tight->GetBinContent(12)+1);
	if(fabs(jeteta->at(0)) < 2.5 and jetpt->at(0) > leadingVBFtight )
	  efficiencyVBFSelections_central->SetBinContent(12,efficiencyVBFSelections_central->GetBinContent(12)+1);

	TLorentzVector jet1, jet2;
	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	if((jet1+jet2).M() < mjj) continue;	
	efficiencyVBFSelections_loose->SetBinContent(13,efficiencyVBFSelections_loose->GetBinContent(13)+1);

	if(jetpt->at(0) > leadingVBFtight and jetpt->at(1) > trailingVBFtight)
	  efficiencyVBFSelections_tight->SetBinContent(13,efficiencyVBFSelections_tight->GetBinContent(13)+1);
	if(fabs(jeteta->at(0)) < 2.5 and jetpt->at(0) > leadingVBFtight )
	  efficiencyVBFSelections_central->SetBinContent(13,efficiencyVBFSelections_central->GetBinContent(13)+1);

	hden_vbf_loose->Fill(*met);
	if(jetpt->at(0) > leadingVBFtight and jetpt->at(1) > trailingVBFtight)
	  hden_vbf_tight->Fill(*met);
	if(fabs(jeteta->at(0)) < 2.5 and jetpt->at(0) > leadingVBFtight )
	  hden_vbf_central->Fill(*met);
	  
	if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300 or *hltmwm100 or *hltmwm110){
	  efficiencyVBFSelections_loose->SetBinContent(14,efficiencyVBFSelections_loose->GetBinContent(14)+1);
	  hnum_vbf_loose->Fill(*met);	
	  if(jetpt->at(0) > leadingVBFtight and jetpt->at(1) > trailingVBFtight){
	    efficiencyVBFSelections_tight->SetBinContent(14,efficiencyVBFSelections_tight->GetBinContent(14)+1);
	    hnum_vbf_tight->Fill(*met);	
	  }
	  if(fabs(jeteta->at(0)) < 2.5 and jetpt->at(0) > leadingVBFtight ){
	    efficiencyVBFSelections_central->SetBinContent(14,efficiencyVBFSelections_central->GetBinContent(14)+1);
	    hnum_vbf_central->Fill(*met);	
	  }
	}
	
	if(*met > 800 and not (*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm170 or *hltmwm300 or *hltmwm100 or *hltmwm110)){
	  cout<<"VBF analysis: run "<<*run<<" lumi "<<*lumi<<" event "<<*event<<" t1pfmet "<<*emet<<" ele pt "<<*el1pt<<" ele eta "<<*el1eta<<" jet pt "<<jetpt->at(0)<<" eta "<<jeteta->at(0)<<" njets "<<*nincjets<<" met "<<*met<<" calomet "<<*metcalo<<endl;
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

  for(int iBin = 0; iBin < efficiencyVBFSelections_loose->GetNbinsX(); iBin++){
    efficiencyVBFSelections_loose->SetBinContent(iBin+1,efficiencyVBFSelections_loose->GetBinContent(iBin+1)/double(nEvents));
    cout<<"VBF selection efficiency: "<<efficiencyVBFSelections_loose->GetXaxis()->GetBinLabel(iBin+1)<<" eff "<<efficiencyVBFSelections_loose->GetBinContent(iBin+1)<<endl;
  }

  for(int iBin = 0; iBin < efficiencyVBFSelections_tight->GetNbinsX(); iBin++){
    efficiencyVBFSelections_tight->SetBinContent(iBin+1,efficiencyVBFSelections_tight->GetBinContent(iBin+1)/double(nEvents));
    cout<<"VBF selection efficiency: "<<efficiencyVBFSelections_tight->GetXaxis()->GetBinLabel(iBin+1)<<" eff "<<efficiencyVBFSelections_tight->GetBinContent(iBin+1)<<endl;
  }

  for(int iBin = 0; iBin < efficiencyVBFSelections_central->GetNbinsX(); iBin++){
    efficiencyVBFSelections_central->SetBinContent(iBin+1,efficiencyVBFSelections_central->GetBinContent(iBin+1)/double(nEvents));
    cout<<"VBF selection efficiency: "<<efficiencyVBFSelections_central->GetXaxis()->GetBinLabel(iBin+1)<<" eff "<<efficiencyVBFSelections_central->GetBinContent(iBin+1)<<endl;
  }

  
  TEfficiency* eff_monojet = new TEfficiency(*hnum_monojet,*hden_monojet);
  eff_monojet->Fit(fitfunc_monojet);
  eff_monojet->SetMarkerColor(kBlack);
  eff_monojet->SetLineColor(kBlack);
  eff_monojet->SetMarkerStyle(20);
  eff_monojet->SetMarkerSize(1);
  fitfunc_monojet->SetLineColor(kBlack);
  fitfunc_monojet->SetLineWidth(2);
  
  TEfficiency* eff_vbf_loose = new TEfficiency(*hnum_vbf_loose,*hden_vbf_loose);
  eff_vbf_loose->Fit(fitfunc_vbf_loose);
  eff_vbf_loose->SetMarkerColor(kRed);
  eff_vbf_loose->SetLineColor(kRed);
  eff_vbf_loose->SetMarkerStyle(20);
  eff_vbf_loose->SetMarkerSize(1);
  fitfunc_vbf_loose->SetLineColor(kRed);
  fitfunc_vbf_loose->SetLineWidth(2);
  
  TEfficiency* eff_vbf_tight = new TEfficiency(*hnum_vbf_tight,*hden_vbf_tight);
  eff_vbf_tight->Fit(fitfunc_vbf_tight);
  eff_vbf_tight->SetMarkerColor(kBlue);
  eff_vbf_tight->SetLineColor(kBlue);
  eff_vbf_tight->SetMarkerStyle(20);
  eff_vbf_tight->SetMarkerSize(1);
  fitfunc_vbf_tight->SetLineColor(kBlue);
  fitfunc_vbf_tight->SetLineWidth(2);

  TEfficiency* eff_vbf_central = new TEfficiency(*hnum_vbf_central,*hden_vbf_central);
  eff_vbf_central->Fit(fitfunc_vbf_central);
  eff_vbf_central->SetMarkerColor(kGreen+1);
  eff_vbf_central->SetLineColor(kGreen+1);
  eff_vbf_central->SetMarkerStyle(20);
  eff_vbf_central->SetMarkerSize(1);
  fitfunc_vbf_central->SetLineColor(kGreen+1);
  fitfunc_vbf_central->SetLineWidth(2);
  
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
  eff_vbf_loose->Draw("E1PSAME");
  eff_vbf_tight->Draw("E1PSAME");
  eff_vbf_central->Draw("E1PSAME");
  fitfunc_monojet->Draw("SAME");
  fitfunc_vbf_loose->Draw("SAME");
  fitfunc_vbf_tight->Draw("SAME");
  fitfunc_vbf_central->Draw("SAME");
  canvas->RedrawAxis();
  CMS_lumi(canvas,string(Form("%.2f",luminosity)),true);

  TLegend leg (0.6,0.2,0.9,0.4);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(eff_monojet,"mono-jet selection","PE");
  leg.AddEntry(eff_vbf_loose,"VBF jet p_{T} [60,50] GeV","PE");
  leg.AddEntry(eff_vbf_tight,"VBF jet p_{T} [80,60] GeV","PE");
  leg.AddEntry(eff_vbf_central,"VBF jet p_{T} [80,50] |#eta| < 2.5 GeV","PE");
  leg.Draw("same");

  if(singleMuon){
    canvas->SaveAs((ouputDIR+"/metTriggerEfficiency.png").c_str(),"png");
    canvas->SaveAs((ouputDIR+"/metTriggerEfficiency.pdf").c_str(),"pdf");
    TFile* outputFile = new TFile((ouputDIR+"/metTriggerEfficiency.root").c_str(),"RECREATE");
    outputFile->cd();
    // efficiency
    eff_monojet->Write("efficiency_monojet");
    eff_vbf_loose->Write("efficiency_vbf_loose");
    eff_vbf_tight->Write("efficiency_vbf_tight");
    eff_vbf_central->Write("efficiency_vbf_central");
    fitfunc_monojet->Write("efficiency_monojet_func");
    fitfunc_vbf_loose->Write("efficiency_vbf_loose_func");
    fitfunc_vbf_tight->Write("efficiency_vbf_tight_func");
    fitfunc_vbf_central->Write("efficiency_vbf_central_func");
    efficiencyMonojetSelections->Write("efficiencyMonojetSelections");
    efficiencyVBFSelections_loose->Write("efficiencyVBFSelections_loose");
    efficiencyVBFSelections_tight->Write("efficiencyVBFSelections_tight");
    efficiencyVBFSelections_central->Write("efficiencyVBFSelections_central");
  }
  else{
    canvas->SaveAs((ouputDIR+"/metTriggerEfficiency_el.png").c_str(),"png");
    canvas->SaveAs((ouputDIR+"/metTriggerEfficiency_el.pdf").c_str(),"pdf");
    TFile* outputFile = new TFile((ouputDIR+"/metTriggerEfficiency_el.root").c_str(),"RECREATE");
    outputFile->cd();
    eff_monojet->Write("efficiency_monojet");
    eff_vbf_loose->Write("efficiency_vbf_loose");
    eff_vbf_tight->Write("efficiency_vbf_tight");
    eff_vbf_central->Write("efficiency_vbf_central");
    fitfunc_monojet->Write("efficiency_monojet_func");
    fitfunc_vbf_loose->Write("efficiency_vbf_func_loose");
    fitfunc_vbf_tight->Write("efficiency_vbf_func_tight");
    fitfunc_vbf_central->Write("efficiency_vbf_func_central");
    efficiencyMonojetSelections->Write("efficiencyMonojetSelections");
    efficiencyVBFSelections_loose->Write("efficiencyVBFSelections_loose");
    efficiencyVBFSelections_tight->Write("efficiencyVBFSelections_tight");
    efficiencyVBFSelections_central->Write("efficiencyVBFSelections_central");
  }  
}

