#include <iostream>
#include <sstream>
#include <cmath>
#include "../triggerUtils.h"
#include "../../CMS_lumi.h"

// Recoil binning
vector <float> bins_recoil_cc_wmn  = {0.,30.,50.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,200.,225.,250.,275.,300.,350.,400.,450.,500.,550.,650.,800.,1000};
vector <float> bins_recoil_cf_wmn  = {0.,30.,50.,70.,90.,100.,110.,120.,130.,140.,150.,160.,180.,200.,225.,250.,275.,300.,400.,550.,1000};
vector <float> bins_recoil_fc_wmn  = {80.,100.,120.,140.,160.,180.,200.,225.,250.,300.,400.,500.,1000.};

vector <float> bins_recoil_cc_zmm  = {0.,30.,50.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,200.,225.,250.,275.,300.,350.,400.,450.,500.,550.,650.,800.,1000};
vector <float> bins_recoil_cf_zmm  = {0.,30.,50.,70.,90.,100.,110.,120.,130.,140.,150.,160.,180.,200.,225.,250.,275.,300.,400.,550.,1000};
vector <float> bins_recoil_fc_zmm  = {80.,100.,120.,140.,160.,180.,200.,225.,250.,300.,400.,500.,1000.};

// HT binning
vector <float> bins_ht_cc_wmn  = {0.,100.,150.,180.,220.,250.,300.,350.,400.,450.,500.,600.,800.,1000};
vector <float> bins_ht_cf_wmn  = {0.,100.,150.,180.,220.,250.,300.,350.,400.,450.,500.,600.,800.,1000};
vector <float> bins_ht_fc_wmn  = {0.,100.,150.,220.,300.,400.,600.,1000};
vector <float> bins_ht_cc_zmm  = {0.,100.,150.,180.,220.,250.,300.,350.,400.,450.,500.,600.,800.,1000};
vector <float> bins_ht_cf_zmm  = {0.,100.,150.,180.,220.,250.,300.,350.,400.,450.,500.,600.,800.,1000};
vector <float> bins_ht_fc_zmm  = {0.,100.,150.,220.,300.,400.,600.,1000};

// Mjj
vector <float> bins_mjj_cc_wmn     = {200.,600,1000.,1500.,2000.,2500,4000};
vector <float> bins_mjj_cc_zmm     = {200.,800,1400.,2200.,4000};
vector <float> bins_mjj_cf_wmn     = {200.,800,1200.,1600.,2000.,2500,4000};
vector <float> bins_mjj_cf_zmm     = {200.,800,1400.,2200.,4000};
vector <float> bins_mjj_fc_wmn     = {200.,800.,1600.,2500.,4000};
vector <float> bins_mjj_fc_zmm     = {200.,800.,1600.,2500.,4000};

// Cuts
vector <float> cuts_detajj_cc_wmn  = {1.0,1.5,2.0,2.5,3.0,3.5,4.0,5.0,10.};
vector <float> cuts_detajj_cf_wmn  = {3.0,3.5,4.0,5.0,6.0,10.};
vector <float> cuts_detajj_fc_wmn  = {3.0,10.};
vector <float> cuts_detajj_cc_zmm  = {1.0,2.0,3.0,4.0,10.};
vector <float> cuts_detajj_cf_zmm  = {3.0,4.0,5.0,10.};
vector <float> cuts_detajj_fc_zmm  = {3.0,10.};

// eras
vector<string> RunEra = {"Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G","Run2016H"};
// sample
enum class Sample {wmn,zmm};

//########## VBF selections
static float leadingVBF  = 80;
static float trailingVBF = 40;
static float detajj      = 1.0; 
static float mjj         = 200; 
static float jetmetdphi  = 0.5; 
static float dphijj      = 1.5;
static float recoil      = 250;
static bool  drawUncertaintyBand   = false;
static bool  useDoubleMuonTriggers = true;
static bool  applyJetSelections    = true;

/// plotting result
void plotTurnOn(TCanvas* canvas, 
		vector<TEfficiency* > eff, 
		vector<TF1*> fitfunc, 
		const string  & axisLabel, 
		const string  & postfix, 
		const string  & outputDIR,
		const vector<float> & binning,
                const string  & binningVar);

void plotTurnOn(TCanvas* canvas, 
		vector<TEfficiency* > eff, 
		vector<TF1*> fitfunc, 
		const string  & axisLabel, 
		const string  & postfix, 
		const string  & outputDIR,
		const vector<string> & label);


/// main function
void makeMETTriggerEfficiencyVBF_Data(string inputDIR, string outputDIR, Sample sample) {

  vector<float> bins_recoil_cc;
  vector<float> bins_recoil_cf;
  vector<float> bins_recoil_fc;
  vector<float> bins_ht_cc;
  vector<float> bins_ht_cf;
  vector<float> bins_ht_fc;
  vector<float> bins_mjj_cc;
  vector<float> bins_mjj_cf;
  vector<float> bins_mjj_fc;
  vector<float> cuts_detajj_cc;
  vector<float> cuts_detajj_cf;
  vector<float> cuts_detajj_fc;


  if(sample == Sample::wmn){
    bins_recoil_cc = bins_recoil_cc_wmn;
    bins_recoil_cf = bins_recoil_cf_wmn;
    bins_recoil_fc = bins_recoil_fc_wmn;
    bins_mjj_cc = bins_mjj_cc_wmn;
    bins_mjj_cf = bins_mjj_cf_wmn;
    bins_mjj_fc = bins_mjj_fc_wmn;
    cuts_detajj_cc = cuts_detajj_cc_wmn;
    cuts_detajj_cf = cuts_detajj_cf_wmn;
    cuts_detajj_fc = cuts_detajj_fc_wmn;
    bins_ht_cc = bins_ht_cc_wmn;
    bins_ht_cf = bins_ht_cf_wmn;
    bins_ht_fc = bins_ht_fc_wmn;
  }
  else{
    bins_recoil_cc = bins_recoil_cc_zmm;
    bins_recoil_cf = bins_recoil_cf_zmm;
    bins_recoil_fc = bins_recoil_fc_zmm;
    bins_mjj_cc = bins_mjj_cc_zmm;
    bins_mjj_cf = bins_mjj_cf_zmm;
    bins_mjj_fc = bins_mjj_fc_zmm;
    cuts_detajj_cc = cuts_detajj_cc_zmm;
    cuts_detajj_cf = cuts_detajj_cf_zmm;
    cuts_detajj_fc = cuts_detajj_fc_zmm;
    bins_ht_cc = bins_ht_cc_zmm;
    bins_ht_cf = bins_ht_cf_zmm;
    bins_ht_fc = bins_ht_fc_zmm;
  }


  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();

  // input tree --> build the list and fill the chain
  TChain* tree = new TChain("tree/tree");
  // use only a subset of directories
  system(("ls "+inputDIR+"  | grep SingleMu > list_dir.txt").c_str());

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

  TFile* outputFile = NULL;
  if(sample == Sample::wmn)
    outputFile = new TFile((outputDIR+"/efficiencyVBF_Wmn.root").c_str(),"RECREATE");
  else if(sample == Sample::zmm)
    outputFile = new TFile((outputDIR+"/efficiencyVBF_Zmm.root").c_str(),"RECREATE");
  outputFile->cd();


  /// CENTRAL JETS vs Mjj
  TH1F* hnum_mjj_cc = new TH1F("hnum_mjj_cc","", bins_mjj_cc.size()-1, &bins_mjj_cc[0]);
  TH1F* hden_mjj_cc = new TH1F("hden_mjj_cc","", bins_mjj_cc.size()-1, &bins_mjj_cc[0]);
  hnum_mjj_cc->Sumw2();
  hden_mjj_cc->Sumw2();

  //// CENTRAL JETS vs Recoil
  TH1F* hnum_recoil_cc = new TH1F("hnum_recoil_cc","", bins_recoil_cc.size()-1, &bins_recoil_cc[0]);
  TH1F* hden_recoil_cc = new TH1F("hden_recoil_cc","", bins_recoil_cc.size()-1, &bins_recoil_cc[0]);
  hnum_recoil_cc->Sumw2();
  hden_recoil_cc->Sumw2();
  TF1* fitfunc_recoil_cc = new TF1("fitfunc_recoil_cc","[0]*1./((1.+[1]*exp(-[2]*(x-[3])))^(1./[4]))",bins_recoil_cc.front(), bins_recoil_cc.back()); 
  fitfunc_recoil_cc->SetParameters(1.,0.01,0.02,50,0.01);
  fitfunc_recoil_cc->SetParLimits(0,0.,1.01);
  fitfunc_recoil_cc->SetParLimits(1,-100.,100.);
  fitfunc_recoil_cc->SetParLimits(3,0.,500.);
  fitfunc_recoil_cc->SetParLimits(4,0.,100);

  // Central jets vs HT
  TH1F* hnum_ht_cc = new TH1F("hnum_ht_cc","", bins_ht_cc.size()-1, &bins_ht_cc[0]);
  TH1F* hden_ht_cc = new TH1F("hden_ht_cc","", bins_ht_cc.size()-1, &bins_ht_cc[0]);
  hnum_ht_cc->Sumw2();
  hden_ht_cc->Sumw2();

  //// CENTRAL JETS --> differential vs recoil in bins of detajj
  vector<TH1F*> hnum_recoil_cc_vs_detajj;
  vector<TH1F*> hden_recoil_cc_vs_detajj;
  vector<TF1*>  fitfunc_recoil_cc_vs_detajj;
  
  for(size_t ibin = 0; ibin < cuts_detajj_cc.size()-1; ibin++){    
    hnum_recoil_cc_vs_detajj.push_back(new TH1F(Form("hnum_recoil_cc_vs_detajj_%.1f_%.1f",cuts_detajj_cc.at(ibin),cuts_detajj_cc.at(ibin+1)),"",bins_recoil_cc.size()-1,&bins_recoil_cc[0]));
    hden_recoil_cc_vs_detajj.push_back(new TH1F(Form("hden_recoil_cc_vs_detajj_%.1f_%.1f",cuts_detajj_cc.at(ibin),cuts_detajj_cc.at(ibin+1)),"",bins_recoil_cc.size()-1,&bins_recoil_cc[0]));
    hnum_recoil_cc_vs_detajj.back()->Sumw2();
    hden_recoil_cc_vs_detajj.back()->Sumw2();
    fitfunc_recoil_cc_vs_detajj.push_back(new TF1(Form("fitfunc_recoil_cc_vs_detajj_%.1f_%.1f",cuts_detajj_cc.at(ibin),cuts_detajj_cc.at(ibin+1)),"[0]*1./((1.+[1]*exp(-[2]*(x-[3])))^(1./[4]))",bins_recoil_cc.front(),bins_recoil_cc.back()));
    fitfunc_recoil_cc_vs_detajj.back()->SetParameters(1.,0.01,0.02,50,0.01);
    fitfunc_recoil_cc_vs_detajj.back()->SetParLimits(0,0.,1.01);
    fitfunc_recoil_cc_vs_detajj.back()->SetParLimits(1,-100.,100.);
    fitfunc_recoil_cc_vs_detajj.back()->SetParLimits(3,-500.,500.);
    fitfunc_recoil_cc_vs_detajj.back()->SetParLimits(4,0.,100);
  }
  
  // Leading central and trailing forwad
  TH1F* hnum_recoil_cf = new TH1F("hnum_recoil_cf","", bins_recoil_cf.size()-1, &bins_recoil_cf[0]);
  TH1F* hden_recoil_cf = new TH1F("hden_recoil_cf","", bins_recoil_cf.size()-1, &bins_recoil_cf[0]);
  hnum_recoil_cf->Sumw2();
  hden_recoil_cf->Sumw2();
  TF1* fitfunc_recoil_cf = new TF1("fitfunc_recoil_cf", "[0]*1./((1.+[1]*exp(-[2]*(x-[3])))^(1./[4]))", bins_recoil_cf.front(), bins_recoil_cf.back());
  fitfunc_recoil_cf->SetParameters(1.,0.01,0.02,50,0.01);
  fitfunc_recoil_cf->SetParLimits(0,0.,1.01);
  fitfunc_recoil_cf->SetParLimits(1,-100.,100.);
  fitfunc_recoil_cf->SetParLimits(3,-500.,500.);
  fitfunc_recoil_cf->SetParLimits(4,0.,100);

  ////
  TH1F* hnum_mjj_cf = new TH1F("hnum_mjj_cf","", bins_mjj_cf.size()-1, &bins_mjj_cf[0]);
  TH1F* hden_mjj_cf = new TH1F("hden_mjj_cf","", bins_mjj_cf.size()-1, &bins_mjj_cf[0]);
  hnum_mjj_cf->Sumw2();
  hden_mjj_cf->Sumw2();

  // Central jets vs HT
  TH1F* hnum_ht_cf = new TH1F("hnum_ht_cf","", bins_ht_cf.size()-1, &bins_ht_cf[0]);
  TH1F* hden_ht_cf = new TH1F("hden_ht_cf","", bins_ht_cf.size()-1, &bins_ht_cf[0]);
  hnum_ht_cf->Sumw2();
  hden_ht_cf->Sumw2();


  ////
  vector<TH1F*> hnum_recoil_cf_vs_detajj;
  vector<TH1F*> hden_recoil_cf_vs_detajj;
  vector<TF1*>  fitfunc_recoil_cf_vs_detajj;
  for(size_t ibin = 0; ibin < cuts_detajj_cf.size()-1; ibin++){
    hnum_recoil_cf_vs_detajj.push_back(new TH1F(Form("hnum_recoil_cf_vs_detajj_%.1f_%.1f",cuts_detajj_cf.at(ibin),cuts_detajj_cf.at(ibin+1)),"",bins_recoil_cf.size()-1,&bins_recoil_cf[0]));
    hden_recoil_cf_vs_detajj.push_back(new TH1F(Form("hden_recoil_cf_vs_detajj_%.1f_%.1f",cuts_detajj_cf.at(ibin),cuts_detajj_cf.at(ibin+1)),"",bins_recoil_cf.size()-1,&bins_recoil_cf[0]));
    hnum_recoil_cf_vs_detajj.back()->Sumw2();
    hden_recoil_cf_vs_detajj.back()->Sumw2();
    fitfunc_recoil_cf_vs_detajj.push_back(new TF1(Form("fitfunc_recoil_cf_vs_detajj_%.1f_%.1f",cuts_detajj_cf.at(ibin),cuts_detajj_cf.at(ibin+1)),"[0]*1./((1.+[1]*exp(-[2]*(x-[3])))^(1./[4]))",bins_recoil_cf.front(),bins_recoil_cf.back()));
    fitfunc_recoil_cf_vs_detajj.back()->SetParameters(1.,0.01,0.02,50,0.01);
    fitfunc_recoil_cf_vs_detajj.back()->SetParLimits(0,0.,1.01);
    fitfunc_recoil_cf_vs_detajj.back()->SetParLimits(1,-100.,100.);
    fitfunc_recoil_cf_vs_detajj.back()->SetParLimits(3,-500.,500.);
    fitfunc_recoil_cf_vs_detajj.back()->SetParLimits(4,0.,100);
  }

  // Leading forward and trailing central
  TH1F* hnum_recoil_fc = new TH1F("hnum_recoil_fc","", bins_recoil_fc.size()-1, &bins_recoil_fc[0]);
  TH1F* hden_recoil_fc = new TH1F("hden_recoil_fc","", bins_recoil_fc.size()-1, &bins_recoil_fc[0]);
  hnum_recoil_fc->Sumw2();
  hden_recoil_fc->Sumw2();
  TF1* fitfunc_recoil_fc = new TF1("fitfunc_recoil_fc", "[0]*1./((1.+[1]*exp(-[2]*(x-[3])))^(1./[4]))", bins_recoil_fc.front(), bins_recoil_fc.back());
  fitfunc_recoil_fc->SetParameters(1.,0.01,0.02,50,0.01);
  fitfunc_recoil_fc->SetParLimits(0,0.,1.01);
  fitfunc_recoil_fc->SetParLimits(1,-100.,100.);
  fitfunc_recoil_fc->SetParLimits(3,-500.,500.);
  fitfunc_recoil_fc->SetParLimits(4,0.,100);

  ////
  TH1F* hnum_mjj_fc = new TH1F("hnum_mjj_fc","", bins_mjj_fc.size()-1, &bins_mjj_fc[0]);
  TH1F* hden_mjj_fc = new TH1F("hden_mjj_fc","", bins_mjj_fc.size()-1, &bins_mjj_fc[0]);
  hnum_mjj_fc->Sumw2();
  hden_mjj_fc->Sumw2();

  TH1F* hnum_ht_fc = new TH1F("hnum_ht_fc","", bins_ht_fc.size()-1, &bins_ht_fc[0]);
  TH1F* hden_ht_fc = new TH1F("hden_ht_fc","", bins_ht_fc.size()-1, &bins_ht_fc[0]);
  hnum_ht_fc->Sumw2();
  hden_ht_fc->Sumw2();

  ////
  vector<TH1F*> hnum_recoil_fc_vs_detajj;
  vector<TH1F*> hden_recoil_fc_vs_detajj;
  vector<TF1*>  fitfunc_recoil_fc_vs_detajj;
  for(size_t ibin = 0; ibin < cuts_detajj_fc.size()-1; ibin++){
    hnum_recoil_fc_vs_detajj.push_back(new TH1F(Form("hnum_recoil_fc_vs_detajj_%.1f_%.1f",cuts_detajj_fc.at(ibin),cuts_detajj_fc.at(ibin+1)),"",bins_recoil_fc.size()-1,&bins_recoil_fc[0]));
    hden_recoil_fc_vs_detajj.push_back(new TH1F(Form("hden_recoil_fc_vs_detajj_%.1f_%.1f",cuts_detajj_fc.at(ibin),cuts_detajj_fc.at(ibin+1)),"",bins_recoil_fc.size()-1,&bins_recoil_fc[0]));
    hnum_recoil_fc_vs_detajj.back()->Sumw2();
    hden_recoil_fc_vs_detajj.back()->Sumw2();
    fitfunc_recoil_fc_vs_detajj.push_back(new TF1(Form("fitfunc_recoil_fc_vs_detajj_%.1f_%.1f",cuts_detajj_fc.at(ibin),cuts_detajj_fc.at(ibin+1)),"[0]*1./((1.+[1]*exp(-[2]*(x-[3])))^(1./[4]))",bins_recoil_fc.front(),bins_recoil_fc.back()));
    fitfunc_recoil_fc_vs_detajj.back()->SetParameters(1.,0.01,0.02,50,0.01);
    fitfunc_recoil_fc_vs_detajj.back()->SetParLimits(0,0.,1.01);
    fitfunc_recoil_fc_vs_detajj.back()->SetParLimits(1,-100.,100.);
    fitfunc_recoil_fc_vs_detajj.back()->SetParLimits(3,-500.,500.);
    fitfunc_recoil_fc_vs_detajj.back()->SetParLimits(4,0.,100);
  }


  /// Reader  
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
  TTreeReaderValue<UChar_t> hltdoublemu (reader,"hltdoublemu");
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
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");                                                                                                                                       
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");                                                                                                                                       
  TTreeReaderValue<unsigned int> ntaus       (reader,"ntaus");
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
  TTreeReaderValue<float> metpf       (reader,"pfmet");
  TTreeReaderValue<float> metcalo     (reader,"calomet");
  TTreeReaderValue<float> jmmdphi (reader,"incjetmumetdphimin4");
  TTreeReaderValue<float> zmass   (reader,"zmass");
 
  //////////////////
  long int nTotal = tree->GetEntries();
  cout<<"Total number of events: "<<nTotal<<endl;
  long int nEvents = 0;

  TH1D* efficiencyVBFSelections = new TH1D("efficiencyVBFSelections","efficiencyVBFSelections",17,0,16);
  efficiencyVBFSelections->Sumw2();
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(1,"Total");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(2,"MET preselection");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(3,"B-veto");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(4,"Tau-veto");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(5,"Photon-veto");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(6,"MET filters");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(7,"HLT single muon");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(8,"Lepton pT and eta");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(9,"Lepton ID");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(10,"Lepton Veto");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(11,"Calo-PF MET");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(12,"Jet met dphi");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(13,"Jet cuts");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(14,"Dphi-jj cut");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(15,"Deta-jj cut");
  efficiencyVBFSelections->GetXaxis()->SetBinLabel(16,"M-jj cut");

  long int nPart = 100000;
  while(reader.Next()){
    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    // define the denominator
    efficiencyVBFSelections->SetBinContent(1,efficiencyVBFSelections->GetBinContent(1)+1);
    if(*mmet < bins_recoil_cc.front()) continue;
    efficiencyVBFSelections->SetBinContent(2,efficiencyVBFSelections->GetBinContent(2)+1);    
    // b-veto
    if(*nbjets != 0) continue;
    efficiencyVBFSelections->SetBinContent(3,efficiencyVBFSelections->GetBinContent(3)+1);
    // tau-veto
    if(*ntaus != 0) continue;
    // photon veto
    efficiencyVBFSelections->SetBinContent(4,efficiencyVBFSelections->GetBinContent(4)+1);
    if(*nphotons  != 0) continue;
    efficiencyVBFSelections->SetBinContent(5,efficiencyVBFSelections->GetBinContent(5)+1);
    // MET filters
    if(not *fcsc)  continue;
    if(not *feeb)  continue;
    if(not *fetp)  continue;
    if(not *fvtx)  continue;
    if(not *fbadmu) continue;
    if(not *fbadch) continue; 
    if(not *fhbhe)  continue;
    if(not *fhbiso) continue;
    efficiencyVBFSelections->SetBinContent(6,efficiencyVBFSelections->GetBinContent(6)+1);
    
    // trigger selection for the denominator
    if(sample == Sample::wmn){
      if(not *hltsinglemu) continue;    
    }
    else if(sample == Sample::zmm){ // impose a di-muon trigger requriement 
      if(not *hltsinglemu) continue;
      if(useDoubleMuonTriggers and not *hltdoublemu) continue; 
    }      
    efficiencyVBFSelections->SetBinContent(7,efficiencyVBFSelections->GetBinContent(7)+1);
    
    // Leading muon
    if(*mu1pt < 20) continue;
    if(fabs(*mu1eta) > 2.4) continue;
    efficiencyVBFSelections->SetBinContent(8,efficiencyVBFSelections->GetBinContent(8)+1);

    // Sample selections
    if(sample == Sample::wmn){
      if(*mu1id != 1) continue;      
      if(*nmuons != 1) continue;      
      float dphi = fabs(*mu1phi-*metphi);
      if(dphi > TMath::Pi())
	dphi = 2*TMath::Pi()-dphi;
      float mtw = sqrt(2*(*mu1pt)*(*met)*(1-cos(dphi)));
      if(mtw > 160) continue;
    }
    else if(sample == Sample::zmm){
      if(*nmuons != 2) continue;      
      if(*mu1pid == *mu2pid) continue; //opposite charge
      if(not ((*mu1pt > 20 and *mu1id == 1) or (*mu2pt > 20 and *mu2id == 1))) continue;
      if(*zmass < 60 or *zmass > 120) continue;
    }
    efficiencyVBFSelections->SetBinContent(9,efficiencyVBFSelections->GetBinContent(9)+1);

    ///
    if(*nelectrons > 0 ) continue;
    efficiencyVBFSelections->SetBinContent(10,efficiencyVBFSelections->GetBinContent(10)+1);
    ///
    if(fabs(*metpf-*metcalo)/(*mmet) > 0.5) continue;
    efficiencyVBFSelections->SetBinContent(11,efficiencyVBFSelections->GetBinContent(11)+1);
    ////
    if(applyJetSelections and *jmmdphi < 0.5) continue;      
    efficiencyVBFSelections->SetBinContent(12,efficiencyVBFSelections->GetBinContent(12)+1);

    /// VBF requirements 
    if(*nincjets <= 1) continue;
    if(jetpt->at(0) < leadingVBF) continue;
    if(jetpt->at(1) < trailingVBF) continue;
    if(fabs(jeteta->at(0)) > 4.7) continue;
    if(fabs(jeteta->at(1)) > 4.7) continue;
    if(fabs(jeteta->at(0)) > 3 and fabs(jeteta->at(1)) > 3) continue;
    // charge fraction cut on the central jet
    if(fabs(jeteta->at(0)) < 2.4 and jetchfrac->at(0) < 0.1) continue;
    if(fabs(jeteta->at(0)) < 2.4 and jetnhfrac->at(0) > 0.8) continue;    
    efficiencyVBFSelections->SetBinContent(13,efficiencyVBFSelections->GetBinContent(13)+1);
    ////
    float deltaPhi = fabs(jetphi->at(0)-jetphi->at(1));
    if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi()-deltaPhi;
    if(sample == Sample::wmn and deltaPhi > dphijj) continue;
    else if(sample == Sample::zmm and deltaPhi > 2.0) continue;
    efficiencyVBFSelections->SetBinContent(14,efficiencyVBFSelections->GetBinContent(12)+1);
	
    /////
    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
    if(fabs(jet1.Eta()-jet2.Eta()) < detajj) continue;
    if(jet1.Eta()*jet2.Eta() > 0) continue;

    efficiencyVBFSelections->SetBinContent(15,efficiencyVBFSelections->GetBinContent(15)+1);

    if((jet1+jet2).M() < mjj) continue;
    efficiencyVBFSelections->SetBinContent(16,efficiencyVBFSelections->GetBinContent(16)+1);

    // value  to be used to fill histograms
    float mjjval = (jet1+jet2).M();
    if(mjjval > bins_mjj_cc.back()) mjjval = bins_mjj_cc.back()-1;

    float metval = *mmet;
    if(metval > bins_recoil_cc.back()) metval =  bins_recoil_cc.back()-1;

    // HT value
    float ht = 0;
    for(size_t ijet = 0 ; ijet < jetpt->size(); ijet++){
      if(fabs(jeteta->at(ijet)) > 3) continue;      
      if(fabs(jetpt->at(ijet)) < 30) continue;      
      ht += jetpt->at(ijet);
    }
    if(ht > bins_ht_cc.back()) ht = bins_ht_cc.back()-1;

    /// Start filling the histograms
    if(fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) < 3){ //central-central case

      // recoil inclusive
      hden_recoil_cc->Fill(metval);
      if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)
	hnum_recoil_cc->Fill(metval);

      /// binned in eta
      for(size_t ihist = 0; ihist < hnum_recoil_cc_vs_detajj.size(); ihist++){
	if(fabs(jeteta->at(0)-jeteta->at(1)) > cuts_detajj_cc.at(ihist) and fabs(jeteta->at(0)-jeteta->at(1)) <= cuts_detajj_cc.at(ihist+1)){
	  hden_recoil_cc_vs_detajj.at(ihist)->Fill(metval);
	  if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)
	    hnum_recoil_cc_vs_detajj.at(ihist)->Fill(metval);
	}
      }
            
      /// fill Mjj histograms
      if(*mmet > recoil){
	hden_mjj_cc->Fill(mjjval);
	if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)
	  hnum_mjj_cc->Fill(mjjval);
	
	hden_ht_cc->Fill(ht);
	if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)
	  hnum_ht_cc->Fill(ht);	  
      }      
    }
    
    // central foward case
    else if(fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) > 3){
	
	// inclusve vs recoil
      hden_recoil_cf->Fill(metval);
      if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)
	hnum_recoil_cf->Fill(metval);

	for(size_t ihist = 0; ihist < hnum_recoil_cf_vs_detajj.size(); ihist++){
	  if(fabs(jeteta->at(0)-jeteta->at(1)) > cuts_detajj_cf.at(ihist) and fabs(jeteta->at(0)-jeteta->at(1)) <= cuts_detajj_cf.at(ihist+1)){
	    hden_recoil_cf_vs_detajj.at(ihist)->Fill(metval);
	    if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)
	      hnum_recoil_cf_vs_detajj.at(ihist)->Fill(metval);
	  }
	}
	
	
	// Mjj histograms
	if(metval > recoil){ 	  
	  hden_mjj_cf->Fill(mjjval);
	  if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)
	    hnum_mjj_cf->Fill(mjjval);

	  hden_ht_cf->Fill(ht);
	  if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)
	    hnum_ht_cf->Fill(ht);	  
	}	
    }

    else if(fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) < 3){
      
      // inclusve vs recoil
      hden_recoil_fc->Fill(metval);
      if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)
	hnum_recoil_fc->Fill(metval);
      
      for(size_t ihist = 0; ihist < hnum_recoil_fc_vs_detajj.size(); ihist++){
	if(fabs(jeteta->at(0)-jeteta->at(1)) > cuts_detajj_fc.at(ihist) and fabs(jeteta->at(0)-jeteta->at(1)) <= cuts_detajj_fc.at(ihist+1)){
	  hden_recoil_fc_vs_detajj.at(ihist)->Fill(metval);
	  if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)
	    hnum_recoil_fc_vs_detajj.at(ihist)->Fill(metval);
	}
      }
      
      // Mjj histograms
      if(metval > recoil){ 	  
	hden_mjj_fc->Fill(mjjval);
	if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)
	  hnum_mjj_fc->Fill(mjjval);

	hden_ht_fc->Fill(ht);
	if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)
	  hnum_ht_fc->Fill(ht);
	
      }	
    }    
  }
  
  cout<<endl;  
  for(int iBin = 0; iBin < efficiencyVBFSelections->GetNbinsX(); iBin++){
    efficiencyVBFSelections->SetBinContent(iBin+1,efficiencyVBFSelections->GetBinContent(iBin+1)/double(nEvents));
    cout<<"Selection efficiency: "<<efficiencyVBFSelections->GetXaxis()->GetBinLabel(iBin+1)<<" eff "<<efficiencyVBFSelections->GetBinContent(iBin+1)<<endl;
  }

  // Build efficiency objects
  TEfficiency* eff_recoil_cc = new TEfficiency(*hnum_recoil_cc,*hden_recoil_cc);
  eff_recoil_cc->SetMarkerColor(kBlack);
  eff_recoil_cc->SetLineColor(kBlack);
  eff_recoil_cc->SetMarkerStyle(20);
  eff_recoil_cc->SetMarkerSize(1);
  fitfunc_recoil_cc->SetLineColor(kBlack);
  fitfunc_recoil_cc->SetLineWidth(2);

  TString name (hnum_recoil_cc->GetName());
  name.ReplaceAll("hnum","eff");
  eff_recoil_cc->SetName(name.Data());

  TEfficiency* eff_recoil_cf = new TEfficiency(*hnum_recoil_cf,*hden_recoil_cf);
  eff_recoil_cf->SetMarkerColor(kBlue);
  eff_recoil_cf->SetLineColor(kBlue);
  eff_recoil_cf->SetMarkerStyle(20);
  eff_recoil_cf->SetMarkerSize(1);
  fitfunc_recoil_cf->SetLineColor(kBlue);
  fitfunc_recoil_cf->SetLineWidth(2);

  name = TString(hnum_recoil_cf->GetName());
  name.ReplaceAll("hnum","eff");
  eff_recoil_cf->SetName(name.Data());

  TEfficiency* eff_recoil_fc = new TEfficiency(*hnum_recoil_fc,*hden_recoil_fc);
  eff_recoil_fc->SetMarkerColor(kRed);
  eff_recoil_fc->SetLineColor(kRed);
  eff_recoil_fc->SetMarkerStyle(20);
  eff_recoil_fc->SetMarkerSize(1);
  fitfunc_recoil_fc->SetLineColor(kRed);
  fitfunc_recoil_fc->SetLineWidth(2);

  name = TString(hnum_recoil_cf->GetName());
  name.ReplaceAll("hnum","eff");
  eff_recoil_fc->SetName(name.Data());

  vector<TEfficiency*> eff_recoil; 
  eff_recoil.push_back(eff_recoil_cc); 
  eff_recoil.push_back(eff_recoil_cf);
  eff_recoil.push_back(eff_recoil_fc);
  vector<TF1*> fit_recoil; 
  fit_recoil.push_back(fitfunc_recoil_cc); 
  fit_recoil.push_back(fitfunc_recoil_cf);
  fit_recoil.push_back(fitfunc_recoil_fc);
  vector<string> label; 
  label.push_back("central-central"); 
  label.push_back("central-forward");
  label.push_back("forward-central");
  if(sample == Sample::wmn)
    plotTurnOn(canvas,eff_recoil,fit_recoil,"Recoil [GeV]","recoil",outputDIR,label);
  else if(sample == Sample::zmm)
    plotTurnOn(canvas,eff_recoil,fit_recoil,"Recoil [GeV]","recoil",outputDIR,label);

  ////////
  TEfficiency* eff_mjj_cc = new TEfficiency(*hnum_mjj_cc,*hden_mjj_cc);
  eff_mjj_cc->SetMarkerColor(kBlack);
  eff_mjj_cc->SetLineColor(kBlack);
  eff_mjj_cc->SetMarkerStyle(20);
  eff_mjj_cc->SetMarkerSize(1);
  TF1* fitfunc_mjj_cc = NULL;

  TEfficiency* eff_mjj_cf = new TEfficiency(*hnum_mjj_cf,*hden_mjj_cf);
  eff_mjj_cf->SetMarkerColor(kBlue);
  eff_mjj_cf->SetLineColor(kBlue);
  eff_mjj_cf->SetMarkerStyle(20);
  eff_mjj_cf->SetMarkerSize(1);
  TF1* fitfunc_mjj_cf = NULL;

  TEfficiency* eff_mjj_fc = new TEfficiency(*hnum_mjj_fc,*hden_mjj_fc);
  eff_mjj_fc->SetMarkerColor(kRed);
  eff_mjj_fc->SetLineColor(kRed);
  eff_mjj_fc->SetMarkerStyle(20);
  eff_mjj_fc->SetMarkerSize(1);
  TF1* fitfunc_mjj_fc = NULL;

  name = TString(hnum_mjj_cc->GetName());
  name.ReplaceAll("hnum","eff");
  eff_mjj_cc->SetName(name.Data());

  name = TString(hnum_mjj_cf->GetName());
  name.ReplaceAll("hnum","eff");
  eff_mjj_cf->SetName(name.Data());

  name = TString(hnum_mjj_fc->GetName());
  name.ReplaceAll("hnum","eff");
  eff_mjj_fc->SetName(name.Data());

  vector<TEfficiency*> eff_mjj; 
  eff_mjj.push_back(eff_mjj_cc); 
  eff_mjj.push_back(eff_mjj_cf);
  eff_mjj.push_back(eff_mjj_fc);
  vector<TF1*> fitfunc_mjj; 
  fitfunc_mjj.push_back(fitfunc_mjj_cc); 
  fitfunc_mjj.push_back(fitfunc_mjj_cf);
  fitfunc_mjj.push_back(fitfunc_mjj_fc);

  if(sample == Sample::wmn)
    plotTurnOn(canvas,eff_mjj,fitfunc_mjj,"M_{jj} [GeV]","mjj",outputDIR,label);
  else if(sample == Sample::zmm)
    plotTurnOn(canvas,eff_mjj,fitfunc_mjj,"M_{jj} [GeV]","mjj",outputDIR,label);

  ///////---    
  vector<TEfficiency*> eff_recoil_cc_vs_detajj;
  int icolor = 1;
  for(size_t ibin = 0; ibin < cuts_detajj_cc.size()-1; ibin++){
    eff_recoil_cc_vs_detajj.push_back(new TEfficiency(*hnum_recoil_cc_vs_detajj.at(ibin),*hden_recoil_cc_vs_detajj.at(ibin)));    
    if(icolor == 5) icolor++;
    if(icolor == 10) icolor++;
    eff_recoil_cc_vs_detajj.back()->SetMarkerColor(icolor);
    eff_recoil_cc_vs_detajj.back()->SetLineColor(icolor);
    eff_recoil_cc_vs_detajj.back()->SetMarkerStyle(20);
    eff_recoil_cc_vs_detajj.back()->SetMarkerSize(1);
    fitfunc_recoil_cc_vs_detajj.at(ibin)->SetLineColor(icolor);
    fitfunc_recoil_cc_vs_detajj.at(ibin)->SetLineWidth(2);

    name = TString(hnum_recoil_cc_vs_detajj.at(ibin)->GetName());
    name.ReplaceAll("hnum","eff");
    eff_recoil_cc_vs_detajj.back()->SetName(name.Data());
    icolor++;
  }
  plotTurnOn(canvas,eff_recoil_cc_vs_detajj,fitfunc_recoil_cc_vs_detajj,"Recoil [GeV]","recoil_cc_vs_detajj",outputDIR,cuts_detajj_cc,"|#Delta#eta_{jj}|");

  ///////---    
  vector<TEfficiency*> eff_recoil_cf_vs_detajj;
  icolor = 1;
  for(size_t ibin = 0; ibin < cuts_detajj_cf.size()-1; ibin++){
    eff_recoil_cf_vs_detajj.push_back(new TEfficiency(*hnum_recoil_cf_vs_detajj.at(ibin),*hden_recoil_cf_vs_detajj.at(ibin)));    
    if(icolor == 3) icolor++;
    if(icolor == 5) icolor++;
    if(icolor == 10) icolor++;
    eff_recoil_cf_vs_detajj.back()->SetMarkerColor(icolor);
    eff_recoil_cf_vs_detajj.back()->SetLineColor(icolor);
    eff_recoil_cf_vs_detajj.back()->SetMarkerStyle(20);
    eff_recoil_cf_vs_detajj.back()->SetMarkerSize(1);
    fitfunc_recoil_cf_vs_detajj.at(ibin)->SetLineColor(icolor);
    fitfunc_recoil_cf_vs_detajj.at(ibin)->SetLineWidth(2);

    name = TString(hnum_recoil_cf_vs_detajj.at(ibin)->GetName());
    name.ReplaceAll("hnum","eff");
    eff_recoil_cf_vs_detajj.back()->SetName(name.Data());
    icolor++;
  }
  plotTurnOn(canvas,eff_recoil_cf_vs_detajj,fitfunc_recoil_cf_vs_detajj,"Recoil [GeV]","recoil_cf_vs_detajj",outputDIR,cuts_detajj_cf,"|#Delta#eta_{jj}|");

  ///////---    
  vector<TEfficiency*> eff_recoil_fc_vs_detajj;
  icolor = 1;
  for(size_t ibin = 0; ibin < cuts_detajj_fc.size()-1; ibin++){
    eff_recoil_fc_vs_detajj.push_back(new TEfficiency(*hnum_recoil_fc_vs_detajj.at(ibin),*hden_recoil_fc_vs_detajj.at(ibin)));    
    if(icolor == 3) icolor++;
    if(icolor == 5) icolor++;
    if(icolor == 10) icolor++;
    eff_recoil_fc_vs_detajj.back()->SetMarkerColor(icolor);
    eff_recoil_fc_vs_detajj.back()->SetLineColor(icolor);
    eff_recoil_fc_vs_detajj.back()->SetMarkerStyle(20);
    eff_recoil_fc_vs_detajj.back()->SetMarkerSize(1);
    fitfunc_recoil_fc_vs_detajj.at(ibin)->SetLineColor(icolor);
    fitfunc_recoil_fc_vs_detajj.at(ibin)->SetLineWidth(2);

    name = TString(hnum_recoil_fc_vs_detajj.at(ibin)->GetName());
    name.ReplaceAll("hnum","eff");
    eff_recoil_fc_vs_detajj.back()->SetName(name.Data());
    icolor++;
  }
  plotTurnOn(canvas,eff_recoil_fc_vs_detajj,fitfunc_recoil_fc_vs_detajj,"Recoil [GeV]","recoil_fc_vs_detajj",outputDIR,cuts_detajj_fc,"|#Delta#eta_{jj}|");

  ////////
  TEfficiency* eff_ht_cc = new TEfficiency(*hnum_ht_cc,*hden_ht_cc);
  eff_ht_cc->SetMarkerColor(kBlack);
  eff_ht_cc->SetLineColor(kBlack);
  eff_ht_cc->SetMarkerStyle(20);
  eff_ht_cc->SetMarkerSize(1);
  TF1* fitfunc_ht_cc = NULL;

  TEfficiency* eff_ht_cf = new TEfficiency(*hnum_ht_cf,*hden_ht_cf);
  eff_ht_cf->SetMarkerColor(kBlue);
  eff_ht_cf->SetLineColor(kBlue);
  eff_ht_cf->SetMarkerStyle(20);
  eff_ht_cf->SetMarkerSize(1);
  TF1* fitfunc_ht_cf = NULL;

  TEfficiency* eff_ht_fc = new TEfficiency(*hnum_ht_fc,*hden_ht_fc);
  eff_ht_fc->SetMarkerColor(kRed);
  eff_ht_fc->SetLineColor(kRed);
  eff_ht_fc->SetMarkerStyle(20);
  eff_ht_fc->SetMarkerSize(1);
  TF1* fitfunc_ht_fc = NULL;

  name = TString(hnum_ht_cc->GetName());
  name.ReplaceAll("hnum","eff");
  eff_ht_cc->SetName(name.Data());
  
  name = TString(hnum_ht_cf->GetName());
  name.ReplaceAll("hnum","eff");
  eff_ht_cf->SetName(name.Data());

  name = TString(hnum_ht_fc->GetName());
  name.ReplaceAll("hnum","eff");
  eff_ht_fc->SetName(name.Data());
  
  vector<TEfficiency*> eff_ht; 
  eff_ht.push_back(eff_ht_cc); 
  eff_ht.push_back(eff_ht_cf);
  eff_ht.push_back(eff_ht_fc);

  vector<TF1*> fitfunc_ht; 
  fitfunc_ht.push_back(fitfunc_ht_cc); 
  fitfunc_ht.push_back(fitfunc_ht_cf);
  fitfunc_ht.push_back(fitfunc_ht_fc);

  if(sample == Sample::wmn)
    plotTurnOn(canvas,eff_ht,fitfunc_ht,"H_{T} [GeV]","ht",outputDIR,label);
  else if(sample == Sample::zmm)
    plotTurnOn(canvas,eff_ht,fitfunc_ht,"H_{T} [GeV]","ht",outputDIR,label);

  // fill output file
  outputFile->cd();
  outputFile->mkdir("efficiency_cc");
  outputFile->cd("efficiency_cc");
  eff_recoil_cc->Write();
  fitfunc_recoil_cc->Write();
  eff_mjj_cc->Write();
  eff_ht_cc->Write();
  for(auto eff : eff_recoil_cc_vs_detajj) eff->Write();
  for(auto func : fitfunc_recoil_cc_vs_detajj) func->Write();

  outputFile->cd();
  outputFile->mkdir("efficiency_cf");
  outputFile->cd("efficiency_cf");
  eff_recoil_cf->Write();
  fitfunc_recoil_cf->Write();
  eff_mjj_cf->Write();
  eff_ht_cf->Write();
  for(auto eff: eff_recoil_cf_vs_detajj) eff->Write();
  for(auto fit: fitfunc_recoil_cf_vs_detajj) fit->Write();

  outputFile->cd();
  outputFile->mkdir("efficiency_fc");
  outputFile->cd("efficiency_fc");
  eff_recoil_fc->Write();
  fitfunc_recoil_fc->Write();
  eff_mjj_fc->Write();
  eff_ht_cf->Write();
  for(auto eff: eff_recoil_fc_vs_detajj) eff->Write();
  for(auto fit: fitfunc_recoil_fc_vs_detajj) fit->Write();

  outputFile->Close();
}

////// -------------
void plotTurnOn(TCanvas* canvas,
                vector<TEfficiency* > eff,
                vector<TF1*> fitfunc,
                const string  & axisLabel,
                const string  & postfix,
                const string  & outputDIR,
		const vector<float> & binning,
		const string  & binningVar){

  const TH1* histo_temp = eff.at(0)->GetPassedHistogram();

  
  TH1* frame = canvas->DrawFrame(histo_temp->GetXaxis()->GetXmin(),0.,histo_temp->GetXaxis()->GetXmax(), 1.1, "");
  frame->GetXaxis()->SetTitle(axisLabel.c_str());  
  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetYaxis()->SetLabelSize(0.85*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetLabelSize(0.85*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetTitleSize(0.95*frame->GetYaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleSize(0.95*frame->GetXaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(1.15);

  vector<TGraphAsymmErrors*> graph;
  for(auto efficiency : eff)
    graph.push_back(efficiency->CreateGraph());
  
  // make the fit
  vector<TH1F*> error_band;
  for(size_t iobj = 0; iobj < graph.size(); iobj++){
    if(fitfunc.size() != 0){
      TFitResultPtr fitResult = graph.at(iobj)->Fit(fitfunc.at(iobj),"RSM");
      int npoints       = 350;                                                                                                                                                                   
      error_band.push_back(new TH1F(Form("%s_error_band",fitfunc.at(iobj)->GetName()),"",npoints,fitfunc.at(iobj)->GetXaxis()->GetXmin(),fitfunc.at(iobj)->GetXaxis()->GetXmax()));
      if(drawUncertaintyBand)
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(error_band.back(),0.68);
    }
  }
  
  canvas->SetRightMargin(0.075);
  canvas->SetTopMargin(0.06);
  canvas->Draw();
  canvas->cd();
  frame->Draw();

  if(drawUncertaintyBand and error_band.size() !=0){
    for(size_t ihist = 0; ihist < error_band.size(); ihist++){
      error_band.at(ihist)->SetFillColor(graph.at(ihist)->GetMarkerColor());
      error_band.at(ihist)->SetFillStyle(3001);
      error_band.at(ihist)->Draw("e4 same");  
    }
  }

  TLegend leg (0.35,0.25,0.65,0.45);
  leg.SetFillColor(0);
  leg.SetFillStyle(1001);
  leg.SetBorderSize(0);
  
  for(size_t iobj = 0; iobj < graph.size(); iobj++){
    graph.at(iobj)->Draw("E1PSAME");
    if(fitfunc.size() != 0 and fitfunc.at(iobj) != NULL and fitfunc.at(iobj) != 0){
      fitfunc.at(iobj)->Draw("SAME");  
    }
    leg.AddEntry(graph.at(iobj),Form("%.1f < %s < %.1f ",binning.at(iobj),binningVar.c_str(),binning.at(iobj+1)),"EP");
  }
  leg.Draw("same");
  
  canvas->RedrawAxis();
  CMS_lumi(canvas,"35.9",true);

  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");
}

////// -------------
void plotTurnOn(TCanvas* canvas,
                vector<TEfficiency* > eff,
                vector<TF1*> fitfunc,
                const string  & axisLabel,
                const string  & postfix,
                const string  & outputDIR,
		const vector<string> & label){

  const TH1* histo_temp = eff.at(0)->GetPassedHistogram();
  
  TH1* frame = canvas->DrawFrame(histo_temp->GetXaxis()->GetXmin(),0.,histo_temp->GetXaxis()->GetXmax(), 1.1, "");
  frame->GetXaxis()->SetTitle(axisLabel.c_str());  
  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetYaxis()->SetLabelSize(0.85*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetLabelSize(0.85*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetTitleSize(0.95*frame->GetYaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleSize(0.95*frame->GetXaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(1.15);

  vector<TGraphAsymmErrors*> graph;
  for(auto efficiency : eff)
    graph.push_back(efficiency->CreateGraph());
  
  // make the fit
  vector<TH1F*> error_band;
  for(size_t iobj = 0; iobj < eff.size(); iobj++){
    if(fitfunc.size() != 0 and fitfunc.at(iobj) != 0 and fitfunc.at(iobj) != NULL){
      TFitResultPtr fitResult = graph.at(iobj)->Fit(fitfunc.at(iobj),"RSM");
      int npoints       = 350;                                                                                                                                                                   
      error_band.push_back(new TH1F(Form("%s_error_band",fitfunc.at(iobj)->GetName()),"",npoints,fitfunc.at(iobj)->GetXaxis()->GetXmin(),fitfunc.at(iobj)->GetXaxis()->GetXmax()));
      if(drawUncertaintyBand)
	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(error_band.back(),0.68);
    }
  }
  
  canvas->SetRightMargin(0.075);
  canvas->SetTopMargin(0.06);
  canvas->Draw();
  canvas->cd();
  frame->Draw();

  if(drawUncertaintyBand and error_band.size() !=0){
    for(size_t ihist = 0; ihist < error_band.size(); ihist++){
      error_band.at(ihist)->SetFillColor(graph.at(ihist)->GetMarkerColor());
      error_band.at(ihist)->SetFillStyle(3001);
      error_band.at(ihist)->Draw("e4 same");  
    }
  }

  TLegend leg (0.35,0.25,0.65,0.45);
  leg.SetFillColor(0);
  leg.SetFillStyle(1001);
  leg.SetBorderSize(0);

  for(size_t iobj = 0; iobj < graph.size(); iobj++){
    graph.at(iobj)->Draw("E1PSAME");
    if(fitfunc.at(iobj) != NULL and fitfunc.at(iobj) != 0){
      fitfunc.at(iobj)->Draw("SAME");  
    }
    leg.AddEntry(graph.at(iobj),label.at(iobj).c_str(),"EP");
  }
  leg.Draw("same");
  
  canvas->RedrawAxis();
  CMS_lumi(canvas,"35.9",true);

  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");

}



