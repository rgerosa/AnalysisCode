#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "../CMS_lumi.h"

static float luminosity = 35.9;
static float jetmetdphi_cut = 0.5;
static float jetpt1_cut  = 80;
static float jetpt2_cut  = 40;
static float jeteta1_cut = 4.7;
static float jeteta2_cut = 4.7;
static float met_cut     = 250;
static float dphijj_relaxed_cut  = 1.5;
static float detajj_relaxed_cut  = 1.0;
static float mjj_relaxed_cut     = 200;
static float dphijj_cut  = 1.5;
static float detajj_cut  = 4.0;
static float mjj_cut     = 1300;

enum class Category {VBFrelaxed,VBF};

void runAnalysisSelection(TChain* chain, const Category & category, const double & xsec, TH1F* passing){

  ////////////////////////  
  TTreeReader reader(chain);
  TTreeReaderValue<float> wgt     (reader,"wgt");
  TTreeReaderValue<double> wgtsum (reader,"wgtsum");
  TTreeReaderValue<float> wgtpileup (reader,"wgtpileup");
  TTreeReaderValue<float> wgtbtag (reader,"wgtbtag");

  /////////////// triggers 
  TTreeReaderValue<UChar_t> hltm90      (reader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100     (reader,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110     (reader,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120     (reader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm90    (reader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hltmwm100   (reader,"hltmetwithmu100");
  TTreeReaderValue<UChar_t> hltmwm110   (reader,"hltmetwithmu110");
  TTreeReaderValue<UChar_t> hltmwm120   (reader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170   (reader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300   (reader,"hltmetwithmu300");

  /////////////// met filters
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");

  // njets
  TTreeReaderValue<unsigned int> nincjets    (reader,"njetsinc");
  TTreeReaderValue<unsigned int> nbjets      (reader,"nbjets");

  // AK4 jets
  TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
  TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
  TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
  TTreeReaderValue<vector<float> > chfrac  (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac  (reader,"combinejetNHfrac");

  // MET
  TTreeReaderValue<float> met  (reader,"t1pfmet");
  TTreeReaderValue<float> metcalo (reader,"calomet");
  TTreeReaderValue<float> metphi  (reader,"t1pfmetphi");

  // Dphi
  TTreeReaderValue<float> jmmdphi (reader,"incjetmumetdphimin4");
  TTreeReaderValue<unsigned int> nmuons     (reader,"nmuons");
  TTreeReaderValue<unsigned int> nelectrons (reader,"nelectrons");
  TTreeReaderValue<unsigned int> ntaus      (reader,"ntausold");
  TTreeReaderValue<unsigned int> nphotons   (reader,"nphotons");

  TFile* triggerfile_MET_wmn = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_VBF_htmiss_Wmn.root");
 
  vector<TF1*> triggermet_func_binned_Wmn;
  triggermet_func_binned_Wmn.push_back((TF1*) triggerfile_MET_wmn->Get("f_eff"));
 
  double event_counter = 0;

  ////////////////////
  while(reader.Next()){    

    // trigger
    Double_t hlt =  *hltm90+*hltm100+*hltm110+*hltm120+*hltmwm90+*hltmwm100+*hltmwm110+*hltmwm120+*hltmwm170+*hltmwm300;
    if(hlt == 0) continue;

    // met filters
    if(*fhbhe == 0 or *fhbiso == 0 or *feeb == 0 or *fetp == 0 or *fvtx == 0 or *fcsc == 0 or *fbadmu == 0 or *fbadch == 0) continue;
    
    // apply vetos
    if(*nmuons != 0) continue;
    if(*nelectrons != 0) continue;
    if(*nphotons != 0) continue;
    if(*ntaus != 0) continue;
    if(*nbjets != 0) continue;

    // jet requirement
    if(jetpt->size() < 2) continue;
    if(*nincjets < 2) continue;

    // recoil
    if(*met < met_cut) continue;    
    if((*met-*metcalo)/(*met) > 0.5) continue;

    //------ Jets requirement
    if(jetpt->at(0) < jetpt1_cut) continue;
    if(jetpt->at(1) < jetpt2_cut) continue;
    if(fabs(jeteta->at(0)) > jeteta1_cut) continue;
    if(fabs(jeteta->at(1)) > jeteta2_cut) continue;
    if(*jmmdphi < jetmetdphi_cut) continue;
    if(fabs(jeteta->at(0)) < 2.4 and chfrac->at(0) < 0.1) continue;
    if(fabs(jeteta->at(0)) < 2.4 and nhfrac->at(0) > 0.8) continue;
    if(fabs(jeteta->at(0)) > 3 and fabs(jeteta->at(1)) > 3) continue;
    if(jeteta->at(0)*jeteta->at(1) > 0) continue;

    //------
    TLorentzVector jet1,jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
    
    if(category == Category::VBF){
      if(fabs(jet1.DeltaPhi(jet2)) > dphijj_cut) continue;
      if(fabs(jet1.Eta()-jet2.Eta()) < detajj_cut) continue;
      if((jet1+jet2).M() < mjj_cut) continue;
    }
    else if(category == Category::VBFrelaxed){
      if(fabs(jet1.DeltaPhi(jet2)) > dphijj_relaxed_cut) continue;
      if(fabs(jet1.Eta()-jet2.Eta()) < detajj_relaxed_cut) continue;
      if((jet1+jet2).M() < mjj_relaxed_cut) continue;
    }

    // pileup weight
    double puwgt = 1;
    if(fabs(*wgtpileup) > 100)  puwgt = 1;
    else if(fabs(*wgtpileup) < 0.01) puwgt = 1;
    else puwgt = *wgtpileup;

    // b-tagging
    double btagw = *wgtbtag;

    // trigger
    double sftrig = 1;
    double htmiss = 0;
    TLorentzVector jet4V_total;
    for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
      if(jetpt->at(ijet) < 30) continue;
      if(fabs(jeteta->at(ijet)) > 3) continue;
      TLorentzVector jet; jet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
      jet4V_total += jet;
    }

    htmiss = jet4V_total.Pt();
    sftrig = triggermet_func_binned_Wmn.at(0)->Eval(htmiss);
    passing->Fill(*met,(xsec)*(luminosity)*(*wgt)*(sftrig)*(btagw)*(puwgt)/(*wgtsum));
  }      
  
}


//////////////////////////
void makeHiggsSignalEfficiency(string inputDirectory, Category category, string outputDIR, bool addGGHEfficiency = false){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  double total_qqH_mh110,total_qqH_mh150, total_qqH_mh200, total_qqH_mh300, total_qqH_mh400, total_qqH_mh500, total_qqH_mh600, total_qqH_mh800, total_qqH_mh1000;
  double total_ggH_mh110,total_ggH_mh150, total_ggH_mh200, total_ggH_mh300, total_ggH_mh400, total_ggH_mh500, total_ggH_mh600, total_ggH_mh800, total_ggH_mh1000;

  double xsec_qqH_mh110 = 4.434E+00;
  double xsec_qqH_mh150 = 3.239E+00;
  double xsec_qqH_mh200 = 2.282E+00;
  double xsec_qqH_mh300 = 1.256E+00;
  double xsec_qqH_mh400 = 7.580E-01;
  double xsec_qqH_mh500 = 4.872E-01;
  double xsec_qqH_mh600 = 3.274E-01;
  double xsec_qqH_mh800 = 1.622E-01;
  double xsec_qqH_mh1000 = 8.732E-02;

  total_qqH_mh110 = 1000*luminosity*4.434E+00;
  total_qqH_mh150 = 1000*luminosity*3.239E+00;
  total_qqH_mh200 = 1000*luminosity*2.282E+00;
  total_qqH_mh300 = 1000*luminosity*1.256E+00;
  total_qqH_mh400 = 1000*luminosity*7.580E-01;
  total_qqH_mh500 = 1000*luminosity*4.872E-01;
  total_qqH_mh600 = 1000*luminosity*3.274E-01;
  total_qqH_mh800 = 1000*luminosity*1.622E-01;
  total_qqH_mh1000 = 1000*luminosity*8.732E-02;

  double xsec_ggH_mh110 = 5.790E+01;
  double xsec_ggH_mh150 = 3.129E+01;
  double xsec_ggH_mh200 = 1.694E+01;
  double xsec_ggH_mh300 = 6.590E+00;
  double xsec_ggH_mh400 = 3.160E+00;
  double xsec_ggH_mh500 = 1.709E+00;
  double xsec_ggH_mh600 = 1.001E+00;
  double xsec_ggH_mh800 = 4.015E-01;
  double xsec_ggH_mh1000 = 1.845E-01;
  
  if(addGGHEfficiency){
    total_ggH_mh110 = 1000*luminosity*5.790E+01;
    total_ggH_mh150 = 1000*luminosity*3.129E+01;
    total_ggH_mh200 = 1000*luminosity*1.694E+01;
    total_ggH_mh300 = 1000*luminosity*6.590E+00;
    total_ggH_mh400 = 1000*luminosity*3.160E+00;
    total_ggH_mh500 = 1000*luminosity*1.709E+00;
    total_ggH_mh600 = 1000*luminosity*1.001E+00;
    total_ggH_mh800 = 1000*luminosity*4.015E-01;
    total_ggH_mh1000 = 1000*luminosity*1.845E-01;
  }
  
  if(category != Category::VBF and category != Category::VBFrelaxed){
    cout<<"Category not adapt for this study --> return"<<endl;
    return;
  }

  ////// ------
  TChain* qqH_mh110 =  new TChain("tree/tree");
  TChain* qqH_mh150 =  new TChain("tree/tree");
  TChain* qqH_mh200 =  new TChain("tree/tree");
  TChain* qqH_mh300 =  new TChain("tree/tree");
  TChain* qqH_mh400 =  new TChain("tree/tree");
  TChain* qqH_mh500 =  new TChain("tree/tree");
  TChain* qqH_mh600 =  new TChain("tree/tree");
  TChain* qqH_mh800 =  new TChain("tree/tree");
  TChain* qqH_mh1000 =  new TChain("tree/tree");
  
  qqH_mh110->Add((inputDirectory+"/HiggsInvisible/sigfilter/*VBF*110*root").c_str());
  qqH_mh150->Add((inputDirectory+"/HiggsInvisible/sigfilter/*VBF*150*root").c_str());
  qqH_mh200->Add((inputDirectory+"/HiggsInvisible/sigfilter/*VBF*200*root").c_str());
  qqH_mh300->Add((inputDirectory+"/HiggsInvisible/sigfilter/*VBF*300*root").c_str());
  qqH_mh400->Add((inputDirectory+"/HiggsInvisible/sigfilter/*VBF*400*root").c_str());
  qqH_mh500->Add((inputDirectory+"/HiggsInvisible/sigfilter/*VBF*500*root").c_str());
  qqH_mh600->Add((inputDirectory+"/HiggsInvisible/sigfilter/*VBF*600*root").c_str());
  qqH_mh800->Add((inputDirectory+"/HiggsInvisible/sigfilter/*VBF*800*root").c_str());
  qqH_mh1000->Add((inputDirectory+"/HiggsInvisible/sigfilter/*VBF*1000*root").c_str());

  TChain* ggH_mh110 = NULL;
  TChain* ggH_mh150 = NULL;
  TChain* ggH_mh200 = NULL;
  TChain* ggH_mh300 = NULL;
  TChain* ggH_mh400 = NULL;
  TChain* ggH_mh500 = NULL;
  TChain* ggH_mh600 = NULL;
  TChain* ggH_mh800 = NULL;
  TChain* ggH_mh1000 = NULL;

  if(addGGHEfficiency){

    ggH_mh110 =  new TChain("tree/tree");
    ggH_mh150 =  new TChain("tree/tree");
    ggH_mh200 =  new TChain("tree/tree");
    ggH_mh300 =  new TChain("tree/tree");
    ggH_mh400 =  new TChain("tree/tree");
    ggH_mh500 =  new TChain("tree/tree");
    ggH_mh600 =  new TChain("tree/tree");
    ggH_mh800 =  new TChain("tree/tree");
    ggH_mh1000 =  new TChain("tree/tree");
    
    ggH_mh110->Add((inputDirectory+"/HiggsInvisible/sigfilter/*GluGlu*110*root").c_str());
    ggH_mh150->Add((inputDirectory+"/HiggsInvisible/sigfilter/*GluGlu*150*root").c_str());
    ggH_mh200->Add((inputDirectory+"/HiggsInvisible/sigfilter/*GluGlu*200*root").c_str());
    ggH_mh300->Add((inputDirectory+"/HiggsInvisible/sigfilter/*GluGlu*300*root").c_str());
    ggH_mh400->Add((inputDirectory+"/HiggsInvisible/sigfilter/*GluGlu*400*root").c_str());
    ggH_mh500->Add((inputDirectory+"/HiggsInvisible/sigfilter/*GluGlu*500*root").c_str());
    ggH_mh600->Add((inputDirectory+"/HiggsInvisible/sigfilter/*GluGlu*600*root").c_str());
    ggH_mh800->Add((inputDirectory+"/HiggsInvisible/sigfilter/*GluGlu*800*root").c_str());
    ggH_mh1000->Add((inputDirectory+"/HiggsInvisible/sigfilter/*GluGlu*1000*root").c_str());
  }

  TH1F* passing_qqH_mh110 = new TH1F("passing_qqH_mh110","",1,0,10000);
  TH1F* passing_qqH_mh150 = new TH1F("passing_qqH_mh150","",1,0,10000);
  TH1F* passing_qqH_mh200 = new TH1F("passing_qqH_mh200","",1,0,10000);
  TH1F* passing_qqH_mh300 = new TH1F("passing_qqH_mh300","",1,0,10000);
  TH1F* passing_qqH_mh400 = new TH1F("passing_qqH_mh400","",1,0,10000);
  TH1F* passing_qqH_mh500 = new TH1F("passing_qqH_mh500","",1,0,10000);
  TH1F* passing_qqH_mh600 = new TH1F("passing_qqH_mh600","",1,0,10000);
  TH1F* passing_qqH_mh800 = new TH1F("passing_qqH_mh800","",1,0,10000);
  TH1F* passing_qqH_mh1000 = new TH1F("passing_qqH_mh1000","",1,0,10000);

  passing_qqH_mh110->Sumw2();
  passing_qqH_mh150->Sumw2();
  passing_qqH_mh200->Sumw2();
  passing_qqH_mh300->Sumw2();
  passing_qqH_mh400->Sumw2();
  passing_qqH_mh500->Sumw2();
  passing_qqH_mh600->Sumw2();
  passing_qqH_mh800->Sumw2();
  passing_qqH_mh1000->Sumw2();
  
  runAnalysisSelection(qqH_mh110,category,xsec_qqH_mh110*1000,passing_qqH_mh110);
  runAnalysisSelection(qqH_mh150,category,xsec_qqH_mh150*1000,passing_qqH_mh150);
  runAnalysisSelection(qqH_mh200,category,xsec_qqH_mh200*1000,passing_qqH_mh200);
  runAnalysisSelection(qqH_mh300,category,xsec_qqH_mh300*1000,passing_qqH_mh300);
  runAnalysisSelection(qqH_mh400,category,xsec_qqH_mh400*1000,passing_qqH_mh400);
  runAnalysisSelection(qqH_mh500,category,xsec_qqH_mh500*1000,passing_qqH_mh500);
  runAnalysisSelection(qqH_mh600,category,xsec_qqH_mh600*1000,passing_qqH_mh600);
  runAnalysisSelection(qqH_mh800,category,xsec_qqH_mh800*1000,passing_qqH_mh800);
  runAnalysisSelection(qqH_mh1000,category,xsec_qqH_mh1000*1000,passing_qqH_mh1000);
  
  TH1F* efficiency_qqH_mh110 = (TH1F*) passing_qqH_mh110->Clone("efficiency_qqH_mh110");
  efficiency_qqH_mh110->Scale(1./total_qqH_mh110);
  TH1F* efficiency_qqH_mh150 = (TH1F*) passing_qqH_mh150->Clone("efficiency_qqH_mh150");
  efficiency_qqH_mh150->Scale(1./total_qqH_mh150);
  TH1F* efficiency_qqH_mh200 = (TH1F*) passing_qqH_mh200->Clone("efficiency_qqH_mh200");
  efficiency_qqH_mh200->Scale(1./total_qqH_mh200);
  TH1F* efficiency_qqH_mh300 = (TH1F*) passing_qqH_mh300->Clone("efficiency_qqH_mh300");
  efficiency_qqH_mh300->Scale(1./total_qqH_mh300);
  TH1F* efficiency_qqH_mh400 = (TH1F*) passing_qqH_mh400->Clone("efficiency_qqH_mh400");
  efficiency_qqH_mh400->Scale(1./total_qqH_mh400);
  TH1F* efficiency_qqH_mh500 = (TH1F*) passing_qqH_mh500->Clone("efficiency_qqH_mh500");
  efficiency_qqH_mh500->Scale(1./total_qqH_mh500);
  TH1F* efficiency_qqH_mh600 = (TH1F*) passing_qqH_mh600->Clone("efficiency_qqH_mh600");
  efficiency_qqH_mh600->Scale(1./total_qqH_mh600);
  TH1F* efficiency_qqH_mh800 = (TH1F*) passing_qqH_mh800->Clone("efficiency_qqH_mh800");
  efficiency_qqH_mh800->Scale(1./total_qqH_mh800);
  TH1F* efficiency_qqH_mh1000 = (TH1F*) passing_qqH_mh1000->Clone("efficiency_qqH_mh1000");
  efficiency_qqH_mh1000->Scale(1./total_qqH_mh1000);

  cout<<"qqH mH = 110 GeV -->  Total events "<<total_qqH_mh110<<" passing "<<passing_qqH_mh110->GetBinContent(1)<<" eff "<<efficiency_qqH_mh110->GetBinContent(1)<<endl;
  cout<<"qqH mH = 150 GeV -->  Total events "<<total_qqH_mh150<<" passing "<<passing_qqH_mh150->GetBinContent(1)<<" eff "<<efficiency_qqH_mh150->GetBinContent(1)<<endl;
  cout<<"qqH mH = 200 GeV -->  Total events "<<total_qqH_mh200<<" passing "<<passing_qqH_mh200->GetBinContent(1)<<" eff "<<efficiency_qqH_mh200->GetBinContent(1)<<endl;
  cout<<"qqH mH = 300 GeV -->  Total events "<<total_qqH_mh300<<" passing "<<passing_qqH_mh300->GetBinContent(1)<<" eff "<<efficiency_qqH_mh300->GetBinContent(1)<<endl;
  cout<<"qqH mH = 400 GeV -->  Total events "<<total_qqH_mh400<<" passing "<<passing_qqH_mh400->GetBinContent(1)<<" eff "<<efficiency_qqH_mh400->GetBinContent(1)<<endl;
  cout<<"qqH mH = 500 GeV -->  Total events "<<total_qqH_mh500<<" passing "<<passing_qqH_mh500->GetBinContent(1)<<" eff "<<efficiency_qqH_mh500->GetBinContent(1)<<endl;
  cout<<"qqH mH = 600 GeV -->  Total events "<<total_qqH_mh600<<" passing "<<passing_qqH_mh600->GetBinContent(1)<<" eff "<<efficiency_qqH_mh600->GetBinContent(1)<<endl;
  cout<<"qqH mH = 800 GeV -->  Total events "<<total_qqH_mh800<<" passing "<<passing_qqH_mh800->GetBinContent(1)<<" eff "<<efficiency_qqH_mh800->GetBinContent(1)<<endl;
  cout<<"qqH mH = 1000 GeV -->  Total events "<<total_qqH_mh1000<<" passing "<<passing_qqH_mh1000->GetBinContent(1)<<" eff "<<efficiency_qqH_mh1000->GetBinContent(1)<<endl;

  /////

  TH1F* passing_ggH_mh110 = new TH1F("passing_ggH_mh110","",1,0,10000);
  TH1F* passing_ggH_mh150 = new TH1F("passing_ggH_mh150","",1,0,10000);
  TH1F* passing_ggH_mh200 = new TH1F("passing_ggH_mh200","",1,0,10000);
  TH1F* passing_ggH_mh300 = new TH1F("passing_ggH_mh300","",1,0,10000);
  TH1F* passing_ggH_mh400 = new TH1F("passing_ggH_mh400","",1,0,10000);
  TH1F* passing_ggH_mh500 = new TH1F("passing_ggH_mh500","",1,0,10000);
  TH1F* passing_ggH_mh600 = new TH1F("passing_ggH_mh600","",1,0,10000);
  TH1F* passing_ggH_mh800 = new TH1F("passing_ggH_mh800","",1,0,10000);
  TH1F* passing_ggH_mh1000 = new TH1F("passing_ggH_mh1000","",1,0,10000);

  passing_ggH_mh110->Sumw2();
  passing_ggH_mh150->Sumw2();
  passing_ggH_mh200->Sumw2();
  passing_ggH_mh300->Sumw2();
  passing_ggH_mh400->Sumw2();
  passing_ggH_mh500->Sumw2();
  passing_ggH_mh600->Sumw2();
  passing_ggH_mh800->Sumw2();
  passing_ggH_mh1000->Sumw2();
  
  if(addGGHEfficiency){
    runAnalysisSelection(ggH_mh110,category,xsec_ggH_mh110*1000,passing_ggH_mh110);
    runAnalysisSelection(ggH_mh150,category,xsec_ggH_mh150*1000,passing_ggH_mh150);
    runAnalysisSelection(ggH_mh200,category,xsec_ggH_mh200*1000,passing_ggH_mh200);
    runAnalysisSelection(ggH_mh300,category,xsec_ggH_mh300*1000,passing_ggH_mh300);
    runAnalysisSelection(ggH_mh400,category,xsec_ggH_mh400*1000,passing_ggH_mh400);
    runAnalysisSelection(ggH_mh500,category,xsec_ggH_mh500*1000,passing_ggH_mh500);
    runAnalysisSelection(ggH_mh600,category,xsec_ggH_mh600*1000,passing_ggH_mh600);
    runAnalysisSelection(ggH_mh800,category,xsec_ggH_mh800*1000,passing_ggH_mh800);
    runAnalysisSelection(ggH_mh1000,category,xsec_ggH_mh1000*1000,passing_ggH_mh1000);
  }

  TH1F* efficiency_ggH_mh110 = (TH1F*) passing_ggH_mh110->Clone("efficiency_ggH_mh110");
  efficiency_ggH_mh110->Scale(1./total_ggH_mh110);
  TH1F* efficiency_ggH_mh150 = (TH1F*) passing_ggH_mh150->Clone("efficiency_ggH_mh150");
  efficiency_ggH_mh150->Scale(1./total_ggH_mh150);
  TH1F* efficiency_ggH_mh200 = (TH1F*) passing_ggH_mh200->Clone("efficiency_ggH_mh200");
  efficiency_ggH_mh200->Scale(1./total_ggH_mh200);
  TH1F* efficiency_ggH_mh300 = (TH1F*) passing_ggH_mh300->Clone("efficiency_ggH_mh300");
  efficiency_ggH_mh300->Scale(1./total_ggH_mh300);
  TH1F* efficiency_ggH_mh400 = (TH1F*) passing_ggH_mh400->Clone("efficiency_ggH_mh400");
  efficiency_ggH_mh400->Scale(1./total_ggH_mh400);
  TH1F* efficiency_ggH_mh500 = (TH1F*) passing_ggH_mh500->Clone("efficiency_ggH_mh500");
  efficiency_ggH_mh500->Scale(1./total_ggH_mh500);
  TH1F* efficiency_ggH_mh600 = (TH1F*) passing_ggH_mh600->Clone("efficiency_ggH_mh600");
  efficiency_ggH_mh600->Scale(1./total_ggH_mh600);
  TH1F* efficiency_ggH_mh800 = (TH1F*) passing_ggH_mh800->Clone("efficiency_ggH_mh800");
  efficiency_ggH_mh800->Scale(1./total_ggH_mh800);
  TH1F* efficiency_ggH_mh1000 = (TH1F*) passing_ggH_mh1000->Clone("efficiency_ggH_mh1000");
  efficiency_ggH_mh1000->Scale(1./total_ggH_mh1000);

  if(addGGHEfficiency){
    cout<<"ggH mH = 110 GeV -->  Total events "<<total_ggH_mh110<<" passing "<<passing_ggH_mh110->GetBinContent(1)<<" eff "<<efficiency_ggH_mh110->GetBinContent(1)<<endl;
    cout<<"ggH mH = 150 GeV -->  Total events "<<total_ggH_mh150<<" passing "<<passing_ggH_mh150->GetBinContent(1)<<" eff "<<efficiency_ggH_mh150->GetBinContent(1)<<endl;
    cout<<"ggH mH = 200 GeV -->  Total events "<<total_ggH_mh200<<" passing "<<passing_ggH_mh200->GetBinContent(1)<<" eff "<<efficiency_ggH_mh200->GetBinContent(1)<<endl;
    cout<<"ggH mH = 300 GeV -->  Total events "<<total_ggH_mh300<<" passing "<<passing_ggH_mh300->GetBinContent(1)<<" eff "<<efficiency_ggH_mh300->GetBinContent(1)<<endl;
    cout<<"ggH mH = 400 GeV -->  Total events "<<total_ggH_mh400<<" passing "<<passing_ggH_mh400->GetBinContent(1)<<" eff "<<efficiency_ggH_mh400->GetBinContent(1)<<endl;
    cout<<"ggH mH = 500 GeV -->  Total events "<<total_ggH_mh500<<" passing "<<passing_ggH_mh500->GetBinContent(1)<<" eff "<<efficiency_ggH_mh500->GetBinContent(1)<<endl;
    cout<<"ggH mH = 600 GeV -->  Total events "<<total_ggH_mh600<<" passing "<<passing_ggH_mh600->GetBinContent(1)<<" eff "<<efficiency_ggH_mh600->GetBinContent(1)<<endl;
    cout<<"ggH mH = 800 GeV -->  Total events "<<total_ggH_mh800<<" passing "<<passing_ggH_mh800->GetBinContent(1)<<" eff "<<efficiency_ggH_mh800->GetBinContent(1)<<endl;
    cout<<"ggH mH = 1000 GeV -->  Total events "<<total_ggH_mh1000<<" passing "<<passing_ggH_mh1000->GetBinContent(1)<<" eff "<<efficiency_ggH_mh1000->GetBinContent(1)<<endl;    
  }


  ///// -------------------
  vector<float> bins = {90,130,175,250,350,450,550,700,900,1000};

  TH1F* efficiency_qqH = new TH1F("efficiency_qqH","",bins.size()-1,&bins[0]);
  efficiency_qqH->Sumw2();
  
  efficiency_qqH->SetBinContent(1,efficiency_qqH_mh110->GetBinContent(1));
  efficiency_qqH->SetBinContent(2,efficiency_qqH_mh150->GetBinContent(1));
  efficiency_qqH->SetBinContent(3,efficiency_qqH_mh200->GetBinContent(1));
  efficiency_qqH->SetBinContent(4,efficiency_qqH_mh300->GetBinContent(1));
  efficiency_qqH->SetBinContent(5,efficiency_qqH_mh400->GetBinContent(1));
  efficiency_qqH->SetBinContent(6,efficiency_qqH_mh500->GetBinContent(1));
  efficiency_qqH->SetBinContent(7,efficiency_qqH_mh600->GetBinContent(1));
  efficiency_qqH->SetBinContent(8,efficiency_qqH_mh800->GetBinContent(1));
  efficiency_qqH->SetBinContent(9,efficiency_qqH_mh1000->GetBinContent(1));

  efficiency_qqH->SetBinError(1,efficiency_qqH_mh110->GetBinError(1));
  efficiency_qqH->SetBinError(2,efficiency_qqH_mh150->GetBinError(1));
  efficiency_qqH->SetBinError(3,efficiency_qqH_mh200->GetBinError(1));
  efficiency_qqH->SetBinError(4,efficiency_qqH_mh300->GetBinError(1));
  efficiency_qqH->SetBinError(5,efficiency_qqH_mh400->GetBinError(1));
  efficiency_qqH->SetBinError(6,efficiency_qqH_mh500->GetBinError(1));
  efficiency_qqH->SetBinError(7,efficiency_qqH_mh600->GetBinError(1));
  efficiency_qqH->SetBinError(8,efficiency_qqH_mh800->GetBinError(1));
  efficiency_qqH->SetBinError(9,efficiency_qqH_mh1000->GetBinError(1));
  
  efficiency_qqH->SetLineColor(kBlack);
  efficiency_qqH->SetMarkerColor(kRed);
  efficiency_qqH->SetMarkerStyle(20);
  efficiency_qqH->SetMarkerSize(1);
  efficiency_qqH->GetYaxis()->SetTitle("efficiency");
  efficiency_qqH->GetXaxis()->SetTitle("m_{H} [GeV]");
  efficiency_qqH->GetYaxis()->SetTitleOffset(1.15);
  efficiency_qqH->GetXaxis()->SetTitleOffset(1.1);

  ///// ------------------- 
  TH1F* efficiency_ggH = NULL;

  if(addGGHEfficiency){

    efficiency_ggH = new TH1F("efficiency_ggH","",bins.size()-1,&bins[0]);
    efficiency_ggH->Sumw2();

    efficiency_ggH->SetBinContent(1,efficiency_ggH_mh110->GetBinContent(1));
    efficiency_ggH->SetBinContent(2,efficiency_ggH_mh150->GetBinContent(1));
    efficiency_ggH->SetBinContent(3,efficiency_ggH_mh200->GetBinContent(1));
    efficiency_ggH->SetBinContent(4,efficiency_ggH_mh300->GetBinContent(1));
    efficiency_ggH->SetBinContent(5,efficiency_ggH_mh400->GetBinContent(1));
    efficiency_ggH->SetBinContent(6,efficiency_ggH_mh500->GetBinContent(1));
    efficiency_ggH->SetBinContent(7,efficiency_ggH_mh600->GetBinContent(1));
    efficiency_ggH->SetBinContent(8,efficiency_ggH_mh800->GetBinContent(1));
    efficiency_ggH->SetBinContent(9,efficiency_ggH_mh1000->GetBinContent(1));
    
    efficiency_ggH->SetBinError(1,efficiency_ggH_mh110->GetBinError(1));
    efficiency_ggH->SetBinError(2,efficiency_ggH_mh150->GetBinError(1));
    efficiency_ggH->SetBinError(3,efficiency_ggH_mh200->GetBinError(1));
    efficiency_ggH->SetBinError(4,efficiency_ggH_mh300->GetBinError(1));
    efficiency_ggH->SetBinError(5,efficiency_ggH_mh400->GetBinError(1));
    efficiency_ggH->SetBinError(6,efficiency_ggH_mh500->GetBinError(1));
    efficiency_ggH->SetBinError(7,efficiency_ggH_mh600->GetBinError(1));
    efficiency_ggH->SetBinError(8,efficiency_ggH_mh800->GetBinError(1));
    efficiency_ggH->SetBinError(9,efficiency_ggH_mh1000->GetBinError(1));
        
    efficiency_ggH->SetLineColor(kBlue);
    efficiency_ggH->SetMarkerColor(kBlue);
    efficiency_qqH->SetMarkerColor(kBlack);
    
    efficiency_ggH->SetMarkerStyle(20);
    efficiency_ggH->SetMarkerSize(1);
    efficiency_ggH->GetYaxis()->SetTitle("efficiency");
    efficiency_ggH->GetXaxis()->SetTitle("m_{H} [GeV]");
  }

  
  /// ---------------
  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();
  efficiency_qqH->Draw("PE");
  if(addGGHEfficiency)
    efficiency_ggH->Draw("PEsame");
  CMS_lumi(canvas,Form("%.1f",luminosity));

  TLegend leg (0.7,0.3,0.9,0.4);
  if(efficiency_ggH){
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.AddEntry(efficiency_qqH,"VBF SM-like","PEL");
    leg.AddEntry(efficiency_ggH,"ggH SM-like","PEL");
    leg.Draw("same");
  }

  if(category == Category::VBFrelaxed)
    efficiency_qqH->GetYaxis()->SetRangeUser(0.0001,0.1);
  else
    efficiency_qqH->GetYaxis()->SetRangeUser(0.00001,0.1);

  canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/efficiency_vs_mass.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/efficiency_vs_mass.pdf").c_str(),"pdf");

}
