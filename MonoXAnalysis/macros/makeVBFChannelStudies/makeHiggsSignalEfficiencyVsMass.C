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

double runAnalysisSelection(TChain* chain, const Category & category, const double & xsec){

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

  TFile* triggerfile_MET_wmn = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_VBF_Wmn.root");
 
  vector<TF1*> triggermet_func_binned_Wmn_cc;
  vector<TF1*> triggermet_func_binned_Wmn_cf;
 
  
  triggermet_func_binned_Wmn_cc.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cc/fitfunc_recoil_cc_vs_detajj_1.0_1.5"));
  triggermet_func_binned_Wmn_cc.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cc/fitfunc_recoil_cc_vs_detajj_1.5_2.0"));
  triggermet_func_binned_Wmn_cc.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cc/fitfunc_recoil_cc_vs_detajj_2.0_2.5"));
  triggermet_func_binned_Wmn_cc.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cc/fitfunc_recoil_cc_vs_detajj_2.5_3.0"));
  triggermet_func_binned_Wmn_cc.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cc/fitfunc_recoil_cc_vs_detajj_3.0_3.5"));
  triggermet_func_binned_Wmn_cc.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cc/fitfunc_recoil_cc_vs_detajj_3.5_4.0"));
  triggermet_func_binned_Wmn_cc.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cc/fitfunc_recoil_cc_vs_detajj_4.0_5.0"));
  triggermet_func_binned_Wmn_cc.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cc/fitfunc_recoil_cc_vs_detajj_5.0_10.0"));

  triggermet_func_binned_Wmn_cf.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cf/fitfunc_recoil_cf_vs_detajj_3.0_3.5"));
  triggermet_func_binned_Wmn_cf.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cf/fitfunc_recoil_cf_vs_detajj_3.5_4.0"));
  triggermet_func_binned_Wmn_cf.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cf/fitfunc_recoil_cf_vs_detajj_4.0_5.0"));
  triggermet_func_binned_Wmn_cf.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cf/fitfunc_recoil_cf_vs_detajj_5.0_6.0"));
  triggermet_func_binned_Wmn_cf.push_back((TF1*) triggerfile_MET_wmn->Get("efficiency_cf/fitfunc_recoil_cf_vs_detajj_6.0_10.0"));

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
    double detajj = fabs(jeteta->at(0)-jeteta->at(1));
    double pfmet = *met;
    if(fabs(jeteta->at(0)) < 3 and fabs(jeteta->at(1)) < 3){ // central-central combination
      if(detajj > 1.0 && detajj < 1.5) sftrig *= triggermet_func_binned_Wmn_cc.at(0)->Eval(min(pfmet,triggermet_func_binned_Wmn_cc.at(0)->GetXaxis()->GetXmax()));
      else if(detajj > 1.5 && detajj < 2.0) sftrig *= triggermet_func_binned_Wmn_cc.at(1)->Eval(min(pfmet,triggermet_func_binned_Wmn_cc.at(1)->GetXaxis()->GetXmax()));
      else if(detajj > 2.0 && detajj < 2.5) sftrig *= triggermet_func_binned_Wmn_cc.at(2)->Eval(min(pfmet,triggermet_func_binned_Wmn_cc.at(2)->GetXaxis()->GetXmax()));
      else if(detajj > 3.0 && detajj < 3.5) sftrig *= triggermet_func_binned_Wmn_cc.at(3)->Eval(min(pfmet,triggermet_func_binned_Wmn_cc.at(3)->GetXaxis()->GetXmax()));
      else if(detajj > 3.5 && detajj < 4.0) sftrig *= triggermet_func_binned_Wmn_cc.at(4)->Eval(min(pfmet,triggermet_func_binned_Wmn_cc.at(4)->GetXaxis()->GetXmax()));
      else if(detajj > 4.0 && detajj < 5.0) sftrig *= triggermet_func_binned_Wmn_cc.at(5)->Eval(min(pfmet,triggermet_func_binned_Wmn_cc.at(5)->GetXaxis()->GetXmax()));
      else if(detajj > 5.0) sftrig *= triggermet_func_binned_Wmn_cc.at(6)->Eval(min(pfmet,triggermet_func_binned_Wmn_cc.at(6)->GetXaxis()->GetXmax()));
    }
    else if((fabs(jeteta->at(0)) < 3 and fabs(jeteta->at(1)) > 3) or (fabs(jeteta->at(0)) > 3 and fabs(jeteta->at(1)) < 3)){ // central-forward combination
      if(detajj > 3.0 && detajj < 3.5) sftrig *= triggermet_func_binned_Wmn_cf.at(0)->Eval(min(pfmet,triggermet_func_binned_Wmn_cf.at(0)->GetXaxis()->GetXmax()));
      else if(detajj > 3.5 && detajj < 4.0) sftrig *= triggermet_func_binned_Wmn_cf.at(1)->Eval(min(pfmet,triggermet_func_binned_Wmn_cf.at(1)->GetXaxis()->GetXmax()));
      else if(detajj > 4.0 && detajj < 5.0) sftrig *= triggermet_func_binned_Wmn_cf.at(2)->Eval(min(pfmet,triggermet_func_binned_Wmn_cf.at(2)->GetXaxis()->GetXmax()));
      else if(detajj > 5.0 && detajj < 6.0) sftrig *= triggermet_func_binned_Wmn_cf.at(3)->Eval(min(pfmet,triggermet_func_binned_Wmn_cf.at(3)->GetXaxis()->GetXmax()));
      else if(detajj > 6.0) sftrig *= triggermet_func_binned_Wmn_cf.at(4)->Eval(min(pfmet,triggermet_func_binned_Wmn_cf.at(4)->GetXaxis()->GetXmax()));
    }

    event_counter += (xsec)*(luminosity)*(*wgt)*(sftrig)*(btagw)*(puwgt)/(*wgtsum);
  }      
  
  return event_counter;
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

  double selected_qqH_mh110 = runAnalysisSelection(qqH_mh110,category,xsec_qqH_mh110*1000);
  double selected_qqH_mh150 = runAnalysisSelection(qqH_mh150,category,xsec_qqH_mh150*1000);
  double selected_qqH_mh200 = runAnalysisSelection(qqH_mh200,category,xsec_qqH_mh200*1000);
  double selected_qqH_mh300 = runAnalysisSelection(qqH_mh300,category,xsec_qqH_mh300*1000);
  double selected_qqH_mh400 = runAnalysisSelection(qqH_mh400,category,xsec_qqH_mh400*1000);
  double selected_qqH_mh500 = runAnalysisSelection(qqH_mh500,category,xsec_qqH_mh500*1000);
  double selected_qqH_mh600 = runAnalysisSelection(qqH_mh600,category,xsec_qqH_mh600*1000);
  double selected_qqH_mh800 = runAnalysisSelection(qqH_mh800,category,xsec_qqH_mh800*1000);
  double selected_qqH_mh1000 = runAnalysisSelection(qqH_mh1000,category,xsec_qqH_mh1000*1000);

  cout<<"qqH mH = 110 GeV -->  Total events "<<total_qqH_mh110<<" selected "<<selected_qqH_mh110<<" eff "<<selected_qqH_mh110/total_qqH_mh110<<endl;
  cout<<"qqH mH = 150 GeV -->  Total events "<<total_qqH_mh150<<" selected "<<selected_qqH_mh150<<" eff "<<selected_qqH_mh150/total_qqH_mh150<<endl;
  cout<<"qqH mH = 200 GeV -->  Total events "<<total_qqH_mh200<<" selected "<<selected_qqH_mh200<<" eff "<<selected_qqH_mh200/total_qqH_mh200<<endl;
  cout<<"qqH mH = 300 GeV -->  Total events "<<total_qqH_mh300<<" selected "<<selected_qqH_mh300<<" eff "<<selected_qqH_mh300/total_qqH_mh300<<endl;
  cout<<"qqH mH = 400 GeV -->  Total events "<<total_qqH_mh400<<" selected "<<selected_qqH_mh400<<" eff "<<selected_qqH_mh400/total_qqH_mh400<<endl;
  cout<<"qqH mH = 500 GeV -->  Total events "<<total_qqH_mh500<<" selected "<<selected_qqH_mh500<<" eff "<<selected_qqH_mh500/total_qqH_mh500<<endl;
  cout<<"qqH mH = 600 GeV -->  Total events "<<total_qqH_mh600<<" selected "<<selected_qqH_mh600<<" eff "<<selected_qqH_mh600/total_qqH_mh600<<endl;
  cout<<"qqH mH = 800 GeV -->  Total events "<<total_qqH_mh800<<" selected "<<selected_qqH_mh800<<" eff "<<selected_qqH_mh800/total_qqH_mh800<<endl;
  cout<<"qqH mH = 1000 GeV -->  Total events "<<total_qqH_mh1000<<" selected "<<selected_qqH_mh1000<<" eff "<<selected_qqH_mh1000/total_qqH_mh1000<<endl;

  double selected_ggH_mh110, selected_ggH_mh150, selected_ggH_mh200, selected_ggH_mh300, selected_ggH_mh400, selected_ggH_mh500, selected_ggH_mh600, selected_ggH_mh800, selected_ggH_mh1000;
  
  if(addGGHEfficiency){
    selected_ggH_mh110 = runAnalysisSelection(ggH_mh110,category,xsec_ggH_mh110*1000);
    selected_ggH_mh150 = runAnalysisSelection(ggH_mh150,category,xsec_ggH_mh150*1000);
    selected_ggH_mh200 = runAnalysisSelection(ggH_mh200,category,xsec_ggH_mh200*1000);
    selected_ggH_mh300 = runAnalysisSelection(ggH_mh300,category,xsec_ggH_mh300*1000);
    selected_ggH_mh400 = runAnalysisSelection(ggH_mh400,category,xsec_ggH_mh400*1000);
    selected_ggH_mh500 = runAnalysisSelection(ggH_mh500,category,xsec_ggH_mh500*1000);
    selected_ggH_mh600 = runAnalysisSelection(ggH_mh600,category,xsec_ggH_mh600*1000);
    selected_ggH_mh800 = runAnalysisSelection(ggH_mh800,category,xsec_ggH_mh800*1000);
    selected_ggH_mh1000 = runAnalysisSelection(ggH_mh1000,category,xsec_ggH_mh1000*1000);
  }

  cout<<"ggH mH = 110 GeV -->  Total events "<<total_ggH_mh110<<" selected "<<selected_ggH_mh110<<" eff "<<selected_ggH_mh110/total_ggH_mh110<<endl;
  cout<<"ggH mH = 150 GeV -->  Total events "<<total_ggH_mh150<<" selected "<<selected_ggH_mh150<<" eff "<<selected_ggH_mh150/total_ggH_mh150<<endl;
  cout<<"ggH mH = 200 GeV -->  Total events "<<total_ggH_mh200<<" selected "<<selected_ggH_mh200<<" eff "<<selected_ggH_mh200/total_ggH_mh200<<endl;
  cout<<"ggH mH = 300 GeV -->  Total events "<<total_ggH_mh300<<" selected "<<selected_ggH_mh300<<" eff "<<selected_ggH_mh300/total_ggH_mh300<<endl;
  cout<<"ggH mH = 400 GeV -->  Total events "<<total_ggH_mh400<<" selected "<<selected_ggH_mh400<<" eff "<<selected_ggH_mh400/total_ggH_mh400<<endl;
  cout<<"ggH mH = 500 GeV -->  Total events "<<total_ggH_mh500<<" selected "<<selected_ggH_mh500<<" eff "<<selected_ggH_mh500/total_ggH_mh500<<endl;
  cout<<"ggH mH = 600 GeV -->  Total events "<<total_ggH_mh600<<" selected "<<selected_ggH_mh600<<" eff "<<selected_ggH_mh600/total_ggH_mh600<<endl;
  cout<<"ggH mH = 800 GeV -->  Total events "<<total_ggH_mh800<<" selected "<<selected_ggH_mh800<<" eff "<<selected_ggH_mh800/total_ggH_mh800<<endl;
  cout<<"ggH mH = 1000 GeV -->  Total events "<<total_ggH_mh1000<<" selected "<<selected_ggH_mh1000<<" eff "<<selected_ggH_mh1000/total_ggH_mh1000<<endl;

  ///// -------------------
  vector<float> bins = {90,130,175,250,350,450,550,700,900,1000};

  TH1F* numerator_qqH = new TH1F("numerator_qqH","",bins.size()-1,&bins[0]);
  TH1F* denominator_qqH = new TH1F("denominator_qqH","",bins.size()-1,&bins[0]);
  numerator_qqH->Sumw2();
  denominator_qqH->Sumw2();
  
  numerator_qqH->SetBinContent(1,selected_qqH_mh110);
  numerator_qqH->SetBinContent(2,selected_qqH_mh150);
  numerator_qqH->SetBinContent(3,selected_qqH_mh200);
  numerator_qqH->SetBinContent(4,selected_qqH_mh300);
  numerator_qqH->SetBinContent(5,selected_qqH_mh400);
  numerator_qqH->SetBinContent(6,selected_qqH_mh500);
  numerator_qqH->SetBinContent(7,selected_qqH_mh600);
  numerator_qqH->SetBinContent(8,selected_qqH_mh800);
  numerator_qqH->SetBinContent(9,selected_qqH_mh1000);

  denominator_qqH->SetBinContent(1,total_qqH_mh110);
  denominator_qqH->SetBinContent(2,total_qqH_mh150);
  denominator_qqH->SetBinContent(3,total_qqH_mh200);
  denominator_qqH->SetBinContent(4,total_qqH_mh300);
  denominator_qqH->SetBinContent(5,total_qqH_mh400);
  denominator_qqH->SetBinContent(6,total_qqH_mh500);
  denominator_qqH->SetBinContent(7,total_qqH_mh600);
  denominator_qqH->SetBinContent(8,total_qqH_mh800);
  denominator_qqH->SetBinContent(9,total_qqH_mh1000);
  
  TGraphAsymmErrors* efficiency = new TGraphAsymmErrors();
  efficiency->BayesDivide(numerator_qqH,denominator_qqH);
  efficiency->SetLineColor(kBlack);
  efficiency->SetMarkerColor(kRed);
  efficiency->SetMarkerStyle(20);
  efficiency->SetMarkerSize(1);
  efficiency->GetYaxis()->SetTitle("efficiency");
  efficiency->GetXaxis()->SetTitle("m_{H} [GeV]");
  efficiency->GetYaxis()->SetTitleOffset(1.15);
  efficiency->GetXaxis()->SetTitleOffset(1.1);

  ///// ------------------- 
  TGraphAsymmErrors* efficiency_ggH = NULL;
  if(addGGHEfficiency){

    TH1F* numerator_ggH = new TH1F("numerator_ggH","",bins.size()-1,&bins[0]);
    TH1F* denominator_ggH = new TH1F("denominator_ggH","",bins.size()-1,&bins[0]);
    numerator_ggH->Sumw2();
    denominator_ggH->Sumw2();
    
    numerator_ggH->SetBinContent(1,selected_ggH_mh110);
    numerator_ggH->SetBinContent(2,selected_ggH_mh150);
    numerator_ggH->SetBinContent(3,selected_ggH_mh200);
    numerator_ggH->SetBinContent(4,selected_ggH_mh300);
    numerator_ggH->SetBinContent(5,selected_ggH_mh400);
    numerator_ggH->SetBinContent(6,selected_ggH_mh500);
    numerator_ggH->SetBinContent(7,selected_ggH_mh600);
    numerator_ggH->SetBinContent(8,selected_ggH_mh800);
    numerator_ggH->SetBinContent(9,selected_ggH_mh1000);
    
    denominator_ggH->SetBinContent(1,total_ggH_mh110);
    denominator_ggH->SetBinContent(2,total_ggH_mh150);
    denominator_ggH->SetBinContent(3,total_ggH_mh200);
    denominator_ggH->SetBinContent(4,total_ggH_mh300);
    denominator_ggH->SetBinContent(5,total_ggH_mh400);
    denominator_ggH->SetBinContent(6,total_ggH_mh500);
    denominator_ggH->SetBinContent(7,total_ggH_mh600);
    denominator_ggH->SetBinContent(8,total_ggH_mh800);
    denominator_ggH->SetBinContent(9,total_ggH_mh1000);
    
    efficiency_ggH = new TGraphAsymmErrors();
    efficiency_ggH->BayesDivide(numerator_ggH,denominator_ggH);
    efficiency_ggH->SetLineColor(kBlue);
    efficiency_ggH->SetMarkerColor(kBlue);
    efficiency->SetMarkerColor(kBlack);
    
    efficiency_ggH->SetMarkerStyle(20);
    efficiency_ggH->SetMarkerSize(1);
    efficiency_ggH->GetYaxis()->SetTitle("efficiency");
    efficiency_ggH->GetXaxis()->SetTitle("m_{H} [GeV]");
  }

  
  /// ---------------
  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();
  efficiency->Draw("APE");
  if(addGGHEfficiency)
    efficiency_ggH->Draw("PEsame");
  CMS_lumi(canvas,Form("%.1f",luminosity));

  TLegend leg (0.7,0.3,0.9,0.4);
  if(efficiency_ggH){
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.AddEntry(efficiency,"VBF SM-like","PEL");
    leg.AddEntry(efficiency_ggH,"ggH SM-like","PEL");
    leg.Draw("same");
  }

  if(category == Category::VBFrelaxed)
    efficiency->GetYaxis()->SetRangeUser(0.0001,0.1);
  else
    efficiency->GetYaxis()->SetRangeUser(0.00001,0.1);

  canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/efficiency_vs_mass.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/efficiency_vs_mass.pdf").c_str(),"pdf");

}
