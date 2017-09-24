#include <iostream>
#include <sstream>
#include <cmath>
#include "../triggerUtils.h"
#include "../../CMS_lumi.h"

// recoil binning for monojet                                                                                                                                                                    
vector <float> bins_monojet_recoil    = {150,160,170,180,190,200,210,220,230,240,250,265,280,300,320,340,360,380,400,430,460,490,520,550,580,610,650,700,740,800,900,1000,1250};
vector <float> bins_VBFrelaxed_recoil = {150,160,170,180,190,200,210,220,230,240,250,265,280,300,320,340,360,380,400,430,460,490,520,550,580,610,650,700,740,800,900,1000,1250};
vector <float> bins_VBF_recoil        = {150,175,200,225,230,250,275,300,350,400,450,500,600,700,850,1000};
vector <float> bins_VBFrelaxed_mjj    = {200,400,600,800,1000,1250,1500,1750,2000,2500,3500,5000};
vector <float> bins_VBF_mjj           = {1300,1500,1750,2000,3500,3500,5000};

vector <float> bins_VBFrelaxed_mjj_2d  = {200,500,800,1200,1600,2000,2750};
vector <float> bins_VBFrelaxed_ptj1_2d = {80,110,140,175,225,300,400};
vector <float> bins_VBFrelaxed_ptj2_2d = {40,60,90,120,160,250};
vector <float> bins_VBFrelaxed_detajj_2d = {1,2.5,4,6,10};

// cut for trigger values
static float L1_ETM_CUT  = 50;
static float caloMET_CUT = 70;
static float caloMETClean_CUT = 60;
static float caloMHT_CUT = 70;
static float PFMHT_CUT   = 90;
static float PFMET_CUT   = 90;
static float PFMHTNoMu_CUT = 90;
static float PFMETNoMu_CUT = 90;
static bool  doJetStudy = false;

// possible samples to be selected
enum class Sample {sig,wmn,zmm};
enum class Category {monojet,VBFrelaxed,VBF};


// k-factor file
static string kfactorFile = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_24bins.root";

// calculate the sum of weights for each file from the gentree
void calculateSumWeight(TChain* gentree, vector<double> & wgtsum){

  TTreeReader reader(gentree);
  TTreeReaderValue<float> wgt (reader,"wgt");
  //////////////////                                                                                                                                                                           
  long int nTotal = gentree->GetEntries();
  long int nEvents = 0;
  long int nPart = 100000;
  string currentFile = "";
  int ifile = 0;    
  while(reader.Next()){

    //check if file name switched or not                                                                                                                                                               
    if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      wgtsum.push_back(0);
      ifile++;
      cout<<currentFile<<" "<<ifile<<endl;
    }
    else if(currentFile == ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      ifile = 0;
      wgtsum.push_back(0);
      cout<<currentFile<<" "<<ifile<<endl;
    }
    nEvents++;
    wgtsum.at(ifile) += *wgt;
  }
}

/////
void createHistograms1D(vector<pair<string,TH1F*> > & histo, const string & postfix, const vector<float> & binning, const string & obs){
  
  TH1F* h_inclusive = new TH1F(Form("h_Inclusive_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  TH1F* h_inclusive_cc = new TH1F(Form("h_Inclusive_cc_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  TH1F* h_inclusive_cf = new TH1F(Form("h_Inclusive_cf_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  TH1F* h_inclusive_ff = new TH1F(Form("h_Inclusive_ff_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  TH1F* h_L1ETM     = new TH1F(Form("h_L1ETM_%s_%s",obs.c_str(),postfix.c_str()),"", binning.size()-1, &binning[0]);
  TH1F* h_caloMET   = new TH1F(Form("h_caloMET_%s_%s",obs.c_str(),postfix.c_str()),"", binning.size()-1, &binning[0]);
  TH1F* h_caloMETClean = new TH1F(Form("h_caloMETClean_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  TH1F* h_caloMHT   = new TH1F(Form("h_caloMHT_%s_%s",obs.c_str(),postfix.c_str()),"", binning.size()-1, &binning[0]);
  TH1F* h_PFMHT     = new TH1F(Form("h_PFMHT_%s_%s",obs.c_str(),postfix.c_str()),"", binning.size()-1, &binning[0]);
  TH1F* h_PFMET     = new TH1F(Form("h_PFMET_%s_%s",obs.c_str(),postfix.c_str()),"", binning.size()-1, &binning[0]);
  TH1F* h_PFMHTNoMu = new TH1F(Form("h_PFMHTNoMu_%s_%s",obs.c_str(),postfix.c_str()),"", binning.size()-1, &binning[0]);
  TH1F* h_PFMETNoMu = new TH1F(Form("h_PFMETNoMu_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  TH1F* h_PFMETPath = new TH1F(Form("h_PFMETPath_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  TH1F* h_PFMETNoMuPath = new TH1F(Form("h_PFMETNoMuPath_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  TH1F* h_TotalPath = new TH1F(Form("h_TotalPath_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  TH1F* h_BigTriggerOR = new TH1F(Form("h_BigTriggerOR_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  TH1F* h_BigTriggerOR_cc = new TH1F(Form("h_BigTriggerOR_cc_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  TH1F* h_BigTriggerOR_cf = new TH1F(Form("h_BigTriggerOR_cf_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  TH1F* h_BigTriggerOR_ff = new TH1F(Form("h_BigTriggerOR_ff_%s_%s",obs.c_str(),postfix.c_str()), "", binning.size()-1, &binning[0]);
  
  h_inclusive->Sumw2();
  h_inclusive_cc->Sumw2();
  h_inclusive_cf->Sumw2();
  h_inclusive_ff->Sumw2();
  h_L1ETM->Sumw2();
  h_caloMET->Sumw2();
  h_caloMETClean->Sumw2();
  h_caloMHT->Sumw2();
  h_PFMHT->Sumw2();
  h_PFMET->Sumw2();
  h_PFMHTNoMu->Sumw2();
  h_PFMETNoMu->Sumw2();
  h_PFMETPath->Sumw2();
  h_PFMETNoMuPath->Sumw2();
  h_TotalPath->Sumw2();
  h_BigTriggerOR->Sumw2();
  h_BigTriggerOR_cc->Sumw2();
  h_BigTriggerOR_cf->Sumw2();
  h_BigTriggerOR_ff->Sumw2();

  h_inclusive->SetDirectory(0);
  h_inclusive_cc->SetDirectory(0);
  h_inclusive_cf->SetDirectory(0);
  h_inclusive_ff->SetDirectory(0);
  h_L1ETM->SetDirectory(0);
  h_caloMET->SetDirectory(0);
  h_caloMETClean->SetDirectory(0);
  h_caloMHT->SetDirectory(0);
  h_PFMHT->SetDirectory(0);
  h_PFMET->SetDirectory(0);
  h_PFMHTNoMu->SetDirectory(0);
  h_PFMETNoMu->SetDirectory(0);
  h_PFMETPath->SetDirectory(0);
  h_PFMETNoMuPath->SetDirectory(0);
  h_TotalPath->SetDirectory(0);
  h_BigTriggerOR->SetDirectory(0);
  h_BigTriggerOR_cc->SetDirectory(0);
  h_BigTriggerOR_cf->SetDirectory(0);
  h_BigTriggerOR_ff->SetDirectory(0);
  
  histo.push_back(pair<string,TH1F*>(obs,h_inclusive));
  if(doJetStudy){
    histo.push_back(pair<string,TH1F*>(obs,h_inclusive_cc));
    histo.push_back(pair<string,TH1F*>(obs,h_inclusive_cf));
    histo.push_back(pair<string,TH1F*>(obs,h_inclusive_ff));    
  }
  histo.push_back(pair<string,TH1F*>(obs,h_L1ETM));
  histo.push_back(pair<string,TH1F*>(obs,h_caloMET));
  histo.push_back(pair<string,TH1F*>(obs,h_caloMETClean));
  histo.push_back(pair<string,TH1F*>(obs,h_caloMHT));
  histo.push_back(pair<string,TH1F*>(obs,h_PFMHTNoMu));
  histo.push_back(pair<string,TH1F*>(obs,h_PFMETNoMu));
  histo.push_back(pair<string,TH1F*>(obs,h_PFMHT));
  histo.push_back(pair<string,TH1F*>(obs,h_PFMET));
  histo.push_back(pair<string,TH1F*>(obs,h_PFMETPath));
  histo.push_back(pair<string,TH1F*>(obs,h_PFMETNoMuPath));
  histo.push_back(pair<string,TH1F*>(obs,h_TotalPath));
  histo.push_back(pair<string,TH1F*>(obs,h_BigTriggerOR));
  if(doJetStudy){
    histo.push_back(pair<string,TH1F*>(obs,h_BigTriggerOR_cc));
    histo.push_back(pair<string,TH1F*>(obs,h_BigTriggerOR_cf));
    histo.push_back(pair<string,TH1F*>(obs,h_BigTriggerOR_ff));
  }
}

/////
void createHistograms2D(vector<pair<string,TH2F*> > & histo, const string & postfix, const vector<float> & binningX, const vector<float> & binningY, const string & obs){
  
  TH2F* h_inclusive = new TH2F(Form("h_Inclusive_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_inclusive_cc = new TH2F(Form("h_Inclusive_cc_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_inclusive_cf = new TH2F(Form("h_Inclusive_cf_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_inclusive_ff = new TH2F(Form("h_Inclusive_ff_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_L1ETM     = new TH2F(Form("h_L1ETM_%s_%s",obs.c_str(),postfix.c_str()),"", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_caloMET   = new TH2F(Form("h_caloMET_%s_%s",obs.c_str(),postfix.c_str()),"", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_caloMETClean = new TH2F(Form("h_caloMETClean_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_caloMHT   = new TH2F(Form("h_caloMHT_%s_%s",obs.c_str(),postfix.c_str()),"", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_PFMHT     = new TH2F(Form("h_PFMHT_%s_%s",obs.c_str(),postfix.c_str()),"", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_PFMET     = new TH2F(Form("h_PFMET_%s_%s",obs.c_str(),postfix.c_str()),"", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_PFMHTNoMu = new TH2F(Form("h_PFMHTNoMu_%s_%s",obs.c_str(),postfix.c_str()),"", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_PFMETNoMu = new TH2F(Form("h_PFMETNoMu_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_PFMETPath = new TH2F(Form("h_PFMETPath_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_PFMETNoMuPath = new TH2F(Form("h_PFMETNoMuPath_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_TotalPath     = new TH2F(Form("h_TotalPath_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_BigTriggerOR  = new TH2F(Form("h_BigTriggerOR_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_BigTriggerOR_cc = new TH2F(Form("h_BigTriggerOR_cc_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_BigTriggerOR_cf = new TH2F(Form("h_BigTriggerOR_cf_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  TH2F* h_BigTriggerOR_ff = new TH2F(Form("h_BigTriggerOR_ff_%s_%s",obs.c_str(),postfix.c_str()), "", binningX.size()-1, &binningX[0], binningY.size()-1, &binningY[0]);
  
  h_inclusive->Sumw2();
  h_inclusive_cc->Sumw2();
  h_inclusive_cf->Sumw2();
  h_inclusive_ff->Sumw2();
  h_L1ETM->Sumw2();
  h_caloMET->Sumw2();
  h_caloMETClean->Sumw2();
  h_caloMHT->Sumw2();
  h_PFMHT->Sumw2();
  h_PFMET->Sumw2();
  h_PFMHTNoMu->Sumw2();
  h_PFMETNoMu->Sumw2();
  h_PFMETPath->Sumw2();
  h_PFMETNoMuPath->Sumw2();
  h_TotalPath->Sumw2();
  h_BigTriggerOR->Sumw2();
  h_BigTriggerOR_cc->Sumw2();
  h_BigTriggerOR_cf->Sumw2();
  h_BigTriggerOR_ff->Sumw2();
  
  h_inclusive->SetDirectory(0);
  h_inclusive_cc->SetDirectory(0);
  h_inclusive_cf->SetDirectory(0);
  h_inclusive_ff->SetDirectory(0);
  h_L1ETM->SetDirectory(0);
  h_caloMET->SetDirectory(0);
  h_caloMETClean->SetDirectory(0);
  h_caloMHT->SetDirectory(0);
  h_PFMHT->SetDirectory(0);
  h_PFMET->SetDirectory(0);
  h_PFMHTNoMu->SetDirectory(0);
  h_PFMETNoMu->SetDirectory(0);
  h_PFMETPath->SetDirectory(0);
  h_PFMETNoMuPath->SetDirectory(0);
  h_BigTriggerOR->SetDirectory(0);
  h_BigTriggerOR_cc->SetDirectory(0);
  h_BigTriggerOR_cf->SetDirectory(0);
  h_BigTriggerOR_ff->SetDirectory(0);
  
  histo.push_back(pair<string,TH2F*>(obs,h_inclusive));
  if(doJetStudy){
    histo.push_back(pair<string,TH2F*>(obs,h_inclusive_cc));
    histo.push_back(pair<string,TH2F*>(obs,h_inclusive_cf));
    histo.push_back(pair<string,TH2F*>(obs,h_inclusive_ff));    
  }
  histo.push_back(pair<string,TH2F*>(obs,h_L1ETM));
  histo.push_back(pair<string,TH2F*>(obs,h_caloMET));
  histo.push_back(pair<string,TH2F*>(obs,h_caloMETClean));
  histo.push_back(pair<string,TH2F*>(obs,h_caloMHT));
  histo.push_back(pair<string,TH2F*>(obs,h_PFMHTNoMu));
  histo.push_back(pair<string,TH2F*>(obs,h_PFMETNoMu));
  histo.push_back(pair<string,TH2F*>(obs,h_PFMHT));
  histo.push_back(pair<string,TH2F*>(obs,h_PFMET));
  histo.push_back(pair<string,TH2F*>(obs,h_PFMETPath));
  histo.push_back(pair<string,TH2F*>(obs,h_PFMETNoMuPath));
  histo.push_back(pair<string,TH2F*>(obs,h_TotalPath));
  histo.push_back(pair<string,TH2F*>(obs,h_BigTriggerOR));
  if(doJetStudy){
    histo.push_back(pair<string,TH2F*>(obs,h_BigTriggerOR_cc));
    histo.push_back(pair<string,TH2F*>(obs,h_BigTriggerOR_cf));
    histo.push_back(pair<string,TH2F*>(obs,h_BigTriggerOR_ff));
  }
    
}

//////////
void makeTriggerAnalysis(TChain* tree, 
			 vector<pair<string,TH1F*> > &  histograms, // depends on the selection one wants to apply
			 vector<pair<string,TH2F*> > &  histograms2D, // depends on the selection one wants to apply
			 const Sample &   sample, // sample to be selected
			 const Category & category,
			 vector<double> & wgtsum, // sum of weights
			 vector<TH1*> &   khists, // NLO k-factors			 
			 const float &    luminosity, // luminosity
			 const TString &  triggerPath // trigger Path to be studied			 
			 ){
  

  string triggerString;
  string triggerStringNoMu;

  // fix the threshold according to the trigger path
  if(triggerPath.Contains("PFMET90")){
    caloMET_CUT = 70;
    caloMETClean_CUT = 60;
    caloMHT_CUT = 70;
    PFMHT_CUT   = 90;
    PFMET_CUT   = 90;
    PFMHTNoMu_CUT = 90;
    PFMETNoMu_CUT = 90;
    triggerString = "hltmetwithmu90";
    triggerStringNoMu = "hltmet90";
  }
  if(triggerPath.Contains("PFMET100")){
    caloMET_CUT = 80;
    caloMETClean_CUT = 70;
    caloMHT_CUT = 80;
    PFMHT_CUT   = 100;
    PFMET_CUT   = 100;
    PFMHTNoMu_CUT = 100;
    PFMETNoMu_CUT = 100;
    triggerString = "hltmetwithmu100";
    triggerStringNoMu = "hltmet100";
  }
  else if(triggerPath.Contains("PFMET110")){
    caloMET_CUT = 80;
    caloMETClean_CUT = 70;
    caloMHT_CUT = 80;
    PFMHT_CUT   = 110;
    PFMET_CUT   = 110;
    PFMHTNoMu_CUT   = 110;
    PFMETNoMu_CUT   = 110;
    triggerString   = "hltmetwithmu110";
    triggerStringNoMu = "hltmet110";
  }
  else if(triggerPath.Contains("PFMET120")){
    caloMET_CUT = 90;
    caloMETClean_CUT = 80;
    caloMHT_CUT = 90;
    PFMHT_CUT   = 120;
    PFMET_CUT   = 120;
    PFMHTNoMu_CUT   = 120;
    PFMETNoMu_CUT   = 120;
    triggerString   = "hltmetwithmu120";
    triggerStringNoMu = "hltmet120";
  }

  /// pileup weight
  TFile* pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt_36.40_summer16.root");
  TH1* puhist = (TH1*)pufile->Get("puhist");
  puhist->SetDirectory(0);

  /// muon scale factors
  TFile* muonSF_file  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/muon_scalefactors.root");
  TH2* msfloose_id  = (TH2*) muonSF_file->Get("scalefactors_MuonLooseId_Muon");
  TH2* msfloose_iso = (TH2*) muonSF_file->Get("scalefactors_Iso_MuonLooseId");
  TH2* msftight_id  = (TH2*) muonSF_file->Get("scalefactors_TightId_Muon");
  TH2* msftight_iso = (TH2*) muonSF_file->Get("scalefactors_Iso_MuonTightId");
  msfloose_id->SetDirectory(0);
  msfloose_iso->SetDirectory(0);
  msftight_id->SetDirectory(0);
  msftight_iso->SetDirectory(0);

  TTreeReader reader(tree);
  TTreeReaderValue<float> xsec     (reader,"xsec");
  TTreeReaderValue<float> wgt      (reader,"wgt");
  TTreeReaderValue<UChar_t> hltmetwithmu (reader,triggerString.c_str());
  TTreeReaderValue<UChar_t> hltmet (reader,triggerStringNoMu.c_str());
  TTreeReaderValue<UChar_t> hltmet90 (reader,"hltmet90");
  TTreeReaderValue<UChar_t> hltmet100 (reader,"hltmet100");
  TTreeReaderValue<UChar_t> hltmet110 (reader,"hltmet110");
  TTreeReaderValue<UChar_t> hltmet120 (reader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmetwithmu90 (reader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hltmetwithmu100 (reader,"hltmetwithmu100");
  TTreeReaderValue<UChar_t> hltmetwithmu110 (reader,"hltmetwithmu110");
  TTreeReaderValue<UChar_t> hltmetwithmu120 (reader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmetwithmu170 (reader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmetwithmu300 (reader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltjetmet (reader,"hltjetmet");
  TTreeReaderValue<UChar_t> hltsinglemu (reader,"hltsinglemu");
  TTreeReaderValue<UChar_t> hltdoublemu (reader,"hltdoublemu");
  TTreeReaderValue<float>   mu1pt  (reader,"mu1pt");
  TTreeReaderValue<float>   mu1eta (reader,"mu1eta");
  TTreeReaderValue<float>   mu1phi (reader,"mu1phi");
  TTreeReaderValue<int>     mu1id  (reader,"mu1id");
  TTreeReaderValue<int>     mu1pid  (reader,"mu1pid");
  TTreeReaderValue<float>   mu2pt  (reader,"mu2pt");
  TTreeReaderValue<float>   mu2eta (reader,"mu2eta");
  TTreeReaderValue<float>   mu2phi (reader,"mu2phi");
  TTreeReaderValue<int>     mu2id  (reader,"mu2id");
  TTreeReaderValue<int>     mu2pid  (reader,"mu2pid");
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");
  TTreeReaderValue<unsigned int> nvtx        (reader,"nvtx");
  TTreeReaderValue<unsigned int> ntaus       (reader,"ntaus");
  TTreeReaderValue<unsigned int> nmuons      (reader,"nmuons");
  TTreeReaderValue<unsigned int> nelectrons  (reader,"nelectrons");
  TTreeReaderValue<unsigned int> nphotons    (reader,"nphotons");
  TTreeReaderValue<unsigned int> nincjets    (reader,"njetsinc");
  TTreeReaderValue<unsigned int> nbjets      (reader,"nbjetslowpt");
  TTreeReaderValue<vector<float> > jetpt     (reader,"combinejetpt");
  TTreeReaderValue<vector<float> > jeteta    (reader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetphi    (reader,"combinejetphi");
  TTreeReaderValue<vector<float> > jetm      (reader,"combinejetm");
  TTreeReaderValue<vector<float> > jetchfrac (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > jetnhfrac (reader,"combinejetNHfrac");
  TTreeReaderValue<float> mmet        (reader,"t1mumet");
  TTreeReaderValue<float> mmetphi     (reader,"t1mumetphi");
  TTreeReaderValue<float> metpf       (reader,"pfmet");
  TTreeReaderValue<float> metcalo     (reader,"calomet");
  TTreeReaderValue<float> jmmdphi (reader,"incjetmumetdphimin4");
  TTreeReaderValue<float> zmass   (reader,"zmass");
  TTreeReaderValue<float> wzpt    (reader,"wzpt");
  TTreeReaderValue<float> trig_L1ETM_pt (reader,"trig_L1ETM_pt");
  TTreeReaderValue<vector<float> > trig_obj_pt (reader,"trig_obj_pt");
  TTreeReaderValue<vector<string> > trig_obj_col (reader,"trig_obj_col");

  string currentFile = "";
  int ifile = 0;    
  
  
  //////////////////                                                                                                                                                                              
  long int nTotal = tree->GetEntries();
  long int nEvents = 0;    
  long int nPart = 100000;

  while(reader.Next()){/////
    
    //check if file name switched or not                                                                                                                                                               
    if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      wgtsum.push_back(0);
      ifile++;
      cout<<currentFile<<" "<<ifile<<endl;
    }
    else if(currentFile == ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      ifile = 0;
      wgtsum.push_back(0);
      cout<<currentFile<<" "<<ifile<<endl;
    }
    
    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    // to speed-up
    if(category == Category::monojet and *mmet < bins_monojet_recoil.front()) continue;
    else if(category == Category::VBFrelaxed and *mmet < bins_VBFrelaxed_recoil.front()) continue;
    else if(category == Category::VBF and *mmet < bins_VBF_recoil.front()) continue;

    // object veto
    if(*nbjets   != 0)    continue;
    if(*ntaus    != 0)    continue;
    if(*nphotons != 0)    continue;
    if(*nelectrons != 0 ) continue;
      
    // muons
    if(sample == Sample::sig and *nmuons != 0) continue;
    else if(sample == Sample::wmn and *nmuons !=1) continue;
    else if(sample == Sample::zmm and *nmuons !=2) continue;
    
    // met filters
    if(not *fcsc)  continue;
    if(not *fcsct) continue;
    if(not *feeb)  continue;
    if(not *fetp)  continue;
    if(not *fvtx)  continue;
    if(not *fbadmu) continue;
    if(not *fbadch) continue;
    if(not *fhbhe)  continue;
    if(not *fhbiso) continue;
    
    // apply calo met cleaning
    if(fabs(*metpf-*metcalo)/(*mmet) > 0.5) continue;
    
    // ask single muon trigger in case of selecting wmn events
    if(sample == Sample::wmn and not *hltsinglemu) continue;
    else if(sample == Sample::zmm and not *hltdoublemu) continue;
    
    // apply standard reco-level cuts for the control region
    if(sample == Sample::wmn){
      if(*mu1pt < 20) continue;
      if(fabs(*mu1eta) > 2.4) continue;
      if(*mu1id  !=1) continue;
      float dphi = fabs(*mu1phi-*mmetphi);
      if(dphi > TMath::Pi())
	dphi = 2*TMath::Pi()-dphi;
      float mtw = sqrt(2*(*mu1pt)*(*mmet)*(1-cos(dphi)));
      if(mtw > 160) continue;
    }
    else if(sample == Sample::zmm){
      if(*mu1pt < 20) continue;
      if(fabs(*mu1eta) > 2.4) continue;
      if(fabs(*mu2eta) > 2.4) continue;
      bool isGoodMuon = false;
      if(*mu1pt > 20 and *mu1id  == 1) isGoodMuon = true;
      if(*mu2pt > 20 and *mu2id  == 1) isGoodMuon = true;
      if(not isGoodMuon) continue;
      if(*mu1pid == *mu2pid) continue;
      if(*zmass > 120 or *zmass < 60) continue;
    }
    
    // apply jet-met dphi--> not for gen level analysis
    if(*jmmdphi < 0.5) continue;
    
    TLorentzVector jet1,jet2;
    
    //// monojet category
    if(category == Category::monojet){
      // apply jet pt selections
      if(*nincjets < 1) continue;
      if(jetpt->size() == 0) continue;
      if(jetpt->at(0) < 100) continue;
      if(fabs(jeteta->at(0)) > 2.5) continue;
      if(jetchfrac->at(0) < 0.1) continue;
      if(jetnhfrac->at(0) > 0.8) continue;
      
    }
    else if(category == Category::VBFrelaxed){
      if(*nincjets < 1) continue;
      if(jetpt->size() <= 1) continue;
      if(jetpt->at(0) < 80) continue;
      if(jetpt->at(1) < 40) continue;
      if(fabs(jeteta->at(0)) > 4.7 or fabs(jeteta->at(1)) > 4.7) continue;
      if(fabs(jeteta->at(0)) < 2.4 and jetchfrac->at(0) < 0.1) continue;
      if(fabs(jeteta->at(0)) < 2.4 and jetnhfrac->at(0) > 0.8) continue;
      if(jeteta->at(0)*jeteta->at(1) > 0) continue;
      if(fabs(jeteta->at(0)-jeteta->at(1)) < 1) continue;
      
      jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
      jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
      
      if((jet1+jet2).M() < 200) continue;
      if(fabs(jet1.DeltaPhi(jet2)) > 1.5) continue;
      if(not doJetStudy and fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) > 3) continue;
    }
    else if(category == Category::VBF){
      if(*nincjets < 1) continue;
      if(jetpt->size() <= 1) continue;
      if(jetpt->at(0) < 80) continue;
      if(jetpt->at(1) < 40) continue;
      if(fabs(jeteta->at(0)) > 4.7 or fabs(jeteta->at(1)) > 4.7) continue;
      if(fabs(jeteta->at(0)) < 2.4 and jetchfrac->at(0) < 0.1) continue;
      if(fabs(jeteta->at(0)) < 2.4 and jetnhfrac->at(0) > 0.8) continue;
      if(jeteta->at(0)*jeteta->at(1) > 0) continue;
      if(fabs(jeteta->at(0)-jeteta->at(1)) < 3) continue;
      
      jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
      jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
      
      if((jet1+jet2).M() < 1300) continue;
      if(fabs(jet1.DeltaPhi(jet2)) > 1.5) continue;
      if(not doJetStudy and fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) > 3) continue;
      }
    
    // pileup re-weight
    double puwgt = 1;
    if(*nvtx < 60)
      puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
    
    // apply muon scale factors
    double sfwgt = 1.;
    if(*mu1pt > 0. and (sample == Sample::wmn or sample == Sample::zmm)){
      float ptValue = *mu1pt;
      if(ptValue < msftight_id->GetYaxis()->GetBinLowEdge(1)) 
	ptValue =  msftight_id->GetYaxis()->GetBinLowEdge(1)+1;
      else if(ptValue > msftight_id->GetYaxis()->GetBinLowEdge(msftight_id->GetNbinsY()+1)) 
	ptValue = msftight_id->GetYaxis()->GetBinLowEdge(msftight_id->GetNbinsY()+1)-1;
	
      if(*mu1id == 1)
	sfwgt *= msftight_id->GetBinContent(msftight_id->FindBin(fabs(*mu1eta),ptValue))*msftight_iso->GetBinContent(msftight_iso->FindBin(fabs(*mu1eta),ptValue));
      else
	sfwgt *= msfloose_id->GetBinContent(msfloose_id->FindBin(fabs(*mu1eta),ptValue))*msfloose_iso->GetBinContent(msfloose_iso->FindBin(fabs(*mu1eta),ptValue));
    }
    if(*mu2pt > 0. and  sample == Sample::zmm){
      float ptValue = *mu2pt;
      if(ptValue < msftight_id->GetYaxis()->GetBinLowEdge(1)) 
	ptValue =  msftight_id->GetYaxis()->GetBinLowEdge(1)+1;
      else if(ptValue > msftight_id->GetYaxis()->GetBinLowEdge(msftight_id->GetNbinsY()+1)) 
	ptValue = msftight_id->GetYaxis()->GetBinLowEdge(msftight_id->GetNbinsY()+1)-1;
      
      if(*mu2id == 1)
	sfwgt *= msftight_id->GetBinContent(msftight_id->FindBin(fabs(*mu2eta),ptValue))*msftight_iso->GetBinContent(msftight_iso->FindBin(fabs(*mu2eta),ptValue));
      else
	sfwgt *= msfloose_id->GetBinContent(msfloose_id->FindBin(fabs(*mu2eta),ptValue))*msfloose_iso->GetBinContent(msfloose_iso->FindBin(fabs(*mu2eta),ptValue));
    }
      
    
    // calculate NLO k-factors
    Double_t kwgt = 1.0;
    double   genpt = *wzpt;
    for (size_t i = 0; i < khists.size(); i++) {
      if (khists[i]) {
	if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
	kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));	  
      }
    }
    
    // found numerical values
    float caloMET = 0;
    float caloMETClean = 0;
    float caloMHT = 0;
    float PFMHT = 0;
    float PFMET = 0;
    float PFMHTNoMu = 0;
    float PFMETNoMu = 0;
      
    for(size_t itrig = 0; itrig < trig_obj_col->size(); itrig++){
      TString name (trig_obj_col->at(itrig));
      if(name.Contains("hltMetClean")) caloMETClean = trig_obj_pt->at(itrig);
      else if(name.Contains("hltMet")) caloMET = trig_obj_pt->at(itrig);
      else if(name.Contains("hltMht")) caloMHT = trig_obj_pt->at(itrig);
      else if(name.Contains("hltPFMHTNoMuTightID")) PFMHTNoMu = trig_obj_pt->at(itrig);
      else if(name.Contains("hltPFMETNoMuProducer")) PFMETNoMu = trig_obj_pt->at(itrig);
      else if(name.Contains("hltPFMHTTightID")) PFMHT = trig_obj_pt->at(itrig);
      else if(name.Contains("hltPFMETProducer")) PFMET = trig_obj_pt->at(itrig);
    }
    
    // to include the overflow
    for(auto obs : histograms){
      
      if(obs.first != "met" and *mmet < 250) continue; // cut harder in met
      
      // values
      float value = 0;
      if(obs.first == "met")
	value = *mmet;
      else if(obs.first == "mjj")
	value = (jet1+jet2).M();
	
      // handling overflow
      if(value > obs.second->GetXaxis()->GetBinLowEdge(obs.second->GetNbinsX()+1))
	value = obs.second->GetXaxis()->GetBinLowEdge(obs.second->GetNbinsX()+1)-1;	
      
      
      TString name (obs.second->GetName());
      if(name.Contains("Inclusive_cc")){
	if(fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) < 3)	    
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      else if(name.Contains("Inclusive_cf")){
	if(((fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) > 3) or (fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) < 3)))
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      else if(name.Contains("Inclusive_ff")){ // fill with all the events i.e. numerator
	if(fabs(jet1.Eta()) > 3.0 and fabs(jet2.Eta()) > 3.)
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      else if(name.Contains("Inclusive")) // fill with all the events i.e. numerator
	obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      
      // L1 efficiency
      else if(name.Contains("L1ETM")){ // fill with all the events passing L1 selection
	if(*trig_L1ETM_pt > L1_ETM_CUT) 
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      
      // L1 + caloMET + caloMETClean
      else if(name.Contains("caloMETClean")){// fill with all the events passing calo-MET + calo-MET clean cut
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT) 
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      
      // L1 + caloMET
      else if(name.Contains("caloMET")){ // fill with all the events passing calo-MET 
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT)
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      
	// L1 + caloMET + caloMETClean + caloMHT
      else if(name.Contains("caloMHT")){// fill with all the events passing calo MHT
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT and
	   caloMHT > caloMHT_CUT)
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      
      // L1 + caloMET + caloMETClean + caloMHT + PFMHTNoMu
      else if(name.Contains("PFMHTNoMu")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and 
	   caloMETClean > caloMETClean_CUT and
	   caloMHT > caloMHT_CUT and
	   PFMHTNoMu > PFMHTNoMu_CUT)
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      
      // L1 + caloMET + caloMETClean + caloMHT + PFMHT
      else if(name.Contains("PFMHT")){  // fill with all the events passing PFMHT
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT and
	   caloMHT > caloMHT_CUT and
	   PFMHT > PFMHT_CUT)
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      
      // L1 + caloMET + caloMETClean + caloMHT + PFMHTNoMu + PFMETNoMu
      else if(name.Contains("PFMETNoMu")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT and
	   caloMHT   > caloMHT_CUT and
	   PFMHTNoMu > PFMHTNoMu_CUT and
	   PFMETNoMu > PFMETNoMu_CUT)
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	    
      }	  
	
      // L1 + caloMET + caloMETClean + caloMHT + PFMHT + PFMET
      else if(name.Contains("PFMET")){ // fill with all the events passing PFMET
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT and
	   caloMHT > caloMHT_CUT and
	   PFMHT > PFMHT_CUT and
	   PFMET > PFMET_CUT) 
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
      
	// PFMET path
      else if(name.Contains("PFMETPath")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   *hltmetwithmu != 0) 
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
      
      // PFMETNoMu path
      else if(name.Contains("PFMETNoMuPath")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   *hltmet != 0) 
	    obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
      
      // PFMETNoMu path
      else if(name.Contains("TotalPath")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   (*hltmet != 0 or *hltmetwithmu != 0)) 
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
      
      else if(name.Contains("h_BigTriggerOR_cc")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   (*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	    *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	    *hltjetmet != 0) and fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) < 3)
	    obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
      
      else if(name.Contains("h_BigTriggerOR_cf")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   (*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	    *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	    *hltjetmet != 0) and ((fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) > 3) or (fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) < 3)))
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
      
      else if(name.Contains("h_BigTriggerOR_ff")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   (*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	    *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	    *hltjetmet != 0) and fabs(jet1.Eta()) > 3.0 and fabs(jet2.Eta()) > 3.)
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
      
      else if(name.Contains("h_BigTriggerOR")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   (*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	    *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	    *hltjetmet != 0))
	  obs.second->Fill(value,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
    }            
    
    // Repeat the same for 2D histograms
    for(auto obs : histograms2D){
      
      if(*mmet < 250) continue; // cut harder in met
      
      // values
      float valueX = 0;
      float valueY = 0;
      if(obs.first == "mjj-ptj1"){
	valueX = (jet1+jet2).M();
	valueY = jet1.Pt();
      }
      else if(obs.first == "mjj-ptj2"){
	valueX = (jet1+jet2).M();
	valueY = jet2.Pt();
      }
      else if(obs.first == "mjj-detajj"){
	valueX = (jet1+jet2).M();
	valueY = fabs(jet1.Eta()-jet2.Eta());
      }
            
      // handling overflow
      if(valueX > obs.second->GetXaxis()->GetBinLowEdge(obs.second->GetNbinsX()+1))
	valueX = obs.second->GetXaxis()->GetBinLowEdge(obs.second->GetNbinsX()+1)-1;	
      // handling overflow
      if(valueY > obs.second->GetYaxis()->GetBinLowEdge(obs.second->GetNbinsY()+1))
	valueY = obs.second->GetYaxis()->GetBinLowEdge(obs.second->GetNbinsY()+1)-1;	

      TString name (obs.second->GetName());
      if(name.Contains("Inclusive_cc")){
	if(fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) < 3)	    
	    obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      else if(name.Contains("Inclusive_cf")){
	if(((fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) > 3) or (fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) < 3)))
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      else if(name.Contains("Inclusive_ff")){ // fill with all the events i.e. numerator
	if(fabs(jet1.Eta()) > 3.0 and fabs(jet2.Eta()) > 3.)
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      else if(name.Contains("Inclusive")) // fill with all the events i.e. numerator
	obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      
      // L1 efficiency
      else if(name.Contains("L1ETM")){ // fill with all the events passing L1 selection
	if(*trig_L1ETM_pt > L1_ETM_CUT) 
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      
      // L1 + caloMET + caloMETClean
      else if(name.Contains("caloMETClean")){// fill with all the events passing calo-MET + calo-MET clean cut
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT) 
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      
      // L1 + caloMET
      else if(name.Contains("caloMET")){ // fill with all the events passing calo-MET 
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT)
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      
      // L1 + caloMET + caloMETClean + caloMHT
      else if(name.Contains("caloMHT")){// fill with all the events passing calo MHT
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT and
	   caloMHT > caloMHT_CUT)
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      
      // L1 + caloMET + caloMETClean + caloMHT + PFMHTNoMu
      else if(name.Contains("PFMHTNoMu")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and 
	   caloMETClean > caloMETClean_CUT and
	   caloMHT > caloMHT_CUT and
	   PFMHTNoMu > PFMHTNoMu_CUT)
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      
      // L1 + caloMET + caloMETClean + caloMHT + PFMHT
      else if(name.Contains("PFMHT")){  // fill with all the events passing PFMHT
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT and
	   caloMHT > caloMHT_CUT and
	   PFMHT > PFMHT_CUT)
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));
      }
      
      // L1 + caloMET + caloMETClean + caloMHT + PFMHTNoMu + PFMETNoMu
      else if(name.Contains("PFMETNoMu")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT and
	   caloMHT   > caloMHT_CUT and
	   PFMHTNoMu > PFMHTNoMu_CUT and
	   PFMETNoMu > PFMETNoMu_CUT)
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	    
      }	  
      
      // L1 + caloMET + caloMETClean + caloMHT + PFMHT + PFMET
      else if(name.Contains("PFMET")){ // fill with all the events passing PFMET
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT and
	   caloMHT > caloMHT_CUT and
	   PFMHT > PFMHT_CUT and
	   PFMET > PFMET_CUT) 
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
      
      // PFMET path
      else if(name.Contains("PFMETPath")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   *hltmetwithmu != 0) 
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
      
      // PFMETNoMu path
      else if(name.Contains("PFMETNoMuPath")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   *hltmet != 0) 
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
      
      else if(name.Contains("TotalPath")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   (*hltmet != 0 or *hltmetwithmu != 0)) 
	   obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      }
      
      /////// ---
      else if(name.Contains("h_BigTriggerOR_cc")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   (*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	    *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	    *hltjetmet != 0) and fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) < 3)
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 

      else if(name.Contains("h_BigTriggerOR_cf")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   (*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	    *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	    *hltjetmet != 0) and ((fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) > 3) or (fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) < 3)))
	    obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
      
      else if(name.Contains("h_BigTriggerOR_ff")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   (*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	    *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	    *hltjetmet != 0) and fabs(jet1.Eta()) > 3.0 and fabs(jet2.Eta()) > 3.)
	    obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
      
      else if(name.Contains("h_BigTriggerOR")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   (*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	    *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	    *hltjetmet != 0))
	  obs.second->Fill(valueX,valueY,luminosity*(*xsec)*(*wgt)*puwgt*sfwgt*kwgt/wgtsum.at(ifile));	  
      } 
    }            
  }
  cout<<endl;
  if(pufile) pufile->Close();
  if(muonSF_file) muonSF_file->Close();
}


//// -----
void fillTreeList(TChain* tree, TChain* gentree, const string & inputDIR, const string & postfix){

  cout<<"Read list of files for "<<postfix<<" process "<<endl;
  system(("ls "+inputDIR+" | grep "+postfix+" > list_dir.txt").c_str());
  ifstream file_dir("list_dir.txt");
  if(file_dir.is_open()){
    string line;
    while(!file_dir.eof()){
      getline(file_dir,line);
      if(line == "") continue;
      system(("find "+inputDIR+"/"+line+" -name  \"*.root\" > list.txt").c_str());
      ifstream file_("list.txt");
      if(file_.is_open()){
	string line2;
	while(!file_.eof()){
	  getline(file_,line2);
	  if(TString(line2).Contains("failed")) continue;	  
	  if(line == "" or not TString(line2).Contains("root")) continue;
	  cout<<"Open "<<postfix<<" file with name: "<<line2<<endl;
	  tree->Add(line2.c_str());
	  gentree->Add(line2.c_str());
	}
      }
      system("rm list.txt");
    }
  }
  system("rm list_dir.txt");
}

/// plotting efficiency
void plotEfficiency(TCanvas* canvas, 
		    TH1* histo_zjet, 
		    TH1* histo_wjet, 
		    TH1* histo_wmn, 
		    TH1* histo_zmm, 
		    const string & outputDIR, 
		    const float &  luminosity, 
		    const string & observable){

  
  histo_zjet->SetMarkerColor(kBlack);
  histo_zjet->SetLineColor(kBlack);
  histo_zjet->SetMarkerStyle(20);
  histo_zjet->SetMarkerSize(0.75);
  histo_zjet->SetLineWidth(1);
  histo_zjet->SetFillStyle(0);
  histo_zjet->SetFillColor(0);

  histo_wjet->SetMarkerColor(kRed);
  histo_wjet->SetLineColor(kRed);
  histo_wjet->SetMarkerStyle(24);
  histo_wjet->SetMarkerSize(0.75);
  histo_wjet->SetLineWidth(1);
  histo_wjet->SetFillStyle(0);
  histo_wjet->SetFillColor(0);

  histo_wmn->SetMarkerColor(kBlue);
  histo_wmn->SetLineColor(kBlue);
  histo_wmn->SetMarkerStyle(20);
  histo_wmn->SetMarkerStyle(26);
  histo_wmn->SetMarkerSize(0.75);
  histo_wmn->SetLineWidth(1);
  histo_wmn->SetFillStyle(0);
  histo_wmn->SetFillColor(0);

  histo_zmm->SetMarkerColor(kGreen+1);
  histo_zmm->SetLineColor(kGreen+1);
  histo_zmm->SetMarkerStyle(20);
  histo_zmm->SetMarkerStyle(26);
  histo_zmm->SetMarkerSize(0.75);
  histo_zmm->SetLineWidth(1);
  histo_zmm->SetFillStyle(0);
  histo_zmm->SetFillColor(0);

  // Plotting final result for MC turn ons
  TH1* frame = NULL;
  if(TString(observable).Contains("mjj"))
    frame = canvas->DrawFrame(histo_zjet->GetXaxis()->GetBinLowEdge(1),0.2,
			      histo_zjet->GetXaxis()->GetBinLowEdge(histo_zjet->GetNbinsX()+1), 1.1, "");
  else
    frame = canvas->DrawFrame(histo_zjet->GetXaxis()->GetBinLowEdge(1),0.6,
			      histo_zjet->GetXaxis()->GetBinLowEdge(histo_zjet->GetNbinsX()+1), 1.1, "");
  
  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(1.30);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetLabelSize(0.03);
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetLabelSize(0.85*frame->GetYaxis()->GetLabelSize());

  canvas->SetTopMargin(0.06);
  canvas->SetBottomMargin(0.3);
  canvas->Draw();
  canvas->cd();
  frame->Draw();
  CMS_lumi(canvas,string(Form("%.2f",luminosity)),true);

  histo_zjet->Draw("EPLsame");
  histo_wjet->Draw("EPLsame");
  histo_wmn->Draw("EPLsame");
  histo_zmm->Draw("EPLsame");

  TLegend* leg = NULL;
  if(TString(observable).Contains("met"))
    leg = new TLegend(0.6,0.4,0.9,0.6);
  else
    leg = new TLegend(0.25,0.4,0.55,0.6);

  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histo_zjet,"Z(#nu#nu)+jets SR","EPL");
  leg->AddEntry(histo_wjet,"W+jets SR","EPL");
  leg->AddEntry(histo_wmn,"W+jets W(#mu#nu)-CR","EPL");
  leg->AddEntry(histo_zmm,"Z+jets Z(#mu#mu)-CR","EPL");
  leg->Draw("same");

  TString plotName(histo_zjet->GetName());
  plotName.ReplaceAll("_zjet","");

  //// ----
  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);
  pad2->Draw();
  pad2->cd();

  // Plotting final result for MC turn ons
  TH1* frame2 = NULL;
  if(TString(observable).Contains("mjj"))
    frame2 = pad2->DrawFrame(histo_zjet->GetXaxis()->GetBinLowEdge(1),0.8,
			     histo_zjet->GetXaxis()->GetBinLowEdge(histo_zjet->GetNbinsX()+1), 1.05, "");
  else
    frame2 = pad2->DrawFrame(histo_zjet->GetXaxis()->GetBinLowEdge(1),0.9,
			     histo_zjet->GetXaxis()->GetBinLowEdge(histo_zjet->GetNbinsX()+1), 1.05, "");

  if(observable == "met")
    frame2->GetXaxis()->SetTitle("Recoil [GeV]");
  else if(observable == "mjj")
    frame2->GetXaxis()->SetTitle("M_{jj} [GeV]");
  else if(observable == "ptj1")
    frame2->GetXaxis()->SetTitle("p_{T}_{j1} [GeV]");
  else if(observable == "ptj2")
    frame2->GetXaxis()->SetTitle("p_{T}_{j2} [GeV]");
  else if(observable == "detajj")
    frame2->GetXaxis()->SetTitle("#Delta#eta_{jj}");

  frame2->GetYaxis()->SetTitle("Ratio");
  frame2->GetXaxis()->SetTitleOffset(1.1);
  frame2->GetYaxis()->SetTitleOffset(1.32);
  frame2->GetYaxis()->SetTitleSize(0.04);
  frame2->GetYaxis()->SetLabelSize(0.03);
  frame2->GetYaxis()->SetNdivisions(504);
  frame2->GetYaxis()->SetLabelSize(0.85*frame2->GetYaxis()->GetLabelSize());
  frame2->Draw();

  /// ------
  TH1F* ratio_wjet_over_zjet = (TH1F*) histo_wjet->Clone();
  ratio_wjet_over_zjet->Divide(histo_zjet);

  TH1F* ratio_wmn_over_zjet = (TH1F*) histo_wmn->Clone();
  ratio_wmn_over_zjet->Divide(histo_zjet);

  TH1F* ratio_zmm_over_zjet = (TH1F*) histo_zmm->Clone();
  ratio_zmm_over_zjet->Divide(histo_zjet);

  ratio_wjet_over_zjet->SetMarkerColor(kRed);
  ratio_wjet_over_zjet->SetLineColor(kRed);
  ratio_wjet_over_zjet->SetMarkerStyle(24);
  ratio_wjet_over_zjet->SetMarkerSize(0.75);
  ratio_wjet_over_zjet->SetLineWidth(1);
  ratio_wjet_over_zjet->Draw("EPsame");

  ratio_wmn_over_zjet->SetMarkerColor(kBlue);
  ratio_wmn_over_zjet->SetLineColor(kBlue);
  ratio_wmn_over_zjet->SetMarkerStyle(20);
  ratio_wmn_over_zjet->SetMarkerStyle(26);
  ratio_wmn_over_zjet->SetMarkerSize(0.75);
  ratio_wmn_over_zjet->SetLineWidth(1);
  ratio_wmn_over_zjet->Draw("EPsame");

  ratio_zmm_over_zjet->SetMarkerColor(kGreen+1);
  ratio_zmm_over_zjet->SetLineColor(kGreen+1);
  ratio_zmm_over_zjet->SetMarkerStyle(20);
  ratio_zmm_over_zjet->SetMarkerStyle(26);
  ratio_zmm_over_zjet->SetMarkerSize(0.75);
  ratio_zmm_over_zjet->SetLineWidth(1);
  ratio_zmm_over_zjet->Draw("EPsame");

  canvas->SaveAs((outputDIR+"/"+string(plotName.Data())+".pdf").c_str(),"pdf");
  if(frame2) delete frame2;
  if(pad2) delete pad2;
  if(leg) delete leg;
}


/// plotting efficiency
void plotEfficiency(TCanvas* canvas, 
		    TEfficiency* histo_zjet, 
		    TEfficiency* histo_wjet, 
		    TEfficiency* histo_wmn, 
		    TEfficiency* histo_zmm, 
		    const string & outputDIR, 
		    const float &  luminosity, 
		    const string & observable){
  
  // Make efficiencies as TGraph
  TGraphAsymmErrors* graph_zjet = histo_zjet->CreateGraph();
  TGraphAsymmErrors* graph_wjet = histo_wjet->CreateGraph();
  TGraphAsymmErrors* graph_wmn  = histo_wmn->CreateGraph();
  TGraphAsymmErrors* graph_zmm  = histo_zmm->CreateGraph();

  graph_zjet->SetMarkerColor(kBlack);
  graph_zjet->SetLineColor(kBlack);
  graph_zjet->SetMarkerStyle(20);
  graph_zjet->SetMarkerSize(0.75);
  graph_zjet->SetLineWidth(1);

  graph_wjet->SetMarkerColor(kRed);
  graph_wjet->SetLineColor(kRed);
  graph_wjet->SetMarkerStyle(24);
  graph_wjet->SetMarkerSize(0.75);
  graph_wjet->SetLineWidth(1);

  graph_wmn->SetMarkerColor(kBlue);
  graph_wmn->SetLineColor(kBlue);
  graph_wmn->SetMarkerStyle(20);
  graph_wmn->SetMarkerStyle(26);
  graph_wmn->SetMarkerSize(0.75);
  graph_wmn->SetLineWidth(1);

  graph_zmm->SetMarkerColor(kGreen+1);
  graph_zmm->SetLineColor(kGreen+1);
  graph_zmm->SetMarkerStyle(20);
  graph_zmm->SetMarkerStyle(26);
  graph_zmm->SetMarkerSize(0.75);
  graph_zmm->SetLineWidth(1);
    
  // Plotting final result for MC turn ons
  // Plotting final result for MC turn ons
  TH1* frame = NULL;
  if(TString(observable).Contains("mjj"))
    frame = canvas->DrawFrame(histo_zjet->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(1),0.2,
			      histo_zjet->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(histo_zjet->GetPassedHistogram()->GetNbinsX()+1), 1.1, "");
  else
    frame = canvas->DrawFrame(histo_zjet->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(1),0.6,
			      histo_zjet->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(histo_zjet->GetPassedHistogram()->GetNbinsX()+1), 1.1, "");
  

  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(1.30);
  frame->GetYaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetLabelSize(0.03);
  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetLabelSize(0.85*frame->GetYaxis()->GetLabelSize());

  canvas->SetTopMargin(0.06);
  canvas->SetBottomMargin(0.3);
  canvas->Draw();
  canvas->cd();
  frame->Draw();

  CMS_lumi(canvas,string(Form("%.2f",luminosity)),true);

  graph_zjet->Draw("EPLsame");
  graph_wjet->Draw("EPLsame");
  graph_wmn->Draw("EPLsame");
  graph_zmm->Draw("EPLsame");
  
  TLegend* leg = NULL;
  if(TString(observable).Contains("met"))
    leg = new TLegend(0.6,0.4,0.9,0.6);
  else
    leg = new TLegend(0.25,0.4,0.55,0.6);

  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(graph_zjet,"Z(#nu#nu)+jets SR","EPL");
  leg->AddEntry(graph_wjet,"W+jets SR","EPL");
  leg->AddEntry(graph_wmn,"W+jets W(#mu#nu)-CR","EPL");
  leg->AddEntry(graph_zmm,"Z+jets Z(#mu#mu)-CR","EPL");
  leg->Draw("same");

  TString plotName(histo_zjet->GetName());
  plotName.ReplaceAll("_zjet","");

  //// ----
  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);
  pad2->Draw();
  pad2->cd();

  // Plotting final result for MC turn ons                                                                                                                                                             
  TH1* frame2 = NULL;
  if(TString(observable).Contains("mjj"))
    frame2 = pad2->DrawFrame(histo_zjet->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(1),0.8,
                             histo_zjet->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(histo_zjet->GetPassedHistogram()->GetNbinsX()+1), 1.05, "");
  else
    frame2 = pad2->DrawFrame(histo_zjet->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(1),0.9,
                             histo_zjet->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(histo_zjet->GetPassedHistogram()->GetNbinsX()+1), 1.05, "");
  
  
  if(observable == "met")
    frame2->GetXaxis()->SetTitle("Recoil [GeV]");
  else if(observable == "mjj")
    frame2->GetXaxis()->SetTitle("M_{jj} [GeV]");
  else if(observable == "ptj1")
    frame2->GetXaxis()->SetTitle("p_{T}_{j1} [GeV]");
  else if(observable == "ptj2")
    frame2->GetXaxis()->SetTitle("p_{T}_{j2} [GeV]");
  else if(observable == "detajj")
    frame2->GetXaxis()->SetTitle("#Delta#eta_{jj}");

  frame2->GetYaxis()->SetTitle("Ratio");
  frame2->GetXaxis()->SetTitleOffset(1.1);
  frame2->GetYaxis()->SetTitleOffset(1.32);
  frame2->GetYaxis()->SetTitleSize(0.04);
  frame2->GetYaxis()->SetLabelSize(0.03);
  frame2->GetYaxis()->SetNdivisions(504);
  frame2->GetYaxis()->SetLabelSize(0.85*frame2->GetYaxis()->GetLabelSize());
  frame2->Draw();

  /// ------
  TH1F* ratio_wjet_over_zjet = (TH1F*) histo_wjet->GetPassedHistogram();
  ratio_wjet_over_zjet->Divide(histo_wjet->GetTotalHistogram());
  TH1F* ratio_temp = (TH1F*) histo_zjet->GetPassedHistogram();
  ratio_temp->Divide(histo_zjet->GetTotalHistogram());
  ratio_wjet_over_zjet->Divide(ratio_temp);

  /// ------
  TH1F* ratio_wmn_over_zjet = (TH1F*) histo_wmn->GetPassedHistogram();
  ratio_wmn_over_zjet->Divide(histo_wmn->GetTotalHistogram());
  ratio_wmn_over_zjet->Divide(ratio_temp);

  /// ------
  TH1F* ratio_zmm_over_zjet = (TH1F*) histo_zmm->GetPassedHistogram();
  ratio_zmm_over_zjet->Divide(histo_zmm->GetTotalHistogram());
  ratio_zmm_over_zjet->Divide(ratio_temp);
  
  ratio_wjet_over_zjet->SetMarkerColor(kRed);
  ratio_wjet_over_zjet->SetLineColor(kRed);
  ratio_wjet_over_zjet->SetMarkerStyle(24);
  ratio_wjet_over_zjet->SetMarkerSize(0.75);
  ratio_wjet_over_zjet->SetLineWidth(1);
  ratio_wjet_over_zjet->Draw("EPsame");

  ratio_wmn_over_zjet->SetMarkerColor(kBlue);
  ratio_wmn_over_zjet->SetLineColor(kBlue);
  ratio_wmn_over_zjet->SetMarkerStyle(20);
  ratio_wmn_over_zjet->SetMarkerStyle(26);
  ratio_wmn_over_zjet->SetMarkerSize(0.75);
  ratio_wmn_over_zjet->SetLineWidth(1);
  ratio_wmn_over_zjet->Draw("EPsame");

  ratio_zmm_over_zjet->SetMarkerColor(kGreen+1);
  ratio_zmm_over_zjet->SetLineColor(kGreen+1);
  ratio_zmm_over_zjet->SetMarkerStyle(20);
  ratio_zmm_over_zjet->SetMarkerStyle(26);
  ratio_zmm_over_zjet->SetMarkerSize(0.75);
  ratio_zmm_over_zjet->SetLineWidth(1);
  ratio_zmm_over_zjet->Draw("EPsame");

  canvas->SaveAs((outputDIR+"/"+string(plotName.Data())+".pdf").c_str(),"pdf");

  if(frame2) delete frame2;
  if(pad2) delete pad2;
  if(leg) delete leg;
}


/// Main function 
void makeMETTriggerEfficiencyMC_ZvsW(string   inputDIR, 
				     string   outputDIR, 
				     Category category,
				     string   triggerPath, // PFMET90, PFMET100, PFMET110, PFMET120
				     vector<string> observable, // study the trigger efficiency vs mjj or recoil
				     vector<string> observable2D,
				     bool     splitJets = false,
				     float    luminosity = 35.9){

  doJetStudy = splitJets;

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  // check trigger path
  if(triggerPath != "PFMET90"  and
     triggerPath != "PFMET100" and 
     triggerPath != "PFMET110" and 
     triggerPath != "PFMET120"){
    cerr<<"######### Trigger path not known --> exit "<<endl;
    return;
  }
     
  //// input tree                                                                                                                                                                                   
  TChain* tree_zjet = new TChain("tree/tree","tree/tree");
  TChain* gentree_zjet = new TChain("gentree/gentree","gentree/gentree");;
  ////
  TChain* tree_wjet = new TChain("tree/tree","tree/tree");
  TChain* gentree_wjet = new TChain("gentree/gentree","gentree/gentree");;
  ////
  TChain* tree_zmm = new TChain("tree/tree","tree/tree");
  TChain* gentree_zmm = new TChain("gentree/gentree","gentree/gentree");;

  fillTreeList(tree_zjet,gentree_zjet,inputDIR,"ZJets");
  fillTreeList(tree_wjet,gentree_wjet,inputDIR,"WJets");
  fillTreeList(tree_zmm,gentree_zmm,inputDIR,"DYJets");
 
  // sum of weights
  vector<double> wgtsum_zjet;
  vector<double> wgtsum_wjet;
  vector<double> wgtsum_zmm;
  cout<<"######### Calculate Weights for : Z+jets "<<endl;
  calculateSumWeight(gentree_zjet,wgtsum_zjet);
  cout<<"######### Calculate Weights for : W+jets "<<endl;
  calculateSumWeight(gentree_zjet,wgtsum_wjet);
  cout<<"######### Calculate Weights for : DY+jets "<<endl;
  calculateSumWeight(gentree_zmm,wgtsum_zmm);
  

  /////// start analysis  
  vector<pair<string,TH1F*> > histo_zjet;
  vector<pair<string,TH1F*> > histo_wjet;
  vector<pair<string,TH1F*> > histo_wmn;
  vector<pair<string,TH1F*> > histo_zmm;

  for(auto obs : observable){
    vector<float> binning;
    if(category == Category::monojet and obs == "met") binning = bins_monojet_recoil;
    else if(category == Category::VBFrelaxed and obs == "met") binning = bins_VBFrelaxed_recoil;
    else if(category == Category::VBF and obs == "met") binning = bins_VBF_recoil;
  
    if(category == Category::VBFrelaxed and obs == "mjj") binning = bins_VBFrelaxed_mjj;
    else if(category == Category::VBF and obs == "mjj") binning = bins_VBF_mjj;

    if(binning.empty()){
      cerr<<"Observable not found --> skip it "<<endl;
      continue;
    }
    
    createHistograms1D(histo_zjet,"zjet",binning,obs);
    createHistograms1D(histo_wjet,"wjet",binning,obs);
    createHistograms1D(histo_wmn,"wmn",binning,obs);
    createHistograms1D(histo_zmm,"zmm",binning,obs);
  }


  /////// start analysis  
  vector<pair<string,TH2F*> > histo_zjet_2d;
  vector<pair<string,TH2F*> > histo_wjet_2d;
  vector<pair<string,TH2F*> > histo_wmn_2d;
  vector<pair<string,TH2F*> > histo_zmm_2d;

  for(auto obs : observable2D){
    vector<float> binningX;
    vector<float> binningY;
    if(category == Category::VBFrelaxed and obs == "mjj-ptj1"){
      binningX = bins_VBFrelaxed_mjj_2d;
      binningY = bins_VBFrelaxed_ptj1_2d;
    }
    else if(category == Category::VBFrelaxed and obs == "mjj-ptj2"){
      binningX = bins_VBFrelaxed_mjj_2d;
      binningY = bins_VBFrelaxed_ptj2_2d;
    }
    else if(category == Category::VBFrelaxed and obs == "mjj-detajj"){
      binningX = bins_VBFrelaxed_mjj_2d;
      binningY = bins_VBFrelaxed_detajj_2d;
    }

    if(binningX.empty() or binningY.empty()){
      cerr<<"Observable not found --> skip it "<<endl;
      continue;
    }
    
    createHistograms2D(histo_zjet_2d,"zjet",binningX,binningY,obs);
    createHistograms2D(histo_wjet_2d,"wjet",binningX,binningY,obs);
    createHistograms2D(histo_wmn_2d,"wmn",binningX,binningY,obs);
    createHistograms2D(histo_zmm_2d,"zmm",binningX,binningY,obs);
  }

  // k-factors
  TFile* kffile = TFile::Open(kfactorFile.c_str());
  TH1* znlohist = (TH1*) kffile->Get("ZJets_012j_NLO/nominal");
  TH1* zlohist  = (TH1*) kffile->Get("ZJets_LO/inv_pt");
  TH1* zewkhist = (TH1*) kffile->Get("EWKcorr/Z");
    
  if(zewkhist)
    zewkhist->Divide(znlohist);
  if(znlohist)
    znlohist->Divide(zlohist);

  TH1* wnlohist = (TH1*) kffile->Get("WJets_012j_NLO/nominal");
  TH1* wlohist  = (TH1*) kffile->Get("WJets_LO/inv_pt");
  TH1* wewkhist = (TH1*) kffile->Get("EWKcorr/W");
    
  if(wewkhist)
    wewkhist->Divide(wnlohist);
  if(wnlohist)
    wnlohist->Divide(wlohist);
 
  vector<TH1*> khists_zjet;
  vector<TH1*> khists_wjet;

  khists_zjet.push_back(zewkhist);
  khists_zjet.push_back(znlohist);
  khists_wjet.push_back(wewkhist);
  khists_wjet.push_back(wnlohist);

  //// ----- 
  cout<<"######### Loop on Z+jets trees for SR selection "<<endl;
  makeTriggerAnalysis(tree_zjet,histo_zjet,histo_zjet_2d,Sample::sig,category,wgtsum_zjet,khists_zjet,luminosity,triggerPath);
  cout<<"######### Loop on W+jets trees for SR selection "<<endl;
  makeTriggerAnalysis(tree_wjet,histo_wjet,histo_wjet_2d,Sample::sig,category,wgtsum_wjet,khists_wjet,luminosity,triggerPath);
  cout<<"######### Loop on W+jets trees for Wmn selection "<<endl;
  makeTriggerAnalysis(tree_wjet,histo_wmn,histo_wmn_2d,Sample::wmn,category,wgtsum_wjet,khists_wjet,luminosity,triggerPath);
  cout<<"######### Loop on Z+jets trees for Zmm selection "<<endl;
  makeTriggerAnalysis(tree_zmm,histo_zmm,histo_zmm_2d,Sample::zmm,category,wgtsum_zmm,khists_zjet,luminosity,triggerPath);

  ////---
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 650);
  canvas->cd();

  TFile* outputFile = new TFile((outputDIR+"/efficiency.root").c_str(),"RECREATE");
  outputFile->cd();
  
  // calculate the efficiency
  vector<pair<string,TEfficiency* > > eff_zjet;
  vector<pair<string,TEfficiency* > > eff_wjet;
  vector<pair<string,TEfficiency* > > eff_wmn;
  vector<pair<string,TEfficiency* > > eff_zmm;
  
  size_t idenom = 0;
  for(size_t ihist = 1; ihist < histo_zjet.size(); ihist++){
    
    if(histo_zjet.at(ihist).first != histo_zjet.at(idenom).first){
      idenom = ihist;
      continue;
    }
    
    TString name (histo_zjet.at(ihist).second->GetName());
    size_t idenom_temp = idenom;
    TString name2;      
    if(name.Contains("_cc") and not name.Contains("Inclusive")){
      for(size_t jhist = 1; jhist < histo_zjet.size(); jhist++){
	name2 = TString(histo_zjet.at(jhist).second->GetName());
	if(name2.Contains("_cc") and name2.Contains("Inclusive")){
	  idenom_temp = jhist;
	  break;
	}
      }
    }
    else if(name.Contains("_cf") and not name.Contains("Inclusive")){
      for(size_t jhist = 1; jhist < histo_zjet.size(); jhist++){
	name2 = TString(histo_zjet.at(jhist).second->GetName());
	if(name2.Contains("_cf") and name2.Contains("Inclusive")){
	  idenom_temp = jhist;
	  break;
	}
      }
    }
    else if(name.Contains("_ff") and not name.Contains("Inclusive")){
      for(size_t jhist = 1; jhist < histo_zjet.size(); jhist++){
	name2 = TString(histo_zjet.at(jhist).second->GetName());
	if(name2.Contains("_ff") and name2.Contains("Inclusive")){
	  idenom_temp = jhist;
	  break;
	}
      }
    }

    name.ReplaceAll("h_","eff_");
    eff_zjet.push_back(pair<string,TEfficiency*>(histo_zjet.at(ihist).first,new TEfficiency(*histo_zjet.at(ihist).second,*histo_zjet.at(idenom_temp).second)));
    eff_zjet.back().second->SetName(name.Data());

    name = TString(histo_wjet.at(ihist).second->GetName());
    name.ReplaceAll("h_","eff_");
    eff_wjet.push_back(pair<string,TEfficiency*>(histo_wjet.at(ihist).first,new TEfficiency(*histo_wjet.at(ihist).second,*histo_wjet.at(idenom_temp).second)));
    eff_wjet.back().second->SetName(name.Data());
    
    name = TString(histo_wmn.at(ihist).second->GetName());
    name.ReplaceAll("h_","eff_");
    eff_wmn.push_back(pair<string,TEfficiency*>(histo_wmn.at(ihist).first,new TEfficiency(*histo_wmn.at(ihist).second,*histo_wmn.at(idenom_temp).second)));
    eff_wmn.back().second->SetName(name.Data());

    name = TString(histo_zmm.at(ihist).second->GetName());
    name.ReplaceAll("h_","eff_");
    eff_zmm.push_back(pair<string,TEfficiency*>(histo_zmm.at(ihist).first,new TEfficiency(*histo_zmm.at(ihist).second,*histo_zmm.at(idenom_temp).second)));
    eff_zmm.back().second->SetName(name.Data());

    eff_zjet.back().second->Write();
    eff_wjet.back().second->Write();
    eff_wmn.back().second->Write();
    eff_zmm.back().second->Write();

  }

  // plot Efficiency
  for(size_t ihist = 0; ihist < eff_zjet.size(); ihist++)
    plotEfficiency(canvas,eff_zjet.at(ihist).second,eff_wjet.at(ihist).second,eff_wmn.at(ihist).second,eff_zmm.at(ihist).second,outputDIR,luminosity,eff_zjet.at(ihist).first);

  // 2D efficiencies
  vector<pair<string,TEfficiency* > > eff_zjet_2d;
  vector<pair<string,TEfficiency* > > eff_wjet_2d;
  vector<pair<string,TEfficiency* > > eff_wmn_2d;
  vector<pair<string,TEfficiency* > > eff_zmm_2d;
  
  idenom = 0;
  for(size_t ihist = 1; ihist < histo_zjet_2d.size(); ihist++){

    if(histo_zjet_2d.at(ihist).first != histo_zjet_2d.at(idenom).first){
      idenom = ihist;
      continue;
    }

    
    TString name (histo_zjet_2d.at(ihist).second->GetName());
    size_t idenom_temp = idenom;
    TString name2;      
    if(name.Contains("_cc") and not name.Contains("Inclusive")){
      for(size_t jhist = 1; jhist < histo_zjet_2d.size(); jhist++){
	name2 = TString(histo_zjet_2d.at(jhist).second->GetName());
	if(name2.Contains("_cc") and name2.Contains("Inclusive")){
	  idenom_temp = jhist;
	  break;
	}
      }
    }
    else if(name.Contains("_cf") and not name.Contains("Inclusive")){
      for(size_t jhist = 1; jhist < histo_zjet_2d.size(); jhist++){
	name2 = TString(histo_zjet_2d.at(jhist).second->GetName());
	if(name2.Contains("_cf") and name2.Contains("Inclusive")){
	  idenom_temp = jhist;
	  break;
	}
      }
    }
    else if(name.Contains("_ff") and not name.Contains("Inclusive")){
      for(size_t jhist = 1; jhist < histo_zjet_2d.size(); jhist++){
	name2 = TString(histo_zjet_2d.at(jhist).second->GetName());
	if(name2.Contains("_ff") and name2.Contains("Inclusive")){
	  idenom_temp = jhist;
	  break;
	}
      }
    }

    name.ReplaceAll("h_","eff_");
    eff_zjet_2d.push_back(pair<string,TEfficiency*>(histo_zjet_2d.at(ihist).first,new TEfficiency(*histo_zjet_2d.at(ihist).second,*histo_zjet_2d.at(idenom_temp).second)));
    eff_zjet_2d.back().second->SetName(name.Data());

    name = TString(histo_wjet_2d.at(ihist).second->GetName());
    name.ReplaceAll("h_","eff_");
    eff_wjet_2d.push_back(pair<string,TEfficiency*>(histo_wjet_2d.at(ihist).first,new TEfficiency(*histo_wjet_2d.at(ihist).second,*histo_wjet_2d.at(idenom_temp).second)));
    eff_wjet_2d.back().second->SetName(name.Data());
    
    name = TString(histo_wmn_2d.at(ihist).second->GetName());
    name.ReplaceAll("h_","eff_");
    eff_wmn_2d.push_back(pair<string,TEfficiency*>(histo_wmn_2d.at(ihist).first,new TEfficiency(*histo_wmn_2d.at(ihist).second,*histo_wmn_2d.at(idenom_temp).second)));
    eff_wmn_2d.back().second->SetName(name.Data());

    name = TString(histo_zmm_2d.at(ihist).second->GetName());
    name.ReplaceAll("h_","eff_");
    eff_zmm_2d.push_back(pair<string,TEfficiency*>(histo_zmm_2d.at(ihist).first,new TEfficiency(*histo_zmm_2d.at(ihist).second,*histo_zmm_2d.at(idenom_temp).second)));
    eff_zmm_2d.back().second->SetName(name.Data());
    
    eff_zjet_2d.back().second->Write();
    eff_wjet_2d.back().second->Write();
    eff_wmn_2d.back().second->Write();
    eff_zmm_2d.back().second->Write();
    
  }

  for(size_t ihist = 0; ihist < eff_zjet_2d.size(); ihist++){
    
    TH2* zjet_2d = eff_zjet_2d.at(ihist).second->CreateHistogram();
    zjet_2d->SetName(eff_zjet_2d.at(ihist).second->GetName());
    TH2* wjet_2d = eff_wjet_2d.at(ihist).second->CreateHistogram();
    wjet_2d->SetName(eff_wjet_2d.at(ihist).second->GetName());
    TH2* wmn_2d = eff_wmn_2d.at(ihist).second->CreateHistogram();
    wmn_2d->SetName(eff_wmn_2d.at(ihist).second->GetName());
    TH2* zmm_2d = eff_zmm_2d.at(ihist).second->CreateHistogram();
    zmm_2d->SetName(eff_zmm_2d.at(ihist).second->GetName());

    if(eff_zjet_2d.at(ihist).first == "mjj-ptj1"){
      
      TProfile* profile_zjet_x = zjet_2d->ProfileX(Form("%s_mjj",zjet_2d->GetName()));
      TProfile* profile_zjet_y = zjet_2d->ProfileY(Form("%s_ptj1",zjet_2d->GetName()));
      TProfile* profile_wjet_x = wjet_2d->ProfileX(Form("%s_mjj",wjet_2d->GetName()));
      TProfile* profile_wjet_y = wjet_2d->ProfileY(Form("%s_ptj1",wjet_2d->GetName()));
      TProfile* profile_wmn_x = wmn_2d->ProfileX(Form("%s_mjj",wmn_2d->GetName()));
      TProfile* profile_wmn_y = wmn_2d->ProfileY(Form("%s_ptj1",wmn_2d->GetName()));
      TProfile* profile_zmm_x = zmm_2d->ProfileX(Form("%s_mjj",zmm_2d->GetName()));
      TProfile* profile_zmm_y = zmm_2d->ProfileY(Form("%s_ptj1",zmm_2d->GetName()));

      plotEfficiency(canvas,profile_zjet_x,profile_wjet_x,profile_wmn_x,profile_zmm_y,outputDIR,luminosity,"mjj");
      plotEfficiency(canvas,profile_zjet_y,profile_wjet_y,profile_wmn_y,profile_zmm_y,outputDIR,luminosity,"ptj1");
      
    }
    else if(eff_wmn_2d.at(ihist).first == "mjj-ptj2"){

      TProfile* profile_zjet_x = zjet_2d->ProfileX(Form("%s_mjj",zjet_2d->GetName()));
      TProfile* profile_zjet_y = zjet_2d->ProfileX(Form("%s_ptj2",zjet_2d->GetName()));
      TProfile* profile_wjet_x = wjet_2d->ProfileX(Form("%s_mjj",wjet_2d->GetName()));
      TProfile* profile_wjet_y = wjet_2d->ProfileX(Form("%s_ptj2",wjet_2d->GetName()));
      TProfile* profile_wmn_x = wmn_2d->ProfileX(Form("%s_mjj",wmn_2d->GetName()));
      TProfile* profile_wmn_y = wmn_2d->ProfileX(Form("%s_ptj2",wmn_2d->GetName()));
      TProfile* profile_zmm_x = wmn_2d->ProfileX(Form("%s_mjj",zmm_2d->GetName()));
      TProfile* profile_zmm_y = wmn_2d->ProfileX(Form("%s_ptj2",zmm_2d->GetName()));

      plotEfficiency(canvas,profile_zjet_x,profile_wjet_x,profile_wmn_x,profile_zmm_y,outputDIR,luminosity,"mjj");
      plotEfficiency(canvas,profile_zjet_y,profile_wjet_y,profile_wmn_y,profile_zmm_y,outputDIR,luminosity,"ptj2");      
    }			 
    else if(eff_wmn_2d.at(ihist).first == "mjj-detajj"){
      
      TProfile* profile_zjet_x = zjet_2d->ProfileX(Form("%s_mjj",zjet_2d->GetName()));
      TProfile* profile_zjet_y = zjet_2d->ProfileY(Form("%s_detajj",zjet_2d->GetName()));
      TProfile* profile_wjet_x = wjet_2d->ProfileX(Form("%s_mjj",wjet_2d->GetName()));
      TProfile* profile_wjet_y = wjet_2d->ProfileY(Form("%s_detajj",wjet_2d->GetName()));
      TProfile* profile_wmn_x = wmn_2d->ProfileX(Form("%s_mjj",wmn_2d->GetName()));
      TProfile* profile_wmn_y = wmn_2d->ProfileY(Form("%s_detajj",wmn_2d->GetName()));
      TProfile* profile_zmm_x = zmm_2d->ProfileX(Form("%s_mjj",zmm_2d->GetName()));
      TProfile* profile_zmm_y = zmm_2d->ProfileY(Form("%s_detajj",zmm_2d->GetName()));
      
      plotEfficiency(canvas,profile_zjet_x,profile_wjet_x,profile_wmn_x,profile_zmm_x,outputDIR,luminosity,"mjj");
      plotEfficiency(canvas,profile_zjet_y,profile_wjet_y,profile_wmn_y,profile_zmm_y,outputDIR,luminosity,"detajj");      
    }
  }
}
