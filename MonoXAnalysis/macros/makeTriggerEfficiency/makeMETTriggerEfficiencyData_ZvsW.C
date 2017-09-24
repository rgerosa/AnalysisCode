#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

// recoil binning for monojet                                                                                                                                                                    
vector <float> bins_monojet_recoil    = {150,160,170,180,190,200,210,220,230,240,250,265,280,300,320,340,360,380,400,430,460,490,520,550,580,610,650,700,740,800,900,1000,1250};
vector <float> bins_VBFrelaxed_recoil = {150,160,170,180,190,200,210,220,230,240,250,265,280,300,320,340,360,380,400,430,460,490,520,550,580,610,650,700,740,800,900,1000,1250};
vector <float> bins_VBF_recoil        = {150,175,200,225,230,250,275,300,350,400,450,500,600,700,850,1000};
vector <float> bins_VBFrelaxed_mjj    = {200,400,600,800,1000,1250,1500,1750,2000,2500,3500,5000};
vector <float> bins_VBF_mjj           = {1300,1500,1750,2000,3500,3500,5000};

vector <float> bins_VBFrelaxed_mjj_2d  = {200,600,1000,1500,2000,3000,5000};
vector <float> bins_VBFrelaxed_ptj1_2d = {80,250,350,600};
vector <float> bins_VBFrelaxed_ptj2_2d = {40,80,150,300};
vector <float> bins_VBFrelaxed_detajj_2d = {1,3,4.5,10};

// cut for trigger values
static float L1_ETM_CUT  = 70;
static float caloMET_CUT = 70;
static float caloMETClean_CUT = 60;
static float caloMHT_CUT = 70;
static float PFMHT_CUT   = 90;
static float PFMET_CUT   = 90;
static float PFMHTNoMu_CUT = 90;
static float PFMETNoMu_CUT = 90;
static bool  doJetStudy = false;

// possible samples to be selected
enum class Sample   {wmn,zmm};
enum class Category {monojet,VBFrelaxed,VBF};


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
void createHistograms2D(vector<pair<string,vector<TH1F* > > > & histo, const string & postfix, const vector<float> & binningX, const vector<float> & binningY, const string & obs){


  vector<TH1F*> binVector;
  
  for(int iBin = 0; iBin < binningY.size()-1; iBin++){
    binVector.push_back(new TH1F(Form("h_Inclusive_%s_%s_bin_%d",obs.c_str(),postfix.c_str(),iBin), "", binningX.size()-1, &binningX[0]));
    binVector.back()->Sumw2();
    binVector.back()->SetDirectory(0);
    binVector.push_back(new TH1F(Form("h_BigTriggerOR_%s_%s_bin_%d",obs.c_str(),postfix.c_str(),iBin), "", binningX.size()-1, &binningX[0]));
    binVector.back()->Sumw2();
    binVector.back()->SetDirectory(0);
  }
  
  histo.push_back(pair<string,vector<TH1F* > >(obs,binVector));
    
}


//////////
void makeTriggerAnalysis(TChain* tree, 
			 vector<pair<string,TH1F*> > &  histograms, // depends on the selection one wants to apply
			 vector<pair<string,vector<TH1F*> > > &  histograms2D, // depends on the selection one wants to apply                                                                                   
 			 const Sample &   sample, // sample to be selected
			 const Category & category,
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

  TTreeReader reader(tree);
  TTreeReaderValue<UChar_t> hltmetwithmu (reader,triggerString.c_str());
  TTreeReaderValue<UChar_t> hltmet (reader,triggerStringNoMu.c_str());
  TTreeReaderValue<UChar_t> hltsinglemu (reader,"hltsinglemu");
  TTreeReaderValue<UChar_t> hltdoublemu (reader,"hltdoublemu");
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
  TTreeReaderValue<float> trig_L1ETM_pt (reader,"trig_L1ETM_pt");
  TTreeReaderValue<vector<float> > trig_obj_pt (reader,"trig_obj_pt");
  TTreeReaderValue<vector<string> > trig_obj_col (reader,"trig_obj_col");
  
  //////////////////                                                                                                                                                                              
  long int nTotal = tree->GetEntries();
  long int nEvents = 0;    
  long int nPart = 100000;
  
  while(reader.Next()){/////
    
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
    if(sample == Sample::wmn and *nmuons !=1) continue;
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
      
      float value = 0;
      if(obs.first == "met")
	value = *mmet;
      else if(obs.first == "mjj")
	value = (jet1+jet2).M();
      
      // fix overflow
      if(value > obs.second->GetXaxis()->GetBinLowEdge(obs.second->GetNbinsX()+1))
	value = obs.second->GetXaxis()->GetBinLowEdge(obs.second->GetNbinsX()+1)-0.001;
      
      TString name (obs.second->GetName());
      if(name.Contains("Inclusive_cc")){
	if(fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) < 3)
            obs.second->Fill(value);
      }
      else if(name.Contains("Inclusive_cf")){
	if(((fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) > 3) or (fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) < 3)))
	  obs.second->Fill(value);
      }
      else if(name.Contains("Inclusive_ff")){ // fill with all the events i.e. numerator                                                                                                              
	if(fabs(jet1.Eta()) > 3.0 and fabs(jet2.Eta()) > 3.)
	  obs.second->Fill(value);
      }
      else if(name.Contains("Inclusive")) // fill with all the events i.e. numerator                                                                                                                 
	obs.second->Fill(value);	
      // L1 efficiency
      else if(name.Contains("L1ETM")){ // fill with all the events passing L1 selection
	if(*trig_L1ETM_pt > L1_ETM_CUT) 
	  obs.second->Fill(value);
      }
      
      // L1 + caloMET + caloMETClean
      else if(name.Contains("caloMETClean")){// fill with all the events passing calo-MET + calo-MET clean cut
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT) 
	  obs.second->Fill(value);
      }
      
      // L1 + caloMET
      else if(name.Contains("caloMET")){ // fill with all the events passing calo-MET 
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT)
	  obs.second->Fill(value);
      }
      
      // L1 + caloMET + caloMETClean + caloMHT
      else if(name.Contains("caloMHT")){// fill with all the events passing calo MHT
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT and
	   caloMHT > caloMHT_CUT)
	  obs.second->Fill(value);
      }
      
      // L1 + caloMET + caloMETClean + caloMHT + PFMHTNoMu
      else if(name.Contains("PFMHTNoMu")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and 
	   caloMETClean > caloMETClean_CUT and
	   caloMHT > caloMHT_CUT and
	   PFMHTNoMu > PFMHTNoMu_CUT)
	  obs.second->Fill(value);
      }
      
      // L1 + caloMET + caloMETClean + caloMHT + PFMHT
      else if(name.Contains("PFMHT")){  // fill with all the events passing PFMHT
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT and
	   caloMHT > caloMHT_CUT and
	   PFMHT > PFMHT_CUT)
	  obs.second->Fill(value);
      }
      
      // L1 + caloMET + caloMETClean + caloMHT + PFMHTNoMu + PFMETNoMu
      else if(name.Contains("PFMETNoMu")){
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT and
	   caloMHT   > caloMHT_CUT and
	   PFMHTNoMu > PFMHTNoMu_CUT and
	   PFMETNoMu > PFMETNoMu_CUT)
	  obs.second->Fill(value);	    
      }	  
      
      // L1 + caloMET + caloMETClean + caloMHT + PFMHT + PFMET
      else if(name.Contains("PFMET")){ // fill with all the events passing PFMET
	if(*trig_L1ETM_pt > L1_ETM_CUT and
	   caloMET > caloMET_CUT and
	   caloMETClean > caloMETClean_CUT and
	   caloMHT > caloMHT_CUT and
	   PFMHT > PFMHT_CUT and
	   PFMET > PFMET_CUT) 
	  obs.second->Fill(value);	  
      } 
      
      // PFMET path
      else if(name.Contains("PFMETPath")){
	if(*hltmetwithmu != 0) 
	  obs.second->Fill(value);	  
      } 
      
      // PFMETNoMu path
      else if(name.Contains("PFMETNoMuPath")){
	if(*hltmet != 0) 
	  obs.second->Fill(value);	  
      } 

      // PFMETNoMu path
      else if(name.Contains("TotalPath")){
	if((*hltmet != 0 or *hltmetwithmu != 0)) 
	  obs.second->Fill(value);	  
      } 
      
      else if(name.Contains("h_BigTriggerOR_cc")){
	if((*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	    *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	    *hltjetmet != 0) and fabs(jet1.Eta()) < 3.0 and fabs(jet2.Eta()) < 3.0)
	  obs.second->Fill(value);
      }
      
      else if(name.Contains("h_BigTriggerOR_cf")){
	if((*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	    *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	    *hltjetmet != 0) and ((fabs(jet1.Eta()) < 3.0 and fabs(jet2.Eta()) > 3.0) or (fabs(jet1.Eta()) > 3.0 and fabs(jet2.Eta()) < 3.0)))
	  obs.second->Fill(value);
        }
      
      else if(name.Contains("h_BigTriggerOR_ff")){
	if((*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	    *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	    *hltjetmet != 0) and fabs(jet1.Eta()) > 3.0 and fabs(jet2.Eta()) > 3.)
	  obs.second->Fill(value);
      }
      
      else if(name.Contains("h_BigTriggerOR")){
	if((*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	    *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	    *hltjetmet != 0))
	  obs.second->Fill(value);
      }	
    }

    // Repeat the same for 2D histograms                                                                                                                                                    
    for(auto obs : histograms2D){

      if(*mmet < 250) continue; // cut harder in met                                                                                                                                
      // values                                                                                                                                                                                       
      float valueX = 0;
      float valueY = 0;

      vector<float> binningY;

      if(obs.first == "mjj-ptj1"){
        valueX = (jet1+jet2).M();
        valueY = jet1.Pt();
	binningY = bins_VBFrelaxed_ptj1_2d;
      }
      else if(obs.first == "mjj-ptj2"){
        valueX = (jet1+jet2).M();
        valueY = jet2.Pt();
	binningY = bins_VBFrelaxed_ptj2_2d;
      }
      else if(obs.first == "mjj-detajj"){
        valueX = (jet1+jet2).M();
        valueY = fabs(jet1.Eta()-jet2.Eta());
	binningY = bins_VBFrelaxed_detajj_2d;
      }
      
      if(valueY > binningY.back())
	valueY = binningY.back()-0.001;
      
      /// decide the histogram
      int binToFill = -1;
      for(int iBinY = 0; iBinY < binningY.size()-1; iBinY++){
	if(valueY >= binningY.at(iBinY) and valueY <= binningY.at(iBinY+1)) binToFill = iBinY;
      }

      TString stringToSearch = Form("_bin_%d",binToFill);

      for(auto histo : obs.second){

	if(not TString(histo->GetName()).Contains(stringToSearch)) continue;

	// handling overflow                                                                                                                                                                           
	if(valueX > histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1))
	  valueX = histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1)-1;

	TString name (histo->GetName());

	if(name.Contains("Inclusive_cc")){
	  if(fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) < 3)
	    histo->Fill(valueX);
	}
	if(name.Contains("Inclusive_cf")){
	  if(((fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) > 3) or (fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) < 3)))
	    histo->Fill(valueX);
	}
	if(name.Contains("Inclusive_ff")){ // fill with all the events i.e. numerator                                                                                                                   
	  if(fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) > 3)
	    histo->Fill(valueX);
	}
	else if(name.Contains("Inclusive")) // fill with all the events i.e. numerator                                                                                                                  
	  histo->Fill(valueX);
	
	// L1 efficiency                                                                                                                                                                           
	else if(name.Contains("L1ETM")){ // fill with all the events passing L1 selection                                                                                                
	  if(*trig_L1ETM_pt > L1_ETM_CUT)
	    histo->Fill(valueX);
	}
      // L1 + caloMET + caloMETClean                                                                                                                                                            
	else if(name.Contains("caloMETClean")){// fill with all the events passing calo-MET + calo-MET clean cut                                                                        
        if(*trig_L1ETM_pt > L1_ETM_CUT and
           caloMET > caloMET_CUT and
           caloMETClean > caloMETClean_CUT)
          histo->Fill(valueX);
	}
	
	// L1 + caloMET                                                                                                                                                                         
	else if(name.Contains("caloMET")){ // fill with all the events passing calo-MET                                                                                                  
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT)
	    histo->Fill(valueX);
	}
	
	// L1 + caloMET + caloMETClean + caloMHT                                                                                                                                                    
	else if(name.Contains("caloMHT")){// fill with all the events passing calo MHT                                                                                                      
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT and
	     caloMETClean > caloMETClean_CUT and
	     caloMHT > caloMHT_CUT)
	    histo->Fill(valueX);
	}
	
	// L1 + caloMET + caloMETClean + caloMHT + PFMHTNoMu                                                                                                                                        
	else if(name.Contains("PFMHTNoMu")){
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT and
	     caloMETClean > caloMETClean_CUT and
	     caloMHT > caloMHT_CUT and
	     PFMHTNoMu > PFMHTNoMu_CUT)
	    histo->Fill(valueX);
	}
	
	// L1 + caloMET + caloMETClean + caloMHT + PFMHT                                                                                                                                           
	else if(name.Contains("PFMHT")){  // fill with all the events passing PFMHT                                                                                                             
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT and
	     caloMETClean > caloMETClean_CUT and
           caloMHT > caloMHT_CUT and
	     PFMHT > PFMHT_CUT)
	    histo->Fill(valueX);
	}
	
	// L1 + caloMET + caloMETClean + caloMHT + PFMHTNoMu + PFMETNoMu                                                                                                                       
	else if(name.Contains("PFMETNoMu")){
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT and
	     caloMETClean > caloMETClean_CUT and
	     caloMHT   > caloMHT_CUT and
	     PFMHTNoMu > PFMHTNoMu_CUT and
	     PFMETNoMu > PFMETNoMu_CUT)
	    histo->Fill(valueX);
	}
	
	// L1 + caloMET + caloMETClean + caloMHT + PFMHT + PFMET                                                                                                                                 
	else if(name.Contains("PFMET")){ // fill with all the events passing PFMET                                                                                                              
	  if(*trig_L1ETM_pt > L1_ETM_CUT and
	     caloMET > caloMET_CUT and
	     caloMETClean > caloMETClean_CUT and
	     caloMHT > caloMHT_CUT and
	     PFMHT > PFMHT_CUT and
	     PFMET > PFMET_CUT)
          histo->Fill(valueX);
	}
	
	// PFMET path                                                                                                                                                                                  
	else if(name.Contains("PFMETPath")){
	  if(*hltmetwithmu != 0)
	    histo->Fill(valueX);
	}
	// PFMETNoMu path                                                                                                                                                                             
	else if(name.Contains("PFMETNoMuPath")){
	  if(*hltmet != 0)
	    histo->Fill(valueX);
	}

	else if(name.Contains("TotalPath")){
	  if(*hltmet != 0 or *hltmetwithmu != 0)
	    histo->Fill(valueX);
	}

	/////// ---                                                                                                                                                                                   
	else if(name.Contains("h_BigTriggerOR_cc")){
	  if((*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	      *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	      *hltjetmet != 0) and fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) < 3)
	    histo->Fill(valueX);
	}
	else if(name.Contains("h_BigTriggerOR_cf")){
	  if((*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	      *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	      *hltjetmet != 0) and ((fabs(jet1.Eta()) < 3 and fabs(jet2.Eta()) > 3) or (fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) < 3)))
	    histo->Fill(valueX);
	}
	
	else if(name.Contains("h_BigTriggerOR_ff")){
	  if((*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	      *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	      *hltjetmet != 0) and fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) > 3.)
	    histo->Fill(valueX);
	}
	
	else if(name.Contains("h_BigTriggerOR")){
	  if((*hltmet90 != 0 or *hltmet100 != 0 or *hltmet110 != 0 or *hltmet120 != 0 or
	      *hltmetwithmu90 != 0 or *hltmetwithmu100 != 0 or *hltmetwithmu110 != 0 or *hltmetwithmu120 != 0 or
	      *hltjetmet != 0))
	    histo->Fill(valueX);
	}
      }            
    }
  }
  cout<<endl;
}


//// -----
void fillTreeList(TChain* tree, const string & inputDIR, const string & postfix){

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
	  cout<<"Fill chain: "<<postfix<<" file with name: "<<line2<<endl;
	  tree->Add(line2.c_str());
	}
      }
      system("rm list.txt");
    }
  }
  system("rm list_dir.txt");
}

/// plotting efficiency
void plotEfficiency(TCanvas* canvas, 
		    TEfficiency* histo_wmn, 
		    TEfficiency* histo_zmm, 
		    const string & outputDIR, 
		    const float  & luminosity, 
		    const string & observable){
  
  // Make efficiencies as TGraph
  TGraphAsymmErrors* graph_wmn  = histo_wmn->CreateGraph();
  TGraphAsymmErrors* graph_zmm  = histo_zmm->CreateGraph();

  graph_wmn->SetMarkerColor(kRed);
  graph_wmn->SetLineColor(kRed);
  graph_wmn->SetMarkerStyle(24);
  graph_wmn->SetMarkerSize(0.75);
  graph_wmn->SetLineWidth(1);

  graph_zmm->SetMarkerColor(kBlue);
  graph_zmm->SetLineColor(kBlue);
  graph_zmm->SetMarkerStyle(20);
  graph_zmm->SetMarkerStyle(26);
  graph_zmm->SetMarkerSize(0.75);
  graph_zmm->SetLineWidth(1);
    
  // Plotting final result for MC turn ons
  TH1* frame = NULL;
  if(TString(observable).Contains("mjj"))
    frame = canvas->DrawFrame(histo_wmn->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(1),0.2,
			      histo_wmn->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(histo_wmn->GetPassedHistogram()->GetNbinsX()+1), 1.1, "");
  else 
    frame = canvas->DrawFrame(histo_wmn->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(1),0.6,
			      histo_wmn->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(histo_wmn->GetPassedHistogram()->GetNbinsX()+1), 1.1, "");


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

  graph_wmn->Draw("EPLsame");
  graph_zmm->Draw("EPLsame");
  
  TLegend* leg = NULL;
  if(TString(observable).Contains("met"))
    leg = new TLegend(0.6,0.4,0.9,0.52);
  else 
    leg = new TLegend(0.25,0.4,0.55,0.52);
  
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(graph_wmn,"W(#mu#nu)-CR","EPL");
  leg->AddEntry(graph_zmm,"Z(#mu#mu)-CR","EPL");
  leg->Draw("same");

  TString plotName(histo_wmn->GetName());
  plotName.ReplaceAll("_wmn","");

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
    frame2 = pad2->DrawFrame(histo_wmn->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(1),0.8,
			     histo_wmn->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(histo_wmn->GetPassedHistogram()->GetNbinsX()+1), 1.05, "");
  else
    frame2 = pad2->DrawFrame(histo_wmn->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(1),0.9,
			     histo_wmn->GetPassedHistogram()->GetXaxis()->GetBinLowEdge(histo_wmn->GetPassedHistogram()->GetNbinsX()+1), 1.05, "");
  
  if(TString(observable).Contains("met"))
    frame2->GetXaxis()->SetTitle("Recoil [GeV]");
  else if(TString(observable).Contains("mjj"))
    frame2->GetXaxis()->SetTitle("M_{jj} [GeV]");

  frame2->GetYaxis()->SetTitle("Ratio");
  frame2->GetXaxis()->SetTitleOffset(1.1);
  frame2->GetYaxis()->SetTitleOffset(1.32);
  frame2->GetYaxis()->SetTitleSize(0.04);
  frame2->GetYaxis()->SetLabelSize(0.03);
  frame2->GetYaxis()->SetNdivisions(504);
  frame2->GetYaxis()->SetLabelSize(0.85*frame2->GetYaxis()->GetLabelSize());
  frame2->Draw();

  /// ------
  TH1F* ratio_zmm_over_wmn = (TH1F*) histo_zmm->GetPassedHistogram();
  ratio_zmm_over_wmn->Divide(histo_zmm->GetTotalHistogram());
  TH1F* ratio_temp = (TH1F*) histo_wmn->GetPassedHistogram();
  ratio_temp->Divide(histo_wmn->GetTotalHistogram());
  ratio_zmm_over_wmn->Divide(ratio_temp);
  
  ratio_zmm_over_wmn->SetMarkerColor(kBlue);
  ratio_zmm_over_wmn->SetLineColor(kBlue);
  ratio_zmm_over_wmn->SetMarkerStyle(24);
  ratio_zmm_over_wmn->SetMarkerSize(0.75);
  ratio_zmm_over_wmn->SetLineWidth(1);
  ratio_zmm_over_wmn->Draw("EPsame");

  canvas->SaveAs((outputDIR+"/"+string(plotName.Data())+".pdf").c_str(),"pdf");

  if(frame2) delete frame2;
  if(pad2) delete pad2;
  if(leg) delete leg;
}



/// Main function 
void makeMETTriggerEfficiencyData_ZvsW(string   inputDIR, 
				       string   outputDIR, 
				       Category category,
				       string   triggerPath, // PFMET90, PFMET100, PFMET110, PFMET120
				       vector<string>   observable, // study the trigger efficiency vs mjj or recoil
				       vector<string>   observable2D, // study the trigger efficiency vs mjj or recoil
				       bool     splitJets = false,
				       float    luminosity = 35.9){

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  doJetStudy = splitJets;

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
  
  // input tree                                                                                                                                                                                       
  TChain* tree_singleMu = new TChain("tree/tree","tree/tree");
  fillTreeList(tree_singleMu,inputDIR,"SingleMuon");
 
  vector<pair<string,TH1F*> > histo_wmn;
  vector<pair<string,TH1F*> > histo_zmm;

  /////// start analysis  
  for(auto obs: observable){
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

    createHistograms1D(histo_wmn,"wmn",binning,obs);
    createHistograms1D(histo_zmm,"zmm",binning,obs);
  }


  ////// start analysis                                                                                                                                                                                
  vector<pair<string,vector<TH1F*> > > histo_wmn_2d;
  vector<pair<string,vector<TH1F*> > > histo_zmm_2d;

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
    createHistograms2D(histo_wmn_2d,"wmn",binningX,binningY,obs);
    createHistograms2D(histo_zmm_2d,"zmm",binningX,binningY,obs);
  }

  //// ----- 
  cout<<"######### Loop on W+jets trees for Wmn selection "<<endl;
  makeTriggerAnalysis(tree_singleMu,histo_wmn,histo_wmn_2d,Sample::wmn,category,triggerPath);
  cout<<"######### Loop on Z+jets trees for Zmm selection "<<endl;
  makeTriggerAnalysis(tree_singleMu,histo_zmm,histo_zmm_2d,Sample::zmm,category,triggerPath);
  
  ////---
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 650);
  canvas->cd();

  TFile* outputFile = new TFile((outputDIR+"/efficiency.root").c_str(),"RECREATE");
  outputFile->cd();

  // calculate the efficiency
  vector<pair<string,TEfficiency* > > eff_wmn;
  vector<pair<string,TEfficiency* > > eff_zmm;

  size_t idenom = 0;
  for(size_t ihist = 1; ihist < histo_wmn.size(); ihist++){

    if(histo_wmn.at(ihist).first != histo_wmn.at(idenom).first){
      idenom = ihist;
      continue;
    }

    TString name (histo_wmn.at(ihist).second->GetName());
    size_t idenom_temp = idenom;
    TString name2;
    if(name.Contains("_cc") and not name.Contains("Inclusive")){
      for(size_t jhist = 1; jhist < histo_wmn.size(); jhist++){
        name2 = TString(histo_wmn.at(jhist).second->GetName());
        if(name2.Contains("_cc") and name2.Contains("Inclusive")){
          idenom_temp = jhist;
          break;
        }
      }
    }
    else if(name.Contains("_cf") and not name.Contains("Inclusive")){
      for(size_t jhist = 1; jhist < histo_wmn.size(); jhist++){
        name2 = TString(histo_wmn.at(jhist).second->GetName());
        if(name2.Contains("_cf") and name2.Contains("Inclusive")){
          idenom_temp = jhist;
          break;
        }
      }
    }
    else if(name.Contains("_ff") and not name.Contains("Inclusive")){
      for(size_t jhist = 1; jhist < histo_wmn.size(); jhist++){
        name2 = TString(histo_wmn.at(jhist).second->GetName());
        if(name2.Contains("_ff") and name2.Contains("Inclusive")){
          idenom_temp = jhist;
          break;
        }
      }
    }
      
    name.ReplaceAll("h_","eff_");
    eff_wmn.push_back(pair<string,TEfficiency*>(histo_wmn.at(ihist).first,new TEfficiency(*histo_wmn.at(ihist).second,*histo_wmn.at(idenom_temp).second)));
    eff_wmn.back().second->SetName(name.Data());
    
    name = TString(histo_zmm.at(ihist).second->GetName());
    name.ReplaceAll("h_","eff_");
    eff_zmm.push_back(pair<string,TEfficiency*>(histo_wmn.at(ihist).first,new TEfficiency(*histo_zmm.at(ihist).second,*histo_zmm.at(idenom_temp).second)));
    eff_zmm.back().second->SetName(name.Data());      

    ////-----
    eff_wmn.back().second->Write();
    eff_zmm.back().second->Write();
    
  }

  ////// ---- plot Efficiency
  for(size_t ihist = 0; ihist < eff_wmn.size(); ihist++)
    plotEfficiency(canvas,eff_wmn.at(ihist).second,eff_zmm.at(ihist).second,outputDIR,luminosity,eff_wmn.at(ihist).first);

  // 2D efficiencies                                                                                                                                                                                
  vector<pair<string,vector<TEfficiency* > > > eff_wmn_2d;
  vector<pair<string,vector<TEfficiency* > > > eff_zmm_2d;

  for(size_t iobs = 0; iobs < histo_wmn_2d.size(); iobs++){ // loop on the observable

    vector<TH1F*> denom_wmn_eff;
    vector<TH1F*> num_wmn_eff;
    
    for(size_t ihist = 0; ihist < histo_wmn_2d.at(iobs).second.size(); ihist++){ // loop on the histogram
      TString name (histo_wmn_2d.at(iobs).second.at(ihist)->GetName());
      if(name.Contains("Inclusive"))
	denom_wmn_eff.push_back(histo_wmn_2d.at(iobs).second.at(ihist));
      else
	num_wmn_eff.push_back(histo_wmn_2d.at(iobs).second.at(ihist));	      
    }

    vector<TH1F*> denom_zmm_eff;
    vector<TH1F*> num_zmm_eff;

    for(size_t ihist = 0; ihist < histo_zmm_2d.at(iobs).second.size(); ihist++){ // loop on the histogram
      TString name (histo_zmm_2d.at(iobs).second.at(ihist)->GetName());
      if(name.Contains("Inclusive")) denom_zmm_eff.push_back(histo_zmm_2d.at(iobs).second.at(ihist));
      else num_zmm_eff.push_back(histo_zmm_2d.at(iobs).second.at(ihist));	      
    }

    // taking into account the order in which histogram vector has been filled
    vector<TEfficiency*> eff_temp_wmn;
    vector<TEfficiency*> eff_temp_zmm;

    for(size_t ihist = 0; ihist < num_wmn_eff.size(); ihist++){
      TString name (num_wmn_eff.at(ihist)->GetName());
      name.ReplaceAll("h_","eff_");
      if(ihist <= denom_wmn_eff.size()-1){
	eff_temp_wmn.push_back(new TEfficiency(*num_wmn_eff.at(ihist),*denom_wmn_eff.at(ihist)));
      }
      else{
	int pos = (ihist % denom_wmn_eff.size());
	eff_temp_wmn.push_back(new TEfficiency(*num_wmn_eff.at(ihist),*denom_wmn_eff.at(pos)));
      }
      eff_temp_wmn.back()->SetName(name.Data());
    }

    for(size_t ihist = 0; ihist < num_zmm_eff.size(); ihist++){
      TString name (num_zmm_eff.at(ihist)->GetName());
      name.ReplaceAll("h_","eff_");
      if(ihist <= denom_zmm_eff.size()-1)
	eff_temp_zmm.push_back(new TEfficiency(*num_zmm_eff.at(ihist),*denom_zmm_eff.at(ihist)));
      else{
	int pos = (ihist % denom_zmm_eff.size());
	eff_temp_zmm.push_back(new TEfficiency(*num_zmm_eff.at(ihist),*denom_zmm_eff.at(pos)));
      }
      eff_temp_zmm.back()->SetName(name.Data());
    }

    ////
    eff_wmn_2d.push_back(pair<string,vector<TEfficiency* > >(histo_wmn_2d.at(iobs).first,eff_temp_wmn));
    eff_zmm_2d.push_back(pair<string,vector<TEfficiency* > >(histo_zmm_2d.at(iobs).first,eff_temp_zmm));
  }

  ///////// ------ plot
  for(size_t iobs = 0; iobs < eff_wmn_2d.size(); iobs++){
    for(size_t ieff = 0; ieff < eff_wmn_2d.at(iobs).second.size(); ieff++){
      plotEfficiency(canvas,eff_wmn_2d.at(iobs).second.at(ieff),eff_zmm_2d.at(iobs).second.at(ieff),outputDIR,luminosity,eff_wmn_2d.at(iobs).first);
    }
  }
  
}
