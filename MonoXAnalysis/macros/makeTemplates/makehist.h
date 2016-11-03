#ifndef MAKEHIST_H
#define MAKEHIST_H

#include <vector>
#include <map>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TString.h"

//#include "histoUtils.h"
#include "histoUtils2D.h"

using namespace std;

// some basic cut values
const float tau2tau1        = 0.6;
const float tau2tau1LP      = 0.75;
const float prunedMassMin   = 65.;
const float prunedMassMax   = 105.;
const float leadingJetPtCut = 100.;
const float leadingJetPtCutVBF  = 70.;
const float trailingJetPtCutVBF = 50.;
const float detajj          = 2.5;
const float mjj             = 500;
const float jetmetdphiVBF   = 0.75;
const float ptJetMinAK8     = 250.;
const float jetEtaAK8       = 2.4;
const float pfMetMonoVLower = 250.;
const float pfMetMonoVUpper = 8000.;
const float pfMetMonoJUpper = 8000.;
const float pfMetMonoJLower = 200.;
const float pfMetVBFLower   = 200.;
const float pfMetVBFUpper   = 8000.;
const float photonPt        = 120;
const int   vBosonCharge    = 0;
const int   nBjets          = 1;
const bool  reweightNVTX    = true;
const bool  reweightPhton   = false;
const bool  applyPhotonScale = true;
const float photonScaleUnc  = 0.020;
const bool  doSmoothing     = false;
const float btagCSVLoose    = 0.460;
const float btagCSVMedium   = 0.800;
const float recoilThresholdTrigger = 350;

string kfactorFile       = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root";
string kfactorFileGJ     = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/photonjets_kfact.root";
string kfactorFileUnc    = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/scalefactors_v4.root";
string baseInputTreePath = "/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/";

VectorSorter jetSorter;

float deltaPhi (float phi1, float phi2){
  if(fabs(phi1-phi2) < TMath::Pi())
    return fabs(phi1-phi2);
  else
    return 2*TMath::Pi()-fabs(phi1-phi2);
}

double reweightTopQuarkPt (double topQuarkPt, double atopQuarkPt){

  double weight_top  = -1.;
  double weight_atop = -1.;
  if(topQuarkPt  != 0)
    weight_top = exp(0.156-0.00137*min(400.,topQuarkPt));
  if(atopQuarkPt != 0)
    weight_atop = exp(0.156-0.00137*min(400.,atopQuarkPt));

  if(weight_top != -1. and weight_atop != -1.)
    return 1.001*sqrt(weight_top*weight_atop);
  else return 1.;
}

double getVtaggingScaleFactor(const double & tau2tau1, const string & sysName){

  double sfwgt = 1;

  if(tau2tau1 == 0.45){
    if(sysName == "VtagUp")
      sfwgt *= (0.97+0.109);
    else if(sysName == "VtagDown")
      sfwgt *= (0.97-0.109);
    else
      sfwgt *= 0.97;
  }
  else if(tau2tau1 == 0.6){
    if(sysName == "VtagUp")
      sfwgt *= (0.97+0.109);
    else if(sysName == "VtagDown")
      sfwgt *= (0.97-0.109);
    else
      sfwgt *= 0.97;
  }
  
  return sfwgt;
}

void makehist4(TTree* tree, /*input tree*/ 
	       vector<TH1*> hist1D, /* set of 1D histogram */ 
	       vector<TH2*> hist2D, /* set of 2D histogram */ 
	       const bool &   isMC, 
	       const Sample & sample, 
	       const Category & category,
	       const bool &   isWJet,
	       const double & scale,
	       const double & lumi,	       
	       vector<TH1*> khists, 
	       const string & sysName,	
	       const bool   & reWeightTopPt = false,
	       const bool   & reweightNVTX  = true,
	       const int    & resonantSelection  = 0,
	       const bool   & isHiggsInvisible   = false, // reject VBF events
	       const bool   & applyPostFitWeight = false,
	       const float  & XSEC = -1.,// fix the cross section from extern
	       TH1*  hhist = NULL,
	       TH2*  ggZHhist = NULL
	       ) {

  if(not tree){
    cout<<" empty tree --> skip process "<<endl;
    return;
  }

  // in case you want to weight the NVTX distribution
  TFile* pufile = NULL;
  TH1*   puhist = NULL;
  if(reweightNVTX){
    pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_12p9fb.root");    
    puhist = (TH1*) pufile->Get("puhist");
  }
  else {
    pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/pileupWeight/puweight_12p9fb.root");
    puhist = (TH1*) pufile->Get("puhist");
  }
  
  // in case one wants to correct the photon pt by the data-to-mc disagreement
  TFile gamPtFile ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/extraNuisance/photon_nuisance.root");
  TH1*  gamPtWeight = (TH1*) gamPtFile.Get("gamma_jets_ratio");
     
  // electron and muon ID scale factor files                                                                                                                                    
  TFile sffile_eleTight ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/scaleFactor_electron_tightid_12p9.root");
  TFile sffile_eleVeto  ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/scaleFactor_electron_vetoid_12p9.root");
  TFile sffile_muTight  ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/scaleFactor_muon_tightid_12p9.root");
  TFile sffile_muLoose  ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/scaleFactor_muon_looseid_12p9.root");

  TH2*  msfloose = (TH2*) sffile_muLoose.Get("scaleFactor_muon_looseid_RooCMSShape");
  TH2*  msftight = (TH2*) sffile_muTight.Get("scaleFactor_muon_tightid_RooCMSShape");
  TH2*  esfveto  = (TH2*) sffile_eleVeto.Get("scaleFactor_electron_vetoid_RooCMSShape");
  TH2*  esftight = (TH2*) sffile_eleTight.Get("scaleFactor_electron_tightid_RooCMSShape");

  // Photon ID scale factor                                                                                                                                                     
  TFile sffile_phoMedium ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF_2016/scaleFactor_photon_mediumid_12p9.root");
  TH2*  psfmedium = (TH2*) sffile_phoMedium.Get("scaleFactor_photon_mediumid_RooCMSShape");

  // Photon Purity                                                                                                                                                              
  TFile purityfile_photon ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF_2016/PhotonSFandEffandPurity_Lumi12p9fb_28072016.root");
  TH2*  purhist = (TH2*) purityfile_photon.Get("purity");

  // Muon track efficiency
  TFile trackingefficiency_muon ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/trackEfficiency/scaleFactor_muon_trackerid.root");
  TH2F* trackingefficiency_muon_lowpu  = (TH2F*) trackingefficiency_muon.Get("scaleFactor_muon_trackerid_Exp_pu_0_16");
  TH2F* trackingefficiency_muon_highpu = (TH2F*) trackingefficiency_muon.Get("scaleFactor_muon_trackerid_Exp_pu_16_50");

  TFile  trackingefficiency_electron ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/trackEfficiency/scaleFactor_electron_recoelectronmatch.root");
  TH2F* trackingefficiency_electron_lowpu  = (TH2F*) trackingefficiency_electron.Get("scaleFactor_electron_recoelectronmatch_RooCMSShape_pu_0_16");
  TH2F* trackingefficiency_electron_highpu = (TH2F*) trackingefficiency_electron.Get("scaleFactor_electron_recoelectronmatch_RooCMSShape_pu_16_50");

  /////////////////////////////////////////
  // trigger files used for 2016                                                                                                                                                
  TFile triggerfile_SinglEle ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/triggerEfficiency_DATA_SingleElectron_12p9fb.root");
  TEfficiency* triggerel_eff = (TEfficiency*) triggerfile_SinglEle.Get("trgeff_ele");
  TH2* triggerelhist         = triggerel_eff->CreateHistogram();
  triggerelhist->SetName("triggerelhist");

  TFile triggerfile_SinglEle_jetHT ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/triggerEfficiency_DATA_SingleElectron_jetHT_12p9.root");
  TEfficiency* triggerel_eff_jetHT = (TEfficiency*) triggerfile_SinglEle_jetHT.Get("efficiency");
  TH2* triggerelhist_ht            = triggerel_eff_jetHT->CreateHistogram();
  triggerelhist_ht->SetName("triggerelhist_ht");

  // Met trigger efficiency
  TFile* triggerfile_MET = NULL;
  if(category != Category::VBF and category != Category::twojet)
    triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/metTriggerEfficiency_12p9.root");
  else
    triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/VBF/metTriggerEfficiency_electron.root");
  //triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/VBF/metTriggerEfficiency_muon.root");

  TEfficiency*  triggermet = (TEfficiency*) triggerfile_MET->Get("trig_eff");
  if(triggermet == 0 or triggermet == NULL)
    triggermet = (TEfficiency*) triggerfile_MET->Get("efficiency_vbf_loose");
  TGraphAsymmErrors* triggermet_graph = triggermet->CreateGraph();
  
  // Photon trigger efficiency measured in jetHT
  TFile triggerfile_SinglePhoton_jetHT ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/VBF/photonTriggerEfficiency_photonpt120_jetHT.root");
  TFile triggerfile_SinglePhoton       ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/VBF/photonTriggerEfficiency_photonpt120.root");
  TFile triggerfile_SinglePhoton_VBF   ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/VBF/photonTriggerEfficiency_photonpt120_VBF.root");
  TEfficiency*  triggerphoton       = NULL;
  TEfficiency*  triggerphoton_jetHT = NULL;

  if(category == Category::twojet){
    triggerphoton = (TEfficiency*) triggerfile_SinglePhoton.Get("eff_vbf_recoil");
    triggerphoton_jetHT = (TEfficiency*) triggerfile_SinglePhoton_jetHT.Get("eff_recover_vbf_recoil");
  }
  else if(category == Category::VBF){
    triggerphoton = (TEfficiency*) triggerfile_SinglePhoton_VBF.Get("eff_vbf_recoil");  
    triggerphoton_jetHT = (TEfficiency*) triggerfile_SinglePhoton_jetHT.Get("eff_recover_vbf_recoil");
  }
  else{
    triggerphoton = (TEfficiency*) triggerfile_SinglePhoton.Get("eff_recoil");  
    triggerphoton_jetHT = (TEfficiency*) triggerfile_SinglePhoton_jetHT.Get("eff_recover_vbf_recoil");
  }
  
  TGraphAsymmErrors* triggerphoton_graph = triggerphoton->CreateGraph();
  TGraphAsymmErrors* triggerphoton_graph_jetHT = triggerphoton_jetHT->CreateGraph();
  /////////////////////////////////////////

  // Post-fit weights
  TFile* postFitFile = NULL;
  if(sample == Sample::sig and applyPostFitWeight and (category == Category::monojet or category == Category::inclusive))
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ_2016_2p6fb/postfit_weights_Sig.root");
  if(sample == Sample::zmm and applyPostFitWeight and (category == Category::monojet or category == Category::inclusive))
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ_2016_2p6fb/postfit_weights_ZM.root");
  else if(sample == Sample::wmn and applyPostFitWeight and (category == Category::monojet or category == Category::inclusive))
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ_2016_2p6fb/postfit_weights_WM.root");
  else if(sample == Sample::zee and applyPostFitWeight and (category == Category::monojet or category == Category::inclusive))
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ_2016_2p6fb/postfit_weights_ZE.root");
  else if(sample == Sample::wen and applyPostFitWeight and (category == Category::monojet or category == Category::inclusive))
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ_2016_2p6fb/postfit_weights_WE.root");
  else if(sample == Sample::qcdgam and applyPostFitWeight and (category == Category::monojet or category == Category::inclusive))
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ_2016_2p6fb/postfit_weights_GJ.root");
  else if(sample == Sample::gam and applyPostFitWeight and (category == Category::monojet or category == Category::inclusive))
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ_2016_2p6fb/postfit_weights_GJ.root");
  else if(sample == Sample::sig and applyPostFitWeight and category == Category::monoV)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoV/postfit_weights_Sig.root");
  else if(sample == Sample::zmm and applyPostFitWeight and category == Category::monoV)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoV/postfit_weights_ZM.root");
  else if(sample == Sample::wmn and applyPostFitWeight and category == Category::monoV)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoV/postfit_weights_WM.root");
  else if(sample == Sample::zee and applyPostFitWeight and category == Category::monoV)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoV/postfit_weights_ZE.root");
  else if(sample == Sample::wen and applyPostFitWeight and category == Category::monoV)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoV/postfit_weights_WE.root");
  else if(sample == Sample::qcdgam and applyPostFitWeight and category == Category::monoV)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoV/postfit_weights_GJ.root");
  else if(sample == Sample::gam and applyPostFitWeight and category == Category::monoV)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoV/postfit_weights_GJ.root");

  TH1* postFitWeight = NULL;
  if(postFitFile != NULL){
    postFitWeight = (TH1*) postFitFile->FindObjectAny("post_fit_post_fit");
    postFitWeight->Divide((TH1*) postFitFile->FindObjectAny("pre_fit_post_fit"));
  }

  // histogram to be filled
  for(size_t ihist  = 0 ; ihist < hist1D.size(); ihist++){
    if(not hist1D.at(ihist)->GetSumw2N())
      hist1D.at(ihist)->Sumw2();
  }
  for(size_t ihist  = 0 ; ihist < hist2D.size(); ihist++){
    if(not hist1D.at(ihist)->GetSumw2N())
      hist2D.at(ihist)->Sumw2();
  }

  // define branches
  TTreeReader myReader(tree);

  // general info
  TTreeReaderValue<unsigned int> run    (myReader,"run");
  TTreeReaderValue<unsigned int> lumisection    (myReader,"lumi");
  TTreeReaderValue<unsigned int> event  (myReader,"event");
  TTreeReaderValue<unsigned int> nvtx   (myReader,"nvtx");
  TTreeReaderValue<int>          putrue (myReader,"putrue");
  TTreeReaderValue<double> xsec         (myReader,"xsec");
  TTreeReaderValue<double> wgt          (myReader,"wgt");

  // take some wgt for MC events
  string wgtname;
  string wgtpileupname;
  string btagname;
  string prescalename;
  string hltphotonname;

  if(isMC){
    // defined branches for re-weigthing only in MC
    wgtname       = "wgtsum";
    wgtpileupname = "wgtpileup";    
    // set b-tag
    if(sysName == "btagUp" or sysName == "bUp")
      btagname =  "wgtbtagUp";
    else if(sysName == "btagDown" or sysName == "btagDw" or sysName == "bDown" or sysName == "bDw")
      btagname =  "wgtbtagDown";
    else
      btagname    = "wgtbtag";
  }
  else if(!isMC){ // in case of data
    wgtname       = "wgt";
    btagname      = "wgt";
    wgtpileupname = "wgt";
  }

  TTreeReaderValue<double> wgtsum       (myReader,wgtname.c_str());
  TTreeReaderValue<double> wgtpileup    (myReader,wgtpileupname.c_str());
  TTreeReaderValue<double> wgtbtag      (myReader,btagname.c_str());
  
  // trigger
  TTreeReaderValue<UChar_t> hltm90     (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100    (myReader,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110    (myReader,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120    (myReader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120  (myReader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170  (myReader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300  (myReader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90   (myReader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hlte       (myReader,"hltsingleel");
  TTreeReaderValue<UChar_t> hltenoiso  (myReader,"hltelnoiso");
  TTreeReaderValue<UChar_t> hltm       (myReader,"hltsinglemu");
  TTreeReaderValue<UChar_t> hltp120    (myReader,"hltphoton120");
  TTreeReaderValue<double>  pswgt_ph120    (myReader,"pswgt_ph120");
  TTreeReaderValue<UChar_t> hltp165    (myReader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175    (myReader,"hltphoton175");
  
  string pfhtname ;
  if(not isMC)
    pfhtname = "hltPFHT800";
  else
    pfhtname = "hltPFHT800";

  TTreeReaderValue<UChar_t> hltpPFHT800    (myReader,pfhtname.c_str());

  TTreeReaderValue<UChar_t> fhbhe  (myReader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (myReader,"flaghbheiso");

  string cscfilter = "flagglobaltighthalo";
  if(isMC) cscfilter = "flagcsctight";
   
  TTreeReaderValue<UChar_t> fcsc   (myReader,cscfilter.c_str());
  TTreeReaderValue<UChar_t> fcsct  (myReader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (myReader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (myReader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (myReader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (myReader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (myReader,"flagbadchpf");

  TTreeReaderValue<unsigned int> njets    (myReader,"njets");
  TTreeReaderValue<unsigned int> ntaus    (myReader,"ntaus");
  TTreeReaderValue<unsigned int> ntausraw    (myReader,"ntausraw");
  TTreeReaderValue<unsigned int> nincjets (myReader,"njetsinc");
  TTreeReaderValue<unsigned int> nbjets   (myReader,"nbjetslowpt");
  TTreeReaderValue<unsigned int> nbjetshigh   (myReader,"nbjets");
  TTreeReaderValue<double> ht   (myReader,"ht");

  TTreeReaderValue<vector<double> > jeteta  (myReader,"combinejeteta");
  TTreeReaderValue<vector<double> > jetpt   (myReader,"combinejetpt");
  TTreeReaderValue<vector<double> > jetphi  (myReader,"combinejetphi");
  TTreeReaderValue<vector<double> > jetbtag (myReader,"combinejetbtag");
  TTreeReaderValue<vector<double> > jetm    (myReader,"combinejetm");
  TTreeReaderValue<vector<double> > jetQGL  (myReader,"combinejetQGL");
  TTreeReaderValue<vector<double> > chfrac  (myReader,"combinejetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac  (myReader,"combinejetNHfrac");
  TTreeReaderValue<vector<double> > emfrac  (myReader,"combinejetEMfrac");
  TTreeReaderValue<vector<double> > jetpuid (myReader,"combinejetPUID");
  TTreeReaderValue<vector<double> > pFlav   (myReader,"combinejetPFlav");

  // AK8 jet
  TTreeReaderValue<vector<double> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<double> > boostedJetQGL   (myReader,"boostedJetQGL");
  TTreeReaderValue<vector<double> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<double> > boostedJetphi   (myReader,"boostedJetphi");
  TTreeReaderValue<vector<double> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<double> > prunedJetm      (myReader,"prunedJetm_v2");
  TTreeReaderValue<vector<double> > prunedJetpt     (myReader,"prunedJetpt_v2");
  TTreeReaderValue<vector<double> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<double> > boostedJettau1  (myReader,"boostedJettau1");
  //  TTreeReaderValue<vector<double> > boostedJetpt    (myReader,"boostedPuppiJetpt");
  //  TTreeReaderValue<vector<double> > boostedJetQGL   (myReader,"boostedPuppiJetQGL");
  //  TTreeReaderValue<vector<double> > boostedJeteta   (myReader,"boostedPuppiJeteta");
  //  TTreeReaderValue<vector<double> > boostedJetphi   (myReader,"boostedPuppiJetphi");
  //  TTreeReaderValue<vector<double> > boostedJetm     (myReader,"boostedPuppiJetm");
  //  TTreeReaderValue<vector<double> > prunedJetm   (myReader,"softDropPuppiJetm");
  //  TTreeReaderValue<vector<double> > prunedJetpt   (myReader,"softDropPuppiJetpt");
  //  TTreeReaderValue<vector<double> > boostedJettau2  (myReader,"boostedPuppiJettau2");
  //  TTreeReaderValue<vector<double> > boostedJettau1  (myReader,"boostedPuppiJettau1");
  TTreeReaderValue<double > hadBosoneta  (myReader,"wzeta_h");
  TTreeReaderValue<double > hadBosonphi  (myReader,"wzphi_h");
  TTreeReaderValue<double > hadBosonpt   (myReader,"wzpt_h");
  TTreeReaderValue<double > hadBosonm    (myReader,"wzmass_h");

  // met systematics
  string metSuffix = "";
  if(sysName == "muUp")
    metSuffix = "MuEnUp";
  else if(sysName == "muDown" or sysName == "muDw")
    metSuffix = "MuEnDown";
  else if(sysName == "elUp")
    metSuffix = "ElEnUp";
  else if(sysName == "elDown" or sysName == "elDw")
    metSuffix = "ElEnDown";
  else if(sysName == "phoUp")
    metSuffix = "PhoEnUp";
  else if(sysName == "phoDown" or sysName == "phoDw")
    metSuffix = "PhoEnDown";
  else if(sysName == "tauUp")
    metSuffix = "TauEnUp";
  else if(sysName == "tauDown" or sysName == "tauDw")
    metSuffix = "TauEnDown";
 else if(sysName == "jesUp")
    metSuffix = "JetEnUp";
  else if(sysName == "jesDown" or sysName == "jesDw")
    metSuffix = "JetEnDown";
 else if(sysName == "jerUp")
    metSuffix = "JetResUp";
  else if(sysName == "jerDown" or sysName == "jerDw")
    metSuffix = "JetResDown";
 else if(sysName == "uncUp")
    metSuffix = "UncEnUp";
  else if(sysName == "uncDown" or sysName == "uncDw")
    metSuffix = "UncEnDown";

  TTreeReaderValue<double> met         (myReader,("t1pfmet"+metSuffix).c_str());
  TTreeReaderValue<double> metOriginal (myReader,"t1pfmet");
  TTreeReaderValue<double> metphi      (myReader,"t1pfmetphi");
  TTreeReaderValue<double> mmet        (myReader,"t1mumet");
  TTreeReaderValue<double> mmetphi     (myReader,"t1mumetphi");
  TTreeReaderValue<double> emet        (myReader,"t1elmet");
  TTreeReaderValue<double> emetphi     (myReader,"t1elmetphi"); 
  TTreeReaderValue<double> pmet        (myReader,"t1phmet");
  TTreeReaderValue<double> pmetphi     (myReader,"t1phmetphi");
  TTreeReaderValue<double> metpf       (myReader,"pfmet");
  TTreeReaderValue<double> metcalo     (myReader,"calomet");
  // dphi
  TTreeReaderValue<double> jmmdphi (myReader,"incjetmumetdphimin4");
  TTreeReaderValue<double> jemdphi (myReader,"incjetelmetdphimin4");
  TTreeReaderValue<double> jpmdphi (myReader,"incjetphmetdphimin4");

  TTreeReaderValue<int>    mu1pid (myReader,"mu1pid");
  TTreeReaderValue<int>    mu2pid (myReader,"mu2pid");
  TTreeReaderValue<int>    mu1id  (myReader,"mu1id");
  TTreeReaderValue<int>    mu2id  (myReader,"mu2id");
  TTreeReaderValue<double> mu1pt  (myReader,"mu1pt");
  TTreeReaderValue<double> mu2pt  (myReader,"mu2pt");
  TTreeReaderValue<double> mu1eta (myReader,"mu1eta");
  TTreeReaderValue<double> mu2eta (myReader,"mu2eta");
  TTreeReaderValue<double> mu1phi (myReader,"mu1phi");
  TTreeReaderValue<double> mu2phi (myReader,"mu2phi");

  TTreeReaderValue<int>    el1pid (myReader,"el1pid");
  TTreeReaderValue<int>    el2pid (myReader,"el2pid");
  TTreeReaderValue<int>    el1id  (myReader,"el1id");
  TTreeReaderValue<int>    el2id  (myReader,"el2id");
  TTreeReaderValue<double> el1pt  (myReader,"el1pt");
  TTreeReaderValue<double> el2pt  (myReader,"el2pt");
  TTreeReaderValue<double> el1eta (myReader,"el1eta");
  TTreeReaderValue<double> el2eta (myReader,"el2eta");
  TTreeReaderValue<double> el1phi (myReader,"el1phi");
  TTreeReaderValue<double> el2phi (myReader,"el2phi");
  
  TTreeReaderValue<int>   phidm  (myReader,"phidm");
  TTreeReaderValue<double> phpt  (myReader,"phpt");
  TTreeReaderValue<double> pheta (myReader,"pheta");
  TTreeReaderValue<double> phphi (myReader,"phphi");
  
  TTreeReaderValue<double> wmt   (myReader,"wmt");
  TTreeReaderValue<double> wemt   (myReader,"wemt");
  TTreeReaderValue<double> wzpt   (myReader,"wzpt");
  TTreeReaderValue<double> wzpt_h (myReader,"wzpt_h");
  TTreeReaderValue<double> wzeta  (myReader,"wzeta");
  TTreeReaderValue<double> zmass  (myReader,"zmass");
  TTreeReaderValue<double> zeemass  (myReader,"zeemass");
  TTreeReaderValue<double> zmmpt  (myReader,"zpt");
  TTreeReaderValue<double> zeept  (myReader,"zeept");
  TTreeReaderValue<double> zeeeta (myReader,"zeeeta");
  TTreeReaderValue<double> zeephi (myReader,"zeephi");
  TTreeReaderValue<double> zmmeta (myReader,"zeta");
  TTreeReaderValue<double> zmmphi (myReader,"zphi");

  TTreeReaderValue<double> dmpt (myReader,"dmpt");

  // other trick to handle the fact that this info is actually only stored for top/s-top samples
  string topptname;
  string atopptname;
  if(reWeightTopPt and isMC){
    topptname  = "toppt";
    atopptname = "atoppt";
  }
  else{
    topptname  = "wgt";
    atopptname = "wgt";
  }

  TTreeReaderValue<double> toppt  (myReader,topptname.c_str());
  TTreeReaderValue<double> atoppt (myReader,atopptname.c_str());

  //  ofstream dump("dumpGAM.txt");
  // loop on events
  while(myReader.Next()){

    //ICHEP dataset
    //    if(not isMC and *run > 274443) continue;
    //    if(not isMC and *run > 274240) continue;
    //    if(not isMC and *run > 275125) continue;
    if(not isMC and *run > 276811) continue;
    //if(sample == Sample::sig and not isMC and (category == Category::VBF or category == Category::twojet) and (*run <  275657 or *run > 276283)) continue;
    // Run E
    //    if(not isMC and *run < 276831) continue;
    //    if(not isMC and *run > 277420) continue;
    // Run F
    //    if(not isMC and *run < 277981) continue;
    //    if(not isMC and *run > 278808) continue;
    //Run G
    //if(not isMC and *run < 278820) continue;
    //if(not isMC and *run > 279116) continue;

    // check trigger depending on the sample
    Double_t hlt   = 0.0;
    Double_t hltw  = 1.0;

    if (sample == Sample::sig || sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu || sample == Sample::qcd)// single and double muon
      hlt = *hltm90+*hltm100+*hltm110+*hltm120+*hltmwm170;
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel) // single and double electron
      hlt = *hlte+*hltenoiso;      
    else if (sample == Sample::qcdgam || sample == Sample::gam){ // single photon
      hlt = *hltp165+*hltp175+*hltpPFHT800;
      hltw = 1;
      //if(*hltp120 and hlt == 0){
      // hltw = *pswgt_ph120;
      // hlt  = *hltp165+*hltp175+*hltpPFHT800+*hltp120;
      //}
    }
    // Trigger Selection
    if (hlt  == 0) continue; // trigger
    
    // MET Filters --> apply on both data and monte-carlo
    if(not isMC and (*fhbhe == 0 || *fhbiso == 0 || *feeb == 0 || *fetp == 0 || *fvtx == 0 || *fcsc == 0 || *fcsct == 0 || *fbadmu == 0 || *fbadch == 0)) continue;

    // check dphi jet-met
    Double_t jmdphi = 0.0;    
    if (sample == Sample::sig || sample == Sample::wmn || sample == Sample::zmm || sample == Sample::topmu || sample == Sample::qcd) jmdphi = fabs(*jmmdphi);
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel) jmdphi = fabs(*jemdphi);
    else if (sample == Sample::qcdgam || sample == Sample::gam) jmdphi = fabs(*jpmdphi);

    //set met
    Double_t pfmet = 0.0;
    Double_t pfmetphi = 0.0;
    if (sample == Sample::sig || sample == Sample::qcd) {pfmet = *mmet; pfmetphi = *mmetphi;}
    else if (sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu){ pfmet = *mmet; pfmetphi = *mmetphi;}
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel){ pfmet = *emet; pfmetphi = *emetphi;}
    else if (sample == Sample::qcdgam || sample == Sample::gam)  { pfmet = *pmet; pfmetphi = *pmetphi;}

    // noise cleaner
    if(fabs(*met-*metcalo)/pfmet > 0.5) continue;

    // propagate met systeamtics on the recoil
    if(metSuffix != "") pfmet += (*met-*metOriginal);

    // set lepton info
    Int_t    id1   = 0;
    Int_t    id2   = 0;
    Double_t pt1   = 0.0;
    Double_t pt2   = 0.0;
    Double_t eta1  = 0.0;
    Double_t eta2  = 0.0;
    Double_t phi1  = 0.0;
    Double_t phi2  = 0.0;
    int pid1  = 0;
    int pid2  = 0;

    if (sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu) {
      id1  = *mu1id;
      id2  = *mu2id;
      pt1  = *mu1pt;
      pt2  = *mu2pt;
      pid1 = *mu1pid;
      pid2 = *mu2pid;
      eta1 = *mu1eta;
      eta2 = *mu2eta;
      phi1 = *el1phi;
      phi2 = *el2phi;
    }
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel) {
      id1  = *el1id;
      id2  = *el2id;
      pt1  = *el1pt;
      pt2  = *el2pt;
      eta1 = *el1eta;
      eta2 = *el2eta;
      phi1 = *el1phi;
      phi2 = *el2phi;
      pid1 = *el1pid;
      pid2 = *el2pid;
    }
    else if (sample == Sample::qcdgam || sample == Sample::gam) {
      id1  = *phidm;
      id2  = 1.0;
      pt1  = *phpt;
      if(applyPhotonScale and not isMC){	
	pt1 += pt1*photonScaleUnc;
      }
      eta1 = *pheta;
    }
    
    // set zpt in case of Zsamples
    Double_t bosonPt  = 0.0;
    Double_t bosonPhi = 0.0;
    Double_t bosonEta = 0.0;
    if (sample == Sample::zmm){
      bosonPt  = *zmmpt; // di-muon CR
      bosonEta = *zmmeta;
      bosonPhi = *zmmphi;
    }
    else if (sample == Sample::zee){
      bosonPt = *zeept; // di-electron CR      
      bosonEta = *zeeeta;
      bosonPhi = *zeephi;
    }
    else if (sample == Sample::qcdgam or sample == Sample::gam){
      if(applyPhotonScale and not isMC)
	bosonPt = *phpt*(1+photonScaleUnc); // gamma+jets
      else
	bosonPt = *phpt;
      bosonEta  = *pheta;
      bosonPhi  = *phphi;
    }
    else if (sample == Sample::sig || sample == Sample::qcd){
      bosonPt = pfmet; // missing energy in case of signal region for Z->nunu
      bosonEta = 0;
      bosonPhi = pfmetphi;
    }
    else if (sample == Sample::wen or sample == Sample::zee){ // single muon or single ele
      TLorentzVector lep4V, met4V;
      lep4V.SetPtEtaPhiM(pt1,eta1,phi1,0.);
      met4V.SetPtEtaPhiM(*met,0.,*metphi,0.);
      bosonPt = (lep4V+met4V).Pt();
      bosonEta = 0;
      bosonPhi = (lep4V+met4V).Phi();
    }

    if (*ntausraw != 0) continue;
    // B-veto, not for top control sample
    if (*nbjets > 0 and sample != Sample::topmu and sample != Sample::topel) continue; 

    // control regions with two leptons --> opposite charge
    if (sample == Sample::zmm && *mu1pid == *mu2pid) continue;
    if (sample == Sample::zee && *el1pid == *el2pid) continue;

    // control regions with two leptons --> one should be tight
    if ((sample == Sample::zmm || sample == Sample::zee)){
      if(sample == Sample::zmm){	
	if(not ((pt1 > 20 and id1 == 1) or (pt2 > 20 and id2 == 1))) continue;
      }
      else if(sample == Sample::zee){
	if(not ((pt1 > 40 and id1 == 1) or (pt2 > 40 and id2 == 1))) continue;
      }
    }
    
    // number of central jets
    if (category != Category::VBF and category != Category::twojet and *njets < 1) continue; 
    else if((category == Category::VBF or category == Category::twojet) and *nincjets < 2) continue;

    // control regions wit one lepton --> tight requirement 
    if ((sample == Sample::wen || sample == Sample::wmn) && id1 != 1) continue;
    if (sample == Sample::wen and *wemt > 160) continue;
    if (sample == Sample::wmn and *wmt > 160) continue;
    if (sample == Sample::wmn and (category == Category::VBF or category == Category::twojet)){
      if(pfMetVBFLower < 200 and *wmt > 100) continue;
      else if(pfMetVBFLower > 200 and *wmt > 160) continue;
    }
    // photon control sample
    if ((sample == Sample::qcdgam || sample == Sample::gam) && pt1 < photonPt) continue;
    if ((sample == Sample::qcdgam || sample == Sample::gam) && fabs(*pheta) > 1.4442) continue;    

    // Wenu kill QCD
    if (sample == Sample::wen && *met < 50.) continue;
    
    // n-bjets cut for unboosted categories
    if ((sample == Sample::topmu || sample == Sample::topel) && (category != Category::monoV and category != Category::boosted and category != Category::prunedMass and category != Category::tau2tau1)  && *nbjets < nBjets) continue;
    if ( sample == Sample::topmu || sample == Sample::topel){ // select only events with one lepton
      // at least one lepton in the plateau region
      if(pt1 <=0 or id1 != 1) continue;
      if(abs(pid1) == 13 && pt1 < 20. ) continue;
      if(abs(pid1) == 11 && pt1 < 40. ) continue;
      // met cut
      if(sample == Sample::topel && *met < 50.) continue;
      // veto di-lepton events
      if(pt2 > 0) continue;
    }
    
    // met selection
    if(category == Category::VBF or category == Category::twojet){
      if (pfmet < pfMetVBFLower) continue;
      if (pfmet > pfMetVBFUpper) continue;
    }
    else  if(category == Category::monojet or category == Category::inclusive){
      if (pfmet < pfMetMonoJLower) continue;
      if (pfmet > pfMetMonoJUpper) continue;
    }
    else{
      if(pfmet < pfMetMonoVLower) continue;
      if(pfmet > pfMetMonoVUpper) continue;
    }

    //apply charge cut to separate W+ from W- in case
    if(sample == Sample::wmn or sample == Sample::wen){ // in case charge is required to be 1 skip events with negative leptons, viceversa
      if(vBosonCharge == 1 and pid1 > 0) continue;
      else if(vBosonCharge == -1 and pid1 < 0) continue;
    }

    // selection on jet --> split them between forward and central
    vector<TLorentzVector> centralJets;
    vector<TLorentzVector> forwardJets;
    int leadingCentralJetPos = -1;
    for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
      TLorentzVector vect;
      if(fabs(jeteta->at(ijet)) > 2.5 and jetpt->at(ijet) > 30){
	vect.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
	forwardJets.push_back(vect);
      }
      else if(fabs(jeteta->at(ijet)) < 2.5 and jetpt->at(ijet) > 30){
      	vect.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
      	centralJets.push_back(vect);
      	if(leadingCentralJetPos == -1)
      	  leadingCentralJetPos = ijet;
      }
    }
   
    if(category != Category::VBF and category != Category::twojet and leadingCentralJetPos < 0)  continue;
    if(category != Category::VBF and category != Category::twojet and leadingCentralJetPos != 0) continue; // asking leading jet to be central for non VBF categories

    // scale factor for leptons
    TH2* sflhist = NULL;
    TH2* sfthist = NULL;
    
    if (sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu) {
	sflhist = msfloose;
	sfthist = msftight;
    }
    if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel) {
	sflhist = msfloose;
	sfthist = esftight;
    }
       
    Double_t sfwgt = 1.0;

    if(isMC && (sample == Sample::zee or sample == Sample::wen)){
      if(pt1 > 0. and *nvtx < 15)
	sfwgt *= trackingefficiency_electron_lowpu->GetBinContent(trackingefficiency_electron_lowpu->FindBin(fabs(eta1),min(pt1,trackingefficiency_electron_lowpu->GetYaxis()->GetBinLowEdge(trackingefficiency_electron_lowpu->GetNbinsY()+1)-1)));
      if(pt2 > 0. and *nvtx < 15)
	sfwgt *= trackingefficiency_electron_lowpu->GetBinContent(trackingefficiency_electron_lowpu->FindBin(fabs(eta2),min(pt2,trackingefficiency_electron_lowpu->GetYaxis()->GetBinLowEdge(trackingefficiency_electron_lowpu->GetNbinsY()+1)-1)));
      if(pt1 > 0. and *nvtx >= 15)
	sfwgt *= trackingefficiency_electron_highpu->GetBinContent(trackingefficiency_electron_highpu->FindBin(fabs(eta1),min(pt1,trackingefficiency_electron_highpu->GetYaxis()->GetBinLowEdge(trackingefficiency_electron_highpu->GetNbinsY()+1)-1)));
      if(pt2 > 0. and *nvtx >= 15)
	sfwgt *= trackingefficiency_electron_highpu->GetBinContent(trackingefficiency_electron_highpu->FindBin(fabs(eta2),min(pt2,trackingefficiency_electron_highpu->GetYaxis()->GetBinLowEdge(trackingefficiency_electron_highpu->GetNbinsY()+1)-1))); 
    }
   
    if (isMC && sflhist && sfthist) {
      if (pt1 > 0.) {
	if (id1 == 1) sfwgt *= sfthist->GetBinContent(sfthist->FindBin(fabs(eta1),min(pt1,sfthist->GetYaxis()->GetBinLowEdge(sfthist->GetNbinsY()+1)-1))); 
	else sfwgt *= sflhist->GetBinContent(sflhist->FindBin(fabs(eta1),min(pt1,sflhist->GetYaxis()->GetBinLowEdge(sflhist->GetNbinsY()+1)-1))); 
      }
      if (pt2 > 0.) {
	if (id2 == 1) sfwgt *= sfthist->GetBinContent(sfthist->FindBin(fabs(eta2),min(pt2,sfthist->GetYaxis()->GetBinLowEdge(sfthist->GetNbinsY()+1)-1))); 
	else sfwgt *= sflhist->GetBinContent(sflhist->FindBin(fabs(eta2),min(pt2,sflhist->GetYaxis()->GetBinLowEdge(sflhist->GetNbinsY()+1)-1))); 
      }
    }

    // reco-muon scale factor    
    if(isMC && (sample == Sample::zmm or sample == Sample::wmn or sample == Sample::topmu)){
      double trackwgt = 1;
      if(pt1 > 0. and *nvtx < 15)
	trackwgt *= trackingefficiency_muon_lowpu->GetBinContent(trackingefficiency_muon_lowpu->FindBin(fabs(eta1),min(pt1,trackingefficiency_muon_lowpu->GetYaxis()->GetBinLowEdge(trackingefficiency_muon_lowpu->GetNbinsY()+1)-1)));
      if(pt2 > 0. and *nvtx < 15)
	trackwgt *= trackingefficiency_muon_lowpu->GetBinContent(trackingefficiency_muon_lowpu->FindBin(fabs(eta2),min(pt2,trackingefficiency_muon_lowpu->GetYaxis()->GetBinLowEdge(trackingefficiency_muon_lowpu->GetNbinsY()+1)-1)));
      if(pt1 > 0. and *nvtx >= 15)
	trackwgt *= trackingefficiency_muon_highpu->GetBinContent(trackingefficiency_muon_highpu->FindBin(fabs(eta1),min(pt1,trackingefficiency_muon_highpu->GetYaxis()->GetBinLowEdge(trackingefficiency_muon_highpu->GetNbinsY()+1)-1)));
      if(pt2 > 0. and *nvtx >= 15)
	trackwgt *= trackingefficiency_muon_highpu->GetBinContent(trackingefficiency_muon_highpu->FindBin(fabs(eta2),min(pt2,trackingefficiency_muon_highpu->GetYaxis()->GetBinLowEdge(trackingefficiency_muon_highpu->GetNbinsY()+1)-1)));      
      sfwgt *= trackwgt;
    }

    // trigger scale factor for electrons
    if (isMC && triggerelhist && triggerelhist_ht && ( sample == Sample::zee || sample == Sample::topel || sample == Sample::wen)) {
      if (pt1 > 40. && id1 == 1 and id2 == 1)
	sfwgt *= 1;
      else if(id1 == 1 and id2 != 1 and pt1 >= 125)
	sfwgt *= triggerelhist_ht->GetBinContent(triggerelhist_ht->FindBin(fabs(eta1),min(pt1,triggerelhist_ht->GetYaxis()->GetBinLowEdge(triggerelhist_ht->GetNbinsY()+1)-1)));
      else if(id1 == 1 and id2 != 1 and pt1 < 125)
	sfwgt *= triggerelhist->GetBinContent(triggerelhist->FindBin(fabs(eta1),min(pt1,triggerelhist->GetYaxis()->GetBinLowEdge(triggerelhist->GetNbinsY()+1)-1)));            
      else if(id2 == 1 and id1 != 1 and pt2 >= 125)
	sfwgt *= triggerelhist_ht->GetBinContent(triggerelhist_ht->FindBin(fabs(eta2),min(pt2,triggerelhist_ht->GetYaxis()->GetBinLowEdge(triggerelhist_ht->GetNbinsY()+1)-1)));
      else if(id2 == 1 and id1 != 1 and pt2 < 125)
	sfwgt *= triggerelhist->GetBinContent(triggerelhist->FindBin(fabs(eta2),min(pt2,triggerelhist->GetYaxis()->GetBinLowEdge(triggerelhist->GetNbinsY()+1)-1)));
    }
    

    // photon id scale factor
    if (isMC && psfmedium && sample == Sample::gam) {
      if (pt1 > 0. && id1 == 1) {
	sfwgt *= psfmedium->GetBinContent(psfmedium->FindBin(fabs(eta1),min(pt1,psfmedium->GetYaxis()->GetBinLowEdge(psfmedium->GetNbinsY()+1)-1)));
      }
    }
    
    // photon purity
    if (!isMC && purhist && sample == Sample::qcdgam) {
      if(pt1 > purhist->GetXaxis()->GetBinLowEdge(1))
	sfwgt *= (1.0 - purhist->GetBinContent(purhist->FindBin(min(pt1,purhist->GetXaxis()->GetBinLowEdge(purhist->GetNbinsX()+1)-1), fabs(eta1))));
      else
	sfwgt *= (1.0 - purhist->GetBinContent(purhist->FindBin(max(pt1,purhist->GetXaxis()->GetBinLowEdge(1)),fabs(eta1))));		  
    }
    
    // met trigger scale factor
    if (isMC && triggermet && (sample == Sample::sig || sample == Sample::wmn || sample == Sample::zmm || sample == Sample::topmu || sample == Sample::qcd)) {
      sfwgt *= triggermet_graph->Eval(min(pfmet,triggermet_graph->GetXaxis()->GetXmax()));
    }
    
    // photon trigger scale factor
    if(isMC && triggerphoton_graph && triggerphoton_graph_jetHT && (sample == Sample::qcdgam || sample == Sample::gam)){ // linear interpolation between graph points            
      if(*pmet < recoilThresholdTrigger)	
	sfwgt *= triggerphoton_graph->Eval(min(*pmet,triggerphoton_graph->GetXaxis()->GetXmax()));
      else
	sfwgt *= triggerphoton_graph_jetHT->Eval(min(*pmet,triggerphoton_graph->GetXaxis()->GetXmax()));
    }
    
    // B10-4 -tag weight
    double btagw = *wgtbtag;    
    if(sample == Sample::topmu or sample == Sample::topel)
      btagw = 0.9;
    else
      btagw = 0.99;
    
    //V-tagging scale factor --> only for mono-V
    if(isMC && category == Category::monoV && isWJet)
      sfwgt *= getVtaggingScaleFactor(tau2tau1,sysName);
    
    //Gen level info --> NLO re-weight    
    Double_t kwgt = 1.0;    
    double genpt = *wzpt;
    for (size_t i = 0; i < khists.size(); i++) {
      if (khists[i]) {
	if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
	kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
      }
    }
    
    // Higgs pT uncertainty
    Double_t hwgt = 1.0;
    if(isHiggsInvisible and hhist and isMC){
      if(*dmpt < hhist->GetBinLowEdge(1))
	*dmpt = hhist->GetBinLowEdge(1)+1;
      else if(*dmpt > hhist->GetBinLowEdge(hhist->GetNbinsX()+1))
	*dmpt =hhist->GetBinLowEdge(hhist->GetNbinsX()+1)-1;
      hwgt *= hhist->GetBinContent(hhist->FindBin(*dmpt));
    }

    // post fit re-weight
    Double_t pfwgt = 1.0;
    if(postFitWeight){
      double met = pfmet;
      if(met < postFitWeight->GetBinLowEdge(1))
	met = postFitWeight->GetBinLowEdge(1)+1;
      else if(met > postFitWeight->GetBinLowEdge(postFitWeight->GetNbinsX()+1))
	met = postFitWeight->GetBinLowEdge(postFitWeight->GetNbinsX()+1)-1;
      pfwgt =  postFitWeight->GetBinContent(postFitWeight->FindBin(met));
    }

    // data based re-weight for gamma+jets
    Double_t gamwgt = 1.0;
    if(sample == Sample::gam and isMC and category == Category::monojet and reweightPhton){
      gamwgt = gamPtWeight->GetBinContent(gamPtWeight->GetXaxis()->FindBin(*phpt));
    }

    // Top quark pt re-weight
    Double_t topptwgt = 1.0;
    if(reWeightTopPt)
      topptwgt = reweightTopQuarkPt(*toppt,*atoppt);

    // ggZH re-weight in case of a non null pointer                                                                                                                             
    Double_t ggZHwgt = 1.0;
    if(isHiggsInvisible and ggZHhist and isMC){
      int binX = ggZHhist->GetXaxis()->FindBin(*dmpt);
      int binY = 0;
      if(*wzpt_h != 0)
        binY = ggZHhist->GetYaxis()->FindBin(*wzpt_h);
      else
        binY = ggZHhist->GetYaxis()->FindBin(*wzpt);

      if(binX == 0) binX = 1;
      if(binY == 0) binY = 1;
      if(binX == ggZHhist->GetNbinsX()+1) binX = ggZHhist->GetNbinsX();
      if(binY == ggZHhist->GetNbinsY()+1) binY = ggZHhist->GetNbinsY();
      ggZHwgt = ggZHhist->GetBinContent(binX,binY);
      if(ggZHwgt == 0)
        ggZHwgt = 1;
    }
    
    /// Start specific analysis selections
    // inclusive mono-jet analysis 
    if(category == Category::inclusive){ 
      if (centralJets.size() == 0) continue;
      if (fabs(jeteta->at(leadingCentralJetPos)) > 2.5) continue;
      if (chfrac->at(leadingCentralJetPos) < 0.1) continue;   // jet id
      if (nhfrac->at(leadingCentralJetPos) > 0.8) continue;   // jet id
      if (jetpt->at(leadingCentralJetPos)  < leadingJetPtCut) continue;  // jet1 > 100 GeV
      if (sample != Sample::qcd and jmdphi < 0.5) continue; // deltaPhi cut
      else if (sample == Sample::qcd and jmdphi > 0.5) continue; // deltaPhi cut
    }
    else{

      // boosted category and monojet
      bool goodMonoJet = false;
      bool goodMonoV   = false;
      
      if(category == Category::monojet){ // mono jet + V-jet veto

	if (centralJets.size() == 0) continue;
	if (fabs(jeteta->at(leadingCentralJetPos)) > 2.5) continue;
	if (chfrac->at(leadingCentralJetPos) < 0.1) continue;   // jet id                                                                                                   
	if (nhfrac->at(leadingCentralJetPos) > 0.8) continue;   // jet id                                                                                                   
	if (jetpt->at(leadingCentralJetPos)  < leadingJetPtCut) continue;  // jet1 > 100 GeV                                                                                          
	if (sample != Sample::qcd and jmdphi < 0.5) continue; 
	else if(sample == Sample::qcd and jmdphi > 0.5) continue;

	if(boostedJetpt->size()  == 0) goodMonoJet = true;
	if(boostedJetpt->size() > 0){ // in case one boosted jet	  
	  if(fabs(boostedJeteta->at(0)) > jetEtaAK8) 
	    goodMonoJet = true;	  

	  if(boostedJetpt->at(0) < ptJetMinAK8) // check pT
	    goodMonoJet = true;	 	  
	  else{ 
	    // pruned mass selection
	    if(prunedJetm->at(0) < prunedMassMin  or prunedJetm->at(0) > prunedMassMax)
	      goodMonoJet= true;
	    // tau2tau1 selection
	    //	    if((boostedJettau2->at(0)/boostedJettau1->at(0)+0.063*log(pow(prunedJetm->at(0),2)/prunedJetpt->at(0))) > tau2tau1)
	    if(boostedJettau2->at(0)/boostedJettau1->at(0) > tau2tau1)
	      goodMonoJet= true;
	  }
	}
      	
	if(not goodMonoJet) continue;
	//	if(pfmet > 800 and not isMC) dump<<"run "<<*run<<" lumi "<<*lumisection<<" event "<<*event<<" met "<<pfmet<<" \n";
      }

      else if(category == Category::monoV or category == Category::boosted or category == Category::prunedMass or category == Category::tau2tau1){
	
	if (centralJets.size() == 0) continue;
	if (boostedJetpt->size() == 0) continue;
	if (boostedJetpt->at(0) < ptJetMinAK8) continue;
	if (fabs(boostedJeteta->at(0)) > jetEtaAK8) continue;
	if (fabs(jeteta->at(leadingCentralJetPos)) > 2.5) continue;	
	//after match apply jetid on leading ak4
	if (chfrac->at(leadingCentralJetPos) < 0.1) continue;   // jet id                                                                                                     
	if (nhfrac->at(leadingCentralJetPos) > 0.8) continue;   // jet id                                                                                                  
	if (jetpt->at(leadingCentralJetPos)  < leadingJetPtCut) continue;  // jet1 > 100 GeV                                                                                           
	if (sample != Sample::qcd and jmdphi < 0.5) continue;
	else if(sample == Sample::qcd and jmdphi > 0.5) continue;

	TLorentzVector jetak4, jetak8;
	jetak4.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jetak8.SetPtEtaPhiM(boostedJetpt->at(0),boostedJeteta->at(0),boostedJetphi->at(0),boostedJetm->at(0));	

	// no overlap between b-jet and v-jet
	if (sample == Sample::topel || sample == Sample::topmu){ 
	  int nbjets = 0;
	  for(size_t ijet = 0 ; ijet < jetbtag->size(); ijet++){
	    jetak4.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
	    if(jetak4.DeltaR(jetak8) < 0.8) continue;
	    if(jetbtag->at(ijet) > 0.8){
	      nbjets++;
	    }
	  }
	  if(nbjets < 1) continue;
	}

	// split among resonant and non resonant wrt gen level
	if(resonantSelection != 0 and isMC){
	  TLorentzVector Wboson4V;
	  Wboson4V.SetPtEtaPhiM(*hadBosonpt,*hadBosoneta,*hadBosonphi,*hadBosonm);
	  if(jetak8.DeltaR(Wboson4V) > 0.4 and resonantSelection == 1)
	    continue;
	  else if(jetak8.DeltaR(Wboson4V) < 0.4 and resonantSelection == 2)
	    continue;
	}

	// category 2 means HP mono-V
	if(category == Category::monoV and (prunedJetm->at(0) > prunedMassMin and prunedJetm->at(0) < prunedMassMax) and boostedJettau2->at(0)/boostedJettau1->at(0) < tau2tau1)
	//	if(category == Category::monoV and (prunedJetm->at(0) > prunedMassMin and prunedJetm->at(0) < prunedMassMax) and (boostedJettau2->at(0)/boostedJettau1->at(0)+0.063*log(pow(prunedJetm->at(0),2)/prunedJetpt->at(0))) < tau2tau1)
	  goodMonoV   = true;
	// category 3 means LP mono-V
	else if(category == Category::prunedMass and (prunedJetm->at(0) > prunedMassMin and prunedJetm->at(0) < prunedMassMax) and 
		(boostedJettau2->at(0)/boostedJettau1->at(0) > tau2tau1 and boostedJettau2->at(0)/boostedJettau1->at(0) < tau2tau1LP))
	  goodMonoV   = true;
	// apply no pruned mass cut --> show full shapes
	else if(category == Category::boosted and (prunedJetm->at(0) > 0 and prunedJetm->at(0) < 200))
	  goodMonoV   = true;
	// apply only n-subjettiness
	else if(category == Category::tau2tau1  and boostedJettau2->at(0)/boostedJettau1->at(0) < tau2tau1)
	  goodMonoV   = true;
	
	if(not goodMonoV) continue;	
	//	if(pfmet > 750. and not isMC) dump<<"run "<<*run<<" lumi "<<*lumisection<<" event "<<*event<<" met "<<pfmet<<" \n";

      }
      else if(category == Category::VBF){
	if(centralJets.size()+forwardJets.size() < 2) continue;
	if(fabs(jeteta->at(0)) > 4.7 or fabs(jeteta->at(1)) > 4.7) continue;
	if(jetpt->at(0) < leadingJetPtCutVBF) continue;
	if(jetpt->at(1) < trailingJetPtCutVBF) continue;

	if (sample != Sample::qcd and jmdphi < jetmetdphiVBF) continue;
	else if(sample == Sample::qcd and jmdphi > jetmetdphiVBF) continue;

	if(fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) < 0.1) continue;
	if(fabs(jeteta->at(0)) < 2.5 and nhfrac->at(0) > 0.8) continue;
	if(jeteta->at(0)*jeteta->at(1) > 0 ) continue;
	if(fabs(jeteta->at(0)-jeteta->at(1)) < detajj) continue;
	if(fabs(jeteta->at(0)) >= 3.0 and nhfrac->at(0) > 0.96) continue;
	if(sample == Sample::qcd and fabs(jeteta->at(1)) >= 2.5 and fabs(jeteta->at(1)) < 3.0 and chfrac->at(1) < 0.05) continue;
	TLorentzVector jet1 ;
	TLorentzVector jet2 ;
	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	if((jet1+jet2).M() < mjj) continue;
      }
      else if(category == Category::twojet){
	if(centralJets.size()+forwardJets.size() < 2) continue;
	if(fabs(jeteta->at(0)) > 4.7 or fabs(jeteta->at(1)) > 4.7) continue;
	if(jetpt->at(0) < leadingJetPtCutVBF) continue;
	if(jetpt->at(1) < trailingJetPtCutVBF) continue;
	if(fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) < 0.1) continue;
	if(fabs(jeteta->at(0)) < 2.5 and nhfrac->at(0) > 0.8) continue;
	if(fabs(jeteta->at(0)) >= 3.0 and nhfrac->at(0) > 0.96) continue;
	if(sample == Sample::qcd and fabs(jeteta->at(1)) >= 2.5 and fabs(jeteta->at(1)) < 3.0 and chfrac->at(1) < 0.05) continue;
	if (sample != Sample::qcd and jmdphi < 0.5) continue;
	else if(sample == Sample::qcd and jmdphi > 0.5) continue;
      }
    }

    // fill 1D histogram
    double fillvar = 0;
    // fill the histograms --> with the right observable
    for(auto hist : hist1D){
      TString name(hist->GetName());      
      if(name.Contains("nvtx")) // number of vertices
	fillvar = *nvtx;
      else if(name.Contains("zmass")) // Z-mass peak to mumu
	fillvar = *zmass;
      else if(name.Contains("zeemass")) // Z-mass peak to ee
	fillvar = *zeemass;
      else if(name.Contains("t1pfmetphi")) // met-phi with type-I correction	
	fillvar = *metphi;      
      else if(name.Contains("t1pfmet")) // met with type-I correction
	fillvar = *met;      
      //Delta phi met leptons: mu1, el1 or ph
      else if(name.Contains("dphi_t1met_lep1"))
	fillvar = deltaPhi(*metphi,phi1);
      //Delta phi met leptons: mu2, el2
      else if(name.Contains("dphi_t1met_lep2"))
	fillvar = deltaPhi(*metphi,phi2);
      //Delta phi met leptons and boson
      else if(name.Contains("dphi_t1met_boson"))
	fillvar = deltaPhi(*metphi,bosonPhi);
      // Lepton Pt and Eta
      else if(name.Contains("elp1pt") and pid1 > 0) // leading e+ pt
	fillvar = pt1;      
      else if(name.Contains("elm1pt") and pid1 < 0) // leading e- pt
	fillvar = pt1;      
      else if(name.Contains("el1pt")) // e pt
	fillvar = pt1;      
      else if(name.Contains("mup1pt") and pid1 > 0) // mu+ pt
	fillvar = pt1;
      else if(name.Contains("mum1pt") and pid1 < 0) // mu- pt
	fillvar = pt1;
      else if(name.Contains("mu1pt")) // mu pt
	fillvar = pt1;
      else if(name.Contains("mup2pt") and pid2 > 0)
	fillvar = pt2;
      else if(name.Contains("mum2pt") and pid2 < 0)
	fillvar = pt2;
      else if(name.Contains("mu2pt"))
	fillvar = pt2;
      else if(name.Contains("elp2pt") and pid2 > 0)
	fillvar = pt2;
      else if(name.Contains("elm2pt") and pid2 < 0)
	fillvar = pt2;
      else if(name.Contains("el2pt"))
	fillvar = pt2;
      else if(name.Contains("mup1eta")){
	if(pid1 > 0)
	  fillvar = eta1;      
	else 
	  fillvar = -99.;
      }
      else if(name.Contains("mum1eta")){
	if(pid1 < 0)
          fillvar = eta1;
	else
          fillvar = -99.;
      }
      else if(name.Contains("mu1eta"))
	fillvar = eta1;      
      else if(name.Contains("mup2eta")){
	if(pid2 > 0)
	  fillvar = eta2;
	else
	  fillvar = -99.;
      } 
      else if(name.Contains("mum2eta")){
	if(pid2 < 0)
	  fillvar = eta2;
	else
	  fillvar = -99.;
      } 
      else if(name.Contains("mu2eta"))
	fillvar = eta2;
      else if(name.Contains("elp1eta")){
	if(pid1 > 0)
	  fillvar = eta1;      
	else 
	  fillvar = -99.;
      }
      else if(name.Contains("elm1eta")){
	if(pid1 < 0)	  
	  fillvar = eta1;      
	else
	  fillvar = -99.;
      }
      else if(name.Contains("el1eta"))
	fillvar = eta1;
      else if(name.Contains("elp2eta")){
	if(pid1 > 0)
	  fillvar = eta2;
	else
	  fillvar = -99.;
      }
      else if(name.Contains("elm2eta")){
	if(pid2 > 0)
	  fillvar = eta2;
	else
	  fillvar = -99.;
      }
      else if(name.Contains("el2eta"))
	fillvar = eta2;
      // W-boson properties
      else if(name.Contains("wmt"))
	fillvar = *wmt;
      else if(name.Contains("wemt"))
	fillvar = *wemt;
      // jet fractions
      else if(name.Contains("nhfrac_v2")){
	if(fabs(jeteta->at(0)) >= 3)
	  fillvar = nhfrac->at(0);
	else 
	  fillvar = -1;
      }
      else if(name.Contains("nhfrac_v3")){
	if(fabs(jeteta->at(0)) >= 2.5 and fabs(jeteta->at(0)) <= 3)
	  fillvar = nhfrac->at(0);
	else 
	  fillvar = -1;
      }
      else if(name.Contains("nhfracj2_v2") and jetpt->size() >= 2){
	if(fabs(jeteta->at(1)) >= 3)
	  fillvar = nhfrac->at(1);
	else 
	  fillvar = -1;
      }
      else if(name.Contains("nhfracj2_v3") and jetpt->size() >= 2){
	if(fabs(jeteta->at(1)) >= 2.5 and fabs(jeteta->at(1)) <= 3)
	  fillvar = nhfrac->at(1);
	else 
	  fillvar = -1;
      }
      else if(name.Contains("emfrac_v2")){
	if(fabs(jeteta->at(0)) >= 3)
	   fillvar = emfrac->at(0);
	else 
	  fillvar = -1;
      }

      else if(name.Contains("chfrac_v2")){
	if(fabs(jeteta->at(0)) >= 2.5 and fabs(jeteta->at(0)) <= 3)
	  fillvar = chfrac->at(0);
      }
      else if(name.Contains("chfracj2_v2") and jetpt->size() >= 2){
	if(fabs(jeteta->at(1)) >= 2.5 and fabs(jeteta->at(1)) <= 3)
	  fillvar = chfrac->at(1);
      }
      else if(name.Contains("chfracj2") and jetpt->size() >= 2){
	if(category == Category::VBF or category == Category::twojet)
	  fillvar = chfrac->at(1);
	else
	  fillvar = chfrac->at(leadingCentralJetPos);
      }      
      else if(name.Contains("nhfrac")){
	if(category == Category::VBF or category == Category::twojet)
	  fillvar = nhfrac->at(0);
	else
	  fillvar = nhfrac->at(leadingCentralJetPos);
      }
      else if(name.Contains("chfrac")){
	if(category == Category::VBF or category == Category::twojet)
	  fillvar = chfrac->at(0);
	else
	  fillvar = chfrac->at(leadingCentralJetPos);
      }
      else if(name.Contains("emfrac")){
	if(category == Category::VBF or category == Category::twojet)
	  fillvar = emfrac->at(0);
	else
	  fillvar = emfrac->at(leadingCentralJetPos);
      } 
      // AK8 jet kinematics
      else if(name.Contains("boostedjetpt") and boostedJetpt->size() > 0)
	  fillvar = boostedJetpt->at(0);
      else if(name.Contains("boostedjeteta") and boostedJeteta->size() > 0)
	  fillvar = boostedJeteta->at(0);
      // Jet met Dphi used for selection
      else if(name.Contains("jetmetdphi"))
	fillvar = jmdphi;
      // calo met
      else if(name.Contains("calomet"))
	fillvar = fabs(*metcalo-*met)/pfmet;
      else if(name.Contains("met"))
	fillvar = pfmet;             
      else if(name.Contains("jetpt2") and jetpt->size() >= 2)
	fillvar = jetpt->at(1);
      else if(name.Contains("jeteta2") and jetpt->size() >= 2)
	fillvar = jeteta->at(1);
      else if(name.Contains("jetpt")){
	if(category == Category::VBF or category == Category::twojet)
	  fillvar = jetpt->at(0);
	else
	  fillvar = jetpt->at(leadingCentralJetPos);
      }
      else if(name.Contains("jetpuid")){
	if(fabs(jeteta->at(0)) >= 3)
	  fillvar = jetpuid->at(0);
	else
	  fillvar = -1;
      }
	  
      else if(name.Contains("jeteta")){
	if(category == Category::VBF or category == Category::twojet)
	  fillvar = jeteta->at(0);
	else
	  fillvar = jeteta->at(leadingCentralJetPos);
      }
      else if(name.Contains("detajj") and jetpt->size() >= 2)
	fillvar = fabs(jeteta->at(0)-jeteta->at(1));
      else if(name.Contains("mjj") and jetpt->size() >= 2){
	TLorentzVector jet1 ;
        TLorentzVector jet2 ;
        jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	fillvar = (jet1+jet2).M();
      }
      else if(name.Contains("mT")){
	float deltaphi = deltaPhi(jetphi->at(0),pfmetphi);
	fillvar    = sqrt(2*jetpt->at(0)*pfmet*(1-cos(deltaphi)));
      }   
      else if(name.Contains("ncjet"))
	fillvar = *njets;      
      else if(name.Contains("njet"))
	fillvar = *nincjets;      
      else if(name.Contains("nbjet_hpt_loose")){
	int nbjet = 0;
	for(size_t iJet = 0; iJet < jetbtag->size(); iJet++){
	  if(jetbtag->at(iJet) > btagCSVLoose and jetpt->at(iJet) > 30)
	    nbjet++;
	}
	fillvar = nbjet;
      }      
      else if(name.Contains("nbjet_hpt")){
	int nbjet = 0;
	for(size_t iJet = 0; iJet < jetbtag->size(); iJet++){
	  if(jetbtag->at(iJet) > btagCSVMedium and jetpt->at(iJet) > 30)
	    nbjet++;
	}
	fillvar = nbjet;
      }      
      else if(name.Contains("nbjet"))
	fillvar = *nbjets;
      else if(name.Contains("bosonpt"))
	fillvar = bosonPt;    	
      else if(name.Contains("jetQGL"))
	fillvar = jetQGL->at(0);
      else if(name.Contains("jet2QGL") and jetQGL->size() > 1)
	fillvar = jetQGL->at(1);
      else if(name.Contains("jetPFlav") and pFlav->size() > 0)
	fillvar = fabs(pFlav->at(0));
      else if(name.Contains("jetQGL_AK8")){
	if(boostedJetpt->size() > 0)
	  fillvar = boostedJetQGL->at(0);
	else
	  fillvar = -1.;
      }
      
      // substructure
      else if(name.Contains("mpruned")){
	if( prunedJetm->size() > 0 and boostedJetpt->at(0) > ptJetMinAK8 )
	  fillvar = prunedJetm->at(0);	
	else fillvar = -1.;
      }
      else if(name.Contains("HT"))
	fillvar = *ht;      
      else if(name.Contains("tau2tau1") and boostedJettau1->size() > 0 and boostedJettau2->size() > 0 and boostedJetpt->at(0) > ptJetMinAK8)
	fillvar = boostedJettau2->at(0)/boostedJettau1->at(0);
      //	fillvar = boostedJettau2->at(0)/boostedJettau1->at(0)+0.063*log(pow(prunedJetm->at(0),2)/prunedJetpt->at(0));
	
      // b-tagging
      else if(name.Contains("btagCSV_max")){
	float btagMax = -10.;
	for(size_t iBjet = 0; iBjet < jetbtag->size(); iBjet++){
	  if(jeteta->at(iBjet) > 2.5) continue;
	  if(jetbtag->at(iBjet) >= btagMax)
	    btagMax = jetbtag->at(iBjet);
	}
	if(fabs(btagMax) != 10.)
	  fillvar = btagMax;
	else
	  fillvar = 0.;
      }
      else if(name.Contains("btagCSV_min")){
	float btagMin = 10.;
	for(size_t iBjet = 0; iBjet < jetbtag->size(); iBjet++){
	  if(jeteta->at(iBjet) >= 2.5) continue;
	  if(jetbtag->at(iBjet) < btagMin and jetbtag->at(iBjet) > 0)
	    btagMin = jetbtag->at(iBjet);
	}
	if(fabs(btagMin) != 10.)
	  fillvar = btagMin;
	else
	  fillvar = 0.;
      }
      else if(name.Contains("btagCSV")){
	if(jetbtag->size() > 0)
	  fillvar = jetbtag->at(0);
	else
	  fillvar = 0;
      }      
      //
      else if(name.Contains("minDphiJJ")){
	if(jetphi->size() <= 1)
	  fillvar = hist->GetXaxis()->GetBinCenter(1);
	else{
	  float minDphi = 2*TMath::Pi();
	  bool  isfound = false;
	  for(size_t ijet = 0 ; ijet < jetphi->size(); ijet++){
	    for(size_t jjet = ijet+1 ; jjet < jetphi->size(); jjet++){
	      if(jetpt->at(ijet) < 30 or jetpt->at(jjet) < 30) continue;
	      float deltaphi = deltaPhi(jetphi->at(ijet),jetphi->at(jjet));
	      if(deltaphi > 0 and deltaphi < minDphi){
		minDphi = deltaphi;
		isfound = true;
	      }
	    }
	  }
	  if(isfound)
	    fillvar = minDphi; 
	  else
	    fillvar = hist->GetXaxis()->GetBinCenter(1);
	}
      }

      else if(name.Contains("minDphiJ1J")){
	float minDphi = 2*TMath::Pi();
	bool  isfound = false;
	if(jetphi->size() <= 1)
	  fillvar = hist->GetXaxis()->GetBinCenter(1);
	else{
	  for(size_t jjet = 1 ; jjet < jetphi->size(); jjet++){
	    float deltaphi = deltaPhi(jetphi->at(0),jetphi->at(jjet));
	    if(jetpt->at(0) < 30 or jetpt->at(jjet) < 30) continue;
	    if(deltaphi > 0 and deltaphi < minDphi){
	      minDphi = deltaphi;
	      isfound = true;
	    }
	  }
	  if(isfound)
	    fillvar = minDphi; 
	  else
	    fillvar = hist->GetXaxis()->GetBinCenter(1);
	}
      }
      else if(name.Contains("dphiJJ")){
	if(jetphi->size() < 2){				
	  cout<<"problem with VBF category "<<endl;
	  cout<<jetpt->at(0)<<" "<<jetpt->at(1)<<endl;
	  fillvar = hist->GetXaxis()->GetBinCenter(1);
	}
	else 
	  fillvar = deltaPhi(jetphi->at(0),jetphi->at(1));
      }
      
      // overflow bin
      if (fillvar >= hist->GetBinLowEdge(hist->GetNbinsX())+hist->GetBinWidth(hist->GetNbinsX())) 
	fillvar = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());

      // total event weight
      double evtwgt  = 1.0;
      Double_t puwgt = 0.;
      if (isMC and not reweightNVTX){
	if (*putrue <= 100)
	  puwgt = puhist->GetBinContent(puhist->FindBin(*putrue));
	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*gamwgt*hwgt*pfwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*gamwgt*hwgt*pfwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)
      }
      else if (isMC and reweightNVTX){
	if (*nvtx <= 60) 
	  puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*topptwgt*sfwgt*kwgt*hwgt*ggZHwgt*pfwgt*gamwgt/(*wgtsum);
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*topptwgt*sfwgt*kwgt*hwgt*ggZHwgt*pfwgt*gamwgt/(*wgtsum);
      }
      
      if (!isMC && sample == Sample::qcdgam) 
	evtwgt = sfwgt*pfwgt*hltw;
      else if (!isMC)
	evtwgt = hltw;      
      hist->Fill(fillvar, evtwgt);
    }

    //fill 2D histo
    double fillvarX = 0;
    double fillvarY = 0;

    for(auto hist: hist2D){
      TString name(hist->GetName());
      if(name.Contains("met_jetpt")){ 
	fillvarX = pfmet;
	fillvarY = jetpt->at(0);
      }
      else if(name.Contains("met_bosonpt")){
	fillvarX = pfmet;
	fillvarY = bosonPt;
      }
      else if(name.Contains("met_mT")){	
	fillvarX = pfmet;
	float deltaPhi = fabs(jetphi->at(0)-pfmetphi);
	if(deltaPhi > TMath::Pi())
	  deltaPhi = 2*TMath::Pi() - deltaPhi;
	fillvarY = sqrt(2*jetpt->at(0)*pfmet*(1-cos(deltaPhi)));
      }
      else if(name.Contains("met_mpruned") and boostedJetpt->size() > 0 and boostedJetpt->at(0) > ptJetMinAK8 ){
	fillvarX = pfmet;
	if(prunedJetm->size() > 0)
	  fillvarY = prunedJetm->at(0);	
      }
      else if(name.Contains("met_tau2tau1") and boostedJetpt->size() > 0 and boostedJetpt->at(0) > ptJetMinAK8){
	fillvarX = pfmet;
	if(boostedJettau2->size() > 0 and boostedJettau1->size() >0)
	  fillvarY = boostedJettau2->at(0)/boostedJettau1->at(0);		
      }
      else if(name.Contains("mpruned_tau2tau1") and boostedJetpt->size() > 0 and boostedJetpt->at(0) > ptJetMinAK8){
	if(prunedJetm->size() > 0)
          fillvarX = prunedJetm->at(0);
	if(boostedJettau2->size() > 0 and boostedJettau1->size() >0)
	  fillvarY = boostedJettau2->at(0)/boostedJettau1->at(0);			
      }

      else if(name.Contains("met_ncjet")){
	fillvarX = pfmet;
	fillvarY = *njets;
      }
      else if(name.Contains("met_njet")){
	fillvarX = pfmet;
	fillvarY = *nincjets;
      }

      else if(name.Contains("mT_ncjet")){
	float deltaPhi = fabs(jetphi->at(0)-pfmetphi);
        if(deltaPhi > TMath::Pi())
          deltaPhi = 2*TMath::Pi() - deltaPhi;
	fillvarX = sqrt(2*jetpt->at(0)*pfmet*(1-cos(deltaPhi)));
	fillvarY = *njets;
      }
      else if(name.Contains("mT_njet")){
	float deltaPhi = fabs(jetphi->at(0)-pfmetphi);
        if(deltaPhi > TMath::Pi())
          deltaPhi = 2*TMath::Pi() - deltaPhi;
	fillvarX = sqrt(2*jetpt->at(0)*pfmet*(1-cos(deltaPhi)));
	fillvarY = *nincjets;
      }

      else if(name.Contains("met_ht")){
        fillvarX = pfmet;
        fillvarY = *ht;
      }
      else if(name.Contains("mT_ht")){
	float deltaPhi = fabs(jetphi->at(0)-pfmetphi);
        if(deltaPhi > TMath::Pi())
          deltaPhi = 2*TMath::Pi() - deltaPhi;
	fillvarX = sqrt(2*jetpt->at(0)*pfmet*(1-cos(deltaPhi)));
        fillvarY = *ht;
      }
      else if(name.Contains("met_QGL")){
	fillvarX = pfmet;
	fillvarY = jetQGL->at(0);
      }
      else if(name.Contains("met_minDphiJJ")){
	fillvarX = pfmet;
	if(jetphi->size() <= 1)
	  fillvarY = hist->GetYaxis()->GetBinCenter(1);	    
	else{
	  float minDphi = 2*TMath::Pi();
	  bool isfound = false;
	  for(size_t ijet = 0 ; ijet < jetphi->size(); ijet++){
	    for(size_t jjet = ijet+1 ; jjet < jetphi->size(); jjet++){
	      if(jetpt->at(ijet) < 30 or jetpt->at(jjet) < 30) continue;
	      float deltaPhi = fabs(jetphi->at(ijet)-jetphi->at(jjet));
	      if(deltaPhi > TMath::Pi())
		deltaPhi = 2*TMath::Pi() - deltaPhi;
	      if(deltaPhi > 0 and deltaPhi < minDphi){
		minDphi = deltaPhi;
	        isfound = true;
	      }
	    }
	  }
	  if(isfound)
	    fillvarY = minDphi;
	  else
	    fillvarY = hist->GetYaxis()->GetBinCenter(1);
	}
      }
      else if(name.Contains("met_minDphiJ1J")){
	fillvarX = pfmet;
	if(jetphi->size() <= 1)
	  fillvarY = hist->GetYaxis()->GetBinCenter(1);	    
	else{
	  float minDphi = 2*TMath::Pi();
	  bool isfound = false;
	  for(size_t jjet = 1 ; jjet < jetphi->size(); jjet++){
	    if(jetpt->at(0) < 30 or jetpt->at(jjet) < 30) continue;
	    float deltaPhi = fabs(jetphi->at(0)-jetphi->at(jjet));
	    if(deltaPhi > TMath::Pi())
	      deltaPhi = 2*TMath::Pi() - deltaPhi;
	    if(deltaPhi > 0 and deltaPhi < minDphi){
	      minDphi = deltaPhi;
	      isfound = true;
	    }
	  }	
	  if(isfound)
	    fillvarY = minDphi;
	  else
	    fillvarY = hist->GetYaxis()->GetBinCenter(1);
	}
      }
      else if(name.Contains("met_dphiJJ")){
	fillvarX = pfmet;
	if(jetphi->size() < 2)
	  fillvarY = hist->GetYaxis()->GetBinCenter(1);
	else {
	  fillvarY = fabs(jetphi->at(0)-jetphi->at(1));
	  if(fillvarY > TMath::Pi())
	    fillvarY = 2*TMath::Pi()-fillvarY;
	}
      }
      
      // overflow bin
      if (fillvarX >= hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX())+hist->GetXaxis()->GetBinWidth(hist->GetNbinsX())) 
	fillvarX = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());
      if (fillvarY >= hist->GetYaxis()->GetBinLowEdge(hist->GetNbinsY())+hist->GetYaxis()->GetBinWidth(hist->GetNbinsY())) 
	fillvarY = hist->GetYaxis()->GetBinCenter(hist->GetNbinsY());

      // total event weight
      double evtwgt  = 1.0;
      Double_t puwgt = 0.;

      if (isMC and not reweightNVTX){

	if (*putrue <= 100)
          puwgt = puhist->GetBinContent(puhist->FindBin(*putrue));
        if(XSEC != -1)
          evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*topptwgt*sfwgt*kwgt*hwgt*ggZHwgt*pfwgt*gamwgt/(*wgtsum);
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*topptwgt*sfwgt*kwgt*hwgt*ggZHwgt*pfwgt*gamwgt/(*wgtsum);
      }
      else if (isMC and reweightNVTX){
	if (*nvtx <= 60) 
	  puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*hwgt*gamwgt/(*wgtsum);
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*hwgt*gamwgt/(*wgtsum);
      }
      if (!isMC && sample == Sample::qcdgam) 
	evtwgt = sfwgt*hltw;
      else if (!isMC)
	evtwgt = hltw;
      
	hist->Fill(fillvarX,fillvarY,evtwgt);            
     }
  }

  //  dump.close();

  sffile_eleTight.Close();
  sffile_eleVeto.Close();
  sffile_muTight.Close();
  sffile_muLoose.Close();
  sffile_phoMedium.Close();
  purityfile_photon.Close();
  triggerfile_SinglEle.Close();
  triggerfile_MET->Close();
  triggerfile_SinglePhoton.Close();
  trackingefficiency_muon.Close();
  trackingefficiency_electron.Close();
  gamPtFile.Close();
}

#endif
