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

// some basic cut values --> Monojet category
const float leadingJetPtCut  = 100.;
const float pfMetMonoJUpper  = 8000.;
const float pfMetMonoJLower  = 250.;
const float btagCSVLoose     = 0.460;
const float btagCSVMedium    = 0.800;
const float numberOfVtxCorrection = 17;
// some basic cut values --> MonoV category
const float tau2tau1        = 0.6;
const float tau2tau1LP      = 0.75;
const float prunedMassMin   = 65.;
const float prunedMassMax   = 105.;
const float ptJetMinAK8     = 250.;
const float jetEtaAK8       = 2.4;
const float pfMetMonoVLower = 250.;
const float pfMetMonoVUpper = 8000.;
// some basic cut values --> VBF category
const float leadingJetPtCutVBF  = 80.;
const float trailingJetPtCutVBF = 40.;
const float detajj          = 4.0;
const float detajjrelaxed   = 1.0;
const float mjj             = 1300;
const float mjjrelaxed      = 200.;
const float jetmetdphiVBF   = 0.5;
const float pfMetVBFLower   = 250.;
const float pfMetVBFUpper   = 8000.;
const float dphijj          = 1.5;
const float dphijjrelaxed   = 1.3;
const bool  removeVBF       = false;
// Additional selections
const float photonPt        = 175;
const int   vBosonCharge    = 0;
const int   nBjets          = 1; // for top-tagged region
const int   njetsMin        = 1;
const int   njetsMax        = 100;
// Re-weight and smoothing
static bool  reweightNVTX    = false;
/// photon scale
const bool  applyPhotonScale = true;
const float photonScaleUnc   = -0.0125;
static bool doSmoothing      = false;
// trigger and object corrections
const float recoilThresholdTrigger = 350; // for photon trigger application
const bool  isSummer16      = true;
const bool  useSingleMuon   = true;
const bool  useMoriondSetup = true;
const bool  usePOGScaleFactors = true;
// other general options
const bool  runOnlyData     = false;
// k-factors
const bool  applyEWKVKfactor = true;

// k-factors
//string kfactorFile       = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root";
string kfactorFile       = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_24bins.root";
string kfactorFileUnc    = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactors_uncertainties.root";
string kfactorFileUNLOPS = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kfactor_gamma_unlops.root";
string kFactorTheoristFile_zvv = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors_theorist_v4/vvj.root";
string kFactorTheoristFile_wln = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors_theorist_v4/evj.root";
string kFactorTheoristFile_zll = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors_theorist_v4/eej.root";
string kFactorTheoristFile_gam = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors_theorist_v4/aj.root";
string kFactorFile_wjetewk = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kFactor_WToLNu_pT_Mjj.root";
string kFactorFile_zjetewk = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/kFactor_ZToNuNu_pT_Mjj.root";

/// basic trees
string baseInputTreePath = "/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_6_06_2017/";
//string baseInputTreePath = "/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_1_02_2017";

VectorSorter jetSorter;

float deltaPhi (float phi1, float phi2){
  if(fabs(phi1-phi2) < TMath::Pi())
    return fabs(phi1-phi2);
  else
    return 2*TMath::Pi()-fabs(phi1-phi2);
}

// top-pt re-weight
double reweightTopQuarkPt (double topQuarkPt, double atopQuarkPt){

  double weight_top  = -1.;
  double weight_atop = -1.;
  if(topQuarkPt  != 0)
    weight_top = exp(0.0615-0.0005*topQuarkPt);
  if(atopQuarkPt != 0)
    weight_atop = exp(0.0615-0.0005*atopQuarkPt);
  
  if(weight_top != -1. and weight_atop != -1.)
    return sqrt(weight_top*weight_atop);
  else return 1.;
}

// v-tagging scale factor
double getVtaggingScaleFactor(const double & tau2tau1, const string & sysName){

  double sfwgt = 1;

  if(tau2tau1 == 0.40 or tau2tau1 == 0.45){
    if(sysName == "VtagUp")
      sfwgt *= (0.92+0.109);
    else if(sysName == "VtagDown")
      sfwgt *= (0.92-0.109);
    else
      sfwgt *= 0.92;
  }
  else if(tau2tau1 == 0.6){
    if(sysName == "VtagUp")
      sfwgt *= (0.92+0.109);
    else if(sysName == "VtagDown")
      sfwgt *= (0.92-0.109);
    else
      sfwgt *= 0.92;
  }
  
  return sfwgt;
}

// main function
void makehist4(TTree* tree,            /*input tree*/ 
	       vector<TH1*> hist1D,    /* set of 1D histogram */ 
	       vector<TH2*> hist2D,    /* set of 2D histogram */ 
	       const bool &   isMC,    // data or MC
	       const Sample   & sample,   // sample to select
	       const Category & category, // category for event
	       const bool &   isWJet,     // is a W-jet sample like ttbar/di-boson
	       const double   & scale,    // overall scale
	       const double   & lumi,     // luminosity        
	       vector<TH1*>   khists,     // NLO k-factors vs boson pT
	       const string   & sysName,  // Sys variation	
	       const bool     & reWeightTopPt      = false,
	       bool             reweightNVTX       = true,	       
	       const int      & resonantSelection  = 0,
	       const bool     & isHiggsInvisible   = false, // reject VBF events
	       const bool     & applyPostFitWeight = false,
	       vector<TH2*>   kVEWKhists = {}, // k-factor for V-EWK process in the VBF case
	       const float    & XSEC = -1.,// fix the cross section from extern	       
	       TH1*  hhist     = NULL,
	       TH1*  higgsNNLO = NULL,
	       TH2*  ggZHhist  = NULL,
	       const string & leptonPID = ""
	       ) {

  if(not tree){
    cout<<" empty tree --> skip process "<<endl;
    return;
  }

  if(runOnlyData and isMC) return;
  // fix by hand for 2016 data
  if((category == Category::monojet or category == Category::monoV) and reweightNVTX == false) reweightNVTX = true;

  // Pileup Weights
  TFile* pufile = NULL;
  TH1*   puhist = NULL;
  if(useMoriondSetup){
    if(reweightNVTX){
      if(not isSummer16){
	pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_35p9fb.root");    
	puhist = (TH1*) pufile->Get("puhist");
      }
      else{
	pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt_36.40_summer16.root");    
	puhist = (TH1*) pufile->Get("puhist");
      }
    }
  }
  else{
    if(reweightNVTX){
      pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_12p9fb.root");    
      puhist = (TH1*) pufile->Get("puhist");
    }
  }
  
  // Lepton ID scale factors
  TFile* sffile_eleTight = NULL;
  TFile* sffile_eleVeto  = NULL;
  TFile* sffile_muTight  = NULL;
  TFile* sffile_muLoose  = NULL;
  
  if(useMoriondSetup){
    if(not isSummer16){
      sffile_eleTight = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/scaleFactor_electron_tightid.root");
      sffile_eleVeto  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/scaleFactor_electron_vetoid.root");
    }
    else{
      if(not usePOGScaleFactors){
	sffile_eleTight = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/scaleFactor_electron_tightid_summer16.root");
	sffile_eleVeto  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/scaleFactor_electron_vetoid_summer16.root");
      }
      else{
	sffile_eleTight = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/scalefactors_80x_egpog_37ifb.root");
	sffile_eleVeto  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/scalefactors_80x_egpog_37ifb.root");	
      }
    }
    if(usePOGScaleFactors){
      sffile_muTight  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/muon_scalefactors_37ifb.root");
      sffile_muLoose  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/muon_scalefactors_37ifb.root");
    }
    else{
      sffile_muTight  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/scaleFactor_muon_tightid.root");
      sffile_muLoose  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/scaleFactor_muon_looseid.root");
    }
  }
  else{
    sffile_eleTight = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_ICHEP/scaleFactor_electron_tightid_12p9.root");
    sffile_eleVeto  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_ICHEP/scaleFactor_electron_vetoid_12p9.root");
    sffile_muTight  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_ICHEP/scaleFactor_muon_tightid_12p9.root");
    sffile_muLoose  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_ICHEP/scaleFactor_muon_looseid_12p9.root");
  }
  

  TH2*  msfloose_lowpu  = NULL;
  TH2*  msftight_lowpu  = NULL;
  TH2*  msfloose_id   = NULL;
  TH2*  msftight_id   = NULL;
  TH2*  msfloose_iso  = NULL;
  TH2*  msftight_iso  = NULL;
  TH2*  esfveto_lowpu   = NULL;
  TH2*  esftight_lowpu  = NULL;
  TH2*  msfloose_highpu = NULL;
  TH2*  msftight_highpu = NULL;
  TH2*  esfveto_highpu  = NULL;
  TH2*  esftight_highpu = NULL;
  TH2*  esfveto  = NULL;
  TH2*  esftight = NULL;
  
  if(useMoriondSetup){
    if(not isSummer16){
      msfloose_lowpu  = (TH2*) sffile_muLoose->Get("scaleFactor_muon_looseid_RooCMSShape_pu_0_17");
      msftight_lowpu  = (TH2*) sffile_muTight->Get("scaleFactor_muon_tightid_RooCMSShape_pu_0_17");
      esfveto_lowpu   = (TH2*) sffile_eleVeto->Get("scaleFactor_electron_vetoid_RooCMSShape_pu_0_17");
      esftight_lowpu  = (TH2*) sffile_eleTight->Get("scaleFactor_electron_tightid_RooCMSShape_pu_0_17");
      msfloose_highpu = (TH2*) sffile_muLoose->Get("scaleFactor_muon_looseid_RooCMSShape_pu_17_50");
      msftight_highpu = (TH2*) sffile_muTight->Get("scaleFactor_muon_tightid_RooCMSShape_pu_17_50");
      esfveto_highpu  = (TH2*) sffile_eleVeto->Get("scaleFactor_electron_vetoid_RooCMSShape_pu_17_50");
      esftight_highpu = (TH2*) sffile_eleTight->Get("scaleFactor_electron_tightid_RooCMSShape_pu_17_50");
    }
    else{
      if(not usePOGScaleFactors){
	msfloose_lowpu  = (TH2*) sffile_muLoose->Get("scaleFactor_muon_looseid_RooCMSShape_pu_0_17");
	msftight_lowpu  = (TH2*) sffile_muTight->Get("scaleFactor_muon_tightid_RooCMSShape_pu_0_17");
	msfloose_highpu = (TH2*) sffile_muLoose->Get("scaleFactor_muon_looseid_RooCMSShape_pu_17_50");
	msftight_highpu = (TH2*) sffile_muTight->Get("scaleFactor_muon_tightid_RooCMSShape_pu_17_50");
      }
      else{
	msfloose_id  = (TH2*) sffile_muLoose->Get("scalefactors_MuonLooseId_Muon");
	msfloose_iso = (TH2*) sffile_muLoose->Get("scalefactors_Iso_MuonLooseId");
	msftight_id  = (TH2*) sffile_muLoose->Get("scalefactors_TightId_Muon");
	msftight_iso = (TH2*) sffile_muLoose->Get("scalefactors_Iso_MuonTightId");
      }
      if(not usePOGScaleFactors){
	esfveto    = (TH2*) sffile_eleVeto->Get("scaleFactor_electron_vetoid_RooCMSShape_pu_0_100");
	esftight   = (TH2*) sffile_eleTight->Get("scaleFactor_electron_tightid_RooCMSShape_pu_0_100");
      }
      else{
	esfveto    = (TH2*) sffile_eleVeto->Get("scalefactors_Veto_Electron");
	esftight   = (TH2*) sffile_eleTight->Get("scalefactors_Tight_Electron");
      }
    }
  }
  else{
    msfloose_lowpu  = (TH2*) sffile_muLoose->Get("scaleFactor_muon_looseid_RooCMSShape");
    msftight_lowpu  = (TH2*) sffile_muTight->Get("scaleFactor_muon_tightid_RooCMSShape");
    esfveto_lowpu   = (TH2*) sffile_eleVeto->Get("scaleFactor_electron_vetoid_RooCMSShape");
    esftight_lowpu  = (TH2*) sffile_eleTight->Get("scaleFactor_electron_tightid_RooCMSShape");
    msfloose_highpu = (TH2*) sffile_muLoose->Get("scaleFactor_muon_looseid_RooCMSShape");
    msftight_highpu = (TH2*) sffile_muTight->Get("scaleFactor_muon_tightid_RooCMSShape");
    esfveto_highpu  = (TH2*) sffile_eleVeto->Get("scaleFactor_electron_vetoid_RooCMSShape");
    esftight_highpu = (TH2*) sffile_eleTight->Get("scaleFactor_electron_tightid_RooCMSShape");
  }

  // Photon ID scale factor                                                                                                                                                     
  TFile* sffile_phoMedium = NULL;
  TH2*  psfmedium_lowpu   = NULL;
  TH2*  psfmedium_highpu  = NULL;
  TH2*  psfmedium         = NULL;

  if(useMoriondSetup){
    if(not isSummer16){
      sffile_phoMedium = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF_2016/photonSF_Moriond/scaleFactor_photon_mediumid.root");
      psfmedium_lowpu  = (TH2*) sffile_phoMedium->Get("scaleFactor_photon_mediumid_RooCMSShape_pu_0_17");
      psfmedium_highpu = (TH2*) sffile_phoMedium->Get("scaleFactor_photon_mediumid_RooCMSShape_pu_17_50");
    }
    else{
      sffile_phoMedium = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF_2016/photonSF_Moriond/photonSF_egamma_80X.root");
      psfmedium = (TH2*) sffile_phoMedium->Get("EGamma_SF2D");
    }
  }
  else{
    sffile_phoMedium = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF_2016/photonSF_ICHEP/scaleFactor_photon_mediumid_12p9.root");
    psfmedium_lowpu  = (TH2*) sffile_phoMedium->Get("scaleFactor_photon_mediumid_RooCMSShape");
    psfmedium_highpu = (TH2*) sffile_phoMedium->Get("scaleFactor_photon_mediumid_RooCMSShape");
  }

  // Photon Purity                                                                                                                                                              
  TFile* purityfile_photon = NULL;
  TH2*   purhist           = NULL; 
  TGraphAsymmErrors*   purgraph  = NULL; 
  if(useMoriondSetup){
    purityfile_photon = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF_2016/photonSF_Moriond/finalPurity_weightedAverage_36fb.root");
    purgraph = (TGraphAsymmErrors*) purityfile_photon->Get("purity_weightedAverage_totUnc");
  }
  else{
    purityfile_photon = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF_2016/photonSF_ICHEP/PhotonSFandEffandPurity_Lumi12p9fb_28072016.root");  
    purhist = (TH2*) purityfile_photon->Get("purity");
  }
  
  // Lepton track efficiency
  TFile* trackingefficiency_muon        = NULL;
  TH2F*  trackingefficiency_muon_lowpu  = NULL;
  TH2F*  trackingefficiency_muon_highpu = NULL;
  TGraphAsymmErrors*  trackingefficiency_pog = NULL;
  if(useMoriondSetup){
    if(not usePOGScaleFactors){
      trackingefficiency_muon = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/trackEfficiency/trackEfficiency_Moriond/scaleFactor_muon_trackerid.root");
      trackingefficiency_muon_lowpu = (TH2F*) trackingefficiency_muon->Get("scaleFactor_muon_trackerid_RooCMSShape_pu_0_17");
      trackingefficiency_muon_highpu = (TH2F*) trackingefficiency_muon->Get("scaleFactor_muon_trackerid_RooCMSShape_pu_17_50");
    }
    else{
      trackingefficiency_muon = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/trackEfficiency/trackEfficiency_Moriond/Tracking_EfficienciesAndSF_BCDEFGH.root");
      trackingefficiency_pog = (TGraphAsymmErrors*) trackingefficiency_muon->Get("ratio_eff_vtx_dr030e030_corr");
    }
  }
  else{
    trackingefficiency_muon = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/trackEfficiency/trackEfficiency_ICHEP/scaleFactor_muon_trackerid.root");
    trackingefficiency_muon_lowpu  = (TH2F*) trackingefficiency_muon->Get("scaleFactor_muon_trackerid_Exp_pu_0_16");
    trackingefficiency_muon_highpu = (TH2F*) trackingefficiency_muon->Get("scaleFactor_muon_trackerid_Exp_pu_16_50");
  }
  
  TFile* trackingefficiencyFile_electron   = NULL;
  TH2F* trackingefficiency_electron_lowpu  = NULL;
  TH2F* trackingefficiency_electron_highpu = NULL;
  TH2F* trackingefficiency_electron        = NULL;
  if(useMoriondSetup){
    if(isSummer16){
      if(usePOGScaleFactors){
	trackingefficiencyFile_electron = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/leptonSF_Moriond/scalefactors_80x_egpog_37ifb.root");
	trackingefficiency_electron = (TH2F*) trackingefficiencyFile_electron->Get("scalefactors_Reco_Electron");
      }
      else{
	trackingefficiencyFile_electron = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/trackEfficiency/trackEfficiency_Moriond/scaleFactor_electron_recoelectronmatch_summer16.root");
	trackingefficiency_electron = (TH2F*) trackingefficiencyFile_electron->Get("scaleFactor_electron_recoelectronmatch_RooCMSShape_pu_0_100");	
      }
    }
    else{
      trackingefficiencyFile_electron = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/trackEfficiency/trackEfficiency_Moriond/scaleFactor_electron_recoelectronmatch.root");
      trackingefficiency_electron_lowpu  = (TH2F*) trackingefficiencyFile_electron->Get("scaleFactor_electron_recoelectronmatch_RooCMSShape_pu_0_17");
      trackingefficiency_electron_highpu = (TH2F*) trackingefficiencyFile_electron->Get("scaleFactor_electron_recoelectronmatch_RooCMSShape_pu_17_50");
    }
  }
  else{
    trackingefficiencyFile_electron = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/trackEfficiency/trackEfficiency_ICHEP/scaleFactor_electron_recoelectronmatch.root");  
    trackingefficiency_electron_lowpu  = (TH2F*) trackingefficiencyFile_electron->Get("scaleFactor_electron_recoelectronmatch_RooCMSShape_pu_0_16");
    trackingefficiency_electron_highpu = (TH2F*) trackingefficiencyFile_electron->Get("scaleFactor_electron_recoelectronmatch_RooCMSShape_pu_16_50");
  } 

  /////////////////////////////////////////
  // trigger files used for 2016                                                                                                                                                
  ////////////////////////////////////////
  TFile* triggerfile_SinglEle       = NULL;
  TFile* triggerfile_SinglEle_jetHT = NULL;
  TEfficiency* triggerel_eff        = NULL;
  TEfficiency* triggerel_eff_jetHT  = NULL;
  TH2* triggerelhist    = NULL;
  TH2* triggerelhist_ht = NULL;


  if(useMoriondSetup){
    if(category != Category::VBF and category != Category::VBFrelaxed){
      triggerfile_SinglEle       = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/Monojet/triggerEfficiency_DATA_SingleElectron.root");
      triggerfile_SinglEle_jetHT = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/Monojet/triggerEfficiency_DATA_SingleElectron.root");
      triggerel_eff        = (TEfficiency*) triggerfile_SinglEle->Get("trgeff_ele");
      triggerel_eff_jetHT  = (TEfficiency*) triggerfile_SinglEle_jetHT->Get("trgeff_ele");
      triggerelhist    = triggerel_eff->CreateHistogram();
      triggerelhist_ht = triggerel_eff_jetHT->CreateHistogram();
      triggerelhist->SetName("triggerelhist");
      triggerelhist_ht->SetName("triggerelhist_ht");
    }
    else{
      triggerfile_SinglEle       = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/eleTrig.root");
      triggerfile_SinglEle_jetHT = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/eleTrig.root");
      triggerelhist     = (TH2*) triggerfile_SinglEle->Get("hEffEtaPt");
      triggerelhist_ht  = (TH2*) triggerfile_SinglEle_jetHT->Get("hEffEtaPt");
    }
  }
  else{
    triggerfile_SinglEle       = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_ICHEP/Monojet/triggerEfficiency_DATA_SingleElectron_12p9fb.root");
    triggerfile_SinglEle_jetHT = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_ICHEP/Monojet/triggerEfficiency_DATA_SingleElectron_jetHT_12p9.root");
    triggerel_eff        = (TEfficiency*) triggerfile_SinglEle->Get("trgeff_ele");
    triggerel_eff_jetHT  = (TEfficiency*) triggerfile_SinglEle_jetHT->Get("efficiency");
  }


  
  // Met trigger efficiency
  TFile* triggerfile_MET     = NULL;
  TFile* triggerfile_MET_zmm = NULL;
  vector<TFile*> triggerfile_MET_binned_Wmn;
  vector<TFile*> triggerfile_MET_binned_Zmm;
  if(useMoriondSetup){
    if(category != Category::VBF and category != Category::twojet and category != Category::VBFrelaxed){ // monojet
      // Use Wmn or Wen for W+jets and SR
      if(useSingleMuon){
	triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/Monojet/metTriggerEfficiency_recoil_monojet.root");
      }
      else
	triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/Monojet/metTriggerEfficiency_ele_recoil_monojet.root");
      // Zmm measurement for Zmm CR
      triggerfile_MET_zmm = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/Monojet/metTriggerEfficiency_recoil_monojet_zmm.root");

    }
    else{
      if(useSingleMuon){
	triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_0.0_800.0.root"));
	triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_800.0_1200.0.root"));
	triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_1200.0_1700.0.root"));
	triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_1700.0_3000.0.root"));
	triggerfile_MET_binned_Zmm.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_zmm_0.0_800.0.root"));
	triggerfile_MET_binned_Zmm.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_zmm_800.0_1200.0.root"));
	triggerfile_MET_binned_Zmm.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_zmm_1200.0_1700.0.root"));
	triggerfile_MET_binned_Zmm.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_zmm_1700.0_3000.0.root"));
      }
      else{
	triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_ele_mjj_vbf_0.0_800.0.root"));
	triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_ele_mjj_vbf_800.0_1200.0.root"));
	triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_ele_mjj_vbf_1200.0_1700.0.root"));
	triggerfile_MET_binned_Wmn.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_ele_mjj_vbf_1700.0_3000.0.root"));
	triggerfile_MET_binned_Zmm.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_zmm_0.0_800.0.root"));
	triggerfile_MET_binned_Zmm.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_zmm_800.0_1200.0.root"));
	triggerfile_MET_binned_Zmm.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_zmm_1200.0_1700.0.root"));
	triggerfile_MET_binned_Zmm.push_back(TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/metTriggerEfficiency_mjj_vbf_zmm_1700.0_3000.0.root"));
      }
    }
  }
  else{
    if(category != Category::VBF and category != Category::twojet and category != Category::VBFrelaxed){ // monojet
      if(useSingleMuon)
	triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_ICHEP/Monojet/metTriggerEfficiency_12p9.root");
      else
	triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_ICHEP/Monojet/metTriggerEfficiency_12p9.root");
    }
    else{ // VBF
      if(useSingleMuon)
	triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_ICHEP/VBF/metTriggerEfficiency_muon.root");      
      else
	triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_ICHEP/VBF/metTriggerEfficiency_electron.root");
    }
  }
  
  // single turn on
  TEfficiency*       triggermet         = NULL;
  TEfficiency*       triggermet_zmm     = NULL;
  ////
  TGraphAsymmErrors* triggermet_graph = NULL;    
  TGraphAsymmErrors* triggermet_graph_zmm = NULL;    
  TGraphAsymmErrors* triggermet_graph_zvv_mc  = NULL;    
  TGraphAsymmErrors* triggermet_graph_wjet_mc = NULL;    
  TGraphAsymmErrors* triggermet_graph_wjet_sf = NULL;    

  if(triggerfile_MET != NULL){
    triggermet = (TEfficiency*) triggerfile_MET->Get("trig_eff");
    if(triggermet == 0 or triggermet == NULL)
      triggermet = (TEfficiency*) triggerfile_MET->Get("efficiency_vbf_loose");
    if(triggermet == 0 or triggermet == NULL)
      triggermet = (TEfficiency*) triggerfile_MET->Get("efficiency");
    triggermet_graph = triggermet->CreateGraph();
  }

  if(triggerfile_MET_zmm != NULL){
    triggermet_zmm = (TEfficiency*) triggerfile_MET_zmm->Get("trig_eff");
    if(triggermet_zmm == 0 or triggermet_zmm == NULL)
      triggermet_zmm = (TEfficiency*) triggerfile_MET_zmm->Get("efficiency_vbf_loose");
    if(triggermet_zmm == 0 or triggermet_zmm == NULL)
      triggermet_zmm = (TEfficiency*) triggerfile_MET_zmm->Get("efficiency");
    triggermet_graph_zmm = triggermet_zmm->CreateGraph();
  }
  
  vector<TF1*> triggermet_func_binned_Wmn;
  vector<TF1*> triggermet_func_binned_Zmm;
  if(triggerfile_MET_binned_Wmn.size() != 0){
    for(auto ifile : triggerfile_MET_binned_Wmn)
      triggermet_func_binned_Wmn.push_back((TF1*) ifile->Get("efficiency_func"));    
  }
  if(triggerfile_MET_binned_Zmm.size() != 0){
    for(auto ifile : triggerfile_MET_binned_Zmm)
      triggermet_func_binned_Zmm.push_back((TF1*) ifile->Get("efficiency_func"));    
  }

  // Photon trigger efficiency measured in jetHT
  TFile* triggerfile_SinglePhoton_jetHT = NULL;
  TFile* triggerfile_SinglePhoton       = NULL;
  TEfficiency* triggerphoton       = NULL;
  TEfficiency* triggerphoton_jetHT = NULL;
  if(useMoriondSetup){
    if(category != Category::VBF and category != Category::VBFrelaxed){
      triggerfile_SinglePhoton_jetHT = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/Monojet/photonTriggerEfficiency_jetHT.root");
      triggerfile_SinglePhoton       = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/Monojet/photonTriggerEfficiency_photon.root");
      triggerphoton = (TEfficiency*) triggerfile_SinglePhoton->Get("eff_photonpt");
      triggerphoton_jetHT = (TEfficiency*) triggerfile_SinglePhoton_jetHT->Get("eff_recover_photonpt");
    }
    else{
      triggerfile_SinglePhoton_jetHT = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/photonTriggerEfficiency_jetHT.root");
      triggerfile_SinglePhoton       = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/VBF/photonTriggerEfficiency_photon.root");
      triggerphoton = (TEfficiency*) triggerfile_SinglePhoton->Get("eff_vbf_recoil");
      triggerphoton_jetHT = (TEfficiency*) triggerfile_SinglePhoton_jetHT->Get("eff_recover_vbf_recoil");
    }
  }
  else{
    triggerfile_SinglePhoton_jetHT = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_ICHEP/VBF/photonTriggerEfficiency_photonpt120_jetHT.root");
    triggerfile_SinglePhoton       = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_ICHEP/VBF/photonTriggerEfficiency_photonpt120.root");
    triggerphoton       = (TEfficiency*) triggerfile_SinglePhoton->Get("eff_recoil");
    triggerphoton_jetHT = (TEfficiency*) triggerfile_SinglePhoton_jetHT->Get("eff_recover_recoil");
  }

  TGraphAsymmErrors* triggerphoton_graph       = triggerphoton->CreateGraph();
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
  TTreeReaderValue<unsigned int> run         (myReader,"run");
  TTreeReaderValue<unsigned int> lumisection (myReader,"lumi");
  TTreeReaderValue<unsigned int> event       (myReader,"event");
  TTreeReaderValue<unsigned int> nvtx        (myReader,"nvtx");
  TTreeReaderValue<int>          putrue      (myReader,"putrue");
  TTreeReaderValue<float> xsec               (myReader,"xsec");
  TTreeReaderValue<float> wgt                (myReader,"wgt");

  TTreeReaderValue<double>* wgtsum = NULL;
  if(isMC){
    wgtsum  = new TTreeReaderValue<double>(myReader,"wgtsum");
  }

  string wgtpuname = "wgtpileup";
  string wgtbtagname = "wgtbtag";
  if(not isMC){
    wgtpuname = "wgt";
    wgtbtagname = "wgt";
  }
  TTreeReaderValue<float> wgtpu              (myReader,wgtpuname.c_str());
  TTreeReaderValue<float> wgtbtag            (myReader,wgtbtagname.c_str());
  
  // take some wgt for MC events
  string prescalename;
  string hltphotonname;
  
  // trigger
  TTreeReaderValue<UChar_t> hltm90      (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100     (myReader,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110     (myReader,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120     (myReader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm90    (myReader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hltmwm120   (myReader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170   (myReader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300   (myReader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hlte        (myReader,"hltsingleel");
  TTreeReaderValue<UChar_t> hltenoiso   (myReader,"hltelnoiso");
  TTreeReaderValue<UChar_t> hltm        (myReader,"hltsinglemu");
  TTreeReaderValue<UChar_t> hltp120     (myReader,"hltphoton120");
  TTreeReaderValue<float>   pswgt_ph120 (myReader,"pswgt_ph120");
  TTreeReaderValue<UChar_t> hltp165     (myReader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175     (myReader,"hltphoton175");
  TTreeReaderValue<UChar_t> hltpPFHT800 (myReader,"hltPFHT800");

  /// met filters
  TTreeReaderValue<UChar_t> fhbhe  (myReader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (myReader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (myReader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (myReader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (myReader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (myReader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (myReader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (myReader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (myReader,"flagbadchpf");

  TTreeReaderValue<unsigned int> njets      (myReader,"njets");
  TTreeReaderValue<unsigned int> ntrigele   (myReader,"ntriggerelectrons");
  TTreeReaderValue<unsigned int> nbjets     (myReader,"nbjetslowpt");
  TTreeReaderValue<float> ht                (myReader,"ht");
  TTreeReaderValue<unsigned int> ntaus      (myReader,"ntausold");

  // AK8 jet
  TTreeReaderValue<vector<float> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<float> > boostedJetQGL   (myReader,"boostedJetQGL");
  TTreeReaderValue<vector<float> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<float> > boostedJetphi   (myReader,"boostedJetphi");
  TTreeReaderValue<vector<float> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<float> > prunedJetm      (myReader,"prunedJetm");
  TTreeReaderValue<vector<float> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<float> > boostedJettau1  (myReader,"boostedJettau1");

  TTreeReaderValue<float> hadBosoneta  (myReader,"wzeta_h");
  TTreeReaderValue<float> hadBosonphi  (myReader,"wzphi_h");
  TTreeReaderValue<float> hadBosonpt   (myReader,"wzpt_h");
  TTreeReaderValue<float> hadBosonm    (myReader,"wzmass_h");

  // met and jet systematics
  string metSuffix = "";
  string jetSuffix = "";
  // muon energy scale --> small effect
  if(sysName == "muUp")
    metSuffix = "MuEnUp";
  else if(sysName == "muDown" or sysName == "muDw")
    metSuffix = "MuEnDown";
  // electron energy scale --> small effect
  else if(sysName == "elUp")
    metSuffix = "ElEnUp";
  else if(sysName == "elDown" or sysName == "elDw")
    metSuffix = "ElEnDown";
  // photon energy scale --> small effect
  else if(sysName == "phoUp")
    metSuffix = "PhoEnUp";
  else if(sysName == "phoDown" or sysName == "phoDw")
    metSuffix = "PhoEnDown";
  // tau energy scale --> small effect
  else if(sysName == "tauUp")
    metSuffix = "TauEnUp";
  else if(sysName == "tauDown" or sysName == "tauDw")
    metSuffix = "TauEnDown";
  // jet energy scale
  else if(sysName == "jesUp"){
    metSuffix = "JetEnUp";
    jetSuffix = "up";
  }
  else if(sysName == "jesDown" or sysName == "jesDw"){
    metSuffix = "JetEnDown";
    jetSuffix = "dw";
  }
  else if(sysName == "jerUp"){
    metSuffix = "JetResUp";
    jetSuffix = "jer";
  }
  else if(sysName == "jerDown" or sysName == "jerDw"){
    metSuffix = "JetResDown";
    jetSuffix = "jer";
  }
  else if(sysName == "uncUp")
    metSuffix = "UncEnUp";
  else if(sysName == "uncDown" or sysName == "uncDw")
    metSuffix = "UncEnDown";
  
  TTreeReaderValue<vector<float> > jeteta  (myReader,("combinejeteta"+jetSuffix).c_str());
  TTreeReaderValue<vector<float> > jetpt   (myReader,("combinejetpt"+jetSuffix).c_str());
  TTreeReaderValue<vector<float> > jetphi  (myReader,("combinejetphi"+jetSuffix).c_str());
  TTreeReaderValue<vector<float> > jetm    (myReader,("combinejetm"+jetSuffix).c_str());
  TTreeReaderValue<vector<float> > jetbtag (myReader,"combinejetbtag");
  TTreeReaderValue<vector<float> > jetQGL  (myReader,"combinejetQGL");
  TTreeReaderValue<vector<float> > chfrac  (myReader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac  (myReader,"combinejetNHfrac");
  TTreeReaderValue<vector<float> > emfrac  (myReader,"combinejetEMfrac");
  TTreeReaderValue<vector<float> > pFlav   (myReader,"combinejetPFlav");
  TTreeReaderValue<unsigned int> nincjets  (myReader,("njetsinc"+jetSuffix).c_str());
  
  TTreeReaderValue<float> met         (myReader,("t1pfmet"+metSuffix).c_str());
  TTreeReaderValue<float> met_orig    (myReader,"t1pfmet");
  TTreeReaderValue<float> metphi      (myReader,"t1pfmetphi");
  TTreeReaderValue<float> mmet        (myReader,"t1mumet");
  TTreeReaderValue<float> mmetphi     (myReader,"t1mumetphi");
  TTreeReaderValue<float> emet        (myReader,"t1elmet");
  TTreeReaderValue<float> emetphi     (myReader,"t1elmetphi"); 
  TTreeReaderValue<float> pmet        (myReader,"t1phmet");
  TTreeReaderValue<float> pmetphi     (myReader,"t1phmetphi");
  TTreeReaderValue<float> tmet        (myReader,"t1pfmet");
  TTreeReaderValue<float> tmetphi     (myReader,"t1pfmetphi");
  TTreeReaderValue<float> metpf       (myReader,"pfmet");
  TTreeReaderValue<float> metcalo     (myReader,"calomet");

  ////////////////// dphi
  TTreeReaderValue<float> jmmdphi (myReader,("incjetmumetdphimin4"+jetSuffix).c_str());
  TTreeReaderValue<float> jemdphi (myReader,("incjetelmetdphimin4"+jetSuffix).c_str());
  TTreeReaderValue<float> jpmdphi (myReader,("incjetphmetdphimin4"+jetSuffix).c_str());

  string hltsafe1 = "el1idt";
  if(not useMoriondSetup)
    hltsafe1 = "el1id";
  string hltsafe2 = "el2idt";
  if(not useMoriondSetup)
    hltsafe2 = "el2id";

  ////////////////
  TTreeReaderValue<int>    mu1pid (myReader,"mu1pid");
  TTreeReaderValue<int>    mu2pid (myReader,"mu2pid");
  TTreeReaderValue<int>    mu1id  (myReader,"mu1id");
  TTreeReaderValue<int>    mu2id  (myReader,"mu2id");
  TTreeReaderValue<float>  mu1pt  (myReader,"mu1pt");
  TTreeReaderValue<float>  mu2pt  (myReader,"mu2pt");
  TTreeReaderValue<float>  mu1eta (myReader,"mu1eta");
  TTreeReaderValue<float>  mu2eta (myReader,"mu2eta");
  TTreeReaderValue<float>  mu1phi (myReader,"mu1phi");
  TTreeReaderValue<float>  mu2phi (myReader,"mu2phi");
  TTreeReaderValue<int>    el1pid (myReader,"el1pid");
  TTreeReaderValue<int>    el2pid (myReader,"el2pid");
  TTreeReaderValue<int>    el1id  (myReader,"el1id");
  TTreeReaderValue<int>    el1idt (myReader,hltsafe1.c_str());
  TTreeReaderValue<int>    el2id  (myReader,"el2id");
  TTreeReaderValue<int>    el2idt (myReader,hltsafe2.c_str());
  TTreeReaderValue<float> el1pt  (myReader,"el1pt");
  TTreeReaderValue<float> el2pt  (myReader,"el2pt");
  TTreeReaderValue<float> el1eta (myReader,"el1eta");
  TTreeReaderValue<float> el2eta (myReader,"el2eta");
  TTreeReaderValue<float> el1phi (myReader,"el1phi");
  TTreeReaderValue<float> el2phi (myReader,"el2phi");
  ////////////////
  TTreeReaderValue<float>  tau1id  (myReader,"tau1idold");
  TTreeReaderValue<float>  tau1pt  (myReader,"tau1pt");
  TTreeReaderValue<float>  tau1eta (myReader,"tau1eta");
  TTreeReaderValue<float>  tau1phi (myReader,"tau1phi");
  TTreeReaderValue<float>  tau1m   (myReader,"tau1m");
  ////////////////  
  TTreeReaderValue<int>   phidm (myReader,"phidm");
  TTreeReaderValue<float> phpt  (myReader,"phpt");
  TTreeReaderValue<float> pheta (myReader,"pheta");
  TTreeReaderValue<float> phphi (myReader,"phphi");
  ////////////////
  TTreeReaderValue<float> wmt    (myReader,"wmt");
  TTreeReaderValue<float> wemt   (myReader,"wemt");
  TTreeReaderValue<float> wtmt   (myReader,"wtmt");
  TTreeReaderValue<float> wzpt   (myReader,"wzpt");
  TTreeReaderValue<float> wzpt_h (myReader,"wzpt_h");
  TTreeReaderValue<float> wzeta  (myReader,"wzeta");
  TTreeReaderValue<float> zmass  (myReader,"zmass");
  TTreeReaderValue<float> zeemass (myReader,"zeemass");
  TTreeReaderValue<float> zmmpt  (myReader,"zpt");
  TTreeReaderValue<float> zeept  (myReader,"zeept");
  TTreeReaderValue<float> zeeeta (myReader,"zeeeta");
  TTreeReaderValue<float> zeephi (myReader,"zeephi");
  TTreeReaderValue<float> zmmeta (myReader,"zeta");
  TTreeReaderValue<float> zmmphi (myReader,"zphi");
  ///////////////
  TTreeReaderValue<float> l1eta (myReader,"l1eta");
  TTreeReaderValue<float> l1phi (myReader,"l1phi");
  TTreeReaderValue<int>   l1pid (myReader,"l1id");
  TTreeReaderValue<float> l2eta (myReader,"l2eta");
  TTreeReaderValue<float> l2phi (myReader,"l2phi");
  TTreeReaderValue<int>   l2pid (myReader,"l2id");
  ////////////////
  TTreeReaderValue<float> dmpt (myReader,"dmpt");

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

  ////////////////
  TTreeReaderValue<float> toppt  (myReader,topptname.c_str());
  TTreeReaderValue<float> atoppt (myReader,atopptname.c_str());

  // loop on events
  while(myReader.Next()){

    //ICHEP dataset
    if(not useMoriondSetup and not isMC and *run > 276811) continue;

    // check trigger depending on the sample
    Double_t hlt   = 0.0;
    Double_t hltw  = 1.0;
    
    if (sample == Sample::sig || sample == Sample::zmm || sample == Sample::wmn || sample == Sample::taun || sample == Sample::topmu || sample == Sample::qcd)// single and double muon
      hlt = *hltm90+*hltm100+*hltm110+*hltm120+*hltmwm90+*hltmwm120+*hltmwm170+*hltmwm300;
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel) // single and double electron
      hlt = *hlte+*hltenoiso;      
    else if (sample == Sample::qcdgam || sample == Sample::gam){ // single photon
      hlt  = *hltp165+*hltp175+*hltpPFHT800;
    }

    // Trigger Selection
    if(hlt  == 0) continue; // trigger    

    // MET Filters --> apply on both data and monte-carlo
    if((category != Category::VBF and category != Category::VBFrelaxed) and 
       (*fhbhe == 0 or *fhbiso == 0 or *feeb == 0 or *fetp == 0 or *fvtx == 0 or *fcsc == 0 or *fbadmu == 0 or *fbadch == 0)) continue;
    if((category == Category::VBF or category == Category::VBFrelaxed) and not isMC and 
       (*fhbhe == 0 or *fhbiso == 0 or *feeb == 0 or *fetp == 0 or *fvtx == 0 or *fcsc == 0 or *fbadmu == 0 or *fbadch == 0)) continue;
    
    // check dphi jet-met
    Double_t jmdphi = 0.0;    
    if (sample == Sample::sig || sample == Sample::wmn || sample == Sample::zmm || sample == Sample::topmu || sample == Sample::qcd || sample == Sample::taun) jmdphi = fabs(*jmmdphi);
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel) jmdphi = fabs(*jemdphi);
    else if (sample == Sample::qcdgam || sample == Sample::gam) jmdphi = fabs(*jpmdphi);
    
    //set met
    Double_t pfmet = 0.0;
    Double_t pfmetphi = 0.0;
    if (sample == Sample::sig || sample == Sample::qcd) {pfmet = *met; pfmetphi = *metphi;}
    else if (sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu){ pfmet = *mmet; pfmetphi = *mmetphi;}
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel){ pfmet = *emet; pfmetphi = *emetphi;}
    else if (sample == Sample::qcdgam || sample == Sample::gam)  { pfmet = *pmet; pfmetphi = *pmetphi;}
    else if (sample == Sample::taun){ pfmet = *tmet; pfmetphi = *tmetphi;}

    // calculate correction on recoil by hand
    if(metSuffix != "" and sample != Sample::sig){
      pfmet += (*met-*met_orig)/2;
    }
    
    // set lepton info
    Int_t    id1   = 0, id2   = 0, id1t  = 0, id2t  = 0;
    Double_t pt1   = 0.0, pt2   = 0.0;
    Double_t eta1  = 0.0, eta2  = 0.0;
    Double_t phi1  = 0.0, phi2  = 0.0;
    int      pid1  = 0, pid2  = 0;
    int leadingjet_notau = -1;
    int trailingjet_notau = -1;

    if (sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu) {
      id1  = *mu1id;  id2  = *mu2id;
      id1t = 1;       id2t = 1;
      pt1  = *mu1pt;  pt2  = *mu2pt;
      pid1 = *mu1pid; pid2 = *mu2pid;
      eta1 = *mu1eta; eta2 = *mu2eta;
      phi1 = *el1phi; phi2 = *el2phi;
    }
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel) {
      id1  = *el1id; id2  = *el2id;
      id1t = *el1idt;
      pt1  = *el1pt;  pt2  = *el2pt;
      eta1 = *el1eta; eta2 = *el2eta;
      phi1 = *el1phi; phi2 = *el2phi;
      pid1 = *el1pid; pid2 = *el2pid;
    }
    else if (sample == Sample::qcdgam || sample == Sample::gam) {
      id1  = *phidm; id2  = 1.0;
      pt1  = *phpt;
      if(applyPhotonScale and not isMC){	
	pt1 += pt1*photonScaleUnc;
      }
      eta1 = *pheta;
    }
    else if (sample == Sample::taun){
      pt1  = *tau1pt;
      eta1 = *tau1eta; phi1 = *tau1phi;
      id1 = *tau1id;
      id2 = 1.0;
      TLorentzVector tau4V; tau4V.SetPtEtaPhiM(*tau1pt,*tau1eta,*tau1phi,*tau1m);
      int ijetmax = 0;
      for(size_t ijet = 0; ijet < jetpt->size(); ijet++){ // clean for jets
	TLorentzVector jet4V;
	jet4V.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
	if(jet4V.DeltaR(tau4V) < 0.4) continue;
	if(leadingjet_notau == -1) leadingjet_notau = ijet;
	else if(leadingjet_notau != -1 and trailingjet_notau == -1) trailingjet_notau = ijet;
	if(jetpt->at(ijet) < 30) continue; // only jets with pt > 30 GeV	
	ijetmax++;
	if(ijetmax > 4) continue;
	TVector2 met2D; met2D.SetMagPhi(pfmet,pfmetphi);
	TVector2 tau2D; tau2D.SetMagPhi(pt1,phi1);
	TVector2 jet2D; jet2D.SetMagPhi(jetpt->at(ijet),jetphi->at(ijet));
	jmdphi = fabs((met2D+tau2D).DeltaPhi(jet2D));	
      }
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
    else if (sample == Sample::wen or sample == Sample::wmn){ // single muon or single ele
      TLorentzVector lep4V, met4V;
      lep4V.SetPtEtaPhiM(pt1,eta1,phi1,0.);
      met4V.SetPtEtaPhiM(*met,0.,*metphi,0.);
      bosonPt = (lep4V+met4V).Pt();
      bosonEta = 0;
      bosonPhi = (lep4V+met4V).Phi();
    }

    // control regions with two leptons --> one should be tight
    if ((sample == Sample::zmm || sample == Sample::zee)){
      if(sample == Sample::zmm){	
	if(not ((pt1 > 20 and id1 == 1) or (pt2 > 20 and id2 == 1))) continue;	
      }
      else if(sample == Sample::zee){
	if(not ((pt1 > 40 and id1 == 1) or (pt2 > 40 and id2 == 1))) continue;
      }
    }

    // control regions with two leptons --> opposite charge
    if (sample == Sample::zmm && *mu1pid == *mu2pid) continue;
    if (sample == Sample::zee && *el1pid == *el2pid) continue;

    /// Tau-veto or tau tag
    if (*ntaus != 0 and sample != Sample::taun) continue;
    else if (*ntaus != 1 and sample == Sample::taun) continue;

    // B-veto, not for top control sample
    if (*nbjets > 0 and sample != Sample::topmu and sample != Sample::topel) continue; 

    // control regions wit one lepton --> tight requirement 
    if ((sample == Sample::wen || sample == Sample::wmn) && id1 !=1) continue;
    if (sample  == Sample::wen and *wemt > 160) continue;
    if (sample  == Sample::wmn and *wmt  > 160) continue;

    // photon control sample
    if ((sample == Sample::qcdgam || sample == Sample::gam) && pt1 < photonPt) continue;
    if ((sample == Sample::qcdgam || sample == Sample::gam) && fabs(*pheta) > 1.4442) continue;    

    // Wenu kill QCD
    if (category != Category::twojet and category != Category::VBF and category != Category::VBFrelaxed and sample == Sample::wen && *met < 50.) continue;
    else if((category  == Category::twojet or category == Category::VBF or category == Category::VBFrelaxed) and sample == Sample::wen && *met < 60.) continue;
    
    // tau-nu control region
    if(sample == Sample::taun){
      if(fabs(eta1) > 2.3) continue;
      if(pt1 < 30) continue;
      if(id1 != 1) continue;
      if(*wtmt > 160) continue;
    }

    // n-bjets cut for unboosted categories
    if ((sample == Sample::topmu || sample == Sample::topel) && 
	(category != Category::monoV and category != Category::boosted and category != Category::prunedMass and category != Category::tau2tau1)  && 
	*nbjets < nBjets) continue;

    if (sample == Sample::topmu || sample == Sample::topel){ // select only events with one lepton
      // at least one lepton in the plateau region
      if(id1 != 1) continue;
      if(abs(pid1) == 13 && pt1 < 20. ) continue;
      if(abs(pid1) == 11 && pt1 < 40. ) continue;
      // met cut
      if(sample == Sample::topel && *met < 50.) continue;
      // veto di-lepton events
      if(pt2 > 0) continue;
    }

    // met selection
    if(category == Category::VBF or category == Category::twojet or category == Category::VBFrelaxed){
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

    // noise cleaner
    if(category != Category::VBF and category != Category::VBFrelaxed and category != Category::twojet){
      if(sample != Sample::gam and sample != Sample::zee and sample != Sample::wen and fabs(*met-*metcalo)/pfmet > 0.5) continue;
    }
    else if((category == Category::VBF or category == Category::VBFrelaxed or category == Category::twojet) and fabs(*met-*metcalo)/pfmet > 0.5) continue;

    // number of central jets
    if (category != Category::VBF and category != Category::twojet and category != Category::VBFrelaxed and *njets < 1) continue; 
    else if((category == Category::VBF or category == Category::twojet or category == Category::VBFrelaxed) and *nincjets < 2) continue;


    //apply charge cut to separate W+ from W- in case
    if(sample == Sample::wmn or sample == Sample::wen){ // in case charge is required to be 1 skip events with negative leptons, viceversa
      if(vBosonCharge == 1 and pid1 > 0) continue;
      else if(vBosonCharge == -1 and pid1 < 0) continue;
    }

    // selection on gen-leptons --> use it only for V+jets sample
    if(leptonPID != ""){
      if((leptonPID == "muon" or leptonPID == "mu") and (fabs(*l1pid) != 13 and fabs(*l1pid) != 14)) continue;
      if((leptonPID == "muon" or leptonPID == "mu") and (fabs(*l2pid) != 13 and fabs(*l2pid) != 14)) continue;
      if((leptonPID == "electron" or leptonPID == "el") and (fabs(*l1pid) != 11 and fabs(*l1pid) != 12)) continue;
      if((leptonPID == "electron" or leptonPID == "el") and (fabs(*l2pid) != 11 and fabs(*l2pid) != 12)) continue;
      if(leptonPID == "tau" and (fabs(*l1pid) != 15 and fabs(*l1pid) != 16)) continue;
      if(leptonPID == "tau" and (fabs(*l2pid) != 15 and fabs(*l2pid) != 16)) continue;
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
      else if(fabs(jeteta->at(ijet)) <= 2.5 and jetpt->at(ijet) > 30){
      	vect.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
      	centralJets.push_back(vect);
      	if(leadingCentralJetPos == -1)
      	  leadingCentralJetPos = ijet;
      }
    }
    
    if(category != Category::VBF and category != Category::twojet and category != Category::VBFrelaxed and leadingCentralJetPos < 0)  continue;
    if(category != Category::VBF and category != Category::twojet and category != Category::VBFrelaxed and leadingCentralJetPos != 0) continue; 
        
    /// re-miniADO specific to adjust lumi
    Double_t sfwgt = 1.0;
    //float sf_track = 1.0;
    // apply tracking efficiency for electrons from POGs / private files
    if(isMC && (sample == Sample::zee or sample == Sample::wen)){
      if(not isSummer16){
	if(pt1 > 0. and *nvtx <= numberOfVtxCorrection){	
	  if(pt1 > trackingefficiency_electron_lowpu->GetYaxis()->GetBinLowEdge(1)) // above minimum
	    sfwgt *= trackingefficiency_electron_lowpu->GetBinContent(trackingefficiency_electron_lowpu->FindBin(eta1,min(pt1,trackingefficiency_electron_lowpu->GetYaxis()->GetBinLowEdge(trackingefficiency_electron_lowpu->GetNbinsY()+1)-1)));
	  else
	    sfwgt *= trackingefficiency_electron_lowpu->GetBinContent(trackingefficiency_electron_lowpu->FindBin(eta1,trackingefficiency_electron_lowpu->GetYaxis()->GetBinLowEdge(1)+1));	
	}
	if(pt2 > 0. and *nvtx <= numberOfVtxCorrection){
	  if(pt2 > trackingefficiency_electron_lowpu->GetYaxis()->GetBinLowEdge(1)) // above minimum
	    sfwgt *= trackingefficiency_electron_lowpu->GetBinContent(trackingefficiency_electron_lowpu->FindBin(eta2,min(pt2,trackingefficiency_electron_lowpu->GetYaxis()->GetBinLowEdge(trackingefficiency_electron_lowpu->GetNbinsY()+1)-1)));
	  else
	    sfwgt *= trackingefficiency_electron_lowpu->GetBinContent(trackingefficiency_electron_lowpu->FindBin(eta2,trackingefficiency_electron_lowpu->GetYaxis()->GetBinLowEdge(1)+1));
	}
	if(pt1 > 0. and *nvtx > numberOfVtxCorrection){
	  if(pt1 > trackingefficiency_electron_highpu->GetYaxis()->GetBinLowEdge(1)) // above minimum
	    sfwgt *= trackingefficiency_electron_highpu->GetBinContent(trackingefficiency_electron_highpu->FindBin(eta1,min(pt1,trackingefficiency_electron_highpu->GetYaxis()->GetBinLowEdge(trackingefficiency_electron_highpu->GetNbinsY()+1)-1)));
	  else
	    sfwgt *= trackingefficiency_electron_highpu->GetBinContent(trackingefficiency_electron_highpu->FindBin(eta1,trackingefficiency_electron_highpu->GetYaxis()->GetBinLowEdge(1)+1));
	}	
	if(pt2 > 0. and *nvtx > numberOfVtxCorrection){ // above minimum
	  if(pt2 > trackingefficiency_electron_highpu->GetYaxis()->GetBinLowEdge(1))
	    sfwgt *= trackingefficiency_electron_highpu->GetBinContent(trackingefficiency_electron_highpu->FindBin(eta2,min(pt2,trackingefficiency_electron_highpu->GetYaxis()->GetBinLowEdge(trackingefficiency_electron_highpu->GetNbinsY()+1)-1))); 
	  else
	    sfwgt *= trackingefficiency_electron_highpu->GetBinContent(trackingefficiency_electron_highpu->FindBin(eta2,trackingefficiency_electron_highpu->GetYaxis()->GetBinLowEdge(1)+1));
	}
      }
      else{// Summer 16 Corrention
	if(pt1 > 0.){	  
	  float ptVal = pt1;
	  if(pt1 < trackingefficiency_electron->GetYaxis()->GetBinLowEdge(1)) ptVal =  trackingefficiency_electron->GetYaxis()->GetBinLowEdge(1)+1;
	  else if(pt1 > trackingefficiency_electron->GetYaxis()->GetBinLowEdge(trackingefficiency_electron->GetNbinsY()+1)) ptVal = trackingefficiency_electron->GetYaxis()->GetBinLowEdge(trackingefficiency_electron->GetNbinsY()+1)-1;
	  sfwgt *= trackingefficiency_electron->GetBinContent(trackingefficiency_electron->FindBin(eta1,ptVal));
	  //sf_track =  trackingefficiency_electron->GetBinContent(trackingefficiency_electron->FindBin(eta1,ptVal));
	}
	if(pt2 > 0.){
	  float ptVal = pt1;
	  if(pt2 < trackingefficiency_electron->GetYaxis()->GetBinLowEdge(1)) ptVal =  trackingefficiency_electron->GetYaxis()->GetBinLowEdge(1)+1;
	  else if(pt2 > trackingefficiency_electron->GetYaxis()->GetBinLowEdge(trackingefficiency_electron->GetNbinsY()+1)) ptVal = trackingefficiency_electron->GetYaxis()->GetBinLowEdge(trackingefficiency_electron->GetNbinsY()+1)-1;
	  sfwgt *= trackingefficiency_electron->GetBinContent(trackingefficiency_electron->FindBin(eta2,ptVal));
	  //sf_track =  trackingefficiency_electron->GetBinContent(trackingefficiency_electron->FindBin(eta1,ptVal));
	}
      }
    }

    // reco-muon scale factor --> private scale factors aleady starting at 10 GeV
    if(isMC && (sample == Sample::zmm or sample == Sample::wmn or sample == Sample::topmu)){
      double trackwgt = 1;
      if(not usePOGScaleFactors){ // private SF binned in nvtx (2bins), eta and pt
	if(pt1 > 0. and *nvtx <= numberOfVtxCorrection)
	  trackwgt *= trackingefficiency_muon_lowpu->GetBinContent(trackingefficiency_muon_lowpu->FindBin(fabs(eta1),min(pt1,trackingefficiency_muon_lowpu->GetYaxis()->GetBinLowEdge(trackingefficiency_muon_lowpu->GetNbinsY()+1)-1)));
	if(pt2 > 0. and *nvtx <= numberOfVtxCorrection)
	  trackwgt *= trackingefficiency_muon_lowpu->GetBinContent(trackingefficiency_muon_lowpu->FindBin(fabs(eta2),min(pt2,trackingefficiency_muon_lowpu->GetYaxis()->GetBinLowEdge(trackingefficiency_muon_lowpu->GetNbinsY()+1)-1)));
	if(pt1 > 0. and *nvtx > numberOfVtxCorrection)
	  trackwgt *= trackingefficiency_muon_highpu->GetBinContent(trackingefficiency_muon_highpu->FindBin(fabs(eta1),min(pt1,trackingefficiency_muon_highpu->GetYaxis()->GetBinLowEdge(trackingefficiency_muon_highpu->GetNbinsY()+1)-1)));
	if(pt2 > 0. and *nvtx > numberOfVtxCorrection)
	  trackwgt *= trackingefficiency_muon_highpu->GetBinContent(trackingefficiency_muon_highpu->FindBin(fabs(eta2),min(pt2,trackingefficiency_muon_highpu->GetYaxis()->GetBinLowEdge(trackingefficiency_muon_highpu->GetNbinsY()+1)-1)));      
      }
      else{
	if(pt1 > 0.)
	  trackwgt *= trackingefficiency_pog->Eval(*nvtx);
	if(pt2 > 0.)
	  trackwgt *= trackingefficiency_pog->Eval(*nvtx);
      }      
      sfwgt *= trackwgt;
    }

    // scale factor for leptons
    TH2* sflhist_lowpu = NULL;
    TH2* sfthist_lowpu = NULL;
    TH2* sflhist_highpu = NULL;
    TH2* sfthist_highpu = NULL;

    //double sf_lepID = 1;

    if (sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu) {

      if(not usePOGScaleFactors){
	if (pt1 > 0.) {	
	  if(*nvtx <= numberOfVtxCorrection){
	    if (id1 == 1)
	      sfwgt *= msftight_lowpu->GetBinContent(msftight_lowpu->FindBin(fabs(eta1),min(pt1,msftight_lowpu->GetYaxis()->GetBinLowEdge(msftight_lowpu->GetNbinsY()+1)-1))); 
	    else
	      sfwgt *= msfloose_lowpu->GetBinContent(msfloose_lowpu->FindBin(fabs(eta1),min(pt1,msfloose_lowpu->GetYaxis()->GetBinLowEdge(msfloose_lowpu->GetNbinsY()+1)-1))); 
	  }
	  else{
	    if (id1 == 1)
	      sfwgt *= msftight_highpu->GetBinContent(msftight_highpu->FindBin(fabs(eta1),min(pt1,msftight_highpu->GetYaxis()->GetBinLowEdge(msftight_highpu->GetNbinsY()+1)-1))); 
	    else
	      sfwgt *= msfloose_highpu->GetBinContent(msfloose_highpu->FindBin(fabs(eta1),min(pt1,msfloose_highpu->GetYaxis()->GetBinLowEdge(msfloose_highpu->GetNbinsY()+1)-1))); 
	  }
	}
	if (pt2 > 0.) {
	  if(*nvtx <= numberOfVtxCorrection){
	    if (id2 == 1)
	      sfwgt *= msftight_lowpu->GetBinContent(msftight_lowpu->FindBin(fabs(eta2),min(pt2,msftight_lowpu->GetYaxis()->GetBinLowEdge(msftight_lowpu->GetNbinsY()+1)-1))); 
	    else
	      sfwgt *= msfloose_lowpu->GetBinContent(msfloose_lowpu->FindBin(fabs(eta2),min(pt2,msfloose_lowpu->GetYaxis()->GetBinLowEdge(msfloose_lowpu->GetNbinsY()+1)-1))); 
	  }
	  else{
	    if (id2 == 1){
	      sfwgt *= msftight_highpu->GetBinContent(msftight_highpu->FindBin(fabs(eta2),min(pt2,msftight_highpu->GetYaxis()->GetBinLowEdge(msftight_highpu->GetNbinsY()+1)-1))); 
	    }
	    else
	      sfwgt *= msfloose_highpu->GetBinContent(msfloose_highpu->FindBin(fabs(eta2),min(pt2,msfloose_highpu->GetYaxis()->GetBinLowEdge(msfloose_highpu->GetNbinsY()+1)-1))); 
	  }
	}      
      }
      else{ // official POG scale factors

	if (pt1 > 0.){
	  float ptValue = pt1;
	  if(ptValue < msftight_id->GetYaxis()->GetBinLowEdge(1)) ptValue =  msftight_id->GetYaxis()->GetBinLowEdge(1)+1;
	  else if(ptValue > msftight_id->GetYaxis()->GetBinLowEdge(msftight_id->GetNbinsY()+1)) ptValue = msftight_id->GetYaxis()->GetBinLowEdge(msftight_id->GetNbinsY()+1)-1; 
	  if(id1 == 1)
	    sfwgt *= msftight_id->GetBinContent(msftight_id->FindBin(fabs(eta1),ptValue))*msftight_iso->GetBinContent(msftight_iso->FindBin(fabs(eta1),ptValue));
	  else 
	    sfwgt *= msfloose_id->GetBinContent(msfloose_id->FindBin(fabs(eta1),ptValue))*msfloose_iso->GetBinContent(msfloose_iso->FindBin(fabs(eta1),ptValue));
	}
	if(pt2 > 0.){
	  float ptValue = pt2;
	  if(ptValue < msftight_id->GetYaxis()->GetBinLowEdge(1)) ptValue =  msftight_id->GetYaxis()->GetBinLowEdge(1)+1;
	  else if(ptValue > msftight_id->GetYaxis()->GetBinLowEdge(msftight_id->GetNbinsY()+1)) ptValue = msftight_id->GetYaxis()->GetBinLowEdge(msftight_id->GetNbinsY()+1)-1; 
	  if(id2 == 1)
	    sfwgt *= msftight_id->GetBinContent(msftight_id->FindBin(fabs(eta2),ptValue))*msftight_iso->GetBinContent(msftight_iso->FindBin(fabs(eta2),ptValue));
	  else 
	    sfwgt *= msfloose_id->GetBinContent(msfloose_id->FindBin(fabs(eta2),ptValue))*msfloose_iso->GetBinContent(msfloose_iso->FindBin(fabs(eta1),ptValue));
	}
      }
    }

    //////////////    
    if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel) {      
      if(isSummer16){ // use e-gamma pog scale factors
	if (pt1 > 0.) {
	  float ptVal = pt1;
	  if(pt1 < esftight->GetYaxis()->GetBinLowEdge(1)) ptVal = esftight->GetYaxis()->GetBinLowEdge(1)+1;
	  else if(pt1 > esftight->GetYaxis()->GetBinLowEdge(esftight->GetNbinsY()+1)) ptVal = esftight->GetYaxis()->GetBinLowEdge(esftight->GetNbinsY()+1)-1;
	  if (id1 == 1) {
	    sfwgt *= esftight->GetBinContent(esftight->FindBin(eta1,ptVal));
	    //sf_lepID *= esftight->GetBinContent(esftight->FindBin(eta1,ptVal));
	  }
	  else{  
	    sfwgt *= esfveto->GetBinContent(esfveto->FindBin(eta1,ptVal));
	    //sf_lepID *=esfveto->GetBinContent(esfveto->FindBin(eta1,ptVal));
	  }
	}
	if (pt2 > 0.) {
	  float ptVal = pt2;
	  if(pt2 < esftight->GetYaxis()->GetBinLowEdge(1)) ptVal = esftight->GetYaxis()->GetBinLowEdge(1)+1;
	  else if(pt2 > esftight->GetYaxis()->GetBinLowEdge(esftight->GetNbinsY()+1)) ptVal = esftight->GetYaxis()->GetBinLowEdge(esftight->GetNbinsY()+1)-1;
	  if (id1 == 1){
	    sfwgt *= esftight->GetBinContent(esftight->FindBin(eta2,ptVal));
	    //sf_lepID *=esftight->GetBinContent(esftight->FindBin(eta2,ptVal));
	  }
	  else{
	    sfwgt *= esfveto->GetBinContent(esfveto->FindBin(eta2,ptVal));
	    //sf_lepID  *= esfveto->GetBinContent(esfveto->FindBin(eta2,ptVal));
	  }
	}
      }
      else{
	// use private SF for Spring 16
	if (pt1 > 0.) {
	  if(*nvtx <= numberOfVtxCorrection){
	    if (id1 == 1) sfwgt *= esftight_lowpu->GetBinContent(esftight_lowpu->FindBin(fabs(eta1),min(pt1,esftight_lowpu->GetYaxis()->GetBinLowEdge(esftight_lowpu->GetNbinsY()+1)-1))); 
	    else sfwgt *= esfveto_lowpu->GetBinContent(esfveto_lowpu->FindBin(fabs(eta1),min(pt1,esfveto_lowpu->GetYaxis()->GetBinLowEdge(esfveto_lowpu->GetNbinsY()+1)-1))); 
	  }
	  else{
	    if (id1 == 1) sfwgt *= esftight_highpu->GetBinContent(esftight_highpu->FindBin(fabs(eta1),min(pt1,esftight_highpu->GetYaxis()->GetBinLowEdge(esftight_highpu->GetNbinsY()+1)-1))); 
	    else sfwgt *= esfveto_highpu->GetBinContent(esfveto_highpu->FindBin(fabs(eta1),min(pt1,esfveto_highpu->GetYaxis()->GetBinLowEdge(esfveto_highpu->GetNbinsY()+1)-1))); 
	  }
	}
	if (pt2 > 0.) {
	  if(*nvtx <= numberOfVtxCorrection){
	    if (id2 == 1) sfwgt *= esftight_lowpu->GetBinContent(esftight_lowpu->FindBin(fabs(eta2),min(pt2,esftight_lowpu->GetYaxis()->GetBinLowEdge(esftight_lowpu->GetNbinsY()+1)-1))); 
	    else sfwgt *= esfveto_lowpu->GetBinContent(esfveto_lowpu->FindBin(fabs(eta2),min(pt2,esfveto_lowpu->GetYaxis()->GetBinLowEdge(esfveto_lowpu->GetNbinsY()+1)-1))); 
	  }
	  else{
	    if (id2 == 1) sfwgt *= esftight_highpu->GetBinContent(esftight_highpu->FindBin(fabs(eta2),min(pt2,esftight_highpu->GetYaxis()->GetBinLowEdge(esftight_highpu->GetNbinsY()+1)-1))); 
	    else sfwgt *= esfveto_highpu->GetBinContent(esfveto_highpu->FindBin(fabs(eta2),min(pt2,esfveto_highpu->GetYaxis()->GetBinLowEdge(esfveto_highpu->GetNbinsY()+1)-1))); 
	  }	  
	}
      }
    }

    // tight tau scale factor--> as suggested by tau pog https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendation13TeV
    if(isMC && sample == Sample::taun){
      sfwgt *= 0.95;
    }
    
    // trigger scale factor for electrons
    float sf_trigger = 1.;
    if (isMC && triggerelhist && triggerelhist_ht && (sample == Sample::zee || sample == Sample::topel || sample == Sample::wen)) {
      float sf1 = 1.;
      float sf2 = 1.;

      if(category != Category::VBF and category != Category::VBFrelaxed){
	sf1 = triggerelhist->GetBinContent(triggerelhist->FindBin(fabs(eta1),min(pt1,triggerelhist->GetYaxis()->GetBinLowEdge(triggerelhist->GetNbinsY()+1)-1)));
	if(pt1 >= 200)
	  sf1 = triggerelhist_ht->GetBinContent(triggerelhist_ht->FindBin(fabs(eta1),min(pt1,triggerelhist_ht->GetYaxis()->GetBinLowEdge(triggerelhist_ht->GetNbinsY()+1)-1)));
	sf2 = triggerelhist->GetBinContent(triggerelhist->FindBin(fabs(eta2),min(pt2,triggerelhist->GetYaxis()->GetBinLowEdge(triggerelhist->GetNbinsY()+1)-1)));
	if(pt2 >= 200)
	  sf2 = triggerelhist_ht->GetBinContent(triggerelhist_ht->FindBin(fabs(eta2),min(pt2,triggerelhist_ht->GetYaxis()->GetBinLowEdge(triggerelhist_ht->GetNbinsY()+1)-1)));	  
      }
      else{

	float ptVal = pt1;
	if(pt1 > triggerelhist->GetYaxis()->GetBinLowEdge(triggerelhist->GetNbinsY()+1)) ptVal = triggerelhist->GetYaxis()->GetBinLowEdge(triggerelhist->GetNbinsY()+1)-1;	
	sf1 = triggerelhist->GetBinContent(triggerelhist->FindBin(eta1,ptVal));
	if(pt1 >= 200){
	  ptVal = pt1;
	  if(ptVal > triggerelhist_ht->GetYaxis()->GetBinLowEdge(triggerelhist_ht->GetNbinsY()+1)) ptVal = triggerelhist_ht->GetYaxis()->GetBinLowEdge(triggerelhist_ht->GetNbinsY()+1)-1;
	  sf1 = triggerelhist_ht->GetBinContent(triggerelhist->FindBin(eta1,ptVal));
	}

	ptVal = pt2;
	if(ptVal > triggerelhist->GetYaxis()->GetBinLowEdge(triggerelhist->GetNbinsY()+1)) ptVal = triggerelhist->GetYaxis()->GetBinLowEdge(triggerelhist->GetNbinsY()+1)-1;	
	sf2 = triggerelhist->GetBinContent(triggerelhist->FindBin(eta2,ptVal));
	if(pt2 >= 200){
	  ptVal = pt2;
	  if(ptVal > triggerelhist_ht->GetYaxis()->GetBinLowEdge(triggerelhist_ht->GetNbinsY()+1)) ptVal = triggerelhist_ht->GetYaxis()->GetBinLowEdge(triggerelhist_ht->GetNbinsY()+1)-1;
	  sf2 = triggerelhist_ht->GetBinContent(triggerelhist->FindBin(eta2,ptVal));
	}	

      }

      //////////////
      if (pt1 > 40. && id1 == 1 and id2 == 1){
	sfwgt *= min(1.,double(sf1+sf2-sf1*sf2));
	sf_trigger *= min(1.,double(sf1+sf2-sf1*sf2));
      }
      else if(pt1 > 40 and id1 == 1 and id2 != 1){	
	sfwgt *= sf1;
	sf_trigger *= sf1;
      }
      else if(pt2 > 40 and id2 == 1 and id1 != 1){
	sfwgt *= sf2;
	sf_trigger *= sf2;
      }
    }
    
    
    // photon id scale factor
    if (isMC && psfmedium_lowpu && psfmedium_highpu && sample == Sample::gam) {
      if(not isSummer16){
	if (pt1 > 0. && id1 == 1) {
	  if(*nvtx <= numberOfVtxCorrection)
	    sfwgt *= psfmedium_lowpu->GetBinContent(psfmedium_lowpu->FindBin(fabs(eta1),min(pt1,psfmedium_lowpu->GetYaxis()->GetBinLowEdge(psfmedium_lowpu->GetNbinsY()+1)-1)));
	  else
	    sfwgt *= psfmedium_highpu->GetBinContent(psfmedium_highpu->FindBin(fabs(eta1),min(pt1,psfmedium_highpu->GetYaxis()->GetBinLowEdge(psfmedium_highpu->GetNbinsY()+1)-1)));
	}
      }
      else{
	if (pt1 > 0. && id1 == 1) {
	  sfwgt *= psfmedium->GetBinContent(psfmedium_lowpu->FindBin(eta1,min(pt1,psfmedium_lowpu->GetYaxis()->GetBinLowEdge(psfmedium_lowpu->GetNbinsY()+1)-1)));
	}
      }
    }

    // photon purity
    if (!isMC && sample == Sample::qcdgam) {
      if(useMoriondSetup){
	double xleft, yleft, xright, yright;
	purgraph->GetPoint(0,xleft,yleft);
	purgraph->GetPoint(purgraph->GetN()-1,xright,yright);	
	if(pt1 >= xleft and pt1 <= xright) // check graph extremes
	  sfwgt *= (1.0 - purgraph->Eval(pt1));
	else if(pt1 < xleft)
	  sfwgt *= (1.0 - purgraph->Eval(xleft));
	else if(pt1 > xright)
	  sfwgt *= (1.0 - purgraph->Eval(xright));
      }      
      else{
	if(pt1 > purhist->GetXaxis()->GetBinLowEdge(1))
	  sfwgt *= (1.0 - purhist->GetBinContent(purhist->FindBin(min(pt1,purhist->GetXaxis()->GetBinLowEdge(purhist->GetNbinsX()+1)-1), fabs(eta1))));
	else
	  sfwgt *= (1.0 - purhist->GetBinContent(purhist->FindBin(max(pt1,purhist->GetXaxis()->GetBinLowEdge(1)),fabs(eta1))));		  
      }
    }

    // met trigger scale factor
    if (isMC && (sample == Sample::sig || sample == Sample::wmn || sample == Sample::zmm || sample == Sample::topmu || sample == Sample::qcd || sample == Sample::taun)) {
      // single trigger turn on to be applied
      if(triggermet_graph and (sample == Sample::wmn ||  sample == Sample::topmu || sample == Sample::qcd || sample == Sample::taun || sample == Sample::sig))
	sfwgt *= triggermet_graph->Eval(min(pfmet,triggermet_graph->GetXaxis()->GetXmax()));
      // for Zmm
      else if(triggermet_graph_zmm and sample == Sample::zmm)
	sfwgt *= triggermet_graph_zmm->Eval(min(pfmet,triggermet_graph_zmm->GetXaxis()->GetXmax()));
      else if(triggermet_graph and sample == Sample::zmm)
	sfwgt *= triggermet_graph->Eval(min(pfmet,triggermet_graph->GetXaxis()->GetXmax()));      
      // for VBF
      else if(category == Category::VBF or category == Category::twojet or category == Category::VBFrelaxed){
	if(centralJets.size()+forwardJets.size() >= 2){
	  TLorentzVector jet1 ;
	  TLorentzVector jet2 ;
	  jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	  jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	  if(sample == Sample::zmm){
	    if((jet1+jet2).M() < 800)
	      sfwgt *= triggermet_func_binned_Zmm.at(0)->Eval(min(pfmet,triggermet_func_binned_Zmm.at(0)->GetXaxis()->GetXmax()));
	    else if((jet1+jet2).M() >= 800 and (jet1+jet2).M() < 1200)
	      sfwgt *= triggermet_func_binned_Zmm.at(1)->Eval(min(pfmet,triggermet_func_binned_Zmm.at(1)->GetXaxis()->GetXmax()));
	    else if((jet1+jet2).M() >= 1200 and (jet1+jet2).M() < 1700)
	      sfwgt *= triggermet_func_binned_Zmm.at(2)->Eval(min(pfmet,triggermet_func_binned_Zmm.at(2)->GetXaxis()->GetXmax()));
	    else if((jet1+jet2).M() >= 1700)
	      sfwgt *= triggermet_func_binned_Zmm.at(3)->Eval(min(pfmet,triggermet_func_binned_Zmm.at(3)->GetXaxis()->GetXmax()));	  
	  }
	  else{
	    if((jet1+jet2).M() < 800)
	      sfwgt *= triggermet_func_binned_Wmn.at(0)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(0)->GetXaxis()->GetXmax()));
	    else if((jet1+jet2).M() >= 800 and (jet1+jet2).M() < 1200)
	      sfwgt *= triggermet_func_binned_Wmn.at(1)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(1)->GetXaxis()->GetXmax()));
	    else if((jet1+jet2).M() >= 1200 and (jet1+jet2).M() < 1700)
	      sfwgt *= triggermet_func_binned_Wmn.at(2)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(2)->GetXaxis()->GetXmax()));
	    else if((jet1+jet2).M() >= 1700)
	      sfwgt *= triggermet_func_binned_Wmn.at(3)->Eval(min(pfmet,triggermet_func_binned_Wmn.at(3)->GetXaxis()->GetXmax()));	  
	  }
	}
      }
    }
    
    // photon trigger scale factor
    if(isMC && triggerphoton_graph && triggerphoton_graph_jetHT and sample == Sample::gam){ // linear interpolation between graph points            
      if(*pmet < recoilThresholdTrigger)	
	sfwgt *= triggerphoton_graph->Eval(min(double(*phpt),triggerphoton_graph->GetXaxis()->GetXmax()));
      else
	sfwgt *= triggerphoton_graph_jetHT->Eval(min(double(*phpt),triggerphoton_graph->GetXaxis()->GetXmax()));
    }
    
    // B-tag weight to be adjusted
    double btagw = 1;
    if(isMC and (sample == Sample::topmu or sample == Sample::topel))
      btagw = 0.920;
    if(isMC and (sample == Sample::sig and category != Category::VBF and category != Category::VBFrelaxed))
      btagw = 1.015;
    if(isMC and (sample == Sample::zee and category != Category::VBF and category != Category::VBFrelaxed))
      btagw = 0.980;
    else
      btagw = *wgtbtag;
    
    //V-tagging scale factor --> only for mono-V
    if(isMC && category == Category::monoV && isWJet)
      sfwgt *= getVtaggingScaleFactor(tau2tau1,sysName);
    
    //Gen level info --> NLO re-weight    
    Double_t kwgt = 1.0;    
    double genpt = *wzpt;
    float sf_ewk = 1.0;
    float sf_qcd_mj  = 1.0;
    float sf_qcd_vbf = 1.0;
    for (size_t i = 0; i < khists.size(); i++) {
      if (khists[i]) {
	if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
	kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));

	if(i == 1)
	  sf_ewk = khists[i]->GetBinContent(khists[i]->FindBin(genpt));
	else if(i == 0)
	  sf_qcd_mj *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
	else if(i == 2)
	  sf_qcd_vbf *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
      }
    }

    //Gen level info --> kfactor NLO for V-EWK processes
    Double_t kewkgt = 1.0;
    if((category == Category::VBF or category == Category::VBFrelaxed or category == Category::twojet) and isMC){
      double genpt = *wzpt;
      if(jetpt->size() >= 2) {
	TLorentzVector jet1 ;
	TLorentzVector jet2 ;
	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	double mjj = (jet1+jet2).M();
	for(size_t i = 0; i < kVEWKhists.size(); i++){
	  if(kVEWKhists[i]){// good histogram
	    if(genpt <= kVEWKhists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = kVEWKhists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	    if(genpt >= kVEWKhists[i]->GetXaxis()->GetBinLowEdge(kVEWKhists[i]->GetNbinsX()+1)) genpt = kVEWKhists[i]->GetXaxis()->GetBinLowEdge(kVEWKhists[i]->GetNbinsX()+1)-1;
	    if(mjj <= kVEWKhists[i]->GetYaxis()->GetBinLowEdge(1)) mjj = kVEWKhists[i]->GetYaxis()->GetBinLowEdge(1) + 1;
	    if(mjj >= kVEWKhists[i]->GetYaxis()->GetBinLowEdge(kVEWKhists[i]->GetNbinsY()+1)) mjj = kVEWKhists[i]->GetYaxis()->GetBinLowEdge(kVEWKhists[i]->GetNbinsY()+1)-1;
	    kewkgt *= kVEWKhists[i]->GetBinContent(kVEWKhists[i]->FindBin(genpt,mjj));
	  }
	}
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

    // NNLO corrections to Higgs pT
    Double_t hnnlowgt = 1.0;
    if(isHiggsInvisible and higgsNNLO and isMC){
      if(*dmpt < higgsNNLO->GetBinLowEdge(1))
	*dmpt = higgsNNLO->GetBinLowEdge(1)+1;
      else if(*dmpt > higgsNNLO->GetBinLowEdge(higgsNNLO->GetNbinsX()+1))
	*dmpt = higgsNNLO->GetBinLowEdge(higgsNNLO->GetNbinsX()+1)-1;
      hnnlowgt *= higgsNNLO->GetBinContent(higgsNNLO->FindBin(*dmpt));      
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
    // boosted category and monojet
    bool goodMonoJet = false;
    bool goodMonoV   = false;
    bool goodVBF     = false;

    if(sample == Sample::taun and leadingjet_notau == -1) continue;
    else if(sample == Sample::taun and leadingjet_notau != -1) leadingCentralJetPos = leadingjet_notau;
    
    if(category == Category::inclusive){ 
      if (centralJets.size() == 0) continue;
      if (fabs(jeteta->at(leadingCentralJetPos)) > 2.5) continue;
      if (chfrac->at(leadingCentralJetPos) < 0.1) continue;   // jet id
      if (nhfrac->at(leadingCentralJetPos) > 0.8) continue;   // jet id
      if (jetpt->at(leadingCentralJetPos)  < leadingJetPtCut) continue;  // jet1 > 100 GeV
      if (sample != Sample::qcd and jmdphi < 0.5) continue; // deltaPhi cut
      else if (sample == Sample::qcd and jmdphi > 0.5) continue; // deltaPhi cut
    }

    else if(category == Category::monojet){ // mono jet + V-jet veto
      if (centralJets.size() == 0) continue;
      if (centralJets.size() < njetsMin) continue;
      if (centralJets.size() > njetsMax) continue;
      if (fabs(jeteta->at(leadingCentralJetPos)) > 2.5) continue;
      if (chfrac->at(leadingCentralJetPos) < 0.1) continue;   // jet id                                                                                                   
      if (nhfrac->at(leadingCentralJetPos) > 0.8) continue;   // jet id                                                                                                   
      if (jetpt->at(leadingCentralJetPos)  < leadingJetPtCut) continue;  // jet1 > 100 GeV                                                                                          
      if (sample != Sample::qcd and jmdphi < 0.5) continue; 
      else if(sample == Sample::qcd and jmdphi > 0.5) continue;	
      goodMonoJet = true;
      
    
      // VBF selection
      if(goodMonoJet and centralJets.size()+forwardJets.size() > 2 and fabs(jeteta->at(0)) < 4.7 and fabs(jeteta->at(1)) < 4.7 and 
	 jetpt->at(0) > leadingJetPtCutVBF and jetpt->at(1) > trailingJetPtCutVBF and
	 jmdphi > jetmetdphiVBF and jeteta->at(0)*jeteta->at(1) < 0 and
	 fabs(jeteta->at(0)-jeteta->at(1)) > detajj){
	
	TLorentzVector jet1 ;
	TLorentzVector jet2 ;
	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	if((jet1+jet2).M() > mjj and fabs(deltaPhi(jetphi->at(0),jetphi->at(1))) < dphijj){
	  if(removeVBF) goodMonoJet = false;
	}
      }
      
      // V-tagging
      if(goodMonoJet and boostedJetpt->size() != 0 and 
	 fabs(boostedJeteta->at(0)) < jetEtaAK8 and
	 boostedJetpt->at(0) > ptJetMinAK8 and
	 prunedJetm->at(0) > prunedMassMin and prunedJetm->at(0) < prunedMassMax and
	 boostedJettau2->at(0)/boostedJettau1->at(0) < tau2tau1)
      goodMonoJet = false;
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


      // remove VBF overlap
      if(goodMonoV and category == Category::monoV and centralJets.size()+forwardJets.size() > 2 and fabs(jeteta->at(0)) < 4.7 and fabs(jeteta->at(1)) < 4.7 and
         jetpt->at(0) > leadingJetPtCutVBF and jetpt->at(1) > trailingJetPtCutVBF and
         (sample != Sample::taun and jmdphi > jetmetdphiVBF) and jeteta->at(0)*jeteta->at(1) < 0 and
         fabs(jeteta->at(0)-jeteta->at(1)) > detajj){
	
        TLorentzVector jet1 ;
        TLorentzVector jet2 ;
        jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
        jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
        if((jet1+jet2).M() > mjj and fabs(deltaPhi(jetphi->at(0),jetphi->at(1))) < dphijj){
          if(removeVBF) goodMonoV = false;
	}
      }             
    }

    ///////////////
    else if(category == Category::VBF){
      if(centralJets.size()+forwardJets.size() < 2) continue;
      if(fabs(jeteta->at(0)) > 4.7 or fabs(jeteta->at(1)) > 4.7) continue;
      if(jetpt->at(0) < leadingJetPtCutVBF)  continue;
      if(jetpt->at(1) < trailingJetPtCutVBF) continue;

      if (sample != Sample::qcd and jmdphi < jetmetdphiVBF) continue;
      else if(sample == Sample::qcd and jmdphi > jetmetdphiVBF) continue;

      if(fabs(jeteta->at(0)) < 2.4 and chfrac->at(0) < 0.1) continue;
      if(fabs(jeteta->at(0)) < 2.4 and nhfrac->at(0) > 0.8) continue;

      if(jeteta->at(0)*jeteta->at(1) > 0 ) continue;
      if(fabs(jeteta->at(0)-jeteta->at(1)) < detajj) continue;
      //if(fabs(jeteta->at(0)) >= 3.0 and fabs(jeteta->at(0)) <= 3.2 and nhfrac->at(0) > 0.96) continue;
      TLorentzVector jet1 ;
      TLorentzVector jet2 ;
      jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
      jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
      if((jet1+jet2).M() < mjj) continue;
      if(fabs(deltaPhi(jetphi->at(0),jetphi->at(1))) > dphijj) continue;
      goodVBF = true;
    }

    ///////////////
    else if(category == Category::VBFrelaxed){
      if(centralJets.size()+forwardJets.size() < 2) continue;
      if(fabs(jeteta->at(0)) > 4.7 or fabs(jeteta->at(1)) > 4.7) continue;
      if(jetpt->at(0) < leadingJetPtCutVBF) continue;
      if(jetpt->at(1) < trailingJetPtCutVBF) continue;
      
      if (sample != Sample::qcd and jmdphi < jetmetdphiVBF) continue;
      else if(sample == Sample::qcd and jmdphi > jetmetdphiVBF) continue;
      
      if(fabs(jeteta->at(0)) < 2.4 and chfrac->at(0) < 0.1) continue;
      if(fabs(jeteta->at(0)) < 2.4 and nhfrac->at(0) > 0.8) continue;
      if(jeteta->at(0)*jeteta->at(1) > 0 ) continue;
      if(fabs(jeteta->at(0)-jeteta->at(1)) < detajjrelaxed) continue;
      //if(fabs(jeteta->at(0)) >= 3.0 and fabs(jeteta->at(0)) <= 3.2 and nhfrac->at(0) > 0.96) continue;
      TLorentzVector jet1 ;
      TLorentzVector jet2 ;
      jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
      jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
      if((jet1+jet2).M() < mjjrelaxed) continue;
      if(fabs(deltaPhi(jetphi->at(0),jetphi->at(1))) > dphijjrelaxed) continue;

      goodVBF = true;
    }
    else if(category == Category::twojet){
      if(centralJets.size()+forwardJets.size() < 2) continue;
      if(fabs(jeteta->at(0)) > 4.7 or fabs(jeteta->at(1)) > 4.7) continue;
      if(jetpt->at(0) < leadingJetPtCutVBF) continue;
      if(jetpt->at(1) < trailingJetPtCutVBF) continue;
      if(fabs(jeteta->at(0)) < 2.4 and chfrac->at(0) < 0.1) continue;
      if(fabs(jeteta->at(0)) < 2.4 and nhfrac->at(0) > 0.8) continue;
      //if(fabs(jeteta->at(0)) >= 3.0 and fabs(jeteta->at(0)) <= 3.2 and nhfrac->at(0) > 0.96) continue;
      if (sample != Sample::qcd and jmdphi < 0.5) continue;
      else if(sample == Sample::qcd and jmdphi > 0.5) continue;	      
    }
  
    
    //////// Remove Category Overlaps    
    if(category == Category::monoV and goodMonoV == false) continue;
    if(category == Category::monojet and goodMonoJet == false) continue;
    

    // fill 1D histogram
    // fill the histograms --> with the right observable
    for(auto hist : hist1D){
      double fillvar = -99;
      TString name(hist->GetName());      
      name.ReplaceAll("metJet","");
      name.ReplaceAll("metRes","");
      name.ReplaceAll("metUnc","");

      // number of vertices
      if(name.Contains("nvtx")) 
	fillvar = *nvtx;
      // number of vertices
      else if(name.Contains("ntaus")) 
	fillvar = *ntaus;
      // Z-mass peak to mumu
      else if(name.Contains("zmass")) 
	fillvar = *zmass;
      // Z-mass peak to ee
      else if(name.Contains("zeemass")) 
	fillvar = *zeemass;
      else if(name.Contains("dRmumu")){
	float dphi = fabs(*mu1phi-*mu2phi);
	if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
	fillvar = sqrt(dphi*dphi+fabs(*mu1eta-*mu2eta)*fabs(*mu1eta-*mu2eta));
      }
      else if(name.Contains("dRee")){
	float dphi = fabs(*el1phi-*el2phi);
	if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
	fillvar = sqrt(dphi*dphi+fabs(*el1eta-*el2eta)*fabs(*el1eta-*el2eta));
      }
      // met-phi with type-I correction 
      else if(name.Contains("t1pfmetphi"))
	fillvar = *metphi;      
      // met with type-I correction 
      else if(name.Contains("t1pfmet")) 
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
      //Delta phi met leptons and boson
      else if(name.Contains("dphi_jet_boson")){
	float mindphi = TMath::Pi();
	for(int ijet = 0; ijet < jetpt->size(); ijet++){
	  if(fabs(jetpt->at(ijet)) < 30) continue;
	  if(deltaPhi(bosonPhi,jetphi->at(ijet)) < mindphi)
	    mindphi = deltaPhi(bosonPhi,jetphi->at(ijet));
	}	
	fillvar = mindphi;
      }
      // Lepton Pt and Eta
      else if(name.Contains("elp1pt")){
	if(pid1 > 0) // leading e+ pt
	fillvar = pt1;      
      }
      else if(name.Contains("elm1pt")){
	if(pid1 < 0) // leading e- pt	  
	fillvar = pt1;      
      }
      else if(name.Contains("el1pt"))	
	fillvar = pt1;      
      else if(name.Contains("tau1pt"))	
	fillvar = pt1;      
      else if(name.Contains("mup1pt")){
	if(pid1 > 0) // mu+ pt
	fillvar = pt1;
      }
      else if(name.Contains("mum1pt")){
	if(pid1 < 0) // mu- pt
	  fillvar = pt1;
      }
      else if(name.Contains("mu1pt")) // mu pt
	fillvar = pt1;
      else if(name.Contains("mup2pt")){
	if(pid2 > 0)
	  fillvar = pt2;
      }
      else if(name.Contains("mum2pt")){
	if(pid2 < 0)
	  fillvar = pt2;
      }
      else if(name.Contains("mu2pt"))
	fillvar = pt2;
      else if(name.Contains("elp2pt")){
	if(pid2 > 0)
	  fillvar = pt2;
      }
      else if(name.Contains("elm2pt")){
	if(pid2 < 0)
	  fillvar = pt2;
      }
      else if(name.Contains("el2pt"))
	fillvar = pt2;
      else if(name.Contains("mup1eta")){
	if(pid1 > 0)
	  fillvar = eta1;      
      }
      else if(name.Contains("mum1eta")){
	if(pid1 < 0)
          fillvar = eta1;
      }
      else if(name.Contains("tau1eta"))
	fillvar = eta1;      
      else if(name.Contains("mu1eta"))
	fillvar = eta1;      
      else if(name.Contains("mup2eta")){
	if(pid2 > 0)
	  fillvar = eta2;
      } 
      else if(name.Contains("mum2eta")){
	if(pid2 < 0)
	  fillvar = eta2;
      } 
      else if(name.Contains("mu2eta"))
	fillvar = eta2;
      else if(name.Contains("elp1eta")){
	if(pid1 > 0)
	  fillvar = eta1;      
      }
      else if(name.Contains("elm1eta")){
	if(pid1 < 0)
	  fillvar = eta1;      
      }
      else if(name.Contains("el1eta"))
	fillvar = eta1;
      else if(name.Contains("elp2eta")){
	if(pid1 > 0)
	  fillvar = eta2;
      }
      else if(name.Contains("elm2eta")){
	if(pid2 < 0)
	  fillvar = eta2;
      }
      else if(name.Contains("el2eta"))
	fillvar = eta2;
      // W-boson properties
      else if(name.Contains("wmt"))
	fillvar = *wmt;
      else if(name.Contains("wemt"))
	fillvar = *wemt;
      // jet neutral hadron fractions: leading jets in hf, he, hb or overall
      else if(name.Contains("jetnhfrac_hf")){	
	if(fabs(jeteta->at(0)) >= 3)
	  fillvar = nhfrac->at(0);	  	
      }
      else if(name.Contains("jetnhfrac_he")){
	if(fabs(jeteta->at(0)) >= 2.5 and fabs(jeteta->at(0)) <= 3)
	  fillvar = nhfrac->at(0);
      }
      else if(name.Contains("jetnhfrac_hb")){
	if(fabs(jeteta->at(0)) <= 2.5)
	  fillvar = nhfrac->at(0);
      }
      else if(name.Contains("jetnhfrac")){
	if((category == Category::VBF or category == Category::twojet or category == Category::VBFrelaxed))
	  fillvar = nhfrac->at(0);
	else
	  fillvar = nhfrac->at(leadingCentralJetPos);
      }
      else if(name.Contains("jet2nhfrac_hf")){
	if(jetpt->size() > 1 and fabs(jeteta->at(1)) >= 3)
	  fillvar = nhfrac->at(1);
      }
      else if(name.Contains("jet2nhfrac_he")){
	if(jetpt->size() > 1 and fabs(jeteta->at(1)) >= 2.5 and fabs(jeteta->at(1)) <= 3)
	  fillvar = nhfrac->at(1);
      }
      else if(name.Contains("jet2nhfrac_hb")){
	if( jetpt->size() > 1 and  fabs(jeteta->at(1)) <= 2.5)
	  fillvar = nhfrac->at(1);
      }
      else if(name.Contains("jet2nhfrac")){
	if(jetpt->size() > 1)
	  fillvar = nhfrac->at(1);
      }
      // jet EM fractions: leading jets in hf, he, hb or overall
      else if(name.Contains("jetemfrac_hf")){
	if(fabs(jeteta->at(0)) >= 3)
	  fillvar = emfrac->at(0);	  
      }
      else if(name.Contains("jetemfrac_he")){
	if(fabs(jeteta->at(0)) >= 2.5 and fabs(jeteta->at(0)) <= 3)
	  fillvar = emfrac->at(0);
      }
      else if(name.Contains("jetemfrac_hb")){
	if(fabs(jeteta->at(0)) <= 2.5)
	  fillvar = emfrac->at(0);
      }
      else if(name.Contains("jetemfrac")){
	if(category == Category::VBF or category == Category::twojet or category == Category::VBFrelaxed)
	fillvar = emfrac->at(0);
	else
	  fillvar = emfrac->at(leadingCentralJetPos);
      }
      else if(name.Contains("jet2emfrac_hf")){
	if(jetpt->size() > 1 and fabs(jeteta->at(1)) >= 3)
	  fillvar = emfrac->at(1);
      }
      else if(name.Contains("jet2emfrac_he")){
	if(jetpt->size() > 1 and fabs(jeteta->at(1)) >= 2.5 and fabs(jeteta->at(1)) <= 3)
	  fillvar = emfrac->at(1);
      }
      else if(name.Contains("jet2emfrac_hb")){
	if(jetpt->size() > 1 and  fabs(jeteta->at(1)) <= 2.5)
	  fillvar = emfrac->at(1);
      }
      else if(name.Contains("jet2emfrac")){
	if(jetpt->size() > 1)
	  fillvar = emfrac->at(1);
      }
      // jet CH fractions: leading jets in hf, he, hb or overall
      else if(name.Contains("jetchfrac_hf")){
	if(fabs(jeteta->at(0)) >= 3)
	fillvar = chfrac->at(0);
      }	  
      else if(name.Contains("jetchfrac_he")){
	if(fabs(jeteta->at(0)) >= 2.5 and fabs(jeteta->at(0)) <= 3)
	  fillvar = chfrac->at(0);
      }
      else if(name.Contains("jetchfrac_hb")){
	if(fabs(jeteta->at(0)) <= 2.5)
	  fillvar = chfrac->at(0);
      }
      else if(name.Contains("jetchfrac")){
	if(category == Category::VBF or category == Category::twojet or category == Category::VBFrelaxed)
	  fillvar = chfrac->at(0);
	else 
	  fillvar = chfrac->at(leadingCentralJetPos);
      }
      else if(name.Contains("jet2chfrac_hf")){
	if(jetpt->size() > 1 and fabs(jeteta->at(1)) >= 3)
	  fillvar = chfrac->at(1);
      }
      else if(name.Contains("jet2chfrac_he")){
	if(jetpt->size() > 1 and fabs(jeteta->at(1)) >= 2.5 and fabs(jeteta->at(1)) <= 3)
	  fillvar = chfrac->at(1);
      }
      else if(name.Contains("jet2chfrac_hb")){
	if(jetpt->size() > 1 and  fabs(jeteta->at(1)) <= 2.5)
	  fillvar = chfrac->at(1);
      }
      else if(name.Contains("jet2chfrac")){
	if(jetpt->size() > 1)
	  fillvar = chfrac->at(1);
      }
      // AK8 jet kinematics
      else if(name.Contains("boostedjetpt")){
	if(boostedJetpt->size() > 0)
	  fillvar = boostedJetpt->at(0);
      }
      else if(name.Contains("boostedjeteta")){
	if(boostedJeteta->size() > 0)
	  fillvar = boostedJeteta->at(0);
      }
      // Jet met Dphi used for selection
      else if(name.Contains("jetmetdphi"))
	fillvar = jmdphi;
      // calo met
      else if(name.Contains("calomet"))
	fillvar = fabs(*metcalo-*met)/pfmet;
      else if(name.Contains("met"))
	fillvar = pfmet;             
      else if(name.Contains("jetpt2")){
	if(jetpt->size() >= 2)
	  fillvar = jetpt->at(1);
      }
      else if(name.Contains("jeteta2")){
	if(jetpt->size() >= 2)
	  fillvar = jeteta->at(1);
      }
      else if(name.Contains("jeteta2jeteta1")){
	if(jetpt->size() >= 2)
	  fillvar = jeteta->at(1)*jeteta->at(0);
      }
      else if(name.Contains("jetpt")){
	if(category == Category::VBF or category == Category::twojet or category == Category::VBFrelaxed)
	  fillvar = jetpt->at(0);
	else
	  fillvar = jetpt->at(leadingCentralJetPos);
      }
      else if(name.Contains("jeteta")){
	if(category == Category::VBF or category == Category::twojet or category == Category::VBFrelaxed)
	  fillvar = jeteta->at(0);
	else
	  fillvar = jeteta->at(leadingCentralJetPos);
      }
      else if(name.Contains("detajj")){
	if(jetpt->size() >= 2)
	  fillvar = fabs(jeteta->at(0)-jeteta->at(1));
      }
      else if(name.Contains("mjj")){
	if(jetpt->size() >= 2){
	  TLorentzVector jet1 ;
	  TLorentzVector jet2 ;
	  jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	  jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	  fillvar = (jet1+jet2).M();
	}
      }
      else if(name.Contains("mT")){
	float deltaphi = deltaPhi(jetphi->at(0),pfmetphi);
	fillvar = sqrt(2*jetpt->at(0)*pfmet*(1-cos(deltaphi)));
      }   
      else if(name.Contains("ncjet"))
	fillvar = *njets;      
      else if(name.Contains("njet"))
	fillvar = *nincjets;      
      else if(name.Contains("nbjet"))
	fillvar = *nbjets;
      else if(name.Contains("bosonpt"))
	fillvar = bosonPt;    	
      else if(name.Contains("bosoneta"))
	fillvar = bosonEta;    	
      else if(name.Contains("bosonphi"))
	fillvar = bosonPhi;    	
      else if(name.Contains("jetQGL"))
	fillvar = jetQGL->at(0);
      else if(name.Contains("jet2QGL")){
	if(jetQGL->size() > 1)
	  fillvar = jetQGL->at(1);
      }
      else if(name.Contains("jetPFlav")){
	if(pFlav->size() > 0)
	  fillvar = fabs(pFlav->at(0));
      }
      else if(name.Contains("jetQGL_AK8")){
	if(boostedJetpt->size() > 0)
	  fillvar = boostedJetQGL->at(0);
      }
      // substructure
      else if(name.Contains("mpruned")){
	if(prunedJetm->size() > 0 and boostedJetpt->at(0) > ptJetMinAK8)	  
	  fillvar = prunedJetm->at(0);	
      }
      else if(name.Contains("HT"))
	fillvar = *ht;      
      else if(name.Contains("tau2tau1")){
	if(boostedJettau1->size() > 0 and boostedJettau2->size() > 0 and boostedJetpt->at(0) > ptJetMinAK8)
	  fillvar = boostedJettau2->at(0)/boostedJettau1->at(0);
	//	fillvar = boostedJettau2->at(0)/boostedJettau1->at(0)+0.063*log(pow(prunedJetm->at(0),2)/boostedJetpt->at(0));
      }
      // b-tagging
      else if(name.Contains("btagCSV_max")){
	float btagMax = -10.;
	for(size_t iBjet = 0; iBjet < jetbtag->size(); iBjet++){
	  if(jeteta->at(iBjet) > 2.4 or jetpt->at(iBjet) < 20) continue;
	  if(jetbtag->at(iBjet) >= btagMax)
	    btagMax = jetbtag->at(iBjet);
	}
	if(fabs(btagMax) != 10.)
	  fillvar = btagMax;
      }
      else if(name.Contains("btagCSV_min")){
	float btagMin = 10.;
	for(size_t iBjet = 0; iBjet < jetbtag->size(); iBjet++){
	  if(jeteta->at(iBjet) >= 2.4 or jetpt->at(iBjet) < 20) continue;
	  if(jetbtag->at(iBjet) < btagMin and jetbtag->at(iBjet) > 0)
	    btagMin = jetbtag->at(iBjet);
	}
	if(fabs(btagMin) != 10.)
	  fillvar = btagMin;
      }
      else if(name.Contains("btagCSV")){
	if(jetbtag->size() > 0)
	  fillvar = jetbtag->at(0);
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

	if(fabs(*wgtpu) > 100)  puwgt = 1;
	else if(fabs(*wgtpu) < 0.01) puwgt = 1;
	else puwgt = *wgtpu;

	cout<<"event "<<*event<<" lumi "<<*lumisection<<" ewk "<<sf_ewk<<" qcd_mj "<<sf_qcd_mj<<" qcd_vbf "<<sf_qcd_vbf<<" wzpt "<<*wzpt<<endl;

	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*kewkgt*hwgt*hnnlowgt*pfwgt/(**wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*kewkgt*hwgt*hnnlowgt*pfwgt/(**wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)
      }
      else if (isMC and reweightNVTX){

	// pu-weight
	if (*nvtx <= 60 and not isSummer16) 
	  puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
	else if(isSummer16 and sample != Sample::sig and sample != Sample::gam and (category == Category::monojet or category == Category::monoV))
	  puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
	else
	  puwgt = 1;
	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*topptwgt*sfwgt*kwgt*kewkgt*hwgt*ggZHwgt*hnnlowgt*pfwgt/(**wgtsum);
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*topptwgt*sfwgt*kwgt*kewkgt*hwgt*ggZHwgt*hnnlowgt*pfwgt/(**wgtsum);
      }
      
      // for data-based events 
      if (!isMC && sample == Sample::qcdgam) 
	evtwgt = sfwgt*pfwgt*hltw;
      else if (!isMC)
	evtwgt = hltw;      
      
      hist->Fill(fillvar, evtwgt);
    }

    //fill 2D histo
    for(auto hist: hist2D){
      double fillvarX = 0;
      double fillvarY = 0;
      TString name(hist->GetName());
      name.ReplaceAll("metJet","");
      name.ReplaceAll("metRes","");
      name.ReplaceAll("metUnc","");

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
      else if(name.Contains("met_mpruned")){
	if(boostedJetpt->size() > 0 and boostedJetpt->at(0) > ptJetMinAK8 ){
	  fillvarX = pfmet;
	  if(prunedJetm->size() > 0)
	    fillvarY = prunedJetm->at(0);	
	}
      }
      else if(name.Contains("met_tau2tau1")){
	if(boostedJetpt->size() > 0 and boostedJetpt->at(0) > ptJetMinAK8){
	  fillvarX = pfmet;
	  if(boostedJettau2->size() > 0 and boostedJettau1->size() >0)
	    fillvarY = boostedJettau2->at(0)/boostedJettau1->at(0);		
	}
      }
      else if(name.Contains("mpruned_tau2tau1")){
	if(boostedJetpt->size() > 0 and boostedJetpt->at(0) > ptJetMinAK8){
	  if(prunedJetm->size() > 0)
	    fillvarX = prunedJetm->at(0);
	  if(boostedJettau2->size() > 0 and boostedJettau1->size() >0)
	    fillvarY = boostedJettau2->at(0)/boostedJettau1->at(0);			
	}
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
      Double_t puwgt = 1.0;

      if (isMC and not reweightNVTX){

	if(fabs(*wgtpu) > 100)  puwgt = 1;
	else if(fabs(*wgtpu) < 0.01) puwgt = 1;
	else puwgt = *wgtpu;

	if(XSEC != -1)
          evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*topptwgt*sfwgt*kwgt*kewkgt*hwgt*ggZHwgt*hnnlowgt*pfwgt/(**wgtsum);
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*topptwgt*sfwgt*kwgt*kewkgt*hwgt*ggZHwgt*hnnlowgt*pfwgt/(**wgtsum);
      }
      else if (isMC and reweightNVTX){
        // pu-weight                                                                                                                                                                                  
        if (*nvtx <= 60 and not isSummer16)
          puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
        else if(isSummer16 and sample != Sample::sig and sample  != Sample::gam and (category == Category::monojet or category == Category::monoV))
          puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
        else if(isSummer16 and (category == Category::VBFrelaxed or category == Category::twojet or category == Category::VBF))
          puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
	else
          puwgt = 1;
	
	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*kewkgt*hwgt*hnnlowgt/(**wgtsum);
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*kewkgt*hwgt*hnnlowgt/(**wgtsum);	
      }

      if (!isMC && sample == Sample::qcdgam) 
	evtwgt = sfwgt*hltw;
      else if (!isMC)
	evtwgt = hltw;
      
	hist->Fill(fillvarX,fillvarY,evtwgt);            
     }
  }

  if(pufile != NULL)
    pufile->Close();
  if(sffile_eleTight != NULL)
    sffile_eleTight->Close();
  if(sffile_eleVeto != NULL)
    sffile_eleVeto->Close();
  if(sffile_muTight != NULL)
    sffile_muTight->Close();
  if(sffile_muLoose != NULL)
    sffile_muLoose->Close();
  if(sffile_phoMedium != NULL)
    sffile_phoMedium->Close();
  if(purityfile_photon != NULL)
    purityfile_photon->Close();
  if(triggerfile_SinglEle != NULL)
    triggerfile_SinglEle->Close();
  if(triggerfile_SinglEle_jetHT != NULL)
    triggerfile_SinglEle_jetHT->Close();
  if(triggerfile_MET != NULL)
    triggerfile_MET->Close();
  if(triggerfile_MET_zmm != NULL)
    triggerfile_MET_zmm->Close();
  if(triggerfile_SinglePhoton != NULL)
    triggerfile_SinglePhoton->Close();
  if(triggerfile_SinglePhoton_jetHT != NULL)
    triggerfile_SinglePhoton_jetHT->Close();
  if(trackingefficiency_muon != NULL)
    trackingefficiency_muon->Close();
  if(trackingefficiencyFile_electron != NULL)
    trackingefficiencyFile_electron->Close();
  if(postFitFile != NULL)
    postFitFile->Close();
  for(auto file : triggerfile_MET_binned_Wmn)
    file->Close();
  for(auto file : triggerfile_MET_binned_Zmm)
    file->Close();
}
#endif
