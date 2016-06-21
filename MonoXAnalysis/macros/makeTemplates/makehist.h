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
const float ptJetMinAK8     = 250.;
const float jetEtaAK8       = 2.4;
const float pfMetMonoVLower = 250.;
const float pfMetMonoVUpper = 8000.;
const int   vBosonCharge    = 0;
const int   nBjets          = 1;
const bool  reweightNVTX    = true;

string kfactorFile       = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root";
string kfactorFileUnc    = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/scalefactors_v4.root";
string baseInputTreePath = "/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_10_06_2016_v2/";

VectorSorter jetSorter;

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
      sfwgt *= (0.692+0.144);
    else if(sysName == "VtagDown")
      sfwgt *= (0.692-0.144);
    else
      sfwgt *= 0.692;
  }
  else if(tau2tau1 == 0.6){
    if(sysName == "VtagUp")
      sfwgt *= (1.031+0.129);
    else if(sysName == "VtagDown")
      sfwgt *= (1.031-0.129);
    else
      sfwgt *= 1.031;
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
	       TH1* hhist = NULL,
	       TH2* ggZHhist = NULL,
	       const bool   & is76Xsample = false
	       ) {

  if(not tree){
    cout<<" empty tree --> skip process "<<endl;
    return;
  }

  cout<<"make resonant selection "<<resonantSelection<<endl;

  // in case you want to weight the NVTX distribution
  TFile* pufile = NULL;
  TH1* puhist = NULL;
  if(not is76Xsample){
    pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt.root");
    puhist = (TH1*) pufile->Get("puhist");
  }
  else{
    pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt_76X.root");
    puhist = (TH1*) pufile->Get("puhist");
  }    
    
  // electron and muon ID scale factor files                                                                                                                                    
  TFile sffile_eleTight("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/scaleFactor_electron_tightid.root");
  TFile sffile_eleVeto("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/scaleFactor_electron_vetoid.root");
  TFile sffile_muTight("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/scaleFactor_muon_tightid.root");
  TFile sffile_muLoose("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/scaleFactor_muon_looseid.root");

  TH2*  msfloose = (TH2*)sffile_muLoose.Get("scaleFactor_muon_looseid_RooCMSShape");
  TH2*  msftight = (TH2*)sffile_muTight.Get("scaleFactor_muon_tightid_RooCMSShape");
  TH2*  esfveto  = (TH2*)sffile_eleVeto.Get("scaleFactor_electron_vetoid_RooCMSShape");
  TH2*  esftight = (TH2*)sffile_eleTight.Get("scaleFactor_electron_tightid_RooCMSShape");

  // Photon ID scale factor                                                                                                                                                     
  TFile sffile_phoLoose("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF_2016/scaleFactor_photon_looseid.root");
  TFile sffile_phoMedium("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF_2016/scaleFactor_photon_mediumid.root");
  TH2*  psfloose  = (TH2*)sffile_phoLoose.Get("scaleFactor_photon_looseid_RooCMSShape");
  TH2*  psfmedium = (TH2*)sffile_phoMedium.Get("scaleFactor_photon_mediumid_RooCMSShape");

  // Photon Purity                                                                                                                                                              
  TFile purityfile_photon ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/PhotonSFandEffandPurity_Lumi2p1fb_0202.root");
  TH2*  purhist = (TH2*) purityfile_photon.Get("PhotonPurity");

  // trigger files used for 2016                                                                                                                                                
  TFile triggerfile_SinglEle("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/triggerEfficiency_DATA_SingleElectron.root");
  TFile triggerfile_SingleMu("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/triggerEfficiency_DATA_SingleMuon.root");
  TFile triggerfile_MET("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/metTriggerEfficiency.root");
  TFile triggerfile_SinglePhoton("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/photonTriggerEfficiency.root");

  TH2*  triggerelhist = (TH2*) triggerfile_SinglEle.Get("trigeff_ele27wptight");
  TH2*  triggermuhist = (TH2*) triggerfile_SingleMu.Get("trigeff_muIso");
  TF1*  triggermet = (TF1*) triggerfile_MET.Get("efficiency_func");
  TEfficiency*  triggerphoton = (TEfficiency*)triggerfile_SinglePhoton.Get("efficiency");
  TGraphAsymmErrors* triggerphoton_graph = triggerphoton->CreateGraph();


  // Post-fit weights
  TFile* postFitFile = NULL;
  if(sample == Sample::sig and applyPostFitWeight and category == Category::monojet)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ/postfit_weights_Sig.root");
  if(sample == Sample::zmm and applyPostFitWeight and category == Category::monojet)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ/postfit_weights_ZM.root");
  else if(sample == Sample::wmn and applyPostFitWeight and category == Category::monojet)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ/postfit_weights_WM.root");
  else if(sample == Sample::zee and applyPostFitWeight and category == Category::monojet)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ/postfit_weights_ZE.root");
  else if(sample == Sample::wen and applyPostFitWeight and category == Category::monojet)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ/postfit_weights_WE.root");
  else if(sample == Sample::qcd and applyPostFitWeight and category == Category::monojet)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ/postfit_weights_GJ.root");
  else if(sample == Sample::gam and applyPostFitWeight and category == Category::monojet)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoJ/postfit_weights_GJ.root");
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
  else if(sample == Sample::qcd and applyPostFitWeight and category == Category::monoV)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoV/postfit_weights_GJ.root");
  else if(sample == Sample::gam and applyPostFitWeight and category == Category::monoV)
    postFitFile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/postFitOverPrefit/monoV/postfit_weights_GJ.root");

  TH1* postFitWeight = NULL;
  if(postFitFile != NULL)
    postFitWeight = (TH1*) postFitFile->FindObjectAny("postfit_over_prefit");

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
  TTreeReaderValue<unsigned int> event  (myReader,"event");
  TTreeReaderValue<unsigned int> nvtx   (myReader,"nvtx");
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
      btagname = "wgtbtagUp";
    else if(sysName == "btagDown" or sysName == "btagDw" or sysName == "bDown" or sysName == "bDw")
      btagname = "wgtbtagDown";
    else
      btagname = "wgtbtag";
    prescalename  = "pswgt_ph120";
    hltphotonname = "hltphoton120";
  }
  else if(!isMC){ // in case of data
    wgtname       = "wgt";
    btagname      = "wgt";
    wgtpileupname = "wgt";
    prescalename  = "pswgt_ph120";
    hltphotonname = "hltphoton120";    
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
  TTreeReaderValue<UChar_t> hltm       (myReader,"hltsinglemu");
  TTreeReaderValue<UChar_t> hltp120    (myReader, hltphotonname.c_str());
  TTreeReaderValue<UChar_t> hltp165    (myReader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175    (myReader,"hltphoton175");
  // prescale
  TTreeReaderValue<double> pswgt(myReader,prescalename.c_str());

  TTreeReaderValue<UChar_t> fhbhe  (myReader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (myReader,"flaghbheiso");

  string cscfilter = "flagglobaltighthalo";
  if(isMC) cscfilter = "flagcsctight";
  TTreeReaderValue<UChar_t> fcsc   (myReader,cscfilter.c_str());
  TTreeReaderValue<UChar_t> feeb   (myReader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (myReader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (myReader,"flaggoodvertices");

  TTreeReaderValue<unsigned int> njets    (myReader,"njets");
  TTreeReaderValue<unsigned int> ntaus    (myReader,"ntaus");
  TTreeReaderValue<unsigned int> ntausraw    (myReader,"ntausraw");
  TTreeReaderValue<unsigned int> nincjets (myReader,"njetsinc");
  TTreeReaderValue<unsigned int> nbjets   (myReader,"nbjetslowpt");
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

  // AK8 jet
  TTreeReaderValue<vector<double> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<double> > boostedJetQGL   (myReader,"boostedJetQGL");
  TTreeReaderValue<vector<double> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<double> > boostedJetphi   (myReader,"boostedJetphi");
  TTreeReaderValue<vector<double> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<double> > prunedJetm      (myReader,"prunedJetm_v2");
  TTreeReaderValue<vector<double> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<double> > boostedJettau1  (myReader,"boostedJettau1");
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
  
  TTreeReaderValue<double> wzpt   (myReader,"wzpt");
  TTreeReaderValue<double> wzpt_h (myReader,"wzpt_h");
  TTreeReaderValue<double> wzeta  (myReader,"wzeta");
  TTreeReaderValue<double> zmass  (myReader,"zmass");
  TTreeReaderValue<double> zmmpt  (myReader,"zpt");
  TTreeReaderValue<double> zeept  (myReader,"zeept");
  TTreeReaderValue<double> zeeeta (myReader,"zeeeta");
  TTreeReaderValue<double> zmmeta (myReader,"zeta");

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

  // loop on events
  while(myReader.Next()){

    if(!isMC and *run > 274240) continue;

    // check trigger depending on the sample
    Double_t hlt   = 0.0;
    Double_t hltw  = 1.0;
    if (sample == Sample::sig || sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu)// single and double muon
      hlt = *hltm90+*hltm100+*hltm110+*hltm120+*hltmwm90+*hltmwm120+*hltmwm170+*hltmwm300;
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel) // single and double electron
      hlt = *hlte+*hltp165+*hltp175;      
    else if (sample == Sample::qcd || sample == Sample::gam) // single photon
      hlt = *hltp165+*hltp175;	        

    // Trigger Selection
    if (hlt  == 0) continue; // trigger

    // MET Filters --> apply on both data and monte-carlo
    if(*fhbhe == 0 || *fhbiso == 0 || *feeb == 0 || *fetp == 0 || *fvtx == 0 || *fcsc == 0) continue;

    // check dphi jet-met
    Double_t jmdphi = 0.0;    
    if (sample == Sample::sig || sample == Sample::wmn || sample == Sample::zmm || sample == Sample::topmu) jmdphi = fabs(*jmmdphi);
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel) jmdphi = fabs(*jemdphi);
    else if (sample == Sample::qcd || sample == Sample::gam) jmdphi = fabs(*jpmdphi);

    //set met
    Double_t pfmet = 0.0;
    Double_t pfmetphi = 0.0;
    if (sample == Sample::sig) {pfmet = *mmet; pfmetphi = *mmetphi;}
    else if (sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu){ pfmet = *mmet; pfmetphi = *mmetphi;}
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel){ pfmet = *emet; pfmetphi = *emetphi;}
    else if (sample == Sample::qcd || sample == Sample::gam)  { pfmet = *pmet; pfmetphi = *pmetphi;}

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
    else if (sample == Sample::qcd || sample == Sample::gam) {
      id1  = *phidm;
      id2  = 1.0;
      pt1  = *phpt;
      eta1 = *pheta;
    }
    

    // set zpt in case of Zsamples
    Double_t bosonPt = 0.0;
    if (sample == Sample::zmm)      bosonPt = *zmmpt; // di-muon CR
    else if (sample == Sample::zee) bosonPt = *zeept; // di-electron CR
    else if (sample == Sample::qcd or sample == Sample::gam) bosonPt = *phpt; // gamma+jets
    else if (sample == Sample::sig) bosonPt = pfmet; // missing energy in case of signal region for Z->nunu
    else if (sample == Sample::wen or sample == Sample::zee){ // single muon or single ele
      TLorentzVector lep4V, met4V;
      lep4V.SetPtEtaPhiM(pt1,eta1,phi1,0.);
      met4V.SetPtEtaPhiM(*met,0.,*metphi,0.);
      bosonPt = (lep4V+met4V).Pt();
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
    if (category != Category::VBF and *njets  < 1) continue; 
    else if(category == Category::VBF and *nincjets < 2) continue;

    // control regions wit one lepton --> tight requirement 
    if ((sample == Sample::wen || sample == Sample::wmn) && id1 != 1) continue;
    
    // photon control sample
    if ((sample == Sample::qcd || sample == Sample::gam) && *phpt < 175.) continue;
    if ((sample == Sample::qcd || sample == Sample::gam) && fabs(*pheta) > 1.4442) continue;

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
    if(category == Category::monojet or category == Category::inclusive){
      if (pfmet < 200.) continue;
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
   
    if(category != Category::VBF and leadingCentralJetPos < 0)  continue;
    if(category != Category::VBF and leadingCentralJetPos != 0) continue;

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
    if (isMC && sflhist && sfthist) {
      if (pt1 > 0.) {
	if (id1 == 1) sfwgt *= sfthist->GetBinContent(sfthist->FindBin(min(pt1,sfthist->GetXaxis()->GetBinLowEdge(sfthist->GetNbinsX()+1)-1),eta1)); 
	else          sfwgt *= sflhist->GetBinContent(sflhist->FindBin(min(pt1,sflhist->GetXaxis()->GetBinLowEdge(sflhist->GetNbinsX()+1)-1),eta1));
      }
      if (pt2 > 0.) {
	if (id2 == 1) sfwgt *= sfthist->GetBinContent(sfthist->FindBin(min(pt2,sfthist->GetXaxis()->GetBinLowEdge(sfthist->GetNbinsX()+1)-1),eta2));
	else          sfwgt *= sflhist->GetBinContent(sflhist->FindBin(min(pt2,sflhist->GetXaxis()->GetBinLowEdge(sflhist->GetNbinsX()+1)-1),eta2));
      }
    }
    
    // trigger scale factor for electrons
    if (isMC && triggerelhist && ( sample == Sample::zee || sample == Sample::topel || sample == Sample::wen)) {
      if (pt1 > 40. && id1 == 1 and id2 == 1)
	  sfwgt *= 1;
      else if(id1 == 1 and id2 != 1)
	sfwgt *= triggerelhist->GetBinContent(triggerelhist->FindBin(min(pt1,triggerelhist->GetXaxis()->GetBinLowEdge(triggerelhist->GetNbinsX()+1)-1),eta1));
      else if(id2 == 1 and id1 != 1)
	sfwgt *= triggerelhist->GetBinContent(triggerelhist->FindBin(min(pt2,triggerelhist->GetXaxis()->GetBinLowEdge(triggerelhist->GetNbinsX()+1)-1),eta2));
    }

    // photon id scale factor
    if (isMC && psfmedium && sample == Sample::gam) {
      if (pt1 > 0. && id1 == 1) {
	//	sfwgt *= psfmedium->GetBinContent(psfmedium->FindBin(min(pt1,psfmedium->GetXaxis()->GetBinLowEdge(psfmedium->GetNbinsX()+1)-1),eta1));
	sfwgt *= 1;
      }
    }

    // photon purity
    if (!isMC && purhist && sample == Sample::qcd) {
      sfwgt *= (1.0 - purhist->GetBinContent(purhist->FindBin(min(pt1,purhist->GetXaxis()->GetBinLowEdge(purhist->GetNbinsX()+1)-1), fabs(eta1))));
    }
    
    // met trigger scale factor
    if (isMC && triggermet && (sample == Sample::sig || sample == Sample::wmn || sample == Sample::zmm || sample == Sample::topmu)) {
      sfwgt *= triggermet->Eval(min(pfmet,triggermet->GetXaxis()->GetXmax()));
    }
        
    // photon trigger scale factor
    if(isMC && triggerphoton_graph && (sample == Sample::qcd || sample == Sample::gam)) // linear interpolation between graph points
      sfwgt *= triggerphoton_graph->Eval(min(pt1,triggerphoton_graph->GetXaxis()->GetXmax()));
    
    
    // b-tag weight
    double btagw = *wgtbtag;
    if( btagw > 2 || btagw <= 0)
      btagw = 1;
    
    //V-tagging scale factor --> only for mono-V
    if(isMC && category == Category::monoV && isWJet)
      sfwgt *= getVtaggingScaleFactor(tau2tau1,sysName);
    
    //Gen level info --> NLO re-weight    
    Double_t kwgt = 1.0;    
    double genpt = *wzpt;
    if (*wzpt < 150. ) genpt = 150.;
    if (*wzpt > 1000.) genpt = 999.;
    for (size_t i = 0; i < khists.size(); i++) {
      if (khists[i]) {
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

    // Top quark pt re-weight
    Double_t topptwgt = 1.0;
    //    if(reWeightTopPt)
      //      topptwgt = reweightTopQuarkPt(*toppt,*atoppt);

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
      if (jetpt->at(leadingCentralJetPos)  < 100.) continue;  // jet1 > 100 GeV
      if (jmdphi < 0.5) continue; // deltaPhi cut
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
	if (jetpt->at(leadingCentralJetPos)  < 100.) continue;  // jet1 > 100 GeV                                                                                          
	if (jmdphi < 0.5) continue; 

	if(boostedJetpt->size()  == 0) goodMonoJet = true;

	if(boostedJetpt->size() > 0){ // in case one boosted jet

	  if(fabs(boostedJeteta->at(0)) > jetEtaAK8) 
	    goodMonoJet = true;
	  
	  if(boostedJetpt->at(0) < ptJetMinAK8) // check pT
	    goodMonoJet = true;
	 	  
	  else{ // if high pT check pruned mass

	    TLorentzVector jetak4, jetak8;
	    jetak4.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	    jetak8.SetPtEtaPhiM(boostedJetpt->at(0),boostedJeteta->at(0),boostedJetphi->at(0),boostedJetm->at(0));

	    if(jetak4.DeltaR(jetak8) > 0.8) continue;	    	    
	    // pruned mass selection
	    if(prunedJetm->at(0) < prunedMassMin  or prunedJetm->at(0) > prunedMassMax)
	      goodMonoJet= true;

	    // tau2tau1 selection
	    if(boostedJettau2->at(0)/boostedJettau1->at(0) > tau2tau1)
	      goodMonoJet= true;
	  }
	}
      	
	if(not goodMonoJet) continue;
      }

      else if(category == Category::monoV or category == Category::boosted or category == Category::prunedMass or category == Category::tau2tau1){
	
	if (centralJets.size() == 0) continue;
	if(boostedJetpt->size() == 0) continue;
	if(boostedJetpt->at(0) < ptJetMinAK8) continue;
	if(fabs(boostedJeteta->at(0)) > jetEtaAK8) continue;

	TLorentzVector jetak4, jetak8;
	jetak4.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jetak8.SetPtEtaPhiM(boostedJetpt->at(0),boostedJeteta->at(0),boostedJetphi->at(0),boostedJetm->at(0));
	
	// match leading ak4 and leading ak8 within 0.8 cone
	if(jetak4.DeltaR(jetak8) > 0.8) continue;
	
	//after match apply jetid on leading ak4
	if (chfrac->at(leadingCentralJetPos) < 0.1) continue;   // jet id                                                                                                     
	if (nhfrac->at(leadingCentralJetPos) > 0.8) continue;   // jet id                                                                                                  
	if (jetpt->at(leadingCentralJetPos)  < 100.) continue;  // jet1 > 100 GeV                                                                                           
	if (jmdphi < 0.5) continue; // deltaPhi cut                                                                                                       

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
	  cout<<"hadBosonpt "<<*hadBosonpt<<" hadBosoneta "<<*hadBosoneta<<" mass "<<*hadBosonm<<endl;
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
	
	if(not goodMonoV) continue;	

      }
      else if(category == Category::VBF){
	if(centralJets.size()+forwardJets.size() < 2) continue;
	if(fabs(jeteta->at(0)) > 4.7 || fabs(jeteta->at(1)) > 4.7) continue;
	if(jetpt->at(0) < 100) continue;
	if(jmdphi < 0.5) continue; // deltaPhi cut                                                                                                                             
	if(fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) < 0.1) continue;
	if(fabs(jeteta->at(0)) < 2.5 and nhfrac->at(0) > 0.8) continue;
	if(fabs(jeteta->at(0)) > 3.0) continue;
	if(jetpt->at(1) < 40) continue;
	if(jeteta->at(0)*jeteta->at(1) > 0 ) continue;
	if(fabs(jeteta->at(0)-jeteta->at(1)) < 3) continue;
	TLorentzVector jet1 ;
	TLorentzVector jet2 ;
	jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
	if((jet1+jet2).M() < 500) continue;
      }
    }

    // fill 1D histogram
    double fillvar = 0;
    // fill the histograms --> with the right observable
    for(auto hist : hist1D){
      TString name(hist->GetName());      
      if(name.Contains("nvtx"))
	fillvar = *nvtx;
      else if(name.Contains("chfrac")){
	if(category == Category::VBF)
	  fillvar = chfrac->at(0);
	else
	  fillvar = chfrac->at(leadingCentralJetPos);
      }
      else if(name.Contains("nhfrac")){
	if(category == Category::VBF)
	  fillvar = nhfrac->at(0);
	else
	  fillvar = nhfrac->at(leadingCentralJetPos);
      }
      else if(name.Contains("emfrac")){
	if(category == Category::VBF)
	  fillvar = emfrac->at(0);
	else
	  fillvar = emfrac->at(leadingCentralJetPos);
      } 
      else if(name.Contains("boostedjetpt")){
	if(boostedJetpt->size() > 0)
	  fillvar = boostedJetpt->at(0);
	else
	  fillvar = 0.;
      }
      else if(name.Contains("boostedjeteta")){
	if(boostedJeteta->size() > 0)
	  fillvar = boostedJeteta->at(0);
	else
	  fillvar = 0.;
      }
      else if(name.Contains("jetmetdphi"))
	fillvar = jmdphi;
      else if(name.Contains("met"))
	fillvar = pfmet;            
      else if(name.Contains("jetpt2") and jetpt->size() >= 2)
	fillvar = jetpt->at(1);
      else if(name.Contains("jeteta2") and jetpt->size() >= 2){
	fillvar = jeteta->at(1);
      }
      else if(name.Contains("jetpt")){
	if(category == Category::VBF)
	  fillvar = jetpt->at(0);
	else
	  fillvar = jetpt->at(leadingCentralJetPos);
      }
      else if(name.Contains("jeteta")){
	if(category == Category::VBF)
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
	float deltaPhi = fabs(jetphi->at(0)-pfmetphi);
	if(deltaPhi > TMath::Pi())
	  deltaPhi = fabs(2*TMath::Pi() - deltaPhi);
	fillvar = sqrt(2*jetpt->at(0)*pfmet*(1-cos(deltaPhi)));
      }   
      else if(name.Contains("ncjet"))
	fillvar = *njets;      
      else if(name.Contains("njet"))
	fillvar = *nincjets;      
      else if(name.Contains("nbjet_hpt_loose")){
	int nbjet = 0;
	for(size_t iJet = 0; iJet < jetbtag->size(); iJet++){
	  if(jetbtag->at(iJet) > 0.460 and jetpt->at(iJet) > 30)
	    nbjet++;
	}
	fillvar = nbjet;
      }      
      else if(name.Contains("nbjet_hpt")){
	int nbjet = 0;
	for(size_t iJet = 0; iJet < jetbtag->size(); iJet++){
	  if(jetbtag->at(iJet) > 0.80 and jetpt->at(iJet) > 30)
	    nbjet++;
	}
	fillvar = nbjet;
      }      
      else if(name.Contains("nbjet"))
	fillvar = *nbjets;
      else if(name.Contains("bosonpt"))
	fillvar = bosonPt;    	
      else if(name.Contains("QGL")){
	fillvar = jetQGL->at(0);
      }
      else if(name.Contains("QGL_AK8")){
	if(boostedJetpt->size() > 0)
	  fillvar = boostedJetQGL->at(0);
	else
	  fillvar = -1.;
      }
      
      // substructure
      else if(name.Contains("mpruned")){
	if( prunedJetm->size() > 0 and boostedJetpt->at(0) > ptJetMinAK8 )
	  fillvar = prunedJetm->at(0);	
	else fillvar = 0.;
      }
      else if(name.Contains("ht"))
	fillvar = *ht;      
      else if(name.Contains("tau2tau1")){
	if( boostedJettau1->size() > 0 and boostedJettau2->size() > 0 and boostedJetpt->at(0) > ptJetMinAK8 )
	  fillvar = boostedJettau2->at(0)/boostedJettau1->at(0);	
	else fillvar = 0.;
      }

      // b-tagging
      else if(name.Contains("btagCSV_max")){
	float btagMax = -10.;
	for(size_t iBjet = 0; iBjet < jetbtag->size(); iBjet++){
	  if(jetbtag->at(iBjet) > btagMax)
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
	  if(jetbtag->at(iBjet) < btagMin)
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
	    float deltaPhi = fabs(jetphi->at(0)-jetphi->at(jjet));
	    if(jetpt->at(0) < 30 or jetpt->at(jjet) < 30) continue;
	    if(deltaPhi > TMath::Pi())
	      deltaPhi = 2*TMath::Pi() - deltaPhi;
	    if(deltaPhi > 0 and deltaPhi < minDphi){
	      minDphi = deltaPhi;
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
	if(jetphi->size() < 2)
	  fillvar = hist->GetXaxis()->GetBinCenter(1);
	else {
	  fillvar = fabs(jetphi->at(0)-jetphi->at(1));
	  if(fillvar > TMath::Pi())
	    fillvar = 2*TMath::Pi()-fillvar;
	}       
      }
      
      // overflow bin
      if (fillvar >= hist->GetBinLowEdge(hist->GetNbinsX())+hist->GetBinWidth(hist->GetNbinsX())) 
	fillvar = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());

      // total event weight
      double evtwgt  = 1.0;
      Double_t puwgt = 0.;
      if (isMC and not reweightNVTX){
	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(*wgtpileup)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*hwgt*pfwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)
	else{
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(*wgtpileup)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*hwgt*pfwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)
	}
      }
      else if (isMC and reweightNVTX){
	if (*nvtx <= 40) 
	  puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*topptwgt*sfwgt*kwgt*hwgt*ggZHwgt*pfwgt/(*wgtsum);
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*topptwgt*sfwgt*kwgt*hwgt*ggZHwgt*pfwgt/(*wgtsum);
      }
      
      if (!isMC && sample == Sample::qcd) 
	evtwgt = sfwgt*pfwgt;
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
	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(*wgtpileup)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*hwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)      
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(*wgtpileup)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*hwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)      
      }
      else if (isMC and reweightNVTX){
	if (*nvtx <= 35) 
	  puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*hwgt/(*wgtsum);
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*sfwgt*topptwgt*ggZHwgt*kwgt*hwgt/(*wgtsum);
      }
      if (!isMC && sample == Sample::qcd) 
	evtwgt = sfwgt;
      else if (!isMC)
	evtwgt = hltw;
      

	hist->Fill(fillvarX,fillvarY,evtwgt);            
     }
  }

  sffile_eleTight.Close();
  sffile_eleVeto.Close();
  sffile_muTight.Close();
  sffile_muLoose.Close();
  sffile_phoLoose.Close();
  sffile_phoMedium.Close();
  purityfile_photon.Close();
  triggerfile_SinglEle.Close();
  triggerfile_SingleMu.Close();
  triggerfile_MET.Close();
  triggerfile_SinglePhoton.Close();
  
}

#endif
