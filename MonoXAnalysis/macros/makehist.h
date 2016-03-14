#ifndef MAKEHIST_H
#define MAKEHIST_H

#include <vector>
#include <fstream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TString.h"

#include "histoUtils.h"

using namespace std;

// some basic cut values
const float tau2tau1      = 0.6;
const float tau2tau1LP    = 0.75;
const float prunedMassMin = 65.;
const float prunedMassMax = 105.;
const float ptJetMinAK8   = 250.;
const float jetEtaAK8     = 2.4;
const float pfMetMonoVLower = 250.;
const float pfMetMonoVUpper = 8000.;
const int   vBosonCharge   = 0;
const int   nBjets         = 1;

string kfactorFile       = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_Wseparated_24bins.root";
//string kfactorFile     = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root";
string kfactorFileUnc    = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/scalefactors_v4.root";
string baseInputTreePath = "/home/rgerosa/MONOJET_ANALYSIS/Production-19-2-2016/";

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


void makehist4(TTree* tree, /*input tree*/ 
	       vector<TH1*> hist1D, /* set of 1D histogram */ 
	       vector<TH2*> hist2D, /* set of 2D histogram */ 
	       bool   isMC, 
	       int    sample, 
	       int    category,
	       bool   isWJet,
	       double scale,
	       double lumi,	       
	       int    QGLweight,
	       vector<TH1*> khists, 
	       string sysName,	
	       bool   reWeightTopPt = false,
	       bool   reweightNVTX  = true,
	       int    resonantSelection = 0,
	       bool   isHiggsInvisible  = false, // reject VBF events
	       float  XSEC = -1.,// fix the cross section from extern
	       TH1*  hhist = NULL
	       ) {

  if(not tree){
    cout<<" empty tree --> skip process "<<endl;
    return;
  }


 
  //  ofstream dump("dump_sample_"+to_string(sample)+".txt");
  
  // in case you want to weight the NVTX distribution
  TFile* pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt.root");
  TH1*   puhist = (TH1*) pufile->Get("puhist");
    
  // Lepton ID scale factor from tag and probe: muons, electrons 
  TFile* sffile  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF/leptonIDsfs.root");
  
  TH2*  msflhist = (TH2*)sffile->Get("muon_loose_SF");
  TH2*  msfthist = (TH2*)sffile->Get("muon_tight_SF");
  
  TH2*  esflhist = (TH2*)sffile->Get("electron_veto_SF");
  TH2*  esfthist = (TH2*)sffile->Get("electron_tight_SF");
  
  // Photon ID scale factor from tag and probe
  TFile* psffile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/PhotonSFandEffandPurity_Lumi2p1fb_0202.root");
  TH2*  psfhist  = (TH2*)psffile->Get("PhotonSF");  
  TH2*  purhist  = (TH2*)psffile->Get("PhotonPurity");
  
  // trigger efficiency correction for single electron trigger
  TFile* trefile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF/leptonTrigsfs.root");
  TH2*   trehist = (TH2*)trefile->Get("hltel27_SF");
  
  // trigger efficiency for met trigger
  TFile* trmfile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF/mettrigSF.root");
  TH1*   trmhist = (TH1*) trmfile->Get("mettrigSF");

  // QGL rewight
  TFile* QGLReweight = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/QGLWeight/QGLWeight.root ");
  TH2* QGLWeightHist = NULL;

  if(QGLweight == 1 and sysName != "QGLup" and sysName != "QGLdw") // zmumu
    QGLWeightHist = (TH2*) QGLReweight->Get("QGL_weight_Z");
  else if(QGLweight == 1 and sysName == "QGLup")
    QGLWeightHist = (TH2*) QGLReweight->Get("QGL_weight_Z_up");
  else if(QGLweight == 1 and sysName == "QGLdw")
    QGLWeightHist = (TH2*) QGLReweight->Get("QGL_weight_Z_dw");

  if(QGLweight == 2 and sysName != "QGLup" and sysName != "QGLdw") // zmumu
    QGLWeightHist = (TH2*) QGLReweight->Get("QGL_weight_W");
  else if(QGLweight == 3 and sysName == "QGLup")
    QGLWeightHist = (TH2*) QGLReweight->Get("QGL_weight_W_up");
  else if(QGLweight == 3 and sysName == "QGLdw")
    QGLWeightHist = (TH2*) QGLReweight->Get("QGL_weight_W_dw");
  
  if(QGLweight == 3 and sysName != "QGLup" and sysName != "QGLdw") // zmumu
    QGLWeightHist = (TH2*) QGLReweight->Get("QGL_weight_G");
  else if(QGLweight == 3 and sysName == "QGLup")
    QGLWeightHist = (TH2*) QGLReweight->Get("QGL_weight_G_up");
  else if(QGLweight == 3 and sysName == "QGLdw")
    QGLWeightHist = (TH2*) QGLReweight->Get("QGL_weight_G_dw");

  if(QGLweight == 4 and sysName != "QGLup" and sysName != "QGLdw") // zmumu
    QGLWeightHist = (TH2*) QGLReweight->Get("QGL_weight_T");
  else if(QGLweight == 4 and sysName == "QGLup")
    QGLWeightHist = (TH2*) QGLReweight->Get("QGL_weight_T_up");
  else if(QGLweight == 4 and sysName == "QGLdw")
    QGLWeightHist = (TH2*) QGLReweight->Get("QGL_weight_T_dw");

  // histogram to be filled
  for(size_t ihist  = 0 ; ihist < hist1D.size(); ihist++)
    hist1D.at(ihist)->Sumw2();
  for(size_t ihist  = 0 ; ihist < hist2D.size(); ihist++)
    hist2D.at(ihist)->Sumw2();
  
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
    if(sysName == "btagUp" or sysName == "bUp")
      btagname = "wgtbtagUp";
    else if(sysName == "btagDown" or sysName == "btagDw" or sysName == "bDown" or sysName == "bDw")
      btagname = "wgtbtagDown";
    else
      btagname = "wgtbtag";

    // tem fix since not all the samples have been re-produced with hltphoton120
    if(sample != 5 and sample != 6){
      prescalename  = "wgt";
      hltphotonname = "hltphoton165";
    }
    else{
      prescalename  = "pswgt";
      hltphotonname = "hltphoton120";
    }
  }
  else if(!isMC){ // in case of data
    wgtname       = "wgt";
    btagname      = "wgt";
    wgtpileupname = "wgt";
    prescalename  = "pswgt";
    hltphotonname = "hltphoton120";    
  }

  TTreeReaderValue<double> wgtsum       (myReader,wgtname.c_str());
  TTreeReaderValue<double> wgtpileup    (myReader,wgtpileupname.c_str());
  TTreeReaderValue<double> wgtbtag      (myReader,btagname.c_str());
  
  // trigger
  TTreeReaderValue<UChar_t> hltm90     (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm120    (myReader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120  (myReader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170  (myReader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300  (myReader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90   (myReader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hlte       (myReader,"hltsingleel");
  TTreeReaderValue<UChar_t> hltp120    (myReader, hltphotonname.c_str());
  TTreeReaderValue<UChar_t> hltp165    (myReader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175    (myReader,"hltphoton175");
  // prescale
  TTreeReaderValue<double> pswgt(myReader,prescalename.c_str());

  TTreeReaderValue<UChar_t> fhbhe  (myReader,"flaghbheloose");
  TTreeReaderValue<UChar_t> fhbiso (myReader,"flaghbheiso");
  
  // MET filters
  string cscname;
  string feebname;
  string fmutrkname;
  string fbtrkname;

  if(isMC){
    cscname    = "flagcsctight";
    feebname   = "flageebadsc";
    fmutrkname = "flagcsctight";
    fbtrkname  = "flagcsctight";
  }
  else{
    cscname    = "flagcscnew";
    feebname   = "flageescnew";
    fmutrkname = "flagmuontrack";
    fbtrkname  = "flagbadtrack";
  }
  
  TTreeReaderValue<UChar_t> fcsc (myReader,cscname.c_str());
  TTreeReaderValue<UChar_t> feeb (myReader,feebname.c_str());

  TTreeReaderValue<UChar_t> fmutrk (myReader,fmutrkname.c_str());
  TTreeReaderValue<UChar_t> fbtrk  (myReader,fbtrkname.c_str());

  TTreeReaderValue<unsigned int> njets  (myReader,"njets");
  TTreeReaderValue<unsigned int> nbjets (myReader,"nbjetslowpt");

  TTreeReaderValue<double> j1pt (myReader,"leadingjetpt");
  TTreeReaderValue<double> ht   (myReader,"ht");
  
  TTreeReaderValue<vector<double> > jetpt   (myReader,"centraljetpt");
  TTreeReaderValue<vector<double> > jetQGL  (myReader,"centraljetQGL");
  TTreeReaderValue<vector<double> > jeteta  (myReader,"centraljeteta");
  TTreeReaderValue<vector<double> > jetphi  (myReader,"centraljetphi");
  TTreeReaderValue<vector<double> > jetbtag (myReader,"centraljetbtag");
  TTreeReaderValue<vector<double> > jetm    (myReader,"centraljetm");
  TTreeReaderValue<vector<double> > chfrac  (myReader,"centraljetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac  (myReader,"centraljetNHfrac");
  TTreeReaderValue<vector<double> > emfrac  (myReader,"centraljetEMfrac");

  TTreeReaderValue<vector<double> > jetpt_fwd   (myReader,"forwardjetpt");
  TTreeReaderValue<vector<double> > jeteta_fwd  (myReader,"forwardjeteta");
  TTreeReaderValue<vector<double> > jetphi_fwd  (myReader,"forwardjetphi");
  TTreeReaderValue<vector<double> > jetm_fwd    (myReader,"forwardjetm");

  // AK8 jet
  TTreeReaderValue<vector<double> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<double> > boostedJetQGL   (myReader,"boostedJetQGL");
  TTreeReaderValue<vector<double> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<double> > boostedJetphi   (myReader,"boostedJetphi");
  TTreeReaderValue<vector<double> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<double> > prunedJetm      (myReader,"prunedJetm");
  TTreeReaderValue<vector<double> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<double> > boostedJettau1  (myReader,"boostedJettau1");
  TTreeReaderValue<vector<double> > boostedJetBosoneta  (myReader,"boostedJetBosoneta");
  TTreeReaderValue<vector<double> > boostedJetBosonphi  (myReader,"boostedJetBosonphi");
  TTreeReaderValue<vector<double> > boostedJetBosonpt   (myReader,"boostedJetBosonpt");
  TTreeReaderValue<vector<double> > boostedJetBosonm    (myReader,"boostedJetBosonm");

  // met
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

  TTreeReaderValue<int> mu1pid (myReader,"mu1pid");
  TTreeReaderValue<int> mu2pid (myReader,"mu2pid");
  TTreeReaderValue<int> mu1id (myReader,"mu1id");
  TTreeReaderValue<int> mu2id (myReader,"mu2id");
  TTreeReaderValue<double> mu1pt (myReader,"mu1pt");
  TTreeReaderValue<double> mu2pt (myReader,"mu2pt");
  TTreeReaderValue<double> mu1eta (myReader,"mu1eta");
  TTreeReaderValue<double> mu2eta (myReader,"mu2eta");
  TTreeReaderValue<double> mu1phi (myReader,"mu1phi");
  TTreeReaderValue<double> mu2phi (myReader,"mu2phi");

  TTreeReaderValue<int> el1pid (myReader,"el1pid");
  TTreeReaderValue<int> el2pid (myReader,"el2pid");
  TTreeReaderValue<int> el1id (myReader,"el1id");
  TTreeReaderValue<int> el2id (myReader,"el2id");
  TTreeReaderValue<double> el1pt (myReader,"el1pt");
  TTreeReaderValue<double> el2pt (myReader,"el2pt");
  TTreeReaderValue<double> el1eta (myReader,"el1eta");
  TTreeReaderValue<double> el2eta (myReader,"el2eta");
  TTreeReaderValue<double> el1phi (myReader,"el1phi");
  TTreeReaderValue<double> el2phi (myReader,"el2phi");
  
  TTreeReaderValue<int> phidm (myReader,"phidm");
  TTreeReaderValue<double> phpt (myReader,"phpt");
  TTreeReaderValue<double> pheta (myReader,"pheta");
  TTreeReaderValue<double> phphi (myReader,"phphi");
  
  TTreeReaderValue<double> wzpt (myReader,"wzpt");
  TTreeReaderValue<double> wzeta (myReader,"wzeta");
  TTreeReaderValue<double> zmass (myReader,"zmass");
  TTreeReaderValue<double> zmmpt (myReader,"zpt");
  TTreeReaderArray<double> zeept (myReader,"zeept.zeeept");
  TTreeReaderValue<double> zeeeta (myReader,"zeeeta");
  TTreeReaderValue<double> zmmeta (myReader,"zeta");

  TTreeReaderValue<double> dmpt (myReader,"dmpt");

  // other trick to handle the fact that this info is actually only stored for top/s-top samples
  string topptname;
  string atopptname;
  if(reWeightTopPt and isMC){
    topptname = "toppt";
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
    // check trigger depending on the sample
    Double_t hlt   = 0.0;
    Double_t hltw  = 1.0;
    if (sample == 0 || sample == 1 || sample == 2 || sample == 7){
      hlt = *hltm90+*hltm120+*hltmwm120+*hltmwm170+*hltmwm300+*hltmwm90;
    }
    else if (sample == 3 || sample == 4 || sample == 8) {
      hlt = *hlte+*hltp165+*hltp175;
    }
    else if (sample == 5 || sample == 6){
      if(isMC){
	  hlt = *hltp165+*hltp175+*hltp120;	  
      }
      else{
	if(*hltp175 || *hltp165)
	  hlt = *hltp165+*hltp175;
	else if(*hltp120 and not *hltp175 and not *hltp165){
	  hlt    = *hltp120;
	  hltw  *= *pswgt;	 
	}
      }
    }
        
    // check dphi jet-met
    Double_t jmdphi = 0.0;
    
    if (sample == 0 || sample == 1 || sample == 2 || sample == 7) jmdphi = fabs(*jmmdphi);
    else if (sample == 3 || sample == 4 || sample == 8)           jmdphi = fabs(*jemdphi);
    else if (sample == 5 || sample == 6)           jmdphi = fabs(*jpmdphi);

    //set met
    Double_t pfmet = 0.0;
    Double_t pfmetphi = 0.0;
    if (sample == 0) {pfmet = *mmet; pfmetphi = *mmetphi;}
    else if (sample == 1 || sample == 2 || sample == 7){ pfmet = *mmet; pfmetphi = *mmetphi;}
    else if (sample == 3 || sample == 4 || sample == 8)          { pfmet = *emet; pfmetphi = *emetphi;}
    else if (sample == 5 || sample == 6)          { pfmet = *pmet; pfmetphi = *pmetphi;}
    else if (sample == 7 and (*hlte or *hltp165 or *hltp175))    { pfmet = *emet; pfmetphi = *emetphi;}
    else if (sample == 7 and not *hlte)           { pfmet = *mmet; pfmetphi = *mmetphi;}

    // propagate met systeamtics on the recoil
    if(metSuffix != ""){
      pfmet += (*met-*metOriginal);
    }
    

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

    if (sample == 1 || sample == 2 || sample == 7) {
      id1  = *mu1id;
      id2  = *mu2id;
      pt1  = *mu1pt;
      pt2  = *mu2pt;
      pid1 = *mu1pid;
      pid2 = *mu2pid;
      eta1 = fabs(*mu1eta);
      eta2 = fabs(*mu2eta);
      phi1 = fabs(*el1phi);
      phi2 = fabs(*el2phi);
    }
    else if (sample == 3 || sample == 4 || sample == 8) {
      id1  = *el1id;
      id2  = *el2id;
      pt1  = *el1pt;
      pt2  = *el2pt;
      eta1 = fabs(*el1eta);
      eta2 = fabs(*el2eta);
      phi1 = fabs(*el1phi);
      phi2 = fabs(*el2phi);
      pid1 = *el1pid;
      pid2 = *el2pid;
    }
    else if (sample == 5 || sample == 6) {
      id1  = 1.0;
      id2  = 1.0;
      pt1  = *phpt;
      eta1 = fabs(*pheta);
    }
    
    if (pt1 >= 1000.) 
      pt1 = 999.0;
    if (pt2 >= 1000.) 
      pt2 = 999.0;


    // set zpt in case of Zsamples
    Double_t bosonPt = 0.0;
    if (sample == 1)      bosonPt = *zmmpt; // di-muon CR
    else if (sample == 3) bosonPt = zeept[0]; // di-electron CR
    else if (sample == 5 or sample == 6) bosonPt = *phpt; // gamma+jets
    else if (sample == 0) bosonPt = pfmet; // missing energy in case of signal region for Z->nunu
    else if (sample == 2 or sample == 4){ // single muon or single ele
      TLorentzVector lep4V, met4V;
      lep4V.SetPtEtaPhiM(pt1,eta1,phi1,0.);
      met4V.SetPtEtaPhiM(*met,0.,*metphi,0.);
      bosonPt = (lep4V+met4V).Pt();
    }

    // scale factor for leptons
    TH2* sflhist = NULL;
    TH2* sfthist = NULL;
    
    if (sample == 1 || sample == 2 || sample == 7) {
	sflhist = msflhist;
	sfthist = msfthist;
    }
    if (sample == 3 || sample == 4 || sample == 8) {
	sflhist = esflhist;
	sfthist = esfthist;
    }
    
    Double_t sfwgt = 1.0;
    if (isMC && sflhist && sfthist) {
      if (pt1 > 0.) {
	if (id1 == 1) sfwgt *= sfthist->GetBinContent(sfthist->FindBin(pt1, eta1)); 
	else          sfwgt *= sflhist->GetBinContent(sflhist->FindBin(pt1, eta1)); 
      }
      if (pt2 > 0.) {
	if (id2 == 1) sfwgt *= sfthist->GetBinContent(sfthist->FindBin(pt2, eta2)); 
	else          sfwgt *= sflhist->GetBinContent(sflhist->FindBin(pt2, eta2)); 
      }
    }
    
    // trigger scale factor
    if (isMC && trehist && ( sample == 4 || sample == 8)) {
      if (pt1 > 0. && id1 == 1) {
	sfwgt *= trehist->GetBinContent(trehist->FindBin(pt1, eta1));
      }
    }

    // photon id scale factor
    if (isMC && psfhist && (sample == 5 || sample == 6)) {
      if (pt1 > 0. && id1 == 1) {
	sfwgt *= psfhist->GetBinContent(psfhist->FindBin(pt1, eta1));
      }
    }

    // photon purity
    if (!isMC && purhist && sample == 6) {
      if(pt1 < purhist->GetXaxis()->GetBinLowEdge(1) and id1 == 1){
	sfwgt *= (1.0 - purhist->GetBinContent(purhist->FindBin(purhist->GetXaxis()->GetBinLowEdge(1)+purhist->GetXaxis()->GetBinWidth(1)/2, eta1)));
      }
      else if (pt1 > purhist->GetXaxis()->GetBinLowEdge(1) && id1 == 1) {
	sfwgt *= (1.0 - purhist->GetBinContent(purhist->FindBin(pt1, eta1)));
      }     
    }

    // met trigger scale factor
    if (isMC && trmhist && (sample == 0 || sample == 1 || sample == 2 || sample == 7)) {
      sfwgt *= trmhist->GetBinContent(trmhist->FindBin(pfmet));
    }

    // QGL weight    
    if(isMC and QGLWeightHist and QGLweight != 0){
      if(jetQGL->size() > 0)
	sfwgt *= QGLWeightHist->GetBinContent(QGLWeightHist->FindBin(jetpt->at(0),jetQGL->at(0)));
      else
	sfwgt *= 1.;
    }

    // b-tag weight
    double btagw = *wgtbtag;
    if( btagw > 2 || btagw < 0)
      btagw = 1;
    
    //V-tagging scale factor --> only for mono-V
    if(isMC && category == 2 && isWJet){
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
    }
    else if(isMC && category == 3 && isWJet){
      if(tau2tau1 == 0.45){
	if(sysName == "VtagUp")
	  sfwgt *= (1.458+0.381);
	else if(sysName == "VtagDown")
	  sfwgt *= (1.458-0.381);
	else
	  sfwgt *= 1.458;
      }
      else if(tau2tau1 == 0.6){
	if(sysName == "VtagUp")
	  sfwgt *= (0.881+0.490);
	else if(sysName == "VtagDown")
	  sfwgt *= (0.881-0.490);
	else
	  sfwgt *= 0.881;
      }
    }

    //Gen level info --> NLO re-weight    
    Double_t kwgt = 1.0;    
    for (unsigned i = 0; i < khists.size(); i++) {      
      // in case of charge independent k-factors
      if(TString(khists[i]->GetName()).Contains("_Wp") and pid1 < 0)
	continue;
      else if(TString(khists[i]->GetName()).Contains("_Wm") and pid1 > 0)
	continue;

      if (isMC && khists[i]){
	if(*wzpt < khists[i]->GetBinLowEdge(1))
	  *wzpt  = khists[i]->GetBinLowEdge(1)+1.;	
	else if(*wzpt > khists[i]->GetBinLowEdge(khists[i]->GetNbinsX()+1))
	  *wzpt = khists[i]->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1.;
	
	kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(*wzpt));
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

    // Top quark pt re-weight
    Double_t topptwgt = 1.0;
    if(reWeightTopPt)
      topptwgt = reweightTopQuarkPt(*toppt,*atoppt);
        
    // Trigger Selection
    if (hlt  == 0) continue; // trigger
    // MET Filters
    if (not isHiggsInvisible and (*fhbhe == 0 || *fhbiso == 0 || *fcsc == 0 || *feeb == 0)) continue;

    // Additional met filters
    if (not isMC){
      if(*fmutrk == 0 || *fbtrk == 0) continue; // met filters
    }

    // N-jets
    if (*njets  < 1) continue; 

    // B-veto, not for top control sample
    if (*nbjets > 0 and sample != 7 and sample != 8) continue; 

    // control regions with two leptons --> opposite charge
    if (sample == 1 && *mu1pid == *mu2pid) continue;
    if (sample == 3 && *el1pid == *el2pid) continue;

    // control regions with two leptons --> one should be tight
    if ((sample == 1 || sample == 3)){
	if(sample == 1 and id1 == 1){
	  if(pt1 < 20) continue;
	}
	else if(sample == 1 and id2 == 1){
	  if(pt2 < 20) continue;
	}

	if(sample == 3 and id1 == 1){
	  if(pt1 < 40) continue;
	}
	else if(sample == 3 and id2 == 1){
	  if(pt2 < 40) continue;
	}	
    }

    // control regions wit one lepton --> tight requirement 
    if ((sample == 2 || sample == 4) && id1 != 1) continue;
    
    // photon control sample
    if ((sample == 5 || sample == 6) && *phpt < 120.) continue;
    if ((sample == 5 || sample == 6) && fabs(*pheta) > 1.4442) continue;

    // Wenu kill QCD
    if (sample == 4 && *met < 50.) continue;

    // n-bjets cut for unboosted categories
    if ((sample == 7 || sample == 8) && (category !=2 and category !=3)  && *nbjets < nBjets) continue;
    if ( sample == 7 || sample == 8){ // select only events with one lepton
      // at least one lepton in the plateau region
      if(pt1 <=0 or id1 != 1) continue;
      if(abs(pid1) == 13 && pt1 < 20. ) continue;
      if(abs(pid1) == 11 && pt1 < 40. ) continue;
      // met cut
      if(sample == 8 && *met < 50.) continue;
      // veto di-lepton events
      if(pt2 > 0) continue;
    }
    
    // met selection
    if(category <= 1){
      if (pfmet < 200.) continue;
    }
    else{
      if(pfmet < pfMetMonoVLower) continue;
      if(pfmet > pfMetMonoVUpper) continue;
    }

    //apply charge cut to separate W+ from W- in case
    if(sample == 2 or sample == 4){ // in case charge is required to be 1 skip events with negative leptons, viceversa
      if(vBosonCharge == 1 and pid1 > 0) continue;
      else if(vBosonCharge == -1 and pid1 < 0) continue;
    }


    if(isHiggsInvisible){ // veto Higgs VBF events

      vector<TLorentzVector> jets;
      for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
	TLorentzVector vec;
	vec.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
	if(fabs(vec.Eta()) > 4.7) continue;
	jets.push_back(vec);
      }
      for(size_t ijet = 0; ijet < jetpt_fwd->size(); ijet++){
	TLorentzVector vec;
	vec.SetPtEtaPhiM(jetpt_fwd->at(ijet),jeteta_fwd->at(ijet),jetphi_fwd->at(ijet),jetm_fwd->at(ijet));
	if(fabs(vec.Eta()) > 4.7) continue;
	jets.push_back(vec);
      }

      if(jets.size() >= 2){
	std::sort(jets.begin(),jets.end(),jetSorter);
	
	if(jets.at(0).Pt() > 80 and jets.at(1).Pt() > 70 and fabs(jets.at(0).Eta()-jets.at(1).Eta()) > 3.6 and (jets.at(0)+jets.at(1)).M() > 1100 and jmdphi > 2.3)
	  continue;
      }     
    }

    // inclusive mono-jet analysis 
    if(category == 0){ 
      if (chfrac->size() == 0 || nhfrac->size() == 0 or jetpt->size() == 0) continue; // at least one leading jet
      if (chfrac->at(0) < 0.1) continue;   // jet id
      if (nhfrac->at(0) > 0.8) continue;   // jet id
      if (jetpt->at(0)  < 100.) continue;  // jet1 > 100 GeV
      if (jetpt->at(0)  < *j1pt) continue; 
      if (sample != 4 and jmdphi < 0.5) continue; // deltaPhi cut
    }
    else{

      // boosted category and monojet
      bool goodMonoJet = false;
      bool goodMonoV   = false;
      
      if(category == 1){ // mono jet + V-jet veto

	if (chfrac->size() == 0 || nhfrac->size() == 0 or jetpt->size() == 0) continue; // at least one leading jet                                                         
	if (chfrac->at(0) < 0.1) continue;   // jet id                                                                                                                       
	if (nhfrac->at(0) > 0.8) continue;   // jet id                                                                                                                        
	if (jetpt->at(0)  < 100.) continue;  // jet1 > 100 GeV                                                                                                             
	if (jetpt->at(0)  < *j1pt) continue;
	if (sample !=4 and jmdphi < 0.5) continue; // deltaPhi cut                                                                                                          

	if(boostedJetpt->size()  == 0)  // no boosted jets (AK8 pT > 200 GeV)
	  goodMonoJet = true;

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
	    	    
	    // jet met dphi     
	    float deltaPhi = 0.;	    
	    if(fabs(pfmetphi-jetak8.Phi()) > TMath::Pi())
	      deltaPhi = fabs(2*TMath::Pi() - fabs(pfmetphi-jetak8.Phi()));
	    else
	      deltaPhi = fabs(pfmetphi-jetak8.Phi());
		
	    if (sample != 4 and deltaPhi < 0.5) continue; // deltaPhi cut                                                                                             

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

      else if(category >= 2){
	
	if(chfrac->size() == 0 || nhfrac->size() == 0 or jetpt->size() == 0) continue; 	
	if(boostedJetpt->size() == 0) continue;
	if(boostedJetpt->at(0) < ptJetMinAK8) continue;
	if(fabs(boostedJeteta->at(0)) > jetEtaAK8) continue;

	TLorentzVector jetak4, jetak8;
	jetak4.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jetak8.SetPtEtaPhiM(boostedJetpt->at(0),boostedJeteta->at(0),boostedJetphi->at(0),boostedJetm->at(0));
	
	// match leading ak4 and leading ak8 within 0.8 cone
	if(jetak4.DeltaR(jetak8) > 0.8) continue;
	
	//after match apply jetid on leading ak4
	if (chfrac->at(0) < 0.1) continue;   // jet id                                                                                                                       
	if (nhfrac->at(0) > 0.8) continue;   // jet id                                                                                                                        
	if (jetpt->at(0)  < 100.) continue;  // jet1 > 100 GeV                                                                                                             
	if (jetpt->at(0)  < *j1pt) continue;
	if (sample != 4 and jmdphi < 0.5) continue; // deltaPhi cut                                                                                                       

	// jet met dphi     
	float deltaPhi = 0.;	    
	if(fabs(pfmetphi-jetak8.Phi()) > TMath::Pi())
	  deltaPhi = fabs(2*TMath::Pi() - fabs(pfmetphi-jetak8.Phi()));
	else
	  deltaPhi = fabs(pfmetphi-jetak8.Phi());
	
	if(sample != 4 and deltaPhi < 0.5) continue; // deltaPhi cut                                                                                               	

	// no overlap between b-jet and v-jet
	if (sample == 7 || sample == 8){ 
	  int nbjets = 0;
	  for(size_t ijet = 0 ; ijet < jetbtag->size(); ijet++){
	    jetak4.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
	    if(jetak4.DeltaR(jetak8) < 0.8) continue;
	    if(jetbtag->at(ijet) > 0.89){
	      nbjets++;
	    }
	  }
	  if(nbjets < 1) continue;
	}

	// split among resonant and non resonant wrt gen level
	if(resonantSelection != 0 and isMC){
	  TLorentzVector Wboson4V;
	  if(boostedJetBosonpt->size() > 0 and boostedJetBosonpt->at(0) > 0){
	    Wboson4V.SetPtEtaPhiM(boostedJetBosonpt->at(0),boostedJetBosoneta->at(0),boostedJetBosonphi->at(0),boostedJetBosonm->at(0));
	    if(jetak8.DeltaR(Wboson4V) > 0.4 and resonantSelection == 1)
	      continue;
	    else if(jetak8.DeltaR(Wboson4V) < 0.4 and resonantSelection == 2)
	      continue;
	  }
	  else if(resonantSelection == 1)
	    continue;
	}

	// category 2 means HP mono-V
	if(category == 2 and (prunedJetm->at(0) > prunedMassMin and prunedJetm->at(0) < prunedMassMax) and boostedJettau2->at(0)/boostedJettau1->at(0) < tau2tau1)
	  goodMonoV   = true;
	// category 3 means LP mono-V
	else if(category == 3 and (prunedJetm->at(0) > prunedMassMin and prunedJetm->at(0) < prunedMassMax) and 
		(boostedJettau2->at(0)/boostedJettau1->at(0) > tau2tau1 and boostedJettau2->at(0)/boostedJettau1->at(0) < tau2tau1LP))
	  goodMonoV   = true;
	// apply no pruned mass cut --> show full shapes
	else if(category == 4 and (prunedJetm->at(0) > 0 and prunedJetm->at(0) < 200))
	  goodMonoV   = true;
	// apply only n-subjettiness
	else if(category == 5 and boostedJettau2->at(0)/boostedJettau1->at(0) < tau2tau1)
	  goodMonoV   = true;
	// apply only pruned mass cut
	else if(category == 6 and (prunedJetm->at(0) > prunedMassMin and prunedJetm->at(0) < prunedMassMax))
	  goodMonoV   = true;
	
	if(not goodMonoV) continue;	
	//if(not isMC)
	// dump<< "event id "<<*event<<" run "<<*run<<"jet pt "<<jetpt->at(0)<<" dphi "<<jmdphi<<" boosted jet pt "<<boostedJetpt->at(0)<<" pruned mass "<<prunedJetm->at(0)<<" tau2tau1 "<<boostedJettau2->at(0)/boostedJettau1->at(0)<<" met "<<pfmet<<" \n";

      }
    }
               
    // fill 1D histogram
    double fillvar = 0;
    
    // fill the histograms
    for(auto hist : hist1D){

      TString name(hist->GetName());      
      if(name.Contains("met"))
	fillvar = pfmet;      
      else if(name.Contains("nvtx"))
	fillvar = *nvtx;
      else if(name.Contains("chfrac"))
	fillvar = chfrac->at(0);
      else if(name.Contains("nhfrac"))
	fillvar = nhfrac->at(0);
      else if(name.Contains("jetPt"))
	fillvar = jetpt->at(0);
      else if(name.Contains("boostedJetPt")){
	if(boostedJetpt->size() > 0)
	  fillvar = boostedJetpt->at(0);
	else
	  fillvar = 0.;
      }
      else if(name.Contains("mT")){
	float deltaPhi = fabs(jetphi->at(0)-pfmetphi);
	if(deltaPhi > TMath::Pi())
	  deltaPhi = fabs(2*TMath::Pi() - deltaPhi);
	fillvar = sqrt(2*jetpt->at(0)*pfmet*(1-cos(deltaPhi)));
      }   
      else if(name.Contains("njet"))
	fillvar = *njets;      
      else if(name.Contains("nbjet_hpt_loose")){
	int nbjet = 0;
	for(size_t iJet = 0; iJet < jetbtag->size(); iJet++){
	  if(jetbtag->at(iJet) > 0.605 and jetpt->at(iJet) > 30)
	    nbjet++;
	}
	fillvar = nbjet;
      }      
      else if(name.Contains("nbjet_hpt")){
	int nbjet = 0;
	for(size_t iJet = 0; iJet < jetbtag->size(); iJet++){
	  if(jetbtag->at(iJet) > 0.89 and jetpt->at(iJet) > 30)
	    nbjet++;
	}
	fillvar = nbjet;
      }      
      else if(name.Contains("nbjet"))
	fillvar = *nbjets;
      else if(name.Contains("bosonPt"))
	fillvar = bosonPt;    	
      else if(name.Contains("QGL_1")){
	if(jetpt->at(0) < 175.)
	  fillvar = jetQGL->at(0);
	else
	  fillvar = -1.;
      }
      else if(name.Contains("QGL_2")){
	if(jetpt->at(0) > 175. and jetpt->at(0) < ptJetMinAK8)
	  fillvar = jetQGL->at(0);
	else
	  fillvar = -1.;
      }
      else if(name.Contains("QGL_3")){
	if(jetpt->at(0) > ptJetMinAK8)
	  fillvar = jetQGL->at(0);
	else
	  fillvar = -1.;
      }
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
      else if(name.Contains("btag_max") or name.Contains("CSV_max")){
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
      else if(name.Contains("btag_min") or name.Contains("CSV_min")){
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
      else if(name.Contains("btag") or name.Contains("CSV")){
	if(jetbtag->size() > 0)
	  fillvar = jetbtag->at(0);
	else
	  fillvar = 0;
      }      
      
      
      // overflow bin
      if (fillvar >= hist->GetBinLowEdge(hist->GetNbinsX())+hist->GetBinWidth(hist->GetNbinsX())) 
	fillvar = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());
      
      // total event weight
      double evtwgt  = 1.0;
      Double_t puwgt = 0.;
      if (isMC and not reweightNVTX){
	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(*wgtpileup)*(btagw)*hltw*sfwgt*topptwgt*kwgt*hwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(*wgtpileup)*(btagw)*hltw*sfwgt*topptwgt*kwgt*hwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)
      }
      else if (isMC and reweightNVTX){
	if (*nvtx <= 35) 
	  puwgt = puhist->GetBinContent(*nvtx);
	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*topptwgt*sfwgt*kwgt*hwgt/(*wgtsum);
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*topptwgt*sfwgt*kwgt*hwgt/(*wgtsum);
      }
      if (!isMC && sample == 6) 
	evtwgt = sfwgt;
      else if (!isMC)
	evtwgt = hltw;

      hist->Fill(fillvar, evtwgt);
    }

    //fill 2D histo
    double fillvarX = 0;
    double fillvarY = 0;

    for(auto hist: hist2D){
      TString name(hist->GetName());
      if(name.Contains("met_jetPt")){ 
	fillvarX = pfmet;
	fillvarY = jetpt->at(0);
      }
      else if(name.Contains("met_bosonPt")){
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

      else if(name.Contains("met_njet")){
	fillvarX = pfmet;
	fillvarY = *njets;
      }

      else if(name.Contains("mT_njet")){
	float deltaPhi = fabs(jetphi->at(0)-pfmetphi);
        if(deltaPhi > TMath::Pi())
          deltaPhi = 2*TMath::Pi() - deltaPhi;
	fillvarX = sqrt(2*jetpt->at(0)*pfmet*(1-cos(deltaPhi)));
	fillvarY = *njets;
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
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(*wgtpileup)*(btagw)*hltw*sfwgt*topptwgt*kwgt*hwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)      
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(*wgtpileup)*(btagw)*hltw*sfwgt*topptwgt*kwgt*hwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)      
      }
      else if (isMC and reweightNVTX){
	if (*nvtx <= 35) 
	  puwgt = puhist->GetBinContent(*nvtx);
	if(XSEC != -1)
	  evtwgt = (XSEC)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*sfwgt*topptwgt*kwgt*hwgt/(*wgtsum);
	else
	  evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*hltw*sfwgt*topptwgt*kwgt*hwgt/(*wgtsum);
      }
      if (!isMC && sample == 6) 
	evtwgt = sfwgt;
      else if (!isMC)
	evtwgt = hltw;
      
      hist->Fill(fillvarX,fillvarY,evtwgt);            
    }
  }

  //  dump.close();
  sffile  ->Close();
  psffile ->Close();
  trefile ->Close();
  pufile  ->Close();
  trmfile ->Close();
  QGLReweight ->Close();
}

#endif
