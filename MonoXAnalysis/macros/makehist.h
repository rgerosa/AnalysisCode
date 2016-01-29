#ifndef MAKEHIST4_H
#define MAKEHIST4_H

#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"
#include "TString.h"

using namespace std;


// define binnings for the different observables
vector<float> bins_monoV = {200., 250., 300., 350., 400., 500., 600., 1000.};
vector<float> bins_monoJ = {200., 250., 300., 350., 400., 500., 600., 1000.};
//vector<float> bins_monoJ = {200, 230, 260, 290, 320, 350, 390, 430, 470, 510, 550, 590, 640, 690, 740, 790, 840, 900, 960, 1020, 1090};

const float ptMax     = 1000.;
const float tau2tau1  = 0.6;

void makehist4(TTree* tree, /*input tree*/ 
	       vector<TH1*> hist1D, /* set of 1D histogram */ 
	       vector<TH2*> hist2D, /* set of 2D histogram */ 
	       bool isMC, 
	       int sample, 
	       int category,
	       bool isWJet,
	       double scale,
	       double lumi,	       
	       vector<TH1*> khists, 
	       string sysName,
	       bool reweightNVTX = true,
	       TH1* rhist = NULL) {

  // in case you want to weight the NVTX distribution
  TFile* pufile;
  TH1*   puhist = NULL;
  if(reweightNVTX){
    pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt.root");
    puhist = (TH1*) pufile->Get("puhist");
  }
    
  // Lepton ID scale factor from tag and probe: muons, electrons 
  TFile* sffile  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF/leptonIDsfs.root");
  
  TH2*  msflhist = (TH2*)sffile->Get("muon_loose_SF");
  TH2*  msfthist = (TH2*)sffile->Get("muon_tight_SF");
  
  TH2*  esflhist = (TH2*)sffile->Get("electron_veto_SF");
  TH2*  esfthist = (TH2*)sffile->Get("electron_tight_SF");
  
  // Photon ID scale factor from tag and probe
  TFile* psffile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/PhotonSFandEffandPurity_Lumi2p6fb_2301.root");
  TH2*  psfhist  = (TH2*)psffile->Get("PhotonSF");  
  TH2*  purhist  = (TH2*)psffile->Get("PhotonPurity");
  
  // trigger efficiency correction for single electron trigger
  TFile* trefile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF/leptonTrigsfs.root");
  TH2*   trehist = (TH2*)trefile->Get("hltel27_SF");
  
  // trigger efficiency for met trigger
  TFile* trmfile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF/mettrigSF.root");
  TH1*   trmhist = (TH1*)trefile->Get("mettrigSF");
  
  // histogram to be filled
  for(size_t ihist  = 0 ; ihist < hist1D.size(); ihist++)
    hist1D.at(ihist)->Sumw2();
  for(size_t ihist  = 0 ; ihist < hist2D.size(); ihist++)
    hist2D.at(ihist)->Sumw2();
  
  // define branches
  TTreeReader myReader(tree);

  // general info
  TTreeReaderValue<unsigned int> run    (myReader,"run");
  TTreeReaderValue<unsigned int> nvtx   (myReader,"nvtx");
  TTreeReaderValue<double> xsec      (myReader,"xsec");
  TTreeReaderValue<double> wgt       (myReader,"wgt");

  // create dummys for data
  string wgtname;
  string wgtpileupname;
  string btagname;
  if(isMC){
    wgtname = "wgtsum";
    wgtpileupname = "wgtpileup";
    if(sysName == "btagUp")
      btagname = "wgtbtagUp";
    else if(sysName == "btagDown")
      btagname = "wgtbtagDown";
    else
      btagname = "wgtbtag";
  }
  else{
    wgtname  = "wgt";
    btagname = "wgt";
    wgtpileupname = "wgt";
  }

  TTreeReaderValue<double> wgtsum    (myReader,wgtname.c_str());
  TTreeReaderValue<double> wgtpileup (myReader,wgtpileupname.c_str());


  TTreeReaderValue<double> wgtbtag   (myReader,btagname.c_str());
  
  // trigger
  TTreeReaderValue<UChar_t> hltm (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hlte (myReader,"hltsingleel");
  TTreeReaderValue<UChar_t> hltp (myReader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp2 (myReader,"hltphoton175");
  TTreeReaderValue<UChar_t> fhbhe (myReader,"flaghbheloose");
  TTreeReaderValue<UChar_t> fhbiso (myReader,"flaghbheiso");
  
  string cscname;
  string feebname;
  string fmutrkname;
  string fbtrkname;

  if(isMC){
    cscname  = "flagcsctight";
    feebname = "flageebadsc";
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
  TTreeReaderValue<UChar_t> fbtrk (myReader,fbtrkname.c_str());

  TTreeReaderValue<unsigned int>    njets  (myReader,"njets");
  TTreeReaderValue<unsigned int>    nbjets (myReader,"nbjetslowpt");
  TTreeReaderValue<double> j1pt   (myReader,"leadingjetpt");
  
  TTreeReaderValue<vector<double> > jetpt  (myReader,"centraljetpt");
  TTreeReaderValue<vector<double> > jeteta (myReader,"centraljeteta");
  TTreeReaderValue<vector<double> > jetphi (myReader,"centraljetphi");
  TTreeReaderValue<vector<double> > jetm   (myReader,"centraljetm");
  TTreeReaderValue<vector<double> > chfrac (myReader,"centraljetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac (myReader,"centraljetNHfrac");
  TTreeReaderValue<vector<double> > emfrac (myReader,"centraljetEMfrac");

  // AK8 jet
  TTreeReaderValue<vector<double> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<double> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<double> > boostedJetphi   (myReader,"boostedJetphi");
  TTreeReaderValue<vector<double> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<double> > prunedJetm      (myReader,"prunedJetm");
  TTreeReaderValue<vector<double> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<double> > boostedJettau1  (myReader,"boostedJettau1");

  // met
  TTreeReaderValue<double> met (myReader,"t1pfmet");
  TTreeReaderValue<double> metphi (myReader,"t1pfmetphi");
  
  TTreeReaderValue<double> mmet (myReader,"t1mumet");
  TTreeReaderValue<double> mmetphi (myReader,"t1mumetphi");
  
  TTreeReaderValue<double> emet (myReader,"t1elmet");
  TTreeReaderValue<double> emetphi (myReader,"t1elmetphi");
  
  TTreeReaderValue<double> pmet (myReader,"t1phmet");
  TTreeReaderValue<double> pmetphi (myReader,"t1phmetphi");
  
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

  // loop on events
  while(myReader.Next()){
    // check trigger depending on the sample
    Double_t hlt = 0.0;
    if (sample == 0 || sample == 1 || sample == 2 || sample == 7)  hlt = *hltm;
    else if (sample == 3 || sample == 4 || sample == 8)  hlt = *hlte;
    else if (sample == 5 || sample == 6)                 hlt = *hltp;

    // check both photon triggers
    if ((sample == 5 || sample == 6) && *hltp2 > 0)  hlt = *hltp2;
    
    // check dphi jet-met
    Double_t jmdphi = 0.0;
    
    if (sample == 0 || sample == 1 || sample == 2 || sample == 7) jmdphi = fabs(*jmmdphi);
    else if (sample == 3 || sample == 4 || sample == 8)           jmdphi = fabs(*jemdphi);
    else if (sample == 5 || sample == 6)           jmdphi = fabs(*jpmdphi);

    //set met
    Double_t pfmet = 0.0;
    Double_t pfmetphi = 0.0;
    if (sample == 0 || sample == 1 || sample == 2 || sample == 7){ pfmet = *mmet; pfmetphi = *mmetphi;}
    else if (sample == 3 || sample == 4 || sample == 8)          { pfmet = *emet; pfmetphi = *emetphi;}
    else if (sample == 5 || sample == 6)          { pfmet = *pmet; pfmetphi = *pmetphi;}
    else if (sample == 7 and *hlte)               { pfmet = *emet; pfmetphi = *emetphi;}
    else if (sample == 7 and not *hlte)           { pfmet = *mmet; pfmetphi = *mmetphi;}

    // set zpt in case of Zsamples
    Double_t zpt = 0.0;
    if (sample == 1)      zpt = *zmmpt;
    else if (sample == 3) zpt = zeept[0];
    else if (sample == 5) zpt = *phpt;

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
    
    if (pt1 >= ptMax) 
      pt1 = 999.0;
    if (pt2 >= ptMax) 
      pt2 = 999.0;

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
      if (pt1 > 175. && id1 == 1) {
	sfwgt *= (1.0 - purhist->GetBinContent(purhist->FindBin(pt1, eta1)));
	}
    }

    // met trigger scale factor
    if (isMC && trmhist && (sample == 0 || sample == 1 || sample == 2 || sample == 7)) {
      sfwgt *= trmhist->GetBinContent(trmhist->FindBin(pfmet));
    }


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
      
    // Gen level info --> NLO re-weight
    Double_t kvar = *wzpt;
    if (*wzpt < 100. ) 
      *wzpt = 100.;
    if (*wzpt > 1000.) 
      *wzpt = 999.;
    
    Double_t rwgt = 1.0;
    if (rhist) 
      rwgt = rhist->GetBinContent(rhist->FindBin(*wzpt));
    
    Double_t kwgt = 1.0;
    for (unsigned i = 0; i < khists.size(); i++) {
      if (isMC && khists[i]) kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(*wzpt));
    }
    
    // common selections (trigger + MET filter + Njets + b-veto)
    if (hlt  == 0) continue; // trigger
    if (*fhbhe == 0 || *fhbiso == 0 || *fcsc == 0 || *feeb == 0) continue;
    if (not isMC){
      if(*fmutrk == 0 || *fbtrk == 0) continue; // met filters
    }
    if (*njets  < 1) continue; //Njets > 1
    if (*nbjets > 0 and sample != 7 and sample != 8) continue; // bjets
    // control regions with leptons or photons (2mu or 2ele, or one high pt photon)
    if (sample == 1 && *mu1pid == *mu2pid) continue;
    if (sample == 3 && *el1pid == *el2pid) continue;
    // take only leading tight leptons in the lepton control sample
    if ((sample == 1 || sample == 3) && not (id1 == 1 or id2 == 1)) continue;
    if ((sample == 2 || sample == 4) && id1 != 1) continue;
    //
    if ((sample == 5 || sample == 6) && *phpt < 175.) continue;
    if ((sample == 5 || sample == 6) && fabs(*pheta) > 1.4442) continue;
    if (sample == 4 && *met < 50.) continue;
    if ((sample == 7 || sample == 8) && *nbjets < 1) continue;
    if (sample == 7 || sample == 8){ // Z-veto in case of more than one lepton    
      if(pt1 <=0) continue;
      if(abs(pid1) == 13 && pt1 < 20. ) continue;
      if(abs(pid1) == 11 && pt1 < 40. ) continue;
      if(pt1 > 0 && pt2 > 0){
	TLorentzVector lep1, lep2;
	lep1.SetPtEtaPhiM(pt1,eta1,phi1,0.);
	lep2.SetPtEtaPhiM(pt2,eta2,phi2,0.);
	if((lep1+lep2).M() > 81. && (lep1+lep2).M() < 101.) continue;
	if(fabs(pid1) != fabs(pid2)) continue;
      }
    }

    // met selection
    if (pfmet < 200.) continue;


    if(category == 0){ // inclusive mono-jet analysis

      if (chfrac->size() == 0 || nhfrac->size() == 0 or jetpt->size() == 0) continue; // at least one leading jet
      if (chfrac->at(0) < 0.1) continue;   // jet id
      if (nhfrac->at(0) > 0.8) continue;   // jet id
      if (jetpt->at(0)  < 100.) continue;  // jet1 > 100 GeV
      if (jetpt->at(0)  < *j1pt) continue; 
      if (jmdphi < 0.5) continue; // deltaPhi cut

    }
    else{

      
      // boosted category and monojet
      // So far MonoJet = jet pt < 200 + jetPt > 200 && pruned mass < 40
      bool goodMonoJet = false;
      bool goodMonoV   = false;

      if(category == 1){ // mono jet + V-jet veto

	if (chfrac->size() == 0 || nhfrac->size() == 0 or jetpt->size() == 0) continue; // at least one leading jet                                                         
	if (chfrac->at(0) < 0.1) continue;   // jet id                                                                                                                       
	if (nhfrac->at(0) > 0.8) continue;   // jet id                                                                                                                        
	if (jetpt->at(0)  < 100.) continue;  // jet1 > 100 GeV                                                                                                             
	if (jetpt->at(0)  < *j1pt) continue;
	if (jmdphi < 0.5) continue; // deltaPhi cut                                                                                                          

	if(boostedJetpt->size()  == 0) // no boosted jets (AK8 pT > 200 GeV)
	  goodMonoJet = true;
	if(boostedJetpt->size() > 0){ // in case one boosted jet
	  if(boostedJetpt->at(0) < 200) // check pT
	    goodMonoJet = true;
	  else{ // if high pT check pruned mass
	    if(prunedJetm->at(0)   < 40)
	      goodMonoJet= true;
	  }
	}      	
	if(not goodMonoJet) continue;
      }

      else if(category >= 2){

	if(chfrac->size() == 0 || nhfrac->size() == 0 or jetpt->size() == 0) continue; 	
	if(boostedJetpt->size() == 0) continue;
	if(boostedJetpt->at(0) < 200) continue;

	TLorentzVector jetak4, jetak8;
	jetak4.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
	jetak8.SetPtEtaPhiM(boostedJetpt->at(0),boostedJeteta->at(0),boostedJetphi->at(0),boostedJetm->at(0));

	// match leading ak4 and leading ak8 within 0.4 cone
	if(jetak4.DeltaR(jetak8) > 0.4) continue;
	
	//after match apply jetid on leading ak4
	if (chfrac->at(0) < 0.1) continue;   // jet id                                                                                                                       
	if (nhfrac->at(0) > 0.8) continue;   // jet id                                                                                                                        
	if (jetpt->at(0)  < 100.) continue;  // jet1 > 100 GeV                                                                                                             
	if (jetpt->at(0)  < *j1pt) continue;

	// jet met dphi
	TLorentzVector met4V;
	met4V.SetPtEtaPhiM(pfmet,0.,pfmetphi,0);
	if (met4V.DeltaPhi(jetak8) < 0.5) continue; // deltaPhi cut                                                                                                          
	
	// category 2 means HP mono-V
	if(category == 2      and (prunedJetm->at(0) > 65 and prunedJetm->at(0) < 105) and boostedJettau2->at(0)/boostedJettau1->at(0) < tau2tau1)
	  goodMonoV   = true;
	// category 3 means LP mono-V
	else if(category == 3 and (prunedJetm->at(0) > 65 and prunedJetm->at(0) < 105) and 
		boostedJettau2->at(0)/boostedJettau1->at(0) > tau2tau1 and boostedJettau2->at(0)/boostedJettau1->at(0) < 0.75)
	  goodMonoV   = true;
	// category 4 means HP mono-V sideband
	else if(category == 4 and (prunedJetm->at(0) > 40 and prunedJetm->at(0) < 65) and 
		boostedJettau2->at(0)/boostedJettau1->at(0) < tau2tau1)
	  goodMonoV   = true;
	// category 5 means LP mono-V sideband
	else if(category == 5 and (prunedJetm->at(0) > 40 and prunedJetm->at(0) < 65) and 
		boostedJettau2->at(0)/boostedJettau1->at(0) > tau2tau1 and boostedJettau2->at(0)/boostedJettau1->at(0) < 0.75)
	  goodMonoV   = true;
	// category 6 means mono-V sideband
	else if(category == 6 and (prunedJetm->at(0) > 40 and prunedJetm->at(0) < 65))
	  goodMonoV   = true;	
	if(not goodMonoV) continue;
      }
    }
               
    // fill 1D histogram
    double fillvar = 0;
    double btagw = *wgtbtag;
    if(*wgtbtag > 2 || *wgtbtag < 0)
      btagw = 1;

    for(auto hist : hist1D){
      TString name(hist->GetName());
      
      if(name.Contains("met") or name.Contains("pfmet") or name.Contains("MET"))
	fillvar = pfmet;
      else if(name.Contains("jetpt") or name.Contains("jetPT"))
	fillvar = jetpt->at(0);
      else if(name.Contains("zpt") or name.Contains("wzpt"))
	fillvar = zpt;
      else if(name.Contains("mT")){
	float deltaPhi = fabs(jetphi->at(0)-pfmetphi);
	if(deltaPhi > TMath::Pi())
	  deltaPhi = 2*TMath::Pi() - deltaPhi;
	fillvar = sqrt(2*jetpt->at(0)*pfmet*(1-cos(deltaPhi)));
      }    
	// overflow bin
      if (fillvar >= hist->GetBinLowEdge(hist->GetNbinsX())+hist->GetBinWidth(hist->GetNbinsX())) 
	fillvar = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());

      // total event weight
      double evtwgt  = 1.0;
      Double_t puwgt = 0.;
      if (isMC and not reweightNVTX) 
	evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(*wgtpileup)*(btagw)*sfwgt*rwgt*kwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)
      else if (isMC and reweightNVTX){
	if (*nvtx <= 35) 
	  puwgt = puhist->GetBinContent(*nvtx);
	evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*sfwgt*rwgt*kwgt/(*wgtsum);
      }
      if (!isMC && sample == 6) 
	evtwgt = sfwgt;

      hist->Fill(fillvar, evtwgt);
    }

    //fill 2D histo
    double fillvarX = 0;
    double fillvarY = 0;

    for(auto hist: hist2D){
      TString name(hist->GetName());

      if(name.Contains("met_jetpt") or name.Contains("pfmet_jetpt")){ 
	fillvarX = pfmet;
	fillvarY = jetpt->at(0);
      }
      else if(name.Contains("met_zpt") or name.Contains("pfmet_zpt") or name.Contains("met_wzpt") or name.Contains("pfmet_wzpt")){
	fillvarX = pfmet;
	fillvarY = zpt;
      }
      else if(name.Contains("met_mT") or name.Contains("pfmet_mT")){	
	fillvarX = pfmet;
	float deltaPhi = fabs(jetphi->at(0)-pfmetphi);
	if(deltaPhi > TMath::Pi())
	  deltaPhi = 2*TMath::Pi() - deltaPhi;
	fillvarY = sqrt(2*jetpt->at(0)*pfmet*(1-cos(deltaPhi)));
      }    
      // overflow bin
      if (fillvarX >= hist->GetXaxis()->GetBinLowEdge(hist->GetNbinsX())+hist->GetXaxis()->GetBinWidth(hist->GetNbinsX())) 
	fillvarX = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());
      if (fillvarY >= hist->GetYaxis()->GetBinLowEdge(hist->GetNbinsY())+hist->GetYaxis()->GetBinWidth(hist->GetNbinsY())) 
	fillvarY = hist->GetYaxis()->GetBinCenter(hist->GetNbinsY());

      // total event weight
      double evtwgt  = 1.0;
      Double_t puwgt = 0.;

      if (isMC and not reweightNVTX)
	evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(*wgtpileup)*(btagw)*sfwgt*rwgt*kwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)      
      else if (isMC and reweightNVTX){
	if (*nvtx <= 35) 
	  puwgt = puhist->GetBinContent(*nvtx);
	evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*(btagw)*sfwgt*rwgt*kwgt/(*wgtsum);
      }
      if (!isMC && sample == 6) 
	evtwgt = sfwgt;
      
      hist->Fill(fillvarX,fillvarY,evtwgt);      
      
    }
  }

  sffile  ->Close();
  psffile ->Close();
  trefile ->Close();
  
}

#endif

