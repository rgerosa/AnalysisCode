#ifndef MAKEHIST4_H
#define MAKEHIST4_H

#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"


// define binnings for the different observables
int   nbins_monoV     = 7;
float bins_monoV[]    = {200., 250., 300., 350., 400., 500., 600., 1000.};
int   nbins_monoJet   = 7;
float bins_monoJet[]  = {200., 250., 300., 350., 400., 500., 600., 1000.};

void makehist4(TTree* tree, /*input tree*/ 
	       TH1* hist, /* MET histogram */ 
	       bool isMC, 
	       int sample, 
	       int category,
	       bool reweightNVTX,
	       double scale, 
	       vector<TH1*> khists, 
	       TH1* rhist=NULL) {

  Double_t lumi = 2.11;
  
  // now pufile not necessary anymore since the info already in the tree
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
  TFile* psffile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/PhotonSFandEffandPurity_Lumi2p1fb_2211.root");
  TH2*  psfhist  = (TH2*)psffile->Get("PhotonSF");  
  TH2*  purhist  = (TH2*)psffile->Get("PhotonPurity");
  
  // trigger efficiency correction for single electron trigger
  TFile* trefile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF/leptonTrigsfs.root");
  TH2*   trehist = (TH2*)trefile->Get("hltel27_SF");
  
  // trigger efficiency for met trigger
  TFile* trmfile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF/mettrigSF.root");
  TH1*   trmhist = (TH1*)trefile->Get("mettrigSF");
  
  // histogram to be filled
  hist->Sumw2();
  
  // define branches
  TTreeReader myReader(tree);

  // general info
  TTreeReaderValue<unsigned int> run          (myReader,"run");
  TTreeReaderValue<unsigned int> nvtx         (myReader,"nvtx");
  TTreeReaderValue<double> xsec      (myReader,"xsec");
  TTreeReaderValue<double> wgt       (myReader,"wgt");
  TTreeReaderValue<double> wgtsum    (myReader,"wgtsum");
  TTreeReaderValue<double> wgtpileup (myReader,"wgtpileup");
  TTreeReaderValue<double> wgtbtag   (myReader,"wgtbtag");
  
  // trigger
  TTreeReaderValue<UChar_t> hltm (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hlte (myReader,"hltsingleel");
  TTreeReaderValue<UChar_t> hltp (myReader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp2 (myReader,"hltphoton175");
  TTreeReaderValue<UChar_t> fhbhe (myReader,"flaghbheloose");
  TTreeReaderValue<UChar_t> fhbiso (myReader,"flaghbheiso");
  
  std::string cscname;
  std::string feebname;
  
  if(isMC){
    cscname  = "flagcsctight";
    feebname = "flageebadsc";
  }
  else{
    cscname  = "flagcscnew";
    feebname = "flageescnew";
  }
  
  TTreeReaderValue<UChar_t> fcsc (myReader,cscname.c_str());
  TTreeReaderValue<UChar_t> feeb (myReader,feebname.c_str());
  
  TTreeReaderValue<unsigned int>    njets  (myReader,"njets");
  TTreeReaderValue<unsigned int>    nbjets (myReader,"nbjetslowpt");
  TTreeReaderValue<double> j1pt   (myReader,"leadingjetpt");
  
  TTreeReaderValue<std::vector<double> > jetpt  (myReader,"centraljetpt");
  TTreeReaderValue<std::vector<double> > chfrac (myReader,"centraljetCHfrac");
  TTreeReaderValue<std::vector<double> > nhfrac (myReader,"centraljetNHfrac");
  TTreeReaderValue<std::vector<double> > emfrac (myReader,"centraljetEMfrac");

  // AK8 jet
  TTreeReaderValue<std::vector<double> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<std::vector<double> > prunedJetm      (myReader,"prunedJetm");
  TTreeReaderValue<std::vector<double> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<std::vector<double> > boostedJettau1  (myReader,"boostedJettau1");

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


  TTreeReaderValue<int> el1pid (myReader,"el1pid");
  TTreeReaderValue<int> el2pid (myReader,"el2pid");
  TTreeReaderValue<int> el1id (myReader,"el1id");
  TTreeReaderValue<int> el2id (myReader,"el2id");
  TTreeReaderValue<double> el1pt (myReader,"el1pt");
  TTreeReaderValue<double> el2pt (myReader,"el2pt");
  TTreeReaderValue<double> el1eta (myReader,"el1eta");
  TTreeReaderValue<double> el2eta (myReader,"el2eta");
  
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

    // check trigger
    Double_t hlt = 0.0;
    if (sample == 0 || sample == 1 || sample == 2) 
      hlt = *hltm;
    else if (sample == 3 || sample == 4)                
      hlt = *hlte;
    else if (sample == 5 || sample == 6)            
      hlt = *hltp;

    // check both photon triggers
    if ((sample == 5 || sample == 6) && *hltp2 > 0)      
      hlt = *hltp2;
    
    // check dphi jet-met
    Double_t jmdphi = 0.0;
    if (sample == 0 || sample == 1 || sample == 2) 
      jmdphi = fabs(*jmmdphi);
    else if (sample == 3 || sample == 4)                
      jmdphi = fabs(*jemdphi);
    else if (sample == 5 || sample == 6)                
      jmdphi = fabs(*jpmdphi);
    
    //set met
    Double_t met = 0.0;
    if (sample == 0 || sample == 1 || sample == 2) 
      met = *mmet;
    else if (sample == 3 || sample == 4)                
      met = *emet;
    else if (sample == 5 || sample == 6)                
      met = *pmet;
    
    // set zpt in case of Zsamples
    Double_t zpt = 0.0;
    if (sample == 1) 
      zpt = *zmmpt;
    else if (sample == 3) 
      zpt = zeept[0];
    else if (sample == 5) 
      zpt = *phpt;

    // set lepton info
    Int_t    id1   = 0;
    Int_t    id2   = 0;
    Double_t pt1   = 0.0;
    Double_t pt2   = 0.0;
    Double_t eta1  = 0.0;
    Double_t eta2  = 0.0;

    if (sample == 1 || sample == 2) {
      id1  = *mu1id;
      id2  = *mu2id;
      pt1  = *mu1pt;
      pt2  = *mu2pt;
      eta1 = fabs(*mu1eta);
      eta2 = fabs(*mu2eta);
    }
    else if (sample == 3 || sample == 4) {
      id1  = *el1id;
      id2  = *el2id;
      pt1  = *el1pt;
      pt2  = *el2pt;
      eta1 = fabs(*el1eta);
      eta2 = fabs(*el2eta);
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

    // scale factor for leptons
    TH2* sflhist = NULL;
    TH2* sfthist = NULL;
    
    if (sample == 1 || sample == 2) {
	sflhist = msflhist;
	sfthist = msfthist;
    }
    if (sample == 3 || sample == 4) {
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
    if (isMC && trehist && sample == 4) {
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
    if (isMC && trmhist && (sample == 0 || sample == 1 || sample == 2)) {
      sfwgt *= trmhist->GetBinContent(trmhist->FindBin(met));
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
     
      // common selections
      if (hlt  == 0) continue; // trigger
      //      if (*fhbhe == 0 || *fhbiso == 0 || *fcsc == 0 || *feeb == 0) continue; // met filters
      if (chfrac->size() == 0 || nhfrac->size() == 0 or jetpt->size() == 0) continue;
      if (chfrac->at(0) < 0.1) continue; // jet id
      if (nhfrac->at(0) > 0.8) continue; // jet id
      if (*njets  < 1) continue; //Njets > 1
      if (*nbjets > 0) continue; // bjets
      if (jetpt->at(0)  < 100.) continue; // jet1 > 100 GeV
      if (jetpt->at(0)  < *j1pt) continue;
      if (jmdphi < 0.5) continue; // deltaPhi cut

      // control regions with leptons or photons
      if (sample == 1 && *mu1pid == *mu2pid) continue;
      if (sample == 3 && *el1pid == *el2pid) continue;
      if ((sample == 5 || sample == 6) && *phpt < 175.) continue;
      if ((sample == 5 || sample == 6) && fabs(*pheta) > 1.4442) continue;
      if (sample == 4 && met < 50.) continue;

      if (met < 200.) continue;

      if(category != 0){
	
	// boosted category and monojet
	// So far MonoJet = jet pt < 200 + jetPt > 200 && pruned mass < 40
	bool goodMonoJet = false;
	bool goodMonoV   = false;

	if(category == 1){
	  if(boostedJetpt->size()  == 0)
	    goodMonoJet = true;
	  if(boostedJetpt->size() > 0){
	    if(boostedJetpt->at(0) < 200)
	      goodMonoJet = true;
	    else{
	      if(prunedJetm->at(0)   < 40)
		goodMonoJet= true;
	    }
	  }
	}

	if(category == 1 and not goodMonoJet) continue;

	if(category == 2){
	  
	  if(boostedJetpt->size()  > 0){
	    if(boostedJetpt->at(0) > 200){
	      if(sample <= 6 and (prunedJetm->at(0) > 65 and prunedJetm->at(0) < 105) and boostedJettau2->at(0)/boostedJettau1->at(0) < 0.45)
		goodMonoV   = true;
	      if(sample >= 7 and (prunedJetm->at(0) > 40 and prunedJetm->at(0) < 65) and boostedJettau2->at(0)/boostedJettau1->at(0) < 0.45)
		goodMonoV   = true;
	    }
	  }	  
	}	      
      
	if(category == 2 and not goodMonoV) continue;
      }
      // fill met histogram
      double fillvar = met;
      if (fillvar >= hist->GetBinLowEdge(hist->GetNbinsX())+hist->GetBinWidth(hist->GetNbinsX())) 
	fillvar = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());
      
      // total event weight
      double evtwgt = 1.0;
      Double_t puwgt = 0.;
      if (isMC and not reweightNVTX) 
	evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(*wgtpileup)*sfwgt*rwgt*kwgt/(*wgtsum); //(xsec, scale, lumi, wgt, pileup, sf, rw, kw, wgtsum)
      else if (isMC and reweightNVTX){
	if (*nvtx <= 35) 
	  puwgt = puhist->GetBinContent(*nvtx);
	evtwgt = (*xsec)*(scale)*(lumi)*(*wgt)*(puwgt)*sfwgt*rwgt*kwgt/(*wgtsum);
      }
      if (!isMC && sample == 6) 
	evtwgt = sfwgt;
      
      hist->Fill(fillvar, evtwgt);
  }
  
  sffile  ->Close();
  psffile ->Close();
  trefile ->Close();
  
}

#endif

