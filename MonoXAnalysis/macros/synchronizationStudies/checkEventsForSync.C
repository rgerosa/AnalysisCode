#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TTreeReader.h"

void checkEventsForSync(string inputFile, string outputDir){

  system(("mkdir -p "+outputDir).c_str());
  system(("rm "+outputDir+"/*").c_str());

  TFile* input = TFile::Open(inputFile.c_str());
  TTree* tree  = (TTree*) input->Get("tree/tree");
  
  // declare branches
  TTreeReader myReader(tree);
  TTreeReaderValue<unsigned int> run    (myReader,"run");
  TTreeReaderValue<unsigned int> lumi   (myReader,"lumi");
  TTreeReaderValue<unsigned int> event  (myReader,"event");

  // triggers 
  TTreeReaderValue<UChar_t> hltm90     (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm120    (myReader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120  (myReader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170  (myReader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300  (myReader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90   (myReader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hlte       (myReader,"hltsingleel");
  TTreeReaderValue<UChar_t> hltp165    (myReader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175    (myReader,"hltphoton175");

  // met filters
  TTreeReaderValue<UChar_t> fhbhe  (myReader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (myReader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (myReader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (myReader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (myReader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (myReader,"flaggoodvertices");

  // njets
  TTreeReaderValue<unsigned int> njets  (myReader,"njets");
  TTreeReaderValue<unsigned int> nbjets (myReader,"nbjets");
  TTreeReaderValue<unsigned int> nbjetslowpt (myReader,"nbjetslowpt");

  // AK4 jets
  TTreeReaderValue<vector<double> > jetpt   (myReader,"centraljetpt");
  TTreeReaderValue<vector<double> > jeteta  (myReader,"centraljeteta");
  TTreeReaderValue<vector<double> > jetphi  (myReader,"centraljetphi");
  TTreeReaderValue<vector<double> > jetm    (myReader,"centraljetm");
  TTreeReaderValue<vector<double> > jetbtag (myReader,"centraljetbtag");
  TTreeReaderValue<vector<double> > chfrac  (myReader,"centraljetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac  (myReader,"centraljetNHfrac");
  TTreeReaderValue<vector<double> > emfrac  (myReader,"centraljetEMfrac");

  TTreeReaderValue<vector<double> > jetpt_f   (myReader,"forwardjetpt");
  TTreeReaderValue<vector<double> > jeteta_f  (myReader,"forwardjeteta");
  TTreeReaderValue<vector<double> > jetphi_f  (myReader,"forwardjetphi");
  TTreeReaderValue<vector<double> > jetm_f    (myReader,"forwardjetm");
  TTreeReaderValue<vector<double> > chfrac_f  (myReader,"forwardjetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac_f  (myReader,"forwardjetNHfrac");
  TTreeReaderValue<vector<double> > emfrac_f  (myReader,"forwardjetEMfrac");
  // Boosted jets
  TTreeReaderValue<vector<double> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<double> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<double> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<double> > prunedJetm      (myReader,"prunedJetm");
  TTreeReaderValue<vector<double> > prunedJetm_v2   (myReader,"prunedJetm_v2");
  TTreeReaderValue<vector<double> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<double> > boostedJettau1  (myReader,"boostedJettau1");

  // MET
  TTreeReaderValue<double> met (myReader,"t1pfmet");
  TTreeReaderValue<double> mmet        (myReader,"t1mumet");
  TTreeReaderValue<double> emet        (myReader,"t1elmet");
  TTreeReaderValue<double> pmet        (myReader,"t1phmet");

  // Dphi
  TTreeReaderValue<double> jmmdphi (myReader,"incjetmumetdphimin4");
  TTreeReaderValue<double> jemdphi (myReader,"incjetelmetdphimin4");
  TTreeReaderValue<double> jpmdphi (myReader,"incjetphmetdphimin4");

  TTreeReaderValue<unsigned int> nmuons     (myReader,"nmuons");
  TTreeReaderValue<unsigned int> nelectrons (myReader,"nelectrons");
  TTreeReaderValue<unsigned int> ntaus      (myReader,"ntaus");
  TTreeReaderValue<unsigned int> nphotons   (myReader,"nphotons");

  // Muon information
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

  // Electron information
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

  // Photon information
  TTreeReaderValue<int> phidm (myReader,"phidm");
  TTreeReaderValue<double> phpt (myReader,"phpt");
  TTreeReaderValue<double> pheta (myReader,"pheta");
  TTreeReaderValue<double> phphi (myReader,"phphi");

  TTreeReaderValue<double> zmass (myReader,"zmass");
  TTreeReaderValue<double> zeemass (myReader,"zeemass");
  TTreeReaderValue<double> zmmpt (myReader,"zpt");
  TTreeReaderValue<double> zeept (myReader,"zeept");
  TTreeReaderValue<double> zeeeta (myReader,"zeeeta");
  TTreeReaderValue<double> zmmeta (myReader,"zeta");

  // output text files for the comparison
  ofstream leptonVeto((outputDir+"/leptonVeto_SR.txt").c_str());
  ofstream leptonPhtonVeto((outputDir+"/leptonPhtonVeto_SR.txt").c_str());
  ofstream leptonPhtonTauVeto((outputDir+"/leptonPhtonTauVeto_SR.txt").c_str());
  ofstream ak4JetSelections((outputDir+"/ak4JetSelections_SR.txt").c_str());
  ofstream metSelections((outputDir+"/metSelections_SR.txt").c_str());
  ofstream btagVetoSelections((outputDir+"/btagVeto_SR.txt").c_str());
  ofstream VtaggingSelections((outputDir+"/Vtagging_SR.txt").c_str());
  ofstream VBFSelections((outputDir+"/VBF_SR.txt").c_str());

  long int n_total        = 0;
  long int n_muonVeto     = 0;
  long int n_electronVeto = 0;
  long int n_photonVeto   = 0;
  long int n_tauVeto      = 0;
  long int n_bjetVeto     = 0;
  long int n_jetpt        = 0;
  long int n_jetid        = 0;
  long int n_jetdphi      = 0;
  long int n_metcut       = 0;
  long int n_ak8pt        = 0;
  long int n_ak8tau2tau1  = 0;
  long int n_ak8mpruned   = 0;
  long int n_metHard      = 0;
  long int n_jetVBF       = 0;
  long int n_emVBF      = 0;
  long int n_detaVBF      = 0;
  long int n_mjjVBF       = 0;
  
  /// event loop
  while(myReader.Next()){
    
    n_total++;

    if(*nmuons != 0) continue;
    n_muonVeto++;

    if(*nelectrons != 0) continue;
    n_electronVeto++;
    leptonVeto << *run << " "<<*lumi<<" "<<*event<<"\n";

    if(*nphotons != 0) continue;    
    n_photonVeto++;
    leptonPhtonVeto << *run << " "<<*lumi<<" "<<*event<<"\n";

    if(*ntaus != 0) continue;
    n_tauVeto++;
    leptonPhtonTauVeto << *run << " "<<*lumi<<" "<<*event<<"\n";
        
    if(*nbjetslowpt > 0) continue;
    n_bjetVeto++;
    btagVetoSelections << *run << " "<<*lumi<<" "<<*event<<"\n";
    
    if(jetpt->size() <= 0) continue;
    if(jetpt->at(0) < 100) continue;
    n_jetpt++;
    
    if(chfrac->at(0) < 0.1) continue;
    if(nhfrac->at(0) > 0.8) continue;
    n_jetid++;

    ak4JetSelections << *run << " "<<*lumi<<" "<<*event<<"\n";

    if(*jmmdphi < 0.5) continue;
    n_jetdphi++;

    if(*met < 200) continue;
    n_metcut++;
    metSelections << *run << " "<<*lumi<<" "<<*event<<"\n";

    /// VBF block
    if(jetpt->size() + jetpt_f->size() >= 2){
      TLorentzVector leadingJet;
      leadingJet.SetPtEtaPhiM(0,0,0,0);
      bool isLeadingCentral;
      TLorentzVector subleadingJet;
      subleadingJet.SetPtEtaPhiM(0,0,0,0);
      bool isSubLeadingCentral;
      
      for(int ijet = 0; ijet<jetpt->size(); ijet++){
	if(jetpt->at(ijet) > leadingJet.Pt()){
	  subleadingJet = leadingJet; 
	  leadingJet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));	  
	  isLeadingCentral = true;
	}
	else if(jetpt->at(ijet) > subleadingJet.Pt()){
	  subleadingJet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
	  isSubLeadingCentral = true;
	}
      }
      
      for(int ijet = 0; ijet<jetpt_f->size(); ijet++){
	if(jetpt_f->at(ijet) > leadingJet.Pt()){
	  subleadingJet = leadingJet; 
	  leadingJet.SetPtEtaPhiM(jetpt_f->at(ijet),jeteta_f->at(ijet),jetphi_f->at(ijet),jetm_f->at(ijet));	  
	  isLeadingCentral = false;
	}
	else if(jetpt_f->at(ijet) > subleadingJet.Pt()){
	  subleadingJet.SetPtEtaPhiM(jetpt_f->at(ijet),jeteta_f->at(ijet),jetphi_f->at(ijet),jetm_f->at(ijet));
	  isSubLeadingCentral = false;
	}
      }
      
      if(leadingJet.Pt() > 100 && fabs(leadingJet.Eta()) < 4.7 && fabs(subleadingJet.Eta()) < 4.7 && subleadingJet.Pt() > 40){
	n_jetVBF++;
	if(leadingJet.Eta()*subleadingJet.Eta() < 0){
	  n_emVBF++;
	  if(fabs(leadingJet.Eta()-subleadingJet.Eta()) > 3.5 and leadingJet.Eta()*subleadingJet.Eta() < 0){
	    n_detaVBF++;
	    TLorentzVector sum;
	    sum = leadingJet+subleadingJet;
	    if(sum.M() > 500){
	      n_mjjVBF++;
	      VBFSelections << *run << "  "<<*lumi<<" "<<*event<<"\n";
	    }
	  }
	}
      }
    }
    if(boostedJetpt->size() <= 0) continue;
    if(boostedJetpt->at(0) < 250 ) continue;
    if(fabs(boostedJeteta->at(0)) > 2.4 ) continue;
    n_ak8pt++;

    if(boostedJettau2->at(0)/boostedJettau1->at(0) > 0.6) continue;
    n_ak8tau2tau1++;

    if(prunedJetm_v2->at(0) < 65 or prunedJetm_v2->at(0) > 105 ) continue;
    n_ak8mpruned++;
    VtaggingSelections << *run << " "<<*lumi<<" "<<*event<<"\n";
   
    if(*met<250) continue;
    n_metHard++;
    
  }
  
  leptonVeto.close();
  leptonPhtonVeto.close();
  leptonPhtonTauVeto.close();
  ak4JetSelections.close();
  metSelections.close();
  btagVetoSelections.close();
  VtaggingSelections.close();
  VBFSelections.close();

  cout<<"############################"<<endl;
  cout<<"## Event report in the SR ##"<<endl;
  cout<<"############################"<<endl;
  cout<<"total event = "<<n_total<<endl;
  cout<<"muon veto event = "<<n_muonVeto<<endl;
  cout<<"electron veto event = "<<n_electronVeto<<endl;
  cout<<"photon veto event = "<<n_photonVeto<<endl;
  cout<<"tau veto event = "<<n_tauVeto<<endl;
  cout<<"bjet veto event = "<<n_bjetVeto<<endl;
  cout<<"jetpt event = "<<n_jetpt<<endl;
  cout<<"jetid event = "<<n_jetid<<endl;
  cout<<"jetphi event = "<<n_jetdphi<<endl;
  cout<<"met cut event = "<<n_metcut<<endl;
  cout<<"ak8 pt event = "<<n_ak8pt<<endl;
  cout<<"ak8 tau2tau1 event = "<<n_ak8tau2tau1<<endl;
  cout<<"ak8 mpruned event = "<<n_ak8mpruned<<endl;
  cout<<"met hard cut = "<<n_metHard<<endl;
  cout<<"njet VBF = "<<n_jetVBF<<endl;
  cout<<"n opposite hem VBF = "<<n_emVBF<<endl;
  cout<<"deta VBF = "<<n_detaVBF<<endl;
  cout<<"mjj VBF = "<<n_mjjVBF<<endl;

  // Double muon CR
  long int n_muons = 0;
  long int n_muonTag     = 0;
  n_electronVeto = 0;
  n_photonVeto   = 0;
  n_tauVeto      = 0;
  n_bjetVeto     = 0;
  n_jetpt        = 0;
  n_jetid        = 0;
  n_jetdphi      = 0;
  n_metcut       = 0;
  n_ak8pt        = 0;
  n_ak8tau2tau1  = 0;
  n_ak8mpruned   = 0;
  n_metHard      = 0;
  n_jetVBF       = 0;
  n_detaVBF      = 0;
  n_mjjVBF       = 0;
  n_emVBF        = 0;

  ofstream leptonTag_ZM((outputDir+"/leptonTag_ZM.txt").c_str());
  ofstream leptonVeto_ZM((outputDir+"/leptonVeto_ZM.txt").c_str());
  ofstream bjetVeto_ZM((outputDir+"/bjetVeto_ZM.txt").c_str());
  ofstream boostedJets_ZM((outputDir+"/boostedJets_ZM.txt").c_str());

  myReader.SetEntry(0);
  while(myReader.Next()){

    if(*nmuons != 2) continue;
    n_muons++;
    if(not (((*mu1pt > 20 and *mu1id == 1) || (*mu2pt > 20 and *mu2id == 1)) and *zmass > 60 and *zmass < 120 and *mu1pid != *mu2pid )) continue; 
    n_muonTag++;

    leptonTag_ZM << *run << " " <<*lumi<<" "<<*event<<"\n";

    if(*nelectrons != 0) continue;
    n_electronVeto++;

    leptonVeto_ZM << *run << " " <<*lumi<<" "<<*event<<"\n";

    if(*nphotons != 0) continue;    
    n_photonVeto++;
    if(*ntaus != 0) continue;
    n_tauVeto++;
        
    if(*nbjetslowpt > 0) continue;
    n_bjetVeto++;
    bjetVeto_ZM << *run << " " <<*lumi<<" "<<*event<<"\n";




    if(jetpt->size() <= 0) continue;
    if(jetpt->at(0) < 100) continue;
    n_jetpt++;
    
    if(chfrac->at(0) < 0.1) continue;
    if(nhfrac->at(0) > 0.8) continue;
    n_jetid++;

    if(*jmmdphi < 0.5) continue;
    n_jetdphi++;
    if(*mmet < 200) continue;
    n_metcut++;


    /// VBF block
    if(jetpt->size() + jetpt_f->size() >= 2){
      TLorentzVector leadingJet;
      leadingJet.SetPtEtaPhiM(0,0,0,0);
      bool isLeadingCentral;
      TLorentzVector subleadingJet;
      subleadingJet.SetPtEtaPhiM(0,0,0,0);
      bool isSubLeadingCentral;
      
      for(int ijet = 0; ijet<jetpt->size(); ijet++){
	if(jetpt->at(ijet) > leadingJet.Pt()){
	  subleadingJet = leadingJet; 
	  leadingJet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));	  
	  isLeadingCentral = true;
	}
	else if(jetpt->at(ijet) > subleadingJet.Pt()){
	  subleadingJet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
	  isSubLeadingCentral = true;
	}
      }
      
      for(int ijet = 0; ijet<jetpt_f->size(); ijet++){
	if(jetpt_f->at(ijet) > leadingJet.Pt()){
	  subleadingJet = leadingJet; 
	  leadingJet.SetPtEtaPhiM(jetpt_f->at(ijet),jeteta_f->at(ijet),jetphi_f->at(ijet),jetm_f->at(ijet));	  
	  isLeadingCentral = false;
	}
	else if(jetpt_f->at(ijet) > subleadingJet.Pt()){
	  subleadingJet.SetPtEtaPhiM(jetpt_f->at(ijet),jeteta_f->at(ijet),jetphi_f->at(ijet),jetm_f->at(ijet));
	  isSubLeadingCentral = false;
	}
      }
      
      if(leadingJet.Pt() > 100 && fabs(leadingJet.Eta()) < 4.7 && fabs(subleadingJet.Eta()) < 4.7 && subleadingJet.Pt() > 40){
	n_jetVBF++;
	if(leadingJet.Eta()*subleadingJet.Eta() < 0){
	  n_emVBF++;
	  if(fabs(leadingJet.Eta()-subleadingJet.Eta()) > 3.5 and leadingJet.Eta()*subleadingJet.Eta() < 0){
	    n_detaVBF++;
	    TLorentzVector sum;
	    sum = leadingJet+subleadingJet;
	    if(sum.M() > 500)
	      n_mjjVBF++;
	  }
	}
      }
    }

    if(boostedJetpt->size() <= 0) continue;
    if(boostedJetpt->at(0) < 250 ) continue;
    if(fabs(boostedJeteta->at(0)) > 2.4 ) continue;
    boostedJets_ZM << *run << " " <<*lumi<<" "<<*event<<"\n";
    n_ak8pt++;
    if(boostedJettau2->at(0)/boostedJettau1->at(0) > 0.6) continue;
    n_ak8tau2tau1++;
    if(prunedJetm_v2->at(0) < 65 or prunedJetm_v2->at(0) > 105 ) continue;
    n_ak8mpruned++;
   
    if(*mmet<250) continue;
      n_metHard++;
    
  }

  leptonTag_ZM.close();
  leptonVeto_ZM.close();
  boostedJets_ZM.close();
  bjetVeto_ZM.close();

  cout<<"################################"<<endl;
  cout<<"## Event report in the DiMuon ##"<<endl;
  cout<<"################################"<<endl;
  cout<<"total event = "<<n_total<<endl;
  cout<<"n muon event = "<<n_muons<<endl;
  cout<<"muon tag event = "<<n_muonTag<<endl;
  cout<<"electron veto event = "<<n_electronVeto<<endl;
  cout<<"photon veto event = "<<n_photonVeto<<endl;
  cout<<"tau veto event = "<<n_tauVeto<<endl;
  cout<<"bjet veto event = "<<n_bjetVeto<<endl;
  cout<<"jetpt event = "<<n_jetpt<<endl;
  cout<<"jetid event = "<<n_jetid<<endl;
  cout<<"jetphi event = "<<n_jetdphi<<endl;
  cout<<"met cut event = "<<n_metcut<<endl;
  cout<<"ak8 pt event = "<<n_ak8pt<<endl;
  cout<<"ak8 tau2tau1 event = "<<n_ak8tau2tau1<<endl;
  cout<<"ak8 mpruned event = "<<n_ak8mpruned<<endl;
  cout<<"met hard cut = "<<n_metHard<<endl;
  cout<<"njet VBF = "<<n_jetVBF<<endl;
  cout<<"n opposite hem VBF = "<<n_emVBF<<endl;
  cout<<"deta VBF = "<<n_detaVBF<<endl;
  cout<<"mjj VBF = "<<n_mjjVBF<<endl;

  // Single muon CR                                                                                                                                                             
  n_muons        = 0;
  n_muonTag      = 0;
  n_electronVeto = 0;
  n_photonVeto   = 0;
  n_tauVeto      = 0;
  n_bjetVeto     = 0;
  n_jetpt        = 0;
  n_jetid        = 0;
  n_jetdphi      = 0;
  n_metcut       = 0;
  n_ak8pt        = 0;
  n_ak8tau2tau1  = 0;
  n_ak8mpruned   = 0;
  n_metHard      = 0;
  n_jetVBF       = 0;
  n_detaVBF      = 0;
  n_mjjVBF       = 0;
  n_emVBF        = 0;

  myReader.SetEntry(0);
  while(myReader.Next()){

    if(*nmuons != 1) continue;
    n_muons++;
    if(not (*mu1pt > 20 and *mu1id == 1)) continue;
    n_muonTag++;

    if(*nelectrons != 0) continue;
    n_electronVeto++;
    if(*nphotons != 0) continue;    
    n_photonVeto++;
    if(*ntaus != 0) continue;
    n_tauVeto++;
        
    if(*nbjetslowpt > 0) continue;
    n_bjetVeto++;

    if(jetpt->size() <= 0) continue;
    if(jetpt->at(0) < 100) continue;
    n_jetpt++;
    
    if(chfrac->at(0) < 0.1) continue;
    if(nhfrac->at(0) > 0.8) continue;
    n_jetid++;

    if(*jmmdphi < 0.5) continue;
    n_jetdphi++;
    if(*mmet < 200) continue;
    n_metcut++;
    /// VBF block
    if(jetpt->size() + jetpt_f->size() >= 2){
      TLorentzVector leadingJet;
      leadingJet.SetPtEtaPhiM(0,0,0,0);
      bool isLeadingCentral;
      TLorentzVector subleadingJet;
      subleadingJet.SetPtEtaPhiM(0,0,0,0);
      bool isSubLeadingCentral;
      
      for(int ijet = 0; ijet<jetpt->size(); ijet++){
	if(jetpt->at(ijet) > leadingJet.Pt()){
	  subleadingJet = leadingJet; 
	  leadingJet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));	  
	  isLeadingCentral = true;
	}
	else if(jetpt->at(ijet) > subleadingJet.Pt()){
	  subleadingJet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
	  isSubLeadingCentral = true;
	}
      }
      
      for(int ijet = 0; ijet<jetpt_f->size(); ijet++){
	if(jetpt_f->at(ijet) > leadingJet.Pt()){
	  subleadingJet = leadingJet; 
	  leadingJet.SetPtEtaPhiM(jetpt_f->at(ijet),jeteta_f->at(ijet),jetphi_f->at(ijet),jetm_f->at(ijet));	  
	  isLeadingCentral = false;
	}
	else if(jetpt_f->at(ijet) > subleadingJet.Pt()){
	  subleadingJet.SetPtEtaPhiM(jetpt_f->at(ijet),jeteta_f->at(ijet),jetphi_f->at(ijet),jetm_f->at(ijet));
	  isSubLeadingCentral = false;
	}
      }
      
      if(leadingJet.Pt() > 100 && fabs(leadingJet.Eta()) < 4.7 && fabs(subleadingJet.Eta()) < 4.7 && subleadingJet.Pt() > 40){
	n_jetVBF++;
	if(leadingJet.Eta()*subleadingJet.Eta() < 0){
	  n_emVBF++;
	  if(fabs(leadingJet.Eta()-subleadingJet.Eta()) > 3.5 and leadingJet.Eta()*subleadingJet.Eta() < 0){
	    n_detaVBF++;
	    TLorentzVector sum;
	    sum = leadingJet+subleadingJet;
	    if(sum.M() > 500)
	      n_mjjVBF++;
	  }
	}
      }
    }



    if(boostedJetpt->size() <= 0) continue;
    if(boostedJetpt->at(0) < 250 ) continue;
    if(fabs(boostedJeteta->at(0)) > 2.4 ) continue;
    n_ak8pt++;
    if(boostedJettau2->at(0)/boostedJettau1->at(0) > 0.6) continue;
    n_ak8tau2tau1++;
    if(prunedJetm_v2->at(0) < 65 or prunedJetm_v2->at(0) > 105 ) continue;
    n_ak8mpruned++;
   
    if(*mmet<250) continue; 
      n_metHard++;
    
  }

  cout<<"################################"<<endl;
  cout<<"## Event report in the SingleMuon ##"<<endl;
  cout<<"################################"<<endl;
  cout<<"total event = "<<n_total<<endl;
  cout<<"n muon event = "<<n_muons<<endl;
  cout<<"muon tag event = "<<n_muonTag<<endl;
  cout<<"electron veto event = "<<n_electronVeto<<endl;
  cout<<"photon veto event = "<<n_photonVeto<<endl;
  cout<<"tau veto event = "<<n_tauVeto<<endl;
  cout<<"bjet veto event = "<<n_bjetVeto<<endl;
  cout<<"jetpt event = "<<n_jetpt<<endl;
  cout<<"jetid event = "<<n_jetid<<endl;
  cout<<"jetphi event = "<<n_jetdphi<<endl;
  cout<<"met cut event = "<<n_metcut<<endl;
  cout<<"ak8 pt event = "<<n_ak8pt<<endl;
  cout<<"ak8 tau2tau1 event = "<<n_ak8tau2tau1<<endl;
  cout<<"ak8 mpruned event = "<<n_ak8mpruned<<endl;
  cout<<"met hard cut = "<<n_metHard<<endl;
  cout<<"njet VBF = "<<n_jetVBF<<endl;
  cout<<"n opposite hem VBF = "<<n_emVBF<<endl;
  cout<<"deta VBF = "<<n_detaVBF<<endl;
  cout<<"mjj VBF = "<<n_mjjVBF<<endl;


  ////////////////////////////////////

  // Double electron CR
  n_muonVeto = 0;
  long int n_electron    = 0;
  long int n_electronTag = 0;
  n_photonVeto   = 0;
  n_tauVeto      = 0;
  n_bjetVeto     = 0;
  n_jetpt        = 0;
  n_jetid        = 0;
  n_jetdphi      = 0;
  n_metcut       = 0;
  n_ak8pt        = 0;
  n_ak8tau2tau1  = 0;
  n_ak8mpruned   = 0;
  n_metHard      = 0;
  n_jetVBF       = 0;
  n_detaVBF      = 0;
  n_mjjVBF       = 0;
  n_emVBF        = 0;

  myReader.SetEntry(0);
  while(myReader.Next()){

    if(*nelectrons != 2) continue;
    n_electron++;
    if(not (((*el1pt > 40 and *el1id == 1)  or (*el2pt > 40 and *el2id == 1)) and *zeemass > 60 and *zeemass < 120 and  *el1pid != *el2pid)) continue; 
    n_electronTag++;

    if(*nmuons != 0) continue;
    n_muonVeto++;
    if(*nphotons != 0) continue;    
    n_photonVeto++;
    if(*ntaus != 0) continue;
    n_tauVeto++;
        
    if(*nbjetslowpt > 0) continue;
    n_bjetVeto++;

    if(jetpt->size() <= 0) continue;
    if(jetpt->at(0) < 100) continue;
    n_jetpt++;
    
    if(chfrac->at(0) < 0.1) continue;
    if(nhfrac->at(0) > 0.8) continue;
    n_jetid++;

    if(*jemdphi < 0.5) continue;
    n_jetdphi++;
    if(*emet < 200) continue;
    n_metcut++;

    /// VBF block
    if(jetpt->size() + jetpt_f->size() >= 2){
      TLorentzVector leadingJet;
      leadingJet.SetPtEtaPhiM(0,0,0,0);
      bool isLeadingCentral;
      TLorentzVector subleadingJet;
      subleadingJet.SetPtEtaPhiM(0,0,0,0);
      bool isSubLeadingCentral;
      
      for(int ijet = 0; ijet<jetpt->size(); ijet++){
	if(jetpt->at(ijet) > leadingJet.Pt()){
	  subleadingJet = leadingJet; 
	  leadingJet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));	  
	  isLeadingCentral = true;
	}
	else if(jetpt->at(ijet) > subleadingJet.Pt()){
	  subleadingJet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
	  isSubLeadingCentral = true;
	}
      }
      
      for(int ijet = 0; ijet<jetpt_f->size(); ijet++){
	if(jetpt_f->at(ijet) > leadingJet.Pt()){
	  subleadingJet = leadingJet; 
	  leadingJet.SetPtEtaPhiM(jetpt_f->at(ijet),jeteta_f->at(ijet),jetphi_f->at(ijet),jetm_f->at(ijet));	  
	  isLeadingCentral = false;
	}
	else if(jetpt_f->at(ijet) > subleadingJet.Pt()){
	  subleadingJet.SetPtEtaPhiM(jetpt_f->at(ijet),jeteta_f->at(ijet),jetphi_f->at(ijet),jetm_f->at(ijet));
	  isSubLeadingCentral = false;
	}
      }
      
      if(leadingJet.Pt() > 100 && fabs(leadingJet.Eta()) < 4.7 && fabs(subleadingJet.Eta()) < 4.7 && subleadingJet.Pt() > 40){
	n_jetVBF++;
	if(leadingJet.Eta()*subleadingJet.Eta() < 0){
	  n_emVBF++;
	  if(fabs(leadingJet.Eta()-subleadingJet.Eta()) > 3.5 and leadingJet.Eta()*subleadingJet.Eta() < 0){
	    n_detaVBF++;
	    TLorentzVector sum;
	    sum = leadingJet+subleadingJet;
	    if(sum.M() > 500)
	      n_mjjVBF++;
	  }
	}
      }
    }


    if(boostedJetpt->size() <= 0) continue;
    if(boostedJetpt->at(0) < 250 ) continue;
    if(fabs(boostedJeteta->at(0)) > 2.4 ) continue;
    n_ak8pt++;
    if(boostedJettau2->at(0)/boostedJettau1->at(0) > 0.6) continue;
    n_ak8tau2tau1++;
    if(prunedJetm_v2->at(0) < 65 or prunedJetm_v2->at(0) > 105 ) continue;
    n_ak8mpruned++;
   
    if(*emet<250) continue;
      n_metHard++;
    
  }

  cout<<"################################"<<endl;
  cout<<"## Event report in the DiElectron ##"<<endl;
  cout<<"################################"<<endl;
  cout<<"total event = "<<n_total<<endl;
  cout<<"n muon event = "<<n_electron<<endl;
  cout<<"muon tag event = "<<n_electronTag<<endl;
  cout<<"electron veto event = "<<n_muonVeto<<endl;
  cout<<"photon veto event = "<<n_photonVeto<<endl;
  cout<<"tau veto event = "<<n_tauVeto<<endl;
  cout<<"bjet veto event = "<<n_bjetVeto<<endl;
  cout<<"jetpt event = "<<n_jetpt<<endl;
  cout<<"jetid event = "<<n_jetid<<endl;
  cout<<"jetphi event = "<<n_jetdphi<<endl;
  cout<<"met cut event = "<<n_metcut<<endl;
  cout<<"ak8 pt event = "<<n_ak8pt<<endl;
  cout<<"ak8 tau2tau1 event = "<<n_ak8tau2tau1<<endl;
  cout<<"ak8 mpruned event = "<<n_ak8mpruned<<endl;
  cout<<"met hard cut = "<<n_metHard<<endl;
  cout<<"njet VBF = "<<n_jetVBF<<endl;
  cout<<"n opposite hem VBF = "<<n_emVBF<<endl;
  cout<<"deta VBF = "<<n_detaVBF<<endl;
  cout<<"mjj VBF = "<<n_mjjVBF<<endl;

  // Single electron CR                                                                                                                                                        
  n_electron     = 0;
  n_electronTag  = 0;
  n_muonVeto     = 0;
  n_photonVeto   = 0;
  n_tauVeto      = 0;
  n_bjetVeto     = 0;
  n_jetpt        = 0;
  n_jetid        = 0;
  n_jetdphi      = 0;
  n_metcut       = 0;
  n_ak8pt        = 0;
  n_ak8tau2tau1  = 0;
  n_ak8mpruned   = 0;
  n_metHard      = 0;
  n_jetVBF       = 0;
  n_detaVBF      = 0;
  n_mjjVBF       = 0;
  n_emVBF        = 0;

  myReader.SetEntry(0);
  while(myReader.Next()){

    if(*nelectrons != 1) continue;
    n_electron++;
    if(not (*el1pt > 40 and *el1id == 1)) continue; 
    n_electronTag++;

    if(*nmuons != 0) continue;
    n_muonVeto++;
    if(*nphotons != 0) continue;    
    n_photonVeto++;
    if(*ntaus != 0) continue;
    n_tauVeto++;
        
    if(*nbjetslowpt > 0) continue;
    n_bjetVeto++;

    if(jetpt->size() <= 0) continue;
    if(jetpt->at(0) < 100) continue;
    n_jetpt++;
    
    if(chfrac->at(0) < 0.1) continue;
    if(nhfrac->at(0) > 0.8) continue;
    n_jetid++;

    if(*jemdphi < 0.5) continue;
    n_jetdphi++;
    if(*emet < 200) continue;
    n_metcut++;

    /// VBF block
    if(jetpt->size() + jetpt_f->size() >= 2){
      TLorentzVector leadingJet;
      leadingJet.SetPtEtaPhiM(0,0,0,0);
      bool isLeadingCentral;
      TLorentzVector subleadingJet;
      subleadingJet.SetPtEtaPhiM(0,0,0,0);
      bool isSubLeadingCentral;
      
      for(int ijet = 0; ijet<jetpt->size(); ijet++){
	if(jetpt->at(ijet) > leadingJet.Pt()){
	  subleadingJet = leadingJet; 
	  leadingJet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));	  
	  isLeadingCentral = true;
	}
	else if(jetpt->at(ijet) > subleadingJet.Pt()){
	  subleadingJet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
	  isSubLeadingCentral = true;
	}
      }
      
      for(int ijet = 0; ijet<jetpt_f->size(); ijet++){
	if(jetpt_f->at(ijet) > leadingJet.Pt()){
	  subleadingJet = leadingJet; 
	  leadingJet.SetPtEtaPhiM(jetpt_f->at(ijet),jeteta_f->at(ijet),jetphi_f->at(ijet),jetm_f->at(ijet));	  
	  isLeadingCentral = false;
	}
	else if(jetpt_f->at(ijet) > subleadingJet.Pt()){
	  subleadingJet.SetPtEtaPhiM(jetpt_f->at(ijet),jeteta_f->at(ijet),jetphi_f->at(ijet),jetm_f->at(ijet));
	  isSubLeadingCentral = false;
	}
      }
      
      if(leadingJet.Pt() > 100 && fabs(leadingJet.Eta()) < 4.7 && fabs(subleadingJet.Eta()) < 4.7 && subleadingJet.Pt() > 40){
	n_jetVBF++;
	if(leadingJet.Eta()*subleadingJet.Eta() < 0){
	  n_emVBF++;
	  if(fabs(leadingJet.Eta()-subleadingJet.Eta()) > 3.5 and leadingJet.Eta()*subleadingJet.Eta() < 0){
	    n_detaVBF++;
	    TLorentzVector sum;
	    sum = leadingJet+subleadingJet;
	    if(sum.M() > 500)
	      n_mjjVBF++;
	  }
	}
      }
    }


    if(boostedJetpt->size() <= 0) continue;
    if(boostedJetpt->at(0) < 250 ) continue;
    if(fabs(boostedJeteta->at(0)) > 2.4 ) continue;
    n_ak8pt++;
    if(boostedJettau2->at(0)/boostedJettau1->at(0) > 0.6) continue;
    n_ak8tau2tau1++;
    if(prunedJetm_v2->at(0) < 65 or prunedJetm_v2->at(0) > 105 ) continue;
    n_ak8mpruned++;
   
    if(*emet<250) continue;
      n_metHard++;
    
  }

  cout<<"################################"<<endl;
  cout<<"## Event report in the SingleElectron ##"<<endl;
  cout<<"################################"<<endl;
  cout<<"total event = "<<n_total<<endl;
  cout<<"n electron event = "<<n_electron<<endl;
  cout<<"electron tag event = "<<n_electronTag<<endl;
  cout<<"electron veto event = "<<n_muonVeto<<endl;
  cout<<"photon veto event = "<<n_photonVeto<<endl;
  cout<<"tau veto event = "<<n_tauVeto<<endl;
  cout<<"bjet veto event = "<<n_bjetVeto<<endl;
  cout<<"jetpt event = "<<n_jetpt<<endl;
  cout<<"jetid event = "<<n_jetid<<endl;
  cout<<"jetphi event = "<<n_jetdphi<<endl;
  cout<<"met cut event = "<<n_metcut<<endl;
  cout<<"ak8 pt event = "<<n_ak8pt<<endl;
  cout<<"ak8 tau2tau1 event = "<<n_ak8tau2tau1<<endl;
  cout<<"ak8 mpruned event = "<<n_ak8mpruned<<endl;
  cout<<"met hard cut = "<<n_metHard<<endl;
  cout<<"njet VBF = "<<n_jetVBF<<endl;
  cout<<"n opposite hem VBF = "<<n_emVBF<<endl;
  cout<<"deta VBF = "<<n_detaVBF<<endl;
  cout<<"mjj VBF = "<<n_mjjVBF<<endl;

  // single photon
  n_electronVeto = 0;
  n_muonVeto     = 0;
  long int n_photonTag = 0;
  n_tauVeto      = 0;
  n_bjetVeto     = 0;
  n_jetpt        = 0;
  n_jetid        = 0;
  n_jetdphi      = 0;
  n_metcut       = 0;
  n_ak8pt        = 0;
  n_ak8tau2tau1  = 0;
  n_ak8mpruned   = 0;
  n_metHard      = 0;
  n_jetVBF       = 0;
  n_detaVBF      = 0;
  n_mjjVBF       = 0;
  n_emVBF        = 0;
  
  ofstream photonTagSelections((outputDir+"/PhotonTag_GJ.txt").c_str());
  ofstream photonBvetoSelections((outputDir+"/PhotonBveto_GJ.txt").c_str());

  myReader.SetEntry(0);
  while(myReader.Next()){

    if(*nphotons != 1) continue;    
    if(*phpt < 175 || fabs(*pheta) > 1.442 or *phidm != 1) continue;    
    n_photonTag++;

    photonTagSelections<< *run << " " <<*lumi<<" "<<*event<<"\n";

    if(*nelectrons != 0) continue;
    n_electronVeto++;
    if(*nmuons != 0) continue;
    n_muonVeto++;
    if(*ntaus != 0) continue;
    n_tauVeto++;
        
    if(*nbjetslowpt > 0) continue;
    n_bjetVeto++;

    photonBvetoSelections << *run << " " <<*lumi<<" "<<*event<<"\n";

    if(jetpt->size() <= 0) continue;
    if(jetpt->at(0) < 100) continue;
    n_jetpt++;
    
    if(chfrac->at(0) < 0.1) continue;
    if(nhfrac->at(0) > 0.8) continue;
    n_jetid++;

    if(*jpmdphi < 0.5) continue;
    n_jetdphi++;
    if(*pmet < 200) continue;
    n_metcut++;

    /// VBF block
    if(jetpt->size() + jetpt_f->size() >= 2){
      TLorentzVector leadingJet;
      leadingJet.SetPtEtaPhiM(0,0,0,0);
      bool isLeadingCentral;
      TLorentzVector subleadingJet;
      subleadingJet.SetPtEtaPhiM(0,0,0,0);
      bool isSubLeadingCentral;
      
      for(int ijet = 0; ijet<jetpt->size(); ijet++){
	if(jetpt->at(ijet) > leadingJet.Pt()){
	  subleadingJet = leadingJet; 
	  leadingJet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));	  
	  isLeadingCentral = true;
	}
	else if(jetpt->at(ijet) > subleadingJet.Pt()){
	  subleadingJet.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetm->at(ijet));
	  isSubLeadingCentral = true;
	}
      }
      
      for(int ijet = 0; ijet<jetpt_f->size(); ijet++){
	if(jetpt_f->at(ijet) > leadingJet.Pt()){
	  subleadingJet = leadingJet; 
	  leadingJet.SetPtEtaPhiM(jetpt_f->at(ijet),jeteta_f->at(ijet),jetphi_f->at(ijet),jetm_f->at(ijet));	  
	  isLeadingCentral = false;
	}
	else if(jetpt_f->at(ijet) > subleadingJet.Pt()){
	  subleadingJet.SetPtEtaPhiM(jetpt_f->at(ijet),jeteta_f->at(ijet),jetphi_f->at(ijet),jetm_f->at(ijet));
	  isSubLeadingCentral = false;
	}
      }
      
      if(leadingJet.Pt() > 100 && fabs(leadingJet.Eta()) < 4.7 && fabs(subleadingJet.Eta()) < 4.7 && subleadingJet.Pt() > 40){
	n_jetVBF++;
	if(leadingJet.Eta()*subleadingJet.Eta() < 0){
	  n_emVBF++;
	  if(fabs(leadingJet.Eta()-subleadingJet.Eta()) > 3.5 and leadingJet.Eta()*subleadingJet.Eta() < 0){
	    n_detaVBF++;
	    TLorentzVector sum;
	    sum = leadingJet+subleadingJet;
	    if(sum.M() > 500)
	      n_mjjVBF++;
	  }
	}
      }
    }

	

    if(boostedJetpt->size() <= 0) continue;
    if(boostedJetpt->at(0) < 250 ) continue;
    if(fabs(boostedJeteta->at(0)) > 2.4 ) continue;
    n_ak8pt++;
    if(boostedJettau2->at(0)/boostedJettau1->at(0) > 0.6) continue;
    n_ak8tau2tau1++;
    if(prunedJetm_v2->at(0) < 65 or prunedJetm_v2->at(0) > 105 ) continue;
    n_ak8mpruned++;
   
    if(*pmet<250) continue;
      n_metHard++;
    
  }

  photonTagSelections.close();
  photonBvetoSelections.close();

  cout<<"################################"<<endl;
  cout<<"## Event report in the SinglePhoton ##"<<endl;
  cout<<"################################"<<endl;
  cout<<"total event = "<<n_total<<endl;
  cout<<"n photon tag = "<<n_photonTag<<endl;
  cout<<"muon veto event = "<<n_muonVeto<<endl;
  cout<<"electron veto event = "<<n_electronVeto<<endl;
  cout<<"tau veto event = "<<n_tauVeto<<endl;
  cout<<"bjet veto event = "<<n_bjetVeto<<endl;
  cout<<"jetpt event = "<<n_jetpt<<endl;
  cout<<"jetid event = "<<n_jetid<<endl;
  cout<<"jetphi event = "<<n_jetdphi<<endl;
  cout<<"met cut event = "<<n_metcut<<endl;
  cout<<"ak8 pt event = "<<n_ak8pt<<endl;
  cout<<"ak8 tau2tau1 event = "<<n_ak8tau2tau1<<endl;
  cout<<"ak8 mpruned event = "<<n_ak8mpruned<<endl;
  cout<<"met hard cut = "<<n_metHard<<endl;
  cout<<"njet VBF = "<<n_jetVBF<<endl;
  cout<<"n opposite hem VBF = "<<n_emVBF<<endl;
  cout<<"deta VBF = "<<n_detaVBF<<endl;
  cout<<"mjj VBF = "<<n_mjjVBF<<endl;


}
