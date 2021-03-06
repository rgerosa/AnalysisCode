#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TTreeReader.h"

#include "../makeTemplates/histoUtils.h"

void checkEventsForSync(string inputFile, string outputDir, Sample controlRegion, Category category){

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
  TTreeReaderValue<vector<float> > jetpt   (myReader,"combinejetpt");
  TTreeReaderValue<vector<float> > jeteta  (myReader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetphi  (myReader,"combinejetphi");
  TTreeReaderValue<vector<float> > jetm    (myReader,"combinejetm");
  TTreeReaderValue<vector<float> > jetbtag (myReader,"combinejetbtag");
  TTreeReaderValue<vector<float> > chfrac  (myReader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac  (myReader,"combinejetNHfrac");
  TTreeReaderValue<vector<float> > emfrac  (myReader,"combinejetEMfrac");

  // Boosted jets
  TTreeReaderValue<vector<float> > boostedJetpt    (myReader,"boostedJetpt");
  TTreeReaderValue<vector<float> > boostedJeteta   (myReader,"boostedJeteta");
  TTreeReaderValue<vector<float> > boostedJetm     (myReader,"boostedJetm");
  TTreeReaderValue<vector<float> > prunedJetm      (myReader,"prunedJetm");
  TTreeReaderValue<vector<float> > prunedJetm_v2   (myReader,"prunedJetm_v2");
  TTreeReaderValue<vector<float> > boostedJettau2  (myReader,"boostedJettau2");
  TTreeReaderValue<vector<float> > boostedJettau1  (myReader,"boostedJettau1");

  // MET
  TTreeReaderValue<float> met (myReader,"t1pfmet");
  TTreeReaderValue<float> mmet        (myReader,"t1mumet");
  TTreeReaderValue<float> emet        (myReader,"t1elmet");
  TTreeReaderValue<float> pmet        (myReader,"t1phmet");

  // Dphi
  TTreeReaderValue<float> jmmdphi (myReader,"incjetmumetdphimin4");
  TTreeReaderValue<float> jemdphi (myReader,"incjetelmetdphimin4");
  TTreeReaderValue<float> jpmdphi (myReader,"incjetphmetdphimin4");
  TTreeReaderValue<unsigned int> nmuons     (myReader,"nmuons");
  TTreeReaderValue<unsigned int> nelectrons (myReader,"nelectrons");
  TTreeReaderValue<unsigned int> ntaus      (myReader,"ntausrawold");
  TTreeReaderValue<unsigned int> nphotons   (myReader,"nphotons");

  // Muon information
  TTreeReaderValue<int> mu1pid (myReader,"mu1pid");
  TTreeReaderValue<int> mu2pid (myReader,"mu2pid");
  TTreeReaderValue<int> mu1id (myReader,"mu1id");
  TTreeReaderValue<int> mu2id (myReader,"mu2id");
  TTreeReaderValue<float> mu1pt (myReader,"mu1pt");
  TTreeReaderValue<float> mu2pt (myReader,"mu2pt");
  TTreeReaderValue<float> mu1eta (myReader,"mu1eta");
  TTreeReaderValue<float> mu2eta (myReader,"mu2eta");
  TTreeReaderValue<float> mu1phi (myReader,"mu1phi");
  TTreeReaderValue<float> mu2phi (myReader,"mu2phi");

  // Electron information
  TTreeReaderValue<int> el1pid (myReader,"el1pid");
  TTreeReaderValue<int> el2pid (myReader,"el2pid");
  TTreeReaderValue<int> el1id (myReader,"el1id");
  TTreeReaderValue<int> el2id (myReader,"el2id");
  TTreeReaderValue<float> el1pt (myReader,"el1pt");
  TTreeReaderValue<float> el2pt (myReader,"el2pt");
  TTreeReaderValue<float> el1eta (myReader,"el1eta");
  TTreeReaderValue<float> el2eta (myReader,"el2eta");
  TTreeReaderValue<float> el1phi (myReader,"el1phi");
  TTreeReaderValue<float> el2phi (myReader,"el2phi");

  // Photon information
  TTreeReaderValue<int> phidm (myReader,"phidm");
  TTreeReaderValue<float> phpt (myReader,"phpt");
  TTreeReaderValue<float> pheta (myReader,"pheta");
  TTreeReaderValue<float> phphi (myReader,"phphi");

  TTreeReaderValue<float> zmass (myReader,"zmass");
  TTreeReaderValue<float> zeemass (myReader,"zeemass");
  TTreeReaderValue<float> zmmpt (myReader,"zpt");
  TTreeReaderValue<float> zeept (myReader,"zeept");
  TTreeReaderValue<float> zeeeta (myReader,"zeeeta");
  TTreeReaderValue<float> zmmeta (myReader,"zeta");

  TTreeReaderValue<float> wmt (myReader,"wmt");
  TTreeReaderValue<float> wemt (myReader,"wemt");

  TTreeReaderValue<float> wgt (myReader,"wgt");

  // output text files for the comparison
  string control = "";
  if(controlRegion == Sample::sig) control = "SR";
  else if(controlRegion == Sample::wmn) control = "wmn";
  else if(controlRegion == Sample::wen) control = "wen";
  else if(controlRegion == Sample::zmm) control = "zmm";
  else if(controlRegion == Sample::zee) control = "zee";
  else if(controlRegion == Sample::gam) control = "gam";
  ofstream leptonSelection((outputDir+"/leptonSelection_"+control+".txt").c_str());
  ofstream photonSelection((outputDir+"/photonSelection_"+control+".txt").c_str());
  ofstream tauSelection((outputDir+"/tauSelection_"+control+".txt").c_str());
  ofstream jetSelection((outputDir+"/jetSelection_"+control+".txt").c_str());
  ofstream metSelection((outputDir+"/metSelection_"+control+".txt").c_str());
  ofstream btagSelection((outputDir+"/btagSelection_"+control+".txt").c_str());
  ofstream VtaggingSelection((outputDir+"/VtaggingSelection_"+control+".txt").c_str());
  ofstream VBFSelection((outputDir+"/VBFSelection_"+control+".txt").c_str());
  
  long int n_total        = 0;
  long int n_muonLooseSelection  = 0;
  long int n_muonSelection     = 0;
  long int n_electronLooseSelection = 0;
  long int n_electronSelection = 0;
  long int n_photonSelection   = 0;
  long int n_tauSelection      = 0;
  long int n_bjetSelection     = 0;
  long int n_metmT        = 0;
  long int n_jetpt        = 0;
  long int n_jetid        = 0;
  long int n_jetdphi      = 0;
  long int n_metcut       = 0;
  long int n_ak8pt        = 0;
  long int n_ak8tau2tau1  = 0;
  long int n_ak8mpruned   = 0;
  long int n_metHard      = 0;
  long int n_emVBF        = 0;
  long int n_detaVBF      = 0;
  long int n_mjjVBF       = 0;

  double nwgt_total        = 0;
  double nwgt_muonLooseSelection  = 0;
  double nwgt_muonSelection     = 0;
  double nwgt_electronLooseSelection = 0;
  double nwgt_electronSelection = 0;
  double nwgt_photonSelection   = 0;
  double nwgt_tauSelection      = 0;
  double nwgt_bjetSelection     = 0;
  double nwgt_metmT        = 0;
  double nwgt_jetpt        = 0;
  double nwgt_jetid        = 0;
  double nwgt_jetdphi      = 0;
  double nwgt_metcut       = 0;
  double nwgt_ak8pt        = 0;
  double nwgt_ak8tau2tau1  = 0;
  double nwgt_ak8mpruned   = 0;
  double nwgt_metHard      = 0;
  double nwgt_emVBF      = 0;
  double nwgt_detaVBF      = 0;
  double nwgt_mjjVBF       = 0;

  /// event loop
  while(myReader.Next()){
    
    n_total++;
    nwgt_total += *wgt;

    // muon selection
    if((controlRegion == Sample::sig or controlRegion == Sample::wen or controlRegion == Sample::zee or controlRegion == Sample::gam) and *nmuons != 0) continue;
    else if(controlRegion == Sample::wmn){
      if(*nmuons != 1) continue;
      n_muonLooseSelection++;
      nwgt_muonLooseSelection += *wgt;
      if(*mu1pt < 20) continue;
      if(fabs(*mu1eta) > 2.4) continue;
      if(*mu1id != 1) continue;
    }
    else if(controlRegion == Sample::zmm){      
      if (*nmuons != 2 ) continue;
      n_muonLooseSelection++;
      nwgt_muonLooseSelection += *wgt;;
      if (fabs(*mu1eta) > 2.4) continue;
      if (fabs(*mu2eta) > 2.4) continue;
      if (*mu1pid == *mu2pid) continue;
      if (not ((*mu1pt > 20 and *mu1id == 1) or (*mu2pt > 20 and *mu2id == 1))) continue;
      if (*zmass < 60 or *zmass > 120) continue;
    }    
    if(controlRegion != Sample::wmn and controlRegion != Sample::zmm){
      n_muonLooseSelection++;
      nwgt_muonLooseSelection += *wgt;;
    }
    n_muonSelection++;
    nwgt_muonSelection += *wgt;

    // electron selection
    if((controlRegion == Sample::sig or controlRegion == Sample::wmn or controlRegion == Sample::zmm or controlRegion == Sample::gam) and *nelectrons != 0) continue;
    if(controlRegion == Sample::wen){
      if (*nelectrons != 1 ) continue;
      n_electronLooseSelection++;
      nwgt_electronLooseSelection += *wgt;;
       if(*el1pt < 40) continue;
      if(fabs(*el1eta) > 2.5) continue;
      if(*el1id != 1) continue;
    }
    else if(controlRegion == Sample::zee){      
      if (*nelectrons != 2 ) continue;
      n_electronLooseSelection++;
      nwgt_electronLooseSelection += *wgt;;
      if (fabs(*el1eta) > 2.5) continue;
      if (fabs(*el2eta) > 2.5) continue;
      if (*el1pid == *el2pid) continue;
      if (not ((*el1pt > 40 and *el1id == 1) or (*el2pt > 40 and *el2id == 1))) continue;
      if (*zeemass < 60 or *zeemass > 120) continue;
    }
    if(controlRegion != Sample::wen and controlRegion != Sample::zee){
      n_electronLooseSelection++;
      nwgt_electronLooseSelection += *wgt;;
    }

    n_electronSelection++;
    nwgt_electronSelection += *wgt;
    leptonSelection << *run << " "<<*lumi<<" "<<*event<<"\n";

    // photon selection
    if((controlRegion == Sample::sig or controlRegion == Sample::wmn or controlRegion == Sample::zmm or controlRegion == Sample::wen or controlRegion == Sample::zee) and *nphotons != 0) continue;
    else if(controlRegion == Sample::gam){
      if(*nphotons != 1) continue;
      if(fabs(*pheta) > 1.444) continue;
      if(*phpt < 175) continue;
    }

    n_photonSelection++;
    nwgt_photonSelection += *wgt;
    photonSelection << *run << " "<<*lumi<<" "<<*event<<"\n";
    
    ///////////
    if(*ntaus != 0) continue;
    n_tauSelection++;
    nwgt_tauSelection += *wgt;
    tauSelection << *run << " "<<*lumi<<" "<<*event<<"\n";
        
    ///////////
    if(*nbjetslowpt > 0) continue;
    n_bjetSelection++;
    nwgt_bjetSelection += *wgt;
    btagSelection << *run << " "<<*lumi<<" "<<*event<<"\n";

    ///////////
    if(controlRegion == Sample::wen){
      if(*met < 50) continue;
      if(*wemt > 160) continue;
    }    
    else if(controlRegion == Sample::wmn){
      if(*wmt > 160) continue;
    }

    n_metmT++;
    nwgt_metmT += *wgt;


    ///////////
    if(category == Category::monojet){

      if(jetpt->size() <= 0) continue;
      if(jetpt->at(0)  < 100) continue;
      if(fabs(jeteta->at(0)) > 2.5) continue;
      n_jetpt++;
      nwgt_jetpt += *wgt;
    
      if(chfrac->at(0) < 0.1) continue;
      if(nhfrac->at(0) > 0.8) continue;
      n_jetid++;
      nwgt_jetid += *wgt;

      jetSelection << *run << " "<<*lumi<<" "<<*event<<"\n";


      if((controlRegion == Sample::sig or controlRegion == Sample::wmn or controlRegion == Sample::zmm) and *mmet < 200) continue;
      else if((controlRegion == Sample::wen or controlRegion == Sample::zee) and *emet < 200) continue;
      else if(controlRegion == Sample::gam and *pmet < 200) continue;
      n_metcut++;
      nwgt_metcut += *wgt;
      metSelection << *run << " "<<*lumi<<" "<<*event<<"\n";

      ///////////
      if((controlRegion == Sample::sig or controlRegion == Sample::wmn or controlRegion == Sample::zmm) and fabs(*jmmdphi) < 0.5) continue;
      else if((controlRegion == Sample::wen or controlRegion == Sample::zee) and fabs(*jemdphi) < 0.5) continue;
      else if(controlRegion == Sample::gam and fabs(*jpmdphi) < 0.5) continue;
      
      n_jetdphi++;
      nwgt_jetdphi += *wgt;
    
    }
    ///////////
    else if(category == Category::monoV){

      // still apply the monojet id on the leading jet
      if(jetpt->size() <= 0) continue;
      if(jetpt->at(0)  < 100) continue;
      if(fabs(jeteta->at(0)) > 2.5) continue;
      n_jetpt++;
      nwgt_jetpt += *wgt;

      if(chfrac->at(0) < 0.1) continue;
      if(nhfrac->at(0) > 0.8) continue;
      n_jetid++;
      nwgt_jetid += *wgt;

      jetSelection << *run << " "<<*lumi<<" "<<*event<<"\n";

      if((controlRegion == Sample::sig or controlRegion == Sample::wmn or controlRegion == Sample::zmm) and *mmet < 200) continue;
      else if((controlRegion == Sample::wen or controlRegion == Sample::zee) and *emet < 200) continue;
      else if(controlRegion == Sample::gam and *pmet < 200) continue;
      n_metcut++;
      nwgt_metcut += *wgt;
      metSelection << *run << " "<<*lumi<<" "<<*event<<"\n";

      ///////////
      if((controlRegion == Sample::sig or controlRegion == Sample::wmn or controlRegion == Sample::zmm) and fabs(*jmmdphi) < 0.5) continue;
      else if((controlRegion == Sample::wen or controlRegion == Sample::zee) and fabs(*jemdphi) < 0.5) continue;
      else if(controlRegion == Sample::gam and fabs(*jpmdphi) < 0.5) continue;
      
      n_jetdphi++;
      nwgt_jetdphi += *wgt;
    
      
      if(boostedJetpt->size() <= 0) continue;
      if(boostedJetpt->at(0)  < 250) continue;
      if(fabs(boostedJeteta->at(0)) > 2.4 ) continue;
      n_ak8pt++;
      nwgt_ak8pt += *wgt;
      
      if(boostedJettau2->at(0)/boostedJettau1->at(0) > 0.6) continue;
      n_ak8tau2tau1++;
      nwgt_ak8tau2tau1 += *wgt;

      if(prunedJetm_v2->at(0) < 65 or prunedJetm_v2->at(0) > 105 ) continue;
      n_ak8mpruned++;
      nwgt_ak8mpruned += *wgt;
      
      VtaggingSelection << *run << " "<<*lumi<<" "<<*event<<"\n";

      if(*met<250) continue;
      n_metHard++;
      nwgt_metHard += *wgt;
    }
    ///////////
    else if(category == Category::VBF){
      
      if(jetpt->size() < 2) continue;

      TLorentzVector leadingJet;
      leadingJet.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
      TLorentzVector subleadingJet;
      subleadingJet.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
      
      if(leadingJet.Pt() < 80 or fabs(leadingJet.Eta()) > 4.7) continue;
      if(subleadingJet.Pt() < 40 or fabs(subleadingJet.Eta()) > 4.7) continue;
      n_jetpt++;
      nwgt_jetpt += *wgt;
      
      //if(fabs(leadingJet.Eta()) < 2.5 and chfrac->at(0) < 0.1) continue;
      //if(fabs(leadingJet.Eta()) < 2.5 and nhfrac->at(0) > 0.8) continue;
      n_jetid++;
      nwgt_jetid += *wgt;

      jetSelection << *run << " "<<*lumi<<" "<<*event<<"\n";

      if((controlRegion == Sample::sig or controlRegion == Sample::wmn or controlRegion == Sample::zmm) and *mmet < 200) continue;
      else if((controlRegion == Sample::wen or controlRegion == Sample::zee) and *emet < 200) continue;
      else if(controlRegion == Sample::gam and *pmet < 200) continue;
      n_metcut++;
      nwgt_metcut += *wgt;
      metSelection << *run << " "<<*lumi<<" "<<*event<<"\n";

      ///////////
      if((controlRegion == Sample::sig or controlRegion == Sample::wmn or controlRegion == Sample::zmm) and fabs(*jmmdphi) < 0.5) continue;
      else if((controlRegion == Sample::wen or controlRegion == Sample::zee) and fabs(*jemdphi) < 0.5) continue;
      else if(controlRegion == Sample::gam and fabs(*jpmdphi) < 0.5) continue;
      
      n_jetdphi++;
      nwgt_jetdphi += *wgt;
      
      if(leadingJet.Eta()*subleadingJet.Eta() > 0) continue;
      n_emVBF++;
      nwgt_emVBF += *wgt;

      if(fabs(leadingJet.Eta()-subleadingJet.Eta()) < 2.5) continue;
      n_detaVBF++;
      nwgt_detaVBF += *wgt;
      
      if((leadingJet+subleadingJet).M() < 500) continue; 
      n_mjjVBF++;
      nwgt_mjjVBF += *wgt;
      VBFSelection << *run << "  "<<*lumi<<" "<<*event<<"\n";
      
    }
  }
  
  leptonSelection.close();
  photonSelection.close();
  tauSelection.close();
  jetSelection.close();
  metSelection.close();
  btagSelection.close();
  VtaggingSelection.close();
  VBFSelection.close();

  cout<<"############################"<<endl;
  cout<<"## Event report in the CR ##"<<endl;
  cout<<"############################"<<endl;
  cout<<"total event = "              <<n_total<<endl;
  cout<<"muon loose selection event = "     <<n_muonLooseSelection<<" efficiency "<<float(n_muonLooseSelection)/n_total<<endl;
  cout<<"muon selection event = "     <<n_muonSelection     <<" efficiency "<<float(n_muonSelection)/n_total<<endl;
  cout<<"electron loose selection event = "     <<n_electronLooseSelection<<" efficiency "<<float(n_electronLooseSelection)/n_total<<endl;
  cout<<"electron selection event = " <<n_electronSelection <<" efficiency "<<float(n_electronSelection)/n_total<<endl;
  cout<<"photon selection event = "   <<n_photonSelection   <<" efficiency "<<float(n_photonSelection)/n_total<<endl;
  cout<<"tau selection event = "      <<n_tauSelection      <<" efficiency "<<float(n_tauSelection)/n_total<<endl;
  cout<<"bjet selection event = "     <<n_bjetSelection     <<" efficiency "<<float(n_bjetSelection)/n_total<<endl;
  cout<<"pfmet+mT selection event = " <<n_metmT             <<" efficiency "<<float(n_metmT)/n_total<<endl;
  cout<<"jetpt event = "              <<n_jetpt             <<" efficiency "<<float(n_jetpt)/n_total<<endl;
  cout<<"jetid event = "              <<n_jetid             <<" efficiency "<<float(n_jetid)/n_total<<endl;
  cout<<"met cut event = "            <<n_metcut            <<" efficiency "<<float(n_metcut)/n_total<<endl;
  cout<<"jetphi event = "             <<n_jetdphi           <<" efficiency "<<float(n_jetdphi)/n_total<<endl;
  if(category == Category::monoV){
    cout<<"ak8 pt event = "           <<n_ak8pt             <<" efficiency "<<float(n_ak8pt)/n_total<<endl;
    cout<<"ak8 tau2tau1 event = "     <<n_ak8tau2tau1       <<" efficiency "<<float(n_ak8tau2tau1)/n_total<<endl;
    cout<<"ak8 mpruned event = "      <<n_ak8mpruned        <<" efficiency "<<float(n_ak8mpruned)/n_total<<endl;
    cout<<"met hard cut = "           <<n_metHard           <<" efficiency "<<float(n_metHard)/n_total<<endl;
  }
  else if(category == Category::VBF){
    cout<<"n opposite hem VBF = "<<n_emVBF   <<" efficiency "<<float(n_emVBF)/n_total<<endl;
    cout<<"deta VBF           = "<<n_detaVBF <<" efficiency "<<float(n_detaVBF)/n_total<<endl;
    cout<<"mjj VBF            = "<<n_mjjVBF  <<" efficiency "<<float(n_mjjVBF)/n_total<<endl;
  }
  

  cout<<"#####################################"<<endl;
  cout<<"## Event report in the CR weighted ##"<<endl;
  cout<<"#####################################"<<endl;
  cout<<"total event = "               <<nwgt_total<<endl;
  cout<<"muon loose selection event = "     <<nwgt_muonLooseSelection<<" efficiency "<<float(nwgt_muonLooseSelection)/nwgt_total<<endl;
  cout<<"muon selection event = "      <<nwgt_muonSelection     <<" efficiency "<<float(nwgt_muonSelection)/nwgt_total<<endl;
  cout<<"electron loose selection event = "     <<nwgt_electronLooseSelection<<" efficiency "<<float(nwgt_electronLooseSelection)/nwgt_total<<endl;
  cout<<"electron selection event = "  <<nwgt_electronSelection <<" efficiency "<<float(nwgt_electronSelection)/nwgt_total<<endl;
  cout<<"photon selection event = "    <<nwgt_photonSelection   <<" efficiency "<<float(nwgt_photonSelection)/nwgt_total<<endl;
  cout<<"tau selection event = "       <<nwgt_tauSelection      <<" efficiency "<<float(nwgt_tauSelection)/nwgt_total<<endl;
  cout<<"bjet selection event = "      <<nwgt_bjetSelection     <<" efficiency "<<float(nwgt_bjetSelection)/nwgt_total<<endl;
  cout<<"pfmet+mT selection event = "  <<nwgt_metmT             <<" efficiency "<<float(nwgt_metmT)/nwgt_total<<endl;
  cout<<"jetpt event = "               <<nwgt_jetpt             <<" efficiency "<<float(nwgt_jetpt)/nwgt_total<<endl;
  cout<<"jetid event = "               <<nwgt_jetid             <<" efficiency "<<float(nwgt_jetid)/nwgt_total<<endl;
  cout<<"met cut event = "             <<nwgt_metcut            <<" efficiency "<<float(nwgt_metcut)/nwgt_total<<endl;
  cout<<"jetphi event = "              <<nwgt_jetdphi           <<" efficiency "<<float(nwgt_jetdphi)/nwgt_total<<endl;
  if(category == Category::monoV){
    cout<<"ak8 pt event = "              <<nwgt_ak8pt           <<" efficiency "<<float(nwgt_ak8pt)/nwgt_total<<endl;
    cout<<"ak8 tau2tau1 event = "        <<nwgt_ak8tau2tau1     <<" efficiency "<<float(nwgt_ak8tau2tau1)/nwgt_total<<endl;
    cout<<"ak8 mpruned event = "         <<nwgt_ak8mpruned      <<" efficiency "<<float(nwgt_ak8mpruned)/nwgt_total<<endl;
    cout<<"met hard cut = "              <<nwgt_metHard         <<" efficiency "<<float(nwgt_metHard)/nwgt_total<<endl;
  }
  else if(category == Category::VBF){
    cout<<"n opposite hem VBF = "  <<nwgt_emVBF   <<" efficiency "<<float(nwgt_emVBF)/nwgt_total<<endl;
    cout<<"deta VBF = "            <<nwgt_detaVBF <<" efficiency "<<float(nwgt_detaVBF)/nwgt_total<<endl;
    cout<<"mjj VBF = "             <<nwgt_mjjVBF  <<" efficiency "<<float(nwgt_mjjVBF)/nwgt_total<<endl;
  }  
}
