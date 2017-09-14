#include "AnalysisCode/MonoXAnalysis/interface/ElectronTreeFiller.h"

ElectronTreeFiller::ElectronTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree):
  rhoTag            (iConfig.getParameter<edm::InputTag>("rho")),
  electronsTag      (iConfig.getParameter<edm::InputTag>("electrons")),
  looseelectronsTag (iConfig.getParameter<edm::InputTag>("looseelectrons")),
  tightelectronsTag (iConfig.getParameter<edm::InputTag>("tightelectrons")),
  triggerelectronsTag (iConfig.getParameter<edm::InputTag>("triggerelectrons")),
  heepelectronsTag    (iConfig.getParameter<edm::InputTag>("heepelectrons")),
  mvalooseelectronsTag  (iConfig.getParameter<edm::InputTag>("mvalooseelectrons")),
  mvatightelectronsTag  (iConfig.getParameter<edm::InputTag>("mvatightelectrons")),
  applyDiElectronFilter (iConfig.existsAs<bool>("applyDiElectronFilter") ? iConfig.getParameter<bool>("applyDiElectronFilter") : false),
  applyDiMuonFilter     (iConfig.existsAs<bool>("applyDiMuonFilter") ? iConfig.getParameter<bool>("applyDiMuonFilter") : false),
  applyPhotonJetsFilter (iConfig.existsAs<bool>("applyPhotonJetsFilter") ? iConfig.getParameter<bool>("applyPhotonJetsFilter") : false),
  isQCDTree       (iConfig.existsAs<bool>("isQCDTree") ? iConfig.getParameter<bool>("isQCDTree") : false),
  isPhotonPurity  (iConfig.existsAs<bool>("isPhotonPurity") ? iConfig.getParameter<bool>("isPhotonPurity") : false),
  isReMiniAOD     (iConfig.existsAs<bool>("isReMiniAOD") ? iConfig.getParameter<bool>("isReMiniAOD") : false),
  isTriggerTree   (iConfig.existsAs<bool>("isTriggerTree") ? iConfig.getParameter<bool>("isTriggerTree") : false),
  addElectronIDVariables (iConfig.existsAs<bool>("addElectronIDVariables") ? iConfig.getParameter<bool>("addElectronIDVariables") : false),
  t1metTag    (iConfig.getParameter<edm::InputTag>("t1met")),
  verticesTag (iConfig.getParameter<edm::InputTag>("vertices")){   

  rhoToken = iC.consumes<double>(rhoTag);
  verticesToken = iC.consumes<std::vector<reco::Vertex> > (verticesTag);
  t1metToken    = iC.consumes<edm::View<pat::MET> > (t1metTag);

  electronsToken       = iC.consumes<pat::ElectronRefVector> (electronsTag);
  looseelectronsToken  = iC.consumes<pat::ElectronRefVector> (looseelectronsTag);
  tightelectronsToken  = iC.consumes<pat::ElectronRefVector> (tightelectronsTag);
  triggerelectronsToken  = iC.consumes<pat::ElectronRefVector> (triggerelectronsTag);
  heepelectronsToken   = iC.consumes<pat::ElectronRefVector> (heepelectronsTag);
  mvalooseelectronsToken   = iC.consumes<pat::ElectronRefVector> (mvalooseelectronsTag);
  mvatightelectronsToken   = iC.consumes<pat::ElectronRefVector> (mvatightelectronsTag);
  
  if(addElectronIDVariables)
    electronIDCollectionToken = iC.consumes<std::vector<pat::Electron>> (iConfig.getParameter<edm::InputTag>("electronIDCollection"));    

  tree_ = tree;

  this->DeclareAndSetBranches();
  this->initBranches();
}

/////
void ElectronTreeFiller::initBranches(){
  
  nelectrons  = 0; nlooseelectrons = 0; ntightelectrons = 0; nheepelectrons = 0; ntriggerelectrons = 0; nmvalooseelectrons = 0; nmvatightelectrons = 0;   
  zeemass     = 0.0; zeept       = 0.0; zeeeta      = 0.0; zeephi      = 0.0;  wemt = 0.0;
  el1pid      = 0; el1pt       = 0.0; el1eta      = 0.0; el1phi      = 0.0; el1id       = 0; el1idl       = 0; el1idt       = 0; el1idmval = 0; el1idmvat = 0;
  el2pid      = 0; el2pt       = 0.0; el2eta      = 0.0; el2phi      = 0.0; el2id       = 0; el2idl       = 0; el2idt       = 0; el2idmval = 0; el2idmvat = 0;
  el1gs       = 0; el2gs       = 0;
  
  electronPt.clear();
  electronEta.clear();
  electronPhi.clear();
  electronE.clear();
  electronSCEta.clear();
  electronSCPhi.clear();
  electronSCEnergy.clear();
  electronSCRawEnergy.clear();
  electronHOverE.clear();
  electronSigmaIetaIeta.clear();    
  electronChargedIso.clear(); 
  electronNeutralIso.clear(); 
  electronEMIso.clear(); 
  electronGsfPt.clear(); 
  electronEOP.clear(); 
  electronDxy.clear();
  electronDz.clear();
  electronDphi.clear(); 
  electronDeta.clear(); 
  electronMissHit.clear(); 
  electronConversion.clear();
  
}

/////
bool ElectronTreeFiller::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace pat;

  this->initBranches();

  Handle<double> rhoH;
  iEvent.getByToken(rhoToken, rhoH);
  rho = *rhoH;

  // VERTEX
  Handle<vector<Vertex> > verticesH;
  iEvent.getByToken(verticesToken, verticesH);

  //MET in the event
  Handle<View<pat::MET> > t1metH;
  iEvent.getByToken(t1metToken, t1metH);
  

  // ELECTRONS
  Handle<pat::ElectronRefVector> electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  pat::ElectronRefVector electrons = *electronsH;

  Handle<pat::ElectronRefVector> looseelectronsH;
  iEvent.getByToken(looseelectronsToken, looseelectronsH);
  pat::ElectronRefVector looseelectrons = *looseelectronsH;

  Handle<pat::ElectronRefVector> tightelectronsH;
  iEvent.getByToken(tightelectronsToken, tightelectronsH);
  pat::ElectronRefVector tightelectrons = *tightelectronsH;

  Handle<pat::ElectronRefVector> triggerelectronsH;
  iEvent.getByToken(triggerelectronsToken, triggerelectronsH);
  pat::ElectronRefVector triggerelectrons = *triggerelectronsH;

  Handle<pat::ElectronRefVector> heepelectronsH;
  iEvent.getByToken(heepelectronsToken, heepelectronsH);
  pat::ElectronRefVector heepelectrons = *heepelectronsH;

  Handle<pat::ElectronRefVector> mvalooseelectronsH;
  iEvent.getByToken(mvalooseelectronsToken, mvalooseelectronsH);
  pat::ElectronRefVector mvalooseelectrons = *mvalooseelectronsH;

  Handle<pat::ElectronRefVector> mvatightelectronsH;
  iEvent.getByToken(mvatightelectronsToken, mvatightelectronsH);
  pat::ElectronRefVector mvatightelectrons = *mvatightelectronsH;

  Handle<vector<pat::Electron> > electronIDH;
  if(addElectronIDVariables)
    iEvent.getByToken(electronIDCollectionToken,electronIDH);

  // electron counters
  vector<pat::ElectronRef> electronvector;
  if(electronsH.isValid()){
    nelectrons      = electronsH->size();
    for (size_t i = 0; i < electrons.size(); i++) 
      electronvector.push_back(electrons[i]);
  }

  // couting objects
  if(looseelectronsH.isValid())
    nlooseelectrons = looseelectronsH->size();
  if(tightelectronsH.isValid())
    ntightelectrons = tightelectronsH->size();
  if(heepelectronsH.isValid())
    nheepelectrons  = heepelectronsH->size();
  if(triggerelectronsH.isValid())
    ntriggerelectrons = triggerelectronsH->size();      
  if(mvalooseelectronsH.isValid())
    nmvalooseelectrons = mvalooseelectronsH->size();      
  if(mvatightelectronsH.isValid())
    nmvatightelectrons = mvatightelectronsH->size();      

  sort(electronvector.begin(), electronvector.end(), electronPtSorter);
  
  // one or two loose electrons
  if (nelectrons == 1 || nelectrons == 2) {
    pat::ElectronRef electron = electronvector[0];
    el1pid = electron->pdgId();
    el1pt  = electron->pt();
    el1eta = electron->eta();
    el1phi = electron->phi();
    if(isReMiniAOD)
      el1gs  = electron->userInt("hasGainSwitchFlag");
    for(std::size_t i = 0; i < looseelectrons.size(); i++) {
      if(electron == looseelectrons[i])
	el1idl = 1;
    }
      
    for (std::size_t i = 0; i < tightelectrons.size(); i++) {
      if (electron == tightelectrons[i]) 
	el1id = 1;
    }

    for (std::size_t i = 0; i < triggerelectrons.size(); i++) {
      if (electron == triggerelectrons[i]) 
	el1idt = 1;
    }
      
    for (std::size_t i = 0; i < heepelectrons.size(); i++) {
      if (electron == heepelectrons[i] and el1id != 1) 
	el1id = 2;
    }

    for (std::size_t i = 0; i < mvalooseelectrons.size(); i++) {
      if (electron == mvalooseelectrons[i]) 
	el1idmval = 1;
    }

    for (std::size_t i = 0; i < mvatightelectrons.size(); i++) {
      if (electron == mvatightelectrons[i]) 
	el1idmvat = 1;
    }
              
    if (electrons.size() == 1) 
      wemt = sqrt(2.0 * el1pt * t1metH->front().corPt() * (1.0 - cos(deltaPhi(el1phi, t1metH->front().corPhi()))));
  }
  
  // two loose electrons
  if (nelectrons == 2) {
    pat::ElectronRef electron = electronvector[1];
    el2pid = electron->pdgId();
    el2pt  = electron->pt();
    el2eta = electron->eta();
    el2phi = electron->phi();
    if(isReMiniAOD)
      el2gs  = electron->userInt("hasGainSwitchFlag");

    for (std::size_t i = 0; i < looseelectrons.size(); i++) {
      if (electron == looseelectrons[i]) el2idl = 1;
    }

    for (std::size_t i = 0; i < tightelectrons.size(); i++) {
      if (electron == tightelectrons[i]) el2id = 1;
    }

    for (std::size_t i = 0; i < triggerelectrons.size(); i++) {
      if (electron == triggerelectrons[i]) el2idt = 1;
    }

    for (std::size_t i = 0; i < heepelectrons.size(); i++) {
      if (electron == heepelectrons[i] and el2id != 1) 
	el2id = 2;
    }
    
    for (std::size_t i = 0; i < mvalooseelectrons.size(); i++) {
      if (electron == mvalooseelectrons[i]) 
	el2idmval = 1;
    }

    for (std::size_t i = 0; i < mvatightelectrons.size(); i++) {
      if (electron == mvatightelectrons[i]) 
	el2idmvat = 1;
    }
    
    TLorentzVector el1vec; el1vec.SetPtEtaPhiE(el1pt, el1eta, el1phi, electronvector[0]->p());
    TLorentzVector el2vec; el2vec.SetPtEtaPhiE(el2pt, el2eta, el2phi, electron->p());

    TLorentzVector zvec(el1vec);
    zvec += el2vec;

    zeemass = zvec.M();
    zeept   = zvec.Pt();
    zeeeta  = zvec.Eta();
    zeephi  = zvec.Phi();
  }

  if(electronIDH.isValid() and addElectronIDVariables and not isTriggerTree and not applyDiMuonFilter and not applyPhotonJetsFilter){

    for(auto electron_iter = electronIDH->begin(); electron_iter != electronIDH->end(); ++electron_iter){
      if(electron_iter->pt() < 35 or fabs(electron_iter->superCluster()->eta()) > 2.5) continue;
      electronPt.push_back(electron_iter->pt());
      electronEta.push_back(electron_iter->eta());
      electronPhi.push_back(electron_iter->phi());
      electronE.push_back(electron_iter->energy());
      electronSCEta.push_back(electron_iter->superCluster()->eta());
      electronSCPhi.push_back(electron_iter->superCluster()->phi());
      electronSCEnergy.push_back(electron_iter->superCluster()->energy());
      electronSCRawEnergy.push_back(electron_iter->superCluster()->rawEnergy());
      electronHOverE.push_back(electron_iter->hadronicOverEm());
      electronSigmaIetaIeta.push_back(electron_iter->full5x5_sigmaIetaIeta());
      electronChargedIso.push_back(electron_iter->pfIsolationVariables().sumChargedHadronPt);
      electronNeutralIso.push_back(electron_iter->pfIsolationVariables().sumNeutralHadronEt);
      electronEMIso.push_back(electron_iter->pfIsolationVariables().sumPhotonEt);
      electronGsfPt.push_back(electron_iter->gsfTrack()->pt());
      electronDxy.push_back(electron_iter->gsfTrack()->dxy(verticesH->begin()->position()));
      electronDz.push_back(electron_iter->gsfTrack()->dz(verticesH->begin()->position()));
      electronEOP.push_back(fabs(1-electron_iter->eSuperClusterOverP())*1./electron_iter->ecalEnergy()); 
      electronDphi.push_back(electron_iter->deltaPhiSuperClusterTrackAtVtx());
      electronDeta.push_back(electron_iter->deltaEtaSuperClusterTrackAtVtx());
      electronMissHit.push_back(electron_iter->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
      electronConversion.push_back(electron_iter->passConversionVeto());
    }      
  }
  
  return true;  
}

void ElectronTreeFiller::DeclareAndSetBranches(){
  
  tree_->Branch("nelectrons"           , &nelectrons           , "nelectrons/i");
  tree_->Branch("nlooseelectrons"      , &nlooseelectrons      , "nlooseelectrons/i");
  tree_->Branch("ntightelectrons"      , &ntightelectrons      , "ntightelectrons/i");
  tree_->Branch("nmvatightelectrons"   , &nmvatightelectrons   , "nmvatightelectrons/i");

  if(not isTriggerTree){
    tree_->Branch("ntriggerelectrons"    , &ntriggerelectrons    , "ntriggerelectrons/i");
    tree_->Branch("nheepelectrons"       , &nheepelectrons       , "nheepelectrons/i");
    tree_->Branch("nmvalooseelectrons"   , &nmvalooseelectrons   , "nmvalooseelectrons/i");
  }

  tree_->Branch("el1pid"               , &el1pid               , "el1pid/I");
  tree_->Branch("el1pt"                , &el1pt                , "el1pt/F");
  tree_->Branch("el1eta"               , &el1eta               , "el1eta/F");
  tree_->Branch("el1phi"               , &el1phi               , "el1phi/F");
  tree_->Branch("el1id"                , &el1id                , "el1id/I");
  tree_->Branch("el1idl"               , &el1idl               , "el1idl/I");
  tree_->Branch("el1idt"               , &el1idt               , "el1idt/I");
  tree_->Branch("el1idmvat"            , &el1idmvat            , "el1idmvat/I");
  tree_->Branch("el1idmval"            , &el1idmval            , "el1idmval/I");
  if(isReMiniAOD)
    tree_->Branch("el1gs"             , &el1gs                 , "el1gs/I");


  tree_->Branch("el2pid"               , &el2pid               , "el2pid/I");
  tree_->Branch("el2pt"                , &el2pt                , "el2pt/F");
  tree_->Branch("el2eta"               , &el2eta               , "el2eta/F");
  tree_->Branch("el2phi"               , &el2phi               , "el2phi/F");
  tree_->Branch("el2id"                , &el2id                , "el2id/I");
  tree_->Branch("el2idl"               , &el2idl               , "el2idl/I");
  tree_->Branch("el2idt"               , &el2idt               , "el2idt/I");
  tree_->Branch("el2idmvat"            , &el2idmvat            , "el2idmvat/I");
  tree_->Branch("el2idmval"            , &el2idmval            , "el2idmval/I");
  if(isReMiniAOD)
    tree_->Branch("el2gs"             , &el2gs                 , "el2gs/I");

  tree_->Branch("zeemass"              , &zeemass              , "zeemass/F");
  if(not isTriggerTree and not isQCDTree and not isPhotonPurity){
    tree_->Branch("zeept"                , &zeept                , "zeept/F");
    tree_->Branch("zeeeta"               , &zeeeta               , "zeeeta/F");
    tree_->Branch("zeephi"               , &zeephi               , "zeephi/F");
  }
  
  if(addElectronIDVariables and not isTriggerTree and not isPhotonPurity and not isQCDTree and not applyDiMuonFilter and not applyPhotonJetsFilter){
    if(not isPhotonPurity)
      tree_->Branch("rho"             , &rho             , "rho/F");
    tree_->Branch("electronPt", "std::vector<float>", &electronPt);
    tree_->Branch("electronEta", "std::vector<float>", &electronEta);
    tree_->Branch("electronPhi", "std::vector<float>", &electronPhi);
    tree_->Branch("electronE", "std::vector<float>", &electronE);
    tree_->Branch("electronSCEta", "std::vector<float>", &electronSCEta);
    tree_->Branch("electronSCPhi", "std::vector<float>", &electronSCPhi);
    tree_->Branch("electronSCEnergy", "std::vector<float>", &electronSCEnergy);
    tree_->Branch("electronSCRawEnergy", "std::vector<float>", &electronSCRawEnergy);
    tree_->Branch("electronHOverE", "std::vector<float>", &electronHOverE);
    tree_->Branch("electronSigmaIetaIeta", "std::vector<float>", &electronSigmaIetaIeta);
    tree_->Branch("electronChargedIso", "std::vector<float>", &electronChargedIso);
    tree_->Branch("electronNeutralIso", "std::vector<float>", &electronNeutralIso);
    tree_->Branch("electronEMIso", "std::vector<float>", &electronEMIso);
    tree_->Branch("electronGsfPt", "std::vector<float>", &electronGsfPt);
    tree_->Branch("electronDphi", "std::vector<float>", &electronDphi);
    tree_->Branch("electronDeta", "std::vector<float>", &electronDeta);
    tree_->Branch("electronEOP", "std::vector<float>", &electronEOP);
    tree_->Branch("electronMissHit", "std::vector<float>", &electronMissHit);
    tree_->Branch("electronConversion", "std::vector<float>", &electronConversion);
    tree_->Branch("electronDz", "std::vector<float>", &electronDz);
    tree_->Branch("electronDxy", "std::vector<float>", &electronDxy);
  }
}

