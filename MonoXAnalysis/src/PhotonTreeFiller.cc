#include "AnalysisCode/MonoXAnalysis/interface/PhotonTreeFiller.h"

PhotonTreeFiller::PhotonTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree):
  rhoTag            (iConfig.getParameter<edm::InputTag>("rho")),
  photonsTag        (iConfig.getParameter<edm::InputTag>("photons")),
  mediumphotonsTag  (iConfig.getParameter<edm::InputTag>("mediumphotons")),
  tightphotonsTag   (iConfig.getParameter<edm::InputTag>("tightphotons")),
  photonHighPtIdTag (iConfig.getParameter<edm::InputTag>("photonHighPtId")),
  mvaloosephotonsTag  (iConfig.getParameter<edm::InputTag>("mvaloosephotons")),
  mvatightphotonsTag  (iConfig.getParameter<edm::InputTag>("mvatightphotons")),
  applyPhotonJetsFilter (iConfig.existsAs<bool>("applyPhotonJetsFilter") ? iConfig.getParameter<bool>("applyPhotonJetsFilter") : false),
  applyDiElectronFilter (iConfig.existsAs<bool>("applyDiElectronFilter") ? iConfig.getParameter<bool>("applyDiElectronFilter") : false),
  applyDiMuonFilter     (iConfig.existsAs<bool>("applyDiMuonFilter") ? iConfig.getParameter<bool>("applyDiMuonFilter") : false),
  isQCDTree       (iConfig.existsAs<bool>("isQCDTree") ? iConfig.getParameter<bool>("isQCDTree") : false),
  isPhotonPurity  (iConfig.existsAs<bool>("isPhotonPurity") ? iConfig.getParameter<bool>("isPhotonPurity") : false),
  isReMiniAOD     (iConfig.existsAs<bool>("isReMiniAOD") ? iConfig.getParameter<bool>("isReMiniAOD") : false),
  isTriggerTree   (iConfig.existsAs<bool>("isTriggerTree") ? iConfig.getParameter<bool>("isTriggerTree") : false),
  addPhotonIDVariables (iConfig.existsAs<bool>("addPhotonIDVariables") ? iConfig.getParameter<bool>("addPhotonIDVariables") : false),
  verticesTag (iConfig.getParameter<edm::InputTag>("vertices")){   

  rhoToken = iC.consumes<double>(rhoTag);
  verticesToken = iC.consumes<std::vector<reco::Vertex> > (verticesTag);

  photonsToken        = iC.consumes<pat::PhotonRefVector> (photonsTag);
  mediumphotonsToken  = iC.consumes<pat::PhotonRefVector> (mediumphotonsTag);
  tightphotonsToken   = iC.consumes<pat::PhotonRefVector> (tightphotonsTag);
  mvaloosephotonsToken   = iC.consumes<pat::PhotonRefVector> (mvaloosephotonsTag);
  mvatightphotonsToken   = iC.consumes<pat::PhotonRefVector> (mvatightphotonsTag);
  photonHighPtIdToken = iC.consumes<edm::ValueMap<bool> > (photonHighPtIdTag);

  if(isPhotonPurity){
    photonsPurityToken      = iC.consumes<pat::PhotonRefVector> (iConfig.getParameter<edm::InputTag>("photonsPurity"));
    photonsieieToken        = iC.consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonsieie")); 
    photonPHisoToken        = iC.consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonPHiso")); 
    photonCHisoToken        = iC.consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonCHiso")); 
    photonNHisoToken        = iC.consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonNHiso")); 
    rndgammaiso04Token      = iC.consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("rndgammaiso04")); 
    rndgammaiso08Token      = iC.consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("rndgammaiso08")); 
    rndchhadiso04Token      = iC.consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("rndchhadiso04"));  
    rndchhadiso08Token      = iC.consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("rndchhadiso08"));  
  }

  if(addPhotonIDVariables){
    photonIDCollectionToken =  iC.consumes<std::vector<pat::Photon> >(iConfig.getParameter<edm::InputTag>("photonIDCollection"));
    if(not isPhotonPurity){
      photonsieieToken        = iC.consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonsieie"));
      photonPHisoToken        = iC.consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonPHiso"));
      photonCHisoToken        = iC.consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonCHiso"));
      photonNHisoToken        = iC.consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonNHiso"));
    }
  }

  tree_ = tree;
  this->DeclareAndSetBranches();
  this->initBranches();
}

/////
void PhotonTreeFiller::initBranches(){

  phidm    = 0; phidt    = 0; phidh    = 0; phgs = 0;
  phpt     = 0; pheta    = 0; phphi    = 0;
  nphotons = 0; nmvaloosephotons = 0; nmvatightphotons = 0;
  phidmval = 0; phidmvat = 0;
 
  phPuritypt     = 0.0; phPurityeta    = 0.0; phPurityphi    = 0.0;
  phPurityPHiso  = 0.0; phPurityNHiso  = 0.0; phPurityCHiso  = 0.0;
  phPurityRND04PHiso  = 0.0; phPurityRND08PHiso  = 0.0;
  phPurityRND08CHiso  = 0.0; phPurityRND04CHiso  = 0.0;
  phNHiso          = 0.0; phCHiso          = 0.0; phPHiso          = 0.0;
  phPuritysieie    = 0.0; phPurityhoe      = 0.0; phPurityElectronVeto = 0;
  phPurityEAEGamma = 0.;

  photonPt.clear();
  photonEta.clear();
  photonPhi.clear();
  photonE.clear(); 
  photonSCEta.clear(); 
  photonSCPhi.clear();
  photonSCEnergy.clear();
  photonSCRawEnergy.clear();
  photonHOverE.clear();
  photonSigmaIetaIeta.clear();
  photonChargedIso.clear();      
  photonNeutralIso.clear();
  photonEMIso.clear(); 
  photonElectronVeto.clear();

}

/////
bool PhotonTreeFiller::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

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


  Handle<pat::PhotonRefVector> photonsH;
  iEvent.getByToken(photonsToken, photonsH);
  pat::PhotonRefVector photons = *photonsH;

  Handle<pat::PhotonRefVector> mediumPhotonsH;
  iEvent.getByToken(mediumphotonsToken, mediumPhotonsH);
  pat::PhotonRefVector mediumphotons = *mediumPhotonsH;

  Handle<pat::PhotonRefVector> tightPhotonsH;
  iEvent.getByToken(tightphotonsToken, tightPhotonsH);
  pat::PhotonRefVector tightphotons = *tightPhotonsH;

  Handle<pat::PhotonRefVector> mvaloosePhotonsH;
  iEvent.getByToken(mvaloosephotonsToken,mvaloosePhotonsH);
  pat::PhotonRefVector mvaloosePhotons = *mvaloosePhotonsH;

  Handle<pat::PhotonRefVector> mvatightPhotonsH;
  iEvent.getByToken(mvatightphotonsToken,mvatightPhotonsH);
  pat::PhotonRefVector mvatightPhotons = *mvatightPhotonsH;


  Handle<ValueMap<bool> > photonHighPtIdH;
  iEvent.getByToken(photonHighPtIdToken, photonHighPtIdH);

  Handle<pat::PhotonRefVector> photonsPurityH;
  pat::PhotonRefVector photonsPurity;
  Handle<edm::ValueMap<float> > photonsieieH;
  Handle<edm::ValueMap<float> > photonPHisoH;
  Handle<edm::ValueMap<float> > photonCHisoH;
  Handle<edm::ValueMap<float> > photonNHisoH;
  Handle<edm::ValueMap<float> > rndgammaiso04H;
  Handle<edm::ValueMap<float> > rndgammaiso08H;
  Handle<edm::ValueMap<float> > rndchhadiso04H;
  Handle<edm::ValueMap<float> > rndchhadiso08H;

  if(isPhotonPurity){
    iEvent.getByToken(photonsPurityToken, photonsPurityH);    
    photonsPurity = *photonsPurityH;
    iEvent.getByToken(photonsieieToken, photonsieieH);
    iEvent.getByToken(photonPHisoToken, photonPHisoH);      
    iEvent.getByToken(photonCHisoToken, photonCHisoH);
    iEvent.getByToken(photonNHisoToken, photonNHisoH);
    iEvent.getByToken(rndgammaiso04Token, rndgammaiso04H);
    iEvent.getByToken(rndgammaiso08Token, rndgammaiso08H);      
    iEvent.getByToken(rndchhadiso04Token, rndchhadiso04H);
    iEvent.getByToken(rndchhadiso08Token, rndchhadiso08H);
  }

  // Photon and electron ID
  Handle<vector<pat::Photon> > photonIDH;
  if(addPhotonIDVariables){
    iEvent.getByToken(photonIDCollectionToken,photonIDH);
    if(not isPhotonPurity){
      iEvent.getByToken(photonsieieToken, photonsieieH);
      iEvent.getByToken(photonPHisoToken, photonPHisoH);      
      iEvent.getByToken(photonCHisoToken, photonCHisoH);
      iEvent.getByToken(photonNHisoToken, photonNHisoH);
    }
  }

  //////
  int hardestPhotonIndex = -1;
  float hardestPhotonPt = 0.0;

  if(photonsH.isValid()){    
    nphotons = photons.size();      
    for (size_t i = 0; i < photons.size(); i++) {
      if (photons[i]->pt() > hardestPhotonPt) {
	hardestPhotonIndex = i;
	hardestPhotonPt = photons[i]->pt();
      }
    }

    if(hardestPhotonIndex >= 0){
      phpt    = photons[hardestPhotonIndex]->pt();
      pheta   = photons[hardestPhotonIndex]->eta();
      phphi   = photons[hardestPhotonIndex]->phi();
      if(isReMiniAOD)
	phgs = photons[hardestPhotonIndex]->userInt("hasGainSwitchFlag");
    }
    if(mediumPhotonsH.isValid()){
      if (hardestPhotonIndex >= 0) {
	for(size_t i = 0; i < mediumphotons.size(); i++){
	  if(photons[hardestPhotonIndex] == mediumphotons[i])
	    phidm = 1;
	}
      }
    }
    if(tightPhotonsH.isValid()){
      if (hardestPhotonIndex >= 0) {
	for(size_t i = 0; i < tightphotons.size(); i++){
	  if(tightphotons[hardestPhotonIndex] == tightphotons[i])
	    phidt = 1;
	}
      }
    }

    if(photonHighPtIdH.isValid()){
      if (hardestPhotonIndex >= 0) 
	phidh   = ((*photonHighPtIdH)[photons[hardestPhotonIndex]] ? 1 : 0);
    }
      
      
    if(mvaloosePhotonsH.isValid()){
      nmvaloosephotons = mvaloosePhotons.size();
      if (hardestPhotonIndex >= 0) {
	for(size_t i = 0; i < mvaloosePhotons.size(); i++){
	  if(mvaloosePhotons[hardestPhotonIndex] == mvaloosePhotons[i])
	    phidmval = 1;
	}
      }
    }

    if(mvatightPhotonsH.isValid()){
      nmvatightphotons = mvatightPhotons.size();
      if (hardestPhotonIndex >= 0) {
	for(size_t i = 0; i < mvatightPhotons.size(); i++){
	  if(mvatightPhotons[hardestPhotonIndex] == mvatightPhotons[i])
	    phidmvat = 1;
	}
      }
    }
  }
    
  int hardestPhotonPurityIndex = -1;
  float hardestPhotonPurityPt = 0.0;


  ///////
  if(isPhotonPurity and not isTriggerTree){

    nphotonsPurity = photonsPurityH->size();      
    for (size_t i = 0; i < photonsPurity.size(); i++) {
      if (photonsPurity[i]->pt() > hardestPhotonPurityPt) {
	hardestPhotonPurityPt = photonsPurity[i]->pt();
	hardestPhotonPurityIndex = i;  
      }
    }
        
    if (hardestPhotonPurityIndex >= 0) {

      phPuritypt    = photonsPurity[hardestPhotonPurityIndex]->pt();
      phPurityeta   = photonsPurity[hardestPhotonPurityIndex]->eta();
      phPurityphi   = photonsPurity[hardestPhotonPurityIndex]->phi();

      phPHiso       = (*photonPHisoH)[photonsPurity[hardestPhotonPurityIndex]];
      phCHiso       = (*photonCHisoH)[photonsPurity[hardestPhotonPurityIndex]];
      phNHiso       = (*photonNHisoH)[photonsPurity[hardestPhotonPurityIndex]];

      phPurityPHiso  = max(0.,double((*photonPHisoH)[photonsPurity[hardestPhotonPurityIndex]]-rho*getGammaEAForPhotonIso(photonsPurity[hardestPhotonPurityIndex]->eta())));
      phPurityCHiso  = max(0.,double((*photonCHisoH)[photonsPurity[hardestPhotonPurityIndex]]-rho*getChargedHadronEAForPhotonIso(photonsPurity[hardestPhotonPurityIndex]->eta())));
      phPurityNHiso  = max(0.,double((*photonNHisoH)[photonsPurity[hardestPhotonPurityIndex]]-rho*getNeutralHadronEAForPhotonIso(photonsPurity[hardestPhotonPurityIndex]->eta()))); 
      
      phPurityRND04CHiso = (*rndchhadiso04H)[photonsPurity[hardestPhotonPurityIndex]];
      phPurityRND04PHiso = (*rndgammaiso04H)[photonsPurity[hardestPhotonPurityIndex]];
      phPurityRND08PHiso = (*rndgammaiso08H)[photonsPurity[hardestPhotonPurityIndex]];
      phPurityRND08CHiso = (*rndchhadiso08H)[photonsPurity[hardestPhotonPurityIndex]];

      phPuritysieie        = (*photonsieieH)[photonsPurity[hardestPhotonPurityIndex]];
      phPurityElectronVeto = photonsPurity[hardestPhotonPurityIndex]->passElectronVeto();
      phPurityhoe          = photonsPurity[hardestPhotonPurityIndex]->hadTowOverEm();
      phPurityEAEGamma     = getGammaEAForPhotonIso(photonsPurity[hardestPhotonPurityIndex]->eta());
    }
  }

  /////
  if(photonIDH.isValid() and addPhotonIDVariables and not isTriggerTree and not applyDiMuonFilter and not applyDiElectronFilter){
   
    for(auto photon_iter = photonIDH->begin(); photon_iter != photonIDH->end(); ++photon_iter){
      if(photon_iter->pt() < 35 or fabs(photon_iter->superCluster()->eta()) > 2.5) continue;
      photonPt.push_back(photon_iter->pt());
      photonEta.push_back(photon_iter->eta());
      photonPhi.push_back(photon_iter->phi());
      photonE.push_back(photon_iter->energy());
      photonElectronVeto.push_back(photon_iter->passElectronVeto());
      photonSCEta.push_back(photon_iter->superCluster()->eta());
      photonSCPhi.push_back(photon_iter->superCluster()->phi());
      photonSCEnergy.push_back(photon_iter->superCluster()->energy());
      photonSCRawEnergy.push_back(photon_iter->superCluster()->rawEnergy());
      photonHOverE.push_back(photon_iter->hadTowOverEm());

      pat::PhotonRef phoref(photonIDH, photon_iter - photonIDH->begin());
      if(phoref.isAvailable() and phoref.isNonnull()){
	if(photonsieieH.isValid())
	  photonSigmaIetaIeta.push_back((*photonsieieH)[phoref]);
	if(photonCHisoH.isValid()){
	  photonChargedIso.push_back((*photonCHisoH)[phoref]);
	}
	if(photonPHisoH.isValid()){
	  photonEMIso.push_back((*photonPHisoH)[phoref]);
	}
	if(photonNHisoH.isValid()){
	  photonNeutralIso.push_back((*photonNHisoH)[phoref]);
	}
      }
    }      
  }

  return true;  
}

void PhotonTreeFiller::DeclareAndSetBranches(){

  tree_->Branch("nphotons"             , &nphotons             , "nphotons/i");
  if(not isTriggerTree)
    tree_->Branch("nmvaloosephotons"     , &nmvaloosephotons     , "nmvaloosephotons/i");
  tree_->Branch("nmvatightphotons"     , &nmvatightphotons     , "nmvatightphotons/i");

  // Photon info
  tree_->Branch("phidm"                , &phidm                , "phidm/I");
  tree_->Branch("phidt"                , &phidt                , "phidt/I");
  tree_->Branch("phidh"                , &phidh                , "phidh/I");
  tree_->Branch("phidmval"             , &phidmval             , "phidmval/I");
  tree_->Branch("phidmvat"             , &phidmvat             , "phidmvat/I");
  tree_->Branch("phpt"                 , &phpt                 , "phpt/F");
  tree_->Branch("pheta"                , &pheta                , "pheta/F");
  tree_->Branch("phphi"                , &phphi                , "phphi/F");
  if(isReMiniAOD){
    tree_->Branch("phgs"               , &phgs                , "phgs/I");
  }
  


  if(isPhotonPurity and not isTriggerTree and not isQCDTree){

    tree_->Branch("nphotonsPurity"  , &nphotonsPurity  , "nphotonsPurity/i");
    tree_->Branch("rho"             , &rho             , "rho/F");
    tree_->Branch("phPHiso"         , &phPHiso         , "phPHiso/F");
    tree_->Branch("phCHiso"         , &phCHiso         , "phCHiso/F");
    tree_->Branch("phNHiso"         , &phNHiso         , "phNHiso/F");

    tree_->Branch("phPuritypt"      , &phPuritypt      , "phPuritypt/F");
    tree_->Branch("phPurityeta"     , &phPurityeta     , "phPurityeta/F");
    tree_->Branch("phPurityphi"     , &phPurityphi     , "phPurityphi/F");

    tree_->Branch("phPurityPHiso"       , &phPurityPHiso   , "phPurityPHiso/F");
    tree_->Branch("phPurityRND04PHiso"  , &phPurityRND04PHiso    , "phPurityRND04PHiso/F");
    tree_->Branch("phPurityRND08PHiso"  , &phPurityRND08PHiso    , "phPurityRND08PHiso/F");
    tree_->Branch("phPurityCHiso"       , &phPurityCHiso         , "phPurityCHiso/F");
    tree_->Branch("phPurityRND04CHiso"  , &phPurityRND04CHiso    , "phPurityRND04CHiso/F");
    tree_->Branch("phPurityRND08CHiso"  , &phPurityRND08CHiso    , "phPurityRND08CHiso/F");
    tree_->Branch("phPurityNHiso"       , &phPurityNHiso         , "phPurityNHiso/F");

    tree_->Branch("phPuritysieie"        , &phPuritysieie         , "phPuritysieie/F");
    tree_->Branch("phPurityhoe"          , &phPurityhoe           , "phPurityhoe/F");
    tree_->Branch("phPurityEAEGamma"     , &phPurityEAEGamma      , "phPurityEAEGamma/F");
    tree_->Branch("phPurityElectronVeto" , &phPurityElectronVeto  , "phPurityElectronVeto/F");

  }
  
  if(addPhotonIDVariables and not isTriggerTree and not isQCDTree and not applyDiMuonFilter and not applyDiElectronFilter){
    if(not isPhotonPurity)
      tree_->Branch("rho"             , &rho             , "rho/F");
    tree_->Branch("photonPt", "std::vector<float>", &photonPt);
    tree_->Branch("photonEta", "std::vector<float>", &photonEta);
    tree_->Branch("photonPhi", "std::vector<float>", &photonPhi);
    tree_->Branch("photonE", "std::vector<float>", &photonE);
    tree_->Branch("photonSCEta", "std::vector<float>", &photonSCEta);
    tree_->Branch("photonSCPhi", "std::vector<float>", &photonSCPhi);
    tree_->Branch("photonSCEnergy", "std::vector<float>", &photonSCEnergy);
    tree_->Branch("photonSCRawEnergy", "std::vector<float>", &photonSCRawEnergy);
    tree_->Branch("photonHOverE", "std::vector<float>", &photonHOverE);
    tree_->Branch("photonSigmaIetaIeta", "std::vector<float>", &photonSigmaIetaIeta);
    tree_->Branch("photonChargedIso", "std::vector<float>", &photonChargedIso);
    tree_->Branch("photonNeutralIso", "std::vector<float>", &photonNeutralIso);
    tree_->Branch("photonEMIso", "std::vector<float>", &photonEMIso);
    tree_->Branch("photonElectronVeto", "std::vector<float>", &photonElectronVeto);
  }  
}

float PhotonTreeFiller::getChargedHadronEAForPhotonIso(float eta) { // 80X values 
  if (fabs(eta) < 1.0) return 0.0360;
  else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return  0.0377;
  else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return  0.0306;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return  0.0283;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return  0.0254;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.4  ) return  0.0217;
  else if (fabs(eta) >= 2.4) return 0.0167 ;
  else return 0.;
}

float PhotonTreeFiller::getNeutralHadronEAForPhotonIso(float eta) { // 80X values  
  if (fabs(eta) < 1.0) return 0.0597;
  else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.0807;
  else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.0629;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.0197;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.0184;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.4  ) return 0.0284;
  else if (fabs(eta) >= 2.4) return 0.0591;
  else return 0.;
}

float PhotonTreeFiller::getGammaEAForPhotonIso(float eta) {// 80X values  
  if (fabs(eta) < 1.0) return 0.1210;
  else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.1107;
  else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.0699;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.1056;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.1457;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.4  ) return 0.1719;
  else if (fabs(eta) >= 2.4) return 0.1998;
  else return 0.;
}

