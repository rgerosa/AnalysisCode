#include "AnalysisCode/MonoXAnalysis/interface/TauTreeFiller.h"

TauTreeFiller::TauTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree):
  muonsTag       (iConfig.getParameter<edm::InputTag>("muons")),
  electronsTag   (iConfig.getParameter<edm::InputTag>("electrons")),
  tausCollection(iConfig.getParameter<edm::InputTag>("tausCollection")),
  tausVLNewTag(iConfig.getParameter<edm::InputTag>("tausVLNew")),
  tausVLOldTag(iConfig.getParameter<edm::InputTag>("tausVLOld")),
  tausRawNewTag(iConfig.getParameter<edm::InputTag>("tausRawNew")),
  tausRawOldTag(iConfig.getParameter<edm::InputTag>("tausRawOld")),
  tausTightNewTag(iConfig.getParameter<edm::InputTag>("tausTightNew")),
  tausTightOldTag(iConfig.getParameter<edm::InputTag>("tausTightOld")),
  applyPhotonJetsFilter (iConfig.existsAs<bool>("applyPhotonJetsFilter") ? iConfig.getParameter<bool>("applyPhotonJetsFilter") : false),
  applyDiElectronFilter (iConfig.existsAs<bool>("applyDiElectronFilter") ? iConfig.getParameter<bool>("applyDiElectronFilter") : false),
  applyDiMuonFilter     (iConfig.existsAs<bool>("applyDiMuonFilter") ? iConfig.getParameter<bool>("applyDiMuonFilter") : false),
  isQCDTree       (iConfig.existsAs<bool>("isQCDTree") ? iConfig.getParameter<bool>("isQCDTree") : false),
  isPhotonPurity  (iConfig.existsAs<bool>("isPhotonPurity") ? iConfig.getParameter<bool>("isPhotonPurity") : false),
  isReMiniAOD     (iConfig.existsAs<bool>("isReMiniAOD") ? iConfig.getParameter<bool>("isReMiniAOD") : false),
  isTriggerTree   (iConfig.existsAs<bool>("isTriggerTree") ? iConfig.getParameter<bool>("isTriggerTree") : false),
  cleanMuonJet     (iConfig.existsAs<bool>("cleanMuonJet") ? iConfig.getParameter<bool>("cleanMuonJet") : false),
  cleanElectronJet (iConfig.existsAs<bool>("cleanElectronJet") ? iConfig.getParameter<bool>("cleanElectronJet") : false),
  dRCleaningAK4    (iConfig.existsAs<double>("dRCleaningAK4") ? iConfig.getParameter<double>("dRCleaningAK4") : 0.4){
  
  muonsToken     = iC.consumes<pat::MuonRefVector> (muonsTag);
  electronsToken = iC.consumes<pat::ElectronRefVector> (electronsTag);
  tausCollectionToken = iC.consumes<pat::TauCollection> (tausCollection);
  tausVLNewToken = iC.consumes<pat::TauRefVector> (tausVLNewTag);
  tausVLOldToken = iC.consumes<pat::TauRefVector> (tausVLOldTag);
  tausRawNewToken = iC.consumes<pat::TauRefVector> (tausRawNewTag);
  tausRawOldToken = iC.consumes<pat::TauRefVector> (tausRawOldTag);
  tausTightNewToken = iC.consumes<pat::TauRefVector> (tausTightNewTag);
  tausTightOldToken = iC.consumes<pat::TauRefVector> (tausTightOldTag);

  tree_ = tree;
  DeclareAndSetBranches();
    
}

/////
void TauTreeFiller::initBranches(){

  combinetaupt.clear(); combinetaueta.clear(); combinetauphi.clear(); combinetaum.clear();
  combinetauGenpt.clear(); combinetauGeneta.clear(); combinetauGenphi.clear(); combinetauGenm.clear();
  combinetaupid.clear(); combinetauidnew.clear(); combinetauidold.clear();
  
  ntausold = 0; ntaus = 0; ntausraw = 0; ntausrawold = 0;

}

/////
bool TauTreeFiller::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace pat;

  this->initBranches();

  Handle<pat::MuonRefVector> muonsH;
  iEvent.getByToken(muonsToken, muonsH);
  pat::MuonRefVector muons = *muonsH;

  Handle<pat::ElectronRefVector> electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  pat::ElectronRefVector electrons = *electronsH;

  Handle<pat::TauCollection> tausH;
  iEvent.getByToken(tausCollectionToken, tausH);
  pat::TauCollection taus = *tausH;

  Handle<pat::TauRefVector > tausVLNewH;
  iEvent.getByToken(tausVLNewToken, tausVLNewH);
  pat::TauRefVector tausVLNew = *tausVLNewH;

  Handle<pat::TauRefVector > tausVLOldH;
  iEvent.getByToken(tausVLOldToken, tausVLOldH);
  pat::TauRefVector tausVLOld = *tausVLOldH;

  Handle<pat::TauRefVector > tausRawNewH;
  iEvent.getByToken(tausRawNewToken, tausRawNewH);
  pat::TauRefVector tausRawNew = *tausRawNewH;

  Handle<pat::TauRefVector > tausRawOldH;
  iEvent.getByToken(tausRawOldToken, tausRawOldH);
  pat::TauRefVector tausRawOld = *tausRawOldH;

  Handle<pat::TauRefVector > tausTightNewH;
  iEvent.getByToken(tausTightNewToken, tausTightNewH);
  pat::TauRefVector tausTightNew = *tausTightNewH;

  Handle<pat::TauRefVector > tausTightOldH;
  iEvent.getByToken(tausTightOldToken, tausTightOldH);
  pat::TauRefVector tausTightOld = *tausTightOldH;

  if(tausVLOldH.isValid()){
    for(std::size_t itau =0 ; itau < tausVLOld.size(); itau++){
      bool skiptau = false;
      for (std::size_t j = 0; j < muons.size(); j++) {
	if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), tausVLOld[itau]->eta(), tausVLOld[itau]->phi()) < dRCleaningAK4) skiptau = true;
      }
      for (std::size_t j = 0; j < electrons.size(); j++) {
	if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), tausVLOld[itau]->eta(), tausVLOld[itau]->phi()) < dRCleaningAK4) skiptau = true;
      }
      if(skiptau) continue;
      ntausold++;
    }
  }

  if(tausVLNewH.isValid()){
    for(std::size_t itau =0 ; itau < tausVLNew.size(); itau++){
      bool skiptau = false;
      for (std::size_t j = 0; j < muons.size(); j++) {
	if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), tausVLNew[itau]->eta(), tausVLNew[itau]->phi()) < dRCleaningAK4) skiptau = true;
      }
      for (std::size_t j = 0; j < electrons.size(); j++) {
	if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), tausVLNew[itau]->eta(), tausVLNew[itau]->phi()) < dRCleaningAK4) skiptau = true;
      }
      if(skiptau) continue;
      ntaus++;
    }
  }

  if(tausRawNewH.isValid()){
    for(std::size_t itau =0 ; itau < tausRawNew.size(); itau++){
      bool skiptau = false;
      for (std::size_t j = 0; j < muons.size(); j++) {
	if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), tausRawNew[itau]->eta(), tausRawNew[itau]->phi()) < dRCleaningAK4) skiptau = true;
      }
      for (std::size_t j = 0; j < electrons.size(); j++) {
	if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), tausRawNew[itau]->eta(), tausRawNew[itau]->phi()) < dRCleaningAK4) skiptau = true;
      }
      if(skiptau) continue;
      ntausraw++;
    }
  }

  if(tausRawOldH.isValid()){
    for(std::size_t itau =0 ; itau < tausRawOld.size(); itau++){
      bool skiptau = false;
      for (std::size_t j = 0; j < muons.size(); j++) {
	if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), tausRawOld[itau]->eta(), tausRawOld[itau]->phi()) < dRCleaningAK4) skiptau = true;
      }
      for (std::size_t j = 0; j < electrons.size(); j++) {
	if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), tausRawOld[itau]->eta(), tausRawOld[itau]->phi()) < dRCleaningAK4) skiptau = true;
      }
      if(skiptau) continue;
      ntausrawold++;
    }
  }


  // store tau-info
  if(tausH.isValid()){ // benchmark collection
    for (auto taus_iter = tausH->begin(); taus_iter != tausH->end(); ++taus_iter) {
      
      //clean from leptons                                                                                                                                                                        
      bool skipjet = false;
      for (std::size_t j = 0; j < muons.size(); j++) {
	if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), taus_iter->eta(), taus_iter->phi()) < dRCleaningAK4)
	  skipjet = true;
      }

      for (std::size_t j = 0; j < electrons.size(); j++) {
	if (cleanMuonJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), taus_iter->eta(), taus_iter->phi()) < dRCleaningAK4)
	  skipjet = true;
      }
      
      if(skipjet) continue;
      
      // dump tau information
      combinetaupt.push_back(taus_iter->pt());
      combinetaueta.push_back(taus_iter->eta());
      combinetauphi.push_back(taus_iter->phi());
      combinetaum.push_back(taus_iter->mass());
      combinetaupid.push_back(taus_iter->pdgId());

      // save the id information
      pat::TauRef tauref = pat::TauRef(tausH, taus_iter - tausH->begin());

      combinetauidold.push_back(0);
      if(std::find(tausVLOld.begin(),tausVLOld.end(),tauref) != tausVLOld.end())
	combinetauidold.back() += 1;
      if(std::find(tausTightOld.begin(),tausTightOld.end(),tauref) != tausTightOld.end())
	combinetauidold.back() += 1;
      
      combinetauidnew.push_back(0);
      if(std::find(tausVLNew.begin(),tausVLNew.end(),tauref) != tausVLNew.end())
	combinetauidnew.back() += 1;
      if(std::find(tausTightNew.begin(),tausTightNew.end(),tauref) != tausTightNew.end())
	combinetauidnew.back() += 1;

      if(taus_iter->genJet()){
	combinetauGenpt.push_back(taus_iter->genJet()->pt()); 
	combinetauGeneta.push_back(taus_iter->genJet()->eta());
	combinetauGenphi.push_back(taus_iter->genJet()->phi());
	combinetauGenm.push_back(taus_iter->genJet()->mass());
      }
      else{
	combinetauGenpt.push_back(0); 
	combinetauGeneta.push_back(0); 
	combinetauGenphi.push_back(0); 
	combinetauGenm.push_back(0);
      }
    }
  }
    
  return true;  
}


////////
void TauTreeFiller::DeclareAndSetBranches(){

  tree_->Branch("ntaus"                , &ntaus                , "ntaus/i");
  tree_->Branch("ntausold"             , &ntausold             , "ntausold/i");
  if(not isTriggerTree){
    tree_->Branch("ntausraw"             , &ntausraw             , "ntausraw/i");
    tree_->Branch("ntausrawold"          , &ntausrawold          , "ntausrawold/i");
  }
  
  if(not isTriggerTree and not isPhotonPurity and not isQCDTree and not applyDiMuonFilter and not applyDiElectronFilter and not applyPhotonJetsFilter){
    tree_->Branch("combinetaupt", "std::vector<float>", &combinetaupt);
    tree_->Branch("combinetaueta", "std::vector<float>", &combinetaueta);
    tree_->Branch("combinetauphi", "std::vector<float>", &combinetauphi);
    tree_->Branch("combinetaum", "std::vector<float>", &combinetaum);
    tree_->Branch("combinetauGenpt", "std::vector<float>", &combinetauGenpt);
    tree_->Branch("combinetauGeneta", "std::vector<float>", &combinetauGeneta);
    tree_->Branch("combinetauGenphi", "std::vector<float>", &combinetauGenphi);
    tree_->Branch("combinetauGenm", "std::vector<float>", &combinetauGenm);
    tree_->Branch("combinetaupid", "std::vector<int>", &combinetaupid);
    tree_->Branch("combinetauidnew", "std::vector<int>", &combinetauidnew);
    tree_->Branch("combinetauidold", "std::vector<int>", &combinetauidold);
  }

}

