#include "AnalysisCode/MonoXAnalysis/interface/JetMetDphiTreeFiller.h"

JetMetDphiTreeFiller::JetMetDphiTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree, const bool & isPuppi):
  muonsTag       (iConfig.getParameter<edm::InputTag>("muons")),
  electronsTag   (iConfig.getParameter<edm::InputTag>("electrons")),
  photonsTag     (iConfig.getParameter<edm::InputTag>("photons")),
  t1metTag       (iConfig.getParameter<edm::InputTag>("t1met")),
  t1mumetTag     (iConfig.getParameter<edm::InputTag>("t1mumet")),
  t1elmetTag     (iConfig.getParameter<edm::InputTag>("t1elmet")),
  t1phmetTag     (iConfig.getParameter<edm::InputTag>("t1phmet")),
  t1taumetTag    (iConfig.getParameter<edm::InputTag>("t1taumet")),
  jetsTag        (iConfig.getParameter<edm::InputTag>("jets")),
  jetsJESUpTag   (iConfig.existsAs<edm::InputTag>("jetsJESUp") ? iConfig.getParameter<edm::InputTag>("jetsJESUp") : edm::InputTag("")),
  jetsJESDwTag   (iConfig.existsAs<edm::InputTag>("jetsJESDw") ? iConfig.getParameter<edm::InputTag>("jetsJESDw") : edm::InputTag("")),
  jetsJERTag     (iConfig.existsAs<edm::InputTag>("jetsJER") ? iConfig.getParameter<edm::InputTag>("jetsJER") : edm::InputTag("")),
  isMC           (iConfig.existsAs<bool>("isMC") ? iConfig.getParameter<bool>("isMC") : false),
  isTriggerTree  (iConfig.existsAs<bool>("isTriggerTree") ? iConfig.getParameter<bool>("isTriggerTree") : false),
  isQCDTree      (iConfig.existsAs<bool>("isQCDTree") ? iConfig.getParameter<bool>("isQCDTree") : false),
  isPhotonPurity (iConfig.existsAs<bool>("isPhotonPurity") ? iConfig.getParameter<bool>("isPhotonPurity") : false),
  cleanMuonJet   (iConfig.existsAs<bool>("cleanMuonJet") ? iConfig.getParameter<bool>("cleanMuonJet") : false),
  cleanElectronJet (iConfig.existsAs<bool>("cleanElectronJet") ? iConfig.getParameter<bool>("cleanElectronJet") : false),
  cleanPhotonJet (iConfig.existsAs<bool>("cleanPhotonJet") ? iConfig.getParameter<bool>("cleanPhotonJet") : false),
  dRCleaningAK4  (iConfig.existsAs<double>("dRCleaningAK4") ? iConfig.getParameter<double>("dRCleaningAK4") : 0.4),
  jetidwp        (iConfig.existsAs<std::string>("jetidwp") ? iConfig.getParameter<std::string>("jetidwp") : "loose"),
  minJetPtCountAK4 (iConfig.existsAs<double>("minJetPtCountAK4") ? iConfig.getParameter<double>("minJetPtCountAK4") : 30),
  addMETSystematics(iConfig.existsAs<bool>("addMETSystematics") ? iConfig.getParameter<bool>("addMETSystematics") : false){
  
  isPuppi_ = isPuppi;

  if(isPuppi_){
    addMETSystematics = (iConfig.existsAs<bool>("addPuppiMETSystematics") ? iConfig.getParameter<bool>("addPuppiMETSystematics") : false);
    t1metTag   = iConfig.getParameter<edm::InputTag>("puppit1met");
    t1mumetTag = iConfig.getParameter<edm::InputTag>("puppit1mumet");
    t1elmetTag = iConfig.getParameter<edm::InputTag>("puppit1elmet");
    t1phmetTag = iConfig.getParameter<edm::InputTag>("puppit1phmet");
    t1taumetTag = iConfig.getParameter<edm::InputTag>("puppit1taumet");
  }

  t1metToken    = iC.consumes<edm::View<pat::MET>  > (t1metTag);
  t1mumetToken  = iC.consumes<edm::View<pat::MET>  > (t1mumetTag);
  t1elmetToken  = iC.consumes<edm::View<pat::MET>   > (t1elmetTag);
  t1phmetToken  = iC.consumes<edm::View<pat::MET>  > (t1phmetTag);

  jetsToken = iC.consumes<std::vector<pat::Jet> > (jetsTag);
  if(jetsJESUpTag.label() != "")
    jetsJESUpToken = iC.consumes<std::vector<pat::Jet> > (jetsJESUpTag);
  if(jetsJESDwTag.label() != "")
    jetsJESDwToken = iC.consumes<std::vector<pat::Jet> > (jetsJESDwTag);
  if(jetsJERTag.label() != "")
    jetsJERToken   = iC.consumes<std::vector<pat::Jet> > (jetsJERTag);
  
  muonsToken     = iC.consumes<pat::MuonRefVector> (muonsTag);
  electronsToken = iC.consumes<pat::ElectronRefVector> (electronsTag);
  photonsToken   = iC.consumes<pat::PhotonRefVector> (photonsTag);

  tree_ = tree;

  this->DeclareAndSetBranches();
  this->initBranches();
}

/////
void JetMetDphiTreeFiller::initBranches(){

  incjetmetdphimin   = 0.0; incjetmumetdphimin  = 0.0; incjetelmetdphimin  = 0.0; incjetphmetdphimin  = 0.0;
  incjetmetdphimin4  = 0.0; incjetmumetdphimin4 = 0.0; incjetelmetdphimin4 = 0.0; incjetphmetdphimin4 = 0.0;
  alljetmetdphimin   = 0.0; alljetmetdphimin4   = 0.0; alljetmumetdphimin  = 0.0; alljetmumetdphimin4 = 0.0;
  alljetelmetdphimin = 0.0; alljetelmetdphimin4 = 0.0; alljetphmetdphimin  = 0.0; alljetphmetdphimin4 = 0.0;

  incjetmetdphimin4up  = 0.0; incjetmumetdphimin4up  = 0.0; incjetelmetdphimin4up  = 0.0; incjetphmetdphimin4up  = 0.0;
  incjetmetdphimin4dw  = 0.0; incjetmumetdphimin4dw  = 0.0; incjetelmetdphimin4dw  = 0.0; incjetphmetdphimin4dw  = 0.0;
  incjetmetdphimin4jer = 0.0; incjetmumetdphimin4jer = 0.0; incjetelmetdphimin4jer = 0.0; incjetphmetdphimin4jer = 0.0;

  incPuppijetmetdphimin   = 0.0; incPuppijetmumetdphimin  = 0.0; incPuppijetelmetdphimin  = 0.0; incPuppijetphmetdphimin  = 0.0;
  incPuppijetmetdphimin4  = 0.0; incPuppijetmumetdphimin4 = 0.0; incPuppijetelmetdphimin4 = 0.0; incPuppijetphmetdphimin4 = 0.0;
  allPuppijetmetdphimin   = 0.0; allPuppijetmetdphimin4   = 0.0; allPuppijetmumetdphimin  = 0.0; allPuppijetmumetdphimin4 = 0.0;
  allPuppijetelmetdphimin = 0.0; allPuppijetelmetdphimin4 = 0.0; allPuppijetphmetdphimin  = 0.0; allPuppijetphmetdphimin4 = 0.0;

  incPuppijetmetdphimin4up  = 0.0; incPuppijetmumetdphimin4up  = 0.0; incPuppijetelmetdphimin4up  = 0.0; incPuppijetphmetdphimin4up  = 0.0;
  incPuppijetmetdphimin4dw  = 0.0; incPuppijetmumetdphimin4dw  = 0.0; incPuppijetelmetdphimin4dw  = 0.0; incPuppijetphmetdphimin4dw  = 0.0;
  incPuppijetmetdphimin4jer = 0.0; incPuppijetmumetdphimin4jer = 0.0; incPuppijetelmetdphimin4jer = 0.0; incPuppijetphmetdphimin4jer = 0.0;

  
}

/////
bool JetMetDphiTreeFiller::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

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

  Handle<pat::PhotonRefVector> photonsH;
  iEvent.getByToken(photonsToken, photonsH);
  pat::PhotonRefVector photons = *photonsH;

  Handle<vector<pat::Jet> > jetsH;
  Handle<vector<pat::Jet> > jetsJESUpH;
  Handle<vector<pat::Jet> > jetsJESDwH;
  Handle<vector<pat::Jet> > jetsJERH;
  iEvent.getByToken(jetsToken, jetsH);
  if(jetsJESUpTag.label() != "")
    iEvent.getByToken(jetsJESUpToken, jetsJESUpH);
  if(jetsJESDwTag.label() != "")
    iEvent.getByToken(jetsJESDwToken, jetsJESDwH);
  if(jetsJERTag.label() != "")
    iEvent.getByToken(jetsJERToken, jetsJERH);

  vector<pat::JetRef> alljets;
  vector<pat::JetRef> incjets;
  vector<pat::JetRef> jets;
  vector<pat::JetRef> incjets_jesup;
  vector<pat::JetRef> incjets_jesdw;
  vector<pat::JetRef> incjets_jer;

  if(addMETSystematics){
    if(jetsJESUpH.isValid())
      fillJetCollections(jetsJESUpH,muons,electrons,photons,incjets_jesup,alljets,isPuppi_);
    alljets.clear();
    if(jetsJESDwH.isValid())
      fillJetCollections(jetsJESDwH,muons,electrons,photons,incjets_jesdw,alljets,isPuppi_);
    alljets.clear();
    if(jetsJERH.isValid())
      fillJetCollections(jetsJERH,muons,electrons,photons,incjets_jer,alljets,isPuppi_);
    alljets.clear();
  }
  if(jetsH.isValid())
    fillJetCollections(jetsH,muons,electrons,photons,incjets,alljets,isPuppi_);

  Handle<View<pat::MET> > t1metH;
  iEvent.getByToken(t1metToken, t1metH);
  
  Handle<View<pat::MET> > t1mumetH;
  iEvent.getByToken(t1mumetToken, t1mumetH);

  Handle<View<pat::MET> > t1elmetH;
  iEvent.getByToken(t1elmetToken, t1elmetH);

  Handle<View<pat::MET> > t1phmetH;
  iEvent.getByToken(t1phmetToken, t1phmetH);

  // start filling info
  std::vector<float> alljetmetdphiminvector;
  std::vector<float> alljetmetdphimin4vector;
  std::vector<float> alljetmumetdphiminvector;
  std::vector<float> alljetmumetdphimin4vector;
  std::vector<float> alljetelmetdphiminvector;
  std::vector<float> alljetelmetdphimin4vector;
  std::vector<float> alljetphmetdphiminvector;
  std::vector<float> alljetphmetdphimin4vector;

  for (size_t i = 0; i < alljets.size(); i++) {
    if (alljets[i]->pt() > minJetPtCountAK4) {  
      float alljetphi = atan2(sin(alljets[i]->phi()), cos(alljets[i]->phi()));
      alljetmetdphiminvector  .push_back(fabs(deltaPhi(alljetphi, t1metH->front().corPhi())));
      alljetmumetdphiminvector.push_back(fabs(deltaPhi(alljetphi, t1mumetH->front().corPhi())));
      alljetelmetdphiminvector.push_back(fabs(deltaPhi(alljetphi, t1elmetH->front().corPhi())));
      alljetphmetdphiminvector.push_back(fabs(deltaPhi(alljetphi, t1phmetH->front().corPhi())));
      if (i < 4) alljetmetdphimin4vector  .push_back(fabs(deltaPhi(alljetphi, t1metH->front().corPhi())));
      if (i < 4) alljetmumetdphimin4vector.push_back(fabs(deltaPhi(alljetphi, t1mumetH->front().corPhi())));
      if (i < 4) alljetelmetdphimin4vector.push_back(fabs(deltaPhi(alljetphi, t1elmetH->front().corPhi())));
      if (i < 4) alljetphmetdphimin4vector.push_back(fabs(deltaPhi(alljetphi, t1phmetH->front().corPhi())));
    }
  }
  
  if(not isPuppi_){ // no ID requirements on the jet --> pnly pt-cut
    if (alljetmetdphiminvector   .size() > 0) alljetmetdphimin    = *min_element(alljetmetdphiminvector   .begin(), alljetmetdphiminvector   .end());
    if (alljetmumetdphiminvector .size() > 0) alljetmumetdphimin  = *min_element(alljetmumetdphiminvector .begin(), alljetmumetdphiminvector .end());
    if (alljetelmetdphiminvector .size() > 0) alljetelmetdphimin  = *min_element(alljetelmetdphiminvector .begin(), alljetelmetdphiminvector .end());
    if (alljetphmetdphiminvector .size() > 0) alljetphmetdphimin  = *min_element(alljetphmetdphiminvector .begin(), alljetphmetdphiminvector .end());
    if (alljetmetdphimin4vector  .size() > 0) alljetmetdphimin4   = *min_element(alljetmetdphimin4vector  .begin(), alljetmetdphimin4vector  .end());
    if (alljetmumetdphimin4vector.size() > 0) alljetmumetdphimin4 = *min_element(alljetmumetdphimin4vector.begin(), alljetmumetdphimin4vector.end());
    if (alljetelmetdphimin4vector.size() > 0) alljetelmetdphimin4 = *min_element(alljetelmetdphimin4vector.begin(), alljetelmetdphimin4vector.end());
    if (alljetphmetdphimin4vector.size() > 0) alljetphmetdphimin4 = *min_element(alljetphmetdphimin4vector.begin(), alljetphmetdphimin4vector.end());
  }
  else{
    if (alljetmetdphiminvector   .size() > 0) allPuppijetmetdphimin    = *min_element(alljetmetdphiminvector   .begin(), alljetmetdphiminvector   .end());
    if (alljetmumetdphiminvector .size() > 0) allPuppijetmumetdphimin  = *min_element(alljetmumetdphiminvector .begin(), alljetmumetdphiminvector .end());
    if (alljetelmetdphiminvector .size() > 0) allPuppijetelmetdphimin  = *min_element(alljetelmetdphiminvector .begin(), alljetelmetdphiminvector .end());
    if (alljetphmetdphiminvector .size() > 0) allPuppijetphmetdphimin  = *min_element(alljetphmetdphiminvector .begin(), alljetphmetdphiminvector .end());
    if (alljetmetdphimin4vector  .size() > 0) allPuppijetmetdphimin4   = *min_element(alljetmetdphimin4vector  .begin(), alljetmetdphimin4vector  .end());
    if (alljetmumetdphimin4vector.size() > 0) allPuppijetmumetdphimin4 = *min_element(alljetmumetdphimin4vector.begin(), alljetmumetdphimin4vector.end());
    if (alljetelmetdphimin4vector.size() > 0) allPuppijetelmetdphimin4 = *min_element(alljetelmetdphimin4vector.begin(), alljetelmetdphimin4vector.end());
    if (alljetphmetdphimin4vector.size() > 0) allPuppijetphmetdphimin4 = *min_element(alljetphmetdphimin4vector.begin(), alljetphmetdphimin4vector.end());

  }
  

  // delta phi jet-met      
  std::vector<float> incjetmetdphiminvector;
  std::vector<float> incjetmetdphimin4vector;
  std::vector<float> incjetmumetdphiminvector;
  std::vector<float> incjetmumetdphimin4vector;
  std::vector<float> incjetelmetdphiminvector;
  std::vector<float> incjetelmetdphimin4vector;
  std::vector<float> incjetphmetdphiminvector;
  std::vector<float> incjetphmetdphimin4vector;

  for (size_t i = 0; i < incjets.size(); i++) {
    if (incjets[i]->pt() > minJetPtCountAK4) {
      float incjetphi = atan2(sin(incjets[i]->phi()), cos(incjets[i]->phi()));
      incjetmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1metH->front().corPhi())));
      incjetmumetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1mumetH->front().corPhi())));
      incjetelmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1elmetH->front().corPhi())));
      incjetphmetdphiminvector.push_back(fabs(deltaPhi(incjetphi, t1phmetH->front().corPhi())));
      if (i < 4) incjetmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1metH->front().corPhi())));
      if (i < 4) incjetmumetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1mumetH->front().corPhi())));
      if (i < 4) incjetelmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1elmetH->front().corPhi())));
      if (i < 4) incjetphmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1phmetH->front().corPhi())));
    }
  }

  if(not isPuppi_){
    
    if (incjetmetdphiminvector .size() > 0) incjetmetdphimin  = *min_element(incjetmetdphiminvector .begin(), incjetmetdphiminvector .end());
    if (incjetmetdphimin4vector.size() > 0) incjetmetdphimin4 = *min_element(incjetmetdphimin4vector.begin(), incjetmetdphimin4vector.end());
    if (incjetmumetdphiminvector .size() > 0) incjetmumetdphimin  = *min_element(incjetmumetdphiminvector .begin(), incjetmumetdphiminvector .end());
    if (incjetmumetdphimin4vector.size() > 0) incjetmumetdphimin4 = *min_element(incjetmumetdphimin4vector.begin(), incjetmumetdphimin4vector.end());
    if (incjetelmetdphiminvector .size() > 0) incjetelmetdphimin  = *min_element(incjetelmetdphiminvector .begin(), incjetelmetdphiminvector .end());
    if (incjetelmetdphimin4vector.size() > 0) incjetelmetdphimin4 = *min_element(incjetelmetdphimin4vector.begin(), incjetelmetdphimin4vector.end());
    if (incjetphmetdphiminvector .size() > 0) incjetphmetdphimin  = *min_element(incjetphmetdphiminvector .begin(), incjetphmetdphiminvector .end());
    if (incjetphmetdphimin4vector.size() > 0) incjetphmetdphimin4 = *min_element(incjetphmetdphimin4vector.begin(), incjetphmetdphimin4vector.end());
  }
  else{
    if (incjetmetdphiminvector .size() > 0) incPuppijetmetdphimin  = *min_element(incjetmetdphiminvector .begin(), incjetmetdphiminvector .end());
    if (incjetmetdphimin4vector.size() > 0) incPuppijetmetdphimin4 = *min_element(incjetmetdphimin4vector.begin(), incjetmetdphimin4vector.end());
    if (incjetmumetdphiminvector .size() > 0) incPuppijetmumetdphimin  = *min_element(incjetmumetdphiminvector .begin(), incjetmumetdphiminvector .end());
    if (incjetmumetdphimin4vector.size() > 0) incPuppijetmumetdphimin4 = *min_element(incjetmumetdphimin4vector.begin(), incjetmumetdphimin4vector.end());
    if (incjetelmetdphiminvector .size() > 0) incPuppijetelmetdphimin  = *min_element(incjetelmetdphiminvector .begin(), incjetelmetdphiminvector .end());
    if (incjetelmetdphimin4vector.size() > 0) incPuppijetelmetdphimin4 = *min_element(incjetelmetdphimin4vector.begin(), incjetelmetdphimin4vector.end());
    if (incjetphmetdphiminvector .size() > 0) incPuppijetphmetdphimin  = *min_element(incjetphmetdphiminvector .begin(), incjetphmetdphiminvector .end());
    if (incjetphmetdphimin4vector.size() > 0) incPuppijetphmetdphimin4 = *min_element(incjetphmetdphimin4vector.begin(), incjetphmetdphimin4vector.end());
  }
  
  // systematics
  // delta phi jet-met      
  if(addMETSystematics){

    incjetmetdphimin4vector.clear();
    incjetmumetdphimin4vector.clear();
    incjetelmetdphimin4vector.clear();
    incjetphmetdphimin4vector.clear();
    
    for (size_t i = 0; i < incjets_jesup.size(); i++) {
      if (incjets_jesup[i]->pt() > minJetPtCountAK4) {
	float incjetphi = atan2(sin(incjets_jesup[i]->phi()), cos(incjets_jesup[i]->phi()));
	if (i < 4) incjetmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetEnUp,   pat::MET::METCorrectionLevel::Type1))));
	if (i < 4) incjetmumetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1mumetH->front().shiftedPhi(pat::MET::METUncertainty::JetEnUp,   pat::MET::METCorrectionLevel::Type1))));
	if (i < 4) incjetelmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1elmetH->front().shiftedPhi(pat::MET::METUncertainty::JetEnUp,   pat::MET::METCorrectionLevel::Type1))));
	if (i < 4) incjetphmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1phmetH->front().shiftedPhi(pat::MET::METUncertainty::JetEnUp,   pat::MET::METCorrectionLevel::Type1))));
      }
    }
    
    if(not isPuppi_){
      if (incjetmetdphimin4vector.size() > 0)   incjetmetdphimin4up   = *min_element(incjetmetdphimin4vector.begin(), incjetmetdphimin4vector.end());
      if (incjetmumetdphimin4vector.size() > 0) incjetmumetdphimin4up = *min_element(incjetmumetdphimin4vector.begin(), incjetmumetdphimin4vector.end());
      if (incjetelmetdphimin4vector.size() > 0) incjetelmetdphimin4up = *min_element(incjetelmetdphimin4vector.begin(), incjetelmetdphimin4vector.end());
      if (incjetphmetdphimin4vector.size() > 0) incjetphmetdphimin4up = *min_element(incjetphmetdphimin4vector.begin(), incjetphmetdphimin4vector.end());
    }
    else{
      if (incjetmetdphimin4vector.size() > 0)   incPuppijetmetdphimin4up   = *min_element(incjetmetdphimin4vector.begin(), incjetmetdphimin4vector.end());
      if (incjetmumetdphimin4vector.size() > 0) incPuppijetmumetdphimin4up = *min_element(incjetmumetdphimin4vector.begin(), incjetmumetdphimin4vector.end());
      if (incjetelmetdphimin4vector.size() > 0) incPuppijetelmetdphimin4up = *min_element(incjetelmetdphimin4vector.begin(), incjetelmetdphimin4vector.end());
      if (incjetphmetdphimin4vector.size() > 0) incPuppijetphmetdphimin4up = *min_element(incjetphmetdphimin4vector.begin(), incjetphmetdphimin4vector.end());
    }

    
    incjetmetdphimin4vector.clear();
    incjetmumetdphimin4vector.clear();
    incjetelmetdphimin4vector.clear();
    incjetphmetdphimin4vector.clear();
    
    for (size_t i = 0; i < incjets_jesdw.size(); i++) {
      if (incjets_jesdw[i]->pt() > minJetPtCountAK4) {
	float incjetphi = atan2(sin(incjets_jesdw[i]->phi()), cos(incjets_jesdw[i]->phi()));
	if (i < 4) incjetmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetEnDown,   pat::MET::METCorrectionLevel::Type1))));
	if (i < 4) incjetmumetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1mumetH->front().shiftedPhi(pat::MET::METUncertainty::JetEnDown,   pat::MET::METCorrectionLevel::Type1))));
	if (i < 4) incjetelmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1elmetH->front().shiftedPhi(pat::MET::METUncertainty::JetEnDown,   pat::MET::METCorrectionLevel::Type1))));
	if (i < 4) incjetphmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1phmetH->front().shiftedPhi(pat::MET::METUncertainty::JetEnDown,   pat::MET::METCorrectionLevel::Type1))));
      }
    }
    
    if(not isPuppi_){
      if (incjetmetdphimin4vector.size() > 0)   incjetmetdphimin4dw   = *min_element(incjetmetdphimin4vector.begin(), incjetmetdphimin4vector.end());
      if (incjetmumetdphimin4vector.size() > 0) incjetmumetdphimin4dw = *min_element(incjetmumetdphimin4vector.begin(), incjetmumetdphimin4vector.end());
      if (incjetelmetdphimin4vector.size() > 0) incjetelmetdphimin4dw = *min_element(incjetelmetdphimin4vector.begin(), incjetelmetdphimin4vector.end());
      if (incjetphmetdphimin4vector.size() > 0) incjetphmetdphimin4dw = *min_element(incjetphmetdphimin4vector.begin(), incjetphmetdphimin4vector.end());
    }
    else{
      if (incjetmetdphimin4vector.size() > 0)   incPuppijetmetdphimin4dw   = *min_element(incjetmetdphimin4vector.begin(), incjetmetdphimin4vector.end());
      if (incjetmumetdphimin4vector.size() > 0) incPuppijetmumetdphimin4dw = *min_element(incjetmumetdphimin4vector.begin(), incjetmumetdphimin4vector.end());
      if (incjetelmetdphimin4vector.size() > 0) incPuppijetelmetdphimin4dw = *min_element(incjetelmetdphimin4vector.begin(), incjetelmetdphimin4vector.end());
      if (incjetphmetdphimin4vector.size() > 0) incPuppijetphmetdphimin4dw = *min_element(incjetphmetdphimin4vector.begin(), incjetphmetdphimin4vector.end());
    }

    //////    
    incjetmetdphimin4vector.clear();
    incjetmumetdphimin4vector.clear();
    incjetelmetdphimin4vector.clear();
    incjetphmetdphimin4vector.clear();
    
    std::cout<<"T1 met "<<t1metH->front().corPt()<<" "<<t1metH->front().corPhi()<<" smear "<<t1metH->front().shiftedPt(pat::MET::METUncertainty::NoShift, pat::MET::METCorrectionLevel::Type1Smear)  <<" "<<t1metH->front().shiftedPhi(pat::MET::METUncertainty::NoShift,   pat::MET::METCorrectionLevel::Type1Smear)<<std::endl;
    
    for (size_t i = 0; i < incjets_jer.size(); i++) {
      if (incjets_jer[i]->pt() > minJetPtCountAK4) {
	float incjetphi = atan2(sin(incjets_jer[i]->phi()), cos(incjets_jer[i]->phi()));
	if (i < 4) incjetmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1metH->front().shiftedPhi(pat::MET::METUncertainty::NoShift,   pat::MET::METCorrectionLevel::Type1Smear))));
	if (i < 4) incjetmumetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1mumetH->front().shiftedPhi(pat::MET::METUncertainty::NoShift,   pat::MET::METCorrectionLevel::Type1Smear))));
	if (i < 4) incjetelmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1elmetH->front().shiftedPhi(pat::MET::METUncertainty::NoShift,   pat::MET::METCorrectionLevel::Type1Smear))));
	if (i < 4) incjetphmetdphimin4vector.push_back(fabs(deltaPhi(incjetphi, t1phmetH->front().shiftedPhi(pat::MET::METUncertainty::NoShift,   pat::MET::METCorrectionLevel::Type1Smear))));
      }
    }
    
    if(not isPuppi_){
      if (incjetmetdphimin4vector.size() > 0)   incjetmetdphimin4jer   = *min_element(incjetmetdphimin4vector.begin(), incjetmetdphimin4vector.end());
      if (incjetmumetdphimin4vector.size() > 0) incjetmumetdphimin4jer = *min_element(incjetmumetdphimin4vector.begin(), incjetmumetdphimin4vector.end());
      if (incjetelmetdphimin4vector.size() > 0) incjetelmetdphimin4jer = *min_element(incjetelmetdphimin4vector.begin(), incjetelmetdphimin4vector.end());
      if (incjetphmetdphimin4vector.size() > 0) incjetphmetdphimin4jer = *min_element(incjetphmetdphimin4vector.begin(), incjetphmetdphimin4vector.end());
    }
    else{
      if (incjetmetdphimin4vector.size() > 0)   incPuppijetmetdphimin4jer   = *min_element(incjetmetdphimin4vector.begin(), incjetmetdphimin4vector.end());
      if (incjetmumetdphimin4vector.size() > 0) incPuppijetmumetdphimin4jer = *min_element(incjetmumetdphimin4vector.begin(), incjetmumetdphimin4vector.end());
      if (incjetelmetdphimin4vector.size() > 0) incPuppijetelmetdphimin4jer = *min_element(incjetelmetdphimin4vector.begin(), incjetelmetdphimin4vector.end());
      if (incjetphmetdphimin4vector.size() > 0) incPuppijetphmetdphimin4jer = *min_element(incjetphmetdphimin4vector.begin(), incjetphmetdphimin4vector.end());      
    }
  }
  return true;
}

//// fill jet collection
void JetMetDphiTreeFiller::fillJetCollections(const edm::Handle<std::vector<pat::Jet> > & jetsH, 
					      const pat::MuonRefVector & muons, 
					      const pat::ElectronRefVector & electrons,
					      const pat::PhotonRefVector & photons, 
					      std::vector<pat::JetRef> & incjets, 
					      std::vector<pat::JetRef> & alljets, 
					      const bool & ispuppi){
  
  if(jetsH.isValid()){      
    for (auto jets_iter = jetsH->begin(); jets_iter != jetsH->end(); ++jets_iter) {
      //clean from leptons
      bool skipjet = false;
      for (std::size_t j = 0; j < muons.size(); j++) {
	if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK4) 
	  skipjet = true;
      }
      for (std::size_t j = 0; j < electrons.size(); j++) {
	if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK4) 
	  skipjet = true;
      }
      for (std::size_t j = 0; j < photons.size(); j++) {
	if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK4) 
	  skipjet = true;
      }
      
      // jet in overlap with lepton
      if (skipjet) continue;
      
      pat::JetRef jetref(jetsH, jets_iter - jetsH->begin());
      if(jetref.isAvailable() and jetref.isNonnull()) alljets.push_back(jetref);
      
      // apply jet id
      bool passjetid = applyJetID(*jets_iter,jetidwp);            
      if (!passjetid) 
	continue;
      
      if(jetref.isAvailable() and jetref.isNonnull())
	incjets.push_back(jetref);
    }
    
    if(incjets.size() > 0) sort(incjets.begin(), incjets.end(), jetPtSorter);
  }  
}


////////
void JetMetDphiTreeFiller::DeclareAndSetBranches(){

  if(not isPuppi_){

    tree_->Branch("incjetmetdphimin4"    , &incjetmetdphimin4    , "incjetmetdphimin4/F");
    tree_->Branch("incjetmumetdphimin4"  , &incjetmumetdphimin4  , "incjetmumetdphimin4/F");
    tree_->Branch("incjetelmetdphimin4"  , &incjetelmetdphimin4  , "incjetelmetdphimin4/F");
    tree_->Branch("incjetphmetdphimin4"  , &incjetphmetdphimin4  , "incjetphmetdphimin4/F");

    tree_->Branch("incjetmetdphimin"     , &incjetmetdphimin     , "incjetmetdphimin/F");
    tree_->Branch("incjetmumetdphimin"   , &incjetmumetdphimin   , "incjetmumetdphimin/F");
    tree_->Branch("incjetelmetdphimin"   , &incjetelmetdphimin   , "incjetelmetdphimin/F");
    tree_->Branch("incjetphmetdphimin"   , &incjetphmetdphimin   , "incjetphmetdphimin/F");
    
    tree_->Branch("alljetmetdphimin"     , &alljetmetdphimin     , "alljetmetdphimin/F");
    tree_->Branch("alljetmumetdphimin"   , &alljetmumetdphimin   , "alljetmumetdphimin/F");
    tree_->Branch("alljetelmetdphimin"   , &alljetelmetdphimin   , "alljetelmetdphimin/F");
    tree_->Branch("alljetphmetdphimin"   , &alljetphmetdphimin   , "alljetphmetdphimin/F");
    
    if(not isQCDTree and not isTriggerTree and jetsJESUpTag.label() != "" and jetsJESDwTag.label() != "" and jetsJERTag.label() != "" and addMETSystematics){
            
      tree_->Branch("incjetmetdphimin4up"    , &incjetmetdphimin4up    , "incjetmetdphimin4up/F");
      tree_->Branch("incjetmumetdphimin4up"  , &incjetmumetdphimin4up  , "incjetmumetdphimin4up/F");
      tree_->Branch("incjetelmetdphimin4up"  , &incjetelmetdphimin4up  , "incjetelmetdphimin4up/F");
      tree_->Branch("incjetphmetdphimin4up"  , &incjetphmetdphimin4up  , "incjetphmetdphimin4up/F");

      tree_->Branch("incjetmetdphimin4dw"    , &incjetmetdphimin4dw    , "incjetmetdphimin4dw/F");
      tree_->Branch("incjetmumetdphimin4dw"  , &incjetmumetdphimin4dw  , "incjetmumetdphimin4dw/F");
      tree_->Branch("incjetelmetdphimin4dw"  , &incjetelmetdphimin4dw  , "incjetelmetdphimin4dw/F");
      tree_->Branch("incjetphmetdphimin4dw"  , &incjetphmetdphimin4dw  , "incjetphmetdphimin4dw/F");

      tree_->Branch("incjetmetdphimin4jer"    , &incjetmetdphimin4jer    , "incjetmetdphimin4jer/F");
      tree_->Branch("incjetmumetdphimin4jer"  , &incjetmumetdphimin4jer  , "incjetmumetdphimin4jer/F");
      tree_->Branch("incjetelmetdphimin4jer"  , &incjetelmetdphimin4jer  , "incjetelmetdphimin4jer/F");
      tree_->Branch("incjetphmetdphimin4jer"  , &incjetphmetdphimin4jer  , "incjetphmetdphimin4jer/F");

      tree_->Branch("alljetmetdphimin4"    , &alljetmetdphimin4    , "alljetmetdphimin4/F");
      tree_->Branch("alljetmumetdphimin4"  , &alljetmumetdphimin4  , "alljetmumetdphimin4/F");
      tree_->Branch("alljetelmetdphimin4"  , &alljetelmetdphimin4  , "alljetelmetdphimin4/F");
      tree_->Branch("alljetphmetdphimin4"  , &alljetphmetdphimin4  , "alljetphmetdphimin4/F");    
    }
  }
  else{

    tree_->Branch("incPuppijetmetdphimin4"    , &incPuppijetmetdphimin4    , "incPuppijetmetdphimin4/F");
    tree_->Branch("incPuppijetmumetdphimin4"  , &incPuppijetmumetdphimin4  , "incPuppijetmumetdphimin4/F");
    tree_->Branch("incPuppijetelmetdphimin4"  , &incPuppijetelmetdphimin4  , "incPuppijetelmetdphimin4/F");
    tree_->Branch("incPuppijetphmetdphimin4"  , &incPuppijetphmetdphimin4  , "incPuppijetphmetdphimin4/F");

    tree_->Branch("incPuppijetmetdphimin"     , &incPuppijetmetdphimin     , "incPuppijetmetdphimin/F");
    tree_->Branch("incPuppijetmumetdphimin"   , &incPuppijetmumetdphimin   , "incPuppijetmumetdphimin/F");
    tree_->Branch("incPuppijetelmetdphimin"   , &incPuppijetelmetdphimin   , "incPuppijetelmetdphimin/F");
    tree_->Branch("incPuppijetphmetdphimin"   , &incPuppijetphmetdphimin   , "incPuppijetphmetdphimin/F");
    
    tree_->Branch("allPuppijetmetdphimin"     , &allPuppijetmetdphimin     , "allPuppijetmetdphimin/F");
    tree_->Branch("allPuppijetmumetdphimin"   , &allPuppijetmumetdphimin   , "allPuppijetmumetdphimin/F");
    tree_->Branch("allPuppijetelmetdphimin"   , &allPuppijetelmetdphimin   , "allPuppijetelmetdphimin/F");
    tree_->Branch("allPuppijetphmetdphimin"   , &allPuppijetphmetdphimin   , "allPuppijetphmetdphimin/F");
    
    if(not isQCDTree and not isTriggerTree and jetsJESUpTag.label() != "" and jetsJESDwTag.label() != "" and jetsJERTag.label() != "" and addMETSystematics){
            
      tree_->Branch("incPuppijetmetdphimin4up"    , &incPuppijetmetdphimin4up    , "incPuppijetmetdphimin4up/F");
      tree_->Branch("incPuppijetmumetdphimin4up"  , &incPuppijetmumetdphimin4up  , "incPuppijetmumetdphimin4up/F");
      tree_->Branch("incPuppijetelmetdphimin4up"  , &incPuppijetelmetdphimin4up  , "incPuppijetelmetdphimin4up/F");
      tree_->Branch("incPuppijetphmetdphimin4up"  , &incPuppijetphmetdphimin4up  , "incPuppijetphmetdphimin4up/F");

      tree_->Branch("incPuppijetmetdphimin4dw"    , &incPuppijetmetdphimin4dw    , "incPuppijetmetdphimin4dw/F");
      tree_->Branch("incPuppijetmumetdphimin4dw"  , &incPuppijetmumetdphimin4dw  , "incPuppijetmumetdphimin4dw/F");
      tree_->Branch("incPuppijetelmetdphimin4dw"  , &incPuppijetelmetdphimin4dw  , "incPuppijetelmetdphimin4dw/F");
      tree_->Branch("incPuppijetphmetdphimin4dw"  , &incPuppijetphmetdphimin4dw  , "incPuppijetphmetdphimin4dw/F");

      tree_->Branch("incPuppijetmetdphimin4jer"    , &incPuppijetmetdphimin4jer    , "incPuppijetmetdphimin4jer/F");
      tree_->Branch("incPuppijetmumetdphimin4jer"  , &incPuppijetmumetdphimin4jer  , "incPuppijetmumetdphimin4jer/F");
      tree_->Branch("incPuppijetelmetdphimin4jer"  , &incPuppijetelmetdphimin4jer  , "incPuppijetelmetdphimin4jer/F");
      tree_->Branch("incPuppijetphmetdphimin4jer"  , &incPuppijetphmetdphimin4jer  , "incPuppijetphmetdphimin4jer/F");

      tree_->Branch("allPuppijetmetdphimin4"    , &allPuppijetmetdphimin4    , "allPuppijetmetdphimin4/F");
      tree_->Branch("allPuppijetmumetdphimin4"  , &allPuppijetmumetdphimin4  , "allPuppijetmumetdphimin4/F");
      tree_->Branch("allPuppijetelmetdphimin4"  , &allPuppijetelmetdphimin4  , "allPuppijetelmetdphimin4/F");
      tree_->Branch("allPuppijetphmetdphimin4"  , &allPuppijetphmetdphimin4  , "allPuppijetphmetdphimin4/F");    
    }

  }
}
