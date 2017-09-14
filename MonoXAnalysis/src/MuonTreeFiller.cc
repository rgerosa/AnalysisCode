#include "AnalysisCode/MonoXAnalysis/interface/MuonTreeFiller.h"

MuonTreeFiller::MuonTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree):
  muonsTag       (iConfig.getParameter<edm::InputTag>("muons")),
  tightmuonsTag  (iConfig.getParameter<edm::InputTag>("tightmuons")),
  highptmuonsTag (iConfig.getParameter<edm::InputTag>("highptmuons")),
  isQCDTree       (iConfig.existsAs<bool>("isQCDTree") ? iConfig.getParameter<bool>("isQCDTree") : false),
  isPhotonPurity  (iConfig.existsAs<bool>("isPhotonPurity") ? iConfig.getParameter<bool>("isPhotonPurity") : false),
  isReMiniAOD     (iConfig.existsAs<bool>("isReMiniAOD") ? iConfig.getParameter<bool>("isReMiniAOD") : false),
  isTriggerTree   (iConfig.existsAs<bool>("isTriggerTree") ? iConfig.getParameter<bool>("isTriggerTree") : false),
  t1metTag    (iConfig.getParameter<edm::InputTag>("t1met")),
  verticesTag (iConfig.getParameter<edm::InputTag>("vertices")){   
  
  t1metToken    = iC.consumes<edm::View<pat::MET> > (t1metTag);
  muonsToken       = iC.consumes<pat::MuonRefVector> (muonsTag);
  tightmuonsToken  = iC.consumes<pat::MuonRefVector> (tightmuonsTag);
  highptmuonsToken = iC.consumes<pat::MuonRefVector> (highptmuonsTag);
  verticesToken    = iC.consumes<std::vector<reco::Vertex> > (verticesTag);
  if(isReMiniAOD)
    fakeMuonCollToken = iC.consumes<edm::View<reco::Candidate> > (iConfig.getParameter<edm::InputTag>("fakeMuonCandidates"));

  tree_ = tree;
  DeclareAndSetBranches();
    
}

/////
void MuonTreeFiller::initBranches(){

  nmuons      = 0;   ntightmuons = 0;   nhighptmuons = 0; 
  mu1pid      = 0;   mu1pt       = 0.0; mu1eta      = 0.0; mu1phi      = 0.0;
  mu1pfpt     = 0.0; mu1pfeta    = 0.0; mu1pfphi    = 0.0; mu1id       = 0;
  mu1idm      = 0;   mu1idt      = 0;   mu1iso      = 0.0;

  mu2pid      = 0;   mu2pt       = 0.0; mu2eta      = 0.0; mu2phi      = 0.0;
  mu2pfpt     = 0.0; mu2pfeta    = 0.0; mu2pfphi    = 0.0; mu2id       = 0;
  mu2idm      = 0;   mu2idt      = 0;   mu2iso      = 0.0;

  zmass = 0.0; zpt = 0.0; zeta = 0.0; zphi = 0.0; wmt = 0.0;
   
}

/////
bool MuonTreeFiller::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace pat;

  this->initBranches();

  //MET in the event
  Handle<View<pat::MET> > t1metH;
  iEvent.getByToken(t1metToken, t1metH);
  
  Handle<vector<Vertex> > verticesH;
  iEvent.getByToken(verticesToken, verticesH);

  Handle<pat::MuonRefVector> muonsH;
  iEvent.getByToken(muonsToken, muonsH);
  pat::MuonRefVector muons = *muonsH;

  Handle<pat::MuonRefVector> tightmuonsH;
  iEvent.getByToken(tightmuonsToken, tightmuonsH);
  pat::MuonRefVector tightmuons = *tightmuonsH;

  Handle<pat::MuonRefVector> highptmuonsH;
  iEvent.getByToken(highptmuonsToken, highptmuonsH);
  pat::MuonRefVector highptmuons = *highptmuonsH;

  Handle<edm::View<reco::Candidate> > fakeMuonsH;
  if(isReMiniAOD)
    iEvent.getByToken(fakeMuonCollToken,fakeMuonsH);

  // muon counters
  vector<pat::MuonRef> muonvector;
  if(muonsH.isValid()){
    nmuons          = muonsH->size();
    for (size_t i = 0; i < muons.size(); i++) 
      muonvector.push_back(muons[i]);
  }
  if(tightmuonsH.isValid())
    ntightmuons     = tightmuonsH->size();
  if(highptmuonsH.isValid())
    nhighptmuons    = highptmuonsH->size(); 
  
  sort(muonvector.begin(), muonvector.end(), muonPtSorter);

  // one or two loose muons
  if (nmuons == 1 || nmuons == 2) {
      
    pat::MuonRef muon = muonvector[0];
    mu1pid   = muon->pdgId(); 
    mu1pt    = muon->pt(); 
    mu1eta   = muon->eta(); 
    mu1phi   = muon->phi();
    mu1pfpt  = muon->pfP4().Pt();
    mu1pfeta = muon->pfP4().Eta();
    mu1pfphi = muon->pfP4().Phi();
    mu1iso   = computeMuonIso(*muon); 
    mu1idm   = (muon::isMediumMuon(*muon) ? 1 : 0);
    if (verticesH->size() > 0) 
      mu1idt = (muon::isTightMuon(*muon, *(verticesH->begin())) ? 1 : 0);
      
    // tight muon
    for (std::size_t i = 0; i < tightmuons.size(); i++) {
      if (muon == tightmuons[i]) 
	mu1id = 1; 
    }
      
    // store high-pt muons that are not tight ones
    for (std::size_t i = 0; i < highptmuons.size(); i++) {
      if (muon == highptmuons[i] and mu1id != 1) 
	mu1id = 2; // high pt muon
    }

    if (nmuons == 1) 
      wmt = sqrt(2.0 * mu1pt * t1metH->front().corPt() * (1.0 - cos(deltaPhi(mu1phi, t1metH->front().corPhi()))));
  }
   
  // two loose muons
  if (nmuons == 2) {        

    pat::MuonRef muon = muonvector[1];
    mu2pid   = muon->pdgId(); 
    mu2pt    = muon->pt(); 
    mu2eta   = muon->eta(); 
    mu2phi   = muon->phi();
    mu2pfpt  = muon->pfP4().Pt();
    mu2pfeta = muon->pfP4().Eta();
    mu2pfphi = muon->pfP4().Phi();
    mu2iso   = computeMuonIso(*muon);       
    mu2idm   = (muon::isMediumMuon(*muon) ? 1 : 0);
    if (verticesH->size() > 0) 
      mu2idt = (muon::isTightMuon(*muon, *(verticesH->begin())) ? 1 : 0);

    // check if belong to the tight / high pt collection
    for (std::size_t i = 0; i < tightmuons.size(); i++) {
      if (muon == tightmuons[i]) 
	mu2id = 1;
    }
      
    // store high-pt muons that are not tight ones
    for (std::size_t i = 0; i < highptmuons.size(); i++) {
      if (muon == highptmuons[i] and mu2id != 1) 
	mu2id = 2;
    }
      
    TLorentzVector mu1vec; 
    mu1vec.SetPtEtaPhiE(mu1pt, mu1eta, mu1phi, muonvector[0]->p());
    TLorentzVector mu2vec; 
    mu2vec.SetPtEtaPhiE(mu2pt, mu2eta, mu2phi, muon->p());
      
    TLorentzVector zvec(mu1vec);
    zvec += mu2vec;
    
    zmass = zvec.M();
    zpt   = zvec.Pt();
    zeta  = zvec.Eta();            
    zphi  = zvec.Phi();
  }

  // fake muons --> reMiniAOD 2016
  nmuonsfake = 0;
  fakemupt.clear();
  fakemueta.clear();
  fakemuphi.clear();
  if(fakeMuonsH.isValid()){
    for(size_t imu = 0; imu < fakeMuonsH->size(); imu++){
      nmuonsfake++;
      fakemupt.push_back(fakeMuonsH->at(imu).pt());
      fakemueta.push_back(fakeMuonsH->at(imu).eta());
      fakemuphi.push_back(fakeMuonsH->at(imu).phi());
    }
  }
  return true;  
}

void MuonTreeFiller::DeclareAndSetBranches(){

  tree_->Branch("nmuons"               , &nmuons               , "nmuons/i");
  tree_->Branch("ntightmuons"          , &ntightmuons          , "ntightmuons/i");
  if(not isTriggerTree)
    tree_->Branch("nhighptmuons"         , &nhighptmuons         , "nhighptmuons/i");

  tree_->Branch("mu1pid"               , &mu1pid               , "mu1pid/I");
  tree_->Branch("mu1pt"                , &mu1pt                , "mu1pt/F");
  tree_->Branch("mu1eta"               , &mu1eta               , "mu1eta/F");
  tree_->Branch("mu1phi"               , &mu1phi               , "mu1phi/F");
  if(not isTriggerTree and not isPhotonPurity and not isQCDTree){
    tree_->Branch("mu1pfpt"              , &mu1pfpt              , "mu1pfpt/F");
    tree_->Branch("mu1pfeta"             , &mu1pfeta             , "mu1pfeta/F");
    tree_->Branch("mu1pfphi"             , &mu1pfphi             , "mu1pfphi/F");
  }
  tree_->Branch("mu1id"                , &mu1id                , "mu1id/I");
  tree_->Branch("mu1idm"               , &mu1idm               , "mu1idm/I");
  tree_->Branch("mu1idt"               , &mu1idt               , "mu1idt/I");
  tree_->Branch("mu1iso"               , &mu1iso               , "mu1iso/F");

  tree_->Branch("mu2pid"               , &mu2pid               , "mu2pid/I");
  tree_->Branch("mu2pt"                , &mu2pt                , "mu2pt/F");
  tree_->Branch("mu2eta"               , &mu2eta               , "mu2eta/F");
  tree_->Branch("mu2phi"               , &mu2phi               , "mu2phi/F");
  if(not isTriggerTree and not isPhotonPurity and not isQCDTree){
    tree_->Branch("mu2pfpt"              , &mu2pfpt              , "mu2pfpt/F");
    tree_->Branch("mu2pfeta"             , &mu2pfeta             , "mu2pfeta/F");
    tree_->Branch("mu2pfphi"             , &mu2pfphi             , "mu2pfphi/F");
  }
  tree_->Branch("mu2id"                , &mu2id                , "mu2id/I");
  tree_->Branch("mu2idm"               , &mu2idm               , "mu2idm/I");
  tree_->Branch("mu2idt"               , &mu2idt               , "mu2idt/I");
  tree_->Branch("mu2iso"               , &mu2iso               , "mu2iso/F");

  if(isReMiniAOD){
    tree_->Branch("nmuonsfake"                , &nmuonsfake                , "nmuonsfake/I");
    tree_->Branch("fakemupt","std::vector<float>",&fakemupt);
    tree_->Branch("fakemueta","std::vector<float>",&fakemueta);
    tree_->Branch("fakemuphi","std::vector<float>",&fakemuphi);
  }  
}

float MuonTreeFiller::computeMuonIso(const reco::Muon& mu) {
  
  float isoval = mu.pfIsolationR04().sumNeutralHadronEt;
  isoval += mu.pfIsolationR04().sumPhotonEt;
  isoval -= 0.5*mu.pfIsolationR04().sumPUPt;
  if (isoval < 0.) isoval = 0.;
  isoval += mu.pfIsolationR04().sumChargedHadronPt;
  isoval /= mu.pt();            

  return isoval;
}
