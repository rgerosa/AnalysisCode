#include "AnalysisCode/MonoXAnalysis/interface/RecoilTreeFiller.h"

RecoilTreeFiller::RecoilTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree, const bool & isPuppi):
  t1metTag(iConfig.getParameter<edm::InputTag>("t1met")),
  t1mumetTag(iConfig.getParameter<edm::InputTag>("t1mumet")),
  t1elmetTag(iConfig.getParameter<edm::InputTag>("t1elmet")),
  t1phmetTag(iConfig.getParameter<edm::InputTag>("t1phmet")),
  t1taumetTag(iConfig.getParameter<edm::InputTag>("t1taumet")),
  addMETBreakDown(iConfig.existsAs<bool>("addMETBreakDown") ? iConfig.getParameter<bool>("addMETBreakDown") : false),
  addMETSystematics(iConfig.existsAs<bool>("addMETSystematics") ? iConfig.getParameter<bool>("addMETSystematics") : false),  
  isMC            (iConfig.existsAs<bool>("isMC") ? iConfig.getParameter<bool>("isMC") : false),
  addBadMuonClean (iConfig.existsAs<bool>("addBadMuonClean") ? iConfig.getParameter<bool>("addBadMuonClean") : false),
  isReMiniAOD     (iConfig.existsAs<bool>("isReMiniAOD") ? iConfig.getParameter<bool>("isReMiniAOD") : false),
  isQCDTree       (iConfig.existsAs<bool>("isQCDTree") ? iConfig.getParameter<bool>("isQCDTree") : false),
  isPhotonPurity  (iConfig.existsAs<bool>("isPhotonPurity") ? iConfig.getParameter<bool>("isPhotonPurity") : false),
  isTriggerTree   (iConfig.existsAs<bool>("isTriggerTree") ? iConfig.getParameter<bool>("isTriggerTree") : false){


  isPuppi_ = isPuppi;
  if(isPuppi_){    
    addMETSystematics = (iConfig.existsAs<bool>("addPuppiMETSystematics") ? iConfig.getParameter<bool>("addPuppiMETSystematics") : false);
    t1metTag   = iConfig.getParameter<edm::InputTag>("puppit1met");
    t1mumetTag = iConfig.getParameter<edm::InputTag>("puppit1mumet");
    t1elmetTag = iConfig.getParameter<edm::InputTag>("puppit1elmet");
    t1phmetTag = iConfig.getParameter<edm::InputTag>("puppit1phmet");
    t1taumetTag = iConfig.getParameter<edm::InputTag>("puppit1taumet");
  }

  t1metToken    = iC.consumes<edm::View<pat::MET> > (t1metTag);
  t1mumetToken  = iC.consumes<edm::View<pat::MET> > (t1mumetTag);
  t1elmetToken  = iC.consumes<edm::View<pat::MET> > (t1elmetTag);
  t1phmetToken  = iC.consumes<edm::View<pat::MET> > (t1phmetTag);
  t1taumetToken = iC.consumes<edm::View<pat::MET> > (t1taumetTag);
  
  if(isReMiniAOD and not isPuppi_){
    
    t1metEGCleanTag = iConfig.getParameter<edm::InputTag>("t1metEGClean");
    t1metMuCleanTag = iConfig.getParameter<edm::InputTag>("t1metMuClean");
    t1metOriginalTag = iConfig.getParameter<edm::InputTag>("t1metOriginal");
    t1mumetEGCleanTag = iConfig.getParameter<edm::InputTag>("t1mumetEGClean");
    t1mumetMuCleanTag = iConfig.getParameter<edm::InputTag>("t1mumetMuClean");
    t1elmetEGCleanTag = iConfig.getParameter<edm::InputTag>("t1elmetEGClean");
    t1elmetMuCleanTag = iConfig.getParameter<edm::InputTag>("t1elmetMuClean");
    t1phmetEGCleanTag = iConfig.getParameter<edm::InputTag>("t1phmetEGClean");
    t1phmetMuCleanTag = iConfig.getParameter<edm::InputTag>("t1phmetMuClean");
    t1taumetEGCleanTag = iConfig.getParameter<edm::InputTag>("t1taumetEGClean");
    t1taumetMuCleanTag = iConfig.getParameter<edm::InputTag>("t1taumetMuClean");

    t1metEGCleanToken    = iC.consumes<edm::View<pat::MET> > (t1metEGCleanTag);
    t1metMuCleanToken    = iC.consumes<edm::View<pat::MET> > (t1metMuCleanTag);
    t1metOriginalToken    = iC.consumes<edm::View<pat::MET> > (t1metOriginalTag);
    t1mumetEGCleanToken    = iC.consumes<edm::View<pat::MET> > (t1mumetEGCleanTag);
    t1mumetMuCleanToken    = iC.consumes<edm::View<pat::MET> > (t1mumetMuCleanTag);
    t1elmetEGCleanToken    = iC.consumes<edm::View<pat::MET> > (t1elmetEGCleanTag);
    t1elmetMuCleanToken    = iC.consumes<edm::View<pat::MET> > (t1elmetMuCleanTag);
    t1phmetEGCleanToken    = iC.consumes<edm::View<pat::MET> > (t1phmetEGCleanTag);
    t1phmetMuCleanToken    = iC.consumes<edm::View<pat::MET> > (t1phmetMuCleanTag);
    t1taumetEGCleanToken    = iC.consumes<edm::View<pat::MET> > (t1taumetEGCleanTag);
    t1taumetMuCleanToken    = iC.consumes<edm::View<pat::MET> > (t1taumetMuCleanTag);    
  }

  if(addBadMuonClean and isMC and not isPuppi_){
    t1metOriginalTag   = iConfig.getParameter<edm::InputTag>("t1metOriginal");
    t1metOriginalToken = iC.consumes<edm::View<pat::MET> > (t1metOriginalTag);
  }
  
  if(addMETBreakDown and not isPuppi_){
    pfMetHadronHFToken = iC.consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetHadronHF"));
    pfMetEgammaHFToken = iC.consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetEgammaHF"));
    pfMetChargedHadronToken = iC.consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetChargedHadron"));
    pfMetNeutralHadronToken = iC.consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetNeutralHadron"));
    pfMetElectronsToken = iC.consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetElectrons"));
    pfMetPhotonsToken   = iC.consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetPhotons"));
    pfMetMuonsToken     = iC.consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetMuons"));
    pfMetUnclusteredToken = iC.consumes<edm::View<pat::MET> > (iConfig.getParameter<edm::InputTag>("pfMetUnclustered"));
  }
 
  
  tree_ = tree;
  this->DeclareAndSetBranches();
}

/////
void RecoilTreeFiller::initBranches(){
  
  genmet     = -99.;      genmetphi = -99.;
  t1pfmet    = -99. ; t1pfmetphi = -99.; pfmet      = -99. ; pfmetphi   = -99.; calomet    = -99. ; calometphi = -99.;
  t1mumet    = -99;   t1mumetphi = -99;  mumet      = -99;   mumetphi   = -99;  t1elmet    = -99;   t1elmetphi = -99;
  elmet      = -99;   elmetphi   = -99;  t1phmet    = -99;   t1phmetphi = -99;  phmet      = -99;   phmetphi   = -99;
  t1taumet   = -99;  t1taumetphi = -99; taumet      = -99;  taumetphi   = -99;

  t1pfmetEGClean    = -99. ; t1pfmetphiEGClean = -99. ; t1pfmetOriginal    = -99. ; t1pfmetphiOriginal = -99. ;
  t1pfmetMuClean    = -99. ; t1pfmetphiMuClean = -99. ; t1mumetEGClean    = -99. ; t1mumetphiEGClean = -99. ;
  t1mumetMuClean    = -99. ; t1mumetphiMuClean = -99. ; t1elmetEGClean    = -99. ; t1elmetphiEGClean = -99. ; 
  t1elmetMuClean    = -99. ; t1elmetphiMuClean = -99. ; t1phmetEGClean    = -99. ; t1phmetphiEGClean = -99. ;
  t1phmetMuClean    = -99. ; t1phmetphiMuClean = -99. ; t1taumetEGClean    = -99. ; t1taumetphiEGClean = -99. ;
  t1taumetMuClean    = -99. ; t1taumetphiMuClean = -99. ;

  t1pfmetMuEnUp  = -99.; t1pfmetMuEnDown  = -99.; t1pfmetElEnUp   = -99.; t1pfmetElEnDown   = -99.;
  t1pfmetPhoEnUp = -99.; t1pfmetPhoEnDown = -99.; t1pfmetTauEnUp  = -99.; t1pfmetTauEnDown  = -99.;
  t1pfmetJetEnUp = -99.; t1pfmetJetEnDown = -99.; t1pfmetJetResUp = -99.; t1pfmetJetResDown = -99.;
  t1pfmetUncEnUp = -99.; t1pfmetUncEnDown = -99.; t1pfmetJetSmear = -99.; t1pfmetXY = -99.;     

  t1pfmetMuEnUpPhi  = -99.; t1pfmetMuEnDownPhi  = -99.; t1pfmetElEnUpPhi   = -99.; t1pfmetElEnDownPhi   = -99.;
  t1pfmetPhoEnUpPhi = -99.; t1pfmetPhoEnDownPhi = -99.; t1pfmetTauEnUpPhi  = -99.; t1pfmetTauEnDownPhi  = -99.;
  t1pfmetJetEnUpPhi = -99.; t1pfmetJetEnDownPhi = -99.; t1pfmetJetResUpPhi = -99.; t1pfmetJetResDownPhi = -99.;
  t1pfmetUncEnUpPhi = -99.; t1pfmetUncEnDownPhi = -99.; t1pfmetJetSmearPhi = -99.; t1pfmetXYPhi = -99.;     

  // MET break down
  pfmethadronHF = -99. ; pfmethadronHFphi = -99. ; 
  pfmetegammaHF = -99. ; pfmetegammaHFphi = -99. ; 
  pfmetchargedhadron = -99. ; pfmetchargedhadronphi = -99.;
  pfmetneutralhadron = -99. ; pfmetneutralhadronphi = -99. ; 
  pfmetelectrons     = -99. ; pfmetelectronsphi   = -99. ; 
  pfmetmuons         = -99. ; pfmetmuonsphi       = -99. ; 
  pfmetphotons       = -99. ; pfmetphotonsphi     = -99. ; 
  pfmetunclustered   = -99. ; pfmetunclusteredphi = -99.;

  // puppi met info    
  puppit1pfmet = -99.; puppit1pfmetphi = -99.;
  puppipfmet   = -99.; puppipfmetphi   = -99.;
  puppit1mumet = -99.; puppit1mumetphi = -99.;
  puppimumet   = -99.; puppimumetphi   = -99.;
  puppit1elmet = -99.; puppit1elmetphi = -99.;
  puppielmet   = -99.; puppielmetphi   = -99.;
  puppit1phmet = -99.; puppit1phmetphi = -99.;
  puppiphmet   = -99.; puppiphmetphi   = -99.;
    
  puppit1pfmetMuEnUp  = -99.; puppit1pfmetMuEnDown  = -99.; puppit1pfmetElEnUp   = -99.; puppit1pfmetElEnDown   = -99.;
  puppit1pfmetPhoEnUp = -99.; puppit1pfmetPhoEnDown = -99.; puppit1pfmetTauEnUp  = -99.; puppit1pfmetTauEnDown  = -99.;
  puppit1pfmetJetEnUp = -99.; puppit1pfmetJetEnDown = -99.; puppit1pfmetJetResUp = -99.; puppit1pfmetJetResDown = -99.;
  puppit1pfmetUncEnUp = -99.; puppit1pfmetUncEnDown = -99.;

  puppit1pfmetMuEnUpPhi  = -99.; puppit1pfmetMuEnDownPhi  = -99.; puppit1pfmetElEnUpPhi   = -99.; puppit1pfmetElEnDownPhi   = -99.;
  puppit1pfmetPhoEnUpPhi = -99.; puppit1pfmetPhoEnDownPhi = -99.; puppit1pfmetTauEnUpPhi  = -99.; puppit1pfmetTauEnDownPhi  = -99.;
  puppit1pfmetJetEnUpPhi = -99.; puppit1pfmetJetEnDownPhi = -99.; puppit1pfmetJetResUpPhi = -99.; puppit1pfmetJetResDownPhi = -99.;
  puppit1pfmetUncEnUpPhi = -99.; puppit1pfmetUncEnDownPhi = -99.;
}

/////
bool RecoilTreeFiller::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace pat;
  
  this->initBranches();
  
  Handle<View<pat::MET> > t1metH;
  iEvent.getByToken(t1metToken, t1metH);
  
  Handle<View<pat::MET> > t1metEGCleanH;
  Handle<View<pat::MET> > t1metMuCleanH;
  Handle<View<pat::MET> > t1metOriginalH;
  if(isReMiniAOD){
    iEvent.getByToken(t1metEGCleanToken, t1metEGCleanH);
    iEvent.getByToken(t1metMuCleanToken, t1metMuCleanH);
    iEvent.getByToken(t1metOriginalToken, t1metOriginalH);
  }
  if(addBadMuonClean and isMC){
    iEvent.getByToken(t1metOriginalToken, t1metOriginalH);
  }
  
  
  Handle<View<pat::MET> > t1mumetH;
  iEvent.getByToken(t1mumetToken, t1mumetH);

  Handle<View<pat::MET> > t1mumetEGCleanH;
  Handle<View<pat::MET> > t1mumetMuCleanH;
  if(isReMiniAOD){
    iEvent.getByToken(t1mumetEGCleanToken, t1mumetEGCleanH);
    iEvent.getByToken(t1mumetMuCleanToken, t1mumetMuCleanH);
  }

  Handle<View<pat::MET> > t1elmetH;
  iEvent.getByToken(t1elmetToken, t1elmetH);
 
  Handle<View<pat::MET> > t1elmetEGCleanH;
  Handle<View<pat::MET> > t1elmetMuCleanH;
  if(isReMiniAOD){
    iEvent.getByToken(t1elmetEGCleanToken, t1elmetEGCleanH);
    iEvent.getByToken(t1elmetMuCleanToken, t1elmetMuCleanH);
  }


  Handle<View<pat::MET> > t1phmetH;
  iEvent.getByToken(t1phmetToken, t1phmetH);
  Handle<View<pat::MET> > t1phmetEGCleanH;
  Handle<View<pat::MET> > t1phmetMuCleanH;
  if(isReMiniAOD){
    iEvent.getByToken(t1phmetEGCleanToken, t1phmetEGCleanH);
    iEvent.getByToken(t1phmetMuCleanToken, t1phmetMuCleanH);
  }
  

  Handle<View<pat::MET> > t1taumetH;
  iEvent.getByToken(t1taumetToken, t1taumetH);
  Handle<View<pat::MET> > t1taumetEGCleanH;
  Handle<View<pat::MET> > t1taumetMuCleanH;
  if(isReMiniAOD){
    iEvent.getByToken(t1taumetEGCleanToken, t1taumetEGCleanH);
    iEvent.getByToken(t1taumetMuCleanToken, t1taumetMuCleanH);
  }

  // MET breakdown 
  Handle<View<pat::MET> > pfMetHadronHFH;
  Handle<View<pat::MET> > pfMetEgammaHFH;
  Handle<View<pat::MET> > pfMetChargedHadronH;
  Handle<View<pat::MET> > pfMetNeutralHadronH;
  Handle<View<pat::MET> > pfMetElectronsH;
  Handle<View<pat::MET> > pfMetPhotonsH;
  Handle<View<pat::MET> > pfMetMuonsH;
  Handle<View<pat::MET> > pfMetUnclusteredH;

  if(addMETBreakDown and not isPuppi_){
    iEvent.getByToken(pfMetHadronHFToken,pfMetHadronHFH);
    iEvent.getByToken(pfMetEgammaHFToken,pfMetEgammaHFH);
    iEvent.getByToken(pfMetChargedHadronToken,pfMetChargedHadronH);
    iEvent.getByToken(pfMetNeutralHadronToken,pfMetNeutralHadronH);
    iEvent.getByToken(pfMetElectronsToken,pfMetElectronsH);
    iEvent.getByToken(pfMetPhotonsToken,pfMetPhotonsH);
    iEvent.getByToken(pfMetMuonsToken,pfMetMuonsH);
    iEvent.getByToken(pfMetUnclusteredToken,pfMetUnclusteredH);
  }

  // Fill info
  if(t1metH.isValid() and not isPuppi_){      
    // dump gen met info
    if(t1metH->front().genMET()){
      genmet    = t1metH->front().genMET()->pt();
      genmetphi = t1metH->front().genMET()->phi();
    }
      
    t1pfmet    = t1metH->front().corPt();    
    t1pfmetphi = t1metH->front().corPhi();
    pfmet      = t1metH->front().uncorPt();
    pfmetphi   = t1metH->front().uncorPhi();
    if(not isReMiniAOD or (isMC and addBadMuonClean)){
      calomet    = t1metH->front().caloMETPt();
      calometphi = t1metH->front().caloMETPhi(); 
    }
    else if(isReMiniAOD and t1metOriginalH.isValid()){
      t1pfmetOriginal = t1metOriginalH->front().corPt();
      t1pfmetphiOriginal = t1metOriginalH->front().corPhi();
    }
  }

  if(t1mumetH.isValid() and not isPuppi_){
    t1mumet    = t1mumetH->front().corPt();
    t1mumetphi = t1mumetH->front().corPhi();
    mumet      = t1mumetH->front().uncorPt();
    mumetphi   = t1mumetH->front().uncorPhi();
  }

  if(t1elmetH.isValid() and not isPuppi_){
    t1elmet    = t1elmetH->front().corPt();
    t1elmetphi = t1elmetH->front().corPhi();
    elmet      = t1elmetH->front().uncorPt();
    elmetphi   = t1elmetH->front().uncorPhi();
  }

  if(t1phmetH.isValid() and not isPuppi_){
    t1phmet    = t1phmetH->front().corPt();
    t1phmetphi = t1phmetH->front().corPhi();
    phmet      = t1phmetH->front().uncorPt();
    phmetphi   = t1phmetH->front().uncorPhi();
  }

  if(t1taumetH.isValid() and not isPuppi_){
    t1taumet    = t1taumetH->front().corPt();
    t1taumetphi = t1taumetH->front().corPhi();
    taumet      = t1taumetH->front().uncorPt();
    taumetphi   = t1taumetH->front().uncorPhi();
  }

  // only for re-miniAOD
  if(isReMiniAOD and not isPuppi_){ 
    if(t1metEGCleanH.isValid()){
      t1pfmetEGClean    = t1metEGCleanH->front().corPt();
      t1pfmetphiEGClean = t1metEGCleanH->front().corPhi();
    }
      
    if(t1metMuCleanH.isValid()){
      t1pfmetMuClean    = t1metMuCleanH->front().corPt();
      t1pfmetphiMuClean = t1metMuCleanH->front().corPhi();
      calomet    = t1metMuCleanH->front().caloMETPt();
      calometphi = t1metMuCleanH->front().caloMETPhi(); 
    }
      
    if(t1mumetEGCleanH.isValid()){
      t1mumetEGClean    = t1mumetEGCleanH->front().corPt();
      t1mumetphiEGClean = t1mumetEGCleanH->front().corPhi();
    }
      
    if(t1mumetMuCleanH.isValid()){
      t1mumetMuClean    = t1mumetMuCleanH->front().corPt();
      t1mumetphiMuClean = t1mumetMuCleanH->front().corPhi();
    }
      
      
    if(t1elmetEGCleanH.isValid()){
      t1elmetEGClean    = t1elmetEGCleanH->front().corPt();
      t1elmetphiEGClean = t1elmetEGCleanH->front().corPhi();
    }
      
    if(t1elmetMuCleanH.isValid()){
      t1elmetMuClean    = t1elmetMuCleanH->front().corPt();
      t1elmetphiMuClean = t1elmetMuCleanH->front().corPhi();
    }

    if(t1phmetEGCleanH.isValid()){
      t1phmetEGClean    = t1phmetEGCleanH->front().corPt();
      t1phmetphiEGClean = t1phmetEGCleanH->front().corPhi();
    }
      
    if(t1phmetMuCleanH.isValid()){
      t1phmetMuClean    = t1phmetMuCleanH->front().corPt();
      t1phmetphiMuClean = t1phmetMuCleanH->front().corPhi();
    }

    if(t1taumetEGCleanH.isValid()){
      t1taumetEGClean    = t1taumetEGCleanH->front().corPt();
      t1taumetphiEGClean = t1taumetEGCleanH->front().corPhi();
    }
      
    if(t1phmetMuCleanH.isValid()){
      t1taumetMuClean    = t1taumetMuCleanH->front().corPt();
      t1taumetphiMuClean = t1taumetMuCleanH->front().corPhi();
    }      
  }
  
  if(addMETSystematics  and not isTriggerTree and not isPhotonPurity){          
    if(t1metH.isValid() and not isPuppi_){
      t1pfmetMuEnUp     = t1metH->front().shiftedPt(pat::MET::METUncertainty::MuonEnUp,   pat::MET::METCorrectionLevel::Type1);
      t1pfmetMuEnDown   = t1metH->front().shiftedPt(pat::MET::METUncertainty::MuonEnDown, pat::MET::METCorrectionLevel::Type1);
      t1pfmetElEnUp     = t1metH->front().shiftedPt(pat::MET::METUncertainty::ElectronEnUp,   pat::MET::METCorrectionLevel::Type1);
      t1pfmetElEnDown   = t1metH->front().shiftedPt(pat::MET::METUncertainty::ElectronEnDown, pat::MET::METCorrectionLevel::Type1);
      t1pfmetPhoEnUp    = t1metH->front().shiftedPt(pat::MET::METUncertainty::PhotonEnUp,   pat::MET::METCorrectionLevel::Type1);
      t1pfmetPhoEnDown  = t1metH->front().shiftedPt(pat::MET::METUncertainty::PhotonEnDown, pat::MET::METCorrectionLevel::Type1);
      t1pfmetTauEnUp    = t1metH->front().shiftedPt(pat::MET::METUncertainty::TauEnUp,   pat::MET::METCorrectionLevel::Type1);
      t1pfmetTauEnDown  = t1metH->front().shiftedPt(pat::MET::METUncertainty::TauEnDown, pat::MET::METCorrectionLevel::Type1);
      t1pfmetJetEnUp    = t1metH->front().shiftedPt(pat::MET::METUncertainty::JetEnUp,   pat::MET::METCorrectionLevel::Type1);
      t1pfmetJetEnDown  = t1metH->front().shiftedPt(pat::MET::METUncertainty::JetEnDown, pat::MET::METCorrectionLevel::Type1);
      t1pfmetJetResUp   = t1metH->front().shiftedPt(pat::MET::METUncertainty::JetResUp,   pat::MET::METCorrectionLevel::Type1Smear);
      t1pfmetJetResDown = t1metH->front().shiftedPt(pat::MET::METUncertainty::JetResDown, pat::MET::METCorrectionLevel::Type1Smear);
      t1pfmetUncEnUp    = t1metH->front().shiftedPt(pat::MET::METUncertainty::UnclusteredEnUp,   pat::MET::METCorrectionLevel::Type1);
      t1pfmetUncEnDown  = t1metH->front().shiftedPt(pat::MET::METUncertainty::UnclusteredEnDown, pat::MET::METCorrectionLevel::Type1);
      t1pfmetJetSmear   = t1metH->front().shiftedPt(pat::MET::NoShift, pat::MET::METCorrectionLevel::Type1Smear);
      t1pfmetXY         = t1metH->front().shiftedPt(pat::MET::NoShift, pat::MET::METCorrectionLevel::Type1XY);
      t1pfmetMuEnUpPhi     = t1metH->front().shiftedPhi(pat::MET::METUncertainty::MuonEnUp,   pat::MET::METCorrectionLevel::Type1);
      t1pfmetMuEnDownPhi   = t1metH->front().shiftedPhi(pat::MET::METUncertainty::MuonEnDown, pat::MET::METCorrectionLevel::Type1);
      t1pfmetElEnUpPhi     = t1metH->front().shiftedPhi(pat::MET::METUncertainty::ElectronEnUp,   pat::MET::METCorrectionLevel::Type1);
      t1pfmetElEnDownPhi   = t1metH->front().shiftedPhi(pat::MET::METUncertainty::ElectronEnDown, pat::MET::METCorrectionLevel::Type1);
      t1pfmetPhoEnUpPhi    = t1metH->front().shiftedPhi(pat::MET::METUncertainty::PhotonEnUp,   pat::MET::METCorrectionLevel::Type1);
      t1pfmetPhoEnDownPhi  = t1metH->front().shiftedPhi(pat::MET::METUncertainty::PhotonEnDown, pat::MET::METCorrectionLevel::Type1);
      t1pfmetTauEnUpPhi    = t1metH->front().shiftedPhi(pat::MET::METUncertainty::TauEnUp,   pat::MET::METCorrectionLevel::Type1);
      t1pfmetTauEnDownPhi  = t1metH->front().shiftedPhi(pat::MET::METUncertainty::TauEnDown, pat::MET::METCorrectionLevel::Type1);
      t1pfmetJetEnUpPhi    = t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetEnUp,   pat::MET::METCorrectionLevel::Type1);
      t1pfmetJetEnDownPhi  = t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetEnDown, pat::MET::METCorrectionLevel::Type1);
      t1pfmetJetResUpPhi   = t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetResUp,   pat::MET::METCorrectionLevel::Type1Smear);
      t1pfmetJetResDownPhi = t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetResDown, pat::MET::METCorrectionLevel::Type1Smear);
      t1pfmetUncEnUpPhi    = t1metH->front().shiftedPhi(pat::MET::METUncertainty::UnclusteredEnUp,   pat::MET::METCorrectionLevel::Type1);
      t1pfmetUncEnDownPhi  = t1metH->front().shiftedPhi(pat::MET::METUncertainty::UnclusteredEnDown, pat::MET::METCorrectionLevel::Type1);
      t1pfmetJetSmearPhi   = t1metH->front().shiftedPhi(pat::MET::NoShift, pat::MET::METCorrectionLevel::Type1Smear);
      t1pfmetXYPhi         = t1metH->front().shiftedPhi(pat::MET::NoShift, pat::MET::METCorrectionLevel::Type1XY);
    }
  }

  ////
  if(addMETBreakDown and not isTriggerTree and not isPhotonPurity and not isPuppi_){
    pfmethadronHF    = pfMetHadronHFH->front().pt();
    pfmethadronHFphi = pfMetHadronHFH->front().phi();
    pfmetegammaHF    = pfMetEgammaHFH->front().pt();
    pfmetegammaHFphi = pfMetEgammaHFH->front().phi();
    pfmetchargedhadron    = pfMetChargedHadronH->front().pt();
    pfmetchargedhadronphi = pfMetChargedHadronH->front().phi();
    pfmetneutralhadron    = pfMetNeutralHadronH->front().pt();
    pfmetneutralhadronphi = pfMetNeutralHadronH->front().phi();
    pfmetelectrons    = pfMetElectronsH->front().pt();
    pfmetelectronsphi = pfMetElectronsH->front().phi();
    pfmetmuons        = pfMetMuonsH->front().pt();
    pfmetmuonsphi     = pfMetMuonsH->front().phi();
    pfmetphotons      = pfMetPhotonsH->front().pt();
    pfmetphotonsphi   = pfMetPhotonsH->front().phi();
    pfmetunclustered  = pfMetUnclusteredH->front().pt();
    pfmetunclusteredphi = pfMetUnclusteredH->front().phi();
  }
    
  // Puppi sector
  if(not isTriggerTree and not isPhotonPurity and isPuppi_){    
    if(t1metH.isValid()){
      puppit1pfmet        = t1metH->front().corPt();
      puppit1pfmetphi     = t1metH->front().corPhi();
      puppipfmet          = t1metH->front().uncorPt();
      puppipfmetphi       = t1metH->front().uncorPhi();
    }

    if(t1mumetH.isValid()){
      puppit1mumet        = t1mumetH->front().corPt();
      puppit1mumetphi     = t1mumetH->front().corPhi();
      puppimumet          = t1mumetH->front().uncorPt();
      puppimumetphi       = t1mumetH->front().uncorPhi();
    }
      
    if(t1elmetH.isValid()){
      puppit1elmet        = t1elmetH->front().corPt();
      puppit1elmetphi     = t1elmetH->front().corPhi();
      puppielmet          = t1elmetH->front().uncorPt();
      puppielmetphi       = t1elmetH->front().uncorPhi();
    }

    if(t1phmetH.isValid()){
      puppit1phmet        = t1phmetH->front().corPt();
      puppit1phmetphi     = t1phmetH->front().corPhi();
      puppiphmet          = t1phmetH->front().uncorPt();
      puppiphmetphi       = t1phmetH->front().uncorPhi();
    }

    if(t1taumetH.isValid()){
      puppit1taumet        = t1taumetH->front().corPt();
      puppit1taumetphi     = t1taumetH->front().corPhi();
      puppitaumet          = t1taumetH->front().uncorPt();
      puppitaumetphi       = t1taumetH->front().uncorPhi();
    }
    
    if(addMETSystematics and not isTriggerTree and not isPhotonPurity and isPuppi_){
      puppit1pfmetMuEnUp     = t1metH->front().shiftedPt(pat::MET::METUncertainty::MuonEnUp, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetMuEnDown   = t1metH->front().shiftedPt(pat::MET::METUncertainty::MuonEnDown, pat::MET::METCorrectionLevel::Type1);       
      puppit1pfmetElEnUp     = t1metH->front().shiftedPt(pat::MET::METUncertainty::ElectronEnUp, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetElEnDown   = t1metH->front().shiftedPt(pat::MET::METUncertainty::ElectronEnDown, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetPhoEnUp    = t1metH->front().shiftedPt(pat::MET::METUncertainty::PhotonEnUp, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetPhoEnDown  = t1metH->front().shiftedPt(pat::MET::METUncertainty::PhotonEnDown, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetTauEnUp    = t1metH->front().shiftedPt(pat::MET::METUncertainty::TauEnUp, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetTauEnDown  = t1metH->front().shiftedPt(pat::MET::METUncertainty::TauEnDown, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetJetEnUp    = t1metH->front().shiftedPt(pat::MET::METUncertainty::JetEnUp, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetJetEnDown  = t1metH->front().shiftedPt(pat::MET::METUncertainty::JetEnDown, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetJetResUp   = t1metH->front().shiftedPt(pat::MET::METUncertainty::JetResUp, pat::MET::METCorrectionLevel::Type1Smear);
      puppit1pfmetJetResDown = t1metH->front().shiftedPt(pat::MET::METUncertainty::JetResDown, pat::MET::METCorrectionLevel::Type1Smear);     
      puppit1pfmetUncEnUp    = t1metH->front().shiftedPt(pat::MET::METUncertainty::UnclusteredEnUp, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetUncEnDown  = t1metH->front().shiftedPt(pat::MET::METUncertainty::UnclusteredEnDown, pat::MET::METCorrectionLevel::Type1);      
      
      puppit1pfmetMuEnUpPhi     = t1metH->front().shiftedPhi(pat::MET::METUncertainty::MuonEnUp, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetMuEnDownPhi   = t1metH->front().shiftedPhi(pat::MET::METUncertainty::MuonEnDown, pat::MET::METCorrectionLevel::Type1);       
      puppit1pfmetElEnUpPhi     = t1metH->front().shiftedPhi(pat::MET::METUncertainty::ElectronEnUp, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetElEnDownPhi   = t1metH->front().shiftedPhi(pat::MET::METUncertainty::ElectronEnDown, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetPhoEnUpPhi    = t1metH->front().shiftedPhi(pat::MET::METUncertainty::PhotonEnUp, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetPhoEnDownPhi  = t1metH->front().shiftedPhi(pat::MET::METUncertainty::PhotonEnDown, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetTauEnUpPhi    = t1metH->front().shiftedPhi(pat::MET::METUncertainty::TauEnUp, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetTauEnDownPhi  = t1metH->front().shiftedPhi(pat::MET::METUncertainty::TauEnDown, pat::MET::METCorrectionLevel::Type1);  
      puppit1pfmetJetEnUpPhi    = t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetEnUp, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetJetEnDownPhi  = t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetEnDown, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetJetResUpPhi   = t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetResUp, pat::MET::METCorrectionLevel::Type1Smear);     
      puppit1pfmetJetResDownPhi = t1metH->front().shiftedPhi(pat::MET::METUncertainty::JetResDown, pat::MET::METCorrectionLevel::Type1Smear);     
      puppit1pfmetUncEnUpPhi    = t1metH->front().shiftedPhi(pat::MET::METUncertainty::UnclusteredEnUp, pat::MET::METCorrectionLevel::Type1);
      puppit1pfmetUncEnDownPhi  = t1metH->front().shiftedPhi(pat::MET::METUncertainty::UnclusteredEnDown, pat::MET::METCorrectionLevel::Type1);      
    }
  }
  
  return true;  
}


////////
void RecoilTreeFiller::DeclareAndSetBranches(){

  if(not isPuppi_){
    
    tree_->Branch("pfmet"                , &pfmet                , "pfmet/F");
    tree_->Branch("pfmetphi"             , &pfmetphi             , "pfmetphi/F");
    tree_->Branch("calomet"              , &calomet              , "calomet/F");   //ND
    tree_->Branch("calometphi"           , &calometphi           , "calometphi/F");//ND

    tree_->Branch("genmet"                , &genmet                , "genmet/F");
    tree_->Branch("genmetphi"             , &genmetphi             , "genmetphi/F");

    tree_->Branch("t1pfmet"              , &t1pfmet              , "t1pfmet/F");
    tree_->Branch("t1pfmetphi"           , &t1pfmetphi           , "t1pfmetphi/F");
    tree_->Branch("t1mumet"              , &t1mumet              , "t1mumet/F");
    tree_->Branch("t1mumetphi"           , &t1mumetphi           , "t1mumetphi/F");
    tree_->Branch("t1elmet"              , &t1elmet              , "t1elmet/F");
    tree_->Branch("t1elmetphi"           , &t1elmetphi           , "t1elmetphi/F");
    tree_->Branch("t1phmet"              , &t1phmet              , "t1phmet/F");
    tree_->Branch("t1phmetphi"           , &t1phmetphi           , "t1phmetphi/F");
    tree_->Branch("t1taumet"              , &t1taumet              , "t1taumet/F");
    tree_->Branch("t1taumetphi"           , &t1taumetphi           , "t1taumetphi/F");

  
    if(isReMiniAOD){
    
      tree_->Branch("t1pfmetEGClean"              , &t1pfmetEGClean              , "t1pfmetEGClean/F");
      tree_->Branch("t1pfmetphiEGClean"           , &t1pfmetphiEGClean           , "t1pfmetphiEGClean/F");
      tree_->Branch("t1mumetEGClean"              , &t1mumetEGClean              , "t1mumetEGClean/F");
      tree_->Branch("t1mumetphiEGClean"           , &t1mumetphiEGClean           , "t1mumetphiEGClean/F");
      tree_->Branch("t1elmetEGClean"              , &t1elmetEGClean              , "t1elmetEGClean/F");
      tree_->Branch("t1elmetphiEGClean"           , &t1elmetphiEGClean           , "t1elmetphiEGClean/F");
      tree_->Branch("t1phmetEGClean"              , &t1phmetEGClean              , "t1phmetEGClean/F");
      tree_->Branch("t1phmetphiEGClean"           , &t1phmetphiEGClean           , "t1phmetphiEGClean/F");
      tree_->Branch("t1taumetEGClean"             , &t1taumetEGClean             , "t1taumetEGClean/F");
      tree_->Branch("t1taumetphiEGClean"          , &t1taumetphiEGClean          , "t1taumetphiEGClean/F");
      
      tree_->Branch("t1pfmetMuClean"              , &t1pfmetMuClean              , "t1pfmetMuClean/F");
      tree_->Branch("t1pfmetphiMuClean"           , &t1pfmetphiMuClean           , "t1pfmetphiMuClean/F");
      tree_->Branch("t1mumetMuClean"              , &t1mumetMuClean              , "t1mumetMuClean/F");
      tree_->Branch("t1mumetphiMuClean"           , &t1mumetphiMuClean           , "t1mumetphiMuClean/F");
      tree_->Branch("t1elmetMuClean"              , &t1elmetMuClean              , "t1elmetMuClean/F");
      tree_->Branch("t1elmetphiMuClean"           , &t1elmetphiMuClean           , "t1elmetphiMuClean/F");
      tree_->Branch("t1phmetMuClean"              , &t1phmetMuClean              , "t1phmetMuClean/F");
      tree_->Branch("t1phmetphiMuClean"           , &t1phmetphiMuClean           , "t1phmetphiMuClean/F");
      tree_->Branch("t1taumetMuClean"              , &t1taumetMuClean              , "t1taumetMuClean/F");
      tree_->Branch("t1taumetphiMuClean"           , &t1taumetphiMuClean           , "t1taumetphiMuClean/F");
      
      tree_->Branch("t1pfmetOriginal"              , &t1pfmetOriginal              , "t1pfmetOriginal/F");
      tree_->Branch("t1pfmetphiOriginal"           , &t1pfmetphiOriginal           , "t1pfmetphiOriginal/F");
      
    }
    
    if(addMETBreakDown and not isTriggerTree and not isPhotonPurity and not isQCDTree and not isPuppi_){    
      tree_->Branch("pfmethadronHF",&pfmethadronHF,"pfmethadronHF/F");
      tree_->Branch("pfmethadronHFphi",&pfmethadronHFphi,"pfmethadronHFphi/F");
      tree_->Branch("pfmetegammaHF",&pfmetegammaHF,"pfmetegammaHF/F");
      tree_->Branch("pfmetegammaHFphi",&pfmetegammaHFphi,"pfmetegammaHFphi/F");
      tree_->Branch("pfmetchargedhadron",&pfmetchargedhadron,"pfmetchargedhadron/F");
      tree_->Branch("pfmetchargedhadronphi",&pfmetchargedhadronphi,"pfmetchargedhadronphi/F");
      tree_->Branch("pfmetneutralhadron",&pfmetneutralhadron,"pfmetneutralhadron/F");
      tree_->Branch("pfmetneutralhadronphi",&pfmetneutralhadronphi,"pfmetneutralhadronphi/F");
      tree_->Branch("pfmetelectrons",&pfmetelectrons,"pfmetelectrons/F");
      tree_->Branch("pfmetelectronsphi",&pfmetelectronsphi,"pfmetelectronsphi/F");
      tree_->Branch("pfmetmuons",&pfmetmuons,"pfmetmuons/F");
      tree_->Branch("pfmetmuonsphi",&pfmetmuonsphi,"pfmetmuonsphi/F");
      tree_->Branch("pfmetphotons",&pfmetphotons,"pfmetphotons/F");
      tree_->Branch("pfmetphotonsphi",&pfmetphotonsphi,"pfmetphotonsphi/F");
      tree_->Branch("pfmetunclustered",&pfmetunclustered,"pfmetunclustered/F");
      tree_->Branch("pfmetunclusteredphi",&pfmetunclusteredphi,"pfmetunclusteredphi/F");
    }
    
    if(addMETSystematics and not isTriggerTree and not isPhotonPurity and not isQCDTree){
      tree_->Branch("t1pfmetMuEnUp"       , &t1pfmetMuEnUp       , "t1pfmetMuEnUp/F");
      tree_->Branch("t1pfmetMuEnDown"     , &t1pfmetMuEnDown     , "t1pfmetMuEnDown/F");
      tree_->Branch("t1pfmetElEnUp"       , &t1pfmetElEnUp       , "t1pfmetElEnUp/F");
      tree_->Branch("t1pfmetElEnDown"     , &t1pfmetElEnDown     , "t1pfmetElEnDown/F");
      tree_->Branch("t1pfmetPhoEnUp"      , &t1pfmetPhoEnUp      , "t1pfmetPhoEnUp/F");
      tree_->Branch("t1pfmetPhoEnDown"    , &t1pfmetPhoEnDown    , "t1pfmetPhoEnDown/F");
      tree_->Branch("t1pfmetTauEnUp"      , &t1pfmetTauEnUp      , "t1pfmetTauEnUp/F");
      tree_->Branch("t1pfmetTauEnDown"    , &t1pfmetTauEnDown    , "t1pfmetTauEnDown/F");
      tree_->Branch("t1pfmetJetEnUp"      , &t1pfmetJetEnUp      , "t1pfmetJetEnUp/F");
      tree_->Branch("t1pfmetJetEnDown"    , &t1pfmetJetEnDown    , "t1pfmetJetEnDown/F");
      tree_->Branch("t1pfmetJetResUp"     , &t1pfmetJetResUp     , "t1pfmetJetResUp/F");
      tree_->Branch("t1pfmetJetResDown"   , &t1pfmetJetResDown   , "t1pfmetJetResDown/F");
      tree_->Branch("t1pfmetUncEnUp"      , &t1pfmetUncEnUp      , "t1pfmetUncEnUp/F");
      tree_->Branch("t1pfmetUncEnDown"    , &t1pfmetUncEnDown    , "t1pfmetUncEnDown/F");
      tree_->Branch("t1pfmetJetSmear"        , &t1pfmetJetSmear        , "t1pfmetJetSmear/F");
      tree_->Branch("t1pfmetXY"           , &t1pfmetXY           , "t1pfmetXY/F");
      
      tree_->Branch("t1pfmetMuEnUpPhi"       , &t1pfmetMuEnUpPhi       , "t1pfmetMuEnUpPhi/F");
      tree_->Branch("t1pfmetMuEnDownPhi"     , &t1pfmetMuEnDownPhi     , "t1pfmetMuEnDownPhi/F");
      tree_->Branch("t1pfmetElEnUpPhi"       , &t1pfmetElEnUpPhi       , "t1pfmetElEnUpPhi/F");
      tree_->Branch("t1pfmetElEnDownPhi"     , &t1pfmetElEnDownPhi     , "t1pfmetElEnDownPhi/F");
      tree_->Branch("t1pfmetPhoEnUpPhi"      , &t1pfmetPhoEnUpPhi      , "t1pfmetPhoEnUpPhi/F");
      tree_->Branch("t1pfmetPhoEnDownPhi"    , &t1pfmetPhoEnDownPhi    , "t1pfmetPhoEnDownPhi/F");
      tree_->Branch("t1pfmetTauEnUpPhi"      , &t1pfmetTauEnUpPhi      , "t1pfmetTauEnUpPhi/F");
      tree_->Branch("t1pfmetTauEnDownPhi"    , &t1pfmetTauEnDownPhi    , "t1pfmetTauEnDownPhi/F");
      tree_->Branch("t1pfmetJetEnUpPhi"      , &t1pfmetJetEnUpPhi      , "t1pfmetJetEnUpPhi/F");
      tree_->Branch("t1pfmetJetEnDownPhi"    , &t1pfmetJetEnDownPhi    , "t1pfmetJetEnDownPhi/F");
      tree_->Branch("t1pfmetJetResUpPhi"     , &t1pfmetJetResUpPhi     , "t1pfmetJetResUpPhi/F");
      tree_->Branch("t1pfmetJetResDownPhi"   , &t1pfmetJetResDownPhi   , "t1pfmetJetResDownPhi/F");
      tree_->Branch("t1pfmetUncEnUpPhi"      , &t1pfmetUncEnUpPhi      , "t1pfmetUncEnUpPhi/F");
      tree_->Branch("t1pfmetUncEnDownPhi"    , &t1pfmetUncEnDownPhi    , "t1pfmetUncEnDownPhi/F");
      tree_->Branch("t1pfmetJetSmearPhi"        , &t1pfmetJetSmearPhi        , "t1pfmetJetSmearPhi/F");
      tree_->Branch("t1pfmetXYPhi"           , &t1pfmetXYPhi           , "t1pfmetXYPhi/F");      
    }
  }
  else{

    if(not isTriggerTree and not isPhotonPurity and not isQCDTree){
      tree_->Branch("puppipfmet"                , &puppipfmet                , "puppipfmet/F");
      tree_->Branch("puppipfmetphi"             , &puppipfmetphi             , "puppipfmetphi/F");
      tree_->Branch("puppit1pfmet"              , &puppit1pfmet              , "puppit1pfmet/F");
      tree_->Branch("puppit1pfmetphi"           , &puppit1pfmetphi           , "puppit1pfmetphi/F");
      tree_->Branch("puppimumet"                , &puppimumet                , "puppimumet/F");
      tree_->Branch("puppimumetphi"             , &puppimumetphi             , "puppimumetphi/F");
      tree_->Branch("puppit1mumet"              , &puppit1mumet              , "puppit1mumet/F");
      tree_->Branch("puppit1mumetphi"           , &puppit1mumetphi           , "puppit1mumetphi/F");
      tree_->Branch("puppielmet"                , &puppielmet                , "puppielmet/F");
      tree_->Branch("puppielmetphi"             , &puppielmetphi             , "elmetphi/F");
      tree_->Branch("puppit1elmet"              , &puppit1elmet              , "puppit1elmet/F");
      tree_->Branch("puppit1elmetphi"           , &puppit1elmetphi           , "puppit1elmetphi/F");
      tree_->Branch("puppiphmet"                , &puppiphmet                , "puppiphmet/F");
      tree_->Branch("puppiphmetphi"             , &puppiphmetphi             , "puppiphmetphi/F");
      tree_->Branch("puppit1phmet"              , &puppit1phmet              , "puppit1phmet/F");
      tree_->Branch("puppit1phmetphi"           , &puppit1phmetphi           , "puppit1phmetphi/F");

      if(addMETSystematics){
	tree_->Branch("puppit1pfmetMuEnUp"       , &puppit1pfmetMuEnUp       , "puppit1pfmetMuEnUp/F");
	tree_->Branch("puppit1pfmetMuEnDown"     , &puppit1pfmetMuEnDown     , "puppit1pfmetMuEnDown/F");
	tree_->Branch("puppit1pfmetElEnUp"       , &puppit1pfmetElEnUp       , "puppit1pfmetElEnUp/F");
	tree_->Branch("puppit1pfmetElEnDown"     , &puppit1pfmetElEnDown     , "puppit1pfmetElEnDown/F");
	tree_->Branch("puppit1pfmetPhoEnUp"      , &puppit1pfmetPhoEnUp      , "puppit1pfmetPhoEnUp/F");
	tree_->Branch("puppit1pfmetPhoEnDown"    , &puppit1pfmetPhoEnDown    , "puppit1pfmetPhoEnDown/F");
	tree_->Branch("puppit1pfmetTauEnUp"      , &puppit1pfmetTauEnUp      , "puppit1pfmetTauEnUp/F");
	tree_->Branch("puppit1pfmetTauEnDown"    , &puppit1pfmetTauEnDown    , "puppit1pfmetTauEnDown/F");
	tree_->Branch("puppit1pfmetJetEnUp"      , &puppit1pfmetJetEnUp      , "puppit1pfmetJetEnUp/F");
	tree_->Branch("puppit1pfmetJetEnDown"    , &puppit1pfmetJetEnDown    , "puppit1pfmetJetEnDown/F");
	tree_->Branch("puppit1pfmetJetResUp"     , &puppit1pfmetJetResUp     , "puppit1pfmetJetResUp/F");
	tree_->Branch("puppit1pfmetJetResDown"   , &puppit1pfmetJetResDown   , "puppit1pfmetJetResDown/F");
	tree_->Branch("puppit1pfmetUncEnUp"      , &puppit1pfmetUncEnUp      , "puppit1pfmetUncEnUp/F");
	tree_->Branch("puppit1pfmetUncEnDown"    , &puppit1pfmetUncEnDown    , "puppit1pfmetUncEnDown/F");

	tree_->Branch("puppit1pfmetMuEnUpPhi"       , &puppit1pfmetMuEnUpPhi       , "puppit1pfmetMuEnUpPhi/F");
	tree_->Branch("puppit1pfmetMuEnDownPhi"     , &puppit1pfmetMuEnDownPhi     , "puppit1pfmetMuEnDownPhi/F");
	tree_->Branch("puppit1pfmetElEnUpPhi"       , &puppit1pfmetElEnUpPhi       , "puppit1pfmetElEnUpPhi/F");
	tree_->Branch("puppit1pfmetElEnDownPhi"     , &puppit1pfmetElEnDownPhi     , "puppit1pfmetElEnDownPhi/F");
	tree_->Branch("puppit1pfmetPhoEnUpPhi"      , &puppit1pfmetPhoEnUpPhi      , "puppit1pfmetPhoEnUpPhi/F");
	tree_->Branch("puppit1pfmetPhoEnDownPhi"    , &puppit1pfmetPhoEnDownPhi    , "puppit1pfmetPhoEnDownPhi/F");
	tree_->Branch("puppit1pfmetTauEnUpPhi"      , &puppit1pfmetTauEnUpPhi      , "puppit1pfmetTauEnUpPhi/F");
	tree_->Branch("puppit1pfmetTauEnDownPhi"    , &puppit1pfmetTauEnDownPhi    , "puppit1pfmetTauEnDownPhi/F");
	tree_->Branch("puppit1pfmetJetEnUpPhi"      , &puppit1pfmetJetEnUpPhi      , "puppit1pfmetJetEnUpPhi/F");
	tree_->Branch("puppit1pfmetJetEnDownPhi"    , &puppit1pfmetJetEnDownPhi    , "puppit1pfmetJetEnDownPhi/F");
	tree_->Branch("puppit1pfmetJetResUpPhi"     , &puppit1pfmetJetResUpPhi     , "puppit1pfmetJetResUpPhi/F");
	tree_->Branch("puppit1pfmetJetResDownPhi"   , &puppit1pfmetJetResDownPhi   , "puppit1pfmetJetResDownPhi/F");
	tree_->Branch("puppit1pfmetUncEnUpPhi"      , &puppit1pfmetUncEnUpPhi      , "puppit1pfmetUncEnUpPhi/F");
	tree_->Branch("puppit1pfmetUncEnDownPhi"    , &puppit1pfmetUncEnDownPhi    , "puppit1pfmetUncEnDownPhi/F");
      }
    }    
  }
}

