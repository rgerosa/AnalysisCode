#ifndef AnalysisCode_MonoXAnalysis_RecoilTreeFiller_h
#define AnalysisCode_MonoXAnalysis_RecoilTreeFiller_h

// basic C++ headers
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

// FWCore
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/MET.h"


#include "TTree.h"
#include "TLorentzVector.h"

class RecoilTreeFiller {

 public:

  ~RecoilTreeFiller(){};

  ////// -------
  RecoilTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree, const bool & isPuppi = false);
  ////// -------
  bool Fill(const edm::Event & iEvent, const edm::EventSetup& iSetup);


 private:

  ////// -------
  void DeclareAndSetBranches();
  void initBranches();


  edm::InputTag t1metTag;
  edm::InputTag t1mumetTag;
  edm::InputTag t1elmetTag;
  edm::InputTag t1phmetTag;
  edm::InputTag t1taumetTag;
  
  // re-miniAOD 2016
  edm::InputTag t1metEGCleanTag;
  edm::InputTag t1metMuCleanTag;
  edm::InputTag t1metOriginalTag;
  edm::InputTag t1mumetEGCleanTag;
  edm::InputTag t1mumetMuCleanTag;
  edm::InputTag t1elmetEGCleanTag;
  edm::InputTag t1elmetMuCleanTag;
  edm::InputTag t1phmetEGCleanTag;
  edm::InputTag t1phmetMuCleanTag;
  edm::InputTag t1taumetEGCleanTag;
  edm::InputTag t1taumetMuCleanTag;

  edm::EDGetTokenT<edm::View<pat::MET> >  t1metToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1metEGCleanToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1metMuCleanToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1metOriginalToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1mumetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1mumetEGCleanToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1mumetMuCleanToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1elmetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1elmetEGCleanToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1elmetMuCleanToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1phmetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1phmetEGCleanToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1phmetMuCleanToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1taumetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1taumetEGCleanToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1taumetMuCleanToken;


  // MET breakdown
  const bool addMETBreakDown;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetHadronHFToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetEgammaHFToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetChargedHadronToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetNeutralHadronToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetElectronsToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetPhotonsToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetMuonsToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  pfMetUnclusteredToken;

  //Systematics
  bool addMETSystematics;

  const bool isMC;
  const bool addBadMuonClean;
  const bool isReMiniAOD;
  const bool isQCDTree;
  const bool isPhotonPurity;
  const bool isTriggerTree;

  bool isPuppi_;
  TTree* tree_;
  
  ///
  float t1pfmet,t1pfmetphi,t1mumet,t1mumetphi,t1elmet,t1elmetphi,t1phmet,t1phmetphi,t1taumet,t1taumetphi;
  float t1pfmetEGClean,t1pfmetphiEGClean,t1mumetEGClean,t1mumetphiEGClean,t1elmetEGClean,t1elmetphiEGClean,t1phmetEGClean,t1phmetphiEGClean,t1taumetEGClean,t1taumetphiEGClean;
  float t1pfmetMuClean,t1pfmetphiMuClean,t1mumetMuClean,t1mumetphiMuClean,t1elmetMuClean,t1elmetphiMuClean,t1phmetMuClean,t1phmetphiMuClean,t1taumetMuClean,t1taumetphiMuClean;
  float t1pfmetOriginal, t1pfmetphiOriginal;
  float pfmet,pfmetphi,mumet,mumetphi,elmet,elmetphi,phmet,phmetphi,taumet,taumetphi;
  float calomet, calometphi;

  // MET break down
  float pfmethadronHF,pfmethadronHFphi,pfmetegammaHF,pfmetegammaHFphi,pfmetchargedhadron,pfmetchargedhadronphi;
  float pfmetneutralhadron,pfmetneutralhadronphi,pfmetelectrons,pfmetelectronsphi,pfmetmuons,pfmetmuonsphi,pfmetphotons,pfmetphotonsphi,pfmetunclustered,pfmetunclusteredphi;

  // Puppi MET info (typeI and Raw)
  float puppipfmet,puppipfmetphi,puppimumet,puppimumetphi,puppielmet,puppielmetphi,puppiphmet,puppiphmetphi,puppitaumet,puppitaumetphi;
  float puppit1pfmet,puppit1pfmetphi,puppit1mumet,puppit1mumetphi,puppit1elmet,puppit1elmetphi,puppit1phmet,puppit1phmetphi,puppit1taumet,puppit1taumetphi;
  
  // gen met
  float genmet,genmetphi;

  // met systematics
  float t1pfmetMuEnUp,t1pfmetMuEnDown,t1pfmetElEnUp,t1pfmetElEnDown,t1pfmetPhoEnUp,t1pfmetPhoEnDown,t1pfmetTauEnUp,t1pfmetTauEnDown;
  float t1pfmetJetEnUp,t1pfmetJetEnDown,t1pfmetJetResUp,t1pfmetJetResDown,t1pfmetUncEnUp,t1pfmetUncEnDown,t1pfmetJetSmear,t1pfmetXY;
  float t1pfmetMuEnUpPhi,t1pfmetMuEnDownPhi,t1pfmetElEnUpPhi,t1pfmetElEnDownPhi,t1pfmetPhoEnUpPhi,t1pfmetPhoEnDownPhi,t1pfmetTauEnUpPhi,t1pfmetTauEnDownPhi;
  float t1pfmetJetEnUpPhi,t1pfmetJetEnDownPhi,t1pfmetJetResUpPhi,t1pfmetJetResDownPhi,t1pfmetUncEnUpPhi,t1pfmetUncEnDownPhi,t1pfmetJetSmearPhi,t1pfmetXYPhi;

  // met systematics puppi
  float puppit1pfmetMuEnUp,puppit1pfmetMuEnDown,puppit1pfmetElEnUp,puppit1pfmetElEnDown,puppit1pfmetPhoEnUp,puppit1pfmetPhoEnDown,puppit1pfmetTauEnUp;
  float puppit1pfmetTauEnDown,puppit1pfmetJetEnUp,puppit1pfmetJetEnDown,puppit1pfmetJetResUp,puppit1pfmetJetResDown,puppit1pfmetUncEnUp,puppit1pfmetUncEnDown;
  float puppit1pfmetMuEnUpPhi,puppit1pfmetMuEnDownPhi,puppit1pfmetElEnUpPhi,puppit1pfmetElEnDownPhi;
  float puppit1pfmetPhoEnUpPhi,puppit1pfmetPhoEnDownPhi,puppit1pfmetTauEnUpPhi,puppit1pfmetTauEnDownPhi;
  float puppit1pfmetJetEnUpPhi,puppit1pfmetJetEnDownPhi,puppit1pfmetJetResUpPhi,puppit1pfmetJetResDownPhi;
  float puppit1pfmetUncEnUpPhi,puppit1pfmetUncEnDownPhi,puppit1pfmetJetSmearPhi,puppit1pfmetXYPhi;
    
};

#endif

