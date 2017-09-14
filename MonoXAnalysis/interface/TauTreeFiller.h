#ifndef AnalysisCode_MonoXAnalysis_TauTreeFiller_h
#define AnalysisCode_MonoXAnalysis_TauTreeFiller_h

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

#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "TTree.h"
#include "TLorentzVector.h"

class TauTreeFiller {

 public:

  ~TauTreeFiller(){};

  ////// -------
  TauTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree);
  ////// -------
  bool Fill(const edm::Event & iEvent, const edm::EventSetup& iSetup);


 private:

  ////// -------
  void DeclareAndSetBranches();
  void initBranches();

  const edm::InputTag muonsTag;
  edm::EDGetTokenT<pat::MuonRefVector> muonsToken;

  const edm::InputTag  electronsTag;
  edm::EDGetTokenT<pat::ElectronRefVector> electronsToken;

  const edm::InputTag tausCollection;
  const edm::InputTag tausVLNewTag;
  const edm::InputTag tausVLOldTag;
  const edm::InputTag tausRawNewTag;
  const edm::InputTag tausRawOldTag;
  const edm::InputTag tausTightNewTag;
  const edm::InputTag tausTightOldTag;

  edm::EDGetTokenT<pat::TauCollection> tausCollectionToken;
  edm::EDGetTokenT<pat::TauRefVector>  tausVLNewToken;
  edm::EDGetTokenT<pat::TauRefVector>  tausVLOldToken;
  edm::EDGetTokenT<pat::TauRefVector>  tausRawNewToken;
  edm::EDGetTokenT<pat::TauRefVector>  tausRawOldToken;
  edm::EDGetTokenT<pat::TauRefVector>  tausTightNewToken;
  edm::EDGetTokenT<pat::TauRefVector>  tausTightOldToken;

  const bool applyPhotonJetsFilter;
  const bool applyDiElectronFilter;
  const bool applyDiMuonFilter;
  const bool isQCDTree;
  const bool isPhotonPurity;
  const bool isReMiniAOD;
  const bool isTriggerTree;

  const bool cleanMuonJet;
  const bool cleanElectronJet;
  const float dRCleaningAK4;

  TTree* tree_;

  uint32_t ntaus,ntausraw,ntausold,ntausrawold;

  std::vector<int>   combinetauidnew, combinetauidold, combinetaupid;
  std::vector<float> combinetaupt,combinetaueta,combinetauphi,combinetaum;
  std::vector<float> combinetauGenpt, combinetauGeneta, combinetauGenphi, combinetauGenm;

  // sorting objects
  template<typename T> 
    class PatPtSorter{
  public:
    bool operator ()(const T & i, const T & j) const {
      return (i->pt() > j->pt());
    }

  };

  PatPtSorter<pat::TauRef> tauPtSorter;
    

};

#endif

