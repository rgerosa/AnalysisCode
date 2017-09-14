#ifndef AnalysisCode_MonoXAnalysis_MuonTreeFiller_h
#define AnalysisCode_MonoXAnalysis_MuonTreeFiller_h

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

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"


#include "TTree.h"
#include "TLorentzVector.h"

class MuonTreeFiller {

 public:

  ~MuonTreeFiller(){};

  ////// -------
  MuonTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree);
  ////// -------
  bool Fill(const edm::Event & iEvent, const edm::EventSetup& iSetup);


 private:

  ////// -------
  void DeclareAndSetBranches();
  void initBranches();

  float computeMuonIso(const reco::Muon &);  

  const edm::InputTag muonsTag;
  const edm::InputTag tightmuonsTag;
  const edm::InputTag highptmuonsTag;

  edm::EDGetTokenT<pat::MuonRefVector> muonsToken;
  edm::EDGetTokenT<pat::MuonRefVector> tightmuonsToken;
  edm::EDGetTokenT<pat::MuonRefVector> highptmuonsToken;
  edm::EDGetTokenT<edm::ValueMap<bool> >    electronLooseIdToken;  
  edm::EDGetTokenT<edm::View<reco::Candidate> > fakeMuonCollToken;

  const bool isQCDTree;
  const bool isPhotonPurity;
  const bool isReMiniAOD;
  const bool isTriggerTree;

  const edm::InputTag t1metTag;
  edm::EDGetTokenT<edm::View<pat::MET> >    t1metToken;

  const edm::InputTag verticesTag;
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;

  TTree* tree_;

  int32_t mu1pid,mu2pid,mu1id,mu2id,mu1idm,mu2idm,mu1idt,mu2idt;

  uint32_t nmuons,ntightmuons,nhighptmuons,nmuonsfake;

  float mu1pt,mu1eta,mu1phi,mu1pfpt,mu1pfeta,mu1pfphi,mu1iso,mu2pt,mu2eta,mu2phi,mu2pfpt,mu2pfeta,mu2pfphi,mu2iso;
  std::vector<float> fakemupt, fakemueta, fakemuphi;   // fake muons in 2016 HIP mitigated data

  float zmass,zpt,zeta,zphi,wmt;

  // sorting objects
  template<typename T> 
    class PatPtSorter{
  public:
    bool operator ()(const T & i, const T & j) const {
      return (i->pt() > j->pt());
    }
  };

  PatPtSorter<pat::MuonRef> muonPtSorter;
  

};

#endif

