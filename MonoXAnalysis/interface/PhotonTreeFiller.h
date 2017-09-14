#ifndef AnalysisCode_MonoXAnalysis_PhotonTreeFiller_h
#define AnalysisCode_MonoXAnalysis_PhotonTreeFiller_h

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

#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"


#include "TTree.h"
#include "TLorentzVector.h"

class PhotonTreeFiller {

 public:

  ~PhotonTreeFiller(){};

  ////// -------
  PhotonTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree);
  ////// -------
  bool Fill(const edm::Event & iEvent, const edm::EventSetup& iSetup);


 private:

  ////// -------
  void DeclareAndSetBranches();
  void initBranches();

  float getGammaEAForPhotonIso(float eta);
  float getChargedHadronEAForPhotonIso(float eta);
  float getNeutralHadronEAForPhotonIso(float eta);

  const edm::InputTag rhoTag;
  edm::EDGetTokenT<double>  rhoToken;
 
  const edm::InputTag  photonsTag;
  const edm::InputTag  mediumphotonsTag;
  const edm::InputTag  tightphotonsTag;
  const edm::InputTag  photonHighPtIdTag;
  const edm::InputTag  mvaloosephotonsTag;
  const edm::InputTag  mvatightphotonsTag;
  const bool applyPhotonJetsFilter;
  const bool applyDiElectronFilter;
  const bool applyDiMuonFilter;
  const bool isQCDTree;
  const bool isPhotonPurity;
  const bool isReMiniAOD;
  const bool isTriggerTree;

  edm::EDGetTokenT<pat::PhotonRefVector>    photonsToken;
  edm::EDGetTokenT<pat::PhotonRefVector>    mediumphotonsToken;
  edm::EDGetTokenT<pat::PhotonRefVector>    tightphotonsToken;
  edm::EDGetTokenT<pat::PhotonRefVector>    mvaloosephotonsToken;
  edm::EDGetTokenT<pat::PhotonRefVector>    mvatightphotonsToken;
  edm::EDGetTokenT<pat::PhotonRefVector>    photonsPurityToken;
  edm::EDGetTokenT<edm::ValueMap<bool> >    photonHighPtIdToken;
  edm::EDGetTokenT<edm::ValueMap<float> >   photonsieieToken;
  edm::EDGetTokenT<edm::ValueMap<float> >   photonPHisoToken;
  edm::EDGetTokenT<edm::ValueMap<float> >   photonCHisoToken;
  edm::EDGetTokenT<edm::ValueMap<float> >   photonNHisoToken;
  edm::EDGetTokenT<edm::ValueMap<float> >   rndgammaiso04Token;
  edm::EDGetTokenT<edm::ValueMap<float> >   rndgammaiso08Token;
  edm::EDGetTokenT<edm::ValueMap<float> >   gammaisoToken;
  edm::EDGetTokenT<edm::ValueMap<float> >   rndchhadiso04Token;
  edm::EDGetTokenT<edm::ValueMap<float> >   rndchhadiso08Token;
 
  const bool addPhotonIDVariables;
  edm::EDGetTokenT<std::vector<pat::Photon> > photonIDCollectionToken;
 
  const edm::InputTag verticesTag;
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;

  TTree* tree_;

  int32_t phidm,phidt,phidh,phidmval, phidmvat, phgs;
  int32_t parid,ancid; 
  uint32_t nphotons,nmvaloosephotons,nmvatightphotons;
  uint32_t nphotonsPurity;

  float rho;
  float phPHiso, phCHiso, phNHiso, phPuritypt, phPurityeta, phPurityphi;
  float phPurityPHiso,phPurityRND04PHiso,phPurityRND08PHiso,phPurityCHiso,phPurityRND04CHiso,phPurityRND08CHiso,phPurityNHiso;
  float    phPuritysieie, phPurityhoe, phPurityElectronVeto, phPurityEAEGamma;
  float phpt,pheta,phphi,phe;

  std::vector<float> photonPt, photonEta, photonPhi, photonE, photonSCEta, photonSCPhi, photonSCEnergy, photonSCRawEnergy;
  std::vector<float> photonHOverE, photonSigmaIetaIeta, photonChargedIso,photonNeutralIso,photonEMIso, photonElectronVeto;

  // sorting objects
  template<typename T> 
    class PatPtSorter{
  public:
    bool operator ()(const T & i, const T & j) const {
      return (i->pt() > j->pt());
    }
  };

  PatPtSorter<pat::PhotonRef> photonPtSorter;
  

};

#endif

