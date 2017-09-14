#ifndef AnalysisCode_MonoXAnalysis_ElectronTreeFiller_h
#define AnalysisCode_MonoXAnalysis_ElectronTreeFiller_h

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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TTree.h"
#include "TLorentzVector.h"

class ElectronTreeFiller {

 public:

  ~ElectronTreeFiller(){};

  ////// -------
  ElectronTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree);
  ////// -------
  bool Fill(const edm::Event & iEvent, const edm::EventSetup& iSetup);


 private:

  ////// -------
  void DeclareAndSetBranches();
  void initBranches();

  const edm::InputTag rhoTag;
  edm::EDGetTokenT<double>  rhoToken;
 
  const edm::InputTag  electronsTag;
  const edm::InputTag  looseelectronsTag;
  const edm::InputTag  tightelectronsTag;
  const edm::InputTag  triggerelectronsTag;
  const edm::InputTag  heepelectronsTag;
  const edm::InputTag  mvalooseelectronsTag;
  const edm::InputTag  mvatightelectronsTag;

  edm::EDGetTokenT<pat::ElectronRefVector>  electronsToken;
  edm::EDGetTokenT<pat::ElectronRefVector>  looseelectronsToken;
  edm::EDGetTokenT<pat::ElectronRefVector>  tightelectronsToken;
  edm::EDGetTokenT<pat::ElectronRefVector>  triggerelectronsToken;
  edm::EDGetTokenT<pat::ElectronRefVector>  heepelectronsToken;
  edm::EDGetTokenT<pat::ElectronRefVector>  mvalooseelectronsToken;
  edm::EDGetTokenT<pat::ElectronRefVector>  mvatightelectronsToken;
  edm::EDGetTokenT<edm::ValueMap<bool> >    electronLooseIdToken;  

  const bool applyDiElectronFilter;
  const bool applyDiMuonFilter;
  const bool applyPhotonJetsFilter;
  const bool isQCDTree;
  const bool isPhotonPurity;
  const bool isReMiniAOD;
  const bool isTriggerTree;

  const bool addElectronIDVariables;
  edm::EDGetTokenT<std::vector<pat::Electron> > electronIDCollectionToken;

  const edm::InputTag t1metTag;
  edm::EDGetTokenT<edm::View<pat::MET> >    t1metToken;

  const edm::InputTag verticesTag;
  edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;

  TTree* tree_;

  int32_t el1pid,el2pid,el1id,el1idl,el1idt,el2id,el2idl,el2idt,el1idmval, el2idmval, el1idmvat, el2idmvat;
  int32_t el1gs, el2gs;
  uint32_t nelectrons,nlooseelectrons,ntightelectrons,nheepelectrons,ntriggerelectrons,nmvalooseelectrons,nmvatightelectrons;
  float rho;
  float el1pt,el1eta,el1phi,ele1e,el2pt,ele2e,el2eta,el2phi,phpt,pheta,phphi,phe;
  float zeemass,zeept,zeeeta,zeephi,wemt;

  std::vector<float> electronPt, electronEta, electronPhi, electronE, electronSCEta, electronSCPhi, electronSCEnergy;
  std::vector<float> electronSCRawEnergy, electronHOverE, electronSigmaIetaIeta, electronChargedIso, electronNeutralIso, electronEMIso, electronGsfPt;
  std::vector<float> electronEOP, electronDphi, electronDeta, electronMissHit, electronConversion, electronDxy, electronDz;

  // sorting objects
  template<typename T> 
    class PatPtSorter{
  public:
    bool operator ()(const T & i, const T & j) const {
      return (i->pt() > j->pt());
    }
  };

  PatPtSorter<pat::ElectronRef> electronPtSorter;
  

};

#endif

