#ifndef AnalysisCode_MonoXAnalysis_JetMetDphiTreeFiller_h
#define AnalysisCode_MonoXAnalysis_JetMetDphiTreeFiller_h

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
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/METReco/interface/MET.h"


#include "TTree.h"
#include "TLorentzVector.h"

class JetMetDphiTreeFiller {

 public:

  ~JetMetDphiTreeFiller(){};

  ////// -------
  JetMetDphiTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree, const bool & isPuppi = false);
  ////// -------
  bool Fill(const edm::Event & iEvent, const edm::EventSetup& iSetup);
  ////// -------
  void fillJetCollections(const edm::Handle<std::vector<pat::Jet> > &,
                          const pat::MuonRefVector &,
                          const pat::ElectronRefVector &,
                          const pat::PhotonRefVector &,
                          std::vector<pat::JetRef> &,
                          std::vector<pat::JetRef> &,
                          const bool & = false);

  // to apply jet ID                                                                                                                                                                                   
  bool applyJetID(const pat::Jet &, const std::string &);


 private:

  ////// -------
  void DeclareAndSetBranches();
  void initBranches();

  const edm::InputTag muonsTag;
  edm::EDGetTokenT<pat::MuonRefVector> muonsToken;

  const edm::InputTag  electronsTag;
  edm::EDGetTokenT<pat::ElectronRefVector> electronsToken;

  const edm::InputTag  photonsTag;
  edm::EDGetTokenT<pat::PhotonRefVector> photonsToken;


  edm::InputTag t1metTag;
  edm::InputTag t1mumetTag;
  edm::InputTag t1elmetTag;
  edm::InputTag t1phmetTag;
  edm::InputTag t1taumetTag;
  
  edm::EDGetTokenT<edm::View<pat::MET> >  t1metToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1mumetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1elmetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1phmetToken;
  edm::EDGetTokenT<edm::View<pat::MET> >  t1taumetToken;

  edm::InputTag jetsTag;
  edm::InputTag jetsJESUpTag;
  edm::InputTag jetsJESDwTag;
  edm::InputTag jetsJERTag;

  edm::EDGetTokenT<std::vector<pat::Jet> >  jetsToken;
  edm::EDGetTokenT<std::vector<pat::Jet> >  jetsJESUpToken;
  edm::EDGetTokenT<std::vector<pat::Jet> >  jetsJESDwToken;
  edm::EDGetTokenT<std::vector<pat::Jet> >  jetsJERToken;


  //Systematics
  const bool isMC;
  const bool isTriggerTree;
  const bool isQCDTree;
  const bool isPhotonPurity;
  const bool cleanMuonJet;
  const bool cleanElectronJet;
  const bool cleanPhotonJet;
  const double dRCleaningAK4;
  const std::string jetidwp;
  const double minJetPtCountAK4;
  const bool addMETSystematics;

  bool isPuppi_;

  TTree* tree_;
  
  ///
  float incjetmetdphimin,incjetmumetdphimin,incjetelmetdphimin,incjetphmetdphimin;
  float incjetmetdphimin4,incjetmumetdphimin4,incjetelmetdphimin4,incjetphmetdphimin4; 
  float alljetmetdphimin,alljetmetdphimin4,alljetmumetdphimin,alljetmumetdphimin4,alljetelmetdphimin,alljetelmetdphimin4,alljetphmetdphimin,alljetphmetdphimin4;
  
  float incjetmetdphimin4up,incjetmumetdphimin4up,incjetelmetdphimin4up,incjetphmetdphimin4up; 
  float incjetmetdphimin4dw,incjetmumetdphimin4dw,incjetelmetdphimin4dw,incjetphmetdphimin4dw; 
  float incjetmetdphimin4jer,incjetmumetdphimin4jer,incjetelmetdphimin4jer,incjetphmetdphimin4jer; 

  float incPuppijetmetdphimin,incPuppijetmumetdphimin,incPuppijetelmetdphimin,incPuppijetphmetdphimin;
  float incPuppijetmetdphimin4,incPuppijetmumetdphimin4,incPuppijetelmetdphimin4,incPuppijetphmetdphimin4; 
  float allPuppijetmetdphimin,allPuppijetmetdphimin4,allPuppijetmumetdphimin,allPuppijetmumetdphimin4,allPuppijetelmetdphimin,allPuppijetelmetdphimin4,allPuppijetphmetdphimin,allPuppijetphmetdphimin4;
  
  float incPuppijetmetdphimin4up,incPuppijetmumetdphimin4up,incPuppijetelmetdphimin4up,incPuppijetphmetdphimin4up; 
  float incPuppijetmetdphimin4dw,incPuppijetmumetdphimin4dw,incPuppijetelmetdphimin4dw,incPuppijetphmetdphimin4dw; 
  float incPuppijetmetdphimin4jer,incPuppijetmumetdphimin4jer,incPuppijetelmetdphimin4jer,incPuppijetphmetdphimin4jer; 

  template<typename T>
    class PatPtSorter{
  public:
    bool operator ()(const T & i, const T & j) const {
      return (i->pt() > j->pt());
    }

  };

  PatPtSorter<pat::JetRef> jetPtSorter;
  
};

#endif

