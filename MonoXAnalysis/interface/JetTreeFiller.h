#ifndef AnalysisCode_MonoXAnalysis_JetTreeFiller_h
#define AnalysisCode_MonoXAnalysis_JetTreeFiller_h

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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "CondTools/BTau/interface/BTagCalibrationReader.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"

#include "TTree.h"
#include "TLorentzVector.h"

class JetTreeFiller {

 public:

  ~JetTreeFiller(){};

  ////// -------
  JetTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree, const bool & isPuppi = false);
  ////// -------
  bool Fill(const edm::Event & iEvent, const edm::EventSetup& iSetup);

 private:

  ///
  void fillJetCollections(const edm::Handle<std::vector<pat::Jet> > &, 
			  const pat::MuonRefVector &, 
			  const pat::ElectronRefVector &,
			  const pat::PhotonRefVector &,
			  std::vector<pat::JetRef> &, 
			  std::vector<pat::JetRef> &, 
			  const bool & = false);

  // to apply jet ID
  bool applyJetID(const pat::Jet &, const std::string &);
  // to apply pileup-jet id
  bool applyPileupJetID(const pat::Jet &, const std::string &, const bool &);
  // fill btag scale factors
  void calculateBtagSF(const pat::Jet &, const std::string &, std::vector<float> &, std::vector<float> &, std::vector<float> &);
 
  ////// -------
  void DeclareAndSetBranches();
  void initBranches();

  const edm::InputTag muonsTag;
  edm::EDGetTokenT<pat::MuonRefVector> muonsToken;

  const edm::InputTag  electronsTag;
  edm::EDGetTokenT<pat::ElectronRefVector> electronsToken;

  const edm::InputTag  photonsTag;
  edm::EDGetTokenT<pat::PhotonRefVector> photonsToken;

  edm::InputTag jetsTag;
  edm::InputTag jetsJESUpTag;
  edm::InputTag jetsJESDwTag;
  edm::InputTag jetsJERTag;

  edm::EDGetTokenT<std::vector<pat::Jet> >  jetsToken;
  edm::EDGetTokenT<std::vector<pat::Jet> >  jetsJESUpToken;
  edm::EDGetTokenT<std::vector<pat::Jet> >  jetsJESDwToken;
  edm::EDGetTokenT<std::vector<pat::Jet> >  jetsJERToken;

  const bool isMC;
  const bool isTriggerTree;
  const bool isQCDTree;
  const bool isPhotonPurity;

  const bool cleanMuonJet;
  const bool cleanElectronJet;
  const bool cleanPhotonJet;   
  const double dRCleaningAK4;

  const std::string jetidwp;
  const std::string pileupjetidwp;

  const bool   applypileupjetid;
  const double btaggingCSVWP;
  const double btaggingMVAWP;
  const double minJetPtCountAK4;
  const double minJetPtBveto;
  const double minJetPtAK4Store;

  // B-tagging SF
  bool addBTagScaleFactor;
  edm::FileInPath bTagScaleFactorFileCSV;
  edm::FileInPath bTagScaleFactorFileMVA;

  BTagCalibration calibCSV;
  BTagCalibration calibMVA;

  std::vector<BTagCalibrationReader> bMediumCSV;
  std::vector<BTagCalibrationReader> bMediumMVA;

  ///////
  TTree* tree_;

  bool isPuppi_;

  uint32_t njets,nbjets,nbjetslowpt,nbjetsMVA,nbjetsMVAlowpt;  
  uint32_t npuppijets,npuppibjets,npuppibjetsMVA,npuppibjetslowpt,npuppibjetsMVAlowpt;
  uint32_t njetsinc,npuppijetsinc;

  uint32_t njetsincup,npuppijetsincup; 
  uint32_t njetsincdw,npuppijetsincdw;
  uint32_t njetsincjer,npuppijetsincjer;

  float ht, htinc, Puppiht, Puppihtinc;

  std::vector<float> combinejetpt,combinejeteta,combinejetphi,combinejetm,combinejetbtag,combinejetbtagMVA,combinejetmetdphi;
  std::vector<float> combinejetCHfrac,combinejetNHfrac,combinejetEMfrac,combinejetCEMfrac,combinejetPHfrac,combinejetELfrac,combinejetMUfrac, combinejetHFHfrac, combinejetHFEMfrac;
  std::vector<unsigned int> combinejetCHmult,combinejetNHmult,combinejetPHmult,combinejetELmult,combinejetMUmult,combinejetHFHmult,combinejetHFEMmult;
  std::vector<float> combinejetHFlav,combinejetPFlav,combinejetQGL,combinejetPUID, combinejetPassPUID; 
  std::vector<float> combinejetGenpt,combinejetGeneta,combinejetGenphi,combinejetGenm;
  std::vector<float> combinejetBtagSF,combinejetBtagSFUp,combinejetBtagSFDown;
  std::vector<float> combinejetBtagMVASF,combinejetBtagMVASFUp,combinejetBtagMVASFDown;
  std::vector<float> combinejetptup,  combinejetptdw,  combinejetptjer,  combinejetetaup, combinejetetadw, combinejetetajer;
  std::vector<float> combinejetphiup, combinejetphidw, combinejetphijer, combinejetmup,   combinejetmdw,   combinejetmjer;

  std::vector<float> combinePuppijetpt,combinePuppijeteta,combinePuppijetphi,combinePuppijetm,combinePuppijetbtag,combinePuppijetbtagMVA;
  std::vector<float> combinePuppijetCHfrac,combinePuppijetNHfrac,combinePuppijetEMfrac,combinePuppijetCEMfrac,combinePuppijetmetdphi;
  std::vector<float> combinePuppijetHFlav,combinePuppijetPFlav,combinePuppijetQGL;
  std::vector<float> combinePuppijetGenpt,combinePuppijetGeneta,combinePuppijetGenphi,combinePuppijetGenm;
  std::vector<float> combinePuppijetBtagSF,combinePuppijetBtagSFUp,combinePuppijetBtagSFDown;
  std::vector<float> combinePuppijetBtagMVASF,combinePuppijetBtagMVASFUp,combinePuppijetBtagMVASFDown;
  std::vector<float> combinePuppijetptup,  combinePuppijetptdw,  combinePuppijetptjer,  combinePuppijetetaup, combinePuppijetetadw, combinePuppijetetajer;
  std::vector<float> combinePuppijetphiup, combinePuppijetphidw, combinePuppijetphijer, combinePuppijetmup,   combinePuppijetmdw,   combinePuppijetmjer;


  // sorting objects
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

