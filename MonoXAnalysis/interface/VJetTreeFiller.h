#ifndef AnalysisCode_MonoXAnalysis_VJetTreeFiller_h
#define AnalysisCode_MonoXAnalysis_VJetTreeFiller_h

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

#include "AnalysisCode/MonoXAnalysis/interface/TreeFillerUtils.h"

class VJetTreeFiller {

 public:

  ~VJetTreeFiller(){};

  ////// -------
  VJetTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree, const bool & isPuppi = false);
  ////// -------
  bool Fill(const edm::Event & iEvent, const edm::EventSetup& iSetup);

 private:

  ///
  void fillJetCollections(const edm::Handle<std::vector<pat::Jet> > &, 
			  const pat::MuonRefVector &, 
			  const pat::ElectronRefVector &,
			  const pat::PhotonRefVector &,
			  std::vector<pat::JetRef> &, 
			  const bool & = false);

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
  edm::EDGetTokenT<std::vector<pat::Jet> >  jetsToken;
  TString jetsLabel;
  

  const bool isMC;
  const bool isTriggerTree;
  const bool isQCDTree;
  const bool isPhotonPurity;
  const bool applyDiMuonFilter;
  const bool applyDiElectronFilter;
  const bool applyPhotonJetsFilter;

  const bool cleanMuonJet;
  const bool cleanElectronJet;
  const bool cleanPhotonJet;   
  const double dRCleaningAK8;

  const std::string jetidwp;
  const double btaggingCSVWP;

  // B-tagging SF
  bool addBTagScaleFactor;
  edm::FileInPath bTagScaleFactorFileSubCSV;
  BTagCalibration calibSubCSV;
  std::vector<BTagCalibrationReader> bMediumSubCSV;

  const bool useMiniAODSubstructure;


  ///////
  TTree* tree_;

  bool isPuppi_;

  // AK8 jet
  std::vector<float> boostedJetpt,boostedJeteta,boostedJetphi,boostedJetm;
  std::vector<float> boostedJetGenpt,boostedJetGenm,boostedJetGeneta,boostedJetGenphi;
  std::vector<float> boostedJettau1,boostedJettau2,boostedJettau3,boostedJettau4;
  std::vector<float> boostedJetGentau1,boostedJetGentau2,boostedJetGentau3,boostedJetGentau4;
  std::vector<float> boostedJetHFlav,boostedJetPFlav,boostedJetQGL,boostedJetBtag,boostedJetDoubleBtag;

  // AK8 pruned-jet
  std::vector<float> prunedJetpt,prunedJetm,prunedJetphi,prunedJeteta;
  std::vector<float> prunedJetm_v2, prunedJetpt_v2, prunedJetphi_v2, prunedJeteta_v2;
  std::vector<float> prunedJetGenpt,prunedJetGenm,prunedJetGeneta,prunedJetGenphi;
  std::vector<float> prunedJetptraw,prunedJetmraw;
  std::vector<float> prunedJetHFlav,prunedJetPFlav,prunedJetQGL,prunedJetBtag,prunedJetDoubleBtag;

  // AK8 pruned subjets
  std::vector<float> prunedSubJetpt_1,prunedSubJetm_1,prunedSubJetphi_1,prunedSubJeteta_1;
  std::vector<float> prunedSubJetHFlav_1,prunedSubJetQGL_1,prunedSubJetBtag_1;
  std::vector<float> prunedSubJetGenpt_1,prunedSubJetGenm_1,prunedSubJetGeneta_1,prunedSubJetGenphi_1,prunedSubJetPFlav_1;
  std::vector<float> prunedSubJetptraw_1,prunedSubJetmraw_1;
  std::vector<float> prunedSubJetBtagSF_1,prunedSubJetBtagSFUp_1,prunedSubJetBtagSFDown_1;

  std::vector<float> prunedSubJetpt_2,prunedSubJetm_2,prunedSubJetphi_2,prunedSubJeteta_2,prunedSubJetHFlav_2,prunedSubJetQGL_2,prunedSubJetBtag_2;
  std::vector<float> prunedSubJetGenpt_2,prunedSubJetGenm_2,prunedSubJetGeneta_2,prunedSubJetGenphi_2,prunedSubJetPFlav_2;
  std::vector<float> prunedSubJetptraw_2,prunedSubJetmraw_2;
  std::vector<float> prunedSubJetBtagSF_2,prunedSubJetBtagSFUp_2,prunedSubJetBtagSFDown_2;

  // AK8 puppi jets
  std::vector<float> boostedPuppiJetpt,boostedPuppiJeteta,boostedPuppiJetphi,boostedPuppiJetm;
  std::vector<float> boostedPuppiJetGenpt,boostedPuppiJetGenm,boostedPuppiJetGeneta,boostedPuppiJetGenphi;
  std::vector<float> boostedPuppiJettau1,boostedPuppiJettau2,boostedPuppiJettau3,boostedPuppiJettau4;
  std::vector<float> boostedPuppiJetGentau1,boostedPuppiJetGentau2,boostedPuppiJetGentau3,boostedPuppiJetGentau4;
  std::vector<float> boostedPuppiJetHFlav,boostedPuppiJetPFlav,boostedPuppiJetQGL,boostedPuppiJetBtag,boostedPuppiJetDoubleBtag;

  // AK8 puppi soft-dropped jets
  std::vector<float> softDropPuppiJetpt,softDropPuppiJetm,softDropPuppiJeteta,softDropPuppiJetphi;
  std::vector<float> softDropPuppiJetm_v2, softDropPuppiJetpt_v2, softDropPuppiJeteta_v2, softDropPuppiJetphi_v2;
  std::vector<float> softDropPuppiJetGenpt,softDropPuppiJetGenm,softDropPuppiJetGeneta,softDropPuppiJetGenphi;
  std::vector<float> softDropPuppiJetHFlav,softDropPuppiJetPFlav,softDropPuppiJetQGL,softDropPuppiJetBtag,softDropPuppiJetDoubleBtag;
  std::vector<float> softDropPuppiJetptraw,softDropPuppiJetmraw;

  // soft drop subjets
  std::vector<float> softDropPuppiSubJetpt_1,softDropPuppiSubJetm_1,softDropPuppiSubJetphi_1,softDropPuppiSubJeteta_1;
  std::vector<float> softDropPuppiSubJetHFlav_1,softDropPuppiSubJetQGL_1,softDropPuppiSubJetBtag_1;
  std::vector<float> softDropPuppiSubJetGenpt_1,softDropPuppiSubJetGenm_1,softDropPuppiSubJetGeneta_1,softDropPuppiSubJetGenphi_1,softDropPuppiSubJetPFlav_1;
  std::vector<float> softDropPuppiSubJetptraw_1,softDropPuppiSubJetmraw_1;
  std::vector<float> softDropPuppiSubJetBtagSF_1,softDropPuppiSubJetBtagSFUp_1,softDropPuppiSubJetBtagSFDown_1;

  std::vector<float> softDropPuppiSubJetpt_2,softDropPuppiSubJetm_2,softDropPuppiSubJetphi_2,softDropPuppiSubJeteta_2,softDropPuppiSubJetHFlav_2; 
  std::vector<float> softDropPuppiSubJetQGL_2,softDropPuppiSubJetBtag_2;
  std::vector<float> softDropPuppiSubJetGenpt_2,softDropPuppiSubJetGenm_2,softDropPuppiSubJetGeneta_2,softDropPuppiSubJetGenphi_2,softDropPuppiSubJetPFlav_2;
  std::vector<float> softDropPuppiSubJetptraw_2,softDropPuppiSubJetmraw_2;
  std::vector<float> softDropPuppiSubJetBtagSF_2,softDropPuppiSubJetBtagSFUp_2,softDropPuppiSubJetBtagSFDown_2;


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

