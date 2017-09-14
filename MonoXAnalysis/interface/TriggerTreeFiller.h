#ifndef AnalysisCode_MonoXAnalysis_TriggerTreeFiller_h
#define AnalysisCode_MonoXAnalysis_TriggerTreeFiller_h

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
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

// HLT info
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

// DataFormats
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/L1Trigger/interface/EGamma.h" 
#include "DataFormats/L1Trigger/interface/Tau.h"    
#include "DataFormats/L1Trigger/interface/Jet.h"    
#include "DataFormats/L1Trigger/interface/Muon.h"   
#include "DataFormats/L1Trigger/interface/EtSum.h"  
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "CondFormats/L1TObjects/interface/L1TUtmAlgorithm.h"
#include "CondFormats/L1TObjects/interface/L1TUtmTriggerMenu.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "TTree.h"
#include "TPRegexp.h"

class TriggerTreeFiller {

 public:

  ~TriggerTreeFiller(){};

  ////// -------
  TriggerTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree);
  ////// -------
  bool Fill(const edm::Event & iEvent, const edm::EventSetup& iSetup);
  ////// -------
  void SetTriggerPaths(edm::Run const& iRun, edm::EventSetup const& iSetup);


 private:
  ////// -------

  bool fillTriggerInfo(const edm::Handle<edm::TriggerResults> &, 
		       const edm::Handle<pat::PackedTriggerPrescales> &,  
		       const bool &, 
		       const std::vector<std::string> &, 
		       const edm::TriggerNames &);

  ////// -------
  void fillTriggerObjects(const edm::Handle<pat::TriggerObjectStandAloneCollection> & triggerObjectsH, 
			  const edm::TriggerNames & trignames);
  
  
  ////// -------
  void fillTriggerL1(const edm::Handle<l1t::EGammaBxCollection> & H_L1EG,  
		     const edm::Handle<l1t::TauBxCollection>  & H_L1Tau,
		     const edm::Handle<l1t::JetBxCollection>    & H_L1Jet, 
		     const edm::Handle<l1t::MuonBxCollection> & H_L1Mu,
		     const edm::Handle<l1t::EtSumBxCollection>  & H_L1Sums);
  
  ////// -------
  void fillAlgosL1(const edm::Event& iEvent, 
		   const edm::EventSetup & eventSetup, 
		   const edm::Handle<GlobalAlgBlkBxCollection> & H_L1Algos);


  void DeclareAndSetBranches();
  void initBranches();

  const bool isTriggerTree;
  const bool addTriggerObjects;
  const bool addNonETMTriggerObjects;

  const edm::InputTag triggerResultsTag;
  const edm::InputTag prescalesTag;
  const edm::InputTag triggerObjectsTag; 
  const edm::InputTag IT_L1_EG;  
  const edm::InputTag IT_L1_Jet; 
  const edm::InputTag IT_L1_Mu;  
  const edm::InputTag IT_L1_Sums;
  const edm::InputTag IT_L1_Algos; 

  edm::EDGetTokenT<edm::TriggerResults>                    triggerResultsToken;
  edm::EDGetTokenT<pat::PackedTriggerPrescales>            triggerPrescalesToken;  
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken; 
  edm::EDGetTokenT<l1t::EGammaBxCollection>  T_L1EG;
  edm::EDGetTokenT<l1t::JetBxCollection>     T_L1Jet;
  edm::EDGetTokenT<l1t::MuonBxCollection>    T_L1Mu;
  edm::EDGetTokenT<l1t::EtSumBxCollection>   T_L1Sums;
  edm::EDGetTokenT<GlobalAlgBlkBxCollection> T_L1Algos;

  const bool isQCDTree;
  const bool isPhotonPurity;
  const bool isMC;
  const bool applyHLTFilter;
  const bool setHLTFilterFlag;
  const float minL1EG;
  const float minL1Jet;
  const float minL1Mu;

  std::vector<std::string>   triggerPathsVector;
  std::map<std::string, int> triggerPathsMap;
  std::unique_ptr<HLTPrescaleProvider> hltPrescaleProvider;

  TTree* tree_;
  
  // MET triggers
  uint8_t hltmet90,hltmet100,hltmet110,hltmet120;
  uint8_t hltmetwithmu90,hltmetwithmu100,hltmetwithmu110,hltmetwithmu120,hltmetwithmu170,hltmetwithmu300;
  uint8_t hltjetmet;

  // photon trigger
  uint8_t hltphoton165,hltphoton175,hltphoton120,hltphoton90,hltphoton120vbf,hltphoton90PFHT;

  // lepton trigger
  uint8_t hltdoublemu,hltsinglemu,hltdoubleel,hltsingleel,hltsingleel27,hltelnoiso;

  // PF HT trigger
  uint8_t hltPFHT125, hltPFHT200, hltPFHT250, hltPFHT300, hltPFHT350;
  uint8_t hltPFHT400, hltPFHT475, hltPFHT600, hltPFHT650, hltPFHT800,hltPFHT900;
  uint8_t hltEcalHT800;

  //pre-scales for PF-HT
  float pswgt_ph120,pswgt_ph90;
  float pswgt_ht125,pswgt_ht200,pswgt_ht250,pswgt_ht300,pswgt_ht350;  
  float pswgt_ht400,pswgt_ht475,pswgt_ht600,pswgt_ht650,pswgt_ht800,pswgt_ht900;
  
  // Trigger objects //ND
  uint32_t                   trig_obj_n; 
  std::vector<float>         trig_obj_pt, trig_obj_eta, trig_obj_phi; 
  std::vector< std::string > trig_obj_col; 

  int trig_L1A_check;
  int trig_L1A_n;
  std::vector< std::string > trig_L1A_list;

  std::vector<float> trig_L1EG_pt  , trig_L1EG_eta  , trig_L1EG_phi  ; 
  std::vector<float> trig_L1Jet_pt , trig_L1Jet_eta , trig_L1Jet_phi ; 
  std::vector<float> trig_L1Mu_pt  , trig_L1Mu_eta  , trig_L1Mu_phi  ; 
  float              trig_L1ETM_pt , trig_L1ETM_phi , trig_L1HTM_pt  , trig_L1HTM_phi;
  float              trig_L1ETT_pt , trig_L1ETT_phi , trig_L1HTT_pt  , trig_L1HTT_phi; 
  

};

#endif

