#include "AnalysisCode/MonoXAnalysis/interface/TriggerTreeFiller.h"

TriggerTreeFiller::TriggerTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree):
  isTriggerTree (iConfig.existsAs<bool>("isTriggerTree") ? iConfig.getParameter<bool>("isTriggerTree") : false),
  addTriggerObjects (iConfig.existsAs<bool>("addTriggerObjects") ? iConfig.getParameter<bool>("addTriggerObjects") : false),
  addNonETMTriggerObjects (iConfig.existsAs<bool>("addNonETMTriggerObjects") ? iConfig.getParameter<bool>("addNonETMTriggerObjects") : false),
  triggerResultsTag  (iConfig.getParameter<edm::InputTag>("triggerResults")),
  prescalesTag       (iConfig.getParameter<edm::InputTag>("prescales")),
  triggerObjectsTag  (iConfig.existsAs<edm::InputTag>("triggerObjects") ? iConfig.getParameter<edm::InputTag>("triggerObjects") : edm::InputTag("")),
  IT_L1_EG           (iConfig.existsAs<edm::InputTag>("triggerL1EG")   ?  iConfig.getParameter<edm::InputTag>("triggerL1EG")   : edm::InputTag("")), 
  IT_L1_Jet          (iConfig.existsAs<edm::InputTag>("triggerL1Jet")  ?  iConfig.getParameter<edm::InputTag>("triggerL1Jet")  : edm::InputTag("")), 
  IT_L1_Mu           (iConfig.existsAs<edm::InputTag>("triggerL1Mu")   ?  iConfig.getParameter<edm::InputTag>("triggerL1Mu")   : edm::InputTag("")),
  IT_L1_Sums         (iConfig.existsAs<edm::InputTag>("triggerL1Sums") ?  iConfig.getParameter<edm::InputTag>("triggerL1Sums") : edm::InputTag("")), 
  IT_L1_Algos        (iConfig.existsAs<edm::InputTag>("triggerL1algos") ? iConfig.getParameter<edm::InputTag>("triggerL1algos") : edm::InputTag("gtStage2Digis")),
  isQCDTree   (iConfig.getParameter<bool>("isQCDTree")),
  isPhotonPurity (iConfig.existsAs<bool>("isPhotonPurity") ? iConfig.getParameter<bool>("isPhotonPurity") : false),
  isMC           (iConfig.existsAs<bool>("isMC") ? iConfig.getParameter<bool>("isMC") : false),
  applyHLTFilter (iConfig.existsAs<bool>("applyHLTFilter") ? iConfig.getParameter<bool>("applyHLTFilter") : false),
  setHLTFilterFlag(iConfig.existsAs<bool>("setHLTFilterFlag") ? iConfig.getParameter<bool>("setHLTFilterFlag") : false),
  minL1EG(iConfig.existsAs<double>("minL1EG") ? iConfig.getParameter<double>("minL1EG") : 20),
  minL1Jet(iConfig.existsAs<double>("minL1Jet") ? iConfig.getParameter<double>("minL1Jet") : 50),
  minL1Mu(iConfig.existsAs<double>("minL1Mu") ? iConfig.getParameter<double>("minL1Mu") : 15){

  // trigger tokens
  triggerResultsToken   = iC.consumes<edm::TriggerResults> (triggerResultsTag);
  triggerPrescalesToken = iC.consumes<pat::PackedTriggerPrescales>(prescalesTag);

  if(isTriggerTree and addTriggerObjects){
    triggerObjectsToken   = iC.consumes<pat::TriggerObjectStandAloneCollection>(triggerObjectsTag); 
    T_L1EG                = iC.consumes<l1t::EGammaBxCollection>(IT_L1_EG  ); //ND
    T_L1Jet               = iC.consumes<l1t::JetBxCollection   >(IT_L1_Jet ); //ND
    T_L1Mu                = iC.consumes<l1t::MuonBxCollection  >(IT_L1_Mu  ); //ND
    T_L1Sums              = iC.consumes<l1t::EtSumBxCollection >(IT_L1_Sums); //ND
    T_L1Algos             = iC.consumes<GlobalAlgBlkBxCollection>(IT_L1_Algos); //ND
  }

  tree_ = tree;
  DeclareAndSetBranches();
    
}

/////
void TriggerTreeFiller::initBranches(){
  
  hltmet90        = 0; hltmet100       = 0; hltmet110       = 0; hltmet120       = 0;
  hltmetwithmu90  = 0; hltmetwithmu100 = 0; hltmetwithmu110 = 0; hltmetwithmu120 = 0;
  hltmetwithmu170 = 0; hltmetwithmu300 = 0;
  hltjetmet       = 0; 
  hltphoton90     = 0; hltphoton120    = 0; hltphoton120vbf = 0;  hltphoton165    = 0; hltphoton175    = 0;
  hltdoublemu     = 0;
  hltsinglemu     = 0;
  hltdoubleel     = 0; 
  hltsingleel     = 0;  hltsingleel27   = 0;  hltelnoiso      = 0;
  hltPFHT125      = 0; hltPFHT200 = 0; hltPFHT250 = 0; hltPFHT300 = 0; hltPFHT350 = 0;
  hltPFHT400      = 0; hltPFHT475 = 0; hltPFHT600 = 0; hltPFHT650 = 0; hltPFHT800 = 0; hltPFHT900 = 0; 
  hltEcalHT800    = 0; 
  hltphoton90PFHT = 0;

  trig_obj_n = 0;
  trig_obj_pt.clear();
  trig_obj_eta.clear();
  trig_obj_phi.clear();
  trig_obj_col.clear();

  trig_L1A_check = 0;
  trig_L1A_n     = 0;
  trig_L1A_list.clear();

  trig_L1EG_pt  .clear(); trig_L1EG_eta  .clear(); trig_L1EG_phi  .clear(); 
  trig_L1Jet_pt .clear(); trig_L1Jet_eta .clear(); trig_L1Jet_phi .clear(); 
  trig_L1Mu_pt  .clear(); trig_L1Mu_eta  .clear(); trig_L1Mu_phi  .clear(); 
  trig_L1ETM_pt = trig_L1ETM_phi = trig_L1HTM_pt  = trig_L1HTM_phi = 0; 
  trig_L1ETT_pt = trig_L1ETT_phi = trig_L1HTT_pt  = trig_L1HTT_phi = 0; 
  

}

/////
bool TriggerTreeFiller::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace pat;

  this->initBranches();

  Handle<TriggerResults> triggerResultsH;
  iEvent.getByToken(triggerResultsToken, triggerResultsH);
  Handle<pat::PackedTriggerPrescales> triggerPrescalesH;
  iEvent.getByToken(triggerPrescalesToken, triggerPrescalesH);
    
  Handle<pat::TriggerObjectStandAloneCollection>   triggerObjectsH; 
  edm::Handle<l1t::EGammaBxCollection> H_L1EG;
  edm::Handle<l1t::TauBxCollection>    H_L1Tau;
  edm::Handle<l1t::JetBxCollection>    H_L1Jet;
  edm::Handle<l1t::MuonBxCollection>   H_L1Mu;
  edm::Handle<l1t::EtSumBxCollection>  H_L1Sums;
  edm::Handle<GlobalAlgBlkBxCollection> H_L1Algos;

  if(isTriggerTree and addTriggerObjects){
    iEvent.getByToken(triggerObjectsToken, triggerObjectsH);
    iEvent.getByToken(T_L1EG  , H_L1EG);
    iEvent.getByToken(T_L1Jet , H_L1Jet);
    iEvent.getByToken(T_L1Mu  , H_L1Mu);
    iEvent.getByToken(T_L1Sums, H_L1Sums);
    iEvent.getByToken(T_L1Algos, H_L1Algos);
  }

  const edm::TriggerNames &trignames = iEvent.triggerNames(*triggerResultsH);
  bool triggered = fillTriggerInfo(triggerResultsH,triggerPrescalesH,setHLTFilterFlag,triggerPathsVector,trignames);
  
  if(isTriggerTree and addTriggerObjects) {
    fillTriggerObjects(triggerObjectsH, trignames); //dump HL trigger objects
    fillTriggerL1(H_L1EG, H_L1Tau, H_L1Jet, H_L1Mu, H_L1Sums); //dump L1 objects
    fillAlgosL1(iEvent, iSetup, H_L1Algos); // dump L1 bits
  }
  
  if (applyHLTFilter && !triggered) return false;
  else return true;

}

//////
bool TriggerTreeFiller::fillTriggerInfo(const edm::Handle<edm::TriggerResults> & triggerResultsH,
				       const edm::Handle<pat::PackedTriggerPrescales> & triggerPrescalesH,
				       const bool & setHLTFilterFlag,
				       const std::vector<std::string> & triggerPathsVector,
				       const edm::TriggerNames & trignames){

  
  if(triggerResultsH.isValid() and setHLTFilterFlag == false){
    for (size_t i = 0; i < triggerPathsVector.size(); i++) {

      if (triggerPathsMap[triggerPathsVector[i]] == -1) continue;
      
      if (i == 0  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90        = 1; // MET trigger
      if (i == 1  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90        = 1; // MET trigger
      if (i == 2  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet90        = 1; // MET trigger
      
      if (i == 3  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet100       = 1; // MET trigger
      if (i == 4  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet110       = 1; // MET trigger
      
      if (i == 5  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120       = 1; // MET trigger
      if (i == 6  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120       = 1; // MET trigger
      if (i == 7  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmet120       = 1; // MET trigger
      
      if (i == 8  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu90   = 1; // MET trigger
      if (i == 9  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu100  = 1; // MET trigger
      if (i == 10 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu110  = 1; // MET trigger
      if (i == 11  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu120 = 1; // MET trigger
      
      if (i == 12  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170 = 1; // MET trigger
      if (i == 13  && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170 = 1; // MET trigger
      if (i == 14 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170  = 1; // MET trigger
      if (i == 15 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu170  = 1; // MET trigger
      
      if (i == 16 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu300 = 1; // MET trigger
      if (i == 17 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu300 = 1; // MET trigger
      if (i == 18 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltmetwithmu300 = 1; // MET trigger
      
      if (i == 19 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
      if (i == 20 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
      if (i == 21 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
      
      if (i == 22 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
      if (i == 23 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
      
      if (i == 24 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
      if (i == 25 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
      if (i == 26 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltjetmet    = 1; // Jet-MET trigger
      
      if (i == 27 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton165    = 1; // Photon trigger
      if (i == 28 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton175    = 1; // Photon trigger
      if (i == 29 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton120    = 1; // Photon trigger
      if (i == 30 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton90     = 1; // Photon trigger
      if (i == 31 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton120vbf = 1; // Photon trigger
    
      if (i == 32 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger
      if (i == 33 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
      if (i == 34 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger
      if (i == 35 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoublemu     = 1; // Double muon trigger
      
      if (i == 36 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
      if (i == 37 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
      if (i == 38 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
      if (i == 39 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
      if (i == 40 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
      if (i == 41 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsinglemu     = 1; // Single muon trigger
      
      if (i == 42 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
      if (i == 43 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltdoubleel     = 1; // Double electron trigger
      
      if (i == 44 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
      if (i == 45 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
      if (i == 46 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
      if (i == 46 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel27   = 1; // Single electron trigger
      if (i == 47 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
      if (i == 48 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltsingleel     = 1; // Single electron trigger
      
      if (i == 49 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltelnoiso      = 1; // Single electron trigger
      if (i == 50 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltelnoiso      = 1; // Single electron trigger
      
      if (i == 51 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT125      = 1; // jet ht
      if (i == 52 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT200      = 1; // jet ht
      if (i == 53 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT250      = 1; // jet ht
      if (i == 54 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT300      = 1; // jet ht
      if (i == 55 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT350      = 1; // jet ht
      if (i == 56 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT400      = 1; // jet ht
      if (i == 57 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT475      = 1; // jet ht
      if (i == 58 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT600      = 1; // jet ht
      if (i == 59 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT650      = 1; // jet ht
      if (i == 60 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT800      = 1; // jet ht
      if (i == 61 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltPFHT900      = 1; // jet ht
      if (i == 62 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltEcalHT800    = 1;
      if (i == 63 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton90PFHT = 1;
      if (i == 64 && triggerResultsH->accept(triggerPathsMap[triggerPathsVector[i]])) hltphoton90PFHT = 1;            
    }
  } 
  else if(setHLTFilterFlag == true){
    /// met no-mu
    hltmet90        = 1; hltmet100       = 1; hltmet110       = 1; hltmet120       = 1;
    hltmetwithmu90  = 1; hltmetwithmu100 = 1; hltmetwithmu110 = 1; hltmetwithmu120 = 1;
    hltmetwithmu170 = 1; hltmetwithmu300 = 1;
    hltjetmet       = 1; 
    /// photon
    hltphoton90     = 1; hltphoton120    = 1; hltphoton120vbf = 1; hltphoton165    = 1; hltphoton175    = 1;
    hltdoublemu     = 1; hltsinglemu     = 1;
    hltdoubleel     = 1; hltsingleel     = 1; hltsingleel27   = 1; hltelnoiso      = 1;
    hltPFHT125   = 1; hltPFHT200 = 1; hltPFHT250 = 1; hltPFHT300 = 1; hltPFHT350 = 1; 
    hltPFHT400   = 1; hltPFHT475 = 1; hltPFHT600 = 1; hltPFHT650 = 1; hltPFHT800 = 1; hltPFHT900 = 1;
    hltEcalHT800 = 1; 
    hltphoton90PFHT = 1;
  }
  
  bool triggered = false;
  if (isPhotonPurity == true){ // for photon purity check only photon triggers
    if (hltphoton165    == 1) triggered = true;
    if (hltphoton175    == 1) triggered = true;
    if (hltphoton120    == 1) triggered = true;
    if (hltphoton120vbf == 1) triggered = true;
    if (hltphoton90    == 1) triggered = true;
  }
  else if(isTriggerTree == true){ // trigger study check met-no-mu, lepton triggers, photon triggers and PF-HT
    if (hltmet90        == 1) triggered = true;
    if (hltmet100       == 1) triggered = true;
    if (hltmet110       == 1) triggered = true;
    if (hltmet120       == 1) triggered = true;
    if (hltmetwithmu90  == 1) triggered = true;
    if (hltmetwithmu100 == 1) triggered = true;
    if (hltmetwithmu110 == 1) triggered = true;
    if (hltmetwithmu120 == 1) triggered = true;
    if (hltmetwithmu170 == 1) triggered = true;
    if (hltmetwithmu300 == 1) triggered = true;
    if (hltjetmet       == 1) triggered = true;  
    if (hltphoton165    == 1) triggered = true;
    if (hltphoton175    == 1) triggered = true;
    if (hltphoton120    == 1) triggered = true;
    if (hltphoton120vbf == 1) triggered = true;
    if (hltphoton90     == 1) triggered = true;
    if (hltdoublemu     == 1) triggered = true;
    if (hltsinglemu     == 1) triggered = true;
    if (hltdoubleel     == 1) triggered = true;
    if (hltsingleel     == 1) triggered = true;
    if (hltelnoiso      == 1) triggered = true;
    if (hltPFHT400      == 1) triggered = true;
    if (hltPFHT400      == 1) triggered = true;
    if (hltPFHT400      == 1) triggered = true;  
    if (hltPFHT400      == 1) triggered = true;
    if (hltPFHT475      == 1) triggered = true;
    if (hltPFHT600      == 1) triggered = true;
    if (hltPFHT650      == 1) triggered = true;
    if (hltPFHT800      == 1) triggered = true;
    if (hltPFHT900      == 1) triggered = true;
    if (hltEcalHT800    == 1) triggered = true;
    if (hltphoton90PFHT == 1) triggered = true;
  }
  else if(isQCDTree == true){ // check only PF HT
    if (hltPFHT125      == 1) triggered = true;
    if (hltPFHT200      == 1) triggered = true;
    if (hltPFHT250      == 1) triggered = true;
    if (hltPFHT300      == 1) triggered = true;
    if (hltPFHT350      == 1) triggered = true;
    if (hltPFHT400      == 1) triggered = true;
    if (hltPFHT475      == 1) triggered = true;
    if (hltPFHT600      == 1) triggered = true;
    if (hltPFHT650      == 1) triggered = true;
    if (hltPFHT800      == 1) triggered = true;
    if (hltPFHT900      == 1) triggered = true;
    if (hltEcalHT800    == 1) triggered = true;
  }
  else{ // normal analysis: check photon triggers, met-nu-mu, single-ele, PF-HT
    if (hltmet90        == 1) triggered = true;
    if (hltmet100       == 1) triggered = true;
    if (hltmet110       == 1) triggered = true;
    if (hltmet120       == 1) triggered = true;
    if (hltmetwithmu90  == 1) triggered = true;
    if (hltmetwithmu100 == 1) triggered = true;
    if (hltmetwithmu110 == 1) triggered = true;
    if (hltmetwithmu120 == 1) triggered = true;
    if (hltmetwithmu170 == 1) triggered = true;
    if (hltmetwithmu300 == 1) triggered = true;
    if (hltjetmet       == 1) triggered = true;  
    if (hltphoton165    == 1) triggered = true;
    if (hltphoton175    == 1) triggered = true;
    if (hltphoton120    == 1) triggered = true;
    if (hltphoton120vbf == 1) triggered = true;
    if (hltsinglemu     == 1) triggered = true;
    if (hltsingleel     == 1) triggered = true;
    if (hltelnoiso      == 1) triggered = true;
    if (hltPFHT800      == 1) triggered = true;
    if (hltPFHT900      == 1) triggered = true;
    if (hltEcalHT800    == 1) triggered = true;
  }

  pswgt_ph120 = 1.0;
  pswgt_ph90  = 1.0;
  pswgt_ht125 = 1.0; 
  pswgt_ht200 = 1.0; 
  pswgt_ht250 = 1.0; 
  pswgt_ht300 = 1.0; 
  pswgt_ht350 = 1.0; 
  pswgt_ht400 = 1.0; 
  pswgt_ht475 = 1.0; 
  pswgt_ht600 = 1.0; 
  pswgt_ht650 = 1.0; 
  pswgt_ht800 = 1.0; 
  pswgt_ht900 = 1.0;
  
  for (size_t i = 0; i < triggerResultsH->size(); i++) {
    TString name (trignames.triggerName(i));
    if (trignames.triggerName(i).find("HLT_Photon90_v") != string::npos) pswgt_ph90 = triggerPrescalesH->getPrescaleForIndex(i);
    if (trignames.triggerName(i).find("HLT_Photon120_v") != string::npos) pswgt_ph120 = triggerPrescalesH->getPrescaleForIndex(i);
    if (trignames.triggerName(i).find("HLT_PFHT125_v") != string::npos) pswgt_ht125 = triggerPrescalesH->getPrescaleForIndex(i);
    if (trignames.triggerName(i).find("HLT_PFHT200_v") != string::npos) pswgt_ht200 = triggerPrescalesH->getPrescaleForIndex(i);
    if (trignames.triggerName(i).find("HLT_PFHT250_v") != string::npos) pswgt_ht250 = triggerPrescalesH->getPrescaleForIndex(i);
    if (trignames.triggerName(i).find("HLT_PFHT300_v") != string::npos) pswgt_ht300 = triggerPrescalesH->getPrescaleForIndex(i);
    if (trignames.triggerName(i).find("HLT_PFHT350_v") != string::npos) pswgt_ht350 = triggerPrescalesH->getPrescaleForIndex(i);
    if (trignames.triggerName(i).find("HLT_PFHT400_v") != string::npos) pswgt_ht400 = triggerPrescalesH->getPrescaleForIndex(i);
    if (trignames.triggerName(i).find("HLT_PFHT475_v") != string::npos) pswgt_ht475 = triggerPrescalesH->getPrescaleForIndex(i);
    if (trignames.triggerName(i).find("HLT_PFHT600_v") != string::npos) pswgt_ht600 = triggerPrescalesH->getPrescaleForIndex(i);
    if (trignames.triggerName(i).find("HLT_PFHT650_v") != string::npos) pswgt_ht650 = triggerPrescalesH->getPrescaleForIndex(i);
    if (trignames.triggerName(i).find("HLT_PFHT800_v") != string::npos) pswgt_ht800 = triggerPrescalesH->getPrescaleForIndex(i);
    if (trignames.triggerName(i).find("HLT_PFHT900_v") != string::npos) pswgt_ht900 = triggerPrescalesH->getPrescaleForIndex(i);
  }
  
  return triggered;
 
}	  

////////
void TriggerTreeFiller::fillTriggerObjects(const edm::Handle<pat::TriggerObjectStandAloneCollection> & triggerObjectsH, 
					   const edm::TriggerNames & trignames) {

  std::string trgColl = "";
  int iObj = 0;
  
  if(triggerObjectsH.isValid()) {
    
    // Loop over trigger objects --> dump all the trigger objects
    for (pat::TriggerObjectStandAlone obj : *triggerObjectsH) {       
      iObj++ ;
      obj.unpackPathNames(trignames);
            
      // collection name
      trgColl = obj.collection();
      // dump only e-gamma candidates / met and PFmuons
      if(TString(trgColl).Contains("hltMet") or 
	 TString(trgColl).Contains("hltMET") or 
	 TString(trgColl).Contains("hltMht") or 
	 TString(trgColl).Contains("hltMHT") or 
	 TString(trgColl).Contains("hltPFMet") or
	 TString(trgColl).Contains("hltPFMHT") or
	 TString(trgColl).Contains("hltPFMET") or 
	 TString(trgColl).Contains("hltL2MuonCandidates") or 
	 TString(trgColl).Contains("hltL3MuonCandidates")){

	trig_obj_col.push_back(trgColl);
	trig_obj_pt.push_back( obj.pt());
	trig_obj_eta.push_back(obj.eta());
	trig_obj_phi.push_back(obj.phi());
	trig_obj_n++ ;
      }
    }
  }  
}

/////////
void TriggerTreeFiller::fillAlgosL1(const edm::Event & iEvent,
				    const edm::EventSetup & eventSetup,
				    const edm::Handle<GlobalAlgBlkBxCollection> & H_L1Algos){
  
  // Check Handle validity
  if(!H_L1Algos.isValid()) {
    trig_L1A_check = -1;
    return;
  }

  // Get L1 Menu and utils
  edm::ESHandle<L1TUtmTriggerMenu> menu;
  eventSetup.get<L1TUtmTriggerMenuRcd>().get(menu);
  int iErrorCode = -1;
  L1GtUtils::TriggerCategory trigCategory = L1GtUtils::AlgorithmTrigger;
  // Get L1 utils for prescales
  L1GtUtils const & l1GtUtils = hltPrescaleProvider->l1GtUtils();
  l1GtUtils.prescaleFactorSetIndex(iEvent, trigCategory, iErrorCode);
  // Get the bit/name association //
  const UInt_t nBits = 512;
  std::string algoBitToName[nBits];
  //
  for (auto const & keyval: menu->getAlgorithmMap()) { 
    std::string const & trigName  = keyval.second.getName(); 
    unsigned int index = keyval.second.getIndex(); 
    if(index < nBits) algoBitToName[index] = trigName;
  } // end algo Map
  
  // Get the L1 decision per algo //
  GlobalAlgBlk const &result = H_L1Algos->at(0,0);
  //
  for (unsigned int itrig = 0; itrig < result.maxPhysicsTriggers; ++itrig) {
    // Check decision for this bit
    bool myflag = result.getAlgoDecisionFinal(itrig) ; 
    if(myflag ) { 
      if(TString(algoBitToName[itrig]).Contains("ETM")){
	trig_L1A_list.push_back(algoBitToName[itrig]);
	trig_L1A_n++ ;
      }
    }
  } // end loop: L1 trigger results  
}

void TriggerTreeFiller::fillTriggerL1(const edm::Handle<l1t::EGammaBxCollection> & H_L1EG,  const edm::Handle<l1t::TauBxCollection>  & H_L1Tau,
				      const edm::Handle<l1t::JetBxCollection>    & H_L1Jet, const edm::Handle<l1t::MuonBxCollection> & H_L1Mu,
				      const edm::Handle<l1t::EtSumBxCollection>  & H_L1Sums) {

  int sumType   = -1;
  int bunchCrossing = 0;
  
  // L1 EG    
  if(H_L1EG.isValid()) {
    for (l1t::EGammaBxCollection::const_iterator it=H_L1EG->begin(bunchCrossing); it!=H_L1EG->end(bunchCrossing); it++){
      if(it->pt() < minL1EG) continue;
      trig_L1EG_pt .push_back( it->pt()  );
      trig_L1EG_eta.push_back( it->eta() );
      trig_L1EG_phi.push_back( it->phi() );
    }
  }
  
  // L1 Jet
  if(H_L1Jet.isValid()) {
    for (l1t::JetBxCollection::const_iterator it=H_L1Jet->begin(bunchCrossing); it!=H_L1Jet->end(bunchCrossing); it++){
      if(it->pt() < minL1Jet) continue;
      trig_L1Jet_pt .push_back( it->pt()  );
      trig_L1Jet_eta.push_back( it->eta() );
      trig_L1Jet_phi.push_back( it->phi() );
    }
  }

  // L1 Mu
  if(H_L1Mu.isValid()) {
    for (l1t::MuonBxCollection::const_iterator it=H_L1Mu->begin(bunchCrossing); it!=H_L1Mu->end(bunchCrossing); it++){
      if(it->pt() < minL1Mu) continue;
      trig_L1Mu_pt .push_back( it->pt()  );
      trig_L1Mu_eta.push_back( it->eta() );
      trig_L1Mu_phi.push_back( it->phi() );
    }
  }
  
  // L1 Sums
  if(H_L1Sums.isValid()) {
    for (l1t::EtSumBxCollection::const_iterator it=H_L1Sums->begin(bunchCrossing); it!=H_L1Sums->end(bunchCrossing); it++){
      
      sumType = static_cast<int>( it->getType() );
      if(sumType == l1t::EtSum::kTotalEt){
	trig_L1ETT_pt  = it->et();
	trig_L1ETT_phi = it->phi();
      }
      else if(sumType == l1t::EtSum::kTotalHt){
	trig_L1HTT_pt  = it->et();
	trig_L1HTT_phi = it->phi();
      }
      else if(sumType == l1t::EtSum::kMissingEt){
	trig_L1ETM_pt  = it->et();
	trig_L1ETM_phi = it->phi();
      }
      else if(sumType == l1t::EtSum::kMissingHt){
	trig_L1HTM_pt  = it->et();
	trig_L1HTM_phi = it->phi();
      }      
    }    
  }
}

	  
void TriggerTreeFiller::SetTriggerPaths(edm::Run const& iRun, edm::EventSetup const& iSetup){

  // triggers for the Analysis
  triggerPathsVector.push_back("HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight"); //0
  triggerPathsVector.push_back("HLT_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight"); //1
  triggerPathsVector.push_back("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight"); //2
  triggerPathsVector.push_back("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v");//3
  triggerPathsVector.push_back("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight_v");//4
  triggerPathsVector.push_back("HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight"); //5
  triggerPathsVector.push_back("HLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight"); //6
  triggerPathsVector.push_back("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight"); //7
  triggerPathsVector.push_back("HLT_PFMET90_PFMHT90_IDTight"); //8
  triggerPathsVector.push_back("HLT_PFMET100_PFMHT100_IDTight_v"); //9 
  triggerPathsVector.push_back("HLT_PFMET110_PFMHT110_IDTight_v");//10 
  triggerPathsVector.push_back("HLT_PFMET120_PFMHT120_IDTight"); //11
  triggerPathsVector.push_back("HLT_PFMET170_NoiseCleaned"); //12
  triggerPathsVector.push_back("HLT_PFMET170_JetIdCleaned"); //13
  triggerPathsVector.push_back("HLT_PFMET170_HBHECleaned"); //14
  triggerPathsVector.push_back("HLT_PFMET170_v"); //15
  triggerPathsVector.push_back("HLT_PFMET300_NoiseCleaned"); //16
  triggerPathsVector.push_back("HLT_PFMET300_JetIdCleaned"); //17
  triggerPathsVector.push_back("HLT_PFMET300_v"); //18
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight"); //19
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_JetIdCleaned_PFMHTNoMu90_IDTight"); //20
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight");    //21
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight_v");//22
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight_v");//23
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight"); //24
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight"); //25
  triggerPathsVector.push_back("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight"); //26
  triggerPathsVector.push_back("HLT_Photon165_HE10"); //27
  triggerPathsVector.push_back("HLT_Photon175");      //28
  triggerPathsVector.push_back("HLT_Photon120_v");    //29
  triggerPathsVector.push_back("HLT_Photon90_v");     //30
  triggerPathsVector.push_back("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40_v");     //31
  triggerPathsVector.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v"); //32
  triggerPathsVector.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v"); //33
  triggerPathsVector.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v"); //34
  triggerPathsVector.push_back("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"); //35
  triggerPathsVector.push_back("HLT_IsoMu20_v"); //36
  triggerPathsVector.push_back("HLT_IsoMu22_v"); //37
  triggerPathsVector.push_back("HLT_IsoMu24_v"); //38
  triggerPathsVector.push_back("HLT_IsoTkMu20"); //39
  triggerPathsVector.push_back("HLT_IsoTkMu22"); //40
  triggerPathsVector.push_back("HLT_IsoTkMu24"); //41
  triggerPathsVector.push_back("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"); //42
  triggerPathsVector.push_back("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"); //43
  triggerPathsVector.push_back("HLT_Ele24_eta2p1_WPLoose_Gsf_v"); //44
  triggerPathsVector.push_back("HLT_Ele25_eta2p1_WPTight_Gsf_v"); //45
  triggerPathsVector.push_back("HLT_Ele27_WPTight_Gsf_v"); //46
  triggerPathsVector.push_back("HLT_Ele27_eta2p1_WPLoose_Gsf_v"); //47 
  triggerPathsVector.push_back("HLT_Ele27_eta2p1_WPTight_Gsf_v"); //48
  triggerPathsVector.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v"); //49
  triggerPathsVector.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v"); //50
  triggerPathsVector.push_back("HLT_PFHT125_v");//51
  triggerPathsVector.push_back("HLT_PFHT200_v");//52
  triggerPathsVector.push_back("HLT_PFHT250_v");//53
  triggerPathsVector.push_back("HLT_PFHT300_v");//54
  triggerPathsVector.push_back("HLT_PFHT350_v");//55  
  triggerPathsVector.push_back("HLT_PFHT400_v");//56
  triggerPathsVector.push_back("HLT_PFHT475_v");//57
  triggerPathsVector.push_back("HLT_PFHT600_v");//58
  triggerPathsVector.push_back("HLT_PFHT650_v");//59
  triggerPathsVector.push_back("HLT_PFHT800_v");//60
  triggerPathsVector.push_back("HLT_PFHT900_v");//61
  triggerPathsVector.push_back("HLT_ECALHT800_v");//62
  triggerPathsVector.push_back("HLT_Photon90_CaloIdL_PFHT500_v");//63
  triggerPathsVector.push_back("HLT_Photon90_CaloIdL_PFHT600_v");//64

  HLTConfigProvider hltConfig;
  bool changedConfig = false;
  hltConfig.init(iRun, iSetup, triggerResultsTag.process(), changedConfig);
  
  for (size_t i = 0; i < triggerPathsVector.size(); i++) {
    triggerPathsMap[triggerPathsVector[i]] = -1;
  }
  
  for(size_t i = 0; i < triggerPathsVector.size(); i++){
    TPRegexp pattern(triggerPathsVector[i]);
    for(size_t j = 0; j < hltConfig.triggerNames().size(); j++){
      std::string pathName = hltConfig.triggerNames()[j];
      if(TString(pathName).Contains(pattern)){
	triggerPathsMap[triggerPathsVector[i]] = j;
      }
    }
  }
  
  bool changedHLTPSP = false;
  hltPrescaleProvider->init(iRun, iSetup, triggerResultsTag.process(), changedHLTPSP); //ND 

}

void TriggerTreeFiller::DeclareAndSetBranches(){
  
  tree_->Branch("hltmet90"             , &hltmet90             , "hltmet90/b");
  tree_->Branch("hltmet100"            , &hltmet100            , "hltmet100/b");
  tree_->Branch("hltmet110"            , &hltmet110            , "hltmet110/b");
  tree_->Branch("hltmet120"            , &hltmet120            , "hltmet120/b");
  tree_->Branch("hltmetwithmu90"       , &hltmetwithmu90       , "hltmetwithmu90/b");
  tree_->Branch("hltmetwithmu100"      , &hltmetwithmu100      , "hltmetwithmu100/b");
  tree_->Branch("hltmetwithmu110"      , &hltmetwithmu110      , "hltmetwithmu110/b");
  tree_->Branch("hltmetwithmu120"      , &hltmetwithmu120      , "hltmetwithmu120/b");
  tree_->Branch("hltmetwithmu170"      , &hltmetwithmu170      , "hltmetwithmu170/b");
  tree_->Branch("hltmetwithmu300"      , &hltmetwithmu300      , "hltmetwithmu300/b");
  tree_->Branch("hltjetmet"            , &hltjetmet            , "hltjetmet/b");

  if(not isQCDTree){
    tree_->Branch("hltphoton90"          , &hltphoton90          , "hltphoton90/b");
    tree_->Branch("hltphoton120"         , &hltphoton120         , "hltphoton120/b");
    tree_->Branch("hltphoton120vbf"      , &hltphoton120vbf      , "hltphoton120vbf/b");
    tree_->Branch("hltphoton165"         , &hltphoton165         , "hltphoton165/b");
    tree_->Branch("hltphoton175"         , &hltphoton175         , "hltphoton175/b");
    tree_->Branch("hltphoton90PFHT"      , &hltphoton90PFHT      , "hltphoton90PFHT/b");
    tree_->Branch("hltEcalHT800"         , &hltEcalHT800         , "hltEcalHT800/b");  

    tree_->Branch("hltdoublemu"          , &hltdoublemu          , "hltdoublemu/b");
    tree_->Branch("hltsinglemu"          , &hltsinglemu          , "hltsinglemu/b");
    tree_->Branch("hltdoubleel"          , &hltdoubleel          , "hltdoubleel/b");
    tree_->Branch("hltsingleel"          , &hltsingleel          , "hltsingleel/b");
    tree_->Branch("hltsingleel27"        , &hltsingleel27        , "hltsingleel27/b");
    tree_->Branch("hltelnoiso"           , &hltelnoiso           , "hltelnoiso/b");
  }

  if(not isPhotonPurity){
    if(isQCDTree){
      tree_->Branch("hltPFHT125"           , &hltPFHT125           , "hltPFHT125/b");
      tree_->Branch("hltPFHT200"           , &hltPFHT200           , "hltPFHT200/b");
      tree_->Branch("hltPFHT250"           , &hltPFHT250           , "hltPFHT250/b");
      tree_->Branch("hltPFHT300"           , &hltPFHT300           , "hltPFHT300/b");
      tree_->Branch("hltPFHT350"           , &hltPFHT350           , "hltPFHT350/b");
    }
    
    if(not isTriggerTree or (isTriggerTree and addNonETMTriggerObjects)){
      tree_->Branch("hltPFHT400"           , &hltPFHT400           , "hltPFHT400/b");
      tree_->Branch("hltPFHT475"           , &hltPFHT475           , "hltPFHT475/b");
      tree_->Branch("hltPFHT600"           , &hltPFHT600           , "hltPFHT600/b");
      tree_->Branch("hltPFHT650"           , &hltPFHT650           , "hltPFHT650/b");
      tree_->Branch("hltPFHT800"           , &hltPFHT800           , "hltPFHT800/b");
      tree_->Branch("hltPFHT900"           , &hltPFHT900           , "hltPFHT900/b");
      tree_->Branch("pswgt_ph120"          , &pswgt_ph120          , "pswgt_ph120/F");
      tree_->Branch("pswgt_ph90"           , &pswgt_ph90           , "pswgt_ph90/F");
    }
    
    if((isTriggerTree and not isMC and addNonETMTriggerObjects) or isQCDTree){
      tree_->Branch("pswgt_ht125"          , &pswgt_ht125          , "pswgt_ht125/F");
      tree_->Branch("pswgt_ht200"          , &pswgt_ht200          , "pswgt_ht200/F");
      tree_->Branch("pswgt_ht250"          , &pswgt_ht250          , "pswgt_ht250/F");
      tree_->Branch("pswgt_ht300"          , &pswgt_ht300          , "pswgt_ht300/F");
      tree_->Branch("pswgt_ht350"          , &pswgt_ht350          , "pswgt_ht350/F");
      tree_->Branch("pswgt_ht400"          , &pswgt_ht400          , "pswgt_ht400/F");
      tree_->Branch("pswgt_ht475"          , &pswgt_ht475          , "pswgt_ht475/F");
      tree_->Branch("pswgt_ht600"          , &pswgt_ht600          , "pswgt_ht600/F");
      tree_->Branch("pswgt_ht650"          , &pswgt_ht650          , "pswgt_ht650/F");
      tree_->Branch("pswgt_ht800"          , &pswgt_ht800          , "pswgt_ht800/F");
      tree_->Branch("pswgt_ht900"          , &pswgt_ht900          , "pswgt_ht900/F");
    }

    if(isTriggerTree and addTriggerObjects){

      tree_->Branch("trig_obj_n"           , &trig_obj_n           , "trig_obj_n/I"); 
      tree_->Branch("trig_obj_pt"          , "std::vector<float>" , &trig_obj_pt);
      tree_->Branch("trig_obj_eta"         , "std::vector<float>" , &trig_obj_eta);
      tree_->Branch("trig_obj_phi"         , "std::vector<float>" , &trig_obj_phi);
      tree_->Branch("trig_obj_col"         , "std::vector<std::string>" , &trig_obj_col);

      if(addNonETMTriggerObjects){
	tree_->Branch("trig_L1A_check"       , &trig_L1A_check            , "trig_L1A_check/I");
	tree_->Branch("trig_L1A_n"           , &trig_L1A_n                , "trig_L1A_n/I");
      }
      tree_->Branch("trig_L1A_list"        , "std::vector<std::string>" , &trig_L1A_list);
      
      if(addNonETMTriggerObjects){ // EG and Jets
	tree_->Branch("trig_L1EG_pt"         , "std::vector<float>" , &trig_L1EG_pt);
	tree_->Branch("trig_L1EG_eta"        , "std::vector<float>" , &trig_L1EG_eta);
	tree_->Branch("trig_L1EG_phi"        , "std::vector<float>" , &trig_L1EG_phi);
	tree_->Branch("trig_L1Jet_pt"        , "std::vector<float>" , &trig_L1Jet_pt);
	tree_->Branch("trig_L1Jet_eta"       , "std::vector<float>" , &trig_L1Jet_eta);
	tree_->Branch("trig_L1Jet_phi"       , "std::vector<float>" , &trig_L1Jet_phi);
      }
   

      tree_->Branch("trig_L1Mu_pt"         , "std::vector<float>" , &trig_L1Mu_pt);
      tree_->Branch("trig_L1Mu_eta"        , "std::vector<float>" , &trig_L1Mu_eta);
      tree_->Branch("trig_L1Mu_phi"        , "std::vector<float>" , &trig_L1Mu_phi);

      tree_->Branch("trig_L1ETM_pt"        , &trig_L1ETM_pt        , "trig_L1ETM_pt/F");
      tree_->Branch("trig_L1ETM_phi"       , &trig_L1ETM_phi       , "trig_L1ETM_phi/F");

      if(addNonETMTriggerObjects){
	tree_->Branch("trig_L1ETT_pt"        , &trig_L1ETT_pt        , "trig_L1ETT_pt/F");
	tree_->Branch("trig_L1ETT_phi"       , &trig_L1ETT_phi       , "trig_L1ETT_phi/F");
      }

      tree_->Branch("trig_L1HTM_pt"        , &trig_L1HTM_pt        , "trig_L1HTM_pt/F");
      tree_->Branch("trig_L1HTM_phi"       , &trig_L1HTM_phi       , "trig_L1HTM_phi/F");
      
      if(addNonETMTriggerObjects){      
	tree_->Branch("trig_L1HTT_pt"        , &trig_L1HTT_pt        , "trig_L1HTT_pt/F");
	tree_->Branch("trig_L1HTT_phi"       , &trig_L1HTT_phi       , "trig_L1HTT_phi/F");
      }
    }
  }  
}

