#include "AnalysisCode/MonoXAnalysis/interface/VJetTreeFiller.h"

VJetTreeFiller::VJetTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree, const bool & isPuppi):
  muonsTag       (iConfig.getParameter<edm::InputTag>("muons")),
  electronsTag   (iConfig.getParameter<edm::InputTag>("electrons")),
  photonsTag     (iConfig.getParameter<edm::InputTag>("photons")),
  jetsTag        (iConfig.getParameter<edm::InputTag>("boostedJetsCHS")),
  isMC           (iConfig.existsAs<bool>("isMC") ? iConfig.getParameter<bool>("isMC") : false),
  isTriggerTree  (iConfig.existsAs<bool>("isTriggerTree") ? iConfig.getParameter<bool>("isTriggerTree") : false),
  isQCDTree      (iConfig.existsAs<bool>("isQCDTree") ? iConfig.getParameter<bool>("isQCDTree") : false),
  isPhotonPurity (iConfig.existsAs<bool>("isPhotonPurity") ? iConfig.getParameter<bool>("isPhotonPurity") : false),
  applyDiMuonFilter (iConfig.existsAs<bool>("applyDiMuonFilter") ? iConfig.getParameter<bool>("applyDiMuonFilter") : false),
  applyDiElectronFilter (iConfig.existsAs<bool>("applyDiElectronFilter") ? iConfig.getParameter<bool>("applyDiElectronFilter") : false),
  applyPhotonJetsFilter (iConfig.existsAs<bool>("applyPhotonJetsFilter") ? iConfig.getParameter<bool>("applyPhotonJetsFilter") : false),
  cleanMuonJet   (iConfig.existsAs<bool>("cleanMuonJet") ? iConfig.getParameter<bool>("cleanMuonJet") : false),
  cleanElectronJet (iConfig.existsAs<bool>("cleanElectronJet") ? iConfig.getParameter<bool>("cleanElectronJet") : false),
  cleanPhotonJet (iConfig.existsAs<bool>("cleanPhotonJet") ? iConfig.getParameter<bool>("cleanPhotonJet") : false),
  dRCleaningAK8  (iConfig.existsAs<double>("dRCleaningAK8") ? iConfig.getParameter<double>("dRCleaningAK8") : 0.4),
  jetidwp        (iConfig.existsAs<std::string>("jetidwp") ? iConfig.getParameter<std::string>("jetidwp") : "loose"),
  btaggingCSVWP  (iConfig.getParameter<double>("btaggingCSVWP")),
  addBTagScaleFactor(iConfig.existsAs<bool>("addBTagScaleFactor") ? iConfig.getParameter<bool>("addBTagScaleFactor") : false),
  useMiniAODSubstructure(iConfig.existsAs<bool>("useMiniAODSubstructure") ? iConfig.getParameter<bool>("useMiniAODSubstructure") : true){
  
  isPuppi_ = isPuppi;

  if(not isPuppi_)
    jetsLabel = TString::Format(iConfig.getParameter<edm::InputTag>("boostedJetsCHS").label().c_str());
  else{
    jetsTag   = iConfig.getParameter<edm::InputTag>("boostedJetsPuppi");
    jetsLabel = TString::Format(iConfig.getParameter<edm::InputTag>("boostedJetsPuppi").label().c_str());
  }
  
  jetsLabel.ReplaceAll("packed","");
  jetsLabel.ReplaceAll("PatJets","");
  jetsLabel.ReplaceAll("patJets","");

  jetsToken      = iC.consumes<std::vector<pat::Jet> > (jetsTag);  
  muonsToken     = iC.consumes<pat::MuonRefVector> (muonsTag);
  electronsToken = iC.consumes<pat::ElectronRefVector> (electronsTag);
  photonsToken   = iC.consumes<pat::PhotonRefVector> (photonsTag);

  //////
  if(addBTagScaleFactor and isMC){
  
    bTagScaleFactorFileSubCSV = iConfig.getParameter<edm::FileInPath>("bTagScaleFactorFileSubCSV");
    if ( bTagScaleFactorFileSubCSV.location()!=edm::FileInPath::Local)
      throw cms::Exception("VJetTreeFiller") << " Failed to find File = " << bTagScaleFactorFileSubCSV << " !!\n";
      
    calibSubCSV = BTagCalibration("CSVv2",bTagScaleFactorFileSubCSV.fullPath());
    bMediumSubCSV.push_back(BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"})); // for light flavor
    bMediumSubCSV.back().load(calibSubCSV,BTagEntry::FLAV_B,"lt");
    bMediumSubCSV.back().load(calibSubCSV,BTagEntry::FLAV_C,"lt");
    bMediumSubCSV.back().load(calibSubCSV,BTagEntry::FLAV_UDSG,"incl");      
  }

  tree_ = tree;
  this->DeclareAndSetBranches();
  this->initBranches();
}

/////
void VJetTreeFiller::initBranches(){

  boostedJetpt    .clear(); boostedJeteta  .clear();  
  boostedJetphi   .clear(); boostedJetm    .clear();
  boostedJetGenpt .clear(); boostedJetGenm .clear();
  boostedJetGeneta.clear(); boostedJetGenphi .clear();      
  boostedJettau1  .clear(); boostedJettau2 .clear();
  boostedJettau3  .clear(); boostedJettau4 .clear();
  boostedJetGentau1.clear(); boostedJetGentau2 .clear();
  boostedJetGentau3.clear(); boostedJetGentau4 .clear();
  boostedJetHFlav.clear(); boostedJetPFlav.clear(); boostedJetQGL.clear(); 
  boostedJetBtag .clear(), boostedJetDoubleBtag .clear();
      
  prunedJetpt     .clear(); prunedJetm    .clear(); prunedJetGenpt.clear(); prunedJetGenm .clear();
  prunedJeteta    .clear(); prunedJetphi  .clear(); prunedJetGeneta.clear(); prunedJetGenphi.clear();
  prunedJetm_v2 .clear();   prunedJetpt_v2 .clear(); prunedJeteta_v2 .clear(); prunedJetphi_v2 .clear(); 
  prunedJetHFlav  .clear(); prunedJetPFlav.clear(); prunedJetQGL  .clear(); prunedJetBtag .clear(); prunedJetDoubleBtag.clear();
  prunedJetptraw  .clear(); prunedJetmraw .clear();
      
  prunedSubJetpt_1 .clear(); prunedSubJetm_1  .clear(); prunedSubJetphi_1 .clear(); prunedSubJeteta_1 .clear();
  prunedSubJetHFlav_1 .clear(); prunedSubJetQGL_1 .clear(); prunedSubJetBtag_1 .clear();
  prunedSubJetBtagSF_1.clear(); prunedSubJetBtagSFUp_1.clear(); prunedSubJetBtagSFDown_1.clear();
  prunedSubJetGenpt_1 .clear(); prunedSubJetGenm_1 .clear(); prunedSubJetPFlav_1 .clear();
  prunedSubJetGenphi_1 .clear(); prunedSubJetGeneta_1 .clear(); 
  prunedSubJetptraw_1 .clear(); prunedSubJetmraw_1  .clear(); 
      
  prunedSubJetpt_2 .clear(); prunedSubJetm_2  .clear(); prunedSubJetphi_2 .clear(); prunedSubJeteta_2 .clear();
  prunedSubJetHFlav_2 .clear(); prunedSubJetQGL_2 .clear(); prunedSubJetBtag_2 .clear();  prunedSubJetPFlav_2 .clear();
  prunedSubJetBtagSF_2.clear(); prunedSubJetBtagSFUp_2.clear(); prunedSubJetBtagSFDown_2.clear();
  prunedSubJetGenpt_2 .clear(); prunedSubJetGenm_2 .clear();
  prunedSubJetGeneta_2 .clear(); prunedSubJetGenphi_2 .clear();
  prunedSubJetptraw_2 .clear(); prunedSubJetmraw_2  .clear(); 

  // Puppi jets
  boostedPuppiJetpt    .clear(); boostedPuppiJeteta  .clear();  boostedPuppiJetphi   .clear(); boostedPuppiJetm    .clear();
  boostedPuppiJetGenpt .clear(); boostedPuppiJetGenm .clear();  boostedPuppiJetGeneta .clear(); boostedPuppiJetGenphi .clear();
  boostedPuppiJettau1  .clear(); boostedPuppiJettau2 .clear();  boostedPuppiJettau3  .clear(); boostedPuppiJettau4 .clear();
  boostedPuppiJetGentau1  .clear(); boostedPuppiJetGentau2 .clear(); boostedPuppiJetGentau3  .clear(); boostedPuppiJetGentau4 .clear();
  boostedPuppiJetHFlav .clear(); boostedPuppiJetPFlav .clear(); boostedPuppiJetQGL  .clear(); 
  boostedPuppiJetBtag .clear(); boostedPuppiJetDoubleBtag .clear();
  
  softDropPuppiJetpt    .clear(); softDropPuppiJetm .clear(); softDropPuppiJetGenpt .clear(); softDropPuppiJetGenm .clear(); 
  softDropPuppiJetm_v2.clear();   softDropPuppiJetpt_v2.clear(); softDropPuppiJeteta_v2.clear(); softDropPuppiJetphi_v2.clear();
  softDropPuppiJeteta    .clear(); softDropPuppiJetphi .clear(); softDropPuppiJetGeneta .clear(); softDropPuppiJetGenphi .clear(); 
  softDropPuppiJetHFlav .clear(); softDropPuppiJetPFlav .clear(); softDropPuppiJetQGL .clear(); softDropPuppiJetBtag .clear();
  softDropPuppiJetDoubleBtag .clear(); softDropPuppiJetmraw .clear(); softDropPuppiJetptraw .clear();
      
  softDropPuppiSubJetpt_1 .clear(); softDropPuppiSubJetm_1  .clear(); softDropPuppiSubJetphi_1 .clear(); softDropPuppiSubJeteta_1 .clear();
  softDropPuppiSubJetHFlav_1 .clear(); softDropPuppiSubJetQGL_1 .clear(); softDropPuppiSubJetBtag_1 .clear(); softDropPuppiSubJetPFlav_1 .clear();
  softDropPuppiSubJetGenpt_1 .clear(); softDropPuppiSubJetGenm_1 .clear(); softDropPuppiSubJetGenphi_1 .clear(); softDropPuppiSubJetGeneta_1 .clear();
  softDropPuppiSubJetptraw_1 .clear(); softDropPuppiSubJetmraw_1  .clear(); 
  softDropPuppiSubJetBtagSF_1.clear(); softDropPuppiSubJetBtagSFUp_1.clear(); softDropPuppiSubJetBtagSFDown_1.clear();

  softDropPuppiSubJetpt_2 .clear(); softDropPuppiSubJetm_2  .clear(); softDropPuppiSubJetphi_2 .clear(); softDropPuppiSubJeteta_2 .clear();
  softDropPuppiSubJetHFlav_2 .clear(); softDropPuppiSubJetQGL_2 .clear(); softDropPuppiSubJetBtag_2 .clear(); softDropPuppiSubJetPFlav_2 .clear();
  softDropPuppiSubJetGenpt_2 .clear(); softDropPuppiSubJetGenm_2 .clear(); softDropPuppiSubJetGenphi_2 .clear(); softDropPuppiSubJetGeneta_2 .clear();
  softDropPuppiSubJetptraw_2 .clear(); softDropPuppiSubJetmraw_2  .clear(); 
  softDropPuppiSubJetBtagSF_2.clear(); softDropPuppiSubJetBtagSFUp_2.clear(); softDropPuppiSubJetBtagSFDown_2.clear();
       
}

/////
bool VJetTreeFiller::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  this->initBranches();

  if(isTriggerTree or isPhotonPurity or applyDiMuonFilter or applyDiElectronFilter or applyPhotonJetsFilter) return true;
  
  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace pat;

  Handle<pat::MuonRefVector> muonsH;
  iEvent.getByToken(muonsToken, muonsH);
  pat::MuonRefVector muons = *muonsH;

  Handle<pat::ElectronRefVector> electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  pat::ElectronRefVector electrons = *electronsH;

  Handle<pat::PhotonRefVector> photonsH;
  iEvent.getByToken(photonsToken, photonsH);
  pat::PhotonRefVector photons = *photonsH;

  Handle<vector<pat::Jet> > jetsH;
  iEvent.getByToken(jetsToken, jetsH);

  vector<pat::JetRef> jetsBoosted;
  
  //sort collection to make sure it is ordered
  if(jetsH.isValid())
    fillJetCollections(jetsH,muons,electrons,photons,jetsBoosted,isPuppi_);

  // Loop on the selected jets
  for(size_t i = 0; i < jetsBoosted.size(); i++){
    
    if(not isPuppi_){ // fill regular jets

      boostedJetpt  .push_back( jetsBoosted[i]->pt());
      boostedJeteta .push_back( jetsBoosted[i]->eta());
      boostedJetphi .push_back( jetsBoosted[i]->phi());
      boostedJetm   .push_back( jetsBoosted[i]->mass());
      
      // gen jets
      if(isMC){
	if(jetsBoosted[i]->genJet()){ // gen AK8 jet
	  boostedJetGenpt .push_back( jetsBoosted[i]->genJet()->pt());
	  boostedJetGenm  .push_back( jetsBoosted[i]->genJet()->mass());
	  boostedJetGeneta  .push_back( jetsBoosted[i]->genJet()->eta());
	  boostedJetGenphi  .push_back( jetsBoosted[i]->genJet()->phi());      
	}  
	else{
	  boostedJetGenpt .push_back(0);
	  boostedJetGenm  .push_back(0);
	  boostedJetGeneta  .push_back(0);
	  boostedJetGenphi  .push_back(0);      
	}
	boostedJetHFlav .push_back( jetsBoosted[i]->hadronFlavour());
	boostedJetPFlav .push_back( jetsBoosted[i]->partonFlavour());  
      }
      
      if(not useMiniAODSubstructure){

	// b-tagging information
	boostedJetBtag .push_back( jetsBoosted[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	boostedJetDoubleBtag .push_back( jetsBoosted[i]->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));
	
	// N-jettiness
	if(jetsBoosted[i]->hasUserFloat("Njettiness"+jetsLabel+":tau1"))
	  boostedJettau1 .push_back( jetsBoosted[i]->userFloat("Njettiness"+jetsLabel+":tau1"));

	if(jetsBoosted[i]->hasUserFloat("Njettiness"+jetsLabel+":tau2"))
	  boostedJettau2 .push_back( jetsBoosted[i]->userFloat("Njettiness"+jetsLabel+":tau2"));
	  
	if(jetsBoosted[i]->hasUserFloat("Njettiness"+jetsLabel+":tau3"))
	  boostedJettau3 .push_back( jetsBoosted[i]->userFloat("Njettiness"+jetsLabel+":tau3"));
	  
	if(jetsBoosted[i]->hasUserFloat("Njettiness"+jetsLabel+":tau4"))
	  boostedJettau4 .push_back( jetsBoosted[i]->userFloat("Njettiness"+jetsLabel+":tau4"));

	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"QGL:qgLikelihood"))
	  boostedJetQGL .push_back( jetsBoosted[i]->userFloat(jetsLabel+"QGL:qgLikelihood"));
	  
	// Gen n-subjettiness
	if(isMC){
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"GenNjettinessMatched:tau1"))
	    boostedJetGentau1 .push_back( jetsBoosted[i]->userFloat(jetsLabel+"GenNjettinessMatched:tau1"));
	      
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"GenNjettinessMatched:tau2"))
	    boostedJetGentau2 .push_back( jetsBoosted[i]->userFloat(jetsLabel+"GenNjettinessMatched:tau2"));
	      
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"GenNjettinessMatched:tau3"))
	    boostedJetGentau3 .push_back( jetsBoosted[i]->userFloat(jetsLabel+"GenNjettinessMatched:tau3"));
	      
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"GenNjettinessMatched:tau4"))
	    boostedJetGentau4 .push_back( jetsBoosted[i]->userFloat(jetsLabel+"GenNjettinessMatched:tau4"));
	}	
      }
      else{

	// N-jettiness
	if(jetsBoosted[i]->hasUserFloat("NjettinessAK8:tau1"))
	  boostedJettau1 .push_back( jetsBoosted[i]->userFloat("NjettinessAK8:tau1"));
	  
	if(jetsBoosted[i]->hasUserFloat("NjettinessAK8:tau2"))
	  boostedJettau2 .push_back( jetsBoosted[i]->userFloat("NjettinessAK8:tau2"));
	  
	if(jetsBoosted[i]->hasUserFloat("NjettinessAK8:tau3"))
	  boostedJettau3 .push_back(jetsBoosted[i]->userFloat("NjettinessAK8:tau3"));

	if(jetsBoosted[i]->hasUserFloat("NjettinessAK8:tau4"))
	  boostedJettau3 .push_back(jetsBoosted[i]->userFloat("NjettinessAK8:tau4"));

      }

      // Info about pruning
      if(not useMiniAODSubstructure){

	// pruned matched jet
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:mass"))
	  prunedJetm .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:mass"));
	  
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:eta"))
	  prunedJeteta .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:eta"));
	
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:phi"))
	  prunedJetphi .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:phi"));
	  
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:pt"))
	  prunedJetpt .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:pt"));
	  
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedQGLMatched:qgLikelihood"))
	  prunedJetQGL .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedQGLMatched:qgLikelihood"));
	  
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"))
	  prunedJetBtag .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"))
	  prunedJetDoubleBtag .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"));
	  
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:rawpt"))
	  prunedJetptraw .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:rawpt"));
	  
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:rawmass")){
	  prunedJetmraw .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:rawmass"));

	  // apply correction by hand from uncorrected variables
	  if(jetsBoosted[i]->availableJECSets().size() > 1 and 
	     jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:raweta") and
	     jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:rawphi") and
	     jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:rawpt")){
	          
	    TLorentzVector correctedP4;
	    correctedP4.SetPtEtaPhiM(jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:rawpt"),
				     jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:raweta"),
				     jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:rawphi"),
				     jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:rawmass")
				     );
	    correctedP4 *= 1./jetsBoosted[i]->jecFactor("Uncorrected","none",jetsBoosted[i]->availableJECSets().at(1)); // apply AK8 corrections
	    prunedJetm_v2 .push_back(correctedP4.M());
	    prunedJetpt_v2 .push_back(correctedP4.Pt());
	    prunedJeteta_v2 .push_back(correctedP4.Eta());
	    prunedJetphi_v2 .push_back(correctedP4.Phi());
	  }
	  else{
	    prunedJetm_v2 .push_back(0.);
	    prunedJetpt_v2 .push_back(0.);
	    prunedJeteta_v2 .push_back(0.);
	    prunedJetphi_v2 .push_back(0.);
	  }
	}
      }
      else{

	if(jetsBoosted[i]->hasUserFloat("ak8PFJetsCHSPrunedMass"))
	  prunedJetmraw .push_back(jetsBoosted[i]->userFloat("ak8PFJetsCHSPrunedMass"));  
	// correct the pruned mass on the fly
	prunedJetm.push_back(jetsBoosted[i]->userFloat("ak8PFJetsCHSPrunedMass")*jetsBoosted[i]->correctedP4(jetsBoosted[i]->availableJECLevels().back()).Pt()/jetsBoosted[i]->correctedP4("L1FastJet").Pt());
      }
            
      if(not useMiniAODSubstructure){  // other pruning info
	if(isMC){
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:hadronFlavour"))
	    prunedJetHFlav  .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:hadronFlavour"));
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:partonFlavour"))
	    prunedJetPFlav .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:partonFlavour"));
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:genMass"))
	    prunedJetGenm  .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:genMass"));
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:genPt"))
	    prunedJetGenpt  .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:genPt"));
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:genEta"))
	    prunedJetGeneta  .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:genEta"));
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"PrunedMatched:genPhi"))
	    prunedJetGenphi  .push_back( jetsBoosted[i]->userFloat(jetsLabel+"PrunedMatched:genPhi"));
	}
      }


      // Subjets after pruning
      if(not useMiniAODSubstructure){

	if(jetsBoosted[i]->hasSubjets("Pruned")){    
	  pat::JetPtrCollection subjets = jetsBoosted[i]->subjets("Pruned");
	  if(subjets.size() > 0 ){
	    prunedSubJetpt_1  .push_back( subjets[0]->pt()); 
	    prunedSubJetm_1   .push_back( subjets[0]->mass()); 
	    prunedSubJetphi_1 .push_back( subjets[0]->phi()); 
	    prunedSubJeteta_1 .push_back( subjets[0]->eta());
	    prunedSubJetBtag_1 .push_back( subjets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	          
	    prunedSubJetptraw_1  .push_back( subjets[0]->correctedP4(0).pt()); 
	    prunedSubJetmraw_1   .push_back( subjets[0]->correctedP4(0).mass()); 
	    
	    if(subjets[0]->hasUserFloat(jetsLabel+"PrunedSubJetsQGL:qgLikelihood"))
	      prunedSubJetQGL_1 .push_back( subjets[0]->userFloat(jetsLabel+"PrunedSubJetsQGL:qgLikelihood"));
	          
	    if(isMC){
	      prunedSubJetHFlav_1 .push_back( subjets[0]->hadronFlavour()); 
	      prunedSubJetPFlav_1 .push_back( subjets[0]->partonFlavour()); 
	      if(subjets[0]->genJet()){
		prunedSubJetGenpt_1 .push_back( subjets[0]->genJet()->pt());
		prunedSubJetGenm_1  .push_back( subjets[0]->genJet()->mass()); 
		prunedSubJetGeneta_1  .push_back( subjets[0]->genJet()->eta()); 
		prunedSubJetGenphi_1  .push_back( subjets[0]->genJet()->phi()); 
	      }
	      else{
		prunedSubJetGenpt_1 .push_back(0);
		prunedSubJetGenm_1  .push_back(0);
		prunedSubJetGeneta_1  .push_back(0);
		prunedSubJetGenphi_1  .push_back(0);
	      }
	      if(addBTagScaleFactor)
		calculateBtagSF(subjets[0],"SubCSV",prunedSubJetBtagSF_1,prunedSubJetBtagSFUp_1,prunedSubJetBtagSFDown_1);
	    }
	  }
	  
	  if(subjets.size() > 1 ){

	    prunedSubJetpt_2  .push_back( subjets.at(1)->pt()); 
	    prunedSubJetm_2   .push_back( subjets.at(1)->mass()); 
	    prunedSubJetphi_2 .push_back( subjets.at(1)->phi()); 
	    prunedSubJeteta_2 .push_back( subjets.at(1)->eta());
	    prunedSubJetBtag_2 .push_back( subjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	          
	    prunedSubJetptraw_2  .push_back( subjets.at(1)->correctedP4(0).pt()); 
	    prunedSubJetmraw_2   .push_back( subjets.at(1)->correctedP4(0).mass()); 
	    
	    if(subjets.at(1)->hasUserFloat(jetsLabel+"PrunedSubJetsQGL:qgLikelihood"))
	      prunedSubJetQGL_2 .push_back( subjets.at(1)->userFloat(jetsLabel+"PrunedSubJetsQGL:qgLikelihood"));
	    
	          
	    if(isMC){
	      prunedSubJetHFlav_2 .push_back( subjets.at(1)->hadronFlavour()); 
	      prunedSubJetPFlav_2 .push_back( subjets.at(1)->partonFlavour()); 
	      if(subjets.at(1)->genJet()){
		prunedSubJetGenpt_2 .push_back( subjets.at(1)->genJet()->pt());
		prunedSubJetGenm_2  .push_back( subjets.at(1)->genJet()->mass());
		prunedSubJetGeneta_2 .push_back( subjets.at(1)->genJet()->eta());
		prunedSubJetGenphi_2 .push_back( subjets.at(1)->genJet()->phi());
	      }
	      else{
		prunedSubJetGenpt_2 .push_back(0);
		prunedSubJetGenm_2  .push_back(0);
		prunedSubJetGeneta_2  .push_back(0);
		prunedSubJetGenphi_2  .push_back(0);

	      }
	      if(addBTagScaleFactor)
		calculateBtagSF(subjets.at(1),"SubCSV",prunedSubJetBtagSF_2,prunedSubJetBtagSFUp_2,prunedSubJetBtagSFDown_2);
	    }
	  }
	}
      }            
    }
    else{ // Puppi case
      
      if(not useMiniAODSubstructure){

	boostedPuppiJetpt  .push_back( jetsBoosted[i]->pt());
	boostedPuppiJeteta .push_back( jetsBoosted[i]->eta());
	boostedPuppiJetphi .push_back( jetsBoosted[i]->phi());
	boostedPuppiJetm   .push_back( jetsBoosted[i]->mass());
	boostedPuppiJetBtag .push_back( jetsBoosted[i]->bDiscriminator("pfCombinedSecondaryVertexV2BJetTags"));
	boostedPuppiJetDoubleBtag .push_back( jetsBoosted[i]->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));
		
	if(jetsBoosted[i]->hasUserFloat("Njettiness"+jetsLabel+":tau1"))
	  boostedPuppiJettau1 .push_back( jetsBoosted[i]->userFloat("Njettiness"+jetsLabel+":tau1"));
	
	if(jetsBoosted[i]->hasUserFloat("Njettiness"+jetsLabel+":tau2"))
	  boostedPuppiJettau2 .push_back( jetsBoosted[i]->userFloat("Njettiness"+jetsLabel+":tau2"));
	
	if(jetsBoosted[i]->hasUserFloat("Njettiness"+jetsLabel+":tau3"))
	    boostedPuppiJettau3 .push_back( jetsBoosted[i]->userFloat("Njettiness"+jetsLabel+":tau3"));
	
	if(jetsBoosted[i]->hasUserFloat("Njettiness"+jetsLabel+":tau4"))
	  boostedPuppiJettau4 .push_back( jetsBoosted[i]->userFloat("Njettiness"+jetsLabel+":tau4"));
	
	if(isMC){
	  if(jetsBoosted[i]->genJet()){ // gen AK8 jet
	    boostedPuppiJetGenpt .push_back( jetsBoosted[i]->genJet()->pt());
	    boostedPuppiJetGenm  .push_back( jetsBoosted[i]->genJet()->mass());
	    boostedPuppiJetGeneta  .push_back( jetsBoosted[i]->genJet()->eta());
	    boostedPuppiJetGenphi  .push_back( jetsBoosted[i]->genJet()->phi());
	  }
	  else{
	    boostedPuppiJetGenpt .push_back(0);
	    boostedPuppiJetGeneta .push_back(0);
	    boostedPuppiJetGenphi .push_back(0);
	    boostedPuppiJetGenm .push_back(0);
	  }
	  boostedPuppiJetHFlav .push_back( jetsBoosted[i]->hadronFlavour());
	  boostedPuppiJetPFlav .push_back( jetsBoosted[i]->partonFlavour());  
	}
	
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"QGL:qgLikelihood"))
	  boostedPuppiJetQGL .push_back( jetsBoosted[i]->userFloat(jetsLabel+"QGL:qgLikelihood"));
	
	if(isMC){
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"GenNjettinessMatched:tau1"))
	    boostedPuppiJetGentau1 .push_back( jetsBoosted[i]->userFloat(jetsLabel+"GenNjettinessMatched:tau1"));
	  
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"GenNjettinessMatched:tau2"))
	      boostedPuppiJetGentau2 .push_back( jetsBoosted[i]->userFloat(jetsLabel+"GenNjettinessMatched:tau2"));
	  
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"GenNjettinessMatched:tau3"))
	    boostedPuppiJetGentau3 .push_back( jetsBoosted[i]->userFloat(jetsLabel+"GenNjettinessMatched:tau3"));
	  
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"GenNjettinessMatched:tau4"))
	    boostedPuppiJetGentau4 .push_back( jetsBoosted[i]->userFloat(jetsLabel+"GenNjettinessMatched:tau4"));
	}
	
	// soft drop matched jet
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:mass"))
	  softDropPuppiJetm .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:mass"));
	
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:pt"))
	  softDropPuppiJetpt .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:pt"));
	
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:eta"))
	  softDropPuppiJeteta .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:eta"));
	
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:phi"))
	  softDropPuppiJetphi .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:phi"));
	
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropQGLMatched:qgLikelihood"))
	  softDropPuppiJetQGL .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropQGLMatched:qgLikelihood"));
	
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"))
	  softDropPuppiJetBtag .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"))
	  softDropPuppiJetDoubleBtag .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:pfBoostedDoubleSecondaryVertexAK8BJetTags"));
	
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:rawpt"))
	  softDropPuppiJetptraw .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:rawpt"));
	
	if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:rawmass")){
	  softDropPuppiJetmraw .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:rawmass"));
	  
	  if(jetsBoosted[i]->availableJECSets().size()>1 and
	     jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:raweta") and
	     jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:rawphi") and
	     jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:rawpt")){
	    
	    TLorentzVector correctedP4;
	    correctedP4.SetPtEtaPhiM(jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:rawpt"),
				     jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:raweta"),
				     jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:rawphi"),
				     jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:rawmass")
				     );
	    correctedP4 *= 1./jetsBoosted[i]->jecFactor("Uncorrected","none",jetsBoosted[i]->availableJECSets().at(1));
	    softDropPuppiJetm_v2 .push_back(correctedP4.M());
	    softDropPuppiJetpt_v2 .push_back(correctedP4.Pt());
	    softDropPuppiJeteta_v2 .push_back(correctedP4.Eta());
	    softDropPuppiJetphi_v2 .push_back(correctedP4.Phi());
	  }
	  else{
	    softDropPuppiJetm_v2 .push_back(0.);
	    softDropPuppiJetpt_v2 .push_back(0.);
	    softDropPuppiJeteta_v2 .push_back(0.);
	    softDropPuppiJetphi_v2 .push_back(0.);
	  }  
	}
	
	if(isMC){
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:hadronFlavour"))
	    softDropPuppiJetHFlav  .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:hadronFlavour"));
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:partonFlavour"))
	    softDropPuppiJetPFlav .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:partonFlavour"));
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:genMass"))
	    softDropPuppiJetGenm  .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:genMass"));
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:genPt"))
	    softDropPuppiJetGenpt  .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:genPt"));
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:genEta"))
	    softDropPuppiJetGeneta  .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:genEta"));
	  if(jetsBoosted[i]->hasUserFloat(jetsLabel+"SoftDropMatched:genPhi"))
	    softDropPuppiJetGenphi  .push_back( jetsBoosted[i]->userFloat(jetsLabel+"SoftDropMatched:genPhi"));    
	}
	
	// sub-jets soft drop 
	if(jetsBoosted[i]->hasSubjets("SoftDrop")){
	  pat::JetPtrCollection subjets = jetsBoosted[i]->subjets("SoftDrop");
	  if(subjets.size() > 0){
	    softDropPuppiSubJetpt_1  .push_back( subjets[0]->pt()); 
	    softDropPuppiSubJetm_1   .push_back( subjets[0]->mass()); 
	    softDropPuppiSubJetphi_1 .push_back( subjets[0]->phi()); 
	    softDropPuppiSubJeteta_1 .push_back( subjets[0]->eta());
	    softDropPuppiSubJetBtag_1 .push_back( subjets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	    
	    softDropPuppiSubJetptraw_1  .push_back( subjets[0]->correctedP4(0).pt()); 
	    softDropPuppiSubJetmraw_1   .push_back( subjets[0]->correctedP4(0).mass()); 
	            
	    if(subjets[0]->hasUserFloat(jetsLabel+"SoftDropSubJetsQGL:qgLikelihood"))
	      softDropPuppiSubJetQGL_1 .push_back( subjets[0]->userFloat(jetsLabel+"SoftDropSubJetsQGL:qgLikelihood"));
	    
	    if(isMC){
	      softDropPuppiSubJetHFlav_1 .push_back( subjets[0]->hadronFlavour()); 
	      if(subjets[0]->genJet()){
		softDropPuppiSubJetGenpt_1 .push_back( subjets[0]->genJet()->pt());
		softDropPuppiSubJetGenm_1  .push_back( subjets[0]->genJet()->mass()); 
		softDropPuppiSubJetGeneta_1  .push_back( subjets[0]->genJet()->eta()); 
		softDropPuppiSubJetGenphi_1  .push_back( subjets[0]->genJet()->phi()); 
	      }
	      else{
		softDropPuppiSubJetGenpt_1 .push_back(0);
		softDropPuppiSubJetGeneta_1 .push_back(0);
		softDropPuppiSubJetGenphi_1 .push_back(0);
		softDropPuppiSubJetGenm_1 .push_back(0);
	      }
	      if(addBTagScaleFactor)
		  calculateBtagSF(subjets.at(0),"SubCSV",softDropPuppiSubJetBtagSF_1,softDropPuppiSubJetBtagSFUp_1,softDropPuppiSubJetBtagSFDown_1);
	    }
	  }
	  
	  if(subjets.size() > 1){
	    softDropPuppiSubJetpt_2  .push_back( subjets.at(1)->pt()); 
	    softDropPuppiSubJetm_2   .push_back( subjets.at(1)->mass()); 
	    softDropPuppiSubJetphi_2 .push_back( subjets.at(1)->phi()); 
	    softDropPuppiSubJeteta_2 .push_back( subjets.at(1)->eta());
	    softDropPuppiSubJetBtag_2 .push_back( subjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	    
	    if(subjets.at(1)->hasUserFloat(jetsLabel+"SoftDropSubJetsQGL:qgLikelihood"))
	      softDropPuppiSubJetQGL_2 .push_back( subjets.at(1)->userFloat(jetsLabel+"SoftDropSubJetsQGL:qgLikelihood"));
	    
	    softDropPuppiSubJetptraw_2  .push_back( subjets.at(1)->correctedP4(0).pt()); 
	    softDropPuppiSubJetmraw_2   .push_back( subjets.at(1)->correctedP4(0).mass()); 
	    
	    if(isMC){
	      softDropPuppiSubJetHFlav_2 .push_back( subjets.at(1)->hadronFlavour()); 
	      if(subjets.at(1)->genJet()){
		softDropPuppiSubJetGenpt_2 .push_back( subjets.at(1)->genJet()->pt());
		softDropPuppiSubJetGenm_2  .push_back( subjets.at(1)->genJet()->mass());
		softDropPuppiSubJetGeneta_2  .push_back( subjets.at(1)->genJet()->eta());
		softDropPuppiSubJetGenphi_2  .push_back( subjets.at(1)->genJet()->phi());
	      }
	      else{
		softDropPuppiSubJetGenpt_2 .push_back(0);
		softDropPuppiSubJetGenm_2 .push_back(0);
		softDropPuppiSubJetGeneta_2 .push_back(0);
		softDropPuppiSubJetGenphi_2 .push_back(0);
	      }
	      if(addBTagScaleFactor)
		calculateBtagSF(subjets.at(1),"SubCSV",softDropPuppiSubJetBtagSF_2,softDropPuppiSubJetBtagSFUp_2,softDropPuppiSubJetBtagSFDown_2);
	    }
	  }
	}
      }
      else{  // Use Mini-AOD info
	
	boostedPuppiJetpt  .push_back( jetsBoosted[i]->userFloat("ak8PFJetsPuppiValueMap:pt"));
	boostedPuppiJeteta  .push_back( jetsBoosted[i]->userFloat("ak8PFJetsPuppiValueMap:eta"));
	boostedPuppiJetphi  .push_back( jetsBoosted[i]->userFloat("ak8PFJetsPuppiValueMap:phi"));
	boostedPuppiJetm  .push_back( jetsBoosted[i]->userFloat("ak8PFJetsPuppiValueMap:mass"));
	
	boostedPuppiJettau1 .push_back( jetsBoosted[i]->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1"));
	boostedPuppiJettau2 .push_back( jetsBoosted[i]->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2"));
	boostedPuppiJettau3 .push_back( jetsBoosted[i]->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3"));

	// sub-jets soft drop 
	if(jetsBoosted[i]->hasSubjets("SoftDropPuppi")){
	  pat::JetPtrCollection subjets = jetsBoosted[i]->subjets("SoftDropPuppi");
	  
	  TLorentzVector puppi_softdrop;
	  TLorentzVector puppi_softdrop_subjet;
	  for(auto jet : subjets){
	    puppi_softdrop_subjet.SetPtEtaPhiM(jet->correctedP4(0).pt(),jet->correctedP4(0).eta(),jet->correctedP4(0).phi(),jet->correctedP4(0).mass());
	    puppi_softdrop += puppi_softdrop_subjet;
	  }
	  
	  softDropPuppiJetmraw.push_back(puppi_softdrop.M()); 
	  softDropPuppiJetptraw.push_back(puppi_softdrop.Pt()); 
	  
	  // not fully correct cause uses AK8CHS JECs
	  softDropPuppiJetm.push_back(puppi_softdrop.M()*jetsBoosted[i]->correctedP4(jetsBoosted[i]->availableJECLevels().back()).Pt()/jetsBoosted[i]->correctedP4("L1FastJet").Pt()); 
	  softDropPuppiJetpt.push_back(puppi_softdrop.Pt()*jetsBoosted[i]->correctedP4(jetsBoosted[i]->availableJECLevels().back()).Pt()/jetsBoosted[i]->correctedP4("L1FastJet").Pt()); 
	  softDropPuppiJeteta.push_back(puppi_softdrop.Eta());
	  softDropPuppiJetphi.push_back(puppi_softdrop.Phi());    
	  
	  if(subjets.size() > 0){
	    
	    softDropPuppiSubJetpt_1  .push_back( subjets[0]->pt()); 
	    softDropPuppiSubJetm_1   .push_back( subjets[0]->mass()); 
	    softDropPuppiSubJetphi_1 .push_back( subjets[0]->phi()); 
	    softDropPuppiSubJeteta_1 .push_back( subjets[0]->eta());
	    softDropPuppiSubJetBtag_1 .push_back( subjets[0]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	    
	    softDropPuppiSubJetptraw_1  .push_back( subjets[0]->correctedP4(0).pt()); 
	    softDropPuppiSubJetmraw_1   .push_back( subjets[0]->correctedP4(0).mass()); 
	    
	    if(isMC)
	      softDropPuppiSubJetHFlav_1 .push_back( subjets[0]->hadronFlavour()); 
	    
	    if(addBTagScaleFactor and isMC)
	      calculateBtagSF(subjets.at(0),"SubCSV",softDropPuppiSubJetBtagSF_1,softDropPuppiSubJetBtagSFUp_1,softDropPuppiSubJetBtagSFDown_1);
	  }
	        
	  if(subjets.size() > 1){
	    softDropPuppiSubJetpt_2  .push_back( subjets.at(1)->pt()); 
	    softDropPuppiSubJetm_2   .push_back( subjets.at(1)->mass()); 
	    softDropPuppiSubJetphi_2 .push_back( subjets.at(1)->phi()); 
	    softDropPuppiSubJeteta_2 .push_back( subjets.at(1)->eta());
	    softDropPuppiSubJetBtag_2 .push_back( subjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	    
	    softDropPuppiSubJetptraw_2  .push_back( subjets.at(1)->correctedP4(0).pt()); 
	    softDropPuppiSubJetmraw_2   .push_back( subjets.at(1)->correctedP4(0).mass()); 
	    
	    if(isMC)
	      softDropPuppiSubJetHFlav_2 .push_back( subjets.at(1)->hadronFlavour()); 

	    if(addBTagScaleFactor and isMC)
	      calculateBtagSF(subjets.at(1),"SubCSV",softDropPuppiSubJetBtagSF_2,softDropPuppiSubJetBtagSFUp_2,softDropPuppiSubJetBtagSFDown_2);
	    
	  }
	}
      }
    }
  }

  return true;
}


// to fill b-tag SF
void VJetTreeFiller::calculateBtagSF(const pat::Jet & jet, 
				    const std::string & algorithm, 
				    std::vector<float> & scalefactor, 
				    std::vector<float> & scalefactorUp, 
				    std::vector<float> & scalefactorDown){

  if(algorithm != "SubCSV") return;

  float jetPt = jet.pt();
  float jetEta = jet.eta();

  if(jet.hadronFlavour() == 5){            
    scalefactor.push_back(bMediumSubCSV.back().eval_auto_bounds("central",BTagEntry::FLAV_B,jetEta,jetPt));      
    scalefactorUp.push_back(bMediumSubCSV.back().eval_auto_bounds("up",BTagEntry::FLAV_B,jetEta,jetPt));      
    scalefactorDown.push_back(bMediumSubCSV.back().eval_auto_bounds("down",BTagEntry::FLAV_B,jetEta,jetPt));      
  }
  else if(jet.hadronFlavour() == 4){
    scalefactor.push_back(bMediumSubCSV.back().eval_auto_bounds("central",BTagEntry::FLAV_C,jetEta,jetPt));      
    scalefactorUp.push_back(bMediumSubCSV.back().eval_auto_bounds("up",BTagEntry::FLAV_C,jetEta,jetPt));      
    scalefactorDown.push_back(bMediumSubCSV.back().eval_auto_bounds("down",BTagEntry::FLAV_C,jetEta,jetPt));      
  }
  else{      
    scalefactor.push_back(bMediumSubCSV.back().eval_auto_bounds("central",BTagEntry::FLAV_UDSG,jetEta,jetPt));      
    scalefactorUp.push_back(bMediumSubCSV.back().eval_auto_bounds("up",BTagEntry::FLAV_UDSG,jetEta,jetPt));      
    scalefactorDown.push_back(bMediumSubCSV.back().eval_auto_bounds("down",BTagEntry::FLAV_UDSG,jetEta,jetPt));      
  }
}


//// fill jet collection
void VJetTreeFiller::fillJetCollections(const edm::Handle<std::vector<pat::Jet> > & jetsH, 
					const pat::MuonRefVector & muons, 
					const pat::ElectronRefVector & electrons,
					const pat::PhotonRefVector & photons, 
					std::vector<pat::JetRef> & jetsBoosted, 
					const bool & ispuppi){
  
  if(jetsH.isValid()){      
    for (auto jets_iter = jetsH->begin(); jets_iter != jetsH->end(); ++jets_iter) {
      //clean from leptons
      bool skipjet = false;
      for (std::size_t j = 0; j < muons.size(); j++) {
	if (cleanMuonJet && deltaR(muons[j]->eta(), muons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK8) 
	  skipjet = true;
      }
      for (std::size_t j = 0; j < electrons.size(); j++) {
	if (cleanElectronJet && deltaR(electrons[j]->eta(), electrons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK8) 
	  skipjet = true;
      }
      for (std::size_t j = 0; j < photons.size(); j++) {
	if (cleanPhotonJet && deltaR(photons[j]->eta(), photons[j]->phi(), jets_iter->eta(), jets_iter->phi()) < dRCleaningAK8) 
	  skipjet = true;
      }
      
      // jet in overlap with lepton
      if (skipjet) continue;
      
      // apply jet id
      bool passjetid = applyJetID(*jets_iter,jetidwp);            
      if (!passjetid) 
	continue;

      pat::JetRef jetref(jetsH, jets_iter - jetsH->begin());
      if(jetref.isAvailable() and jetref.isNonnull())
	jetsBoosted.push_back(jetref);
    }
    
    if(jetsBoosted.size() > 0) sort(jetsBoosted.begin(), jetsBoosted.end(), jetPtSorter);
  }  
}


////////
void VJetTreeFiller::DeclareAndSetBranches(){
  
  // AK8 Puppi jets                                                                                                                                                             
  if(not isTriggerTree and not isPhotonPurity and not applyDiMuonFilter and not applyDiElectronFilter and not applyPhotonJetsFilter and not isPuppi_){
    
    tree_->Branch("boostedJetpt",  "std::vector<float>", &boostedJetpt);
    tree_->Branch("boostedJeteta", "std::vector<float>", &boostedJeteta);
    tree_->Branch("boostedJetphi", "std::vector<float>", &boostedJetphi);
    tree_->Branch("boostedJetm",   "std::vector<float>", &boostedJetm);
    
    tree_->Branch("boostedJetGenpt",  "std::vector<float>", &boostedJetGenpt);
    tree_->Branch("boostedJetGeneta", "std::vector<float>", &boostedJetGeneta);
    tree_->Branch("boostedJetGenphi", "std::vector<float>", &boostedJetGenphi);
    tree_->Branch("boostedJetGenm",   "std::vector<float>", &boostedJetGenm);

    tree_->Branch("boostedJetHFlav", "std::vector<float>", &boostedJetHFlav);
    tree_->Branch("boostedJetPFlav", "std::vector<float>", &boostedJetPFlav);
    tree_->Branch("boostedJetQGL",   "std::vector<float>", &boostedJetQGL);

    tree_->Branch("boostedJettau1", "std::vector<float>", &boostedJettau1);
    tree_->Branch("boostedJettau2", "std::vector<float>", &boostedJettau2);
    tree_->Branch("boostedJettau3", "std::vector<float>", &boostedJettau3);


    if(not useMiniAODSubstructure){
      tree_->Branch("boostedJetBtag",  "std::vector<float>", &boostedJetBtag);
      tree_->Branch("boostedJetDoubleBtag", "std::vector<float>", &boostedJetDoubleBtag);
      tree_->Branch("boostedJettau4", "std::vector<float>", &boostedJettau4);
      tree_->Branch("boostedJetGentau1", "std::vector<float>", &boostedJetGentau1);
      tree_->Branch("boostedJetGentau2", "std::vector<float>", &boostedJetGentau2);
      tree_->Branch("boostedJetGentau3", "std::vector<float>", &boostedJetGentau3);
      tree_->Branch("boostedJetGentau4", "std::vector<float>", &boostedJetGentau4);
      
    }

    tree_->Branch("prunedJetm",     "std::vector<float>", &prunedJetm);
    tree_->Branch("prunedJetmraw",  "std::vector<float>", &prunedJetmraw);
    
    if(not useMiniAODSubstructure){

      tree_->Branch("prunedJetpt",    "std::vector<float>", &prunedJetpt);
      tree_->Branch("prunedJeteta",   "std::vector<float>", &prunedJeteta);
      tree_->Branch("prunedJetphi",   "std::vector<float>", &prunedJetphi);

      tree_->Branch("prunedJetptraw", "std::vector<float>", &prunedJetptraw);

      tree_->Branch("prunedJetm_v2",   "std::vector<float>", &prunedJetm_v2);
      tree_->Branch("prunedJetpt_v2",  "std::vector<float>", &prunedJetpt_v2);
      tree_->Branch("prunedJeteta_v2", "std::vector<float>", &prunedJeteta_v2);
      tree_->Branch("prunedJetphi_v2", "std::vector<float>", &prunedJetphi_v2);
      
      tree_->Branch("prunedJetGenpt",  "std::vector<float>", &prunedJetGenpt);
      tree_->Branch("prunedJetGeneta", "std::vector<float>", &prunedJetGeneta);
      tree_->Branch("prunedJetGenphi", "std::vector<float>", &prunedJetGenphi);
      tree_->Branch("prunedJetGenm",   "std::vector<float>", &prunedJetGenm);
      
      tree_->Branch("prunedJetHFlav", "std::vector<float>", &prunedJetHFlav);
      tree_->Branch("prunedJetPFlav", "std::vector<float>", &prunedJetPFlav);
      tree_->Branch("prunedJetQGL",   "std::vector<float>", &prunedJetQGL);


      tree_->Branch("prunedSubJetpt_1",  "std::vector<float>",  &prunedSubJetpt_1);
      tree_->Branch("prunedSubJeteta_1", "std::vector<float>",  &prunedSubJeteta_1);
      tree_->Branch("prunedSubJetphi_1", "std::vector<float>",  &prunedSubJetphi_1);
      tree_->Branch("prunedSubJetm_1",   "std::vector<float>", &prunedSubJetm_1);
      tree_->Branch("prunedSubJetGenpt_1","std::vector<float>",  &prunedSubJetGenpt_1);
      tree_->Branch("prunedSubJetGenm_1", "std::vector<float>", &prunedSubJetGenm_1);
      tree_->Branch("prunedSubJetGeneta_1", "std::vector<float>", &prunedSubJetGeneta_1);
      tree_->Branch("prunedSubJetGenphi_1", "std::vector<float>", &prunedSubJetGenphi_1);
      tree_->Branch("prunedSubJetHFlav_1",  "std::vector<float>", &prunedSubJetHFlav_1);
      tree_->Branch("prunedSubJetPFlav_1",  "std::vector<float>", &prunedSubJetPFlav_1);
      tree_->Branch("prunedSubJetQGL_1",    "std::vector<float>", &prunedSubJetQGL_1);
      tree_->Branch("prunedSubJetBtag_1",   "std::vector<float>", &prunedSubJetBtag_1);
      tree_->Branch("prunedSubJetptraw_1",  "std::vector<float>", &prunedSubJetptraw_1);
      tree_->Branch("prunedSubJetmraw_1",   "std::vector<float>", &prunedSubJetmraw_1);
      tree_->Branch("prunedSubJetBtagSF_1",   "std::vector<float>", &prunedSubJetBtagSF_1);
      tree_->Branch("prunedSubJetBtagSFUp_1",   "std::vector<float>", &prunedSubJetBtagSFUp_1);
      tree_->Branch("prunedSubJetBtagSFDown_1",   "std::vector<float>", &prunedSubJetBtagSFDown_1);
      
      tree_->Branch("prunedSubJetpt_2",  "std::vector<float>",  &prunedSubJetpt_2);
      tree_->Branch("prunedSubJeteta_2", "std::vector<float>",  &prunedSubJeteta_2);
      tree_->Branch("prunedSubJetphi_2", "std::vector<float>",  &prunedSubJetphi_2);
      tree_->Branch("prunedSubJetm_2",   "std::vector<float>", &prunedSubJetm_2);
      tree_->Branch("prunedSubJetGenpt_2","std::vector<float>",  &prunedSubJetGenpt_2);
      tree_->Branch("prunedSubJetGenm_2", "std::vector<float>", &prunedSubJetGenm_2);
      tree_->Branch("prunedSubJetGeneta_2", "std::vector<float>", &prunedSubJetGeneta_2);
      tree_->Branch("prunedSubJetGenphi_2", "std::vector<float>", &prunedSubJetGenphi_2);
      tree_->Branch("prunedSubJetHFlav_2",  "std::vector<float>", &prunedSubJetHFlav_2);
      tree_->Branch("prunedSubJetPFlav_2", "std::vector<float>", &prunedSubJetPFlav_2);
      tree_->Branch("prunedSubJetQGL_2",   "std::vector<float>", &prunedSubJetQGL_2);
      tree_->Branch("prunedSubJetBtag_2",  "std::vector<float>", &prunedSubJetBtag_2);
      tree_->Branch("prunedSubJetptraw_2", "std::vector<float>", &prunedSubJetptraw_2);
      tree_->Branch("prunedSubJetmraw_2",  "std::vector<float>", &prunedSubJetmraw_2);
      tree_->Branch("prunedSubJetBtagSF_2",   "std::vector<float>", &prunedSubJetBtagSF_2);
      tree_->Branch("prunedSubJetBtagSFUp_2",   "std::vector<float>", &prunedSubJetBtagSFUp_2);
      tree_->Branch("prunedSubJetBtagSFDown_2",   "std::vector<float>", &prunedSubJetBtagSFDown_2);
    }
  }
  else if(not isTriggerTree and not isPhotonPurity and not applyDiMuonFilter and not applyDiElectronFilter and not applyPhotonJetsFilter and isPuppi_){

    tree_->Branch("boostedPuppiJetpt", "std::vector<float>", &boostedPuppiJetpt);
    tree_->Branch("boostedPuppiJeteta", "std::vector<float>", &boostedPuppiJeteta);
    tree_->Branch("boostedPuppiJetphi", "std::vector<float>", &boostedPuppiJetphi);
    tree_->Branch("boostedPuppiJetm", "std::vector<float>", &boostedPuppiJetm);
    tree_->Branch("boostedPuppiJettau1", "std::vector<float>", &boostedPuppiJettau1);
    tree_->Branch("boostedPuppiJettau2", "std::vector<float>", &boostedPuppiJettau2);
    tree_->Branch("boostedPuppiJettau3", "std::vector<float>", &boostedPuppiJettau3);

    if(not useMiniAODSubstructure){

      tree_->Branch("boostedPuppiJetGenpt", "std::vector<float>", &boostedPuppiJetGenpt);
      tree_->Branch("boostedPuppiJetGeneta", "std::vector<float>", &boostedPuppiJetGeneta);
      tree_->Branch("boostedPuppiJetGenphi", "std::vector<float>", &boostedPuppiJetGenphi);
      tree_->Branch("boostedPuppiJetGenm", "std::vector<float>", &boostedPuppiJetGenm);
      
      tree_->Branch("boostedPuppiJetHFlav", "std::vector<float>", &boostedPuppiJetHFlav);
      tree_->Branch("boostedPuppiJetPFlav", "std::vector<float>", &boostedPuppiJetPFlav);    
      tree_->Branch("boostedPuppiJetQGL", "std::vector<float>", &boostedPuppiJetQGL);
      tree_->Branch("boostedPuppiJetBtag", "std::vector<float>", &boostedPuppiJetBtag);
      tree_->Branch("boostedPuppiJetDoubleBtag", "std::vector<float>", &boostedPuppiJetDoubleBtag);

      tree_->Branch("boostedPuppiJettau4", "std::vector<float>", &boostedPuppiJettau4);

    }

    tree_->Branch("softDropPuppiJetpt", "std::vector<float>", &softDropPuppiJetpt);
    tree_->Branch("softDropPuppiJeteta", "std::vector<float>", &softDropPuppiJeteta);
    tree_->Branch("softDropPuppiJetphi", "std::vector<float>", &softDropPuppiJetphi);
    tree_->Branch("softDropPuppiJetm", "std::vector<float>", &softDropPuppiJetm);
    tree_->Branch("softDropPuppiJetptraw", "std::vector<float>", &softDropPuppiJetptraw);
    tree_->Branch("softDropPuppiJetmraw", "std::vector<float>", &softDropPuppiJetmraw);

    tree_->Branch("softDropPuppiJetm_v2", "std::vector<float>", &softDropPuppiJetm_v2);
    tree_->Branch("softDropPuppiJetpt_v2", "std::vector<float>", &softDropPuppiJetpt_v2);
    tree_->Branch("softDropPuppiJeteta_v2", "std::vector<float>", &softDropPuppiJeteta_v2);
    tree_->Branch("softDropPuppiJetphi_v2", "std::vector<float>", &softDropPuppiJetphi_v2);
      
    if(not useMiniAODSubstructure){

      tree_->Branch("softDropPuppiJetGenpt", "std::vector<float>", &softDropPuppiJetGenpt);
      tree_->Branch("softDropPuppiJetGeneta", "std::vector<float>", &softDropPuppiJetGeneta);
      tree_->Branch("softDropPuppiJetGenphi", "std::vector<float>", &softDropPuppiJetGenphi);
      tree_->Branch("softDropPuppiJetGenm", "std::vector<float>", &softDropPuppiJetGenm);
      
      tree_->Branch("softDropPuppiJetHFlav", "std::vector<float>", &softDropPuppiJetHFlav);
      tree_->Branch("softDropPuppiJetPFlav", "std::vector<float>", &softDropPuppiJetPFlav);
      tree_->Branch("softDropPuppiJetQGL", "std::vector<float>", &softDropPuppiJetQGL);
      tree_->Branch("softDropPuppiJetBtag", "std::vector<float>", &softDropPuppiJetBtag);
      tree_->Branch("softDropPuppiJetDoubleBtag", "std::vector<float>", &softDropPuppiJetDoubleBtag);
      
    }
    
    tree_->Branch("softDropPuppiSubJetpt_1","std::vector<float>",  &softDropPuppiSubJetpt_1);
    tree_->Branch("softDropPuppiSubJeteta_1","std::vector<float>",  &softDropPuppiSubJeteta_1);
    tree_->Branch("softDropPuppiSubJetphi_1","std::vector<float>",  &softDropPuppiSubJetphi_1);
    tree_->Branch("softDropPuppiSubJetm_1", "std::vector<float>", &softDropPuppiSubJetm_1);
    tree_->Branch("softDropPuppiSubJetBtagSF_1", "std::vector<float>", &softDropPuppiSubJetBtagSF_1);
    tree_->Branch("softDropPuppiSubJetBtagSFUp_1", "std::vector<float>", &softDropPuppiSubJetBtagSFUp_1);
    tree_->Branch("softDropPuppiSubJetBtagSFDown_1", "std::vector<float>", &softDropPuppiSubJetBtagSFDown_1);
    tree_->Branch("softDropPuppiSubJetHFlav_1", "std::vector<float>", &softDropPuppiSubJetHFlav_1);

    tree_->Branch("softDropPuppiSubJetpt_2","std::vector<float>",  &softDropPuppiSubJetpt_2);
    tree_->Branch("softDropPuppiSubJeteta_2","std::vector<float>",  &softDropPuppiSubJeteta_2);
    tree_->Branch("softDropPuppiSubJetphi_2","std::vector<float>",  &softDropPuppiSubJetphi_2);
    tree_->Branch("softDropPuppiSubJetm_2", "std::vector<float>", &softDropPuppiSubJetm_2);
    tree_->Branch("softDropPuppiSubJetBtagSF_2", "std::vector<float>", &softDropPuppiSubJetBtagSF_2);
    tree_->Branch("softDropPuppiSubJetBtagSFUp_2", "std::vector<float>", &softDropPuppiSubJetBtagSFUp_2);
    tree_->Branch("softDropPuppiSubJetBtagSFDown_2", "std::vector<float>", &softDropPuppiSubJetBtagSFDown_2);
    tree_->Branch("softDropPuppiSubJetHFlav_2", "std::vector<float>", &softDropPuppiSubJetHFlav_2);

    if(not useMiniAODSubstructure){

      tree_->Branch("softDropPuppiSubJetGenpt_1","std::vector<float>",  &softDropPuppiSubJetGenpt_1);
      tree_->Branch("softDropPuppiSubJetGenm_1", "std::vector<float>", &softDropPuppiSubJetGenm_1);
      tree_->Branch("softDropPuppiSubJetGeneta_1", "std::vector<float>", &softDropPuppiSubJetGeneta_1);
      tree_->Branch("softDropPuppiSubJetGenphi_1", "std::vector<float>", &softDropPuppiSubJetGenphi_1);
      tree_->Branch("softDropPuppiSubJetPFlav_1", "std::vector<float>", &softDropPuppiSubJetPFlav_1);
      tree_->Branch("softDropPuppiSubJetQGL_1", "std::vector<float>", &softDropPuppiSubJetQGL_1);
      tree_->Branch("softDropPuppiSubJetBtag_1", "std::vector<float>", &softDropPuppiSubJetBtag_1);
      tree_->Branch("softDropPuppiSubJetptraw_1", "std::vector<float>", &softDropPuppiSubJetptraw_1);
      tree_->Branch("softDropPuppiSubJetmraw_1", "std::vector<float>", &softDropPuppiSubJetmraw_1);

      tree_->Branch("softDropPuppiSubJetGenpt_2","std::vector<float>",  &softDropPuppiSubJetGenpt_2);
      tree_->Branch("softDropPuppiSubJetGenm_2", "std::vector<float>", &softDropPuppiSubJetGenm_2);
      tree_->Branch("softDropPuppiSubJetGeneta_2", "std::vector<float>", &softDropPuppiSubJetGeneta_2);
      tree_->Branch("softDropPuppiSubJetGenphi_2", "std::vector<float>", &softDropPuppiSubJetGenphi_2);
      tree_->Branch("softDropPuppiSubJetPFlav_2", "std::vector<float>", &softDropPuppiSubJetPFlav_2);
      tree_->Branch("softDropPuppiSubJetQGL_2", "std::vector<float>", &softDropPuppiSubJetQGL_2);
      tree_->Branch("softDropPuppiSubJetBtag_2", "std::vector<float>", &softDropPuppiSubJetBtag_2);
      tree_->Branch("softDropPuppiSubJetptraw_2", "std::vector<float>", &softDropPuppiSubJetptraw_2);
      tree_->Branch("softDropPuppiSubJetmraw_2", "std::vector<float>", &softDropPuppiSubJetmraw_2);
    }
  }
}
