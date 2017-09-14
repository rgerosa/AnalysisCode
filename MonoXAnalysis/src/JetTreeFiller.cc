#include "AnalysisCode/MonoXAnalysis/interface/JetTreeFiller.h"

JetTreeFiller::JetTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree):
  muonsTag       (iConfig.getParameter<edm::InputTag>("muons")),
  electronsTag   (iConfig.getParameter<edm::InputTag>("electrons")),
  photonsTag     (iConfig.getParameter<edm::InputTag>("photons")),
  jetsTag        (iConfig.getParameter<edm::InputTag>("jets")),
  jetsJESUpTag   (iConfig.existsAs<edm::InputTag>("jetsJESUp") ? iConfig.getParameter<edm::InputTag>("jetsJESUp") : edm::InputTag("")),
  jetsJESDwTag   (iConfig.existsAs<edm::InputTag>("jetsJESDw") ? iConfig.getParameter<edm::InputTag>("jetsJESDw") : edm::InputTag("")),
  jetsJERTag     (iConfig.existsAs<edm::InputTag>("jetsJER") ? iConfig.getParameter<edm::InputTag>("jetsJER") : edm::InputTag("")),
  isMC           (iConfig.existsAs<bool>("isMC") ? iConfig.getParameter<bool>("isMC") : false),
  isTriggerTree  (iConfig.existsAs<bool>("isTriggerTree") ? iConfig.getParameter<bool>("isTriggerTree") : false),
  isQCDTree      (iConfig.existsAs<bool>("isQCDTree") ? iConfig.getParameter<bool>("isQCDTree") : false),
  isPhotonPurity (iConfig.existsAs<bool>("isPhotonPurity") ? iConfig.getParameter<bool>("isPhotonPurity") : false),
  cleanMuonJet   (iConfig.existsAs<bool>("cleanMuonJet") ? iConfig.getParameter<bool>("cleanMuonJet") : false),
  cleanElectronJet (iConfig.existsAs<bool>("cleanElectronJet") ? iConfig.getParameter<bool>("cleanElectronJet") : false),
  cleanPhotonJet (iConfig.existsAs<bool>("cleanPhotonJet") ? iConfig.getParameter<bool>("cleanPhotonJet") : false),
  dRCleaningAK4  (iConfig.existsAs<double>("dRCleaningAK4") ? iConfig.getParameter<double>("dRCleaningAK4") : 0.4),
  jetidwp        (iConfig.existsAs<std::string>("jetidwp") ? iConfig.getParameter<std::string>("jetidwp") : "loose"),
  pileupjetidwp  (iConfig.existsAs<std::string>("pileupjetidwp") ? iConfig.getParameter<std::string>("pileupjetidwp") : "medium"),
  applypileupjetid(iConfig.existsAs<bool>("applypileupjetid") ? iConfig.getParameter<bool>("applypileupjetid") : false),
  btaggingCSVWP  (iConfig.getParameter<double>("btaggingCSVWP")),
  btaggingMVAWP  (iConfig.getParameter<double>("btaggingMVAWP")),
  minJetPtCountAK4 (iConfig.existsAs<double>("minJetPtCountAK4") ? iConfig.getParameter<double>("minJetPtCountAK4") : 30),
  minJetPtBveto    (iConfig.existsAs<double>("minJetPtBveto") ? iConfig.getParameter<double>("minJetPtBveto") : 20),
  minJetPtAK4Store (iConfig.existsAs<double>("minJetPtAK4Store") ? iConfig.getParameter<double>("minJetPtAK4Store") : 20),
  addBTagScaleFactor(iConfig.existsAs<bool>("addBTagScaleFactor") ? iConfig.getParameter<bool>("addBTagScaleFactor") : false){

  isPuppi_ = isPuppi;

  jetsToken = iC.consumes<std::vector<pat::Jet> > (jetsTag);
  if(jetsJESUpTag.label() != "")
    jetsJESUpToken = iC.consumes<std::vector<pat::Jet> > (jetsJESUpTag);
  if(jetsJESDwTag.label() != "")
    jetsJESDwToken = iC.consumes<std::vector<pat::Jet> > (jetsJESDwTag);
  if(jetsJERTag.label() != "")
    jetsJERToken   = iC.consumes<std::vector<pat::Jet> > (jetsJERTag);
  
  muonsToken     = iC.consumes<pat::MuonRefVector> (muonsTag);
  electronsToken = iC.consumes<pat::ElectronRefVector> (electronsTag);
  photonsToken   = iC.consumes<pat::PhotonRefVector> (photonsTag);

  //////
  if(addBTagScaleFactor and isMC){
  
    bTagScaleFactorFileCSV = iConfig.getParameter<edm::FileInPath>("bTagScaleFactorFileCSV");

    if ( bTagScaleFactorFileCSV.location()!=edm::FileInPath::Local)
      throw cms::Exception("JetTreeFiller") << " Failed to find File = " << bTagScaleFactorFileCSV << " !!\n";

    calibCSV = BTagCalibration("CSVv2",bTagScaleFactorFileCSV.fullPath());
    bMediumCSV.push_back(BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"})); // for light flavor
    bMediumCSV.back().load(calibCSV,BTagEntry::FLAV_B,"comb");
    bMediumCSV.back().load(calibCSV,BTagEntry::FLAV_C,"comb");
    bMediumCSV.back().load(calibCSV,BTagEntry::FLAV_UDSG,"incl");

    bTagScaleFactorFileMVA = iConfig.getParameter<edm::FileInPath>("bTagScaleFactorFileMVA");
    if ( bTagScaleFactorFileMVA.location()!=edm::FileInPath::Local)
      throw cms::Exception("JetTreeFiller") << " Failed to find File = " << bTagScaleFactorFileMVA << " !!\n";

    calibMVA = BTagCalibration("CMVAv2",bTagScaleFactorFileMVA.fullPath());
    bMediumMVA.push_back(BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central",{"up","down"})); // for light flavor
    bMediumMVA.back().load(calibMVA,BTagEntry::FLAV_B,"ttbar");
    bMediumMVA.back().load(calibMVA,BTagEntry::FLAV_C,"ttbar");
    bMediumMVA.back().load(calibMVA,BTagEntry::FLAV_UDSG,"incl");
    
  }

  tree_ = tree;
  DeclareAndSetBranches();
    
}

/////
void JetTreeFiller::initBranches(){

  njets       = 0; njetsinc    = 0; njetsincup    = 0; njetsincdw    = 0; njetsincjer    = 0;
  nbjets      = 0; nbjetslowpt   = 0; nbjetsMVA     = 0; nbjetsMVAlowpt = 0;
    
  combinejetpt        .clear(); combinejeteta       .clear(); combinejetphi       .clear(); combinejetbtag      .clear(); 
  combinejetCHfrac    .clear(); combinejetNHfrac    .clear(); combinejetEMfrac    .clear(); combinejetCEMfrac   .clear(); combinejetmetdphi  .clear();
  combinejetPHfrac    .clear(); combinejetELfrac    .clear(); combinejetMUfrac    .clear(); combinejetHFHfrac   .clear(); combinejetHFEMfrac .clear();
  combinejetCHmult    .clear(); combinejetNHmult    .clear(); combinejetPHmult    .clear(); combinejetMUmult    .clear(); combinejetHFHmult  .clear(); combinejetHFEMmult .clear();
 
  combinejetHFlav     .clear(); combinejetPFlav     .clear(); combinejetQGL       .clear(); combinejetPUID      .clear();
  combinejetGenpt     .clear(); combinejetGeneta    .clear(); combinejetGenphi    .clear(); combinejetGenm      .clear(); 
  combinejetm         .clear(); combinejetbtagMVA   .clear();
  combinejetBtagSF .clear(); combinejetBtagSFUp .clear(); combinejetBtagSFDown .clear();
  combinejetBtagMVASF .clear(); combinejetBtagMVASFUp .clear(); combinejetBtagMVASFDown .clear();
  combinejetPassPUID.clear();

  combinejetptup        .clear(); combinejetetaup        .clear(); combinejetphiup        .clear(); combinejetmup        .clear();
  combinejetptdw        .clear(); combinejetetadw        .clear(); combinejetphidw        .clear(); combinejetmdw        .clear();
  combinejetptjer       .clear(); combinejetetajer       .clear(); combinejetphijer       .clear(); combinejetmjer       .clear();

  jetjetdphi = 0.0; ht = 0.; htinc  = 0.; ht30 = 0.;

  npuppijets       = 0; npuppijetsinc    = 0; npuppijetsincup = 0;
  npuppijetsincdw  = 0; npuppijetsincjer = 0; npuppibjets     = 0; npuppibjetslowpt = 0;
  npuppibjetsMVA   = 0; npuppibjetsMVAlowpt = 0;
      
  combinePuppijetpt        .clear(); combinePuppijeteta       .clear(); combinePuppijetphi       .clear(); combinePuppijetbtag      .clear(); 
  combinePuppijetCHfrac    .clear(); combinePuppijetbtagMVA   .clear();
  combinePuppijetNHfrac    .clear(); combinePuppijetEMfrac    .clear(); combinePuppijetCEMfrac   .clear(); combinePuppijetmetdphi   .clear();
  combinePuppijetHFlav     .clear(); combinePuppijetPFlav     .clear(); combinePuppijetQGL       .clear(); 
  combinePuppijetGenpt     .clear(); combinePuppijetGeneta    .clear(); combinePuppijetGenphi    .clear(); combinePuppijetGenm      .clear();
  combinePuppijetm         .clear(); 
  combinePuppijetBtagSF    .clear(); combinePuppijetBtagSFUp .clear(); combinePuppijetBtagSFDown .clear();
  combinePuppijetBtagMVASF    .clear(); combinePuppijetBtagMVASFUp .clear(); combinePuppijetBtagMVASFDown .clear();

  combinePuppijetptup        .clear(); combinePuppijetetaup       .clear(); combinePuppijetphiup       .clear(); combinePuppijetmup     .clear(); 
  combinePuppijetptdw        .clear(); combinePuppijetetadw       .clear(); combinePuppijetphidw       .clear(); combinePuppijetmdw     .clear(); 
  combinePuppijetptjer       .clear(); combinePuppijetetajer      .clear(); combinePuppijetphijer      .clear(); combinePuppijetmjer    .clear(); 

  PuppijetPuppijetdphi = 0.0;    Puppiht = 0.;
  
}

/////
bool JetTreeFiller::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace pat;

  this->initBranches();

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
  Handle<vector<pat::Jet> > jetsJESUpH;
  Handle<vector<pat::Jet> > jetsJESDwH;
  Handle<vector<pat::Jet> > jetsJERH;
  iEvent.getByToken(jetsToken, jetsH);
  if(jetsJESUpTag.label() != "")
    iEvent.getByToken(jetsJESUpToken, jetsJESUpH);
  if(jetsJESDwTag.label() != "")
    iEvent.getByToken(jetsJESDwToken, jetsJESDwH);
  if(jetsJERTag.label() != "")
    iEvent.getByToken(jetsJERToken, jetsJERH);

  vector<pat::JetRef> alljets;
  vector<pat::JetRef> incjets;
  vector<pat::JetRef> jets;
  vector<pat::JetRef> incjets_jesup;
  vector<pat::JetRef> incjets_jesdw;
  vector<pat::JetRef> incjets_jer;

  if(jetsJESUpH.isValid())
    fillJetCollections(jetsJESUpH,muons,electrons,photons,incjets_jesup,alljets);
  alljets.clear();
  if(jetsJESDwH.isValid())
    fillJetCollections(jetsJESDwH,muons,electrons,photons,incjets_jesdw,alljets);
  alljets.clear();
  if(jetsJERH.isValid())
    fillJetCollections(jetsJERH,muons,electrons,photons,incjets_jer,alljets);
  alljets.clear();
  if(jetsH.isValid())
    fillJetCollections(jetsH,muons,electrons,photons,incjets,alljets);

  
  // only central jets for nominal scale
  for (size_t i = 0; i < incjets.size(); i++) {
    if (fabs(incjets[i]->eta()) <= 2.5) 
      jets.push_back(incjets[i]);
  }        

  // sort them in pt
  if(jets.size() > 0)  sort(jets.begin(), jets.end(), jetPtSorter);

  // all jets
  for(size_t i = 0; i < incjets.size(); i++){
    if(incjets[i]->pt() > minJetPtCountAK4 and not isPuppi_)
      njetsinc++;
    else if(incjets[i]->pt() > minJetPtCountAK4 and  isPuppi_)
      npuppijetsinc++;
  }
    
  for(size_t i = 0; i < incjets_jesup.size(); i++){
    if(incjets_jesup[i]->pt() > minJetPtCountAK4 and not isPuppi_) 
      njetsincup++;
    else if(incjets_jesup[i]->pt() > minJetPtCountAK4 and isPuppi_) 
      npuppijetsincup++;
  }

  for(size_t i = 0; i < incjets_jesdw.size(); i++){
    if(incjets_jesdw[i]->pt() > minJetPtCountAK4 and not isPuppi_) 
      njetsincdw++;
    else if(incjets_jesdw[i]->pt() > minJetPtCountAK4 and isPuppi_) 
      npuppijetsincdw++;
  }
    
  for(size_t i = 0; i < incjets_jer.size(); i++){
    if(incjets_jer[i]->pt() > minJetPtCountAK4 and not isPuppi_) 
      njetsincjer++;
    else if(incjets_jer[i]->pt() > minJetPtCountAK4 and isPuppi_) 
      npuppijetsincjer++;    
  }
    
  // only central jets
  for (size_t i = 0; i < jets.size(); i++) {
      
    if (jets[i]->pt() > minJetPtCountAK4 and not isPuppi_) njets++;
    else (jets[i]->pt() > minJetPtCountAK4 and isPuppi_) npuppijets++;

    // btagging
    if (jets[i]->pt() > minJetPtCountAK4 && fabs(jets[i]->eta()) < 2.4 && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btaggingCSVWP and not isPuppi_) nbjets++;
    else if (jets[i]->pt() > minJetPtCountAK4 && fabs(jets[i]->eta()) < 2.4 && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btaggingCSVWP and isPuppi_) npuppibjets++;

    if (jets[i]->pt() > minJetPtBveto    && fabs(jets[i]->eta()) < 2.4 && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btaggingCSVWP and not isPuppi_) nbjetslowpt++;
    else if (jets[i]->pt() > minJetPtBveto    && fabs(jets[i]->eta()) < 2.4 && jets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") > btaggingCSVWP and isPuppi) npuppibjetslowpt++;
      
    if (jets[i]->pt() > minJetPtCountAK4 && fabs(jets[i]->eta()) < 2.4 && jets[i]->bDiscriminator("pfCombinedMVAV2BJetTags") > btaggingMVAWP and not isPuppi_) nbjetsMVA++;
    else if (jets[i]->pt() > minJetPtCountAK4 && fabs(jets[i]->eta()) < 2.4 && jets[i]->bDiscriminator("pfCombinedMVAV2BJetTags") > btaggingMVAWP and isPuppi_)npuppibjetsMVA++;

    if (jets[i]->pt() > minJetPtBveto && fabs(jets[i]->eta()) < 2.4 && jets[i]->bDiscriminator("pfCombinedMVAV2BJetTags") > btaggingMVAWP and not isPuppi_) nbjetsMVAlowpt++;
    else if (jets[i]->pt() > minJetPtBveto && fabs(jets[i]->eta()) < 2.4 && jets[i]->bDiscriminator("pfCombinedMVAV2BJetTags") > btaggingMVAWP and isPuppi_) npuppibjetsMVAlowpt++;
  }
    

  // fill collections
  for(size_t i = 0; i < incjets.size(); i++){

    if (incjets[i]->pt() > minJetPtAK4Store){ //      
      if(not isPuppi_){ // standard jets

	combinejetpt. push_back(incjets[i]->pt());
	combinejeteta.push_back(incjets[i]->eta());
	combinejetphi.push_back(incjets[i]->phi());
	combinejetm.  push_back(incjets[i]->mass());
	combinejetbtag.push_back(incjets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	combinejetbtagMVA.push_back(incjets[i]->bDiscriminator("pfCombinedMVAV2BJetTags"));
	combinejetCHfrac .push_back(incjets[i]->chargedHadronEnergyFraction());
	combinejetNHfrac .push_back(incjets[i]->neutralHadronEnergyFraction());
	combinejetEMfrac .push_back(incjets[i]->neutralEmEnergyFraction());
	combinejetCEMfrac.push_back(incjets[i]->chargedEmEnergyFraction());
	combinejetPHfrac .push_back(incjets[i]->photonEnergyFraction());
	combinejetELfrac .push_back(incjets[i]->electronEnergyFraction());
	combinejetMUfrac .push_back(incjets[i]->muonEnergyFraction());
	combinejetHFHfrac .push_back(incjets[i]->HFHadronEnergyFraction());
	combinejetHFEMfrac.push_back(incjets[i]->HFEMEnergyFraction());
	combinejetCHmult  .push_back(incjets[i]->chargedHadronMultiplicity());
	combinejetNHmult  .push_back(incjets[i]->neutralHadronMultiplicity());
	combinejetPHmult  .push_back(incjets[i]->photonMultiplicity());
	combinejetELmult  .push_back(incjets[i]->electronMultiplicity());
	combinejetMUmult  .push_back(incjets[i]->muonMultiplicity());
	combinejetHFHmult .push_back(incjets[i]->HFHadronMultiplicity());
	combinejetHFEMmult.push_back(incjets[i]->HFEMMultiplicity());

	if(incjets[i]->hasUserFloat("QGTagger:qgLikelihood"))
	  combinejetQGL.push_back(incjets[i]->userFloat("QGTagger:qgLikelihood")); 
	// pileup jet id
	if(incjets[i]->hasUserFloat("puid:fullDiscriminant"))
	  combinejetPUID.push_back(incjets[i]->userFloat("puid:fullDiscriminant"));
	else
	  combinejetPUID.push_back(incjets[i]->userFloat("pileupJetId:fullDiscriminant"));	
	combinejetPassPUID.push_back(applyPileupJetID(*incjets[i],pileupjetidwp,false));  
	
	// MC based info
	if(isMC){
	  combinejetHFlav.push_back(incjets[i]->hadronFlavour()); 
	  combinejetPFlav.push_back(incjets[i]->partonFlavour()); 
	  if(incjets[i]->genJet()){
	    combinejetGenpt.push_back(incjets[i]->genJet()->pt()); 
	    combinejetGeneta.push_back(incjets[i]->genJet()->eta()); 
	    combinejetGenphi.push_back(incjets[i]->genJet()->phi()); 
	    combinejetGenm.push_back(incjets[i]->genJet()->mass()); 
	  }
	  else{
	    combinejetGenpt.push_back(0.); 
	    combinejetGeneta.push_back(0.); 
	    combinejetGenphi.push_back(0.); 
	    combinejetGenm.push_back(0.); 
	  }
	  // b-tag SF for jets
	  if(addBTagScaleFactor){
	    calculateBtagSF(*incjets[i],"CSV",combinejetBtagSF,combinejetBtagSFUp,combinejetBtagSFDown);
	    calculateBtagSF(*incjets[i],"MVA",combinejetBtagMVASF,combinejetBtagMVASFUp,combinejetBtagMVASFDown);
	  }
	}
      }
    }
    
    // systematics
    for(size_t i = 0; i < incjets_jesup.size(); i++){
      if (incjets_jesup[i]->pt() > minJetPtAK4Store){
	combinejetptup.push_back(incjets_jesup[i]->pt());
	combinejetetaup.push_back(incjets_jesup[i]->eta());
	combinejetphiup.push_back(incjets_jesup[i]->phi());
	combinejetmup.push_back(incjets_jesup[i]->mass());
      }
    }

    for(size_t i = 0; i < incjets_jesdw.size(); i++){
      if (incjets_jesdw[i]->pt() > minJetPtAK4Store){
	combinejetptdw.push_back(incjets_jesdw[i]->pt());
	combinejetetadw.push_back(incjets_jesdw[i]->eta());
	combinejetphidw.push_back(incjets_jesdw[i]->phi());
	combinejetmdw.push_back(incjets_jesdw[i]->mass());
      }
    }
    
    for(size_t i = 0; i < incjets_jer.size(); i++){
      if (incjets_jer[i]->pt() > minJetPtAK4Store){
	combinejetptjer.push_back(incjets_jer[i]->pt());
	combinejetetajer.push_back(incjets_jer[i]->eta());
	combinejetphijer.push_back(incjets_jer[i]->phi());
	combinejetmjer.push_back(incjets_jer[i]->mass());
      }
    }
    
    // delta phi between jets
    if (combinejetphi.size() > 1)
      jetjetdphi = deltaPhi(combinejetphi[0], combinejetphi[1]);
    
    for (size_t i = 0; i < incjets.size(); i++) {
      if (incjets[i]->pt() > minJetPtCountAK4) {
	htinc += incjets[i]->pt(); 
      if (fabs(incjets[i]->eta()) < 3.0) 
	ht30 += incjets[i]->pt();
      }
    }  
    for (size_t i = 0; i < jets.size(); i++) {
      if (jets[i]->pt() > minJetPtCountAK4)
	ht += jets[i]->pt();      
    }
  }
  else{ // Puppi jets 


    combinePuppijetpt.push_back(incPuppijets[i]->pt());
    combinePuppijeteta.push_back(incPuppijets[i]->eta());
    combinePuppijetphi.push_back(incPuppijets[i]->phi());
    combinePuppijetm.push_back(incPuppijets[i]->mass());
    combinePuppijetbtag.push_back(incPuppijets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    combinePuppijetbtagMVA.push_back(incPuppijets[i]->bDiscriminator("pfCombinedMVAV2BJetTags"));
    combinePuppijetCHfrac  .push_back(incPuppijets[i]->chargedHadronEnergyFraction());
    combinePuppijetNHfrac  .push_back(incPuppijets[i]->neutralHadronEnergyFraction());
    combinePuppijetEMfrac  .push_back(incPuppijets[i]->neutralEmEnergyFraction());
    combinePuppijetCEMfrac .push_back(incPuppijets[i]->chargedEmEnergyFraction());
    combinePuppijetmetdphi.push_back(deltaPhi(incPuppijets[i]->phi(), puppit1pfmetphi));
      
    if(incPuppijets[i]->hasUserFloat("QGTaggerPuppi:qgLikelihood"))
      combinePuppijetQGL   .push_back(incPuppijets[i]->userFloat("QGTaggerPuppi:qgLikelihood")); 
      
    // MC based info
    if(isMC){
      combinePuppijetHFlav.push_back(incPuppijets[i]->hadronFlavour()); 
      combinePuppijetPFlav.push_back(incPuppijets[i]->partonFlavour()); 
      if(incPuppijets[i]->genJet()){
	combinePuppijetGenpt.push_back(incPuppijets[i]->genJet()->pt()); 
	combinePuppijetGeneta.push_back(incPuppijets[i]->genJet()->eta()); 
	combinePuppijetGenphi.push_back(incPuppijets[i]->genJet()->phi()); 
	combinePuppijetGenm.push_back(incPuppijets[i]->genJet()->mass()); 
      }  
      else{
	combinePuppijetGenpt.push_back(0.); 
	combinePuppijetGeneta.push_back(0.); 
	combinePuppijetGenphi.push_back(0.); 
	combinePuppijetGenm.push_back(0.); 
      }
          
      // b-tag SF for Puppijets
      if(addBTagScaleFactor){
	calculateBtagSF(*incPuppijets[i],"CSV",combinePuppijetBtagSF,combinePuppijetBtagSFUp,combinePuppijetBtagSFDown);
	calculateBtagSF(*incPuppijets[i],"MVA",combinePuppijetBtagMVASF,combinePuppijetBtagMVASFUp,combinePuppijetBtagMVASFDown);
      }    
    }
  }
}

///////////
for(size_t i = 0; i < incPuppijets_jesup.size(); i++){
  if (incPuppijets_jesup[i]->pt() > minJetPtCountAK4){
      
    combinePuppijetptup.push_back(incPuppijets_jesup[i]->pt());
    combinePuppijetetaup.push_back(incPuppijets_jesup[i]->eta());
    combinePuppijetphiup.push_back(incPuppijets_jesup[i]->phi());
    combinePuppijetmup.push_back(incPuppijets_jesup[i]->mass());
  }
 }

///////////
for(size_t i = 0; i < incPuppijets_jesdw.size(); i++){
  if (incPuppijets_jesdw[i]->pt() > minJetPtCountAK4){
      
    combinePuppijetptdw.push_back(incPuppijets_jesdw[i]->pt());
    combinePuppijetetadw.push_back(incPuppijets_jesdw[i]->eta());
    combinePuppijetphidw.push_back(incPuppijets_jesdw[i]->phi());
    combinePuppijetmdw.push_back(incPuppijets_jesdw[i]->mass());
  }
 }

///////////
for(size_t i = 0; i < incPuppijets_jer.size(); i++){
  if (incPuppijets_jer[i]->pt() > minJetPtCountAK4){
      
    combinePuppijetptjer.push_back(incPuppijets_jer[i]->pt());
    combinePuppijetetajer.push_back(incPuppijets_jer[i]->eta());
    combinePuppijetphijer.push_back(incPuppijets_jer[i]->phi());
    combinePuppijetmjer.push_back(incPuppijets_jer[i]->mass());
  }
 }


  }
  return true;  
}


////////
void JetTreeFiller::DeclareAndSetBranches(){

  if(not isPuppi_){

    tree_->Branch("njets"                , &njets                , "njets/i");
    tree_->Branch("njetsinc"             , &njetsinc             , "njetsinc/i");
    tree_->Branch("nbjets"               , &nbjets               , "nbjets/i");
    tree_->Branch("nbjetslowpt"          , &nbjetslowpt          , "nbjetslowpt/i");

    if(addMETSystematics){
      tree_->Branch("njetsincup"         , &njetsincup           , "njetsincup/i");
      tree_->Branch("njetsincdw"         , &njetsincdw           , "njetsincdw/i");
      tree_->Branch("njetsincjer"        , &njetsincjer          , "njetsincjer/i");
    }

    if(not isTriggerTree and not isPhotonPurity){
      tree_->Branch("nbjetsMVA"            , &nbjetsMVA            , "nbjetsMVA/i");
      tree_->Branch("nbjetsMVAlowpt"       , &nbjetsMVAlowpt       , "nbjetsMVAlowpt/i");
    }
    
    tree_->Branch("combinejetpt",      "std::vector<float>", &combinejetpt);
    tree_->Branch("combinejeteta",     "std::vector<float>", &combinejeteta);
    tree_->Branch("combinejetphi",     "std::vector<float>", &combinejetphi);
    tree_->Branch("combinejetm",       "std::vector<float>", &combinejetm);

    if(not isTriggerTree){
      tree_->Branch("combinejetbtag",    "std::vector<float>", &combinejetbtag);
      tree_->Branch("combinejetbtagMVA", "std::vector<float>", &combinejetbtagMVA);
    }

    tree_->Branch("combinejetCHfrac",  "std::vector<float>", &combinejetCHfrac);
    tree_->Branch("combinejetNHfrac",  "std::vector<float>", &combinejetNHfrac);
    
    if(not isTriggerTree){

      tree_->Branch("combinejetEMfrac",  "std::vector<float>", &combinejetEMfrac);
      tree_->Branch("combinejetCEMfrac", "std::vector<float>", &combinejetCEMfrac);
      tree_->Branch("combinejetPHfrac", "std::vector<float>", &combinejetPHfrac);
      tree_->Branch("combinejetELfrac", "std::vector<float>", &combinejetELfrac);
      tree_->Branch("combinejetMUfrac", "std::vector<float>", &combinejetMUfrac);
      tree_->Branch("combinejetHFHfrac", "std::vector<float>", &combinejetHFHfrac);
      tree_->Branch("combinejetHFEMfrac", "std::vector<float>", &combinejetHFEMfrac);

      tree_->Branch("combinejetCHmult",  "std::vector<unsigned int>", &combinejetCHmult);
      tree_->Branch("combinejetNHmult",  "std::vector<unsigned int>", &combinejetNHmult);
      tree_->Branch("combinejetPHmult", "std::vector<unsigned int>", &combinejetPHmult);
      tree_->Branch("combinejetELmult", "std::vector<unsigned int>", &combinejetELmult);
      tree_->Branch("combinejetMUmult", "std::vector<unsigned int>", &combinejetMUmult);
      tree_->Branch("combinejetHFHmult", "std::vector<unsigned int>", &combinejetHFHmult);
      tree_->Branch("combinejetHFEMmult", "std::vector<unsigned int>", &combinejetHFEMmult);

      tree_->Branch("combinejetHFlav",   "std::vector<float>", &combinejetHFlav);
      tree_->Branch("combinejetPFlav",   "std::vector<float>", &combinejetPFlav);
      tree_->Branch("combinejetQGL",     "std::vector<float>", &combinejetQGL);
      tree_->Branch("combinejetPUID",    "std::vector<float>", &combinejetPUID);
      tree_->Branch("combinejetPassPUIF",    "std::vector<float>", &combinejetPassPUID);
      tree_->Branch("combinejetGenpt",   "std::vector<float>", &combinejetGenpt);
      tree_->Branch("combinejetGeneta",  "std::vector<float>", &combinejetGeneta);
      tree_->Branch("combinejetGenphi",  "std::vector<float>", &combinejetGenphi);
      tree_->Branch("combinejetGenm",    "std::vector<float>", &combinejetGenm);
      tree_->Branch("combinejetBtagSF",  "std::vector<float>", &combinejetBtagSF);
      tree_->Branch("combinejetBtagSFUp", "std::vector<float>", &combinejetBtagSFUp);
      tree_->Branch("combinejetBtagSFDown",    "std::vector<float>", &combinejetBtagSFDown);
      tree_->Branch("combinejetBtagMVASF",     "std::vector<float>", &combinejetBtagMVASF);
      tree_->Branch("combinejetBtagMVASFUp",   "std::vector<float>", &combinejetBtagMVASFUp);
      tree_->Branch("combinejetBtagMVASFDown", "std::vector<float>", &combinejetBtagMVASFDown);    

      tree_->Branch("ht"                   , &ht                   , "ht/F");
      tree_->Branch("htinc"                , &htinc                , "htinc/F");
      tree_->Branch("ht30"                 , &ht30                 , "ht30/F");
      tree_->Branch("jetjetdphi"           , &jetjetdphi           , "jetjetdphi/F");

    }

    if(jetsJESUpTag.label() != "" and jetsJESDwTag.label() != "" and jetsJERTag.label() != ""){
      tree_->Branch("combinejetptup",     "std::vector<float>", &combinejetptup);
      tree_->Branch("combinejetetaup",    "std::vector<float>", &combinejetetaup);
      tree_->Branch("combinejetphiup",    "std::vector<float>", &combinejetphiup);
      tree_->Branch("combinejetmup",      "std::vector<float>", &combinejetmup);

      tree_->Branch("combinejetptdw",     "std::vector<float>", &combinejetptdw);
      tree_->Branch("combinejetetadw",    "std::vector<float>", &combinejetetadw);
      tree_->Branch("combinejetphidw",    "std::vector<float>", &combinejetphidw);
      tree_->Branch("combinejetmdw",      "std::vector<float>", &combinejetmdw);

      tree_->Branch("combinejetptjer",     "std::vector<float>", &combinejetptjer);
      tree_->Branch("combinejetetajer",    "std::vector<float>", &combinejetetajer);
      tree_->Branch("combinejetphijer",    "std::vector<float>", &combinejetphijer);
      tree_->Branch("combinejetmjer",      "std::vector<float>", &combinejetmjer);

    }
  }
  else{

    if(isPuppi_ and not isTriggerTree and not isPhotonPurity and not isQCDTree){

      tree_->Branch("npuppijets"                , &npuppijets                , "npuppijets/i");
      tree_->Branch("npuppijetsinc"             , &npuppijetsinc             , "npuppijetsinc/i");

      if(addMETSystematics){
	tree_->Branch("npuppijetsincup"             , &npuppijetsincup             , "npuppijetsincup/i");
	tree_->Branch("npuppijetsincdw"             , &npuppijetsincdw             , "npuppijetsincdw/i");
	tree_->Branch("npuppijetsincjer"            , &npuppijetsincjer            , "npuppijetsincjer/i");
      }
      
      tree_->Branch("npuppibjets"               , &npuppibjets               , "npuppibjets/i");
      tree_->Branch("npuppibjetslowpt"          , &npuppibjetslowpt          , "npuppibjetslowpt/i");
      tree_->Branch("npuppibjetsMVA"            , &npuppibjetsMVA            , "npuppibjetsMVA/i");
      tree_->Branch("npuppibjetsMVAlowpt"       , &npuppibjetsMVAlowpt       , "npuppibjetsMVAlowpt/i");

      tree_->Branch("combinePuppijetpt",  "std::vector<float>", &combinePuppijetpt);
      tree_->Branch("combinePuppijeteta", "std::vector<float>", &combinePuppijeteta);
      tree_->Branch("combinePuppijetphi", "std::vector<float>", &combinePuppijetphi);
      tree_->Branch("combinePuppijetm",   "std::vector<float>", &combinePuppijetm);
      tree_->Branch("combinePuppijetbtag", "std::vector<float>", &combinePuppijetbtag);
      tree_->Branch("combinePuppijetbtagMVA", "std::vector<float>", &combinePuppijetbtagMVA);
      tree_->Branch("combinePuppijetCHfrac", "std::vector<float>", &combinePuppijetCHfrac);
      tree_->Branch("combinePuppijetNHfrac", "std::vector<float>", &combinePuppijetNHfrac);
      tree_->Branch("combinePuppijetEMfrac", "std::vector<float>", &combinePuppijetEMfrac);
      tree_->Branch("combinePuppijetCEMfrac", "std::vector<float>", &combinePuppijetCEMfrac);
      tree_->Branch("combinePuppijetHFlav", "std::vector<float>", &combinePuppijetHFlav);
      tree_->Branch("combinePuppijetPFlav", "std::vector<float>", &combinePuppijetPFlav);
      tree_->Branch("combinePuppijetQGL",   "std::vector<float>", &combinePuppijetQGL);
      tree_->Branch("combinePuppijetGenpt", "std::vector<float>", &combinePuppijetGenpt);
      tree_->Branch("combinePuppijetGeneta", "std::vector<float>", &combinePuppijetGeneta);
      tree_->Branch("combinePuppijetGenphi", "std::vector<float>", &combinePuppijetGenphi);
      tree_->Branch("combinePuppijetGenm",   "std::vector<float>", &combinePuppijetGenm);
      tree_->Branch("combinePuppijetBtagSF", "std::vector<float>", &combinePuppijetBtagSF);
      tree_->Branch("combinePuppijetBtagSFUp", "std::vector<float>", &combinePuppijetBtagSFUp);
      tree_->Branch("combinePuppijetBtagSFDown", "std::vector<float>", &combinePuppijetBtagSFDown);
      tree_->Branch("combinePuppijetBtagMVASF", "std::vector<float>", &combinePuppijetBtagMVASF);
      tree_->Branch("combinePuppijetBtagMVASFUp", "std::vector<float>", &combinePuppijetBtagMVASFUp);
      tree_->Branch("combinePuppijetBtagMVASFDown", "std::vector<float>", &combinePuppijetBtagMVASFDown);
      tree_->Branch("PuppijetPuppijetdphi"      , &PuppijetPuppijetdphi      , "PuppijetPuppijetdphi/F");
      tree_->Branch("Puppiht"                   , &Puppiht                   , "Puppiht/F");

      if(jetsJESUpTag.label() != "" and jetsJESDwTag.label() != "" and jetsJERTag.label() != ""){
	
	tree_->Branch("combinePuppijetptup",  "std::vector<float>", &combinePuppijetptup);
	tree_->Branch("combinePuppijetetaup", "std::vector<float>", &combinePuppijetetaup);
	tree_->Branch("combinePuppijetphiup", "std::vector<float>", &combinePuppijetphiup);
	tree_->Branch("combinePuppijetmup",   "std::vector<float>", &combinePuppijetmup);

	tree_->Branch("combinePuppijetptdw",  "std::vector<float>", &combinePuppijetptdw);
	tree_->Branch("combinePuppijetetadw", "std::vector<float>", &combinePuppijetetadw);
	tree_->Branch("combinePuppijetphidw", "std::vector<float>", &combinePuppijetphidw);
	tree_->Branch("combinePuppijetmdw",   "std::vector<float>", &combinePuppijetmdw);

	tree_->Branch("combinePuppijetptjer",  "std::vector<float>", &combinePuppijetptjer);
	tree_->Branch("combinePuppijetetajer", "std::vector<float>", &combinePuppijetetajer);
	tree_->Branch("combinePuppijetphijer", "std::vector<float>", &combinePuppijetphijer);
	tree_->Branch("combinePuppijetmjer",   "std::vector<float>", &combinePuppijetmjer);

	tree_->Branch("incPuppijetmetdphimin4up"    , &incPuppijetmetdphimin4up    , "incPuppijetmetdphimin4up/F");
	tree_->Branch("incPuppijetmumetdphimin4up"    , &incPuppijetmumetdphimin4up    , "incPuppijetmumetdphimin4up/F");
	tree_->Branch("incPuppijetelmetdphimin4up"  , &incPuppijetelmetdphimin4up  , "incPuppijetelmetdphimin4up/F");
	tree_->Branch("incPuppijetphmetdphimin4up"  , &incPuppijetphmetdphimin4up  , "incPuppijetphmetdphimin4up/F");
  
	tree_->Branch("incPuppijetmetdphimin4dw"    , &incPuppijetmetdphimin4dw    , "incPuppijetmetdphimin4dw/F");
	tree_->Branch("incPuppijetmumetdphimin4dw"    , &incPuppijetmumetdphimin4dw    , "incPuppijetmumetdphimin4dw/F");
	tree_->Branch("incPuppijetelmetdphimin4dw"  , &incPuppijetelmetdphimin4dw  , "incPuppijetelmetdphimin4dw/F");
	tree_->Branch("incPuppijetphmetdphimin4dw"  , &incPuppijetphmetdphimin4dw  , "incPuppijetphmetdphimin4dw/F");
  
	tree_->Branch("incPuppijetmetdphimin4jer"    , &incPuppijetmetdphimin4jer    , "incPuppijetmetdphimin4jer/F");
	tree_->Branch("incPuppijetmumetdphimin4jer"    , &incPuppijetmumetdphimin4jer    , "incPuppijetmumetdphimin4jer/F");
	tree_->Branch("incPuppijetelmetdphimin4jer"  , &incPuppijetelmetdphimin4jer  , "incPuppijetelmetdphimin4jer/F");
	tree_->Branch("incPuppijetphmetdphimin4jer"  , &incPuppijetphmetdphimin4jer  , "incPuppijetphmetdphimin4jer/F");
 
      }
    }  
  }
}
