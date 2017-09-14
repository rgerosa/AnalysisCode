#include "AnalysisCode/MonoXAnalysis/interface/GenParticleTreeFiller.h"

GenParticleTreeFiller::GenParticleTreeFiller(const edm::ParameterSet & iConfig, edm::ConsumesCollector & iC, TTree* tree):
  isMC          (iConfig.existsAs<bool>("isMC") ? iConfig.getParameter<bool>("isMC") : false),
  useLHEWeights (iConfig.existsAs<bool>("useLHEWeights") ? iConfig.getParameter<bool>("useLHEWeights") : false),
  isSignalSample   (iConfig.existsAs<bool>("isSignalSample") ? iConfig.getParameter<bool>("isSignalSample") : false),
  addGenParticles  (iConfig.existsAs<bool>("addGenParticles") ? iConfig.getParameter<bool>("addGenParticles") : false),
  lheEventTag   (iConfig.getParameter<edm::InputTag>("lheinfo")),
  lheRunTag     (iConfig.getParameter<edm::InputTag>("lheRuninfo")),
  isQCDTree  (iConfig.existsAs<bool>("isQCDTree") ? iConfig.getParameter<bool>("isQCDTree") : false),
  applyDiMuonFilter  (iConfig.existsAs<bool>("applyDiMuonFilter") ? iConfig.getParameter<bool>("applyDiMuonFilter") : false),
  applyDiElectronFilter  (iConfig.existsAs<bool>("applyDiElectronFilter") ? iConfig.getParameter<bool>("applyDiElectronFilter") : false),
  applyPhotonJetsFilter  (iConfig.existsAs<bool>("applyPhotonJetsFilter") ? iConfig.getParameter<bool>("applyPhotonJetsFilter") : false),
  isTriggerTree  (iConfig.existsAs<bool>("isTriggerTree") ? iConfig.getParameter<bool>("isTriggerTree") : false),
  isPhotonPurity(iConfig.existsAs<bool>("isPhotonPurity") ? iConfig.getParameter<bool>("isPhotonPurity") : false){

  if(isMC){
    genevtInfoToken = iC.consumes<GenEventInfoProduct> (iConfig.getParameter<edm::InputTag>("genevt"));
    lheInfoToken    = iC.consumes<LHEEventProduct> (lheEventTag);
    lheRunInfoToken = iC.consumes<LHERunInfoProduct,edm::InRun> (lheRunTag);
    gensToken       = iC.consumes<edm::View<reco::GenParticle> > (iConfig.getParameter<edm::InputTag>("gens"));   
  }

  if(isPhotonPurity)
    photonsPurityToken = iC.consumes<pat::PhotonRefVector> (iConfig.getParameter<edm::InputTag>("photonsPurity"));

  readDMFromGenParticle_ = false;
  tree_ = tree;
  DeclareAndSetBranches();
    
}

/////
void GenParticleTreeFiller::initBranches(){

  wgt = 0;
  wzid          = 0; wzmass        = 0.0; wzpt          = 0.0; wzeta         = 0.0; wzphi         = 0.0;
  l1id          = 0; l1pt          = 0.0; l1eta         = 0.0; l1phi         = 0.0;
  l2id          = 0; l2pt          = 0.0; l2eta         = 0.0; l2phi         = 0.0;
  wzid_h        = 0; wzmass_h      = 0.0; wzpt_h        = 0.0; wzeta_h       = 0.0; wzphi_h       = 0.0;
  topmass       = 0; toppt         = 0.0; topeta        = 0.0; topphi        = 0.0;
  atopmass      = 0; atoppt        = 0.0; atopeta       = 0.0; atopphi       = 0.0;
  q1id          = 0; q1pt          = 0.0; q1eta         = 0.0; q1phi         = 0.0;
  q2id          = 0; q2pt          = 0.0; q2eta         = 0.0; q2phi         = 0.0;
  parid         = 0; parpt         = 0.0; pareta        = 0.0; parphi        = 0.0; parmass       = 0;
  ancid         = 0; ancpt         = 0.0; anceta        = 0.0; ancphi        = 0.0; ancmass       = 0;
  wzmothid      = 0.0;
  isdirect      = 0;
  ismatch       = 0;

  dmmass   = 0.; dmphi   = 0.; dmeta   = 0.; dmpt   = 0.; dmid   = 0;
  dmX1mass = 0.; dmX1phi = 0.; dmX1eta = 0.; dmX1pt = 0.; dmX1id = 0;
  dmX2mass = 0.; dmX2phi = 0.; dmX2eta = 0.; dmX2pt = 0.; dmX2id = 0;

  qcdscalewgt.clear();
  couplingwgt.clear();
  gDMV.clear(); gDMA.clear(); gV.clear(); gA.clear(); couplingwgt.clear(); gTheta.clear();
  

}

/////
bool GenParticleTreeFiller::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace pat;
  using namespace boost::algorithm;

  this->initBranches();

  Handle<GenEventInfoProduct>        genevtInfoH;
  Handle<LHEEventProduct>            lheInfoH;
  Handle<View<GenParticle> >         gensH;
  
  Handle<pat::PhotonRefVector> photonsPurityH;
  pat::PhotonRefVector photonsPurity;
   
  if(isPhotonPurity){
    iEvent.getByToken(photonsPurityToken, photonsPurityH);    
    photonsPurity = *photonsPurityH;
  }

  int hardestPhotonPurityIndex = -1;
  float hardestPhotonPurityPt = 0.0;

  for (size_t i = 0; i < photonsPurity.size(); i++) {
    if (photonsPurity[i]->pt() > hardestPhotonPurityPt) {
      hardestPhotonPurityPt = photonsPurity[i]->pt();
      hardestPhotonPurityIndex = i;  
    }
  }
  

  if(isMC){
    if (useLHEWeights){
      iEvent.getByToken(genevtInfoToken, genevtInfoH);
      iEvent.getByToken(lheInfoToken, lheInfoH);
    }
    if (addGenParticles or isSignalSample)
      iEvent.getByToken(gensToken, gensH);
  }

  ////-- event weight
  if (useLHEWeights && genevtInfoH.isValid())
    wgt = genevtInfoH->weight();    
  else wgt = 1.0;

  ///-----
  // add weights for QCD scale, PDF and couplings
  if(useLHEWeights && lheInfoH.isValid()){

    vector<gen::WeightsInfo> weights = lheInfoH->weights();
    std::vector<std::string> tokens;
    for (size_t i = 0; i < weights.size(); i++) {
      TString weight_name (weights[i].id);
      split(tokens, weights[i].id, is_any_of("_"));
      tokens.erase(std::remove(tokens.begin(), tokens.end(),""), tokens.end());

      if(weight_name.Contains("gdms") and weight_name.Contains("gdmp") and weight_name.Contains("gs") and weight_name.Contains("gp")){ // DMsimp Scalar-PS
	gDMV.push_back(std::stod(std::string(TString(tokens.at(1)).ReplaceAll("p","."))));
	gDMA.push_back(std::stod(std::string(TString(tokens.at(3)).ReplaceAll("p","."))));
	gV.push_back(std::stod(std::string(TString(tokens.at(5)).ReplaceAll("p","."))));
	gA.push_back(std::stod(std::string(TString(tokens.at(7)).ReplaceAll("p","."))));
	couplingwgt.push_back(weights[i].wgt);
      }
      else if(weight_name.Contains("gdmv") and weight_name.Contains("gdma") and weight_name.Contains("gv") and weight_name.Contains("ga")){ // DMSimp V/AV
	gDMV.push_back(std::stod(std::string(TString(tokens.at(1)).ReplaceAll("p","."))));
	gDMA.push_back(std::stod(std::string(TString(tokens.at(3)).ReplaceAll("p","."))));
	gV.push_back(std::stod(std::string(TString(tokens.at(5)).ReplaceAll("p","."))));
	gA.push_back(std::stod(std::string(TString(tokens.at(7)).ReplaceAll("p","."))));
	couplingwgt.push_back(weights[i].wgt);
      }
      else if(weight_name.Contains("sin") and weight_name.Contains("gDM") and weight_name.Contains("gH")){ //SMM  
	gTheta.push_back(std::stod(std::string(TString(tokens.at(1)).ReplaceAll("p","."))));
	gDMV.push_back(std::stod(std::string(TString(tokens.at(3)).ReplaceAll("p","."))));
	gV.push_back(std::stod(std::string(TString(tokens.at(5)).ReplaceAll("p","."))));
	couplingwgt.push_back(weights[i].wgt);
      }
      else if(weight_name.Contains("rwgt")) continue;
      else if(qcdscale.size() != 0){ // qcd scale variations
	if(find(qcdscale.begin(),qcdscale.end(),std::stoi(weights[i].id)) != qcdscale.end())
	  qcdscalewgt.push_back(weights[i].wgt);
      }
      else if(qcdscale.size() == 0 and ((std::stoi(weights[i].id) >=1 and std::stoi(weights[i].id) <= 9) or (std::stoi(weights[i].id) >= 1000 and std::stoi(weights[i].id) <= 1009)))
	qcdscalewgt.push_back(weights[i].wgt);
    }      
  }


  ///--- add gen paritcles
  if (isSignalSample and gensH.isValid() and isMC) {
      
    TLorentzVector dm1vec; 
    TLorentzVector dm2vec; 
    bool foundfirst = false;
    
    // loop on gen particles looking for the DM particles --> then to the mediator
    for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) { 

      bool goodParticle = false;
      if (abs(gens_iter->pdgId()) >= 1000001 and abs(gens_iter->pdgId()) <= 1000039)
	goodParticle = true;
      else if(abs(gens_iter->pdgId()) >= 2000001 and abs(gens_iter->pdgId()) <= 2000015)
	goodParticle = true;
      else if(abs(gens_iter->pdgId()) == 9100012)
	goodParticle = true;
      else if(abs(gens_iter->pdgId()) == 18) // DM particles in MG DMSimp and SMM
	goodParticle = true;

      if(not goodParticle)
	continue;
      
      if(!foundfirst) { // first DM particle
	dm1vec.SetPtEtaPhiM(gens_iter->pt(), gens_iter->eta(), gens_iter->phi(), gens_iter->mass());
	dmX1id = gens_iter->pdgId();
	foundfirst = true;	
	if(readDMFromGenParticle_)
	  sampledmM = gens_iter->mass();
      }
      else{
	dm2vec.SetPtEtaPhiM(gens_iter->pt(), gens_iter->eta(), gens_iter->phi(), gens_iter->mass());
	dmX2id = gens_iter->pdgId();
	break;
      }
    }
    
    dmX1pt   = dm1vec.Pt();
    dmX1eta  = dm1vec.Eta();
    dmX1phi  = dm1vec.Phi();
    dmX1mass = dm1vec.M();

    dmX2pt   = dm2vec.Pt();
    dmX2eta  = dm2vec.Eta();
    dmX2phi  = dm2vec.Phi();
    dmX2mass = dm2vec.M();
      
    TLorentzVector medvec(dm1vec);
    medvec += dm2vec;
    dmpt  = medvec.Pt();
    dmeta = medvec.Eta();
    dmphi = medvec.Phi();
    dmmass = medvec.M();

    if(foundfirst == false){ //not found the DM particles and mediator --> look for Higgs invisible
      
      for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter){
	if(gens_iter->pdgId() != 25) continue;
	if(gens_iter->numberOfDaughters() <= 1) continue;

	dmpt   = gens_iter->pt();
	dmeta  = gens_iter->eta();
	dmphi  = gens_iter->phi();
	dmmass = gens_iter->mass();
	dmid   = gens_iter->pdgId();

	dmX1pt   = gens_iter->daughter(0)->pt();
	dmX1eta  = gens_iter->daughter(0)->eta();
	dmX1phi  = gens_iter->daughter(0)->phi();
	dmX1mass = gens_iter->daughter(0)->mass();
	dmX1id   = gens_iter->daughter(0)->pdgId();
	dmX2pt   = gens_iter->daughter(1)->pt();
	dmX2eta  = gens_iter->daughter(1)->eta();
	dmX2phi  = gens_iter->daughter(1)->phi();
	dmX2mass = gens_iter->daughter(1)->mass();
	dmX2id   = gens_iter->daughter(1)->pdgId();
      }
    }
  }
  
  // dump inportant gen particles
  if(addGenParticles and gensH.isValid() and isMC){
      
    // loop on genParticles (prunedGenParticles) trying to find W/Z decying leptonically or hadronically, top and anti-top quarks
    for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {

      if ( (gens_iter->pdgId() == 23 || abs(gens_iter->pdgId()) == 24) && // Z or W-boson
	   gens_iter->numberOfDaughters() > 1 && // before the decay (more than one daughter)
	   abs(gens_iter->daughter(0)->pdgId()) > 10 && 
	   abs(gens_iter->daughter(0)->pdgId()) < 17)  { // decays into leptons, neutrinos 

	wzid   = gens_iter->pdgId();
	wzmass = gens_iter->mass();
	wzpt   = gens_iter->pt();
	wzeta  = gens_iter->eta();
	wzphi  = gens_iter->phi();

	l1id   = gens_iter->daughter(0)->pdgId();
	l1pt   = gens_iter->daughter(0)->pt();
	l1eta  = gens_iter->daughter(0)->eta();
	l1phi  = gens_iter->daughter(0)->phi();
	
	l2id   = gens_iter->daughter(1)->pdgId();
	l2pt   = gens_iter->daughter(1)->pt();
	l2eta  = gens_iter->daughter(1)->eta();
	l2phi  = gens_iter->daughter(1)->phi();

	// look for the tau --> to be more specific
	if(abs(gens_iter->daughter(0)->pdgId()) == 15){
	  for(size_t ipart = 0; ipart < gens_iter->daughter(0)->numberOfDaughters(); ipart++){
	    if(abs(gens_iter->daughter(0)->daughter(ipart)->pdgId()) == 11 or abs(gens_iter->daughter(0)->daughter(ipart)->pdgId()) == 13) // if there is a muon / electron -> leptonic decay
	      l1id = gens_iter->daughter(0)->daughter(ipart)->pdgId();
	    // otherwise it remains a tau --> hadronic decays
	  }
	}
	  
	if(abs(gens_iter->daughter(1)->pdgId()) == 15){
	  for(size_t ipart = 0; ipart < gens_iter->daughter(1)->numberOfDaughters(); ipart++){
	    if(abs(gens_iter->daughter(1)->daughter(ipart)->pdgId()) == 11 or abs(gens_iter->daughter(1)->daughter(ipart)->pdgId()) == 13) // if there is a muon / electron -> leptonic decay
	      l1id = gens_iter->daughter(1)->daughter(ipart)->pdgId();
	    // otherwise it remains a tau --> hadronic decays
	  }
	}
	wzmt   = sqrt(2.0 * l1pt * l2pt * (1.0 - cos(deltaPhi(l1phi, l2phi))));	  
      }
      else if ( (gens_iter->pdgId() == 23 || abs(gens_iter->pdgId()) == 24) && // Z or W-boson
		gens_iter->numberOfDaughters() > 1 && // before the decay (more than one daughter)
		( (abs(gens_iter->daughter(0)->pdgId()) > 0 && abs(gens_iter->daughter(0)->pdgId()) <= 5) or
		  (abs(gens_iter->daughter(1)->pdgId()) > 0 && abs(gens_iter->daughter(1)->pdgId()) <= 5)))  { // decays into quarks
	  
	wzid_h   = gens_iter->pdgId();
	wzmass_h = gens_iter->mass();
	wzpt_h   = gens_iter->pt();
	wzeta_h  = gens_iter->eta();
	wzphi_h  = gens_iter->phi();
	  
	q1id   = gens_iter->daughter(0)->pdgId();
	q1pt   = gens_iter->daughter(0)->pt();
	q1eta  = gens_iter->daughter(0)->eta();
	q1phi  = gens_iter->daughter(0)->phi();
	  
	q2id   = gens_iter->daughter(1)->pdgId();
	q2pt   = gens_iter->daughter(1)->pt();
	q2eta  = gens_iter->daughter(1)->eta();
	q2phi  = gens_iter->daughter(1)->phi();
	wzmt_h   = sqrt(2.0 * l1pt * l2pt * (1.0 - cos(deltaPhi(l1phi, l2phi))));
      }  
      
      else if(gens_iter->pdgId() == 6 and gens_iter->numberOfDaughters() > 1 and 
	      (((abs(gens_iter->daughter(0)->pdgId()) > 0 and  abs(gens_iter->daughter(0)->pdgId()) <= 5) and abs(gens_iter->daughter(1)->pdgId()) == 24) or
	       ((abs(gens_iter->daughter(1)->pdgId()) > 0 and  abs(gens_iter->daughter(1)->pdgId()) <= 5) and abs(gens_iter->daughter(0)->pdgId()) == 24))){
	  
	topmass = gens_iter->mass();
	toppt   = gens_iter->pt();
	topeta  = gens_iter->eta();
	topphi  = gens_iter->eta();
      }
      else if(gens_iter->pdgId() == -6 and gens_iter->numberOfDaughters() > 1 and 
	      (((abs(gens_iter->daughter(0)->pdgId()) > 0 and  abs(gens_iter->daughter(0)->pdgId()) <= 5) and abs(gens_iter->daughter(1)->pdgId()) == 24) or
	       ((abs(gens_iter->daughter(1)->pdgId()) > 0 and  abs(gens_iter->daughter(1)->pdgId()) <= 5) and abs(gens_iter->daughter(0)->pdgId()) == 24))){
	  
	atopmass = gens_iter->mass();
	atoppt   = gens_iter->pt();
	atopeta  = gens_iter->eta();
	atopphi  = gens_iter->eta();  
      }
    }
    // if a Z/W is not found look for a pair of lepton .. this way with the pdgId is not guaranteed that you catch a Z/W boson and also recover DY production
    if (wzid == 0) {
      for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
	if (gens_iter->isPromptFinalState() || gens_iter->isPromptDecayed()) {
	  if (gens_iter->pdgId() >  10 && gens_iter->pdgId() <  17) {
	    l1id   = gens_iter->pdgId();
	    l1pt   = gens_iter->pt();
	    l1eta  = gens_iter->eta();
	    l1phi  = gens_iter->phi();
	  }
	  if (gens_iter->pdgId() < -10 && gens_iter->pdgId() > -17) {
	    l2id   = gens_iter->pdgId();
	    l2pt   = gens_iter->pt();
	    l2eta  = gens_iter->eta();
	    l2phi  = gens_iter->phi();
	  }
	}
      }
      if (l1id > 0) {
	TLorentzVector l1vec;
	TLorentzVector l2vec;
	l1vec.SetPtEtaPhiM(l1pt, l1eta, l1phi, 0.);
	l2vec.SetPtEtaPhiM(l2pt, l2eta, l2phi, 0.);
	TLorentzVector wzvec(l1vec);
	wzvec += l2vec;
	wzmass = wzvec.M();
	wzpt   = wzvec.Pt();
	wzeta  = wzvec.Eta();
	wzphi  = wzvec.Phi();
	wzmt   = sqrt(2.0 * l1pt * l2pt * (1.0 - cos(deltaPhi(l1phi, l2phi))));
	if (l1id+l2id == 0) wzid = 23;
	else                wzid = 24;
      }
    }

    // no W or Z decay leptonically
    if (wzid == 0) {

      for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) { // loop on prunedGenParticles                                                          
	if (gens_iter->pdgId() == 22 && // photons                                                                                                                           
	    gens_iter->status() == 1 && // final state                                                                                                                       
	    gens_iter->isPromptFinalState() &&
	    gens_iter->pt() > wzpt) {

	  wzid   = gens_iter->pdgId();
	  wzpt   = gens_iter->pt();
	  wzeta  = gens_iter->eta();
	  wzphi  = gens_iter->phi();
	  wzmass = gens_iter->mass();
	  findFirstNonPhotonMother(&(*gens_iter), ancid, ancpt, anceta, ancphi, ancmass);
	  findMother(&(*gens_iter), parid, parpt, pareta, parphi, parmass);        
	}

	if(isPhotonPurity){
	  if (gens_iter->pdgId() == 22 && // photons
	      gens_iter->status() == 1 && // final state
	      gens_iter->isPromptFinalState() &&
	      gens_iter->pt() >= wzpt) {
	          
	    findFirstNonPhotonMother(&(*gens_iter), ancid, ancpt, anceta, ancphi,ancmass);
	    findMother(&(*gens_iter), parid, parpt, pareta, parphi, parmass);
	          
	    if( (abs(ancid) <= 5 || abs(ancid) == 2212) and hardestPhotonPurityIndex >= 0){ 
	      float dR = computeDR(&(*gens_iter),photonsPurity[hardestPhotonPurityIndex] );
	      wzid   = gens_iter->pdgId();
	      wzpt   = gens_iter->pt();
	      wzeta  = gens_iter->eta();
	      wzphi  = gens_iter->phi();
	      wzmass = gens_iter->mass();
	      wzmothid = gens_iter->mother(0)->pdgId();
	      if(dR < 0.3 && fabs((photonsPurity[hardestPhotonPurityIndex]->pt()-gens_iter->pt())/photonsPurity[hardestPhotonPurityIndex]->pt()) < 0.5){
		ismatch=1;
		float dRFrag = sqrt(fabs(anceta-wzeta)*fabs(anceta-wzeta)+deltaPhi(wzphi,ancphi)*deltaPhi(wzphi,ancphi));
		if(dRFrag > 0.4) isdirect = 1;
	      }
	    }
	  }
	}          
      }
    }
  }
  
  return true;  
}

void GenParticleTreeFiller::ReadLHERunProduct(edm::Run const& iRun, float & xsec){

  // info about MC cross section in case the xsec parsed has a dummy value
  if(isMC and useLHEWeights){    
    edm::Handle<LHERunInfoProduct> run;
    iRun.getByLabel(lheRunTag,run);
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());
    if(xsec < 0)
      xsec = myLHERunInfoProduct.heprup().XSECUP.at(0);
    
    using namespace boost::algorithm;

    if(isSignalSample){ // in case of DM signal

      for (auto iter = myLHERunInfoProduct.headers_begin(); iter != myLHERunInfoProduct.headers_end(); iter++){
	std::vector<std::string> lines = iter->lines();
	for (unsigned int iLine = 0; iLine<lines.size(); iLine++) { 
	  std::vector<std::string> tokens;
	  if(lines.at(iLine).find("DMmass") !=std::string::npos){// powheg mono-j
	    split(tokens, lines.at(iLine), is_any_of(" "));
	    tokens.erase(std::remove(tokens.begin(), tokens.end(),""), tokens.end());
	    sampledmM = std::stod(tokens.at(1));
	  }
	  else if(lines.at(iLine).find("DMVmass") !=std::string::npos){// powheg mono-j
	    split(tokens, lines.at(iLine), is_any_of(" "));
	    tokens.erase(std::remove(tokens.begin(), tokens.end(),""), tokens.end());
	    samplemedM = std::stod(tokens.at(1));
	  }
	  else if(lines.at(iLine).find("import model") !=std::string::npos){ // madgraph mono-V                                                                               
	    split(tokens, lines.at(iLine), is_any_of(" "));
	    tokens.erase(std::remove(tokens.begin(), tokens.end(),""), tokens.end());
	    std::vector<std::string> subtokens;
	    split(subtokens,tokens.at(2),is_any_of("_"));
	    if(subtokens.size() >= 5){
	      samplemedM = std::stod(subtokens.at(3));
	      sampledmM = std::stod(subtokens.at(4));
	    }
	    else{
	      samplemedM = std::stod(subtokens.at(1));
	      sampledmM = std::stod(subtokens.at(2));
	    }
	  }
	  else if(lines.at(iLine).find("Resonance:") != std::string::npos){ // JHUGen --> only resonance mass (mediator) .. dM fixed in the event loop                     
	    split(tokens, lines.at(iLine), is_any_of(" "));
	    tokens.erase(std::remove(tokens.begin(), tokens.end(),""), tokens.end());
	    samplemedM = std::stod(tokens.at(3));
	    sampledmM  = -1.;
	    readDMFromGenParticle_ = true;
	  }
	    
	  // read-weights for scale variation                                                                                                                                                         
	  if(lines.at(iLine).find("Central scale variation") != std::string::npos or lines.at(iLine).find("scale_variation") != std::string::npos){
	    for(unsigned int iLine2 = iLine+1; iLine2 < lines.size(); iLine2++){
	      TString line_string (lines.at(iLine2));
	      if(lines.at(iLine2) != "" and line_string.Contains("id=") and not line_string.Contains("</weightgroup>")){
		split(tokens, lines.at(iLine2), is_any_of("\""));
		tokens.erase(std::remove(tokens.begin(), tokens.end(),""), tokens.end());
		qcdscale.push_back(std::stoi(tokens.at(1)));
	      }
	      else if(lines.at(iLine2) != "" and line_string.Contains("</weightgroup>"))
		break;
	    }
	  }  
	}
      }
    }
  }
}


////////
void GenParticleTreeFiller::DeclareAndSetBranches(){

  // W/Z gen-level info: leptonic and hadronic
  if(not isQCDTree){    

    if(not isTriggerTree or (isTriggerTree and isMC)){
      tree_->Branch("wzid"                 , &wzid                 , "wzid/I");
      tree_->Branch("wzmass"               , &wzmass               , "wzmass/F");
      tree_->Branch("wzmt"                 , &wzmt                 , "wzmt/F");
      tree_->Branch("wzpt"                 , &wzpt                 , "wzpt/F");
      tree_->Branch("wzeta"                , &wzeta                , "wzeta/F");
      tree_->Branch("wzphi"                , &wzphi                , "wzphi/F");
      
      tree_->Branch("l1id"                 , &l1id                 , "l1id/I");
      tree_->Branch("l1pt"                 , &l1pt                 , "l1pt/F");
      tree_->Branch("l1eta"                , &l1eta                , "l1eta/F");
      tree_->Branch("l1phi"                , &l1phi                , "l1phi/F");
      tree_->Branch("l2id"                 , &l2id                 , "l2id/I");
      tree_->Branch("l2pt"                 , &l2pt                 , "l2pt/F");
      tree_->Branch("l2eta"                , &l2eta                , "l2eta/F");
      tree_->Branch("l2phi"                , &l2phi                , "l2phi/F");
    }

    if(not isTriggerTree){

      tree_->Branch("wzmothid"             , &wzmothid             , "wzmothid/F");
      tree_->Branch("ismatch"              , &ismatch              , "ismatch/I");
      tree_->Branch("isdirect"             , &isdirect             , "isdirect/I");          
      tree_->Branch("wzid_h"               , &wzid_h               , "wzid_h/I");
      tree_->Branch("wzmass_h"             , &wzmass_h             , "wzmass_h/F");
      tree_->Branch("wzmt_h"               , &wzmt_h               , "wzmt_h/F");
      tree_->Branch("wzpt_h"               , &wzpt_h               , "wzpt_h/F");
      tree_->Branch("wzeta_h"              , &wzeta_h              , "wzeta_h/F");
      tree_->Branch("wzphi_h"              , &wzphi_h              , "wzphi_h/F");
      tree_->Branch("q1id"                 , &q1id                 , "q1id/I");
      tree_->Branch("q1pt"                 , &q1pt                 , "q1pt/F");
      tree_->Branch("q1eta"                , &q1eta                , "q1eta/F");
      tree_->Branch("q1phi"                , &q1phi                , "q1phi/F");
      tree_->Branch("q2id"                 , &q2id                 , "q2id/I");
      tree_->Branch("q2pt"                 , &q2pt                 , "q2pt/F");
      tree_->Branch("q2eta"                , &q2eta                , "q2eta/F");
      tree_->Branch("q2phi"                , &q2phi                , "q2phi/F");
  
      // Top info
      if(not applyDiMuonFilter and not applyDiElectronFilter and not applyPhotonJetsFilter){

	tree_->Branch("topmass"               , &topmass               , "topmass/F");
	tree_->Branch("toppt"                 , &toppt                 , "toppt/F");
	tree_->Branch("topeta"                , &topeta                , "topeta/F");
	tree_->Branch("topphi"                , &topphi                , "topphi/F");
	tree_->Branch("atopmass"               , &atopmass               , "atopmass/F");
	tree_->Branch("atoppt"                 , &atoppt                 , "atoppt/F");
	tree_->Branch("atopeta"                , &atopeta                , "atopeta/F");
	tree_->Branch("atopphi"                , &atopphi                , "atopphi/F");
	
	// photon gen info
	tree_->Branch("parid"                , &parid                , "parid/I");
	tree_->Branch("parpt"                , &parpt                , "parpt/F");
	tree_->Branch("pareta"               , &pareta               , "pareta/F");
	tree_->Branch("parphi"               , &parphi               , "parphi/F");
	tree_->Branch("parmass"              , &parmass              , "parmass/F");
	tree_->Branch("ancid"                , &ancid                , "ancid/I");
	tree_->Branch("ancpt"                , &ancpt                , "ancpt/F");
	tree_->Branch("anceta"               , &anceta               , "anceta/F");
	tree_->Branch("ancphi"               , &ancphi               , "ancphi/F");
	tree_->Branch("ancmass"              , &ancmass              , "ancmass/F");
      }
      
      if(not applyDiMuonFilter and not applyDiElectronFilter and not applyPhotonJetsFilter){

	// DM mediator
	tree_->Branch("dmmass",&dmmass,"dmmass/F");
	tree_->Branch("dmpt",&dmpt,"dmpt/F");
	tree_->Branch("dmeta",&dmeta,"dmeta/F");
	tree_->Branch("dmphi",&dmphi,"dmphi/F");
	tree_->Branch("dmid",&dmid,"dmid/I");
	
	// DM particles
	tree_->Branch("dmX1id",&dmX1id,"dmX1id/I");
	tree_->Branch("dmX1pt",&dmX1pt,"dmX1pt/F");
	tree_->Branch("dmX1eta",&dmX1eta,"dmX1eta/F");
	tree_->Branch("dmX1phi",&dmX1phi,"dmX1phi/F");
	tree_->Branch("dmX1mass",&dmX1mass,"dmX1mass/F");
	
	tree_->Branch("dmX2id",&dmX2id,"dmX2id/I");
	tree_->Branch("dmX2pt",&dmX2pt,"dmX2pt/F");
	tree_->Branch("dmX2eta",&dmX2eta,"dmX2eta/F");
	tree_->Branch("dmX2phi",&dmX2phi,"dmX2phi/F");
	tree_->Branch("dmX2mass",&dmX2mass,"dmX2mass/F");
	
	if(useLHEWeights){
	  tree_->Branch("qcdscalewgt","std::vector<float>",&qcdscalewgt);
	  if(isSignalSample){
	    tree_->Branch("couplingwgt","std::vector<float>",&couplingwgt);
	    tree_->Branch("gDMV","std::vector<float>",&gDMV);
	    tree_->Branch("gTheta","std::vector<float>",&gTheta);
	    tree_->Branch("gDMA","std::vector<float>",&gDMA);
	    tree_->Branch("gV","std::vector<float>",&gV);
	    tree_->Branch("gA","std::vector<float>",&gA);
	  }
	}
	
	// sample info: mediator and DM mass, useful for fast sim                                                                                                                     
	tree_->Branch("samplemedM",   &samplemedM, "samplemedM/F");
	tree_->Branch("sampledmM",    &sampledmM, "sampledmM/F");
      }
    }
  }

}


void GenParticleTreeFiller::findFirstNonPhotonMother(const reco::Candidate *particle, 
						     int & ancestorid, 
						     float & ancestorpt, 
						     float & ancestoreta, 
						     float & ancestorphi, 
						     float & ancestormass) {
  
  if (particle == 0)
    return;
  
  if (abs(particle->pdgId()) == 22) 
    findFirstNonPhotonMother(particle->mother(0), ancestorid, ancestorpt, ancestoreta, ancestorphi,ancestormass);
  else {
    if(particle->pt() <= 0) return;
    ancestorid  = particle->pdgId();
    ancestorpt  = particle->pt();
    ancestoreta = particle->eta();
    ancestorphi = particle->phi();
    ancestormass = particle->mass();
  }
  return;
}

// for photons
void GenParticleTreeFiller::findMother(const reco::Candidate *particle, 
				       int & ancestorid, 
				       float & ancestorpt, 
				       float & ancestoreta, 
				       float & ancestorphi, 
				       float & ancestormass) {
  
  if (particle == 0) 
    return;
  
  if (abs(particle->pdgId()) == 22) {
    if(particle->pt() <= 0) return;
    ancestorid  = particle->pdgId();
    ancestorpt  = particle->pt();
    ancestoreta = particle->eta();
    ancestorphi = particle->phi();
    ancestormass = particle->mass();
  }
  return;
}

float GenParticleTreeFiller::computeDR(const reco::Candidate *genPart,pat::PhotonRef phot){
  float dR = 999.;
  TLorentzVector phop4;
  phop4.SetPtEtaPhiM(phot->pt(),phot->eta(), phot->phi(),phot->mass());
  TLorentzVector p4;
  p4.SetPtEtaPhiM(genPart->pt(),genPart->eta(),genPart->phi(),genPart->mass());
  if(phot->pt() != 0 and genPart->pt() != 0)
    dR = phop4.DeltaR(p4);
  return dR;
}
