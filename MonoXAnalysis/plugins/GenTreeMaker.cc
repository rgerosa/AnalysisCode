// basic C++ headers
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>

// FWCore
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 

// Gen Info
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"

// DataFormats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"

class GenTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns> {

public:
  explicit GenTreeMaker(const edm::ParameterSet&);
  ~GenTreeMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  
  int  sample;
  bool isMiniAOD;
  edm::InputTag lheruntag;  
  edm::EDGetTokenT<LHEEventProduct> lheInfoToken;
  edm::EDGetTokenT<LHERunInfoProduct> lheRunInfoToken;
  edm::EDGetTokenT<GenEventInfoProduct>            genevtInfoToken;
  edm::EDGetTokenT<edm::View<reco::GenParticle> >  gensToken;
  edm::EDGetTokenT<std::vector<reco::GenJet> >     jetsToken;
  edm::EDGetTokenT<edm::View<pat::MET> >           metToken;
  edm::EDGetTokenT<edm::View<reco::GenMET> >       metTrueToken;
  
  TTree* tree;
  
  // event information
  bool     readDMFromGenParticle;
  bool     addEvtInfo;
  bool     isSignalSample;
  edm::InputTag lheRunTag;
  double   minBosonPt;

  uint32_t event, run, lumi;  
  float    xsec, wgt, wgtoriginal;

  std::vector<float>  wgtpdf; 
  std::vector<float>  qcdscalewgt;
  std::vector<int>    qcdscale;

  // boson pdgid, and leptons ids
  int32_t  wzid, mvid, l1id, l2id;
  // information about V-bosons
  float   wzmass, wzpt, wzeta, wzphi;
  float   mvmass, mvpt, mveta, mvphi;
  float   l1pt, l1eta, l1phi, l2pt, l2eta, l2phi;
  // store gen met info
  float   met, metphi;
  // number of jets
  uint32_t njets, njetsinc;
  // store gen jet informations
  std::vector<float> jetpt, jeteta, jetphi, jetmass;
  // DM mediator info
  float sampledmM, samplemedM;
  float dmmass,dmpt,dmeta,dmphi,dmX1pt,dmX1eta,dmX1phi,dmX1mass,dmX2pt,dmX2eta,dmX2phi,dmX2mass;
  int   dmid,dmX1id,dmX2id;
  std::vector<float> couplingwgt, gDMV, gDMA, gV, gA, gTheta;


  template<typename T> 
  class PtSorter{
  public:
    bool operator ()(const T & i, const T & j) const {
      return (i->pt() > j->pt());
    }
    
  };
  
  PtSorter<reco::GenJetRef> jetSorter;
  
};


GenTreeMaker::GenTreeMaker(const edm::ParameterSet& iConfig) {
  isMiniAOD       = iConfig.getParameter<bool>("isMiniAOD");
  lheruntag       = iConfig.getParameter<edm::InputTag>("lherun");
  addEvtInfo      = (iConfig.existsAs<bool>("addEventInfo")    ? iConfig.getParameter<bool>("addEventInfo")  : false);
  isSignalSample  = (iConfig.existsAs<bool>("isSignalSample")  ? iConfig.getParameter<bool>("isSignalSample")  : false);
  xsec            = (iConfig.existsAs<double>("xsec")        ? iConfig.getParameter<double>("xsec") * 1000.0 : -1000.);
  sample          = (iConfig.existsAs<int>("sample")         ? iConfig.getParameter<int>("sample")           : -1);    
  jetsToken       = consumes<std::vector<reco::GenJet> >   (iConfig.getParameter<edm::InputTag>("jets"));
  if(isMiniAOD)
    metToken        = consumes<edm::View<pat::MET> >         (iConfig.getParameter<edm::InputTag>("met"));
  else
    metTrueToken    = consumes<edm::View<reco::GenMET> >     (iConfig.getParameter<edm::InputTag>("met"));

  lheInfoToken    = consumes<LHEEventProduct>              (iConfig.getParameter<edm::InputTag>("lheevt")); 
  lheRunTag       = iConfig.getParameter<edm::InputTag>("lherun");
  lheRunInfoToken = consumes<LHERunInfoProduct,edm::InRun> (iConfig.getParameter<edm::InputTag>("lherun"));
  genevtInfoToken = consumes<GenEventInfoProduct>          (iConfig.getParameter<edm::InputTag>("genevt"));
  gensToken       = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("gens"));   
  minBosonPt      = (iConfig.existsAs<double>("minBosonPt")? iConfig.getParameter<double>("minBosonPt") : 100.);
  usesResource();
  usesResource("TFileService");    
}


GenTreeMaker::~GenTreeMaker() {}

void GenTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;
  using namespace boost::algorithm;
    
  Handle<LHEEventProduct> lheInfoH;
  Handle<GenEventInfoProduct>  genevtInfoH;
  Handle<View<GenParticle> >   gensH;
  Handle<vector<GenJet> >      jetsH;
  Handle<View<pat::MET> >      metH;
  Handle<View<reco::GenMET> >  metTrueH;
  
  iEvent.getByToken(lheInfoToken, lheInfoH);
  iEvent.getByToken(genevtInfoToken, genevtInfoH);
  iEvent.getByToken(gensToken      , gensH);
  iEvent.getByToken(jetsToken      , jetsH);
  if(isMiniAOD)
    iEvent.getByToken(metToken       , metH);
  else
    iEvent.getByToken(metTrueToken   , metTrueH);
  
  event = iEvent.id().event();
  run   = iEvent.id().run();
  lumi  = iEvent.luminosityBlock();
  
  if(genevtInfoH.isValid()) 
    wgt = genevtInfoH->weight();    
  else 
    wgt = 1.0;
  if(lheInfoH.isValid()) 
    wgtoriginal = lheInfoH->originalXWGTUP();
  else 
    wgtoriginal = 1.0;
  
  if(metH.isValid()) {      
    met    = metH->front().genMET()->pt();
    metphi = metH->front().genMET()->phi();
  }
  else {
    met    = -99.;
    metphi = -99.;
  }

  if(metTrueH.isValid()){
    met = metTrueH->front().pt();
    metphi = metTrueH->front().phi();
  }
  else {
    met    = -99.;
    metphi = -99.;
  }
  

  // PDF and scale weights
  wgtpdf.clear();
  qcdscalewgt.clear();
  couplingwgt.clear();
  gDMV.clear(); 
  gDMA.clear(); 
  gV.clear(); 
  gA.clear(); 
  couplingwgt.clear();

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
    else if(qcdscale.size() != 0){ // qcd scale variations                                                                                                                                         
      if(weight_name.Contains("rwgt")) continue;
      if(find(qcdscale.begin(),qcdscale.end(),std::stoi(weights[i].id)) != qcdscale.end())
	qcdscalewgt.push_back(weights[i].wgt);
    }
    else if(qcdscale.size() == 0){
      if(weight_name.Contains("rwgt")) continue;
      else if((std::stoi(weights[i].id) >=1 and std::stoi(weights[i].id) <= 9) or (std::stoi(weights[i].id) >= 1000 and std::stoi(weights[i].id) <= 1009))
	qcdscalewgt.push_back(weights[i].wgt);
    }
  }
  
  wzid = 0; wzmass = -99; wzpt  = -99; wzeta = -99; wzphi = -99;
  mvid = 0; mvmass = -99; mvpt  = -99; mveta = -99; mvphi = -99;
  l1id = 0; l1pt   = -99; l1eta = -99; l1phi = -99;
  l2id = 0; l2pt   = -99; l2eta = -99; l2phi = -99;

  dmmass   = 0.; dmphi   = 0.; dmeta   = 0.; dmpt   = 0.; dmid   = 0;
  dmX1mass = 0.; dmX1phi = 0.; dmX1eta = 0.; dmX1pt = 0.; dmX1id = 0;
  dmX2mass = 0.; dmX2phi = 0.; dmX2eta = 0.; dmX2pt = 0.; dmX2id = 0;

  if(gensH.isValid() && isSignalSample){

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

	if(readDMFromGenParticle)
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
  
  else  if (gensH.isValid() && (sample == 23 || sample == 24)) { // Z+jets or W+jets
    for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
      
      if (abs(gens_iter->pdgId()) == sample && gens_iter->status() == 22) { // take some specific status particles
	mvid   = gens_iter->pdgId();
	mvmass = gens_iter->mass();
	mvpt   = gens_iter->pt();
	mveta  = gens_iter->eta();
	mvphi  = gens_iter->phi();
      }
	  
      if (abs(gens_iter->pdgId()) == sample && gens_iter->numberOfDaughters() > 1 && abs(gens_iter->daughter(0)->pdgId()) > 10 && abs(gens_iter->daughter(0)->pdgId()) < 17) { // take the gen particle with that specific pdgId, decayed in more than 1 daugheter, one of the two to be a lepton/parton 
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
      } 
    }
    
    // if id is not valid --> do a trick for madgraph. A V-boson is not a boson with the right pdgID if outside the BW cutoff
    if (wzid == 0) {      
      float l1mass = 0.;
      float l2mass = 0.;
      for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
	if (gens_iter->isPromptFinalState() || gens_iter->isPromptDecayed()) {
	  if (gens_iter->pdgId() >  10 && gens_iter->pdgId() <  17) { //lepton
	    l1id   = gens_iter->pdgId();
	    l1pt   = gens_iter->pt();
	    l1eta  = gens_iter->eta();
	    l1phi  = gens_iter->phi();
	    l1mass = gens_iter->mass();
	  }
	  if (gens_iter->pdgId() < -10 && gens_iter->pdgId() > -17) { //anti-lepton
	    l2id   = gens_iter->pdgId();
	    l2pt   = gens_iter->pt();
	    l2eta  = gens_iter->eta();
	    l2phi  = gens_iter->phi();
	    l2mass = gens_iter->mass();
	  }
	}
      }
      if (l1id > 0 && ( (l2id == -l1id) || abs(abs(l1id) - abs(l2id)) == 1)) {
	TLorentzVector l1vec;
	TLorentzVector l2vec;
	l1vec.SetPtEtaPhiM(l1pt, l1eta, l1phi, l1mass);
	l2vec.SetPtEtaPhiM(l2pt, l2eta, l2phi, l2mass);
	TLorentzVector wzvec(l1vec);
	wzvec += l2vec;
	wzmass = wzvec.M();
	wzpt   = wzvec.Pt();
	wzeta  = wzvec.Eta();
	wzphi  = wzvec.Phi();
	if (sample == 23 and l2id == -l1id) wzid = 23;
	else if (sample == 24 and (l1id == 11 || l1id == 13 || l1id == 15)) wzid = 24;
	else if (sample == 24 and (l1id == -11 || l1id == -13 || l1id == -15)) wzid = -24;
      }
    }
  }
    
  // photon +jets
  if (gensH.isValid() && sample == 22) {
    for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {
      if (gens_iter->pdgId() == sample && gens_iter->status() == 22) {
	mvid   = gens_iter->pdgId();
	mvmass = gens_iter->mass();
	mvpt   = gens_iter->pt();
	mveta  = gens_iter->eta();
	mvphi  = gens_iter->phi();
      }
      
      if (gens_iter->pdgId() == sample && gens_iter->status() == 1 && gens_iter->isPromptFinalState() && gens_iter->pt() > wzpt) {
	wzid   = gens_iter->pdgId();
	wzpt   = gens_iter->pt();
	wzeta  = gens_iter->eta();
	wzphi  = gens_iter->phi();
      }
    }
  }
  

  // make bosonPt slections
  if(isSignalSample and dmpt < minBosonPt) return;
  else if(not isSignalSample and wzpt < minBosonPt) return;

  // gen jets
  vector<GenJetRef> jets;
  if (jetsH.isValid()) {
    for (auto jets_iter = jetsH->begin(); jets_iter != jetsH->end(); ++jets_iter) {
      GenJetRef jetref(jetsH, jets_iter - jetsH->begin());
      if (jetref.isAvailable() and jetref.isNonnull()) {
	// remove overlap with boson
	if (sample == 22 && wzid == 22 && deltaR(jets_iter->eta(), jets_iter->phi(), wzeta, wzphi) < 0.4) continue;
	// remove overlpas with leptons
	if ((sample == 23 || sample == 24) && (abs(l1id) == 11 || abs(l1id) == 13 || abs(l1id) == 15) && deltaR(jets_iter->eta(), jets_iter->phi(), l1eta, l1phi) < 0.4) continue;
	if ((sample == 23 || sample == 24) && (abs(l2id) == 11 || abs(l2id) == 13 || abs(l2id) == 15) && deltaR(jets_iter->eta(), jets_iter->phi(), l2eta, l2phi) < 0.4) continue;
	jets.push_back(jetref);
      }
    }
  }

  // sort in pt
  if(jets.size() > 0) sort(jets.begin(), jets.end(), jetSorter);
    
  jetpt  .clear();
  jeteta .clear();
  jetphi .clear();
  jetmass.clear();
  
  njets    = 0;
  njetsinc = 0;
  
  for(size_t i = 0; i < jets.size(); i++){
    if(jets[i]->pt() > 30.) {
      njetsinc++;
      if (fabs(jets[i]->eta()) < 2.5) njets++;
      jetpt  .push_back(jets[i]->pt()  );
      jeteta .push_back(jets[i]->eta() );
      jetphi .push_back(jets[i]->phi() );
      jetmass.push_back(jets[i]->mass());
    }
  }
  tree->Fill();    
}    

void GenTreeMaker::beginJob() {

  ///
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree","tree");
  
  if (addEvtInfo) { 
    tree->Branch("event"                , &event                , "event/i");
    tree->Branch("run"                  , &run                  , "run/i");
    tree->Branch("lumi"                 , &lumi                 , "lumi/i");
  }

  tree->Branch("xsec"                 , &xsec                 , "xsec/F");
  tree->Branch("wgt"                  , &wgt                  , "wgt/F");
  tree->Branch("wgtoriginal"          , &wgtoriginal          , "wgtoriginal/F");
  tree->Branch("wgtpdf"               , "std::vector<float>",  &wgtpdf);
  tree->Branch("wgtqcd"               , "std::vector<float>",  &qcdscalewgt);
  
  tree->Branch("njets"                , &njets                , "njets/i");
  tree->Branch("njetsinc"             , &njetsinc             , "njetsinc/i");
  tree->Branch("met"                  , &met                  , "met/F");
  tree->Branch("metphi"               , &metphi               , "metphi/F");

  if(not isSignalSample){
    tree->Branch("mvid"                 , &mvid                 , "mvid/I");
    tree->Branch("mvmass"               , &mvmass               , "mvmass/F");
    tree->Branch("mvpt"                 , &mvpt                 , "mvpt/F");
    tree->Branch("mveta"                , &mveta                , "mveta/F");
    tree->Branch("mvphi"                , &mvphi                , "mvphi/F");
    tree->Branch("wzid"                 , &wzid                 , "wzid/I");
    tree->Branch("wzmass"               , &wzmass               , "wzmass/F");
    tree->Branch("wzpt"                 , &wzpt                 , "wzpt/F");
    tree->Branch("wzeta"                , &wzeta                , "wzeta/F");
    tree->Branch("wzphi"                , &wzphi                , "wzphi/F");
    tree->Branch("l1id"                 , &l1id                 , "l1id/I");
    tree->Branch("l1pt"                 , &l1pt                 , "l1pt/F");
    tree->Branch("l1eta"                , &l1eta                , "l1eta/F");
    tree->Branch("l1phi"                , &l1phi                , "l1phi/F");
    tree->Branch("l2id"                 , &l2id                 , "l2id/I");
    tree->Branch("l2pt"                 , &l2pt                 , "l2pt/F");
    tree->Branch("l2eta"                , &l2eta                , "l2eta/F");
    tree->Branch("l2phi"                , &l2phi                , "l2phi/F");
  }

  tree->Branch("jetpt"                , "std::vector<float>" , &jetpt);
  tree->Branch("jeteta"               , "std::vector<float>" , &jeteta);
  tree->Branch("jetphi"               , "std::vector<float>" , &jetphi);
  tree->Branch("jetmass"              , "std::vector<float>" , &jetmass);

  if(isSignalSample){
    tree->Branch("sampledmM",&sampledmM,"sampledmM/F");
    tree->Branch("samplemedM",&samplemedM,"samplemedM/F");
    tree->Branch("couplingwgt","std::vector<float>",&couplingwgt);
    tree->Branch("gDMV","std::vector<float>",&gDMV);
    tree->Branch("gTheta","std::vector<float>",&gTheta);
    tree->Branch("gDMA","std::vector<float>",&gDMA);
    tree->Branch("gV","std::vector<float>",&gV);
    tree->Branch("gA","std::vector<float>",&gA);     
    tree->Branch("dmmass",&dmmass,"dmmass/F");
    tree->Branch("dmpt",&dmpt,"dmpt/F");
    tree->Branch("dmeta",&dmeta,"dmeta/F");
    tree->Branch("dmphi",&dmphi,"dmphi/F");
    tree->Branch("dmid",&dmid,"dmid/I");

    // DM particles                                                                                                                                                                                   
    tree->Branch("dmX1id",&dmX1id,"dmX1id/I");
    tree->Branch("dmX1pt",&dmX1pt,"dmX1pt/F");
    tree->Branch("dmX1eta",&dmX1eta,"dmX1eta/F");
    tree->Branch("dmX1phi",&dmX1phi,"dmX1phi/F");
    tree->Branch("dmX1mass",&dmX1mass,"dmX1mass/F");

    tree->Branch("dmX2id",&dmX2id,"dmX2id/I");
    tree->Branch("dmX2pt",&dmX2pt,"dmX2pt/F");
    tree->Branch("dmX2eta",&dmX2eta,"dmX2eta/F");
    tree->Branch("dmX2phi",&dmX2phi,"dmX2phi/F");
    tree->Branch("dmX2mass",&dmX2mass,"dmX2mass/F");

  }

}

void GenTreeMaker::endJob() {
}

void GenTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {

  using namespace boost::algorithm;
  using namespace edm;

  edm::Handle<LHERunInfoProduct> run;
  iRun.getByLabel(lheRunTag,run);
  LHERunInfoProduct myLHERunInfoProduct = *(run.product());
  
  if(isSignalSample){
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
	else if(lines.at(iLine).find("DMSmass") !=std::string::npos){// powheg mono-j                                                                                                                
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
	  readDMFromGenParticle = true;
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

void GenTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}


void GenTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(GenTreeMaker);

