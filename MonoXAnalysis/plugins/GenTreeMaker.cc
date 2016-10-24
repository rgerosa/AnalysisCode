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
  
  edm::EDGetTokenT<GenEventInfoProduct>            genevtInfoToken;
  edm::EDGetTokenT<edm::View<reco::GenParticle> >  gensToken;
  edm::EDGetTokenT<std::vector<reco::GenJet> >     jetsToken;
  edm::EDGetTokenT<edm::View<reco::GenMET> >       metToken;
  
  TTree* tree;
  
  // event information
  uint32_t event, run, lumi;  
  bool     addEvtInfo;
  double   wgt, xsec;
  
  // boson pdgid, and leptons ids
  int32_t  wzid, mvid, l1id, l2id;
  // information about V-bosons
  double   wzmass, wzpt, wzeta, wzphi, mvmass, mvpt, mveta, mvphi, l1pt, l1eta, l1phi, l2pt, l2eta, l2phi;
  // store gen met info
  double    met, metphi;
  // store gen jet informations
  std::vector<double> jetpt, jeteta, jetphi, jetmass;
  // number of jets
  int njets, njetsinc;
    
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
    addEvtInfo      = (iConfig.existsAs<bool>("addEventInfo")  ? iConfig.getParameter<double>("addEventInfo")  : false);
    xsec            = (iConfig.existsAs<double>("xsec")        ? iConfig.getParameter<double>("xsec") * 1000.0 : -1000.);
    sample          = (iConfig.existsAs<int>("sample")         ? iConfig.getParameter<int>("sample")           : -1);    
    usesResource();
    usesResource("TFileService");    
    jetsToken       = consumes<std::vector<reco::GenJet> >   (iConfig.getParameter<edm::InputTag>("jets"));
    metToken        = consumes<edm::View<reco::GenMET> >     (iConfig.getParameter<edm::InputTag>("met"));
    genevtInfoToken = consumes<GenEventInfoProduct>          (iConfig.getParameter<edm::InputTag>("genevt"));
    gensToken       = consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("gens"));   
}


GenTreeMaker::~GenTreeMaker() {}

void GenTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    using namespace edm;
    using namespace reco;
    using namespace std;
    
    Handle<GenEventInfoProduct>  genevtInfoH;
    Handle<View<GenParticle> >   gensH;
    Handle<vector<GenJet> >      jetsH;
    Handle<View<GenMET> >        metH;
    
    iEvent.getByToken(genevtInfoToken, genevtInfoH);
    iEvent.getByToken(gensToken      , gensH);
    iEvent.getByToken(jetsToken      , jetsH);
    iEvent.getByToken(metToken       , metH);

    event = iEvent.id().event();
    run   = iEvent.id().run();
    lumi  = iEvent.luminosityBlock();
    
    if (genevtInfoH.isValid()) wgt = genevtInfoH->weight();    
    else wgt = 1.0;
    

    if(metH.isValid()) {      
        met    = metH->front().pt();
        metphi = metH->front().phi();
    }
    else {
        met = -1.;
        metphi = -1.;
    }
      
    wzid = 0; wzmass = 0.0; wzpt  = 0.0; wzeta = 0.0; wzphi = 0.0;
    mvid = 0; mvmass = 0.0; mvpt  = 0.0; mveta = 0.0; mvphi = 0.0;
    l1id = 0; l1pt   = 0.0; l1eta = 0.0; l1phi = 0.0;
    l2id = 0; l2pt   = 0.0; l2eta = 0.0; l2phi = 0.0;

    if (gensH.isValid() && (sample == 23 || sample == 24)) { // Z+jets or W+jets
        for (auto gens_iter = gensH->begin(); gens_iter != gensH->end(); ++gens_iter) {

            /*
            if (abs(gens_iter->pdgId()) == sample) {
                cout << "Status : " << gens_iter->status() << endl; 
                cout << "pT     : " << gens_iter->pt()     << endl; 
                cout << "eta    : " << gens_iter->eta()    << endl; 
                cout << "phi    : " << gens_iter->phi()    << endl; 
                cout << "-----------------------------"    << endl; 
            }
            */
	  
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

	  double l1mass = 0.;
	  double l2mass = 0.;
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
	    if (l2id == -l1id) wzid = 23;
	    else if (l1id == 11 || l1id == 13 || l1id == 15) wzid = 24;
	    else wzid = -24;
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
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("tree"       , "tree");
  
  if (addEvtInfo) { 
    tree->Branch("event"                , &event                , "event/i");
    tree->Branch("run"                  , &run                  , "run/i");
    tree->Branch("lumi"                 , &lumi                 , "lumi/i");
  }
  tree->Branch("xsec"                 , &xsec                 , "xsec/D");
  tree->Branch("wgt"                  , &wgt                  , "wgt/D");
  tree->Branch("njets"                , &njets                , "njets/i");
  tree->Branch("njetsinc"             , &njetsinc             , "njetsinc/i");
  tree->Branch("met"                  , &met                  , "met/D");
  tree->Branch("metphi"               , &metphi               , "metphi/D");
  tree->Branch("mvid"                 , &mvid                 , "mvid/I");
  tree->Branch("mvmass"               , &mvmass               , "mvmass/D");
  tree->Branch("mvpt"                 , &mvpt                 , "mvpt/D");
  tree->Branch("mveta"                , &mveta                , "mveta/D");
  tree->Branch("mvphi"                , &mvphi                , "mvphi/D");
  tree->Branch("wzid"                 , &wzid                 , "wzid/I");
  tree->Branch("wzmass"               , &wzmass               , "wzmass/D");
  tree->Branch("wzpt"                 , &wzpt                 , "wzpt/D");
  tree->Branch("wzeta"                , &wzeta                , "wzeta/D");
  tree->Branch("wzphi"                , &wzphi                , "wzphi/D");
  tree->Branch("l1id"                 , &l1id                 , "l1id/I");
  tree->Branch("l1pt"                 , &l1pt                 , "l1pt/D");
  tree->Branch("l1eta"                , &l1eta                , "l1eta/D");
  tree->Branch("l1phi"                , &l1phi                , "l1phi/D");
  tree->Branch("l2id"                 , &l2id                 , "l2id/I");
  tree->Branch("l2pt"                 , &l2pt                 , "l2pt/D");
  tree->Branch("l2eta"                , &l2eta                , "l2eta/D");
  tree->Branch("l2phi"                , &l2phi                , "l2phi/D");
  tree->Branch("jetpt"                , "std::vector<double>" , &jetpt);
  tree->Branch("jeteta"               , "std::vector<double>" , &jeteta);
  tree->Branch("jetphi"               , "std::vector<double>" , &jetphi);
  tree->Branch("jetmass"              , "std::vector<double>" , &jetmass);
}

void GenTreeMaker::endJob() {
}

void GenTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
}

void GenTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}


void GenTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(GenTreeMaker);

