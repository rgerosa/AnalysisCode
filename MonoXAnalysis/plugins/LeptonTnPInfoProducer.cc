#include <memory>
#include <vector>
#include <iostream>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

class LeptonTnPInfoProducer : public edm::stream::EDProducer<> {
public:
  explicit LeptonTnPInfoProducer(const edm::ParameterSet&);
  ~LeptonTnPInfoProducer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&) ;
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  const edm::EDGetTokenT<GenEventInfoProduct>         geninfoToken;
  const edm::EDGetTokenT<std::vector<reco::Vertex> >  verticesToken;
  const edm::EDGetTokenT<std::vector<pat::Muon> >     muonsToken;
  const edm::EDGetTokenT<std::vector<pat::Electron> > electronsToken;
  const edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken;
  const edm::EDGetTokenT<edm::TriggerResults>  triggerResultsToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken;

  const double loosemuisocut;
  const double tightmuisocut;
  const double tagmuonptcut;
  const double tagmuonetacut;
  const double tagmuontrigmatchdR;
  const bool   requiremuonhlt;
  const std::vector<std::string> tagmuontriggers;
  const double tagelectronptcut;
  const double tagelectronetacut;
  const double tagelectrontrigmatchdR;
  const bool   requireelectronhlt;
  const std::vector<std::string> tagelectrontriggers;
};

LeptonTnPInfoProducer::LeptonTnPInfoProducer(const edm::ParameterSet& iConfig): 
    geninfoToken(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("geninfo"))),
    verticesToken(consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonsToken(consumes<std::vector<pat::Muon> >(iConfig.getParameter<edm::InputTag>("muons"))),
    electronsToken(consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("electrons"))),
    triggerObjectsToken(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects"))),
    triggerResultsToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
    electronVetoIdMapToken(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronvetoid"))),
    electronLooseIdMapToken(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronlooseid"))),
    electronMediumIdMapToken(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronmediumid"))),
    electronTightIdMapToken(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electrontightid"))),
    loosemuisocut(iConfig.getParameter<double>("loosemuisocut")),
    tightmuisocut(iConfig.getParameter<double>("tightmuisocut")),
    tagmuonptcut(iConfig.getParameter<double>("tagmuonptcut")),
    tagmuonetacut(iConfig.getParameter<double>("tagmuonetacut")),
    tagmuontrigmatchdR(iConfig.getParameter<double>("tagmuontrigmatchdR")),
    requiremuonhlt(iConfig.getParameter<bool>("requiremuonhlt")),
    tagmuontriggers(iConfig.getParameter<std::vector<std::string> >("tagmuontriggers")),
    tagelectronptcut(iConfig.getParameter<double>("tagelectronptcut")),
    tagelectronetacut(iConfig.getParameter<double>("tagelectronetacut")),
    tagelectrontrigmatchdR(iConfig.getParameter<double>("tagelectrontrigmatchdR")),
    requireelectronhlt(iConfig.getParameter<bool>("requireelectronhlt")),
    tagelectrontriggers(iConfig.getParameter<std::vector<std::string> >("tagelectrontriggers"))
{
  // produce a map with nvtx for a given muon
  produces<edm::ValueMap<float> >("munvtxmap");
  // produce a map with the generator weight
  produces<edm::ValueMap<float> >("muwgtmap");

  // produce a ref vector with tag muon trigger info
  produces<pat::MuonRefVector>("hltmu20muonrefs");
  produces<pat::MuonRefVector>("hlttkmu20muonrefs");
  produces<pat::MuonRefVector>("hltmu22muonrefs");
  produces<pat::MuonRefVector>("hlttkmu22muonrefs");
  produces<pat::MuonRefVector>("hltmu24muonrefs");
  produces<pat::MuonRefVector>("hlttkmu24muonrefs");
  produces<pat::MuonRefVector>("hltmumuonrefs");
  produces<pat::MuonRefVector>("hlttkmumuonrefs");

  // produces collection for loose and tight muons
  produces<pat::MuonRefVector>("loosemuonrefs");
  produces<pat::MuonRefVector>("tightmuonrefs");
  produces<pat::MuonCollection>("tightmuons");

  /// same logic also for electrons
  produces<edm::ValueMap<float> >("elnvtxmap");
  produces<edm::ValueMap<float> >("elwgtmap");

  // Map for electron trigger
  produces<pat::ElectronRefVector>("hltele24eta2p1wplooseelectronrefs");
  produces<pat::ElectronRefVector>("hltele25eta2p1wptightelectronrefs");
  produces<pat::ElectronRefVector>("hltele27eta2p1wplooseelectronrefs");
  produces<pat::ElectronRefVector>("hltele27eta2p1wptightelectronrefs");
  produces<pat::ElectronRefVector>("hltele27wptightelectronrefs");
  produces<pat::ElectronRefVector>("hltele105electronrefs");
  produces<pat::ElectronRefVector>("hltele115electronrefs");
  produces<pat::ElectronRefVector>("hltelelectronrefs");

  /// Map for electron ID studies
  produces<pat::ElectronRefVector>("vetoelectronrefs");
  produces<pat::ElectronRefVector>("looseelectronrefs");
  produces<pat::ElectronRefVector>("mediumelectronrefs");
  produces<pat::ElectronRefVector>("tightelectronrefs");
  produces<pat::ElectronCollection>("tightelectrons");
}


LeptonTnPInfoProducer::~LeptonTnPInfoProducer() {
}

void LeptonTnPInfoProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;
  
  Handle<GenEventInfoProduct> geninfoH;
  iEvent.getByToken(geninfoToken, geninfoH);
  
  Handle<std::vector<reco::Vertex> > verticesH;
  iEvent.getByToken(verticesToken, verticesH);
  
  Handle<std::vector<pat::Muon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);
  
  Handle<std::vector<pat::Electron> > electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  
  Handle<pat::TriggerObjectStandAloneCollection> triggerObjectsH;
  iEvent.getByToken(triggerObjectsToken, triggerObjectsH);
  
  Handle<edm::TriggerResults> triggerResultsH;
  iEvent.getByToken(triggerResultsToken, triggerResultsH);
  
  Handle<edm::ValueMap<bool> > electronVetoIdH;
  iEvent.getByToken(electronVetoIdMapToken, electronVetoIdH);
  
  Handle<edm::ValueMap<bool> > electronLooseIdH;
  iEvent.getByToken(electronLooseIdMapToken, electronLooseIdH);
  
  Handle<edm::ValueMap<bool> > electronMediumIdH;
  iEvent.getByToken(electronMediumIdMapToken, electronMediumIdH);
  
  Handle<edm::ValueMap<bool> > electronTightIdH;
  iEvent.getByToken(electronTightIdMapToken, electronTightIdH);

  // output collection
  std::auto_ptr<edm::ValueMap<float> > outputmunvtxmap(new ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > outputmuwgtmap(new ValueMap<float>());
  // trigger muons
  std::auto_ptr<pat::MuonRefVector> outputhltmu20muonrefs(new pat::MuonRefVector);
  std::auto_ptr<pat::MuonRefVector> outputhlttkmu20muonrefs(new pat::MuonRefVector);
  std::auto_ptr<pat::MuonRefVector> outputhltmu22muonrefs(new pat::MuonRefVector);
  std::auto_ptr<pat::MuonRefVector> outputhlttkmu22muonrefs(new pat::MuonRefVector);
  std::auto_ptr<pat::MuonRefVector> outputhltmu24muonrefs(new pat::MuonRefVector);
  std::auto_ptr<pat::MuonRefVector> outputhlttkmu24muonrefs(new pat::MuonRefVector);
  std::auto_ptr<pat::MuonRefVector> outputhltmumuonrefs(new pat::MuonRefVector);
  std::auto_ptr<pat::MuonRefVector> outputhlttkmumuonrefs(new pat::MuonRefVector);
  // muon ID
  std::auto_ptr<pat::MuonRefVector> outputloosemuonrefs(new pat::MuonRefVector);
  std::auto_ptr<pat::MuonRefVector> outputtightmuonrefs(new pat::MuonRefVector);
  std::auto_ptr<pat::MuonCollection> outputtightmuons(new pat::MuonCollection);

  std::auto_ptr<edm::ValueMap<float> > outputelnvtxmap(new ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > outputelwgtmap(new ValueMap<float>());
  // single ele trigger info
  std::auto_ptr<pat::ElectronRefVector> outputhltele24eta2p1wplooseelectronrefs (new pat::ElectronRefVector);
  std::auto_ptr<pat::ElectronRefVector> outputhltele25eta2p1wptightelectronrefs (new pat::ElectronRefVector);
  std::auto_ptr<pat::ElectronRefVector> outputhltele27eta2p1wplooseelectronrefs (new pat::ElectronRefVector);
  std::auto_ptr<pat::ElectronRefVector> outputhltele27eta2p1wptightelectronrefs (new pat::ElectronRefVector);
  std::auto_ptr<pat::ElectronRefVector> outputhltele27wptightelectronrefs (new pat::ElectronRefVector);
  std::auto_ptr<pat::ElectronRefVector> outputhltele105electronrefs (new pat::ElectronRefVector);
  std::auto_ptr<pat::ElectronRefVector> outputhltele115electronrefs (new pat::ElectronRefVector);
  std::auto_ptr<pat::ElectronRefVector> outputhlthltelelectronrefs (new pat::ElectronRefVector);
  //electron id
  std::auto_ptr<pat::ElectronRefVector> outputvetoelectronrefs(new pat::ElectronRefVector);
  std::auto_ptr<pat::ElectronRefVector> outputlooseelectronrefs(new pat::ElectronRefVector);
  std::auto_ptr<pat::ElectronRefVector> outputmediumelectronrefs(new pat::ElectronRefVector);
  std::auto_ptr<pat::ElectronRefVector> outputtightelectronrefs(new pat::ElectronRefVector);
  std::auto_ptr<pat::ElectronCollection> outputtightelectrons(new pat::ElectronCollection);
  
  // event weight
  float wgt = 1.0;
  if (geninfoH.isValid()) wgt = geninfoH->weight(); 
  // take all the trigger name for the event
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerResultsH);

  vector<float> munvtxvector;
  vector<float> muwgtvector;
  // loop on the probe muon collection
  for (vector<pat::Muon>::const_iterator muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
    // calculate isolation
    float isoval = muons_iter->pfIsolationR04().sumNeutralHadronEt;
    isoval += muons_iter->pfIsolationR04().sumPhotonEt;
    isoval -= 0.5*muons_iter->pfIsolationR04().sumPUPt;
    if (isoval < 0.) isoval = 0.;
    isoval += muons_iter->pfIsolationR04().sumChargedHadronPt;
    isoval /= muons_iter->pt();

    bool triggermatched = false;
    bool hltisomu20matched = false;
    bool hltisotkmu20matched = false;
    bool hltisomu22matched = false;
    bool hltisotkmu22matched = false;
    bool hltisomu24matched = false;
    bool hltisotkmu24matched = false;
    bool hltisomumatched = false;
    bool hltisotkmumatched = false;
    // loop on the trigger result objects and upack single objects
    for (pat::TriggerObjectStandAlone trgobj : *triggerObjectsH) {
      trgobj.unpackPathNames(trigNames);
      if(not (deltaR(trgobj.eta(), trgobj.phi(), muons_iter->eta(), muons_iter->phi()) < tagmuontrigmatchdR)) continue;
      for (std::string trigpath : tagmuontriggers) { 
	// loop on the list of tag muon triggers and check whether the trigger object belongs to the path and matched the offilen muon
	if (trgobj.hasPathName(trigpath, false, true)) triggermatched = true; 
      }
    
      // match specifically only IsoMu20
      if (trgobj.hasPathName("HLT_IsoMu20_v*"  , false, true) and not hltisomu20matched)
	hltisomu20matched   = true;
      
      // match specifically on IsoTk20
      if (trgobj.hasPathName("HLT_IsoTkMu20_v*", false, true) and not hltisotkmu20matched)
	hltisotkmu20matched = true;
      
      // match specifically only IsoMu22
      if (trgobj.hasPathName("HLT_IsoMu22_v*"  , false, true) and not hltisomu22matched)
	hltisomu22matched   = true;
      
      // match specifically on IsoTk22
      if (trgobj.hasPathName("HLT_IsoTkMu22_v*", false, true) and not hltisotkmu22matched)
	hltisotkmu22matched = true;
      // match specifically only IsoMu24
      if (trgobj.hasPathName("HLT_IsoMu24_v*"  , false, true) and not hltisomu24matched)	
	hltisomu24matched   = true;
      
      // match specifically on IsoTk24
      if (trgobj.hasPathName("HLT_IsoTkMu24_v*", false, true) and not hltisotkmu24matched) 
	hltisotkmu24matched = true;
      
      if(hltisomu20matched || hltisomu22matched || hltisomu24matched) hltisomumatched = true;
      if(hltisotkmu20matched || hltisotkmu22matched || hltisotkmu24matched) hltisotkmumatched = true;      
    }
    
    if (!requiremuonhlt) triggermatched = true;    
    // matched to mu20
    if (verticesH->size() != 0 && hltisomu20matched) 
      outputhltmu20muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to IsoMu20
    if (verticesH->size() != 0 && hltisotkmu20matched) 
      outputhlttkmu20muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to mu22
    if (verticesH->size() != 0 && hltisomu22matched) 
      outputhltmu22muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to IsoMu22
    if (verticesH->size() != 0 && hltisotkmu22matched) 
    outputhlttkmu22muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to mu24
    if (verticesH->size() != 0 && hltisomu24matched) 
      outputhltmu24muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to IsoMu24
    if (verticesH->size() != 0 && hltisotkmu24matched) 
      outputhlttkmu24muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
  // matched to mu20 || mu22 || mu24
    if (verticesH->size() != 0 && hltisomumatched) 
      outputhltmumuonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to IsoMu20 || IsoMu22 || IsoMu24
    if (verticesH->size() != 0 && hltisotkmumatched) 
      outputhlttkmumuonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    

    // Loose muons
    if (verticesH->size() != 0 && muon::isLooseMuon(*muons_iter) && isoval <= loosemuisocut) 
      outputloosemuonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // Tight muons
    if (verticesH->size() != 0 && muon::isTightMuon(*muons_iter, *(verticesH->begin())) && isoval <= tightmuisocut) 
      outputtightmuonrefs    ->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // real tught matched with trigger and passing eta and pt cuts
    if (verticesH->size() != 0 && muon::isTightMuon(*muons_iter, *(verticesH->begin())) && isoval <= tightmuisocut) {
      if (triggermatched && muons_iter->pt() > tagmuonptcut && fabs(muons_iter->eta()) < tagmuonetacut) 
	outputtightmuons->push_back(*muons_iter);
    }
  
    munvtxvector.push_back(float(verticesH->size()));
    muwgtvector.push_back(wgt);
  }

  // electron part  
  vector<float> elnvtxvector;
  vector<float> elwgtvector;
  for (vector<pat::Electron>::const_iterator electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) {
    const Ptr<pat::Electron> electronPtr(electronsH, electrons_iter - electronsH->begin());

    bool triggermatched = false;
    bool hltele24eta2p1wploosematched = false;
    bool hltele25eta2p1wptightmatched = false;
    bool hltele27eta2p1wploosematched = false;
    bool hltele27eta2p1wptightmatched = false;
    bool hltele27wptightmatched = false;
    bool hltele105matched = false;
    bool hltele115matched = false;
    bool hltelematched = false;

    for (pat::TriggerObjectStandAlone trgobj : *triggerObjectsH) {
      trgobj.unpackPathNames(trigNames);
      if(not (deltaR(trgobj.eta(), trgobj.phi(), electrons_iter->eta(), electrons_iter->phi()) < tagelectrontrigmatchdR)) continue;
      for (std::string trigpath : tagelectrontriggers) {
	if (trgobj.hasPathName(trigpath, false, true)) triggermatched = true;
      }

      if (trgobj.hasPathName("HLT_Ele24_eta2p1_WPLoose_Gsf_v*", false, true) and not hltele24eta2p1wploosematched)
        hltele24eta2p1wploosematched = true;

      if (trgobj.hasPathName("HLT_Ele25_eta2p1_WPTight_Gsf_v*"  , false, true) and not hltele25eta2p1wptightmatched)
        hltele25eta2p1wptightmatched = true;

      if (trgobj.hasPathName("HLT_Ele27_eta2p1_WPLoose_Gsf_v*"  , false, true) and not hltele27eta2p1wploosematched)
        hltele27eta2p1wploosematched = true;

      if (trgobj.hasPathName("HLT_Ele27_eta2p1_WPTight_Gsf_v*"  , false, true) and not hltele27eta2p1wptightmatched)
        hltele27eta2p1wptightmatched = true;

      if (trgobj.hasPathName("HLT_Ele27_WPTight_Gsf_v*"  , false, true) and not hltele27wptightmatched)
        hltele27wptightmatched = true;
      
      if (trgobj.hasPathName("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*"  , false, true) and not hltele105matched)
        hltele105matched = true;

      if (trgobj.hasPathName("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"  , false, true) and not hltele115matched)
        hltele115matched = true;
      
      if(hltele24eta2p1wploosematched || hltele25eta2p1wptightmatched || hltele27eta2p1wploosematched || hltele27eta2p1wptightmatched || hltele27wptightmatched || hltele105matched || hltele115matched) hltelematched = true;
    }

     if (!requireelectronhlt) triggermatched = true;    
    if (verticesH->size() != 0 && hltele24eta2p1wploosematched) 
      outputhltele24eta2p1wplooseelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 && hltele25eta2p1wptightmatched) 
      outputhltele25eta2p1wptightelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 && hltele27eta2p1wploosematched) 
      outputhltele27eta2p1wplooseelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 && hltele27eta2p1wptightmatched) 
      outputhltele27eta2p1wptightelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 && hltele27wptightmatched) 
      outputhltele27wptightelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 && hltele105matched) 
      outputhltele105electronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 && hltele115matched) 
      outputhltele115electronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 && hltelematched) 
      outputhlthltelelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin())); 

    // veto electrons
    if (verticesH->size() != 0 && (*electronVetoIdH)  [electronPtr]) 
      outputvetoelectronrefs  ->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    // loose electrons
    if (verticesH->size() != 0 && (*electronLooseIdH) [electronPtr]) 
      outputlooseelectronrefs ->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    // medium electrons
    if (verticesH->size() != 0 && (*electronMediumIdH)[electronPtr]) 
      outputmediumelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    // tight electrons
    if (verticesH->size() != 0 && (*electronTightIdH) [electronPtr]) 
      outputtightelectronrefs ->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    //
    if (verticesH->size() != 0 && (*electronTightIdH) [electronPtr]) {
      if (triggermatched && electrons_iter->pt() > tagelectronptcut && fabs(electrons_iter->eta()) < tagelectronetacut) 
	outputtightelectrons->push_back(*electrons_iter);
    }
    elnvtxvector.push_back(float(verticesH->size()));
    elwgtvector.push_back(wgt);
  }
  
  edm::ValueMap<float>::Filler munvtxfiller(*outputmunvtxmap);
  munvtxfiller.insert(muonsH, munvtxvector.begin(), munvtxvector.end());
  munvtxfiller.fill();
  
  edm::ValueMap<float>::Filler muwgtfiller(*outputmuwgtmap);
  muwgtfiller.insert(muonsH, muwgtvector.begin(), muwgtvector.end());
  muwgtfiller.fill();
  
  edm::ValueMap<float>::Filler elnvtxfiller(*outputelnvtxmap);
  elnvtxfiller.insert(electronsH, elnvtxvector.begin(), elnvtxvector.end());
  elnvtxfiller.fill();
  
  edm::ValueMap<float>::Filler elwgtfiller(*outputelwgtmap);
  elwgtfiller.insert(electronsH, elwgtvector.begin(), elwgtvector.end());
  elwgtfiller.fill();
  
  iEvent.put(outputmunvtxmap,"munvtxmap");
  iEvent.put(outputmuwgtmap, "muwgtmap");

  iEvent.put(outputhltmu20muonrefs,  "hltmu20muonrefs");
  iEvent.put(outputhlttkmu20muonrefs,"hlttkmu20muonrefs");
  iEvent.put(outputhltmu22muonrefs,  "hltmu22muonrefs");
  iEvent.put(outputhlttkmu22muonrefs,"hlttkmu22muonrefs");
  iEvent.put(outputhltmu24muonrefs,  "hltmu24muonrefs");
  iEvent.put(outputhlttkmu24muonrefs,"hlttkmu24muonrefs");
  iEvent.put(outputhltmumuonrefs,    "hltmumuonrefs");
  iEvent.put(outputhlttkmumuonrefs,  "hlttkmumuonrefs");

  iEvent.put(outputloosemuonrefs, "loosemuonrefs");
  iEvent.put(outputtightmuonrefs, "tightmuonrefs");
  iEvent.put(outputtightmuons,    "tightmuons");

  iEvent.put(outputelnvtxmap, "elnvtxmap");
  iEvent.put(outputelwgtmap,  "elwgtmap");

  // single ele trigger info                                                                                                                                                   
  iEvent.put(outputhltele24eta2p1wplooseelectronrefs,"hltele24eta2p1wplooseelectronrefs");
  iEvent.put(outputhltele25eta2p1wptightelectronrefs,"hltele25eta2p1wptightelectronrefs");
  iEvent.put(outputhltele27eta2p1wplooseelectronrefs,"hltele27eta2p1wplooseelectronrefs");
  iEvent.put(outputhltele27eta2p1wptightelectronrefs,"hltele27eta2p1wptightelectronrefs");
  iEvent.put(outputhltele27wptightelectronrefs,"hltele27wptightelectronrefs");
  iEvent.put(outputhltele105electronrefs,      "hltele105electronrefs");
  iEvent.put(outputhltele115electronrefs,      "hltele115electronrefs");
  iEvent.put(outputhlthltelelectronrefs,       "hltelelectronrefs");
  //electron id                                                                                                                                                                
  iEvent.put(outputvetoelectronrefs,  "vetoelectronrefs");
  iEvent.put(outputlooseelectronrefs, "looseelectronrefs");
  iEvent.put(outputmediumelectronrefs,"mediumelectronrefs");
  iEvent.put(outputtightelectronrefs, "tightelectronrefs");
  iEvent.put(outputtightelectrons,    "tightelectrons");
}

void LeptonTnPInfoProducer::beginJob() {
}

void LeptonTnPInfoProducer::endJob() {
}

void LeptonTnPInfoProducer::beginRun(edm::Run const&, edm::EventSetup const&) {
}

void LeptonTnPInfoProducer::endRun(edm::Run const&, edm::EventSetup const&) {
}

void LeptonTnPInfoProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void LeptonTnPInfoProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void LeptonTnPInfoProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(LeptonTnPInfoProducer);
