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
#include "DataFormats/PatCandidates/interface/Photon.h"
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

  // Gen info
  const edm::EDGetTokenT<GenEventInfoProduct>        geninfoToken;
  // Vertices
  const edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;
  //Leptons
  const edm::EDGetTokenT<pat::MuonCollection>     muonsToken;
  const edm::EDGetTokenT<pat::ElectronCollection> electronsToken;
  const edm::EDGetTokenT<pat::ElectronCollection> electronsFullCollectionToken;
  const edm::EDGetTokenT<pat::PhotonCollection>   photonsToken;

  // Triggers
  const edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken;
  const edm::EDGetTokenT<edm::TriggerResults>  triggerResultsToken;

  // ID value maps
  const edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > electronLooseIdMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > electronTightIdMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > electronHLTSafeIdMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > electronMvaLooseMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > electronMvaTightMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > photonLooseIdMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > photonMediumIdMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > photonTightIdMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > photonMvaLooseMapToken;
  const edm::EDGetTokenT<edm::ValueMap<bool> > photonMvaTightMapToken;

  // Additional selections for muons
  const edm::ParameterSet tagloosemuons;
  const edm::ParameterSet tagtightmuons;
  const edm::ParameterSet tagmuontriggermatch;
  // Additional electron selections
  const edm::ParameterSet tagelectrons;
  const edm::ParameterSet tagelectrontriggermatch;
};

LeptonTnPInfoProducer::LeptonTnPInfoProducer(const edm::ParameterSet& iConfig):
    geninfoToken   (consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("geninfo"))),
    verticesToken  (consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("vertices"))),
    // leptons
    muonsToken     (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronsToken (consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    electronsFullCollectionToken(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronsFullCollection"))),
    photonsToken   (consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    // triggers
    triggerObjectsToken (consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerobjects"))),
    triggerResultsToken (consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
    // id and isolation maps
    electronVetoIdMapToken    (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronvetoid"))),
    electronLooseIdMapToken   (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronlooseid"))),
    electronMediumIdMapToken  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronmediumid"))),
    electronTightIdMapToken   (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electrontightid"))),
    electronHLTSafeIdMapToken (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronhltsafeid"))),
    electronMvaLooseMapToken (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronmvalooseid"))),
    electronMvaTightMapToken (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronmvatightid"))),
    photonLooseIdMapToken  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonlooseid"))),
    photonMediumIdMapToken (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonmediumid"))),
    photonTightIdMapToken  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photontightid"))),
    photonMvaLooseMapToken  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonmvalooseid"))),
    photonMvaTightMapToken  (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonmvatightid"))),
    // addition
    tagloosemuons (iConfig.getParameter<edm::ParameterSet>("tagloosemuons")),
    tagtightmuons (iConfig.getParameter<edm::ParameterSet>("tagtightmuons")),
    tagmuontriggermatch (iConfig.getParameter<edm::ParameterSet>("tagmuontriggermatch")),
    tagelectrons  (iConfig.getParameter<edm::ParameterSet>("tagelectrons")),
    tagelectrontriggermatch  (iConfig.getParameter<edm::ParameterSet>("tagelectrontriggermatch"))

{
  // produce a map with the generator weight
  produces<edm::ValueMap<float> >("muwgtmap");
  // produce a map with nvtx for a given muon
  produces<edm::ValueMap<float> >("munvtxmap");
  produces<edm::ValueMap<float> >("mudxymap");
  produces<edm::ValueMap<float> >("mudzmap");
  produces<edm::ValueMap<float> >("muchi2map");
  produces<edm::ValueMap<float> >("munvalidhitmap");
  produces<edm::ValueMap<float> >("munpixelhitmap");
  produces<edm::ValueMap<float> >("muntrackerlayermap");


  // produce a ref vector with tag muon trigger info
  produces<pat::MuonRefVector>("hltmu20muonrefs");
  produces<pat::MuonRefVector>("hlttkmu20muonrefs");
  produces<pat::MuonRefVector>("hltmu22muonrefs");
  produces<pat::MuonRefVector>("hlttkmu22muonrefs");
  produces<pat::MuonRefVector>("hltmu24muonrefs");
  produces<pat::MuonRefVector>("hlttkmu24muonrefs");
  produces<pat::MuonRefVector>("hltmu50muonrefs");
  produces<pat::MuonRefVector>("hlttkmu50muonrefs");
  produces<pat::MuonRefVector>("hltmumuonrefs");
  produces<pat::MuonRefVector>("hlttkmumuonrefs");

  produces<pat::MuonRefVector>("hltmu17mu8Leg17");
  produces<pat::MuonRefVector>("hltmu17mu8Leg8");
  produces<pat::MuonRefVector>("hltmu17tkmu8Leg17");
  produces<pat::MuonRefVector>("hltmu17tkmu8Leg8");
  produces<pat::MuonRefVector>("hltmu17mu8dzLeg17");
  produces<pat::MuonRefVector>("hltmu17mu8dzLeg8");
  produces<pat::MuonRefVector>("hltmu17tkmu8dzLeg17");
  produces<pat::MuonRefVector>("hltmu17tkmu8dzLeg8");

  // produces collection for loose and tight muons
  produces<pat::MuonRefVector>("loosemuonrefs");
  produces<pat::MuonRefVector>("tightmuonrefs");
  produces<pat::MuonCollection>("tightmuons");

  /// same logic also for electrons
  produces<edm::ValueMap<float> >("elnvtxmap");
  produces<edm::ValueMap<float> >("elwgtmap");
  produces<edm::ValueMap<float> >("eldxymap");
  produces<edm::ValueMap<float> >("eldzmap");
  produces<edm::ValueMap<float> >("phnvtxmap");
  produces<edm::ValueMap<float> >("phwgtmap");


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
  produces<pat::ElectronRefVector>("hltsafeelectronrefs");
  produces<pat::ElectronRefVector>("mvalooseelectronrefs");
  produces<pat::ElectronRefVector>("mvatightelectronrefs");
  produces<pat::ElectronCollection>("tightelectrons");

  produces<pat::PhotonRefVector>("loosephotonrefs");
  produces<pat::PhotonRefVector>("mediumphotonrefs");
  produces<pat::PhotonRefVector>("tightphotonrefs");
  produces<pat::PhotonRefVector>("mvaloosephotonrefs");
  produces<pat::PhotonRefVector>("mvatightphotonrefs");
  produces<pat::PhotonRefVector>("recoelectronmatch");

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

  Handle<pat::MuonCollection> muonsH;
  iEvent.getByToken(muonsToken, muonsH);

  Handle<pat::ElectronCollection> electronsH;
  iEvent.getByToken(electronsToken, electronsH);

  Handle<pat::ElectronCollection> electronsFullCollectionH;
  iEvent.getByToken(electronsFullCollectionToken, electronsFullCollectionH);

  Handle<pat::PhotonCollection> photonsH;
  iEvent.getByToken(photonsToken, photonsH);

  Handle<pat::TriggerObjectStandAloneCollection> triggerObjectsH;
  iEvent.getByToken(triggerObjectsToken, triggerObjectsH);
  const pat::TriggerObjectStandAloneCollection triggerObjects = *triggerObjectsH;

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
  Handle<edm::ValueMap<bool> > electronHLTSafeIdH;
  iEvent.getByToken(electronHLTSafeIdMapToken, electronHLTSafeIdH);
  Handle<edm::ValueMap<bool> > electronMvaLooseIdH;
  iEvent.getByToken(electronMvaLooseMapToken, electronMvaLooseIdH);
  Handle<edm::ValueMap<bool> > electronMvaTightIdH;
  iEvent.getByToken(electronMvaTightMapToken, electronMvaTightIdH);
  Handle<edm::ValueMap<bool> > photonLooseIdH;
  iEvent.getByToken(photonLooseIdMapToken, photonLooseIdH);
  Handle<edm::ValueMap<bool> > photonMediumIdH;
  iEvent.getByToken(photonMediumIdMapToken, photonMediumIdH);
  Handle<edm::ValueMap<bool> > photonTightIdH;
  iEvent.getByToken(photonTightIdMapToken, photonTightIdH);
  Handle<edm::ValueMap<bool> > photonMvaLooseIdH;
  iEvent.getByToken(photonMvaLooseMapToken, photonMvaLooseIdH);
  Handle<edm::ValueMap<bool> > photonMvaTightIdH;
  iEvent.getByToken(photonMvaTightMapToken, photonMvaTightIdH);

  // output collection
  std::unique_ptr<edm::ValueMap<float> > outputmunvtxmap(new ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> > outputmuwgtmap(new ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> > outputmudxymap(new ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> > outputmudzmap(new ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> > outputmuchi2map(new ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> > outputmunvalidhitmap(new ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> > outputmunpixelhitmap(new ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> > outputmuntrackerlayermap(new ValueMap<float>());
  // trigger muons
  std::unique_ptr<pat::MuonRefVector> outputhltmu20muonrefs(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhlttkmu20muonrefs(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhltmu22muonrefs(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhlttkmu22muonrefs(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhltmu24muonrefs(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhlttkmu24muonrefs(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhltmu50muonrefs(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhlttkmu50muonrefs(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhltmumuonrefs(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhlttkmumuonrefs(new pat::MuonRefVector);

  std::unique_ptr<pat::MuonRefVector> outputhltmu17mu8muonrefs_leg17(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhltmu17mu8muonrefs_leg8(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhltmu17tkmu8muonrefs_leg17(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhltmu17tkmu8muonrefs_leg8(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhltmu17mu8dzmuonrefs_leg17(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhltmu17mu8dzmuonrefs_leg8(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhltmu17tkmu8dzmuonrefs_leg17(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputhltmu17tkmu8dzmuonrefs_leg8(new pat::MuonRefVector);

  // muon ID
  std::unique_ptr<pat::MuonRefVector> outputloosemuonrefs(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonRefVector> outputtightmuonrefs(new pat::MuonRefVector);
  std::unique_ptr<pat::MuonCollection> outputtightmuons(new pat::MuonCollection);

  std::unique_ptr<edm::ValueMap<float> > outputelnvtxmap(new ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> > outputelwgtmap(new ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> > outputeldxymap(new ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> > outputeldzmap(new ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> > outputphnvtxmap(new ValueMap<float>());
  std::unique_ptr<edm::ValueMap<float> > outputphwgtmap(new ValueMap<float>());
  // single ele trigger info
  std::unique_ptr<pat::ElectronRefVector> outputhltele24eta2p1wplooseelectronrefs (new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputhltele25eta2p1wptightelectronrefs (new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputhltele27eta2p1wplooseelectronrefs (new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputhltele27eta2p1wptightelectronrefs (new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputhltele27wptightelectronrefs (new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputhltele105electronrefs (new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputhltele115electronrefs (new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputhlthltelelectronrefs (new pat::ElectronRefVector);
  //electron id
  std::unique_ptr<pat::ElectronRefVector> outputvetoelectronrefs(new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputlooseelectronrefs(new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputmediumelectronrefs(new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputtightelectronrefs(new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputhltsafeelectronrefs(new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputmvalooseelectronrefs(new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronRefVector> outputmvatightelectronrefs(new pat::ElectronRefVector);
  std::unique_ptr<pat::ElectronCollection> outputtightelectrons(new pat::ElectronCollection);
  //photon id
  std::unique_ptr<pat::PhotonRefVector> outputloosephotonrefs(new pat::PhotonRefVector);
  std::unique_ptr<pat::PhotonRefVector> outputmediumphotonrefs(new pat::PhotonRefVector);
  std::unique_ptr<pat::PhotonRefVector> outputtightphotonrefs(new pat::PhotonRefVector);
  std::unique_ptr<pat::PhotonRefVector> outputmvaloosephotonrefs(new pat::PhotonRefVector);
  std::unique_ptr<pat::PhotonRefVector> outputmvatightphotonrefs(new pat::PhotonRefVector);
  std::unique_ptr<pat::PhotonCollection> outputtightphotons(new pat::PhotonCollection);
  std::unique_ptr<pat::PhotonRefVector> outputrecoelectronmatchrefs (new pat::PhotonRefVector);

  // event weight
  float wgt = 1.0;
  if (geninfoH.isValid()) wgt = geninfoH->weight();
  // take all the trigger name for the event
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*triggerResultsH);

  vector<float> munvtxvector;
  vector<float> muwgtvector;
  vector<float> mudxyvector;
  vector<float> mudzvector;
  vector<float> muchi2vector;
  vector<float> munvalidhitvector;
  vector<float> munpixelhitvector;
  vector<float> muntrackerlayervector;

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
    bool hltmu50matched = false;
    bool hlttkmu50matched = false;
    bool hltmumatched = false;
    bool hlttkmumatched = false;

    bool hltmu17mu8matched_leg17 = false;
    bool hltmu17mu8matched_leg8 = false;
    bool hltmu17tkmu8matched_leg17 = false;
    bool hltmu17tkmu8matched_leg8   = false;
    bool hltmu17mu8dzmatched_leg17 = false;
    bool hltmu17mu8dzmatched_leg8  = false;
    bool hltmu17tkmu8dzmatched_leg17 = false;
    bool hltmu17tkmu8dzmatched_leg8  = false;

    // loop on the whole trigger object collection
    for (pat::TriggerObjectStandAlone trgobj : *triggerObjectsH) {
      trgobj.unpackPathNames(trigNames); // un-pack names
      if(not (deltaR(trgobj.eta(), trgobj.phi(), muons_iter->eta(), muons_iter->phi()) < tagmuontriggermatch.getParameter<double>("tagmuontrigmatchdR"))) continue; //check dR matching

      for (std::string trigpath : tagmuontriggermatch.getParameter<std::vector<std::string> >("tagmuontriggers")) {
	// loop on the list of tag muon triggers and check whether the trigger object belongs to the path and matched the offilen muon
	if (trgobj.hasPathName(trigpath, true, false) or trgobj.hasPathName(trigpath, true, true) ) triggermatched = true;
      }

      // loop on the list of tag muon triggers and check whether the trigger object belongs to the path and matched the offilen muon
      if (trgobj.hasPathName("HLT_IsoMu20_v*" , true, false) or trgobj.hasPathName("HLT_IsoMu20_v*", true, true)) hltisomu20matched = true;
      if (trgobj.hasPathName("HLT_IsoMu22_v*" , true, false) or trgobj.hasPathName("HLT_IsoMu22_v*", true, true)) hltisomu22matched = true;
      if (trgobj.hasPathName("HLT_IsoMu24_v*" , true, false) or trgobj.hasPathName("HLT_IsoMu24_v*", true, true)) hltisomu24matched = true;
      if (trgobj.hasPathName("HLT_Mu50_v*" , true, false) or trgobj.hasPathName("HLT_Mu50_v*", true, true)) hltmu50matched = true;

      if (trgobj.hasPathName("HLT_IsoTkMu20_v*" , true, false) or trgobj.hasPathName("HLT_IsoTkMu20_v*", true, true)) hltisotkmu20matched = true;
      if (trgobj.hasPathName("HLT_IsoTkMu22_v*" , true, false) or trgobj.hasPathName("HLT_IsoTkMu22_v*", true, true)) hltisotkmu22matched = true;
      if (trgobj.hasPathName("HLT_IsoTkMu24_v*" , true, false) or trgobj.hasPathName("HLT_IsoTkMu24_v*", true, true)) hltisotkmu24matched = true;
      if (trgobj.hasPathName("HLT_TkMu50_v*" , true, false) or trgobj.hasPathName("HLT_TkMu50_v*", true, true)) hlttkmu50matched = true;

      // Double Muon leg 17
      if (trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*" , true, false) or trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", true, true)){
	if(trgobj.hasFilterLabel("hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17"))
	  hltmu17mu8dzmatched_leg17 = true;
      }
      if (trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*" , true, false) or trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", true, true)){
	if(trgobj.hasFilterLabel("hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17"))
	  hltmu17mu8matched_leg17 = true;
      }
      if (trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*" , true, false) or trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*", true, true)){
	if(trgobj.hasFilterLabel("hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17") or trgobj.hasFilterLabel("hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17"))
	  hltmu17tkmu8dzmatched_leg17 = true;
      }
      if (trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*" , true, false) or trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*", true, true)){
	if(trgobj.hasFilterLabel("hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17") or trgobj.hasFilterLabel("hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17"))
	  hltmu17tkmu8matched_leg17 = true;
      }
      // Double Muon leg 8
      if (trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*" , true, false) or trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*", true, true)){
	if(trgobj.hasFilterLabel("hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8") or trgobj.hasFilterLabel("hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8"))
	  hltmu17mu8dzmatched_leg8 = true;
      }
      if (trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*" , true, false) or trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*", true, true)){
	if(trgobj.hasFilterLabel("hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8") or trgobj.hasFilterLabel("hltL3pfL1sDoubleMu114ORDoubleMu125L1f0L2pf0L3PreFiltered8"))
	  hltmu17mu8matched_leg8 = true;
      }
      if (trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*" , true, false) or trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*", true, true)){
	if(trgobj.hasFilterLabel("hltDiMuonGlbFiltered17TrkFiltered8"))
	  hltmu17tkmu8dzmatched_leg8 = true;
      }
      if (trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*" , true, false) or trgobj.hasPathName("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*", true, true)){
	if(trgobj.hasFilterLabel("hltDiMuonGlbFiltered17TrkFiltered8"))
	  hltmu17tkmu8matched_leg8 = true;
      }
    }

    if(hltisomu24matched or hltmu50matched) hltmumatched = true;
    if(hltisotkmu24matched or hlttkmu50matched) hlttkmumatched = true;
    if(not tagmuontriggermatch.getParameter<bool>("requiremuonhlt")) triggermatched = true;

    if(not triggermatched and (hltisomu20matched || hltisomu22matched || hltisomu24matched || hltisotkmu20matched || hltisotkmu22matched || hltisotkmu24matched || hltmu50matched || hlttkmu50matched))
      std::cout<<"Problem with the trigger matching for muons --> triggermathc should be always >= than the big or "<<std::endl;

    // matched to mu20
    if (verticesH->size() != 0 and hltisomu20matched)
      outputhltmu20muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to IsoMu20
    if (verticesH->size() != 0 and hltisotkmu20matched)
      outputhlttkmu20muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to mu22
    if (verticesH->size() != 0 and hltisomu22matched)
      outputhltmu22muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to IsoMu22
    if (verticesH->size() != 0 and hltisotkmu22matched)
      outputhlttkmu22muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to mu24
    if (verticesH->size() != 0 and hltisomu24matched)
      outputhltmu24muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to IsoMu24
    if (verticesH->size() != 0 and hltisotkmu24matched)
      outputhlttkmu24muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to mu50
    if (verticesH->size() != 0 and hltmu50matched)
      outputhltmu50muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to Mu50
    if (verticesH->size() != 0 and hlttkmu50matched)
      outputhlttkmu50muonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to mu20 || mu22 || mu24
    if (verticesH->size() != 0 and hltmumatched)
      outputhltmumuonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // matched to IsoMu20 || IsoMu22 || IsoMu24
    if (verticesH->size() != 0 and hlttkmumatched)
      outputhlttkmumuonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));

    if (verticesH->size() != 0 and hltmu17mu8matched_leg17)
      outputhltmu17mu8muonrefs_leg17->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    if (verticesH->size() != 0 and hltmu17mu8matched_leg8)
      outputhltmu17mu8muonrefs_leg8->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    if (verticesH->size() != 0 and hltmu17tkmu8matched_leg17)
      outputhltmu17tkmu8muonrefs_leg17->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    if (verticesH->size() != 0 and hltmu17tkmu8matched_leg8)
      outputhltmu17tkmu8muonrefs_leg8->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));

    if (verticesH->size() != 0 and hltmu17mu8dzmatched_leg17)
      outputhltmu17mu8dzmuonrefs_leg17->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    if (verticesH->size() != 0 and hltmu17mu8dzmatched_leg8)
      outputhltmu17mu8dzmuonrefs_leg8->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    if (verticesH->size() != 0 and hltmu17tkmu8dzmatched_leg17)
      outputhltmu17tkmu8dzmuonrefs_leg17->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    if (verticesH->size() != 0 and hltmu17tkmu8dzmatched_leg8)
      outputhltmu17tkmu8dzmuonrefs_leg8->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));

    // Loose muons
    if (verticesH->size() != 0 and
	muon::isLooseMuon(*muons_iter) and
	isoval <= tagloosemuons.getParameter<double>("isocut"))
      outputloosemuonrefs->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
    // Tight muons
    if (verticesH->size() != 0 and
	muon::isTightMuon(*muons_iter, *(verticesH->begin())) and
	isoval <= tagtightmuons.getParameter<double>("isocut"))
      outputtightmuonrefs    ->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));

    // real tight matched with trigger and passing eta and pt cuts
    if (verticesH->size() != 0 and
	muon::isTightMuon(*muons_iter, *(verticesH->begin())) and
	isoval <= tagtightmuons.getParameter<double>("isocut")) {
      if (triggermatched and
	  muons_iter->pt() > tagtightmuons.getParameter<double>("ptcut") and
	  fabs(muons_iter->eta()) < tagtightmuons.getParameter<double>("etacut"))
	outputtightmuons->push_back(*muons_iter);
    }
    ///
    munvtxvector.push_back(float(verticesH->size()));
    muwgtvector.push_back(wgt);

    if(muons_iter->isPFMuon() or muons_iter->isGlobalMuon()){
      if(muons_iter->isGlobalMuon() and muons_iter->globalTrack().isNonnull() and muons_iter->globalTrack().isAvailable()){
	muchi2vector.push_back(muons_iter->globalTrack()->normalizedChi2());
	munvalidhitvector.push_back(muons_iter->globalTrack()->hitPattern().numberOfValidMuonHits());
      }else{
	muchi2vector.push_back(-99.);
	munvalidhitvector.push_back(-99.);
      }
      if(muons_iter->innerTrack().isNonnull() and muons_iter->innerTrack().isAvailable()){
	munpixelhitvector.push_back(muons_iter->innerTrack()->hitPattern().numberOfValidPixelHits());
	muntrackerlayervector.push_back(muons_iter->innerTrack()->hitPattern().trackerLayersWithMeasurement());
      }
      else{
	munpixelhitvector.push_back(-99.);
	muntrackerlayervector.push_back(-99.);
      }
      if(muons_iter->muonBestTrack().isNonnull() and muons_iter->muonBestTrack().isAvailable() and verticesH->begin() != verticesH->end()){
	mudxyvector.push_back(muons_iter->muonBestTrack()->dxy(verticesH->begin()->position()));
	mudzvector.push_back(muons_iter->muonBestTrack()->dz(verticesH->begin()->position()));
      }
      else{
	mudxyvector.push_back(-99.);
	mudzvector.push_back(-99.);
      }
    }
    else{
      muchi2vector.push_back(-99.);
      munvalidhitvector.push_back(-99.);
      munpixelhitvector.push_back(-99.);
      muntrackerlayervector.push_back(-99.);
      mudxyvector.push_back(-99.);
      mudzvector.push_back(-99.);
    }
  }
  // electron part
  vector<float> elnvtxvector;
  vector<float> elwgtvector;
  vector<float> eldxyvector;
  vector<float> eldzvector;

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

    // loop on the trigger result objects and upack single objects
    for (pat::TriggerObjectStandAlone trgobj : *triggerObjectsH) {
      trgobj.unpackPathNames(trigNames);
      if(not (deltaR(trgobj.eta(), trgobj.phi(), electrons_iter->eta(), electrons_iter->phi()) < tagelectrontriggermatch.getParameter<double>("tagelectrontrigmatchdR"))) continue; //check dR matching

      for (std::string trigpath : tagelectrontriggermatch.getParameter<std::vector<std::string> >("tagelectrontriggers")) {
	if (trgobj.hasPathName(trigpath, true, false) or trgobj.hasPathName(trigpath, true, true)) triggermatched = true;
      }

      if (trgobj.hasPathName("HLT_Ele24_eta2p1_WPLoose_Gsf_v*", true, false) or trgobj.hasPathName("HLT_Ele24_eta2p1_WPLoose_Gsf_v*", true, true))
	hltele24eta2p1wploosematched = true;

      if (trgobj.hasPathName("HLT_Ele25_eta2p1_WPTight_Gsf_v*", true, false) or trgobj.hasPathName("HLT_Ele25_eta2p1_WPTight_Gsf_v*", true, true))
	hltele25eta2p1wptightmatched = true;

      if (trgobj.hasPathName("HLT_Ele27_eta2p1_WPLoose_Gsf_v*", true, false) or trgobj.hasPathName("HLT_Ele27_eta2p1_WPLoose_Gsf_v*", true, true))
	hltele27eta2p1wploosematched = true;

      if (trgobj.hasPathName("HLT_Ele27_eta2p1_WPTight_Gsf_v*", true, false) or trgobj.hasPathName("HLT_Ele27_eta2p1_WPTight_Gsf_v*", true, true))
	hltele27eta2p1wptightmatched = true;

      if (trgobj.hasPathName("HLT_Ele27_WPTight_Gsf_v*", true, false) or trgobj.hasPathName("HLT_Ele27_WPTight_Gsf_v*", true, true))
	hltele27wptightmatched = true;

      if (trgobj.hasPathName("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*", true, false) or trgobj.hasPathName("HLT_Ele105_CaloIdVT_GsfTrkIdT_v*", true, true))
	hltele105matched = true;

      if (trgobj.hasPathName("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*", true, false) or trgobj.hasPathName("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*", true, true))
	hltele115matched = true;
    }

    if(triggermatched or
       hltele24eta2p1wploosematched or
       hltele25eta2p1wptightmatched or
       hltele27eta2p1wploosematched or
       hltele27eta2p1wptightmatched or
       hltele27wptightmatched or
       hltele105matched or
       hltele115matched) hltelematched = true;

    if (not tagelectrontriggermatch.getParameter<bool>("requireelectronhlt")) triggermatched = true;


    if(not triggermatched and
       (hltele24eta2p1wploosematched or hltele25eta2p1wptightmatched or
	hltele27eta2p1wploosematched or hltele27eta2p1wptightmatched or
	hltele27wptightmatched or hltele105matched or hltele115matched))
      std::cout<<"Problem with the trigger matching for electrons --> triggermathcing should be always >= than the big or "<<std::endl;

    if(fabs(electrons_iter->eta()) < 1 and hltele27wptightmatched and not hltele27eta2p1wploosematched)
      std::cout<<"Problem with electorns "<<endl;

    if (verticesH->size() != 0 and hltele24eta2p1wploosematched)
      outputhltele24eta2p1wplooseelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 and hltele25eta2p1wptightmatched)
      outputhltele25eta2p1wptightelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 and hltele27eta2p1wploosematched)
      outputhltele27eta2p1wplooseelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 and hltele27eta2p1wptightmatched)
      outputhltele27eta2p1wptightelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 and hltele27wptightmatched)
      outputhltele27wptightelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 and hltele105matched)
      outputhltele105electronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 and hltele115matched)
      outputhltele115electronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    if (verticesH->size() != 0 and hltelematched)
      outputhlthltelelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));

    //
    bool passPVconstraints = true;
    if(tagelectrons.getParameter<bool>("applyPVSelection") and verticesH->size() != 0){
      if(fabs(electrons_iter->superCluster()->eta()) < 1.479 and fabs(electrons_iter->gsfTrack()->dxy(verticesH->at(0).position())) > tagelectrons.getParameter<double>("d0Barrel"))
	passPVconstraints = false;
      else if(fabs(electrons_iter->superCluster()->eta()) > 1.479 and fabs(electrons_iter->gsfTrack()->dxy(verticesH->at(0).position())) > tagelectrons.getParameter<double>("d0Endcap"))
	passPVconstraints = false;
      if(fabs(electrons_iter->superCluster()->eta()) < 1.479 and fabs(electrons_iter->gsfTrack()->dz(verticesH->at(0).position())) > tagelectrons.getParameter<double>("dzBarrel"))
	passPVconstraints = false;
      else if(fabs(electrons_iter->superCluster()->eta()) > 1.479 and fabs(electrons_iter->gsfTrack()->dz(verticesH->at(0).position())) > tagelectrons.getParameter<double>("dzEndcap"))
	passPVconstraints = false;
    }

    // veto electrons
    if (verticesH->size() != 0 and (*electronVetoIdH)  [electronPtr] and passPVconstraints) {
      outputvetoelectronrefs  ->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    }
    // loose electrons
    if (verticesH->size() != 0 and (*electronLooseIdH) [electronPtr] and passPVconstraints)
      outputlooseelectronrefs ->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    // medium electrons
    if (verticesH->size() != 0 and (*electronMediumIdH)[electronPtr] and passPVconstraints)
      outputmediumelectronrefs->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    // tight electrons
    if (verticesH->size() != 0 and (*electronTightIdH) [electronPtr] and passPVconstraints) {
      outputtightelectronrefs ->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    }
    // mva-loose electrons
    if (verticesH->size() != 0 and (*electronMvaLooseIdH) [electronPtr] and passPVconstraints) {
      outputmvalooseelectronrefs ->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    }
    // mva-tight electrons
    if (verticesH->size() != 0 and (*electronMvaTightIdH) [electronPtr] and passPVconstraints) {
      outputmvatightelectronrefs ->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    }
    // HLT safe electrons
    if (verticesH->size() != 0 and (*electronHLTSafeIdH) [electronPtr] and passPVconstraints) {
      outputhltsafeelectronrefs ->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    }
    //
    if (verticesH->size() != 0 and (*electronTightIdH) [electronPtr] and passPVconstraints) {
      if (triggermatched and
	  electrons_iter->pt() > tagelectrons.getParameter<double>("tagelectronptcut") and
	  fabs(electrons_iter->superCluster()->eta()) < tagelectrons.getParameter<double>("tagelectronetacut")) {
	outputtightelectrons->push_back(*electrons_iter);
      }
    }

    elnvtxvector.push_back(float(verticesH->size()));
    elwgtvector.push_back(wgt);

    if(electrons_iter->gsfTrack().isNonnull() and electrons_iter->gsfTrack().isAvailable() and verticesH->begin() != verticesH->end()){
      eldxyvector.push_back(electrons_iter->gsfTrack()->dxy(verticesH->begin()->position()));
      eldzvector.push_back(electrons_iter->gsfTrack()->dz(verticesH->begin()->position()));
    }

    else{
      eldxyvector.push_back(-99.);
      eldzvector.push_back(-99.);
    }

  }

  // photon part
  vector<float> phnvtxvector;
  vector<float> phwgtvector;
  for (vector<pat::Photon>::const_iterator photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
    const Ptr<pat::Photon> photonPtr(photonsH, photons_iter - photonsH->begin());

    phnvtxvector.push_back(float(verticesH->size()));
    phwgtvector.push_back(wgt);

    bool recoelectronmatch =  false;

    // match probe with full collection
    float minDR = 100.;
    for(vector<pat::Electron>::const_iterator itEle = electronsFullCollectionH->begin(); itEle != electronsFullCollectionH->end(); ++itEle) {
      float DR = deltaR(photons_iter->eta(), photons_iter->phi(),itEle->eta(),itEle->phi());
      if(DR < tagelectrontriggermatch.getParameter<double>("tagelectrontrigmatchdR") and DR < minDR){
	DR = deltaR(photons_iter->eta(), photons_iter->phi(),itEle->eta(),itEle->phi());
	recoelectronmatch = true;
      }
    }
    // loose photons --> NOT ADDING ELECTRON VETO otherwise Zee events cannot be used for photon id
    if (verticesH->size() != 0 and (*photonLooseIdH) [photonPtr])
      outputloosephotonrefs ->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
    // medium photons
    if (verticesH->size() != 0 and (*photonMediumIdH)[photonPtr])
      outputmediumphotonrefs->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
    // tight photons
    if (verticesH->size() != 0 and (*photonTightIdH) [photonPtr])
      outputtightphotonrefs ->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
    // tight photons
    if (verticesH->size() != 0 and (*photonMvaLooseIdH) [photonPtr])
      outputmvaloosephotonrefs ->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
    // tight photons
    if (verticesH->size() != 0 and (*photonMvaTightIdH) [photonPtr])
      outputmvatightphotonrefs ->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
    // mathcing with reco electrons
    if (verticesH->size() != 0 and recoelectronmatch)
      outputrecoelectronmatchrefs->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));

  }

  edm::ValueMap<float>::Filler munvtxfiller(*outputmunvtxmap);
  munvtxfiller.insert(muonsH, munvtxvector.begin(), munvtxvector.end());
  munvtxfiller.fill();

  edm::ValueMap<float>::Filler muwgtfiller(*outputmuwgtmap);
  muwgtfiller.insert(muonsH, muwgtvector.begin(), muwgtvector.end());
  muwgtfiller.fill();

  edm::ValueMap<float>::Filler mudxyfiller(*outputmudxymap);
  mudxyfiller.insert(muonsH, mudxyvector.begin(), mudxyvector.end());
  mudxyfiller.fill();

  edm::ValueMap<float>::Filler mudzfiller(*outputmudzmap);
  mudzfiller.insert(muonsH, mudzvector.begin(), mudzvector.end());
  mudzfiller.fill();

  edm::ValueMap<float>::Filler muchi2filler(*outputmuchi2map);
  muchi2filler.insert(muonsH, muchi2vector.begin(), muchi2vector.end());
  muchi2filler.fill();

  edm::ValueMap<float>::Filler munvalidhitfiller(*outputmunvalidhitmap);
  munvalidhitfiller.insert(muonsH, munvalidhitvector.begin(), munvalidhitvector.end());
  munvalidhitfiller.fill();

  edm::ValueMap<float>::Filler munpixelhitfiller(*outputmunpixelhitmap);
  munpixelhitfiller.insert(muonsH, munpixelhitvector.begin(), munpixelhitvector.end());
  munpixelhitfiller.fill();

  edm::ValueMap<float>::Filler muntrackerlayerfiller(*outputmuntrackerlayermap);
  muntrackerlayerfiller.insert(muonsH, muntrackerlayervector.begin(), muntrackerlayervector.end());
  muntrackerlayerfiller.fill();

  edm::ValueMap<float>::Filler elnvtxfiller(*outputelnvtxmap);
  elnvtxfiller.insert(electronsH, elnvtxvector.begin(), elnvtxvector.end());
  elnvtxfiller.fill();

  edm::ValueMap<float>::Filler elwgtfiller(*outputelwgtmap);
  elwgtfiller.insert(electronsH, elwgtvector.begin(), elwgtvector.end());
  elwgtfiller.fill();

  edm::ValueMap<float>::Filler eldxyfiller(*outputeldxymap);
  eldxyfiller.insert(electronsH, eldxyvector.begin(), eldxyvector.end());
  eldxyfiller.fill();

  edm::ValueMap<float>::Filler eldzfiller(*outputeldzmap);
  eldzfiller.insert(electronsH, eldzvector.begin(), eldzvector.end());
  eldzfiller.fill();

  edm::ValueMap<float>::Filler phnvtxfiller(*outputphnvtxmap);
  phnvtxfiller.insert(photonsH, phnvtxvector.begin(), phnvtxvector.end());
  phnvtxfiller.fill();

  edm::ValueMap<float>::Filler phwgtfiller(*outputphwgtmap);
  phwgtfiller.insert(photonsH, phwgtvector.begin(), phwgtvector.end());
  phwgtfiller.fill();

  iEvent.put(std::move(outputmunvtxmap),"munvtxmap");
  iEvent.put(std::move(outputmuwgtmap), "muwgtmap");
  iEvent.put(std::move(outputmudxymap), "mudxymap");
  iEvent.put(std::move(outputmudzmap),  "mudzmap");
  iEvent.put(std::move(outputmuchi2map),  "muchi2map");
  iEvent.put(std::move(outputmunvalidhitmap),  "munvalidhitmap");
  iEvent.put(std::move(outputmunpixelhitmap),  "munpixelhitmap");
  iEvent.put(std::move(outputmuntrackerlayermap),  "muntrackerlayermap");

  iEvent.put(std::move(outputhltmu20muonrefs),  "hltmu20muonrefs");
  iEvent.put(std::move(outputhlttkmu20muonrefs),"hlttkmu20muonrefs");
  iEvent.put(std::move(outputhltmu22muonrefs),  "hltmu22muonrefs");
  iEvent.put(std::move(outputhlttkmu22muonrefs),"hlttkmu22muonrefs");
  iEvent.put(std::move(outputhltmu24muonrefs),  "hltmu24muonrefs");
  iEvent.put(std::move(outputhlttkmu24muonrefs),"hlttkmu24muonrefs");
  iEvent.put(std::move(outputhltmu50muonrefs),  "hltmu50muonrefs");
  iEvent.put(std::move(outputhlttkmu50muonrefs),"hlttkmu50muonrefs");
  iEvent.put(std::move(outputhltmumuonrefs),    "hltmumuonrefs");
  iEvent.put(std::move(outputhlttkmumuonrefs),  "hlttkmumuonrefs");

  iEvent.put(std::move(outputhltmu17mu8muonrefs_leg17),    "hltmu17mu8Leg17");
  iEvent.put(std::move(outputhltmu17mu8muonrefs_leg8),    "hltmu17mu8Leg8");
  iEvent.put(std::move(outputhltmu17mu8dzmuonrefs_leg17),    "hltmu17mu8dzLeg17");
  iEvent.put(std::move(outputhltmu17mu8dzmuonrefs_leg8),    "hltmu17mu8dzLeg8");
  iEvent.put(std::move(outputhltmu17tkmu8muonrefs_leg17),    "hltmu17tkmu8Leg17");
  iEvent.put(std::move(outputhltmu17tkmu8muonrefs_leg8),    "hltmu17tkmu8Leg8");
  iEvent.put(std::move(outputhltmu17tkmu8dzmuonrefs_leg17),    "hltmu17tkmu8dzLeg17");
  iEvent.put(std::move(outputhltmu17tkmu8dzmuonrefs_leg8),    "hltmu17tkmu8dzLeg8");

  iEvent.put(std::move(outputloosemuonrefs), "loosemuonrefs");
  iEvent.put(std::move(outputtightmuonrefs), "tightmuonrefs");
  iEvent.put(std::move(outputtightmuons),    "tightmuons");

  iEvent.put(std::move(outputelnvtxmap), "elnvtxmap");
  iEvent.put(std::move(outputelwgtmap),  "elwgtmap");
  iEvent.put(std::move(outputeldxymap),  "eldxymap");
  iEvent.put(std::move(outputeldzmap),   "eldzmap");
  iEvent.put(std::move(outputphnvtxmap), "phnvtxmap");
  iEvent.put(std::move(outputphwgtmap),  "phwgtmap");

  // single ele trigger info
  iEvent.put(std::move(outputhltele24eta2p1wplooseelectronrefs),"hltele24eta2p1wplooseelectronrefs");
  iEvent.put(std::move(outputhltele25eta2p1wptightelectronrefs),"hltele25eta2p1wptightelectronrefs");
  iEvent.put(std::move(outputhltele27eta2p1wplooseelectronrefs),"hltele27eta2p1wplooseelectronrefs");
  iEvent.put(std::move(outputhltele27eta2p1wptightelectronrefs),"hltele27eta2p1wptightelectronrefs");
  iEvent.put(std::move(outputhltele27wptightelectronrefs),"hltele27wptightelectronrefs");
  iEvent.put(std::move(outputhltele105electronrefs),      "hltele105electronrefs");
  iEvent.put(std::move(outputhltele115electronrefs),      "hltele115electronrefs");
  iEvent.put(std::move(outputhlthltelelectronrefs),       "hltelelectronrefs");

  //electron id
  iEvent.put(std::move(outputvetoelectronrefs),  "vetoelectronrefs");
  iEvent.put(std::move(outputlooseelectronrefs), "looseelectronrefs");
  iEvent.put(std::move(outputmediumelectronrefs),"mediumelectronrefs");
  iEvent.put(std::move(outputtightelectronrefs), "tightelectronrefs");
  iEvent.put(std::move(outputhltsafeelectronrefs), "hltsafeelectronrefs");
  iEvent.put(std::move(outputmvalooseelectronrefs), "mvalooseelectronrefs");
  iEvent.put(std::move(outputmvatightelectronrefs), "mvatightelectronrefs");
  iEvent.put(std::move(outputtightelectrons),    "tightelectrons");

  //photon id
  iEvent.put(std::move(outputloosephotonrefs), "loosephotonrefs");
  iEvent.put(std::move(outputmediumphotonrefs),"mediumphotonrefs");
  iEvent.put(std::move(outputtightphotonrefs), "tightphotonrefs");
  iEvent.put(std::move(outputmvatightphotonrefs), "mvaloosephotonrefs");
  iEvent.put(std::move(outputmvaloosephotonrefs), "mvatightphotonrefs");
  iEvent.put(std::move(outputrecoelectronmatchrefs),      "recoelectronmatch");

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
