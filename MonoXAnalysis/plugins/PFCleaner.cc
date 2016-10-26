#include <memory>
#include <vector>
#include <iostream>

#include <TRandom3.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

class PFCleaner : public edm::stream::EDProducer<> {
public:
  explicit PFCleaner(const edm::ParameterSet&);
  ~PFCleaner();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob();
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  
  //livia for purity studies  
  bool   randomConeOverlaps(double, double, double, std::vector<pat::Jet>, double);
  double computePhotonIso( edm::Handle<edm::View<reco::Candidate>> ,double, double , double);
  bool   testPhotonIsolationLoose(std::vector<pat::Photon>::const_iterator  & photon, float chargedHadronIsolation, float neutralHadronIsolation,float rhoval, float hoe);
  bool   testPhotonIsolationMedium(std::vector<pat::Photon>::const_iterator & photon, float chargedHadronIsolation, float rhoval, float hoe);
  double getGammaEAForPhotonIso(double eta);
  double getNeutralHadronEAForPhotonIso(double eta);
  double getChargedHadronEAForPhotonIso(double eta);
  bool   isPassingPhotonHighPtID(double, double, double, double, double, double, double);

  // rho value
  const edm::EDGetTokenT<double>  rhoToken;
  // vertex collection
  const edm::EDGetTokenT<std::vector<reco::Vertex> >  verticesToken;
  // pf candidate collection
  const edm::EDGetTokenT<edm::View<reco::Candidate> > pfcandsToken;
  // jets input collection
  const edm::EDGetTokenT<std::vector<pat::Jet> >  jetsToken;
  // muon input collection
  const edm::EDGetTokenT<std::vector<pat::Muon> > muonsToken;
  const std::vector<edm::ParameterSet > muonSelection;
  // tau input collection
  const edm::EDGetTokenT<std::vector<pat::Tau> >  tausToken;
  const std::vector<edm::ParameterSet > tauSelection;

  // electron input collection
  const edm::EDGetTokenT<std::vector<pat::Electron> > electronsToken;
  const std::vector<edm::ParameterSet > electronSelection;
  std::vector<edm::EDGetTokenT<edm::ValueMap<bool> > > electronIdMapToken;
  // calibrated Pat electrons
  edm::EDGetTokenT<std::vector<pat::Electron> > calibratedElectronsToken;
  const bool useCalibratedElectrons;

  // photon input collection
  const edm::EDGetTokenT<std::vector<pat::Photon> > photonsToken;
  std::vector<edm::EDGetTokenT<edm::ValueMap<bool> > > photonIdMapToken;
  const std::vector<edm::ParameterSet > photonSelection;
  edm::EDGetTokenT<std::vector<pat::Photon> > calibratedPhotonsToken;
  const bool useCalibratedPhotons;

  const edm::EDGetTokenT<edm::ValueMap<float> > photonsieieToken;
  const edm::EDGetTokenT<edm::ValueMap<float> > photonPHisoToken;
  const edm::EDGetTokenT<edm::ValueMap<float> > photonCHisoToken;
  // calibrated Pat photons
  const edm::EDGetTokenT<std::vector<pat::Photon> > calibratedPhotonToken;
  const edm::ParameterSet  highPtPhotonID;
  const edm::ParameterSet  loosePhotonID;

  const bool userandomphi;
  const bool addPhotonPurity;
  int looseMuonPosition;
  int vetoElectronPosition;
  TRandom3 rand;

};

PFCleaner::PFCleaner(const edm::ParameterSet& iConfig): 
  rhoToken                 (consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  verticesToken            (consumes<std::vector<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("vertices"))),
  pfcandsToken             (consumes<edm::View<reco::Candidate> > (iConfig.getParameter<edm::InputTag>("pfcands"))),
  jetsToken                (consumes<std::vector<pat::Jet> > (iConfig.getParameter<edm::InputTag>("jets"))),
  muonsToken               (consumes<std::vector<pat::Muon> > (iConfig.getParameter<edm::InputTag>("muons"))), 
  muonSelection            (iConfig.getParameter<std::vector<edm::ParameterSet > >("muonSelection")),
  tausToken                (consumes<std::vector<pat::Tau> >  (iConfig.getParameter<edm::InputTag>("taus"))),
  tauSelection             (iConfig.getParameter<std::vector<edm::ParameterSet> >("tauSelection")),
  electronsToken           (consumes<std::vector<pat::Electron> > (iConfig.getParameter<edm::InputTag>("electrons"))),
  electronSelection        (iConfig.getParameter<std::vector<edm::ParameterSet> >("electronSelection")),
  useCalibratedElectrons   (iConfig.existsAs<bool>("useCalibratedElectrons") ? iConfig.getParameter<bool>("useCalibratedElectrons") : false),
  photonsToken             (consumes<std::vector<pat::Photon> > (iConfig.getParameter<edm::InputTag>("photons"))),
  photonSelection          (iConfig.getParameter<std::vector<edm::ParameterSet> >("photonSelection")),
  useCalibratedPhotons     (iConfig.existsAs<bool>("useCalibratedPhotons") ? iConfig.getParameter<bool>("useCalibratedPhotons") : false),
  photonsieieToken         (consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonsieie"))), 
  photonPHisoToken         (consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonphiso"))), 
  photonCHisoToken         (consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonchiso"))),   
  highPtPhotonID           (iConfig.getParameter<edm::ParameterSet>("highPtPhotonID")),
  loosePhotonID            (iConfig.getParameter<edm::ParameterSet>("loosePhotonID")),
  userandomphi             (iConfig.existsAs<bool>("userandomphiforRC") ? iConfig.getParameter<bool>("userandomphiforRC") : true),
  addPhotonPurity          (iConfig.existsAs<bool>("addPhotonPurity") ? iConfig.getParameter<bool>("addPhotonPurity") : false)

{

  rand.SetSeed(0);
  
  if(useCalibratedElectrons)
    calibratedElectronsToken = consumes<std::vector<pat::Electron> > (iConfig.getParameter<edm::InputTag>("calibratedElectrons"));
  
  if(useCalibratedPhotons)
    calibratedPhotonsToken   = consumes<std::vector<pat::Photon> >   (iConfig.getParameter<edm::InputTag>("calibratedPhotons"));

  // produces muon output
  for(auto imuon : muonSelection)
    produces<pat::MuonRefVector>(imuon.getParameter<std::string>("muonCollectionName"));
  // produces tau output
  for(auto itau : tauSelection)
    produces<pat::TauRefVector>(itau.getParameter<std::string>("tauCollectionName"));
  
  // produces ele + get value map
  for (auto iele : electronSelection){
    electronIdMapToken.push_back(consumes<edm::ValueMap<bool> > (iele.getParameter<edm::InputTag>("eleValueMap")));
    produces<pat::ElectronRefVector>(iele.getParameter<std::string>("electronCollectionName"));
  }

  // produces photon + get value map
  for (auto ipho : photonSelection){
    photonIdMapToken.push_back(consumes<edm::ValueMap<bool> > (ipho.getParameter<edm::InputTag>("photonValueMap")));
    produces<pat::PhotonRefVector>(ipho.getParameter<std::string>("photonCollectionName"));
  }
  
  // random cone information
  if(addPhotonPurity){
    produces<edm::ValueMap<float> >("rndchhadiso04");
    produces<edm::ValueMap<float> >("rndchhadiso08");
    produces<edm::ValueMap<float> >("rndgammaiso04");
    produces<edm::ValueMap<float> >("rndgammaiso08");
    produces<pat::PhotonRefVector> ("photonsPurity");
    produces<pat::PhotonRefVector> ("tightphotonsPurity");
  }
  // value map for highptId
  produces<edm::ValueMap<bool> > ("photonHighPtId");
  produces<edm::ValueMap<bool> > ("photonLooseId");
  produces<edm::ValueMap<float> >("gammaiso");
  produces<edm::ValueMap<float> >("chhadiso");

  // useful for taus
  looseMuonPosition     = -1;
  vetoElectronPosition  = -1;
  for(size_t imuon = 0; imuon < muonSelection.size(); imuon++){
    if(muonSelection.at(imuon).getParameter<std::string>("idType") == "loose" or muonSelection.at(imuon).getParameter<std::string>("idType") == "Loose")
      looseMuonPosition = imuon;
  }
  
  if(looseMuonPosition < 0 or looseMuonPosition == int(muonSelection.size()))
    throw cms::Exception("PFCleaner") <<" no loose muons found  --> check \n";

  for(size_t iele = 0; iele < electronSelection.size(); iele++){
    if(electronSelection.at(iele).getParameter<std::string>("idType") == "veto" or electronSelection.at(iele).getParameter<std::string>("idType") == "Veto")
      vetoElectronPosition = iele;
  }
  
  if(vetoElectronPosition < 0 or vetoElectronPosition == int(electronSelection.size()))
    throw cms::Exception("PFCleaner") <<" no veto electrons found  --> check \n";

}


PFCleaner::~PFCleaner() {}

void PFCleaner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;
  
  Handle<std::vector<reco::Vertex> > verticesH;
  iEvent.getByToken(verticesToken, verticesH);

  Handle<edm::View<reco::Candidate> > pfcandsH;
  iEvent.getByToken(pfcandsToken, pfcandsH);
  
  Handle<std::vector<pat::Jet> > jetsH;
  iEvent.getByToken(jetsToken, jetsH);
  
  // muons
  Handle<std::vector<pat::Muon> > muonsH;
  iEvent.getByToken(muonsToken, muonsH);
  
  // electrons
  Handle<std::vector<pat::Electron> > electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  Handle<std::vector<pat::Electron> > calibratedElectronsH;
  if(useCalibratedElectrons)
    iEvent.getByToken(calibratedElectronsToken, calibratedElectronsH);
  std::vector<Handle<edm::ValueMap<bool> > > electronMapH;
  electronMapH.reserve(electronIdMapToken.size());
  for(auto eleToken :  electronIdMapToken){
    electronMapH.push_back(Handle<edm::ValueMap<bool> > ());
    iEvent.getByToken(eleToken,electronMapH.back());
  }
  
  // taus
  Handle<std::vector<pat::Tau> > tausH;
  iEvent.getByToken(tausToken, tausH);

  // photons
  Handle<std::vector<pat::Photon> > photonsH;
  iEvent.getByToken(photonsToken, photonsH);
  Handle<std::vector<pat::Photon> > calibratedPhotonsH;
  if(useCalibratedPhotons)
    iEvent.getByToken(calibratedPhotonsToken, calibratedPhotonsH);
  std::vector<Handle<edm::ValueMap<bool> > > photonsMapH;
  for(auto phoToken :  photonIdMapToken){
    photonsMapH.push_back(Handle<edm::ValueMap<bool> > ());
    iEvent.getByToken(phoToken,photonsMapH.back());
  }

  Handle<edm::ValueMap<float> > photonsieieH;
  iEvent.getByToken(photonsieieToken, photonsieieH);  
  Handle<edm::ValueMap<float> > photonPHisoH;
  iEvent.getByToken(photonPHisoToken, photonPHisoH);
  Handle<edm::ValueMap<float> > photonCHisoH;
  iEvent.getByToken(photonCHisoToken, photonCHisoH);

  Handle<double> rhoH;
  iEvent.getByToken(rhoToken, rhoH);
  double rho = *rhoH;

  // output
  std::vector<std::auto_ptr<pat::MuonRefVector> > outputmuons;
  for(auto imuon : muonSelection)
    outputmuons.push_back(std::auto_ptr<pat::MuonRefVector>(new pat::MuonRefVector));
  std::vector<std::auto_ptr<pat::ElectronRefVector> > outputelectrons;
  for(auto iele  : electronSelection)
    outputelectrons.push_back(std::auto_ptr<pat::ElectronRefVector>(new pat::ElectronRefVector));
  std::vector<std::auto_ptr<pat::TauRefVector> > outputtaus;
  for(auto itau  : tauSelection)
    outputtaus.push_back(std::auto_ptr<pat::TauRefVector>(new pat::TauRefVector));
  std::vector<std::auto_ptr<pat::PhotonRefVector> > outputphotons;
  for(auto ipho  : photonSelection)
    outputphotons.push_back(std::auto_ptr<pat::PhotonRefVector>(new pat::PhotonRefVector));
  
  std::auto_ptr<edm::ValueMap<bool>  > outputphotonhighptidmap(new ValueMap<bool>());
  std::auto_ptr<edm::ValueMap<bool>  > outputloosephotonmap(new ValueMap<bool>());

  std::auto_ptr<edm::ValueMap<float> > outputgammaisomap(new ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > outputchhadisomap(new ValueMap<float>());

  std::auto_ptr<edm::ValueMap<float> > outputrndgammaiso04map(new ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > outputrndgammaiso08map(new ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > outputrndchhadiso04map(new ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > outputrndchhadiso08map(new ValueMap<float>());

  std::auto_ptr<pat::PhotonRefVector>  outputphotonsPurity(new pat::PhotonRefVector);
  std::auto_ptr<pat::PhotonRefVector>  outputtightphotonsPurity(new pat::PhotonRefVector);
  
  //muon info https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
  for (vector<pat::Muon>::const_iterator muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
    if (verticesH->size() == 0) continue;

    // loop on the muon selection definition
    size_t ipos = 0;
    for(auto imuon : muonSelection){      
      bool passeskincuts = (muons_iter->pt() > imuon.getParameter<double>("ptMin") && 
			    fabs(muons_iter->eta()) < imuon.getParameter<double>("absEta"));
      float isoval       = std::max(0.,muons_iter->pfIsolationR04().sumNeutralHadronEt+muons_iter->pfIsolationR04().sumPhotonEt 
				    -imuon.getParameter<double>("deltaBeta")*muons_iter->pfIsolationR04().sumPUPt);
      isoval += muons_iter->pfIsolationR04().sumChargedHadronPt;
      isoval /= muons_iter->pt();

      if (passeskincuts) {
	if(imuon.getParameter<std::string>("idType") == "loose" or imuon.getParameter<std::string>("idType") =="Loose"){
	  if (muon::isLooseMuon(*muons_iter) && isoval < imuon.getParameter<double>("isolation")) 
	    outputmuons.at(ipos)->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
	}
	else if(imuon.getParameter<std::string>("idType") == "tight" or imuon.getParameter<std::string>("idType") =="Tight"){
	  if (muon::isTightMuon(*muons_iter,verticesH->at(0)) && isoval < imuon.getParameter<double>("isolation")) 
	    outputmuons.at(ipos)->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
	}
	else if(imuon.getParameter<std::string>("idType") == "highPt" or imuon.getParameter<std::string>("idType") =="HighPt"){
	  if (muon::isHighPtMuon(*muons_iter,verticesH->at(0)) && isoval < imuon.getParameter<double>("isolation")) 
	    outputmuons.at(ipos)->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
	}
      }
      ipos++;
    }
  }
  
  // electrons --> id applied on standard slimmed, in case of calibrated use them to fill the output collection
  //electron info https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2
  for (vector<pat::Electron>::const_iterator electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) {
    const Ptr<pat::Electron> electronPtr(electronsH, electrons_iter - electronsH->begin());
    size_t ipos = 0;
    for(auto iele : electronSelection){

      // bool for the electron id
      bool passesid = (*electronMapH.at(ipos))[electronPtr];      
      // bool for dxy and dz cut that are taken out from standard electron VID
      bool pass_dxy = false;
      bool pass_dz  = false;
      if(verticesH->size() > 0){
	const reco::Vertex & vtx  = verticesH->at(0);
	// barrel
	if(fabs(electrons_iter->superCluster()->eta()) < 1.5 and 
	   fabs(electrons_iter->gsfTrack()->dxy(vtx.position())) < iele.getParameter<edm::ParameterSet>("PVSelection").getParameter<double>("d0Barrel"))
	  pass_dxy = true;
	// endcalp
	else if(fabs(electrons_iter->superCluster()->eta()) > 1.5 and 
		fabs(electrons_iter->gsfTrack()->dxy(vtx.position())) < iele.getParameter<edm::ParameterSet>("PVSelection").getParameter<double>("d0Endcap"))
	  pass_dxy = true;
	// barrel
	if(fabs(electrons_iter->superCluster()->eta()) < 1.5 and 
	   fabs(electrons_iter->gsfTrack()->dz(vtx.position())) < iele.getParameter<edm::ParameterSet>("PVSelection").getParameter<double>("dzBarrel"))
	  pass_dz = true;
	//endcap
	else if(fabs(electrons_iter->superCluster()->eta()) > 1.5 and 
		fabs(electrons_iter->gsfTrack()->dz(vtx.position())) < iele.getParameter<edm::ParameterSet>("PVSelection").getParameter<double>("dzEndcap"))
	  pass_dz = true;
      }
      else{
	//barrel
	if(fabs(electrons_iter->superCluster()->eta()) < 1.5 and 
	   fabs(electrons_iter->gsfTrack()->dxy()) < iele.getParameter<edm::ParameterSet>("PVSelection").getParameter<double>("d0Barrel"))
	  pass_dxy = true;
	//endcap
	else if(fabs(electrons_iter->superCluster()->eta()) > 1.5 and 
		fabs(electrons_iter->gsfTrack()->dxy()) < iele.getParameter<edm::ParameterSet>("PVSelection").getParameter<double>("d0Endcap"))
	  pass_dxy = true;
	
	//barrel
	if(fabs(electrons_iter->superCluster()->eta()) < 1.5 and 
	   fabs(electrons_iter->gsfTrack()->dz()) < iele.getParameter<edm::ParameterSet>("PVSelection").getParameter<double>("dzBarrel"))
	  pass_dz = true;
	//endcap
	else if(fabs(electrons_iter->superCluster()->eta()) > 1.5 and 
		fabs(electrons_iter->gsfTrack()->dz()) < iele.getParameter<edm::ParameterSet>("PVSelection").getParameter<double>("dzEndcap"))
	  pass_dz = true;
      }



      if(passesid and pass_dxy and pass_dz){
	// match with calibrate electrons if possible
	if(calibratedElectronsH.isValid()){ // check if it exists
	  for(vector<pat::Electron>::const_iterator calibele_iter = calibratedElectronsH->begin();
	      calibele_iter != calibratedElectronsH->end(); ++calibele_iter){
	    // match by reference
	    if(calibele_iter->pfCandidateRef()     == electrons_iter->pfCandidateRef() or 
	       calibele_iter->core()               == electrons_iter->core() or
	       calibele_iter->originalObjectRef()  == electrons_iter->originalObjectRef()){
	      
	      bool passeskincuts  = (calibele_iter->pt() > iele.getParameter<double>("ptMin") &&
				     fabs(calibele_iter->superCluster()->eta()) < iele.getParameter<double>("absEta"));
	      
	      if(not passeskincuts) continue;
	      outputelectrons.at(ipos)->push_back(pat::ElectronRef(calibratedElectronsH,calibele_iter-calibratedElectronsH->begin()));
	    }
	  }
	}
	else{	  
	  // apply basic kinematic cuts
	  bool passeskincuts  = (electrons_iter->pt() > iele.getParameter<double>("ptMin") &&
				 fabs(electrons_iter->superCluster()->eta()) < iele.getParameter<double>("absEta"));
	  if(not passeskincuts) continue;
	  outputelectrons.at(ipos)->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
	}
      }
      ipos++;
    }
  }
  
  // taus
  for (vector<pat::Tau>::const_iterator taus_iter = tausH->begin(); taus_iter != tausH->end(); ++taus_iter) {
    if (verticesH->size() == 0) continue;
        
    size_t ipos = 0;
    for(auto itau : tauSelection){
      // clean tau candidates from identified muons and electrons (loosely identified)
      bool skiptau = false;

      for (std::size_t j = 0; j < outputmuons.at(looseMuonPosition)->size(); j++) {
	if (deltaR(outputmuons.at(looseMuonPosition)->at(j)->eta(), outputmuons.at(looseMuonPosition)->at(j)->phi(), 
		   taus_iter->eta(), taus_iter->phi()) < itau.getParameter<double>("dRCleaning")) skiptau = true;
      }
      for (std::size_t j = 0; j < outputelectrons.at(vetoElectronPosition)->size(); j++) {
	if (deltaR(outputelectrons.at(vetoElectronPosition)->at(j)->eta(), outputelectrons.at(vetoElectronPosition)->at(j)->phi(), 
		   taus_iter->eta(), taus_iter->phi()) < itau.getParameter<double>("dRCleaning")) skiptau = true;
      }

      // apply loose id and store vector
      string decayMode = "decayModeFinding";
      if(itau.getParameter<bool>("useNewDecayMode"))
	decayMode = "decayModeFindingNewDMs";
      
      if(itau.getParameter<bool>("graterThan") == true){
	if (taus_iter->pt() > itau.getParameter<double>("ptMin") &&
	    fabs(taus_iter->eta()) < itau.getParameter<double>("absEta") &&
	    taus_iter->tauID(decayMode) > itau.getParameter<double>("decayModeFinding") && 	    
	    taus_iter->tauID(itau.getParameter<std::string>("tauIDName")) > itau.getParameter<double>("isolation") && !skiptau){	
	  outputtaus.at(ipos)->push_back(pat::TauRef(tausH, taus_iter - tausH->begin()));
	}
      }
      else{
	if (taus_iter->pt() > itau.getParameter<double>("ptMin") &&
	    fabs(taus_iter->eta()) < itau.getParameter<double>("absEta") &&
	    taus_iter->tauID(decayMode) > itau.getParameter<double>("decayModeFinding") &&
	    taus_iter->tauID(itau.getParameter<std::string>("tauIDName")) < itau.getParameter<double>("isolation") && !skiptau){	
	  outputtaus.at(ipos)->push_back(pat::TauRef(tausH, taus_iter - tausH->begin()));
	}
      }
      ipos++;      
    }
  }

  // photon https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2
  std::vector<float> rndgammaiso04;
  std::vector<float> rndgammaiso08;
  std::vector<float> gammaiso;
  std::vector<float> chhadiso;
  std::vector<float> rndchhadiso04;
  std::vector<float> rndchhadiso08;
  std::vector<bool>  photonidhighpt;
  std::vector<bool>  photonloose;
  
  // loop on the photon colection
  for (vector<pat::Photon>::const_iterator photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {

    float gaisoval = computePhotonIso(pfcandsH,photons_iter->eta(),  photons_iter->phi(),0.3);     
    gammaiso.push_back(gaisoval);
    float chisoval = 0;
    for(size_t i = 0; i < pfcandsH->size(); i++) {
      const auto& pfcand = pfcandsH->ptrAt(i);
      if (abs(pfcand->pdgId()) == 211 && deltaR(photons_iter->eta(),photons_iter->phi(), pfcand->eta(), pfcand->phi()) <= 0.3)
	chisoval += pfcand->pt();
    }
    chhadiso.push_back(chisoval);
    
    if(addPhotonPurity){
      double rndphi04 = photons_iter->phi() + M_PI/2.0;
      double rndphi08 = photons_iter->phi() + M_PI/2.0;
      float  chisoval04 = 0.;
      float  chisoval08 = 0.;
      unsigned int counterR04 = 1;
      unsigned int counterR08 = 1;
      if (userandomphi) {
	rndphi04 = rand.Uniform(-M_PI, M_PI);
	while (randomConeOverlaps(rndphi04, photons_iter->eta(), photons_iter->phi(), *jetsH, 0.4)) {
	  rndphi04 = rand.Uniform(-M_PI, M_PI);
	  counterR04++;
	  if(counterR04 > 5000){
	    rndphi04 = -99;
	    continue;
	  }
	}
      }
      if (userandomphi){
	rndphi08 = rand.Uniform(-M_PI, M_PI);
	while (randomConeOverlaps(rndphi08, photons_iter->eta(), photons_iter->phi(), *jetsH, 0.8)){
	  rndphi08 = rand.Uniform(-M_PI, M_PI);
	  counterR08++;
	  if(counterR08 > 5000){
	    rndphi08 = -99;
	    continue;
	  }
	}
      }

      if(rndphi04 != -99){
	for(size_t i = 0; i < pfcandsH->size(); i++) {
	  const auto& pfcand = pfcandsH->ptrAt(i);
	  if (abs(pfcand->pdgId()) == 211 && deltaR(photons_iter->eta(),rndphi04, pfcand->eta(), pfcand->phi()) <= 0.3) 
	    chisoval04 += pfcand->pt();
	}
      }
      else
	chisoval04 = -99;
      
      if(rndphi08 != -99){
	for(size_t i = 0; i < pfcandsH->size(); i++) {
	  const auto& pfcand = pfcandsH->ptrAt(i);
	  if (abs(pfcand->pdgId()) == 211 && deltaR(photons_iter->eta(),rndphi08, pfcand->eta(), pfcand->phi()) <= 0.3) 
	    chisoval08 += pfcand->pt();
	}      
      }
      else chisoval08 = -99;

      float rndgaisoval04 = 0;
      float rndgaisoval08 = 0;
      if(rndphi04 != -99)
	rndgaisoval04 = computePhotonIso(pfcandsH,photons_iter->eta(), rndphi04, 0.3);      
      else 
	rndgaisoval04 = -99;
      if(rndphi08 != -99)
	rndgaisoval08 = computePhotonIso(pfcandsH,photons_iter->eta(), rndphi08, 0.3);      
      else
	rndgaisoval08 = -99;
      rndgammaiso04.push_back(rndgaisoval04);
      rndgammaiso08.push_back(rndgaisoval08);
      rndchhadiso04.push_back(chisoval04);
      rndchhadiso08.push_back(chisoval08);
    }

    //livia
    const Ptr<pat::Photon> photonPtr(photonsH, photons_iter - photonsH->begin());
    double photonsieie = (*photonsieieH)[photonPtr];
    double photonphiso = (*photonPHisoH)[photonPtr];
    double photonchiso = (*photonCHisoH)[photonPtr];
    double photonhoe   =  photons_iter->hadTowOverEm();
    // high pt photon ID as described in AN (2016/079)
    bool  passesphotonidhighpt = false;
    passesphotonidhighpt = isPassingPhotonHighPtID(photons_iter->pt(), photons_iter->superCluster()->eta(), photonhoe, photonchiso, photonphiso, photonsieie, rho);
    photonidhighpt.push_back(passesphotonidhighpt);

    // additional loose photon definition
    if (fabs(photons_iter->superCluster()->eta()) < loosePhotonID.getParameter<double>("absEta") and
	photons_iter->pt() > loosePhotonID.getParameter<double>("ptMin") and
	(photons_iter->r9() > loosePhotonID.getParameter<double>("R9min") or photons_iter->chargedHadronIso() < loosePhotonID.getParameter<double>("chIso") or photons_iter->chargedHadronIso() < photons_iter->pt()*loosePhotonID.getParameter<double>("chIsoFrac")))    
      photonloose.push_back(true);
    else
      photonloose.push_back(false);

    // Purity studies
    if(addPhotonPurity){
      bool passeslooseidPurity = testPhotonIsolationLoose(photons_iter, photons_iter->chargedHadronIso(), photons_iter->neutralHadronIso(), rho, photons_iter->hadTowOverEm());
      if (passeslooseidPurity) {
	outputphotonsPurity->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
	if (photons_iter->pt() > 150) outputtightphotonsPurity->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
      }
    }

    size_t ipos = 0;
    for(auto ipho : photonSelection){
      // check only photons with the following proprierties
      // standard ID
      bool passesid = (*photonsMapH.at(ipos))[photonPtr];
      if(passesid){	
	// match with calibrate photons if possible                                                                                                                      
	if(calibratedPhotonsH.isValid()){ // check if it exists                                                                                                          
	  if(calibratedPhotonsH->size() != photonsH->size())
	    throw cms::Exception("PFCleaner") <<" different size between calibrated and un-calibrated photons --> check \n";
	  bool isMatched = false;
	  for(vector<pat::Photon>::const_iterator calibpho_iter = calibratedPhotonsH->begin();
	      calibpho_iter != calibratedPhotonsH->end(); ++calibpho_iter){
	    // match by reference                                                                                                                                        
	    if(calibpho_iter->superCluster()       == photons_iter->superCluster() or
	       calibpho_iter->originalObjectRef()  == photons_iter->originalObjectRef()){
	      isMatched = true;
	      if (fabs(calibpho_iter->superCluster()->eta()) > ipho.getParameter<double>("absEta") or 
		  calibpho_iter->pt() < ipho.getParameter<double>("ptMin") or 
		  not calibpho_iter->passElectronVeto()) continue;	      
	      outputphotons.at(ipos)->push_back(pat::PhotonRef(calibratedPhotonsH,calibpho_iter-calibratedPhotonsH->begin()));
	    }
	  }
	  if(not isMatched)
	    throw cms::Exception("PFCleaner") <<" missing matching for one photons between calib and un-calib collections --> check \n";
	}
	else{
	  if (fabs(photons_iter->superCluster()->eta()) > ipho.getParameter<double>("absEta") or 
	      photons_iter->pt() < ipho.getParameter<double>("ptMin") or 
	      not photons_iter->passElectronVeto()) continue;
	  outputphotons.at(ipos)->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
	}
      }
      ipos++;
    }    
  }

  if(addPhotonPurity){
    edm::ValueMap<float>::Filler rnd04gafiller(*outputrndgammaiso04map);
    rnd04gafiller.insert(photonsH, rndgammaiso04.begin(), rndgammaiso04.end());
    rnd04gafiller.fill();
    edm::ValueMap<float>::Filler rnd08gafiller(*outputrndgammaiso08map);
    rnd08gafiller.insert(photonsH, rndgammaiso08.begin(), rndgammaiso08.end());
    rnd08gafiller.fill();
    
    edm::ValueMap<float>::Filler rnd04chfiller(*outputrndchhadiso04map);
    rnd04chfiller.insert(photonsH, rndchhadiso04.begin(), rndchhadiso04.end());
    rnd04chfiller.fill();

    edm::ValueMap<float>::Filler rnd08chfiller(*outputrndchhadiso08map);
    rnd08chfiller.insert(photonsH, rndchhadiso08.begin(), rndchhadiso08.end());
    rnd08chfiller.fill();
  }
  
  edm::ValueMap<float>::Filler gafiller(*outputgammaisomap);
  gafiller.insert(photonsH, gammaiso.begin(), gammaiso.end());
  gafiller.fill();

  edm::ValueMap<float>::Filler chfiller(*outputchhadisomap);
  chfiller.insert(photonsH, chhadiso.begin(), chhadiso.end());
  chfiller.fill();


  //livia    
  edm::ValueMap<bool>::Filler photonhighptidfiller(*outputphotonhighptidmap);
  if(calibratedPhotonsH.isValid())
    photonhighptidfiller.insert(calibratedPhotonsH, photonidhighpt.begin(), photonidhighpt.end());
  else
    photonhighptidfiller.insert(photonsH, photonidhighpt.begin(), photonidhighpt.end());
  photonhighptidfiller.fill();

  //livia    
  edm::ValueMap<bool>::Filler photonlooseidfiller(*outputloosephotonmap);
  if(calibratedPhotonsH.isValid())
    photonlooseidfiller.insert(calibratedPhotonsH, photonloose.begin(), photonloose.end());
  else
    photonlooseidfiller.insert(photonsH, photonloose.begin(), photonloose.end());
  photonlooseidfiller.fill();

  iEvent.put(outputphotonhighptidmap,"photonHighPtId");
  iEvent.put(outputloosephotonmap,   "photonLooseId");
  iEvent.put(outputgammaisomap,      "gammaiso");  
  iEvent.put(outputchhadisomap,      "chhadiso");  
  
  if(addPhotonPurity){
    iEvent.put(outputphotonsPurity,    "photonsPurity");
    iEvent.put(outputtightphotonsPurity,"tightphotonsPurity");
    iEvent.put(outputrndgammaiso04map, "rndgammaiso04");
    iEvent.put(outputrndgammaiso08map, "rndgammaiso08");
    iEvent.put(outputrndchhadiso04map, "rndchhadiso04");
    iEvent.put(outputrndchhadiso08map, "rndchhadiso08");
  }

  for(size_t imuon = 0; imuon < outputmuons.size(); imuon++)
    iEvent.put(outputmuons.at(imuon), muonSelection.at(imuon).getParameter<std::string>("muonCollectionName"));

  for(size_t ielectron = 0; ielectron < outputelectrons.size(); ielectron++)
    iEvent.put(outputelectrons.at(ielectron), electronSelection.at(ielectron).getParameter<std::string>("electronCollectionName"));

  for(size_t itau = 0; itau < outputtaus.size(); itau++)
    iEvent.put(outputtaus.at(itau), tauSelection.at(itau).getParameter<std::string>("tauCollectionName"));

  for(size_t iphoton = 0; iphoton < outputphotons.size(); iphoton++)
    iEvent.put(outputphotons.at(iphoton), photonSelection.at(iphoton).getParameter<std::string>("photonCollectionName"));

}

void PFCleaner::beginJob() {
}

void PFCleaner::endJob() {
}

void PFCleaner::beginRun(edm::Run const&, edm::EventSetup const&) {
}

void PFCleaner::endRun(edm::Run const&, edm::EventSetup const&) {
}

void PFCleaner::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void PFCleaner::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void PFCleaner::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// Photon purity studies
bool PFCleaner::randomConeOverlaps(double randomphi, double photoneta, double photonphi, std::vector<pat::Jet> jets, double vetoCone) {
  if (reco::deltaR(photoneta, randomphi, photoneta, photonphi) < vetoCone) return true;
  for (std::size_t i = 0; i < jets.size(); i++) {    
    if (reco::deltaR(photoneta,photonphi,jets[i].eta(),jets[i].phi()) < 0.4) continue;
    if (jets[i].pt() > 30. && reco::deltaR(photoneta, randomphi, jets[i].eta(), jets[i].phi()) < vetoCone){
      return true;
    }
  }
  return false;
}

double PFCleaner::computePhotonIso(edm::Handle<edm::View<reco::Candidate> >  pfcandsH, double eta, double phi, double isoCone){
  double isoval=0;
  for(size_t i = 0; i < pfcandsH->size(); i++) {
    const auto& pfcand = pfcandsH->ptrAt(i);
    if (pfcand->pdgId()  ==  22 && deltaR(eta, phi, pfcand->eta(), pfcand->phi()) <= isoCone) isoval += pfcand->pt();
  }
  return isoval;
}

bool PFCleaner::testPhotonIsolationMedium(std::vector<pat::Photon>::const_iterator& photon, float chargedHadronIsolation, float rhoval, float hoe) {
  if (photon->isEB()) {
    if (chargedHadronIsolation<10 && hoe < 0.05) return true;
    else return false;
  } 
  else return false;
}

bool PFCleaner::testPhotonIsolationLoose(std::vector<pat::Photon>::const_iterator& photon, float chargedHadronIsolation,float neutralHadronIsolation, float rhoval, float hoe) { //cut only on NH and HoE. Cuts on CH and PH and sieie will be applied later for purity studies
  double corrNHIso = neutralHadronIsolation - rhoval * getNeutralHadronEAForPhotonIso(photon->eta());
  if (photon->isEB()) {
    if (corrNHIso < 1.92 + 0.014 * photon->pt() + 0.000019 * pow ( photon->pt(),2) && hoe < 0.05) return true;
    else return false;
  }
  else if (photon->isEE()) {//we care only about EB
    if (chargedHadronIsolation<10 && hoe<0.05) return true;
    else return false;
  }
  else return false;
}


// check the high pt id
bool PFCleaner::isPassingPhotonHighPtID(double pt, double eta, double hoe, double chiso, double phiso, double sieie, double rho ){

    double chisoCUT = 0;
    double phisoCUT = 0;
    double sieieCUT = 0;
    double hoeCUT   = 0; 
    double alpha    = 0;
    double k        = 0;
    double EA       = 0;
    double newphiso = 0;
    int isPassing = false;
    if(abs(eta)<highPtPhotonID.getParameter<double>("absEta") && pt > highPtPhotonID.getParameter<double>("ptMin")){
      chisoCUT = highPtPhotonID.getParameter<double>("chIso");
      sieieCUT = highPtPhotonID.getParameter<double>("sigmaIetaIeta");
      hoeCUT   = highPtPhotonID.getParameter<double>("HOverE");
      phisoCUT = highPtPhotonID.getParameter<double>("isolation");
      alpha    = highPtPhotonID.getParameter<double>("alpha");
      k        = highPtPhotonID.getParameter<double>("k");
      if(abs(eta)<0.9) EA=0.17;
      if(abs(eta)>0.9 && abs(eta)<1.4442) EA=0.14;
    }
    newphiso = alpha+phiso-rho*EA-k*pt;
    if(newphiso < phisoCUT && chiso < chisoCUT && sieie < sieieCUT && hoe < hoeCUT)
      isPassing =true;    
    return isPassing;
}

double PFCleaner::getChargedHadronEAForPhotonIso(double eta) {
  if (fabs(eta) < 1.0) return 0.000000000001;
  else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return  0.000000000001;
  else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return  0.000000000001;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return  0.000000000001;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return  0.000000000001;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.3  ) return  0.000000000001;
  else if (fabs(eta) >= 2.4) return 0.000000000001 ;
  else return 0.;
}

double PFCleaner::getNeutralHadronEAForPhotonIso(double eta) {
  if (fabs(eta) < 1.0) return 0.0599;
  else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.0819;
  else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.0696;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.036;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.036;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.3  ) return 0.0462;
  else if (fabs(eta) >= 2.4) return 0.0656;
  else return 0.;
}

double PFCleaner::getGammaEAForPhotonIso(double eta) {
  if (fabs(eta) < 1.0) return 0.1271;
  else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.1101;
  else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.0756;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.1175;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.1498;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.3  ) return 0.1857;
  else if (fabs(eta) >= 2.4) return 0.2183;
  else return 0.;
}



DEFINE_FWK_MODULE(PFCleaner);
