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
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {
    if( itr->key() == theCandidate.key() ) return true;
  }
  return false;
}

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

  
  //for purity
  bool   randomConeOverlaps(const double &, const double &, const double &, const std::vector<pat::Jet> &, const double &);
  // for isolation
  double computePhotonIso( const pat::Photon &, edm::Handle<edm::View<reco::Candidate>> , const double &, const double &, const double &);
  double computeCHhadronIso(const pat::Photon &, edm::Handle<edm::View<reco::Candidate> > , const double &, const double &, const double &, const reco::Vertex &);
  double computeNHhadronIso(const pat::Photon &, edm::Handle<edm::View<reco::Candidate> > , const double &, const double &, const double &);  
  // effective area according to 80X Spring 16 ID
  double getGammaEAForPhotonIso(const double & eta);
  double getNeutralHadronEAForPhotonIso(const double & eta);
  double getChargedHadronEAForPhotonIso(const double & eta);

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
  const edm::EDGetTokenT<edm::ValueMap<float> > photonNHisoToken;
  // calibrated Pat photons
  const edm::EDGetTokenT<std::vector<pat::Photon> > calibratedPhotonToken;
  const edm::ParameterSet  photonPurityID;

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
  photonNHisoToken         (consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonnhiso"))),   
  photonPurityID           (iConfig.getParameter<edm::ParameterSet>("photonPurityID")),
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
  }

  // value map for isolation
  produces<edm::ValueMap<float> > ("gammaiso");
  produces<edm::ValueMap<float> > ("chhadiso");
  produces<edm::ValueMap<float> > ("nhhadiso");
  produces<edm::ValueMap<float> > ("sigmaietaieta");

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
  Handle<edm::ValueMap<float> > photonNHisoH;
  iEvent.getByToken(photonNHisoToken, photonNHisoH);

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
  
  std::auto_ptr<edm::ValueMap<float> > outputgammaisomap(new ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > outputchhadisomap(new ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > outputnhhadisomap(new ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > outputsigietaietamap(new ValueMap<float>());

  std::auto_ptr<edm::ValueMap<float> > outputrndgammaiso04map(new ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > outputrndgammaiso08map(new ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > outputrndchhadiso04map(new ValueMap<float>());
  std::auto_ptr<edm::ValueMap<float> > outputrndchhadiso08map(new ValueMap<float>());
  // for purity studies
  std::auto_ptr<pat::PhotonRefVector>  outputphotonsPurity(new pat::PhotonRefVector);
  
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

      // in case no IP selections needs to be applied
      if(iele.getParameter<bool>("applyPVSelection") == false){
	pass_dxy = true;
	pass_dz  = true;
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
  std::vector<float> gammaiso;
  std::vector<float> chhadiso;
  std::vector<float> nhhadiso;
  std::vector<float> sigietaieta;
  std::vector<float> rndgammaiso04;
  std::vector<float> rndgammaiso08;
  std::vector<float> rndchhadiso04;
  std::vector<float> rndchhadiso08;
  
  // loop on the photon colection
  for (vector<pat::Photon>::const_iterator photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {

    // if a calibrate photon collection exists --> match them --> if no matching skip the candidate (calibrated are photns after scale/smearing corrections)
    vector<pat::Photon>::const_iterator calibpho_iter ;
    bool isValid = false;
    bool isMatched = false;

    if(calibratedPhotonsH.isValid()){ // check if it exists                                                                                                                                        
      isValid = true;
      calibpho_iter = calibratedPhotonsH->begin();
      if(calibratedPhotonsH->size() != photonsH->size())
	throw cms::Exception("PFCleaner") <<" different size between calibrated and un-calibrated photons --> check \n";
      for( ; calibpho_iter != calibratedPhotonsH->end(); ++calibpho_iter){
	// match by reference                                                                                                                                                                      
	if(calibpho_iter->superCluster()       == photons_iter->superCluster() and
	   calibpho_iter->originalObjectRef()  == photons_iter->originalObjectRef()){
	  isMatched = true;
	  break;
	}
      }
      if(not isMatched) continue;
    }
    
    // take the photon candidate
    pat::Photon photonCand;
    if(isValid)
      photonCand = *calibpho_iter;
    else
      photonCand = *photons_iter;

    // photon purity case
    if(addPhotonPurity){
      double rndphi04 = photonCand.phi() + M_PI/2.0;
      double rndphi08 = photonCand.phi() + M_PI/2.0;
      float  rndchisoval04 = 0.;
      float  rndchisoval08 = 0.;
      float  rndphisoval04 = 0.;
      float  rndphisoval08 = 0.;
      unsigned int counterR04 = 1;
      unsigned int counterR08 = 1;
      // select the random phi directions at the same eta of the photon candidate
      if(userandomphi) {
	rndphi04 = rand.Uniform(-M_PI, M_PI);
	while (randomConeOverlaps(rndphi04, photonCand.eta(), photonCand.phi(), *jetsH, 0.4)) {
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
	while (randomConeOverlaps(rndphi08, photonCand.eta(), photonCand.phi(), *jetsH, 0.8)){
	  rndphi08 = rand.Uniform(-M_PI, M_PI);
	  counterR08++;
	  if(counterR08 > 5000){
	    rndphi08 = -99;
	    continue;
	  }
	}
      }
    
      // if the random cone is found --> charged hadron isolation from pfCandidates with dR = 0.3 smaller than random cones
      if(rndphi04 != -99)
	rndchisoval04 = computeCHhadronIso(photonCand,pfcandsH,photonCand.eta(),rndphi04,0.3, verticesH->at(0));
      else
	rndchisoval04 = -99;
      
      if(rndphi08 != -99)
	rndchisoval08 = computeCHhadronIso(photonCand,pfcandsH,photonCand.eta(),rndphi08,0.3, verticesH->at(0));
      else 
	rndchisoval08 = -99;

      // compute the photon isolation in the random cone with dR = 0.3 smaller than random cones
      if(rndphi04 != -99)
	rndphisoval04 = computePhotonIso(photonCand,pfcandsH,photonCand.eta(), rndphi04, 0.3);      
      else 
	rndphisoval04 = -99;
      if(rndphi08 != -99)
	rndphisoval08 = computePhotonIso(photonCand,pfcandsH,photonCand.eta(), rndphi08, 0.3);      
      else
	rndphisoval08 = -99;
      
      rndgammaiso04.push_back(rndphisoval04);
      rndgammaiso08.push_back(rndphisoval08);
      rndchhadiso04.push_back(rndchisoval04);
      rndchhadiso08.push_back(rndchisoval08);
    }

    /// read the super-cluster and iso variables
    const Ptr<pat::Photon> photonPtr(photonsH, photons_iter - photonsH->begin());
    double photonsieie = (*photonsieieH)[photonPtr];
    double photonphiso = (*photonPHisoH)[photonPtr];
    double photonchiso = (*photonCHisoH)[photonPtr];
    double photonnhiso = (*photonNHisoH)[photonPtr];
    double photonhoe   =  photons_iter->hadTowOverEm();

    sigietaieta.push_back(photonsieie);
    // photon component from PF-candidate pdgId
    float gaisoval = computePhotonIso(photonCand,pfcandsH,photonCand.eta(),photonCand.phi(),0.3);     
    if(gaisoval != photonphiso)
      std::cout<<"PFCleaner::Photon isolation: difference between re-calculated and original one: "<<gaisoval<<" vs "<<photonphiso<<std::endl;
    gammaiso.push_back(gaisoval);
    
    // charged hadron isolation
    float chisoval = computeCHhadronIso(photonCand,pfcandsH,photonCand.eta(),photonCand.phi(),0.3, verticesH->at(0));
    if(chisoval != photonchiso)
      std::cout<<"PFCleaner::Photon Charged Hadron isolation: difference between re-calculated and original one: "<<chisoval<<" vs "<<photonchiso<<std::endl;
    chhadiso.push_back(chisoval);
    
    // charged hadron isolation
    float nhisoval = computeNHhadronIso(photonCand,pfcandsH,photonCand.eta(),photonCand.phi(),0.3);
    if(nhisoval != photonnhiso)
      std::cout<<"PFCleaner::Photon Neutral Hadron isolation: difference between re-calculated and original one: "<<nhisoval<<" vs "<<photonnhiso<<std::endl;
    nhhadiso.push_back(nhisoval);
    
    // Purity studies
    if(addPhotonPurity){
      if( photonCand.pt() > photonPurityID.getParameter<double>("ptMin") and
	  fabs(photonCand.eta()) < photonPurityID.getParameter<double>("etaMax") and
	  photonhoe < photonPurityID.getParameter<double>("HOverE") and
	  photonsieie <  photonPurityID.getParameter<double>("sigmaIetaIeta")  and
	  max(0.,photonchiso - rho*getChargedHadronEAForPhotonIso(photonCand.eta())) < photonPurityID.getParameter<double>("chIso") and
	  max(0.,photonnhiso - rho*getNeutralHadronEAForPhotonIso(photonCand.eta())) < photonPurityID.getParameter<double>("nhIso") ){

	if(isValid)
	  outputphotonsPurity->push_back(pat::PhotonRef(calibratedPhotonsH, calibpho_iter - calibratedPhotonsH->begin()));
	else
	  outputphotonsPurity->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));

      }
    }
    size_t ipos = 0;
    for(auto ipho : photonSelection){
      // check only photons with the following proprierties
      // standard ID
      bool passesid = (*photonsMapH.at(ipos))[photonPtr];
      if(passesid and isValid and isMatched){	
	if (fabs(calibpho_iter->superCluster()->eta()) > ipho.getParameter<double>("absEta") ) continue;
	if (calibpho_iter->pt() < ipho.getParameter<double>("ptMin")) continue;
	if (not calibpho_iter->passElectronVeto()) continue;
	outputphotons.at(ipos)->push_back(pat::PhotonRef(calibratedPhotonsH,calibpho_iter-calibratedPhotonsH->begin()));
      }
      else if(passesid and not isValid){
	if (fabs(photons_iter->superCluster()->eta()) > ipho.getParameter<double>("absEta") ) continue;
	if (photons_iter->pt() < ipho.getParameter<double>("ptMin")) continue;
	if (not photons_iter->passElectronVeto()) continue;
	outputphotons.at(ipos)->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
      }
      ipos++;
    }    
  }

  // Value maps
  if(addPhotonPurity){
    if(calibratedPhotonsH.isValid()){
      edm::ValueMap<float>::Filler rnd04gafiller(*outputrndgammaiso04map);
      rnd04gafiller.insert(calibratedPhotonsH, rndgammaiso04.begin(), rndgammaiso04.end());
      rnd04gafiller.fill();
      edm::ValueMap<float>::Filler rnd08gafiller(*outputrndgammaiso08map);
      rnd08gafiller.insert(calibratedPhotonsH, rndgammaiso08.begin(), rndgammaiso08.end());
      rnd08gafiller.fill();
      
      edm::ValueMap<float>::Filler rnd04chfiller(*outputrndchhadiso04map);
      rnd04chfiller.insert(calibratedPhotonsH, rndchhadiso04.begin(), rndchhadiso04.end());
      rnd04chfiller.fill();
      
      edm::ValueMap<float>::Filler rnd08chfiller(*outputrndchhadiso08map);
      rnd08chfiller.insert(calibratedPhotonsH, rndchhadiso08.begin(), rndchhadiso08.end());
      rnd08chfiller.fill();
    }
    else{
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
  }

  if(calibratedPhotonsH.isValid()){
    edm::ValueMap<float>::Filler gafiller(*outputgammaisomap);
    gafiller.insert(calibratedPhotonsH, gammaiso.begin(), gammaiso.end());
    gafiller.fill();
    
    edm::ValueMap<float>::Filler chfiller(*outputchhadisomap);
    chfiller.insert(calibratedPhotonsH, chhadiso.begin(), chhadiso.end());
    chfiller.fill();
    
    edm::ValueMap<float>::Filler nhfiller(*outputnhhadisomap);
    nhfiller.insert(calibratedPhotonsH, nhhadiso.begin(), nhhadiso.end());
    nhfiller.fill();

    edm::ValueMap<float>::Filler sigietaietafiller(*outputsigietaietamap);
    sigietaietafiller.insert(calibratedPhotonsH, sigietaieta.begin(), sigietaieta.end());
    sigietaietafiller.fill();

  }
  else{
    edm::ValueMap<float>::Filler gafiller(*outputgammaisomap);
    gafiller.insert(photonsH, gammaiso.begin(), gammaiso.end());
    gafiller.fill();
    
    edm::ValueMap<float>::Filler chfiller(*outputchhadisomap);
    chfiller.insert(photonsH, chhadiso.begin(), chhadiso.end());
    chfiller.fill();
    
    edm::ValueMap<float>::Filler nhfiller(*outputnhhadisomap);
    nhfiller.insert(photonsH, nhhadiso.begin(), nhhadiso.end());
    nhfiller.fill();

    edm::ValueMap<float>::Filler sigietaietafiller(*outputsigietaietamap);
    sigietaietafiller.insert(photonsH, sigietaieta.begin(), sigietaieta.end());
    sigietaietafiller.fill();

  }


  iEvent.put(outputgammaisomap,      "gammaiso");  
  iEvent.put(outputchhadisomap,      "chhadiso");  
  iEvent.put(outputnhhadisomap,      "nhhadiso");  
  iEvent.put(outputsigietaietamap,   "sigmaietaieta");  

  if(addPhotonPurity){
    iEvent.put(outputrndgammaiso04map, "rndgammaiso04");
    iEvent.put(outputrndgammaiso08map, "rndgammaiso08");
    iEvent.put(outputrndchhadiso04map, "rndchhadiso04");
    iEvent.put(outputrndchhadiso08map, "rndchhadiso08");
    iEvent.put(outputphotonsPurity,    "photonsPurity");
  }


  ///////////
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
bool PFCleaner::randomConeOverlaps(const double & randomphi, const double & photoneta, const double & photonphi, const std::vector<pat::Jet> & jets, const double & vetoCone) {
  if (reco::deltaR(photoneta, randomphi, photoneta, photonphi) < vetoCone) return true;
  for (std::size_t i = 0; i < jets.size(); i++) {    
    if (reco::deltaR(photoneta,photonphi,jets[i].eta(),jets[i].phi()) < 0.4) continue;
    if (jets[i].pt() > 30. && reco::deltaR(photoneta, randomphi, jets[i].eta(), jets[i].phi()) < vetoCone){
      return true;
    }
  }
  return false;
}


/// to calculate photon isolation
double PFCleaner::computePhotonIso(const pat::Photon & photon, edm::Handle<edm::View<reco::Candidate> >  pfcandsH, const double & eta, const double & phi, const double & isoCone){

  double isoval=0;
  for(size_t i = 0; i < pfcandsH->size(); i++) {
    const auto& pfcand = pfcandsH->ptrAt(i);
    if (pfcand->pdgId()  ==  22 && // pf photon candidate
	deltaR(eta, phi, pfcand->eta(), pfcand->phi()) <= isoCone){
      // check footprint
      bool skipCandidate = isInFootprint(photon.associatedPackedPFCandidates(),pfcand);
      if(skipCandidate) continue;
      else
	isoval += pfcand->pt();
    }
  }
  return isoval;
}


// to compute charged hadronic isolation
double PFCleaner::computeCHhadronIso(const pat::Photon & photon, edm::Handle<edm::View<reco::Candidate> >  pfcandsH, const double & eta, const double & phi, const double & isoCone, const reco::Vertex & pv){

  double isoval=0;
  for(size_t i = 0; i < pfcandsH->size(); i++) {
    const auto& pfcand = pfcandsH->ptrAt(i);
    if (abs(pfcand->pdgId()) == 211 && deltaR(eta,phi,pfcand->eta(),pfcand->phi()) <= isoCone){
      // check footprint
      bool skipCandidate = isInFootprint(photon.associatedPackedPFCandidates(),pfcand);
      if(skipCandidate) continue;
      
      // only charged hadrons from PV
      const reco::Track *theTrack = &( ((const edm::Ptr<pat::PackedCandidate>) pfcand)->pseudoTrack());      
      float dxy = theTrack->dxy(pv.position());
      if(fabs(dxy) > 0.1) continue;
      
      float dz  = theTrack->dz(pv.position());
      if (fabs(dz) > 0.2) continue;
      
      isoval += pfcand->pt();
    }
  }
  return isoval;
}


// to compute charged hadronic isolation
double PFCleaner::computeNHhadronIso(const pat::Photon & photon, edm::Handle<edm::View<reco::Candidate> >  pfcandsH, const double & eta, const double & phi, const double & isoCone){

  double isoval=0;
  for(size_t i = 0; i < pfcandsH->size(); i++) {
    const auto& pfcand = pfcandsH->ptrAt(i);
    if (abs(pfcand->pdgId()) == 130 && deltaR(eta,phi, pfcand->eta(), pfcand->phi()) <= isoCone){
      // check footprint
      bool skipCandidate = isInFootprint(photon.associatedPackedPFCandidates(),pfcand);
      if(skipCandidate) continue;
      else
	isoval += pfcand->pt();
    }
  }
  return isoval;
}


double PFCleaner::getChargedHadronEAForPhotonIso(const double & eta) {
  if (fabs(eta) < 1.0) return 0.0360;
  else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return  0.0377;
  else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return  0.0306;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return  0.0283;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return  0.0254;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.4  ) return  0.0217;
  else if (fabs(eta) >= 2.4) return 0.0167 ;
  else return 0.;
}

double PFCleaner::getNeutralHadronEAForPhotonIso(const double & eta) {
  if (fabs(eta) < 1.0) return 0.0597;
  else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.0807;
  else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.0629;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.0197;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.0184;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.4  ) return 0.0284;
  else if (fabs(eta) >= 2.4) return 0.0591;
  else return 0.;
}

double PFCleaner::getGammaEAForPhotonIso(const double & eta) {
  if (fabs(eta) < 1.0) return 0.1210;
  else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.1107;
  else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.0699;
  else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.1056;
  else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.1457;
  else if (fabs(eta) >= 2.3   && fabs(eta) < 2.3  ) return 0.1719;
  else if (fabs(eta) >= 2.4) return 0.1998;
  else return 0.;
}



DEFINE_FWK_MODULE(PFCleaner);
