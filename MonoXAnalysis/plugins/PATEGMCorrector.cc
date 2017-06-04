#include <memory>
#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"


template< class T>
class PATEGMCorrectorT : public edm::stream::EDProducer<> {

public:
  explicit PATEGMCorrectorT(const edm::ParameterSet&);
  ~PATEGMCorrectorT();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void beginJob();
  virtual void endJob();

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  const edm::InputTag    objectTag;
  const std::vector<edm::ParameterSet> correctionBin;
  const bool isMC;
  const edm::InputTag rechitEBTag;
  const edm::InputTag rechitEETag;

  edm::EDGetTokenT<std::vector<T> > objectToken;
  edm::EDGetTokenT<EcalRecHitCollection> recHitEBToken;
  edm::EDGetTokenT<EcalRecHitCollection> recHitEEToken;
};

template< class T>
PATEGMCorrectorT<T>::PATEGMCorrectorT(const edm::ParameterSet& iConfig):
  objectTag(iConfig.getParameter<edm::InputTag>("src")),
  correctionBin(iConfig.getParameter<std::vector<edm::ParameterSet> >("correction")),
  isMC(iConfig.getParameter<bool>("isMC")),
  rechitEBTag(iConfig.getParameter<edm::InputTag>("recHitEB")),
  rechitEETag(iConfig.getParameter<edm::InputTag>("recHitEE")){

  recHitEBToken = consumes<EcalRecHitCollection>(rechitEBTag);
  recHitEEToken = consumes<EcalRecHitCollection>(rechitEETag);
  objectToken   = consumes<std::vector<T> >(objectTag);

  produces<std::vector<T> >();
}

template< class T>
PATEGMCorrectorT<T>::~PATEGMCorrectorT() {}


template<>
void PATEGMCorrectorT<pat::Electron>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace reco;
  using namespace std;

  Handle<std::vector<pat::Electron> > egammaCandidatesH;
  iEvent.getByToken(objectToken,egammaCandidatesH);

  Handle<EcalRecHitCollection> barrelRecHitH;
  iEvent.getByToken(recHitEBToken,barrelRecHitH);

  Handle<EcalRecHitCollection> endcapRecHitH;
  iEvent.getByToken(recHitEEToken,endcapRecHitH);

  std::unique_ptr<std::vector<pat::Electron>> egammaoutput (new std::vector<pat::Electron>());
  for (auto egamma_iter : *egammaCandidatesH){
    if(isMC){// just make a close
      egammaoutput->push_back(egamma_iter);
    }
    else{
      DetId detid = egamma_iter.superCluster()->seed()->seed();
      const EcalRecHit * rh = NULL;
      float Ecorr = 1;
      if (detid.subdetId() == EcalBarrel) {
	auto rh_i =  barrelRecHitH->find(detid);
	if(rh_i !=  barrelRecHitH->end())
	  rh = &(*rh_i);
	else
	  rh = NULL;
      } else {
	rh = NULL;
      }

      if(rh!=NULL){
	for(auto pset : correctionBin){
	  if(rh->energy() > pset.getParameter<double>("eMin") and
	     rh->energy() < pset.getParameter<double>("eMax") )
	    Ecorr = pset.getParameter<double>("value");
	}
      }

      float newEcalEnergy = egamma_iter.correctedEcalEnergy()*Ecorr;
      float newEcalEnergyError = egamma_iter.correctedEcalEnergyError()*Ecorr;
      // clone the electron as it was --> no need for Eecal and P combination
      if(Ecorr == 1)
	egammaoutput->push_back(egamma_iter);
      else{
	math::XYZTLorentzVector oldMomentum = egamma_iter.p4();
	math::XYZTLorentzVector newMomentum = math::XYZTLorentzVector(oldMomentum.x()*newEcalEnergy/egamma_iter.correctedEcalEnergy(),
								      oldMomentum.y()*newEcalEnergy/egamma_iter.correctedEcalEnergy(),
								    oldMomentum.z()*newEcalEnergy/egamma_iter.correctedEcalEnergy(),
								      newEcalEnergy);
	egamma_iter.setCorrectedEcalEnergy(newEcalEnergy);
	egamma_iter.setCorrectedEcalEnergyError(newEcalEnergyError);
	egamma_iter.correctMomentum(newMomentum,egamma_iter.trackMomentumError(),egamma_iter.p4Error(reco::GsfElectron::P4_COMBINATION));
	egammaoutput->push_back(egamma_iter);
      }
    }
  }
  iEvent.put(std::move(egammaoutput));  
}

template<>
void PATEGMCorrectorT<pat::Photon>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace reco;
  using namespace std;

  Handle<std::vector<pat::Photon> > egammaCandidatesH;
  iEvent.getByToken(objectToken,egammaCandidatesH);

  Handle<EcalRecHitCollection> barrelRecHitH;
  iEvent.getByToken(recHitEBToken,barrelRecHitH);

  Handle<EcalRecHitCollection> endcapRecHitH;
  iEvent.getByToken(recHitEEToken,endcapRecHitH);

  std::unique_ptr<std::vector<pat::Photon>> egammaoutput (new std::vector<pat::Photon>());

  for (auto egamma_iter : *egammaCandidatesH){
    if(isMC)
      egammaoutput->push_back(egamma_iter);
    else{
      DetId detid = egamma_iter.superCluster()->seed()->seed();
      const EcalRecHit * rh = NULL;
      float Ecorr = 1;
      if (detid.subdetId() == EcalBarrel) {
	auto rh_i =  barrelRecHitH->find(detid);
	if(rh_i !=  barrelRecHitH->end())
	rh = &(*rh_i);
	else
	  rh = NULL;
      } else {
	rh = NULL;
      }

      if(rh!=NULL){
	for(auto pset : correctionBin){
	  if(rh->energy() > pset.getParameter<double>("eMin") and
	     rh->energy() < pset.getParameter<double>("eMax") )
	    Ecorr = pset.getParameter<double>("value");
	}
      }

      float newEcalEnergy = egamma_iter.getCorrectedEnergy(reco::Photon::P4type::regression2)*Ecorr;
      float newEcalEnergyError = egamma_iter.getCorrectedEnergyError(reco::Photon::P4type::regression2)*Ecorr;
      // clone the electron as it was --> no need for Eecal and P combination
      if(Ecorr == 1)
	egammaoutput->push_back(egamma_iter);
      else{
	egamma_iter.setCorrectedEnergy(reco::Photon::P4type::regression2,newEcalEnergy,newEcalEnergyError, true);
	egammaoutput->push_back(egamma_iter);
      }
    }
  }
  iEvent.put(std::move(egammaoutput));

}


template< class T>
void PATEGMCorrectorT<T>::beginJob() {
}

template< class T>
void PATEGMCorrectorT<T>::endJob() {
}

template< class T>
void PATEGMCorrectorT<T>::beginRun(edm::Run const&, edm::EventSetup const&) {
}

template< class T>
void PATEGMCorrectorT<T>::endRun(edm::Run const&, edm::EventSetup const&) {
}

template< class T>
void PATEGMCorrectorT<T>::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

template< class T>
void PATEGMCorrectorT<T>::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

template< class T>
void PATEGMCorrectorT<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

typedef PATEGMCorrectorT<pat::Electron> PATElectronCorrector ;
DEFINE_FWK_MODULE(PATElectronCorrector);

typedef PATEGMCorrectorT<pat::Photon> PATPhotonCorrector;
DEFINE_FWK_MODULE(PATPhotonCorrector);
