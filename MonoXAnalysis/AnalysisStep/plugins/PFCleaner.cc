#include <memory>
#include <vector>
#include <iostream>

#include <TRandom3.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
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

class PFCleaner : public edm::EDProducer {
    public:
        explicit PFCleaner(const edm::ParameterSet&);
        ~PFCleaner();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    private:
        virtual void beginJob() override;
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        
        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

        bool randomConeOverlaps(double, double, double, std::vector<pat::Jet>);
  //livia 
  bool IsPassingPhotonIsoNew(double,double,double, double,double,double,double);
        edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken;
        edm::EDGetTokenT<edm::View<reco::Candidate> > pfcandsToken;
        edm::EDGetTokenT<std::vector<pat::Jet> > jetsToken;
        edm::EDGetTokenT<std::vector<pat::Muon> > muonsToken;
        edm::EDGetTokenT<std::vector<pat::Electron> > electronsToken;
        edm::EDGetTokenT<std::vector<pat::Photon> > photonsToken;
        edm::EDGetTokenT<edm::ValueMap<bool> > electronVetoIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool> > electronMediumIdMapToken;
        edm::EDGetTokenT<edm::ValueMap<bool> > photonLooseIdMapToken;
  //livia 
  edm::InputTag rhoTag;
   
  edm::EDGetTokenT<edm::ValueMap<float> > photonsieieToken;
  edm::EDGetTokenT<edm::ValueMap<float> > photonPHisoToken;
  edm::EDGetTokenT<edm::ValueMap<float> > photonCHisoToken;
 
  
        bool userandomphi;
};

PFCleaner::PFCleaner(const edm::ParameterSet& iConfig): 

    verticesToken            (consumes<std::vector<reco::Vertex> > (iConfig.getParameter<edm::InputTag>("vertices"))),
    pfcandsToken             (consumes<edm::View<reco::Candidate> > (iConfig.getParameter<edm::InputTag>("pfcands"))),
    jetsToken                (consumes<std::vector<pat::Jet> > (iConfig.getParameter<edm::InputTag>("jets"))),
    muonsToken               (consumes<std::vector<pat::Muon> > (iConfig.getParameter<edm::InputTag>("muons"))), 
    electronsToken           (consumes<std::vector<pat::Electron> > (iConfig.getParameter<edm::InputTag>("electrons"))),
    photonsToken             (consumes<std::vector<pat::Photon> > (iConfig.getParameter<edm::InputTag>("photons"))),
    electronVetoIdMapToken   (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronidveto"))),
    electronMediumIdMapToken (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electronidmedium"))),
    photonLooseIdMapToken    (consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("photonidloose"))),
    //livia
    rhoTag(iConfig.getParameter<edm::InputTag>("rho")),
    photonsieieToken(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonsieie"))), 
    photonPHisoToken(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonPHiso"))), 
    photonCHisoToken(consumes<edm::ValueMap<float> > (iConfig.getParameter<edm::InputTag>("photonCHiso"))),   
    userandomphi             (iConfig.existsAs<bool>("userandomphiforRC") ? iConfig.getParameter<bool>("userandomphiforRC") : false)
{
    produces<pat::MuonRefVector>("muons");
    produces<pat::ElectronRefVector>("electrons");
    produces<pat::PhotonRefVector>("photons");
    produces<pat::MuonRefVector>("tightmuons");
    produces<pat::ElectronRefVector>("tightelectrons");
    produces<pat::PhotonRefVector>("tightphotons");
    produces<pat::PhotonRefVector>("loosephotons");
    produces<edm::ValueMap<float> >("rndgammaiso");
    produces<edm::ValueMap<float> >("rndchhadiso");
    produces<edm::ValueMap<bool> >("photonidNew");
    produces<pat::PhotonRefVector>("photonsNew");
    produces<pat::PhotonRefVector>("tightphotonsNew");
}


PFCleaner::~PFCleaner() {
}

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

    Handle<std::vector<pat::Muon> > muonsH;
    iEvent.getByToken(muonsToken, muonsH);

    Handle<std::vector<pat::Electron> > electronsH;
    iEvent.getByToken(electronsToken, electronsH);

    Handle<std::vector<pat::Photon> > photonsH;
    iEvent.getByToken(photonsToken, photonsH);

    Handle<edm::ValueMap<bool> > electronVetoIdH;
    iEvent.getByToken(electronVetoIdMapToken, electronVetoIdH);

    Handle<edm::ValueMap<bool> > electronMediumIdH;
    iEvent.getByToken(electronMediumIdMapToken, electronMediumIdH);

    Handle<edm::ValueMap<bool> > photonLooseIdH;
    iEvent.getByToken(photonLooseIdMapToken, photonLooseIdH);
    
    Handle<edm::ValueMap<float> > photonsieieH;
    iEvent.getByToken(photonsieieToken, photonsieieH);

    Handle<edm::ValueMap<float> > photonPHisoH;
    iEvent.getByToken(photonPHisoToken, photonPHisoH);

    Handle<edm::ValueMap<float> > photonCHisoH;
    iEvent.getByToken(photonCHisoToken, photonCHisoH);

    
    Handle<double> rhoH;
    iEvent.getByLabel(rhoTag, rhoH);
    double rho = *rhoH;
    
    

    std::auto_ptr<pat::MuonRefVector> outputmuons(new pat::MuonRefVector);
    std::auto_ptr<pat::ElectronRefVector> outputelectrons(new pat::ElectronRefVector);
    std::auto_ptr<pat::PhotonRefVector> outputphotons(new pat::PhotonRefVector);
    std::auto_ptr<pat::MuonRefVector> outputtightmuons(new pat::MuonRefVector);
    std::auto_ptr<pat::ElectronRefVector> outputtightelectrons(new pat::ElectronRefVector);
    std::auto_ptr<pat::PhotonRefVector> outputtightphotons(new pat::PhotonRefVector);
    std::auto_ptr<pat::PhotonRefVector> outputloosephotons(new pat::PhotonRefVector);
    std::auto_ptr<edm::ValueMap<float> > outputgammaisomap(new ValueMap<float>());
    std::auto_ptr<edm::ValueMap<float> > outputchhadisomap(new ValueMap<float>());
    //livia
    std::auto_ptr<edm::ValueMap<bool> > outputphotonisoNewmap(new ValueMap<bool>());
    std::auto_ptr<pat::PhotonRefVector> outputphotonsNew(new pat::PhotonRefVector);
    std::auto_ptr<pat::PhotonRefVector> outputtightphotonsNew(new pat::PhotonRefVector);
 

    TRandom3 rand;

    for (vector<pat::Muon>::const_iterator muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
        if (verticesH->size() == 0) continue;
        bool passeskincuts = (muons_iter->pt() > 10 && fabs(muons_iter->eta()) < 2.4);
        float isoval = muons_iter->pfIsolationR04().sumNeutralHadronEt;
        isoval += muons_iter->pfIsolationR04().sumPhotonEt;
        isoval -= 0.5*muons_iter->pfIsolationR04().sumPUPt;
        if (isoval < 0.) isoval = 0.;
        isoval += muons_iter->pfIsolationR04().sumChargedHadronPt;
        isoval /= muons_iter->pt();

        if (passeskincuts) {
            if (muon::isLooseMuon(*muons_iter) && isoval < 0.2) outputmuons->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
            if (muon::isTightMuon(*muons_iter, *(verticesH->begin())) && isoval < 0.12) outputtightmuons->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
        }
    }

    for (vector<pat::Electron>::const_iterator electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) {
        const Ptr<pat::Electron> electronPtr(electronsH, electrons_iter - electronsH->begin());
        bool passeskincuts = (electrons_iter->pt() > 10 && fabs(electrons_iter->superCluster()->eta()) < 2.5);
        bool passesvetoid = (*electronVetoIdH)[electronPtr];
        bool passesmediumid = (*electronMediumIdH)[electronPtr];

        if (passeskincuts && passesvetoid) outputelectrons->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
        if (passeskincuts && passesmediumid) outputtightelectrons->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
    }

    std::vector<float> rndgammaiso;
    std::vector<float> rndchhadiso;
    std::vector<bool> photonidNew;//livia

    for (vector<pat::Photon>::const_iterator photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
        float gaisoval = 0.;
        float chisoval = 0.;
	bool 	passesphotonisoNew = false;
        double rndphi = photons_iter->phi() + M_PI/2.0;
        if (userandomphi) {
            rndphi = rand.Uniform(-M_PI, M_PI);
            while (randomConeOverlaps(rndphi, photons_iter->eta(), photons_iter->phi(), *jetsH)) rndphi = rand.Uniform(-M_PI, M_PI);
        }

        for(size_t i = 0; i < pfcandsH->size(); i++) {
            const auto& pfcand = pfcandsH->ptrAt(i);
            if (    pfcand->pdgId()  ==  22 && deltaR(photons_iter->eta(), rndphi, pfcand->eta(), pfcand->phi()) <= 0.3) gaisoval += pfcand->pt();
            if (abs(pfcand->pdgId()) == 211 && deltaR(photons_iter->eta(), rndphi, pfcand->eta(), pfcand->phi()) <= 0.3) chisoval += pfcand->pt();
        }
        rndgammaiso.push_back(gaisoval);
        rndchhadiso.push_back(chisoval);


	//livia
	const Ptr<pat::Photon> photonPtr(photonsH, photons_iter - photonsH->begin());
	double photonsieie=(*photonsieieH)[photonPtr];
	double photonphiso=(*photonPHisoH)[photonPtr];
	double photonchiso=(*photonCHisoH)[photonPtr];
	double photonhoe =  photons_iter->hadTowOverEm();
	passesphotonisoNew = IsPassingPhotonIsoNew(photons_iter->pt(),photons_iter->superCluster()->eta(),photonhoe,photonchiso,photonphiso,photonsieie,rho );
	photonidNew.push_back(passesphotonisoNew);


        if (fabs(photons_iter->superCluster()->eta()) > 2.5 || photons_iter->pt() < 15) continue;
        if (photons_iter->r9() > 0.8 || photons_iter->chargedHadronIso() < 20. || photons_iter->chargedHadronIso() < photons_iter->pt()*0.3) outputloosephotons->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin())); 

        
        bool passeslooseid = (*photonLooseIdH)[photonPtr];
        if (passeslooseid && photons_iter->passElectronVeto()) {
            outputphotons->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
            if (photons_iter->pt() > 175) outputtightphotons->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
        }

	//livia
	if (passesphotonisoNew && photons_iter->passElectronVeto()){
	  std::cout<<"pass"<<std::endl;
	  outputphotonsNew->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
            if (photons_iter->pt() > 175) outputtightphotonsNew->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
	}
    }


   

    edm::ValueMap<float>::Filler gafiller(*outputgammaisomap);
    gafiller.insert(photonsH, rndgammaiso.begin(), rndgammaiso.end());
    gafiller.fill();

    edm::ValueMap<float>::Filler chfiller(*outputchhadisomap);
    chfiller.insert(photonsH, rndchhadiso.begin(), rndchhadiso.end());
    chfiller.fill();

    //livia
    
    edm::ValueMap<bool>::Filler photonisoNewfiller(*outputphotonisoNewmap);
    photonisoNewfiller.insert(photonsH, photonidNew.begin(), photonidNew.end());
    photonisoNewfiller.fill();
    
    iEvent.put(outputmuons, "muons");
    iEvent.put(outputelectrons, "electrons");
    iEvent.put(outputphotons, "photons");
    iEvent.put(outputtightmuons, "tightmuons");
    iEvent.put(outputtightelectrons, "tightelectrons");
    iEvent.put(outputtightphotons, "tightphotons");
    iEvent.put(outputloosephotons, "loosephotons");
    iEvent.put(outputgammaisomap, "rndgammaiso");
    iEvent.put(outputchhadisomap, "rndchhadiso");
    //  iEvent.put(outputphotonisoNewmap, "photonidNew");

    //livia
    iEvent.put(outputphotonsNew, "photonsNew");
    iEvent.put(outputtightphotonsNew, "tightphotonsNew");
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

bool PFCleaner::randomConeOverlaps(double randomphi, double photoneta, double photonphi, std::vector<pat::Jet> jets) {
    if (reco::deltaR(photoneta, randomphi, photoneta, photonphi) < 0.4) return true;
    for (std::size_t i = 0; i < jets.size(); i++) {
        if (jets[i].pt() > 30. && reco::deltaR(photoneta, randomphi, jets[i].eta(), jets[i].phi()) < 0.4) return true;
    }
    return false;
}

bool PFCleaner::IsPassingPhotonIsoNew(double pt,double eta,double hoe, double chiso,double phiso,double sieie,double rho ){
  double chisoCUT;
  double phisoCUT;
  double sieieCUT;
  double hoeCUT;
  double alpha;
  double k;
  double EA=0;
  double newphiso;
  int isPassing = false;
  if(abs(eta)<1.4442){
    chisoCUT= 5;
    sieieCUT=0.0105;
    hoeCUT=0.05;
    phisoCUT=2.75;
    alpha=2.5;
    k=0.0045;
    if(abs(eta)<0.9)EA=0.17;
    if(abs(eta)>0.9&&abs(eta)<1.4442)EA=0.14;
  }
  newphiso = alpha+phiso-rho*EA-k*pt;
  if(newphiso<phisoCUT&&chiso<chisoCUT&&sieie<sieieCUT&&hoe<hoeCUT)isPassing =true;

  return isPassing;
}

DEFINE_FWK_MODULE(PFCleaner);
