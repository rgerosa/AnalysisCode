#include <memory>
#include <vector>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
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

        double getChargedHadronEAForPhotonIso(double);
        double getNeutralHadronEAForPhotonIso(double);
        double getGammaEAForPhotonIso(double);
        bool testPhotonIsolation(const reco::Photon&, double, double, double, double);

        edm::InputTag vertices;
        edm::InputTag rhoTag;
        edm::InputTag muons;
        edm::InputTag electrons;
        edm::InputTag photons;
        edm::InputTag electronVetoIdMap;
        edm::InputTag electronMediumIdMap;
        edm::InputTag photonSIEIEMap;
        edm::InputTag photonChargedIsoMap;
        edm::InputTag photonNeutralIsoMap;
        edm::InputTag photonGammaIsoMap;
};

PFCleaner::PFCleaner(const edm::ParameterSet& iConfig): 
    vertices(iConfig.getParameter<edm::InputTag>("vertices")),
    rhoTag(iConfig.getParameter<edm::InputTag>("rho")),
    muons(iConfig.getParameter<edm::InputTag>("muons")),
    electrons(iConfig.getParameter<edm::InputTag>("electrons")),
    photons(iConfig.getParameter<edm::InputTag>("photons")),
    electronVetoIdMap(iConfig.getParameter<edm::InputTag>("electronidveto")),
    electronMediumIdMap(iConfig.getParameter<edm::InputTag>("electronidmedium")),
    photonSIEIEMap(iConfig.getParameter<edm::InputTag>("photonsigmaietaieta")),
    photonChargedIsoMap(iConfig.getParameter<edm::InputTag>("photonchargediso")),
    photonNeutralIsoMap(iConfig.getParameter<edm::InputTag>("photonneutraliso")),
    photonGammaIsoMap(iConfig.getParameter<edm::InputTag>("photongammaiso"))
{
    produces<pat::MuonRefVector>("muons");
    produces<pat::ElectronRefVector>("electrons");
    produces<pat::PhotonRefVector>("photons");
    produces<pat::MuonRefVector>("tightmuons");
    produces<pat::ElectronRefVector>("tightelectrons");
    produces<pat::PhotonRefVector>("tightphotons");
}


PFCleaner::~PFCleaner() {
}

void PFCleaner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    Handle<vector<Vertex> > verticesH;
    iEvent.getByLabel(vertices, verticesH);

    Handle<vector<pat::Muon> > muonsH;
    iEvent.getByLabel(muons, muonsH);

    Handle<vector<pat::Electron> > electronsH;
    iEvent.getByLabel(electrons, electronsH);

    Handle<vector<pat::Photon> > photonsH;
    iEvent.getByLabel(photons, photonsH);

    Handle<ValueMap<bool> > electronVetoIdH;
    iEvent.getByLabel(electronVetoIdMap, electronVetoIdH);

    Handle<ValueMap<bool> > electronMediumIdH;
    iEvent.getByLabel(electronMediumIdMap, electronMediumIdH);

    Handle<ValueMap<float> > photonSIEIEH;
    iEvent.getByLabel(photonSIEIEMap, photonSIEIEH);

    Handle<ValueMap<float> > photonChargedIsoH;
    iEvent.getByLabel(photonChargedIsoMap, photonChargedIsoH);

    Handle<ValueMap<float> > photonNeutralIsoH;
    iEvent.getByLabel(photonNeutralIsoMap, photonNeutralIsoH);

    Handle<ValueMap<float> > photonGammaIsoH;
    iEvent.getByLabel(photonGammaIsoMap, photonGammaIsoH);

    Handle<double> rhoH;
    iEvent.getByLabel(rhoTag, rhoH);
    double rho = *rhoH;

    std::auto_ptr<pat::MuonRefVector> outputmuons(new pat::MuonRefVector);
    std::auto_ptr<pat::ElectronRefVector> outputelectrons(new pat::ElectronRefVector);
    std::auto_ptr<pat::PhotonRefVector> outputphotons(new pat::PhotonRefVector);
    std::auto_ptr<pat::MuonRefVector> outputtightmuons(new pat::MuonRefVector);
    std::auto_ptr<pat::ElectronRefVector> outputtightelectrons(new pat::ElectronRefVector);
    std::auto_ptr<pat::PhotonRefVector> outputtightphotons(new pat::PhotonRefVector);

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

    for (vector<pat::Photon>::const_iterator photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
        const Ptr<pat::Photon> photonPtr(photonsH, photons_iter - photonsH->begin());
        bool passeskincuts = (photons_iter->pt() > 15 && fabs(photons_iter->superCluster()->eta()) < 2.5);
        bool haspixelseed = photons_iter->hasPixelSeed();
        double hovere = photons_iter->hadTowOverEm();
        double sieie = (*photonSIEIEH)[photonPtr];
        double chiso = (*photonChargedIsoH)[photonPtr];
        double nhiso = (*photonNeutralIsoH)[photonPtr];
        double phiso = (*photonGammaIsoH)[photonPtr];

        bool passesiso = testPhotonIsolation(*photons_iter, chiso, nhiso, phiso, rho);
        bool passesselection = false;

        if (passesiso && !haspixelseed && photons_iter->isEB() && sieie < 0.0101 && hovere < 0.0320) passesselection = true;
        if (passesiso && !haspixelseed && photons_iter->isEE() && sieie < 0.0264 && hovere < 0.0166) passesselection = true;

        if (passeskincuts && passesselection) {
            outputphotons->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
            if (photons_iter->pt() > 160) outputtightphotons->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
        }
    }

    iEvent.put(outputmuons, "muons");
    iEvent.put(outputelectrons, "electrons");
    iEvent.put(outputphotons, "photons");
    iEvent.put(outputtightmuons, "tightmuons");
    iEvent.put(outputtightelectrons, "tightelectrons");
    iEvent.put(outputtightphotons, "tightphotons");
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

double PFCleaner::getChargedHadronEAForPhotonIso(double eta) {
    if      (fabs(eta) < 1.0)                         return 0.0130;
    else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.0096;
    else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.0107;
    else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.0077;
    else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.0088;
    else if (fabs(eta) >= 2.3   && fabs(eta) < 2.4  ) return 0.0065;
    else if (fabs(eta) >= 2.4)                        return 0.0030;
    else return 0.;
}

double PFCleaner::getNeutralHadronEAForPhotonIso(double eta) {
    if      (fabs(eta) < 1.0)                         return 0.0056;
    else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.0107;
    else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.0019;
    else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.0011;
    else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.0077;
    else if (fabs(eta) >= 2.3   && fabs(eta) < 2.4  ) return 0.0178;
    else if (fabs(eta) >= 2.4)                        return 0.1675;
    else return 0.;
}

double PFCleaner::getGammaEAForPhotonIso(double eta) {
    if      (fabs(eta) < 1.0)                         return 0.0896;
    else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.0762;
    else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.0383;
    else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.0534;
    else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.0846;
    else if (fabs(eta) >= 2.3   && fabs(eta) < 2.4  ) return 0.1032;
    else if (fabs(eta) >= 2.4)                        return 0.1598;
    else return 0.;
}

bool PFCleaner::testPhotonIsolation(const reco::Photon& photon, double chargedHadronIsolation, double neutralHadronIsolation, double gammaIsolation, double rhoval) {
    double corrCHIso = chargedHadronIsolation - rhoval * getChargedHadronEAForPhotonIso(photon.eta());
    double corrNHIso = neutralHadronIsolation - rhoval * getNeutralHadronEAForPhotonIso(photon.eta());
    double corrPHIso = gammaIsolation - rhoval * getGammaEAForPhotonIso(photon.eta());

    if (corrCHIso < 0.) corrCHIso = 0.;
    if (corrNHIso < 0.) corrNHIso = 0.;
    if (corrPHIso < 0.) corrPHIso = 0.;

    if (photon.isEB()) {
        if (corrCHIso < 1.90 && corrNHIso < 2.96 + 0.0025*photon.pt() && corrPHIso < 1.39 + 0.0010*photon.pt()) return true;
        else return false;
    }
    else if (photon.isEE()) {
        if (corrCHIso < 1.95 && corrNHIso < 4.42 + 0.0118*photon.pt() && corrPHIso < 1.89 + 0.0059*photon.pt()) return true;
        else return false;
    }
    else return false;
}

DEFINE_FWK_MODULE(PFCleaner);
