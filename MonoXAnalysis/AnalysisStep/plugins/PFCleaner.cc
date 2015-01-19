/*

Notes:

Muon ID and isolation : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
Photon ID and isolation : https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012
Electron ID and isolation : https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification

Need to switch to Run-II electron ID : https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2

*/


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
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "MonoXAnalysis/AnalysisStep/interface/EGammaCutBasedEleId.h"
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

        edm::InputTag beamspot;
        edm::InputTag vertices;
        edm::InputTag rhoTag;
        edm::InputTag muons;
        edm::InputTag electrons;
        edm::InputTag conversions;
        edm::InputTag photons;

        double d0cut, dzcut;

        bool debug;
};

PFCleaner::PFCleaner(const edm::ParameterSet& iConfig): 
    beamspot(iConfig.getParameter<edm::InputTag>("beamspot")),
    vertices(iConfig.getParameter<edm::InputTag>("vertices")),
    rhoTag(iConfig.getParameter<edm::InputTag>("rho")),
    muons(iConfig.getParameter<edm::InputTag>("muons")),
    electrons(iConfig.getParameter<edm::InputTag>("electrons")),
    conversions(iConfig.getParameter<edm::InputTag>("conversions")),
    photons(iConfig.getParameter<edm::InputTag>("photons")),
    d0cut(iConfig.getParameter<double>("d0cut")),
    dzcut(iConfig.getParameter<double>("dzcut")),
    debug(iConfig.existsAs<bool>("debug") ? iConfig.getParameter<bool>("debug") : false)
{
    produces<pat::MuonRefVector>("muons");
    produces<pat::ElectronRefVector>("electrons");
    produces<pat::ElectronRefVector>("electronsnew");
    produces<pat::MuonRefVector>("tightmuons");
    produces<pat::ElectronRefVector>("tightelectrons");
    produces<pat::PhotonRefVector>("photons");
}


PFCleaner::~PFCleaner() {
}

void PFCleaner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    Handle<BeamSpot> beamspotH;
    iEvent.getByLabel(beamspot, beamspotH);

    Handle<vector<Vertex> > verticesH;
    iEvent.getByLabel(vertices, verticesH);

    Handle<vector<pat::Muon> > muonsH;
    iEvent.getByLabel(muons, muonsH);

    Handle<vector<pat::Electron> > electronsH;
    iEvent.getByLabel(electrons, electronsH);

    Handle<vector<Conversion> > conversionsH;
    iEvent.getByLabel(conversions, conversionsH);

    Handle<vector<pat::Photon> > photonsH;
    iEvent.getByLabel(photons, photonsH);

    Handle<double> rhoH;
    iEvent.getByLabel(rhoTag, rhoH);
    double rho = *rhoH;

    std::auto_ptr<pat::MuonRefVector> outputmuons(new pat::MuonRefVector);
    std::auto_ptr<pat::ElectronRefVector> outputelectrons(new pat::ElectronRefVector);
    std::auto_ptr<pat::ElectronRefVector> outputelectronsnew(new pat::ElectronRefVector);
    std::auto_ptr<pat::MuonRefVector> outputtightmuons(new pat::MuonRefVector);
    std::auto_ptr<pat::ElectronRefVector> outputtightelectrons(new pat::ElectronRefVector);
    std::auto_ptr<pat::PhotonRefVector> outputphotons(new pat::PhotonRefVector);

    for (vector<pat::Muon>::const_iterator muons_iter = muonsH->begin(); muons_iter != muonsH->end(); ++muons_iter) {
        if (verticesH->size() == 0) continue;
        bool passeskincuts = (muons_iter->pt() > 10 && fabs(muons_iter->eta()) < 2.4);
        bool passestipcut = false;
        bool passeslipcut = false;
        if (muons_iter->innerTrack().isAvailable() && fabs(muons_iter->innerTrack()->dxy(verticesH->begin()->position())) < d0cut) passestipcut = true;
        if (muons_iter->innerTrack().isAvailable() && fabs(muons_iter->innerTrack()->dz (verticesH->begin()->position())) < dzcut) passeslipcut = true;
        float isoval = muons_iter->pfIsolationR04().sumNeutralHadronEt;
        isoval += muons_iter->pfIsolationR04().sumPhotonEt;
        isoval -= 0.5*muons_iter->pfIsolationR04().sumPUPt;
        if (isoval < 0.) isoval = 0.;
        isoval += muons_iter->pfIsolationR04().sumChargedHadronPt;
        isoval /= muons_iter->pt();
        bool passesselection = (muon::isLooseMuon(*muons_iter) && isoval < 0.2);

        if (debug && muons_iter->pt() > 10) {
            std::cout << "**** Muon Info ****" << std::endl;
            std::cout << "pT : " << muons_iter->pt() << ", eta : " << muons_iter->eta() << ", phi : " << muons_iter->phi() << std::endl;
            std::cout << "chiso  : " << muons_iter->pfIsolationR04().sumNeutralHadronEt       << std::endl;
            std::cout << "nhiso  : " << muons_iter->pfIsolationR04().sumPhotonEt              << std::endl;
            std::cout << "phiso  : " << muons_iter->pfIsolationR04().sumChargedHadronPt       << std::endl;
            std::cout << "puiso  : " << muons_iter->pfIsolationR04().sumPUPt                  << std::endl;
            std::cout << "isoval : " << isoval                                                << std::endl;
            std::cout << "ID (V) : " << muon::isLooseMuon(*muons_iter)                        << std::endl;
            std::cout << "ID (T) : " << muon::isTightMuon(*muons_iter, *(verticesH->begin())) << std::endl;
            std::cout << std::endl;
        }

        if (passeskincuts && passestipcut && passeslipcut && passesselection) {
            outputmuons->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
            if (muon::isTightMuon(*muons_iter, *(verticesH->begin())) && isoval < 0.12) outputtightmuons->push_back(pat::MuonRef(muonsH, muons_iter - muonsH->begin()));
        }
    }

    for (vector<pat::Electron>::const_iterator electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) {
        if (verticesH->size() == 0) continue;
        bool passeskincuts = (electrons_iter->pt() > 10 && fabs(electrons_iter->superCluster()->eta()) < 2.5);
        bool passestipcut = false;
        bool passeslipcut = false;
        if (electrons_iter->gsfTrack().isAvailable() && fabs(electrons_iter->gsfTrack()->dxy(verticesH->begin()->position())) < d0cut) passestipcut = true;
        if (electrons_iter->gsfTrack().isAvailable() && fabs(electrons_iter->gsfTrack()->dz (verticesH->begin()->position())) < dzcut) passeslipcut = true;
        double chiso = electrons_iter->pfIsolationVariables().sumChargedHadronPt;
        double nhiso = electrons_iter->pfIsolationVariables().sumNeutralHadronEt;
        double phiso = electrons_iter->pfIsolationVariables().sumPhotonEt;
        bool passesselection = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::VETO, *electrons_iter, conversionsH, *(beamspotH.product()), verticesH, chiso, phiso, nhiso, rho);
        if (passeskincuts && passestipcut && passeslipcut && passesselection) {
            outputelectrons->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
            bool passestightid = EgammaCutBasedEleId::PassWP(EgammaCutBasedEleId::TIGHT, *electrons_iter, conversionsH, *(beamspotH.product()), verticesH, chiso, phiso, nhiso, rho);
            bool passestighttrig = EgammaCutBasedEleId::PassTriggerCuts(EgammaCutBasedEleId::TRIGGERTIGHT, *electrons_iter);
            if (passestightid && passestighttrig) outputtightelectrons->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
        }
    }

    for (vector<pat::Electron>::const_iterator electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) {
        if (verticesH->size() == 0) continue;
        bool passeskincuts = (electrons_iter->pt() > 10 && fabs(electrons_iter->superCluster()->eta()) < 2.5);

        bool isBarrel = (fabs(electrons_iter->superCluster()->eta()) <= 1.479 ? true : false);

        bool passestipcut = false;
        bool passeslipcut = false;
        float d0val = fabs(electrons_iter->gsfTrack()->dxy(verticesH->begin()->position()));
        float dzval = fabs(electrons_iter->gsfTrack()->dz (verticesH->begin()->position()));
        if ((isBarrel && d0val < 0.0250) || (!isBarrel && d0val < 0.2232)) passestipcut = true;
        if ((isBarrel && dzval < 0.5863) || (!isBarrel && dzval < 0.9513)) passeslipcut = true;

        bool passesselection = false;
        bool passesid = false;
        bool passesiso = false;
        float dEtaIn = fabs(electrons_iter->deltaEtaSuperClusterTrackAtVtx());
        float dPhiIn = fabs(electrons_iter->deltaPhiSuperClusterTrackAtVtx());
        float sieie  = electrons_iter->full5x5_sigmaIetaIeta();
        float hoe    = electrons_iter->hcalOverEcal();
        float iemip  = 1e30;
        if (electrons_iter->ecalEnergy() != 0 && std::isfinite(electrons_iter->ecalEnergy())) iemip =  fabs(1.0/electrons_iter->ecalEnergy() - electrons_iter->eSuperClusterOverP()/electrons_iter->ecalEnergy());
        int mhits = electrons_iter->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS); 
        bool passesconversionveto = !ConversionTools::hasMatchedConversion(*electrons_iter, conversionsH, beamspotH->position());

        if (isBarrel) {
            if (dEtaIn < 0.0200 && dPhiIn < 0.2579 && sieie < 0.0125 && hoe < 0.2564 && iemip < 0.1508 && passesconversionveto && mhits <= 2) passesid = true;
        }
        else {
            if (dEtaIn < 0.0141 && dPhiIn < 0.2591 && sieie < 0.0371 && hoe < 0.1335 && iemip < 0.1542 && passesconversionveto && mhits <= 3) passesid = true;
        }

        double chiso = electrons_iter->pfIsolationVariables().sumChargedHadronPt;
        double nhiso = electrons_iter->pfIsolationVariables().sumNeutralHadronEt;
        double phiso = electrons_iter->pfIsolationVariables().sumPhotonEt;
        double puiso = electrons_iter->pfIsolationVariables().sumPUPt;
        double iso   = nhiso + phiso - 0.5*puiso;
        if (iso < 0.0) iso = 0.0;
        iso += chiso;
        iso /= electrons_iter->pt();
        if (isBarrel) {
            if (iso <= 0.3313) passesiso = true;
        }            
        else {
            if (iso <= 0.3816) passesiso = true;
        }
       
        passesselection = passesid && passesiso; 

        if (debug && electrons_iter->pt() > 10) {
            std::cout << "**** Electron Info ****" << std::endl;
            std::cout << "pT : " << electrons_iter->pt() << ", eta : " << electrons_iter->eta() << ", phi : " << electrons_iter->phi() << std::endl;
            std::cout << "dEtaIn : " << dEtaIn               << std::endl;
            std::cout << "dPhiIn : " << dPhiIn               << std::endl;
            std::cout << "sieie  : " << sieie                << std::endl;
            std::cout << "hoe    : " << hoe                  << std::endl;
            std::cout << "iemip  : " << iemip                << std::endl;
            std::cout << "mhits  : " << mhits                << std::endl;
            std::cout << "conv   : " << passesconversionveto << std::endl;
            std::cout << "d0     : " << d0val                << std::endl;
            std::cout << "dz     : " << dzval                << std::endl;
            std::cout << "chiso  : " << chiso                << std::endl;
            std::cout << "nhiso  : " << nhiso                << std::endl;
            std::cout << "phiso  : " << phiso                << std::endl;
            std::cout << "puiso  : " << puiso                << std::endl;
            std::cout << "isoval : " << iso                  << std::endl;
            std::cout << "ID (V) : " << passesid             << std::endl;
            std::cout << "Iso    : " << passesiso            << std::endl;
            std::cout << std::endl;
        }

        if (passeskincuts && passestipcut && passeslipcut && passesselection) {
            outputelectronsnew->push_back(pat::ElectronRef(electronsH, electrons_iter - electronsH->begin()));
        }
    }

    for (vector<pat::Photon>::const_iterator photons_iter = photonsH->begin(); photons_iter != photonsH->end(); ++photons_iter) {
        bool passeskincuts = (photons_iter->pt() > 30 && fabs(photons_iter->eta()) < 2.5);
        bool isconversionssafe = true;
        for (vector<pat::Electron>::const_iterator electrons_iter = electronsH->begin(); electrons_iter != electronsH->end(); ++electrons_iter) {
            if (electrons_iter->superCluster() != photons_iter->superCluster()) continue;
            if (electrons_iter->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) > 0) continue;
            if (ConversionTools::hasMatchedConversion(*electrons_iter, conversionsH, beamspotH->position())) continue;
            isconversionssafe = false;
        }
        Ref<vector<pat::Photon> > photonR(photonsH, photons_iter - photonsH->begin());
        double chiso = photons_iter->chargedHadronIso();
        double nhiso = photons_iter->neutralHadronIso();
        double phiso = photons_iter->photonIso();
        bool passesiso = testPhotonIsolation(*photons_iter, chiso, nhiso, phiso, rho);
        bool passesselection = (passesiso && photons_iter->hadTowOverEm() < 0.05 && photons_iter->r9() > 0.9 &&  
                               ((photons_iter->isEB() && photons_iter->sigmaIetaIeta() < 0.011) || (photons_iter->isEE() && photons_iter->sigmaIetaIeta() < 0.033)));

        if (debug && photons_iter->pt() > 30) {
            std::cout << "**** Photon Info ****" << std::endl;
            std::cout << "pT : " << photons_iter->pt() << ", eta : " << photons_iter->eta() << ", phi : " << photons_iter->phi() << std::endl;
            std::cout << "R9     : " << photons_iter->r9()            << std::endl;
            std::cout << "hoe    : " << photons_iter->hadTowOverEm()  << std::endl;
            std::cout << "sieie  : " << photons_iter->sigmaIetaIeta() << std::endl;
            std::cout << "conv   : " << isconversionssafe             << std::endl;
            std::cout << "chiso  : " << chiso                << std::endl;
            std::cout << "nhiso  : " << nhiso                << std::endl;
            std::cout << "phiso  : " << phiso                << std::endl;
            std::cout << "rho    : " << rho                  << std::endl;
            std::cout << "ID,Iso : " << passesselection      << std::endl;
            std::cout << std::endl;
        }

        if (passeskincuts && isconversionssafe && passesselection) {
            outputphotons->push_back(pat::PhotonRef(photonsH, photons_iter - photonsH->begin()));
        }
    }

    iEvent.put(outputmuons, "muons");
    iEvent.put(outputelectrons, "electrons");
    iEvent.put(outputelectronsnew, "electronsnew");
    iEvent.put(outputtightmuons, "tightmuons");
    iEvent.put(outputtightelectrons, "tightelectrons");
    iEvent.put(outputphotons, "photons");
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
    if (fabs(eta) < 1.0) return 0.012;
    else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.010;
    else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.014;
    else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.012;
    else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.016;
    else if (fabs(eta) >= 2.3   && fabs(eta) < 2.3  ) return 0.020;
    else if (fabs(eta) >= 2.4) return 0.012;
    else return 0.;
}

double PFCleaner::getNeutralHadronEAForPhotonIso(double eta) {
    if (fabs(eta) < 1.0) return 0.030;
    else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.057;
    else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.039;
    else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.015;
    else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.024;
    else if (fabs(eta) >= 2.3   && fabs(eta) < 2.3  ) return 0.039;
    else if (fabs(eta) >= 2.4) return 0.072;
    else return 0.;
}

double PFCleaner::getGammaEAForPhotonIso(double eta) {
    if (fabs(eta) < 1.0) return 0.148;
    else if (fabs(eta) >= 1.0   && fabs(eta) < 1.479) return 0.130;
    else if (fabs(eta) >= 1.479 && fabs(eta) < 2.0  ) return 0.112;
    else if (fabs(eta) >= 2.0   && fabs(eta) < 2.2  ) return 0.216;
    else if (fabs(eta) >= 2.2   && fabs(eta) < 2.3  ) return 0.262;
    else if (fabs(eta) >= 2.3   && fabs(eta) < 2.3  ) return 0.260;
    else if (fabs(eta) >= 2.4) return 0.266;
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
        if (corrCHIso < 1.5 && corrNHIso < 1.0 + 0.04*photon.pt() && corrPHIso < 0.7 + 0.005*photon.pt()) return true;
        else return false;
    }
    else if (photon.isEE()) {
        if (corrCHIso < 1.2 && corrNHIso < 1.5 + 0.04*photon.pt() && corrPHIso < 1.0 + 0.005*photon.pt()) return true;
        else return false;
    }
    else return false;
}

DEFINE_FWK_MODULE(PFCleaner);
