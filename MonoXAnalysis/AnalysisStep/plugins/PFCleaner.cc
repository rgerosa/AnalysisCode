#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

class PFCleaner : public edm::EDProducer {
    public:
        explicit PFCleaner(const edm::ParameterSet&);
        ~PFCleaner();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
        
        virtual void beginRun(edm::Run&, edm::EventSetup const&);
        virtual void endRun(edm::Run&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

        edm::InputTag src;
        edm::InputTag vertices;
        edm::InputTag pfpu;
        edm::InputTag rhoTag;
        double d0cut, dzcut;
        double muisocut, eleisocut;

        StringCutObjectSelector<reco::Muon> muselector;
        StringCutObjectSelector<reco::GsfElectron> eleselector;
        StringCutObjectSelector<reco::Photon> photonselector;

        bool vetophotons;
};

PFCleaner::PFCleaner(const edm::ParameterSet& iConfig): 
    src(iConfig.getParameter<edm::InputTag>("src")),
    vertices(iConfig.getParameter<edm::InputTag>("vertices")),
    pfpu(iConfig.getParameter<edm::InputTag>("pfpileup")),
    rhoTag(iConfig.getParameter<edm::InputTag>("rho")),
    d0cut(iConfig.getParameter<double>("d0cut")),
    dzcut(iConfig.getParameter<double>("dzcut")),
    muisocut(iConfig.getParameter<double>("muisocut")),
    eleisocut(iConfig.getParameter<double>("eleisocut")),
    muselector(iConfig.getParameter<std::string>("muselection")),
    eleselector(iConfig.getParameter<std::string>("eleselection")),
    photonselector(iConfig.getParameter<std::string>("photonselection")),
    vetophotons(iConfig.getParameter<bool>("vetophotons"))
{
    produces<reco::PFCandidateCollection>("pfcands");
    produces<reco::MuonCollection>("muons");
    produces<reco::GsfElectronCollection>("electrons");
    produces<reco::PhotonCollection>("photons");
}


PFCleaner::~PFCleaner() {
}

void PFCleaner::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    Handle<PFCandidateCollection> inputH;
    iEvent.getByLabel(src, inputH);
    PFCandidateCollection input = *inputH;    

    Handle<View<PFCandidate> > pfpileupH;
    iEvent.getByLabel(pfpu, pfpileupH);

    Handle<vector<Vertex> > verticesH;
    iEvent.getByLabel(vertices, verticesH);

    Handle<double> rhoH;
    iEvent.getByLabel(rhoTag, rhoH);
    double rho = *rhoH;

    std::auto_ptr<PFCandidateCollection> output(new PFCandidateCollection);
    std::auto_ptr<MuonCollection> outputmuons(new MuonCollection);
    std::auto_ptr<GsfElectronCollection> outputelectrons(new GsfElectronCollection);
    std::auto_ptr<PhotonCollection> outputphotons(new PhotonCollection);

    for (size_t i = 0; i  < input.size(); i++) {
        bool veto = false;

        for (View<PFCandidate>::const_iterator pfpileup_iter = pfpileupH->begin(); pfpileup_iter != pfpileupH->end(); ++pfpileup_iter) {
            if (input[i].pdgId() == pfpileup_iter->pdgId() && input[i].pt() == pfpileup_iter->pt() && input[i].eta() == pfpileup_iter->eta() && input[i].phi() == pfpileup_iter->phi()) veto = true;
        }

        bool isMu = (abs(input[i].pdgId()) == 13 && input[i].muonRef().isAvailable());
        if (isMu && verticesH->size() > 0) {
            bool passestipcut = false;
            bool passeslipcut = false;
            if (input[i].muonRef()->innerTrack().isAvailable() && fabs(input[i].muonRef()->innerTrack()->dxy(verticesH->begin()->position())) < d0cut) passestipcut = true;
            if (input[i].muonRef()->innerTrack().isAvailable() && fabs(input[i].muonRef()->innerTrack()->dz (verticesH->begin()->position())) < dzcut) passeslipcut = true;
            bool passesselection = muselector(*(input[i].muonRef()));
            float isoval = input[i].muonRef()->pfIsolationR04().sumNeutralHadronEt;
            isoval += input[i].muonRef()->pfIsolationR04().sumPhotonEt;
            isoval -= 0.5*input[i].muonRef()->pfIsolationR04().sumPUPt;
            if (isoval < 0.) isoval = 0.;
            isoval += input[i].muonRef()->pfIsolationR04().sumChargedHadronPt;
            isoval /= input[i].pt();
            bool passesisolation = (isoval < muisocut);
            if (passestipcut && passeslipcut && passesselection && passesisolation) {
                veto = true;
                outputmuons->push_back(*(input[i].muonRef()));
            }
        }

        bool isEl = (abs(input[i].pdgId()) == 11 && input[i].gsfElectronRef().isAvailable());
        if (isEl && verticesH->size() > 0) {
            bool passestipcut = false;
            bool passeslipcut = false;
            if (input[i].gsfElectronRef()->gsfTrack().isAvailable() && fabs(input[i].gsfElectronRef()->gsfTrack()->dxy(verticesH->begin()->position())) < d0cut) passestipcut = true;
            if (input[i].gsfElectronRef()->gsfTrack().isAvailable() && fabs(input[i].gsfElectronRef()->gsfTrack()->dz (verticesH->begin()->position())) < dzcut) passeslipcut = true;
            bool passesselection = eleselector(*(input[i].gsfElectronRef()));
            float isoval = input[i].gsfElectronRef()->pfIsolationVariables().neutralHadronIso;
            isoval += input[i].gsfElectronRef()->pfIsolationVariables().photonIso;
            for (View<PFCandidate>::const_iterator pfpileup_iter = pfpileupH->begin(); pfpileup_iter != pfpileupH->end(); ++pfpileup_iter) {
                if (pfpileup_iter->charge() == 0) continue;
                if (deltaR(pfpileup_iter->eta(), pfpileup_iter->phi(), input[i].eta(), input[i].phi()) > 0.4) continue;
                isoval -= 0.5*(pfpileup_iter->pt());
            }
            if (isoval < 0.) isoval = 0.;
            isoval += input[i].gsfElectronRef()->pfIsolationVariables().chargedHadronIso;
            isoval /= input[i].pt();
            bool passesisolation = (isoval < eleisocut);
            if (passestipcut && passeslipcut && passesselection && passesisolation) {
                veto = true;
                outputelectrons->push_back(*(input[i].gsfElectronRef()));
            }
        }

        bool isPh = (abs(input[i].pdgId()) == 22 && input[i].photonRef().isAvailable());
        if (isPh && verticesH->size() > 0) {
            bool passesselection = photonselector(*(input[i].photonRef()));
            float tkisoval       = input[i].photonRef()->trkSumPtHollowConeDR04();
            float ecalisoval     = input[i].photonRef()->ecalRecHitSumEtConeDR04();
            float hcalisoval     = input[i].photonRef()->hcalTowerSumEtConeDR04();
            bool passestkiso     = ((input[i].photonRef()->isEB()&& tkisoval   < 2.0 + input[i].photonRef()->et()*0.0010 + rho*0.0167) || (input[i].photonRef()->isEE()&& tkisoval   < 2.0 + input[i].photonRef()->et()*0.0010 + rho*0.032));
            bool passesecaliso   = ((input[i].photonRef()->isEB()&& ecalisoval < 4.2 + input[i].photonRef()->et()*0.0060 + rho*0.1830) || (input[i].photonRef()->isEE()&& ecalisoval < 4.2 + input[i].photonRef()->et()*0.0060 + rho*0.090));
            bool passeshcaliso   = ((input[i].photonRef()->isEB()&& hcalisoval < 2.2 + input[i].photonRef()->et()*0.0025 + rho*0.0620) || (input[i].photonRef()->isEE()&& hcalisoval < 2.2 + input[i].photonRef()->et()*0.0025 + rho*0.180));
            if (passesselection && passestkiso && passesecaliso && passeshcaliso) {
                if (vetophotons) veto = true;
                outputphotons->push_back(*(input[i].photonRef()));
            }
        }
        if (!veto) {
            PFCandidatePtr ptrToMother(inputH, i);
            output->push_back(input[i]);
            output->back().setSourceCandidatePtr(ptrToMother);
        }
    }

    iEvent.put(output, "pfcands");
    iEvent.put(outputmuons, "muons");
    iEvent.put(outputelectrons, "electrons");
    iEvent.put(outputphotons, "photons");
}

void PFCleaner::beginJob() {
}

void PFCleaner::endJob() {
}

void PFCleaner::beginRun(edm::Run&, edm::EventSetup const&) {
}

void PFCleaner::endRun(edm::Run&, edm::EventSetup const&) {
}

void PFCleaner::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}

void PFCleaner::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}

void PFCleaner::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PFCleaner);
