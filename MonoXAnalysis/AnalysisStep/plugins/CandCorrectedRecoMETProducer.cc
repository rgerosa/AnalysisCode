#include <memory>
#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"

class CandCorrectedRecoMETProducer : public edm::EDProducer {
    public:
        explicit CandCorrectedRecoMETProducer(const edm::ParameterSet&);
        ~CandCorrectedRecoMETProducer();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    private:
        virtual void beginJob() override;
        virtual void produce(edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        
        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

        edm::InputTag    metTag;
        std::vector<edm::InputTag> candTags;

        edm::EDGetTokenT<edm::View<reco::MET> > metToken;
        edm::EDGetTokenT<edm::View<reco::Candidate> > theCandToken;
        std::vector<edm::EDGetTokenT<edm::View<reco::Candidate> > > candTokens;
};

CandCorrectedRecoMETProducer::CandCorrectedRecoMETProducer(const edm::ParameterSet& iConfig): 
    metTag(iConfig.getParameter<edm::InputTag>("met")),
    candTags(iConfig.getParameter<std::vector<edm::InputTag> >("cands"))
{
    produces<reco::METCollection>();
    
    metToken = consumes<edm::View<reco::MET> > (metTag);
    for (std::size_t i = 0; i < candTags.size(); i++) {
        theCandToken = consumes<edm::View<reco::Candidate> > (candTags[i]);
        candTokens.push_back(theCandToken);
    }
}


CandCorrectedRecoMETProducer::~CandCorrectedRecoMETProducer() {
}

void CandCorrectedRecoMETProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    Handle<View<reco::MET> > metH;
    iEvent.getByToken(metToken, metH);

    std::vector<Handle<View<Candidate> > > candHs;
    for (std::size_t i = 0; i < candTokens.size(); i++) {
        Handle<View<reco::Candidate> > candH;
        iEvent.getByToken(candTokens[i], candH);
        candHs.push_back(candH);
    }

    std::auto_ptr<METCollection> output(new METCollection);

    double met    = metH->front().et();
    double metphi = metH->front().phi();
    
    double ccmetx = met * cos(metphi);
    double ccmety = met * sin(metphi);

    for (size_t i = 0; i < candHs.size(); i++) {
        for (View<Candidate>::const_iterator cands_iter = candHs[i]->begin(); cands_iter != candHs[i]->end(); ++cands_iter) {
            ccmetx += cands_iter->pt() * cos(cands_iter->phi());
            ccmety += cands_iter->pt() * sin(cands_iter->phi());
        }            
    }
    double ccmet = sqrt(ccmetx*ccmetx + ccmety*ccmety);

    MET* ccmetcand = metH->front().clone();
    ccmetcand->setP4(reco::Candidate::LorentzVector(ccmetx, ccmety, 0., ccmet));
    output->push_back(*ccmetcand);

    iEvent.put(output);
}

void CandCorrectedRecoMETProducer::beginJob() {
}

void CandCorrectedRecoMETProducer::endJob() {
}

void CandCorrectedRecoMETProducer::beginRun(edm::Run const&, edm::EventSetup const&) {
}

void CandCorrectedRecoMETProducer::endRun(edm::Run const&, edm::EventSetup const&) {
}

void CandCorrectedRecoMETProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void CandCorrectedRecoMETProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void CandCorrectedRecoMETProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(CandCorrectedRecoMETProducer);
