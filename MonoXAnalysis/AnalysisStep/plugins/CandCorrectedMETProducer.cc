/*

Notes:

Muon ID and isolation : https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId
Photon ID and isolation : https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonID2012
Electron ID and isolation : https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgammaCutBasedIdentification

*/


#include <memory>
#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"

class CandCorrectedMETProducer : public edm::EDProducer {
    public:
        explicit CandCorrectedMETProducer(const edm::ParameterSet&);
        ~CandCorrectedMETProducer();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
        
        virtual void beginRun(edm::Run&, edm::EventSetup const&);
        virtual void endRun(edm::Run&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

        edm::InputTag metTag;
        std::vector<edm::InputTag> candTags;
};

CandCorrectedMETProducer::CandCorrectedMETProducer(const edm::ParameterSet& iConfig): 
    metTag(iConfig.getParameter<edm::InputTag>("met")),
    candTags(iConfig.getParameter<std::vector<edm::InputTag> >("cands"))
{
    produces<reco::METCollection>();
}


CandCorrectedMETProducer::~CandCorrectedMETProducer() {
}

void CandCorrectedMETProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    Handle<View<MET> > metH;
    iEvent.getByLabel(metTag, metH);

    std::vector<Handle<View<Candidate> > > candHs;
    for (std::size_t i = 0; i < candTags.size(); i++) {
        Handle<View<Candidate> > candH;
        iEvent.getByLabel(candTags[i], candH);
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

void CandCorrectedMETProducer::beginJob() {
}

void CandCorrectedMETProducer::endJob() {
}

void CandCorrectedMETProducer::beginRun(edm::Run&, edm::EventSetup const&) {
}

void CandCorrectedMETProducer::endRun(edm::Run&, edm::EventSetup const&) {
}

void CandCorrectedMETProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}

void CandCorrectedMETProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}

void CandCorrectedMETProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(CandCorrectedMETProducer);