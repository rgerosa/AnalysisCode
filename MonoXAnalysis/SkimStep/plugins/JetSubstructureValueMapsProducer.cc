#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"

#include "MonoXAnalysis/SkimStep/interface/Njettiness.hh"

class JetSubstructureValueMapsProducer : public edm::EDProducer {
    public:
        explicit JetSubstructureValueMapsProducer(const edm::ParameterSet&);
        ~JetSubstructureValueMapsProducer();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&);
        virtual void endJob() ;
        
        virtual void beginRun(edm::Run&, edm::EventSetup const&);
        virtual void endRun(edm::Run&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

        edm::InputTag jetsTag;
        double jetRadius;
};

JetSubstructureValueMapsProducer::JetSubstructureValueMapsProducer(const edm::ParameterSet& iConfig):
    jetsTag(iConfig.getParameter<edm::InputTag>("src")),
    jetRadius(iConfig.getParameter<double>("jetRadius"))
{
    produces<edm::ValueMap<float> >("tau3");
    produces<edm::ValueMap<float> >("tau2");
    produces<edm::ValueMap<float> >("tau1");
}


JetSubstructureValueMapsProducer::~JetSubstructureValueMapsProducer() {
}

void JetSubstructureValueMapsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    Handle<PFJetCollection> jetsH;
    iEvent.getByLabel(jetsTag, jetsH);
    PFJetCollection jets = *jetsH;

    std::auto_ptr< edm::ValueMap<float> > tau3ValueMap(new edm::ValueMap<float>());
    std::auto_ptr< edm::ValueMap<float> > tau2ValueMap(new edm::ValueMap<float>());
    std::auto_ptr< edm::ValueMap<float> > tau1ValueMap(new edm::ValueMap<float>());

    edm::ValueMap<float>::Filler tau3ValueMapFiller(*tau3ValueMap);
    edm::ValueMap<float>::Filler tau2ValueMapFiller(*tau2ValueMap);
    edm::ValueMap<float>::Filler tau1ValueMapFiller(*tau1ValueMap);

    std::vector<float> tau3Values(jetsH->size(), -1.0);
    std::vector<float> tau2Values(jetsH->size(), -1.0);
    std::vector<float> tau1Values(jetsH->size(), -1.0);

    std::vector<fastjet::PseudoJet> FJConstituents;
    
    for (size_t i = 0; i < jets.size(); i++) {
        FJConstituents.clear();
        vector<Ptr<PFCandidate> > jetConstituents = jets[i].getPFConstituents();
        for (size_t j = 0; j < jetConstituents.size(); j++) FJConstituents.push_back(fastjet::PseudoJet(jetConstituents[j]->px(), jetConstituents[j]->py(), jetConstituents[j]->pz(), jetConstituents[j]->energy()));

        NsubParameters paramsNjettiness = NsubParameters(1.0, jetRadius);
        Njettiness moduleNjettiness(Njettiness::onepass_kt_axes, paramsNjettiness);
        tau3Values[i] = moduleNjettiness.getTau(3, FJConstituents);
        tau2Values[i] = moduleNjettiness.getTau(2, FJConstituents);
        tau1Values[i] = moduleNjettiness.getTau(1, FJConstituents);
    }

    tau3ValueMapFiller.insert(jetsH, tau3Values.begin(), tau3Values.end());
    tau2ValueMapFiller.insert(jetsH, tau2Values.begin(), tau2Values.end());
    tau1ValueMapFiller.insert(jetsH, tau1Values.begin(), tau1Values.end());

    tau3ValueMapFiller.fill();
    tau2ValueMapFiller.fill();
    tau1ValueMapFiller.fill();

    iEvent.put(tau3ValueMap, "tau3");
    iEvent.put(tau2ValueMap, "tau2");
    iEvent.put(tau1ValueMap, "tau1");
}

void JetSubstructureValueMapsProducer::beginJob() {
}

void JetSubstructureValueMapsProducer::endJob() {
}

void JetSubstructureValueMapsProducer::beginRun(edm::Run&, edm::EventSetup const&) {
}

void JetSubstructureValueMapsProducer::endRun(edm::Run&, edm::EventSetup const&) {
}

void JetSubstructureValueMapsProducer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}

void JetSubstructureValueMapsProducer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {
}

void JetSubstructureValueMapsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(JetSubstructureValueMapsProducer);
