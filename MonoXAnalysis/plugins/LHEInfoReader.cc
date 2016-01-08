#include <memory>
#include <iostream>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

class LHEInfoReader : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
    public:
        explicit LHEInfoReader(const edm::ParameterSet&);
        ~LHEInfoReader();
        
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
    
    private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;
        
        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
};

LHEInfoReader::LHEInfoReader(const edm::ParameterSet& iConfig) {
}


LHEInfoReader::~LHEInfoReader() {
}

void LHEInfoReader::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {}

void LHEInfoReader::beginJob() {}
void LHEInfoReader::endJob() {}
void LHEInfoReader::beginRun(edm::Run const& iRun, edm::EventSetup const&) {}

void LHEInfoReader::endRun(edm::Run const& iRun, edm::EventSetup const&) {

    edm::Handle<LHERunInfoProduct> run;    
    iRun.getByLabel("externalLHEProducer", run);
    LHERunInfoProduct myLHERunInfoProduct = *(run.product());
    
    for (auto iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
        std::cout << iter->tag() << std::endl;
        std::vector<std::string> lines = iter->lines();
        for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
            std::cout << lines.at(iLine);
        }
    }
}

void LHEInfoReader::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(LHEInfoReader);
