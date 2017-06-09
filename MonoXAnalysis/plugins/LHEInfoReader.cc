#include <memory>
#include <iostream>
#include <fstream>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Run.h" 
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

class LHEInfoReader : public edm::one::EDAnalyzer<edm::one::WatchRuns> {
public:
  explicit LHEInfoReader(const edm::ParameterSet &);
  ~LHEInfoReader();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  std::string outputLHEFileName_ ;

};

LHEInfoReader::LHEInfoReader(const edm::ParameterSet& iConfig) {
  outputLHEFileName_ = iConfig.getParameter<std::string>("outputLHEFileName");

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

    std::ofstream outfile(outputLHEFileName_.c_str());

    for (auto iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
      outfile << iter->tag();
      std::vector<std::string> lines = iter->lines();
      for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
	outfile << lines.at(iLine);
      }
    }

    outfile.close();
}

void LHEInfoReader::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(LHEInfoReader);
