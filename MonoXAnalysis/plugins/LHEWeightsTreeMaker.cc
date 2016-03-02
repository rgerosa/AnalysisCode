// basic C++ headers
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>

// FWCore headers
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// additional
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// ROOT headers
#include <TTree.h>

class LHEWeightsTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources,edm::one::WatchRuns> {

public:
  explicit LHEWeightsTreeMaker(const edm::ParameterSet&);
  ~LHEWeightsTreeMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  
  // InputTags
  const edm::InputTag lheInfoTag;
  const edm::InputTag lheRunInfoTag;
  const edm::InputTag genInfoTag;
  const edm::InputTag pileupInfoTag;
  
  // Tokens
  edm::EDGetTokenT<LHEEventProduct> lheInfoToken;
  edm::EDGetTokenT<GenEventInfoProduct> genInfoToken;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfoToken;
  
  const bool uselheweights, addqcdpdfweights;

  TTree* tree;
  
  uint32_t event, run, lumi;
  int      puobs, putrue;
  double   wgtsign, wgtxsec, wgtpdf1, wgtpdf2, wgtpdf3, wgtpdf4, wgtpdf5;
  double   lheXSEC;
  std::auto_ptr<double>  wgtpdf;
  std::auto_ptr<double>  wgtqcd;
};

LHEWeightsTreeMaker::LHEWeightsTreeMaker(const edm::ParameterSet& iConfig): 
  lheInfoTag(iConfig.getParameter<edm::InputTag>("lheinfo")),
  lheRunInfoTag(iConfig.getParameter<edm::InputTag>("lheRuninfo")),
  genInfoTag(iConfig.getParameter<edm::InputTag>("geninfo")),
  pileupInfoTag(iConfig.getParameter<edm::InputTag>("pileupinfo")),
  uselheweights(iConfig.getParameter<bool>("uselheweights")),
  addqcdpdfweights(iConfig.getParameter<bool>("addqcdpdfweights"))
{
  // Token consumes instructions
  lheInfoToken = consumes<LHEEventProduct>(lheInfoTag);
  genInfoToken = consumes<GenEventInfoProduct>(genInfoTag);
  pileupInfoToken = consumes<std::vector<PileupSummaryInfo> >(pileupInfoTag);

  wgtqcd = std::auto_ptr<double> (new double[8]);
  wgtpdf = std::auto_ptr<double> (new double[100]);

  for (size_t i = 0; i < 8  ; i++) wgtqcd.get()[i] = 0.;
  for (size_t i = 0; i < 100; i++) wgtpdf.get()[i] = 0.;
  
  // state that TFileService is used
  usesResource();
  usesResource("TFileService");

}


LHEWeightsTreeMaker::~LHEWeightsTreeMaker() {}

void LHEWeightsTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  using namespace edm;
  using namespace std;

  // Get handles to all the requisite collections
  Handle<LHEEventProduct> lheInfoH;
  if (uselheweights) 
    iEvent.getByToken(lheInfoToken, lheInfoH);
    
  Handle<GenEventInfoProduct> genInfoH;
  if (uselheweights) 
    iEvent.getByToken(genInfoToken, genInfoH);

  Handle<std::vector<PileupSummaryInfo> > pileupInfoH;
  iEvent.getByToken(pileupInfoToken,pileupInfoH);

  // Event, lumi, run info
  event = iEvent.id().event();
  run   = iEvent.id().run();
  lumi  = iEvent.luminosityBlock();

  puobs = 0;
  putrue = 0;

  if(pileupInfoH.isValid()) {
    for (auto pileupInfo_iter = pileupInfoH->begin(); pileupInfo_iter != pileupInfoH->end(); ++pileupInfo_iter) {
      if (pileupInfo_iter->getBunchCrossing() == 0) {
	puobs  = pileupInfo_iter->getPU_NumInteractions();
	putrue = pileupInfo_iter->getTrueNumInteractions();
      }
    }
  }
  

  // Weights info
  wgtsign = 1.0;
  wgtxsec = 1.0;
  if (uselheweights) {
    wgtsign = genInfoH->weight();
    wgtxsec = lheInfoH->originalXWGTUP();
  }
  wgtpdf1 = 0.0;
  wgtpdf2 = 0.0;
  wgtpdf3 = 0.0;
  wgtpdf4 = 0.0;
  wgtpdf5 = 0.0;
  
  if (addqcdpdfweights) {

      vector<gen::WeightsInfo> weights = lheInfoH->weights();
      for (size_t i = 0; i < weights.size(); i++) {
	
	if (weights[i].id == "315")      wgtpdf1 = weights[i].wgt; // cteq6l1
	else if (weights[i].id == "316") wgtpdf2 = weights[i].wgt; // MMHT2014lo68cl
	else if (weights[i].id == "370") wgtpdf3 = weights[i].wgt; // HERAPDF15LO
	else if (weights[i].id == "393") wgtpdf4 = weights[i].wgt; // CT10nlo
	else if (weights[i].id == "446") wgtpdf5 = weights[i].wgt; // MMHT2014nlo68cl

	else if(std::stoi(weights[i].id) >= 2 and std::stoi(weights[i].id) <=9)
	  wgtqcd.get()[std::stoi(weights[i].id)-2] = weights[i].wgt;
	else if(std::stoi(weights[i].id) >= 11 and std::stoi(weights[i].id) <=110)
	  wgtpdf.get()[std::stoi(weights[i].id)-11] = weights[i].wgt;
	else
	  continue;
      }
  }
  
  tree->Fill();

}


void LHEWeightsTreeMaker::beginJob() {

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("gentree"    , "gentree");
  // Run, Lumi, Event info
  tree->Branch("event"                , &event                , "event/i");
  tree->Branch("run"                  , &run                  , "run/i");
  tree->Branch("lumi"                 , &lumi                 , "lumi/i");
  // Event weights
  tree->Branch("wgtsign"              , &wgtsign              , "wgtsign/D");
  tree->Branch("wgtxsec"              , &wgtxsec              , "wgtxsec/D");
  // pileup info
  tree->Branch("puobs"                , &puobs                , "puobs/I");
  tree->Branch("putrue"               , &putrue               , "putrue/I");

  if (addqcdpdfweights) {
    tree->Branch("wgtpdf1"              , &wgtpdf1              , "wgtpdf1/D");
    tree->Branch("wgtpdf2"              , &wgtpdf2              , "wgtpdf2/D");
    tree->Branch("wgtpdf3"              , &wgtpdf3              , "wgtpdf3/D");
    tree->Branch("wgtpdf4"              , &wgtpdf4              , "wgtpdf4/D");
    tree->Branch("wgtpdf5"              , &wgtpdf5              , "wgtpdf5/D");
    tree->Branch("wgtpdf"               ,  wgtpdf.get()         , "wgtpdf[100]/D");
    tree->Branch("wgtqcd"               ,  wgtqcd.get()         , "wgtqcd[8]/D");
  }

  tree->Branch("lheXSEC",   &lheXSEC, "lheXSEC/D");
}

void LHEWeightsTreeMaker::endJob() {}

void LHEWeightsTreeMaker::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  // in cae of MC store XS value                                                                                                                                               
  edm::Handle<LHERunInfoProduct> run;
  iRun.getByLabel(lheRunInfoTag,run);
  LHERunInfoProduct myLHERunInfoProduct = *(run.product());
  lheXSEC = myLHERunInfoProduct.heprup().XSECUP.at(0);  
}

void LHEWeightsTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {}

void LHEWeightsTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(LHEWeightsTreeMaker);

