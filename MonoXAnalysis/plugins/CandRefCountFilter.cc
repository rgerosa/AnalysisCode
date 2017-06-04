#include <memory>
#include <vector>
#include <string>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

/////////////////////////////////////////////
class MuonRefCountFilter : public edm::stream::EDFilter<> {

public:
  explicit MuonRefCountFilter(const edm::ParameterSet&);
  ~MuonRefCountFilter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual bool filter(edm::Event&, const edm::EventSetup&) ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  const bool filterEvents_;
  const bool produceOutputCollection_;
  edm::InputTag src_;
  int minNumber_;
  int maxNumber_;
  std::string selection_ ;
  edm::EDGetTokenT<pat::MuonRefVector>  srcToken_;
  StringCutObjectSelector<pat::Muon>  selectionObj_;
};

MuonRefCountFilter::MuonRefCountFilter(const edm::ParameterSet& iConfig):
  filterEvents_(iConfig.existsAs<bool>("filterEvents") ? iConfig.getParameter<bool>("filterEvents") : true),
  produceOutputCollection_(iConfig.existsAs<bool>("produceOutputCollection") ? iConfig.getParameter<bool>("produceOutputCollection") : false),
  src_(iConfig.getParameter<edm::InputTag>("src")),
  minNumber_(iConfig.getParameter<int>("minNumber")),
  maxNumber_(iConfig.getParameter<int>("maxNumber")),
  selection_(iConfig.existsAs<std::string>("selection") ? iConfig.getParameter<std::string>("selection") : "abs(eta) < 5.0"),
  selectionObj_(StringCutObjectSelector<pat::Muon>(selection_))
{
  srcToken_    = consumes<pat::MuonRefVector> (src_);
  if(produceOutputCollection_)
    produces<pat::MuonCollection>();

}

MuonRefCountFilter::~MuonRefCountFilter() {}

bool MuonRefCountFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // if you don't want to filter just exit
  if(filterEvents_ == false) return true;

  edm::Handle<pat::MuonRefVector> candCollectionH;
  iEvent.getByToken(srcToken_,candCollectionH);
  pat::MuonRefVector candCollection = *candCollectionH;

  std::unique_ptr<pat::MuonCollection> filteredObjects (new pat::MuonCollection);

  int size = 0;
  for(size_t itMuon = 0; itMuon < candCollection.size(); itMuon++){
    pat::MuonRef candRef = candCollection[itMuon];
    const pat::Muon* cand = dynamic_cast<const pat::Muon*> ((*candRef).clone());
    if(not selectionObj_(*cand)) continue;
    size++;
    if(produceOutputCollection_){
      filteredObjects->push_back(*cand);
    }
  }

  bool accept = false;
  if(size >= minNumber_ and size <= maxNumber_)
    accept = true;

  if(produceOutputCollection_)
    iEvent.put(std::move(filteredObjects));

  return accept;
}

void MuonRefCountFilter::beginRun(edm::Run const& iRun , edm::EventSetup const& iSetup) {}
void MuonRefCountFilter::endRun(edm::Run const&, edm::EventSetup const&) {}
void MuonRefCountFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
void MuonRefCountFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
void MuonRefCountFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MuonRefCountFilter);


/////////////////////////////////////////////
class ElectronRefCountFilter : public edm::stream::EDFilter<> {

public:
  explicit ElectronRefCountFilter(const edm::ParameterSet&);
  ~ElectronRefCountFilter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual bool filter(edm::Event&, const edm::EventSetup&) ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  const bool filterEvents_;
  const bool produceOutputCollection_;
  edm::InputTag src_;
  int minNumber_;
  int maxNumber_;
  std::string selection_ ;
  edm::EDGetTokenT<pat::ElectronRefVector>  srcToken_;
  StringCutObjectSelector<pat::Electron>  selectionObj_;
};

ElectronRefCountFilter::ElectronRefCountFilter(const edm::ParameterSet& iConfig):
  filterEvents_(iConfig.existsAs<bool>("filterEvents") ? iConfig.getParameter<bool>("filterEvents") : true),
  produceOutputCollection_(iConfig.existsAs<bool>("produceOutputCollection") ? iConfig.getParameter<bool>("produceOutputCollection") : true),
  src_(iConfig.getParameter<edm::InputTag>("src")),
  minNumber_(iConfig.getParameter<int>("minNumber")),
  maxNumber_(iConfig.getParameter<int>("maxNumber")),
  selection_(iConfig.existsAs<std::string>("selection") ? iConfig.getParameter<std::string>("selection") : "abs(eta) < 5.0"),
  selectionObj_(StringCutObjectSelector<pat::Electron>(selection_))
{
  srcToken_    = consumes<pat::ElectronRefVector> (src_);
  if(produceOutputCollection_)
    produces<pat::ElectronCollection>();

}

ElectronRefCountFilter::~ElectronRefCountFilter() {}

bool ElectronRefCountFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // if you don't want to filter just exit
  if(filterEvents_ == false) return true;

  edm::Handle<pat::ElectronRefVector> candCollectionH;
  iEvent.getByToken(srcToken_,candCollectionH);
  pat::ElectronRefVector candCollection = *candCollectionH;

  std::unique_ptr<pat::ElectronCollection> filteredObjects (new pat::ElectronCollection);

  int size = 0;
  for(size_t itElectron = 0; itElectron < candCollection.size(); itElectron++){
    pat::ElectronRef candRef = candCollection[itElectron];
    const pat::Electron* cand = dynamic_cast<const pat::Electron*> ((*candRef).clone());
    if(not selectionObj_(*cand)) continue;
    size++;
    if(produceOutputCollection_){
      filteredObjects->push_back(*cand);
    }
  }

  bool accept = false;
  if(size >= minNumber_ and size <= maxNumber_)
    accept = true;

  if(produceOutputCollection_)
    iEvent.put(std::move(filteredObjects));

  return accept;
}

void ElectronRefCountFilter::beginRun(edm::Run const& iRun , edm::EventSetup const& iSetup) {}
void ElectronRefCountFilter::endRun(edm::Run const&, edm::EventSetup const&) {}
void ElectronRefCountFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
void ElectronRefCountFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
void ElectronRefCountFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ElectronRefCountFilter);

/////////////////////////////////////////////

class PhotonRefCountFilter : public edm::stream::EDFilter<> {

public:
  explicit PhotonRefCountFilter(const edm::ParameterSet&);
  ~PhotonRefCountFilter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual bool filter(edm::Event&, const edm::EventSetup&) ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  const bool filterEvents_;
  const bool produceOutputCollection_;
  edm::InputTag src_;
  int minNumber_;
  int maxNumber_;
  std::string selection_ ;
  edm::EDGetTokenT<pat::PhotonRefVector>  srcToken_;
  StringCutObjectSelector<pat::Photon>  selectionObj_;
};

PhotonRefCountFilter::PhotonRefCountFilter(const edm::ParameterSet& iConfig):
  filterEvents_(iConfig.existsAs<bool>("filterEvents") ? iConfig.getParameter<bool>("filterEvents") : true),
  produceOutputCollection_(iConfig.existsAs<bool>("produceOutputCollection") ? iConfig.getParameter<bool>("produceOutputCollection") : true),
  src_(iConfig.getParameter<edm::InputTag>("src")),
  minNumber_(iConfig.getParameter<int>("minNumber")),
  maxNumber_(iConfig.getParameter<int>("maxNumber")),
  selection_(iConfig.existsAs<std::string>("selection") ? iConfig.getParameter<std::string>("selection") : "abs(eta) < 5.0"),
  selectionObj_(StringCutObjectSelector<pat::Photon>(selection_))
{
  srcToken_    = consumes<pat::PhotonRefVector> (src_);
  if(produceOutputCollection_)
    produces<pat::PhotonCollection>();

}

PhotonRefCountFilter::~PhotonRefCountFilter() {}

bool PhotonRefCountFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // if you don't want to filter just exit
  if(filterEvents_ == false) return true;

  edm::Handle<pat::PhotonRefVector> candCollectionH;
  iEvent.getByToken(srcToken_,candCollectionH);
  pat::PhotonRefVector candCollection = *candCollectionH;

  std::unique_ptr<pat::PhotonCollection> filteredObjects (new pat::PhotonCollection);
  int size = 0;
  for(size_t itPhoton = 0; itPhoton < candCollection.size(); itPhoton++){
    pat::PhotonRef candRef = candCollection[itPhoton];
    const pat::Photon* cand = dynamic_cast<const pat::Photon*> ((*candRef).clone());
    if(not selectionObj_(*cand)) continue;
    size++;
    if(produceOutputCollection_){
      filteredObjects->push_back(*cand);
    }
  }

  bool accept = false;
  if(size >= minNumber_ and size <= maxNumber_)
    accept = true;

  if(produceOutputCollection_)
    iEvent.put(std::move(filteredObjects));

  return accept;
}

void PhotonRefCountFilter::beginRun(edm::Run const& iRun , edm::EventSetup const& iSetup) {}
void PhotonRefCountFilter::endRun(edm::Run const&, edm::EventSetup const&) {}
void PhotonRefCountFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
void PhotonRefCountFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
void PhotonRefCountFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(PhotonRefCountFilter);


/////////////////////////////////////////////
class TauRefCountFilter : public edm::stream::EDFilter<> {

public:
  explicit TauRefCountFilter(const edm::ParameterSet&);
  ~TauRefCountFilter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual bool filter(edm::Event&, const edm::EventSetup&) ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  const bool filterEvents_;
  const bool produceOutputCollection_;
  edm::InputTag src_;
  int minNumber_;
  int maxNumber_;
  std::string selection_ ;
  edm::EDGetTokenT<pat::TauRefVector>  srcToken_;
  StringCutObjectSelector<pat::Tau>  selectionObj_;
};

TauRefCountFilter::TauRefCountFilter(const edm::ParameterSet& iConfig):
  filterEvents_(iConfig.existsAs<bool>("filterEvents") ? iConfig.getParameter<bool>("filterEvents") : true),
  produceOutputCollection_(iConfig.existsAs<bool>("produceOutputCollection") ? iConfig.getParameter<bool>("produceOutputCollection") : true),
  src_(iConfig.getParameter<edm::InputTag>("src")),
  minNumber_(iConfig.getParameter<int>("minNumber")),
  maxNumber_(iConfig.getParameter<int>("maxNumber")),
  selection_(iConfig.existsAs<std::string>("selection") ? iConfig.getParameter<std::string>("selection") : "abs(eta) < 5.0"),
  selectionObj_(StringCutObjectSelector<pat::Tau>(selection_))
{
  srcToken_    = consumes<pat::TauRefVector> (src_);
  if(produceOutputCollection_)
    produces<pat::TauCollection>();

}

TauRefCountFilter::~TauRefCountFilter() {}

bool TauRefCountFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // if you don't want to filter just exit
  if(filterEvents_ == false) return true;

  edm::Handle<pat::TauRefVector> candCollectionH;
  iEvent.getByToken(srcToken_,candCollectionH);
  pat::TauRefVector candCollection = *candCollectionH;

  std::unique_ptr<pat::TauCollection> filteredObjects (new pat::TauCollection);

  int size = 0;
  for(size_t itTau = 0; itTau < candCollection.size(); itTau++){
    pat::TauRef candRef = candCollection[itTau];
    const pat::Tau* cand = dynamic_cast<const pat::Tau*> ((*candRef).clone());
    if(not selectionObj_(*cand)) continue;
    size++;
    if(produceOutputCollection_){
      filteredObjects->push_back(*cand);
    }
  }

  bool accept = false;
  if(size >= minNumber_ and size <= maxNumber_)
    accept = true;

  if(produceOutputCollection_)
    iEvent.put(std::move(filteredObjects));

  return accept;
}

void TauRefCountFilter::beginRun(edm::Run const& iRun , edm::EventSetup const& iSetup) {}
void TauRefCountFilter::endRun(edm::Run const&, edm::EventSetup const&) {}
void TauRefCountFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
void TauRefCountFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
void TauRefCountFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(TauRefCountFilter);


/////////////////////////////////////////////
class JetRefCountFilter : public edm::stream::EDFilter<> {

public:
  explicit JetRefCountFilter(const edm::ParameterSet&);
  ~JetRefCountFilter();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  virtual bool filter(edm::Event&, const edm::EventSetup&) ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  const bool filterEvents_;
  const bool produceOutputCollection_;
  edm::InputTag src_;
  int minNumber_;
  int maxNumber_;
  std::string selection_ ;
  edm::EDGetTokenT<pat::JetRefVector>  srcToken_;
  StringCutObjectSelector<pat::Jet>  selectionObj_;
};

JetRefCountFilter::JetRefCountFilter(const edm::ParameterSet& iConfig):
  filterEvents_(iConfig.existsAs<bool>("filterEvents") ? iConfig.getParameter<bool>("filterEvents") : true),
  produceOutputCollection_(iConfig.existsAs<bool>("produceOutputCollection") ? iConfig.getParameter<bool>("produceOutputCollection") : true),
  src_(iConfig.getParameter<edm::InputTag>("src")),
  minNumber_(iConfig.getParameter<int>("minNumber")),
  maxNumber_(iConfig.getParameter<int>("maxNumber")),
  selection_(iConfig.existsAs<std::string>("selection") ? iConfig.getParameter<std::string>("selection") : "abs(eta) < 5.0"),
  selectionObj_(StringCutObjectSelector<pat::Jet>(selection_))
{
  srcToken_    = consumes<pat::JetRefVector> (src_);
  if(produceOutputCollection_)
    produces<pat::JetCollection>();

}

JetRefCountFilter::~JetRefCountFilter() {}

bool JetRefCountFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // if you don't want to filter just exit
  if(filterEvents_ == false) return true;

  edm::Handle<pat::JetRefVector> candCollectionH;
  iEvent.getByToken(srcToken_,candCollectionH);
  pat::JetRefVector candCollection = *candCollectionH;

  std::unique_ptr<pat::JetCollection> filteredObjects (new pat::JetCollection);

  int size = 0;
  for(size_t itJet = 0; itJet < candCollection.size(); itJet++){
    pat::JetRef candRef = candCollection[itJet];
    const pat::Jet* cand = dynamic_cast<const pat::Jet*> ((*candRef).clone());
    if(not selectionObj_(*cand)) continue;
    size++;
    if(produceOutputCollection_){
      filteredObjects->push_back(*cand);
    }
  }

  bool accept = false;
  if(size >= minNumber_ and size <= maxNumber_)
    accept = true;

  if(produceOutputCollection_)
    iEvent.put(std::move(filteredObjects));

  return accept;
}

void JetRefCountFilter::beginRun(edm::Run const& iRun , edm::EventSetup const& iSetup) {}
void JetRefCountFilter::endRun(edm::Run const&, edm::EventSetup const&) {}
void JetRefCountFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
void JetRefCountFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}
void JetRefCountFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(JetRefCountFilter);
