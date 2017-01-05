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

template<class T>
class CandRefCountFilterT : public edm::stream::EDFilter<> {

public:
  explicit CandRefCountFilterT(const edm::ParameterSet&);
  ~CandRefCountFilterT();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  
  virtual bool filter(edm::Event&, const edm::EventSetup&) ;        
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override; 
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override; 
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  std::string   selection_;
  const bool filterEvents_;  
  edm::InputTag src_;
  int minNumber_;
  int maxNumber_;
  edm::EDGetTokenT<T>  srcToken_;
};

template<class T>
CandRefCountFilterT<T>::CandRefCountFilterT(const edm::ParameterSet& iConfig):
  filterEvents_(iConfig.existsAs<bool>("filterEvents") ? iConfig.getParameter<bool>("filterEvents") : true),
  src_(iConfig.getParameter<edm::InputTag>("src")),
  minNumber_(iConfig.getParameter<int>("minNumber")),
  maxNumber_(iConfig.getParameter<int>("maxNumber"))  
{
  srcToken_ = consumes<T>(src_);
}

template<class T>
CandRefCountFilterT<T>::~CandRefCountFilterT() {}

template<class T>
bool CandRefCountFilterT<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  // if you don't want to filter just exit
  if(filterEvents_ == false) return true;

  edm::Handle<T> candCollectionH;
  iEvent.getByToken(srcToken_,candCollectionH);

  int size = 0;
  for(typename T::const_iterator itCand = candCollectionH->begin(); itCand != candCollectionH->end(); itCand++){
    size++;
  }
  
  bool accept = false;
  if(size >= minNumber_ and size <= maxNumber_)
    accept = true;
  
  return accept;
}

template<class T>
void CandRefCountFilterT<T>::beginRun(edm::Run const& iRun , edm::EventSetup const& iSetup) {}

template<class T>
void CandRefCountFilterT<T>::endRun(edm::Run const&, edm::EventSetup const&) {}

template<class T>
void CandRefCountFilterT<T>::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

template<class T>
void CandRefCountFilterT<T>::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

template<class T>
void CandRefCountFilterT<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

// specification
typedef CandRefCountFilterT<pat::MuonRefVector> MuonRefCountFilter;
DEFINE_FWK_MODULE(MuonRefCountFilter);
typedef CandRefCountFilterT<pat::ElectronRefVector> ElectronRefCountFilter;
DEFINE_FWK_MODULE(ElectronRefCountFilter);
typedef CandRefCountFilterT<pat::TauRefVector> TauRefCountFilter;
DEFINE_FWK_MODULE(TauRefCountFilter);
typedef CandRefCountFilterT<pat::PhotonRefVector> PhotonRefCountFilter;
DEFINE_FWK_MODULE(PhotonRefCountFilter);
typedef CandRefCountFilterT<pat::JetRefVector> JetRefCountFilter;
DEFINE_FWK_MODULE(JetRefCountFilter);
