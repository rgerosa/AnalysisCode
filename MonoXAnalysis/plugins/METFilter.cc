#include <memory>
#include <vector>
#include <string>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"


template<class T>
class METFilter : public edm::stream::EDFilter<> {

public:
  explicit METFilter(const edm::ParameterSet&);
  ~METFilter();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  
  virtual bool filter(edm::Event&, const edm::EventSetup&) ;        
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override; 
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override; 
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  const std::vector<edm::ParameterSet> metCollectionInfo_;
  std::vector<edm::InputTag> metCollectionTag_;
  const bool filterEvents_;
  const bool graterThan_;
  const bool applyAndInsteadOfOr_;
  std::vector<double> metTreshold_;
  
  std::vector<edm::EDGetTokenT<T> > metCollectionToken_;
};

template<class T>
METFilter<T>::METFilter(const edm::ParameterSet& iConfig):
  metCollectionInfo_(iConfig.getParameter<std::vector<edm::ParameterSet> >("metCollections")),
  filterEvents_(iConfig.existsAs<bool>("filterEvents") ? iConfig.getParameter<bool>("filterEvents") : true),
  graterThan_(iConfig.existsAs<bool>("graterThan") ? iConfig.getParameter<bool>("graterThan") : true),
  applyAndInsteadOfOr_(iConfig.existsAs<bool>("applyAndInsteadOfOr") ? iConfig.getParameter<bool>("applyAndInsteadOfOr") : false)
{

  // loop on the met collection info to extract the collection tag and thresholds
  for(auto it = metCollectionInfo_.begin(); it != metCollectionInfo_.end(); ++it){   
    metCollectionTag_.push_back((*it).getParameter<edm::InputTag> ("srcMet"));
    metCollectionToken_.push_back(consumes<T>(metCollectionTag_.back()));
    metTreshold_.push_back((*it).getParameter<double> ("metCut"));
  }
}

template<class T>
METFilter<T>::~METFilter() {
  metCollectionTag_.clear();
  metCollectionToken_.clear();
  metTreshold_.clear();
}

template<class T>
bool METFilter<T>::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // if you don't want to filter just exit
  if(filterEvents_ == false) return true;

  bool eventAccepted = false;
  // get by token
  edm::Handle<T> metCollectionH;
  size_t iMetColl = 0;
  for(auto token : metCollectionToken_){ // loop on the met collection
    iEvent.getByToken(token,metCollectionH); // take the met
    
    for(auto met : *metCollectionH)
      if(graterThan_){
	if(not applyAndInsteadOfOr_ and met.pt() > metTreshold_.at(iMetColl))
	  eventAccepted += true;
	else if(applyAndInsteadOfOr_ and met.pt() > metTreshold_.at(iMetColl))
	  eventAccepted *= true;
      }
      else{
	if(not applyAndInsteadOfOr_ and met.pt() < metTreshold_.at(iMetColl))
	  eventAccepted += true;
	else if(applyAndInsteadOfOr_ and met.pt() < metTreshold_.at(iMetColl))
	  eventAccepted *= true;
      }
  }

  return eventAccepted;

}

template<class T>
void METFilter<T>::beginRun(edm::Run const& iRun , edm::EventSetup const& iSetup) {}

template<class T>
void METFilter<T>::endRun(edm::Run const&, edm::EventSetup const&) {}

template<class T>
void METFilter<T>::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

template<class T>
void METFilter<T>::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

template<class T>
void METFilter<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// specification
typedef METFilter<pat::METCollection> PATMETFilter;
DEFINE_FWK_MODULE(PATMETFilter);

typedef METFilter<reco::METCollection> RecoMETFilter;
DEFINE_FWK_MODULE(RecoMETFilter);
