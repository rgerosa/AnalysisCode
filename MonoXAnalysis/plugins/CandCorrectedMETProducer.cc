#include <memory>
#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"

template< class T>
class CandCorrectedMETProducerT : public edm::stream::EDProducer<> {

public:
  explicit CandCorrectedMETProducerT(const edm::ParameterSet&);
  ~CandCorrectedMETProducerT();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void produce(edm::Event&, const edm::EventSetup&); 
  virtual void beginJob();
  virtual void endJob();

  reco::Candidate::LorentzVector findParticle(const T & particle, const edm::View<reco::Candidate> & pfCandCollection);
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  
  const edm::InputTag    metTag;
  std::vector<edm::InputTag> candTags;
  const edm::InputTag    pfCandidatesTag;
  const bool isPuppiTag;
  const bool useuncorrmet;
  
  edm::EDGetTokenT<edm::View<pat::MET> > metToken;
  edm::EDGetTokenT<edm::View<reco::Candidate> > pfCandidateToken; 
  edm::EDGetTokenT<edm::View<T> > theCandToken;
  std::vector<edm::EDGetTokenT<edm::View<T> > > candTokens;
};

template< class T>
CandCorrectedMETProducerT<T>::CandCorrectedMETProducerT(const edm::ParameterSet& iConfig): 
  metTag(iConfig.getParameter<edm::InputTag>("met")),
  candTags(iConfig.getParameter<std::vector<edm::InputTag> >("cands")),
  pfCandidatesTag(iConfig.existsAs<edm::InputTag>("pfCandidates") ? iConfig.getParameter<edm::InputTag>("pfCandidates") : edm::InputTag("packedPFCandidates")), 
  isPuppiTag(iConfig.existsAs<bool>("isPuppiTag") ? iConfig.getParameter<bool>("isPuppi") : false),
  useuncorrmet(iConfig.existsAs<bool>("useuncorrmet") ? iConfig.getParameter<bool>("useuncorrmet") : false){

  produces<pat::METCollection>();
  
  metToken = consumes<edm::View<pat::MET> > (metTag);
  pfCandidateToken = consumes<edm::View<reco::Candidate> >(pfCandidatesTag);
  
  for (std::size_t i = 0; i < candTags.size(); i++) {
    theCandToken = consumes<edm::View<T> > (candTags[i]);
    candTokens.push_back(theCandToken);
  }
}

template< class T>
CandCorrectedMETProducerT<T>::~CandCorrectedMETProducerT() {
  candTags.clear();
  candTokens.clear();
}

template< class T>
void CandCorrectedMETProducerT<T>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    Handle<View<pat::MET> > metH;
    iEvent.getByToken(metToken, metH);

    Handle<View<reco::Candidate> > pfCandH;
    iEvent.getByToken(pfCandidateToken, pfCandH);

    std::vector<Handle<View<T> > > candHs;
    for (std::size_t i = 0; i < candTokens.size(); i++) {
        Handle<View<T> > candH;
        iEvent.getByToken(candTokens[i], candH);
        candHs.push_back(candH);
    }

    std::auto_ptr<pat::METCollection> output(new pat::METCollection);

    double met    = (useuncorrmet ? metH->front().uncorPt()  : metH->front().corPt());
    double metphi = (useuncorrmet ? metH->front().uncorPhi() : metH->front().corPhi());
    
    double ccmetx = met * cos(metphi);
    double ccmety = met * sin(metphi);

    for (size_t i = 0; i < candHs.size(); i++) {
        for (auto cands_iter = candHs[i]->begin(); cands_iter != candHs[i]->end(); ++cands_iter) {
	  reco::Candidate::LorentzVector total4V;
	  if(isPuppiTag)
	    total4V = findParticle(*cands_iter,*pfCandH);
	  else 
	    total4V = cands_iter->p4();
	  ccmetx += total4V.pt() * cos(total4V.phi());
	  ccmety += total4V.pt() * sin(total4V.phi());
        }            
    }
    double ccmet = sqrt(ccmetx*ccmetx + ccmety*ccmety);

    pat::MET* ccmetcand = metH->front().clone();
    ccmetcand->setP4(reco::Candidate::LorentzVector(ccmetx, ccmety, 0., ccmet));
    output->push_back(*ccmetcand);

    iEvent.put(output);

}

template< class T>
reco::Candidate::LorentzVector CandCorrectedMETProducerT<T>::findParticle(const T & particle, const edm::View<reco::Candidate> & pfCandCollection){

  reco::Candidate::LorentzVector total4V;
  std::vector<reco::CandidatePtr> particles;
  for(size_t ipart = 0 ; ipart < particle.numberOfSourceCandidatePtrs(); ipart++){
    if(particle.sourceCandidatePtr(ipart).isNonnull() and particle.sourceCandidatePtr(ipart).isAvailable()){
      particles.push_back(particle.sourceCandidatePtr(ipart));
    }
  }

  if(particles.size() == 0)
    return particle.p4();

  for(unsigned int icand = 0; icand <  pfCandCollection.size(); icand++){
    
    reco::CandidatePtr ptrCand = pfCandCollection.ptrAt(icand);
    const pat::PackedCandidate* lPack = dynamic_cast<const pat::PackedCandidate *>(&(pfCandCollection.at(icand)));
    if(lPack->puppiWeightNoLep() == 0) // in case of puppi, only the one used for MET calculations are useful                                                                   
      continue;

    for(auto ipart : particles){
      if(ipart == ptrCand){
        total4V += ptrCand->p4()*lPack->puppiWeightNoLep();
        break;
      }
    }
  }

  particles.clear();
  return total4V;
}

template< class T>
void CandCorrectedMETProducerT<T>::beginJob() {
}

template< class T>
void CandCorrectedMETProducerT<T>::endJob() {
}

template< class T>
void CandCorrectedMETProducerT<T>::beginRun(edm::Run const&, edm::EventSetup const&) {
}

template< class T>
void CandCorrectedMETProducerT<T>::endRun(edm::Run const&, edm::EventSetup const&) {
}

template< class T>
void CandCorrectedMETProducerT<T>::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

template< class T>
void CandCorrectedMETProducerT<T>::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

template< class T>
void CandCorrectedMETProducerT<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

// specialization
// specialization for muons
template<>
void CandCorrectedMETProducerT<pat::Muon>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    using namespace reco;
    using namespace std;

    Handle<View<pat::MET> > metH;
    iEvent.getByToken(metToken, metH);

    Handle<View<reco::Candidate> > pfCandH;
    iEvent.getByToken(pfCandidateToken, pfCandH);

    std::vector<Handle<View<pat::Muon> > > candHs;
    for (std::size_t i = 0; i < candTokens.size(); i++) {
      Handle<View<pat::Muon> > candH;
      iEvent.getByToken(candTokens[i], candH);
      candHs.push_back(candH);
    }

    std::auto_ptr<pat::METCollection> output(new pat::METCollection);

    double met    = (useuncorrmet ? metH->front().uncorPt()  : metH->front().corPt());
    double metphi = (useuncorrmet ? metH->front().uncorPhi() : metH->front().corPhi());
    
    double ccmetx = met * cos(metphi);
    double ccmety = met * sin(metphi);

    for (size_t i = 0; i < candHs.size(); i++) {
      for (auto cands_iter = candHs[i]->begin(); cands_iter != candHs[i]->end(); ++cands_iter) {

	  reco::Candidate::LorentzVector total4V;
	  if(isPuppiTag)
	    total4V = findParticle(*cands_iter,*pfCandH);
	  else 
	    total4V = cands_iter->pfP4();

	  ccmetx += total4V.pt() * cos(total4V.phi());
	  ccmety += total4V.pt() * sin(total4V.phi());
      }            
    }
    double ccmet = sqrt(ccmetx*ccmetx + ccmety*ccmety);

    pat::MET* ccmetcand = metH->front().clone();
    ccmetcand->setP4(reco::Candidate::LorentzVector(ccmetx, ccmety, 0., ccmet));
    output->push_back(*ccmetcand);

    iEvent.put(output);
}




typedef CandCorrectedMETProducerT<pat::Muon> MuonCorrectedMETProducer;
DEFINE_FWK_MODULE(MuonCorrectedMETProducer);

typedef CandCorrectedMETProducerT<pat::Electron> ElectronCorrectedMETProducer;
DEFINE_FWK_MODULE(ElectronCorrectedMETProducer);

typedef CandCorrectedMETProducerT<pat::Photon> PhotonCorrectedMETProducer;
DEFINE_FWK_MODULE(PhotonCorrectedMETProducer);

typedef CandCorrectedMETProducerT<pat::Tau> TauCorrectedMETProducer;
DEFINE_FWK_MODULE(TauCorrectedMETProducer);

typedef CandCorrectedMETProducerT<pat::Jet> JetCorrectedMETProducer;
DEFINE_FWK_MODULE(JetCorrectedMETProducer);
