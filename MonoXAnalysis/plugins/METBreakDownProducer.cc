#include <memory>
#include <vector>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/LorentzVector.h"

template <class T>
class METBreakDownProducer : public edm::stream::EDProducer<> {
    public:
        explicit METBreakDownProducer(const edm::ParameterSet&);
  ~METBreakDownProducer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
        virtual void beginJob() ;
        virtual void produce(edm::Event&, const edm::EventSetup&) ;
        virtual void endJob() ;

        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

        const edm::InputTag candsTag;
        edm::EDGetTokenT<edm::View<reco::Candidate> > candsToken;
};

template <class T>
METBreakDownProducer<T>::METBreakDownProducer(const edm::ParameterSet& iConfig):
  candsTag(iConfig.getParameter<edm::InputTag>("pfcands")){
  produces<std::vector<T> >("pfMetHadronHF");
  produces<std::vector<T> >("pfMetEgammaHF");
  produces<std::vector<T> >("pfMetChargedHadron");
  produces<std::vector<T> >("pfMetNeutralHadron");
  produces<std::vector<T> >("pfMetElectrons");
  produces<std::vector<T> >("pfMetMuons");
  produces<std::vector<T> >("pfMetPhotons");
  produces<std::vector<T> >("pfMetUnclustered");
  produces<std::vector<T> >("pfMet");
  candsToken = consumes<edm::View<reco::Candidate> >(candsTag);
}

template <class T>
METBreakDownProducer<T>::~METBreakDownProducer() {
}

template <class T>
void METBreakDownProducer<T>::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace reco;
  using namespace std;

  Handle<View<Candidate> > candsH;
  iEvent.getByToken(candsToken, candsH);

  std::unique_ptr<std::vector<T> > pfMet(new std::vector<T>);
  std::unique_ptr<std::vector<T> > pfMetHadronHF(new std::vector<T>);
  std::unique_ptr<std::vector<T> > pfMetEgammaHF(new std::vector<T>);
  std::unique_ptr<std::vector<T> > pfMetChargedHadron(new std::vector<T>);
  std::unique_ptr<std::vector<T> > pfMetNeutralHadron(new std::vector<T>);
  std::unique_ptr<std::vector<T> > pfMetElectrons(new std::vector<T>);
  std::unique_ptr<std::vector<T> > pfMetMuons(new std::vector<T>);
  std::unique_ptr<std::vector<T> > pfMetPhotons(new std::vector<T>);
  std::unique_ptr<std::vector<T> > pfMetUnclustered(new std::vector<T>);

  // full met
  double metx  = 0.0, mety  = 0.0;
  // neutral hadron met
  double hmetx = 0.0, hmety = 0.0;
  // charged hadron met
  double cmetx = 0.0, cmety = 0.0;
  // e-gamma forward region
  double ametx = 0.0, amety = 0.0;
  // hadron forward region
  double bmetx = 0.0, bmety = 0.0;
  // electrons
  double emetx = 0.0, emety = 0.0;
  // muons
  double mmetx = 0.0, mmety = 0.0;
  // photons
  double pmetx = 0.0, pmety = 0.0;
  // others
  double ometx = 0.0, omety = 0.0;

  for (View<Candidate>::const_iterator cands_iter = candsH->begin(); cands_iter != candsH->end(); ++cands_iter) {
    int absid = abs(cands_iter->pdgId()); // read info from  http://cmslxr.fnal.gov/source/DataFormats/ParticleFlowCandidate/src/PFCandidate.cc#0224
    if (absid == 211) { // charged hadrons in the tracker covered region
      cmetx -= cands_iter->pt() * cos(cands_iter->phi());
      cmety -= cands_iter->pt() * sin(cands_iter->phi());
    }
    else if (absid == 130) { // neutral hadrons in the tracker covered region -->
      hmetx -= cands_iter->pt() * cos(cands_iter->phi());
      hmety -= cands_iter->pt() * sin(cands_iter->phi());
    }
    else if (absid == 1) { // hadrons in the HF region
      ametx -= cands_iter->pt() * cos(cands_iter->phi());
      amety -= cands_iter->pt() * sin(cands_iter->phi());
    }
    else if (absid == 2) { // egamma in the HF region
      bmetx -= cands_iter->pt() * cos(cands_iter->phi());
      bmety -= cands_iter->pt() * sin(cands_iter->phi());
    }
    else if (absid == 11) { // electrons
      emetx -= cands_iter->pt() * cos(cands_iter->phi());
      emety -= cands_iter->pt() * sin(cands_iter->phi());
    }
    else if (absid == 13) { // muons
      mmetx -= cands_iter->pt() * cos(cands_iter->phi());
      mmety -= cands_iter->pt() * sin(cands_iter->phi());
    }
    else if (absid == 22) { // photons
      pmetx -= cands_iter->pt() * cos(cands_iter->phi());
      pmety -= cands_iter->pt() * sin(cands_iter->phi());
    }
    else {
      ometx -= cands_iter->pt() * cos(cands_iter->phi());
      omety -= cands_iter->pt() * sin(cands_iter->phi());
    }
    metx -= cands_iter->pt() * cos(cands_iter->phi());
    mety -= cands_iter->pt() * sin(cands_iter->phi());
  }
  // calculate met values
  double  met = sqrt( metx* metx +  mety* mety);
  double hmet = sqrt(hmetx*hmetx + hmety*hmety);
  double amet = sqrt(ametx*ametx + amety*amety);
  double bmet = sqrt(bmetx*bmetx + bmety*bmety);
  double cmet = sqrt(cmetx*cmetx + cmety*cmety);
  double emet = sqrt(emetx*emetx + emety*emety);
  double mmet = sqrt(mmetx*mmetx + mmety*mmety);
  double pmet = sqrt(pmetx*pmetx + pmety*pmety);
  double omet = sqrt(ometx*ometx + omety*omety);

  T metcand; metcand.setP4(reco::Candidate::LorentzVector( metx,  mety, 0., met));
  pfMet->push_back(metcand);

  T cmetcand; cmetcand.setP4(reco::Candidate::LorentzVector( cmetx,  cmety, 0., cmet));
  pfMetChargedHadron->push_back(cmetcand);

  T hmetcand; hmetcand.setP4(reco::Candidate::LorentzVector( hmetx,  hmety, 0., hmet));
  pfMetNeutralHadron->push_back(hmetcand);

  T ametcand; ametcand.setP4(reco::Candidate::LorentzVector( ametx,  amety, 0., amet));
  pfMetHadronHF->push_back(ametcand);

  T bmetcand; bmetcand.setP4(reco::Candidate::LorentzVector( bmetx,  bmety, 0., bmet));
  pfMetEgammaHF->push_back(bmetcand);

  T emetcand; emetcand.setP4(reco::Candidate::LorentzVector( emetx,  emety, 0., emet));
  pfMetElectrons->push_back(emetcand);

  T mmetcand; mmetcand.setP4(reco::Candidate::LorentzVector( mmetx,  mmety, 0., mmet));
  pfMetMuons->push_back(mmetcand);

  T pmetcand; pmetcand.setP4(reco::Candidate::LorentzVector( pmetx,  pmety, 0., pmet));
  pfMetPhotons->push_back(pmetcand);

  T ometcand; ometcand.setP4(reco::Candidate::LorentzVector( ometx,  omety, 0., omet));
  pfMetUnclustered->push_back(ometcand);

  iEvent.put(std::move(pfMet),"pfMet");
  iEvent.put(std::move(pfMetChargedHadron),"pfMetChargedHadron");
  iEvent.put(std::move(pfMetNeutralHadron),"pfMetNeutralHadron");
  iEvent.put(std::move(pfMetHadronHF),"pfMetHadronHF");
  iEvent.put(std::move(pfMetEgammaHF),"pfMetEgammaHF");
  iEvent.put(std::move(pfMetElectrons),"pfMetElectrons");
  iEvent.put(std::move(pfMetMuons),"pfMetMuons");
  iEvent.put(std::move(pfMetPhotons),"pfMetPhotons");
  iEvent.put(std::move(pfMetUnclustered),"pfMetUnclustered");

}

template<class T>
void METBreakDownProducer<T>::beginJob() {}

template<class T>
void METBreakDownProducer<T>::endJob() {}

template<class T>
void METBreakDownProducer<T>::beginRun(edm::Run const&, edm::EventSetup const&) {}

template<class T>
void METBreakDownProducer<T>::endRun(edm::Run const&, edm::EventSetup const&) {}

template<class T>
void METBreakDownProducer<T>::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

template<class T>
void METBreakDownProducer<T>::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {}

template<class T>
void METBreakDownProducer<T>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

typedef METBreakDownProducer<reco::MET> RecoMETBreakDownProducer;
DEFINE_FWK_MODULE(RecoMETBreakDownProducer);

typedef METBreakDownProducer<pat::MET> PATMETBreakDownProducer;
DEFINE_FWK_MODULE(PATMETBreakDownProducer);
