#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"


template < class T>
class SelectorJetConstituent : public edm::stream::EDProducer<> {
public:

  typedef std::vector<T> JetsOutput;

  SelectorJetConstituent ( edm::ParameterSet const & params ) :
    srcToken_( consumes< typename edm::View<T> >( params.getParameter<edm::InputTag>("src") ) ),
    cut_( params.getParameter<std::string>("cut") ),
    selector_( cut_ )
  {
   produces< JetsOutput >();
   produces< edm::PtrVector<reco::Candidate> > ("constituents");
  }
 
  virtual ~SelectorJetConstituent() {}
  
  virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override {
    
    std::auto_ptr< JetsOutput > jets ( new std::vector<T>() );
    std::auto_ptr< edm::PtrVector<reco::Candidate> > candsOut( new edm::PtrVector<reco::Candidate>  );
    
    edm::Handle< typename edm::View<T> > h_jets;
    iEvent.getByToken( srcToken_, h_jets );
    // Now set the Ptrs with the orphan handles.
    for ( typename edm::View<T>::const_iterator ibegin = h_jets->begin(),
	    iend = h_jets->end(), ijet = ibegin;
	  ijet != iend; ++ijet ) {
     // Check the selection
      if ( selector_(*ijet) ) {
	// Add the jets that pass to the output collection
	jets->push_back( *ijet );
	for ( unsigned int ida = 0; ida < ijet->numberOfDaughters(); ++ida ) {
	  candsOut->push_back((*ijet).daughterPtr(ida));
	}
      }
    }
    
    // put  in Event
    iEvent.put(jets);
    iEvent.put(candsOut, "constituents");
    
  }

protected:
  edm::EDGetTokenT< typename edm::View<T> > srcToken_;
  std::string                cut_;
  StringCutObjectSelector<T> selector_;
 
 };

typedef SelectorJetConstituent<pat::Jet>       PatSelectorJetConstituent;
typedef SelectorJetConstituent<reco::PFJet>    PFSelectorJetConstituent;
typedef SelectorJetConstituent<reco::GenJet>   GenSelectorJetConstituent;

DEFINE_FWK_MODULE( PatSelectorJetConstituent );
DEFINE_FWK_MODULE( PFSelectorJetConstituent );
DEFINE_FWK_MODULE( GenSelectorJetConstituent );
