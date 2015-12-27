#include "DataFormats/PatCandidates/interface/Jet.h"
#include "CommonTools/RecoAlgos/plugins/JetDeltaRValueMapProducer.cc"

typedef JetDeltaRValueMapProducer<pat::Jet> PATJetDeltaRValueMapProducer;
DEFINE_FWK_MODULE( PATJetDeltaRValueMapProducer);
