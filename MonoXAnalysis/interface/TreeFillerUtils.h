#ifndef AnalysisCode_MonoXAnalysis_TreeFillerUtils_h
#define AnalysisCode_MonoXAnalysis_TreeFillerUtils_h

// basic C++ headers
#include <memory>
#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

// to apply jet ID                                                                                                                                                                                    
bool applyJetID(const pat::Jet &, const std::string &);
// to apply pileup-jet id                                                                                                                                                                             
bool applyPileupJetID(const pat::Jet &, const std::string &, const bool &);


#endif
