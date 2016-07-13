import os
import FWCore.ParameterSet.Config as cms

# Select good primary vertices
goodVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True))

## Bad muon filter
from RecoMET.METFilters.BadPFMuonFilter_cfi import BadPFMuonFilter
BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

##Bad charged tracks
from RecoMET.METFilters.BadChargedCandidateFilter_cfi import BadChargedCandidateFilter
BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

metFilters = cms.Sequence(goodVertices*BadPFMuonFilter*BadChargedCandidateFilter)
