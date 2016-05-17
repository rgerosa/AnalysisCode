import os
import FWCore.ParameterSet.Config as cms

# Select good primary vertices
goodVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True))

metFilters = cms.Sequence(goodVertices)
