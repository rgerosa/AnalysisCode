import os
import FWCore.ParameterSet.Config as cms

# Re-run the HBHE Noise filter
from CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi import *
HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False) 
HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")

# Select good primary vertices
goodVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True))

metFilters = cms.Sequence(#HBHENoiseFilterResultProducer+
                          goodVertices)
