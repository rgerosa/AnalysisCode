import os
import FWCore.ParameterSet.Config as cms
from CondCore.DBCommon.CondDBSetup_cfi import *

## setup JEC on PAT jets from miniAOD                                                                                                                                    
def JERConfiguration(process,usePrivateSQlite,JEREra,isMC):

    process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
    
    process.JetResolutionESProducer_AK8PFchs = process.JetResolutionESProducer_AK4PFchs.clone(
        label = cms.string('AK8PFchs'))

    process.JetResolutionESProducer_SF_AK8PFchs = process.JetResolutionESProducer_SF_AK4PFchs.clone(
        label = cms.string('AK8PFchs'))

    process.JetResolutionESProducer_AK4PFPuppi = process.JetResolutionESProducer_AK4PFchs.clone(
        label = cms.string('AK4PFPuppi'))

    process.JetResolutionESProducer_SF_AK4PFPuppi = process.JetResolutionESProducer_SF_AK4PFchs.clone(
        label = cms.string('AK4PFPuppi'))

    process.JetResolutionESProducer_AK8PFPuppi = process.JetResolutionESProducer_AK8PFchs.clone(
        label = cms.string('AK8PFPuppi'))

    process.JetResolutionESProducer_SF_AK8PFPuppi = process.JetResolutionESProducer_SF_AK8PFchs.clone(
        label = cms.string('AK8PFPuppi'))

    era = JEREra
    if isMC :
        era += "_MC"
    else :
        era += "_DATA"

    dBFile = os.path.expandvars(era+".db")
    
    
    if usePrivateSQlite:

        process.jer = cms.ESSource("PoolDBESSource",
                                   CondDBSetup,
                                   connect = cms.string("sqlite_file:../data/JER/"+dBFile),
                                   toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetResolutionRcd'),
                    tag    = cms.string('JR_'+era+'_PtResolution_AK4PFchs'),
                    label  = cms.untracked.string('AK4PFchs_pt')
                    ),
                cms.PSet(
                    record = cms.string('JetResolutionRcd'),
                    tag    = cms.string('JR_'+era+'_PtResolution_AK4PFchs'), ## to be fixed in future
                    label  = cms.untracked.string('AK4PFPuppi_pt')
                    ),
        
                cms.PSet(
                    record = cms.string('JetResolutionRcd'),
                    tag    = cms.string('JR_'+era+'_PtResolution_AK8PFchs'),
                    label  = cms.untracked.string('AK8PFchs_pt')
                    ),
                cms.PSet(
                    record = cms.string('JetResolutionRcd'),
                    tag    = cms.string('JR_'+era+'_PtResolution_AK8PFchs'), ## to be fixed in future
                    label  = cms.untracked.string('AK8PFPuppi_pt')
                    ),
                
                # Scale factors
                cms.PSet(
                    record = cms.string('JetResolutionScaleFactorRcd'),
                    tag    = cms.string('JR_'+era+'_SF_AK4PFchs'),
                    label  = cms.untracked.string('AK4PFchs')
                    ),

                cms.PSet(
                    record = cms.string('JetResolutionScaleFactorRcd'),
                    tag    = cms.string('JR_'+era+'_SF_AK4PFchs'), ## to be fixed in future
                    label  = cms.untracked.string('AK4PFPuppi')
                    ),

                cms.PSet(
                    record = cms.string('JetResolutionScaleFactorRcd'),
                    tag    = cms.string('JR_'+era+'_SF_AK4PFchs'),
                    label  = cms.untracked.string('AK8PFchs')
                    ),

               cms.PSet(
                    record = cms.string('JetResolutionScaleFactorRcd'),
                    tag    = cms.string('JR_'+era+'_SF_AK4PFchs'), ## to be fixed in future
                    label  = cms.untracked.string('AK8PFPuppi')
                    ),
                ),
                                   )
        
        process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')
        
