import os, sys
from multiprocessing import Process
from WMCore.Configuration import Configuration

#### Basic crab config structure

config = Configuration()

pyCfgParams = ['isMC=True',
               'filterOnHLT=False',
               'setHLTFilterFlag=True',
               'filterHighMETEvents=True',
               'metCut=0',
               'applyL2L3Residuals=False',
               'addQGLikelihood=True',
               'addPileupJetID=False',
               'addPuppiJets=True',
               'addPuppiMET=True',
               'addEGMSmear=True',
               'addMETSystematics=True',
               'useOfficialMETSystematics=True',
               'addMETBreakDown=True',
               'addSubstructureCHS=False',
               'addSubstructurePuppi=False',
               'isPhotonPurity=True',
               'miniAODProcess=PAT',
               'globalTag=80X_mcRun2_asymptotic_2016_miniAODv2_v1',
               'outputFileName=tree.root',
               'usePrivateSQliteJEC=True',       
               'nThreads=3',
               'isCrab=True']

config.section_('General')
config.General.transferLogs = False
config.General.workArea     = 'crab_projects_MC_80X_PhotonPurity'  # Make sure you set this parameter

config.section_('JobType')
config.JobType.psetName         = '../tree.py'
config.JobType.pluginName       = 'Analysis'
config.JobType.outputFiles      = ['tree.root']
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB      = 2450
config.JobType.numCores         = 3


config.section_('Data')    
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'EventAwareLumiBased'
config.Data.unitsPerJob   = 60000
config.Data.outLFNDirBase = '/store/user/rgerosa/MONOJET_ANALYSIS/ProductionMC-05-01-2017_PhotonPurity/'
config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.JobType.pyCfgParams = list(pyCfgParams)
