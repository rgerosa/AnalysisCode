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
               'addEGMSmear=False',
               'addMETSystematics=True',
               'useOfficialMETSystematics=True',
               'addMETBreakDown=True',
               'addSubstructureCHS=True',
               'addSubstructurePuppi=False',
               'miniAODProcess=PAT',
               'globalTag=80X_mcRun2_asymptotic_2016_miniAODv2_v1',
               'outputFileName=tree.root',
               'nThreads=3',
               'usePrivateSQliteJEC=True',
               'isCrab=True']

config.section_('General')
config.General.transferLogs = False
config.General.workArea     = 'crab_projects_MC_80X'  # Make sure you set this parameter

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
config.Data.unitsPerJob   = 30000
config.Data.outLFNDirBase = '/store/group/upgrade/delphes/VBS_SS/Production-30-05-2016-80X-MC_Inclusive/'
config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

config.JobType.pyCfgParams = list(pyCfgParams)
