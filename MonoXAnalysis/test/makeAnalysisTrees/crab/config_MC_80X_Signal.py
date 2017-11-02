import os, sys
from multiprocessing import Process
from WMCore.Configuration import Configuration

#### Basic crab config structure

config = Configuration()

pyCfgParams = ['isMC=True',
               'filterOnHLT=False',
               'setHLTFilterFlag=True',
               'filterHighMETEvents=False',
               'metCut=0',
               'applyL2L3Residuals=False',
               'addQGLikelihood=False',
               'addPileupJetID=False',
               'addPuppiJets=False',
               'addPuppiMET=False',
               'addEGMSmear=True',
               'addMETSystematics=True',
               'useOfficialMETSystematics=True',
               'addPuppiMETSystematics=False',
               'addMETBreakDown=False',
               'addSubstructureCHS=True',
               'useMiniAODSubstructure=False',
               'addSubstructurePuppi=False',
               'miniAODProcess=PAT',
               'globalTag=80X_mcRun2_asymptotic_2016_TrancheIV_v8',
               'outputFileName=tree.root',
               'nThreads=4',
               'isCrab=True']

config.section_('General')
config.General.transferLogs = False
config.General.workArea     = 'crab_projects_MC_80X_Signal'  # Make sure you set this parameter

config.section_('JobType')
config.JobType.psetName         = '../tree.py'
config.JobType.pluginName       = 'Analysis'
config.JobType.outputFiles      = ['tree.root']
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB      = 2450
config.JobType.numCores         = 4


config.section_('Data')    
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'EventAwareLumiBased'
config.Data.unitsPerJob   = 10000
config.Data.outLFNDirBase = '/store/group/phys_exotica/monojet/rgerosa/ProductionMC_Signal_MonoJ_LL'
config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

config.JobType.pyCfgParams = list(pyCfgParams)
