import os, sys
from multiprocessing import Process
from WMCore.Configuration import Configuration

#### Basic crab config structure

config = Configuration()

pyCfgParams = ['isMC=True',
               'filterOnHLT=False',
               'setHLTFilterFlag=True',
               'filterHighMETEvents=True',
               'metCut=190',
               'applyL2L3Residuals=False',
               'addQGLikelihood=True',
               'addPileupJetID=False',
               'addPuppiJets=True',
               'addPuppiMET=True',
               'addEGMSmear=True',
               'addMETSystematics=True',
               'addPuppiMETSystematics=False',
               'useOfficialMETSystematics=True',
               'addMETBreakDown=False',
               'addSubstructureCHS=False',
               'addSubstructurePuppi=False',
               'useMiniAODSubstructure=True',
               'miniAODProcess=PAT',
               'outputFileName=tree.root',
               'globalTag=80X_mcRun2_asymptotic_2016_TrancheIV_v8',
               'nThreads=4',
               'isCrab=True']

config.section_('General')
config.General.transferLogs = False
config.General.workArea     = 'crab_projects_MC_80X'  # Make sure you set this parameter

config.section_('JobType')
config.JobType.psetName         = '../tree.py'
config.JobType.pluginName       = 'Analysis'
config.JobType.outputFiles      = ['tree.root']
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB      = 2480
config.JobType.numCores         = 4


config.section_('Data')    
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'EventAwareLumiBased'
config.Data.unitsPerJob   = 25000
config.Data.outLFNDirBase = '/store/group/phys_exotica/monojet/rgerosa/ProductionMC_03_06_2017/'
config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

config.JobType.pyCfgParams = list(pyCfgParams)
