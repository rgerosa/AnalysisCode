import os, sys
from multiprocessing import Process
from WMCore.Configuration import Configuration

#### Basic crab config structure

config = Configuration()

pyCfgParams = ['isMC=True',
               'filterOnHLT=False',
               'filterHighMETEvents=True',
               'metCut=80',
               'applyL2L3Residuals=False',
               'addQGLikelihood=False',
               'addPileupJetID=False',
               'addPuppiJets=False',
               'addPuppiMET=False',
               'addEGMSmear=False',
               'addMVAMet=False',
               'addMETSystematics=False',
               'useOfficialMETSystematics=False',
               'addMETBreakDown=False',
               'addSubstructureCHS=False',
               'addSubstructurePuppi=False',
               'miniAODProcess=PAT'
               'outputFileName=tree.root',
               'isTriggerTree=True',
               'addTriggerObjects=True',
               'nThreads=4',
               'isCrab=True']

config.section_('General')
config.General.transferLogs = False
config.General.workArea     = 'crab_projects_MC_80X_trigger'  # Make sure you set this parameter

config.section_('JobType')
config.JobType.psetName         = '../tree.py'
config.JobType.pluginName       = 'Analysis'
config.JobType.outputFiles      = ['tree.root']
config.JobType.allowUndistributedCMSSW = True
config.JobType.numCores         = 4
config.JobType.maxMemoryMB      = 2500


config.section_('Data')    
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'EventAwareLumiBased'
config.Data.unitsPerJob   = 50000
config.Data.outLFNDirBase = '/store/group/phys_exotica/monojet/rgerosa/ProductionMC_13_03_2017_trigger/'
config.Data.allowNonValidInputDataset = True
config.Data.publication   = False


config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

config.JobType.pyCfgParams = list(pyCfgParams)
