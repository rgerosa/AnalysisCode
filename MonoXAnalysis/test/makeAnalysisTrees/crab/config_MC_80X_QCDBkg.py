import os, sys
from multiprocessing import Process
from WMCore.Configuration import Configuration
from CRABClient.UserUtilities import config, getUsernameFromSiteDB

#### Basic crab config structure
config = Configuration()

pyCfgParams = ['isMC=True',
               'filterOnHLT=False', ## don't filter event according to the trigger bit
               'setHLTFilterFlag=False',
               'filterHighMETEvents=True', ## apply a low recoil cut at 50 GeV
               'metCut=100',                ## 50 GeV
               'applyL2L3Residuals=False', 
               'addQGLikelihood=False', 
               'addPileupJetID=False',
               'addPuppiJets=True',
               'addPuppiMET=True',
               'addEGMSmear=True',
               'addMETSystematics=True',
               'addPuppiMETSystematics=False',
               'useOfficialMETSystematics=True',
               'addMETBreakDown=False',
               'addSubstructureCHS=False',
               'addSubstructurePuppi=False', ## to spead up
               'miniAODProcess=PAT',
               'globalTag=80X_mcRun2_asymptotic_2016_TrancheIV_v8',
               'outputFileName=tree.root',
               'usePrivateSQliteJEC=False',
               'isQCDTree=True',
               'nThreads=4',
               'isCrab=True']

config.section_('General')
config.General.transferLogs = False
config.General.workArea     = 'crab_projects_MC_80X_QCDBkg'  # Make sure you set this parameter
config.General.requestName  = ''

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
config.Data.outLFNDirBase = '/store/group/phys_exotica/monojet/rgerosa/ProductionMC_21_05_2017_QCDBkg'
config.Data.allowNonValidInputDataset = True

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

config.JobType.pyCfgParams = list(pyCfgParams)
 
 
