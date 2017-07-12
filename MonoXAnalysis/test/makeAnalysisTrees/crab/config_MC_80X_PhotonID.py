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
               'usePrivateSQliteJEC=False',
               'applyL2L3Residuals=False',
               'addQGLikelihood=False',
               'addPileupJetID=False',
               'addPuppiJets=False',
               'addPuppiMET=False',
               'addEGMSmear=False',
               'addMETSystematics=False',
               'useOfficialMETSystematics=False',
               'addMETBreakDown=False',
               'addSubstructureCHS=False',
               'addSubstructurePuppi=False',
               'addPhotonPurity=False',
               'miniAODProcess=PAT',
               'globalTag=80X_mcRun2_asymptotic_2016_miniAODv2_v1',
               'outputFileName=tree.root',
               'addPhotonIDVariables=True',
               'addElectronIDVariables=True',
               'nThreads=3',
               'isCrab=True']

config.section_('General')
config.General.transferLogs = False
config.General.workArea     = 'crab_projects_MC_80X_PhotonID'  # Make sure you set this parameter

config.section_('JobType')
config.JobType.psetName         = '../tree.py'
config.JobType.pluginName       = 'Analysis'
config.JobType.outputFiles      = ['tree.root']
config.JobType.allowUndistributedCMSSW = True
config.JobType.numCores         = 3
config.JobType.maxMemoryMB      = 2500


config.section_('Data')    
config.Data.inputDBS      = 'global'
config.Data.splitting     = 'EventAwareLumiBased'
config.Data.unitsPerJob   = 40000
config.Data.outLFNDirBase = '/store/user/rgerosa/MONOJET_ANALYSIS/Production-28-08-2016_PhotonID/'
config.Data.publication   = False


config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
config.JobType.pyCfgParams = list(pyCfgParams)
