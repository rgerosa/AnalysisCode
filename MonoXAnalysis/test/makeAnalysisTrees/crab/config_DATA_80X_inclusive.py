import os, sys
from multiprocessing import Process
from WMCore.Configuration import Configuration

#### Basic crab config structure

config = Configuration()

pyCfgParams = ['isMC=False',
               'filterOnHLT=True',
               'filterHighMETEvents=False',
               'applyL2L3Residuals=True',
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
               'miniAODProcess=RECO',
               'globalTag=80X_dataRun2_Prompt_v8',
               'outputFileName=tree.root',
               'nThreads=3',
               'isCrab=True']

config.section_('General')
config.General.transferLogs = False
config.General.workArea     = 'crab_projects_DATA_80X'  # Make sure you set this parameter

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
config.Data.unitsPerJob   = 30000
config.Data.outLFNDirBase = '/store/group/upgrade/delphes/VBS_SS/Production-06-06-2016_80X_Data7p65fb-1_Inclusive/'
config.Data.lumiMask      = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-276097_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'  
config.Data.publication   = False


config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

config.JobType.pyCfgParams = list(pyCfgParams)
