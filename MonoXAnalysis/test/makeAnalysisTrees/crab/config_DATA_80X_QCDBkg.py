import os, sys
from multiprocessing import Process
from WMCore.Configuration import Configuration

#### Basic crab config structure

config = Configuration()

pyCfgParams = ['isMC=False',
               'filterOnHLT=False',
               'filterHighMETEvents=True',
               'metCut=50',
               'applyL2L3Residuals=True',
               'addQGLikelihood=False',
               'addPileupJetID=False',
               'addPuppiJets=True',
               'addPuppiMET=True',
               'addEGMSmear=True',
               'addMETSystematics=True',
               'useOfficialMETSystematics=True',
               'addMETBreakDown=False',
               'addSubstructureCHS=True',
               'addSubstructurePuppi=False',
               'miniAODProcess=RECO',
               'outputFileName=tree.root',
               'isTriggerTree=False',
               'addTriggerObjects=False',
               'isQCDTree=True',
               'usePrivateSQliteJEC=False',
               'JECEra=Spring16_23Sep2016V1',
               'nThreads=4',
               'isCrab=True']

config.section_('General')
config.General.transferLogs = False
config.General.workArea     = 'crab_projects_DATA_80X_QCDBkg'  # Make sure you set this parameter

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
config.Data.outLFNDirBase = '/store/group/phys_exotica/monojet/rgerosa/ProductionMC_21_05_2017_QCDBkg/'
config.Data.lumiMask      = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.runRange
config.Data.publication   = False


config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

config.JobType.pyCfgParams = list(pyCfgParams)
