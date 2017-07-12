import os, sys
from multiprocessing import Process
from WMCore.Configuration import Configuration

#### Basic crab config structure

config = Configuration()

pyCfgParams = ['isMC=False',
               'filterOnHLT=True', ## apply HLT path filters
               'filterHighMETEvents=True', ## apply met filter for high recoil events
               'metCut=190',
               'applyL2L3Residuals=True',
               'addQGLikelihood=True',
               'addPileupJetID=False',
               'addPuppiJets=False', ## store puppi jets
               'addPuppiMET=False', ## store puppi met
               'addEGMSmear=True',
               'addEGMRegression=False',
               'addMETSystematics=True', ## add sys in the tree
               'useOfficialMETSystematics=True', ## use official met tool
               'addMETBreakDown=False',
               'addSubstructureCHS=True',
               'addSubstructurePuppi=False',
               'usePrivateSQliteJEC=False',
               'JECEra=Summer16_23Sep2016V3',
               'outputFileName=tree.root',
               'nThreads=4',
               'isReMiniAOD=True', ## re-miniaod option for data
               'useMiniAODPuppiMet=True', ## use puppi met reading from miniAOD directly
               'useMiniAODMet=False', ## MUEGclean redone on the flu
               'isCrab=True']

config.section_('General')
config.General.transferLogs = False
config.General.workArea     = 'crab_projects_DATA_80X_reminiAOD'  # Make sure you set this parameter

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
config.Data.unitsPerJob   = 45000
config.Data.outLFNDirBase = '/store/group/phys_exotica/monojet/rgerosa/ProductionData_ReReco_36fb-1_reMiniAOD/'
config.Data.lumiMask      = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.runRange
config.Data.publication   = False

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'

config.JobType.pyCfgParams = list(pyCfgParams)
