import os, sys
from multiprocessing import Process
from WMCore.Configuration import Configuration

#### Basic crab config structure

config = Configuration()

pyCfgParams = ['isMC=True',
               'filterOnHLT=False', ## don't filter event according to the trigger bit
               'setHLTFilterFlag=False',
               'filterHighMETEvents=True', ## apply a low recoil cut at 50 GeV
               'metCut=50',                ## 50 GeV
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
               'addSubstructureCHS=True',
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

## multicrab section
if __name__ == '__main__':

    print "################################"
    print "#### Begin multicrab script ####"
    print "################################"

    from CRABAPI.RawCommand import crabCommand
    
    ## to submit jobs
    def submit(config):
        print " to do: ",config
        res = crabCommand('submit', config = config)

    ## python multicrab.py <sample file>
    print sys.argv
    if len(sys.argv) <= 1 :
        print "no arguments?"
        print "Usage to submit:     python multicrab.py samples_file.py <additional info>"
        print "Usage to get status: python multicrab.py folder <additional info>"
        exit()


    samples = {}
    SamplesFile = sys.argv[1]
    print "SamplesFile = ", SamplesFile
    
    additionalConfiguration = ''
    if len(sys.argv) == 4 :
        additionalConfiguration = sys.argv[3]
    print "AdditionalConfiguration = ", additionalConfiguration


    # submit--> if found sample file --> open it
    if os.path.exists(SamplesFile) and not os.path.isdir(SamplesFile) :
        handle = open(SamplesFile,'r')
        exec(handle)
        handle.close()

        # samples to be analysed                   
        for key, value in samples.iteritems():
           print key, ' -> ', value

           ## set name
           config.General.requestName = key
           ## set dataset
           config.Data.inputDataset   = value[0]
           ## set list of python cfg pset parameters extending it
           config.JobType.pyCfgParams = list(pyCfgParams)
           config.JobType.pyCfgParams.extend(value[1])
           ## declare submitter
           p = Process(target=submit, args=(config,))           
           ## start application
           # see https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3FAQ#Multiple_submission_fails_with_a
           p.start()
           p.join()

    # status and resubmit
    else :
        if len(sys.argv) >= 3 :
            if sys.argv[2] == 'report' :
                os.system("ls " + SamplesFile + " | awk '{print \" crab report "   + SamplesFile + "/\"$1" + "\" " + additionalConfiguration + "\"}' | /bin/sh")
            elif sys.argv[2] == 'status' :
                os.system("ls " + SamplesFile + " | awk '{print \" crab status "   + SamplesFile + "/\"$1" + "\" " + additionalConfiguration + "\"}' | /bin/sh")
            elif sys.argv[2] == 'resubmit' :
                os.system("ls " + SamplesFile + " | awk '{print \" crab resubmit " + SamplesFile + "/\"$1" + "\" " + additionalConfiguration + "\"}' | /bin/sh") 
            elif sys.argv[2] == 'kill' :
                os.system("ls " + SamplesFile + " | awk '{print \" crab kill " + SamplesFile + "/\"$1" + "\" " + additionalConfiguration + "\"}' | /bin/sh") 
            else :
                os.system("ls " + SamplesFile + " | awk '{print \" crab status " + SamplesFile + "/\"$1}' | /bin/sh")
