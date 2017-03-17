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

## multicrab section
if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand

    print "################################"
    print "#### Begin multicrab script ####"
    print "################################"
    
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
        
           config.General.requestName = key
           config.Data.inputDataset = value[0]
           config.JobType.pyCfgParams = list(pyCfgParams)
           config.JobType.pyCfgParams.extend(value[1])
        
           # see https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3FAQ#Multiple_submission_fails_with_a
           p = Process(target=submit, args=(config,))
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

