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

