#python makeTagAndProbeFits.py --inputDIR /home/rgerosa/MONOJET_ANALYSIS/TagAndProbe_600pb-1/SingleMuon/ --outputDIR testFit --typeID looseid --leptonType muon 
import os
import glob
import math
from array import array
import sys
import time
import subprocess
import ROOT

from optparse import OptionParser
from subprocess import Popen

############################################                                                                                                                                 
#            Job steering                  #                                                                                                                                 
############################################                                                                                                                                   

##################
muonPtBinning       = [10.,20.,30.,40.,50.,60.,80.,100.,500.];
muonEtaBinning      = [0.,0.9,1.2,2.1,2.4];
muonNvtxBinning     = [0.,17.,50.];
muonPtBinningReco   = [10.,20.,30.,40.,150.];
muonEtaBinningReco  = [0.,0.4,0.9,1.2,1.8,2.4];
muonNvtxBinningReco = [0.,17.,50.];

electronPtBinning       = [10.,20.,30.,40.,50.,60.,80.,100.,125.,500.];
electronEtaBinning      = [0.,0.8,1.5,2.,2.5];
electronNvtxBinning     = [0.,17.,50.];
electronPtBinningReco   = [10.,20.,30.,50.,150.];
electronEtaBinningReco  = [0.,0.3,0.8,1.5,2.0,2.5];
electronNvtxBinningReco = [0.,17.,50.];

photonPtBinning   = [10.,20.,30.,40.,50.,60.,80.,100.,125.,150.,500.];
photonEtaBinning  = [0.,0.5,1.0,1.5];
photonNvtxBinning = [0.,17.,50.];

##################

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

## parse files                                                                                                                                                                 
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",     default="",   help="input directory where files are contained")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",    default="",   help="output DIR")
parser.add_option('--isMC',         action="store_true",           dest="isMC",                       help="isMC")
parser.add_option('--typeID',       action="store", type="string", dest="typeID",       default="",   help="lepton id type: looseid or tightid for muons, vetoid or tightid for electrons")
parser.add_option('--leptonType',   action="store", type="string", dest="leptonType",   default="",   help="lepton type: muon or electron or photon")
parser.add_option('--doAlternativeBkg',  action="store_true",      dest="doAlternativeBkg",           help="run the fits with a exponential background shape")
parser.add_option('--doAlternativeSig',  action="store_true",      dest="doAlternativeSig",           help="run the fits with alternative signal template")
parser.add_option('--doAnalyticalFit',   action="store_true",      dest="doAnalyticalFit",            help="run analytical fits signal template")
parser.add_option('--absetaBin',         action="store_true",      dest="absetaBin",                  help="bin in abs eta instead of eta")
parser.add_option('--isRecoEff',         action="store_true",      dest="isRecoEff",                  help="isRecoEff")
##  for submitting jobs in lxbatch
parser.add_option('--batchMode',    action="store_true",           dest="batchMode",                  help="batchMode")
parser.add_option('--jobDIR',       action="store", type="string", dest="jobDIR",  default="",        help="directory for job")
parser.add_option('--queque',       action="store", type="string", dest="queque",  default="",        help="queque for LSF")
parser.add_option('--submit',       action="store_true",           dest="submit",                     help="submit")

(options, args) = parser.parse_args()
                   
##################################
########### Main Code ############
##################################

if __name__ == '__main__':
                       
    
    #### parsing some options
    isMC = False;
    if options.isMC == True:
        isMC = True;

    isRecoEff = False;
    if options.isRecoEff == True:
        isRecoEff = True;

    os.system("mkdir -p "+options.outputDIR);
    currentDIR = os.getcwd();

    absetaBin = False;
    if options.absetaBin:
        absetaBin = True
    

    #### select the right binning
    leptonPtBinning = [];
    leptonEtaBinning = [];
    leptonNvtxBinning = [];

    if options.leptonType == "muon" and not isRecoEff:
        leptonPtBinning = muonPtBinning;
        leptonEtaBinning = muonEtaBinning;    
        leptonNvtxBinning = muonNvtxBinning;    
    elif options.leptonType == "muon" and isRecoEff:
        leptonPtBinning = muonPtBinningReco;
        leptonEtaBinning = muonEtaBinningReco;    
        leptonNvtxBinning = muonNvtxBinningReco;    
    elif options.leptonType == "electron" and not isRecoEff:
        leptonPtBinning = electronPtBinning;
        leptonEtaBinning = electronEtaBinning;    
        leptonNvtxBinning = electronNvtxBinning;    
    elif options.leptonType == "electron" and isRecoEff:
        leptonPtBinning = electronPtBinningReco;
        leptonEtaBinning = electronEtaBinningReco;    
        leptonNvtxBinning = electronNvtxBinningReco;    
    elif options.leptonType == "photon" and not isRecoEff:
        leptonPtBinning = photonPtBinning;
        leptonEtaBinning = photonEtaBinning;    
        leptonNvtxBinning = photonNvtxBinning;    


    leptonPID = 0;
    if options.leptonType == "muon":
        leptonPID = 13;
    elif options.leptonType == "electron":
        leptonPID = 11;
    elif options.leptonType == "photon":
        leptonPID = 22;
    else:
        sys.exit('Problem with lepton type --> muon or electron or photon are the recognized options --> exit');
        
                
    ###### start loop for jobs 
    for pt in range(len(leptonPtBinning)-1):
        for eta in range(len(leptonEtaBinning)-1):
            for nvtx in range(len(leptonNvtxBinning)-1):
                #### do the fit interactively
                if not options.batchMode:
                    ### look for the template file
                    if not options.isRecoEff :
                        templatePath = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeTagAndProbe/TemplateNominal/')
                    else:
                        templatePath = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeTagAndProbe/TemplateNominalReco/')

                    os.system("ls "+templatePath+" | grep root | grep "+options.leptonType+" | grep "+options.typeID+" | grep pt_"+str(leptonPtBinning[pt])+"_"+str(leptonPtBinning[pt+1])+"_eta_"+str(leptonEtaBinning[eta])+"_"+str(leptonEtaBinning[eta+1])+"_pu_"+str(leptonNvtxBinning[nvtx])+"_"+str(leptonNvtxBinning[nvtx+1])+" > file_temp_"+options.leptonType+"_"+options.typeID);
                    file = open("file_temp_"+options.leptonType+"_"+options.typeID,"r");
                    listOffile = [];
                    for line in file:
                        if line == "" or line =="\n": continue;
                        listOffile.append(line.replace("\n",""));
                        if len(listOffile) > 1:
                            sys.exit('Problem more than one template file for a single configuration --> return');

                    command = "root -l -b -q makeTagAndProbeFits.C\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",%d,TagAndProbeBin\(%f,%f,%f,%f,%f,%f\),\\\"%s\\\",\\\"%s\\\",%d,%d,%d,%d\)"%(options.inputDIR,options.outputDIR,templatePath+"/"+listOffile[0],options.typeID,leptonPID,leptonPtBinning[pt],leptonPtBinning[pt+1],leptonEtaBinning[eta],leptonEtaBinning[eta+1],leptonNvtxBinning[nvtx],leptonNvtxBinning[nvtx+1],"RooCMSShape","",False,isRecoEff,absetaBin,False);
                    os.system(command);

                    if options.doAlternativeBkg:
                        command = "root -l -b -q makeTagAndProbeFits.C\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",%d,TagAndProbeBin\(%f,%f,%f,%f,%f,%f\),\\\"%s\\\",\\\"%s\\\",%d,%d,%d,%d\)"%(options.inputDIR,options.outputDIR,templatePath+"/"+listOffile[0],options.typeID,leptonPID,leptonPtBinning[pt],leptonPtBinning[pt+1],leptonEtaBinning[eta],leptonEtaBinning[eta+1],leptonNvtxBinning[nvtx],leptonNvtxBinning[nvtx+1],"Exponential","",False,isRecoEff,absetaBin,False);
                        os.system(command);
                        #### do the analytical fit
                    if options.doAnalyticalFit:
                        command = "root -l -b -q makeTagAndProbeFits.C\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",%d,TagAndProbeBin\(%f,%f,%f,%f,%f,%f\),\\\"%s\\\",\\\"%s\\\",%d,%d,%d,%d\)"%(options.inputDIR,options.outputDIR,templatePath+"/"+listOffile[0],options.typeID,leptonPID,leptonPtBinning[pt],leptonPtBinning[pt+1],leptonEtaBinning[eta],leptonEtaBinning[eta+1],leptonNvtxBinning[nvtx],leptonNvtxBinning[nvtx+1],"RooCMSShape","",False,isRecoEff,absetaBin,True);
                        os.system(command);

                    ## do alternative signal template
                    if options.doAlternativeSig:
                        if not options.isRecoEff:
                            templatePathAlt = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeTagAndProbe/TemplateAlternative/')
                        else:
                            templatePathAlt = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeTagAndProbe/TemplateAlternativeReco/')

                            command = "root -l -b -q makeTagAndProbeFits.C(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",%d,TagAndProbeBin\(%f,%f,%f,%f,%f,%f\),\\\"%s\\\",\\\"%s\\\",%d,%d,%d,%d\)"%(options.inputDIR,options.outputDIR,templatePath+"/"+listOffile[0],options.typeID,leptonPID,leptonPtBinning[pt],leptonPtBinning[pt+1],leptonEtaBinning[eta],leptonEtaBinning[eta+1],leptonNvtxBinning[nvtx],leptonNvtxBinning[nvtx+1],"RooCMSShape","Alternative",False,isRecoEff,absetaBin,False);
                            command = command.replace(templatePath,templatePathAlt);
                            os.system(command);

                ##### batch mode
                else:

                    if not options.isRecoEff:
                        templatePath = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeTagAndProbe/TemplatesNominal/')
                    else:
                        templatePath = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeTagAndProbe/TemplatesNominalReco/')

                        
                    os.system("ls "+templatePath+" | grep root | grep "+options.leptonType+" | grep "+options.typeID+" | grep pt_"+str(leptonPtBinning[pt])+"_"+str(leptonPtBinning[pt+1])+"_eta_"+str(leptonEtaBinning[eta])+"_"+str(leptonEtaBinning[eta+1])+"_pu_"+str(leptonNvtxBinning[nvtx])+"_"+str(leptonNvtxBinning[nvtx+1])+" > file_temp_"+options.leptonType+"_"+options.typeID);
                    file = open("file_temp_"+options.leptonType+"_"+options.typeID,"r");
                    listOffile = [];
                    for line in file:
                        if line == "" or line =="\n": continue;
                        listOffile.append(line.replace("\n",""));
                        if len(listOffile) > 1:
                            sys.exit('Problem more than one template file for a single configuration --> return');
                            
                    command = "root -l -b -q %s/makeTagAndProbeFits.C\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",%d,TagAndProbeBin\(%f,%f,%f,%f,%f,%f\),\\\"%s\\\",\\\"%s\\\",%d,%d,%d,%d\)"%(currentDIR,options.inputDIR,"./",templatePath+"/"+listOffile[0],options.typeID,leptonPID,leptonPtBinning[pt],leptonPtBinning[pt+1],leptonEtaBinning[eta],leptonEtaBinning[eta+1],leptonNvtxBinning[nvtx],leptonNvtxBinning[nvtx+1],"RooCMSShape","",True,isRecoEff,absetaBin,False);


                    ### submit jobs
                    os.system("mkdir -p "+options.jobDIR);                    
                    jobName = 'job_%s_%s_pt_%.1f_%.1f_eta_%.1f_%.1f_pu_%.1f_%.1f'%(options.leptonType,options.typeID,leptonPtBinning[pt],leptonPtBinning[pt+1],leptonEtaBinning[eta],leptonEtaBinning[eta+1],leptonNvtxBinning[nvtx],leptonNvtxBinning[nvtx+1])

                    jobscript = open('%s/%s_RooCMSShape.sh'%(options.jobDIR,jobName),'w');
                    jobscript.write('cd %s \n'%currentDIR)
                    jobscript.write('eval ` scramv1 runtime -sh ` \n')
                    jobscript.write('cd - \n')
                    jobscript.write(command+'\n');
                    jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*pu*root "+currentDIR+'/'+options.outputDIR+' \n')                    
                    os.system('chmod a+x %s/%s_RooCMSShape.sh'%(options.jobDIR,jobName))

                    if options.submit:
                        os.system('bsub -q %s -o %s/%s_RooCMSShape.log -e %s/%s_RooCMSShape.err %s/%s_RooCMSShape.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));

                    ### alternative background
                    if options.doAlternativeBkg:
                        
                        jobscript = open('%s/%s_Exp.sh'%(options.jobDIR,jobName),'w');
                        jobscript.write('cd %s \n'%currentDIR)
                        jobscript.write('eval ` scramv1 runtime -sh ` \n')
                        jobscript.write('cd - \n')

                        command = "root -l -b -q %s/makeTagAndProbeFits.C\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",%d,TagAndProbeBin\(%f,%f,%f,%f,%f,%f\),\\\"%s\\\",\\\"%s\\\",%d,%d,%d,%d\)"%(currentDIR,options.inputDIR,"./",templatePath+"/"+listOffile[0],options.typeID,leptonPID,leptonPtBinning[pt],leptonPtBinning[pt+1],leptonEtaBinning[eta],leptonEtaBinning[eta+1],leptonNvtxBinning[nvtx],leptonNvtxBinning[nvtx+1],"Exponential","",True,isRecoEff,absetaBin,False);
                        jobscript.write(command+'\n');
                        jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*root "+currentDIR+'/'+options.outputDIR+' \n')                    
                        os.system('chmod a+x %s/%s_Exp.sh'%(options.jobDIR,jobName))
                        
                        if options.submit:
                            os.system('bsub -q %s -o %s/%s_Exp.log -e %s/%s_Exp.err %s/%s_Exp.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));

                        ### analytical fit
                        if options.doAnalyticalFit:
                            jobscript = open('%s/%s_Analytical.sh'%(options.jobDIR,jobName),'w');
                            jobscript.write('cd %s \n'%currentDIR)
                            jobscript.write('eval ` scramv1 runtime -sh ` \n')
                            jobscript.write('cd - \n')
                            command = "root -l -b -q %s/makeTagAndProbeFits.C\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",%d,TagAndProbeBin\(%f,%f,%f,%f,%f,%f\),\\\"%s\\\",\\\"%s\\\",%d,%d,%d,%d\)"%(currentDIR,options.inputDIR,"./",templatePath+"/"+listOffile[0],options.typeID,leptonPID,leptonPtBinning[pt],leptonPtBinning[pt+1],leptonEtaBinning[eta],leptonEtaBinning[eta+1],leptonNvtxBinning[nvtx],leptonNvtxBinning[nvtx+1],"RooCMSShape","",True,isRecoEff,absetaBin,True);
                            jobscript.write(command+'\n');
                            jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*pu*Analytical*root "+currentDIR+'/'+options.outputDIR+' \n')
                            os.system('chmod a+x %s/%s_Analytical.sh'%(options.jobDIR,jobName))

                            if options.submit:
                                os.system('bsub -q %s -o %s/%s_Analytical.log -e %s/%s_Analytical.err %s/%s_Analytical.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));


                        #### alternative signal
                        if options.doAlternativeSig:
                            if not options.isRecoEff:
                                templatePathAlt = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeTagAndProbe/TemplateAlternative/')
                            else:
                                templatePathAlt = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeTagAndProbe/TemplateAlternativeReco/')


                            command = "root -l -b -q %s/makeTagAndProbeFits.C\(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",%d,TagAndProbeBin\(%f,%f,%f,%f,%f,%f\),\\\"%s\\\",\\\"%s\\\",%d,%d,%d,%d\)"%(currentDIR,options.inputDIR,"./",templatePath+"/"+listOffile[0],options.typeID,leptonPID,leptonPtBinning[pt],leptonPtBinning[pt+1],leptonEtaBinning[eta],leptonEtaBinning[eta+1],leptonNvtxBinning[nvtx],leptonNvtxBinning[nvtx+1],"RooCMSShape","Alternative",True,isRecoEff,absetaBin,True);
                            command = command.replace(templatePath,templatePathAlt);
                            jobscript = open('%s/%s_Alternative.sh'%(options.jobDIR,jobName),'w');
                            jobscript.write('cd %s \n'%currentDIR)
                            jobscript.write('eval ` scramv1 runtime -sh ` \n')
                            jobscript.write('cd - \n')
                            jobscript.write(command+'\n');
                            jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*root "+currentDIR+'/'+options.outputDIR+' \n')
                            os.system('chmod a+x %s/%s_Alternative.sh'%(options.jobDIR,jobName))

                            if options.submit:
                                os.system('bsub -q %s -o %s/%s_Alternative.log -e %s/%s_Alternative.err %s/%s_Alternative.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));

    os.system("rm file_temp_"+options.leptonType+"_"+options.typeID);
