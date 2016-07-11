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

muonPtBinning  = [10.,20.,30.,40.,50.,70.,100.,135.,500.];
muonEtaBinning = [-2.4,-1.2,0.0,1.2,2.4];

electronPtBinning  = [10.,20.,30.,40.,50.,70.,100.,135.,500.];
electronEtaBinning = [-2.5,-1.4,0.0,1.4,2.5];

photonPtBinning  = [10.,20.,30.,40.,50.,70.,100.,135.,500.];
photonEtaBinning = [-2.5,-1.4,0.0,1.4,2.5];

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
parser.add_option('--absetaBin',         action="store_true",      dest="absetaBin",                  help="bin in abs eta instead of eta")

##  for submitting jobs in lxbatch
parser.add_option('--batchMode',    action="store_true",           dest="batchMode",                  help="batchMode")
parser.add_option('--jobDIR',       action="store", type="string", dest="jobDIR",  default="",        help="directory for job")
parser.add_option('--queque',       action="store", type="string", dest="queque",  default="",        help="queque for LSF")
parser.add_option('--submit',       action="store_true",           dest="submit",                     help="submit")

(options, args) = parser.parse_args()

if __name__ == '__main__':

    isMC = False;
    if options.isMC == True:
        isMC = True;

    os.system("mkdir -p "+options.outputDIR);

    currentDIR = os.getcwd();

    add_option = "";
    if options.absetaBin:
        add_option = "absEta=True";

    if options.leptonType == "muon":
        for pt in range(len(muonPtBinning)-1):
            for eta in range(len(muonEtaBinning)-1):
                if not options.batchMode:
                    command = "cmsRun tnpanalysis.py isMC="+str(isMC)+" inputDIR="+options.inputDIR+" outputDIR="+options.outputDIR+" typeID="+options.typeID+" leptonPID="+str(13)+" ptMin="+str(muonPtBinning[pt])+" ptMax="+str(muonPtBinning[pt+1])+" etaMin="+str(muonEtaBinning[eta])+" etaMax="+str(muonEtaBinning[eta+1])+" "+add_option;
                    ### look for the template file
                    templatePath = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/TagAndProbeTemplates/TemplateNominal/')
                    os.system("ls "+templatePath+" | grep root | grep "+options.leptonType+" | grep "+options.typeID+" | grep pt_"+str(muonPtBinning[pt])+"_"+str(muonPtBinning[pt+1])+"_eta_"+str(muonEtaBinning[eta])+"_"+str(muonEtaBinning[eta+1])+" > file_temp_"+options.leptonType+"_"+options.typeID);
                    file = open("file_temp_"+options.leptonType+"_"+options.typeID,"r");
                    listOffile = [];
                    for line in file:
                        if line == "" or line =="\n": continue;
                        listOffile.append(line.replace("\n",""));
                    if len(listOffile) > 1:
                        sys.exit('Problem more than one template file for a single configuration --> return');
                    command += " templateFile="+templatePath+"/"+listOffile[0];
                    os.system(command+" backgroundType=RooCMSShape");
                    ## do alternative background
                    if options.doAlternativeBkg:
                        os.system(command+" backgroundType=Exponential");
                    ## do alternative signal template
                    if options.doAlternativeSig:
                        templatePathAlt = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/TagAndProbeTemplates/TemplateAlternative/')
                        command = command.replace(templatePath,templatePathAlt);
                        os.system(command+" backgroundType=RooCMSShape signalType=Alternative");

                else:
                    command = "cmsRun tnpanalysis.py isMC="+str(isMC)+" inputDIR="+options.inputDIR+" isEOSDIR=True outputDIR=./ typeID="+options.typeID+" leptonPID="+str(13)+" ptMin="+str(muonPtBinning[pt])+" ptMax="+str(muonPtBinning[pt+1])+" etaMin="+str(muonEtaBinning[eta])+" etaMax="+str(muonEtaBinning[eta+1])+" "+add_option;
                    ### look for the template file
                    templatePath = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/TagAndProbeTemplates/TemplateNominal/')
                    os.system("ls "+templatePath+" | grep root | grep "+options.leptonType+" | grep "+options.typeID+" | grep pt_"+str(muonPtBinning[pt])+"_"+str(muonPtBinning[pt+1])+"_eta_"+str(muonEtaBinning[eta])+"_"+str(muonEtaBinning[eta+1])+" > file_temp_"+options.leptonType+"_"+options.typeID);
                    file = open("file_temp_"+options.leptonType+"_"+options.typeID,"r");
                    listOffile = [];
                    for line in file:
                        if line == "" or line =="\n": continue;
                        listOffile.append(line.replace("\n",""));
                    if len(listOffile) > 1:
                        sys.exit('Problem more than one template file for a single configuration --> return');
                    command += " templateFile="+templatePath+"/"+listOffile[0];
                    ### submit jobs
                    os.system("mkdir -p "+options.jobDIR);
                    jobName = 'job_%s_%s_pt_%.1f_%.1f_eta_%.1f_%.1f'%(options.leptonType,options.typeID,muonPtBinning[pt],muonPtBinning[pt+1],muonEtaBinning[eta],muonEtaBinning[eta+1])
                    jobscript = open('%s/%s_RooCMSShape.sh'%(options.jobDIR,jobName),'w');
                    jobscript.write('cd %s \n'%currentDIR)
                    jobscript.write('eval ` scramv1 runtime -sh ` \n')
                    jobscript.write('cd - \n')
                    jobscript.write('scp '+currentDIR+'/tnpanalysis.py ./ \n')                    
                    jobscript.write(command+' backgroundType=RooCMSShape\n');
                    jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*root "+currentDIR+'/'+options.outputDIR+' \n')                    
                    os.system('chmod a+x %s/%s_RooCMSShape.sh'%(options.jobDIR,jobName))

                    if options.submit:
                        os.system('bsub -q %s -o %s/%s_RooCMSShape.log -e %s/%s_RooCMSShape.err %s/%s_RooCMSShape.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));

                    if options.doAlternativeBkg:

                        jobscript = open('%s/%s_Exp.sh'%(options.jobDIR,jobName),'w');
                        jobscript.write('cd %s \n'%currentDIR)
                        jobscript.write('eval ` scramv1 runtime -sh ` \n')
                        jobscript.write('cd - \n')
                        jobscript.write('scp '+currentDIR+'/tnpanalysis.py ./ \n')                    
                        jobscript.write(command+' backgroundType=Exponential\n');
                        jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*root "+currentDIR+'/'+options.outputDIR+' \n')                    
                        os.system('chmod a+x %s/%s_Exp.sh'%(options.jobDIR,jobName))

                        if options.submit:
                            os.system('bsub -q %s -o %s/%s_Exp.log -e %s/%s_Exp.err %s/%s_Exp.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));

                    if options.doAlternativeSig:
                        templatePathAlt = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/TagAndProbeTemplates/TemplateAlternative/')
                        command = command.replace(templatePath,templatePathAlt);
                        jobscript = open('%s/%s_Alternative.sh'%(options.jobDIR,jobName),'w');
                        jobscript.write('cd %s \n'%currentDIR)
                        jobscript.write('eval ` scramv1 runtime -sh ` \n')
                        jobscript.write('cd - \n')
                        jobscript.write('scp '+currentDIR+'/tnpanalysis.py ./ \n')
                        jobscript.write(command+' backgroundType=RooCMSShape signalType=Alternative\n');
                        jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*root "+currentDIR+'/'+options.outputDIR+' \n')
                        os.system('chmod a+x %s/%s_Alternative.sh'%(options.jobDIR,jobName))

                        if options.submit:
                            os.system('bsub -q %s -o %s/%s_Alternative.log -e %s/%s_Alternative.err %s/%s_Alternative.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));
                        

    elif options.leptonType == "electron":
        for pt in range(len(electronPtBinning)-1):
            for eta in range(len(electronEtaBinning)-1):
                if not options.batchMode:
                    command = "cmsRun tnpanalysis.py isMC="+str(isMC)+" inputDIR="+options.inputDIR+" outputDIR="+options.outputDIR+" typeID="+options.typeID+" leptonPID="+str(11)+" ptMin="+str(electronPtBinning[pt])+" ptMax="+str(electronPtBinning[pt+1])+" etaMin="+str(electronEtaBinning[eta])+" etaMax="+str(electronEtaBinning[eta+1])+" "+add_option

                    ### look for the template file
                    templatePath = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/TagAndProbeTemplates/TemplateNominal/')
                    os.system("ls "+templatePath+" | grep root | grep "+options.leptonType+" | grep "+options.typeID+" | grep pt_"+str(electronPtBinning[pt])+"_"+str(electronPtBinning[pt+1])+"_eta_"+str(electronEtaBinning[eta])+"_"+str(electronEtaBinning[eta+1])+" > file_temp_"+options.leptonType+"_"+options.typeID);
                    file = open("file_temp_"+options.leptonType+"_"+options.typeID,"r");
                    listOffile = [];
                    for line in file:
                        if line == "" or line =="\n": continue;
                        listOffile.append(line.replace("\n",""));
                    if len(listOffile) > 1:
                        sys.exit('Problem more than one template file for a single configuration --> return');
                    command += " templateFile="+templatePath+"/"+listOffile[0];
                    ## run first fit with RooCMSShape
                    os.system(command+" backgroundType=RooCMSShape");

                    if options.doAlternativeBkg:
                        os.system(command+" backgroundType=Exponential");

                    ## do alternative signal template
                    if options.doAlternativeSig:
                        templatePathAlt = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/TagAndProbeTemplates/TemplateAlternative/')
                        command = command.replace(templatePath,templatePathAlt);
                        os.system(command+" backgroundType=RooCMSShape signalType=Alternative");

                else:
                    command = "cmsRun tnpanalysis.py isMC="+str(isMC)+" inputDIR="+options.inputDIR+" isEOSDIR=True outputDIR=./"+" typeID="+options.typeID+" leptonPID="+str(11)+" ptMin="+str(electronPtBinning[pt])+" ptMax="+str(electronPtBinning[pt+1])+" etaMin="+str(electronEtaBinning[eta])+" etaMax="+str(electronEtaBinning[eta+1])+ " "+add_option
                    ### look for the template file
                    templatePath = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/TagAndProbeTemplates/TemplateNominal/')
                    os.system("ls "+templatePath+" | grep root | grep "+options.leptonType+" | grep "+options.typeID+" | grep pt_"+str(electronPtBinning[pt])+"_"+str(electronPtBinning[pt+1])+"_eta_"+str(electronEtaBinning[eta])+"_"+str(electronEtaBinning[eta+1])+" > file_temp_"+options.leptonType+"_"+options.typeID);
                    file = open("file_temp_"+options.leptonType+"_"+options.typeID,"r");
                    listOffile = [];
                    for line in file:
                        if line == "" or line =="\n": continue;
                        listOffile.append(line.replace("\n",""));
                    if len(listOffile) > 1:
                        sys.exit('Problem more than one template file for a single configuration --> return');
                    command += " templateFile="+templatePath+"/"+listOffile[0];
                    ### submit jobs
                    os.system("mkdir -p "+options.jobDIR);
                    jobName = 'job_%s_%s_pt_%.1f_%.1f_eta_%.1f_%.1f'%(options.leptonType,options.typeID,electronPtBinning[pt],electronPtBinning[pt+1],electronEtaBinning[eta],electronEtaBinning[eta+1])
                    jobscript = open('%s/%s_RooCMSShape.sh'%(options.jobDIR,jobName),'w');
                    jobscript.write('cd %s \n'%currentDIR)
                    jobscript.write('eval ` scramv1 runtime -sh ` \n')
                    jobscript.write('cd - \n')
                    jobscript.write('scp '+currentDIR+'/tnpanalysis.py ./\n')                    
                    jobscript.write(command+' backgroundType=RooCMSShape\n');
                    jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*root "+currentDIR+'/'+options.outputDIR+' \n')                    
                    os.system('chmod a+x %s/%s_RooCMSShape.sh'%(options.jobDIR,jobName))

                    if options.submit:
                        os.system('bsub -q %s -o %s/%s_RooCMSShape.log -e %s/%s_RooCMSShape.err %s/%s_RooCMSShape.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));

                    if options.doAlternativeBkg:
                        jobscript = open('%s/%s_Exp.sh'%(options.jobDIR,jobName),'w');
                        jobscript.write('cd %s \n'%currentDIR)
                        jobscript.write('eval ` scramv1 runtime -sh ` \n')
                        jobscript.write('cd - \n')
                        jobscript.write('scp '+currentDIR+'/tnpanalysis.py ./\n')                    
                        jobscript.write(command+' backgroundType=Exponential\n');
                        jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*root "+currentDIR+'/'+options.outputDIR+' \n')                    
                        os.system('chmod a+x %s/%s_Exp.sh'%(options.jobDIR,jobName))

                        if options.submit:
                            os.system('bsub -q %s -o %s/%s_Exp.log -e %s/%s_Exp.err %s/%s_Exp.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));

                    if options.doAlternativeSig:
                        templatePathAlt = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/TagAndProbeTemplates/TemplateAlternative/')
                        command = command.replace(templatePath,templatePathAlt);
                        jobscript = open('%s/%s_Alternative.sh'%(options.jobDIR,jobName),'w');
                        jobscript.write('cd %s \n'%currentDIR)
                        jobscript.write('eval ` scramv1 runtime -sh ` \n')
                        jobscript.write('cd - \n')
                        jobscript.write('scp '+currentDIR+'/tnpanalysis.py ./ \n')
                        jobscript.write(command+' backgroundType=RooCMSShape signalType=Alternative\n');
                        jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*root "+currentDIR+'/'+options.outputDIR+' \n')
                        os.system('chmod a+x %s/%s_Alternative.sh'%(options.jobDIR,jobName))

                        if options.submit:
                            os.system('bsub -q %s -o %s/%s_Alternative.log -e %s/%s_Alternative.err %s/%s_Alternative.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));


    elif options.leptonType == "photon":
        for pt in range(len(photonPtBinning)-1):
            for eta in range(len(photonEtaBinning)-1):
                if not options.batchMode:
                    command = "cmsRun tnpanalysis.py isMC="+str(isMC)+" inputDIR="+options.inputDIR+" outputDIR="+options.outputDIR+" typeID="+options.typeID+" leptonPID="+str(22)+" ptMin="+str(photonPtBinning[pt])+" ptMax="+str(photonPtBinning[pt+1])+" etaMin="+str(photonEtaBinning[eta])+" etaMax="+str(photonEtaBinning[eta+1])+ " "+add_option

                    ### look for the template file
                    templatePath = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/TagAndProbeTemplates/TemplateNominal/')
                    os.system("ls "+templatePath+" | grep root | grep "+options.leptonType+" | grep "+options.typeID+" | grep pt_"+str(photonPtBinning[pt])+"_"+str(photonPtBinning[pt+1])+"_eta_"+str(photonEtaBinning[eta])+"_"+str(photonEtaBinning[eta+1])+" > file_temp_"+options.leptonType+"_"+options.typeID);
                    file = open("file_temp_"+options.leptonType+"_"+options.typeID,"r");
                    listOffile = [];
                    for line in file:
                        if line == "" or line =="\n": continue;
                        listOffile.append(line.replace("\n",""));
                    if len(listOffile) > 1:
                        sys.exit('Problem more than one template file for a single configuration --> return');
                    command += " templateFile="+templatePath+"/"+listOffile[0];
                    ## run first fit with RooCMSShape
                    os.system(command+" backgroundType=RooCMSShape");

                    if options.doAlternativeBkg:
                        os.system(command+" backgroundType=Exponential");

                    ## do alternative signal template
                    if options.doAlternativeSig:
                        templatePathAlt = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/TagAndProbeTemplates/TemplateAlternative/')
                        command = command.replace(templatePath,templatePathAlt);
                        os.system(command+" backgroundType=RooCMSShape signalType=Alternative");

                else:
                    command = "cmsRun tnpanalysis.py isMC="+str(isMC)+" inputDIR="+options.inputDIR+" isEOSDIR=True outputDIR=./ typeID="+options.typeID+" leptonPID="+str(22)+" ptMin="+str(photonPtBinning[pt])+" ptMax="+str(photonPtBinning[pt+1])+" etaMin="+str(photonEtaBinning[eta])+" etaMax="+str(photonEtaBinning[eta+1])+" "+add_option
                    ### look for the template file
                    templatePath = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/TagAndProbeTemplates/TemplateNominal/')
                    os.system("ls "+templatePath+" | grep root | grep "+options.leptonType+" | grep "+options.typeID+" | grep pt_"+str(photonPtBinning[pt])+"_"+str(photonPtBinning[pt+1])+"_eta_"+str(photonEtaBinning[eta])+"_"+str(photonEtaBinning[eta+1])+" > file_temp_"+options.leptonType+"_"+options.typeID);
                    file = open("file_temp_"+options.leptonType+"_"+options.typeID,"r");
                    listOffile = [];
                    for line in file:
                        if line == "" or line =="\n": continue;
                        listOffile.append(line.replace("\n",""));
                    if len(listOffile) > 1:
                        sys.exit('Problem more than one template file for a single configuration --> return');
                    command += " templateFile="+templatePath+"/"+listOffile[0];
                    ## submit jobs
                    os.system("mkdir -p "+options.jobDIR);
                    jobName = 'job_%s_%s_pt_%.1f_%.1f_eta_%.1f_%.1f'%(options.leptonType,options.typeID,photonPtBinning[pt],photonPtBinning[pt+1],photonEtaBinning[eta],photonEtaBinning[eta+1])
                    jobscript = open('%s/%s_RooCMSShape.sh'%(options.jobDIR,jobName),'w');
                    jobscript.write('cd %s \n'%currentDIR)
                    jobscript.write('eval ` scramv1 runtime -sh ` \n')
                    jobscript.write('cd - \n')
                    jobscript.write('scp '+currentDIR+'/tnpanalysis.py ./\n')                    
                    jobscript.write(command+' backgroundType=RooCMSShape\n');
                    jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*root "+currentDIR+'/'+options.outputDIR+' \n')                    
                    os.system('chmod a+x %s/%s_RooCMSShape.sh'%(options.jobDIR,jobName))

                    if options.submit:
                        os.system('bsub -q %s -o %s/%s_RooCMSShape.log -e %s/%s_RooCMSShape.err %s/%s_RooCMSShape.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));

                    if options.doAlternativeBkg:
                        jobscript = open('%s/%s_Exp.sh'%(options.jobDIR,jobName),'w');
                        jobscript.write('cd %s \n'%currentDIR)
                        jobscript.write('eval ` scramv1 runtime -sh ` \n')
                        jobscript.write('cd - \n')
                        jobscript.write('scp '+currentDIR+'/tnpanalysis.py ./\n')                    
                        jobscript.write(command+' backgroundType=Exponential\n');
                        jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*root "+currentDIR+'/'+options.outputDIR+' \n')                    
                        os.system('chmod a+x %s/%s_Exp.sh'%(options.jobDIR,jobName))

                        if options.submit:
                            os.system('bsub -q %s -o %s/%s_Exp.log -e %s/%s_Exp.err %s/%s_Exp.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));

                    if options.doAlternativeSig:
                        templatePathAlt = os.path.expandvars('$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/TagAndProbeTemplates/TemplateAlternative/')
                        command = command.replace(templatePath,templatePathAlt);
                        jobscript = open('%s/%s_Alternative.sh'%(options.jobDIR,jobName),'w');
                        jobscript.write('cd %s \n'%currentDIR)
                        jobscript.write('eval ` scramv1 runtime -sh ` \n')
                        jobscript.write('cd - \n')
                        jobscript.write('scp '+currentDIR+'/tnpanalysis.py ./ \n')
                        jobscript.write(command+' backgroundType=RooCMSShape signalType=Alternative\n');
                        jobscript.write('scp efficiency*'+options.leptonType+"*pt*eta*root "+currentDIR+'/'+options.outputDIR+' \n')
                        os.system('chmod a+x %s/%s_Alternative.sh'%(options.jobDIR,jobName))

                        if options.submit:
                            os.system('bsub -q %s -o %s/%s_Alternative.log -e %s/%s_Alternative.err %s/%s_Alternative.sh'%(options.queque,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName,currentDIR+"/"+options.jobDIR,jobName));


    else:
        sys.exit('Problem with lepton type --> muon or electron or photon are the recognized options --> exit');

    os.system("rm file_temp_"+options.leptonType+"_"+options.typeID);
