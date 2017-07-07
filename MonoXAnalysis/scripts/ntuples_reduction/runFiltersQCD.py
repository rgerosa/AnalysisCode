#python scripts/runFilters.py --inputDIR /store/caf/user/rgerosa/MONOJET_ANALYSIS/Production-14-1-2016/DMS_Scalar/ --filterName all --isMC --applyBTagSF --batchMode --jobDIR JOB  --submit --queque 1nd
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

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

## parse files                                                                                                                                                                 
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",     default="",   help="input directory where files are contained")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",    default="",   help="output DIR")
parser.add_option('--calculateXSfromSW', action="store_true",      dest="calculateXSfromSW",          help="calculateXSfromSW means use event weight for XS")
parser.add_option('--calculateXSfromLHE',action="store_true",      dest="calculateXSfromLHE",         help="calculateXSfromLHE means take the values in the LHE")
parser.add_option('--isMC',         action="store_true",           dest="isMC",                       help="isMC")
parser.add_option('--applyBTagSF',  action="store_true",           dest="applyBTagSF",                help="applyBTagSF")
parser.add_option('--storeGenTree', action="store_true",           dest="storeGenTree",               help="storeGenTree")
parser.add_option('--isCrabDirectory', action="store_true",        dest="isCrabDirectory",            help="isCrabDirectory: when the input directory has been created by crab with many files")
parser.add_option('--isOnEOS',      action="store_true",           dest="isOnEOS",                    help="isOnEOS when the input directory is located in EOS")
parser.add_option('--dropHLTFilter', action="store_true",          dest="dropHLTFilter",              help="drop HLT filter requirement")
##  for submitting jobs in lxbatch
parser.add_option('--batchMode',    action="store_true",           dest="batchMode",                  help="batchMode")
parser.add_option('--jobDIR',       action="store", type="string", dest="jobDIR",  default="",        help="directory for job")
parser.add_option('--queque',       action="store", type="string", dest="queque",  default="",        help="queque for LSF")
parser.add_option('--submit',       action="store_true",           dest="submit",                     help="submit")

(options, args) = parser.parse_args()

if __name__ == '__main__':

    print "################################";    
    print "##### Start job submission #####";
    print "################################";    

    currentDIR = os.getcwd();
    ## generate binary file
    ROOT.gROOT.ProcessLine(".L $CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeQCDBackgroundEstimation/makeReducedTrees.C+");
    ## not in batchMode
    if not options.batchMode:        
        os.chdir(options.inputDIR);
        os.system("mkdir -p "+options.outputDIR);

    fileList = [];
    ## if the input is a directory --> make the list of files
    if not options.isCrabDirectory and not options.isOnEOS and not options.batchMode:
        os.system("ls | grep -v txt | grep root | grep -v failed   > file_temp.txt");
        fs = open("file_temp.txt","r");
        for line in fs:
            line = line.replace('\n','');
            fileList.append(line);
        os.system("rm file_temp.txt");

    elif not options.isCrabDirectory and options.isOnEOS: ## not a crab directory do the same with eos ls command
        os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select ls "+options.inputDIR+" | grep -v txt | grep -v root | grep -v failed  > file_temp.txt");
        fs = open("file_temp.txt","r");
        for line in fs:
            line = line.replace('\n','');
            fileList.append(line);
        os.system("rm file_temp.txt");
    elif options.isCrabDirectory and options.isOnEOS:
        fileList.append(options.inputDIR);
    else:
        print "Problem in parsing the following job informations: batchMode = ",options.batchMode," isOnEOS = ",options.isOnEOS," isCrabDir ",options.isCrabDirectory;
                
    ## fix options
    storeGenTree = 0;
    if options.storeGenTree:
        storeGenTree = 1;

    isMC = 0;
    if options.isMC:
        isMC = 1;

    applyBTagSF = 0;
    if options.applyBTagSF:
        applyBTagSF = 1;

    isCrabDirectory = 0;
    if options.isCrabDirectory:
        isCrabDirectory = 1;

    isOnEOS = 0;
    if options.isOnEOS:
        isOnEOS = 1;

    dropHLTFilter = 0;
    if options.dropHLTFilter:
        dropHLTFilter = 1;

    xsType = 0;
    if options.calculateXSfromSW:
        xsType = 1;
    if options.calculateXSfromLHE:
        xsType = 2;
    if options.calculateXSfromSW and options.calculateXSfromLHE:
        sys.exit("decide to fix the cross section as sum of weights or taking LHE value");

    ####################### loop on the file list
    for ifile in fileList:        

        fileName    = ifile;
        outFileName = ifile;
        ## in case of running over unmerged crab files: fix input file path and output file name
        if isCrabDirectory:
           nameList = ifile.split("/");
           while True:                 
               if nameList[len(nameList)-1] != '':
                   break
               else:
                   nameList.pop();
           fileName    = nameList[len(nameList)-1]; 
           outFileName = nameList[len(nameList)-1]+".root";
           options.outputDIR += "/"+fileName ## modify this for a better organization of output files

           if not options.batchMode:
                
               command = ROOT.TString("qcdfilter(\"%s\",\"%s\",%i,%i,%i,%i,%i,%i,%i,%i,%i)"%(ifile,"qcd_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropHLTFilter))
               print command
               
               ROOT.gROOT.ProcessLine(command.Data());
               os.system("mkdir -p "+options.outputDIR);
               os.system("mkdir -p "+options.outputDIR+"/qcdfilter")
               os.system("mv qcd_"+outFileName+" "+options.outputDIR+"/qcdfilter/")
               
           else:
               
               ## make job directory
               subdirName = fileName.replace(".root","");
               os.system("mkdir -p "+options.jobDIR)
               os.system("mkdir -p "+options.jobDIR+"/"+"JOB_qcd_"+subdirName);
               
                ## write job sh file
               jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_qcd_"+subdirName),'w')
               jobmacro.write("{\n");
               jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/makeQCDBackgroundEstimation/makeReducedTrees.C+\");\n");
               jobmacro.write("gROOT->ProcessLine(\""+"qcdfilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i,%i,%i,%i,%i)"%(ifile,"qcd_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropHLTFilter)+"\");\n");
                    
               jobmacro.write("}\n");
               jobmacro.close();

               jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_qcd_"+subdirName),'w')
               jobscript.write('cd %s \n'%currentDIR)
               jobscript.write('eval ` scramv1 runtime -sh ` \n')
               jobscript.write('cd - \n')
               if not isCrabDirectory:
                   jobscript.write("xrdcp -f root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+fileName+" ./\n")
               jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_qcd_"+subdirName))
               jobscript.write('root -l -b -q job.C\n');
               jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+options.outputDIR+"\n");
               jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+options.outputDIR+"/qcdfilter/\n");
               jobscript.write("xrdcp -f qcd_"+outFileName+" root://eoscms.cern.ch//eos/cms"+options.outputDIR+"/qcdfilter/");

               os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_qcd_"+subdirName))
                
               if options.submit:
                   os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_qcd_"+subdirName,currentDIR,options.jobDIR,"JOB_qcd_"+subdirName,currentDIR,options.jobDIR,"JOB_qcd_"+subdirName))
                   
