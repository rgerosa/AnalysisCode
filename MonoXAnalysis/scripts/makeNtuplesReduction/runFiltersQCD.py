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
parser.add_option('--splitPerFile', action="store_true",           dest="splitPerFile",               help="splitPerFile")

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
    elif options.isCrabDirectory and options.isOnEOS and not options.splitPerFile:
        fileList.append(options.inputDIR);
    elif options.splitPerFile and options.isOnEOS and options.isCrabDirectory:
        os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select find "+options.inputDIR+" | grep -v txt | grep root | grep -v failed  > file_temp.txt");
        fs = open("file_temp.txt","r");
        for line in fs:
            line = line.replace('\n','');
            fileList.append(line);
        os.system("rm file_temp.txt");
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

    splitPerFile = 0;
    if options.splitPerFile:
        splitPerFile = 1;

    xsType = 0;
    if options.calculateXSfromSW:
        xsType = 1;
    if options.calculateXSfromLHE:
        xsType = 2;
    if options.calculateXSfromSW and options.calculateXSfromLHE:
        sys.exit("decide to fix the cross section as sum of weights or taking LHE value");

    ####################### loop on the file list
    print "Number of files or directories:", len(fileList)
    nFile = 0;
    for ifile in fileList:        
        nFile = nFile+1;
        fileName    = ifile;
        outFileName = ifile;
        fileDir = ifile;
        outputDIR = options.outputDIR
        ## in case of running over unmerged crab files: fix input file path and output file name
        if options.isCrabDirectory:
           nameList = ifile.split("/");
           while True:                 
               if nameList[len(nameList)-1] != '':
                   break
               else:
                   nameList.pop();
           if not splitPerFile:
               fileName    = nameList[len(nameList)-1];                
               outFileName = nameList[len(nameList)-1]+".root";
               outputDIR += "/"+fileName ## modify this for a better organization of output files
           else:
               fileDir = nameList[len(nameList)-4];
               outFileName = nameList[len(nameList)-1];
               outputDIR += "/"+fileDir ## modify this for a better organization of output files
               
        if not options.batchMode:
                
            command = ROOT.TString("qcdfilter(\"%s\",\"%s\",%i,%i,%i,%i,%i,%i,%i,%i,%i)"%(ifile,"qcd_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropHLTFilter))
            print command
            
            ROOT.gROOT.ProcessLine(command.Data());
            os.system("mkdir -p "+outputDIR);
            os.system("mkdir -p "+outputDIR+"/qcdfilter")
            os.system("mv qcd_"+outFileName+" "+outputDIR+"/qcdfilter/")
               
        else:

            if splitPerFile and isCrabDirectory :
                isCrabDirectory = 0;

            ## make job directory
            if splitPerFile:
                subdirName = fileDir+"/job_%d"%(nFile)
            else:
                subdirName = fileName.replace(".root","");

            os.system("mkdir -p "+options.jobDIR)
            os.system("mkdir -p "+options.jobDIR+"/"+"JOB_qcd_"+subdirName);

            print subdirName;

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
            jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_qcd_"+subdirName))
            jobscript.write('root -l -b -q job.C\n');
            jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+outputDIR+"\n");
            jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+outputDIR+"/qcdfilter/\n");
            jobscript.write("xrdcp -f qcd_"+outFileName+" root://eoscms.cern.ch//eos/cms"+outputDIR+"/qcdfilter/");

            os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_qcd_"+subdirName))
                
            if options.submit:
                os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_qcd_"+subdirName,currentDIR,options.jobDIR,"JOB_qcd_"+subdirName,currentDIR,options.jobDIR,"JOB_qcd_"+subdirName))                   
