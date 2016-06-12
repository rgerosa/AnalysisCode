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
parser.add_option('--filterName',   action="store", type="string", dest="filterName",   default="",   help="filter name to be run .. all means all filters")
parser.add_option('--calculateXSfromSW', action="store_true",      dest="calculateXSfromSW",          help="calculateXSfromSW means use event weight for XS")
parser.add_option('--calculateXSfromLHE',action="store_true",      dest="calculateXSfromLHE",         help="calculateXSfromLHE means take the values in the LHE")
parser.add_option('--isMC',         action="store_true",           dest="isMC",                       help="isMC")
parser.add_option('--applyBTagSF',  action="store_true",           dest="applyBTagSF",                help="applyBTagSF")
parser.add_option('--storeGenTree', action="store_true",           dest="storeGenTree",               help="storeGenTree")
parser.add_option('--isSinglePhoton', action="store_true",         dest="isSinglePhoton",             help="isSinglePhoton")
parser.add_option('--isCrabDirectory', action="store_true",        dest="isCrabDirectory",            help="isCrabDirectory: when the input directory has been created by crab with many files")
parser.add_option('--isOnEOS',      action="store_true",           dest="isOnEOS",                    help="isOnEOS when the input directory is located in EOS")
parser.add_option('--dropPuppiBranches',    action="store_true",   dest="dropPuppiBranches",          help="drop all puppi branches for ak4 and met")
parser.add_option('--dropPuppiBoostedJets', action="store_true",   dest="dropPuppiBoostedJets",       help="drop all puppi branches for boosted jets")
parser.add_option('--dropSubJetsBranches',  action="store_true",   dest="dropSubJetsBranches",        help="drop all subjet branches")
parser.add_option('--dropHLTFilter', action="store_true",          dest="dropHLTFilter",              help="drop HLT filter requirement")
parser.add_option('--metCut',        action="store", type="string",dest="metCut",                     help="met/Recoil threshold to be applied")
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
    ROOT.gROOT.ProcessLine(".L macros/filters.C+");
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
        os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select ls "+options.inputDIR+" | grep -v txt | grep root | grep -v failed  > file_temp.txt");
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

    isSinglePhoton = 0;
    if options.isSinglePhoton:
        isSinglePhoton = 1;

    isCrabDirectory = 0;
    if options.isCrabDirectory:
        isCrabDirectory = 1;

    isOnEOS = 0;
    if options.isOnEOS:
        isOnEOS = 1;

    dropPuppiBranches = 0
    if options.dropPuppiBranches:
        dropPuppiBranches = 1;

    dropPuppiBoostedJets = 0;
    if options.dropPuppiBoostedJets:
        dropPuppiBoostedJets = 1;

    dropSubJetsBranches = 0;
    if options.dropSubJetsBranches:
        dropSubJetsBranches = 1;

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
        ########
        if options.filterName == "sigfilter" or options.filterName == "all":

            if not options.batchMode:
                
                command = ROOT.TString("sigfilter(\"%s\",\"%s\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\"%s\")"%(ifile,"sig_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut))
                print command
                
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p "+options.outputDIR);
                os.system("mkdir -p "+options.outputDIR+"/sigfilter")
                os.system("mv sig_"+outFileName+" "+options.outputDIR+"/sigfilter/")
                
            else:

                ## make job directory
                subdirName = fileName.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_sig_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_sig_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"sigfilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\\\"%s\\\")"%(ifile,"sig_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut)+"\");\n");
                    
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_sig_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                if not isCrabDirectory:
                    jobscript.write("xrdcp -f root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+fileName+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_sig_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+options.outputDIR+"\n");
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+options.outputDIR+"/sigfilter/\n");
                jobscript.write("xrdcp -f sig_"+outFileName+" root://eoscms.cern.ch//eos/cms"+options.outputDIR+"/sigfilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_sig_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_sig_"+subdirName,currentDIR,options.jobDIR,"JOB_sig_"+subdirName,currentDIR,options.jobDIR,"JOB_sig_"+subdirName))
                   

        ########
        if options.filterName == "zmmfilter"  or options.filterName == "all":

            
            if not options.batchMode:

                command = ROOT.TString("zmmfilter(\"%s\",\"%s\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\"%s\")"%(ifile,"zmm_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p "+options.outputDIR);
                os.system("mkdir -p "+options.outputDIR+"/zmmfilter")
                os.system("mv zmm_"+outFileName+" "+options.outputDIR+"/zmmfilter/")

            else:

                ## make job directory
                subdirName = fileName.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_zmm_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_zmm_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"zmmfilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\\\"%s\\\")"%(ifile,"zmm_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_zmm_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                if not isCrabDirectory:
                    jobscript.write("xrdcp -f root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+fileName+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_zmm_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"\n");
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"/zmmfilter/\n");
                jobscript.write("xrdcp -f zmm_"+outFileName+" root://eoscms.cern.ch//eos/cms"+options.outputDIR+"/zmmfilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_zmm_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_zmm_"+subdirName,currentDIR,options.jobDIR,"JOB_zmm_"+subdirName,currentDIR,options.jobDIR,"JOB_zmm_"+subdirName))

        ########
        if options.filterName == "zeefilter" or  options.filterName == "all":

            if not options.batchMode:

                command = ROOT.TString("zeefilter(\"%s\",\"%s\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\"%s\")"%(ifile,"zee_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,isSinglePhoton,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p "+options.outputDIR);
                os.system("mkdir -p "+options.outputDIR+"/zeefilter")
                os.system("mv zee_"+outFileName+" "+options.outputDIR+"/zeefilter/")

            else:

                ## make job directory
                subdirName = fileName.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_zee_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_zee_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"zeefilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\\\"%s\\\")"%(ifile,"zee_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,isSinglePhoton,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_zee_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                if not isCrabDirectory:
                    jobscript.write("xrdcp -f root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+fileName+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_zee_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"\n");
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"/zeefilter/\n");
                jobscript.write("xrdcp -f zee_"+outFileName+" root://eoscms.cern.ch//eos/cms"+options.outputDIR+"/zeefilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_zee_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_zee_"+subdirName,currentDIR,options.jobDIR,"JOB_zee_"+subdirName,currentDIR,options.jobDIR,"JOB_zee_"+subdirName))

        ########
        if options.filterName == "wmnfilter" or  options.filterName == "all":

            if not options.batchMode:

                command = ROOT.TString("wmnfilter(\"%s\",\"%s\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\"%s\")"%(ifile,"wmn_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p "+options.outputDIR);
                os.system("mkdir -p "+options.outputDIR+"/wmnfilter")
                os.system("mv wmn_"+outFileName+" "+options.outputDIR+"/wmnfilter/")

            else:

                ## make job directory
                subdirName = fileName.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_wmn_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_wmn_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"wmnfilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\\\"%s\\\")"%(ifile,"wmn_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_wmn_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                if not isCrabDirectory:
                    jobscript.write("xrdcp -f root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+fileName+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_wmn_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"\n");
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"/wmnfilter/\n");
                jobscript.write("xrdcp -f wmn_"+outFileName+" root://eoscms.cern.ch//eos/cms"+options.outputDIR+"/wmnfilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_wmn_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_wmn_"+subdirName,currentDIR,options.jobDIR,"JOB_wmn_"+subdirName,currentDIR,options.jobDIR,"JOB_wmn_"+subdirName))

        ########
        if options.filterName == "wenfilter" or  options.filterName == "all":

            if not options.batchMode:

                command = ROOT.TString("wenfilter(\"%s\",\"%s\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\"%s\")"%(ifile,"wen_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,isSinglePhoton,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p "+options.outputDIR);
                os.system("mkdir -p "+options.outputDIR+"/wenfilter")
                os.system("mv wen_"+outFileName+" "+options.outputDIR+"/wenfilter/")

            else:

                ## make job directory
                subdirName = fileName.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_wen_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_wen_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"wenfilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\\\"%s\\\")"%(ifile,"wen_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,isSinglePhoton,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_wen_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                if not isCrabDirectory:
                    jobscript.write("xrdcp -f root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+fileName+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_wen_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"\n");
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"/wenfilter/\n");
                jobscript.write("xrdcp -f wen_"+outFileName+" root://eoscms.cern.ch//eos/cms"+options.outputDIR+"/wenfilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_wen_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_wen_"+subdirName,currentDIR,options.jobDIR,"JOB_wen_"+subdirName,currentDIR,options.jobDIR,"JOB_wen_"+subdirName))

        ########
        if options.filterName == "gamfilter" or  options.filterName == "all":

            if not options.batchMode:

                command = ROOT.TString("gamfilter(\"%s\",\"%s\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\"%s\")"%(ifile,"gam_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p "+options.outputDIR);
                os.system("mkdir -p "+options.outputDIR+"/gamfilter")
                os.system("mv gam_"+outFileName+" "+options.outputDIR+"/gamfilter/")

            else:

                ## make job directory
                subdirName = fileName.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_gam_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_gam_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"gamfilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\\\"%s\\\")"%(ifile,"gam_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_gam_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                if not isCrabDirectory:
                    jobscript.write("xrdcp -f root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+fileName+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_gam_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"\n");
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"/gamfilter/\n");
                jobscript.write("xrdcp -f gam_"+outFileName+" root://eoscms.cern.ch//eos/cms"+options.outputDIR+"/gamfilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_gam_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_gam_"+subdirName,currentDIR,options.jobDIR,"JOB_gam_"+subdirName,currentDIR,options.jobDIR,"JOB_gam_"+subdirName))


        ########
        if options.filterName == "topmufilter" or  options.filterName == "all":

            if not options.batchMode:

                command = ROOT.TString("topmufilter(\"%s\",\"%s\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\"%s\")"%(ifile,"topmu_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p "+options.outputDIR);
                os.system("mkdir -p "+options.outputDIR+"/topmufilter")
                os.system("mv top_"+fileName+" "+options.outputDIR+"/topmufilter/")

            else:

                ## make job directory
                subdirName = fileName.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_topmu_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_topmu_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"topmufilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\\\"%s\\\")"%(ifile,"topmu_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_topmu_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                if not isCrabDirectory:
                    jobscript.write("xrdcp -f root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+fileName+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_topmu_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"\n");
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"/topmufilter/\n");
                jobscript.write("xrdcp -f topmu_"+outFileName+" root://eoscms.cern.ch//eos/cms"+options.outputDIR+"/topmufilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_topmu_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_topmu_"+subdirName,currentDIR,options.jobDIR,"JOB_topmu_"+subdirName,currentDIR,options.jobDIR,"JOB_topmu_"+subdirName))


        ########
        if options.filterName == "topelfilter" or  options.filterName == "all":

            if not options.batchMode:

                command = ROOT.TString("topelfilter(\"%s\",\"%s\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\"%s\")"%(ifile,"topel_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p "+options.outputDIR);
                os.system("mkdir -p "+options.outputDIR+"/topelfilter")
                os.system("mv top_"+fileName+" "+options.outputDIR+"/topelfilter/")

            else:

                ## make job directory
                subdirName = fileName.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_topel_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_topel_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"topelfilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,\\\"%s\\\")"%(ifile,"topel_"+outFileName,isMC,applyBTagSF,isCrabDirectory,isOnEOS,xsType,storeGenTree,dropPuppiBranches,dropPuppiBoostedJets,dropSubJetsBranches,dropHLTFilter,options.metCut)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_topel_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                if not isCrabDirectory:
                    jobscript.write("xrdcp -f root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+fileName+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_topel_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"\n");
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.outputDIR+"/topelfilter/\n");
                jobscript.write("xrdcp -f topel_"+outFileName+" root://eoscms.cern.ch//eos/cms"+options.outputDIR+"/topelfilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_topel_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_topel_"+subdirName,currentDIR,options.jobDIR,"JOB_topel_"+subdirName,currentDIR,options.jobDIR,"JOB_topel_"+subdirName))
