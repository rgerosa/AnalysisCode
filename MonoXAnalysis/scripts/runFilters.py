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
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",     default="",   help="input DIR")
parser.add_option('--filterName',   action="store", type="string", dest="filterName",   default="",   help="filter name")
parser.add_option('--isMC',         action="store_true", dest="isMC",         help="isMC")
parser.add_option('--applyBTagSF',  action="store_true", dest="applyBTagSF",  help="applyBTagSF")
parser.add_option('--storeGenTree', action="store_true", dest="storeGenTree", help="storeGenTree")

##  for submitting jobs in lxbatch
parser.add_option('--batchMode',    action="store_true", dest="batchMode",   help="batchMode")
parser.add_option('--jobDIR',       action="store", type="string", dest="jobDIR", default="",   help="directory for job")
parser.add_option('--queque',        action="store", type="string", dest="queque", default="",   help="queque for LSF")
parser.add_option('--submit',       action="store_true", dest="submit", help="submit")

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

    fileList = [];

    ## if the input is a directory
    if not options.batchMode:
        os.system("ls | grep -v txt | grep root  > file_temp.txt");
        fs = open("file_temp.txt","r");
        for line in fs:
            line = line.replace('\n','');
            fileList.append(line);
        os.system("rm file_temp.txt");

    else:
        os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select ls "+options.inputDIR+" | grep -v txt | grep root  > file_temp.txt");
        fs = open("file_temp.txt","r");
        for line in fs:
            line = line.replace('\n','');
            fileList.append(line);
        os.system("rm file_temp.txt");

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


    #######################
    for ifile in fileList:        

        ########
        if options.filterName == "sigfilter" or options.filterName == "all":

            if not options.batchMode:
                
                command = ROOT.TString("sigfilter(\"%s\",\"%s\",%i,%i,%i)"%(ifile,"sig_"+ifile,isMC,applyBTagSF,storeGenTree))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p sigfilter")
                os.system("mv sig_"+ifile+" sigfilter/")
                
            else:

                ## make job directory
                subdirName = ifile.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_sig_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_sig_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"sigfilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i)"%(ifile,"sig_"+ifile,isMC,applyBTagSF,storeGenTree)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_sig_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                jobscript.write("xrdcp root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+ifile+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_sig_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.inputDIR+"/sigfilter/\n");
                jobscript.write("xrdcp sig_"+ifile+" root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/sigfilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_sig_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_sig_"+subdirName,currentDIR,options.jobDIR,"JOB_sig_"+subdirName,currentDIR,options.jobDIR,"JOB_sig_"+subdirName))
                   

        ########
        if options.filterName == "zmmfilter"  or options.filterName == "all":

            
            if not options.batchMode:

                command = ROOT.TString("zmmfilter(\"%s\",\"%s\",%i,%i,%i)"%(ifile,"zmm_"+ifile,isMC,applyBTagSF,storeGenTree))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p zmmfilter")
                os.system("mv zmm_"+ifile+" zmmfilter/")

            else:

                ## make job directory
                subdirName = ifile.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_zmm_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_zmm_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"zmmfilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i)"%(ifile,"zmm_"+ifile,isMC,applyBTagSF,storeGenTree)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_zmm_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                jobscript.write("xrdcp root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+ifile+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_zmm_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.inputDIR+"/zmmfilter/\n");
                jobscript.write("xrdcp zmm_"+ifile+" root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/zmmfilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_zmm_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_zmm_"+subdirName,currentDIR,options.jobDIR,"JOB_zmm_"+subdirName,currentDIR,options.jobDIR,"JOB_zmm_"+subdirName))

        ########
        if options.filterName == "zeefilter" or  options.filterName == "all":

            if not options.batchMode:

                command = ROOT.TString("zeefilter(\"%s\",\"%s\",%i,%i,%i)"%(ifile,"zee_"+ifile,isMC,applyBTagSF,storeGenTree))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p zeefilter")
                os.system("mv zee_"+ifile+" zeefilter/")

            else:

                ## make job directory
                subdirName = ifile.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_zee_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_zee_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"zeefilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i)"%(ifile,"zee_"+ifile,isMC,applyBTagSF,storeGenTree)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_zee_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                jobscript.write("xrdcp root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+ifile+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_zee_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.inputDIR+"/zeefilter/\n");
                jobscript.write("xrdcp zee_"+ifile+" root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/zeefilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_zee_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_zee_"+subdirName,currentDIR,options.jobDIR,"JOB_zee_"+subdirName,currentDIR,options.jobDIR,"JOB_zee_"+subdirName))

        ########
        if options.filterName == "wmnfilter" or  options.filterName == "all":

            if not options.batchMode:

                command = ROOT.TString("wmnfilter(\"%s\",\"%s\",%i,%i,%i)"%(ifile,"wmn_"+ifile,isMC,applyBTagSF,storeGenTree))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p wmnfilter")
                os.system("mv wmn_"+ifile+" wmnfilter/")

            else:

                ## make job directory
                subdirName = ifile.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_wmn_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_wmn_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"wmnfilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i)"%(ifile,"wmn_"+ifile,isMC,applyBTagSF,storeGenTree)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_wmn_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                jobscript.write("xrdcp root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+ifile+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_wmn_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.inputDIR+"/wmnfilter/\n");
                jobscript.write("xrdcp wmn_"+ifile+" root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/wmnfilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_wmn_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_wmn_"+subdirName,currentDIR,options.jobDIR,"JOB_wmn_"+subdirName,currentDIR,options.jobDIR,"JOB_wmn_"+subdirName))

        ########
        if options.filterName == "wenfilter" or  options.filterName == "all":

            if not options.batchMode:

                command = ROOT.TString("wenfilter(\"%s\",\"%s\",%i,%i,%i)"%(ifile,"wen_"+ifile,isMC,applyBTagSF,storeGenTree))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p wenfilter")
                os.system("mv wen_"+ifile+" wenfilter/")

            else:

                ## make job directory
                subdirName = ifile.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_wen_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_wen_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"wenfilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i)"%(ifile,"wen_"+ifile,isMC,applyBTagSF,storeGenTree)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_wen_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                jobscript.write("xrdcp root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+ifile+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_wen_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.inputDIR+"/wenfilter/\n");
                jobscript.write("xrdcp wen_"+ifile+" root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/wenfilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_wen_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_wen_"+subdirName,currentDIR,options.jobDIR,"JOB_wen_"+subdirName,currentDIR,options.jobDIR,"JOB_wen_"+subdirName))

        ########
        if options.filterName == "gamfilter" or  options.filterName == "all":

            if not options.batchMode:

                command = ROOT.TString("gamfilter(\"%s\",\"%s\",%i,%i,%i)"%(ifile,"gam_"+ifile,isMC,applyBTagSF,storeGenTree))
                print command
                ROOT.gROOT.ProcessLine(command.Data());
                os.system("mkdir -p gamfilter")
                os.system("mv gam_"+ifile+" gamfilter/")

            else:

                ## make job directory
                subdirName = ifile.replace(".root","");
                os.system("mkdir -p "+options.jobDIR)
                os.system("mkdir -p "+options.jobDIR+"/"+"JOB_gam_"+subdirName);

                ## write job sh file
                jobmacro = open('%s/%s/job.C'%(options.jobDIR,"JOB_gam_"+subdirName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/filters.C+\");\n");
                jobmacro.write("gROOT->ProcessLine(\""+"gamfilter(\\\"%s\\\",\\\"%s\\\",%i,%i,%i)"%(ifile,"gam_"+ifile,isMC,applyBTagSF,storeGenTree)+"\");\n");
                jobmacro.write("}\n");
                jobmacro.close();

                jobscript = open('%s/%s/job.sh'%(options.jobDIR,"JOB_gam_"+subdirName),'w')
                jobscript.write('cd %s \n'%currentDIR)
                jobscript.write('eval ` scramv1 runtime -sh ` \n')
                jobscript.write('cd - \n')
                jobscript.write("xrdcp root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/"+ifile+" ./\n")
                jobscript.write('scp '+currentDIR+'/%s/%s/job.C ./ \n'%(options.jobDIR,"JOB_gam_"+subdirName))
                jobscript.write('root -l -b -q job.C\n');
                jobscript.write("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir  -p "+options.inputDIR+"/gamfilter/\n");
                jobscript.write("xrdcp gam_"+ifile+" root://eoscms.cern.ch//eos/cms"+options.inputDIR+"/gamfilter/");

                os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,"JOB_gam_"+subdirName))
                
                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,"JOB_gam_"+subdirName,currentDIR,options.jobDIR,"JOB_gam_"+subdirName,currentDIR,options.jobDIR,"JOB_gam_"+subdirName))
