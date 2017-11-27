import os
import glob
import math
from array import array
import sys
import time
import subprocess
import ROOT
import random

from optparse import OptionParser
from subprocess import Popen

############################################                                                                                                                                              
#            Job steering                  #                                                                                                                                                     
############################################                                                                                                                                                      

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

## parse files                                                                                                                                                                                        
parser.add_option('--inputDIRGridPack', action="store", type="string", dest="inputDIRGridPack",   default="",   help="input directory with all tarballs containing the gridpack --> must be on EOS")
parser.add_option('--outputDIR',        action="store", type="string", dest="outputDIR",          default="",   help="output directory on eos where storing lhe file --> must be on EOS")
parser.add_option('--numberOfEvents',   action="store", type=int,      dest="numberOfEvents",     default=-1,   help="total number of events to generate")
parser.add_option('--eventsPerJob',     action="store", type=int,      dest="eventsPerJob",       default=-1,   help="max number of events per job")
parser.add_option('--jobDIR',           action="store", type="string", dest="jobDIR",             default="",   help="directory for creating jobs")
parser.add_option('--queque',           action="store", type="string", dest="queque",             default="",   help="queque for LSF")
parser.add_option('--submit',           action="store_true",           dest="submit",                           help="submit")

(options, args) = parser.parse_args()

if __name__ == '__main__':

    print "######################################";
    print "### Submit jobs for LHE generation ###";
    print "######################################";

    os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select ls "+options.inputDIRGridPack+" | grep tar > file.list");
    currentDIR = os.getcwd();

    inputfile = open("file.list","r");
    
    for filename in inputfile:

        jobname = file.replace(".tar.gz","");

        ### create output directory
        print "/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+options.outputDIR;
        os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+options.outputDIR);
        os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+options.outputDIR+"/"+jobname);
        
        ### split in the right number of jobs
        if options.numberOfEvents/options.eventsPerJob - int(options.numberOfEvents/options.eventsPerJob) != 0:
            njobs =  int(options.numberOfEvents/options.eventsPerJob)+1
        else:
            njobs =  int(options.numberOfEvents/options.eventsPerJob)

            print "options.numberOfEvents ",options.numberOfEvents," options.eventsPerJob ",options.eventsPerJob," njobs ",njobs;

        ### create general job directory
        os.system("mkdir -p "+options.jobDIR+"_"+jobname);
        
        for ijob in range(njobs):
            jobscript = open('%s/job_%d.sh'%(options.jobDIR+"/"+jobname,ijob),'w')
            jobscript.write("cd "+currentDIR+" \n");
            jobscript.write("eval `scramv1 runtime -sh` \n");
            jobscript.write("xrdcp root://eoscms.cern.ch//eos/cms/"+options.inputDIRGridPack+"/"+filename+" ./ \n");
            jobscript.write("tar -xvf "+filename+" \n");
            jobscript.write("./runcmsgrid.sh %d %d 1 \n"%(options.eventsPerJob,random.randint(0,1000000)));
            jobscript.write("gunzip events.lhe.gz \n");
            jobscript.write("scp events.lhe events_%d.lhe \n"%(ijob));
            jobscript.write("xrdcp events_%d.lhe root://eoscms.cern.ch//eos/cms/%s \n"%(ijob,options.outputDIR+"/"+jobname));

            os.system("chmod a+x %s/job_%d.sh"%(options.jobDIR+"/"+jobname,ijob));
            
            if options.submit:
                os.system('bsub -q %s -o %s/%s/job_%d.log -e %s/%s/job_%d.err %s/%s/job_%d.sh'%(options.queque,currentDIR,options.jobDIR+"/"+jobname,ijob,currentDIR,options.jobDIR+"/"+jobname,ijob,currentDIR,options.jobDIR+"/"+jobname,ijob));
