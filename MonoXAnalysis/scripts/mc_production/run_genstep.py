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
parser.add_option('--inputDIR',      action="store", type="string", dest="inputDIR",       default="",   help="Input directory where LHE files are located")
parser.add_option('--outputDIR',     action="store", type="string", dest="outputDIR",      default="",   help="Output directory where GEN files should be copied")
parser.add_option('--pythonScript',  action="store", type="string", dest="pythonScript",   default="",   help="Location of the python configuration for GEN step")
parser.add_option('--isAMCNLO',      action="store_true", dest="isAMCNLO",  help="isAMCNLO")
parser.add_option('--outputFileName', action="store", type="string", dest="outputFileName", default="gentree", help="outputFileName : name of the output root file")
parser.add_option('--partonMultiplicity', action="store", type=int,      dest="partonMultiplicity", default=0, help="partonMultiplicity: to match with parton shower")
parser.add_option('--jobDIR',        action="store", type="string", dest="jobDIR",           default="",   help="directory for creating jobs")
parser.add_option('--queque',        action="store", type="string", dest="queque",           default="",   help="queque for LSF")
parser.add_option('--submit',        action="store_true",           dest="submit",                         help="submit")

(options, args) = parser.parse_args()

if __name__ == '__main__':

    print "################################";
    print "### Submit jobs for GEN step ###";
    print "################################";

    currentDIR = os.getcwd();

    isAMCNLO = False;
    if options.isAMCNLO:
        isAMCNLO = True;

    ### create output directory
    os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+options.outputDIR);
    print "/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir -p "+options.outputDIR;

    ### create general job directory
    os.system("mkdir -p "+options.jobDIR);

    ### make list of files
    listFiles = [];
    os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select ls "+options.inputDIR+" > file.list");
    list = open("file.list","r");
    for iline in list:
        if iline != "" and iline.find(".root"):
            iline = iline.replace("\n","");
            listFiles.append(options.inputDIR+"/"+iline);

    os.system("rm file.list");
    ijob = 0;

    cmsswConfig = options.pythonScript.split("/");
    cmsswConfig = cmsswConfig[len(cmsswConfig)-1];
    
    for ifile in listFiles:
        
        lheFile = ifile.split("/");
        lheFile = lheFile[len(lheFile)-1];
        jobscript = open('%s/job_%d.sh'%(options.jobDIR,ijob),'w')
        jobscript.write("cd "+currentDIR+" \n");
        jobscript.write("eval `scramv1 runtime -sh` \n");
        jobscript.write("cd - \n");
        jobscript.write("scp "+currentDIR+"/"+options.pythonScript+"  ./ \n");
        jobscript.write("xrdcp root://eoscms.cern.ch//eos/cms/"+ifile+" ./ \n");
        jobscript.write("sed -i -- 's/<custom_folder>/custom_folder/g' "+lheFile+" \n"); 
        jobscript.write("cmsRun %s partonMultiplicity=%d isAMCNLO=%d outputFileName=%s_%d.root inputFiles=file:%s \n"%(cmsswConfig,options.partonMultiplicity,isAMCNLO,options.outputFileName,ijob,lheFile));
        jobscript.write("xrdcp -f %s_%d.root root://eoscms.cern.ch//eos/cms/%s"%(options.outputFileName,ijob,options.outputDIR));        
        os.system("chmod a+x %s/job_%d.sh"%(options.jobDIR,ijob));
        ijob = ijob+1;
        if options.submit:
            os.system('bsub -q %s -o %s/%s/job_%d.log -e %s/%s/job_%d.err %s/%s/job_%d.sh'%(options.queque,currentDIR,options.jobDIR,ijob,currentDIR,options.jobDIR,ijob,currentDIR,options.jobDIR,ijob));

            
