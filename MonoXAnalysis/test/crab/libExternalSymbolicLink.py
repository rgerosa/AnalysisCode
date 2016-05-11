import os
import glob
import math
from array import array
import sys
import time
import subprocess

### make symbolic links to all external libraries into the lib one before submitting crab jobs
if __name__ == '__main__':

    ### set workig directory in the lib dir
    libDIR = os.path.expandvars('$CMSSW_BASE/lib/$SCRAM_ARCH/')
    os.chdir(libDIR);
    externalDIR = os.path.expandvars('$CMSSW_BASE/external/$SCRAM_ARCH/lib');
    os.system("ls "+externalDIR+" | grep .a > list.temp");
    os.system("ls "+externalDIR+" | grep .so >> list.temp");
    os.system("ls "+externalDIR+" | grep .la >> list.temp");
    
    fs = open("list.temp","r");
    libList = [];
    for line in fs:
        line = line.replace('\n','');
        if line == "" : continue;
        libList.append(line);

    os.system("rm list.temp");
    for lib in libList:
        command = "ln -sf $CMSSW_BASE/external/$SCRAM_ARCH/lib/"+lib+" "+lib;
        print command
        os.system("ln -sf $CMSSW_BASE/external/$SCRAM_ARCH/lib/"+lib+" "+lib);
        

