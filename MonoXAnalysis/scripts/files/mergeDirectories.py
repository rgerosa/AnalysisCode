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

def foo_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

## parse files                                                                                                                                                                
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",      default="",   help="input directory where all the directories to be merged are located")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",     default="",   help="output common directory")
parser.add_option('--grepName',     action="callback", type="string", dest="grepName", default="", callback=foo_callback, help="grep a set of names in the directory")
parser.add_option('--skipName',     action="callback", type="string", dest="skipName", default="", callback=foo_callback, help="drop a set of names in the directory")
(options, args) = parser.parse_args()

filterList = ["sigfilter","wmnfilter","wenfilter","gamfilter","zmmfilter","zeefilter","topmufilter","topelfilter"];

if __name__ == '__main__':

    dirList = [];
    command = "ls "+options.inputDIR+" | grep -v txt | grep -v root ";
    for name in options.grepName:
        command += " | grep "+name;
    for name in options.skipName:
        command += " | grep -v "+name;
    command += " > dir_list.txt";
    os.system(command);
    fs = open("dir_list.txt","r");
    for line in fs:
        line = line.replace('\n','');
        dirList.append(line);
    os.system("rm dir_list.txt");

    ### having the dir list to be merged
    os.system("mkdir -p "+options.outputDIR);

    for dir in dirList:
        for filter in filterList:
            os.system("mkdir -p "+options.outputDIR+"/"+filter);
            command = "mv -u "+options.inputDIR+"/"+dir+"/"+filter+"/*root "+options.outputDIR+"/"+filter;
            os.system(command);

    ## remove all directories
    for dir in dirList:
      os.system("rm -r "+options.inputDIR+"/"+dir);
