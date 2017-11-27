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

def foo_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--inputDIR',  action="store", type="string",   dest="inputDIR",     default="",   help="input directory where files are contained")
parser.add_option('--grepName',  action="callback", type="string", dest="grepName",    default="", callback=foo_callback, help="grep a set of names in the directory") 
parser.add_option('--maxFiles',  action="store", type=int     , dest="maxFiles",   default=1, help="max files per mass point") 

(options, args) = parser.parse_args()

if __name__ == '__main__':

    command = "/afs/cern.ch/project/eos/installation/cms/bin/eos.select ls "+options.inputDIR; 
    for name in options.grepName:
        command += " | grep "+name;
    command =  command +" > dirList.tmp ";
    os.system(command);
    file = open("dirList.tmp","r");
    for ifile in file:
      ifile = ifile.replace("\n","");
      print ifile;

      os.system("/afs/cern.ch/project/eos/installation/cms/bin/eos.select find "+options.inputDIR+"/"+ifile+" -name \"*root\" > fileList.tmp");
      fileList = open("fileList.tmp","r");
      inputFiles = [];
      command = "cmsRun test/makeGenTrees/genAnalyzer.py inputFiles="      
      nfile = 0;
      for jfile in fileList:
        if not ".root" in jfile:
          continue;
        if nfile >= options.maxFiles:
          break;
        jfile = jfile.replace("\n","");
        jfile = jfile.replace("/eos/cms","");
        command = command+jfile+","        
        nfile = nfile +1;
      command = command[:-1];
      command = command + " debugName="+ifile+".log useDebug=True"
      print command
      os.system(command)
      
      
    os.system("rm fileList.tmp");
    os.system("rm dirList.tmp");
        
