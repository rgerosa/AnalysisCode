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
parser.add_option('--inputFile',    action="store", type="string", dest="inputFile",    default="",   help="input file")
parser.add_option('--filterName',   action="store", type="string", dest="filterName",   default="",   help="filter name")
parser.add_option('--isMC',         action="store_true", dest="isMC",         help="isMC")
parser.add_option('--applyBTagSF',  action="store_true", dest="applyBTagSF",  help="applyBTagSF")
parser.add_option('--storeGenTree', action="store_true", dest="storeGenTree", help="storeGenTree")
parser.add_option('--inputIsDIR',   action="store_true", dest="inputIsDIR",   help="inputIsDIR")

(options, args) = parser.parse_args()

if __name__ == '__main__':

    print "################################";    
    print "##### Start job submission #####";
    print "################################";    
    
    ROOT.gROOT.ProcessLine(".L macros/filters.C");

    os.chdir(options.inputDIR);
    fileList = [];

    ## if the input is a directory
    if options.inputIsDIR:

         os.system("ls "+options.inputFile+" | grep -v txt | grep root  > file_temp.txt");
         fs = open("file_temp.txt","r");
         for line in fs:
             line = line.replace('\n','');
             fileList.append(line);

         os.system("rm file_temp.txt");
    
    else:
        fileList.append(options.inputFile);
        os.system("mkdir -p "+options.outputDIR);


    ## fix options
    storeGenTree = 0;
    if options.storeGenTree:
        storeGenTree = 1;

    isMC = 0;
    if not options.isMC:
        isMC = 1;

    applyBTagSF = 0;
    if not options.applyBTagSF:
        applyBTagSF = 1;

    #######################
    for ifile in fileList:        

        ########
        if options.filterName == "sigfilter" or options.filterName == "all":
            
            command = ROOT.TString("sigfilter(\"%s\",\"%s\",%i,%i,%i)"%(ifile,"sig_"+ifile,isMC,applyBTagSF,storeGenTree))
            print command
            ROOT.gROOT.ProcessLine(command.Data());
            os.system("mkdir -p sigfilter")
            os.system("mv sig_"+ifile+" sigfilter/")

        ########
        if options.filterName == "zmmfilter"  or options.filterName == "all":

            command = ROOT.TString("zmmfilter(\"%s\",\"%s\",%i,%i,%i)"%(ifile,"zmm_"+ifile,isMC,applyBTagSF,storeGenTree))
            print command
            ROOT.gROOT.ProcessLine(command.Data());
            os.system("mkdir -p zmmfilter")
            os.system("mv zmm_"+ifile+" zmmfilter/")

        ########
        if options.filterName == "zeefilter" or  options.filterName == "all":

            command = ROOT.TString("zeefilter(\"%s\",\"%s\",%i,%i,%i)"%(ifile,"zee_"+ifile,isMC,applyBTagSF,storeGenTree))
            print command
            ROOT.gROOT.ProcessLine(command.Data());
            os.system("mkdir -p zeefilter")
            os.system("mv zee_"+ifile+" zeefilter/")

        ########
        if options.filterName == "wmnfilter" or  options.filterName == "all":

            command = ROOT.TString("wmnfilter(\"%s\",\"%s\",%i,%i,%i)"%(ifile,"wmn_"+ifile,isMC,applyBTagSF,storeGenTree))
            print command
            ROOT.gROOT.ProcessLine(command.Data());
            os.system("mkdir -p wmnfilter")
            os.system("mv wmn_"+ifile+" wmnfilter/")

        ########
        if options.filterName == "wenfilter" or  options.filterName == "all":

            command = ROOT.TString("wenfilter(\"%s\",\"%s\",%i,%i,%i)"%(ifile,"wen_"+ifile,isMC,applyBTagSF,storeGenTree))
            print command
            ROOT.gROOT.ProcessLine(command.Data());
            os.system("mkdir -p wenfilter")
            os.system("mv wen_"+ifile+" wenfilter/")

        ########
        if options.filterName == "gamfilter" or  options.filterName == "all":

            command = ROOT.TString("gamfilter(\"%s\",\"%s\",%i,%i,%i)"%(ifile,"gam_"+ifile,isMC,applyBTagSF,storeGenTree))
            print command
            ROOT.gROOT.ProcessLine(command.Data());
            os.system("mkdir -p gamfilter")
            os.system("mv gam_"+ifile+" gamfilter/")
