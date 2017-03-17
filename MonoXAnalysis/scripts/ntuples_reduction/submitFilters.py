### script to submit runFilters.py on all the sub-directories found inside a mother one located on eos
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
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",     default="",   help="input directory where files are contained")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",    default="",   help="output DIR")
parser.add_option('--filterName',   action="store", type="string", dest="filterName",   default="",   help="filter name to be run .. all means all filters")
parser.add_option('--calculateXSfromSW', action="store_true",      dest="calculateXSfromSW",          help="calculateXSfromSW means use event weight for XS")
parser.add_option('--calculateXSfromLHE',action="store_true",      dest="calculateXSfromLHE",         help="calculateXSfromLHE means take the values in the LHE")
parser.add_option('--isMC',         action="store_true",           dest="isMC",                       help="isMC")
parser.add_option('--applyBTagSF',  action="store_true",           dest="applyBTagSF",                help="applyBTagSF")
parser.add_option('--storeGenTree', action="store_true",           dest="storeGenTree",               help="storeGenTree")
parser.add_option('--isSinglePhoton', action="store_true",         dest="isSinglePhoton",             help="isSinglePhoton")
parser.add_option('--isJetHT',      action="store_true",           dest="isJetHT",                    help="isJetHT")
parser.add_option('--isDoubleEG',   action="store_true",           dest="isDoubleEG",                 help="isDoubleEG")
parser.add_option('--isCrabDirectory', action="store_true",        dest="isCrabDirectory",            help="isCrabDirectory: when the input directory has been created by crab with many files")
parser.add_option('--dropPuppiBranches',   action="store_true",    dest="dropPuppiBranches",          help="drop all puppi branches")
parser.add_option('--dropPuppiBoostedJets', action="store_true",   dest="dropPuppiBoostedJets",       help="drop all puppi branches for boosted jets")
parser.add_option('--dropSubJetsBranches', action="store_true",    dest="dropSubJetsBranches",        help="drop all subjet branches")
parser.add_option('--dropHLTFilter',       action="store_true",    dest="dropHLTFilter",              help="drop HLT filter requirement")
parser.add_option('--metCut',              action="store", type="string", dest="metCut",        help="apply MET/Recoil threshold")

parser.add_option('--grepName', action="callback", type="string", dest="grepName", default="", callback=foo_callback, help="grep a set of names in the directory") ## useful when submitting only signal or data
parser.add_option('--skipName', action="callback", type="string", dest="skipName", default="", callback=foo_callback, help="drop a set of names in the directory")

##  for submitting jobs in lxbatch                                                                                                                                              
parser.add_option('--jobDIR',       action="store", type="string", dest="jobDIR",  default="",        help="directory for job")
parser.add_option('--queque',       action="store", type="string", dest="queque",  default="",        help="queque for LSF")

(options, args) = parser.parse_args()

if __name__ == '__main__':

    ###### make the directory list
    command = "/afs/cern.ch/project/eos/installation/cms/bin/eos.select ls "+options.inputDIR+" | grep -v txt | grep -v root ";
    for name in options.grepName:
        command += " | grep "+name;
    for name in options.skipName:
        command += " | grep -v "+name;
    command += " > dir_list.txt";
    os.system(command);
    fs = open("dir_list.txt","r");
    dirList = [];
    for line in fs:
        line = line.replace('\n','');
        dirList.append(line);
    os.system("rm dir_list.txt");

    print dirList
    for dir in dirList:
        command = "python scripts/ntuples_reduction/runFilters.py --inputDIR "+options.inputDIR+"/"+dir+" --outputDIR "+options.outputDIR+" --filterName "+options.filterName+" --jobDIR "+options.jobDIR+" --queque "+options.queque+" --batchMode --submit --isOnEOS --metCut "+options.metCut+" ";
        if options.calculateXSfromSW:
          command += "--calculateXSfromSW ";
        if options.calculateXSfromLHE:
          command += "--calculateXSfromLHE ";
        if options.isMC:
          command += "--isMC ";
        if options.applyBTagSF:
          command += "--applyBTagSF ";
        if options.storeGenTree:
          command += "--storeGenTree ";
        if options.isSinglePhoton:
          command += "--isSinglePhoton ";
        if options.isJetHT:
          command += "--isJetHT ";
        if options.isDoubleEG:
          command += "--isDoubleEG ";
        if options.isCrabDirectory:
          command += "--isCrabDirectory ";
        if options.dropPuppiBranches:
          command += "--dropPuppiBranches ";
        if options.dropSubJetsBranches:
          command += "--dropSubJetsBranches ";
        if options.dropPuppiBoostedJets:
          command += "--dropPuppiBoostedJets ";
        if options.dropHLTFilter:
          command += "--dropHLTFilter";        

        print command    
        os.system(command)
