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

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

## parse files                                                                                                                                                                                         
parser.add_option('--inputFile1',   action="store", type="string", dest="inputFile1",  default="",   help="input directory where files are located on eos")
parser.add_option('--inputFile2',   action="store", type="string", dest="inputFile2",  default="",   help="number of parallel streams")
parser.add_option('--modeltype',    action="store", type="string", dest="modeltype",  default="",   help="number of parallel streams")


(options, args) = parser.parse_args()

if __name__ == '__main__':

    model_code = -1;
    if options.modeltype == "vector":
        model_code = 800;
    elif options.modeltype == "axial":
        model_code = 801;
    elif options.modeltype == "scalar":
        model_code = 805;
    elif options.modeltype == "pseudo":
        model_code = 806;


    file1 = ROOT.TFile(options.inputFile1,"READ")
    file2 = ROOT.TFile(options.inputFile2,"READ")

    workspace1 = file1.Get("combinedws");
    workspace2 = file2.Get("combinedws");

    filelist1 = workspace1.allData ()  
    filelist2 = workspace2.allData () 

    massPoint_1 = []
    massPoint_2 = []
    noffshellpoints = 0;
    for key1 in filelist1:
        name1 = str(key1.GetName())
        list1 = name1.split("_")
        massPoint_1.append('%d'%(int(list1[3])))

    for key2 in filelist2:
        name2 = str(key2.GetName())
        list2 = name2.split("_")
        massPoint_2.append('%d'%(int(list2[3])))
        
    common_point = [];
    non_common_point = [];

    for point in massPoint_1:
        if point in massPoint_2:

            mh = int(point);
            mmed = -1;
            mdm  = -1;
            if model_code == 800: mmed = ((int)(mh-80000000000))/10000;
            if model_code == 801: mmed = ((int)(mh-80100000000))/10000;
            if model_code == 805: mmed = ((int)(mh-80500000000))/10000;
            if model_code == 806: mmed = ((int)(mh-80600000000))/10000;

            if model_code == 800: mdm =  (mh-80000000000)  - ( ((int)(mh-80000000000))/10000 )*10000
            if model_code == 801: mdm =  (mh-80100000000)  - ( ((int)(mh-80100000000))/10000 )*10000
            if model_code == 805: mdm =  (mh-80500000000)  - ( ((int)(mh-80500000000))/10000 )*10000
            if model_code == 806: mdm =  (mh-80600000000)  - ( ((int)(mh-80600000000))/10000 )*10000

            if mdm*2 > mmed :
                noffshellpoints = noffshellpoints +1;

            common_point.append(point)
        else:
            non_common_point.append(point);

    for point in massPoint_2:
        if not point in massPoint_1:
            non_common_point.append(point);

    print "common point ",common_point;
    print "non common point ",non_common_point;

    print "number of common point ",len(common_point)
    print "number of non common point ",len(non_common_point)
    print "number of off-shell points ",noffshellpoints
