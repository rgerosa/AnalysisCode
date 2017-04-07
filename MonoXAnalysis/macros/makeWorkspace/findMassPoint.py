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

    filelist1 = file1.GetListOfKeys() 
    filelist2 = file2.GetListOfKeys() 

    massPoint_1 = []
    massPoint_2 = []

    for key1 in filelist1:
        name1 = str(key1.GetName())
        list1 = name1.split("_")
        if "MonoZ" in name1 or "MonoW" in name1 or "MonoJ" in name1:
            if (model_code == 800 or model_code == 801) and list1[4] != "0.25" : continue;         
            if (model_code == 805 or model_code == 806) and list1[4] != "1.0" : continue; 
            massPoint_1.append('%d%04d%04d'%(model_code,int(list1[2]),int(list1[3])))
        else:
            if (model_code == 800 or model_code == 801) and list1[3] != "0.25" : continue;         
            if (model_code == 805 or model_code == 806) and list1[3] != "1.0" : continue; 
            massPoint_1.append('%d%04d%04d'%(model_code,int(list1[1]),int(list1[2])))

    for key2 in filelist2:
        name2 = str(key2.GetName())
        list2 = name2.split("_")
        if "MonoZ" in name2 or "MonoW" in name2 or "MonoJ" in name2:
            if (model_code == 800 or model_code == 801) and list2[4] != "0.25" : continue;         
            if (model_code == 805 or model_code == 806) and list2[4] != "1.0" : continue; 
            massPoint_2.append('%d%04d%04d'%(model_code,int(list2[2]),int(list2[3])))
        else:
            if (model_code == 800 or model_code == 801) and list2[3] != "0.25" : continue;         
            if (model_code == 805 or model_code == 806) and list2[3] != "1.0" : continue; 
            massPoint_2.append('%d%04d%04d'%(model_code,int(list2[1]),int(list2[2])))

    common_point = [];
    non_common_point = [];

    for point in massPoint_1:
        if point in massPoint_2:
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
