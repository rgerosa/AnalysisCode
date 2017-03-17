
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
parser.add_option('--inputFile',    action="store", type="string", dest="inputFile",     default="",   help="input directory where files are located on eos")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",     default="",   help="input directory where files are located on eos")
(options, args) = parser.parse_args()

if __name__ == '__main__':


    os.system("mkdir -p "+options.outputDIR);
    input_file = open(options.inputFile,'r');

    histograms = [];
    bins = [];
    start_histo = False;

    for line in input_file:
        line = line.replace('\n','');
        if "# xlow" in line:
            start_histo = True;
            continue;
        elif "# END HISTO1D" in line:
            start_histo = False;
            break;
        elif start_histo == True:
            if(len(bins) == 0):
                bins.append(float(line.split(" ")[0]));
                bins.append(float(line.split(" ")[1]));
            else:
                if float(line.split(" ")[1]) != 13000:
                    bins.append(float(line.split(" ")[1]));
                else:
                    bins.append(bins[len(bins)-1]+500);

    bins.sort();
    binsArray = array('d',bins);

    start_histo = False;
    iBin = 0;
    for line in input_file:
        if "# BEGIN HISTO1D" in line:
            line = line.replace('\n','');
            print "Create histogram = ",line.split(" ")[3]
            histograms.append(ROOT.TH1F(line.split(" ")[3],"",len(binsArray)-1,binsArray));
        elif "# xlow" in line:
            start_histo= True;
            continue;
        elif "# END HISTO1D" in line:
            start_histo = False;
            iBin = 0;
            continue;
        elif start_histo == True:
            iBin = iBin+1;
            histograms[-1].SetBinContent(iBin,(float(line.split(" ")[2])));

    outputFile = ROOT.TFile(options.outputDIR+"/"+options.inputFile.replace(".dat",".root"),"RECREATE")
    for histo in histograms:
        histo.Write();
    outputFile.Close();
