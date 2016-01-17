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
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",      default="",   help="input directory where files are located on eos")
parser.add_option('--nStreams',     action="store", type=int,      dest="nStreams",      default=3,    help="number of parallel streams")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",     default="",   help="director to be copied")
(options, args) = parser.parse_args()

if __name__ == '__main__':

    os.system('source /afs/cern.ch/project/eos/installation/cms/etc/setup.sh')
    os.system('mkdir -p '+options.outputDIR);

    if 'eos/cms/' in options.inputDIR or '/eos/cms/' in options.inputDIR:
        options.inputDIR = options.inputDIR.replace('/eos/cms/','');
        options.inputDIR = options.inputDIR.replace('eos/cms/','');
    
    os.system('xrdcp -r -S '+str(options.nStreams)+' root://eoscms.cern.ch//eos/cms'+options.inputDIR+' '+options.outputDIR);
        

