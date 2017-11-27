import os
import glob
import math
from array import array
import sys
import time
import subprocess
import random

from optparse import OptionParser
from subprocess import Popen

############################################                                                                                                                                                           
#            Job steering                  #                                                                                                                                                           
############################################                                                                                                                                                           

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--inputProcessDIR', action="store", type="string", dest="inputProcessDIR",      default="",   help="Input directory where the template file are stored")
parser.add_option('--queque',          action="store", type="string", dest="queque",      default="8nh",   help="queque")
(options, args) = parser.parse_args()

if __name__ == '__main__':

    os.system("ls "+options.inputProcessDIR+" > file.list\n");
    file = open("file.list","r");    
    for lines in file:
        entry = lines.replace("\n","");
        path = entry.split("/")        
        command = "./submit_gridpack_generation.sh 1200 1200 "+options.queque+" "+path[len(path)-1]+" "+options.inputProcessDIR+"/"+entry+" "+options.queque;
        os.system(command);
            
    os.system("rm file.list");
