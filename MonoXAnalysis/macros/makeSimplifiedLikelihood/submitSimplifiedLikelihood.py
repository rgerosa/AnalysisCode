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
from itertools import product

############################################                                                                                                                                                
#            Job steering                  #                                                                                                                                                      
############################################                                                                                                                                                       

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

## parse files                                                                                                                                                                                     
parser.add_option('--inputWorkspace',  action="store", type="string", dest="inputWorkspace",  default="",   help="input ROOT file with workspace with data")
parser.add_option('--inputMLFitFile',  action="store", type="string", dest="inputMLFitFile",  default="",   help="input ROOT file with post-fit backgrond and covariance matrix")
parser.add_option('--inputSignalDIR',  action="store", type="string", dest="inputSignalDIR",  default="",   help="input director where signal templates are located")
parser.add_option('--mediatorType',    action="store", type="string", dest="mediatorType",    default="",   help="mediator type: vector, axial, scalar, pseudoscalar")
parser.add_option('--category',        action="store", type="string", dest="category",        default="",   help="mediator type: category for limit -> monojet or monov")
parser.add_option('--outputDIR',       action="store", type="string", dest="outputDIR",       default="",   help="output DIR")
##  for submitting jobs in lxbatch                                                                                                                                                                    
parser.add_option('--batchMode',    action="store_true",           dest="batchMode",                  help="batchMode")
parser.add_option('--jobDIR',       action="store", type="string", dest="jobDIR",  default="",        help="directory for job")
parser.add_option('--queque',       action="store", type="string", dest="queque",  default="",        help="queque for LSF")
parser.add_option('--submit',       action="store_true",           dest="submit",                     help="submit")

(options, args) = parser.parse_args()

mMED = [50,60,70,80,90,100,125,150,175,200,300,325,400,525,600,725,800,925,1000,1125,1200,1325,1400,1525,1600,1725,1800,1925,2000,2500,3000,3500,4000,5000];
mDM  = [1,5,10,25,50,100,150,200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000];

if __name__ == '__main__':

    monojetFile = "";
    monowFile = "";
    monozFile = "";
    index = "";

    if options.mediatorType == "Vector" or options.mediatorType == "vector":
        monojetFile = options.inputSignalDIR+"/MonoJ_800_0.25_cat"+options.category+"_13TeV_v1.root";
        monowFile = options.inputSignalDIR+"/MonoW_800_0.25_cat"+options.category+"_13TeV_v1.root";
        monozFile = options.inputSignalDIR+"/MonoZ_800_0.25_cat"+options.category+"_13TeV_v1.root";
        index = "800";
    if options.mediatorType == "Axial"  or options.mediatorType == "axial":
        monojetFile = options.inputSignalDIR+"/MonoJ_801_0.25_cat"+options.category+"_13TeV_v1.root";
        monowFile = options.inputSignalDIR+"/MonoW_801_0.25_cat"+options.category+"_13TeV_v1.root";
        monozFile = options.inputSignalDIR+"/MonoZ_801_0.25_cat"+options.category+"_13TeV_v1.root";
        index = "801";
    if options.mediatorType == "Scalar" or options.mediatorType == "scalar":
        monojetFile = options.inputSignalDIR+"/MonoJ_805_1.0_cat"+options.category+"_13TeV_v1.root";
        monowFile = options.inputSignalDIR+"/MonoW_805_1.0_cat"+options.category+"_13TeV_v1.root";
        monozFile = options.inputSignalDIR+"/MonoZ_805_1.0_cat"+options.category+"_13TeV_v1.root";
        index = "805";
    if options.mediatorType == "Pseudoscalar" or options.mediatorType == "pseudoscalar":
        monojetFile = options.inputSignalDIR+"/MonoJ_805_1.0_cat"+options.category+"_13TeV_v1.root";
        monowFile = options.inputSignalDIR+"/MonoW_805_1.0_cat"+options.category+"_13TeV_v1.root";
        monozFile = options.inputSignalDIR+"/MonoZ_805_1.0_cat"+options.category+"_13TeV_v1.root";
        index = "806";

    
    currentDIR = os.getcwd();
    os.system("mkdir -p "+options.outputDIR);
    os.system("mkdir -p "+options.jobDIR);
    os.system("rm -r "+options.jobDIR+"/*");

    workspacedir = options.inputWorkspace.split("/")[0];
    mlfitdir = options.inputMLFitFile.split("/")[0];


    ### loop to submit jobs for different mass point
    monojetfile = ROOT.TFile(monojetFile,"READ");
    monojetws   = monojetfile.Get("combinedws");
    monowfile   = ROOT.TFile(monowFile,"READ");
    monowws   = monowfile.Get("combinedws");
    monozfile = ROOT.TFile(monozFile,"READ");
    monozws   = monozfile.Get("combinedws");

    for mmed in mMED:
        for mdm in mDM:

            expMassPoint = index+"%04d%04d"%(mmed,mdm);

            if options.category == "monojet":    
                if not monojetws.obj("monojet_signal_signal_"+expMassPoint) or not monowws.obj("monojet_signal_signal_"+expMassPoint) or not monozws.obj("monojet_signal_signal_"+expMassPoint):
                    continue; ### signal histogram not found
            elif options.category == "monov":
                if not monojetws.obj("monov_signal_signal_"+expMassPoint) or not monowws.obj("monov_signal_signal_"+expMassPoint) or not monozws.obj("monov_signal_signal_"+expMassPoint):
                    continue; ### signal histogram not found


            jobname = "job_"+options.mediatorType+"_cat"+options.category+"_"+expMassPoint;
            os.system("mkdir -p "+options.jobDIR+"/"+jobname);
            
            if options.category == "monojet":
                category = 1;
            elif options.category == "monov":
                category = 2;
            
    
            jobmacro = open('%s/%s/job.C'%(options.jobDIR,jobname),'w')
            jobmacro.write("{\n");
            jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/limitSimplifiedLikelihood.C\");\n");
            jobmacro.write("gROOT->ProcessLine(\""+"limitSimplifiedLikelihood(\\\"%s\\\",\\\"%s\\\",{\\\"%s\\\",\\\"%s\\\",\\\"%s\\\"},\\\"%s\\\",%d,\\\"%s\\\")"%(options.inputWorkspace,options.inputMLFitFile,monojetfile.GetName(),monowfile.GetName(),monozfile.GetName(),expMassPoint,category,options.outputDIR)+"\");\n");
            jobmacro.write("}\n");
            jobmacro.close();

            
            jobscript = open('%s/%s/job.sh'%(options.jobDIR,jobname),'w')
            jobscript.write('cd %s \n'%currentDIR)
            jobscript.write('eval ` scramv1 runtime -sh ` \n')
            jobscript.write('cd - \n')            

            jobscript.write('scp '+currentDIR+'/'+options.jobDIR+'/'+jobname+'/job.C ./ \n');
            jobscript.write('scp -r '+currentDIR+'/'+workspacedir+' ./ \n');
            jobscript.write('scp -r '+currentDIR+'/'+mlfitdir+' ./ \n');
            jobscript.write('scp -r '+currentDIR+'/'+options.inputSignalDIR+' ./ \n');
            jobscript.write('root -l -b -q job.C\n');
            jobscript.write('scp '+options.outputDIR+"/*"+expMassPoint+"*root "+currentDIR+"/"+options.outputDIR);
            os.system('chmod a+x %s/%s/job.sh'%(options.jobDIR,jobname))

            if options.submit:
                os.system('bsub -q %s -o %s/%s/%s/job.log -e %s/%s/%s/job.err %s/%s/%s/job.sh'%(options.queque,currentDIR,options.jobDIR,jobname,currentDIR,options.jobDIR,jobname,currentDIR,options.jobDIR,jobname))
