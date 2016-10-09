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
parser.add_option('--calculateExpSigma', action="store", type=int,    dest="calculateExpSigma", default=1,   help="to calculate expected for quantiles = 0.025,0.16,0.84,0.975")
parser.add_option('--doLikelihoodScan', action="store", type=int,     dest="doLikelihoodScan",  default=1,   help="in order to perform a likelihood scan")
parser.add_option('--skipCorrelations', action="store", type=int,     dest="skipCorrelations",  default=0,   help="to skip correlation across bins")

##  for submitting jobs in lxbatch                                                                                                                                                                
parser.add_option('--batchMode',    action="store_true",           dest="batchMode",                  help="batchMode")
parser.add_option('--jobDIR',       action="store", type="string", dest="jobDIR",  default="",        help="directory for job")
parser.add_option('--queque',       action="store", type="string", dest="queque",  default="",        help="queque for LSF")
parser.add_option('--submit',       action="store_true",           dest="submit",                     help="submit")

(options, args) = parser.parse_args()

mMED_v = [50,60,70,80,90,100,125,150,175,200,300,325,400,525,600,725,800,925,1000,1125,1200,1325,1400,1525,1600,1725,1800,1925,2000,2500,3000,3500,4000,5000];
mDM_v  = [1,5,10,25,50,100,150,200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000];

mMED_av = [50,60,70,80,90,100,125,150,175,200,300,325,400,525,600,725,800,925,1000,1125,1200,1325,1400,1525,1600,1725,1800,1925,2000,2500,3000,3500,4000,5000];
mDM_av  = [1,5,10,25,50,100,150,200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000];

mMED_s  = [10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,325,400,525,600,725,800,925,1000,1125,1200,1325,1400,1525,1600,1725,1800,1925,2000,2500,3000,3500,4000,5000]
mDM_s   = [1,5,10,25,50,100,150,200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000]

mMED_ps  = [10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,325,400,525,600,725,800,925,1000,1125,1200,1325,1400,1525,1600,1725,1800,1925,2000,2500,3000,3500,4000,5000]
mDM_ps   = [1,5,10,25,50,100,150,200,300,400,500,600,700,800,900,1000,1250,1500,1750,2000]


if __name__ == '__main__':

    monojetFile_mj = "";
    monowFile_mj = "";
    monozFile_mj = "";
    monojetFile_mv = "";
    monowFile_mv = "";
    monozFile_mv = "";
    index = "";

    ### signal file list
    if options.mediatorType == "Vector" or options.mediatorType == "vector" :

        monojetFile_mj = options.inputSignalDIR+"/MonoJ_800_0.25_catmonojet_13TeV_v1.root";
        monowFile_mj = options.inputSignalDIR+"/MonoW_800_0.25_catmonojet_13TeV_v1.root";
        monozFile_mj = options.inputSignalDIR+"/MonoZ_800_0.25_catmonojet_13TeV_v1.root";

        monojetFile_mv = options.inputSignalDIR+"/MonoJ_800_0.25_catmonov_13TeV_v1.root";
        monowFile_mv   = options.inputSignalDIR+"/MonoW_800_0.25_catmonov_13TeV_v1.root";
        monozFile_mv   = options.inputSignalDIR+"/MonoZ_800_0.25_catmonov_13TeV_v1.root";

        index = "800";

    elif options.mediatorType == "Axial"  or options.mediatorType == "axial" :

        monojetFile_mj = options.inputSignalDIR+"/MonoJ_801_0.25_catmonojet_13TeV_v1.root";
        monowFile_mj = options.inputSignalDIR+"/MonoW_801_0.25_catmonojet_13TeV_v1.root";
        monozFile_mj = options.inputSignalDIR+"/MonoZ_801_0.25_catmonojet_13TeV_v1.root";

        monojetFile_mv = options.inputSignalDIR+"/MonoJ_801_0.25_catmonov_13TeV_v1.root";
        monowFile_mv   = options.inputSignalDIR+"/MonoW_801_0.25_catmonov_13TeV_v1.root";
        monozFile_mv   = options.inputSignalDIR+"/MonoZ_801_0.25_catmonov_13TeV_v1.root";

        index = "801";

    elif options.mediatorType == "Scalar" or options.mediatorType == "scalar" :

        monojetFile_mj = options.inputSignalDIR+"/MonoJ_805_1.0_catmonojet_13TeV_v1.root";
        monowFile_mj   = options.inputSignalDIR+"/MonoW_805_1.0_catmonojet_13TeV_v1.root";
        monozFile_mj   = options.inputSignalDIR+"/MonoZ_805_1.0_catmonojet_13TeV_v1.root";

        monojetFile_mv = options.inputSignalDIR+"/MonoJ_805_1.0_catmonov_13TeV_v1.root";
        monowFile_mv   = options.inputSignalDIR+"/MonoW_805_1.0_catmonov_13TeV_v1.root";
        monozFile_mv   = options.inputSignalDIR+"/MonoZ_805_1.0_catmonov_13TeV_v1.root";

        index = "805";

    elif options.mediatorType == "Pseudoscalar" or options.mediatorType == "pseudoscalar" :

        monojetFile_mj  = options.inputSignalDIR+"/MonoJ_806_1.0_catmonojet_13TeV_v1.root";
        monowFile_mj    = options.inputSignalDIR+"/MonoW_806_1.0_catmonojet_13TeV_v1.root";
        monozFile_mj    = options.inputSignalDIR+"/MonoZ_806_1.0_catmonojet_13TeV_v1.root";

        monojetFile_mv = options.inputSignalDIR+"/MonoJ_806_1.0_catmonov_13TeV_v1.root";
        monowFile_mv   = options.inputSignalDIR+"/MonoW_806_1.0_catmonov_13TeV_v1.root";
        monozFile_mv   = options.inputSignalDIR+"/MonoZ_806_1.0_catmonov_13TeV_v1.root";

        index = "806";

    
    currentDIR = os.getcwd();
    os.system("mkdir -p "+options.outputDIR);
    os.system("mkdir -p "+options.jobDIR);
    os.system("rm -r "+options.jobDIR+"/*"+options.category+"*");
    
    workspacedir = options.inputWorkspace.split("/")[0];
    mlfitdir = options.inputMLFitFile.split("/")[0];

    ### loop to submit jobs for different mass point
    if monojetFile_mj != "":
        monojetfile_mj = ROOT.TFile(monojetFile_mj,"READ");
        monojetws_mj   = monojetfile_mj.Get("combinedws");
    if monowFile_mj != "":
        monowfile_mj   = ROOT.TFile(monowFile_mj,"READ");
        monowws_mj     = monowfile_mj.Get("combinedws");
    if monozFile_mj != "":
        monozfile_mj = ROOT.TFile(monozFile_mj,"READ");
        monozws_mj   = monozfile_mj.Get("combinedws");

    if monojetFile_mv != "":
        monojetfile_mv = ROOT.TFile(monojetFile_mv,"READ");
        monojetws_mv   = monojetfile_mv.Get("combinedws");
    if monowFile_mv != "":
        monowfile_mv   = ROOT.TFile(monowFile_mv,"READ");
        monowws_mv   = monowfile_mv.Get("combinedws");
    if monozFile_mv != "":
        monozfile_mv = ROOT.TFile(monozFile_mv,"READ");
        monozws_mv   = monozfile_mv.Get("combinedws");

    ###############
    mMED = [];
    mDM  = [];

    if options.mediatorType == "vector":
        mMED = mMED_v;
        mDM  = mDM_v;
    elif options.mediatorType == "axial":
        mMED = mMED_av;
        mDM  = mDM_av;
    elif options.mediatorType == "scalar":
        mMED = mMED_s;
        mDM  = mDM_s;
    elif options.mediatorType == "pseudoscalar":
        mMED = mMED_ps;
        mDM  = mDM_ps;

    ##########################
    for mmed in mMED:
        for mdm in mDM:

            expMassPoint = index+"%04d%04d"%(mmed,mdm);
            if options.category == "monojet" or options.category == "combined" or options.category == "supercombo" and options.mediatorType != "pseudoscalar":    
                if not monojetws_mj.obj("monojet_signal_signal_"+expMassPoint) or not monowws_mj.obj("monojet_signal_signal_"+expMassPoint) or not monozws_mj.obj("monojet_signal_signal_"+expMassPoint):
                    continue; ### signal histogram not found
            elif options.category == "monov" or options.category == "combined" or options.category == "supercombo" and options.mediatorType != "pseudoscalar":
                if not monojetws_mv.obj("monov_signal_signal_"+expMassPoint) or not monowws_mv.obj("monov_signal_signal_"+expMassPoint) or not monozws_mv.obj("monov_signal_signal_"+expMassPoint):
                    continue; ### signal histogram not found
            elif options.category == "monojet" or options.category == "combined" or options.category == "supercombo" and options.mediatorType == "pseudoscalar":
                if not monojetws_mj.obj("monojet_signal_signal_"+expMassPoint):
                    continue


            jobname = "job_"+options.mediatorType+"_cat"+options.category+"_"+expMassPoint;
            os.system("mkdir -p "+options.jobDIR+"/"+jobname);
            
            if options.category == "monojet":
                category = 1;
            elif options.category == "monov":
                category = 2;
            elif options.category == "combined":
                category = 0;
            elif options.category == "supercombo":
                category = -1;
            
    
            jobmacro = open('%s/%s/job.C'%(options.jobDIR,jobname),'w')
            jobmacro.write("{\n");
            jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/limitSimplifiedLikelihood.C\");\n");
            if options.mediatorType != "pseudoscalar":
                jobmacro.write("gROOT->ProcessLine(\""+"limitSimplifiedLikelihood(\\\"%s\\\",\\\"%s\\\",{\\\"%s\\\",\\\"%s\\\",\\\"%s\\\"},{\\\"%s\\\",\\\"%s\\\",\\\"%s\\\"},\\\"%s\\\",%d,\\\"%s\\\",%d,%d,%d)"%(options.inputWorkspace,options.inputMLFitFile,monojetfile_mj.GetName(),monowfile_mj.GetName(),monozfile_mj.GetName(),monojetfile_mv.GetName(),monowfile_mv.GetName(),monozfile_mv.GetName(),expMassPoint,category,options.outputDIR,options.calculateExpSigma,options.doLikelihoodScan,options.skipCorrelations)+"\");\n");
            else:                
                jobmacro.write("gROOT->ProcessLine(\""+"limitSimplifiedLikelihood(\\\"%s\\\",\\\"%s\\\",{\\\"%s\\\"},{\\\"%s\\\"},\\\"%s\\\",%d,\\\"%s\\\",%d,%d,%d)"%(options.inputWorkspace,options.inputMLFitFile,monojetfile_mj.GetName(),monojetfile_mv.GetName(),expMassPoint,category,options.outputDIR,options.calculateExpSigma,options.doLikelihoodScan,options.skipCorrelations)+"\");\n");

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
