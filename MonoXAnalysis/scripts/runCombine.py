# python scripts/runCombine.py --computeAsymptotic --batchMode --queque cmscaf1nh --datacardDIR cards/monoV/MonoJMonoV/Shapes/  --workspaceDIR macros/monoV/workspaces/ --outputDIR macros/monoV/workspaces/AsymptoticLimit --jobDIR JOB --category 2 --Interaction Vector
#! /usr/bin/env python
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
from ROOT import gROOT, gStyle, gSystem
from collections import defaultdict

############################################
#            Job steering                  #
############################################

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

#### basic allowed methods
parser.add_option('--computeAsymptotic',        action='store_true', dest='computeAsymptotic',        default=False, help='run asymptotic limits')
parser.add_option('--computeProfileLikelihood', action='store_true', dest='computeProfileLikelihood', default=False, help='run profile likelihood')
parser.add_option('--maximumLikelihoodFit',     action='store_true', dest='maximumLikelihoodFit',     default=False, help='run max Likelihood fit inside combination tool')
parser.add_option('--generateOnly',             action='store_true', dest='generateOnly',             default=False, help='only generation with combiner')
parser.add_option('--makeLikelihoodScan',       action='store_true', dest='makeLikelihoodScan',       default=False, help='run the likelihood scan on mu (1D scan)')

##### submit jobs to condor, lxbatch and hercules 
parser.add_option('--batchMode',       action='store_true', dest='batchMode',      default=False, help='to run jobs on lxbatch')
parser.add_option('--queque',          action="store",      type="string",         dest="queque", default="cmscaf1nh")
parser.add_option('--submit',          action="store_true", dest='submit',         default=False, help='submit jobs')

##### other basci options for all the methods 
parser.add_option('--datacardDIR',  action="store", type="string", dest="datacardDIR",  default="", help='directory where to find the basic datacard templates')
parser.add_option('--workspaceDIR', action="store", type="string", dest="workspaceDIR", default="", help='directory with the workspaces .. submuit jobs for all of them')
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",    default="", help='output dir where to store .root files from combine')
parser.add_option('--jobDIR',       action="store", type="string", dest="jobDIR",       default="JOB", help='directory where to create jobs')

###### options for Bias test in the combination tool
parser.add_option('--nToys',                 action="store", type="int", dest="nToys", default=0,  help="number of toys to generate")
parser.add_option('--submitOneJobPerToy',    action="store_true",  dest="submitOneJobPerToy",        default=False,  help="one job one toy")
parser.add_option('--injectSignal',          action="store", type=float,    dest="injectSignal",  default=0., help='inject a signal in the toy generation')
parser.add_option('--runBlind',              action="store", type=int,      dest="runBlind",      default=0,  help='runBlind')
parser.add_option('--toyGeneration',         action="store", type=int,      dest="toyGeneration", default=0,  help='toyGeneration: 0 means toysNoSystematics, 1 means toys in a frequentis way, 2 standard toys')

###### input set of variables : 
parser.add_option('--rMin',          action="store", type=float, dest="rMin", default=0)
parser.add_option('--rMax',          action="store", type=float, dest="rMax", default=10)
parser.add_option('--category',      action="store", type=int,   dest="category", default=2, help="category = 0 means mono-j+mono-v combo, 1 means mono-j and 2 means mono-V")
parser.add_option('--Interaction',   action="store", type="string", dest="Interaction", default="all")
parser.add_option('--DMMass',        action="store", type="string", dest="DMMass"     , default="")
parser.add_option('--MedMass',       action="store", type="string", dest="MedMass"    , default="")


(options, args) = parser.parse_args()

##########################################
###### Submit batch job for combine ######
##########################################

def submitBatchJobCombine(command, fn, subname, outputDIR):

    currentDir = os.getcwd();    
    
    # create a dummy bash/csh
    outScript = open(fn,"w");

    outScript.write('#!/bin/bash \n');
    outScript.write('cd '+currentDir+'\n');
    outScript.write('eval `scram runtime -sh`'+'\n');
    outScript.write('cd - \n');
    outScript.write('scp '+currentDir+"/datacard*txt ./\n");
    outScript.write('scp '+currentDir+"/workspace*root ./\n");
    outScript.write(command+'\n');
    outScript.write("scp higgsCombine*"+subname+"* "+outputDIR+'\n');
    outScript.write("scp mlfit*"+subname+"* "+outputDIR+'\n');    
    outScript.write("rm rootstats* "+'\n');
    outScript.close();

    os.system("chmod 777 "+currentDir+"/"+fn);

    name = fn.replace(".sh","");

    if options.submit:
        if options.queque!="" :
            os.system("bsub -q "+options.queque+" -o "+currentDir+"/"+name+".log -e "+currentDir+"/"+name+".err "+fn);
        else: 
            os.system("bsub -q 1nh -o "+currentDir+"/"+name+".log -e "+currentDir+"/"+name+".err "+fn);


##################################
########### Main Code ############
##################################    

if __name__ == '__main__':

    print "###### start combine analysis ########";
    originalDIR = os.getcwd();

    if options.datacardDIR == "":
        sys.exit("provide a datacard directory");
    if options.workspaceDIR == "":
        sys.exit("provide a workspace directory");
        
    print "mkdir -p "+options.outputDIR;
    os.system("mkdir -p "+options.outputDIR);

    if options.batchMode:            
        print "mkdir -p "+options.jobDIR
        os.system("mkdir -p "+options.jobDIR);

    os.system("eval `scramv1 runtime -sh`"); ## to set the right root for combine
    
    rMin = options.rMin ; 
    rMax = options.rMax ;


    ##### make the workspace list
    command = "ls "+options.workspaceDIR+" | grep workspace | grep root ";
    if options.Interaction != "all":
        command += " | grep "+options.Interaction;
    if options.MedMass != "" and options.DMMass != "":
        command += " | grep _"+options.MedMass+"_"+options.DMMass;

    suffix = ""
    if options.category == 1:
        suffix =  "_MJ"
    elif options.category == 2:
        suffix =  "_MV"
        
    os.system(command+" > file_temp"+suffix);

    ## open the file
    fileWS = open("file_temp"+suffix,'r')
    for line_tmp in fileWS:
        os.chdir(originalDIR);
        line = line_tmp.rstrip()
        if options.category == 0 and "_MV_" in line : continue; ## already done
        if line == "" or line == '\n': continue;
        if options.category == 1 and "_MV_" in line : continue;
        if options.category == 2 and "_MJ_" in line : continue;
        if options.category != 0:
            dirName = line.replace("workspace_","").replace(".root","")
        else:
            dirName = line.replace("workspace_","").replace(".root","").replace("MJ_","").replace("MV_","")
            
        ## create single job dir
        os.system("mkdir -p "+options.jobDIR+"/job_"+dirName)
        ## copy the workspace
        os.system("scp "+options.workspaceDIR+"/"+line+" "+options.jobDIR+"/job_"+dirName);
        if options.category == 0:
            os.system("scp "+options.workspaceDIR+"/"+line.replace("MJ","MV")+" "+options.jobDIR+"/job_"+dirName);

        ## copy the datacards
        os.system("scp "+options.datacardDIR+"/*.txt "+options.jobDIR+"/job_"+dirName);
        ## set working dir in the job dir
        os.chdir(options.jobDIR+"/job_"+dirName);
        ## sobstitute the workspace name in the cards        
        os.system("ls | grep datacard | grep txt > datacard_list"+suffix)
        fileDC = open("datacard_list"+suffix,'r')
        dcList = [];
        for dc_temp in fileDC:            
            dc = dc_temp.rstrip();
            dcList.append(dc);
            if options.category == 1 and "_monoV_" in dc : continue;
            if options.category == 2 and "_monoJ_" in dc : continue;

            os.system("sed -i -e 's/workspace_MJ.root/"+line+"/g' "+dc);
            os.system("sed -i -e 's/workspace_MV.root/"+line.replace("MJ","MV")+"/g' "+dc);
                    
            ### combine cards
            if options.category == 1:
                command = "combineCards.py datacard_SR_monoJ.txt datacard_ZM_monoJ.txt datacard_WM_monoJ.txt datacard_GJ_monoJ.txt datacard_ZE_monoJ.txt datacard_WE_monoJ.txt" 
                if "datacard_TM_monoJ.txt" in dcList:
                    command += "datacard_TM_monoJ.txt"
                if "datacard_TE_monoJ.txt" in dcList:
                    command += "datacard_TE_monoJ.txt"
                
                command += " > datacard_COMB_monoJ.txt"
                dcList.append("datacard_COMB_monoJ.txt");
                os.system(command);

            elif options.category == 2:
                command = "combineCards.py datacard_SR_monoV.txt datacard_ZM_monoV.txt datacard_WM_monoV.txt datacard_GJ_monoV.txt datacard_ZE_monoV.txt datacard_WE_monoV.txt"
                if "datacard_TM_monoV.txt" in dcList:
                    command += "datacard_TM_monoV.txt"
                if "datacard_TE_monoV.txt" in dcList:
                    command += "datacard_TE_monoV.txt"
                
                command += " > datacard_COMB_monoV.txt"
                dcList.append("datacard_COMB_monoV.txt");
                os.system(command);

            else:

                command = "combineCards.py datacard_SR_monoJ.txt datacard_SR_monoV.txt datacard_ZM_monoJ.txt datacard_ZM_monoV.txt datacard_WM_monoJ.txt datacard_WM_monoV.txt datacard_GJ_monoJ.txt datacard_GJ_monoV.txt datacard_ZE_monoJ.txt datacard_ZE_monoV.txt datacard_WE_monoJ.txt datacard_WE_monoV.txt"
                if "datacard_TM_monoJ.txt" in dcList:
                    command += "datacard_TM_monoJ.txt"
                if "datacard_TM_monoV.txt" in dcList:
                    command += "datacard_TM_monoV.txt"
                if "datacard_TE_monoJ.txt" in dcList:
                    command += "datacard_TE_monoJ.txt"
                if "datacard_TE_monoV.txt" in dcList:
                    command += "datacard_TE_monoV.txt"

                command += " > datacard_COMB.txt"
                dcList.append("datacard_COMB.txt");
                os.system(command);

        os.system("rm datacard_list"+suffix)

        option = "";
        if options.category == 0:
            name = line.replace("workspace_","").replace(".root","").replace("MJ_","")
        else:
            name = line.replace("workspace_","").replace(".root","")

        ######## start to run combine options
        if options.generateOnly:

           if options.toyGeneration == 0:
               option += "--toysNoSystematics";
           elif options.toyGeneration == 1:
               option += "--toysFrequentist";

           command = "combine -M GenerateOnly --saveToys -s -1 -n %s -d %s -t %d --expectSignal=%d %s "%(name,dcList[len(dcList)-1],options.nToys,options.injectSignal,option);
           print "runCmmd ",runCmmd;
           if options.batchMode:
               fn = "job_Combine_generateOnly_%s_nToys_%d.sh"%(name,options.nToys);               
               submitBatchJobCombine(runCmmd,fn,name,originalDIR+"/"+options.outputDIR);
           else: 
               os.system(runCmmd);
               os.system("mv higgsCombine* "+originalDIR+"/"+options.outputDIR);   
           continue ;

        ########################################
        elif options.maximumLikelihoodFit == 1:
                                    
            if options.nToys == 0 : 
                runCmmd = "combine -M MaxLikelihoodFit --rMin %f --rMax %f --saveNormalizations --saveWithUncertainties --saveShapes -n %s -d  %s --saveNLL"%(rMin,rMax,name,dcList[len(dcList)-1]);
                print "runCmmd ",runCmmd;
                if options.batchMode:
                    fn = "job_Combine_MaxLikelihoodFit_%s.sh"%(name);
                    submitBatchJobCombine(runCmmd,fn,name,originalDIR+"/"+options.outputDIR);
                else:   
                    os.system(runCmmd);
                    os.system("mv higgsCombine* "+originalDIR+"/"+options.outputDIR);   
                    os.system("mv mlfit* "+originalDIR+"/"+options.outputDIR);   
                   
                continue;

            elif options.nToys != 0  :

                if options.toyGeneration == 0:
                    option += "--toysNoSystematics";
                elif options.toyGeneration == 1:
                    option += "--toysFrequentist";
                    
                if not options.submitOneJobPerToy:
                    runCmmd =  "combine -M MaxLikelihoodFit --rMin %f --rMax %f --saveNormalizations --saveWithUncertainties --saveShapes --saveToys -s -1 -n %s -d %s -t %d --expectSignal=%d --saveNLL %s"%(rMin,rMax,name,dcList[len(dcList)-1],options.nToys,options.injectSignal,option);
                    print "runCmmd ",runCmmd;
                    if options.batchMode:
                        fn = "job_Combine_MaxLikelihoodFit_%s_nToys_%d.sh"%(name,options.nToys);
                        submitBatchJobCombine(runCmmd,fn,name,originalDIR+"/"+options.outputDIR);
                    else: 
                        os.system(runCmmd);
                        os.system("mv higgsCombine* "+originalDIR+"/"+options.outputDIR);   
                        os.system("mv mlfit* "+originalDIR+"/"+options.outputDIR);   
                        
                    continue ;

                else:

                    for itoy in range(options.nToys):
                        outname = "iToy_%d"%(itoy);
                        runCmmd =  "combine -M MaxLikelihoodFit  --rMin %f --rMax %f --saveNormalizations --saveWithUncertainties  --saveShapes --saveToys -s -1 -n %s  -d %s  -t 1 --expectSignal=%d %s"%(rMin,rMax,name,dcList[len(dcList)-1],options.injectSignal,option);
                        print "runCmmd ",runCmmd;
                        if options.batchMode:
                            fn = "job_Combine_MaxLikelihoodFit_%s.sh"%(name);
                            submitBatchJobCombine(runCmmd,fn,name,originalDIR+"/"+options.outputDIR);
                        else: 
                            os.system(runCmmd);
                            os.system("mv higgsCombine* "+originalDIR+"/"+options.outputDIR);   
                            os.system("mv mlfit* "+originalDIR+"/"+options.outputDIR);   

                        continue ;
                      
       ########################################
        elif options.computeAsymptotic == 1 :

            if options.runBlind == 1 : option += "--run blind"            
            if options.nToys == 0 : ## run on data               
                runCmmd = "combine -M Asymptotic  --minosAlgo stepping -n %s -d %s  -H ProfileLikelihood %s "%(name,dcList[len(dcList)-1],option)
            else:
                if options.toyGeneration == 0:
                    option += "--toysNoSystematics";
                elif options.toyGeneration == 1:
                    option += "--toysFrequentist";

                runCmmd = "combine -M Asymptotic  --minosAlgo stepping -n %s -d %s -s -1 --expectSignal=%d -t %d  -H ProfileLikelihood %s"%(name,dcList[len(dcList)-1],options.injectSignal,options.nToys,option)

            print "runCmmd ",runCmmd ;
            if options.batchMode:
                fn = "job_Combine_Asymptotic_%s.sh"%(name);
                submitBatchJobCombine(runCmmd,fn,name,originalDIR+"/"+options.outputDIR);
            else: 
                os.system(runCmmd);
                os.system("mv higgsCombine* "+originalDIR+"/"+options.outputDIR);
                continue;
                                                
         ########################################
        elif options.computeProfileLikelihood == 1 : 

           if options.nToys == 0:
               runCmmd = "combine -M ProfileLikelihood --signif  -n %s  -d %s"%(name,dcList[len(dcList)-1]);
           else:    
               if options.toyGeneration == 0:
                   option += "--toysNoSystematics";
               elif options.toyGeneration == 1:
                   option += "--toysFrequentist";

               runCmmd = "combine -M ProfileLikelihood --signif  -n %s  -d %s -t %d --expectSignal=%d -s -1 %s"%(name,dcList[len(dcList)-1],options.nToys,options.injectSignal,option);

           print "runCmmd ",runCmmd;
                          
           if options.batchMode:
               fn = "job_Combine_ProfileLikelihood_%s.sh"%(name);
               submitBatchJobCombine(runCmmd,fn,name,originalDIR+"/"+options.outputDIR);
           else:
               os.system(runCmmd);
               os.system("mv higgsCombine* "+originalDIR+"/"+options.outputDIR);

        ########################################
        elif options.makeLikelihoodScan == 1:

            if options.category == 2 or "10000" in name or "5000" in name or "2000" in name:
                rMax = 3.;
            else:
                rMax = 1.;

            if options.nToys == 0:          
                runCmmd = "combine -M MultiDimFit -n %s  -d %s --algo=grid --points=100 --setPhysicsModelParameterRanges r=%f,%f"%(name,dcList[len(dcList)-1],rMin,rMax);
                
            else:
                if options.toyGeneration == 0:
                    option += "--toysNoSystematics";
                elif options.toyGeneration == 1:
                    option += "--toysFrequentist";
                    
                runCmmd = "combine -M MultiDimFit -n %s  -d %s --algo=grid --points=100 --setPhysicsModelParameterRanges r=%f,%f -s -1 --expectSignal=%d -t %d %s"%(name,dcList[len(dcList)-1],rMin,rMax,options.injectSignal,options.nToys,option);

                
            print "runCmmd ",runCmmd;
            if options.batchMode:
                fn = "job_Combine_LikelihoodScan_%s.sh"%(name);
                submitBatchJobCombine(runCmmd,fn,name,originalDIR+"/"+options.outputDIR);
            else:
                os.system(runCmmd);
                os.system("mv higgsCombine*MultiDimFit* "+originalDIR+"/"+options.outputDIR);


    os.chdir(originalDIR);
    os.system("rm file_temp"+suffix)
