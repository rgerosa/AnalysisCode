#! /usr/bin/env python
import os
import glob
import math
import sys
import time
import subprocess
import ROOT

from array import array
from optparse import OptionParser
from subprocess import Popen
from ROOT import gROOT, gStyle, gSystem
from collections import defaultdict

############################################
#            Job steering                  #
############################################

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
##### Basic options for to effectively run combine
parser.add_option('--datacardDIR',  action="store", type="string", dest="datacardDIR",  default="",    help='directory where to find the basic datacard templates')
parser.add_option('--workspaceDIR', action="store", type="string", dest="workspaceDIR", default="",    help='directory with the workspaces: one for each mass point')
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",    default="",    help='output dir where to store .root files from combine')
parser.add_option('--jobDIR',       action="store", type="string", dest="jobDIR",       default="JOB", help='directory where to create jobs')
##### option to run limits on a specific set of workspaces
parser.add_option('--category',     action="store", type=int,      dest="category",     default=1,     help="category = 0 means mono-j+mono-v combo, 1 means mono-j and 2 means mono-V")
parser.add_option('--interaction',   action="store", type="string", dest="Interaction", default="all", help="interaction type for DM searches")
parser.add_option('--dmMass',        action="store", type="string", dest="dmMass"     , default="",    help="DM mass value")
parser.add_option('--medMass',       action="store", type="string", dest="medMass"    , default="",    help="Mediator mass value")

##### submit jobs to condor, lxbatch and hercules 
parser.add_option('--batchMode',  action='store_true', dest='batchMode', default=False, help='to run jobs on lxbatch')
parser.add_option('--queque',     action="store",      dest="queque",    default="cmscaf1nh", type="string",  help='name of the lxbatch queque for jobs')
parser.add_option('--submit',     action="store_true", dest='submit',    default=False, help='submit jobs')

##### toy infomation
parser.add_option('--ntoys',              action="store",       type="int", dest="ntoys",         default=0,     help="number of toys to generate for each single mass point, -1 means Asimov Toy")
parser.add_option('--submitOneJobPerToy', action="store_true",  dest="submitOneJobPerToy",        default=False, help="one job one toy --> useful in case of bias studies")
parser.add_option('--expectSignal',       action="store",       type=float, dest="expectSignal",  default=0.,    help='inject a signal in the toy generation')
parser.add_option('--toyGeneration',      action="store",       type=int,   dest="toyGeneration", default=-1,    help='toyGeneration: 0 mean frequentist toys, 1 means a priori toy without sys, other')
parser.add_option('--crossedToys',        action="store_true",  dest="crossedToys", default=False, help='Generate dataset with a model and fit with another oen')
parser.add_option('--inputGeneratedDataset', action="store",    type="string", dest="inputGeneratedDataset", default="", help='directory with root files with toys')

#### basic allowed methods and options 
parser.add_option('--Asymptotic',  action='store_true', dest='Asymptotic',  default=False, help='run asymptotic limits')
parser.add_option('--runBlind',    action="store_true", dest="runBlind",    default=False, help='runBlind mode')
parser.add_option('--runExpected', action="store_true", dest="runExpected", default=False, help='run expected limit only')

#####
parser.add_option('--generateOnly', action='store_true', dest='generateOnly',  default=False, help='only generation with combiner')

####
parser.add_option('--MaxLikelihoodFit', action='store_true', dest='MaxLikelihoodFit', default=False, help='run max Likelihood fit inside combination tool: by default --saveNLL, --saveShapes, --saveNormalizations, --saveWithUncertainties')

parser.add_option('--MultiDimFit',  action='store_true', dest='MaxLikelihoodFit', default=False, help='run max Likelihood fit inside combination tool')
parser.add_option('--rMin',  action='store', dest='rMin', type=float, default=0, help='rMin to be used by all the methods')
parser.add_option('--rMax',  action='store', dest='rMax', type=float, default=100, help='rMax to be used by all the methods')

#####
parser.add_option('--ProfileLikelihood',  action='store_true', dest='ProfileLikelihood', default=False, help='run profile likelihood method in combine')
parser.add_option('--pvalue',             action='store_true', dest='--pvalue',   default=False, help='run profile likelihood method in combine')

parser.add_option('--calculateStatOnly',  action='store_true', dest='calculateStatOnly', default=False, help='to calculate stat only limits or significance, needs a snapshot')
parser.add_option('--inputSnapShotDIR',   action='store', type="string", dest='inputSnapShotDIR', default=False, help='directory with snapshot files')

#####
parser.add_option('--makeLikelihoodScan', action='store_true', dest='makeLikelihoodScan',  default=False, help='run the likelihood scan on mu (1D scan)')

(options, args) = parser.parse_args()

##########################################
###### Submit batch job for combine ######
##########################################

def submitBatchJobCombine(command, fn, subname, outputDIR):

    currentDir = os.getcwd();    
    
    # create a dummy bash/csh file with the instruction for the job
    outScript = open(fn,"w");

    outScript.write('#!/bin/bash \n');
    outScript.write('cd '+currentDir+'\n');
    outScript.write('eval `scram runtime -sh`'+'\n');
    outScript.write('cd - \n');
    ### copy datacards and workspaces at the node
    outScript.write('scp '+currentDir+"/datacard*txt ./\n");
    outScript.write('scp '+currentDir+"/workspace*root ./\n");
    ### execute combine command
    outScript.write(command+'\n');
    ### copyt output files
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

    #### dump location of basic directories 
    print "###### start combine analysis ########";
    originalDIR  = os.getcwd();
    os.chdir(options.workspaceDIR)
    workspaceDIR = os.getcwd();
    os.chdir(originalDIR)

    #### basic checks
    if options.datacardDIR == "":
        sys.exit("provide a datacard directory");
    if options.workspaceDIR == "":
        sys.exit("provide a workspace directory");
        
    ### make datacard dir
    print "mkdir -p "+options.outputDIR;
    os.system("mkdir -p "+options.outputDIR);

    ### make job dir in case of batch mode
    if options.batchMode:            
        print "mkdir -p "+options.jobDIR
        os.system("mkdir -p "+options.jobDIR);

    os.system("eval `scramv1 runtime -sh`"); ## to set the right root for combine
    
    ## set rMin and rMax
    rMin = options.rMin ; 
    rMax = options.rMax ;

    ##### make the workspace list
    command = "ls "+options.workspaceDIR+" | grep workspace | grep root ";
    if options.Interaction != "all":
        command += " | grep "+options.interaction;
    if options.medMass != "" and options.dmMass != "":
        command += " | grep _"+options.medMass+"_"+options.dmMass;

    #### to handle workspace correctly
    suffix = ""
    if options.category == 1:
        suffix =  "_MJ"
    elif options.category == 2:
        suffix =  "_MV"
        
    os.system(command+" > file_temp"+suffix);

    ## open the file --> make the workspace list
    fileWS = open("file_temp"+suffix,'r')
    for line_tmp in fileWS:
        os.chdir(originalDIR);
        line = line_tmp.rstrip()
        if options.category == 0 and "_MV_" in line : continue; ## already done .. since there are two workspaces, avoid to submit twice the same job
        if line == "" or line == '\n': continue;
        if options.category == 1 and "_MV_" in line : continue; ## consider the right workspace for the chosen category
        if options.category == 2 and "_MJ_" in line : continue; ## consider the right workspace for the chosen category
        if options.category != 0:
            dirName = line.replace("workspace_","").replace(".root","")
        else:
            dirName = line.replace("workspace_","").replace(".root","").replace("MJ_","").replace("MV_","")
            
        ## create single job dir
        os.system("mkdir -p "+options.jobDIR+"/job_"+dirName)
        ## copy the workspace into the job dir
        os.system("scp "+options.workspaceDIR+"/"+line+" "+options.jobDIR+"/job_"+dirName);
        if options.category == 0:
            os.system("scp "+options.workspaceDIR+"/"+line.replace("MJ","MV")+" "+options.jobDIR+"/job_"+dirName); ## to ensure to copy both WS in combo case

        ## copy the datacards
        os.system("scp "+options.datacardDIR+"/*.txt "+options.jobDIR+"/job_"+dirName);
        ## set working dir in the job dir
        os.chdir(options.jobDIR+"/job_"+dirName);
        ## sobstitute the workspace name in the cards        
        os.system("ls | grep datacard | grep txt > datacard_list"+suffix)
        fileDC = open("datacard_list"+suffix,'r')
        dcList = [];
        command = "combineCards.py ";
        for dc_temp in fileDC:            
            dc = dc_temp.rstrip();
            dcList.append(dc);
            if options.category == 1 and "_monoV_" in dc : continue;
            if options.category == 2 and "_monoJ_" in dc : continue;
            ### put the right workspace name into the datacard
            os.system("sed -i -e 's/workspace_MJ.root/"+line+"/g' "+dc);
            os.system("sed -i -e 's/workspace_MV.root/"+line.replace("MJ","MV")+"/g' "+dc);
            command += dc+" "; ## to run the combine cards command

        ### execute the combine card command
        if options.category == 1:                
            command += " > datacard_COMB_monoJ.txt"
            dcList.append("datacard_COMB_monoJ.txt");
            os.system(command);

        elif options.category == 2:
            command += " > datacard_COMB_monoV.txt"
            dcList.append("datacard_COMB_monoV.txt");
            os.system(command);

        else:
            command += " > datacard_COMB.txt"
            dcList.append("datacard_COMB.txt");
            os.system(command);

        os.system("rm datacard_list"+suffix)

        ###### name to be used for the output combine file
        if options.category == 0:
            name = line.replace("workspace_","").replace(".root","").replace("MJ_","").replace("MV_","")
        else:
            name = line.replace("workspace_","").replace(".root","")

        ### options for toy generation
        toyOption = "";
        if options.toyGeneration == 0:
            toyOption += " --toysFrequentist ";
        elif options.toyGeneration == 1:
            toyOption += " --toysNoSystematics ";

        ### options for asym limits
        asymOpt = "";
        if options.runBlind and options.runExpected: sys.exit("Cannot run the asymptotic limits with --run blind and --run expected at the same time")
        if options.runBlind: asymOpt += " --run blind "            
        if options.runExpected: asymOpt += " --run expected "            
            
        ### options for profile likelihood
        profOpt = ""
        if options.pvalue :
            profOpt += " --pvalue ";
            
        ####### for a string to be used as mH value in the limit tree
        massPoint = "";
        if "Vector" in name:   massPoint += "1";
        elif "Axial" in name:  massPoint += "2";
        elif "Scalar" in name: massPoint += "3";
        elif "Pseudoscalar" in name or "PseudoScalar" in name: massPoint += "4";
        else: massPoint += "5";

        point = name.split("_");
        mMED  = "";
        mDM   = "";
        for substring in point:
            if substring == "MJ" or substring == "MV": continue;
            if mMED == "" :
                mMED = substring;
                continue;
            if mDM == "":
                mDM = substring;
                continue;
            
        if mMED == "":       sys.exit("mediator mass information not found in workspace name --> bad ");
        elif len(mMED) > 4:  sys.exit("mediator mass information too large --> bad ");
        elif len(mMED) == 3: mMED = "0"+mMED;
        elif len(mMED) == 2: mMED = "00"+mMED;
        elif len(mMED) == 1: mMED = "000"+mMED;
        
        massPoint += mMED;
        
        if mDM == "": mDM = "0000";
        elif len(mDM) > 4:  sys.exit("dark matter mass information too large --> bad ");
        elif len(mDM) == 3: mDM = "0"+mDM;
        elif len(mDM) == 2: mDM = "00"+mDM;
        elif len(mDM) == 3: mDM = "000"+mDM;

        massPoint += mDM;

        ######## performing only toys: saveToys allow to store the generated dataset in the output root file
        if options.generateOnly:
            if not options.submitOneJobPerToy or options.ntoys == -1: ### one job, many toys
                command = "combine -M GenerateOnly --saveToys -s -1 -n %s -d %s -t %d --expectSignal=%f %s -m %d"%(name,dcList[len(dcList)-1],options.ntoys,options.expectSignal,toyOption,int(massPoint));
                print "runCmmd ",runCmmd;
                if options.batchMode:
                    if options.ntoys != -1:
                        fn = "job_Combine_generateOnly_%s_ntoys_%d.sh"%(name,options.ntoys);               
                    else:
                        fn = "job_Combine_generateOnly_%s_ntoys_Asimov.sh"%(name);               
                    submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                else: 
                    os.system(runCmmd);
                    os.system("mv higgsCombine* "+workspaceDIR+"/"+options.outputDIR);   
                continue ;
            else: ## one job per toy
                for itoy in range(options.ntoys):
                    command = "combine -M GenerateOnly --saveToys -s -1 -n %s_%d -d %s -t 1 --expectSignal=%f %s -m %d"%(name+"_itoy_",itoy,dcList[len(dcList)-1],options.expectSignal,toyOption,int(massPoint));
                print "runCmmd ",runCmmd;
                if options.batchMode:
                    fn = "job_Combine_generateOnly_%s_itoy_%d.sh"%(name,itoy);               
                    submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                else: 
                    os.system(runCmmd);
                    os.system("mv higgsCombine* "+workspaceDIR+"/"+options.outputDIR);   
                continue ;            

        #########################
        elif options.Asymptotic and not options.calculateStatOnly:
            
            if options.ntoys == 0 : ## run on data --> option contains type of run and type of toys             
                runCmmd = "combine -M Asymptotic -n %s -d %s --rMin %f --rMax %f -m %d %s"%(name,dcList[len(dcList)-1],rMin,rMax,int(massPoint),asymOpt)
                print "runCmmd ",runCmmd ;
                if options.batchMode:
                    fn = "job_Combine_Asymptotic_%s.sh"%(name);
                    submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                else: 
                    os.system(runCmmd);
                    os.system("mv higgsCombine* "+workspaceDIR+"/"+options.outputDIR);
                    continue;

            #########################
            elif not options.submitOneJobPerToy and options.ntoys != 0: ## make an asimov toy or a set of toys in a single job
                runCmmd = "combine -M Asymptotic -n %s -d %s -s -1 --expectSignal=%f -t %d %s -m %d --rMin %f --rMax %f %s"%(name,dcList[len(dcList)-1],options.expectSignal,options.ntoys,toyOption,int(massPoint),rMin,rMax,asymOpt)

                print "runCmmd ",runCmmd ;
                if options.batchMode:
                    if options.ntoys != -1:
                        fn = "job_Combine_Asymptotic_%s_ntoys_%d.sh"%(name,options.ntoys);
                    else:
                        fn = "job_Combine_Asymptotic_%s_ntoys_Asimov.sh"%(name);
                    submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                else: 
                    os.system(runCmmd);
                    os.system("mv higgsCombine* "+workspaceDIR+"/"+options.outputDIR);
                    continue;

            #########################
            elif options.submitOneJobPerToy and options.ntoys > 0:
                for itoy in range(options.ntoys):
                    runCmmd = "combine -M Asymptotic -n %s_%d -d %s -s -1 --expectSignal=%f -t 1 %s -m %d --rMin %f --rMax %f %s"%(name+"_itoy_",itoy,dcList[len(dcList)-1],options.expectSignal,options.ntoys,toyOption,int(massPoint),rMin,rMax,asymOpt)
                    print "runCmmd ",runCmmd ;
                    if options.batchMode:
                        fn = "job_Combine_Asymptotic_%s_itoy_%d.sh"%(name,itoy);
                        submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                    else: 
                        os.system(runCmmd);
                        os.system("mv higgsCombine* "+workspaceDIR+"/"+options.outputDIR);
                        continue;

            else:
                sys.exit("Check options used for asymptotic limit calculation");
        ########################################
        elif options.ProfileLikelihood and not options.calculateStatOnly: 
            
            if options.ntoys == 0: ## compute the significance on real data --> observed one
                runCmmd = "combine -M ProfileLikelihood --signif -n %s -d %s %s -m %d --rMin %f --rMax %f %s"%(name,dcList[len(dcList)-1],options.pvalue,int(massPoint),rMin,rMax,profOpt);
                print "runCmmd ",runCmmd;               
                if options.batchMode:
                    fn = "job_Combine_ProfileLikelihood_%s.sh"%(name);
                    submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                else:
                    os.system(runCmmd);
                    os.system("mv higgsCombine* "+workspaceDIR+"/"+options.outputDIR);

            elif options.ntoys != 0 and not options.submitOneJobPerToy:    
                runCmmd = "combine -M ProfileLikelihood --signif  -n %s  -d %s -t %d --expectSignal=%f -s -1 %s -m %d --rMin %f --rMax %f %s"%(name,dcList[len(dcList)-1],options.ntoys,options.expectSignal,toyOption,int(massPoint),rMin,rMax,profOpt);
                print "runCmmd ",runCmmd;               
                if options.batchMode:
                    if options.ntoys != -1:
                        fn = "job_Combine_ProfileLikelihood_%s_ntoys_%d.sh"%(name,options.ntoys);
                    else:
                        fn = "job_Combine_ProfileLikelihood_%s_ntoys_Asimov.sh"%(name);

                    submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                else:
                    os.system(runCmmd);
                    os.system("mv higgsCombine* "+workspaceDIR+"/"+options.outputDIR);

            elif options.ntoys > 0 and options.submitOneJobPerToy:
                for itoy in range(options.ntoys):
                    runCmmd = "combine -M ProfileLikelihood --signif  -n %s_%d  -d %s -t 1 --expectSignal=%f -s -1 %s -m %d --rMin %f --rMax %f %s"%(name+"_itoy_",itoy,dcList[len(dcList)-1],options.expectSignal,toyOption,int(massPoint),rMin,rMax,profOpt);
                    print "runCmmd ",runCmmd;               
                    if options.batchMode:
                        fn = "job_Combine_ProfileLikelihood_%s_itoy_%d.sh"%(name,itoy);
                        submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                    else:
                        os.system(runCmmd);
                        os.system("mv higgsCombine* "+workspaceDIR+"/"+options.outputDIR);
                    
            else:
                sys.exit("Check options used for profile likelihood calculation");

        ########################################
        elif options.MaxLikelihoodFit:
                                    
            if options.ntoys == 0 : ### fit real data
                runCmmd = "combine -M MaxLikelihoodFit --rMin %f --rMax %f --saveNormalizations --saveWithUncertainties --saveShapes -n %s -d --saveNLL %s -m %d"%(rMin,rMax,name,dcList[len(dcList)-1],int(massPoint));
                print "runCmmd ",runCmmd;
                if options.batchMode:
                    fn = "job_Combine_MaxLikelihoodFit_%s.sh"%(name);
                    submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                else:   
                    os.system(runCmmd);
                    os.system("mv higgsCombine* "+workspaceDIR+"/"+options.outputDIR);   
                    os.system("mv mlfit* "+workspaceDIR+"/"+options.outputDIR);   
                   
                continue;
            
            ## toys with the same mode
            elif options.ntoys != 0  and not options.crossedToys:
                if not options.submitOneJobPerToy:
                    runCmmd =  "combine -M MaxLikelihoodFit --rMin %f --rMax %f --saveNormalizations --saveWithUncertainties --saveShapes --saveToys -s -1 -n %s -d %s -t %d --expectSignal=%f --saveNLL %s -m %d"%(rMin,rMax,name,dcList[len(dcList)-1],options.nToys,options.expectSignal,toyOption,int(massPoint));
                    print "runCmmd ",runCmmd;
                    if options.batchMode:
                        if options.ntoys != -1:
                            fn = "job_Combine_MaxLikelihoodFit_%s_ntoys_%d.sh"%(name,options.ntoys);
                        else:
                            fn = "job_Combine_MaxLikelihoodFit_%s_ntoys_Asimov.sh"%(name,options.ntoys);

                        submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                    else: 
                        os.system(runCmmd);
                        os.system("mv higgsCombine* "+workspaceDIR+"/"+options.outputDIR);   
                        os.system("mv mlfit* "+workspaceDIR+"/"+options.outputDIR);   
                        
                    continue ;
                elif options.submitOneJobPerToy and options.ntoys > 0: 
                    for itoy in range(options.ntoys):
                        runCmmd =  "combine -M MaxLikelihoodFit  --rMin %f --rMax %f --saveNormalizations --saveWithUncertainties  --saveShapes --saveToys -s -1 -n %s_%d -d %s  -t 1 --expectSignal=%f %s -m %d"%(rMin,rMax,name+"_itoy_",itoy,dcList[len(dcList)-1],options.expectSignal,toyOption,int(massPoint));
                        print "runCmmd ",runCmmd;
                        if options.batchMode:
                            fn = "job_Combine_MaxLikelihoodFit_%s_itoy_%d.sh"%(name,itoy);
                            submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                        else: 
                            os.system(runCmmd);
                            os.system("mv higgsCombine* "+workspaceDIR+"/"+options.outputDIR);   
                            os.system("mv mlfit* "+workspaceDIR+"/"+options.outputDIR);   
                        continue ;
                    
            ###### crossed toys for bias studies
            elif  options.ntoys != 0 and options.crossedToys:
                ## make a list of all generated dataset for the massPoint
                os.system("ls "+options.inputGeneratedDataset+" | grep root | grep higgsCombine | grep "+name+" > list_temp"+suffix);
                itoy = 0;
                list = open("list_temp"+suffix,"r")
                for line in list:
                    itoy += 1;
                    if line == "": continue;
                    ### -t in this case refers to how many toys should be read from the file
                    if not options.submitOneJobPerToy:
                        runCmmd =  "combine -M MaxLikelihoodFit  --rMin %f --rMax %f --saveNormalizations --saveWithUncertainties  --saveShapes -n %s_%d -d %s  -t %d -m %d --toysFile %s/%s"%(rMin,rMax,name,itoy,dcList[len(dcList)-1],options.ntoys,int(massPoint),options.inputGeneratedDataset,line);
                    else:
                        runCmmd =  "combine -M MaxLikelihoodFit  --rMin %f --rMax %f --saveNormalizations --saveWithUncertainties  --saveShapes -n %s_%d -d %s  -t 1 -m %d --toysFile %s/%s"%(rMin,rMax,name,itoy,dcList[len(dcList)-1],options.ntoys,int(massPoint),options.inputGeneratedDataset,line);
                    print "runCmmd ",runCmmd;
                    if options.batchMode:
                        if not options.submitOneJobPerToy:
                            fn = "job_Combine_MaxLikelihoodFit_%s_ntoys_%d.sh"%(name,options.ntoys);
                        else:
                            fn = "job_Combine_MaxLikelihoodFit_%s_itoy_%d.sh"%(name,itoy);
                        submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                    else: 
                        os.system(runCmmd);
                        os.system("mv higgsCombine* "+workspaceDIR+"/"+options.outputDIR);   
                        os.system("mv mlfit* "+workspaceDIR+"/"+options.outputDIR);   
                    continue ;
                continue;

        ########################################
        elif options.makeLikelihoodScan:

            if options.ntoys == 0: ## on rela data     
                runCmmd = "combine -M MultiDimFit -n %s -d %s --algo=grid --points=100 --setPhysicsModelParameterRanges r=%f,%f -m"%(name,dcList[len(dcList)-1],rMin,rMax,int(massPoint));                
                print "runCmmd ",runCmmd;
                if options.batchMode:
                    fn = "job_Combine_LikelihoodScan_%s.sh"%(name);
                    submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                else:
                    os.system(runCmmd);
                    os.system("mv higgsCombine*MultiDimFit* "+workspaceDIR+"/"+options.outputDIR);
            else:
                runCmmd = "combine -M MultiDimFit -n %s -d %s --algo=grid --points=100 --setPhysicsModelParameterRanges r=%f,%f -s -1 --expectSignal=%f -t %d %s"%(name,dcList[len(dcList)-1],rMin,rMax,options.expectSignal,options.ntoys,option);
                print "runCmmd ",runCmmd;
                if options.batchMode:
                    if options.ntoys != -1:
                        fn = "job_Combine_LikelihoodScan_%s_ntoys.sh"%(name,options.ntoys);
                    else:
                        fn = "job_Combine_LikelihoodScan_%s_ntoys_Asimov.sh"%(name);
                    submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                else:
                    os.system(runCmmd);
                    os.system("mv higgsCombine*MultiDimFit* "+workspaceDIR+"/"+options.outputDIR);
            continue;
                
        ########################################
        elif options.MultiDimFit:
            if options.ntoys == 0: ## fit real data
                runCmmd = "combine -M MultiDimFit -n %s -d %s --setPhysicsModelParameterRanges r=%f,%f --saveWorkspace -m %d",(name,dcList[len(dcList)-1],rMin,rMax,int(massPoint));
                print "runCmmd ",runCmmd;
                if options.batchMode:
                    fn = "job_Combine_MultiDimFit_%s.sh"%(name);
                    submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                else:
                    os.system(runCmmd);
                    os.system("mv higgsCombine*MultiDimFit* "+workspaceDIR+"/"+options.outputDIR);
                continue;

            else:                
                runCmmd = "combine -M MultiDimFit -n %s -d %s --setPhysicsModelParameterRanges r=%f,%f --saveWorkspace -m %d -s -1 -t %d --expectSignal=%f",(name,dcList[len(dcList)-1],rMin,rMax,int(massPoint),options.ntoys,options.expectSignal);
                print "runCmmd ",runCmmd;
                if options.batchMode:
                    if options.ntoys != -1:
                        fn = "job_Combine_MultiDimFit_%s_ntoys_%d.sh"%(name,options.ntoys);
                    else:
                        fn = "job_Combine_MultiDimFit_%s_ntoys_Asimov.sh"%(name);
                    submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                else:
                    os.system(runCmmd);
                    os.system("mv higgsCombine*MultiDimFit* "+workspaceDIR+"/"+options.outputDIR);
                continue;

        ########################################
        elif (options.Asymptotic or options.ProfileLikelihood) and options.calculateStatOnly:

             ## make a list of all generated dataset for the massPoint                                                                                                        
            os.system("ls "+options.inputSnapShotDIR+" | grep root | grep higgsCombine | grep "+name+" > list_temp"+suffix);
            list = open("list_temp"+suffix,"r")
            for line in list:
                if options.ntoys == 0: ## fit real data
                    if options.Asymptotic:
                        runCmmd = "combine -M Asymptotic -n %s -d %s --rMin %f --rMax %f --snapshotName MultiDimFit -m %d %s -S 0",(name,line,rMin,rMax,int(massPoint),asymOpt);
                        print "runCmmd ",runCmmd;
                        if options.batchMode:
                            fn = "job_Combine_Asymptotic_%s_statOnly.sh"%(name);
                            submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                        else:
                            os.system(runCmmd);
                            os.system("mv higgsCombine*Asymptotic* "+workspaceDIR+"/"+options.outputDIR);
                        continue;

                    elif options.ProfileLikelihood:
                        runCmmd = "combine -M ProfileLikelihood -n %s -d %s --rMin %f --rMax %f --snapshotName MultiDimFit -m %d %s -S 0",(name,line,rMin,rMax,int(massPoint),profOpt);
                        print "runCmmd ",runCmmd;
                        if options.batchMode:
                            fn = "job_Combine_ProfileLikelihood_%s_statOnly.sh"%(name);
                            submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                        else:
                            os.system(runCmmd);
                            os.system("mv higgsCombine*ProfileLikelihood* "+workspaceDIR+"/"+options.outputDIR);
                        continue;
                
                elif options.ntoys !=0:
                    if options.Asymptotic:
                        runCmmd = "combine -M Asymptotic -n %s -d %s --rMin %f --rMax %f --snapshotName MultiDimFit -m %d %s -S 0 -t %d --expectSignal=%f %s",(name,line,rMin,rMax,int(massPoint),asymOpt,options.ntoys,options.expectSignal,option);
                        print "runCmmd ",runCmmd;
                        if options.batchMode:
                            if options.ntoys != -1:
                                fn = "job_Combine_Asymptotic_%s_ntoys_%d.sh"%(name,options.ntoys);
                            else:
                                fn = "job_Combine_Asymptotic_%s_ntoys_Asimov.sh"%(name);
                            submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                        else:
                            os.system(runCmmd);
                            os.system("mv higgsCombine*Asymptotic* "+workspaceDIR+"/"+options.outputDIR);
                        continue;
                    elif options.ProfileLikelihood:
                        runCmmd = "combine -M ProfileLikelihood -n %s -d %s --rMin %f --rMax %f --snapshotName MultiDimFit -m %d %s -S 0 -t %d --expectSignal=%f %s",(name,line,rMin,rMax,int(massPoint),profOpt,options.ntoys,options.expectSignal,option);
                        print "runCmmd ",runCmmd;
                        if options.batchMode:
                            if options.ntoys != -1:
                                fn = "job_Combine_ProfileLikelihood_%s_ntoys_%d.sh"%(name,options.ntoys);
                            else:
                                fn = "job_Combine_ProfileLikelihood_%s_ntoys_Asimov.sh"%(name);
                            submitBatchJobCombine(runCmmd,fn,name,workspaceDIR+"/"+options.outputDIR);
                        else:
                            os.system(runCmmd);
                            os.system("mv higgsCombine*ProfileLikelihood* "+workspaceDIR+"/"+options.outputDIR);
                        continue;
                
                                                                                                                               
    os.chdir(originalDIR);
    os.system("rm file_temp"+suffix)
    os.system("rm list_temp"+suffix)
