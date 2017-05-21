import ROOT 
import sys,os 

from itertools import product

mmed  =   [10,20,30,40,50,60,70,80,90,100,125,150,200,300,325,400,500,525,600,725,800,925,1000]
mdm   =   [1]
r_range = [2,2,2,2,2,2,2,2,3,3,3,4,4,5,5,6,7,8,10,10,10,10,10]

pos = 0;
for med in mmed:
    for dm in mdm:
        expAV  = '805%04d%04d'%(med,dm) 
        
        cmd = 'combineTool.py -M MultiDimFit datacard_COMB_monoJ.txt --algo=grid --points 40  -m %s --setPhysicsModelParameterRanges r=-%f,%f --job-mode lxbatch --task-name "likelihoodScan_MoriondScalar_MJ_%s" -n _scan_monoJ_%s --sub-opts "-q 1nh " --freezeNuisances lumiScale_MonoJ_MJ,lumiScale_MonoW_MJ,lumiScale_MonoZ_MJ,lumiScale_MonoJ_MV,lumiScale_MonoW_MV,lumiScale_MonoZ_MV --setPhysicsModelParameters lumiScale_MonoJ_MJ=2.78,lumiScale_MonoW_MJ=2.78,lumiScale_MonoZ_MJ=2.78,lumiScale_MonoJ_MV=2.78,lumiScale_MonoW_MV=2.78,lumiScale_MonoZ_MV=2.78 '%(expAV,r_range[pos],r_range[pos],expAV,expAV)
        #print cmd
        os.system(cmd)


        cmd = 'combineTool.py -M MultiDimFit datacard_COMB_monoV.txt --algo=grid --points 40  -m %s --setPhysicsModelParameterRanges r=-%f,%f --job-mode lxbatch --task-name "likelihoodScan_MoriondScalar_MV_%s" -n _scan_monoV_%s --sub-opts "-q 1nh " --freezeNuisances lumiScale_MonoJ_MJ,lumiScale_MonoW_MJ,lumiScale_MonoZ_MJ,lumiScale_MonoJ_MV,lumiScale_MonoW_MV,lumiScale_MonoZ_MV --setPhysicsModelParameters lumiScale_MonoJ_MJ=2.78,lumiScale_MonoW_MJ=2.78,lumiScale_MonoZ_MJ=2.78,lumiScale_MonoJ_MV=2.78,lumiScale_MonoW_MV=2.78,lumiScale_MonoZ_MV=2.78 '%(expAV,r_range[pos],r_range[pos],expAV,expAV)
        #os.system(cmd)

        cmd = 'combineTool.py -M MultiDimFit datacard_COMB.txt --algo=grid --points 40  -m %s --setPhysicsModelParameterRanges r=-%f,%f --job-mode lxbatch --task-name "likelihoodScan_MoriondScalar_COMB_%s" -n _scan_COMB_%s --sub-opts "-q 1nh " --freezeNuisances lumiScale_MonoJ_MJ,lumiScale_MonoW_MJ,lumiScale_MonoZ_MJ,lumiScale_MonoJ_MV,lumiScale_MonoW_MV,lumiScale_MonoZ_MV --setPhysicsModelParameters lumiScale_MonoJ_MJ=2.78,lumiScale_MonoW_MJ=2.78,lumiScale_MonoZ_MJ=2.78,lumiScale_MonoJ_MV=2.78,lumiScale_MonoW_MV=2.78,lumiScale_MonoZ_MV=2.78 '%(expAV,r_range[pos],r_range[pos],expAV,expAV)
        #print cmd
        #os.system(cmd)

    pos = pos+1;
