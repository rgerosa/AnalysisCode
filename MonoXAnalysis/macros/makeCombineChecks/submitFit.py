import ROOT 
import sys,os 

from itertools import product

mmed  =   [10,20,30,40,50,60,70,80,90,100,125,150,200,300,325,400,500,525,600,725,800,925,1000]
mdm   =   [1]
r_range = [1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,2,2,2,3,3,4,4,5,6,8,10,10,10,10,10]

pos = 0;
for med in mmed:
    for dm in mdm:
        expAV  = '805%04d%04d'%(med,dm) 
        
        cmd = 'combineTool.py -M MaxLikelihoodFit datacard_COMB_monoJ.txt --rMin -%f --rMax %f --saveShapes --saveNormalizations -m %s --job-mode lxbatch --task-name "fit_MoriondScalar_MJ_%s" -n _fit_monoJ_%s --sub-opts "-q 1nh " --freezeNuisances lumiScale_MonoJ_MJ,lumiScale_MonoW_MJ,lumiScale_MonoZ_MJ,lumiScale_MonoJ_MV,lumiScale_MonoW_MV,lumiScale_MonoZ_MV --setPhysicsModelParameters lumiScale_MonoJ_MJ=2.78,lumiScale_MonoW_MJ=2.78,lumiScale_MonoZ_MJ=2.78,lumiScale_MonoJ_MV=2.78,lumiScale_MonoW_MV=2.78,lumiScale_MonoZ_MV=2.78 '%(r_range[pos],r_range[pos],expAV,expAV,expAV)
        #print cmd
        #os.system(cmd)

        cmd = 'combineTool.py -M MaxLikelihoodFit datacard_COMB_monoV.txt --rMin -%f --rMax %f --saveShapes --saveNormalizations -m %s --job-mode lxbatch --task-name "fit_MoriondScalar_MV_%s" -n _fit_monoV_%s --sub-opts "-q 1nh " --freezeNuisances lumiScale_MonoJ_MJ,lumiScale_MonoW_MJ,lumiScale_MonoZ_MJ,lumiScale_MonoJ_MV,lumiScale_MonoW_MV,lumiScale_MonoZ_MV --setPhysicsModelParameters lumiScale_MonoJ_MJ=2.78,lumiScale_MonoW_MJ=2.78,lumiScale_MonoZ_MJ=2.78,lumiScale_MonoJ_MV=2.78,lumiScale_MonoW_MV=2.78,lumiScale_MonoZ_MV=2.78 '%(r_range[pos],r_range[pos],expAV,expAV,expAV)
        #os.system(cmd)

        cmd = 'combineTool.py -M MaxLikelihoodFit datacard_COMB.txt --rMin -%f --rMax %f --saveShapes --saveNormalizations -m %s --job-mode lxbatch --task-name "fit_MoriondScalar_COMB_%s" -n _fit_COMB_%s --sub-opts "-q 1nh " --freezeNuisances lumiScale_MonoJ_MJ,lumiScale_MonoW_MJ,lumiScale_MonoZ_MJ,lumiScale_MonoJ_MV,lumiScale_MonoW_MV,lumiScale_MonoZ_MV --setPhysicsModelParameters lumiScale_MonoJ_MJ=2.78,lumiScale_MonoW_MJ=2.78,lumiScale_MonoZ_MJ=2.78,lumiScale_MonoJ_MV=2.78,lumiScale_MonoW_MV=2.78,lumiScale_MonoZ_MV=2.78 '%(r_range[pos],r_range[pos],expAV,expAV,expAV)
        os.system(cmd)

    pos = pos+1;
