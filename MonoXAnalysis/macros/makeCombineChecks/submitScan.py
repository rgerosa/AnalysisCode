import ROOT 
import sys,os 

from itertools import product

mmed  = [10,20,30,40,50,60,70,80,90,100,125,150,200,300,325,400,500,525,600,725,800,925,1000]
mdm   = [1]

expAV  = [ '805%04d%04d'%(i,j) for i,j in product(mmed,mdm) ] 

# make a crazy string
# This string has the number for coupling, mediator mass, and darkmattermass 

striAV = ",".join(expAV)

cmd = 'combineTool.py -M Asymptotic datacard_COMB_monoJ.txt --cl=0.95 -m %s --job-mode lxbatch --task-name "MoriondScalar_MJ" -n _monoJ --sub-opts "-q 1nh " --freezeNuisances lumiScale_MonoJ_MJ,lumiScale_MonoW_MJ,lumiScale_MonoZ_MJ,lumiScale_MonoJ_MV,lumiScale_MonoW_MV,lumiScale_MonoZ_MV --setPhysicsModelParameters lumiScale_MonoJ_MJ=2.78,lumiScale_MonoW_MJ=2.78,lumiScale_MonoZ_MJ=2.78,lumiScale_MonoJ_MV=2.78,lumiScale_MonoW_MV=2.78,lumiScale_MonoZ_MV=2.78 '%(striAV)

#print cmd 
#os.system(cmd)

cmd = 'combineTool.py -M Asymptotic datacard_COMB_monoV.txt --cl=0.95 -m %s --job-mode lxbatch --task-name "MoriondScalar_MV" -n _monoV --sub-opts "-q 1nh " --freezeNuisances lumiScale_MonoJ_MJ,lumiScale_MonoW_MJ,lumiScale_MonoZ_MJ,lumiScale_MonoJ_MV,lumiScale_MonoW_MV,lumiScale_MonoZ_MV --setPhysicsModelParameters lumiScale_MonoJ_MJ=2.78,lumiScale_MonoW_MJ=2.78,lumiScale_MonoZ_MJ=2.78,lumiScale_MonoJ_MV=2.78,lumiScale_MonoW_MV=2.78,lumiScale_MonoZ_MV=2.78 '%(striAV)

#print cmd 
#os.system(cmd)

cmd = 'combineTool.py -M Asymptotic datacard_COMB.txt --cl=0.95 -m %s --job-mode lxbatch --task-name "MoriondScalar" -n _COMB --sub-opts "-q 1nh " --freezeNuisances lumiScale_MonoJ_MJ,lumiScale_MonoW_MJ,lumiScale_MonoZ_MJ,lumiScale_MonoJ_MV,lumiScale_MonoW_MV,lumiScale_MonoZ_MV --setPhysicsModelParameters lumiScale_MonoJ_MJ=2.78,lumiScale_MonoW_MJ=2.78,lumiScale_MonoZ_MJ=2.78,lumiScale_MonoJ_MV=2.78,lumiScale_MonoW_MV=2.78,lumiScale_MonoZ_MV=2.78 '%(striAV)

#print cmd 
os.system(cmd)
