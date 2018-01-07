import os
import glob
import math
from array import array
import sys
import time
import subprocess
import ROOT
import random

from optparse import OptionParser
from subprocess import Popen

############################################                                                                                                                                                           
#            Job steering                  #                                                                                                                                                           
############################################                                                                                                                                                           

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

## parse files                                                                                                                                                                                         
parser.add_option('--inputTemplateDIR', action="store", type="string", dest="inputTemplateDIR",      default="",   help="Input directory where the template file are stored")
parser.add_option('--outputBaseDIR',    action="store", type="string", dest="outputBaseDIR", default="",   help="Output base directory where madgraph card will be placed")
parser.add_option('--processType',      action="store", type="string", dest="processType",   default="",   help="Type of process to be generated: MonoJ, MonoW or MonZ")

(options, args) = parser.parse_args()

if __name__ == '__main__':


    if options.processType != "MonoJ" and options.processType != "MonoW" and options.processType != "MonoZ":
        sys.exit("Unknown process type --> please fix");

    os.system("mkdir -p "+options.outputBaseDIR);
    

    mass_scan = { 'mh4' : [100,200,300,400,500,600,700,800,900],
                  'mh3' : [100,200,300,400,500,600,700,800,900,1000,1100,1200,1300],
                  'tanbeta': 1,
                  'sinp': 0.35,
                  'lam3' : 3,
                  'laP1' : 3,
                  'laP2' : 3,
                  'gPXd' : 1,
                  'sinbma' : 1,
                  'Mxd' : 10,
                  'mh1' : 125
                  };

    dm_scan   = { 'mh4' : [250],
                  'mh3' : [600],
                  'mh1' : 125,
                  'tanbeta': 1,
                  'sinp': 0.35,
                  'lam3' : 3,
                  'laP1' : 3,
                  'laP2' : 3,
                  'gPXd' : 1,
                  'sinbma' : 1,
                  'Mxd' : [1,50,100,150,200,250,300,350,400,450,500]
                  };

    tanb_scan   = { 'mh4' : [100,150,200,250,300,350,400,450,500],
                    'mh3' : [600],
                    'mh1' : 125,
                    'tanbeta': [0.5,1,1.5,2,2.5,3,3.5,4,5,6,7,8,9,10],
                    'sinp': 0.35,
                    'lam3' : 3,
                    'laP1' : 3,
                    'laP2' : 3,
                    'gPXd' : 1,
                    'sinbma' : 1,
                    'Mxd' : 10,
                    };


    if len(dm_scan['mh4']) != len(dm_scan['mh3']) :
        sys.exit("Problem with the DM mass scan lenght --> exit");


    #### make the output task name
    for mh3 in mass_scan['mh3']:
        for mh4 in mass_scan['mh4']:
        
            name = options.processType+"_mh3_"+str(mh3)+"_mh4_"+str(mh4);
            os.system("mkdir -p "+options.outputBaseDIR+"/"+name);
            os.system("cp "+options.inputTemplateDIR+"/customizecards.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            
            ### make customized card
            os.system("sed -i -e 's/@gPXd/"+str(mass_scan['gPXd'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@Mxd/"+str(mass_scan['Mxd'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@sinbma/"+str(mass_scan['sinbma'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@laP1/"+str(mass_scan['laP1'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@laP2/"+str(mass_scan['laP2'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@lam3/"+str(mass_scan['lam3'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@tanbeta/"+str(mass_scan['tanbeta'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@sinp/"+str(mass_scan['sinp'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@mh1/"+str(mass_scan['mh1'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@mh2/"+str(mh3)+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@mh3/"+str(mh3)+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@mhc/"+str(mh3)+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@mh4/"+str(mh4)+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");

            ### make extramodel
            os.system("cp "+options.inputTemplateDIR+"/extramodels.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_extramodels.dat");
            os.system("cp "+options.inputTemplateDIR+"/*tgz "+options.outputBaseDIR+"/"+name+"/");
            os.system("cp "+options.inputTemplateDIR+"/*zip "+options.outputBaseDIR+"/"+name+"/");
            ### make run card
            os.system("cp "+options.inputTemplateDIR+"/run_card.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_run_card.dat");
            ### make madsping card
            os.system("cp "+options.inputTemplateDIR+"/madspin_card.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_madspin_card.dat");
            ### make proc card
            os.system("cp "+options.inputTemplateDIR+"/proc_card.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
            
            if options.processType == "MonoZ":
                os.system("sed -i -e 's/@MODEL/Pseudoscalar_2HDM-bbMET_5FS/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
                process = "define p =  g u c d s u~ c~ d~ s~ b b~\\n";
                process += "define j = g u c d s b u~ c~ d~ s~ b~ \\n";
                process += "generate g g  > xd xd~ z [QCD] \\n";
                os.system("sed -i -e 's/@PROCESS/"+process+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
                os.system("sed -i -e 's/@OUTPUT/"+name+" -nojpeg/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");

            elif options.processType == "MonoW":
                os.system("sed -i -e 's/@MODEL/Pseudoscalar_2HDM-WMET/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
                process = "define p = g u c d s u~ c~ d~ s~ b b~ \\n";
                process += "define j =  g u c d s b u~ c~ d~ s~ b~ \\n";
                process += "generate p p  > xd xd~ w+ [QCD]  \\n";
                process += "add process p p  > xd xd~ w- [QCD]  \\n";
                os.system("sed -i -e 's/@PROCESS/"+process+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
                os.system("sed -i -e 's/@OUTPUT/"+name+" -nojpeg/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");

            elif options.processType == "MonoJ":
                os.system("sed -i -e 's/@MODEL/Pseudoscalar_2HDM-bbMET_5FS/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
                process = "define p = g u c d s u~ c~ d~ s~ b b~ \\n";
                process += "define j =  g u c d s b u~ c~ d~ s~ b~ \\n";
                process += "generate p p  > xd xd~ j  [QCD] \\n";
                os.system("sed -i -e 's/@PROCESS/"+process+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
                os.system("sed -i -e 's/@OUTPUT/"+name+" -nojpeg/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");


    #### make the output task name
    for mxd in dm_scan['Mxd']:
        
        name = options.processType+"_mxd_"+str(mxd);
        os.system("mkdir -p "+options.outputBaseDIR+"/"+name);
        os.system("cp "+options.inputTemplateDIR+"/customizecards.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            
        ### make customized card
        os.system("sed -i -e 's/@gPXd/"+str(dm_scan['gPXd'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        os.system("sed -i -e 's/@Mxd/"+str(mxd)+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        os.system("sed -i -e 's/@sinbma/"+str(dm_scan['sinbma'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        os.system("sed -i -e 's/@laP1/"+str(dm_scan['laP1'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        os.system("sed -i -e 's/@laP2/"+str(dm_scan['laP2'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        os.system("sed -i -e 's/@lam3/"+str(dm_scan['lam3'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        os.system("sed -i -e 's/@tanbeta/"+str(dm_scan['tanbeta'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        os.system("sed -i -e 's/@sinp/"+str(dm_scan['sinp'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        os.system("sed -i -e 's/@mh1/"+str(dm_scan['mh1'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        os.system("sed -i -e 's/@mh2/"+str(dm_scan['mh3'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        os.system("sed -i -e 's/@mh3/"+str(dm_scan['mh3'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        os.system("sed -i -e 's/@mhc/"+str(dm_scan['mh3'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        os.system("sed -i -e 's/@mh4/"+str(dm_scan['mh4'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
        
        ### make extramodel
        os.system("cp "+options.inputTemplateDIR+"/extramodels.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_extramodels.dat");
        os.system("cp "+options.inputTemplateDIR+"/*tgz "+options.outputBaseDIR+"/"+name+"/");
        os.system("cp "+options.inputTemplateDIR+"/*zip "+options.outputBaseDIR+"/"+name+"/");

        ### make run card
        os.system("cp "+options.inputTemplateDIR+"/run_card.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_run_card.dat");

        ### make madsping card
        os.system("cp "+options.inputTemplateDIR+"/madspin_card.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_madspin_card.dat");
            
        ### make proc card
        os.system("cp "+options.inputTemplateDIR+"/proc_card.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
        
        if options.processType == "MonoZ":
            os.system("sed -i -e 's/@MODEL/Pseudoscalar_2HDM-bbMET_5FS/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
            process = "define p = g u c d s u~ c~ d~ s~ b b~ \\n";
            process += "define j =  g u c d s b u~ c~ d~ s~ b~ \\n";
            process += "generate g g  > xd xd~ z  [QCD] \\n";
            os.system("sed -i -e 's/@PROCESS/"+process+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
            os.system("sed -i -e 's/@OUTPUT/"+name+" -nojpeg/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");

        elif options.processType == "MonoW":
            os.system("sed -i -e 's/@MODEL/Pseudoscalar_2HDM-WMET/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
            process = "define p =  g u c d s u~ c~ d~ s~ b b~ \\n";
            process += "define j =  g u c d s b u~ c~ d~ s~ b~ \\n";
            process += "generate p p  > xd xd~ w+  [QCD] \\n";
            process += "add process p p  > xd xd~ w-  [QCD] \\n";
            os.system("sed -i -e 's/@PROCESS/"+process+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
            os.system("sed -i -e 's/@OUTPUT/"+name+" -nojpeg/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");

        elif options.processType == "MonoJ":
            os.system("sed -i -e 's/@MODEL/Pseudoscalar_2HDM-bbMET_5FS/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
            process = "define p = g u c d s u~ c~ d~ s~ b b~ \\n";
            process += "define j = g u c d s b u~ c~ d~ s~ b~ \\n";
            process += "generate p p  > xd xd~ j  [QCD] \\n";
            os.system("sed -i -e 's/@PROCESS/"+process+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
            os.system("sed -i -e 's/@OUTPUT/"+name+" -nojpeg/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");


    #### make the output task name
    for mh4 in tanb_scan['mh4']:
        for tanb in tanb_scan['tanbeta']:
        
            name = options.processType+"_mh4_"+str(mh4)+"_tanb_"+str(tanb);
            os.system("mkdir -p "+options.outputBaseDIR+"/"+name);
            os.system("cp "+options.inputTemplateDIR+"/customizecards.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            
            ### make customized card
            os.system("sed -i -e 's/@gPXd/"+str(tanb_scan['gPXd'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@Mxd/"+str(tanb_scan['Mxd'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@sinbma/"+str(tanb_scan['sinbma'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@laP1/"+str(tanb_scan['laP1'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@laP2/"+str(tanb_scan['laP2'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@lam3/"+str(tanb_scan['lam3'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@tanbeta/"+str(tanb)+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@sinp/"+str(tanb_scan['sinp'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@mh1/"+str(tanb_scan['mh1'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@mh2/"+str(tanb_scan['mh3'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@mh3/"+str(tanb_scan['mh3'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@mhc/"+str(tanb_scan['mh3'])+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");
            os.system("sed -i -e 's/@mh4/"+str(mh4)+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_customizecards.dat");

            ### make extramodel
            os.system("cp "+options.inputTemplateDIR+"/extramodels.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_extramodels.dat");
            os.system("cp "+options.inputTemplateDIR+"/*tgz "+options.outputBaseDIR+"/"+name+"/");
            os.system("cp "+options.inputTemplateDIR+"/*zip "+options.outputBaseDIR+"/"+name+"/");
            ### make run card
            os.system("cp "+options.inputTemplateDIR+"/run_card.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_run_card.dat");
            ### make madsping card
            os.system("cp "+options.inputTemplateDIR+"/madspin_card.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_madspin_card.dat");
            ### make proc card
            os.system("cp "+options.inputTemplateDIR+"/proc_card.dat "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
            
            if options.processType == "MonoZ":
                os.system("sed -i -e 's/@MODEL/Pseudoscalar_2HDM-bbMET_5FS/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
                process = "define p =  g u c d s u~ c~ d~ s~ b b~\\n";
                process += "define j = g u c d s b u~ c~ d~ s~ b~ \\n";
                process += "generate g g  > xd xd~ z [QCD] \\n";
                os.system("sed -i -e 's/@PROCESS/"+process+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
                os.system("sed -i -e 's/@OUTPUT/"+name+" -nojpeg/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");

            elif options.processType == "MonoW":
                os.system("sed -i -e 's/@MODEL/Pseudoscalar_2HDM-WMET/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
                process = "define p = g u c d s u~ c~ d~ s~ b b~ \\n";
                process += "define j =  g u c d s b u~ c~ d~ s~ b~ \\n";
                process += "generate p p  > xd xd~ w+ [QCD]  \\n";
                process += "add process p p  > xd xd~ w- [QCD]  \\n";
                os.system("sed -i -e 's/@PROCESS/"+process+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
                os.system("sed -i -e 's/@OUTPUT/"+name+" -nojpeg/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");

            elif options.processType == "MonoJ":
                os.system("sed -i -e 's/@MODEL/Pseudoscalar_2HDM-bbMET_5FS/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
                process = "define p = g u c d s u~ c~ d~ s~ b b~ \\n";
                process += "define j =  g u c d s b u~ c~ d~ s~ b~ \\n";
                process += "generate p p  > xd xd~ j  [QCD] \\n";
                os.system("sed -i -e 's/@PROCESS/"+process+"/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
                os.system("sed -i -e 's/@OUTPUT/"+name+" -nojpeg/g' "+options.outputBaseDIR+"/"+name+"/"+name+"_proc_card.dat");
