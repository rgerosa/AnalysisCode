--batchMode#python makeWorkspace.py --inputDIR output_exclusive_23bins --observable met --category 1 --templateFile templates_met.root 
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
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",      default="",             help="input directory where files are located")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",     default="workspace",    help="input directory where files are located")
parser.add_option('--templateFile', action="store", type="string", dest="templateFile",  default="",    help="input directory where files are located")
parser.add_option('--category',     action="store", type=int     , dest="category",      default=0,     help="input directory where files are located")
parser.add_option('--observable',   action="store", type="string", dest="observable",    default="met", help="input directory where files are located")
parser.add_option('--shapeSys',     action="store", type=int     , dest="shapeSys",      default=0,     help="input directory where files are located")
parser.add_option('--connectWZ',    action="store", type=int     , dest="connectWZ",     default=1,     help="input directory where files are located")
parser.add_option('--connectTop',   action="store", type=int     , dest="connectTop",    default=0,     help="input directory where files are located")
parser.add_option('--scaleQCD',     action="store", type=float   , dest="scaleQCD",      default=2,     help="input directory where files are located")

### batch submission
parser.add_option('--batchMode',    action="store_true", dest="batchMode",   help="batchMode")
parser.add_option('--submit',       action="store_true", dest="submit", help="submit")
parser.add_option('--jobDIR',       action="store",      type="string", dest="jobDIR", default="",  help="directory for job")
parser.add_option('--queque',       action="store",      type="string", dest="queque", default="",  help="queque for LSF")


(options, args) = parser.parse_args()

if __name__ == '__main__':

    ## load the create workspace macro
    ROOT.gROOT.ProcessLine(".L ./macros/createWorkspace.C");

    startingDIR = os.getcwd()

    ## set the working dir in the directory where the template file is located
    if not options.batchMode:
        os.chdir(options.inputDIR)
    print os.getcwd()

    currentDIR = os.getcwd();

    ## open the template file and make the list of signal points
    if currentDIR == startingDIR+"/"+options.inputDIR:
        inputTemplate = ROOT.TFile(options.templateFile,"READ")
    else:
        inputTemplate = ROOT.TFile(options.inputDIR+"/"+options.templateFile,"READ")
    inputTemplate.cd();

    monoJInteraction = [];
    monoJMediatorMass = [];
    monoJdmMass = [];

    monoWInteraction = [];
    monoWMediatorMass = [];
    monoWdmMass = [];

    monoZInteraction = [];
    monoZMediatorMass = [];
    monoZdmMass = [];

    for key in inputTemplate.GetListOfKeys():
        if(key.GetClassName() == "TH1" or key.GetClassName() == "TH1F" or key.GetClassName() == "TH1D"):
            h1 = key.ReadObj()
            if not ROOT.TString(h1.GetName()).Contains("monoJ") and not ROOT.TString(h1.GetName()).Contains("monoW") and not ROOT.TString(h1.GetName()).Contains("monoZ"):
                continue;
            if ROOT.TString(h1.GetName()).Contains("Up") or ROOT.TString(h1.GetName()).Contains("Down") or ROOT.TString(h1.GetName()).Contains("Dw"):
                continue;
            list = str(h1.GetName()).split("_")
            if(len(list) < 4):
                sys.exit("Problem in the signal template name convention ---> exit ");
            
            if(ROOT.TString(h1.GetName()).Contains("monoJ")):                
                monoJInteraction.append(list [1]);
                monoJMediatorMass.append(list[2]);
                monoJdmMass.append(list[3]);
            elif (ROOT.TString(h1.GetName()).Contains("monoW")):
                monoWInteraction.append(list [1]);
                monoWMediatorMass.append(list[2]);
                monoWdmMass.append(list[3]);
            elif (ROOT.TString(h1.GetName()).Contains("monoZ")):
                monoZInteraction.append(list [1]);
                monoZMediatorMass.append(list[2]);
                monoZdmMass.append(list[3]);

    
    if ((len(monoJInteraction) != len(monoWInteraction)) or  (len(monoWInteraction) != len(monoZInteraction)) or (len(monoJInteraction) != len(monoZInteraction))):
        sys.exit("Problem with the mass point size between monoJ, monoW and monoZ  ---> exit ");

        
    for isig in range(len(monoJInteraction)):
        
        if options.category <= 1:
            cat = "MJ"
        else:
            cat = "MV"

        if monoJInteraction[isig] != monoWInteraction[isig] or monoJInteraction[isig] != monoZInteraction[isig]:
             sys.exit("Problem with the mass point interaction among monoJ, monoW and monoZ  ---> exit ");
        if monoJMediatorMass[isig] != monoWMediatorMass[isig] or monoJMediatorMass[isig] != monoZMediatorMass[isig]:
             sys.exit("Problem with the mass point mediator among monoJ, monoW and monoZ  ---> exit ");
        if monoJdmMass[isig] != monoWdmMass[isig] or monoJdmMass[isig] != monoZdmMass[isig]:
             sys.exit("Problem with the mass point DM among monoJ, monoW and monoZ  ---> exit ");

        if not options.batchMode:


            os.system("mkdir -p "+options.outputDIR);
    
            command = ROOT.TString("createWorkspace(\"%s\",%d,\"%s/workspace_%s_%s_%s_%s.root\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d,%d,%d)"%(options.templateFile,options.category,options.outputDIR,cat,monoJInteraction[isig],monoJMediatorMass[isig],monoJdmMass[isig],options.observable,monoJInteraction[isig],monoJMediatorMass[isig],monoJdmMass[isig],options.scaleQCD,options.connectWZ,options.connectTop,options.shapeSys))
            ROOT.gROOT.ProcessLine(command.Data());
        else:
            
            os.system("mkdir -p "+options.jobDIR);
            os.system("mkdir -p "+currentDIR+"/"+options.inputDIR+"/"+options.outputDIR);

            jobName = "job_%s_%s_%s_%s"%(cat,monoJInteraction[isig],monoJMediatorMass[isig],monoJdmMass[isig]);
            
            ## write job sh file                                                                                                                                             
            jobmacro = open('%s/%s.C'%(options.jobDIR,jobName),'w')
            jobmacro.write("{\n");
            jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/createWorkspace.C\");\n");
            command = ROOT.TString("\"createWorkspace(\\\"%s\\\",%d,\\\"workspace_%s_%s_%s_%s.root\\\",\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",%f,%d,%d,%d)\""%(options.templateFile,options.category,cat,monoJInteraction[isig],monoJMediatorMass[isig],monoJdmMass[isig],options.observable,monoJInteraction[isig],monoJMediatorMass[isig],monoJdmMass[isig],options.scaleQCD,options.connectWZ,options.connectTop,options.shapeSys))
            jobmacro.write("gROOT->ProcessLine("+command.Data()+");\n");
            jobmacro.write("}\n");
            jobmacro.close();


            jobscript = open('%s/%s.sh'%(options.jobDIR,jobName),'w')
            jobscript.write('cd %s \n'%currentDIR)
            jobscript.write('eval ` scramv1 runtime -sh ` \n')
            jobscript.write('cd - \n')
            jobscript.write('scp %s/%s ./\n'%(currentDIR+"/"+options.inputDIR,options.templateFile))
            jobscript.write('scp '+currentDIR+'/%s/%s.C ./ \n'%(options.jobDIR,jobName))
            jobscript.write('root -l -b -q %s.C \n'%(jobName))
            jobscript.write('scp workspace_%s_%s_%s_%s.root %s/%s/%s\n'%(cat,monoJInteraction[isig],monoJMediatorMass[isig],monoJdmMass[isig],currentDIR,options.inputDIR,options.outputDIR))
            
            os.system('chmod a+x %s/%s.sh'%(options.jobDIR,jobName))
            
            if options.submit:
                os.system('bsub -q %s -o %s/%s/%s.log -e %s/%s/%s.err %s/%s/%s.sh'%(options.queque,currentDIR,options.jobDIR,jobName,currentDIR,options.jobDIR,jobName,currentDIR,options.jobDIR,jobName))
        
