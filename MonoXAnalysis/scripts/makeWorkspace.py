#python makeWorkspace.py --inputDIR output_exclusive_23bins --observable met --category 1 --templateFile templates_met.root 
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

############################################                                                                                                                                 
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('--inputDIR',     action="store", type="string", dest="inputDIR",      default="",             help="to be used to correctly set the working area")
parser.add_option('--outputDIR',    action="store", type="string", dest="outputDIR",     default="workspace",    help="output dir where workspaces are copied")
parser.add_option('--templateFile', action="store", type="string", dest="templateFile",  default="",             help="file with all tempaltes: signal, bkg and data")
parser.add_option('--category',     action="store", type=int     , dest="category",      default=0,              help="analysis category to pick up right binning")
parser.add_option('--observable',   action="store", type="string", dest="observable",    default="met",          help="observable name")

############################################                                                                                                                                 
parser.add_option('--addShapeSystematics',     action="store_true",  dest="addShapeSystematics", default=False,  help="add shape sys for bin-by-bin and experimental sources")
parser.add_option('--connectWZ',    action="store_true",  dest="connectWZ",   default=False,  help="connect W+jets and Z->nunu SR")
parser.add_option('--connectTop',   action="store_true",  dest="connectTop",  default=False,  help="extract ttbar+stop from control sample")
parser.add_option('--isHiggsInvisible', action="store_true",  dest="isHiggsInvisible",  default=False,  help="create workspace for Higgs Invisible different mH")
parser.add_option('--scaleQCD',     action="store", type=int, dest="scaleQCD",          default=1,      help="scale factor for QCD backgrond in SR")
parser.add_option('--mergeLeptons', action="store_true",  dest="mergeLeptons",          default=False,  help="merge muon and electron datasets")
parser.add_option('--isCombination',action="store_true",  dest="isCombination",         default=False,  help="naming convention for Higgs invisible combination")
parser.add_option('--interaction',  action="store", type="string", dest="interaction",  default="",     help="Inteaction type")

############################################                                                                                                                                 
parser.add_option('--batchMode',    action="store_true", dest="batchMode",   help="batchMode")
parser.add_option('--submit',       action="store_true", dest="submit", help="submit")
parser.add_option('--jobDIR',       action="store",      type="string", dest="jobDIR", default="",  help="directory for job")
parser.add_option('--queque',       action="store",      type="string", dest="queque", default="",  help="queque for LSF")


(options, args) = parser.parse_args()

############################################                                                                                                                                 
if __name__ == '__main__':


    ## transform some inputs
    connectWZ = 0;
    if options.connectWZ:
        connectWZ = 1;

    connectTop = 0;
    if options.connectTop:
        connectTop = 1;

    isHiggsInvisible = 0;
    if options.isHiggsInvisible:
        isHiggsInvisible = 1;

    addShapeSystematics = 0;
    if options.addShapeSystematics:
        addShapeSystematics = 1;
    
    mergeLeptons = 0;
    if options.mergeLeptons:
        mergeLeptons = 1;

    isCombination = 0;
    if options.isCombination:
        isCombination = 1;

    ## load the create workspace macro
    ROOT.gROOT.ProcessLine(".L $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so");
    ROOT.gROOT.ProcessLine(".L ./macros/createWorkspace.C");

    ## set starting and working dir
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

    ggHMass  = [];
    vbfHMass = [];
    WHMass   = [];
    ZHMass   = [];
    ggZHMass = [];

    ## make signal mass point list in case of DM analysis
    ## loop on template file list of keys
    for key in inputTemplate.GetListOfKeys():
        inputTemplate.cd();
        ## enter inside directory
        if(key.GetClassName() == "TDirectory" or key.GetClassName() == "TDirectoryFile"):
            dir1 = key.ReadObj();
            dir1.cd();
            for key2 in dir1.GetListOfKeys():                
                if(key2.GetClassName() == "TH1" or key2.GetClassName() == "TH1F" or key2.GetClassName() == "TH1D"):
                    h1 = key2.ReadObj()
                    if not options.isHiggsInvisible:
                        if not ROOT.TString(h1.GetName()).Contains("monoJ") and not ROOT.TString(h1.GetName()).Contains("monoW") and not ROOT.TString(h1.GetName()).Contains("monoZ"):
                            continue;
                        if ROOT.TString(h1.GetName()).Contains("Up") or ROOT.TString(h1.GetName()).Contains("Down") or ROOT.TString(h1.GetName()).Contains("Dw"):
                            continue;
                        if not ROOT.TString(h1.GetName()).Contains(options.interaction):
                            continue;
                        list = str(h1.GetName()).split("_")
                        if(len(list) < 4):
                            sys.exit("Problem in the signal template name convention ---> exit ");
                        else:                            
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
                    else:

                        if not ROOT.TString(h1.GetName()).Contains("ggH") and not ROOT.TString(h1.GetName()).Contains("vbfH") and not ROOT.TString(h1.GetName()).Contains("wH") and not ROOT.TString(h1.GetName()).Contains("zH") and not ROOT.TString(h1.GetName()).Contains("ggZH"):
                            continue;
                        if ROOT.TString(h1.GetName()).Contains("Up") or ROOT.TString(h1.GetName()).Contains("Down") or ROOT.TString(h1.GetName()).Contains("Dw"):
                            continue;
                        list = str(h1.GetName()).split("_")
                        if(len(list) < 3):
                            sys.exit("Problem in the signal template name convention ---> exit ");
                        else:
                            if(ROOT.TString(h1.GetName()).Contains("ggH")):
                                ggHMass.append(list[2]);
                            elif (ROOT.TString(h1.GetName()).Contains("vbfH")):
                                vbfHMass.append(list[2]);
                            elif (ROOT.TString(h1.GetName()).Contains("wH")):
                                WHMass.append(list[2]);
                            elif (ROOT.TString(h1.GetName()).Contains("zH")):
                                ZHMass.append(list[2]);
                            elif (ROOT.TString(h1.GetName()).Contains("ggZH")):
                                ggZHMass.append(list[2]);


            if not isHiggsInvisible:                    
                if ((len(monoJInteraction) != len(monoWInteraction)) or  (len(monoWInteraction) != len(monoZInteraction)) or (len(monoJInteraction) != len(monoZInteraction))):
                    sys.exit("Problem with the mass point size between monoJ, monoW and monoZ  ---> exit ");
            else:
                if ((len(ggHMass) != len(vbfHMass)) or  (len(ggHMass) != len(WHMass)) or (len(ggHMass) != len(ZHMass)) or (len(ggHMass) != len(ggZHMass)) ):
                    sys.exit("Problem with the mass point size between ggH, qqH, wH, zH and ggZH  ---> exit ");
            
        else:                    
            print"Problem .. in the base file found a key which is not a TDirectory : ",key.GetClassName()," ",key.GetName();
                
            
    if options.category <= 1:
        cat = "MJ"
    else:
        cat = "MV"
                
    #### create jobs for DM
    if not isHiggsInvisible:
        for isig in range(len(monoJInteraction)):            

            if monoJInteraction[isig] != monoWInteraction[isig] or monoJInteraction[isig] != monoZInteraction[isig]:
                sys.exit("Problem with the mass point interaction among monoJ, monoW and monoZ  ---> exit ");
            if monoJMediatorMass[isig] != monoWMediatorMass[isig] or monoJMediatorMass[isig] != monoZMediatorMass[isig]:
                sys.exit("Problem with the mass point mediator among monoJ, monoW and monoZ  ---> exit ");
            if monoJdmMass[isig] != monoWdmMass[isig] or monoJdmMass[isig] != monoZdmMass[isig]:
                sys.exit("Problem with the mass point DM among monoJ, monoW and monoZ  ---> exit ");

            if not options.batchMode:

                os.system("mkdir -p "+options.outputDIR);
    
                command = ROOT.TString("createWorkspace(\"%s\",%d,\"%s/workspace_%s_%s_%s_%s.root\",\"%s\",%d,%f,%d,%d,%d,%d,%d,\"%s\",\"%s\",\"%s\")"%(options.templateFile,options.category,options.outputDIR,cat,monoJInteraction[isig],monoJMediatorMass[isig],monoJdmMass[isig],options.observable,isHiggsInvisible,options.scaleQCD,connectWZ,connectTop,addShapeSystematics,mergeLeptons,isCombination,monoJInteraction[isig],monoJMediatorMass[isig],monoJdmMass[isig]))
                ROOT.gROOT.ProcessLine(command.Data());
            else:
            
                os.system("mkdir -p "+options.jobDIR);
                os.system("mkdir -p "+currentDIR+"/"+options.inputDIR+"/"+options.outputDIR);

                jobName = "job_%s_%s_%s_%s"%(cat,monoJInteraction[isig],monoJMediatorMass[isig],monoJdmMass[isig]);
            
                ## write job sh file                                                                                                                                      
                jobmacro = open('%s/%s.C'%(options.jobDIR,jobName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so\");\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/createWorkspace.C\");\n");
                command = ROOT.TString("\"createWorkspace(\\\"%s\\\",%d,\\\"workspace_%s_%s_%s_%s.root\\\",\\\"%s\\\",%d,%f,%d,%d,%d,%d,%d,\\\"%s\\\",\\\"%s\\\",\\\"%s\\\")\""%(options.templateFile,options.category,cat,monoJInteraction[isig],monoJMediatorMass[isig],monoJdmMass[isig],options.observable,isHiggsInvisible,options.scaleQCD,connectWZ,connectTop,addShapeSystematics,mergeLeptons,isCombination,monoJInteraction[isig],monoJMediatorMass[isig],monoJdmMass[isig]))
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


    ## for Higgs invisible                    
    else: 

        for mass in ggHMass:

            if not options.batchMode:

                command = ROOT.TString("createWorkspace(\"%s\",%d,\"%s/workspace_%s_%s.root\",\"%s\",%d,%f,%d,%d,%d,%d,%d,\"%s\",\"%s\",\"%s\")"%(options.templateFile,options.category,options.outputDIR,cat,mass,options.observable,isHiggsInvisible,options.scaleQCD,connectWZ,connectTop,addShapeSystematics,mergeLeptons,isCombination,"",mass,""))
            
                os.system("mkdir -p "+options.outputDIR);
                ROOT.gROOT.ProcessLine(command.Data());
            else:
                
                os.system("mkdir -p "+options.jobDIR);
                os.system("mkdir -p "+currentDIR+"/"+options.inputDIR+"/"+options.outputDIR);
                
                jobName = "job_%s_%s_hinv"%(cat,mass);

                ## write job sh file                                                                                                                                      
                jobmacro = open('%s/%s.C'%(options.jobDIR,jobName),'w')
                jobmacro.write("{\n");
                jobmacro.write("gROOT->ProcessLine(\".L $CMSSW_BASE/lib/$SCRAM_ARCH/libHiggsAnalysisCombinedLimit.so\");\n");
                jobmacro.write("gROOT->ProcessLine(\".L "+currentDIR+"/macros/createWorkspace.C\");\n");
                command = ROOT.TString("\"createWorkspace(\\\"%s\\\",%d,\\\"workspace_%s_%s.root\\\",\\\"%s\\\",%d,%f,%d,%d,%d,%d,%d,\\\"%s\\\",\\\"%s\\\",\\\"%s\\\")\""%(options.templateFile,options.category,cat,mass,options.observable,isHiggsInvisible,options.scaleQCD,connectWZ,connectTop,addShapeSystematics,mergeLeptons,isCombination,"",mass,""))
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
                jobscript.write('scp workspace_%s_%s.root %s/%s/%s\n'%(cat,mass));
                os.system('chmod a+x %s/%s.sh'%(options.jobDIR,jobName))

                if options.submit:
                    os.system('bsub -q %s -o %s/%s/%s.log -e %s/%s/%s.err %s/%s/%s.sh'%(options.queque,currentDIR,options.jobDIR,jobName,currentDIR,options.jobDIR,jobName,currentDIR,options.jobDIR,jobName))

