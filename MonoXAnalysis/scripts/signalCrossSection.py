#! /usr/bin/env python

import ROOT as r
import argparse,commands,os

aparser = argparse.ArgumentParser(description='Process benchmarks.')
aparser.add_argument('-proc' ,'--proc'  ,action='store' ,dest='proc',default='800',  help='800,801,805,806 => V,A,S,P')
aparser.add_argument('-med'  ,'--med'   ,action='store' ,dest='med' ,default='500',  help='med mass')
aparser.add_argument('-dm'   ,'--dm'    ,action='store' ,dest='dm'  ,default='10',    help='dm mass')
aparser.add_argument('-vid'  ,'--vid'   ,action='store' ,dest='vid' ,default=0,      help='boson id (0,23,24)')
aparser.add_argument('-gSM'  ,'--gSM'   ,action='store' ,dest='gSM' ,default=1,      help='gSM')

aparser.add_argument('-list' ,'--list'  ,action='store_true',dest='list'  ,  help='list everything')
options = aparser.parse_args()

entries_sum = 0
entries_inc = 0
entries_ex1 = 0
entries_ex2 = 0

eos='/afs/cern.ch/project/eos/installation/cms/bin/eos.select'
basedir='eos/cms/store/cmst3/group/monojet/mc/model3/'

def getXS(iMed,iId,basedir='/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_7_3_3/src/genproductions/bin/JHUGen/'):
    lFile  = r.TFile(basedir+'/patches/WZXS.root')
    label="ZH"
    if int(iId) == 24:
       label="WH"
    lG     = lFile.Get(label)
    lScale = lFile.Get("scaleUp")
    lBR    = lFile.Get("BRbb")
    if iMed > 500:
       lBR    = lFile.Get("BRtt")
    
    scale=int(iMed)+91+15 # 15 is an approximation of the extra energy based on matching xsections at 125                                                                  
    if int(iId) == 24:
       scale=int(iMed)+80+15
    #Correct for the BR to fermions assuming Scalar decays to bosons                                                                                                       
    lBaseMass=4.2
    if iMed > 500:
       lBaseMass=172.5
    BRCorr = min(lBR.Eval(iMed)*246.*246./lBaseMass/lBaseMass,1.)
    return lG.Eval(iMed)*lScale.Eval(scale)*BRCorr

def computeXS(med,dm,proc,gSM):
    global entries_sum,entries_inc,entries_ex1,entries_ex2 

    infile='%s/MonoJ_%s_%s_%s_%s.root' % (basedir,med,dm,gSM,proc)
    lFile  = r.TFile().Open("root://eoscms//%s" % (infile))
    lTree  = lFile.Get("Events")
    lWHist = r.TH1F("W","W",1,-100000,10000000)
    lZHist = r.TH1F("Z","Z",1,-100000,10000000)
    lXHist = r.TH1F("X","X",1,-100000,10000000)
    lTree.Draw("v_pt>>Z","evtweight*(jpt_1 > 0)","goff")
    lTree.Draw("v_pt>>W","evtweight*(v_pt > 200)","goff")
    lTree.Draw("v_pt>>X","evtweight*(v_pt > 500)","goff")
    #print sumxs,sumweights,sumentries
    entries_inc += lZHist.Integral()
    entries_ex1 += lWHist.Integral()
    entries_ex2 += lXHist.Integral()
    entries_sum += lTree.GetEntries()

def computeXSV(med,dm,proc,iId,gSM):
    global entries_sum,entries_inc,entries_ex1,entries_ex2 
    infile='%s/MonoV_%s_%s_%s_%s.root' % (basedir,med,dm,gSM,proc)
    lFile  = r.TFile().Open("root://eoscms//%s" % (infile))
    lTree  = lFile.Get("Events")
    lZHist = r.TH1F("Z","W",1,-100000,10000000)
    lWHist = r.TH1F("W","W",1,-100000,10000000)
    lXHist = r.TH1F("X","X",1,-100000,10000000)
    bosonid="(abs(v_id) == "+str(iId)+")"
    weight1="xs"
    if int(proc) < 802:
         weight1=weight1+"2"
    weight1+="*"
    lTree.Draw("v_pt>>Z",weight1+"(jpt_1 > 40)*"+bosonid,"goff")
    lTree.Draw("v_pt>>W",weight1+"(v_pt > 200)*"+bosonid,"goff")
    lTree.Draw("v_pt>>X",weight1+"(v_pt > 500)*"+bosonid,"goff")
    #print sumxs,sumweights,sumentries
    xs=1
    if int(proc) > 802:
        xs=getXS(float(med),iId)
    elif int(iId) == 23:
        xs=float(lTree.GetEntries(bosonid))/float(lTree.GetEntries())
    entries_ex1 += lWHist.Integral()*xs
    entries_ex2 += lXHist.Integral()*xs
    entries_sum += lTree.GetEntries(bosonid)
    entries_inc += lZHist.Integral()*xs

if __name__ == '__main__':

    if options.list: 
        command = '%s ls eos/cms/store/cmst3/group/monojet/mc/model3/ | sed "s@_@ @g" | awk \'{print "mMed="$2" mDM="$3}\' | uniq' % eos
        if options.vid > 0:
            command = '%s ls eos/cms/store/cmst3/group/monojet/mc/model3/ | grep MonoV | sed "s@_@ @g" | awk \'{print "mMed="$2" mDM="$3}\' | uniq' % eos
        exists = commands.getoutput(command)
        for line in exists.splitlines():
            print line
        quit()

    if int(options.vid) == 0:
        computeXS(options.med,options.dm,options.proc,options.gSM)

    if int(options.vid) > 0:
        computeXSV(options.med,options.dm,options.proc,options.vid,options.gSM)
    
    print options.proc,"-",options.vid,"-",options.med,"-",options.dm,"XS (inclusive) ",entries_inc/entries_sum," XS (Met > 200):",entries_ex1/entries_sum,"XS (Met > 500):",entries_ex2/entries_sum
