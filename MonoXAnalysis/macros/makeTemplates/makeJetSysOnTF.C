#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#include "makehist.h"
#include "makeCorrHistograms.C"
#include "makeTemplatesUtils.h"

static SamplesNLO nloSamples (false,false,false,false);

void makeJetSysOnTF(Category category, // define selections
		    float    lumi, // define lumi
		    string   outDir, // outputdirectory
		    vector<string> observables, // observable
		    vector<string> observables_2D // observable 2D		    
		    ){

  system(("mkdir -p "+outDir).c_str());
  // to initialize the binning map                                                                                                                                                                    
  initializeBinning();

  //////////////////////////                                                                                                                                                                           
  nloSamples.WJetsDIR = "WJets";
  nloSamples.ZJetsDIR = "ZJets";
  nloSamples.DYJetsDIR = "DYJets";
  nloSamples.PhotonJetsDIR = "PhotonJets";

  if(nloSamples.useWJetsNLO)
    nloSamples.WJetsDIR = "WJetsNLO";
  if(nloSamples.useZJetsNLO)
    nloSamples.ZJetsDIR = "ZJetsNLO";
  if(nloSamples.useDYJetsNLO)
    nloSamples.DYJetsDIR = "DYJetsNLO";
  if(nloSamples.usePhotonJetsNLO)
    nloSamples.PhotonJetsDIR = "PhotonJetsNLO";

  
  cout<<"make correction histogram for Zmm to Znn"<<endl;
  makezmmcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		 baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/",
		 category,nloSamples,observables,observables_2D,lumi,outDir,"",true,false);

  cout<<"make correction histogram for Zmm to Znn --> jesUp"<<endl;
  makezmmcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		 baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/",
		 category,nloSamples,observables,observables_2D,lumi,outDir,"jesUp",true,false,"jesUp");

  cout<<"make correction histogram for Zmm to Znn --> jesDw"<<endl;
  makezmmcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		 baseInputTreePath+"/"+nloSamples.DYJetsDIR+"/zmmfilter/",
		 category,nloSamples,observables,observables_2D,lumi,outDir,"jesDw",true,false,"jesDw");

  cout<<"make correction histogram for Wmn to Wln"<<endl;
  makewmncorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR +"/sigfilter/",
		 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/",
		 category,nloSamples,observables,observables_2D,lumi,outDir,"",true,false);

  cout<<"make correction histogram for Wmn to Wln --> jesUp"<<endl;
  makewmncorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR +"/sigfilter/",
		 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/",
		 category,nloSamples,observables,observables_2D,lumi,outDir,"jesUp",true,false,"jesUp");

  cout<<"make correction histogram for Wmn to Wln --> jesDw"<<endl;
  makewmncorhist(baseInputTreePath+"/"+nloSamples.WJetsDIR +"/sigfilter/",
		 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/wmnfilter/",
		 category,nloSamples,observables,observables_2D,lumi,outDir,"jesDw",true,false,"jesDw");

  cout<<"make correction histogram for Znn to Wln"<<endl;
  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		 category,nloSamples,observables,observables_2D,lumi,outDir,"",true,false,"",3);

  cout<<"make correction histogram for Znn to Wln --> jesUp"<<endl;
  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		 category,nloSamples,observables,observables_2D,lumi,outDir,"jesUp",true,false,"jesUp",3);
  
  cout<<"make correction histogram for Znn to Wln --> jesDw"<<endl;
  makezwjcorhist(baseInputTreePath+"/"+nloSamples.ZJetsDIR +"/sigfilter/",
		 baseInputTreePath+"/"+nloSamples.WJetsDIR+"/sigfilter/",
		 category,nloSamples,observables,observables_2D,lumi,outDir,"jesDw",true,false,"jesDw",3);
  
}
