//monoV_signalAnalysis("/store/ciMassaf/user/rgerosa/MONOJET_ANALYSIS/Production-14-1-2016/",1,"Pseudoscalar","Plots_MonoW_Pseudoscalar")
#include <iostream>
#include <vector>
#include <string>
#include <string>
#include <fstream>
#include <unordered_map>

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TChain.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

using namespace std;

vector<string> Mphi_Vector        = {"100","200","500","1000"};
vector<string> Mphi_Axial         = {"100","200","500","1000"};
vector<string> Mphi_PseudoScalar  = {"100","200","500","1000"};
vector<string> Mphi_Scalar        = {"100","125","1000","2000"};

vector<string> Mchi_Vector       = {"1","10","10","50"};
vector<string> Mchi_Axial        = {"1","10","50","50"};
vector<string> Mchi_PseudoScalar = {"1","10","10","50"};
vector<string> Mchi_Scalar       = {"50","10","10","10"};

vector<pair<string,string> > selectionsSig;
vector<pair<string,string> > selectionsWmn;
vector<pair<string,string> > selectionsWen;
vector<pair<string,string> > selectionsZmm;
vector<pair<string,string> > selectionsZee;


void makeEfficiencyPlots(TCanvas* cCanvas, unordered_map<string,TH1F*> eff_sig, 
			 vector<string> Mphi, vector<string> Mchi, 
			 string bosonMode, string interactionModel, 
			 TLegend* legend, string outputDirectory, 
			 string plotName, bool isRelative);

void makeShapePlots(TCanvas* cCanvas, vector<TH1F*> histoList, TLegend* legend, string outputDirectory, string plotName, string xAxisTile);
void makeShapePlots(TCanvas* cCanvas, vector<TH1F*> histoList, TH1F* histBkg, TLegend* legend, string outputDirectory, string plotName, string xAxisTile);

// make signal efficiency studies
void makeMonoVSignalEfficiency(string baseInputPath,  // base path with all the directories on the local pc
			       bool isW, // to pick up MonoW or MonoZ directories 
			       string interactionModel, // interaction 
			       string outputDirectory // output directory for plots
			       ){

  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetErrorX(0.5);
  gStyle->SetOptTitle(0);

  // take input files
  string bosonMode;
  if(isW)
    bosonMode = "MonoW";
  else 
    bosonMode = "MonoZ";
  

  // take vectors for Mphi and Mchi
  vector<string> Mphi;
  vector<string> Mchi;
  
  if(interactionModel == "Vector"){
    Mphi = Mphi_Vector;
    Mchi = Mchi_Vector;
  }
  else if(interactionModel == "Axial"){
    Mphi = Mphi_Axial;
    Mchi = Mchi_Axial;
  }
  else if(interactionModel == "Pseudoscalar"){
    Mphi = Mphi_PseudoScalar;
    Mchi = Mchi_PseudoScalar;
  }
  else if(interactionModel == "Scalar"){
    Mphi = Mphi_Scalar;
    Mchi = Mchi_Scalar;
  }

  //create output directories
  system(("mkdir -p "+outputDirectory).c_str());
  system(("mkdir -p "+outputDirectory+"/efficiency/").c_str());
  system(("mkdir -p "+outputDirectory+"/shapes/").c_str());

  //////
  TCanvas *cCanvas = new TCanvas("cCanvas","",180,52,550,550);
  cCanvas->SetTicks();
  cCanvas->SetFillColor(0);
  cCanvas->SetBorderMode(0);
  cCanvas->SetBorderSize(2);
  cCanvas->SetRightMargin(0.05);
  cCanvas->SetBottomMargin(0.12);
  cCanvas->SetFrameBorderMode(0);
  
  TLegend* legend = new TLegend(0.25,0.68,0.7,0.89);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.031);
  legend->SetTextFont(42);

  // start taking input files before the skim
  vector<TFile*> fileMonoV;  
  vector<pair<string,TTree*> > treeMonoV;  

  system(("eos ls /store/caf/user/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/"+bosonMode+"_"+interactionModel+" > file_list.txt").c_str());
  ifstream infile;
  string line;
  infile.open("file_list.txt",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);	
      if(line == "") continue;
      // find matching with Mphi, Mchi
      for(size_t iMass = 0; iMass < Mphi.size(); iMass++){	  
	if(TString(line).Contains(("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_").c_str())){
	  fileMonoV.push_back(TFile::Open(("root://eoscms.cern.ch///store/caf/user/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/"+bosonMode+"_"+interactionModel+"/"+line).c_str()));
	  if(fileMonoV.back() and not fileMonoV.back()->IsZombie())
	    treeMonoV.push_back(make_pair("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass],(TTree*) fileMonoV.back()->Get("tree/tree")));
	}
      }
    }
    infile.close();
  }
  
  // take signal region files for signal
  vector<TFile*> fileMonoV_Sig;  
  vector<pair<string,TTree*> > treeMonoV_Sig;  
  system(("ls "+baseInputPath+"/"+bosonMode+"_"+interactionModel+"/sigfilter/ > file_list.txt").c_str());
  infile.open("file_list.txt",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "") continue;
      // find matching with Mphi, Mchi                                                                                                                                         
      for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
	if(TString(line).Contains(("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_").c_str())){
	  cout<<baseInputPath+"/"+bosonMode+"_"+interactionModel+"/sigfilter/"+line<<endl;
	  fileMonoV_Sig.push_back(TFile::Open((baseInputPath+"/"+bosonMode+"_"+interactionModel+"/sigfilter/"+line).c_str()));
	  if(fileMonoV_Sig.back() and not fileMonoV_Sig.back()->IsZombie())
	    treeMonoV_Sig.push_back(make_pair("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass],(TTree*) fileMonoV_Sig.back()->Get("tree/tree")));
	}
      }
    }
  }

  infile.close();

  // start taking input files in the signal region  
  vector<TFile*> fileMonoV_CRm;  
  vector<pair<string,TTree*> > treeMonoV_CRm;  
  vector<TFile*> fileMonoV_CRe;  
  vector<pair<string,TTree*> > treeMonoV_CRe;  

  if(isW){	
    system(("eos ls /store/caf/user/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/"+bosonMode+"_"+interactionModel+"/wmnfilter/ > file_list.txt").c_str());
    infile.open("file_list.txt",ifstream::in);
    if(infile.is_open()){
      while(!infile.eof()){
	getline(infile,line);	  
	if(line == "") continue;
	for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
	  if(TString(line).Contains(("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_").c_str())){
	    fileMonoV_CRm.push_back(TFile::Open(("root://eoscms.cern.ch///store/caf/user/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/"+bosonMode+"_"+interactionModel+"/wmnfilter/"+line).c_str()));
	    if(fileMonoV_CRm.back() and not fileMonoV_CRm.back()->IsZombie())
	      treeMonoV_CRm.push_back(make_pair("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass],(TTree*) fileMonoV_CRm.back()->Get("tree/tree")));      
	  }
	}
      }
    }
    infile.close();
    
    system(("eos ls /store/caf/user/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/"+bosonMode+"_"+interactionModel+"/wenfilter/ > file_list.txt").c_str());
    infile.open("file_list.txt",ifstream::in);
    if(infile.is_open()){
      while(!infile.eof()){
	getline(infile,line);	  
	if(line == "") continue;
	for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
	  if(TString(line).Contains(("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_").c_str())){
	    fileMonoV_CRe.push_back(TFile::Open(("root://eoscms.cern.ch///store/caf/user/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/"+bosonMode+"_"+interactionModel+"/wenfilter/"+line).c_str()));
	    if(fileMonoV_CRe.back() and not fileMonoV_CRe.back()->IsZombie())
	      treeMonoV_CRe.push_back(make_pair("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass],(TTree*) fileMonoV_CRe.back()->Get("tree/tree")));      
	  }
	}
      }
    }
    infile.close();
  }	
  else{
    
    system(("eos ls /store/caf/user/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/"+bosonMode+"_"+interactionModel+"/zmmfilter/ > file_list.txt").c_str());
    infile.open("file_list.txt",ifstream::in);
    if(infile.is_open()){
      while(!infile.eof()){
	getline(infile,line);	  
	if(line == "") continue;
	for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
	  if(TString(line).Contains(("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_").c_str())){
	    fileMonoV_CRm.push_back(TFile::Open(("root://eoscms.cern.ch///store/caf/user/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/"+bosonMode+"_"+interactionModel+"/zmmfilter/"+line).c_str()));
	    if(fileMonoV_CRm.back() and not fileMonoV_CRm.back()->IsZombie())
	      treeMonoV_CRm.push_back(make_pair("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass],(TTree*) fileMonoV_CRm.back()->Get("tree/tree")));      
	  }
	}
      }
    }
    infile.close();
    system(("eos ls "+baseInputPath+"/"+bosonMode+"_"+interactionModel+"/zeefilter/ > file_list.txt").c_str());
    infile.open("file_list.txt",ifstream::in);
    if(infile.is_open()){
      while(!infile.eof()){
	getline(infile,line);	  
	if(line == "") continue;
	for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
	  if(TString(line).Contains(("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_").c_str())){
	    fileMonoV_CRe.push_back(TFile::Open(("root://eoscms.cern.ch///store/caf/user/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/"+bosonMode+"_"+interactionModel+"/zeefilter/"+line).c_str()));
	    if(fileMonoV_CRe.back() and not fileMonoV_CRe.back()->IsZombie())
	      treeMonoV_CRe.push_back(make_pair("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass],(TTree*) fileMonoV_CRe.back()->Get("tree/tree")));      
	  }
	}
      }
    }
    infile.close();      
  }
  
  system("rm file_list.txt");
  
  // define signal region cuts for efficiency                                                                                                                                
  selectionsSig.push_back(make_pair("trigger","hltmet90"));
  selectionsSig.push_back(make_pair("njets","hltmet90 && njets >= 1"));
  selectionsSig.push_back(make_pair("jetpt","hltmet90 && njets >= 1 && centraljetpt[0] > 100"));
  selectionsSig.push_back(make_pair("bveto","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0"));
  selectionsSig.push_back(make_pair("jetid","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8"));
  selectionsSig.push_back(make_pair("dphijemet","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5"));
  selectionsSig.push_back(make_pair("met","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1pfmet > 200"));
  selectionsSig.push_back(make_pair("ak8pt","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1pfmet > 250 && boostedJetpt[0] > 250"));

  // define signal region cuts for efficiency
  selectionsSig.push_back(make_pair("trigger","hltmet90"));
  selectionsSig.push_back(make_pair("njets","hltmet90 && njets >= 1"));
  selectionsSig.push_back(make_pair("jetpt","hltmet90 && njets >= 1 && centraljetpt[0] > 100"));
  selectionsSig.push_back(make_pair("bveto","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0"));
  selectionsSig.push_back(make_pair("jetid","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8"));
  selectionsSig.push_back(make_pair("dphijemet","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5"));
  selectionsSig.push_back(make_pair("met","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1pfmet > 200"));
  selectionsSig.push_back(make_pair("ak8pt","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1pfmet > 250 && boostedJetpt[0] > 250"));
  
  // define Wmn cuts
  selectionsWmn.push_back(make_pair("trigger", "hltmet90"));
  selectionsWmn.push_back(make_pair("njets", "hltmet90 && njets >= 1"));
  selectionsWmn.push_back(make_pair("jetpt", "hltmet90 && njets >= 1 && centraljetpt[0] > 100"));
  selectionsWmn.push_back(make_pair("bveto", "hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0"));
  selectionsWmn.push_back(make_pair("jetid", "hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8"));
  selectionsWmn.push_back(make_pair("dphijemet", "hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5"));
  selectionsWmn.push_back(make_pair("met", "hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1mumet > 200"));
  
  // define Wen cuts
  selectionsWen.push_back(make_pair("trigger","hltsingleel"));
  selectionsWen.push_back(make_pair("njets","hltsingleel && njets >= 1"));
  selectionsWen.push_back(make_pair("jetpt", "hltsingleel && njets >= 1 && centraljetpt[0] > 100"));
  selectionsWen.push_back(make_pair("bveto", "hltsingleel && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0"));
  selectionsWen.push_back(make_pair("jetid", "hltsingleel && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8"));
  selectionsWen.push_back(make_pair("dphijemet", "hltsingleel && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5"));
  selectionsWen.push_back(make_pair("met","hltsingleel && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetelmetdphimin4 > 0.5 && t1elmet > 200"));
  
  // define Zmm cuts
  selectionsZmm.push_back(make_pair("trigger","hltmet90"));
  selectionsZmm.push_back(make_pair("njets","hltmet90 && njets >= 1"));
  selectionsZmm.push_back(make_pair("jetpt","hltmet90 && njets >= 1 && centraljetpt[0] > 100"));
  selectionsZmm.push_back(make_pair("bveto","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0"));
  selectionsZmm.push_back(make_pair("jetid","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8"));
  selectionsZmm.push_back(make_pair("dphijemet","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5"));
  selectionsZmm.push_back(make_pair("met","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1mumet > 200"));
  selectionsZmm.push_back(make_pair("lepcharge","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1mumet > 200 && mu1pid != mu2pid"));
  
  // define Zee cuts
  selectionsZee.push_back(make_pair("trigger","hltsingleel"));
  selectionsZee.push_back(make_pair("njets", "hltsingleel && njets >= 1"));
  selectionsZee.push_back(make_pair("jetpt", "hltsingleel && njets >= 1 && centraljetpt[0] > 100"));
  selectionsZee.push_back(make_pair("bveto", "hltsingleel && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0"));
  selectionsZee.push_back(make_pair("jetid", "hltsingleel && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8"));
  selectionsZee.push_back(make_pair("dphijemet", "hltsingleel && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5"));
  selectionsZee.push_back(make_pair("met","hltsingleel && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetelmetdphimin4 > 0.5 && t1elmet > 200"));
  selectionsZee.push_back(make_pair("lepcharge","hltsingleel && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1elmet > 200 && el1pid != el2pid"));

  // calculate efficiency of basic control region selection
  unordered_map<string,TH1F*> eff_sig;
  unordered_map<string,TH1F*> eff_CRm;
  unordered_map<string,TH1F*> eff_CRe;
  
  unordered_map<string,TH1F*> eff_sig_relative;
  unordered_map<string,TH1F*> eff_CRm_relative;
  unordered_map<string,TH1F*> eff_CRe_relative;
  
  TH1F* temp_CRm;
  TH1F* temp_CRe;
  TH1F* temp_CRm_relative;
  TH1F* temp_CRe_relative;
  
  for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
    
    ////////// declare histograms for efficiency ////////
    eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] = 
      new TH1F(("eff_sig_"+bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"",selectionsSig.size()+1,0,selectionsSig.size()+1);	
    eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> SetBinContent(1,float(treeMonoV_Sig.at(iMass).second->GetEntries())/treeMonoV.at(iMass).second->GetEntries());
    eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> GetXaxis()->SetBinLabel(1,"Preselection");
    
    eff_sig_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] = 
      new TH1F(("eff_sig_rel_"+bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"",selectionsSig.size()+1,0,selectionsSig.size()+1);	
    eff_sig_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> SetBinContent(1,float(treeMonoV_Sig.at(iMass).second->GetEntries())/treeMonoV.at(iMass).second->GetEntries());
    eff_sig_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> GetXaxis()->SetBinLabel(1,"Preselection");
    
    /////////////////////
    if(isW){
      temp_CRm = new TH1F(("eff_wmn_"+bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"",selectionsWmn.size()+1,0,selectionsWmn.size()+1);
      temp_CRe = new TH1F(("eff_wen_"+bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"",selectionsWen.size()+1,0,selectionsWen.size()+1);
      
      temp_CRm_relative = new TH1F(("eff_wmn_rel_"+bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"",selectionsWmn.size()+1,0,selectionsWmn.size()+1);
      temp_CRe_relative = new TH1F(("eff_wen_rel_"+bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"",selectionsWen.size()+1,0,selectionsWen.size()+1);
      
    }    
    else{
      temp_CRm = new TH1F(("eff_zmm_"+bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"",selectionsZmm.size()+1,0,selectionsZmm.size()+1);
      temp_CRe = new TH1F(("eff_zee_"+bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"",selectionsZee.size()+1,0,selectionsZee.size()+1);
      
      temp_CRm_relative = new TH1F(("eff_wmn_rel_"+bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"",selectionsZmm.size()+1,0,selectionsZmm.size()+1);
      temp_CRe_relative = new TH1F(("eff_wen_rel_"+bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"",selectionsZee.size()+1,0,selectionsZee.size()+1);
    }   
    
    eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] = temp_CRm;
    eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] = temp_CRe;
    
    eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] = temp_CRm_relative;
    eff_CRe_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] = temp_CRe_relative;
    
    // set first bin    
    eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> SetBinContent(1,float(treeMonoV_CRm.at(iMass).second->GetEntries())/treeMonoV.at(iMass).second->GetEntries());
    eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> GetXaxis()->SetBinLabel(1,"Preselection");
    
    eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> SetBinContent(1,float(treeMonoV_CRe.at(iMass).second->GetEntries())/treeMonoV.at(iMass).second->GetEntries());
    eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> GetXaxis()->SetBinLabel(1,"Preselection");
    
    eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> SetBinContent(1,float(treeMonoV_CRm.at(iMass).second->GetEntries())/treeMonoV.at(iMass).second->GetEntries());
    eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> GetXaxis()->SetBinLabel(1,"Preselection");
    
    eff_CRe_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> SetBinContent(1,float(treeMonoV_CRe.at(iMass).second->GetEntries())/treeMonoV.at(iMass).second->GetEntries());
    eff_CRe_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> GetXaxis()->SetBinLabel(1,"Preselection");
    
    // apply signal selections
    int iselection = 1;
    for(auto iSel = selectionsSig.begin(); iSel != selectionsSig.end(); iSel++){
      iselection++;
      eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_Sig.at(iMass).second->GetEntries(iSel->second.c_str()))/treeMonoV.at(iMass).second->GetEntries());
      eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
      
      eff_sig_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_Sig.at(iMass).second->GetEntries(iSel->second.c_str()))/(treeMonoV.at(iMass).second->GetEntries())*1/eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetBinContent(iselection-1));
      eff_sig_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
    }
    
    if(isW){
      iselection = 1;
      for(auto iSel = selectionsWmn.begin(); iSel != selectionsWmn.end(); iSel++){
	iselection ++;
	eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRm.at(iMass).second->GetEntries(iSel->second.c_str()))/treeMonoV.at(iMass).second->GetEntries());
	eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
	  
	eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRm.at(iMass).second->GetEntries(iSel->second.c_str()))/(treeMonoV.at(iMass).second->GetEntries())*1/eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetBinContent(iselection-1));
	eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
      }
      
      iselection = 1;
      for(auto iSel = selectionsWen.begin(); iSel != selectionsWen.end(); iSel++){
	iselection ++;
	eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRe.at(iMass).second->GetEntries(iSel->second.c_str()))/treeMonoV.at(iMass).second->GetEntries());
	eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
	
	eff_CRe_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRe.at(iMass).second->GetEntries(iSel->second.c_str()))/(treeMonoV.at(iMass).second->GetEntries())*1/eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetBinContent(iselection-1));
	eff_CRe_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
      }      
    }
    else{	
      iselection = 1;
      for(auto iSel = selectionsZmm.begin(); iSel != selectionsZmm.end(); iSel++){
	iselection ++;
	eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRm.at(iMass).second->GetEntries(iSel->second.c_str()))/treeMonoV.at(iMass).second->GetEntries());
	eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
	
	eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRm.at(iMass).second->GetEntries(iSel->second.c_str()))/(treeMonoV.at(iMass).second->GetEntries())*1/eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetBinContent(iselection-1));
	eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());	    
      }
      
      iselection = 1;
      for(auto iSel = selectionsZee.begin(); iSel != selectionsZee.end(); iSel++){
	iselection ++;
	eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRe.at(iMass).second->GetEntries(iSel->second.c_str()))/treeMonoV.at(iMass).second->GetEntries());
	eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
	
	eff_CRe_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRe.at(iMass).second->GetEntries(iSel->second.c_str()))/(treeMonoV.at(iMass).second->GetEntries())*1/eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetBinContent(iselection-1));
	eff_CRe_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
      }      
    } 
  }
  // plotting efficiencies: signal region
  TFile* outputEfficiency = new TFile((outputDirectory+"/efficiency.root").c_str(),"RECREATE");
  outputEfficiency->cd();
  
  makeEfficiencyPlots(cCanvas,eff_sig_relative,Mphi,Mchi,bosonMode,interactionModel,legend,outputDirectory,"SR_efficiency_rel",true);
  makeEfficiencyPlots(cCanvas,eff_CRm_relative,Mphi,Mchi,bosonMode,interactionModel,legend,outputDirectory,"CRm_efficiency_rel",true);
  makeEfficiencyPlots(cCanvas,eff_CRe_relative,Mphi,Mchi,bosonMode,interactionModel,legend,outputDirectory,"CRe_efficiency_rel",true);
  makeEfficiencyPlots(cCanvas,eff_sig,Mphi,Mchi,bosonMode,interactionModel,legend,outputDirectory,"SR_efficiency",false);
  makeEfficiencyPlots(cCanvas,eff_CRm,Mphi,Mchi,bosonMode,interactionModel,legend,outputDirectory,"CRm_efficiency",false);
  makeEfficiencyPlots(cCanvas,eff_CRe,Mphi,Mchi,bosonMode,interactionModel,legend,outputDirectory,"CRe_efficiency",false);

}
  		       		 
void makeMonoVSignalShape(string baseInputPath,  // base path with all the directories on the local pc
			  bool isW, // to pick up MonoW or MonoZ directories 
			  string interactionModel, // interaction 
			  string outputDirectory, // output directory for plots
			  bool skipBackground = false
			  ){
  
  gROOT->SetBatch(1);
  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetErrorX(0.5);
  gStyle->SetOptTitle(0);

  // take input files
  string bosonMode;
  if(isW)
    bosonMode = "MonoW";
  else 
    bosonMode = "MonoZ";
  

  // take vectors for Mphi and Mchi
  vector<string> Mphi;
  vector<string> Mchi;
  
  if(interactionModel == "Vector"){
    Mphi = Mphi_Vector;
    Mchi = Mchi_Vector;
  }
  else if(interactionModel == "Axial"){
    Mphi = Mphi_Axial;
    Mchi = Mchi_Axial;
  }
  else if(interactionModel == "Pseudoscalar"){
    Mphi = Mphi_PseudoScalar;
    Mchi = Mchi_PseudoScalar;
  }
  else if(interactionModel == "Scalar"){
    Mphi = Mphi_Scalar;
    Mchi = Mchi_Scalar;
  }

  //create output directories
  system(("mkdir -p "+outputDirectory).c_str());
  system(("mkdir -p "+outputDirectory+"/efficiency/").c_str());
  system(("mkdir -p "+outputDirectory+"/shapes/").c_str());

  //////
  TCanvas *cCanvas = new TCanvas("cCanvas","",180,52,550,550);
  cCanvas->SetTicks();
  cCanvas->SetFillColor(0);
  cCanvas->SetBorderMode(0);
  cCanvas->SetBorderSize(2);
  cCanvas->SetRightMargin(0.05);
  cCanvas->SetBottomMargin(0.12);
  cCanvas->SetFrameBorderMode(0);
  
  TLegend* legend = new TLegend(0.25,0.68,0.7,0.89);
  legend->SetBorderSize(0);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.031);
  legend->SetTextFont(42);

  // take signal region files for signal
  vector<TFile*> fileMonoV_Sig;  
  vector<pair<string,TTree*> > treeMonoV_Sig;  
  system(("ls "+baseInputPath+"/"+bosonMode+"_"+interactionModel+"/sigfilter/ > file_list.txt").c_str());
  ifstream infile;
  string line;
  infile.open("file_list.txt",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "") continue;
      // find matching with Mphi, Mchi                                                                                                                                         
      for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
	if(TString(line).Contains(("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_").c_str())){
	  cout<<baseInputPath+"/"+bosonMode+"_"+interactionModel+"/sigfilter/"+line<<endl;
	  fileMonoV_Sig.push_back(TFile::Open((baseInputPath+"/"+bosonMode+"_"+interactionModel+"/sigfilter/"+line).c_str()));
	  if(fileMonoV_Sig.back() and not fileMonoV_Sig.back()->IsZombie())
	    treeMonoV_Sig.push_back(make_pair("Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass],(TTree*) fileMonoV_Sig.back()->Get("tree/tree")));
	}
      }
    }
  }

  infile.close();
  system("rm file_list.txt");

  //Background Znunu
  TChain* backgroundZnunu = new TChain("tree/tree");
  TFile* file1 = TFile::Open((baseInputPath+"/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root").c_str());
  backgroundZnunu->Add(file1->GetName());

  string cut = "hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1pfmet > 250 && boostedJetpt[0] > 250";

  /// now look at substructure in the signal region only
  vector<TH1F*> leadingJetPt;
  vector<TH1F*> leadingJetEta;
  vector<TH1F*> leadingJetTau2tau1;
  vector<TH1F*> leadingPrunedJetpt;
  vector<TH1F*> leadingPrunedJetm;
  vector<TH1F*> leadingPrunedJetmraw;
  vector<TH1F*> leadingSoftDropJetpt;
  vector<TH1F*> leadingSoftDropJetm;
  vector<TH1F*> leadingSoftDropJetmraw;
  vector<TH1F*> dRGenJet;
  vector<TH1F*> dRVboson;
  vector<TH1F*> subjetPtRatio;
  vector<TH1F*> responseLeadJetPt;
  vector<TH1F*> responseLeadPrunedJetMass;
  vector<TH1F*> responseLeadPrunedJetMassRaw;
  vector<TH1F*> responseLeadSoftDropJetMass;
  vector<TH1F*> responseLeadSoftDropJetMassRaw;
  vector<TH1F*> boostedNJet;
  vector<TH1F*> subjetPtRatio_pr;
  vector<TH1F*> responseLeadPrunedJetMass_pr;
  vector<TH1F*> responseLeadPrunedJetMassRaw_pr;
  vector<TH1F*> leadingJetTau2Tau1_pr;
  vector<TH1F*> leadingJetTau2Tau1_sb;

  cout<<"Run signal analysis "<<endl;
  for(size_t iMass = 0; iMass < treeMonoV_Sig.size(); iMass++){

    string name = treeMonoV_Sig.at(iMass).first;

    TTreeFormula* formula = new TTreeFormula((bosonMode+"_"+interactionModel+"_"+name+"_formula").c_str(),
					     cut.c_str(),treeMonoV_Sig.at(iMass).second);
    
    // allocate all the histograms
    cout<<" histogram name "<<bosonMode+"_"+interactionModel+"_"+name<<endl;
    leadingJetPt.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingJetPt").c_str(),
				    (bosonMode+"_"+interactionModel+"_"+name).c_str(),25,200,650));

    leadingJetEta.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingJetEta").c_str(),
				     (bosonMode+"_"+interactionModel+"_"+name).c_str(),20,-2.5,2.5));
    leadingJetTau2tau1.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingJetTau2tau1").c_str(),
					  (bosonMode+"_"+interactionModel+"_"+name).c_str(),20,0,1));    
    leadingPrunedJetpt.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingPrunedJetpt").c_str(),
					  (bosonMode+"_"+interactionModel+"_"+name).c_str(),25,200,650));
    leadingPrunedJetm.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingPrunedJetm").c_str(),
					 (bosonMode+"_"+interactionModel+"_"+name).c_str(),30,0,120));
    leadingPrunedJetmraw.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingPrunedJetmraw").c_str(),
					    (bosonMode+"_"+interactionModel+"_"+name).c_str(),30,0,120));
    leadingSoftDropJetpt.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingSoftDropJetpt").c_str(),
					    (bosonMode+"_"+interactionModel+"_"+name).c_str(),25,200,650));
    leadingSoftDropJetm.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingSoftDropJetm").c_str(),
					   (bosonMode+"_"+interactionModel+"_"+name).c_str(),30,0,120));
    leadingSoftDropJetmraw.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingSoftDropJetmraw").c_str(),
					      (bosonMode+"_"+interactionModel+"_"+name).c_str(),30,0,120));
    boostedNJet.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_boostedNJet").c_str(),
				   (bosonMode+"_"+interactionModel+"_"+name).c_str(),5,0,5));
    dRGenJet.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_dRGenJet").c_str(),
				(bosonMode+"_"+interactionModel+"_"+name).c_str(),50,0,0.8));
    dRVboson.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_dRVboson").c_str(),
				(bosonMode+"_"+interactionModel+"_"+name).c_str(),50,0,0.8));
    responseLeadJetPt.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_responseLeadJetPt").c_str(),
					 (bosonMode+"_"+interactionModel+"_"+name).c_str(),30,0.5,1.5));
    responseLeadPrunedJetMass.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_responseLeadPrunedJetMass").c_str(),
						 (bosonMode+"_"+interactionModel+"_"+name).c_str(),30,0.5,1.5));
    responseLeadPrunedJetMassRaw.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_responseLeadPrunedJetMassRaw").c_str(),
						    (bosonMode+"_"+interactionModel+"_"+name).c_str(),30,0.5,1.5));
    responseLeadSoftDropJetMass.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_responseLeadSoftDropJetMass").c_str(),
						   (bosonMode+"_"+interactionModel+"_"+name).c_str(),30,0.5,1.5));
    responseLeadSoftDropJetMassRaw.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_responseLeadSoftDropJetMassRaw").c_str(),
						      (bosonMode+"_"+interactionModel+"_"+name).c_str(),30,0.5,1.5));
    subjetPtRatio.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_subjetPtRatio").c_str(),
				     (bosonMode+"_"+interactionModel+"_"+name).c_str(),20,0,1));
    // after pruned mass cut
    responseLeadPrunedJetMass_pr.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_responseLeadPrunedJetMass_pr").c_str(),
						    (bosonMode+"_"+interactionModel+"_"+name).c_str(),30,0.5,1.5));
    responseLeadPrunedJetMassRaw_pr.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_responseLeadPrunedJetMassRaw_pr").c_str(),
						       (bosonMode+"_"+interactionModel+"_"+name).c_str(),30,0.5,1.5));
    subjetPtRatio_pr.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_subjetPtRatio_PR").c_str(),
					(bosonMode+"_"+interactionModel+"_"+name).c_str(),30,0,2));
    leadingJetTau2Tau1_pr.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingJetTau2Tau1_pr").c_str(),
					     (bosonMode+"_"+interactionModel+"_"+name).c_str(),20,0,1));

    //sideband region
    leadingJetTau2Tau1_sb.push_back(new TH1F((bosonMode+"_"+interactionModel+"_"+name+"_leadingJetTau2Tau1_sb").c_str(),
					     (bosonMode+"_"+interactionModel+"_"+name).c_str(),20,0,1));
    // make cut formula
    string cut_pr = cut + " && prunedJetm[0] > 65 && prunedJetm[0] < 105";
    TTreeFormula* formula_pr = new TTreeFormula((bosonMode+"_"+interactionModel+"_"+name+"_formula_pr").c_str(),
						cut_pr.c_str(),treeMonoV_Sig.at(iMass).second);
    string cut_sb = cut + " && prunedJetm[0] > 40 && prunedJetm[0] < 65";
    TTreeFormula* formula_sb = new TTreeFormula((bosonMode+"_"+interactionModel+"_"+name+"_formula_sb").c_str(),
						cut_pr.c_str(),treeMonoV_Sig.at(iMass).second);
    
    TTreeReader myReader(treeMonoV_Sig.at(iMass).second);
    TTreeReaderValue<vector<double> > jetpt    (myReader,"boostedJetpt");
    TTreeReaderValue<vector<double> > jeteta   (myReader,"boostedJeteta");
    TTreeReaderValue<vector<double> > jetphi   (myReader,"boostedJetphi");
    TTreeReaderValue<vector<double> > jetm     (myReader,"boostedJetm");
    TTreeReaderValue<vector<double> > jetpt_gen    (myReader,"boostedJetGenpt");
    TTreeReaderValue<vector<double> > jeteta_gen   (myReader,"boostedJetGeneta");
    TTreeReaderValue<vector<double> > jetphi_gen   (myReader,"boostedJetGenphi");
    TTreeReaderValue<vector<double> > jetm_gen   (myReader,"boostedJetGenm");

    TTreeReaderValue<vector<double> > jettau2  (myReader,"boostedJettau2");
    TTreeReaderValue<vector<double> > jettau1  (myReader,"boostedJettau1");

    TTreeReaderValue<vector<double> > jetpt_pr (myReader,"prunedJetpt");
    TTreeReaderValue<vector<double> > jetm_pr  (myReader,"prunedJetm");
    TTreeReaderValue<vector<double> > jetmraw_pr  (myReader,"prunedJetmraw");
    TTreeReaderValue<vector<double> > jetsub1pt_pr  (myReader,"prunedSubJetpt_1");
    TTreeReaderValue<vector<double> > jetsub2pt_pr  (myReader,"prunedSubJetpt_2");

    TTreeReaderValue<vector<double> > jetpt_pr_gen (myReader,"prunedJetGenpt");
    TTreeReaderValue<vector<double> > jetm_pr_gen  (myReader,"prunedJetGenm");

    TTreeReaderValue<vector<double> > jetpt_sd (myReader,"softDropJetpt");
    TTreeReaderValue<vector<double> > jetm_sd  (myReader,"softDropJetm");
    TTreeReaderValue<vector<double> > jetpt_sd_gen (myReader,"softDropJetGenpt");
    TTreeReaderValue<vector<double> > jetm_sd_gen  (myReader,"softDropJetGenm");
    TTreeReaderValue<vector<double> > jetmraw_sd  (myReader,"softDropJetmraw");

    TTreeReaderValue<double> wzpt_h (myReader,"wzpt_h");
    TTreeReaderValue<double> wzeta_h (myReader,"wzeta_h");
    TTreeReaderValue<double> wzphi_h (myReader,"wzphi_h");
    TTreeReaderValue<double> wzmass_h (myReader,"wzmass_h");

    int iEvent = 0;

    while(myReader.Next()){
      float Njet = 0;
      treeMonoV_Sig.at(iMass).second->GetEntry(iEvent);
      iEvent++;

      if(not formula->EvalInstance()) continue;
      
      // AK8 jets with pT > 250 GeV
      for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
	if(jetpt->at(ijet) > 250)
	  Njet++;
      }
      
      // Fill histos
      boostedNJet.back()->Fill(Njet);

      if(jetpt->size()>0)
	leadingJetPt.back()->Fill(jetpt->at(0));
      if(jeteta->size()>0)
	leadingJetEta.back()->Fill(jeteta->at(0));
      if(jettau2->size()>0 and jettau1->size() > 0)
	leadingJetTau2tau1.back()->Fill(jettau2->at(0)/jettau1->at(0));
      if(jetpt_pr->size()>0)
	leadingPrunedJetpt.back()->Fill(jetpt_pr->at(0));
      if(jetm_pr->size()>0)
	leadingPrunedJetm.back()->Fill(jetm_pr->at(0));
      if(jetmraw_pr->size()>0)
	leadingPrunedJetmraw.back()->Fill(jetmraw_pr->at(0));
      if(jetpt_sd->size()>0)
	leadingSoftDropJetpt.back()->Fill(jetpt_sd->at(0));
      if(jetm_sd->size()>0)
	leadingSoftDropJetm.back()->Fill(jetm_sd->at(0));
      if(jetmraw_sd->size()>0)
	leadingSoftDropJetmraw.back()->Fill(jetmraw_sd->at(0));

      TLorentzVector ak8, ak8_gen, vBoson;
      if(jetpt->size() > 0 )
	ak8.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
      if(jeteta_gen->size() > 0)
	ak8_gen.SetPtEtaPhiM(jetpt_gen->at(0),jeteta_gen->at(0),jetphi_gen->at(0),jetm_gen->at(0));
      vBoson.SetPtEtaPhiM(*wzpt_h,*wzeta_h,*wzphi_h,*wzmass_h);
      
      dRGenJet.back()->Fill(ak8.DeltaR(ak8_gen));
      dRVboson.back()->Fill(ak8.DeltaR(vBoson));
      
      if(jetpt_gen->size() > 0 and jetpt->size() > 0)
	responseLeadJetPt.back()->Fill(jetpt->at(0)/jetpt_gen->at(0));
      if(jetm_pr->size() > 0 and jetm_pr_gen->size() > 0)
	responseLeadPrunedJetMass.back()->Fill(jetm_pr->at(0)/jetm_pr_gen->at(0));
      if(jetmraw_pr->size() > 0 and jetm_pr_gen->size() > 0)
	responseLeadPrunedJetMassRaw.back()->Fill(jetmraw_pr->at(0)/jetm_pr_gen->at(0));
      if(jetm_sd->size() > 0 and jetm_sd_gen -> size() > 0)
	responseLeadSoftDropJetMass.back()->Fill(jetm_sd->at(0)/jetm_sd_gen->at(0));
      if(jetmraw_sd->size() > 0 and jetm_sd_gen -> size() > 0)
	responseLeadSoftDropJetMassRaw.back()->Fill(jetmraw_sd->at(0)/jetm_sd_gen->at(0));

      if(jetsub2pt_pr->size() > 0 and jetsub1pt_pr->size() > 0)
	subjetPtRatio.back()->Fill(jetsub2pt_pr->at(0)/jetsub1pt_pr->at(0));
      
      if(formula_pr->EvalInstance()){
	if(jetm_pr_gen->size() > 0 and jetm_pr->size() > 0)
	  responseLeadPrunedJetMass_pr.back()->Fill(jetm_pr->at(0)/jetm_pr_gen->at(0));
	if(jetmraw_pr->size() > 0 and jetm_pr_gen->size() > 0)
	  responseLeadPrunedJetMassRaw_pr.back()->Fill(jetmraw_pr->at(0)/jetm_pr_gen->at(0));
	if(jetsub2pt_pr->size() > 0 and jetsub1pt_pr->size() > 0)
	  subjetPtRatio_pr.back()->Fill(jetsub2pt_pr->at(0)/jetsub1pt_pr->at(0));
	if(jettau1->size() > 0 and jettau2->size() > 0)
	  leadingJetTau2Tau1_pr.back()->Fill(jettau2->at(0)/jettau1->at(0));
      }
     
      if(formula_sb->EvalInstance()){
	if(jettau2->size() > 0 and jettau1->size() > 0)
	  leadingJetTau2Tau1_sb.back()->Fill(jettau2->at(0)/jettau1->at(0));
      }
    }
  }

  // background 
  if(not skipBackground){

    cout<<"Backgroud analysis"<<endl;
    TH1F* leadingJetPt_bkg       = new TH1F("leadingJetPt_bkg","Z #rightarrow #nu #nu",20,250,650);
    TH1F* leadingJetEta_bkg      = new TH1F("leadingJetEta_bkg","Z #rightarrow #nu #nu",20,-2.5,2.5);
    TH1F* leadingJetTau2tau1_bkg = new TH1F("leadingJetTau2tau1_bkg","Z #rightarrow #nu #nu",20,0,1);
    TH1F* leadingPrunedJetpt_bkg = new TH1F("leadingPrunedJetpt_bkg","Z #rightarrow #nu #nu",20,250,650);
    TH1F* leadingPrunedJetm_bkg  = new TH1F("leadingPrunedJetm_bkg","Z #rightarrow #nu #nu",30,0,120);
    TH1F* leadingPrunedJetmraw_bkg  = new TH1F("leadingPrunedJetmraw_bkg","Z #rightarrow #nu #nu",30,0,120);
    TH1F* leadingSoftDropJetpt_bkg = new TH1F("leadingSoftDropJetpt_bkg","Z #rightarrow #nu #nu",20,250,650);
    TH1F* leadingSoftDropJetm_bkg  = new TH1F("leadingSoftDropJetm_bkg","Z #rightarrow #nu #nu",30,0,120);
    TH1F* leadingSoftDropJetmraw_bkg  = new TH1F("leadingSoftDropJetmraw_bkg","Z #rightarrow #nu #nu",30,0,120);
    TH1F* boostedNJet_bkg = new TH1F("boostedNJet_bkg","Z #rightarrow #nu #nu",5,0,5);
    TH1F* responseLeadJetPt_bkg = new TH1F("responseLeadJetPt_bkg","Z #rightarrow #nu #nu",30,0.5,1.5);
    TH1F* responseLeadPrunedJetMass_bkg = new TH1F("responseLeadPrunedJetMass_bkg","Z #rightarrow #nu #nu",30,0.5,1.5);
    TH1F* responseLeadPrunedJetMassRaw_bkg = new TH1F("responseLeadPrunedJetMassRaw_bkg","Z #rightarrow #nu #nu",30,0.5,1.5);
    TH1F* responseLeadSoftDropJetMass_bkg = new TH1F("responseLeadSoftDropJetMass_bkg","Z #rightarrow #nu #nu",30,0.5,1.5);
    TH1F* responseLeadSoftDropJetMassRaw_bkg = new TH1F("responseLeadSoftDropJetMassRaw_bkg","Z #rightarrow #nu #nu",30,0.5,1.5);
    TH1F* subjetPtRatio_bkg = new TH1F("subjetPtRatio_bkg","Z #rightarrow #nu #nu",20,0,1);
    TH1F* responseLeadPrunedJetMass_pr_bkg = new TH1F("responseLeadPrunedJetMass_pr_bkg","Z #rightarrow #nu #nu",30,0.5,1.5);
    TH1F* responseLeadPrunedJetMassRaw_pr_bkg = new TH1F("responseLeadPrunedJetMassRaw_pr_bkg","Z #rightarrow #nu #nu",30,0.5,1.5);
    TH1F* subjetPtRatio_pr_bkg = new TH1F("subjetPtRatio_pr_bkg","Z #rightarrow #nu #nu",20,0,1);
    TH1F* leadingJetTau2Tau1_pr_bkg = new TH1F("leadingJetTau2Tau1_pr_bkg","Z #rightarrow #nu #nu",20,0,1);
    TH1F* leadingJetTau2Tau1_sb_bkg = new TH1F("leadingJetTau2Tau1_sb_bkg","Z #rightarrow #nu #nu",20,0,1);
    
    TTreeReader myReader(backgroundZnunu);

    TTreeReaderValue<UChar_t> hltm (myReader,"hltmet90");
    TTreeReaderValue<unsigned int> nbjets (myReader,"nbjetslowpt");
    TTreeReaderValue<unsigned int> njets (myReader,"njets");
    TTreeReaderValue<double> t1pfmet (myReader,"t1pfmet");
    TTreeReaderValue<double> jmdphi (myReader,"incjetmumetdphimin4");    
    TTreeReaderValue<vector<double> > cjetCH    (myReader,"centraljetCHfrac");
    TTreeReaderValue<vector<double> > cjetNH    (myReader,"centraljetNHfrac");
    TTreeReaderValue<vector<double> > cjetpt    (myReader,"centraljetpt");
    
    TTreeReaderValue<vector<double> > jetpt    (myReader,"boostedJetpt");
    TTreeReaderValue<vector<double> > jeteta   (myReader,"boostedJeteta");
    TTreeReaderValue<vector<double> > jetphi   (myReader,"boostedJetphi");
    TTreeReaderValue<vector<double> > jetm     (myReader,"boostedJetm");
    TTreeReaderValue<vector<double> > jetpt_gen    (myReader,"boostedJetGenpt");
    TTreeReaderValue<vector<double> > jeteta_gen   (myReader,"boostedJetGeneta");
    TTreeReaderValue<vector<double> > jetphi_gen   (myReader,"boostedJetGenphi");
    TTreeReaderValue<vector<double> > jetm_gen   (myReader,"boostedJetGenm");
    
    TTreeReaderValue<vector<double> > jettau2  (myReader,"boostedJettau2");
    TTreeReaderValue<vector<double> > jettau1  (myReader,"boostedJettau1");
    
    TTreeReaderValue<vector<double> > jetpt_pr (myReader,"prunedJetpt");
    TTreeReaderValue<vector<double> > jetm_pr  (myReader,"prunedJetm");
    TTreeReaderValue<vector<double> > jetmraw_pr  (myReader,"prunedJetmraw");
    TTreeReaderValue<vector<double> > jetsub1pt_pr  (myReader,"prunedSubJetpt_1");
    TTreeReaderValue<vector<double> > jetsub2pt_pr  (myReader,"prunedSubJetpt_2");
    
    TTreeReaderValue<vector<double> > jetpt_pr_gen (myReader,"prunedJetGenpt");
    TTreeReaderValue<vector<double> > jetm_pr_gen  (myReader,"prunedJetGenm");
    
    TTreeReaderValue<vector<double> > jetpt_sd (myReader,"softDropJetpt");
    TTreeReaderValue<vector<double> > jetm_sd  (myReader,"softDropJetm");
    TTreeReaderValue<vector<double> > jetpt_sd_gen (myReader,"softDropJetGenpt");
    TTreeReaderValue<vector<double> > jetm_sd_gen  (myReader,"softDropJetGenm");
    TTreeReaderValue<vector<double> > jetmraw_sd  (myReader,"softDropJetmraw");
    
    int iEvent = 0;
    cout<<"starting loop on events "<<backgroundZnunu->GetEntries()<<endl;
    while(myReader.Next()){
      float Njet = 0;
      iEvent++;    

      if(*hltm == 0) continue;
      if(*nbjets > 0) continue;
      if(*njets < 1) continue;
      if(*t1pfmet < 250) continue;
      if(*jmdphi < 0.5) continue;
      if(cjetCH->size() == 0 or cjetpt->size() == 0 or cjetNH->size() == 0) continue;
      if(cjetCH->at(0) < 0.1) continue;
      if(cjetNH->at(0) > 0.8) continue;
      if(cjetpt->at(0) < 100) continue;
      
      // AK8 jets with pT > 250 GeV
      for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
	if(jetpt->at(ijet) > 250)
	  Njet++;
      }
      
      // Fill histos
      boostedNJet_bkg->Fill(Njet);
      if(jetpt->size() > 0)
	leadingJetPt_bkg->Fill(jetpt->at(0));
      if(jeteta->size() > 0)
	leadingJetEta_bkg->Fill(jeteta->at(0));
      if(jettau2->size() >  0 and jettau1->size() > 0)
	leadingJetTau2tau1_bkg->Fill(jettau2->at(0)/jettau1->at(0));
      if(jetpt_pr->size() > 0)
	leadingPrunedJetpt_bkg->Fill(jetpt_pr->at(0));
      if(jetm_pr->size() > 0)
	leadingPrunedJetm_bkg->Fill(jetm_pr->at(0));
      if(jetmraw_pr->size() > 0)
	leadingPrunedJetmraw_bkg->Fill(jetmraw_pr->at(0));
      if(jetpt_sd->size() > 0)
	leadingSoftDropJetpt_bkg->Fill(jetpt_sd->at(0));
      if(jetm_sd->size() > 0)
	leadingSoftDropJetm_bkg->Fill(jetm_sd->at(0));
      if(jetmraw_sd->size() > 0)
	leadingSoftDropJetmraw_bkg->Fill(jetmraw_sd->at(0));
      
      if(jetpt->size() > 0 and jetpt_gen->size() > 0)
	responseLeadJetPt_bkg->Fill(jetpt->at(0)/jetpt_gen->at(0));
      if(jetm_pr->size() > 0 and jetm_pr_gen->size()> 0)
	responseLeadPrunedJetMass_bkg->Fill(jetm_pr->at(0)/jetm_pr_gen->at(0));
      if(jetmraw_pr->size() > 0 and jetm_pr_gen->size() > 0)
	responseLeadPrunedJetMassRaw_bkg->Fill(jetmraw_pr->at(0)/jetm_pr_gen->at(0));
      if(jetm_sd->size() > 0 and jetm_sd_gen->size() > 0)
	responseLeadSoftDropJetMass_bkg->Fill(jetm_sd->at(0)/jetm_sd_gen->at(0));
      if(jetmraw_sd->size() > 0 and jetm_sd_gen->size() > 0)
	responseLeadSoftDropJetMassRaw_bkg->Fill(jetmraw_sd->at(0)/jetm_sd_gen->at(0));
      if(jetsub2pt_pr->size() > 0 and jetsub1pt_pr->size() > 0)
	subjetPtRatio_bkg->Fill(jetsub2pt_pr->at(0)/jetsub1pt_pr->at(0));
      
      if(jetm_pr->size() > 0 and (jetm_pr->at(0) > 65 and jetm_pr->at(0) < 105)){
	if(jetm_pr->size() > 0 and jetm_pr_gen->size() > 0)
	  responseLeadPrunedJetMass_pr_bkg->Fill(jetm_pr->at(0)/jetm_pr_gen->at(0));
	if(jetmraw_pr->size() > 0 and jetm_pr_gen->size() > 0)
	  responseLeadPrunedJetMassRaw_pr_bkg->Fill(jetmraw_pr->at(0)/jetm_pr_gen->at(0));
	if(jetsub2pt_pr->size() > 0 and jetsub1pt_pr->size() > 0)
	  subjetPtRatio_pr_bkg->Fill(jetsub2pt_pr->at(0)/jetsub1pt_pr->at(0));
	if(jettau2->size() > 0 and jettau1->size() > 0)
	  leadingJetTau2Tau1_pr_bkg->Fill(jettau2->at(0)/jettau1->at(0));
      }
      else if(jetm_pr->size() > 0 and (jetm_pr->at(0) > 40 and jetm_pr->at(0) < 65)){
	leadingJetTau2Tau1_sb_bkg->Fill(jettau2->at(0)/jettau1->at(0));
      }
      else 
	continue;      
    }
    
    // make single plots
    makeShapePlots(cCanvas,leadingJetPt,leadingJetPt_bkg,
		   legend,outputDirectory,"jetPt","p_{T} (GeV)");
    makeShapePlots(cCanvas,leadingJetEta,leadingJetEta_bkg,
		   legend,outputDirectory,"jetEta","#eta");
    makeShapePlots(cCanvas,leadingJetTau2tau1,leadingJetTau2tau1_bkg,
		   legend,outputDirectory,"jetTau2tau1","#tau_{2}/#tau_{1}");
    makeShapePlots(cCanvas,leadingPrunedJetpt,leadingPrunedJetpt_bkg,
		   legend,outputDirectory,"jetPrunedPt","p_{T}^{pr} (GeV)");
    makeShapePlots(cCanvas,leadingPrunedJetm,leadingPrunedJetm_bkg,
		   legend,outputDirectory,"jetPrunedM","m_{pr} (GeV)");
    makeShapePlots(cCanvas,leadingPrunedJetmraw,leadingPrunedJetmraw_bkg,
		   legend,outputDirectory,"jetPrunedMraw","m_{pr} raw (GeV)");
    makeShapePlots(cCanvas,leadingSoftDropJetpt,leadingSoftDropJetpt_bkg,
		   legend,outputDirectory,"jetSoftDropPt","p_{T}^{sd} (GeV)");
    makeShapePlots(cCanvas,leadingSoftDropJetm,leadingSoftDropJetm_bkg,
		   legend,outputDirectory,"jetSoftDropM","m_{sd} (GeV)");
    makeShapePlots(cCanvas,leadingSoftDropJetmraw,leadingSoftDropJetmraw_bkg,
		   legend,outputDirectory,"jetSoftDropMraw","m_{sd} raw (GeV)");
    makeShapePlots(cCanvas,boostedNJet,boostedNJet_bkg,
		   legend,outputDirectory,"NJet","N_{jet}^{AK8} (GeV)");
    makeShapePlots(cCanvas,dRGenJet,
		   legend,outputDirectory,"dRGenJet","#DeltaR(reco,gen)");
    makeShapePlots(cCanvas,dRVboson,
		   legend,outputDirectory,"dRVboson","#DeltaR(reco,V)");
    makeShapePlots(cCanvas,subjetPtRatio,subjetPtRatio_bkg,
		 legend,outputDirectory,"SubJetPtRatio","#Pt ratio");
    
    makeShapePlots(cCanvas,responseLeadJetPt,responseLeadJetPt_bkg,
		   legend,outputDirectory,"ResponsePt","p_{T}^{reco}-p_{T}^{gen} (GeV)");
    makeShapePlots(cCanvas,responseLeadPrunedJetMass,responseLeadPrunedJetMass_bkg,
		   legend,outputDirectory,"ResponsePrunedJetMass","m_{pr}^{reco}-m_{pr}^{gen} (GeV)");
    makeShapePlots(cCanvas,responseLeadPrunedJetMassRaw,responseLeadPrunedJetMassRaw_bkg,
		   legend,outputDirectory,"ResponsePrunedJetMassRaw","m_{pr}^{reco}-m_{pr}^{gen} raw (GeV)");
    makeShapePlots(cCanvas,responseLeadSoftDropJetMass,responseLeadSoftDropJetMass_bkg,
		 legend,outputDirectory,"ResponseSoftDropJetMass","m_{sd}^{reco}-m_{sd}^{gen} (GeV)");
    makeShapePlots(cCanvas,responseLeadSoftDropJetMassRaw,responseLeadSoftDropJetMassRaw_bkg,
		   legend,outputDirectory,"ResponseSoftDropJetMassRaw","m_{sd}^{reco}-m_{sd}^{gen} raw (GeV)");
    
    makeShapePlots(cCanvas,responseLeadPrunedJetMass_pr,responseLeadPrunedJetMass_pr_bkg,
		   legend,outputDirectory,"ResponsePrunedJetMass_pr","m_{pr}^{reco}-m_{pr}^{gen} (GeV)");
    makeShapePlots(cCanvas,responseLeadPrunedJetMassRaw_pr,responseLeadPrunedJetMassRaw_pr_bkg,
		   legend,outputDirectory,"ResponsePrunedJetMassRaw_pr","m_{pr}^{reco}-m_{pr}^{gen} raw (GeV)");
    makeShapePlots(cCanvas,subjetPtRatio_pr,subjetPtRatio_pr_bkg,
		   legend,outputDirectory,"SubJetPtRatio_pr","p_{T} ratio");
    makeShapePlots(cCanvas,leadingJetTau2Tau1_pr,leadingJetTau2Tau1_pr_bkg,
		   legend,outputDirectory,"jetTau2tau1_pr","#tau_{2}/#tau_{1}");
    makeShapePlots(cCanvas,leadingJetTau2Tau1_sb,leadingJetTau2Tau1_sb_bkg,
		   legend,outputDirectory,"jetTau2tau1_sb","#tau_{2}/#tau_{1}");
  }
  else{

    makeShapePlots(cCanvas,leadingJetPt,
		   legend,outputDirectory,"jetPt","p_{T} (GeV)");
    makeShapePlots(cCanvas,leadingJetEta,
		   legend,outputDirectory,"jetEta","#eta");
    makeShapePlots(cCanvas,leadingJetTau2tau1,
		   legend,outputDirectory,"jetTau2tau1","#tau_{2}/#tau_{1}");
    makeShapePlots(cCanvas,leadingPrunedJetpt,
		   legend,outputDirectory,"jetPrunedPt","p_{T}^{pr} (GeV)");
    makeShapePlots(cCanvas,leadingPrunedJetm,
		   legend,outputDirectory,"jetPrunedM","m_{pr} (GeV)");
    makeShapePlots(cCanvas,leadingPrunedJetmraw,
		   legend,outputDirectory,"jetPrunedMraw","m_{pr} raw (GeV)");
    makeShapePlots(cCanvas,leadingSoftDropJetpt,
		   legend,outputDirectory,"jetSoftDropPt","p_{T}^{sd} (GeV)");
    makeShapePlots(cCanvas,leadingSoftDropJetm,
		   legend,outputDirectory,"jetSoftDropM","m_{sd} (GeV)");
    makeShapePlots(cCanvas,leadingSoftDropJetmraw,
		   legend,outputDirectory,"jetSoftDropMraw","m_{sd} raw (GeV)");
    makeShapePlots(cCanvas,boostedNJet,
		   legend,outputDirectory,"NJet","N_{jet}^{AK8} (GeV)");
    makeShapePlots(cCanvas,dRGenJet,
		   legend,outputDirectory,"dRGenJet","#DeltaR(reco,gen)");
    makeShapePlots(cCanvas,dRVboson,
		   legend,outputDirectory,"dRVboson","#DeltaR(reco,V)");
    makeShapePlots(cCanvas,subjetPtRatio,
		 legend,outputDirectory,"SubJetPtRatio","#Pt ratio");
    
    makeShapePlots(cCanvas,responseLeadJetPt,
		   legend,outputDirectory,"ResponsePt","p_{T}^{reco}-p_{T}^{gen} (GeV)");
    makeShapePlots(cCanvas,responseLeadPrunedJetMass,
		   legend,outputDirectory,"ResponsePrunedJetMass","m_{pr}^{reco}-m_{pr}^{gen} (GeV)");
    makeShapePlots(cCanvas,responseLeadPrunedJetMassRaw,
		   legend,outputDirectory,"ResponsePrunedJetMassRaw","m_{pr}^{reco}-m_{pr}^{gen} raw (GeV)");
    makeShapePlots(cCanvas,responseLeadSoftDropJetMass,
		 legend,outputDirectory,"ResponseSoftDropJetMass","m_{sd}^{reco}-m_{sd}^{gen} (GeV)");
    makeShapePlots(cCanvas,responseLeadSoftDropJetMassRaw,
		   legend,outputDirectory,"ResponseSoftDropJetMassRaw","m_{sd}^{reco}-m_{sd}^{gen} raw (GeV)");
    
    makeShapePlots(cCanvas,responseLeadPrunedJetMass_pr,
		   legend,outputDirectory,"ResponsePrunedJetMass_pr","m_{pr}^{reco}-m_{pr}^{gen} (GeV)");
    makeShapePlots(cCanvas,responseLeadPrunedJetMassRaw_pr,
		   legend,outputDirectory,"ResponsePrunedJetMassRaw_pr","m_{pr}^{reco}-m_{pr}^{gen} raw (GeV)");
    makeShapePlots(cCanvas,subjetPtRatio_pr,
		   legend,outputDirectory,"SubJetPtRatio_pr","p_{T} ratio");
    makeShapePlots(cCanvas,leadingJetTau2Tau1_pr,
		   legend,outputDirectory,"jetTau2tau1_pr","#tau_{2}/#tau_{1}");
    makeShapePlots(cCanvas,leadingJetTau2Tau1_sb,
		   legend,outputDirectory,"jetTau2tau1_sb","#tau_{2}/#tau_{1}");


  }
  return;
}

void makeShapePlots(TCanvas* cCanvas, vector<TH1F*> histoList, TLegend* legend, string outputDirectory, string plotName, string xAxisTitle){

  legend->Clear();

  TLatex * tex = new TLatex(0.94,0.92," 13 TeV");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  TLatex * tex2 = new TLatex(0.14,0.92,"CMS");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.04);
  tex2->SetLineWidth(2);
  TLatex * tex3 = new TLatex(0.286,0.92,"Preliminary Simulation");
  tex3->SetNDC();
  tex3->SetTextFont(52);
  tex3->SetTextSize(0.035);
  tex3->SetLineWidth(2);

  float maxYScale = 0;
  for(size_t iHisto = 0; iHisto < histoList.size(); iHisto++){

    float integral = histoList.at(iHisto)->Integral();
    histoList.at(iHisto)->Scale(1./integral);

    if(iHisto == 0){
      maxYScale = histoList.at(iHisto)->GetMaximum();
    }
    else{
      if(histoList.at(iHisto)->GetMaximum() > maxYScale)
	maxYScale = histoList.at(iHisto)->GetMaximum();
    }
  }    
  
  for(size_t iHisto = 0; iHisto < histoList.size(); iHisto++){

    histoList.at(iHisto)->SetLineWidth(2);
    histoList.at(iHisto)->SetLineColor(iHisto+1);

    if(iHisto%2==0)
      histoList.at(iHisto)->SetLineStyle(1);
    else
      histoList.at(iHisto)->SetLineStyle(2);

    histoList.at(iHisto)->SetMarkerStyle(20);
    histoList.at(iHisto)->SetMarkerSize(1.5);

    legend->AddEntry(histoList.at(iHisto),histoList.at(iHisto)->GetTitle(),"l");

    if(iHisto == 0){
      histoList.at(iHisto)->GetXaxis()->SetTitle(xAxisTitle.c_str());
      histoList.at(iHisto)->GetYaxis()->SetRangeUser(0,maxYScale*1.5);
      histoList.at(iHisto)->Draw("hist");
    }
    else
      histoList.at(iHisto)->Draw("histsame");

  }


  legend->Draw("same");
  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");

  cCanvas->SaveAs((outputDirectory+"/shapes/"+plotName+".png").c_str(),"png");
  cCanvas->SaveAs((outputDirectory+"/shapes/"+plotName+".pdf").c_str(),"pdf");
};



void makeShapePlots(TCanvas* cCanvas, vector<TH1F*> histoList, TH1F* histoBkg, TLegend* legend, string outputDirectory, string plotName, string xAxisTitle){

  legend->Clear();

  TLatex * tex = new TLatex(0.94,0.92," 13 TeV");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  TLatex * tex2 = new TLatex(0.14,0.92,"CMS");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.04);
  tex2->SetLineWidth(2);
  TLatex * tex3 = new TLatex(0.286,0.92,"Preliminary Simulation");
  tex3->SetNDC();
  tex3->SetTextFont(52);
  tex3->SetTextSize(0.035);
  tex3->SetLineWidth(2);

  float maxYScale = 0;

  float integral = histoBkg->Integral();
  histoBkg->Scale(1./integral);
  maxYScale = histoBkg->GetMaximum();

  for(size_t iHisto = 0; iHisto < histoList.size(); iHisto++){

    float integral = histoList.at(iHisto)->Integral();
    histoList.at(iHisto)->Scale(1./integral);

    if(histoList.at(iHisto)->GetMaximum() > maxYScale)
      maxYScale = histoList.at(iHisto)->GetMaximum();
  }

  histoBkg->SetLineColor(kBlack);
  histoBkg->SetLineWidth(2);
  histoBkg->SetFillColor(kGray);
  histoBkg->SetFillStyle(3001);
  histoBkg->GetXaxis()->SetTitle(xAxisTitle.c_str());
  histoBkg->GetYaxis()->SetRangeUser(0,maxYScale*1.5);
  histoBkg->Draw("hist");

  legend->AddEntry(histoBkg,histoBkg->GetTitle(),"FL");

  for(size_t iHisto = 0; iHisto < histoList.size(); iHisto++){
    
    histoList.at(iHisto)->SetLineWidth(2);
    histoList.at(iHisto)->SetLineColor(iHisto+1);

    if(iHisto%2==0)
      histoList.at(iHisto)->SetLineStyle(1);
    else
      histoList.at(iHisto)->SetLineStyle(2);

    histoList.at(iHisto)->SetMarkerStyle(20);
    histoList.at(iHisto)->SetMarkerSize(1.5);

    legend->AddEntry(histoList.at(iHisto),histoList.at(iHisto)->GetTitle(),"l");

    histoList.at(iHisto)->Draw("histsame");
    
  }
  

  legend->Draw("same");
  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");

  cCanvas->SaveAs((outputDirectory+"/shapes/"+plotName+"_SandB.png").c_str(),"png");
  cCanvas->SaveAs((outputDirectory+"/shapes/"+plotName+"_SandB.pdf").c_str(),"pdf");
};



void makeEfficiencyPlots(TCanvas* cCanvas, unordered_map<string,TH1F*> eff_sig, vector<string> Mphi, vector<string> Mchi,string bosonMode, string interactionModel, 
			 TLegend* legend, string outputDirectory, string plotName, bool isRelative){

  legend->Clear();

  TLatex * tex = new TLatex(0.94,0.92," 13 TeV");
  tex->SetNDC();
  tex->SetTextAlign(31);
  tex->SetTextFont(42);
  tex->SetTextSize(0.04);
  tex->SetLineWidth(2);
  TLatex * tex2 = new TLatex(0.14,0.92,"CMS");
  tex2->SetNDC();
  tex2->SetTextFont(61);
  tex2->SetTextSize(0.04);
  tex2->SetLineWidth(2);
  TLatex * tex3 = new TLatex(0.286,0.92,"Preliminary Simulation");
  tex3->SetNDC();
  tex3->SetTextFont(52);
  tex3->SetTextSize(0.035);
  tex3->SetLineWidth(2);

  for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
    
    eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetLineWidth(2);
    eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetLineColor(iMass+1);
    eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetMarkerColor(iMass+1);

    if(iMass %2 == 0)
      eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetLineStyle(1);
    else
      eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetLineStyle(2);

    eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetMarkerStyle(20);
    eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetMarkerSize(1.5);

    if(not isRelative)
      eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetYaxis()->SetRangeUser(1e-3,0.7);
    else
      eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetYaxis()->SetRangeUser(1e-3,1.0);

    if(iMass == 0){
      eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetYaxis()->SetTitle("efficiency");
      if(not isRelative)
	eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->Draw("hist");
      else
	eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->Draw("p");
    }
    else{
      if(not isRelative)
	eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->Draw("histsame");
      else
	eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->Draw("psame");
    }

    if(not isRelative)
      legend->AddEntry(eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]],(bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"l");
    else
      legend->AddEntry(eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]],(bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"p");
    
  }

  legend->Draw("same");
  tex->Draw("same");
  tex2->Draw("same");
  tex3->Draw("same");
  
  cCanvas->SaveAs((outputDirectory+"/efficiency/"+plotName+".png").c_str(),"png");
  cCanvas->SaveAs((outputDirectory+"/efficiency/"+plotName+".pdf").c_str(),"pdf");
  cCanvas->SetLogy(1);
  cCanvas->SaveAs((outputDirectory+"/efficiency/"+plotName+"_Log.png").c_str(),"png");
  cCanvas->SaveAs((outputDirectory+"/efficiency/"+plotName+"_Log.pdf").c_str(),"pdf");
  cCanvas->SetLogy(0);
}
