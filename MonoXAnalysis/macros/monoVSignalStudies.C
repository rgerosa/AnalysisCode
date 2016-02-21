//monoV_signalAnalysis("/store/caf/user/rgerosa/MONOJET_ANALYSIS/Production-14-1-2016/",1,"Pseudoscalar","Plots_MonoW_Pseudoscalar")
#include <iostream>
#include <vector>
#include <string>
#include <string>
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

using namespace std;

vector<string> Mphi_Vector = {"100","200","500","1000"};
vector<string> Mphi_Axial  = {"100","200","500","1000"};
vector<string> Mphi_PseudoScalar  = {"100","200","500","1000"};
vector<string> Mphi_Scalar = {"100","125","1000","2000"};

vector<string> Mchi_Vector = {"50","10","50"};
vector<string> Mchi_Axial = {"50","10","50"};
vector<string> Mchi_PseudoScalar = {"50","10","50"};
vector<string> Mchi_Scalar = {"50","10","10","10"};

vector<pair<string,string> > selectionsSig;
vector<pair<string,string> > selectionsWmn;
vector<pair<string,string> > selectionsWen;
vector<pair<string,string> > selectionsZmm;
vector<pair<string,string> > selectionsZee;

void makeEfficiencyPlots(TCanvas* cCanvas, unordered_map<string,TH1F*> eff_sig, vector<string> Mphi, vector<string> Mchi, string bosonMode, string interactionModel, 
			 TLegend* legend, string outputDirectory, string plotName, bool isRelative);

void makeShapePlots(TCanvas* cCanvas, vector<TH1F*> histoList, TLegend* legend, string outputDirectory, string plotName, string xAxisTile);
void makeShapePlots(TCanvas* cCanvas, vector<TH1F*> histoList, TH1F* histBkg, TLegend* legend, string outputDirectory, string plotName, string xAxisTile);
  		       		 
void monoV_signalAnalysis(string baseInputPath, bool isW, string interactionModel, string generator, string outputDirectory, bool makeEfficiency = false){

  gROOT->SetBatch(1);

  gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadTopMargin(0.09);
  gStyle->SetErrorX(0.5);
  gStyle->SetOptTitle(0);

  // define signal region cuts
  selectionsSig.push_back(make_pair("trigger","hltmet90"));
  selectionsSig.push_back(make_pair("njets","hltmet90 && njets >= 1"));
  selectionsSig.push_back(make_pair("jetpt","hltmet90 && njets >= 1 && centraljetpt[0] > 100"));
  selectionsSig.push_back(make_pair("bveto","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0"));
  selectionsSig.push_back(make_pair("jetid","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8"));
  selectionsSig.push_back(make_pair("dphijemet","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5"));
  selectionsSig.push_back(make_pair("met","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1pfmet > 200"));
  selectionsSig.push_back(make_pair("ak8pt","hltmet90 && njets >= 1 && centraljetpt[0] > 100 && nbjetslowpt == 0 && centraljetCHfrac[0] > 0.1 && centraljetNHfrac[0] < 0.8 && incjetmumetdphimin4 > 0.5 && t1pfmet > 200 && boostedJetpt[0] > 200"));

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

  // take input files
  string bosonMode;
  if(isW)
    bosonMode = "MonoW";
  else
    bosonMode = "MonoZ";

  // start taking input files before the skim
  vector<TFile*> fileMonoV;  
  vector<TTree*> treeMonoV;  

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

  system(("mkdir "+outputDirectory).c_str());
  system(("mkdir "+outputDirectory+"/efficiency/").c_str());
  system(("mkdir "+outputDirectory+"/shapes/").c_str());
 

  if(makeEfficiency){
    for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
      fileMonoV.push_back(TFile::Open(("root://eoscms.cern.ch//"+baseInputPath+"/"+bosonMode+"_"+interactionModel+"/sigfilter/tree_"+interactionModel+bosonMode+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_gSM-1p0_gDM-1p0_13TeV-"+generator+".root").c_str()));
      if(fileMonoV.back() and not fileMonoV.back()->IsZombie()){
	treeMonoV.push_back((TTree*) fileMonoV.back()->Get("tree/tree"));
      }
    }
  }
  // start taking input files in the signal region
  // "nmuons == 0 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && hltmet90 > 0"
  vector<TFile*> fileMonoV_Sig;  
  vector<TTree*> treeMonoV_Sig;  
    
  for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
    fileMonoV_Sig.push_back(TFile::Open(("root://eoscms.cern.ch//"+baseInputPath+"/"+bosonMode+"_"+interactionModel+"/sigfilter/sig_tree_"+interactionModel+bosonMode+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_gSM-1p0_gDM-1p0_13TeV-"+generator+".root").c_str()));
    if(fileMonoV_Sig.back() and not fileMonoV_Sig.back()->IsZombie())
      treeMonoV_Sig.push_back((TTree*) fileMonoV_Sig.back()->Get("tree/tree"));
  }



  //Background Znunu
  TChain* backgroundZnunu = new TChain("tree/tree");
  TFile* file1 = TFile::Open(("root://eoscms.cern.ch//"+baseInputPath+"/ZJets/sigfilter/sig_tree_ZJetsToNuNu_HT-100To200_13TeV-madgraph.root").c_str());
  TFile* file2 = TFile::Open(("root://eoscms.cern.ch//"+baseInputPath+"/ZJets/sigfilter/sig_tree_ZJetsToNuNu_HT-200To400_13TeV-madgraph.root").c_str());
  TFile* file3 = TFile::Open(("root://eoscms.cern.ch//"+baseInputPath+"/ZJets/sigfilter/sig_tree_ZJetsToNuNu_HT-400To600_13TeV-madgraph.root").c_str());
  TFile* file4 = TFile::Open(("root://eoscms.cern.ch//"+baseInputPath+"/ZJets/sigfilter/sig_tree_ZJetsToNuNu_HT-600ToInf_13TeV-madgraph.root").c_str());
  backgroundZnunu->Add(file1->GetName());
  backgroundZnunu->Add(file2->GetName());
  backgroundZnunu->Add(file3->GetName());
  backgroundZnunu->Add(file4->GetName());

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
      
  // start taking input files in the signal region  
  vector<TFile*> fileMonoV_CRm;  
  vector<TTree*> treeMonoV_CRm;  
  vector<TFile*> fileMonoV_CRe;  
  vector<TTree*> treeMonoV_CRe;  

  if(makeEfficiency){

      if(isW){
	
	//"nmuons == 1 && nelectrons == 0 && ntaus == 0 && nphotons == 0 && mu1pt > 20 && mu1id >= 1"
	for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
	  fileMonoV_CRm.push_back(TFile::Open(("root://eoscms.cern.ch//"+baseInputPath+"/"+bosonMode+"_"+interactionModel+"/wmnfilter/wmn_tree_"+interactionModel+bosonMode+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_gSM-1p0_gDM-1p0_13TeV-"+generator+".root").c_str()));
	  if(fileMonoV_CRm.back() and not fileMonoV_CRm.back()->IsZombie())
	    treeMonoV_CRm.push_back((TTree*) fileMonoV_CRm.back()->Get("tree/tree"));      
	}
	
	for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
	  fileMonoV_CRe.push_back(TFile::Open(("root://eoscms.cern.ch//"+baseInputPath+"/"+bosonMode+"_"+interactionModel+"/wenfilter/wen_tree_"+interactionModel+bosonMode+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_gSM-1p0_gDM-1p0_13TeV-"+generator+".root").c_str()));
	  if(fileMonoV_CRe.back() and not fileMonoV_CRe.back()->IsZombie())
	    treeMonoV_CRe.push_back((TTree*) fileMonoV_CRe.back()->Get("tree/tree"));
    }
	
      }
      else{
	
	//"nmuons == 0 && nelectrons == 1 && ntaus == 0 && nphotons == 0 && el1pt > 40 && el1id >= 1"
	for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
	  fileMonoV_CRm.push_back(TFile::Open(("root://eoscms.cern.ch//"+baseInputPath+"/"+bosonMode+"_"+interactionModel+"/zmmfilter/zmm_tree_"+interactionModel+bosonMode+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_gSM-1p0_gDM-1p0_13TeV-"+generator+".root").c_str()));
	  if(fileMonoV_CRm.back() and not fileMonoV_CRm.back()->IsZombie())
	    treeMonoV_CRm.push_back((TTree*) fileMonoV_CRm.back()->Get("tree/tree"));
	}
	
	for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
	  fileMonoV_CRe.push_back(TFile::Open(("root://eoscms.cern.ch//"+baseInputPath+"/"+bosonMode+"_"+interactionModel+"/zeefilter/zee_tree_"+interactionModel+bosonMode+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_gSM-1p0_gDM-1p0_13TeV-"+generator+".root").c_str()));
	  if(fileMonoV_CRe.back() and not fileMonoV_CRe.back()->IsZombie())
	    treeMonoV_CRe.push_back((TTree*) fileMonoV_CRe.back()->Get("tree/tree"));
	}
      }

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
	eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] = new TH1F(("eff_sig_"+bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"",selectionsSig.size()+1,0,selectionsSig.size()+1);
	
	eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> SetBinContent(1,float(treeMonoV_Sig.at(iMass)->GetEntries())/treeMonoV.at(iMass)->GetEntries());
	eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> GetXaxis()->SetBinLabel(1,"Preselection");
	
	
	eff_sig_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] = new TH1F(("eff_sig_rel_"+bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),"",selectionsSig.size()+1,0,selectionsSig.size()+1);
	
	eff_sig_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> SetBinContent(1,float(treeMonoV_Sig.at(iMass)->GetEntries())/treeMonoV.at(iMass)->GetEntries());
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
	eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> SetBinContent(1,float(treeMonoV_CRm.at(iMass)->GetEntries())/treeMonoV.at(iMass)->GetEntries());
	eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> GetXaxis()->SetBinLabel(1,"Preselection");
	
	eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> SetBinContent(1,float(treeMonoV_CRe.at(iMass)->GetEntries())/treeMonoV.at(iMass)->GetEntries());
	eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> GetXaxis()->SetBinLabel(1,"Preselection");
	
	eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> SetBinContent(1,float(treeMonoV_CRm.at(iMass)->GetEntries())/treeMonoV.at(iMass)->GetEntries());
	eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> GetXaxis()->SetBinLabel(1,"Preselection");
	
	eff_CRe_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> SetBinContent(1,float(treeMonoV_CRe.at(iMass)->GetEntries())/treeMonoV.at(iMass)->GetEntries());
	eff_CRe_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]] -> GetXaxis()->SetBinLabel(1,"Preselection");
	
	// apply signal selections
	int iselection = 1;
	for(auto iSel = selectionsSig.begin(); iSel != selectionsSig.end(); iSel++){
	  iselection++;
	  eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_Sig.at(iMass)->GetEntries(iSel->second.c_str()))/treeMonoV.at(iMass)->GetEntries());
	  eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
	  
	  eff_sig_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_Sig.at(iMass)->GetEntries(iSel->second.c_str()))/(treeMonoV.at(iMass)->GetEntries())*1/eff_sig[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetBinContent(iselection-1));
	  eff_sig_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
	}
	
	if(isW){
	  iselection = 1;
	  for(auto iSel = selectionsWmn.begin(); iSel != selectionsWmn.end(); iSel++){
	    iselection ++;
	    eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRm.at(iMass)->GetEntries(iSel->second.c_str()))/treeMonoV.at(iMass)->GetEntries());
	    eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
	    
	    eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRm.at(iMass)->GetEntries(iSel->second.c_str()))/(treeMonoV.at(iMass)->GetEntries())*1/eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetBinContent(iselection-1));
	    eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
	  }
	  
	  iselection = 1;
	  for(auto iSel = selectionsWen.begin(); iSel != selectionsWen.end(); iSel++){
	    iselection ++;
	    eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRe.at(iMass)->GetEntries(iSel->second.c_str()))/treeMonoV.at(iMass)->GetEntries());
	    eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
	    
	    eff_CRe_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRe.at(iMass)->GetEntries(iSel->second.c_str()))/(treeMonoV.at(iMass)->GetEntries())*1/eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetBinContent(iselection-1));
	    eff_CRe_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
	  }      
	}
	else{
	  
	  iselection = 1;
	  for(auto iSel = selectionsZmm.begin(); iSel != selectionsZmm.end(); iSel++){
	    iselection ++;
	    eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRm.at(iMass)->GetEntries(iSel->second.c_str()))/treeMonoV.at(iMass)->GetEntries());
	    eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
	    
	    eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRm.at(iMass)->GetEntries(iSel->second.c_str()))/(treeMonoV.at(iMass)->GetEntries())*1/eff_CRm[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetBinContent(iselection-1));
	    eff_CRm_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());
	    
	  }
	  
	  iselection = 1;
	  for(auto iSel = selectionsZee.begin(); iSel != selectionsZee.end(); iSel++){
	    iselection ++;
	    eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRe.at(iMass)->GetEntries(iSel->second.c_str()))/treeMonoV.at(iMass)->GetEntries());
	    eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetXaxis()->SetBinLabel(iselection,iSel->first.c_str());

	    eff_CRe_relative[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->SetBinContent(iselection,float(treeMonoV_CRe.at(iMass)->GetEntries(iSel->second.c_str()))/(treeMonoV.at(iMass)->GetEntries())*1/eff_CRe[bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]]->GetBinContent(iselection-1));
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

  /// now look at substructure in the signal region only
  vector<TH1F*> leadingJetPt;
  vector<TH1F*> leadingJetEta;
  vector<TH1F*> leadingJetTau2tau1;
  vector<TH1F*> leadingPrunedJetpt;
  vector<TH1F*> leadingPrunedJetm;
  vector<TH1F*> leadingPrunedJetptraw;
  vector<TH1F*> leadingPrunedJetmraw;
  vector<TH1F*> leadingSoftDropJetpt;
  vector<TH1F*> leadingSoftDropJetm;
  vector<TH1F*> leadingSoftDropJetptraw;
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

  
  for(size_t iMass = 0; iMass < Mphi.size(); iMass++){
    string cut;
    for(size_t iSelection = 0; iSelection < selectionsSig.size(); iSelection++){
      if(selectionsSig.at(iSelection).first == "ak8pt")
	cut = selectionsSig.at(iSelection).second;
    }

    leadingJetPt.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingJetPt").c_str(),
				    (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,150,650));
    treeMonoV_Sig.at(iMass)->Draw(("boostedJetpt[0] >> "+string(leadingJetPt.back()->GetName())).c_str(),cut.c_str(),"goff");

    leadingJetEta.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingJetEta").c_str(),
				     (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),20,-2.5,2.5));
    treeMonoV_Sig.at(iMass)->Draw(("boostedJeteta[0] >> "+string(leadingJetEta.back()->GetName())).c_str(),cut.c_str(),"goff");

    leadingJetTau2tau1.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingJetTau2tau1").c_str(),
					  (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),20,0,1));
    treeMonoV_Sig.at(iMass)->Draw(("boostedJettau2[0]/boostedJettau1[0] >> "+string(leadingJetTau2tau1.back()->GetName())).c_str(),cut.c_str(),"goff");

    leadingPrunedJetpt.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingPrunedJetpt").c_str(),
					  (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,150,650));
    treeMonoV_Sig.at(iMass)->Draw(("prunedJetpt[0] >> "+string(leadingPrunedJetpt.back()->GetName())).c_str(),cut.c_str(),"goff");

    leadingPrunedJetm.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingPrunedJetm").c_str(),
					 (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0,120));
    treeMonoV_Sig.at(iMass)->Draw(("prunedJetm[0] >> "+string(leadingPrunedJetm.back()->GetName())).c_str(),cut.c_str(),"goff");

    leadingPrunedJetptraw.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingPrunedJetptraw").c_str(),
					     (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,150,650));
    treeMonoV_Sig.at(iMass)->Draw(("prunedJetptraw[0] >> "+string(leadingPrunedJetptraw.back()->GetName())).c_str(),cut.c_str(),"goff");

    leadingPrunedJetmraw.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingPrunedJetmraw").c_str(),
					    (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0,120));
    treeMonoV_Sig.at(iMass)->Draw(("prunedJetmraw[0] >> "+string(leadingPrunedJetmraw.back()->GetName())).c_str(),cut.c_str(),"goff");

    leadingSoftDropJetpt.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingSoftDropJetpt").c_str(),
					    (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,150,650));
    treeMonoV_Sig.at(iMass)->Draw(("prunedJetpt[0] >> "+string(leadingSoftDropJetpt.back()->GetName())).c_str(),cut.c_str(),"goff");

    leadingSoftDropJetm.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingSoftDropJetm").c_str(),
					   (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0,120));
    treeMonoV_Sig.at(iMass)->Draw(("prunedJetm[0] >> "+string(leadingSoftDropJetm.back()->GetName())).c_str(),cut.c_str(),"goff");

    leadingSoftDropJetptraw.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingSoftDropJetptraw").c_str(),
					       (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,150,650));
    treeMonoV_Sig.at(iMass)->Draw(("prunedJetptraw[0] >> "+string(leadingSoftDropJetptraw.back()->GetName())).c_str(),cut.c_str(),"goff");

    leadingSoftDropJetmraw.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingSoftDropJetmraw").c_str(),
					      (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0,120));
    treeMonoV_Sig.at(iMass)->Draw(("prunedJetmraw[0] >> "+string(leadingSoftDropJetmraw.back()->GetName())).c_str(),cut.c_str(),"goff");    

    boostedNJet.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_boostedNJet").c_str(),
				   (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),5,0,5));
    TTreeReader myReader(treeMonoV_Sig.at(iMass));
    TTreeReaderValue<vector<double> > jetpt (myReader,"boostedJetpt");
    TTreeFormula sel ("sel",cut.c_str(),treeMonoV_Sig.at(iMass));
    int iEvent = 0;
    while(myReader.Next()){
      float Njet = 0;
      treeMonoV_Sig.at(iMass)->GetEntry(iEvent);
      iEvent++;
      if(not sel.EvalInstance())
	continue;
      for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
	if(jetpt->at(ijet) > 200)
	  Njet++;
      }
      
      boostedNJet.back()->Fill(Njet);
    }

    dRGenJet.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_dRGenJet").c_str(),
				(bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),50,0,0.8));
    treeMonoV_Sig.at(iMass)->Draw(("sqrt((boostedJeteta[0]-boostedJetGeneta[0])*(boostedJeteta[0]-boostedJetGeneta[0]) + (boostedJetphi[0]-boostedJetGenphi[0])*(boostedJetphi[0]-boostedJetGenphi[0])) >> "+string(dRGenJet.back()->GetName())).c_str(),cut.c_str(),"goff");    

    dRVboson.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_dRVboson").c_str(),
				(bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),50,0,0.8));
    treeMonoV_Sig.at(iMass)->Draw(("sqrt((boostedJeteta[0]-boostedJetBosoneta[0])*(boostedJeteta[0]-boostedJetBosoneta[0]) + (boostedJetphi[0]-boostedJetBosonphi[0])*(boostedJetphi[0]-boostedJetBosonphi[0])) >> "+string(dRVboson.back()->GetName())).c_str(),cut.c_str(),"goff");    

    responseLeadJetPt.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_responseLeadJetPt").c_str(),
					 (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0.5,1.5));
    treeMonoV_Sig.at(iMass)->Draw(("boostedJetpt[0]/boostedJetGenpt[0] >> "+string(responseLeadJetPt.back()->GetName())).c_str(),cut.c_str(),"goff");    

    responseLeadPrunedJetMass.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_responseLeadPrunedJetMass").c_str(),
						 (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0.5,1.5));
    treeMonoV_Sig.at(iMass)->Draw(("prunedJetm[0]/prunedJetGenm[0] >> "+string(responseLeadPrunedJetMass.back()->GetName())).c_str(),cut.c_str(),"goff");    

    responseLeadPrunedJetMassRaw.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_responseLeadPrunedJetMassRaw").c_str(),
						    (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0.5,1.5));
    treeMonoV_Sig.at(iMass)->Draw(("(prunedJetmraw[0]/prunedJetGenm[0]) >> "+string(responseLeadPrunedJetMassRaw.back()->GetName())).c_str(),cut.c_str(),"goff");    


    responseLeadSoftDropJetMass.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_responseLeadSoftDropJetMass").c_str(),
						   (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0.5,1.5));
    treeMonoV_Sig.at(iMass)->Draw(("(softDropJetm[0]/softDropJetGenm[0]) >> "+string(responseLeadSoftDropJetMass.back()->GetName())).c_str(),cut.c_str(),"goff");    

    responseLeadSoftDropJetMassRaw.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_responseLeadSoftDropJetMassRaw").c_str(),
						      (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0.5,1.5));
    treeMonoV_Sig.at(iMass)->Draw(("(softDropJetmraw[0]/softDropJetGenm[0]) >> "+string(responseLeadSoftDropJetMassRaw.back()->GetName())).c_str(),cut.c_str(),"goff");    

    subjetPtRatio.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_subjetPtRatio").c_str(),
				     (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0,2));
    treeMonoV_Sig.at(iMass)->Draw(("(prunedSubJetpt_2[0]/prunedSubJetpt_1[0]) >> "+string(subjetPtRatio.back()->GetName())).c_str(),cut.c_str(),"goff");

    //
    string cut_pr = cut + " && prunedJetm[0] > 65 && prunedJetm[0] < 105";

    responseLeadPrunedJetMass_pr.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_responseLeadPrunedJetMass_pr").c_str(),
						    (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0.5,1.5));
    treeMonoV_Sig.at(iMass)->Draw(("(prunedJetm[0]/prunedJetGenm[0]) >> "+string(responseLeadPrunedJetMass_pr.back()->GetName())).c_str(),cut_pr.c_str(),"goff");    

    responseLeadPrunedJetMassRaw_pr.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_responseLeadPrunedJetMassRaw_pr").c_str(),
						       (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0.5,1.5));
    treeMonoV_Sig.at(iMass)->Draw(("(prunedJetmraw[0]/prunedJetGenm[0]) >> "+string(responseLeadPrunedJetMassRaw_pr.back()->GetName())).c_str(),cut_pr.c_str(),"goff");    

    subjetPtRatio_pr.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_subjetPtRatio_PR").c_str(),
					(bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),30,0,2));
    treeMonoV_Sig.at(iMass)->Draw(("(prunedSubJetpt_2[0]/prunedSubJetpt_1[0]) >> "+string(subjetPtRatio_pr.back()->GetName())).c_str(),cut_pr.c_str(),"goff");

    leadingJetTau2Tau1_pr.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingJetTau2Tau1_pr").c_str(),
					     (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),20,0,1));
    treeMonoV_Sig.at(iMass)->Draw(("boostedJettau2[0]/boostedJettau1[0] >> "+string(leadingJetTau2Tau1_pr.back()->GetName())).c_str(),cut_pr.c_str(),"goff");    

    string cut_sb = cut + " && prunedJetm[0] > 40 && prunedJetm[0] <65";

    leadingJetTau2Tau1_sb.push_back(new TH1F((bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]+"_leadingJetTau2Tau1_sb").c_str(),
					     (bosonMode+"_"+interactionModel+"_Mphi-"+Mphi[iMass]+"_Mchi-"+Mchi[iMass]).c_str(),20,0,1));
    treeMonoV_Sig.at(iMass)->Draw(("boostedJettau2[0]/boostedJettau1[0] >> "+string(leadingJetTau2Tau1_sb.back()->GetName())).c_str(),cut_sb.c_str(),"goff");    


  }

  // background 
  TH1F* prunedMass_Znunu  = new TH1F("prunedMass_Znunu","Z #rightarrow #nu #nu",30,0,120);
  TH1F* tau2tau1_Znunu    = new TH1F("tau2tau1_Znunu","Z #rightarrow #nu #nu",20,0,1);

  string cut;
  for(size_t iSelection = 0; iSelection < selectionsSig.size(); iSelection++){
    if(selectionsSig.at(iSelection).first == "ak8pt")
      cut = selectionsSig.at(iSelection).second;
  }

  backgroundZnunu->Draw(("prunedJetm[0] >> "+string(prunedMass_Znunu->GetName())).c_str(),cut.c_str(),"goff");
  backgroundZnunu->Draw(("boostedJettau2[0]/boostedJettau1[0] >> "+string(tau2tau1_Znunu->GetName())).c_str(),cut.c_str(),"goff");

  // make single plots
  makeShapePlots(cCanvas,leadingJetPt,legend,outputDirectory,"jetPt","p_{T} (GeV)");
  makeShapePlots(cCanvas,leadingJetEta,legend,outputDirectory,"jetEta","#eta");
  makeShapePlots(cCanvas,leadingJetTau2tau1,legend,outputDirectory,"jetTau2tau1","#tau_{2}/#tau_{1}");
  makeShapePlots(cCanvas,leadingPrunedJetpt,legend,outputDirectory,"jetPrunedPt","p_{T}^{pr} (GeV)");
  makeShapePlots(cCanvas,leadingPrunedJetm,legend,outputDirectory,"jetPrunedM","m_{pr} (GeV)");
  makeShapePlots(cCanvas,leadingPrunedJetptraw,legend,outputDirectory,"jetPrunedPtraw","p_{T}^{pr} raw (GeV)");
  makeShapePlots(cCanvas,leadingPrunedJetmraw,legend,outputDirectory,"jetPrunedMraw","m_{pr} raw (GeV)");
  makeShapePlots(cCanvas,leadingSoftDropJetpt,legend,outputDirectory,"jetSoftDropPt","p_{T}^{sd} (GeV)");
  makeShapePlots(cCanvas,leadingSoftDropJetm,legend,outputDirectory,"jetSoftDropM","m_{sd} (GeV)");
  makeShapePlots(cCanvas,leadingSoftDropJetptraw,legend,outputDirectory,"jetSoftDropPtraw","p_{T}^{sd} raw (GeV)");
  makeShapePlots(cCanvas,leadingSoftDropJetmraw,legend,outputDirectory,"jetSoftDropMraw","m_{sd} raw (GeV)");
  makeShapePlots(cCanvas,boostedNJet,legend,outputDirectory,"NJet","N_{jet}^{AK8} (GeV)");
  makeShapePlots(cCanvas,dRGenJet,legend,outputDirectory,"dRGenJet","#DeltaR(reco,gen)");
  makeShapePlots(cCanvas,dRVboson,legend,outputDirectory,"dRVboson","#DeltaR(reco,V)");
  makeShapePlots(cCanvas,subjetPtRatio,legend,outputDirectory,"SubJetPtRatio","#Pt ratio");

  makeShapePlots(cCanvas,responseLeadJetPt,legend,outputDirectory,"ResponsePt","p_{T}^{reco}-p_{T}^{gen} (GeV)");
  makeShapePlots(cCanvas,responseLeadPrunedJetMass,legend,outputDirectory,"ResponsePrunedJetMass","m_{pr}^{reco}-m_{pr}^{gen} (GeV)");
  makeShapePlots(cCanvas,responseLeadPrunedJetMassRaw,legend,outputDirectory,"ResponsePrunedJetMassRaw","m_{pr}^{reco}-m_{pr}^{gen} raw (GeV)");
  makeShapePlots(cCanvas,responseLeadSoftDropJetMass,legend,outputDirectory,"ResponseSoftDropJetMass","m_{sd}^{reco}-m_{sd}^{gen} (GeV)");
  makeShapePlots(cCanvas,responseLeadSoftDropJetMassRaw,legend,outputDirectory,"ResponseSoftDropJetMassRaw","m_{sd}^{reco}-m_{sd}^{gen} raw (GeV)");

  makeShapePlots(cCanvas,responseLeadPrunedJetMass_pr,legend,outputDirectory,"ResponsePrunedJetMass_pr","m_{pr}^{reco}-m_{pr}^{gen} (GeV)");
  makeShapePlots(cCanvas,responseLeadPrunedJetMassRaw_pr,legend,outputDirectory,"ResponsePrunedJetMassRaw_pr","m_{pr}^{reco}-m_{pr}^{gen} raw (GeV)");
  makeShapePlots(cCanvas,subjetPtRatio_pr,legend,outputDirectory,"SubJetPtRatio_pr","p_{T} ratio");
  makeShapePlots(cCanvas,leadingJetTau2Tau1_pr,legend,outputDirectory,"jetTau2tau1_pr","#tau_{2}/#tau_{1}");
  makeShapePlots(cCanvas,leadingJetTau2Tau1_sb,legend,outputDirectory,"jetTau2tau1_sb","#tau_{2}/#tau_{1}");

  makeShapePlots(cCanvas,leadingPrunedJetm,prunedMass_Znunu,legend,outputDirectory,"jetPrunedM","m_{pr} (GeV)");
  makeShapePlots(cCanvas,leadingJetTau2tau1,tau2tau1_Znunu,legend,outputDirectory,"jetTau2tau1","#tau_{2}/#tau_{1}");

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
  histoBkg->SetFillColor(kGray);
  histoBkg->SetFillStyle(3001);
  histoBkg->GetXaxis()->SetTitle(xAxisTitle.c_str());
  histoBkg->GetYaxis()->SetRangeUser(0,maxYScale*1.5);
  histoBkg->Draw("hist");

  legend->AddEntry(histoBkg,histoBkg->GetTitle(),"l");

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
